import json
import sys
import types

import numpy as np
import pandas as pd
from biotite.structure import AtomArray

try:
    from rdkit.Chem import GetPeriodicTable  # noqa: F401
except ImportError:
    rdkit = types.ModuleType("rdkit")
    chem = types.ModuleType("rdkit.Chem")

    class _PeriodicTable:
        def GetElementSymbol(self, idx):
            return "C"

    chem.GetPeriodicTable = lambda: _PeriodicTable()
    chem.Mol = type("Mol", (), {})
    rdkit.Chem = chem
    sys.modules["rdkit"] = rdkit
    sys.modules["rdkit.Chem"] = chem

from protenix.data.constants import DENSE_ATOM
from protenix.data.constraint.constraint_featurizer import (
    ConstraintFeatureGenerator,
    ConstraintSourceManager,
)
from protenix.data.pipeline.dataset import BaseSingleDataset
from protenix.data.tokenizer import AtomArrayTokenizer


def _build_minimal_rna_atom_array(sequence: str, asym_id: int = 2) -> AtomArray:
    atom_names = DENSE_ATOM["A"]
    atoms_per_res = len(atom_names)
    atom_array = AtomArray(len(sequence) * atoms_per_res)

    coords = []
    chain_ids = []
    res_ids = []
    res_names = []
    atom_name_ann = []
    elements = []
    hetero = []
    centre_atom_mask = []
    mol_type = []
    asym_id_int = []
    is_rna = []
    is_dna = []
    is_ligand = []
    is_protein = []
    is_resolved = []
    cano_seq_resname = []

    centre_atom_name = "C4" if "C4" in atom_names else atom_names[0]

    for res_index, _ in enumerate(sequence, start=1):
        base = float(res_index * 4)
        for atom_index, atom_name in enumerate(atom_names):
            coords.append(
                [
                    base + atom_index * 0.1,
                    float(atom_index % 4) * 0.25,
                    float(atom_index % 5) * 0.2,
                ]
            )
            chain_ids.append("A")
            res_ids.append(res_index)
            res_names.append("A")
            atom_name_ann.append(atom_name)
            elements.append(next(ch for ch in atom_name if ch.isalpha()).upper())
            hetero.append(False)
            centre_atom_mask.append(atom_name == centre_atom_name)
            mol_type.append("rna")
            asym_id_int.append(asym_id)
            is_rna.append(True)
            is_dna.append(False)
            is_ligand.append(False)
            is_protein.append(False)
            is_resolved.append(True)
            cano_seq_resname.append("A")

    atom_array.coord = np.asarray(coords, dtype=np.float32)
    atom_array.chain_id = np.asarray(chain_ids)
    atom_array.res_id = np.asarray(res_ids)
    atom_array.res_name = np.asarray(res_names)
    atom_array.atom_name = np.asarray(atom_name_ann)
    atom_array.element = np.asarray(elements)
    atom_array.hetero = np.asarray(hetero)
    atom_array.set_annotation("centre_atom_mask", np.asarray(centre_atom_mask, dtype=int))
    atom_array.set_annotation("mol_type", np.asarray(mol_type, dtype=object))
    atom_array.set_annotation("asym_id_int", np.asarray(asym_id_int, dtype=int))
    atom_array.set_annotation("is_rna", np.asarray(is_rna, dtype=bool))
    atom_array.set_annotation("is_dna", np.asarray(is_dna, dtype=bool))
    atom_array.set_annotation("is_ligand", np.asarray(is_ligand, dtype=bool))
    atom_array.set_annotation("is_protein", np.asarray(is_protein, dtype=bool))
    atom_array.set_annotation("is_resolved", np.asarray(is_resolved, dtype=bool))
    atom_array.set_annotation("cano_seq_resname", np.asarray(cano_seq_resname))
    return atom_array


def test_constraint_source_manager_supports_sequence_uid_mapping(tmp_path):
    constraint_path = tmp_path / "1abc_2.json"
    constraint_path.write_text(json.dumps({"contact": []}))

    mapping_path = tmp_path / "mapping.json"
    mapping_path.write_text(json.dumps({"1abc_2": str(constraint_path)}))

    source_mgr = ConstraintSourceManager(
        raw_paths=[""],
        seq_or_filename_to_constraint_jsons=[str(mapping_path)],
        indexing_methods=["sequence_uid"],
    )

    resolved = source_mgr.fetch_constraint_path(
        query_sequence="AAAA",
        sequence_uid="1abc_2",
    )
    assert resolved == str(constraint_path)


def test_constraint_generator_supports_external_same_chain_rna_contacts():
    atom_array = _build_minimal_rna_atom_array("AAAA")
    token_array = AtomArrayTokenizer(atom_array).get_token_array()
    generator = ConstraintFeatureGenerator(
        constraint_config={
            "contact": {"prob": 0.0, "feature_type": "continuous"},
            "pocket": {"prob": 0.0},
            "contact_atom": {"prob": 0.0},
            "substructure": {"prob": 0.0},
        },
        ab_top2_clusters=set(),
    )

    (
        _,
        _,
        _,
        constraint_feature,
        _,
        _,
        _,
    ) = generator.generate(
        atom_array=atom_array,
        token_array=token_array,
        sample_indice=pd.Series({"pdb_id": "1abc"}),
        pdb_indice=pd.Series(dtype=object),
        msa_features={},
        max_entity_mol_id=None,
        full_atom_array=atom_array,
        external_rna_contact_payload={
            "constraint": {
                "contact": [
                    {
                        "entity1": 1,
                        "copy1": 1,
                        "position1": 1,
                        "entity2": 1,
                        "copy2": 1,
                        "position2": 4,
                        "max_distance": 10.0,
                        "min_distance": 0.0,
                    },
                    {
                        "entity1": 1,
                        "copy1": 1,
                        "position1": 2,
                        "entity2": 1,
                        "copy2": 1,
                        "position2": 3,
                        "max_distance": 8.0,
                        "min_distance": 0.0,
                    },
                ]
            }
        },
    )

    contact_feature = constraint_feature["contact"]
    assert tuple(contact_feature.shape) == (4, 4, 2)
    assert float(contact_feature[0, 3, 1]) == 10.0
    assert float(contact_feature[3, 0, 1]) == 10.0
    assert float(contact_feature[1, 2, 1]) == 8.0
    assert float(contact_feature[2, 1, 1]) == 8.0


def test_dataset_loader_uses_sequence_uid_for_external_rna_ss(tmp_path):
    constraint_path = tmp_path / "constraint.json"
    constraint_payload = {
        "contact": [
            {
                "entity1": 1,
                "copy1": 1,
                "position1": 1,
                "entity2": 1,
                "copy2": 1,
                "position2": 4,
                "max_distance": 10.0,
                "min_distance": 0.0,
            }
        ]
    }
    constraint_path.write_text(json.dumps(constraint_payload))

    mapping_path = tmp_path / "mapping.json"
    mapping_path.write_text(json.dumps({"1abc_2": str(constraint_path)}))

    dataset = BaseSingleDataset.__new__(BaseSingleDataset)
    dataset.constraint = {"rna_ss": {"strict": False}}
    dataset.rna_ss_source_mgr = ConstraintSourceManager(
        raw_paths=[""],
        seq_or_filename_to_constraint_jsons=[str(mapping_path)],
        indexing_methods=["sequence_uid"],
    )

    atom_array = _build_minimal_rna_atom_array("AAAA", asym_id=2)
    token_array = AtomArrayTokenizer(atom_array).get_token_array()

    loaded = dataset._load_external_rna_ss_constraint(
        sample_indice=pd.Series({"pdb_id": "1abc"}),
        token_array=token_array,
        atom_array=atom_array,
    )
    assert loaded == constraint_payload
