import sys
import types

import numpy as np
import pytest
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
from protenix.data.constraint.constraint_featurizer import ConstraintFeatureGenerator
from protenix.data.tokenizer import AtomArrayTokenizer


def _build_minimal_rna_atom_array(sequence: str) -> AtomArray:
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
    label_entity_id = []
    copy_id = []
    is_ligand = []
    is_protein = []
    is_rna = []
    is_dna = []

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
            asym_id_int.append(0)
            label_entity_id.append("1")
            copy_id.append(1)
            is_ligand.append(False)
            is_protein.append(False)
            is_rna.append(True)
            is_dna.append(False)

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
    atom_array.set_annotation("label_entity_id", np.asarray(label_entity_id))
    atom_array.set_annotation("copy_id", np.asarray(copy_id, dtype=int))
    atom_array.set_annotation("is_ligand", np.asarray(is_ligand, dtype=bool))
    atom_array.set_annotation("is_protein", np.asarray(is_protein, dtype=bool))
    atom_array.set_annotation("is_rna", np.asarray(is_rna, dtype=bool))
    atom_array.set_annotation("is_dna", np.asarray(is_dna, dtype=bool))
    return atom_array


def _build_minimal_protein_atom_array(sequence: str) -> AtomArray:
    atom_names = ["N", "CA", "C", "O", "CB"]
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
    label_entity_id = []
    copy_id = []
    is_ligand = []
    is_protein = []
    is_rna = []
    is_dna = []

    for res_index, _ in enumerate(sequence, start=1):
        base = float(res_index * 3)
        atom_coords = {
            "N": (base, 0.0, 0.0, "N"),
            "CA": (base + 1.0, 0.2, 0.0, "C"),
            "C": (base + 2.0, 0.0, 0.0, "C"),
            "O": (base + 2.6, 0.5, 0.0, "O"),
            "CB": (base + 1.0, 1.2, 0.0, "C"),
        }
        for atom_name in atom_names:
            x, y, z, elem = atom_coords[atom_name]
            coords.append([x, y, z])
            chain_ids.append("A")
            res_ids.append(res_index)
            res_names.append("ALA")
            atom_name_ann.append(atom_name)
            elements.append(elem)
            hetero.append(False)
            centre_atom_mask.append(atom_name == "CA")
            mol_type.append("protein")
            asym_id_int.append(0)
            label_entity_id.append("1")
            copy_id.append(1)
            is_ligand.append(False)
            is_protein.append(True)
            is_rna.append(False)
            is_dna.append(False)

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
    atom_array.set_annotation("label_entity_id", np.asarray(label_entity_id))
    atom_array.set_annotation("copy_id", np.asarray(copy_id, dtype=int))
    atom_array.set_annotation("is_ligand", np.asarray(is_ligand, dtype=bool))
    atom_array.set_annotation("is_protein", np.asarray(is_protein, dtype=bool))
    atom_array.set_annotation("is_rna", np.asarray(is_rna, dtype=bool))
    atom_array.set_annotation("is_dna", np.asarray(is_dna, dtype=bool))
    return atom_array


def test_inference_constraint_allows_same_chain_rna_contact():
    atom_array = _build_minimal_rna_atom_array("AAAA")
    token_array = AtomArrayTokenizer(atom_array).get_token_array()

    feature_dict, _, _ = ConstraintFeatureGenerator.generate_from_json(
        token_array=token_array,
        atom_array=atom_array,
        sequences=[{"rnaSequence": {"sequence": "AAAA", "count": 1}}],
        constraint_param={
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
        },
    )

    assert tuple(feature_dict["contact"].shape) == (4, 4, 2)
    assert float(feature_dict["contact"][0, 3, 1]) == 10.0
    assert float(feature_dict["contact"][3, 0, 1]) == 10.0


def test_inference_constraint_keeps_same_chain_protein_rejected():
    atom_array = _build_minimal_protein_atom_array("AAAA")
    token_array = AtomArrayTokenizer(atom_array).get_token_array()

    with pytest.raises(ValueError, match="same chain"):
        ConstraintFeatureGenerator.generate_from_json(
            token_array=token_array,
            atom_array=atom_array,
            sequences=[{"proteinChain": {"sequence": "AAAA", "count": 1}}],
            constraint_param={
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
            },
        )
