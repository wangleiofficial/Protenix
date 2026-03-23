import sys
import types

import numpy as np
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

from protenix.data.core.parser import MMCIFParser
from protenix.data.pipeline.dataset import BaseSingleDataset
from protenix.data.tokenizer import Token, TokenArray
from scripts import prepare_training_data as prep_script


def _chain_row(
    *,
    chain_id: str,
    mol_type: str,
    sub_mol_type: str,
) -> dict:
    return {
        "type": "chain",
        "chain_1_id": chain_id,
        "chain_2_id": "",
        "mol_1_type": mol_type,
        "mol_2_type": "",
        "sub_mol_1_type": sub_mol_type,
        "sub_mol_2_type": "",
    }


def _interface_row(
    *,
    chain_1_id: str,
    mol_1_type: str,
    sub_mol_1_type: str,
    chain_2_id: str,
    mol_2_type: str,
    sub_mol_2_type: str,
) -> dict:
    return {
        "type": "interface",
        "chain_1_id": chain_1_id,
        "chain_2_id": chain_2_id,
        "mol_1_type": mol_1_type,
        "mol_2_type": mol_2_type,
        "sub_mol_1_type": sub_mol_1_type,
        "sub_mol_2_type": sub_mol_2_type,
    }


def _make_single_atom_token(atom_index: int) -> Token:
    token = Token(1)
    token.atom_indices = [atom_index]
    token.atom_names = ["C1'"]
    token.centre_atom_index = atom_index
    return token


def _build_mock_rna_bioassembly_dict() -> dict:
    atom_array = AtomArray(4)
    atom_array.coord = np.zeros((4, 3), dtype=float)
    atom_array.chain_id = np.array(["A", "B", "C", "D"], dtype=object)
    atom_array.res_id = np.array([1, 1, 1, 1], dtype=int)
    atom_array.res_name = np.array(["A", "U", "DA", "A"], dtype=object)
    atom_array.atom_name = np.array(["C1'", "C1'", "C1'", "C1'"], dtype=object)
    atom_array.element = np.array(["C", "C", "C", "C"], dtype=object)
    atom_array.hetero = np.array([False, False, False, False], dtype=bool)
    atom_array.set_annotation("centre_atom_mask", np.array([1, 1, 1, 1], dtype=int))
    atom_array.set_annotation("asym_id_int", np.array([5, 9, 11, 13], dtype=int))
    atom_array.set_annotation(
        "label_entity_id", np.array(["1", "2", "3", "4"], dtype=object)
    )

    token_array = TokenArray(
        [
            _make_single_atom_token(0),
            _make_single_atom_token(1),
            _make_single_atom_token(2),
            _make_single_atom_token(3),
        ]
    )
    return {
        "pdb_id": "1abc",
        "assembly_id": "1",
        "release_date": "2024-01-01",
        "num_assembly_polymer_chains": 4,
        "num_prot_chains": 0,
        "entity_poly_type": {
            "1": "polyribonucleotide",
            "2": "polyribonucleotide",
            "3": "polydeoxyribonucleotide",
            "4": "polydeoxyribonucleotide/polyribonucleotide hybrid",
        },
        "sequences": {"1": "A", "2": "U", "3": "T", "4": "A"},
        "resolution": 2.0,
        "num_tokens": 4,
        "atom_array": atom_array,
        "token_array": token_array,
        "msa_features": None,
        "template_features": None,
    }


def test_rna_monomer_collects_rna_chain_rows_from_complexes():
    rows = [
        _chain_row(chain_id="A", mol_type="nuc", sub_mol_type="rna"),
        _chain_row(chain_id="B", mol_type="nuc", sub_mol_type="modified_rna"),
        _chain_row(chain_id="C", mol_type="nuc", sub_mol_type="dna"),
        _chain_row(chain_id="D", mol_type="nuc", sub_mol_type="dna_rna_hybrid"),
        _chain_row(chain_id="P", mol_type="prot", sub_mol_type="prot"),
        _interface_row(
            chain_1_id="A",
            mol_1_type="nuc",
            sub_mol_1_type="rna",
            chain_2_id="P",
            mol_2_type="prot",
            sub_mol_2_type="prot",
        ),
    ]

    filtered = prep_script.filter_sample_indices_by_rna_mode(rows, rna_mode="rna_monomer")

    assert filtered == [
        _chain_row(chain_id="A", mol_type="nuc", sub_mol_type="rna"),
        _chain_row(chain_id="B", mol_type="nuc", sub_mol_type="modified_rna"),
    ]


def test_build_rna_monomer_samples_splits_each_rna_chain_and_preserves_original_pdb_id():
    rows = [
        _chain_row(chain_id="A", mol_type="nuc", sub_mol_type="rna"),
        _chain_row(chain_id="B", mol_type="nuc", sub_mol_type="modified_rna"),
        _chain_row(chain_id="C", mol_type="nuc", sub_mol_type="dna"),
        _chain_row(chain_id="D", mol_type="nuc", sub_mol_type="dna_rna_hybrid"),
    ]

    samples = prep_script.build_rna_monomer_samples(
        bioassembly_dict=_build_mock_rna_bioassembly_dict(),
        sample_indices_list=rows,
    )

    assert [sample[0] for sample in samples] == ["1abc_5", "1abc_9"]
    assert [sample[1]["bioassembly_name"] for sample in samples] == ["1abc_5", "1abc_9"]
    assert [sample[1]["pdb_id"] for sample in samples] == ["1abc", "1abc"]
    assert [sample[2]["num_tokens"] for sample in samples] == [1, 1]
    assert [sample[2]["num_prot_chains"] for sample in samples] == [0, 0]
    assert [sample[2]["atom_array"].asym_id_int.tolist() for sample in samples] == [[5], [9]]


def test_rna_only_keeps_rna_chains_and_rna_rna_interfaces():
    rows = [
        _chain_row(chain_id="A", mol_type="nuc", sub_mol_type="rna"),
        _chain_row(chain_id="B", mol_type="nuc", sub_mol_type="modified_rna"),
        _interface_row(
            chain_1_id="A",
            mol_1_type="nuc",
            sub_mol_1_type="rna",
            chain_2_id="B",
            mol_2_type="nuc",
            sub_mol_2_type="modified_rna",
        ),
    ]

    filtered = prep_script.filter_sample_indices_by_rna_mode(rows, rna_mode="rna_only")

    assert filtered == rows


def test_rna_only_excludes_rna_protein_complexes():
    rows = [
        _chain_row(chain_id="A", mol_type="nuc", sub_mol_type="rna"),
        _chain_row(chain_id="P", mol_type="prot", sub_mol_type="prot"),
        _interface_row(
            chain_1_id="A",
            mol_1_type="nuc",
            sub_mol_1_type="rna",
            chain_2_id="P",
            mol_2_type="prot",
            sub_mol_2_type="prot",
        ),
    ]

    filtered = prep_script.filter_sample_indices_by_rna_mode(rows, rna_mode="rna_only")

    assert filtered == []


def test_rna_all_keeps_rna_and_rna_protein_interfaces_but_drops_protein_chain_rows():
    rna_chain = _chain_row(chain_id="A", mol_type="nuc", sub_mol_type="rna")
    protein_chain = _chain_row(chain_id="P", mol_type="prot", sub_mol_type="prot")
    rna_protein = _interface_row(
        chain_1_id="A",
        mol_1_type="nuc",
        sub_mol_1_type="rna",
        chain_2_id="P",
        mol_2_type="prot",
        sub_mol_2_type="prot",
    )
    rows = [rna_chain, protein_chain, rna_protein]

    filtered = prep_script.filter_sample_indices_by_rna_mode(rows, rna_mode="rna_all")

    assert filtered == [rna_chain, rna_protein]


def test_rna_modes_exclude_dna_or_hybrid_entries():
    rows = [
        _chain_row(chain_id="A", mol_type="nuc", sub_mol_type="dna"),
        _chain_row(chain_id="B", mol_type="nuc", sub_mol_type="rna"),
    ]

    filtered = prep_script.filter_sample_indices_by_rna_mode(rows, rna_mode="rna_all")

    assert filtered == []


def test_make_chain_indices_classifies_dna_and_hybrid_as_nuc(monkeypatch):
    atom_array = AtomArray(3)
    atom_array.coord = np.zeros((3, 3), dtype=float)
    atom_array.chain_id = np.array(["A", "B", "C"], dtype=object)
    atom_array.label_entity_id = np.array(["1", "2", "3"], dtype=object)
    atom_array.res_id = np.array([1, 1, 1], dtype=int)
    atom_array.res_name = np.array(["A", "DA", "A"], dtype=object)
    atom_array.atom_name = np.array(["C1'", "C1'", "C1'"], dtype=object)
    atom_array.element = np.array(["C", "C", "C"], dtype=object)
    atom_array.hetero = np.array([False, False, False], dtype=bool)
    atom_array.set_annotation("is_resolved", np.array([True, True, True], dtype=bool))
    atom_array.set_annotation("centre_atom_mask", np.array([1, 1, 1], dtype=int))

    class _DummyParser:
        pdb_id = "1abc"
        entity_poly_type = {
            "1": "polyribonucleotide",
            "2": "polydeoxyribonucleotide",
            "3": "polydeoxyribonucleotide/polyribonucleotide hybrid",
        }

        @staticmethod
        def get_poly_res_names(_atom_array):
            return {"1": ["A"], "2": ["DA"], "3": ["A"]}

    monkeypatch.setattr(
        "protenix.data.core.parser.ccd.res_names_to_sequence",
        lambda res_names: "".join(res_names),
    )

    chain_indices = MMCIFParser.make_chain_indices(_DummyParser(), atom_array)

    assert [row["mol_type"] for row in chain_indices] == ["nuc", "nuc", "nuc"]


def test_bioassembly_pickle_stem_prefers_bioassembly_name():
    class _SampleIndice:
        pdb_id = "1abc"
        bioassembly_name = "1abc_9"

    assert BaseSingleDataset.get_bioassembly_pickle_stem(_SampleIndice()) == "1abc_9"


def test_bioassembly_pickle_stem_falls_back_to_pdb_id():
    class _SampleIndice:
        pdb_id = "1abc"
        bioassembly_name = ""

    assert BaseSingleDataset.get_bioassembly_pickle_stem(_SampleIndice()) == "1abc"
