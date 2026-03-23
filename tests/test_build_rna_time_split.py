import gzip
import pickle

import pandas as pd

from scripts import build_rna_time_split as split_script


def _dump_pickle(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as f:
        pickle.dump(data, f)


def test_build_time_split_filters_test_sequences_seen_in_train(tmp_path):
    bioassembly_dir = tmp_path / "bioassembly"
    indices_path = tmp_path / "indices.csv"

    _dump_pickle(
        bioassembly_dir / "1abc_2.pkl.gz",
        {"sequences": {"1": "AUGC"}},
    )
    _dump_pickle(
        bioassembly_dir / "2xyz_5.pkl.gz",
        {"sequences": {"7": "AUGC"}},
    )
    _dump_pickle(
        bioassembly_dir / "3def_1.pkl.gz",
        {"sequences": {"3": "GGAA"}},
    )

    pd.DataFrame(
        [
            {
                "type": "chain",
                "pdb_id": "1abc",
                "bioassembly_name": "1abc_2",
                "entity_1_id": "1",
                "chain_1_id": "A",
                "sub_mol_1_type": "rna",
                "release_date": "2020-01-01",
            },
            {
                "type": "chain",
                "pdb_id": "2xyz",
                "bioassembly_name": "2xyz_5",
                "entity_1_id": "7",
                "chain_1_id": "B",
                "sub_mol_1_type": "modified_rna",
                "release_date": "2022-01-01",
            },
            {
                "type": "chain",
                "pdb_id": "3def",
                "bioassembly_name": "3def_1",
                "entity_1_id": "3",
                "chain_1_id": "C",
                "sub_mol_1_type": "rna",
                "release_date": "2022-02-01",
            },
            {
                "type": "chain",
                "pdb_id": "4dna",
                "bioassembly_name": "4dna_1",
                "entity_1_id": "4",
                "chain_1_id": "D",
                "sub_mol_1_type": "dna",
                "release_date": "2022-03-01",
            },
        ]
    ).to_csv(indices_path, index=False)

    train_df, test_df, dropped_test_df = split_script.build_time_split(
        indices_path=indices_path,
        bioassembly_dir=bioassembly_dir,
        test_release_date_from="2021-09-30",
    )

    assert train_df["pdb_id"].tolist() == ["1abc"]
    assert test_df["pdb_id"].tolist() == ["3def"]
    assert dropped_test_df["pdb_id"].tolist() == ["2xyz"]


def test_write_pdb_id_list_outputs_unique_sorted_ids(tmp_path):
    output_path = tmp_path / "test_pdb_ids.txt"
    df = pd.DataFrame({"pdb_id": ["3def", "1abc", "3def"]})

    split_script.write_pdb_id_list(df, output_path)

    assert output_path.read_text() == "1abc\n3def\n"
