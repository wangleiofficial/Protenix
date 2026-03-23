import json
from pathlib import Path

from scripts.split_rna_template_hits import split_template_hits


def test_split_template_hits_with_query_rename_and_relative_mapping(tmp_path):
    input_csv = tmp_path / "combined.csv"
    input_csv.write_text(
        "\n".join(
            [
                "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
                "9v44_1_6_F,9v44_1_6_F,1.0,24,0,0,1,24,1,24,7.832e-08,45",
                "9v44_1_6_F,9v48_1_3_C,1.0,24,0,0,1,24,1,24,7.832e-08,45",
                "9wnq_1_32_GA,4v4q_1_22_V,1.0,120,0,0,1,120,1,120,2.928e-62,216",
            ]
        )
        + "\n"
    )

    query_map = tmp_path / "query_map.json"
    query_map.write_text(
        json.dumps(
            {
                "9v44_1_6_F": "sample_a",
                "9wnq_1_32_GA": "sample_b",
            }
        )
    )

    output_dir = tmp_path / "split"
    mapping_out = tmp_path / "mapping.json"

    mapping = split_template_hits(
        input_path=input_csv,
        output_dir=output_dir,
        mapping_out=mapping_out,
        query_map=json.loads(query_map.read_text()),
        output_format="csv",
        layout="flat",
        path_mode="relative",
    )

    assert mapping == {
        "sample_a": "sample_a.csv",
        "sample_b": "sample_b.csv",
    }

    saved_mapping = json.loads(mapping_out.read_text())
    assert saved_mapping == mapping

    sample_a = output_dir / "sample_a.csv"
    sample_b = output_dir / "sample_b.csv"
    assert sample_a.exists()
    assert sample_b.exists()

    sample_a_lines = sample_a.read_text().strip().splitlines()
    assert len(sample_a_lines) == 3
    assert sample_a_lines[0].startswith("query,target,fident")
    assert sample_a_lines[1].startswith("9v44_1_6_F,9v44_1_6_F")
    assert sample_a_lines[2].startswith("9v44_1_6_F,9v48_1_3_C")

    sample_b_lines = sample_b.read_text().strip().splitlines()
    assert len(sample_b_lines) == 2
    assert sample_b_lines[1].startswith("9wnq_1_32_GA,4v4q_1_22_V")


def test_split_template_hits_dir_layout_for_m8(tmp_path):
    input_m8 = tmp_path / "combined.m8"
    input_m8.write_text(
        "\n".join(
            [
                "q1 t1 1.0 24 0 0 1 24 1 24 1e-5 40",
                "q2 t2 1.0 30 0 0 1 30 1 30 1e-6 50",
            ]
        )
        + "\n"
    )

    output_dir = tmp_path / "split_dirs"
    mapping_out = tmp_path / "mapping.json"

    mapping = split_template_hits(
        input_path=input_m8,
        output_dir=output_dir,
        mapping_out=mapping_out,
        input_format="m8",
        output_format="same",
        layout="dir",
        path_mode="relative",
    )

    assert mapping == {
        "q1": "q1",
        "q2": "q2",
    }
    assert (output_dir / "q1" / "result.m8").exists()
    assert (output_dir / "q2" / "result.m8").exists()


def test_split_template_hits_can_drop_self_hits_and_same_pdb(tmp_path):
    input_csv = tmp_path / "combined.csv"
    input_csv.write_text(
        "\n".join(
            [
                "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
                "9v44_1_6_F,9v44_1_6_F,1.0,24,0,0,1,24,1,24,7.832e-08,45",
                "9v44_1_6_F,9v48_1_3_C,1.0,24,0,0,1,24,1,24,7.832e-08,45",
                "9wnq_1_32_GA,9wnq_2_10_A,0.9,110,1,0,1,110,2,111,3.2e-50,180",
                "9wnq_1_32_GA,4v4q_1_22_V,1.0,120,0,0,1,120,1,120,2.928e-62,216",
            ]
        )
        + "\n"
    )

    output_dir = tmp_path / "split"
    mapping_out = tmp_path / "mapping.json"

    mapping = split_template_hits(
        input_path=input_csv,
        output_dir=output_dir,
        mapping_out=mapping_out,
        output_format="csv",
        layout="flat",
        path_mode="relative",
        drop_self_hit=True,
        drop_same_pdb=True,
    )

    assert mapping == {
        "9v44_1_6_F": "9v44_1_6_F.csv",
        "9wnq_1_32_GA": "9wnq_1_32_GA.csv",
    }

    sample_a_lines = (output_dir / "9v44_1_6_F.csv").read_text().strip().splitlines()
    assert len(sample_a_lines) == 2
    assert sample_a_lines[1].startswith("9v44_1_6_F,9v48_1_3_C")

    sample_b_lines = (output_dir / "9wnq_1_32_GA.csv").read_text().strip().splitlines()
    assert len(sample_b_lines) == 2
    assert sample_b_lines[1].startswith("9wnq_1_32_GA,4v4q_1_22_V")
