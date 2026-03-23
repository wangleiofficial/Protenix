from scripts.database import download_rna_pdb_dataset as dataset_script


def test_build_rna_entry_query_includes_expected_filters():
    query = dataset_script.build_rna_entry_query(
        max_resolution=4.5,
        release_date_before="2021-09-30",
        limit=25,
        start=50,
        rows=100,
    )

    assert query["return_type"] == "entry"
    assert query["request_options"]["paginate"] == {"start": 50, "rows": 25}

    nodes = query["query"]["nodes"]
    assert nodes[0]["parameters"]["attribute"] == "rcsb_entry_info.polymer_entity_count_RNA"
    assert nodes[0]["parameters"]["operator"] == "greater"
    assert nodes[0]["parameters"]["value"] == 0
    assert nodes[1]["parameters"]["attribute"] == "rcsb_entry_info.resolution_combined"
    assert nodes[1]["parameters"]["value"] == 4.5
    assert (
        nodes[2]["parameters"]["attribute"]
        == "rcsb_accession_info.initial_release_date"
    )
    assert nodes[2]["parameters"]["value"] == "2021-09-30"


def test_search_rna_entry_ids_handles_pagination(monkeypatch):
    responses = [
        {
            "total_count": 3,
            "result_set": [
                {"identifier": "1ABC"},
                {"identifier": "2XYZ"},
            ],
        },
        {
            "total_count": 3,
            "result_set": [
                {"identifier": "3DEF"},
            ],
        },
    ]

    def fake_post(payload, *, timeout=60, user_agent=dataset_script.DEFAULT_USER_AGENT):
        del timeout, user_agent
        start = payload["request_options"]["paginate"]["start"]
        if start == 0:
            return responses[0]
        if start == 2:
            return responses[1]
        raise AssertionError(f"Unexpected start offset: {start}")

    monkeypatch.setattr(dataset_script, "_post_search_query", fake_post)

    entry_ids = dataset_script.search_rna_entry_ids(
        max_resolution=4.5,
        limit=None,
        rows_per_page=2,
    )

    assert entry_ids == ["1abc", "2xyz", "3def"]


def test_write_lines_writes_one_entry_per_line(tmp_path):
    output_path = tmp_path / "pdb_ids.txt"

    dataset_script.write_lines(["1abc", "2xyz"], output_path)

    assert output_path.read_text() == "1abc\n2xyz\n"
