import importlib


def test_rna_config_autoload_prefers_standard_mapping_names(tmp_path, monkeypatch):
    root = tmp_path
    (root / "rna_msa").mkdir()
    (root / "rna_msa" / "mapping_sequence.json").write_text("{}\n")
    (root / "rna_template_hits").mkdir()
    (root / "rna_template_hits" / "mapping_sequence_uid.json").write_text("{}\n")
    (root / "rna_ss_constraints").mkdir()
    (root / "rna_ss_constraints" / "mapping_sequence_uid.json").write_text("{}\n")

    monkeypatch.setenv("PROTENIX_ROOT_DIR", str(root))

    import configs.configs_data as configs_data_module

    configs_data_module = importlib.reload(configs_data_module)

    assert (
        configs_data_module.data_configs["rna_monomer_train"]["base_info"]["indices_fpath"]
        == str(root / "indices" / "rna_monomer_train.csv")
    )
    assert configs_data_module.data_configs["msa"][
        "rna_seq_or_filename_to_msadir_jsons"
    ].value == [str(root / "rna_msa" / "mapping_sequence.json")]
    assert configs_data_module.data_configs["msa"]["rna_msadir_raw_paths"].value == [
        str(root / "rna_msa" / "msas")
    ]
    assert configs_data_module.data_configs["template"][
        "rna_seq_or_filename_to_templatedir_jsons"
    ].value == [str(root / "rna_template_hits" / "mapping_sequence_uid.json")]
    assert configs_data_module.default_weighted_pdb_configs["constraint"]["rna_ss"][
        "seq_or_filename_to_ss_jsons"
    ].value == [str(root / "rna_ss_constraints" / "mapping_sequence_uid.json")]


def test_rna_config_autoload_falls_back_to_legacy_msa_mapping_name(tmp_path, monkeypatch):
    root = tmp_path
    (root / "rna_msa").mkdir()
    (root / "rna_msa" / "rna_sequence_to_pdb_chains.json").write_text("{}\n")

    monkeypatch.setenv("PROTENIX_ROOT_DIR", str(root))

    import configs.configs_data as configs_data_module

    configs_data_module = importlib.reload(configs_data_module)

    assert configs_data_module.data_configs["msa"][
        "rna_seq_or_filename_to_msadir_jsons"
    ].value == [str(root / "rna_msa" / "rna_sequence_to_pdb_chains.json")]
