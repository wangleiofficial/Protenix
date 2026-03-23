import json
import sys
import types

import numpy as np

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

from protenix.data.constants import DENSE_ATOM, PROTEIN_CHAIN, RNA_CHAIN
from protenix.data.template.template_parser import TemplateHit
from protenix.data.template.template_featurizer import TemplateFeaturizer
from protenix.data.template.template_utils import (
    DistogramFeaturesConfig,
    TemplateFeatures,
    TemplateHitFeaturizer,
)


def _build_minimal_protein_mmcif(sequence: str) -> str:
    res_names = ["ALA"] * len(sequence)
    atom_names = ["N", "CA", "C", "O", "CB"]

    lines = [
        "data_1abc",
        "#",
        "loop_",
        "_entity_poly_seq.entity_id",
        "_entity_poly_seq.num",
        "_entity_poly_seq.mon_id",
    ]
    for i, res_name in enumerate(res_names, start=1):
        lines.append(f"1 {i} {res_name}")

    lines += [
        "#",
        "loop_",
        "_entity_poly.entity_id",
        "_entity_poly.type",
        "1 'polypeptide(L)'",
        "#",
        "loop_",
        "_struct_asym.id",
        "_struct_asym.entity_id",
        "A 1",
        "#",
        "loop_",
        "_exptl.entry_id",
        "_exptl.method",
        "1ABC 'X-RAY DIFFRACTION'",
        "#",
        "loop_",
        "_pdbx_audit_revision_history.revision_date",
        "2020-01-01",
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_asym_id",
        "_atom_site.pdbx_PDB_model_num",
    ]

    atom_id = 1
    for res_id, res_name in enumerate(res_names, start=1):
        base = float(res_id * 3)
        coords = {
            "N": (base, 0.0, 0.0, "N"),
            "CA": (base + 1.0, 0.2, 0.0, "C"),
            "C": (base + 2.0, 0.0, 0.0, "C"),
            "O": (base + 2.6, 0.5, 0.0, "O"),
            "CB": (base + 1.0, 1.2, 0.0, "C"),
        }
        for atom_name in atom_names:
            x, y, z, elem = coords[atom_name]
            lines.append(
                f"ATOM {atom_id} {elem} {atom_name} . {res_name} A 1 {res_id} ? "
                f"{x:.3f} {y:.3f} {z:.3f} 1.00 10.00 {res_id} A 1"
            )
            atom_id += 1

    lines += ["#", ""]
    return "\n".join(lines)


def _build_minimal_rna_mmcif(sequence: str) -> str:
    atom_names = DENSE_ATOM["A"]

    lines = [
        "data_1abc",
        "#",
        "loop_",
        "_entity_poly_seq.entity_id",
        "_entity_poly_seq.num",
        "_entity_poly_seq.mon_id",
    ]
    for i, _ in enumerate(sequence, start=1):
        lines.append(f"2 {i} A")

    lines += [
        "#",
        "loop_",
        "_entity_poly.entity_id",
        "_entity_poly.type",
        "2 'polyribonucleotide'",
        "#",
        "loop_",
        "_struct_asym.id",
        "_struct_asym.entity_id",
        "C 2",
        "#",
        "loop_",
        "_exptl.entry_id",
        "_exptl.method",
        "1ABC 'X-RAY DIFFRACTION'",
        "#",
        "loop_",
        "_pdbx_audit_revision_history.revision_date",
        "2020-01-01",
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.type_symbol",
        "_atom_site.label_atom_id",
        "_atom_site.label_alt_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.label_entity_id",
        "_atom_site.label_seq_id",
        "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv",
        "_atom_site.auth_seq_id",
        "_atom_site.auth_asym_id",
        "_atom_site.pdbx_PDB_model_num",
    ]

    atom_id = 1
    for res_id, _ in enumerate(sequence, start=1):
        base = float(res_id * 4)
        for atom_idx, atom_name in enumerate(atom_names):
            element = next(ch for ch in atom_name if ch.isalpha()).upper()
            x = base + atom_idx * 0.05
            y = (atom_idx % 4) * 0.3
            z = (atom_idx % 5) * 0.2
            lines.append(
                f"ATOM {atom_id} {element} {atom_name} . A C 2 {res_id} ? "
                f"{x:.3f} {y:.3f} {z:.3f} 1.00 10.00 {res_id} T 1"
            )
            atom_id += 1

    lines += ["#", ""]
    return "\n".join(lines)


def test_protein_template_regression_keeps_legacy_feature_path(tmp_path):
    sequence = "AAAAAAAAAA"
    mmcif_dir = tmp_path / "mmcif"
    mmcif_dir.mkdir()
    (mmcif_dir / "1abc.cif").write_text(_build_minimal_protein_mmcif(sequence))

    release_dates_path = tmp_path / "release_dates.json"
    release_dates_path.write_text(json.dumps({"1abc": {"release_date": "2020-01-01"}}))

    # Keep the hit sequence non-identical so the protein duplicate prefilter does not
    # reject it; the actual template sequence still comes from the mmCIF chain.
    hit = TemplateHit(
        index=1,
        name="1abc_A",
        aligned_cols=len(sequence),
        sum_probs=99.0,
        query=sequence,
        hit_sequence="AAAAAAAAAB",
        indices_query=list(range(len(sequence))),
        indices_hit=list(range(len(sequence))),
    )

    featurizer = TemplateHitFeaturizer(
        mmcif_dir=str(mmcif_dir),
        kalign_binary_path="kalign",
        max_hits=1,
        max_template_date="2021-01-01",
        release_dates_path=str(release_dates_path),
    )

    result, _ = featurizer.get_templates(
        sequence_uid="protein_regression",
        query_sequence=sequence,
        hits=[hit],
        query_chain_type=PROTEIN_CHAIN,
    )

    assert result.errors == []
    assert result.warnings == []
    assert len(result.features) == 1

    packed = TemplateFeatures.package_template_features(hit_features=result.features)
    fixed = TemplateFeatures.fix_template_features(packed, num_res=len(sequence))

    assert fixed["template_aatype"].shape == (1, len(sequence))
    assert fixed["template_atom_positions"].shape == (1, len(sequence), 24, 3)
    assert fixed["template_atom_mask"].shape == (1, len(sequence), 24)
    assert fixed["template_atom_mask"].sum() > 0

    aatype = fixed["template_aatype"][0]
    atom_positions = fixed["template_atom_positions"][0]
    atom_mask = fixed["template_atom_mask"][0].astype(np.float32)

    pseudo_beta, _ = TemplateFeatures.pseudo_beta_fn(
        aatype=aatype,
        dense_atom_positions=atom_positions,
        dense_atom_masks=atom_mask,
    )
    dgram = TemplateFeatures.dgram_from_positions(
        pseudo_beta, DistogramFeaturesConfig()
    )
    unit_vector, backbone_frame_mask = TemplateFeatures.compute_template_unit_vector(
        atom_positions=atom_positions,
        atom_mask=atom_mask.astype(bool),
        aatype=aatype,
    )

    assert dgram.shape == (len(sequence), len(sequence), 39)
    assert float(dgram.sum()) > 0.0
    assert unit_vector.shape == (len(sequence), len(sequence), 3)
    assert float(np.abs(unit_vector).sum()) > 0.0
    assert float(backbone_frame_mask.sum()) > 0.0


def test_training_rna_template_supports_sequence_uid_and_direct_csv_mapping(tmp_path):
    sequence = "AAAAAAAAAA"
    mmcif_dir = tmp_path / "rna_mmcif"
    mmcif_dir.mkdir()
    (mmcif_dir / "1abc.cif").write_text(_build_minimal_rna_mmcif(sequence))

    release_dates_path = tmp_path / "release_dates.json"
    release_dates_path.write_text(json.dumps({"1abc": {"release_date": "2020-01-01"}}))
    obsolete_path = tmp_path / "obsolete.json"
    obsolete_path.write_text("{}")

    csv_path = tmp_path / "rna_templates.csv"
    csv_path.write_text(
        "\n".join(
            [
                "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
                f"rna_chain_sample,1abc_1_2_C,1.0,{len(sequence)},0,0,1,{len(sequence)},1,{len(sequence)},0.0,100",
            ]
        )
    )
    mapping_path = tmp_path / "rna_mapping.json"
    mapping_path.write_text(json.dumps({"rna_chain_sample": str(csv_path)}))

    featurizer = TemplateFeaturizer(
        dataset_name="WeightedPDB",
        enable_prot_template=False,
        enable_rna_template=True,
        prot_template_mmcif_dir=str(mmcif_dir),
        rna_template_mmcif_dir=str(mmcif_dir),
        rna_template_raw_paths=[""],
        rna_seq_or_filename_to_templatedir_jsons=[str(mapping_path)],
        rna_indexing_methods=["sequence_uid"],
        kalign_binary_path="kalign",
        release_dates_path=str(release_dates_path),
        obsolete_pdbs_path=str(obsolete_path),
    )

    templates, raw_hit_count, _ = featurizer.get_template(
        pdb_id="query",
        query_sequence=sequence,
        sequence_uid="rna_chain_sample",
        query_release_date=None,
        chain_entity_type=RNA_CHAIN,
    )

    assert raw_hit_count == 1
    assert len(templates) == 1

    packed = TemplateFeatures.package_template_features(hit_features=templates)
    fixed = TemplateFeatures.fix_template_features(packed, num_res=len(sequence))

    assert fixed["template_aatype"].shape == (1, len(sequence))
    assert fixed["template_atom_positions"].shape == (1, len(sequence), 24, 3)
    assert fixed["template_atom_mask"].shape == (1, len(sequence), 24)
    assert fixed["template_atom_mask"].sum() > 0

    aatype = fixed["template_aatype"][0]
    atom_positions = fixed["template_atom_positions"][0]
    atom_mask = fixed["template_atom_mask"][0].astype(np.float32)

    pseudo_beta, _ = TemplateFeatures.pseudo_beta_fn(
        aatype=aatype,
        dense_atom_positions=atom_positions,
        dense_atom_masks=atom_mask,
    )
    dgram = TemplateFeatures.dgram_from_positions(
        pseudo_beta, DistogramFeaturesConfig()
    )
    unit_vector, backbone_frame_mask = TemplateFeatures.compute_template_unit_vector(
        atom_positions=atom_positions,
        atom_mask=atom_mask.astype(bool),
        aatype=aatype,
    )

    assert dgram.shape == (len(sequence), len(sequence), 39)
    assert float(dgram.sum()) > 0.0
    assert unit_vector.shape == (len(sequence), len(sequence), 3)
    assert float(np.abs(unit_vector).sum()) > 0.0
    assert float(backbone_frame_mask.sum()) > 0.0
