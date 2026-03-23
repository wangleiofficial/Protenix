# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import csv
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from joblib import delayed, Parallel

from tqdm import tqdm

RNA_MODE_CHOICES = ("all", "rna_monomer", "rna_only", "rna_all")


def _is_rna_sub_mol_type(sub_mol_type: str) -> bool:
    """Returns True when a sub_mol_type denotes RNA but not DNA/RNA hybrid."""
    subtype = str(sub_mol_type or "")
    return subtype != "dna_rna_hybrid" and subtype.split("_")[-1] == "rna"


def _is_rna_chain_row(row: dict) -> bool:
    """Returns True for chain rows representing RNA."""
    return (
        row.get("type") == "chain"
        and row.get("mol_1_type") == "nuc"
        and _is_rna_sub_mol_type(row.get("sub_mol_1_type", ""))
    )


def _is_rna_rna_interface_row(row: dict) -> bool:
    """Returns True for interfaces where both sides are RNA."""
    return (
        row.get("type") == "interface"
        and row.get("mol_1_type") == "nuc"
        and row.get("mol_2_type") == "nuc"
        and _is_rna_sub_mol_type(row.get("sub_mol_1_type", ""))
        and _is_rna_sub_mol_type(row.get("sub_mol_2_type", ""))
    )


def _is_rna_protein_interface_row(row: dict) -> bool:
    """Returns True for interfaces between RNA and protein."""
    if row.get("type") != "interface":
        return False

    left_is_rna = row.get("mol_1_type") == "nuc" and _is_rna_sub_mol_type(
        row.get("sub_mol_1_type", "")
    )
    right_is_rna = row.get("mol_2_type") == "nuc" and _is_rna_sub_mol_type(
        row.get("sub_mol_2_type", "")
    )
    left_is_protein = row.get("mol_1_type") == "prot"
    right_is_protein = row.get("mol_2_type") == "prot"
    return (left_is_rna and right_is_protein) or (right_is_rna and left_is_protein)


def make_bioassembly_name(pdb_id: str, asym_id: int) -> str:
    """Builds the on-disk name for a split single-chain bioassembly."""
    return f"{pdb_id}_{int(asym_id)}"


def build_rna_monomer_samples(
    bioassembly_dict: dict,
    sample_indices_list: list[dict],
) -> list[tuple[str, dict, dict]]:
    """
    Splits every RNA chain row into an independent single-chain bioassembly sample.

    Returns:
        list of (bioassembly_name, sample_row, single_chain_bioassembly_dict)
    """
    rna_chain_rows = [row for row in sample_indices_list if _is_rna_chain_row(row)]
    if not rna_chain_rows:
        return []

    import torch

    from protenix.utils.cropping import CropData

    atom_array = bioassembly_dict["atom_array"]
    token_array = bioassembly_dict["token_array"]
    centre_atom_indices = np.asarray(token_array.get_annotation("centre_atom_index"))
    token_chain_ids = atom_array[centre_atom_indices].chain_id

    samples = []
    for row in rna_chain_rows:
        chain_id = row["chain_1_id"]
        selected_token_indices = np.where(token_chain_ids == chain_id)[0]
        if selected_token_indices.size == 0:
            continue

        cropped_token_array, cropped_atom_array = CropData.select_by_token_indices(
            token_array=token_array,
            atom_array=atom_array,
            selected_token_indices=torch.as_tensor(selected_token_indices, dtype=torch.long),
        )

        asym_ids = np.unique(cropped_atom_array.asym_id_int)
        if len(asym_ids) != 1:
            raise ValueError(
                f"Expected one asym_id_int for split RNA chain {chain_id}, got {asym_ids.tolist()}."
            )
        entity_ids = np.unique(cropped_atom_array.label_entity_id)
        if len(entity_ids) != 1:
            raise ValueError(
                f"Expected one label_entity_id for split RNA chain {chain_id}, got {entity_ids.tolist()}."
            )

        asym_id = int(asym_ids[0])
        entity_id = str(entity_ids[0])
        pdb_id = bioassembly_dict["pdb_id"]
        bioassembly_name = make_bioassembly_name(pdb_id, asym_id)
        num_tokens = int(len(cropped_token_array))

        sample_row = dict(row)
        sample_row["pdb_id"] = pdb_id
        sample_row["assembly_id"] = bioassembly_dict.get("assembly_id", "1")
        sample_row["release_date"] = bioassembly_dict.get("release_date", "")
        sample_row["resolution"] = bioassembly_dict.get("resolution", -1)
        sample_row["bioassembly_name"] = bioassembly_name
        sample_row["num_tokens"] = num_tokens
        sample_row["num_prot_chains"] = 0

        single_chain_bioassembly = {
            "pdb_id": pdb_id,
            "bioassembly_name": bioassembly_name,
            "assembly_id": bioassembly_dict.get("assembly_id", "1"),
            "sequences": {entity_id: bioassembly_dict["sequences"][entity_id]},
            "release_date": bioassembly_dict.get("release_date"),
            "num_assembly_polymer_chains": 1,
            "num_prot_chains": 0,
            "entity_poly_type": {
                entity_id: bioassembly_dict["entity_poly_type"][entity_id]
            },
            "resolution": bioassembly_dict.get("resolution", -1),
            "num_tokens": num_tokens,
            "atom_array": cropped_atom_array,
            "token_array": cropped_token_array,
            "msa_features": None,
            "template_features": None,
        }
        if "num_asym_chains" in bioassembly_dict:
            single_chain_bioassembly["num_asym_chains"] = 1
        if "max_res_num_per_chain" in bioassembly_dict:
            single_chain_bioassembly["max_res_num_per_chain"] = num_tokens

        samples.append((bioassembly_name, sample_row, single_chain_bioassembly))
    return samples


def filter_sample_indices_by_rna_mode(
    sample_indices_list: list[dict],
    rna_mode: str = "all",
) -> list[dict]:
    """
    Filters chain/interface rows into RNA-centric training subsets.

    Modes:
        - all: keep original behavior.
        - rna_monomer: keep RNA chain rows and split each RNA chain into its own
          single-chain sample, regardless of whether the source assembly was a
          monomer, RNA complex, RNA-DNA complex, or RNA-protein complex.
        - rna_only: keep RNA-only assemblies (RNA monomer or RNA-RNA complex); output
          RNA chain rows and RNA-RNA interface rows.
        - rna_all: keep RNA-containing assemblies without DNA/hybrid polymer chains;
          output RNA chain rows, RNA-RNA interfaces, and RNA-protein interfaces.
    """
    if rna_mode == "all" or not sample_indices_list:
        return sample_indices_list
    if rna_mode not in RNA_MODE_CHOICES:
        raise ValueError(f"Unsupported rna_mode={rna_mode!r}")

    chain_rows = [row for row in sample_indices_list if row.get("type") == "chain"]
    if rna_mode == "rna_monomer":
        return [row for row in chain_rows if _is_rna_chain_row(row)]

    rna_chain_count = sum(_is_rna_chain_row(row) for row in chain_rows)
    protein_chain_count = sum(row.get("mol_1_type") == "prot" for row in chain_rows)
    non_rna_nuc_chain_count = sum(
        row.get("mol_1_type") == "nuc" and not _is_rna_sub_mol_type(row.get("sub_mol_1_type", ""))
        for row in chain_rows
    )

    if non_rna_nuc_chain_count > 0 or rna_chain_count == 0:
        return []

    if rna_mode == "rna_only":
        if protein_chain_count != 0:
            return []
        return [
            row
            for row in sample_indices_list
            if _is_rna_chain_row(row) or _is_rna_rna_interface_row(row)
        ]

    return [
        row
        for row in sample_indices_list
        if _is_rna_chain_row(row)
        or _is_rna_rna_interface_row(row)
        or _is_rna_protein_interface_row(row)
    ]


def gen_a_bioassembly_data(
    mmcif: Path,
    bioassembly_output_dir: Path,
    cluster_file: Optional[Path],
    distillation: bool = False,
    rna_mode: str = "all",
) -> Optional[list[dict]]:
    """
    Generates bioassembly data from an mmCIF file and saves it to the specified output directory.

    Args:
        mmcif (Path): Path to the mmCIF file.
        bioassembly_output_dir (Path): Directory where the bioassembly data will be saved.
        cluster_file (Optional[Path]): Path to the cluster file, if available.
        distillation (bool, optional): Flag indicating whether to use the 'Distillation' setting. Defaults to False.

    Returns:
        Optional[list[dict]]: A list of sample indices if data is successfully generated, otherwise None.
    """
    if distillation:
        dataset = "Distillation"
    else:
        dataset = "WeightedPDB"

    from protenix.data.pipeline.data_pipeline import DataPipeline

    sample_indices_list, bioassembly_dict = DataPipeline.get_data_from_mmcif(
        mmcif, cluster_file, dataset
    )
    if not sample_indices_list or not bioassembly_dict:
        return None

    if rna_mode == "rna_monomer":
        rna_monomer_samples = build_rna_monomer_samples(
            bioassembly_dict=bioassembly_dict,
            sample_indices_list=sample_indices_list,
        )
        for bioassembly_name, _sample_row, single_chain_bioassembly in rna_monomer_samples:
            from protenix.utils.file_io import dump_gzip_pickle

            dump_gzip_pickle(
                single_chain_bioassembly,
                bioassembly_output_dir / f"{bioassembly_name}.pkl.gz",
            )
        return [sample_row for _name, sample_row, _dict in rna_monomer_samples]

    sample_indices_list = filter_sample_indices_by_rna_mode(
        sample_indices_list, rna_mode=rna_mode
    )

    if sample_indices_list and bioassembly_dict:
        from protenix.utils.file_io import dump_gzip_pickle

        pdb_id = bioassembly_dict["pdb_id"]
        # save to output dir
        dump_gzip_pickle(bioassembly_dict, bioassembly_output_dir / f"{pdb_id}.pkl.gz")
        return sample_indices_list


def gen_data_from_mmcifs(
    mmcif_list: list[Path],
    output_indices_csv: Path,
    bioassembly_output_dir: Path,
    cluster_file: Optional[Path],
    distillation: bool = False,
    rna_mode: str = "all",
    num_workers: int = 1,
):
    """
    Generates training data from a list of mmCIF files and saves the results to a CSV file.

    Args:
        mmcif_list (list[Path]): List of paths to mmCIF files.
        output_indices_csv (Path): Path to the output CSV file where the indices will be saved.
        bioassembly_output_dir (Path): Directory where the bioassembly output will be stored.
        cluster_file (Optional[Path]): Path to the cluster file. If None, clustering is not performed.
        distillation (bool, optional): Flag indicating whether to use the 'Distillation' setting. Defaults to False.
        num_workers (int, optional): Number of parallel workers to use. Defaults to 1.
    """

    all_sample_indices_list = [
        r
        for r in tqdm(
            Parallel(n_jobs=num_workers, return_as="generator_unordered")(
                delayed(gen_a_bioassembly_data)(
                    mmcif,
                    bioassembly_output_dir,
                    cluster_file,
                    distillation,
                    rna_mode,
                )
                for mmcif in mmcif_list
            ),
            total=len(mmcif_list),
        )
    ]

    merged_results = []
    for sample_indices_list in all_sample_indices_list:
        if sample_indices_list:
            merged_results += sample_indices_list
    df = pd.DataFrame(merged_results)

    df.to_csv(output_indices_csv, index=False, quoting=csv.QUOTE_NONNUMERIC)


def run_gen_data(
    input_path: Path,
    output_indices_csv: Path,
    bioassembly_output_dir: Path,
    cluster_file: Optional[Path],
    distillation: bool = False,
    rna_mode: str = "all",
    num_workers: int = 1,
):
    """
    Generates data from MMCIF files and saves the output to specified locations.

    Args:
        input_path (str): Path to the input directory containing MMCIF files or a text file listing MMCIF file paths.
        output_indices_csv (str): Path to the output CSV file where indices will be saved.
        bioassembly_output_dir (str): Directory where bioassembly outputs will be saved.
        cluster_file (Optional[str]): Path to the cluster file, if any.
        distillation (bool, optional): Flag indicating whether to use the 'Distillation' setting. Defaults to False.
        num_workers (int, optional): Number of worker processes to use. Defaults to 1.

    Raises:
        NotImplementedError: If the input path is not a directory or a text file.
    """

    input_path = Path(input_path)
    bioassembly_output_dir = Path(bioassembly_output_dir)
    output_indices_csv = Path(output_indices_csv)

    # create directory for output
    output_indices_csv.parent.mkdir(parents=True, exist_ok=True)
    bioassembly_output_dir.mkdir(parents=True, exist_ok=True)

    if input_path.is_dir():
        mmcif_list = list(input_path.glob("*.cif")) + list(input_path.glob("*.cif.gz"))
    elif input_path.suffix == ".txt":
        with open(input_path) as f:
            mmcif_list = [i.strip() for i in f.readlines()]
    else:
        raise NotImplementedError(f"Unsupported input path: {input_path}")

    gen_data_from_mmcifs(
        mmcif_list,
        output_indices_csv,
        bioassembly_output_dir,
        cluster_file,
        distillation,
        rna_mode,
        num_workers,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_path",
        type=Path,
        default=None,
        help="Path to the input directory containing MMCIF files or a .txt file listing MMCIF file paths.",
    )
    parser.add_argument(
        "-o",
        "--output_csv",
        type=Path,
        default=None,
        help="Path to the output CSV file where indices will be saved.",
    )
    parser.add_argument(
        "-b",
        "--bio_output_dir",
        type=Path,
        default=None,
        help="Directory where bioassembly outputs will be saved.",
    )
    parser.add_argument(
        "-c",
        "--cluster_file",
        type=Path,
        default=None,
        help="Path to the cluster txt file, if any",
    )

    parser.add_argument(
        "-d",
        "--distillation",
        action="store_true",
        help="Whether to use the 'Distillation' setting",
    )

    parser.add_argument(
        "-n",
        "--n_cpu",
        type=int,
        default=1,
        help="Number of worker processes to use. Defaults to 1.",
    )
    parser.add_argument(
        "--rna-mode",
        choices=RNA_MODE_CHOICES,
        default="all",
        help=(
            "Optional RNA-specific filtering mode. "
            "'rna_monomer' keeps only RNA monomers; "
            "'rna_only' keeps RNA-only monomers/complexes; "
            "'rna_all' keeps RNA monomers, RNA-RNA complexes, and RNA-protein complexes. "
            "Default: all."
        ),
    )

    args = parser.parse_args()

    run_gen_data(
        input_path=args.input_path,
        output_indices_csv=args.output_csv,
        bioassembly_output_dir=args.bio_output_dir,
        cluster_file=args.cluster_file,
        distillation=args.distillation,
        rna_mode=args.rna_mode,
        num_workers=args.n_cpu,
    )
