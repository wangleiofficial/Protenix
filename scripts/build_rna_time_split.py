#!/usr/bin/env python3
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
import gzip
import pickle
from pathlib import Path
from typing import Any

import pandas as pd

RNA_SUBTYPES = {"rna", "modified_rna"}


def load_gzip_pickle(path: Path) -> Any:
    """Loads a gzip pickle file without importing the full training stack."""
    with gzip.open(path, "rb") as f:
        return pickle.load(f)


def read_indices_csv(path: Path) -> pd.DataFrame:
    """Reads an indices CSV as strings."""
    return pd.read_csv(path, keep_default_na=False, dtype=str)


def write_indices_csv(df: pd.DataFrame, path: Path) -> None:
    """Writes an indices DataFrame with stable quoting."""
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, quoting=csv.QUOTE_NONNUMERIC)


def filter_to_rna_chain_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Keeps only RNA chain rows suitable for monomer RNA evaluation."""
    mask = (df["type"] == "chain") & (df["sub_mol_1_type"].isin(RNA_SUBTYPES))
    return df[mask].copy()


def get_bioassembly_stem(row: pd.Series) -> str:
    """Returns the file stem used to locate the bioassembly pickle."""
    bioassembly_name = row.get("bioassembly_name", "")
    return bioassembly_name or row["pdb_id"]


def resolve_sequence_for_row(
    row: pd.Series,
    bioassembly_dir: Path,
    bioassembly_cache: dict[str, dict[str, Any]],
) -> str:
    """Resolves the canonical sequence for one indices row."""
    stem = get_bioassembly_stem(row)
    if stem not in bioassembly_cache:
        bioassembly_cache[stem] = load_gzip_pickle(bioassembly_dir / f"{stem}.pkl.gz")

    bioassembly = bioassembly_cache[stem]
    entity_id = str(row["entity_1_id"])
    sequences = bioassembly.get("sequences", {})
    if entity_id in sequences:
        return str(sequences[entity_id])

    if len(sequences) == 1:
        return str(next(iter(sequences.values())))

    raise KeyError(
        f"Unable to resolve sequence for entity_1_id={entity_id} in bioassembly={stem}."
    )


def attach_canonical_sequences(
    df: pd.DataFrame,
    bioassembly_dir: Path,
) -> pd.DataFrame:
    """Adds a canonical_sequence column by reading saved bioassembly pickles."""
    bioassembly_cache: dict[str, dict[str, Any]] = {}
    out = df.copy()
    out["canonical_sequence"] = out.apply(
        lambda row: resolve_sequence_for_row(row, bioassembly_dir, bioassembly_cache),
        axis=1,
    )
    return out


def split_by_release_date(
    df: pd.DataFrame,
    test_release_date_from: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Splits rows into train/test by an inclusive test start date."""
    release_dates = pd.to_datetime(df["release_date"], format="%Y-%m-%d")
    cutoff = pd.Timestamp(test_release_date_from)
    train_df = df[release_dates < cutoff].copy()
    test_df = df[release_dates >= cutoff].copy()
    return train_df, test_df


def remove_test_sequences_seen_in_train(
    train_df: pd.DataFrame,
    test_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Removes test rows whose canonical sequence is identical to any training row."""
    train_sequences = set(train_df["canonical_sequence"].tolist())
    duplicate_mask = test_df["canonical_sequence"].isin(train_sequences)
    dropped_test_df = test_df[duplicate_mask].copy()
    kept_test_df = test_df[~duplicate_mask].copy()
    return kept_test_df, dropped_test_df


def build_time_split(
    *,
    indices_path: Path,
    bioassembly_dir: Path,
    test_release_date_from: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Builds train/test RNA splits with exact-sequence filtering on the test set."""
    df = read_indices_csv(indices_path)
    df = filter_to_rna_chain_rows(df)
    if len(df) == 0:
        raise ValueError("No RNA chain rows were found in the input indices CSV.")

    if "release_date" not in df.columns:
        raise ValueError("The input indices CSV must contain a release_date column.")

    df = attach_canonical_sequences(df, bioassembly_dir=bioassembly_dir)
    train_df, test_df = split_by_release_date(
        df, test_release_date_from=test_release_date_from
    )
    if len(train_df) == 0:
        raise ValueError("The training split is empty after applying the release-date cutoff.")
    if len(test_df) == 0:
        raise ValueError("The test split is empty after applying the release-date cutoff.")

    filtered_test_df, dropped_test_df = remove_test_sequences_seen_in_train(
        train_df=train_df, test_df=test_df
    )
    if len(filtered_test_df) == 0:
        raise ValueError(
            "The test split became empty after removing sequences identical to the training split."
        )
    return train_df, filtered_test_df, dropped_test_df


def write_pdb_id_list(df: pd.DataFrame, output_path: Path) -> None:
    """Writes sorted unique PDB IDs, one per line."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    unique_pdb_ids = sorted(set(df["pdb_id"].tolist()))
    output_path.write_text("".join(f"{pdb_id}\n" for pdb_id in unique_pdb_ids))


def build_argparser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a monomer RNA train/test split by release date and remove "
            "test samples whose canonical sequence is identical to any training sample."
        )
    )
    parser.add_argument(
        "--indices",
        type=Path,
        required=True,
        help="Input indices CSV produced by prepare_training_data.py.",
    )
    parser.add_argument(
        "--bioassembly-dir",
        type=Path,
        required=True,
        help="Directory containing saved Bioassembly Dict .pkl.gz files.",
    )
    parser.add_argument(
        "--test-release-date-from",
        type=str,
        required=True,
        help="Inclusive release-date threshold (YYYY-MM-DD) for the test split.",
    )
    parser.add_argument(
        "--train-out",
        type=Path,
        required=True,
        help="Output CSV path for the training split.",
    )
    parser.add_argument(
        "--test-out",
        type=Path,
        required=True,
        help="Output CSV path for the filtered test split.",
    )
    parser.add_argument(
        "--dropped-test-out",
        type=Path,
        default=None,
        help="Optional CSV path for test rows dropped due to exact sequence matches in train.",
    )
    parser.add_argument(
        "--test-pdb-list-out",
        type=Path,
        default=None,
        help="Optional text file containing one unique test PDB ID per line.",
    )
    return parser


def main() -> None:
    args = build_argparser().parse_args()
    train_df, test_df, dropped_test_df = build_time_split(
        indices_path=args.indices,
        bioassembly_dir=args.bioassembly_dir,
        test_release_date_from=args.test_release_date_from,
    )

    write_indices_csv(train_df, args.train_out)
    write_indices_csv(test_df, args.test_out)
    if args.dropped_test_out is not None:
        write_indices_csv(dropped_test_df, args.dropped_test_out)
    if args.test_pdb_list_out is not None:
        write_pdb_id_list(test_df, args.test_pdb_list_out)

    print(
        f"train_rows={len(train_df)} test_rows={len(test_df)} "
        f"dropped_test_rows={len(dropped_test_df)}"
    )


if __name__ == "__main__":
    main()
