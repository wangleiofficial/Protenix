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
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence


M8_COLUMNS = [
    "query",
    "target",
    "fident",
    "alnlen",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "evalue",
    "bits",
]


def detect_input_format(path: Path, explicit_format: str) -> str:
    """Returns the effective input format."""
    if explicit_format != "auto":
        return explicit_format
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return "csv"
    if suffix == ".m8":
        return "m8"
    raise ValueError(
        f"Cannot infer input format from suffix '{path.suffix}'. "
        "Use --input-format csv or --input-format m8."
    )


def load_query_map(path: Optional[Path]) -> Optional[Dict[str, str]]:
    """Loads a query -> output key mapping from JSON."""
    if path is None:
        return None
    with open(path, "r") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("query map JSON must be an object: {source_query: output_key}")
    return {str(k): str(v) for k, v in data.items()}


def read_rows(input_path: Path, input_format: str) -> List[Dict[str, str]]:
    """Reads RNA template hit rows from csv or m8."""
    rows: List[Dict[str, str]] = []
    with open(input_path, "r", newline="") as f:
        if input_format == "csv":
            reader = csv.DictReader(f)
            missing = [col for col in M8_COLUMNS if col not in (reader.fieldnames or [])]
            if missing:
                raise ValueError(
                    f"CSV file is missing required columns: {', '.join(missing)}"
                )
            for row in reader:
                rows.append({col: row[col] for col in M8_COLUMNS})
            return rows

        if input_format == "m8":
            for line_number, line in enumerate(f, start=1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) < len(M8_COLUMNS):
                    raise ValueError(
                        f"Invalid m8 line {line_number}: expected at least "
                        f"{len(M8_COLUMNS)} columns, got {len(fields)}"
                    )
                rows.append(dict(zip(M8_COLUMNS, fields[: len(M8_COLUMNS)])))
            return rows

    raise ValueError(f"Unsupported input format: {input_format}")


def sanitize_filename(name: str) -> str:
    """Converts an output key into a filesystem-safe stem."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("._")
    return stem or "query"


def unique_stem(stem: str, original_key: str, used: Dict[str, str]) -> str:
    """Avoids collisions when multiple keys sanitize to the same filename."""
    current = used.get(stem)
    if current is None or current == original_key:
        used[stem] = original_key
        return stem
    digest = hashlib.sha1(original_key.encode("utf-8")).hexdigest()[:8]
    candidate = f"{stem}_{digest}"
    used[candidate] = original_key
    return candidate


def extract_pdb_id(identifier: str) -> Optional[str]:
    """Extracts the leading 4-character PDB ID from a query/target identifier."""
    token = identifier.strip()
    if len(token) < 4:
        return None
    pdb_id = token[:4].lower()
    if not re.fullmatch(r"[a-zA-Z0-9]{4}", pdb_id):
        return None
    return pdb_id


def filter_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    drop_self_hit: bool,
    drop_same_pdb: bool,
) -> tuple[List[Dict[str, str]], Dict[str, int]]:
    """Filters rows before splitting them into per-query files."""
    filtered_rows: List[Dict[str, str]] = []
    stats = {
        "drop_self_hit": 0,
        "drop_same_pdb": 0,
    }

    for row in rows:
        query = row["query"].strip()
        target = row["target"].strip()

        if drop_self_hit and query == target:
            stats["drop_self_hit"] += 1
            continue

        if drop_same_pdb:
            query_pdb = extract_pdb_id(query)
            target_pdb = extract_pdb_id(target)
            if query_pdb is not None and target_pdb is not None and query_pdb == target_pdb:
                stats["drop_same_pdb"] += 1
                continue

        filtered_rows.append(dict(row))

    return filtered_rows, stats


def group_rows_by_output_key(
    rows: Sequence[Mapping[str, str]],
    query_map: Optional[Mapping[str, str]],
    skip_unmapped: bool,
) -> tuple[Dict[str, List[Dict[str, str]]], List[str]]:
    """Groups rows by the final output key used for writing split files."""
    grouped: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    skipped_queries = set()
    for row in rows:
        query = row["query"]
        if query_map is None:
            output_key = query
        else:
            output_key = query_map.get(query)
            if output_key is None:
                if skip_unmapped:
                    skipped_queries.add(query)
                    continue
                raise KeyError(
                    f"Query '{query}' was not found in the provided query map."
                )
        grouped[output_key].append(dict(row))
    return dict(grouped), sorted(skipped_queries)


def write_rows(rows: Iterable[Mapping[str, str]], output_path: Path, output_format: str) -> None:
    """Writes a query-specific hit table to csv or m8."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="") as f:
        if output_format == "csv":
            writer = csv.DictWriter(f, fieldnames=M8_COLUMNS)
            writer.writeheader()
            for row in rows:
                writer.writerow({col: row[col] for col in M8_COLUMNS})
            return

        if output_format == "m8":
            for row in rows:
                f.write("\t".join(row[col] for col in M8_COLUMNS) + "\n")
            return

    raise ValueError(f"Unsupported output format: {output_format}")


def split_template_hits(
    input_path: Path,
    output_dir: Path,
    mapping_out: Path,
    *,
    input_format: str = "auto",
    output_format: str = "same",
    query_map: Optional[Mapping[str, str]] = None,
    skip_unmapped: bool = False,
    layout: str = "flat",
    path_mode: str = "absolute",
    drop_self_hit: bool = False,
    drop_same_pdb: bool = False,
) -> Dict[str, str]:
    """Splits a combined RNA template hit table into per-query files and a mapping JSON."""
    effective_input_format = detect_input_format(input_path, input_format)
    effective_output_format = (
        effective_input_format if output_format == "same" else output_format
    )

    rows = read_rows(input_path, effective_input_format)
    rows, filter_stats = filter_rows(
        rows,
        drop_self_hit=drop_self_hit,
        drop_same_pdb=drop_same_pdb,
    )
    grouped, skipped_queries = group_rows_by_output_key(rows, query_map, skip_unmapped)

    output_dir.mkdir(parents=True, exist_ok=True)
    mapping: Dict[str, str] = {}
    used_stems: Dict[str, str] = {}

    for output_key, group_rows in sorted(grouped.items()):
        stem = unique_stem(sanitize_filename(output_key), output_key, used_stems)
        if layout == "flat":
            output_path = output_dir / f"{stem}.{effective_output_format}"
            mapping_value_path = output_path
        elif layout == "dir":
            output_path = output_dir / stem / f"result.{effective_output_format}"
            mapping_value_path = output_path.parent
        else:
            raise ValueError(f"Unsupported layout: {layout}")

        write_rows(group_rows, output_path, effective_output_format)
        if path_mode == "absolute":
            mapping[output_key] = str(mapping_value_path.resolve())
        elif path_mode == "relative":
            mapping[output_key] = str(mapping_value_path.relative_to(output_dir))
        else:
            raise ValueError(f"Unsupported path_mode: {path_mode}")

    mapping_out.parent.mkdir(parents=True, exist_ok=True)
    with open(mapping_out, "w") as f:
        json.dump(mapping, f, indent=2, sort_keys=True)
        f.write("\n")

    print(
        "Wrote "
        f"{len(mapping)} split files to {output_dir} "
        f"from {len(rows)} rows across {len(grouped)} output keys."
    )
    if filter_stats["drop_self_hit"]:
        print(f"Dropped {filter_stats['drop_self_hit']} self-hit rows.")
    if filter_stats["drop_same_pdb"]:
        print(f"Dropped {filter_stats['drop_same_pdb']} same-PDB rows.")
    if skipped_queries:
        print(f"Skipped {len(skipped_queries)} unmapped queries.")
    print(f"Wrote mapping JSON to {mapping_out}")
    return mapping


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Split a combined RNA template hit table (csv or m8) into one file per "
            "query/output key and write a mapping JSON for Protenix training."
        )
    )
    parser.add_argument("--input", required=True, help="Input combined csv/m8 file.")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for split query files.",
    )
    parser.add_argument(
        "--mapping-out",
        required=True,
        help="Output JSON path for {output_key: file_or_dir}.",
    )
    parser.add_argument(
        "--input-format",
        choices=["auto", "csv", "m8"],
        default="auto",
        help="Input table format. Default infers from file suffix.",
    )
    parser.add_argument(
        "--output-format",
        choices=["same", "csv", "m8"],
        default="same",
        help="Output file format. Default preserves the input format.",
    )
    parser.add_argument(
        "--query-map-json",
        default=None,
        help=(
            "Optional JSON mapping {source_query: output_key}. "
            "Use this to rename query IDs to training sequence_uid values."
        ),
    )
    parser.add_argument(
        "--skip-unmapped",
        action="store_true",
        help="Skip queries that are missing from --query-map-json instead of failing.",
    )
    parser.add_argument(
        "--layout",
        choices=["flat", "dir"],
        default="flat",
        help=(
            "flat: mapping points to split files. "
            "dir: mapping points to per-query directories containing result.csv/m8."
        ),
    )
    parser.add_argument(
        "--path-mode",
        choices=["absolute", "relative"],
        default="absolute",
        help=(
            "Whether mapping JSON values should be absolute paths or paths "
            "relative to --output-dir."
        ),
    )
    parser.add_argument(
        "--drop-self-hit",
        action="store_true",
        help="Drop rows where query and target are exactly the same identifier.",
    )
    parser.add_argument(
        "--drop-same-pdb",
        action="store_true",
        help="Drop rows where query and target share the same 4-character PDB ID.",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    split_template_hits(
        input_path=Path(args.input),
        output_dir=Path(args.output_dir),
        mapping_out=Path(args.mapping_out),
        input_format=args.input_format,
        output_format=args.output_format,
        query_map=load_query_map(Path(args.query_map_json))
        if args.query_map_json
        else None,
        skip_unmapped=args.skip_unmapped,
        layout=args.layout,
        path_mode=args.path_mode,
        drop_self_hit=args.drop_self_hit,
        drop_same_pdb=args.drop_same_pdb,
    )


if __name__ == "__main__":
    main()
