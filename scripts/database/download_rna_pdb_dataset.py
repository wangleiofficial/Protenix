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
import json
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable, List, Optional

SEARCH_API_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DEFAULT_ROWS_PER_PAGE = 1000
DEFAULT_USER_AGENT = "ProtenixRNAPDBDataset/1.0"
DEFAULT_CIF_BASE_URL = "https://files.rcsb.org/download"


def _make_terminal_node(attribute: str, operator: str, value: Any) -> dict[str, Any]:
    """Builds one RCSB Search API terminal node."""
    return {
        "type": "terminal",
        "service": "text",
        "parameters": {
            "attribute": attribute,
            "operator": operator,
            "value": value,
        },
    }


def build_rna_entry_query(
    *,
    max_resolution: Optional[float] = 4.5,
    release_date_before: Optional[str] = None,
    limit: Optional[int] = None,
    start: int = 0,
    rows: int = DEFAULT_ROWS_PER_PAGE,
) -> dict[str, Any]:
    """Builds an RCSB Search API query for RNA-containing entries."""
    nodes = [
        _make_terminal_node(
            "rcsb_entry_info.polymer_entity_count_RNA",
            "greater",
            0,
        )
    ]
    if max_resolution is not None:
        nodes.append(
            _make_terminal_node(
                "rcsb_entry_info.resolution_combined",
                "less_or_equal",
                max_resolution,
            )
        )
    if release_date_before is not None:
        nodes.append(
            _make_terminal_node(
                "rcsb_accession_info.initial_release_date",
                "less_or_equal",
                release_date_before,
            )
        )

    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": nodes,
        },
        "return_type": "entry",
        "request_options": {
            "results_verbosity": "compact",
            "paginate": {
                "start": start,
                "rows": rows if limit is None else min(rows, limit),
            },
        },
    }


def _post_search_query(
    payload: dict[str, Any],
    *,
    timeout: int = 60,
    user_agent: str = DEFAULT_USER_AGENT,
) -> dict[str, Any]:
    """Posts a query to the RCSB Search API."""
    request = urllib.request.Request(
        SEARCH_API_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Content-Type": "application/json",
            "User-Agent": user_agent,
        },
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.load(response)


def search_rna_entry_ids(
    *,
    max_resolution: Optional[float] = 4.5,
    release_date_before: Optional[str] = None,
    limit: Optional[int] = None,
    rows_per_page: int = DEFAULT_ROWS_PER_PAGE,
    timeout: int = 60,
) -> List[str]:
    """Searches RCSB for RNA-containing entry IDs."""
    start = 0
    remaining = limit
    entry_ids: List[str] = []

    while True:
        payload = build_rna_entry_query(
            max_resolution=max_resolution,
            release_date_before=release_date_before,
            limit=remaining,
            start=start,
            rows=rows_per_page,
        )
        result = _post_search_query(payload, timeout=timeout)
        result_set = result.get("result_set", [])
        batch = [row["identifier"].lower() for row in result_set if "identifier" in row]
        entry_ids.extend(batch)

        if limit is not None:
            remaining = limit - len(entry_ids)
            if remaining <= 0:
                break

        if not batch:
            break

        total_count = int(result.get("total_count", len(entry_ids)))
        start += len(batch)
        if start >= total_count:
            break

    if limit is not None:
        return entry_ids[:limit]
    return entry_ids


def write_lines(lines: Iterable[str], output_path: Path) -> None:
    """Writes one item per line."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        for line in lines:
            f.write(f"{line}\n")


def download_one_pdb_cif(
    pdb_id: str,
    output_dir: Path,
    *,
    base_url: str = DEFAULT_CIF_BASE_URL,
    timeout: int = 60,
    skip_existing: bool = True,
    user_agent: str = DEFAULT_USER_AGENT,
) -> Path:
    """Downloads one mmCIF file from RCSB."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{pdb_id}.cif"
    if skip_existing and output_path.exists():
        return output_path

    request = urllib.request.Request(
        f"{base_url.rstrip('/')}/{pdb_id}.cif",
        headers={"User-Agent": user_agent},
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        output_path.write_bytes(response.read())
    return output_path


def download_required_pdb_cifs(
    pdb_ids: Iterable[str],
    output_dir: Path,
    *,
    base_url: str = DEFAULT_CIF_BASE_URL,
    timeout: int = 60,
    skip_existing: bool = True,
    workers: int = 4,
) -> tuple[List[Path], List[str]]:
    """Downloads all required mmCIF files and returns (successes, failed_ids)."""
    unique_pdb_ids = sorted({pdb_id.lower() for pdb_id in pdb_ids})
    if workers <= 1:
        successes: List[Path] = []
        failures: List[str] = []
        for pdb_id in unique_pdb_ids:
            try:
                successes.append(
                    download_one_pdb_cif(
                        pdb_id,
                        output_dir,
                        base_url=base_url,
                        timeout=timeout,
                        skip_existing=skip_existing,
                    )
                )
            except Exception:
                failures.append(pdb_id)
        return successes, failures

    successes = []
    failures = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_pdb = {
            executor.submit(
                download_one_pdb_cif,
                pdb_id,
                output_dir,
                base_url=base_url,
                timeout=timeout,
                skip_existing=skip_existing,
            ): pdb_id
            for pdb_id in unique_pdb_ids
        }
        for future in as_completed(future_to_pdb):
            pdb_id = future_to_pdb[future]
            try:
                successes.append(future.result())
            except Exception:
                failures.append(pdb_id)
    return sorted(successes), sorted(failures)


def build_argparser() -> argparse.ArgumentParser:
    """Builds the CLI parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Search RNA-containing entries from RCSB PDB and download their mmCIF files."
        )
    )
    parser.add_argument(
        "--mmcif-dir",
        type=Path,
        required=True,
        help="Directory where downloaded mmCIF files will be stored.",
    )
    parser.add_argument(
        "--entry-id-list-out",
        type=Path,
        default=None,
        help="Optional text file for saving searched PDB IDs, one per line.",
    )
    parser.add_argument(
        "--mmcif-path-list-out",
        type=Path,
        default=None,
        help="Optional text file for saving downloaded mmCIF paths, one per line.",
    )
    parser.add_argument(
        "--max-resolution",
        type=float,
        default=4.5,
        help=(
            "Maximum allowed resolution for RCSB search. "
            "Use a negative value to disable this filter."
        ),
    )
    parser.add_argument(
        "--release-date-before",
        type=str,
        default=None,
        help="Optional release-date cutoff in YYYY-MM-DD format.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional maximum number of RNA-containing entries to retrieve.",
    )
    parser.add_argument(
        "--rows-per-page",
        type=int,
        default=DEFAULT_ROWS_PER_PAGE,
        help="Pagination size for the RCSB Search API.",
    )
    parser.add_argument(
        "--download-workers",
        type=int,
        default=8,
        help="Number of concurrent workers for mmCIF downloads.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=60,
        help="Network timeout in seconds.",
    )
    parser.add_argument(
        "--no-skip-existing",
        action="store_true",
        help="Re-download mmCIF files even if they already exist locally.",
    )
    return parser


def main() -> None:
    """CLI entrypoint."""
    args = build_argparser().parse_args()

    max_resolution = args.max_resolution
    if max_resolution is not None and max_resolution < 0:
        max_resolution = None

    entry_ids = search_rna_entry_ids(
        max_resolution=max_resolution,
        release_date_before=args.release_date_before,
        limit=args.limit,
        rows_per_page=args.rows_per_page,
        timeout=args.timeout,
    )
    if not entry_ids:
        raise SystemExit("No RNA-containing PDB entries were found.")

    print(f"Found {len(entry_ids)} RNA-containing entries from RCSB.")
    if args.entry_id_list_out is not None:
        write_lines(entry_ids, args.entry_id_list_out)

    downloaded_paths, failed_ids = download_required_pdb_cifs(
        entry_ids,
        args.mmcif_dir,
        timeout=args.timeout,
        skip_existing=not args.no_skip_existing,
        workers=args.download_workers,
    )
    print(f"Downloaded or verified {len(downloaded_paths)} mmCIF files.")
    if args.mmcif_path_list_out is not None:
        write_lines((str(path) for path in downloaded_paths), args.mmcif_path_list_out)

    if failed_ids:
        print("Failed downloads:")
        for pdb_id in failed_ids:
            print(f"  - {pdb_id}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
