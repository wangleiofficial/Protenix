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
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


UNPAIRED_DOT_BRACKET = set(".-_:,")
BRACKET_PAIRS = {
    "(": ")",
    "[": "]",
    "{": "}",
    "<": ">",
}
BRACKET_PAIRS.update({chr(ord("A") + i): chr(ord("a") + i) for i in range(26)})
CLOSE_TO_OPEN = {close: open_ for open_, close in BRACKET_PAIRS.items()}


def read_nonempty_lines(text: str) -> List[str]:
    """Returns non-empty, non-comment lines."""
    return [
        line.strip()
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def is_dot_bracket_line(line: str) -> bool:
    """Checks whether a line looks like dot-bracket notation."""
    line = line.replace(" ", "")
    valid = UNPAIRED_DOT_BRACKET | set(BRACKET_PAIRS) | set(CLOSE_TO_OPEN)
    return bool(line) and all(char in valid for char in line)


def detect_input_format(path: Path, explicit_format: str) -> str:
    """Infers the input format from suffix or file content."""
    if explicit_format != "auto":
        return explicit_format

    suffix = path.suffix.lower()
    if suffix in {".dbn", ".dot", ".dotbracket"}:
        return "dotbracket"
    if suffix == ".bpseq":
        return "bpseq"
    if suffix in {".csv", ".tsv", ".pairs", ".lst"}:
        return "pairs"

    text = path.read_text()
    lines = read_nonempty_lines(text)
    if not lines:
        raise ValueError("Input file is empty.")
    if lines[0].startswith(">"):
        return "dotbracket"

    first_fields = re.split(r"[\s,;]+", lines[0])
    if is_dot_bracket_line(lines[-1]):
        return "dotbracket"
    if (
        len(first_fields) >= 3
        and re.fullmatch(r"\d+", first_fields[0])
        and re.fullmatch(r"[A-Za-z]", first_fields[1])
        and re.fullmatch(r"-?\d+", first_fields[2])
    ):
        return "bpseq"
    if len(first_fields) >= 2:
        return "pairs"
    raise ValueError(
        f"Cannot infer input format from '{path}'. Use --input-format explicitly."
    )


def normalize_pairs(
    pairs: Iterable[Tuple[int, int]],
    *,
    pair_index_base: int = 1,
) -> List[Tuple[int, int]]:
    """Normalizes, deduplicates, and sorts base pairs in 1-based indexing."""
    normalized = set()
    for left, right in pairs:
        if pair_index_base == 0:
            left += 1
            right += 1
        if left <= 0 or right <= 0:
            raise ValueError(f"Base-pair positions must be positive. Got ({left}, {right}).")
        if left == right:
            raise ValueError(f"Self-pair is invalid: ({left}, {right}).")
        if left > right:
            left, right = right, left
        normalized.add((left, right))
    return sorted(normalized)


def parse_dot_bracket(text: str) -> List[Tuple[int, int]]:
    """Parses raw dot-bracket text into 1-based base pairs."""
    lines = read_nonempty_lines(text)
    if not lines:
        raise ValueError("Dot-bracket input is empty.")

    if lines[0].startswith(">"):
        lines = lines[1:]
    structure_lines = [line.replace(" ", "") for line in lines if is_dot_bracket_line(line)]
    if not structure_lines:
        raise ValueError("No dot-bracket structure line found.")

    structure = structure_lines[-1]
    stacks: Dict[str, List[int]] = {open_: [] for open_ in BRACKET_PAIRS}
    pairs: List[Tuple[int, int]] = []

    for idx, char in enumerate(structure, start=1):
        if char in BRACKET_PAIRS:
            stacks[char].append(idx)
        elif char in CLOSE_TO_OPEN:
            open_ = CLOSE_TO_OPEN[char]
            if not stacks[open_]:
                raise ValueError(f"Unmatched closing bracket '{char}' at position {idx}.")
            pairs.append((stacks[open_].pop(), idx))
        elif char not in UNPAIRED_DOT_BRACKET:
            raise ValueError(f"Unsupported dot-bracket character '{char}' at position {idx}.")

    unclosed = [(open_, positions) for open_, positions in stacks.items() if positions]
    if unclosed:
        open_, positions = unclosed[0]
        raise ValueError(f"Unmatched opening bracket '{open_}' at positions {positions}.")

    return normalize_pairs(pairs)


def parse_bpseq(text: str) -> List[Tuple[int, int]]:
    """Parses BPSEQ text into 1-based base pairs."""
    pairs: List[Tuple[int, int]] = []
    for line in read_nonempty_lines(text):
        fields = line.split()
        if len(fields) < 3:
            raise ValueError(f"Invalid BPSEQ line: {line}")
        if not re.fullmatch(r"\d+", fields[0]) or not re.fullmatch(r"-?\d+", fields[2]):
            raise ValueError(f"Invalid BPSEQ line: {line}")
        position = int(fields[0])
        partner = int(fields[2])
        if partner > position:
            pairs.append((position, partner))
    return normalize_pairs(pairs)


def parse_pair_list(text: str, *, pair_index_base: int = 1) -> List[Tuple[int, int]]:
    """Parses a simple base-pair list from csv/tsv/whitespace text."""
    pairs: List[Tuple[int, int]] = []
    saw_pair = False
    for line in read_nonempty_lines(text):
        fields = [field for field in re.split(r"[\s,;]+", line) if field]
        if len(fields) < 2:
            continue
        if not re.fullmatch(r"-?\d+", fields[0]) or not re.fullmatch(r"-?\d+", fields[1]):
            continue
        saw_pair = True
        pairs.append((int(fields[0]), int(fields[1])))
    if not saw_pair:
        raise ValueError("No base-pair rows found in pair-list input.")
    return normalize_pairs(pairs, pair_index_base=pair_index_base)


def parse_input(
    input_path: Path,
    *,
    input_format: str = "auto",
    pair_index_base: int = 1,
) -> List[Tuple[int, int]]:
    """Parses a secondary-structure file into base pairs."""
    effective_format = detect_input_format(input_path, input_format)
    text = input_path.read_text()
    if effective_format == "dotbracket":
        return parse_dot_bracket(text)
    if effective_format == "bpseq":
        return parse_bpseq(text)
    if effective_format == "pairs":
        return parse_pair_list(text, pair_index_base=pair_index_base)
    raise ValueError(f"Unsupported input format: {effective_format}")


def build_contact_json(
    pairs: Sequence[Tuple[int, int]],
    *,
    entity1: int,
    copy1: int,
    entity2: int,
    copy2: int,
    max_distance: float,
    min_distance: float,
    atom1: Optional[str] = None,
    atom2: Optional[str] = None,
    wrap: str = "contact",
) -> Dict[str, object]:
    """Converts base pairs into Protenix contact constraints."""
    contacts = []
    for left, right in pairs:
        contact = {
            "entity1": entity1,
            "copy1": copy1,
            "position1": left,
            "entity2": entity2,
            "copy2": copy2,
            "position2": right,
            "max_distance": max_distance,
            "min_distance": min_distance,
        }
        if atom1:
            contact["atom1"] = atom1
        if atom2:
            contact["atom2"] = atom2
        contacts.append(contact)

    if wrap == "contact":
        return {"contact": contacts}
    if wrap == "constraint":
        return {"constraint": {"contact": contacts}}
    raise ValueError(f"Unsupported wrap mode: {wrap}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert RNA secondary structure formats (dot-bracket, bpseq, or "
            "base-pair list) into Protenix constraint.contact JSON."
        )
    )
    parser.add_argument("--input", required=True, help="Input secondary-structure file.")
    parser.add_argument(
        "--input-format",
        choices=["auto", "dotbracket", "bpseq", "pairs"],
        default="auto",
        help="Input format. Default infers from suffix or content.",
    )
    parser.add_argument(
        "--output",
        default="-",
        help="Output JSON path. Use '-' to write to stdout.",
    )
    parser.add_argument(
        "--wrap",
        choices=["contact", "constraint"],
        default="contact",
        help="Output {'contact': [...]} or {'constraint': {'contact': [...]}}.",
    )
    parser.add_argument(
        "--entity1",
        type=int,
        default=1,
        help="Entity ID for the first residue in each pair.",
    )
    parser.add_argument(
        "--copy1",
        type=int,
        default=1,
        help="Copy ID for the first residue in each pair.",
    )
    parser.add_argument(
        "--entity2",
        type=int,
        default=1,
        help="Entity ID for the second residue in each pair.",
    )
    parser.add_argument(
        "--copy2",
        type=int,
        default=1,
        help="Copy ID for the second residue in each pair.",
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=10.0,
        help="Constraint max_distance value in Angstrom.",
    )
    parser.add_argument(
        "--min-distance",
        type=float,
        default=0.0,
        help="Constraint min_distance value in Angstrom.",
    )
    parser.add_argument(
        "--atom1",
        default=None,
        help="Optional atom name for residue 1, e.g. C4 or C3'.",
    )
    parser.add_argument(
        "--atom2",
        default=None,
        help="Optional atom name for residue 2, e.g. C4 or C3'.",
    )
    parser.add_argument(
        "--pair-index-base",
        choices=[0, 1],
        type=int,
        default=1,
        help="Index base for pair-list inputs only. Dot-bracket and BPSEQ are always 1-based.",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    pairs = parse_input(
        Path(args.input),
        input_format=args.input_format,
        pair_index_base=args.pair_index_base,
    )
    payload = build_contact_json(
        pairs,
        entity1=args.entity1,
        copy1=args.copy1,
        entity2=args.entity2,
        copy2=args.copy2,
        max_distance=args.max_distance,
        min_distance=args.min_distance,
        atom1=args.atom1,
        atom2=args.atom2,
        wrap=args.wrap,
    )

    output = json.dumps(payload, indent=2, sort_keys=False) + "\n"
    if args.output == "-":
        print(output, end="")
    else:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(output)


if __name__ == "__main__":
    main()
