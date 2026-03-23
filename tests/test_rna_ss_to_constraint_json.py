import json

from scripts.rna_ss_to_constraint_json import (
    build_contact_json,
    parse_bpseq,
    parse_dot_bracket,
    parse_pair_list,
)


def test_parse_dot_bracket_to_contacts():
    pairs = parse_dot_bracket(">rna\nGGGAAACCC\n(((...)))\n")
    assert pairs == [(1, 9), (2, 8), (3, 7)]

    payload = build_contact_json(
        pairs,
        entity1=1,
        copy1=1,
        entity2=1,
        copy2=1,
        max_distance=10.0,
        min_distance=0.0,
    )
    assert payload["contact"][0] == {
        "entity1": 1,
        "copy1": 1,
        "position1": 1,
        "entity2": 1,
        "copy2": 1,
        "position2": 9,
        "max_distance": 10.0,
        "min_distance": 0.0,
    }


def test_parse_bpseq_deduplicates_symmetric_pairs():
    pairs = parse_bpseq(
        "\n".join(
            [
                "1 G 6",
                "2 C 5",
                "3 A 0",
                "4 U 0",
                "5 G 2",
                "6 C 1",
            ]
        )
        + "\n"
    )
    assert pairs == [(1, 6), (2, 5)]


def test_parse_pair_list_supports_csv_and_zero_based_indices():
    pairs = parse_pair_list(
        "i,j,prob\n0,8,0.99\n1,7,0.95\n",
        pair_index_base=0,
    )
    assert pairs == [(1, 9), (2, 8)]

    payload = build_contact_json(
        pairs,
        entity1=2,
        copy1=1,
        entity2=2,
        copy2=1,
        max_distance=8.0,
        min_distance=2.0,
        atom1="C4",
        atom2="C4",
        wrap="constraint",
    )
    assert json.loads(json.dumps(payload)) == {
        "constraint": {
            "contact": [
                {
                    "entity1": 2,
                    "copy1": 1,
                    "position1": 1,
                    "entity2": 2,
                    "copy2": 1,
                    "position2": 9,
                    "max_distance": 8.0,
                    "min_distance": 2.0,
                    "atom1": "C4",
                    "atom2": "C4",
                },
                {
                    "entity1": 2,
                    "copy1": 1,
                    "position1": 2,
                    "entity2": 2,
                    "copy2": 1,
                    "position2": 8,
                    "max_distance": 8.0,
                    "min_distance": 2.0,
                    "atom1": "C4",
                    "atom2": "C4",
                },
            ]
        }
    }
