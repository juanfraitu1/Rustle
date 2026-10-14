import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import hier_match as H  # noqa: E402
import posthoc_checks as P  # noqa: E402

NODES = {
    ("family", "F1"): frozenset("ab"),
    ("expansion", "F1.x1"): frozenset("ab"),      # a twin of its family
    ("superfamily", "S1"): frozenset("abc"),
    ("sd_complex", "C1"): frozenset("abc"),       # a twin of its superfamily
    ("family", "F2"): frozenset("de"),
}
MEMBERS = sorted("abcde")


def test_distinct_set_columns_collapse_twins_onto_the_lowest_rank_node():
    assert P.distinct_nodes(NODES) == {("family", "F1"): frozenset("ab"), ("superfamily", "S1"): frozenset("abc"),
                                       ("family", "F2"): frozenset("de")}


def test_duplicate_truth_families_lose_the_twin_match_over_distinct_sets():
    truth = {"T1": frozenset("ab"), "T2": frozenset("ab")}  # two truth families with one member set
    as_scored = {m.label: m for m in H.Scorer(NODES, MEMBERS).run(truth)}
    assert [(m.node_id, m.jaccard) for m in as_scored.values()] == [("F1", 1.0), ("F1.x1", 1.0)]
    distinct = H.Scorer(P.distinct_nodes(NODES), MEMBERS).run(truth)
    # one set, one column: the second twin falls to the next set (which twin gets which is scipy's choice)
    assert sorted((m.node_id, m.jaccard) for m in distinct) == [("F1", 1.0), ("S1", 2 / 3)]
    assert P.duplicate_truth_sets(truth) == [["T1", "T2"]]


def test_shared_set_matches_count_the_matches_beyond_the_first_and_separate_twins():
    nodes = {("family", "F1"): frozenset("ab"), ("expansion", "F1.x1"): frozenset("ab"),
             ("family", "F2"): frozenset("de"), ("expansion", "F2.x1"): frozenset("de"),
             ("superfamily", "S2"): frozenset("de")}
    truth = {"T1": frozenset("ab"), "T2": frozenset("ab"),                    # twins of one Soto set: 1 extra
             "T3": frozenset("de"), "T4": frozenset("d"), "T5": frozenset("e")}  # three different families: 2 extra
    ms = H.Scorer(nodes, sorted("abde")).run(truth)
    assert {m.label: m.node_id for m in ms} == {"T1": "F1", "T2": "F1.x1", "T3": "F2", "T4": "F2.x1", "T5": "S2"}
    assert P.shared_sets(ms, nodes, truth) == {"shared_sets": 2, "matches_on_shared_sets": 5, "extra_matches": 3,
                                               "extra_matches_different_truth": 2, "extra_matches_twins": 1}


def test_strata_pool_the_tail():
    assert [P.stratum(s) for s in (2, 3, 4, 5, 43)] == ["2", "3", "4", "5+", "5+"]
    assert [P.stratum(s, top=None) for s in (2, 5, 43)] == ["2", "5", "43"]


def test_reweighting_is_the_dev_stratum_rates_under_the_held_out_mix():
    dev = {"2": (3, 4), "3": (1, 4), "5+": (0, 2)}  # stratum: (exact, families)
    held_out = {"2": 6, "3": 2, "5+": 2}             # stratum: families
    # (6 x 3/4 + 2 x 1/4 + 2 x 0) / 10 = (4.5 + 0.5) / 10
    assert P.reweight(dev, held_out) == pytest.approx(0.5)
    assert P.reweight(dev, {"2": 1, "4": 1}) is None  # a held-out stratum the dev catalog lacks: undefined


def test_rates_by_stratum_count_exact_and_evaluable_per_stratum():
    sizes = {"T1": 2, "T2": 2, "T3": 3, "T4": 7}
    exact = {"T1": True, "T2": False, "T3": True, "T4": False}
    assert P.rates_by(sizes, exact, P.stratum) == {"2": (1, 2), "3": (1, 1), "5+": (0, 1)}


def test_near_end_sides_keep_same_length_sides_within_the_window_of_a_v1_end():
    lengths = {"chr1": 10_000}
    regions = [
        (("chr1", 100, 300), ("chr1", 100, 300)),        # near the start: kept
        (("chr1", 9_500, 9_900), ("chr1", 9_400, 9_800)),  # near the end (v1 coordinates): kept
        (("chr1", 5_000, 5_100), ("chr1", 5_000, 5_100)),  # mid-chromosome: dropped
        (("chr1", 0, 300), ("chr1", 0, 290)),            # spans an indel (lengths differ): dropped
    ]
    assert P.near_end_sides(regions, lengths, 1000) == regions[:2]


@pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools not on PATH")
def test_differing_sides_compare_v1_and_v2_sequence(tmp_path):
    v1, v2 = tmp_path / "v1.fa", tmp_path / "v2.fa"
    v1.write_text(">chr1\nACGTACGTAC\n")
    v2.write_text(">chr1\nACGTTCGTAC\n")  # one substitution at 0-based 4
    for fa in (v1, v2):
        subprocess.run(["samtools", "faidx", str(fa)], check=True)
    sides = [(("chr1", 0, 4), ("chr1", 0, 4)), (("chr1", 2, 6), ("chr1", 2, 6))]
    assert P.differing_sides(str(v1), str(v2), sides) == [sides[1]]
