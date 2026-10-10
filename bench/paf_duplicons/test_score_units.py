"""Tests for score_units.py (docs/PREREG_paf_duplicon_units_2026-10-10.md)."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import score_units as S  # noqa: E402


def test_cover_reads_expansions_and_unique():
    rows = [("F", "a", "F.x1"), ("F", "b", "F.x1,F.x2"), ("F", "c", "F.x2"), ("F", "d", "F.unique"), ("G", "e", "G.x1")]
    cov = S.covers(rows)
    assert cov["F"] == {"F.x1": {"a", "b"}, "F.x2": {"b", "c"}, "unique": {"d"}}
    assert cov["G"] == {"G.x1": {"e"}}


def test_family_score_identical_and_relabelled_is_one():
    d = {"F.x1": {"a", "b"}, "unique": {"c"}}
    p = {"F.x9": {"a", "b"}, "unique": {"c"}}
    assert S.family_score(d, p) == 1.0


def test_family_score_penalises_unmatched_blocks():
    d = {"x1": {"a", "b", "c", "d"}}
    p = {"y1": {"a", "b"}, "y2": {"c", "d"}}
    assert abs(S.family_score(d, p) - 0.25) < 1e-9          # best match 0.5, over max(1, 2) blocks


def test_unique_only_matches_unique():
    d = {"unique": {"a", "b"}}
    p = {"y1": {"a", "b"}}
    assert S.family_score(d, p) == 0.0


def test_b1_is_one_expansion_holding_every_member():
    d = {"x1": {"a", "b"}, "x2": {"b", "c"}, "unique": {"d"}}
    assert S.b1(d) == {"all": {"a", "b", "c", "d"}}


def test_structured_means_more_than_one_block():
    assert S.structured({"x1": {"a", "b"}, "unique": {"c"}})
    assert not S.structured({"x1": {"a", "b", "c"}})
    assert not S.structured({"unique": {"a", "b"}})


def test_half_is_a_stable_hash_split():
    import hashlib
    for f in ("MCL9", "MCL28", "MCL286"):
        exp = "dev" if int(hashlib.sha256(f.encode()).hexdigest(), 16) % 2 == 0 else "heldout"
        assert S.half(f) == exp


def test_exact_sensitivity_precision_ignore_unique():
    d = {"F": {"x1": {"a", "b"}, "x2": {"c"}, "unique": {"d"}}}
    p = {"F": {"y1": {"a", "b"}, "y2": {"c", "d"}}}
    sens, prec = S.exact_sens_prec(d, p, ["F"])
    assert sens == (1, 2) and prec == (1, 2)


def test_sign_test_one_sided():
    assert abs(S.sign_test(5, 0) - 1 / 32) < 1e-12
    assert abs(S.sign_test(1, 1) - 0.75) < 1e-12
    assert S.sign_test(0, 0) == 1.0


def test_covered_bp_counts_overlap_once():
    bed = S.merge([(0, 100), (50, 150), (300, 400)])
    assert bed == [[0, 150], [300, 400]]
    assert S.covered_bp([[100, 350]], bed) == 50 + 50
