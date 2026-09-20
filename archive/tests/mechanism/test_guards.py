import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism"))
from v2_v3_guards import classify_guard, soto_family_sizes, measure_v2


def test_guard_inert_when_observed_below():
    assert classify_guard(guard=6000, observed_max=812)["fires"] is False


def test_guard_fires_when_observed_meets():
    r = classify_guard(guard=30, observed_max=44)
    assert r["fires"] is True


def test_classify_guard_shape():
    r = classify_guard(guard=30, observed_max=25)
    assert r == {"guard": 30, "observed_max": 25, "fires": False}


def test_soto_family_sizes_nonempty_and_matches_known_max():
    sizes = soto_family_sizes()
    assert len(sizes) == 66
    # Real data: max n_copies across the 66 Soto families is 25 (see brief).
    assert max(sizes) == 25


def test_measure_v2_both_guards_inert_on_real_data():
    v2 = measure_v2()
    assert set(v2.keys()) == {"MAX_MEMBERS(30)", "MAX_LOCI(60)"}
    for k, v in v2.items():
        assert v["observed_max"] == 25
        assert v["fires"] is False
