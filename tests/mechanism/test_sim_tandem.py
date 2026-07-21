import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism"))
from sim_tandem import build_tandem, set_match

def test_build_tandem_plants_divergence():
    ref, copies = build_tandem("ACGTACGTACGTACGT"*20, n_copies=3, divergence=0.02, seed_rng=7)
    assert len(copies) == 3
    # divergent copies must differ pairwise (at least one PSV)
    assert copies[0] != copies[1] and copies[1] != copies[2]

def test_build_tandem_identical_when_zero_divergence():
    ref, copies = build_tandem("ACGT"*100, n_copies=3, divergence=0.0, seed_rng=7)
    assert copies[0] == copies[1] == copies[2]

def test_set_match_is_permutation_invariant():
    truth = ["AAAA", "CCCC", "GGGG"]
    recovered = ["GGGG", "AAAA", "CCCC"]   # reordered
    m = set_match(recovered, truth)
    assert len(m) == 3
    assert all(identity == 1.0 for _, _, identity in m)
