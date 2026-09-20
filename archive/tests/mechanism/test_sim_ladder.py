import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2] / "bench/mechanism"))
from sim_ladder import (
    _find_divergence_for_identity,
    _mean_pairwise_identity,
    _pick_host_hidden_partition,
)
from sim_tandem import build_tandem


def test_mean_pairwise_identity_all_identical_is_one():
    bodies = ["ACGTACGT"] * 4
    assert _mean_pairwise_identity(bodies) == 1.0


def test_mean_pairwise_identity_detects_divergence():
    bodies = ["AAAA", "AAAT", "AATT", "ATTT"]
    # 6 pairs, hamming distances: (0,1)=1 (0,2)=2 (0,3)=3 (1,2)=1 (1,3)=2 (2,3)=1
    # identities: 3/4 2/4 1/4 3/4 2/4 3/4 -> mean = (3+2+1+3+2+3)/(6*4) = 14/24
    assert abs(_mean_pairwise_identity(bodies) - 14 / 24) < 1e-9


def test_find_divergence_for_identity_hits_target_closely():
    seed_seq = "ACGT" * 350  # 1400bp, matches EXON_LEN*2 in sim_ladder
    for target in (0.90, 0.96, 0.995):
        d = _find_divergence_for_identity(seed_seq, 4, target, seed_rng=123)
        _, bodies = build_tandem(seed_seq, 4, d, seed_rng=123)
        achieved = _mean_pairwise_identity(bodies)
        assert abs(achieved - target) < 0.01, (target, achieved, d)


def test_find_divergence_for_identity_zero_at_target_one():
    assert _find_divergence_for_identity("ACGT" * 10, 4, 1.0, seed_rng=1) == 0.0


def test_pick_host_hidden_partition_finds_the_planted_tight_cluster():
    # A is far from everyone; host/hidden1/hidden2 are mutually close.
    a = "AAAA" * 25
    host = "CCCC" * 25
    hidden1 = host[:-1] + "G"       # 1 base off host
    hidden2 = host[:-2] + "GG"      # 2 bases off host
    bodies = [a, host, hidden1, hidden2]
    score, a_idx, host_idx, h1_idx, h2_idx = _pick_host_hidden_partition(bodies)
    assert a_idx == 0
    assert {host_idx, h1_idx, h2_idx} == {1, 2, 3}
    assert score > 0
