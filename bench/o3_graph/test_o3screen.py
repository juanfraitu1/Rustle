"""EXPLORATORY (spike, 2026-10-22): one sanity test of o3screen's tier-1 core on a toy copy."""
import importlib.util
import os

spec = importlib.util.spec_from_file_location("o3screen_t", os.path.join(os.path.dirname(os.path.abspath(__file__)), "o3screen.py"))
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)


def test_three_haplotypes_on_a_diploid_copy_exceed_capacity_and_one_site_errors_do_not():
    h = {"h1": "AAAAAA", "h2": "CCAAAA", "h3": "AACCAA"}           # pairwise 2 differences: three haplotypes
    pats = {f"{k}_{i}": {p: b for p, b in enumerate(s)} for k, s in h.items() for i in range(2)}
    pats["err"] = {p: b for p, b in enumerate("AAAAAT")}           # h1 with one error: joins h1, opens no class
    assert not M.conflict(pats["err"], pats["h1_0"]) and M.conflict(pats["h2_0"], pats["h3_0"])
    for pair in (False, True):
        cls = M.classes(pats)
        assert sorted(len(r) for r, _, _ in cls) == [2, 2, 3]
        assert len(M.greedy_clique(cls, pair)) == 3 > M.capacity("chr1", 0, 10, "male")
    assert M.capacity("chrX", 5_000_000, 5_001_000, "male") == 1 and M.capacity("chrX", 1000, 2000, "male") == 2


def test_deficiency_counts_classes_no_slot_can_host():
    assert M.deficiency([[("c", 0), ("c", 1)]] * 3) == 1
    assert M.deficiency([[("c", 0), ("c", 1)], [("c", 0), ("c", 1)], [("c", 0), ("d", 0)]]) == 0
