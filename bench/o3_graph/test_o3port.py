"""Tests for o3port.py's pure parts: the significance test, informative columns, the zero-disagreement clustering with
refinement, the conflicting set, the slot matching and the consensus."""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import o3port as P  # noqa: E402


def obs_from(haps, counts, q=30, cols=None):
    """Reads drawn from haplotypes: haps = {name: {col: base}}, counts = {name: n}; every read covers every column of its
    haplotype unless `cols` restricts a read range."""
    obs = {}
    reads = set()
    for name, h in haps.items():
        for i in range(counts[name]):
            r = f"{name}{i}"
            reads.add(r)
            for c, b in h.items():
                if cols and c not in cols.get(r, h):
                    continue
                obs.setdefault(c, []).append((r, b, q))
    return obs, reads


def test_poisson_tail_small_and_large():
    assert abs(P.poisson_tail(1, 0.5) - (1 - 2.718281828 ** -0.5)) < 1e-6
    assert 0 < P.poisson_tail(5, 1e-6) < 1e-28            # leading-term branch, no underflow to 0 too early
    assert P.poisson_tail(400, 1e-3) == 0.0 or P.poisson_tail(400, 1e-3) < 1e-300   # log space: no OverflowError


def test_informative_needs_two_carriers_and_significance():
    # 10 A / 10 G at Q30 -> informative; a single G -> not; 2 G among 1,000 reads at Q10 -> explained by error
    obs = {1: [(f"r{i}", "A", 30) for i in range(10)] + [(f"s{i}", "G", 30) for i in range(10)],
           2: [(f"r{i}", "A", 30) for i in range(19)] + [("s0", "G", 30)],
           3: [(f"r{i}", "A", 10) for i in range(998)] + [("x", "G", 10), ("y", "G", 10)]}
    cols = P.informative(obs)
    assert cols == {1: {"A", "G"}}


def test_three_haplotypes_on_one_pile_including_a_minor_one():
    # two alleles of the copy (A/B, 45 reads each) and an extra copy at 9% (C, 9 reads): C shares A's allele at column 1
    # and B's at column 2, and is separated only by the combination plus its private column 3
    haps = {"A": {1: "A", 2: "C", 3: "G", 4: "A"}, "B": {1: "T", 2: "T", 3: "G", 4: "A"}, "C": {1: "A", 2: "T", 3: "C", 4: "G"}}
    obs, reads = obs_from(haps, {"A": 45, "B": 45, "C": 9})
    clusters, cols, top = P.cluster_pile(reads, obs)
    got = sorted(sorted({r[0] for r in c}) for c in clusters)
    assert got == [["A"], ["B"], ["C"]]
    hs = [P.haplotype(c, obs, top) for c in clusters]
    assert len(P.conflicting_set(hs)) == 3


def test_bridge_reads_do_not_merge_haplotypes():
    # reads covering only column 1 (where A and C agree) touch both; zero disagreement + label propagation keeps A and C apart
    haps = {"A": {1: "A", 2: "C"}, "C": {1: "A", 2: "T"}}
    obs, reads = obs_from(haps, {"A": 20, "C": 20})
    for i in range(6):
        obs[1].append((f"bridge{i}", "A", 30))
        reads.add(f"bridge{i}")
    clusters, _, _ = P.cluster_pile(reads, obs)
    big = [c for c in clusters if len(c) >= 10]
    assert len(big) == 2
    assert all(len({r[0] for r in c if not r.startswith("bridge")}) == 1 for c in big)


def test_fragments_of_one_haplotype_do_not_conflict():
    # one haplotype read in two non-overlapping halves: two clusters, but no shared column -> one slot suffices
    hs = [{1: "A", 2: "C"}, {7: "G", 8: "T"}]
    assert len(P.conflicting_set(hs)) == 1


def test_b_matching_leaves_the_last_in_order_unmatched():
    options = {0: ["c1"], 1: ["c1"], 2: ["c1"]}
    m = P.b_matching(options, {"c1": 2}, order=[0, 1, 2])
    assert m[0] == "c1" and m[1] == "c1" and m[2] is None
    # augmenting: haplotype 0 can move to c2, so 2 fits in c1
    options = {0: ["c1", "c2"], 1: ["c1"], 2: ["c1"]}
    m = P.b_matching(options, {"c1": 2, "c2": 2}, order=[0, 1, 2])
    assert all(v is not None for v in m.values())


def test_consensus_majority_with_deletion_and_insertion():
    al = {"r1": [(10, "A", ""), (11, "C", "GG"), (12, "-", "")],
          "r2": [(10, "A", ""), (11, "C", "GG"), (12, "-", "")],
          "r3": [(10, "T", ""), (11, "C", ""), (12, "G", "")]}
    assert P.consensus(al) == "ACGG"


def test_divergence():
    assert P.divergence({1: "A", 2: "G"}, {1: "A", 2: "T"}) == 0.5
    assert P.divergence({}, {}) == 0.0


def test_a_single_differing_column_is_not_a_conflict():
    # a piece that differs from a big haplotype at one column only (one observation) can share its slot
    assert len(P.conflicting_set([{1: "A", 2: "C", 3: "G"}, {1: "A", 2: "C", 3: "T"}])) == 1
    assert len(P.conflicting_set([{1: "A", 2: "C", 3: "G"}, {1: "A", 2: "T", 3: "T"}])) == 2


def test_private_alleles_mosaic_versus_third_pattern():
    merged = [{"hap": {"c": {1: "A", 2: "C", 3: "G", 4: "T"}}},     # slot holder A
              {"hap": {"c": {1: "T", 2: "T", 3: "C", 4: "A"}}},     # slot holder B
              {"hap": {"c": {1: "A", 2: "C", 3: "C", 4: "A"}}},     # mosaic of A and B: no private allele
              {"hap": {"c": {1: "G", 2: "C", 3: "G", 4: "G"}}}]     # third pattern: private G at 1 and 4
    match = {0: "c", 1: "c", 2: None, 3: None}
    options = {0: ["c"], 1: ["c"], 2: ["c"], 3: ["c"]}
    assert P.private_alleles(2, merged, match, options) == 0
    assert P.private_alleles(3, merged, match, options) == 2
    assert P.private_detail(3, merged, match, options) == ("c", [(1, "G"), (4, "G")])


def test_explanation_dada2_rule():
    hp = {1: "A", 2: "C", 3: "G", 4: "T"}
    left = ("p1", {1: "A", 2: "C", 3: "A", 4: "A"})
    right = ("p2", {1: "G", 2: "G", 3: "G", 4: "T"})
    assert P.explanation(hp, [left, right]) == ("join", ("p1", "p2"))
    assert P.explanation(hp, [("same", dict(hp))]) == ("single", "same")
    assert P.explanation(hp, [left]) is None                      # one parent cannot make a join
    assert P.explanation(hp, [("x", {1: "A", 4: "T"})]) is None    # uncovered columns are mismatches


def test_isocon_keep_rules():
    members = {"a", "b", "c", "d"}
    assert P.isocon_keep(members, 0, {"a", "b"}, "no-pile")           # half of the members in no phased pile
    assert not P.isocon_keep(members, 0, {"a", "b", "c"}, "no-pile")
    assert P.isocon_keep(members, 2, {"a", "b", "c", "d"}, "via-b")   # v4: half joined by attribution
    assert not P.isocon_keep(members, 1, set(), "via-b")


def test_absorb_joins_a_fragment_to_its_one_compatible_haplotype():
    hs = [{1: "A", 2: "C", 3: "G"},          # kept
          {1: "T", 2: "T", 3: "C"},          # kept
          {2: "C", 3: "G"},                  # fragment of 0
          {1: "A", 2: "T"},                  # one difference from each kept haplotype: fits both -> stays out
          {9: "G"}]                          # shares no column -> stays out
    assert P.absorb(hs, [0, 1]) == {2: 0}
