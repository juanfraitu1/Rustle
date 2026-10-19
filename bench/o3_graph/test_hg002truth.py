import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_hg002truth", os.path.join(HERE, "hg002truth.py"))
H = importlib.util.module_from_spec(spec)
spec.loader.exec_module(H)

IDX = {"c1": [(100, 300, "F.0"), (500, 700, "F.0"), (5100, 5300, "F.1"), (5500, 5700, "F.1")], "c2": [(100, 400, "G.0")]}


def test_a_read_is_tied_when_its_best_alignment_score_is_shared_by_two_copies():
    alns = [("c1", [(150, 300), (500, 650)], 900), ("c1", [(5150, 5300), (5500, 5650)], 900), ("c2", [(100, 300)], 400)]
    assert H.tied_copies(alns, IDX) == (900, {"F.0", "F.1"})
    alns[1] = ("c1", [(5150, 5300), (5500, 5650)], 899)                  # one point lower: a unique best
    assert H.tied_copies(alns, IDX) == (900, {"F.0"})
    assert H.tied_copies([("c1", [(310, 490)], 900)], IDX) == (900, set())   # on no copy's exons


def test_the_donor_s_own_genome_settles_a_tie_only_when_its_best_placements_lift_to_one_copy():
    lift = {("chr1_MATERNAL", 1000): ("c1", 150), ("chr1_PATERNAL", 2000): ("c1", 160),   # both alleles of copy F.0
            ("chr1_MATERNAL", 9000): ("c1", 5150), ("chr7_MATERNAL", 500): None}
    lf = lambda contig, s, e: [(lift[(contig, s)][0], lift[(contig, s)][1], lift[(contig, s)][1] + (e - s))] if lift.get((contig, s)) else []  # noqa: E731
    own = [("chr1_MATERNAL", [(1000, 1150)], 950), ("chr1_PATERNAL", [(2000, 2150)], 950), ("chr1_MATERNAL", [(9000, 9150)], 940)]
    assert H.own_truth(own, lf, IDX) == ("unique", "F.0")                  # one copy, both haplotypes
    own[2] = ("chr1_MATERNAL", [(9000, 9150)], 950)
    assert H.own_truth(own, lf, IDX) == ("own_tie", None)                  # the donor's copies tie too: undecidable
    assert H.own_truth([("chr7_MATERNAL", [(500, 650)], 950)], lf, IDX) == ("no_counterpart", None)   # a donor copy CHM13 lacks
    assert H.own_truth([], lf, IDX) == ("unmapped", None)
