import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2correct", os.path.join(HERE, "o2correct.py"))
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)


def test_the_corrector_takes_the_tied_copy_with_fewest_substitutions_and_keeps_minimap2_s_pick_on_a_tie():
    per = {"F.0": dict(AS=1000, sub=3), "F.1": dict(AS=995, sub=1)}        # minimap2 took F.0 for its indel-free AS
    assert C.pick(per, "F.0") == ("F.1", "fewest_substitutions")
    assert C.pick({"F.0": dict(AS=1000, sub=1), "F.1": dict(AS=995, sub=1)}, "F.0") == ("F.0", "minimap2")
    assert C.pick({"F.0": dict(AS=1000, sub=0), "F.1": dict(AS=995, sub=2)}, "F.0") == ("F.0", "fewest_substitutions")


def test_a_read_s_tied_copies_are_its_best_alignment_per_copy_within_ten_as_of_its_best():
    alns = [dict(AS=1000, copy="F.0", sub=3, prim=True), dict(AS=995, copy="F.1", sub=1, prim=False),
            dict(AS=990, copy="F.1", sub=0, prim=False), dict(AS=980, copy="F.2", sub=0, prim=False)]
    per, prim = C.tied(alns)
    assert prim == "F.0" and per == {"F.0": dict(AS=1000, copy="F.0", sub=3, prim=True),
                                      "F.1": dict(AS=995, copy="F.1", sub=1, prim=False)}
    assert C.tied(alns + [dict(AS=999, copy=None, sub=0, prim=False)]) == (None, "F.0")      # tied off the family
    assert C.tied(alns[:1]) == (None, "F.0")                                                # one copy: minimap2's


def test_fixes_and_breaks_count_where_the_corrector_and_minimap2_disagree_against_the_truth():
    n = C.tally([("F.1", "F.0", "F.1"), ("F.1", "F.0", "F.0"), ("F.0", "F.0", "F.0"), ("F.2", "F.0", "F.1")])
    assert (n["reads"], n["corrector right"], n["minimap2 right"], n["fixes"], n["breaks"], n["both wrong"]) == (4, 2, 2, 1, 1, 1)


def _fetch(ref):
    return lambda lo, hi: ref[max(0, lo):hi]


def test_an_alignment_s_events_are_its_substitutions_and_left_aligned_indels_and_its_blocks_skip_introns():
    ev, bl = C.alignment_events(0, [(0, 4)], "AGGT", _fetch("ACGT"))
    assert ev == {("S", 1, "G")} and bl == [(0, 4)]
    ev, _ = C.alignment_events(0, [(0, 4), (2, 1), (0, 4)], "ACGTTGCA", _fetch("ACGTTTGCA"))     # a T dropped from TTT
    assert ev == {("D", 3, 1)}
    ev, _ = C.alignment_events(0, [(0, 5), (1, 1), (0, 3)], "ACGTTTGCA", _fetch("ACGTTGCA"))     # a T added to TT
    assert ev == {("I", 3, "T")}
    ev, bl = C.alignment_events(10, [(4, 2), (0, 2), (3, 3), (0, 2)], "NNACGT", _fetch("x" * 10 + "ACgggGT"))
    assert ev == set() and bl == [(10, 12), (15, 17)]


def test_the_sample_explains_an_event_two_or_more_of_its_confident_reads_share_and_holds_the_ones_nearly_all_carry():
    d = C.Donor()
    for evs in ({("I", 5, "T")}, {("I", 5, "T")}, {("I", 5, "T"), ("S", 7, "A")}, set()):
        d.add("F.0", evs, [(0, 20)])
    d.add("F.1", set(), [(0, 20)])
    d.add("F.1", set(), [(0, 4)])
    assert d.depth("F.0", 5) == 4 and d.depth("F.1", 5) == 1 and d.depth("F.0", 20) == 0
    assert d.explained("F.0", ("I", 5, "T")) and not d.explained("F.0", ("S", 7, "A"))        # one read: an error
    assert not d.explained("F.1", ("I", 5, "T"))
    d.add("F.0", {("I", 5, "T")}, [(0, 20)])
    assert d.hom("F.0") == [("I", 5, "T")]                                                   # 4 of 5 reads


def test_the_sample_aware_corrector_takes_the_copy_whose_own_reads_explain_the_read_and_keeps_minimap2_on_a_tie():
    d = C.Donor()
    for _ in range(4):
        d.add("F.0", {("I", 5, "T")}, [(0, 20)])                                              # the donor's F.0 has the T
        d.add("F.1", set(), [(0, 20)])
    per = {"F.0": ({("I", 5, "T")}, [(0, 20)]), "F.1": ({("S", 9, "C")}, [(0, 20)])}         # minimap2 took F.1 (AS: 5 < 8)
    assert C.donor_pick(per, "F.1", d) == ("F.0", "sample")
    per = {"F.0": ({("I", 5, "T"), ("S", 2, "A")}, [(0, 20)]), "F.1": ({("S", 9, "C")}, [(0, 20)])}
    assert C.donor_pick(per, "F.1", d) == ("F.1", "tie")                                       # one unexplained each
    per = {"F.0": (set(), [(0, 20)]), "F.1": ({("S", 9, "C")}, [(0, 20)])}                    # the read lacks F.0's T
    assert C.donor_pick(per, "F.1", d) == ("F.1", "tie")


def test_a_fragment_s_own_source_read_is_left_out_of_the_sample_that_judges_it():
    d = C.Donor()
    src = ("F.0", {("I", 5, "T")}, [(0, 20)])
    d.add(*src)
    d.add("F.0", {("I", 5, "T")}, [(0, 20)])
    d.add("F.0", set(), [(0, 20)])
    d.add("F.0", set(), [(0, 20)])
    assert d.explained("F.0", ("I", 5, "T")) and not d.explained("F.0", ("I", 5, "T"), own=src)   # 2 of 4; 1 of 3
    assert d.depth("F.0", 5, own=src) == 3 and d.count("F.0", ("I", 5, "T"), own=src) == 1
    h = C.Donor()
    h.add(*src)
    for _ in range(4):
        h.add("F.0", {("I", 5, "T")}, [(0, 20)])
    h.add("F.0", set(), [(0, 20)])                                                             # 5 of 6 carry it
    assert h.hom("F.0") == [("I", 5, "T")] and h.hom("F.0", own=src) == [("I", 5, "T")]        # 4 of 5
    h.add("F.0", set(), [(0, 20)])                                                             # 5 of 7; 5 of 6 without one
    assert h.hom("F.0") == [] and h.hom("F.0", own=("F.0", set(), [(0, 20)])) == [("I", 5, "T")]


def test_minimap2_s_copy_is_overruled_only_when_the_alternative_wins_with_every_uncovered_event_granted_to_minimap2():
    d = C.Donor()
    for _ in range(4):
        d.add("F.0", {("I", 5, "T")}, [(0, 20)])
        d.add("F.1", set(), [(0, 10)])                                                            # F.1's reads stop at 10
    covered = {"F.0": ({("I", 5, "T")}, [(0, 20)]), "F.1": ({("S", 9, "C")}, [(0, 20)])}
    assert C.donor_pick(covered, "F.1", d) == ("F.0", "sample")                                  # F.1's 9C refuted
    uncovered = {"F.0": ({("I", 5, "T")}, [(0, 20)]), "F.1": ({("S", 15, "C")}, [(0, 20)])}
    assert C.donor_pick(uncovered, "F.1", d) == ("F.1", "undominated")                           # F.1's 15C unknown
    assert C.donor_pick(uncovered, "F.0", d) == ("F.0", "sample")


def test_a_read_s_aligned_span_and_its_events_within_a_read_window():
    assert C.read_span([(4, 2), (0, 2), (3, 3), (0, 2), (4, 1)]) == (2, 6)
    ev, bl = C.alignment_events(0, [(0, 8)], "AGGTACGA", _fetch("ACGTACGT"), window=(0, 4))
    assert ev == {("S", 1, "G")} and bl == [(0, 4)]
    ev, bl = C.alignment_events(0, [(0, 3), (3, 4), (0, 5)], "AGGTACGA", _fetch("AGGxxxxTACGT"), window=(2, 8))
    assert ev == {("S", 11, "A")} and bl == [(2, 3), (7, 12)]


def test_tied_alignments_are_compared_over_the_read_bases_every_one_of_them_aligns():
    # the read's first 3 bases (an error, 1D) align to F.0 only: F.1 soft-clips them, so they cannot count against F.0
    assert C.common_window({"F.0": ((0, 10), False), "F.1": ((3, 10), False)}, 10) == {"F.0": (3, 10), "F.1": (3, 10)}
    # F.1 aligns the reverse strand: its read coordinates run from the read's other end
    assert C.common_window({"F.0": ((0, 8), False), "F.1": ((0, 8), True)}, 10) == {"F.0": (2, 8), "F.1": (2, 8)}


def test_the_copies_an_alignment_can_lie_on_are_found_by_span_in_their_table_order():
    cands = [("F.2", "+", [(500, 600)], (500, 900)), ("F.0", "+", [(0, 100)], (0, 300)), ("F.1", "-", [(250, 280)], (250, 260))]
    idx = C.span_index({"chr1": cands})
    assert C.overlapping(idx, "chr1", 255, 520) == [cands[0], cands[1], cands[2]]                 # table order kept
    assert C.overlapping(idx, "chr1", 300, 500) == []
    assert C.overlapping(idx, "chr1", 899, 2000) == [cands[0]] and C.overlapping(idx, "chr2", 0, 10) == []
