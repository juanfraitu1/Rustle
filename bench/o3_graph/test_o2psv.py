import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2psv", os.path.join(HERE, "o2psv.py"))
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)


def test_an_alignment_s_substitutions_come_from_x_ops_or_from_nm_less_its_indel_bases():
    assert P.substitutions([(7, 100), (8, 2), (1, 3), (7, 50)], nm=None) == (2, 152)        # =/X: X counted, 152 aligned
    assert P.substitutions([(0, 150), (2, 4), (0, 10)], nm=7) == (3, 160)                  # M: NM 7 less 4 deleted bases


def test_a_read_is_placed_when_two_psvs_separate_its_copy_and_that_copy_explains_it():
    D = P.psv_decide
    assert D({"F.0": (0, 2000), "F.1": (2, 2000)}) == ("F.0", "psv_certificate", ["F.0"])    # two PSVs, a perfect fit
    assert D({"F.0": (0, 2000), "F.1": (1, 2000)}) == ("-", "psv_set", ["F.0", "F.1"])       # one PSV: a donor allele can do that
    assert D({"F.0": (3, 2000), "F.1": (5, 2000)}) == ("-", "psv_unexplained", ["F.0", "F.1"])   # 0.15% > the HiFi floor
    assert D({"F.0": (2, 2000), "F.1": (4, 2000), "F.2": (2, 2000)}) == ("-", "psv_set", ["F.0", "F.2"])
    assert D({"F.0": (1, 1000), "F.1": (3, 1000), "F.2": (9, 1000)}) == ("F.0", "psv_certificate", ["F.0"])


def test_the_tied_copies_are_the_alignments_within_ten_as_of_the_best_one_per_copy_off_family_ones_marked():
    alns = [dict(AS=1000, copy="F.0", mism=0, alen=900, prim=True), dict(AS=995, copy="F.1", mism=1, alen=900, prim=False),
            dict(AS=995, copy="F.1", mism=3, alen=880, prim=False), dict(AS=980, copy="F.2", mism=4, alen=900, prim=False)]
    assert P.tied_copies(alns) == {"F.0": (0, 900), "F.1": (1, 900)}                         # F.2 is 20 AS down
    assert P.tied_copies(alns + [dict(AS=999, copy=None, mism=0, alen=900, prim=False)]) is None   # tied off the family
    assert P.tied_copies(alns[:1]) == {"F.0": (0, 900)}


def test_fates_judge_a_call_or_a_set_against_a_truth_copy_and_against_a_donor_class():
    F = P.fate
    assert F("F.0", ["F.0"], "F.0") == ("right", True)
    assert F("-", ["F.0", "F.1"], "F.1") == ("abstain", True)
    assert F("F.1", ["F.1"], "F.0") == ("wrong", False)
    assert F("F.1", ["F.1"], None, cls="own_tie") == ("unsupported", None)                  # the donor ties: no truth copy
    assert F("F.1", ["F.1"], None, cls="elsewhere") == ("wrong", False)                     # its copy is off the family
    assert F("F.1", ["F.1"], None, cls="poor_match") == ("unscored", None)                  # the donor fits it badly: no truth


def test_the_score_counts_o2_s_calls_and_sets_and_minimap2_s_primary_against_the_truth():
    rows = [dict(read="a", call="F.0", candidates="F.0,F.1", set="F.0", primary="F.1"),
            dict(read="b", call="-", candidates="F.0,F.1,F.2", set="F.1,F.2", primary="F.2"),
            dict(read="c", call="F.1", candidates="F.0,F.1", set="F.1", primary="F.1"),
            dict(read="z", call="F.1", candidates="F.0,F.1", set="F.1", primary="F.1")]
    n = P.score(rows, {"a": "F.0", "b": "F.2", "c": "F.0"})
    assert (n["reads"], n["right"], n["wrong"], n["abstain"], n["set holds truth"], n["judged sets"]) == (3, 1, 1, 1, 2, 3)
    assert (n["candidates"], n["set size"], n["minimap2 right"], n["minimap2 wrong"]) == (7, 4, 1, 2)
    d = P.score(rows, {"a": ("unique_elsewhere", "F.0"), "b": ("own_tie", None), "c": ("elsewhere", None)}, donor=True)
    assert (d["right"], d["wrong"], d["abstain"], d["judged sets"], d["minimap2 unsupported"]) == (1, 1, 1, 2, 1)


def test_aligned_pairs_count_soft_clips_and_a_leading_hard_clip_in_the_read_s_positions():
    assert P.aligned_pairs(100, [(4, 2), (7, 3), (1, 1), (8, 1), (3, 50), (7, 2)]) == [(2, 100), (3, 101), (4, 102),
                                                                                       (6, 103), (7, 154), (8, 155)]
    assert P.aligned_pairs(10, [(5, 4), (0, 2)]) == [(4, 10), (5, 11)]


def test_psv_columns_are_read_positions_every_tied_copy_aligns_with_differing_reference_bases():
    ref = {"c": "AAAAACAAAAGAAAA"}
    fetch = lambda ch, s, e: ref[ch][s:e]                                      # noqa: E731
    al = {"F.0": ("c", 0, [(7, 5)], False), "F.1": ("c", 5, [(7, 5)], False), "F.2": ("c", 10, [(7, 3)], False)}
    cols = P.psv_columns(5, False, al, fetch)
    assert cols == {0: {"F.0": ("c", 0, "A", False), "F.1": ("c", 5, "C", False), "F.2": ("c", 10, "G", False)}}
    rc = P.psv_columns(5, False, {"F.0": ("c", 0, [(7, 5)], True), "F.1": ("c", 5, [(7, 5)], False)}, fetch)
    assert rc[4]["F.0"] == ("c", 0, "T", True) and rc[4]["F.1"] == ("c", 9, "A", False)   # opposite strand: complemented


def test_a_copy_is_excluded_when_the_read_carries_another_copy_s_base_that_its_own_confident_reads_never_show():
    cols = {3: {"F.0": ("c", 103, "A", False), "F.1": ("c", 503, "G", False)}}
    read = "TTTAT"
    six = {("F.1", "c", 503): {f"t{i}": "G" for i in range(6)}}
    assert P.excluded(cols, read, six, k=3) == {"F.1"}                         # 6 F.1 reads, none with A
    one_a = {("F.1", "c", 503): {**{f"t{i}": "G" for i in range(5)}, "u": "A"}}
    assert P.excluded(cols, read, one_a, k=3) == set()                          # the sample's F.1 shows A too
    assert P.excluded(cols, read, one_a, k=3, own="u") == {"F.1"}               # ... only through the read itself
    assert P.excluded(cols, read, {("F.1", "c", 503): {"t1": "G", "t2": "G"}}, k=3) == set()   # too few to tell
    assert P.excluded(cols, "TTTCT", six, k=3) == set()                         # the read carries neither copy's base
    flip = {3: {"F.0": ("c", 103, "A", False), "F.1": ("c", 503, "G", True)}}
    assert P.excluded(flip, read, {("F.1", "c", 503): {f"t{i}": "C" for i in range(4)}}, k=3) == {"F.1"}   # genome strand


def test_the_validated_rule_places_the_one_copy_left_if_it_explains_the_read_else_keeps_the_set():
    V = P.svp_decide
    cands = {"F.0": (0, 2000), "F.1": (1, 2000), "F.2": (1, 2000)}
    assert V(cands, {"F.1", "F.2"}) == ("F.0", "svp_certificate", ["F.0"])
    assert V(cands, {"F.1"}) == ("-", "svp_set", ["F.0", "F.2"])
    assert V(cands, set(cands)) == ("-", "svp_contradiction", ["F.0", "F.1", "F.2"])    # every copy excluded: a mosaic
    assert V({"F.0": (5, 2000), "F.1": (6, 2000)}, {"F.1"}) == ("-", "psv_unexplained", ["F.0", "F.1"])


def test_the_fit_floor_is_the_sample_s_own_substitution_rate_quantile_among_its_uniquely_placed_reads():
    rates = [0.0] * 50 + [0.001] * 40 + [0.002] * 9 + [0.01]
    assert P.quantile_rate(rates, 0.99) == 0.002 and P.quantile_rate(rates, 0.5) == 0.0
    assert P.svp_decide({"F.0": (3, 1000), "F.1": (5, 1000)}, {"F.1"}, fit=0.004) == ("F.0", "svp_certificate", ["F.0"])
    assert P.svp_decide({"F.0": (3, 1000), "F.1": (5, 1000)}, {"F.1"}, fit=0.002) == ("-", "psv_unexplained", ["F.0", "F.1"])


def test_a_copy_the_read_s_own_bases_contradict_somewhere_is_never_placed():
    cols = {3: {"F.0": ("c", 103, "A", False), "F.1": ("c", 503, "G", False)},
            7: {"F.0": ("c", 107, "C", False), "F.1": ("c", 507, "T", False)}}
    assert P.contradicted(cols, "TTTATTTT") == {"F.0", "F.1"}       # A (F.0's) at 3, T (F.1's) at 7: a mosaic
    assert P.contradicted(cols, "TTTATTTC") == {"F.1"}
    assert P.contradicted(cols, "TTTCTTTG") == set()               # neither copy's base anywhere


def test_the_regimen_places_by_two_psvs_or_by_one_the_sample_validates_and_never_against_the_read_s_own_bases():
    R = P.regimen_decide
    two = {"F.0": (0, 2000), "F.1": (2, 2000)}
    one = {"F.0": (0, 2000), "F.1": (1, 2000)}
    assert R(two, set(), set()) == ("F.0", "psv_two", ["F.0"])                       # two PSVs: no validation needed
    assert R(one, set(), set()) == ("-", "psv_set", ["F.0", "F.1"])                  # one PSV the sample cannot vouch for
    assert R(one, {"F.1"}, set()) == ("F.0", "psv_validated", ["F.0"])               # ... or can
    assert R(one, {"F.1"}, {"F.0"}) == ("-", "psv_set", ["F.0", "F.1"])              # the read also carries F.1's base
    assert R(one, {"F.0"}, set()) == ("F.1", "psv_validated", ["F.1"])
    assert R(two, {"F.0"}, set()) == ("-", "psv_conflict", ["F.0", "F.1"])          # the two rules disagree
    assert R({"F.0": (5, 2000), "F.1": (7, 2000)}, set(), set(), fit=0.002) == ("-", "psv_unexplained", ["F.0", "F.1"])
    assert R({"F.0": (0, 900), "F.1": (1, 900), "F.2": (3, 900)}, set(), set()) == ("-", "psv_set", ["F.0", "F.1"])


def test_exclusions_shape_the_set_only_when_the_read_does_not_contradict_every_copy_they_leave():
    one = {"F.0": (0, 2000), "F.1": (1, 2000)}
    assert P.regimen_decide(one, {"F.0"}, {"F.1"}) == ("-", "psv_set", ["F.0", "F.1"])   # MAGEA9B: an error made a mosaic


def test_with_validated_only_every_placement_and_every_narrowing_needs_the_sample_s_validation():
    R = P.regimen_decide
    two = {"F.0": (0, 2000), "F.1": (2, 2000)}
    assert R(two, set(), set(), validated_only=True) == ("-", "psv_set", ["F.0", "F.1"])     # two PSVs alone: a set
    assert R(two, {"F.1"}, set(), validated_only=True) == ("F.0", "psv_validated", ["F.0"])
    three = {"F.0": (0, 900), "F.1": (1, 900), "F.2": (3, 900)}
    assert R(three, set(), set(), validated_only=True) == ("-", "psv_set", ["F.0", "F.1", "F.2"])   # no unvalidated cut
    assert R(three, {"F.2"}, set(), validated_only=True) == ("-", "psv_set", ["F.0", "F.1"])


def test_a_fragment_whose_truth_is_none_of_its_tied_copies_can_be_scored_apart():
    rows = [dict(read="a", call="F.1", candidates="F.1,F.2", set="F.1", primary="F.1"),     # truth F.0: not a candidate
            dict(read="b", call="F.0", candidates="F.0,F.1", set="F.0", primary="F.1")]
    truth = {"a": "F.0", "b": "F.0"}
    assert P.score(rows, truth)["wrong"] == 1
    n = P.score(rows, truth, among_candidates=True)
    assert (n["reads"], n["right"], n["wrong"], n["outside candidates"]) == (1, 1, 0, 1)


def test_a_fragment_s_own_source_read_is_the_read_left_out_of_its_validation():
    assert P.own_name("m64/1/ccs|f2|130|700") == "m64/1/ccs" and P.own_name("m64/1/ccs") == "m64/1/ccs"
