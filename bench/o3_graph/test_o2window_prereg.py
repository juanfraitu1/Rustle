"""PREREG_o2_window_set_2026-10-20: the generalized prepare / truth converters / scorer of o2window.py."""
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2window_p", os.path.join(HERE, "o2window.py"))
OW = importlib.util.module_from_spec(spec)
spec.loader.exec_module(OW)


def test_match_columns_compare_bases_under_an_m_op_and_trust_eqx_ops():
    ref = "AAAACCCCGGGG"
    fetch = lambda a, b: ref[a:b]
    cig = [(4, 1), (0, 4), (3, 4), (7, 2), (8, 1)]                                  # 1S 4M 4N 2= 1X
    eq, x = OW.match_columns_ref(cig, False, 8, "TAATACCG", 0, fetch)              # read AATA vs ref AAAA: X at read 3
    assert eq == {1, 2, 4, 5, 6} and x == {3, 7}
    eq, x = OW.match_columns_ref(cig, True, 8, "TAATACCG", 0, fetch)
    assert eq == {6, 5, 3, 2, 1} and x == {4, 0}


def test_a_simulated_truth_row_becomes_determined_only_when_its_copy_is_in_the_reference():
    rows = [dict(read="r1", copy="F.0", kind="present"), dict(read="r2", copy="F.9", kind="present")]
    out = OW.truth_sim(rows, {"F.0"})
    assert out == {"r1": ("determined", "F.0"), "r2": ("off-copy", "")}


def test_an_hg002_truth_row_maps_its_class_to_a_status():
    rows = [dict(read="a", **{"class": "unique"}, truth="F.0"), dict(read="b", **{"class": "own_tie"}, truth=""),
            dict(read="c", **{"class": "no_counterpart"}, truth=""), dict(read="d", **{"class": "poor_match"}, truth=""),
            dict(read="e", **{"class": "elsewhere"}, truth=""), dict(read="f", **{"class": "unique_elsewhere"}, truth="F.1")]
    out = OW.truth_hg002(rows)
    assert out["a"] == ("determined", "F.0") and out["b"] == ("undetermined", "") and out["c"] == ("off-copy", "")
    assert out["d"] == ("poor_match", "") and out["e"] == ("off-copy", "") and out["f"] == ("determined", "F.1")


def _call(read, fam, copies, mm, o2, st, w=30, stitched="0"):
    return dict(read=read, family=fam, window=str(w), copies=",".join(copies), minimap2=mm, o2=o2, o2_set=",".join(st),
                primary_stitched=stitched)


def test_the_bars_count_truth_held_singletons_and_singleton_precision():
    calls = [_call("r1", "F", ["a", "b"], "a", "a", ["a"]), _call("r2", "F", ["a", "b"], "b", "a", ["a", "b"]),
             _call("r3", "G", ["c", "d"], "c", "c", ["d"]), _call("r4", "G", ["c", "d"], "d", "d", ["d"])]
    truth = {"r1": ("determined", "a"), "r2": ("determined", "a"), "r3": ("determined", "c"), "r4": ("undetermined", "")}
    b = OW.bars(calls, truth)
    assert (b["scored"], b["families"], b["truth_held"], b["singletons"], b["singletons_right"]) == (3, 2, 2, 2, 1)
    assert b["minimap2_right"] == 2 and b["minimap2_right_on_singletons"] == 2
    assert b["H"] is False and b["R"] is True and b["Q"] is False and b["P"] is False
    assert b["verdict"] == "UNDERPOWERED"


def test_pooled_worlds_key_reads_by_world_so_repeated_simulated_names_do_not_collide():
    assert OW.world_tag("/x/heldout24/A119b_s1401/CROSS") == "A119b_s1401_CROSS"
    assert OW.calls_tag("/y/s1/A119b_s1401_CROSS.calls.tsv") == "A119b_s1401_CROSS"
    assert OW.keyed("s1", "A119b_s1401_CROSS", pooled=True) == "A119b_s1401_CROSS/s1"
    assert OW.keyed("s1", "A119b_s1401_CROSS", pooled=False) == "s1"


def test_the_band_scorer_scores_the_set_and_both_nulls_on_the_same_reads():
    calls = [_call("r1", "F", ["a", "b"], "a", "a", ["a"]), _call("r2", "F", ["a", "b"], "b", "a", ["a", "b"]),
             _call("r3", "G", ["c", "d"], "c", "c", ["d"]), _call("r4", "G", ["c", "d"], "d", "d", ["c"])]
    truth = {q: ("determined", t) for q, t in (("r1", "a"), ("r2", "a"), ("r3", "c"), ("r4", "d"))}
    b = OW.band_bars(calls, truth)
    assert b["set"]["scored"] == 4 and b["set"]["truth_held"] == 2 and b["set"]["singletons"] == 3 and b["set"]["singletons_right"] == 1
    assert b["keep_all"]["truth_held"] == 4 and b["keep_all"]["singletons"] == 0
    assert b["primary"]["singletons"] == 4 and b["primary"]["singletons_right"] == 3
    assert b["bars"] == dict(P=False, H=False, N=True, Q=False, verdict="UNDERPOWERED")


def test_window_calls_become_o2loci_calls_with_an_empty_copy_for_an_abstained_read():
    calls = [_call("r1", "F", ["F.0", "F.1"], "F.0", "F.1", ["F.1"]), _call("r2", "F", ["F.0", "F.1"], "F.0", "F.0", ["F.0", "F.1"]),
             _call("r3", "G", ["G.0", "G.1"], "G.0", "G.0", ["G.0"])]
    assert OW.loci_calls(calls, "F") == {"r1": "F.1", "r2": ""}
