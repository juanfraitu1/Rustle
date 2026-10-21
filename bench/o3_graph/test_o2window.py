import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("o3_graph_o2window", os.path.join(HERE, "o2window.py"))
OW = importlib.util.module_from_spec(spec)
spec.loader.exec_module(OW)


def test_an_alignment_is_stitched_when_it_lays_bases_on_two_copies_exons():
    exons = {"F.0": [(0, 2500)], "F.1": [(14000, 15000)]}
    assert OW.stitched([(100, 2100), (14320, 14795)], exons) is True
    assert OW.stitched([(100, 2100)], exons) is False
    assert OW.stitched([(100, 2100), (14000, 14040)], exons) is False                # 40 bp on the second copy


def _a(ms, copy, prim=False):
    return dict(ms=ms, copy=copy, prim=prim)


def test_the_window_holds_the_copies_within_w_of_the_read_s_best_intron_free_score():
    al = [_a(1000, "F.0", True), _a(995, "F.1"), _a(980, "F.2")]
    assert OW.window_copies(al, 10) == ({"F.0", "F.1"}, "F.0")
    assert OW.window_copies(al, 30) == ({"F.0", "F.1", "F.2"}, "F.0")
    assert OW.window_copies(al, 5) == (None, "F.0")                                   # one copy is no tie
    assert OW.window_copies(al + [_a(985, None)], 20) == (None, "F.0")                # an off-copy locus inside the window
    assert OW.window_copies([_a(1000, "F.0", True), _a(995, "G.1")], 10) == (None, "F.0")   # two families
    assert OW.window_copies([_a(1000, "F.0"), _a(999, "F.1"), _a(900, "F.2", True)], 10) == (None, "F.2")   # primary outside


def test_the_set_drops_only_a_copy_another_beats_at_every_column_where_they_differ():
    win = (0, 100)
    clean = {"a": ({20, 30}, {5}), "b": (set(), {5, 20, 30})}                       # the read is a's: b loses at 20 and 30
    assert OW.pareto_set(clean, win) == {"a"}
    mosaic = {"a": ({30}, {20}), "b": ({20}, {30})}                                  # each wins one column
    assert OW.pareto_set(mosaic, win) == {"a", "b"}
    same = {"a": ({1, 2}, set()), "b": ({1, 2}, set())}                              # identical over the read
    assert OW.pareto_set(same, win) == {"a", "b"}
    assert OW.pareto_set(clean, (0, 10)) == {"a", "b"}                               # the PSVs fall outside the window


def test_the_ms_truth_is_the_one_copy_within_one_substitution_of_the_best_intron_free_score():
    ref, alt = "ref", "alt"
    assert OW.truth_ms([(214, "a", ref), (208, "b", ref)]) == ("determined", "a", "refcons")     # 0 mismatches + 4-bp clip
    assert OW.truth_ms([(214, "a", alt), (209, "a", ref), (200, "b", ref)]) == ("determined", "a", "altonly")
    assert OW.truth_ms([(214, "a", ref), (210, "b", alt)]) == ("undetermined", None, None)      # under one substitution apart
    assert OW.truth_ms([(214, None, alt), (200, "a", ref)]) == ("off-copy", None, None)
    assert OW.truth_ms([]) == ("unmapped", None, None)
