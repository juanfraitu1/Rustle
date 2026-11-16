"""tools/eq_classes.py: ranges hold the EM estimate, no read is lost, groups need two shared reads, repeats link nothing."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pytest
from eq_classes import summarise, em, load_ties, floor_fn, OUTSIDE

META = {l: ["c1", i * 100, i * 100 + 50, "+", "."] for i, l in enumerate("ABCDE")}

def test_ranges_and_groups():
    placed = {"A": 90, "B": 10, "E": 5}
    classes = {frozenset("AB"): 100, frozenset("CD"): 8, frozenset(["D", "E"]): 1, frozenset(["C", OUTSIDE]): 4}
    loci, cls = summarise(placed, classes, META)
    by = {r[0]: r for r in loci}
    assert (by["A"][6], by["A"][7], by["A"][8]) == ("mixed", 90, 190) and by["C"][6] == "tied_only" and by["E"][7:9] == [5, 6]
    assert abs(float(by["A"][9]) - 180) < 0.5 and abs(float(by["B"][9]) - 20) < 0.5  # EM follows the placed reads 9:1
    assert by["A"][10] == by["B"][10] != "." and by["C"][10] == by["D"][10] != by["A"][10]
    assert by["E"][10] == "."            # one shared read links nothing
    assert by["C"][13] == 0              # the C/D group has no placed read: its EM is the starting split
    assert sum(n for n, *_ in cls) == 113

def test_repeat_and_conservation():
    loci, _ = summarise({"A": 3}, {frozenset("AB"): 10}, META, repeat=frozenset("B"))
    by = {r[0]: r for r in loci}
    assert by["B"][6] == "repeat" and by["A"][10] == "."
    t = em({frozenset("AB"): 10, frozenset("BC"): 6}, {"A": 3})
    assert abs(sum(t.values()) - 19) < 1e-6 and t["C"] < 3  # B is explained by A's placed reads, C keeps a remainder

def test_load_ties_reads_the_named_score(tmp_path):
    row = lambda n, a, a2, pa, m, m2, pm: "\t".join(map(str, [n, a, a2, 2, "c1", pa, m, m2, pm])) + "\n"
    t = tmp_path / "t.tsv"
    # as_only: AS tie, ms apart (gene vs processed pseudogene); ms_only: the reverse; swap: primary below its best record
    t.write_text("#as_table\n" + row("as_only", 100, 100, 100, 300, 100, 300) + row("ms_only", 100, 50, 100, 300, 300, 300) + row("swap", 100, 50, 50, 300, 100, 100))
    assert load_ties(t, floor_fn()) == {"ms_only": 300, "swap": 300}
    assert load_ties(t, floor_fn(ratio=0.98), "AS") == {"as_only": 100, "swap": 100}
    t.write_text("#as_table\n" + row("near", 100, 100, 100, 2000, 1985, 2000))  # 15 below: inside 0.98 x best, outside the window of 10
    assert load_ties(t, floor_fn()) == {} and load_ties(t, floor_fn(window=20)) == {"near": 2000} == load_ties(t, floor_fn(ratio=0.98))
    t.write_text("r\t100\t100\t2\tc1\t100\n")
    with pytest.raises(SystemExit):
        load_ties(t, floor_fn())
