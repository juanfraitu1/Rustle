import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("loop_home", os.path.join(HERE, "loop_home.py"))
LH = importlib.util.module_from_spec(spec)
spec.loader.exec_module(LH)

SPANS = {"F.0": ("chr5", 100, 200), "F.1": ("chr5", 900, 1000)}


def test_o2_homes_the_reads_the_corrector_decided_at_its_copy_and_the_oracle_every_read_at_its_true_copy():
    rows = [dict(read="a", corrector="F.1", minimap2="F.0", reason="sample"),     # moved by O2
            dict(read="b", corrector="F.0", minimap2="F.0", reason="sample"),     # decided, as minimap2 had it
            dict(read="c", corrector="F.0", minimap2="F.0", reason="tie")]        # undecided: keeps the base rule
    assert LH.o2_homes(rows, SPANS) == {"a": ("chr5", 900, 1000, "o2", "F.1"), "b": ("chr5", 100, 200, "o2", "F.0")}
    truth = {"a": ("unique", "F.1"), "b": ("poor_match", ""), "c": ("unique_elsewhere", "F.1")}
    assert LH.o2_homes(rows, SPANS, truth) == {"a": ("chr5", 900, 1000, "oracle", "F.1"),
                                               "c": ("chr5", 900, 1000, "oracle", "F.1")}


def test_pass2_homes_every_family_read_once_except_the_ties_o2_could_not_break():
    prim = {"a": "F.0", "b": "F.0", "c": "F.0", "d": "F.1"}                       # every read's minimap2 primary copy
    rows = [dict(read="a", corrector="F.1", minimap2="F.0", reason="sample"),
            dict(read="b", corrector="F.0", minimap2="F.0", reason="tie")]
    h = LH.pass2_homes(prim, rows, SPANS)
    assert h == {"a": ("chr5", 900, 1000, "o2", "F.1"),                           # O2 decided: its copy
                 "c": ("chr5", 100, 200, "minimap2", "F.0"),                      # not tied within 10 AS: minimap2's
                 "d": ("chr5", 900, 1000, "minimap2", "F.1")}                     # b, a tie O2 left: no home (echoes)
    truth = {"b": ("unique", "F.1"), "d": ("poor_match", "")}
    o = LH.pass2_homes(prim, rows, SPANS, truth)
    assert o["b"] == ("chr5", 900, 1000, "oracle", "F.1") and o["d"][3] == "minimap2"
    assert "a" not in o                                                           # a tie without a truth: no home
