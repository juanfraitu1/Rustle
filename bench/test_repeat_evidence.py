import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import repeat_evidence as re_  # noqa: E402

RMSK = """   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  311  32.7  3.2  3.9  chrA         11      20   (14400) +  AluY           SINE/Alu              1  155 (1128) 1
  650   5.1  0.0  0.0  chrB          1      78    (13839) +  (CA)n          Simple_repeat          1   78    (0) 2
"""


def test_parse_rmsk_plain_and_gz(tmp_path):
    p = tmp_path / "x.out"
    p.write_text(RMSK)
    g = tmp_path / "x.out.gz"
    with gzip.open(g, "wt") as fh:
        fh.write(RMSK)
    want = [("chrA", 10, 20, "SINE/Alu"), ("chrB", 0, 78, "Simple_repeat")]
    assert re_.parse_rmsk(p) == want and re_.parse_rmsk(g) == want and re_.parse_rmsk(p, {"chrA"}) == want[:1]


def test_lowercase_runs(tmp_path):
    fa = tmp_path / "t.fa"
    fa.write_text(">c1\nACgtaCG\nttA\n>c2\nACGT\n")
    import pysam
    pysam.faidx(str(fa))
    assert re_.lowercase_runs(fa) == [("c1", 2, 5, "."), ("c1", 7, 9, ".")]


def test_masker_interval_parse(tmp_path):
    p = tmp_path / "wm.txt"
    p.write_text(">c1 some description\n0 - 9\n20 - 20\n>c2\n5 - 6\n")
    assert re_.parse_masker_intervals(p) == [("c1", 0, 10, "."), ("c1", 20, 21, "."), ("c2", 5, 7, ".")]


def test_merge_and_masked_bases():
    m = re_.merge([("c", 0, 10, "."), ("c", 5, 20, "."), ("c", 30, 40, ".")])
    assert m == {"c": [(0, 20), (30, 40)]}
    assert re_.masked_bases(m, "c", 15, 35) == 10 and re_.masked_bases(m, "x", 0, 5) == 0


def test_bed_roundtrip(tmp_path):
    ivs = [("c", 0, 10, "SINE/Alu")]
    re_.write_bed(ivs, tmp_path / "r.bed", "R1:rmsk")
    assert re_.read_bed(tmp_path / "r.bed") == ivs
