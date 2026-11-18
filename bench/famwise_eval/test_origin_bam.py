"""Tests of origin_bam.py and region_chains.py (bench/famwise_eval)."""
import os, subprocess, sys
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_records import A, B, C, D, F, HERE, bam, loci, truth

def run(tmp_path, script, *args):
    p = subprocess.run([sys.executable, os.path.join(HERE, script), *map(str, args)], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return p.stdout

def recs_of(path):
    assert os.path.exists(str(path) + ".bai")
    with pysam.AlignmentFile(str(path)) as fh:
        return sorted((r.query_name, r.reference_name, r.reference_start, r.is_secondary) for r in fh)

def test_three_inputs(tmp_path):
    rows = [("mv", "decidable", "1", loci(A, B)), ("sv", "tie", "0,1", loci(A, B, C)), ("nh", "decidable", "1", loci(D, F))]
    sel = bam(tmp_path, [("mv", "c1", 0, 0, 200, 200), ("mv", "c1", 1000, 256, 200, 200),
                         ("sv", "c1", 0, 0, 200, 200), ("sv", "c1", 1000, 256, 200, 200), ("sv", "c2", 0, 256, 200, 150),
                         ("nt", "c2", 0, 0, 200, 200), ("nt", "c1", 4000, 256, 200, 200),
                         ("nh", "c2", 1000, 0, 200, 200), ("nh", "c1", 2000, 256, 190, 190)])
    out = tmp_path / "o"
    run(tmp_path, "origin_bam.py", "--truth", truth(tmp_path, rows), "--bam", sel, "--out", out)
    assert recs_of(out / "D.bam") == [("mv", "c1", 0, False), ("mv", "c1", 1000, True), ("nh", "c2", 1000, False), ("nt", "c2", 0, False),
                                     ("sv", "c1", 0, False), ("sv", "c1", 1000, True), ("sv", "c2", 0, True)]
    # moved: the secondary at the origin is promoted and the old primary is gone; several: two records, one primary
    assert recs_of(out / "O.bam") == [("mv", "c1", 1000, False), ("nh", "c2", 1000, False), ("nt", "c2", 0, False),
                                     ("sv", "c1", 0, False), ("sv", "c1", 1000, True)]
    assert [l.split() for l in open(out / "regions.bed")] == [["c1", "0", "250"], ["c1", "1000", "1250"]]
    t = dict(l.rstrip("\n").split("\t") for l in open(out / "inputs.tsv") if not l.startswith("#"))
    assert (t["kept_as_is"], t["restricted"], t["moved"], t["no_home"], t["primary_only"]) == ("0", "1", "1", "1", "1")
    assert t["D_records_tied_on_ms_below_the_AS_band"] == "1"      # sv's record on c2: ms 200 = best, AS 150 < 0.98 x 200

def test_region_chains(tmp_path):
    from test_records import gtf
    (tmp_path / "r.bed").write_text("c1\t0\t250\nc1\t1000\t1250\n")
    out = run(tmp_path, "region_chains.py", gtf(tmp_path, "ref.gtf", "ABC"), tmp_path / "r.bed",
              f"P={gtf(tmp_path, 'p.gtf', 'A')}", f"D={gtf(tmp_path, 'd.gtf', 'ABF')}", f"O={gtf(tmp_path, 'o.gtf', 'AB')}")
    rows = [l.split("\t") for l in out.splitlines()]
    assert rows[0] == ["scope", "input", "chains", "refseq_exact"]
    assert rows[1:7] == [["regions", "P", "1", "1"], ["regions", "D", "2", "2"], ["regions", "O", "2", "2"],
                         ["genome", "P", "1", "1"], ["genome", "D", "3", "2"], ["genome", "O", "2", "2"]]
    assert rows[7:] == [["R", "regions", "2.00"], ["R", "genome", "2.00"], ["D_not_O", "regions", "0", "0"], ["D_not_O", "genome", "1", "0"]]

def test_regions_are_the_tied_or_moved_reads(tmp_path):
    """PREREG r1794: regions come from reads that are tied or whose primary is off its origin. A one-origin read at home
    with a secondary elsewhere adds none (O = P for it); a tie whose records all lie on its origin copies adds both."""
    X, Y, Z, W = "c2:2000-2250:+", "c2:3000-3250:+", "c2:4000-4250:+", "c1:3000-3250:+"
    rows = [("o1", "decidable", "0", loci(X, W)), ("t2", "tie", "0,1", loci(Y, Z))]
    sel = bam(tmp_path, [("o1", "c2", 2000, 0, 200, 200), ("o1", "c1", 3000, 256, 190, 190),
                         ("t2", "c2", 3000, 0, 200, 200), ("t2", "c2", 4000, 256, 200, 200)])
    out = tmp_path / "o"
    run(tmp_path, "origin_bam.py", "--truth", truth(tmp_path, rows), "--bam", sel, "--out", out)
    assert [l.split() for l in open(out / "regions.bed")] == [["c2", "3000", "3250"], ["c2", "4000", "4250"]]
