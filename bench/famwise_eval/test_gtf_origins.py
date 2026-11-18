"""Tests of gtf_origins.py (PREREG r1799): simulated transcripts and their reads in the form sim_chains.py takes."""
import csv, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gtf_origins as G

GTF = '''c1\tx\ttranscript\t1\t300\t.\t+\t.\tgene_id "g1"; transcript_id "FBtr1"; gene_name "A";
c1\tx\texon\t201\t300\t.\t+\t.\tgene_id "g1"; transcript_id "FBtr1";
c1\tx\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "FBtr1";
c1\tx\texon\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "FBtr2";
c2\tx\texon\t11\t50\t.\t-\t.\tgene_id "g2"; transcript_id "FBti9";
c2\tx\texon\t11\t50\t.\t-\t.\tgene_id "g3"; transcript_id "FBtr7";
'''
EXACT = "read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\tsource\thome\n" + "".join(
    f"{q}\tS0\tNone\tdecidable\tprimary=0\t0\t\t{t}\t1\n" for q, t in (("r1", "FBtr1"), ("r2", "FBtr1"), ("r3", "FBti9")))

def run(tmp_path, *extra):
    (tmp_path / "t.gtf").write_text(GTF); (tmp_path / "e.tsv").write_text(EXACT)
    G.main([str(tmp_path / "t.gtf"), str(tmp_path / "e.tsv"), str(tmp_path / "o"), *extra])
    fam = json.load(open(tmp_path / "o.fam.json"))["sets"]["sim"]
    return fam, list(csv.DictReader(open(tmp_path / "o.origin.tsv"), delimiter="\t"))

def test_only_transcripts_with_reads_and_every_read_is_full(tmp_path):
    fam, rows = run(tmp_path)
    assert [f["family"] for f in fam] == ["g1", "g2"]                                   # FBtr2 and g3 have no read
    c = fam[0]["copies"][0]
    assert c["transcripts"] == {"FBtr1": [[0, 100], [200, 300]]} and (c["chrom"], c["strand"], c["start"], c["end"], c["gene_id"]) == ("c1", "+", 0, 300, "g1")
    assert [(r["read"], r["member"], r["family"], r["isoform"], r["full"]) for r in rows] == [("r1", "g1", "g1", "FBtr1", "1"), ("r2", "g1", "g1", "FBtr1", "1"), ("r3", "g2", "g2", "FBti9", "1")]

def test_prefix_keeps_one_source_type(tmp_path):
    fam, rows = run(tmp_path, "--tid-prefix", "FBti")
    assert [f["family"] for f in fam] == ["g2"] and [r["read"] for r in rows] == ["r3"]
