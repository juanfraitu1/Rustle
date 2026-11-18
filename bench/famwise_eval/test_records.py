"""Tests of bench/famwise_eval (records.py, origin_score.py, retoss.py): hand-made truth rows and tiny unsorted BAMs."""
import os, subprocess, sys
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R

HERE = os.path.dirname(os.path.abspath(__file__))
HDR = "read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\n"
A, B, C = "c1:0-250:+", "c1:1000-1250:+", "c2:0-250:+"   # truth spans of three copies

def loci(*spans, AS=200):
    return ";".join(f"{i}:{s}:AS{AS}" for i, s in enumerate(spans))

def truth(tmp_path, rows, name="reads.tsv"):
    """rows = (read, status, origin_loci, loci); -> path"""
    (tmp_path / name).write_text(HDR + "".join(f"{r}\tS2\t0\t{st}\tprimary=0\t{o}\t{l}\n" for r, st, o, l in rows))
    return str(tmp_path / name)

def bam(tmp_path, recs, name="sel.bam"):
    """recs = (read, contig, start, flag, ms, AS[, cigar]); spliced 100=50N100= by default; unsorted, no index -> path"""
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 5000, "SN": "c1"}, {"LN": 5000, "SN": "c2"}]}
    with pysam.AlignmentFile(str(tmp_path / name), "wb", header=hdr) as w:
        for q, contig, start, flag, ms, AS, *cig in recs:
            a = pysam.AlignedSegment()
            a.query_name, a.flag, a.reference_id, a.reference_start, a.mapping_quality = q, flag, int(contig[1]) - 1, start, 0 if flag & 256 else 60
            a.cigarstring = cig[0] if cig else "100=50N100="
            a.query_sequence = "A" * a.infer_query_length()
            a.set_tag("AS", AS); a.set_tag("ms", ms)
            w.write(a)
    return str(tmp_path / name)

def test_kinds(tmp_path):
    twelve = loci(*[f"c1:{i * 1000}-{i * 1000 + 250}:+" for i in range(12)])
    split = tmp_path / "split.tsv"
    split.write_text("read\tfamily\tclass\nsev\tF\tSEVERAL COPIES\nonep\tF\tONE PLACE\nuns\tF\tUNSPLIT\n")
    t = R.load_truth(truth(tmp_path, [
        ("one", "decidable", "0", loci(A, B)), ("sev", "tie", "0,1", loci(A, B)), ("onep", "tie", "0,1", loci(A, B)),
        ("nosplit", "tie", "0,1", loci(A, B)), ("uns", "tie", "0,1", loci(A, B)), ("out", "outside", "", loci(A, B)),
        ("nohit", "no near-perfect hit", "", loci(A, B)), ("far", "decidable", "12", twelve)]), str(split))
    assert {q: v.kind for q, v in t.items()} == {"one": "one", "sev": "several", "onep": "one_place", "nosplit": "tie_unsplit",
                                                 "uns": "tie_unsplit", "out": "outside", "nohit": "unjudged", "far": "unjudged"}
    assert t["sev"].origin == frozenset({0, 1}) and t["one"].loci[1] == ("c1", 1000, 1250, "+", 200)
    # without a split table every tie is tie_unsplit
    assert R.load_truth(truth(tmp_path, [("sev", "tie", "0,1", loci(A, B))]))["sev"].kind == "tie_unsplit"

def test_locus_index(tmp_path):
    t = R.load_truth(truth(tmp_path, [("r", "decidable", "0", "0:c1:0-250:+:AS200;1:c1:1000-1250:+:AS199;2:c2:0-250:+:AS100")]))["r"]
    rec = lambda chrom, start, strand: R.Rec(chrom, start, start + 250, strand, (), 200, 200, True)
    assert R.locus_index(rec("c1", 1100, "+"), t) == 1
    assert R.locus_index(rec("c1", 1100, "-"), t) is None      # other strand
    assert R.locus_index(rec("c1", 3000, "+"), t) is None      # no truth locus there
    assert R.locus_index(rec("c2", 0, "+"), t) is None         # AS 100 < 0.8 x 200: listed, never judged

def test_statement_both_scores(tmp_path):
    b = bam(tmp_path, [("r", "c1", 0, 0, 200, 190), ("r", "c1", 1000, 256, 195, 196), ("r", "c2", 0, 256, 185, 150),
                       ("x", "c1", 0, 4, 0, 0), ("r", "c2", 3000, 2048, 300, 300)])
    recs = R.load_records(b)
    assert set(recs) == {"r"} and len(recs["r"]) == 3          # unmapped and supplementary records are skipped
    r = recs["r"]
    assert r[0].chain == ((100, 151),) and r[0].primary and r[0].strand == "+" and (r[0].start, r[0].end) == (0, 250)
    at = lambda s: sorted((x.chrom, x.start) for x in s)
    assert at(R.statement(r, "ms", 10)) == [("c1", 0), ("c1", 1000)]
    assert at(R.statement(r, "AS", 10)) == [("c1", 0), ("c1", 1000)]     # 196 and the primary (190)
    assert at(R.statement(r, "AS", 5)) == [("c1", 0), ("c1", 1000)]      # the primary is always in, here below the band
    assert at(R.statement(r, "ms", 20)) == [("c1", 0), ("c1", 1000), ("c2", 0)]
    assert at(R.statement(r, "ms", None)) == [("c1", 0)]

# ---- origin_score.py ----
SCORE = os.path.join(HERE, "origin_score.py")
D, E, F, G = "c2:1000-1250:+", "c1:2000-2250:+", "c1:3000-3250:+", "c1:4000-4250:+"
START = {"A": ("c1", 0), "B": ("c1", 1000), "C": ("c2", 0), "D": ("c2", 1000), "E": ("c1", 2000), "F": ("c1", 3000), "G": ("c1", 4000)}

def gtf(tmp_path, name, where):
    """one 2-exon transcript T<letter> per letter of `where`, with the chain the default test records carry"""
    out = ""
    for x in where:
        c, s = START[x]
        out += "".join(f'{c}\tx\texon\t{a}\t{b}\t.\t+\t.\tgene_id "g{x}"; transcript_id "T{x}";\n' for a, b in ((s + 1, s + 100), (s + 151, s + 250)))
    (tmp_path / name).write_text(out)
    return str(tmp_path / name)

def table(path):
    ls = [l.rstrip("\n").split("\t") for l in open(path) if not l.startswith("#")]
    return [dict(zip(ls[0], l)) for l in ls[1:]]

def score(tmp_path, rows, recs, arms=(), extra=()):
    out = tmp_path / "out"
    a = [x for n, w in arms for x in ("--arm", f"{n}={gtf(tmp_path, n + '.gtf', w)}")]
    p = subprocess.run([sys.executable, SCORE, "--truth", truth(tmp_path, rows), "--bam", bam(tmp_path, recs), "--out", str(out), *a, *extra], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return out

def cell(out, **want):
    r = [x for x in table(out / "molecule.tsv") if all(x[k] == v for k, v in want.items())]
    assert len(r) == 1, r
    return {k: int(r[0][k]) for k in ("n", "equals", "sound", "no_verdict")}

def test_molecule_table(tmp_path):
    out = score(tmp_path, [("m1", "decidable", "1", loci(A, B))], [("m1", "c1", 0, 0, 200, 200), ("m1", "c1", 1000, 256, 200, 200)])
    assert open(out / "molecule.tsv").readline().startswith("# limits:")
    assert cell(out, score="ms", window="primary", near_tied="yes", kind="one") == dict(n=1, equals=0, sound=0, no_verdict=0)
    assert cell(out, score="ms", window="0", near_tied="yes", kind="one") == dict(n=1, equals=0, sound=1, no_verdict=0)
    assert cell(out, score="AS", window="0", near_tied="yes", kind="one") == dict(n=1, equals=0, sound=1, no_verdict=0)

def test_claims(tmp_path):
    (tmp_path / "split.tsv").write_text("read\tfamily\tclass\ns1\tF\tSEVERAL COPIES\n")
    tx = "band\tnode\tlocus\tchrom\tstart\tend\tstrand\ttranscripts\tft\tlo\thi\tstatus\n" + "".join(
        f"wide\t{t}\tg\tc1\t1\t250\t+\t{t}\tFT1\t0\t1\t{s}\n" for t, s in (("TA", "certified"), ("TB", "tied_only"), ("TC", "tied_only"), ("TF", "no_band_record")))
    (tmp_path / "fw.tx.tsv").write_text(tx)
    out = score(tmp_path, [("p1", "decidable", "0", loci(A)), ("s1", "tie", "0,1", loci(B, E)), ("e1", "decidable", "1", loci(C, D))],
                [("p1", "c1", 0, 0, 200, 200), ("s1", "c1", 1000, 0, 200, 200), ("s1", "c1", 2000, 256, 200, 200),
                 ("e1", "c2", 0, 0, 200, 200), ("e1", "c2", 1000, 256, 200, 200)],
                arms=[("a1", "ABCF"), ("a2", "A")], extra=["--tie-split", str(tmp_path / "split.tsv"), "--tx", str(tmp_path / "fw.tx.tsv"), "--tx-arm", "a1"])
    got = {(r["arm"], r["status"]): tuple(int(r[k]) for k in ("claims", "TRUE", "SHARED", "PHANTOM", "UNJUDGED")) for r in table(out / "claims.tsv")}
    assert got == {("a1", "certified"): (1, 1, 0, 0, 0), ("a1", "tied_only"): (2, 0, 1, 1, 0), ("a1", "no_band_record"): (1, 0, 0, 0, 1),
                   ("a2", "."): (1, 1, 0, 0, 0)}
    bars = {r["quantity"]: r for r in table(out / "bars.tsv")}
    assert bars["B2i tied_only PHANTOM+SHARED"]["k"] == "2" and bars["B2i verdict"]["verdict"] == "undecidable"   # 2 judged tied_only claims < 30
    assert bars["B2ii certified+best_only of TRUE"]["k"] == "1" and bars["B2ii certified+best_only of TRUE"]["n"] == "1"
    # range: e1's exact band holds C and D, its origin is D: node TC has lo 0 <= true 0 <= hi 1; p1: TA 1 <= 1 <= 1
    rng = {r["band"]: r for r in table(out / "range.tsv")}
    assert (rng["exact"]["nodes"], rng["exact"]["held"], rng["exact"]["a_holds"], rng["exact"]["a_n"]) == ("2", "2", "2", "2")

def test_no_verdict_is_not_wrong(tmp_path):
    out = score(tmp_path, [("m1", "decidable", "0", loci(A))], [("m1", "c1", 0, 0, 200, 200), ("m1", "c1", 4000, 256, 200, 200)], arms=[("a1", "G")])
    assert cell(out, score="ms", window="primary", near_tied="yes", kind="one") == dict(n=1, equals=1, sound=1, no_verdict=0)
    assert cell(out, score="ms", window="0", near_tied="yes", kind="one") == dict(n=1, equals=0, sound=0, no_verdict=1)
    assert [(r["UNJUDGED"], r["PHANTOM"]) for r in table(out / "claims.tsv")] == [("1", "0")]   # its only supporter has no verdict there

def test_one_score_only(tmp_path):
    out = score(tmp_path, [("m1", "decidable", "0", loci(A, B))], [("m1", "c1", 0, 0, 200, 200), ("m1", "c1", 1000, 256, 185, 195)])
    assert table(out / "one_score_only.tsv") == [dict(read="m1", near_tied_on="AS", lead_ms="15", lead_AS="5")]
    assert cell(out, score="AS", window="primary", near_tied="yes", kind="one")["n"] == 1
    assert cell(out, score="ms", window="primary", near_tied="no", kind="one")["n"] == 1

def test_two_records_at_one_place_are_no_near_tie(tmp_path):
    # a second alignment at the same place (other splice) does not make the molecule near-tied
    out = score(tmp_path, [("m1", "decidable", "0", loci(A, B))],
                [("m1", "c1", 0, 0, 200, 200), ("m1", "c1", 0, 256, 198, 198, "100=60N90="), ("m1", "c1", 1000, 256, 150, 150)])
    assert cell(out, score="ms", window="primary", near_tied="no", kind="one")["n"] == 1

def test_unmatched_reads_counted(tmp_path):
    out = score(tmp_path, [("m1", "decidable", "0", loci(A)), ("gone", "decidable", "0", loci(A))],
                [("m1", "c1", 0, 0, 200, 200), ("new", "c1", 0, 0, 200, 200)])
    head = [l for l in open(out / "molecule.tsv") if l.startswith("#")]
    assert any("truth 2" in l and "BAM 2" in l and "both 1" in l for l in head), head

def test_family_table_and_b1(tmp_path):
    # 20 `several` molecules of family F1 tied between A and B: the exact band equals the truth, the primary alone never
    rows = [(f"m{i}", "tie", "0,1", loci(A, B)) for i in range(20)]
    recs = [x for i in range(20) for x in ((f"m{i}", "c1", 0, 0, 200, 200), (f"m{i}", "c1", 1000, 256, 200, 200))]
    (tmp_path / "split.tsv").write_text("read\tfamily\tclass\n" + "".join(f"m{i}\tF\tSEVERAL COPIES\n" for i in range(20)))
    (tmp_path / "fam.tsv").write_text("".join(f"m{i}\tF1\n" for i in range(20)))
    out = score(tmp_path, rows, recs, extra=["--tie-split", str(tmp_path / "split.tsv"), "--families", str(tmp_path / "fam.tsv")])
    f = [r for r in table(out / "molecule_family.tsv") if r["family"] == "F1" and r["window"] == "0" and r["kind"] == "several"]
    assert len(f) == 1 and (f[0]["n"], f[0]["equals"]) == ("20", "20")
    bars = {r["quantity"]: r for r in table(out / "bars.tsv")}
    assert (bars["B1 exact band equals"]["k"], bars["B1 primary equals"]["k"], bars["B1 exact band equals"]["n"]) == ("20", "0", "20")
    assert bars["B1 families exact > primary"]["k"] == "1" and bars["B1 families exact > primary"]["n"] == "1"
    assert bars["B1 verdict"]["verdict"] == "PASS"

# ---- retoss.py ----
def retoss(*args):
    p = subprocess.run([sys.executable, os.path.join(HERE, "retoss.py"), *map(str, args)], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return [l.split("\t") for l in p.stdout.splitlines()]

def test_retoss(tmp_path):
    a = bam(tmp_path, [("r1", "c1", 0, 0, 200, 200), ("r1", "c1", 1000, 256, 200, 200), ("r2", "c1", 0, 0, 200, 200),
                       ("r3", "c1", 0, 0, 200, 200), ("r3", "c2", 0, 256, 190, 190)], "a.bam")
    b = bam(tmp_path, [("r1", "c1", 0, 256, 200, 200), ("r1", "c1", 1000, 0, 200, 200), ("r2", "c1", 0, 0, 200, 200),
                       ("r3", "c1", 0, 0, 200, 200)], "b.bam")
    assert retoss("records", a, b) == [["records_in_both", "4"], ["only_A", "1"], ["only_B", "0"], ["reads_primary_differs", "1"], ["reads", "3"]]
    assert retoss("records", a, a)[1:4] == [["only_A", "0"], ["only_B", "0"], ["reads_primary_differs", "0"]]
    g = lambda n, w: gtf(tmp_path, n, w)
    assert retoss("claims", f"x={g('1.gtf', 'AB')},{g('2.gtf', 'A')}", f"y={g('3.gtf', 'C')},{g('4.gtf', 'C')}") == [
        ["arm", "gtfs", "chains_in_any", "chains_in_every", "not_in_every", "share"], ["x", "2", "2", "1", "1", "0.5000"], ["y", "2", "1", "1", "0", "0.0000"]]
