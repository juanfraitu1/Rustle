"""Tests of tools/famwise.py and of the --famwise switch of tools/eq_classes.py (driver stage `famwise`)."""
import collections, csv, os, random, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famwise as fw
from gtf_names import chain

GTF = "\n".join("\t".join(r) for r in [
    ("c1", "x", "exon", "1", "100", ".", "+", ".", 'gene_id "A"; transcript_id "A.2";'),
    ("c1", "x", "exon", "151", "250", ".", "+", ".", 'gene_id "A"; transcript_id "A.2";'),
    ("c1", "x", "exon", "1", "100", ".", "+", ".", 'gene_id "A"; transcript_id "A.1";'),
    ("c1", "x", "exon", "151", "300", ".", "+", ".", 'gene_id "A"; transcript_id "A.1";'),
    ("c1", "x", "exon", "500", "700", ".", ".", ".", 'gene_id "S"; transcript_id "S.1";'),
]) + "\n"

def test_record_chain():
    assert fw.record_chain(0, [(7, 100), (3, 50), (8, 1), (7, 99)]) == ((100, 151),) == chain([(1, 100), (151, 250)])
    assert fw.record_chain(0, [(7, 40), (2, 3), (7, 57)]) == ()          # a deletion is not an intron
    assert fw.record_chain(9, [(4, 5), (7, 10), (3, 20), (7, 10)]) == ((19, 40),)   # soft clip consumes no reference

def test_load_nodes(tmp_path):
    p = tmp_path / "a.gtf"; p.write_text(GTF)
    nodes, tx = fw.load_nodes(str(p))
    assert nodes == {"A": {((100, 151),): "A.1"}} and tx == {"A.1": ["A.1", "A.2"]}   # same chain: one node; S is single-exon
    assert fw.member_of("A", ((100, 151),), nodes) == "A.1"
    assert fw.member_of("A", ((100, 160),), nodes) == fw.OTHER
    assert fw.member_of("A", (), nodes) == fw.UNSPLICED and fw.member_of(fw.OUTSIDE, (), nodes) == "."

A1, A2, B1 = ("A", "a1"), ("A", "a2"), ("B", "b1")

def L(records, base=(), **kw):
    s = fw.bands(records); b = collections.Counter(dict(base))
    return fw.layer(s["exact"], b, **kw), fw.layer(s["wide"], b, **kw), s, b

def test_two_tied_emitted_chains():
    e, w, s, b = L({"r1": [("A", "a1", 0, 1), ("B", "b1", 0, 0)], "r2": [("A", "a1", 0, 0), ("B", "b1", 0, 1)]})
    assert e["ft"][A1] == e["ft"][B1] and e["inside"][e["ft"][A1]] == 2 and not e["open"]
    assert (e["lo"][A1], e["hi"][A1], e["hi"][B1]) == (0, 2, 2)
    assert fw.statuses(e, w, [A1, B1]) == {A1: "tied_only", B1: "tied_only"}
    assert fw.min_copies(set(), [frozenset("AB")]) == (1, True)

def test_one_emitted_one_not():
    e, w, s, b = L({"r1": [("A", "a1", 0, 1), ("B", fw.OTHER, 0, 0)]})
    k = e["ft"][A1]
    assert e["open"][k] == 1 and e["inside"][k] == 0 and e["hi"][A1] == 1 and e["lo"][A1] == 0
    assert sum(e["placed"].values()) + sum(e["classes"].values()) == 1            # P1: nobody lost or doubled

def test_two_records_at_one_locus():
    e, w, s, b = L({"r1": [("A", "a1", 0, 1), ("A", "a2", 0, 0), ("B", "b1", 0, 0)], "r2": [("A", "a1", 0, 1), ("A", "a2", 0, 0)]})
    assert e["ft"][A1] == e["ft"][A2] == e["ft"][B1] and e["ft"][A1] in e["collision"]
    pl, mx = fw.locus_counts(s["exact"], b)
    assert (pl["A"], mx["A"], pl["B"], mx["B"]) == (1, 2, 0, 1) and e["lo"][A1] == 0   # r2: placed on the locus, on no node

def test_statuses_and_primary_below_best():
    recs = {"r1": [("A", "a1", 0, 1), ("B", "b1", 5, 0)],          # best at A, B within 10
            "r2": [("C", "c1", 3, 1), ("D", "d1", 0, 0)]}          # primary 3 below the best record: tied between both
    e, w, s, b = L(recs, base={("E", "e1"): 1})
    st = fw.statuses(e, w, [A1, B1, ("C", "c1"), ("D", "d1"), ("E", "e1"), ("F", "f1")])
    assert st == {A1: "best_only", B1: "near_only", ("C", "c1"): "tied_only", ("D", "d1"): "tied_only",
                  ("E", "e1"): "certified", ("F", "f1"): "no_band_record"}
    assert s["exact"]["r2"] == frozenset({("C", "c1"), ("D", "d1")})

def test_min_link_repeat_and_unseen():
    recs = {"r1": [("A", "a1", 0, 1), ("B", "b1", 0, 0)]}
    assert L(recs)[0]["inside"] and not L(recs, min_link=2)[0]["inside"] and L(recs, min_link=2)[0]["open"]
    assert not L(recs, repeat=frozenset({"B"}))[0]["inside"]
    e, w, s, b = L({"r1": [("A", "a1", 0, 1), (fw.UNSEEN, ".", 0, 0)]})
    assert e["lo"][A1] == 0 and e["open"][e["ft"][A1]] == 1 and fw.locus_counts(s["exact"], b)[0]["A"] == 1

def test_min_copies():
    assert fw.min_copies(set(), [frozenset("AB"), frozenset("BC")]) == (1, True)
    assert fw.min_copies({"A"}, [frozenset("BC")]) == (2, True)
    assert fw.min_copies(set(), [frozenset("AB"), frozenset("CD"), frozenset("BC")], cap=0) == (2, False)

def test_invariance_and_bounds_random():
    rng = random.Random(7); nodes = [(l, l.lower() + str(i)) for l in "ABCD" for i in (1, 2)]
    recs, origin = {}, {}
    for i in range(300):
        k = rng.sample(nodes, rng.choice((1, 1, 2, 3))); origin[i] = k[0]
        recs[i] = [(l, m, 0, int(j == 0)) for j, (l, m) in enumerate(k)]
    e, w, s, b = L(recs)
    flipped = {i: [(l, m, g, int(j == len(r) - 1)) for j, (l, m, g, p) in enumerate(r)] for i, r in recs.items()}
    assert fw.bands(flipped) == s                                                  # P2: the coin does not matter
    x = collections.Counter(origin.values())
    assert all(e["lo"][n] <= x[n] <= e["hi"][n] for n in nodes)                    # P3 per node
    for key in set(e["ft"].values()):                                              # P5: exact total without open molecules
        comp = [n for n in nodes if e["ft"].get(n) == key]
        if not e["open"][key]:
            assert sum(x[n] for n in comp) == sum(e["placed"][n] for n in comp) + e["inside"][key]

EQ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "eq_classes.py")
FWPY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "famwise.py")

def tiny(tmp_path, bcontig="c1", more=(), more_rows=()):
    """c1: locus A (exons 1-100, 151-250, +); locus B (1001-1100, 1151-1250, +) on bcontig; one 2-exon transcript each
    (A.1, B.1). Every record carries AS equal to its ms. t1, t2: primary at A, secondary at B, ms 200 / 200.
    u1: one record at A. w1: primary at A ms 200, secondary at B ms 185. s1: unspliced primary at A.
    -> (gtf, bam, score table) paths."""
    import pysam
    ex = lambda c, g, s: "".join(f'{c}\tx\texon\t{a}\t{b}\t.\t+\t.\tgene_id "{g}"; transcript_id "{g}.1";\n' for a, b in ((s + 1, s + 100), (s + 151, s + 250)))
    (tmp_path / "t.gtf").write_text(ex("c1", "A", 0) + ex(bcontig, "B", 1000))
    recs = [("t1", "c1", 0, 0, 200), ("t1", bcontig, 1000, 256, 200), ("t2", "c1", 0, 0, 200), ("t2", bcontig, 1000, 256, 200),
            ("u1", "c1", 0, 0, 200), ("w1", "c1", 0, 0, 200), ("w1", bcontig, 1000, 256, 185), ("s1", "c1", 0, 0, 100), *more]
    hdr = {"HD": {"VN": "1.6"}, "SQ": [{"LN": 5000, "SN": "c1"}, {"LN": 5000, "SN": "c2"}]}
    with pysam.AlignmentFile(str(tmp_path / "u.bam"), "wb", header=hdr) as w:
        for name, contig, start, flag, ms in recs:
            a = pysam.AlignedSegment()
            a.query_name, a.flag, a.reference_id, a.reference_start, a.mapping_quality = name, flag, int(contig[1]) - 1, start, 60
            a.cigarstring = "100=" if name == "s1" else "100=50N100="
            a.query_sequence = "A" * (100 if name == "s1" else 200)
            a.set_tag("AS", ms); a.set_tag("ms", ms)
            w.write(a)
    pysam.sort("-o", str(tmp_path / "t.bam"), str(tmp_path / "u.bam")); pysam.index(str(tmp_path / "t.bam"))
    row = lambda n, best, second, prim: "\t".join(map(str, [n, best, second, 2 if second > 0 else 1, "c1", prim, best, second, prim, 60, 1, "NA", "NA", "NA", "NA", 0, "NA"])) + "\n"
    (tmp_path / "t.tab").write_text("#as_table\ttest\n" + row("t1", 200, 200, 200) + row("t2", 200, 200, 200) + row("u1", 200, -1, 200) + row("w1", 200, 185, 200) + row("s1", 100, -1, 100) + "".join(row(*r) for r in more_rows))
    return str(tmp_path / "t.gtf"), str(tmp_path / "t.bam"), str(tmp_path / "t.tab")

def run_eq(tmp_path, out, *x, **paths):
    gtf, bam, tab = paths.get("inputs") or tiny(tmp_path)
    return subprocess.run([sys.executable, EQ, gtf, bam, tab, str(tmp_path / out), "--threads", "1", *x], capture_output=True, text=True)

def test_switch_keeps_eq_tables_and_writes_reads(tmp_path):
    inp = tiny(tmp_path)
    a, b = run_eq(tmp_path, "off", inputs=inp), run_eq(tmp_path, "on", "--famwise", str(tmp_path / "fw"), inputs=inp)
    assert a.returncode == b.returncode == 0, a.stderr + b.stderr
    assert a.stdout == b.stdout and len(a.stdout.splitlines()) == 2
    for t in ("loci", "classes"):
        assert (tmp_path / f"off.{t}.tsv").read_bytes() == (tmp_path / f"on.{t}.tsv").read_bytes()
    assert "A,B" in (tmp_path / "off.classes.tsv").read_text()                    # the fixture does tie t1, t2 over A and B
    rows = [l.split("\t") for l in (tmp_path / "fw.reads.tsv").read_text().splitlines()]
    assert rows[0] == "read locus member gap ms AS primary".split()
    assert ["w1", "B", "B.1", "15", "185", "185", "0"] in rows and sum(r[0] == "t1" for r in rows) == 2
    assert not any(r[1] == fw.UNSEEN for r in rows)
    assert (tmp_path / "fw.placed.tsv").read_text().splitlines() == ["locus\tmember\tn", "A\tA.1\t1", "A\tU\t1"]

def test_ratio_with_switch_is_refused(tmp_path):
    r = run_eq(tmp_path, "bad", "--famwise", str(tmp_path / "x"), "--ratio", "0.98")
    assert r.returncode == 2 and "--famwise needs the window band" in r.stderr

def test_unseen_rival(tmp_path):
    inp = tiny(tmp_path, bcontig="c2")                      # the tied record of t1 lies on a contig the panel does not scan
    r = run_eq(tmp_path, "p", "--contigs", "c1", "--famwise", str(tmp_path / "fw"), inputs=inp)
    assert r.returncode == 0, r.stderr
    rows = [l.split("\t") for l in (tmp_path / "fw.reads.tsv").read_text().splitlines()]
    assert ["t1", "unseen", ".", "0", "200", "NA", "0"] in rows and ["w1", "unseen", ".", "15", "185", "NA", "0"] in rows

def run_fw(gtf, prefix, *x):
    return subprocess.run([sys.executable, FWPY, str(gtf), str(prefix), *x], capture_output=True, text=True)

def table(path):
    return list(csv.DictReader(open(path), delimiter="\t"))

def test_end_to_end_tiny(tmp_path):
    inp = tiny(tmp_path)
    assert run_eq(tmp_path, "eq", "--famwise", str(tmp_path / "fw"), inputs=inp).returncode == 0
    r = run_fw(inp[0], tmp_path / "fw", "--eq-loci", str(tmp_path / "eq.loci.tsv"))
    assert r.returncode == 0, r.stderr                                             # P1 holds against the eq table
    top = next(f for f in table(tmp_path / "fw.families.tsv") if f["band"] == "exact" and f["ft"] == "FT1")
    assert (top["nodes"], top["copies"], top["inside"], top["open"], top["min_copies"], top["collision"], top["rep"]) == ("2", "A,B", "2", "0", "1", "no", "A.1")
    tx = {(t["band"], t["node"]): t for t in table(tmp_path / "fw.tx.tsv")}
    assert (tx["wide", "A.1"]["status"], tx["wide", "A.1"]["lo"], tx["wide", "A.1"]["hi"]) == ("certified", "2", "4")   # u1 and w1 fit nowhere else
    assert (tx["wide", "B.1"]["status"], tx["wide", "B.1"]["lo"], tx["wide", "B.1"]["hi"]) == ("tied_only", "0", "2")
    cl = table(tmp_path / "fw.classes.tsv")
    assert [(c["band"], c["reads"], c["members"], c["ft"], c["kind"]) for c in cl] == [(b, "2", "A|A.1,B|B.1", "FT1", "inside") for b in ("exact", "wide")]
    assert len(r.stdout.splitlines()) == 3

def test_p1_mismatch_exits_1(tmp_path):
    inp = tiny(tmp_path)
    run_eq(tmp_path, "eq", "--famwise", str(tmp_path / "fw"), inputs=inp)
    rows = (tmp_path / "eq.loci.tsv").read_text().splitlines()
    f = rows[1].split("\t"); f[7] = str(int(f[7]) + 1); rows[1] = "\t".join(f)       # one `placed` value off by one
    (tmp_path / "bad.loci.tsv").write_text("\n".join(rows) + "\n")
    r = run_fw(inp[0], tmp_path / "fw", "--eq-loci", str(tmp_path / "bad.loci.tsv"))
    assert r.returncode == 1 and "P1 fails at locus" in r.stderr

TX = lambda g, t, a, b, s="+": f'c1\tx\ttranscript\t{a}\t{b}\t.\t{s}\t.\tgene_id "{g}"; transcript_id "{t}";\n'

def hand(tmp_path, gtf, reads="", placed=""):
    (tmp_path / "g.gtf").write_text(gtf)
    (tmp_path / "fw.reads.tsv").write_text("read\tlocus\tmember\tgap\tms\tAS\tprimary\n" + reads)
    (tmp_path / "fw.placed.tsv").write_text("locus\tmember\tn\n" + placed)
    return run_fw(tmp_path / "g.gtf", tmp_path / "fw")

def test_gtf_attributes(tmp_path):
    gtf = TX("A", "A.2", 1, 250) + TX("A", "A.1", 1, 300) + TX("S", "S.1", 500, 700, ".") + GTF     # A.1 and A.2 share a chain
    r = hand(tmp_path, gtf, placed="A\tA.1\t3\n")
    assert r.returncode == 0, r.stderr
    out = (tmp_path / "fw.gtf").read_text().splitlines()
    add = ' fw_tx "FT1"; fw_status "certified"; fw_lo "3"; fw_hi "3"; fw_copies "1"; fw_rep "yes";'
    assert out[0] == TX("A", "A.2", 1, 250).rstrip("\n") + add and out[1] == TX("A", "A.1", 1, 300).rstrip("\n") + add
    assert out[2:] == gtf.splitlines()[2:]                                          # single-exon S.1 and every exon line untouched
    assert [t["transcripts"] for t in table(tmp_path / "fw.tx.tsv")] == ["A.1,A.2", "A.1,A.2"]

def test_empty(tmp_path):
    gtf = TX("S", "S.1", 500, 700, ".") + GTF.splitlines(True)[-1]
    r = hand(tmp_path, gtf)
    assert r.returncode == 0, r.stderr
    assert all(len((tmp_path / f"fw.{t}.tsv").read_text().splitlines()) == 1 for t in ("tx", "classes", "families"))
    assert (tmp_path / "fw.gtf").read_text() == gtf

def test_reads_of_another_gtf_are_refused(tmp_path):
    r = hand(tmp_path, GTF, placed="A\tZ.9\t1\n")
    assert r.returncode == 1 and "Z.9" in r.stderr

def test_unscanned_primary_certifies_nothing(tmp_path):
    # p1: best record a secondary at A (ms 200), primary at B (ms 170) on a contig the panel does not scan
    inp = tiny(tmp_path, bcontig="c2", more=[("p1", "c1", 0, 256, 200), ("p1", "c2", 1000, 0, 170)], more_rows=[("p1", 200, 170, 170)])
    assert run_eq(tmp_path, "eq", "--contigs", "c1", "--famwise", str(tmp_path / "fw"), inputs=inp).returncode == 0
    rows = [l.split("\t") for l in (tmp_path / "fw.reads.tsv").read_text().splitlines()]
    assert ["p1", "unseen", ".", "30", "170", "NA", "1"] in rows
    r = run_fw(inp[0], tmp_path / "fw", "--eq-loci", str(tmp_path / "eq.loci.tsv"))
    assert r.returncode == 0, r.stderr                                             # the panel's eq table still agrees per locus
    tx = {(t["band"], t["node"]): t for t in table(tmp_path / "fw.tx.tsv")}
    assert (tx["wide", "A.1"]["lo"], tx["wide", "A.1"]["hi"]) == ("2", "5")        # u1, w1 placed; t1, t2, p1 open, as on the whole sample
