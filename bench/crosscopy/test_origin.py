"""bench/crosscopy/origin.py: do the reads behind a cross-copy call come from the copy they are counted at?"""
import os, sys
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from origin import placed_introns, read_verdict, score, main
from test_calls import write, HEADER
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "famwise_eval"))
from records import Truth, Rec

T = lambda kind, origin, loci: Truth("S1", "x", kind, frozenset(origin), loci)
R = lambda chain, ms, primary, chrom="c1", strand="+": Rec(chrom, 0, 1000, strand, chain, ms, ms, primary)

def test_read_verdicts():
    ex = [("c1", "+", 101, 200), ("c1", "+", 301, 400)]                         # copy X's column exons, 1-based closed
    here, both = [("c1", 150, 320, "+", 900)], [("c1", 150, 320, "+", 900), ("c2", 10, 500, "+", 890)]
    assert read_verdict(T("one", {0}, here), ex) == "SUPPORTED" and read_verdict(T("one", {1}, both), ex) == "AGAINST"
    assert read_verdict(T("one", {0}, [("c1", 150, 320, "-", 900)]), ex) == "AGAINST"            # antisense
    assert read_verdict(T("one", {0}, [("c1", 200, 300, "+", 900)]), ex) == "AGAINST"            # same contig, between the exons
    assert read_verdict(T("several", {0, 1}, both), ex) == "AGAINST"
    assert read_verdict(T("outside", (), here), ex) == "OUTSIDE"
    assert read_verdict(T("unjudged", (), here), ex) is None and read_verdict(T("one_place", {0}, here), ex) is None
    assert read_verdict(None, ex) is None                                                        # a read the truth run never scored

def test_placed_needs_every_record_in_the_band():
    a, b, far = R(((200, 300),), 900, True), R(((200, 300), (400, 500)), 895, False), R((), 700, False)
    assert placed_introns([a, b, far]) == {("c1", "+", 201, 299)}                                # 1-based first / last intron base
    assert placed_introns([a, R(((200, 300),), 900, False, chrom="c2")]) == set()                # the tied record is on another contig

def score_fixture(tmp_path, u_key, u_other, s_key, s_other, s_outside_key=0, s_from_u=0):
    """one family, copies c1 (U) and c2 (S): exons 61-160, 361-420, 621-720; the call = the inclusion junction 161-360 at the
    donor, stated 0.5 at U and 0.0 at S. Reads: N on the junction and N skipping it (161-620) per copy, each with its truth."""
    d = tmp_path / f"f{u_key}_{s_key}_{s_outside_key}_{s_from_u}"; d.mkdir()
    m = {"c1": "c1:1-800", "c2": "c2:1-800"}
    psi = write(d, [("F", m["c1"], "donor", "C0", "C1", "c1", "+", 161, 360, 15, 30, 0.5, 0.5, 0.5, 0),
                    ("F", m["c2"], "donor", "C0", "C1", "c2", "+", 161, 360, 0, 30, 0.0, 0.0, 0.0, 1)])
    (d / "cols.tsv").write_text("family\tmembers\tcolumn\tclass\tn_loci\tn_exons\tcollision\tlen_ratio\texons\n" + "".join(
        f"F\t2\t{c}\tcore\t2\t2\t0\t1.0\t{m['c1']}#0:{s - 1}-{e},{m['c2']}#0:{s - 1}-{e}\n" for c, s, e in (("C0", 61, 160), ("C1", 361, 420), ("C2", 621, 720))))
    (d / "loci.gff3").write_text("".join(f"{c}\t.\tgene\t1\t800\t.\t+\t.\tID=g{c}\n" for c in m))
    reads, n = [], 0
    def add(chrom, key, count, status, origin_chrom=None):
        nonlocal n
        for _ in range(count):
            n += 1; reads.append((f"r{n}", chrom, key, status, origin_chrom or chrom))
    add("c1", True, u_key, "decidable"); add("c1", False, u_other, "decidable")
    add("c2", True, s_key, "decidable"); add("c2", False, s_other, "decidable")
    add("c2", True, s_outside_key, "outside"); add("c2", False, s_from_u, "decidable", origin_chrom="c1")
    cnt = {}
    with pysam.AlignmentFile(str(d / "sel.bam"), "wb", header={"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": 1000} for c in m]}) as bam:
        for name, chrom, key, _, _ in reads:
            a = pysam.AlignedSegment(); a.query_name, a.flag, a.reference_id, a.reference_start, a.mapping_quality = name, 0, list(m).index(chrom), 60, 60
            a.cigartuples = [(0, 100), (3, 200), (0, 60)] if key else [(0, 100), (3, 460), (0, 100)]
            a.query_sequence = "A" * (160 if key else 200); a.set_tag("ms", 300); a.set_tag("AS", 300); bam.write(a)
            k = (chrom, 161, 360 if key else 620); cnt[k] = cnt.get(k, 0) + 1
    (d / "dump.tsv").write_text("chrom\tstrand\tstart\tend\tplaced\tmax\tplaced_w\tmax_w\n" + "".join(f"{c}\t+\t{s}\t{e}\t{v}\t{v}\t{v}\t{v}\n" for (c, s, e), v in cnt.items()))
    (d / "truth.tsv").write_text("read\tstratum\tlead\tstatus\tpicks\torigin_loci\tloci\n" + "".join(
        f"{name}\tS1\t50\t{status}\t.\t{'' if status == 'outside' else 0}\t0:{origin}:60-{420 if key else 720}:+:AS300\n"
        for name, chrom, key, status, origin in reads))
    (d / "tie.tsv").write_text("read\tfamily\tclass\n")
    return d, score(psi, str(d / "cols.tsv"), str(d / "loci.gff3"), str(d / "dump.tsv"), str(d / "sel.bam"), str(d / "truth.tsv"), str(d / "tie.tsv"), str(d / "o"))

def test_call_holds_on_origin_reads_only(tmp_path):
    d, out = score_fixture(tmp_path, u_key=15, u_other=15, s_key=0, s_other=30, s_outside_key=10)   # 10 reads of a copy the reference lacks, on the junction at S
    assert out.calls[0]["verdict"] == "HOLDS" and out.c1 == (60, 60) and out.outside == 10
    assert (out.calls[0]["sup_u"], out.calls[0]["sup_s"], out.calls[0]["origin_gap"]) == (30, 30, 0.5)
    assert len(open(d / "o.reads.tsv").read().splitlines()) == 71 and len(open(d / "o.calls.tsv").read().splitlines()) == 2
    _, out = score_fixture(tmp_path, u_key=15, u_other=15, s_key=14, s_other=16)                    # the stated gap 0.5 shrinks to 0.03
    assert out.calls[0]["verdict"] == "FAILS"
    assert score_fixture(tmp_path, u_key=9, u_other=9, s_key=0, s_other=30)[1].calls[0]["verdict"] == "UNJUDGED"

def test_reads_of_the_other_copy_count_against(tmp_path):
    _, out = score_fixture(tmp_path, u_key=15, u_other=15, s_key=0, s_other=20, s_from_u=10)        # 10 skip reads at S whose origin is U
    assert out.c1 == (50, 60) and out.calls[0]["sup_s"] == 20 and out.calls[0]["verdict"] == "HOLDS"

def test_loci_bed(tmp_path):
    d, _ = score_fixture(tmp_path, u_key=15, u_other=15, s_key=0, s_other=30)
    main(["loci", str(d / "f.tsv"), str(d / "loci.bed")])
    assert open(d / "loci.bed").read() == "c1\t0\t800\tF\nc2\t0\t800\tF\n"

def test_origin_psi_sums_the_copys_junctions_under_the_key(tmp_path):
    # U uses the key through two acceptors 3 bases apart (the table sums them: 0.5); the origin PSI must sum them too
    d, _ = score_fixture(tmp_path, u_key=15, u_other=15, s_key=0, s_other=30)
    psi = open(d / "f.tsv").read().splitlines()
    u = psi[1].split("\t"); v = list(u); v[HEADER.index("end")] = "357"; u[HEADER.index("placed_w")] = "8"; v[HEADER.index("placed_w")] = "7"
    (d / "f.tsv").write_text("\n".join([psi[0], "\t".join(u), "\t".join(v), psi[2]]) + "\n")
    bam = pysam.AlignmentFile(str(d / "sel.bam")); recs = list(bam)
    with pysam.AlignmentFile(str(d / "sel2.bam"), "wb", header=bam.header) as out:
        for i, a in enumerate(recs):
            if a.reference_name == "c1" and a.cigartuples[1][1] == 200 and i % 2:       # 7 of U's 15 key reads take the shifted acceptor
                a.cigartuples = [(0, 100), (3, 197), (0, 63)]; a.query_sequence = "A" * 163
            out.write(a)
    open(d / "dump.tsv", "a").write("c1\t+\t161\t357\t7\t7\t7\t7\n")
    res = score(str(d / "f.tsv"), str(d / "cols.tsv"), str(d / "loci.gff3"), str(d / "dump.tsv"), str(d / "sel2.bam"), str(d / "truth.tsv"), str(d / "tie.tsv"), str(d / "o2"))
    assert (res.calls[0]["sup_u"], res.calls[0]["key_u"], res.calls[0]["origin_gap"], res.calls[0]["verdict"]) == (30, 15, 0.5, "HOLDS")

def test_development_and_miskeyed_calls_are_left_out(tmp_path):
    d, out = score_fixture(tmp_path, u_key=15, u_other=15, s_key=0, s_other=30)
    args = [str(d / x) for x in ("f.tsv", "cols.tsv", "loci.gff3", "dump.tsv", "sel.bam", "truth.tsv", "tie.tsv")]
    sites = lambda name, verdict: (d / name).write_text("family\tanchor\tcol_anchor\tcol_other\tverdict_1\tunit\nF\tdonor\tC0\tC1\t%s\tc1:+:161-360|c2:+:161\n" % verdict) and str(d / name)
    assert len(score(*args, str(d / "a"), sites=sites("s.tsv", "AGREES")).calls) == 1
    assert score(*args, str(d / "b"), sites=sites("m.tsv", "MISKEYED")).calls == []
    r = score(*args, str(d / "c"), sites=sites("s.tsv", "AGREES"), dev=sites("dev.tsv", "AGREES"))
    assert r.calls == [] and r.c1 == (0, 0)
