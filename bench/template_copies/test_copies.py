import json
import os
import re
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import copies as K  # noqa: E402
from templates import Transcript  # noqa: E402


def paf(q, qlen, strand, chrom, ts, nm, cg, qs=0, qe=None):
    """One PAF line; the target end and block length follow from the CIGAR."""
    ops = re.findall(r"(\d+)([MIDN])", cg)
    te = ts + sum(int(n) for n, op in ops if op in "MDN")
    bl = sum(int(n) for n, op in ops if op in "MID")
    qe = qlen if qe is None else qe
    return "\t".join(map(str, [q, qlen, qs, qe, strand, chrom, 10 ** 9, ts, te, nm, bl, 60, "tp:A:P", "cg:Z:" + cg])) + "\n"


def lrow(src, cls, contig, s, e, short=False):
    return {"source_id": src, "cls": cls, "contig": contig, "start0": str(s), "end": str(e), "iv": [(s, e)],
            "short": short, "seqid_f": 1.0}


def brow(locus, contig, s, e):
    return {"contig": contig, "locus": locus, "iv": [(s, e)]}


def tmpl(tid, contig, exons):
    return Transcript(tid, tid, contig, "+", exons, 5)


def test_metrics_follow_liftoff():
    # 10 template bases unaligned, 5 target-only bases: coverage 150/160, identity 145/(160+5)
    assert K.metrics(160, 145, "100M5D50M") == (150 / 160, 145 / 165)
    assert K.metrics(100, 99, "50M1000N50M") == (1.0, 0.99)       # an intron is neither
    assert K.metrics(100, 94, "50M5I45M") == (0.95, 0.94)         # template-only bases are not aligned


# Liftoff: gene S1 placed on chr1 with copies on chr2, chr3, chr6; gene S2 placed on chr1 (not assembled), copy on chr4.
LOCI = [lrow("S1", "in_place", "chr1", 1000, 2000), lrow("S1", "extra_copy", "chr2", 5000, 6000),
        lrow("S1", "extra_copy", "chr3", 5000, 6000), lrow("S1", "extra_copy", "chr6", 0, 1000),
        lrow("S2", "in_place", "chr1", 50000, 51000), lrow("S2", "extra_copy", "chr4", 100, 1100),
        lrow("S3", "extra_copy", "chr9", 0, 1000, short=True)]
BASE = [brow("t1", "chr1", 1000, 2000), brow("t2", "chr5", 0, 1000), brow("t3", "chr7", 0, 1200)]
TEMPLATES = {"t1": tmpl("t1", "chr1", [(1000, 2000)]), "t2": tmpl("t2", "chr5", [(0, 1000)]),
             "t3": tmpl("t3", "chr7", [(0, 400), (600, 1200)])}


def lines():
    return [paf("t1", 1000, "+", "chr1", 1000, 1000, "1000M"),        # t1 self
            paf("t1", 1000, "+", "chr2", 5000, 990, "1000M"),         # t1 at copy chr2: identity 0.99
            paf("t1", 1000, "+", "chr3", 5000, 900, "1000M"),         # t1 at copy chr3: identity 0.90, below -sc
            paf("t2", 1000, "+", "chr5", 0, 1000, "1000M"),           # t2 self
            paf("t2", 1000, "+", "chr3", 5000, 1000, "1000M"),        # t2 at copy chr3 (another gene's template)
            paf("t2", 1000, "-", "chr1", 50000, 1000, "1000M"),       # t2 at S2's placement (a Liftoff gene)
            paf("t1", 1000, "+", "chr5", 0, 1000, "1000M"),           # t1 at t2's locus (assembled, no Liftoff)
            paf("t3", 1000, "+", "chr7", 0, 1000, "400M200N600M"),    # t3 self (spliced)
            paf("t3", 1000, "-", "chr2", 5000, 980, "1000M"),         # t3 at chr2, other strand: same site as t1's
            paf("t3", 1000, "+", "chr8", 0, 1000, "1000M")]           # t3 unspliced at no annotation: retro-like


def write_run(tmp_path, paf_lines, done=True):
    run = tmp_path / "run"
    (run / "batches").mkdir(parents=True)
    with open(run / "templates.tsv", "w") as fh:
        fh.write("tid\tgene\tcontig\tstrand\texons\treads\n")
        for t in TEMPLATES.values():
            fh.write(f"{t.tid}\t{t.gene}\t{t.contig}\t{t.strand}\t{','.join(f'{s}-{e}' for s, e in t.exons)}\t5\n")
    (run / "batches" / "b0000.fa").write_text(">t1\nA\n")
    (run / "batches" / "b0000.paf").write_text("".join(paf_lines))
    if done:
        (run / "batches" / "b0000.done").write_text("")
    return run


def test_scan_keeps_liftoff_passing_hits_and_counts(tmp_path):
    run = write_run(tmp_path, lines())
    sc = K.scan([str(run / "batches" / "b0000.paf")], TEMPLATES, keep={"t1"})
    assert sorted((h["template"], h["contig"]) for h in sc["projected"]) == [
        ("t1", "chr1"), ("t1", "chr2"), ("t1", "chr5"), ("t2", "chr1"), ("t2", "chr3"), ("t2", "chr5"),
        ("t3", "chr2"), ("t3", "chr7"), ("t3", "chr8")]
    assert len(sc["kept"]) == 4                       # every t1 alignment, whatever its identity
    assert sc["alignments"] == {"t1": 4, "t2": 3, "t3": 3}
    assert sc["self"] == {"t1", "t2", "t3"}


def test_recall_own_any_and_misses(tmp_path):
    run = write_run(tmp_path, lines())
    sc = K.scan([str(run / "batches" / "b0000.paf")], TEMPLATES, keep=K.own_templates(LOCI, BASE, TEMPLATES)[1])
    rec = K.recall(LOCI, BASE, TEMPLATES, sc)
    # S1's copies are reachable (its placement is assembled: t1); S2's is not; S3's copy is short (never scored)
    assert (rec["copies"], rec["reachable"], rec["unreachable_unassembled"], rec["unreachable_unplaced"]) == (4, 3, 1, 0)
    assert (rec["found_own"], rec["found_any"]) == (1, 2)          # chr2 by t1; chr3 only by t2's template
    assert rec["misses"] == {"below_threshold": 1, "capped": 0, "no_alignment": 1}    # chr3 at 0.90; chr6 nothing
    assert rec["by_copy_count"] == {"1": [0, 0], "2-5": [1, 3], ">=6": [0, 0]}


def test_unplaced_source_is_counted_not_scored():
    loci = [lrow("S9", "extra_copy", "chr2", 0, 1000), lrow("S9", "dropped_M1", "chr1", 0, 1000)]
    rec = K.recall(loci, BASE, TEMPLATES, {"projected": [], "kept": [], "alignments": {}, "self": set()})
    assert (rec["copies"], rec["reachable"], rec["unreachable_unplaced"]) == (1, 0, 1)


def test_capped_miss(tmp_path):
    many = [paf("t1", 1000, "+", "chr1", 1000, 1000, "1000M")] + [
        paf("t1", 1000, "+", "chr12", 10000 * i, 1000, "1000M") for i in range(K.CAP - 1)]
    run = write_run(tmp_path, many)
    sc = K.scan([str(run / "batches" / "b0000.paf")], TEMPLATES, keep={"t1"})
    rec = K.recall(LOCI, BASE, TEMPLATES, sc)
    assert rec["misses"] == {"below_threshold": 0, "capped": 3, "no_alignment": 0}


def test_sites_off_own_locus_merge_across_strands_and_classes(tmp_path):
    run = write_run(tmp_path, lines())
    sc = K.scan([str(run / "batches" / "b0000.paf")], TEMPLATES, keep=set())
    st = K.sites(LOCI, BASE, TEMPLATES, sc["projected"])
    # self hits drop out; chr2 (t1 +, t3 -) is one site
    assert st["n"] == 5
    assert st["classes"] == {"liftoff_copy": [2, 0], "liftoff_gene": [1, 0], "assembled": [1, 0], "none": [1, 1]}
    assert st["unexplained"] == 0 and st["templates_with_sites"] == 3


def test_self_rate_gate(tmp_path):
    run = write_run(tmp_path, [paf("t1", 1000, "+", "chr1", 1000, 1000, "1000M")])
    sc = K.scan([str(run / "batches" / "b0000.paf")], TEMPLATES, keep=set())
    assert K.self_rate(sc, TEMPLATES) == pytest.approx(1 / 3)


def test_cli_end_to_end(tmp_path):
    run = write_run(tmp_path, lines())
    lo = tmp_path / "lift" / "human"
    lo.mkdir(parents=True)
    with open(lo / "liftoff_loci.tsv", "w") as fh:
        fh.write("source_id\tcls\tcontig\tstart0\tend\tshort\tsequence_id\texons\texonic_len\n")
        for r in LOCI:
            s, e = r["iv"][0]
            fh.write(f"{r['source_id']}\t{r['cls']}\t{r['contig']}\t{s}\t{e}\t{int(r['short'])}\t1.0\t{s}-{e}\t{e - s}\n")
    gtf = tmp_path / "base.gtf"
    with open(gtf, "w") as fh:
        for b in BASE:
            s, e = b["iv"][0]
            fh.write(f'{b["contig"]}\tx\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "{b["locus"]}"; transcript_id "{b["locus"]}.1";\n')
    out = tmp_path / "res" / "copies.json"
    assert K.main(["--species", "human", "--sample", "toy", "--base", str(gtf), "--run", str(run),
                   "--out", str(out), "--liftoff", str(tmp_path / "lift")]) == 0
    res = json.load(open(out))
    assert res["recall"]["found_own"] == 1 and res["sites"]["n"] == 5
    assert res["c1"] == {"bar": 0.8, "recall_own": 1 / 3, "reachable": 3, "powered": False, "passes": None}
    assert res["c2"]["passes"] is True


def test_cli_refuses_unprojected_batches_and_low_self_rate(tmp_path):
    run = write_run(tmp_path, lines(), done=False)
    with pytest.raises(SystemExit):
        K.batch_pafs(str(run))
    run2 = write_run(tmp_path / "b", [paf("t1", 1000, "+", "chr1", 1000, 1000, "1000M")])
    sc = K.scan(K.batch_pafs(str(run2)), TEMPLATES, keep=set())
    with pytest.raises(SystemExit):
        K.check_self(sc, TEMPLATES)
