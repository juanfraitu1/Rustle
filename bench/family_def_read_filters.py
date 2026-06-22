#!/usr/bin/env python3
"""
family_def_read_filters.py — does adding two READ-level requirements to the de-tie
conflict predicate sharpen the family definition?
    (i)  high query coverage : each conflicting read aligns over >= QCOV of its OWN
         length at BOTH placements (not a partial / clipped / local hit), and
    (ii) >= 1 intron         : the read is SPLICED (has an N in its CIGAR) at BOTH
         placements -> excludes intronless processed-pseudogene / retrocopy / repeat
         cross-mapping, which is the dominant RNA-only "false edge" mode.
(Gene-level multi-exon does NOT work: OCLN/SEPTIN7 and the domain-sharers are all
multi-exon genes; the retro/repeat cross-map is a per-read property.)

Compares the BASELINE de-tie graph to the FILTERED one, re-derives precision against
the DNA cDNA-homology truth (rep_ava.tsv), and audits exactly which edges are removed
(homology-supported real paralogs = a COST; no-homology cross-chrom bridges /
retrocopies = the intended GAIN).

Run: /home/juanfra/miniforge3/bin/python bench/family_def_read_filters.py
"""
import collections
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_genomewide import best_gene, GENES_BED, BAM, DELTA, DE_MAX, MIN_READS

AVA = "/home/juanfra/winloci_scratch/rep_ava.tsv"
EXONCOUNT = "/home/juanfra/winloci_scratch/gene_exoncount.json"
QCOV = 0.80          # high query-coverage floor (per placement)
MIN_INTRON = 1       # require the read to be spliced at the placement
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def cig_stats(cigar):
    aln_q = full_q = ref = nintron = 0
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "MI=X":
            aln_q += n; full_q += n
        elif op in "SH":
            full_q += n
        if op in "MDN=X":
            ref += n
        if op == "N":
            nintron += 1
    qcov = aln_q / full_q if full_q else 0.0
    return ref, qcov, nintron


def de_of(fields):
    for t in fields:
        if t.startswith("de:f:"):
            return float(t[5:])
    return None


def load_genes():
    by_chrom = collections.defaultdict(list)
    g2chrom = {}
    with open(GENES_BED) as f:
        for line in f:
            c, s, e, name = line.rstrip("\n").split("\t")
            by_chrom[c].append((int(s), int(e), name))
            g2chrom.setdefault(name, c)
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom, g2chrom


def scan(by_chrom):
    """per read -> {gene: (min_de, qcov_at_that_de, nintron_at_that_de)}."""
    mm = collections.defaultdict(dict)
    for flag in ("0x100", None):  # secondaries, then primaries
        args = ["samtools", "view", "-f", "0x100", BAM] if flag else ["samtools", "view", "-F", "0x900", BAM]
        p = subprocess.Popen(args, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.split("\t")
            if len(f) < 9:
                continue
            de = de_of(f[11:])
            if de is None or de > DE_MAX:
                continue
            ref, qcov, nintron = cig_stats(f[5])
            g = best_gene(by_chrom, f[2], int(f[3]) - 1, int(f[3]) - 1 + ref)
            if g is None:
                continue
            d = mm[f[0]]
            if g not in d or de < d[g][0]:
                d[g] = (de, qcov, nintron)
        p.wait()
    return mm


def build_edges(mm, qcov_min=0.0, intron_min=0):
    ev = collections.defaultdict(int)
    for genes in mm.values():
        if len(genes) < 2:
            continue
        items = sorted(genes.items())
        for a in range(len(items)):
            ga, (da, qa, ia) = items[a]
            for b in range(a + 1, len(items)):
                gb, (db, qb, ib) = items[b]
                if abs(da - db) > DELTA or max(da, db) > DE_MAX:
                    continue
                if not (qa >= qcov_min and qb >= qcov_min and ia >= intron_min and ib >= intron_min):
                    continue
                ev[(ga, gb)] += 1
    return {k for k, n in ev.items() if n >= MIN_READS}


def dna_homology():
    H = {}
    with open(AVA) as f:
        for line in f:
            q, t, ql, tl, m, bl, qs, qe, strand, mapq, de = line.rstrip("\n").split("\t")
            if q == t:
                continue
            ql = int(ql); m = int(m); bl = int(bl); qs = int(qs); qe = int(qe)
            if bl == 0 or ql == 0:
                continue
            ident = m / bl; cov = (qe - qs) / ql
            a, b = (q, t) if q < t else (t, q)
            r = H.setdefault((a, b), {"id": 0.0, "cov_a": 0.0, "cov_b": 0.0})
            r["id"] = max(r["id"], ident)
            if q == a: r["cov_a"] = max(r["cov_a"], cov)
            else: r["cov_b"] = max(r["cov_b"], cov)
    dna = {k for k, r in H.items() if r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30}
    return H, dna


def main():
    by_chrom, g2chrom = load_genes()
    print("[scan] reading GGO.bam (tracking per-read qcov + intron) ...", flush=True)
    mm = scan(by_chrom)
    H, dna = dna_homology()

    def homol_cat(k):
        r = H.get(k)
        if r is None or r["id"] == 0:
            return "no_homology"
        if r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
            return "DNA_paralog"
        if r["id"] >= 0.80:
            return "sub_bar_homology"
        return "weak"

    base = build_edges(mm)
    base_real = {k for k in base if homol_cat(k) in ("DNA_paralog", "sub_bar_homology")}  # genuine
    base_junk = {k for k in base if homol_cat(k) in ("no_homology", "weak")}  # true FP candidates

    CONFIGS = [
        ("baseline (no filter)", 0.0, 0),
        (">=1 intron, both", 0.0, 1),
        ("qcov>=0.5, both", 0.5, 0),
        ("qcov>=0.8, both", 0.8, 0),
        ("qcov>=0.8 & >=1 intron", 0.8, 1),
        ("qcov>=0.5 & >=1 intron", 0.5, 1),
    ]
    print(f"\n  {'filter':28} {'edges':>6} {'TP':>5} {'realFP*':>7} {'P':>6} {'effP':>6} "
          f"{'realLost':>8} {'junkLost':>8} {'good:bad':>9}")
    rows = []
    for name, qc, im in CONFIGS:
        e = build_edges(mm, qc, im)
        tp = len(e & dna)
        real = {k for k in e if homol_cat(k) in ("DNA_paralog", "sub_bar_homology")}
        junk = {k for k in e if homol_cat(k) in ("no_homology", "weak")}
        P = tp / max(len(e), 1)
        effP = len(real) / max(len(e), 1)
        real_lost = len(base_real - real)
        junk_lost = len(base_junk - junk)
        ratio = real_lost / junk_lost if junk_lost else float("inf")
        print(f"  {name:28} {len(e):6d} {tp:5d} {len(junk):7d} {P:6.3f} {effP:6.3f} "
              f"{real_lost:8d} {junk_lost:8d} {ratio:9.2f}")
        rows.append(dict(filter=name, qcov=qc, intron=im, edges=len(e), tp=tp,
                         precision=P, eff_precision=effP, real_lost=real_lost,
                         junk_lost=junk_lost, good_to_bad_removed=ratio))
    print("  * realFP = no-homology/weak edges (the true precision targets); "
          "good:bad = real paralog edges lost per junk edge lost (want <1; >1 = net-harmful).")

    print("\n--- named cases under (qcov>=0.8 & >=1 intron) ---")
    strict = build_edges(mm, 0.8, 1)
    for a, b in [("RABL2A", "RABL2B"), ("OCLN", "SEPTIN7"), ("BCAS4", "CCDC30"),
                 ("TBC1D1", "LOC134756953"), ("RABL2B", "LOC134756389")]:
        k = (a, b) if a < b else (b, a)
        tag = "KEPT" if k in strict else ("REMOVED" if k in base else "not-an-edge")
        print(f"  {a:10}~{b:14} {tag}  (homology={homol_cat(k)})")

    json.dump(dict(configs=rows), open(os.path.join(HERE, "family_def_read_filters.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_read_filters.json')}")


if __name__ == "__main__":
    main()
