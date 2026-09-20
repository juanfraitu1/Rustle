#!/usr/bin/env python3
"""
family_def_junction_concordance.py — exploit the SPLICED structure: require a
conflicting read to be spliced the SAME WAY at both loci.

Idea (sharper than "is it spliced", immune to UTR divergence unlike cDNA coverage):
a real paralog read recovers the SAME exon-intron architecture at both copies
(homologous junctions => matching exon-length signature); a retrocopy/processed
pseudogene maps the read CONTIGUOUSLY (intronless) or with a different junction
pattern at the copy => signatures DISAGREE. So "junction-concordant" conflicts should
keep true (multi-exon) paralogs and drop retrocopies/repeat-bridges.

signature(read at a locus) = sorted tuple of aligned exon-block lengths (M/=/X runs
split by N introns). Concordant = both multi-exon, same #exons, each length within
tolerance. We build the graph from concordant conflicts only and score it with the
good:bad metric vs the DNA homology truth; we also report the concordance RATE per
edge class (the diagnostic: are no-homology bridges less concordant than real paralogs?).

Run: /home/juanfra/miniforge3/bin/python bench/family_def_junction_concordance.py
"""
import collections
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_genomewide import best_gene, ref_span, de_of, GENES_BED, BAM, DELTA, DE_MAX, MIN_READS
from family_def_read_filters import dna_homology

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
ABS_TOL = 20      # bp
REL_TOL = 0.15


def load_genes():
    by_chrom = collections.defaultdict(list)
    with open(GENES_BED) as f:
        for line in f:
            c, s, e, name = line.rstrip("\n").split("\t")
            by_chrom[c].append((int(s), int(e), name))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def signature(cigar):
    """sorted tuple of aligned exon-block lengths (M/=/X runs split by N)."""
    segs = []; cur = 0
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "M=X":
            cur += n
        elif op == "N":
            if cur > 0:
                segs.append(cur); cur = 0
    if cur > 0:
        segs.append(cur)
    return tuple(sorted(segs))


def concordant(a, b):
    if len(a) < 2 or len(b) < 2 or len(a) != len(b):
        return False
    for x, y in zip(a, b):
        if abs(x - y) > max(ABS_TOL, REL_TOL * max(x, y)):
            return False
    return True


def scan(by_chrom):
    mm = collections.defaultdict(dict)
    for flag in ("0x100", None):
        args = ["samtools", "view", "-f", "0x100", BAM] if flag else ["samtools", "view", "-F", "0x900", BAM]
        p = subprocess.Popen(args, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.split("\t")
            if len(f) < 9:
                continue
            if flag is None and f[0] not in mm:
                continue
            de = de_of(f[11:])
            if de is None or de > DE_MAX:
                continue
            pos0 = int(f[3]) - 1
            g = best_gene(by_chrom, f[2], pos0, pos0 + ref_span(f[5]))
            if g is None:
                continue
            d = mm[f[0]]
            if g not in d or de < d[g][0]:
                d[g] = (de, signature(f[5]))
        p.wait()
    return mm


def main():
    by_chrom = load_genes()
    print("[scan] de-tie graph over genes, storing per-placement exon signature ...", flush=True)
    mm = scan(by_chrom)
    H, dna = dna_homology()

    def cat(k):
        r = H.get(k)
        if r is None or r["id"] == 0:
            return "no_homology"
        if r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
            return "DNA_paralog"
        return "sub_bar" if r["id"] >= 0.80 else "weak"

    # baseline edges and concordant edges + per-conflict concordance bookkeeping
    ev_all = collections.defaultdict(int)
    ev_conc = collections.defaultdict(int)
    conc_hits = collections.defaultdict(lambda: [0, 0])  # pair -> [concordant, total]
    for genes in mm.values():
        items = sorted(genes.items())
        for a in range(len(items)):
            ga, (da, sa) = items[a]
            for b in range(a + 1, len(items)):
                gb, (db, sb) = items[b]
                if abs(da - db) > DELTA or max(da, db) > DE_MAX:
                    continue
                ev_all[(ga, gb)] += 1
                conc_hits[(ga, gb)][1] += 1
                if concordant(sa, sb):
                    ev_conc[(ga, gb)] += 1
                    conc_hits[(ga, gb)][0] += 1
    base = {k for k, n in ev_all.items() if n >= MIN_READS}
    strict = {k for k, n in ev_conc.items() if n >= MIN_READS}

    base_real = {k for k in base if cat(k) in ("DNA_paralog", "sub_bar")}
    base_junk = {k for k in base if cat(k) in ("no_homology", "weak")}
    real2 = {k for k in strict if cat(k) in ("DNA_paralog", "sub_bar")}
    junk2 = {k for k in strict if cat(k) in ("no_homology", "weak")}
    rl, jl = len(base_real - real2), len(base_junk - junk2)
    print(f"\n  baseline edges        : {len(base)}  (TP={len(base&dna)} P={len(base&dna)/max(len(base),1):.3f})")
    print(f"  junction-concordant   : {len(strict)}  (TP={len(strict&dna)} P={len(strict&dna)/max(len(strict),1):.3f})")
    print(f"  removed real={rl}  junk={jl}  good:bad={rl/jl if jl else float('inf'):.2f}  "
          f"(want < 1)")

    print("\n--- diagnostic: mean per-edge junction-concordance RATE by class (baseline edges) ---")
    by_cls = collections.defaultdict(list)
    for k in base:
        c, t = conc_hits[k]
        by_cls[cat(k)].append(c / t if t else 0)
    for c in ("DNA_paralog", "sub_bar", "no_homology", "weak"):
        v = by_cls.get(c, [])
        if v:
            print(f"  {c:14} n={len(v):5d}  mean concordance-rate={sum(v)/len(v):.3f}")

    json.dump(dict(baseline_edges=len(base), concordant_edges=len(strict),
                   tp_base=len(base & dna), tp_strict=len(strict & dna),
                   precision_base=len(base & dna) / max(len(base), 1),
                   precision_strict=len(strict & dna) / max(len(strict), 1),
                   real_lost=rl, junk_lost=jl,
                   good_to_bad=rl / jl if jl else None),
              open(os.path.join(HERE, "family_def_junction_concordance.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_junction_concordance.json')}")


if __name__ == "__main__":
    main()
