#!/usr/bin/env python3
"""
family_def_repeat_mask.py — does masking REPEATS clean the family definition?

Unlike read-coverage / intron / multi-exon filters (all net-harmful, because the
contamination is read-indistinguishable from real paralogy), repeats are a
mechanism-specific signal: spurious cross-locus bridges (OCLN~SEPTIN7) are driven by
shared transposable elements / satellites, while real paralog cross-mapping is over
GENIC (coding) sequence. The NCBI gorilla assembly is already soft-masked (lowercase =
TE/satellite/low-complexity; segmental duplications are NOT masked, so real segdup
paralogs are preserved).

A placement is "repeat-driven" if its ALIGNED (exonic, M/=/X) reference blocks are
mostly lowercase. We drop repeat-driven placements before building the de-tie graph and
score the result against the DNA homology truth with the same good:bad metric as the
rejected filters (good:bad = real-paralog edges lost per junk edge lost; want < 1).

Run: /home/juanfra/miniforge3/bin/python bench/family_def_repeat_mask.py
"""
import collections
import json
import os
import re
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_genomewide import best_gene, ref_span, de_of, GENES_BED, BAM, DELTA, DE_MAX, MIN_READS
from family_def_read_filters import dna_homology

GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
THRESHOLDS = [0.3, 0.5, 0.7]


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


def aligned_blocks(pos0, cigar):
    """list of (start,end) reference blocks consuming query identity (M/=/X)."""
    blocks = []; ref = pos0
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "M=X":
            blocks.append((ref, ref + n)); ref += n
        elif op in "DN":
            ref += n
    return blocks, ref


def repeat_frac(genome, chrom, pos0, cigar):
    blocks, ref_end = aligned_blocks(pos0, cigar)
    aligned = sum(e - s for s, e in blocks)
    if aligned == 0:
        return 0.0
    try:
        seq = genome.fetch(chrom, pos0, ref_end)
    except (KeyError, ValueError):
        return 0.0
    lc = 0
    for s, e in blocks:
        sub = seq[s - pos0:e - pos0]
        lc += sub.count("a") + sub.count("c") + sub.count("g") + sub.count("t")
    return lc / aligned


def scan(by_chrom):
    """mm[qname][gene] = (min_de, repeat_frac_at_that_placement)."""
    genome = pysam.FastaFile(GENOME)
    mm = collections.defaultdict(dict)
    for flag in ("0x100", None):
        args = ["samtools", "view", "-f", "0x100", BAM] if flag else ["samtools", "view", "-F", "0x900", BAM]
        p = subprocess.Popen(args, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.split("\t")
            if len(f) < 9:
                continue
            if flag is None and f[0] not in mm:   # primaries: only multimapper reads
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
                d[g] = (de, repeat_frac(genome, f[2], pos0, f[5]))
        p.wait()
    genome.close()
    return mm


def build_edges(mm, rmax=1.0):
    ev = collections.defaultdict(int)
    for genes in mm.values():
        if len(genes) < 2:
            continue
        items = sorted(genes.items())
        for a in range(len(items)):
            ga, (da, ra) = items[a]
            for b in range(a + 1, len(items)):
                gb, (db, rb) = items[b]
                if abs(da - db) > DELTA or max(da, db) > DE_MAX:
                    continue
                if ra > rmax or rb > rmax:        # drop repeat-driven placements
                    continue
                ev[(ga, gb)] += 1
    return {k for k, n in ev.items() if n >= MIN_READS}


def main():
    by_chrom, g2chrom = load_genes()
    print("[scan] reading GGO.bam + soft-mask (per-placement repeat fraction) ...", flush=True)
    mm = scan(by_chrom)
    H, dna = dna_homology()

    def cat(k):
        r = H.get(k)
        if r is None or r["id"] == 0:
            return "no_homology"
        if r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
            return "DNA_paralog"
        return "sub_bar" if r["id"] >= 0.80 else "weak"

    base = build_edges(mm, 1.0)
    base_real = {k for k in base if cat(k) in ("DNA_paralog", "sub_bar")}
    base_junk = {k for k in base if cat(k) in ("no_homology", "weak")}
    print(f"\n  {'mask threshold':22} {'edges':>6} {'TP':>5} {'junk':>5} {'P':>6} "
          f"{'realLost':>8} {'junkLost':>8} {'good:bad':>9}")
    print(f"  {'baseline (no mask)':22} {len(base):6d} {len(base&dna):5d} {len(base_junk):5d} "
          f"{len(base&dna)/max(len(base),1):6.3f} {0:8d} {0:8d} {'—':>9}")
    rows = []
    for T in THRESHOLDS:
        e = build_edges(mm, T)
        real = {k for k in e if cat(k) in ("DNA_paralog", "sub_bar")}
        junk = {k for k in e if cat(k) in ("no_homology", "weak")}
        rl = len(base_real - real); jl = len(base_junk - junk)
        ratio = rl / jl if jl else float("inf")
        P = len(e & dna) / max(len(e), 1)
        print(f"  {'drop placement rep>'+str(T):22} {len(e):6d} {len(e&dna):5d} {len(junk):5d} "
              f"{P:6.3f} {rl:8d} {jl:8d} {ratio:9.2f}")
        rows.append(dict(threshold=T, edges=len(e), tp=len(e & dna), precision=P,
                         real_lost=rl, junk_lost=jl, good_to_bad=ratio))
    print("  good:bad = real paralog edges lost per junk edge lost (want < 1 — the first filter that would).")

    # is repeat-fraction actually higher on the junk (no-homology) edges' reads?
    print("\n--- mean placement repeat-fraction by edge class (baseline edges) ---")
    cls_fracs = collections.defaultdict(list)
    # recover per-edge mean repeat frac from the reads supporting it
    edge_fracs = collections.defaultdict(list)
    for genes in mm.values():
        items = sorted(genes.items())
        for a in range(len(items)):
            ga, (da, ra) = items[a]
            for b in range(a + 1, len(items)):
                gb, (db, rb) = items[b]
                if abs(da - db) <= DELTA and max(da, db) <= DE_MAX:
                    edge_fracs[(ga, gb)].append(max(ra, rb))
    for k in base:
        m = sum(edge_fracs[k]) / len(edge_fracs[k]) if edge_fracs[k] else 0
        cls_fracs[cat(k)].append(m)
    for c in ("DNA_paralog", "sub_bar", "no_homology", "weak"):
        v = cls_fracs.get(c, [])
        if v:
            print(f"  {c:14} n={len(v):5d}  mean max-placement repeat-frac={sum(v)/len(v):.3f}")

    json.dump(dict(thresholds=rows), open(os.path.join(HERE, "family_def_repeat_mask.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_repeat_mask.json')}")


if __name__ == "__main__":
    main()
