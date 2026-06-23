#!/usr/bin/env python3
"""family_def_masked_b.py — lever (2): build ~B on REPEAT-MASKED copy models S(v).

The 9 false survivors pass ~B because S(v) (read-block exon-union) includes intronic/repeat
sequence: two genes sharing a TE/satellite align >=30% despite ZERO spliced homology. The
NCBI gorilla genome is soft-masked (lowercase = TE/satellite/low-complexity; segmental dups
stay UPPERCASE, so real recent paralogs are uppercase). So we rebuild S(v) keeping only
UPPERCASE (non-repeat) bases and recompute ~B. Prediction: false survivors (share a repeat)
drop below tau; true paralogs (share unique genic sequence) survive.
Run: python bench/family_def_masked_b.py NC_073240.2 NC_073244.2 NC_073248.2
"""
import collections
import json
import os
import subprocess
import sys

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import GENES_BED, DELTA, DE_MAX, MIN_READS
from family_def_feature_discovery import scan_feats
from family_def_newbam_validate import GENOME, NEW, TMP, COV_MIN, PAD, MAX_READS, MIN_LEN, merge
from family_def_read_filters import dna_homology

HERE = os.path.dirname(os.path.abspath(__file__))


def load_genes(chroms):
    by, coord = collections.defaultdict(list), {}
    with open(GENES_BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] not in chroms:
                continue
            by[p[0]].append((int(p[1]), int(p[2]), p[3]))
            coord[p[3]] = (p[0], int(p[1]), int(p[2]))
    for c in by:
        by[c].sort()
    return by, coord


def _subtract(A, B):
    """merged interval lists A minus B."""
    res = []
    for s, e in A:
        cur = s
        for bs, be in B:
            if be <= cur or bs >= e:
                continue
            if bs > cur:
                res.append((cur, min(bs, e)))
            cur = max(cur, be)
            if cur >= e:
                break
        if cur < e:
            res.append((cur, e))
    return res


def build_model(bam, genome, chrom, start, end, mode):
    """exon-union of primary reads.
    mode='raw'    : current S(v) (read blocks, uppercased)
    mode='mask'   : keep only UPPERCASE (non-repeat) bases
    mode='exonic' : exon-union MINUS any region some read splices over (strict exonic; drops
                    retained-intron / intronic sequence that creates spurious ~B)."""
    exon_iv, intron_iv, n = [], [], 0
    try:
        it = bam.fetch(chrom, max(0, start - PAD), end + PAD)
    except (ValueError, KeyError):
        return ""
    for r in it:
        if r.is_unmapped or r.is_supplementary or r.is_secondary:
            continue
        if r.reference_end is None or r.reference_end < start or r.reference_start > end:
            continue
        pos = r.reference_start
        for op, length in (r.cigartuples or []):
            if op in (0, 2, 7, 8):          # M,D,=,X  -> exon (ref-consuming, aligned)
                exon_iv.append((pos, pos + length)); pos += length
            elif op == 3:                    # N -> intron
                intron_iv.append((pos, pos + length)); pos += length
        n += 1
        if n >= MAX_READS:
            break
    if not exon_iv:
        return ""
    ex = merge([list(x) for x in exon_iv])
    if mode == "exonic":
        ex = _subtract(ex, merge([list(x) for x in intron_iv]) if intron_iv else [])
    parts = []
    for s, e in ex:
        seq = genome.fetch(chrom, s, e)
        parts.append("".join(c for c in seq if c.isupper()) if mode == "mask" else seq.upper())
    seq = "".join(parts)
    return seq if len(seq) >= MIN_LEN else ""


def recip(a, b):
    if not a or not b:
        return 0.0
    os.makedirs(TMP, exist_ok=True)
    open(f"{TMP}/a.fa", "w").write(f">a\n{a}\n"); open(f"{TMP}/b.fa", "w").write(f">b\n{b}\n")
    out = subprocess.run(f"minimap2 -c -x asm20 {TMP}/a.fa {TMP}/b.fa 2>/dev/null",
                         shell=True, capture_output=True, text=True).stdout
    best = 0.0
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 11:
            continue
        ql, qs, qe, tl, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8])
        best = max(best, min((qe - qs) / max(ql, 1), (te - ts) / max(tl, 1)))
    return best


def main():
    chroms = sys.argv[1:] or ["NC_073240.2", "NC_073244.2", "NC_073248.2"]
    by_chrom, coord = load_genes(set(chroms))
    Hd, _ = dna_homology()
    bam = pysam.AlignmentFile(NEW, "rb"); genome = pysam.FastaFile(GENOME)
    cache = {}

    def models(g):
        if g not in cache:
            c, s, e = coord[g]
            cache[g] = (build_model(bam, genome, c, s, e, "raw"),
                        build_model(bam, genome, c, s, e, "mask"),
                        build_model(bam, genome, c, s, e, "exonic"))
        return cache[g]

    def dna_label(ga, gb):
        r = Hd.get((ga, gb) if ga < gb else (gb, ga))
        if r is None or r.get("id", 0) == 0:
            return 1
        if r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
            return 0
        return None

    edges = {}
    for region in chroms:
        mm = scan_feats(by_chrom, region)
        ev = collections.defaultdict(int)
        for genes in mm.values():
            if len(genes) < 2:
                continue
            it = sorted(genes.items())
            for i in range(len(it)):
                for j in range(i + 1, len(it)):
                    ga, pa = it[i]; gb, pb = it[j]
                    if ga != gb and abs(pa["de"] - pb["de"]) <= DELTA and max(pa["de"], pb["de"]) <= DE_MAX:
                        ev[(ga, gb)] += 1
        for (ga, gb), n in ev.items():
            if n < MIN_READS or ga not in coord or gb not in coord:
                continue
            lab = dna_label(ga, gb)
            if lab is None:
                continue
            ra, ma, ea = models(ga); rb, mb, eb = models(gb)
            edges[(ga, gb)] = (lab, recip(ra, rb), recip(ma, mb), recip(ea, eb))
    bam.close(); genome.close()

    surv = {k: v for k, v in edges.items() if v[1] >= COV_MIN}   # survivors on RAW ~B
    true = [(k, v) for k, v in surv.items() if v[0] == 0]
    false = [(k, v) for k, v in surv.items() if v[0] == 1]
    print(f"=== lever (2): cleaner S(v) — {len(surv)} RAW ~B survivors "
          f"({len(true)} true paralog / {len(false)} false survivor) ===")
    for idx, name in [(2, "REPEAT-MASKED (uppercase only)"), (3, "STRICT-EXONIC (minus introns)")]:
        f_drop = sum(1 for _, v in false if v[idx] < COV_MIN)
        t_lost = sum(1 for _, v in true if v[idx] < COV_MIN)
        print(f"\n  ~B on {name}:")
        print(f"    FALSE survivors DROPPED: {f_drop}/{len(false)}  ({100*f_drop/max(len(false),1):.0f}%)  <- precision gain")
        print(f"    TRUE paralogs LOST     : {t_lost}/{len(true)}  ({100*t_lost/max(len(true),1):.0f}%)  <- recall cost")
        print(f"    net (false dropped - true lost) = {f_drop - t_lost:+d}")
    print(f"\n  the 9 false survivors  (raw -> masked -> exonic recip-cov):")
    for (ga, gb), v in sorted(false, key=lambda x: x[1][3]):
        print(f"    {ga}~{gb:22} {v[1]:.3f} -> {v[2]:.3f} -> {v[3]:.3f}"
              f"   {'EXONIC DROPS ✓' if v[3] < COV_MIN else 'still survives'}")

    json.dump(dict(survivors=len(surv), n_true=len(true), n_false=len(false),
                   detail=[dict(ga=k[0], gb=k[1], dna_paralog=(v[0] == 0),
                                recip_raw=round(v[1], 3), recip_masked=round(v[2], 3),
                                recip_exonic=round(v[3], 3)) for k, v in surv.items()]),
              open(os.path.join(HERE, "family_def_masked_b.json"), "w"), indent=2)
    print("\n[+] wrote family_def_masked_b.json")


if __name__ == "__main__":
    main()
