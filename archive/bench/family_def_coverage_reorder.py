#!/usr/bin/env python3
"""family_def_coverage_reorder.py — ~C in a DIFFERENT ORDER (per-read full-length, not
coverage-of-tied-reads).

Order A (family_def_coverage_condition.py) = cov_frac: fraction of the LOCUS covered by the
DE-TIED crossmap reads. The de-tie restricts a diverged paralog to its conserved region ->
low cov_frac -> wrongly pruned GBA1/GBAP1, SHLD2.

Order B (here) = qaln: per-READ, does the molecule align FULL-LENGTH at BOTH loci? A bridge
read maps fragment at the 2nd locus (low qaln); a diverged paralog read maps full-length
(high qaln) regardless of de. Test on the 6 families Order A mishandled: does qaln KEEP the
real paralogs while still cutting the genuine bridge?
Run: python bench/family_def_coverage_reorder.py
"""
import collections
import os
import sys

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import best_gene, GENES_BED, DELTA, DE_MAX, MIN_READS
from family_def_read_filters import dna_homology

BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
PAD = 2000
OVERLAP_MERGE = 0.5
AFFECTED = [24, 64, 124, 138, 150, 153]


def merge(iv):
    iv = sorted(iv); out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def covered(iv):
    return sum(e - s for s, e in merge(iv))


def dedup(loci):
    by = collections.defaultdict(list)
    for c, s, e in loci:
        by[c].append((s, e))
    regs = []
    for c, it in by.items():
        it.sort(); cur = None
        for s, e in it:
            if cur and s < cur[1] and (min(cur[1], e) - s) > OVERLAP_MERGE * min(cur[1] - cur[0], e - s):
                cur = (cur[0], max(cur[1], e)); continue
            if cur:
                regs.append((c, cur[0], cur[1]))
            cur = (s, e)
        if cur:
            regs.append((c, cur[0], cur[1]))
    return regs


def reads_at(bam, chrom, s, e):
    """qname -> (de, blocks, qaln_frac)."""
    out = {}
    try:
        it = bam.fetch(chrom, max(0, s - PAD), e + PAD)
    except (ValueError, KeyError):
        return out
    for r in it:
        if r.is_unmapped or r.is_supplementary:
            continue
        if r.reference_end is None or r.reference_end < s or r.reference_start > e:
            continue
        de = dict(r.get_tags()).get("de")
        if de is None or de > DE_MAX:
            continue
        ql = r.query_length or r.infer_read_length() or 0
        qaln = (r.query_alignment_length or 0) / ql if ql else 0
        b = r.get_blocks()
        if b and (r.query_name not in out or de < out[r.query_name][0]):
            out[r.query_name] = (de, b, qaln)
    return out


def main():
    fams = collections.defaultdict(list)
    for line in open(FAM_TSV).read().splitlines()[1:]:
        fi, lid, c, s, e, nr = line.split("\t")
        fams[int(fi)].append((c, int(s), int(e)))
    by = collections.defaultdict(list)
    for line in open(GENES_BED):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            by[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c in by:
        by[c].sort()
    Hd, _ = dna_homology()
    bam = pysam.AlignmentFile(BAM, "rb")

    def dverdict(genes):
        v = []
        gs = sorted(set(genes) - {None})
        for i in range(len(gs)):
            for j in range(i + 1, len(gs)):
                ga, gb = gs[i], gs[j]
                r = Hd.get((ga, gb) if ga < gb else (gb, ga))
                if r is None or r.get("id", 0) == 0:
                    v.append("bridge")
                elif r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
                    v.append(f"PARALOG{r['id']:.2f}")
                else:
                    v.append(f"sub{r['id']:.2f}")
        return gs, v

    print("=== ~C reorder: cov_frac (Order A, tied-read coverage) vs qaln (Order B, per-read full-length) ===")
    print(f"  {'fam':>5} {'genes':28} {'DNA':16} {'edge':>4} {'covfrac':>8} {'qaln_full%':>11} {'n_tied':>6}")
    for fi in AFFECTED:
        regs = dedup(fams[fi])
        reads = {r: reads_at(bam, *r) for r in regs}
        gns, ver = dverdict([best_gene(by, *r) for r in regs])
        for i in range(len(regs)):
            for j in range(i + 1, len(regs)):
                ri, rj = reads[regs[i]], reads[regs[j]]
                eu_i = covered([list(b) for _, bl, _ in ri.values() for b in bl])
                eu_j = covered([list(b) for _, bl, _ in rj.values() for b in bl])
                if eu_i == 0 or eu_j == 0:
                    continue
                shared = set(ri) & set(rj)
                tied = [q for q in shared if abs(ri[q][0] - rj[q][0]) <= DELTA and max(ri[q][0], rj[q][0]) <= DE_MAX]
                if len(tied) < MIN_READS:
                    continue
                cov_i = covered([list(b) for q in tied for b in ri[q][1]])
                cov_j = covered([list(b) for q in tied for b in rj[q][1]])
                covfrac = min(cov_i / eu_i, cov_j / eu_j)
                # Order B: per-read full-length at BOTH (qaln>=0.9), over ALL shared crossmap reads
                fulllen = [q for q in shared if min(ri[q][2], rj[q][2]) >= 0.90]
                qaln_full = len(fulllen) / max(len(shared), 1)
                print(f"  {fi:>5} {'~'.join(gns)[:28]:28} {','.join(ver)[:16]:16} "
                      f"{i}-{j:>1} {covfrac:>8.2f} {100*qaln_full:>10.0f}% {len(tied):>6}")
    bam.close()
    print("\n  read: if qaln_full% is HIGH for the PARALOG rows and LOW for the bridge row,")
    print("  the per-read full-length order separates them where cov_frac conflated them.")


if __name__ == "__main__":
    main()
