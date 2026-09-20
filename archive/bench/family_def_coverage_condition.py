#!/usr/bin/env python3
"""family_def_coverage_condition.py — add the FLNC full-length condition ~C to the definition.

Definition becomes a CONJUNCTION of THREE relations:
  ~R  read-confusability   (>=k reads cross-map at de-tied divergence)
  ~B  backbone-homology    (reciprocal coverage of copy models >= tau)
  ~C  full-length coverage (the cross-mapping reads cover >= tau_C of BOTH loci's expressed
      exon-union, not just a shared fragment) -- leverages FLNC = full molecules.

For each ~B-validated family (the 196), dedup members to copy-regions, recompute ~C per
copy-region pair from the multimapper BAM, keep ~C-passing edges, re-derive components, and
report how many families ~C sharpens (prunes as fragment-bridges / splits).
Run: python bench/family_def_coverage_condition.py
"""
import collections
import json
import os
import sys

import networkx as nx
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import DELTA, DE_MAX, MIN_READS

FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
PAD = 2000
TAU_C = 0.20
OVERLAP_MERGE = 0.5


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def covered(iv):
    return sum(e - s for s, e in merge(iv))


def dedup(loci):
    """member loci -> copy-regions (same chrom, reciprocal overlap > OVERLAP_MERGE merged)."""
    by = collections.defaultdict(list)
    for c, s, e in loci:
        by[c].append((s, e))
    regs = []
    for c, it in by.items():
        it.sort()
        cur = None
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
        b = r.get_blocks()
        if b and (r.query_name not in out or de < out[r.query_name][0]):
            out[r.query_name] = (de, b)
    return out


def coverage_C(bam, ri, rj):
    """~C: fraction of each locus's expressed exon-union covered by the DE-TIED crossmap reads."""
    eu_i = covered([list(b) for _, bl in ri.values() for b in bl])
    eu_j = covered([list(b) for _, bl in rj.values() for b in bl])
    if eu_i == 0 or eu_j == 0:
        return None
    xmap = [(ri[q], rj[q]) for q in set(ri) & set(rj)
            if abs(ri[q][0] - rj[q][0]) <= DELTA and max(ri[q][0], rj[q][0]) <= DE_MAX]
    if len(xmap) < MIN_READS:
        return (0.0, len(xmap))   # ~R no longer holds on this BAM
    cov_i = covered([list(b) for (de, bl), _ in xmap for b in bl])
    cov_j = covered([list(b) for _, (de, bl) in xmap for b in bl])
    return (min(cov_i / eu_i, cov_j / eu_j), len(xmap))


def main():
    fams = collections.defaultdict(list)
    for line in open(FAM_TSV).read().splitlines()[1:]:
        fi, lid, c, s, e, nr = line.split("\t")
        fams[int(fi)].append((c, int(s), int(e)))
    bam = pysam.AlignmentFile(BAM, "rb")

    kept_B = kept_C = pruned = split = 0
    sharpened = []
    cov_cache = {}
    for fi in sorted(fams):
        regs = dedup(fams[fi])
        if len(regs) < 2:
            continue
        kept_B += 1
        reads = {r: reads_at(bam, *r) for r in regs}
        G = nx.Graph(); G.add_nodes_from(range(len(regs)))
        edge_cov = []
        for i in range(len(regs)):
            for j in range(i + 1, len(regs)):
                rc = coverage_C(bam, reads[regs[i]], reads[regs[j]])
                if rc is None:
                    continue
                cf, nx_reads = rc
                edge_cov.append(cf)
                if cf >= TAU_C and nx_reads >= MIN_READS:
                    G.add_edge(i, j)
        comps = [c for c in nx.connected_components(G) if len(c) >= 2]
        if not comps:
            pruned += 1
            sharpened.append((fi, len(regs), "PRUNED (all edges fragment-only)",
                              round(max(edge_cov) if edge_cov else 0, 2)))
        else:
            kept_C += len(comps)
            if len(comps) > 1 or sum(len(c) for c in comps) < len(regs):
                split += 1
                sharpened.append((fi, len(regs), f"SPLIT -> {len(comps)} cores, "
                                  f"{len(regs)-sum(len(c) for c in comps)} members dropped",
                                  round(min(edge_cov) if edge_cov else 0, 2)))
    bam.close()

    print("=== ~C (full-length coverage) applied to the ~B-validated families ===")
    print(f"  multi-copy ~B families (deduped >=2 copies) : {kept_B}")
    print(f"  ~C-validated families (cores)               : {kept_C}")
    print(f"  families PRUNED entirely (fragment-bridges)  : {pruned}")
    print(f"  families SHARPENED (split / member dropped)  : {split}")
    print(f"  NET families ~R∩~B  -> ~R∩~B∩~C              : {kept_B} -> {kept_C}")
    print(f"\n  sharpened families (id, #copies, effect, edge cov_frac):")
    for fi, n, eff, cf in sorted(sharpened)[:25]:
        print(f"    fam{fi:>3} ({n}cp) {eff}  [cov_frac={cf}]")
    json.dump(dict(tau_C=TAU_C, multicopy_B=kept_B, validated_C=kept_C,
                   pruned=pruned, sharpened=split,
                   detail=[dict(family=fi, copies=n, effect=eff, cov_frac=cf) for fi, n, eff, cf in sharpened]),
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "family_def_coverage_condition.json"), "w"), indent=2)
    print("\n[+] wrote family_def_coverage_condition.json")


if __name__ == "__main__":
    main()
