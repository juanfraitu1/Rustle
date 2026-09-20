#!/usr/bin/env python3
"""family_def_coverage_v2.py — ~C reordered: PER-READ full-length, applied genome-wide.

Order A (cov_frac of de-tied reads) was net-negative (pruned diverged real paralogs
GBA1/GBAP1, SHLD2). Order B (this): a family edge passes ~C iff a majority of its
cross-mapping molecules align FULL-LENGTH (qaln>=0.9) at BOTH loci -- which keeps diverged
paralogs (full-length, just divergent) and cuts fragment-bridges. Re-run on the 196 and
check the net + that the affected families now resolve correctly.
Run: python bench/family_def_coverage_v2.py
"""
import collections
import json
import os
import sys

import networkx as nx
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import best_gene, GENES_BED, DELTA, DE_MAX, MIN_READS
from family_def_read_filters import dna_homology
from family_def_coverage_reorder import dedup, reads_at   # reuse

BAM = "/home/juanfra/winloci_scratch/GGO_mm.bam"
FAM_TSV = "/home/juanfra/winloci_scratch/validated_families.tsv"
QALN_MIN = 0.90        # a molecule maps "full-length" if >=90% of it aligns at the locus
C_FRAC = 0.50          # ~C: >= this fraction of crossmap molecules must be full-length at both


def edge_qaln(ri, rj):
    """(fraction of crossmap molecules full-length at both, n_detied)."""
    shared = set(ri) & set(rj)
    if not shared:
        return None
    tied = [q for q in shared if abs(ri[q][0] - rj[q][0]) <= DELTA and max(ri[q][0], rj[q][0]) <= DE_MAX]
    if len(tied) < MIN_READS:
        return (0.0, len(tied))
    full = sum(1 for q in shared if min(ri[q][2], rj[q][2]) >= QALN_MIN)
    return (full / len(shared), len(tied))


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

    def is_paralog(ga, gb):
        r = Hd.get((ga, gb) if ga < gb else (gb, ga))
        return r is not None and r.get("id", 0) >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30

    # collect every ~B edge's (full-length fraction, DNA truth) once, then sweep the threshold
    kept_B = 0
    edges_data = []                       # (frac, is_paralog, is_bridge)
    for fi in sorted(fams):
        regs = dedup(fams[fi])
        if len(regs) < 2:
            continue
        kept_B += 1
        reads = {r: reads_at(bam, *r) for r in regs}
        for i in range(len(regs)):
            for j in range(i + 1, len(regs)):
                rc = edge_qaln(reads[regs[i]], reads[regs[j]])
                if rc is None or rc[1] < MIN_READS:
                    continue
                gi, gj = best_gene(by, *regs[i]), best_gene(by, *regs[j])
                par = bool(gi and gj and is_paralog(gi, gj))
                brg = bool(gi and gj and gi != gj and not par and
                           (lambda r: r is None or r.get("id", 0) == 0)(Hd.get((gi, gj) if gi < gj else (gj, gi))))
                edges_data.append((rc[0], par, brg))
    bam.close()

    n_par = sum(1 for _, p, _ in edges_data if p)
    n_brg = sum(1 for _, _, b in edges_data if b)
    print("=== ~C reordered (per-read full-length), genome-wide ===")
    print(f"  multi-copy ~B families: {kept_B}; ~B edges with DNA truth: "
          f"{n_par} PARALOG / {n_brg} bridge / {len(edges_data)-n_par-n_brg} unannotated")
    print(f"\n  threshold sweep (edge passes ~C iff full-length-frac >= T):")
    print(f"  {'T':>5} {'bridge_cut':>10} {'paralog_cut':>11} {'bridge_recall':>13} {'paralog_cost':>12}")
    for T in [0.30, 0.40, 0.50, 0.60, 0.70, 0.80]:
        bc = sum(1 for f, p, b in edges_data if b and f < T)
        pc = sum(1 for f, p, b in edges_data if p and f < T)
        print(f"  {T:>5.2f} {bc:>10} {pc:>11} {bc/max(n_brg,1):>13.2f} {pc/max(n_par,1):>12.2f}")
    pruned_edges_bridge = sum(1 for f, p, b in edges_data if b and f < C_FRAC)
    pruned_edges_paralog = sum(1 for f, p, b in edges_data if p and f < C_FRAC)
    kept_C = pruned = split = 0
    affected = []
    print(f"\n  at T={C_FRAC}: bridge edges cut={pruned_edges_bridge}/{n_brg}, "
          f"paralog edges cut={pruned_edges_paralog}/{n_par} (recall cost)")
    json.dump(dict(qaln_min=QALN_MIN, multicopy_B=kept_B, n_paralog_edges=n_par, n_bridge_edges=n_brg,
                   sweep=[dict(T=T, bridge_cut=sum(1 for f, p, b in edges_data if b and f < T),
                               paralog_cut=sum(1 for f, p, b in edges_data if p and f < T)) for T in
                          [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]]),
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "family_def_coverage_v2.json"), "w"), indent=2)
    print("\n[+] wrote family_def_coverage_v2.json")


if __name__ == "__main__":
    main()
