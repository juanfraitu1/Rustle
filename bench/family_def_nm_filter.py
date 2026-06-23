#!/usr/bin/env python3
"""family_def_nm_filter.py — (b) can an NM read-level condition remove bridges BEFORE ~B,
and at what cost to real copies?

de does not separate ~B-bridge from ~B-copy cross-mappings (0.9x); NM-rate does (2.6x,
family_def_bam_signals.md). Here we test NM as an added ~R criterion: form ~R edges
(de-tie quorum), classify each edge ~B-copy vs ~B-bridge (the ground truth for this test),
then ask how well an NM-rate cut separates the two:
  - EDGE level: sweep T_nm; an edge is "NM-bridge" if its reads' median NM-rate > T_nm.
    Report bridge recall (frac of ~B-bridges flagged) vs copy FP (frac of ~B-copies flagged).
  - READ level: require NM-rate <= T_nm for a read to count toward the k=3 tie quorum; report
    how many ~B-bridge vs ~B-copy edges fall below quorum (i.e. dissolve) after the filter.
Aggregated over chromosomes. Run: python bench/family_def_nm_filter.py NC_073240.2 NC_073244.2 NC_073248.2
"""
import collections
import json
import os
import sys
import statistics as st

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import DELTA, DE_MAX, MIN_READS
from family_def_newbam_validate import load_denovo, build_model, recip_cov, load_models, GENOME, OLD, COV_MIN
from family_def_bam_signals import scan_sig

HERE = os.path.dirname(os.path.abspath(__file__))


def edges_with_reads(mm):
    ev = collections.defaultdict(list)
    for genes in mm.values():
        if len(genes) < 2:
            continue
        items = sorted(genes.items())
        for a in range(len(items)):
            ga, pa = items[a]
            for b in range(a + 1, len(items)):
                gb, pb = items[b]
                if abs(pa[0] - pb[0]) <= DELTA and max(pa[0], pb[0]) <= DE_MAX:
                    ev[(ga, gb)].append((pa, pb))   # pa=(de,nm,rl,aln)
    return ev


def main():
    chroms = sys.argv[1:] or ["NC_073240.2", "NC_073244.2", "NC_073248.2"]
    models = load_models()
    bam_old = pysam.AlignmentFile(OLD, "rb"); genome = pysam.FastaFile(GENOME)
    all_edges = []   # (is_copy, [read NM-rates], nreads)
    for region in chroms:
        by_chrom, info = load_denovo({region})
        mm = scan_sig(by_chrom, region)
        ev = edges_with_reads(mm)
        loci = {g for (ga, gb) in ev for g in (ga, gb)}
        for v in loci:
            if not models.get(v) and v in info:
                c, s, e, _ = info[v]
                models[v] = build_model(bam_old, genome, c, s, e)
        for (ga, gb), recs in ev.items():
            if len(recs) < MIN_READS or not models.get(ga) or not models.get(gb):
                continue
            is_copy = recip_cov(models[ga], models[gb]) >= COV_MIN
            nmrates = [(pa[1] / max(pa[3], 1) + pb[1] / max(pb[3], 1)) / 2 for (pa, pb) in recs]
            nmdiffs = [abs(pa[1] / max(pa[3], 1) - pb[1] / max(pb[3], 1)) for (pa, pb) in recs]
            all_edges.append((is_copy, nmrates, nmdiffs, len(recs)))
    bam_old.close(); genome.close()

    copies = [e for e in all_edges if e[0]]
    bridges = [e for e in all_edges if not e[0]]
    print(f"=== (b) NM as ~R filter — {len(copies)} ~B-copy edges, {len(bridges)} ~B-bridge edges ===")

    # EDGE-level sweep: edge NM signal = median read NM-rate
    print("\n  EDGE-level (flag edge as bridge if median read NM-rate > T):")
    print(f"  {'T_nm':>7} {'bridge_recall':>13} {'copy_FP':>9} {'kept_copies':>12} {'killed_bridges':>15}")
    best = None
    for T in [0.0015, 0.0020, 0.0025, 0.0030, 0.0040, 0.0050]:
        br = sum(1 for e in bridges if st.median(e[1]) > T)
        cf = sum(1 for e in copies if st.median(e[1]) > T)
        recall = br / max(len(bridges), 1)
        fp = cf / max(len(copies), 1)
        print(f"  {T:>7.4f} {recall:>13.2f} {fp:>9.2f} {len(copies)-cf:>12} {br:>15}")
        # operating point: maximize bridge recall at copy FP <= 0.05
        if fp <= 0.05 and (best is None or recall > best[1]):
            best = (T, recall, fp, br, cf)

    # NM-TIE: |NM_rate_i - NM_rate_j| per read (the "DE plus NM" reading: read must be equally
    # well-explained by BOTH copies). Edge NM-tie signal = median read |diff|.
    print("\n  NM-TIE edge-level (flag bridge if median read |NM_rate_i - NM_rate_j| > T):")
    print(f"  {'T_tie':>7} {'bridge_recall':>13} {'copy_FP':>9} {'kept_copies':>12} {'killed_bridges':>15}")
    for T in [0.0010, 0.0015, 0.0020, 0.0025, 0.0030, 0.0040]:
        br = sum(1 for e in bridges if st.median(e[2]) > T)
        cf = sum(1 for e in copies if st.median(e[2]) > T)
        print(f"  {T:>7.4f} {br/max(len(bridges),1):>13.2f} {cf/max(len(copies),1):>9.2f} "
              f"{len(copies)-cf:>12} {br:>15}")

    # READ-level: drop reads with NM-rate > T from the quorum; does the edge keep k>=3?
    print("\n  READ-level (count only reads with NM-rate <= T toward the k=3 quorum):")
    print(f"  {'T_nm':>7} {'bridges_dissolved':>18} {'copies_dissolved':>17}")
    for T in [0.0020, 0.0025, 0.0030]:
        bd = sum(1 for e in bridges if sum(1 for r in e[1] if r <= T) < MIN_READS)
        cd = sum(1 for e in copies if sum(1 for r in e[1] if r <= T) < MIN_READS)
        print(f"  {T:>7.4f} {bd:>6}/{len(bridges):<3} ({100*bd/max(len(bridges),1):>4.0f}%) "
              f"      {cd:>6}/{len(copies):<3} ({100*cd/max(len(copies),1):>4.0f}%)")

    out = dict(n_copy=len(copies), n_bridge=len(bridges),
               best_edge_op=dict(zip(["T_nm", "bridge_recall", "copy_fp", "killed_bridges", "killed_copies"],
                                     best)) if best else None)
    json.dump(out, open(os.path.join(HERE, "family_def_nm_filter.json"), "w"), indent=2)
    if best:
        print(f"\n  BEST edge operating point (copy-FP<=5%): T_nm={best[0]:.4f} -> removes "
              f"{best[3]}/{len(bridges)} bridges ({100*best[1]:.0f}%) at {best[4]} copy false-drops "
              f"({100*best[2]:.0f}%).")
    print("\n[+] wrote family_def_nm_filter.json")


if __name__ == "__main__":
    main()
