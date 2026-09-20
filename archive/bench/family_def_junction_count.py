#!/usr/bin/env python3
"""
family_def_junction_count.py — refined angle A2.

The panel probe showed intron-LENGTH concordance (Carch) is the WRONG invariant: real
paralogs preserve intron COUNT / exon STRUCTURE but NOT intron lengths (introns drift in
length after duplication; RABL2A vs RABL2B = same 9 introns, very different lengths). And
"both spliced" (Cspliced) keeps everything (domain-sharers are multi-exon too).

The duplication-preserved, locus-independent invariant is the INTRON COUNT (number of
junctions). Real whole-gene-duplication copies have COMMENSURATE junction counts; retrocopies
have 0; many domain-bridges join genes of very different exon counts. We test a ladder of
count criteria and report the panel scorecard + genome-wide mega-family breakup, plus a
ROC-style sweep of the count-ratio threshold to find a principled operating point (and the
intron-count AUC separating real-panel from domain/retro counter pairs).

Run: /home/juanfra/miniforge3/bin/python bench/family_def_junction_count.py
"""
import collections
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology
from family_def_junction_genomewide import components, ID_MIN, COV_MIN

INDEX = "/home/juanfra/winloci_scratch/gene_intron_index.json"

PANEL_REAL = {"RABL2": ["RABL2A", "RABL2B"],
              "APOBEC3": ["APOBEC3D", "APOBEC3F"],
              "RFPL": ["RFPL1", "RFPL2", "RFPL3"]}
REAL_PAIRS = [("RABL2A", "RABL2B"), ("APOBEC3D", "APOBEC3F"),
              ("RFPL1", "RFPL2"), ("RFPL2", "RFPL3"), ("RFPL1", "RFPL3")]
COUNTER_DOMAIN = [("CASP8", "FLACC1"), ("ASDURF", "ASNSD1"),
                  ("CREB1", "METTL21A"), ("GPR39", "LYPD1")]


def count_crit(ra, rb, kind):
    """kind:
       both_spliced : na>=1 and nb>=1
       ratio        : both spliced AND min(na,nb)/max(na,nb) >= R  (we sweep R)
    """
    if ra is None or rb is None:
        return False
    na, nb = ra["n_intron"], rb["n_intron"]
    if na < 1 or nb < 1:
        return False
    return True  # placeholder; ratio handled separately


def ratio_keep(ra, rb, R):
    if ra is None or rb is None:
        return False
    na, nb = ra["n_intron"], rb["n_intron"]
    if na < 1 or nb < 1:
        return False
    return min(na, nb) / max(na, nb) >= R


def main():
    H, dna = dna_homology()
    idx = json.load(open(INDEX))
    base_edges = sorted(dna)
    base_comp = components(base_edges)
    base_max = max(len(c) for c in base_comp)
    mega = max(base_comp, key=len)
    mega_edges = {(a, b) for (a, b) in base_edges if a in mega and b in mega}
    print(f"[baseline] edges={len(base_edges)} families={len(base_comp)} maxcomp={base_max} "
          f"mega={len(mega)} megaEdges={len(mega_edges)}")

    # ---- intron-count AUC: panel real vs counter (domain+retro) ----
    # feature = min/max intron-count ratio (1.0 = identical architecture count)
    def feat(a, b):
        ra, rb = idx.get(a), idx.get(b)
        if ra is None or rb is None:
            return None
        na, nb = ra["n_intron"], rb["n_intron"]
        if max(na, nb) == 0:
            return 0.0
        return min(na, nb) / max(na, nb)
    pos = [feat(a, b) for a, b in REAL_PAIRS]
    # add retro pairs (intronless partner) as negatives
    retro_neg = []
    for a in ("EEF1A1", "CNN2"):
        # best intronless partner
        best = None
        for (x, y), r in H.items():
            if r["id"] < 0.90 or max(r["cov_a"], r["cov_b"]) < 0.30:
                continue
            p = y if x == a else (x if y == a else None)
            if p is None:
                continue
            rp = idx.get(p)
            if rp and rp["n_intron"] == 0:
                retro_neg.append(feat(a, p)); break
    neg = [feat(a, b) for a, b in COUNTER_DOMAIN] + retro_neg
    pos = [p for p in pos if p is not None]; neg = [n for n in neg if n is not None]
    # AUC (ratio higher => more "real")
    auc = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))
    print(f"\n[panel] intron-count-ratio  REAL={sorted(round(p,2) for p in pos)}  "
          f"COUNTER={sorted(round(n,2) for n in neg)}")
    print(f"[panel] AUC(intron-count-ratio, real>counter) = {auc:.3f}  "
          f"(n_real={len(pos)} n_counter={len(neg)})")

    # ---- criteria ladder ----
    def panel_conn(edge_set):
        comps = components(edge_set); g2c = {}
        for i, c in enumerate(comps):
            for g in c:
                g2c[g] = i
        out = {}
        for fam, mem in PANEL_REAL.items():
            pres = [m for m in mem if m in g2c]
            out[fam] = "CONN" if len(pres) >= 2 and len({g2c[m] for m in pres}) == 1 else \
                       ("SPLIT" if len(pres) >= 2 else "N/A")
        return out

    print("\n=== INTRON-COUNT CRITERIA LADDER ===")
    print(f"  {'mode':14} {'edges':>6} {'fams':>5} {'maxC':>5} {'megaDrop':>8} "
          f"{'panel':>20} {'realKept':>8} {'domKept':>7}")
    modes = [("both_spliced", None),
             ("ratio>=0.50", 0.50), ("ratio>=0.60", 0.60),
             ("ratio>=0.67", 0.67), ("ratio>=0.75", 0.75), ("ratio>=1.0", 1.0)]
    results = {}
    for name, R in modes:
        if R is None:
            keep = lambda a, b: (idx.get(a) and idx.get(b) and
                                 idx[a]["n_intron"] >= 1 and idx[b]["n_intron"] >= 1)
        else:
            keep = lambda a, b, R=R: ratio_keep(idx.get(a), idx.get(b), R)
        kept = [(a, b) for (a, b) in base_edges if keep(a, b)]
        comp = components(kept)
        maxc = max((len(c) for c in comp), default=0)
        mdrop = sum(1 for (a, b) in mega_edges if not keep(a, b))
        pc = panel_conn(set(kept))
        realk = sum(1 for a, b in REAL_PAIRS if keep(a, b))
        domk = sum(1 for a, b in COUNTER_DOMAIN if keep(a, b))
        pstr = "/".join(pc[f] for f in ("RABL2", "APOBEC3", "RFPL"))
        print(f"  {name:14} {len(kept):6d} {len(comp):5d} {maxc:5d} {mdrop:8d} "
              f"{pstr:>20} {realk:6d}/5 {domk:5d}/4")
        results[name] = dict(edges=len(kept), families=len(comp), maxcomp=maxc,
                             mega_drop=mdrop, panel=pc, real_kept=realk, dom_kept=domk)

    # ---- the retrocopy axis: how many baseline edges are retro-edges (one side intronless)? ----
    il = {g for g in idx if idx[g]["n_intron"] == 0}
    retro_edges = [(a, b) for (a, b) in base_edges if a in il or b in il]
    retro_in_mega = [(a, b) for (a, b) in mega_edges if a in il or b in il]
    print(f"\n=== RETROCOPY AXIS (one side intronless, cDNA id high) ===")
    print(f"  baseline homology edges touching an intronless gene: {len(retro_edges)} "
          f"({100*len(retro_edges)/len(base_edges):.1f}% of all edges)")
    print(f"  of which inside the 389-mega: {len(retro_in_mega)} "
          f"({100*len(retro_in_mega)/len(mega_edges):.1f}% of mega edges) — cut by 'both_spliced'")
    # genome-wide: how many WHOLE families are pure retro/intronless artifacts?
    retro_only_comp = 0
    for c in base_comp:
        e_in = [(a, b) for (a, b) in base_edges if a in c and b in c]
        if e_in and all(a in il or b in il for (a, b) in e_in):
            retro_only_comp += 1
    print(f"  baseline families whose every edge is retro (all-intronless-incident): {retro_only_comp}")

    json.dump(dict(baseline=dict(edges=len(base_edges), families=len(base_comp),
                                 maxcomp=base_max, mega=len(mega)),
                   auc_intron_count_ratio=auc, pos=pos, neg=neg,
                   ladder=results,
                   retro_edges=len(retro_edges), retro_in_mega=len(retro_in_mega),
                   retro_only_families=retro_only_comp),
              open(os.path.join(HERE, "family_def_junction_count.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_junction_count.json')}")


if __name__ == "__main__":
    main()
