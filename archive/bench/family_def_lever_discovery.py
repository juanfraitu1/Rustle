#!/usr/bin/env python3
"""family_def_lever_discovery.py — what closes the definition FURTHER (precision residual)?

~B (reciprocal coverage >= 0.30) prunes most bridges, but some non-paralog edges survive it
(DNA-bridges that pass ~B = FALSE SURVIVORS). This searches for a SECOND-order lever: among
the ~B-survivors, which feature best separates true DNA-paralogs from the false survivors?
Adds copy-level / distributional features (cross-mapping SPATIAL SPREAD, gene-length symmetry,
homology-block shape) to the read-level set, ranks by AUC vs the EXTERNAL DNA label among
survivors, with a permutation null. The winner (if any) is the lever to add after ~B.
Run: python bench/family_def_lever_discovery.py NC_073240.2 NC_073244.2 NC_073248.2
"""
import collections
import json
import os
import re
import subprocess
import sys

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import GENES_BED, ref_span, best_gene, DELTA, DE_MAX, MIN_READS
from family_def_feature_discovery import scan_feats, edge_features, auc
from family_def_newbam_validate import build_model, GENOME, NEW, TMP, COV_MIN
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


def recip_and_core(a, b):
    """one minimap2 alignment -> (reciprocal coverage, contiguous-core fraction)."""
    os.makedirs(TMP, exist_ok=True)
    if not a or not b:
        return 0.0, 0.0
    open(f"{TMP}/a.fa", "w").write(f">a\n{a}\n"); open(f"{TMP}/b.fa", "w").write(f">b\n{b}\n")
    out = subprocess.run(f"minimap2 -c --eqx -x asm20 {TMP}/a.fa {TMP}/b.fa 2>/dev/null",
                         shell=True, capture_output=True, text=True).stdout
    recip = core = 0.0
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        ql, qs, qe, tl, ts, te = int(f[1]), int(f[2]), int(f[3]), int(f[6]), int(f[7]), int(f[8])
        rc = min((qe - qs) / max(ql, 1), (te - ts) / max(tl, 1))
        cg = next((t[5:] for t in f[12:] if t.startswith("cg:Z:")), "")
        run = longest = 0
        for n, op in re.findall(r"(\d+)([=XID])", cg):
            n = int(n)
            if op in "=X":
                run += n; longest = max(longest, run)
            else:
                run = 0
        if rc > recip:
            recip, core = rc, longest / max(min(ql, tl), 1)
    return recip, core


def main():
    chroms = sys.argv[1:] or ["NC_073240.2", "NC_073244.2", "NC_073248.2"]
    by_chrom, coord = load_genes(set(chroms))
    Hd, _ = dna_homology()
    bam = pysam.AlignmentFile(NEW, "rb"); genome = pysam.FastaFile(GENOME)
    model_cache = {}

    def model(g):
        if g not in model_cache:
            c, s, e = coord[g]
            model_cache[g] = build_model(bam, genome, c, s, e)
        return model_cache[g]

    def dna_label(ga, gb):
        r = Hd.get((ga, gb) if ga < gb else (gb, ga))
        if r is None or r.get("id", 0) == 0:
            return 1
        if r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
            return 0
        return None

    rows = []
    for region in chroms:
        mm = scan_feats(by_chrom, region)
        ev = collections.defaultdict(list)
        for genes in mm.values():
            if len(genes) < 2:
                continue
            it = sorted(genes.items())
            for i in range(len(it)):
                for j in range(i + 1, len(it)):
                    ga, pa = it[i]; gb, pb = it[j]
                    if ga != gb and abs(pa["de"] - pb["de"]) <= DELTA and max(pa["de"], pb["de"]) <= DE_MAX:
                        ev[(ga, gb)].append((pa, pb))
        for (ga, gb), recs in ev.items():
            if len(recs) < MIN_READS:
                continue
            lab = dna_label(ga, gb)
            if lab is None or ga not in coord or gb not in coord:
                continue
            feat = {k: float(np.median([edge_features(pa, pb)[k] for (pa, pb) in recs]))
                    for k in edge_features(*recs[0])}
            # NEW distributional / copy-level features
            posA = [pa["refspan"] for (pa, pb) in recs]; posB = [pb["refspan"] for (pa, pb) in recs]
            lenA = coord[ga][2] - coord[ga][1]; lenB = coord[gb][2] - coord[gb][1]
            feat["xmap_spread_min"] = min((max(posA) - min(posA)) / max(lenA, 1),
                                          (max(posB) - min(posB)) / max(lenB, 1))
            feat["len_ratio"] = min(lenA, lenB) / max(lenA, lenB)
            feat["nreads"] = len(recs)
            recip, core = recip_and_core(model(ga), model(gb))
            feat["recip_cov"] = recip          # ~B itself
            feat["core_frac"] = core           # copy-level contiguous homology block
            rows.append((feat, lab, recip, ga, gb, region))
    bam.close(); genome.close()

    labels = np.array([l for _, l, *_ in rows])
    surv = [(f, l) for (f, l, rc, *_) in rows if rc >= COV_MIN]    # ~B survivors
    surv_full = [(f, l, ga, gb, reg) for (f, l, rc, ga, gb, reg) in rows if rc >= COV_MIN]
    sl = np.array([l for _, l in surv])
    n_all, n_par, n_brg = len(rows), int((labels == 0).sum()), int((labels == 1).sum())
    print(f"=== precision residual: {n_all} genuine-crossmap gene edges, DNA {n_par} paralog / {n_brg} bridge ===")
    print(f"  ~B (recip-cov>=0.30): {sum(1 for _,_,rc,*_ in rows if rc>=COV_MIN)} survive; "
          f"of survivors {int((sl==0).sum())} DNA-paralog (TRUE), {int((sl==1).sum())} DNA-bridge (FALSE SURVIVOR)")
    if int((sl == 1).sum()) == 0:
        print("  -> ~B already removes every DNA-bridge here; no residual lever needed on this set.")
    feats = [k for k in surv[0][0] if k != "nreads"]
    print(f"\n  among the {len(surv)} ~B-survivors, feature separation of TRUE vs FALSE-SURVIVOR (AUC vs DNA):")
    res = []
    for fk in feats:
        a = auc([f[fk] for f, _ in surv], sl)
        res.append((fk, a))
    res.sort(key=lambda t: -t[1])
    for fk, a in res:
        print(f"    {fk:20} AUC={a:.3f}{'  <== candidate lever' if a >= 0.75 else ''}")

    rng = np.random.default_rng(1729)
    real_max = max(a for _, a in res)
    null = [max(auc([f[fk] for f, _ in surv], rng.permutation(sl)) for fk in feats) for _ in range(2000)]
    p = float((np.array(null) >= real_max).mean())
    print(f"\n  best survivor-lever: {res[0][0]} AUC={res[0][1]:.3f}; permutation p={p:.4f} "
          f"-> {'real lever' if p < 0.05 else 'NOT past chance (residual too small / no lever)'}")
    # --- adversarial verification of the top lever ---
    best_fk = res[0][0]
    print(f"\n  === VERIFY lever '{best_fk}' ===")
    # 1) the 9 false survivors: identity + lever value vs the true-paralog distribution
    fs = [(ga, gb, reg, f[best_fk]) for (f, l, ga, gb, reg) in surv_full if l == 1]
    tp_vals = sorted(f[best_fk] for (f, l, *_ ) in surv_full if l == 0)
    print(f"  TRUE paralogs {best_fk}: median={np.median(tp_vals):.4f} p90={np.percentile(tp_vals,90):.4f}")
    print(f"  FALSE survivors (n={len(fs)}):")
    for ga, gb, reg, v in sorted(fs, key=lambda x: -x[3]):
        print(f"    {ga}~{gb} [{reg}] {best_fk}={v:.4f}")
    # 2) per-chrom: is the separation driven by one chromosome?
    print("  per-chrom false-survivor count:",
          dict(collections.Counter(reg for _, _, reg, _ in fs)))
    # 3) operating point: a cut on the lever that removes the false survivors, at what TP cost?
    import numpy as _np
    fvals = _np.array([v for *_, v in fs]); tvals = _np.array(tp_vals)
    hi = fvals.mean() > tvals.mean()   # do false survivors sit HIGH on this lever?
    print("  operating points (cut to drop false survivors; cost = true paralogs also dropped):")
    for q in [50, 75, 90, 95]:
        thr = _np.percentile(fvals, 100 - q) if hi else _np.percentile(fvals, q)
        dropped_false = int((fvals >= thr).sum()) if hi else int((fvals <= thr).sum())
        cost_true = int((tvals >= thr).sum()) if hi else int((tvals <= thr).sum())
        print(f"    thr={thr:.4f}: drops {dropped_false}/{len(fs)} false survivors, "
              f"costs {cost_true}/{len(tvals)} true paralogs")

    json.dump(dict(n_edges=n_all, n_paralog=n_par, n_bridge=n_brg,
                   survivors=len(surv), false_survivors=int((sl == 1).sum()),
                   survivor_lever_auc=[dict(feature=f, auc=a) for f, a in res], perm_p=p,
                   best_lever=best_fk,
                   false_survivor_detail=[dict(ga=ga, gb=gb, chrom=reg, val=round(v, 4)) for ga, gb, reg, v in fs],
                   per_chrom_false=dict(collections.Counter(reg for _, _, reg, _ in fs))),
              open(os.path.join(HERE, "family_def_lever_discovery.json"), "w"), indent=2)
    print("\n[+] wrote family_def_lever_discovery.json")


if __name__ == "__main__":
    main()
