#!/usr/bin/env python3
"""family_def_feature_discovery.py — UNBIASED distinguishing-feature discovery.

Instead of hand-picking a feature (de, NM, rl) and testing it, give the data a RICH feature
set and let it pick. Principle: a genuine copy-vs-bridge distinguisher reflects two real
populations, so its distribution over the LARGE unlabeled cross-mapping set is BIMODAL,
label-free. Rank every feature by intrinsic bimodality (1-D two-cluster separation); the
handful of labels only ORIENT (which mode = copies) and VALIDATE — zero boundary-fitting,
so nothing overfits the handful. A permutation null guards against spurious bimodality
under multiple comparisons.

Key move: --eqx lets us DECOMPOSE de into mismatch-rate vs indel-rate (and their per-read
ties between the two placements), testing whether INDELS or MISMATCHES carry the signal
without asserting it.

Run: python bench/family_def_feature_discovery.py NC_073240.2 NC_073244.2 NC_073248.2
"""
import collections
import json
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import best_gene, ref_span, DELTA, DE_MAX, MIN_READS
from family_def_newbam_validate import load_denovo, build_model, recip_cov, load_models, GENOME, OLD, COV_MIN

HERE = os.path.dirname(os.path.abspath(__file__))
NEW = "/home/juanfra/winloci_scratch/GGO_mm.bam"
SAM = "/home/juanfra/miniforge3/bin/samtools"
CIG = re.compile(r"(\d+)([=XIDNMS])")


def placement_feats(f):
    """per-alignment features from an --eqx SAM line; None if unusable."""
    de = nm = as_ = rl = None
    for t in f[11:]:
        if t.startswith("de:f:"): de = float(t[5:])
        elif t.startswith("NM:i:"): nm = int(t[5:])
        elif t.startswith("AS:i:"): as_ = int(t[5:])
        elif t.startswith("rl:i:"): rl = int(t[5:])
    if de is None or de > DE_MAX:
        return None
    eq = mm = ins_b = del_b = ins_r = del_r = 0
    for n, op in CIG.findall(f[5]):
        n = int(n)
        if op == "=": eq += n
        elif op == "X": mm += n
        elif op == "I": ins_b += n; ins_r += 1
        elif op == "D": del_b += n; del_r += 1
    aln = eq + mm
    if aln < 100:
        return None
    return dict(de=de, mm_rate=mm / aln, indel_rate=(ins_b + del_b) / aln,
                indelcnt_rate=(ins_r + del_r) / aln, nm_rate=(nm or 0) / aln,
                rl_frac=(rl or 0) / aln, aln=aln, as_pb=(as_ or 0) / aln,
                refspan=int(f[3]))


def scan_feats(by_chrom, region):
    """qname -> {gene: placement_feats(best by de)}."""
    mm = collections.defaultdict(dict)
    for flag in ("-f", "0x100"), ("-F", "0x900"):
        p = subprocess.Popen([SAM, "view", flag[0], flag[1], NEW, region],
                             stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.split("\t")
            if len(f) < 11:
                continue
            pf = placement_feats(f)
            if pf is None:
                continue
            g = best_gene(by_chrom, f[2], int(f[3]) - 1, int(f[3]) - 1 + ref_span(f[5]))
            if g is None:
                continue
            d = mm[f[0]]
            if g not in d or pf["de"] < d[g]["de"]:
                d[g] = pf
        p.wait()
    return mm


def edge_features(pa, pb):
    """magnitude (mean) + tie (|diff|) features for one cross-mapping read over two placements."""
    feat = {}
    for k in ("de", "mm_rate", "indel_rate", "indelcnt_rate", "nm_rate", "rl_frac"):
        feat[f"{k}__mag"] = (pa[k] + pb[k]) / 2
        feat[f"{k}__tie"] = abs(pa[k] - pb[k])
    feat["aln_ratio"] = min(pa["aln"], pb["aln"]) / max(pa["aln"], pb["aln"])
    feat["as_ratio"] = min(pa["as_pb"], pb["as_pb"]) / max(pa["as_pb"], pb["as_pb"], 1e-9)
    return feat


def bimodality(x):
    """label-free: standardized two-cluster separation via 1-D 2-means (Lloyd), + Sarle BC."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 8 or np.ptp(x) == 0:
        return 0.0, 0.0
    c = np.percentile(x, [25, 75]).astype(float)
    for _ in range(50):
        a = np.abs(x[:, None] - c[None, :]).argmin(1)
        nc = np.array([x[a == k].mean() if (a == k).any() else c[k] for k in (0, 1)])
        if np.allclose(nc, c):
            break
        c = nc
    x0, x1 = x[a == 0], x[a == 1]
    if len(x0) < 2 or len(x1) < 2:
        return 0.0, 0.0
    sep = abs(x1.mean() - x0.mean()) / np.sqrt((x0.var() + x1.var()) / 2 + 1e-12)
    # Sarle's bimodality coefficient
    n = len(x); m3 = ((x - x.mean()) ** 3).mean(); m2 = x.var()
    g1 = m3 / (m2 ** 1.5 + 1e-12)
    k = ((x - x.mean()) ** 4).mean() / (m2 ** 2 + 1e-12) - 3
    bc = (g1 ** 2 + 1) / (k + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3) + 1e-9) + 1e-9)
    return round(float(sep), 3), round(float(bc), 3)


def auc(score, label):
    """rank AUC of score vs binary label (1 = bridge)."""
    score = np.asarray(score, float); label = np.asarray(label)
    pos, neg = score[label == 1], score[label == 0]
    if len(pos) == 0 or len(neg) == 0:
        return 0.5
    order = score.argsort()
    ranks = np.empty(len(score)); ranks[order] = np.arange(1, len(score) + 1)
    a = (ranks[label == 1].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))
    return max(a, 1 - a)   # direction-agnostic


def main():
    chroms = sys.argv[1:] or ["NC_073240.2", "NC_073244.2", "NC_073248.2"]
    models = load_models()
    import pysam
    bam_old = pysam.AlignmentFile(OLD, "rb"); genome = pysam.FastaFile(GENOME)
    rows, labels = [], []
    for region in chroms:
        by_chrom, info = load_denovo({region})
        mm = scan_feats(by_chrom, region)
        ev = collections.defaultdict(list)
        for genes in mm.values():
            if len(genes) < 2:
                continue
            it = sorted(genes.items())
            for i in range(len(it)):
                for j in range(i + 1, len(it)):
                    ga, pa = it[i]; gb, pb = it[j]
                    if abs(pa["de"] - pb["de"]) <= DELTA and max(pa["de"], pb["de"]) <= DE_MAX:
                        ev[(ga, gb)].append(edge_features(pa, pb))
        loci = {g for e in ev for g in e}
        for v in loci:
            if not models.get(v) and v in info:
                c, s, e, _ = info[v]; models[v] = build_model(bam_old, genome, c, s, e)
        for (ga, gb), recs in ev.items():
            if len(recs) < MIN_READS or not models.get(ga) or not models.get(gb):
                continue
            is_bridge = 0 if recip_cov(models[ga], models[gb]) >= COV_MIN else 1
            agg = {k: float(np.median([r[k] for r in recs])) for k in recs[0]}
            rows.append(agg); labels.append(is_bridge)
    bam_old.close(); genome.close()

    feats = list(rows[0].keys())
    labels = np.array(labels)
    n_copy, n_bridge = int((labels == 0).sum()), int((labels == 1).sum())
    print(f"=== UNBIASED feature discovery: {len(rows)} edges ({n_copy} copy / {n_bridge} bridge), "
          f"{len(feats)} features ===")

    # 1) LABEL-FREE bimodality ranking
    ranking = []
    for fk in feats:
        vals = [r[fk] for r in rows]
        sep, bc = bimodality(vals)
        a = auc(vals, labels)              # validation only (handful orients)
        ranking.append((fk, sep, bc, a))
    ranking.sort(key=lambda t: -t[1])
    print("\n  ranked by INTRINSIC BIMODALITY (label-free 2-cluster separation); AUC = handful validation")
    print(f"  {'feature':22} {'bimod_sep':>9} {'sarle_BC':>9} {'AUC_vs_~B':>10}")
    for fk, sep, bc, a in ranking:
        star = " <== bimodal+discriminative" if (sep > 1.0 and a > 0.75) else ""
        print(f"  {fk:22} {sep:>9.3f} {bc:>9.3f} {a:>10.3f}{star}")

    # 2) PERMUTATION NULL on the best feature's AUC (controls multiple comparisons / overfit)
    rng = np.random.default_rng(1729)
    best_fk = max(ranking, key=lambda t: t[3])[0]
    best_vals = np.array([r[best_fk] for r in rows])
    real_auc = auc(best_vals, labels)
    null = []
    for _ in range(2000):
        # max AUC over ALL features under shuffled labels (multiple-comparison null)
        perm = rng.permutation(labels)
        null.append(max(auc([r[fk] for r in rows], perm) for fk in feats))
    null = np.array(null)
    pval = float((null >= real_auc).mean())
    print(f"\n  PERMUTATION NULL (max-AUC-over-features, 2000 shuffles): "
          f"real best={best_fk} AUC={real_auc:.3f}; null max-AUC 95th pct={np.percentile(null,95):.3f}; "
          f"p={pval:.4f} -> {'SIGNAL beats multiple-comparison chance' if pval < 0.05 else 'NOT past chance'}")

    json.dump(dict(n_copy=n_copy, n_bridge=n_bridge,
                   ranking=[dict(feature=f, bimod_sep=s, sarle_bc=b, auc=a) for f, s, b, a in ranking],
                   best_feature=best_fk, best_auc=round(real_auc, 3),
                   perm_null_p=pval, perm_null_95=round(float(np.percentile(null, 95)), 3),
                   note="Bimodality ranking is LABEL-FREE; ~B (recip-cov>=0.30) used only to ORIENT/validate "
                        "(a proxy label; external Compara labels would make it fully unbiased). Permutation "
                        "null uses max-AUC-over-features to control the small-n multiple-comparison overfit."),
              open(os.path.join(HERE, "family_def_feature_discovery.json"), "w"), indent=2)
    print("\n[+] wrote family_def_feature_discovery.json")


if __name__ == "__main__":
    main()
