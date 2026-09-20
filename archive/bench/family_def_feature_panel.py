#!/usr/bin/env python3
"""family_def_feature_panel.py — DE-CIRCULARIZED feature validation on the Compara panel.

family_def_feature_discovery.py ranked features by LABEL-FREE bimodality on the large
unlabeled cross-mapping set, then validated against ~B (recip-cov) — but the winner
aln_ratio IS read-level recip-cov, so that AUC was semi-circular. Here we replace the
~B label with EXTERNAL Compara labels on the handful of panel families (copies that share
a backbone vs domain-sharer bridges), keeping the SAME feature extraction. The bimodality
ranking (label-free, from the large set) is unchanged; only the orientation/validation is
now external -> a fully unbiased per-feature claim.
Run: python bench/family_def_feature_panel.py
"""
import collections
import json
import os
import sys

import numpy as np
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_feature_discovery import edge_features, auc

NEW = "/home/juanfra/winloci_scratch/GGO_mm.bam"
HERE = os.path.dirname(os.path.abspath(__file__))
DE_MAX = 0.05

GENES = {
    "RABL2A": ("NC_073235.2", 15131653, 15147533), "RABL2B": ("NC_086018.1", 48818440, 48832011),
    "APOBEC3D": ("NC_086018.1", 37023493, 37035944), "APOBEC3F": ("NC_086018.1", 37043278, 37058621),
    "RFPL1": ("NC_086018.1", 27445867, 27472024), "RFPL2": ("NC_086018.1", 30208643, 30221673),
    "RFPL3": ("NC_086018.1", 30368560, 30376055),
    "AK6": ("NC_073243.2", 40875017, 40892984), "LOC115934278": ("NC_073227.2", 47415194, 47415981),
    "CCDC196": ("NC_073239.2", 83376335, 83388628), "LOC129526440": ("NC_073238.2", 33024134, 33036450),
    "ASDURF": ("NC_073235.2", 92106380, 92111318), "ASNSD1": ("NC_073235.2", 92111059, 92115776),
    "CASP8": ("NC_073235.2", 103850223, 103901051), "FLACC1": ("NC_073235.2", 103898904, 103977183),
    "CREB1": ("NC_073235.2", 110133100, 110220981), "METTL21A": ("NC_073235.2", 110198842, 110241136),
    "GCA": ("NC_073235.2", 64785386, 64835227), "KCNH7": ("NC_073235.2", 64832333, 65309180),
    "GPR39": ("NC_073235.2", 34964263, 35195860), "LYPD1": ("NC_073235.2", 35194087, 35224540),
}
# external Compara labels: 0 = real copy (shares backbone), 1 = domain-sharer bridge
PAIRS = [
    ("RABL2", "RABL2A", "RABL2B", 0), ("AK6", "AK6", "LOC115934278", 0),
    ("CCDC196", "CCDC196", "LOC129526440", 0), ("APOBEC3", "APOBEC3D", "APOBEC3F", 0),
    ("RFPL1_2", "RFPL1", "RFPL2", 0), ("RFPL1_3", "RFPL1", "RFPL3", 0), ("RFPL2_3", "RFPL2", "RFPL3", 0),
    ("ASDURF_ds", "ASDURF", "ASNSD1", 1), ("CASP8_ds", "CASP8", "FLACC1", 1),
    ("CREB1_ds", "CREB1", "METTL21A", 1), ("GCA_ds", "GCA", "KCNH7", 1), ("GPR39_ds", "GPR39", "LYPD1", 1),
]


def feats_for_gene(bam, gene):
    """qname -> placement features for reads overlapping the gene (best by de)."""
    c, s, e = GENES[gene]
    out = {}
    for r in bam.fetch(c, max(0, s - 2000), e + 2000):
        if r.is_unmapped or r.is_supplementary:
            continue
        if r.reference_end is None or r.reference_end < s or r.reference_start > e:
            continue
        de = r.get_tag("de") if r.has_tag("de") else None
        if de is None or de > DE_MAX:
            continue
        eq = mm = ins_b = dl_b = ins_r = dl_r = 0
        for op, n in (r.cigartuples or []):
            if op == 7: eq += n
            elif op == 8: mm += n
            elif op == 1: ins_b += n; ins_r += 1
            elif op == 2: dl_b += n; dl_r += 1
        aln = eq + mm
        if aln < 100:
            continue
        nm = r.get_tag("NM") if r.has_tag("NM") else 0
        as_ = r.get_tag("AS") if r.has_tag("AS") else 0
        rl = r.get_tag("rl") if r.has_tag("rl") else 0
        pf = dict(de=de, mm_rate=mm / aln, indel_rate=(ins_b + dl_b) / aln,
                  indelcnt_rate=(ins_r + dl_r) / aln, nm_rate=nm / aln,
                  rl_frac=rl / aln, aln=aln, as_pb=as_ / aln)
        if r.query_name not in out or de < out[r.query_name]["de"]:
            out[r.query_name] = pf
    return out


def main():
    bam = pysam.AlignmentFile(NEW, "rb")
    rows, labels, names = [], [], []
    for name, ga, gb, lab in PAIRS:
        fa, fb = feats_for_gene(bam, ga), feats_for_gene(bam, gb)
        shared = set(fa) & set(fb)
        if len(shared) < 3:
            print(f"  {name:12} ({'bridge' if lab else 'copy '}): {len(shared)} cross-mapping reads — skip")
            continue
        ef = [edge_features(fa[q], fb[q]) for q in shared]
        agg = {k: float(np.median([e[k] for e in ef])) for k in ef[0]}
        rows.append(agg); labels.append(lab); names.append(name)
        print(f"  {name:12} ({'bridge' if lab else 'copy '}): {len(shared)} reads  "
              f"aln_ratio={agg['aln_ratio']:.3f} nm_rate__tie={agg['nm_rate__tie']:.4f}")
    bam.close()

    labels = np.array(labels)
    n_copy, n_bridge = int((labels == 0).sum()), int((labels == 1).sum())
    feats = list(rows[0].keys())
    print(f"\n=== EXTERNAL (Compara) validation: {len(rows)} panel pairs "
          f"({n_copy} copy / {n_bridge} bridge) ===")

    # bimodality ranking from the large unlabeled set (unchanged, label-free)
    disc = json.load(open(os.path.join(HERE, "family_def_feature_discovery.json")))
    bimod = {r["feature"]: r["bimod_sep"] for r in disc["ranking"]}

    res = []
    for fk in feats:
        a = auc([r[fk] for r in rows], labels)   # AUC vs EXTERNAL Compara label
        res.append((fk, bimod.get(fk, 0.0), a))
    res.sort(key=lambda t: -t[1])   # sort by intrinsic bimodality (the discovery axis)
    print(f"\n  {'feature':22} {'bimod_sep(large)':>16} {'AUC_vs_COMPARA':>15}")
    for fk, b, a in res:
        star = " <==" if (b > 1.0 and a >= 0.80) else ""
        print(f"  {fk:22} {b:>16.3f} {a:>15.3f}{star}")

    # permutation null on the panel (max-AUC over features under shuffled external labels)
    rng = np.random.default_rng(1729)
    real_max = max(a for _, _, a in res)
    null = []
    for _ in range(5000):
        perm = rng.permutation(labels)
        null.append(max(auc([r[fk] for r in rows], perm) for fk in feats))
    null = np.array(null)
    p = float((null >= real_max).mean())
    top_bimodal = res[0]
    print(f"\n  top-by-bimodality feature: {top_bimodal[0]} -> Compara AUC {top_bimodal[2]:.3f}")
    print(f"  PERMUTATION NULL (external labels, 5000 shuffles): real max AUC={real_max:.3f}; "
          f"null 95th={np.percentile(null,95):.3f}; p={p:.4f} -> "
          f"{'beats chance' if p < 0.05 else 'NOT past chance (handful too small)'}")
    json.dump(dict(n_copy=n_copy, n_bridge=n_bridge, pairs=names,
                   feature_auc=[dict(feature=f, bimod_sep=b, compara_auc=a) for f, b, a in res],
                   perm_p=p, note="External Compara labels (not ~B). Bimodality from the large unlabeled "
                   "set; AUC on the panel handful. n is tiny — read directionally, with the permutation null."),
              open(os.path.join(HERE, "family_def_feature_panel.json"), "w"), indent=2)
    print("\n[+] wrote family_def_feature_panel.json")


if __name__ == "__main__":
    main()
