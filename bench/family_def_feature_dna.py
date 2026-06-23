#!/usr/bin/env python3
"""family_def_feature_dna.py — DE-CIRCULARIZED feature validation with EXTERNAL DNA labels.

The discovery validated against ~B (recip-cov), and the winner aln_ratio IS read-level
recip-cov => semi-circular. Here we keep the LABEL-FREE bimodality ranking but replace the
label with an EXTERNAL, DNA-based one: scan GENUINE cross-mapping gene-gene edges (per-record
best-overlap attribution, so no nesting artifact), and label each edge by all-vs-all cDNA
DNA homology (paralog = id>=0.90 AND max reciprocal cDNA coverage >=0.30) vs bridge (reads
cross-map but the genes are NOT DNA paralogs). DNA homology is independent of RNA reciprocal
coverage, so a feature that separates the DNA classes is validated non-circularly.
Run: python bench/family_def_feature_dna.py NC_073240.2 NC_073244.2 NC_073248.2
"""
import collections
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import GENES_BED, DELTA, DE_MAX, MIN_READS
from family_def_feature_discovery import scan_feats, edge_features, bimodality, auc
from family_def_read_filters import dna_homology

HERE = os.path.dirname(os.path.abspath(__file__))


def load_genes(chroms):
    by = collections.defaultdict(list)
    with open(GENES_BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4 or p[0] not in chroms:
                continue
            by[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c in by:
        by[c].sort()
    return by


def main():
    chroms = sys.argv[1:] or ["NC_073240.2", "NC_073244.2", "NC_073248.2"]
    by_chrom = load_genes(set(chroms))
    Hd, _ = dna_homology()

    def dna_label(ga, gb):
        r = Hd.get((ga, gb) if ga < gb else (gb, ga))
        if r is None or r.get("id", 0) == 0:
            return 1                      # no DNA homology -> bridge
        if r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30:
            return 0                      # DNA paralog -> real copy
        return None                       # sub-bar / ambiguous -> exclude

    rows, labels = [], []
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
                    if ga == gb:
                        continue
                    if abs(pa["de"] - pb["de"]) <= DELTA and max(pa["de"], pb["de"]) <= DE_MAX:
                        ev[(ga, gb)].append(edge_features(pa, pb))
        for (ga, gb), recs in ev.items():
            if len(recs) < MIN_READS:
                continue
            lab = dna_label(ga, gb)
            if lab is None:
                continue
            agg = {k: float(np.median([r[k] for r in recs])) for k in recs[0]}
            rows.append(agg); labels.append(lab)

    labels = np.array(labels)
    n_copy, n_bridge = int((labels == 0).sum()), int((labels == 1).sum())
    feats = list(rows[0].keys())
    print(f"=== DE-CIRCULARIZED: {len(rows)} gene-gene crossmap edges with DNA labels "
          f"({n_copy} DNA-paralog / {n_bridge} non-paralog bridge), {len(feats)} features ===")

    res = []
    for fk in feats:
        vals = [r[fk] for r in rows]
        sep, bc = bimodality(vals)
        a = auc(vals, labels)
        res.append((fk, sep, a))
    res.sort(key=lambda t: -t[1])
    print(f"\n  ranked by LABEL-FREE bimodality; AUC vs EXTERNAL DNA label")
    print(f"  {'feature':22} {'bimod_sep':>9} {'AUC_vs_DNA':>11}")
    for fk, sep, a in res:
        star = " <== bimodal+discriminative(DNA)" if (sep > 1.0 and a >= 0.75) else ""
        print(f"  {fk:22} {sep:>9.3f} {a:>11.3f}{star}")

    rng = np.random.default_rng(1729)
    real_max = max(a for _, _, a in res)
    null = [max(auc([r[fk] for r in rows], rng.permutation(labels)) for fk in feats) for _ in range(2000)]
    null = np.array(null); p = float((null >= real_max).mean())
    top = res[0]
    print(f"\n  top-by-bimodality: {top[0]} -> DNA AUC {top[2]:.3f}")
    print(f"  PERMUTATION NULL (DNA labels, 2000 shuffles): real max AUC={real_max:.3f}; "
          f"null 95th={np.percentile(null,95):.3f}; p={p:.4f} -> "
          f"{'beats chance' if p < 0.05 else 'NOT past chance'}")
    json.dump(dict(n_copy=n_copy, n_bridge=n_bridge,
                   ranking=[dict(feature=f, bimod_sep=s, dna_auc=a) for f, s, a in res],
                   perm_p=p, note="EXTERNAL DNA-paralogy label (cDNA all-vs-all, NOT recip-cov). Genuine "
                   "per-record crossmap edges. De-circularizes the ~B validation."),
              open(os.path.join(HERE, "family_def_feature_dna.json"), "w"), indent=2)
    print("\n[+] wrote family_def_feature_dna.json")


if __name__ == "__main__":
    main()
