#!/usr/bin/env python3
"""family_def_joint.py — the HIGHER-DIMENSIONAL test: do the orthogonal levers COMBINE?

Each single lever is soft as a hard cut, but they fail on DIFFERENT subsets (repeat-rich,
diverged, size-2) => orthogonal. Combine all axes -- sequence (recip_cov, aln_ratio),
coverage/FLNC (cov_frac, qaln), topology (jaccard_nbr, adamic_adar), repeat (rl_frac),
divergence (de, nm_tie) -- in an L2-logistic, evaluated by LEAVE-ONE-FAMILY-OUT grouped CV
(group = connected component of the ~R graph, so correlated within-family edges never split
across train/test = no leakage). Does the joint CV-AUC beat the best single feature, and is
there a clean operating point (bridges out, real families in)?
Run: python bench/family_def_joint.py NC_073240.2 NC_073244.2 NC_073248.2 [more chroms]
"""
import collections
import json
import math
import os
import re
import sys

import numpy as np
import networkx as nx
import pysam
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_def_genomewide import GENES_BED, best_gene, ref_span, pair_evidence, components, DELTA, DE_MAX, MIN_READS
from family_def_newbam_validate import build_model, recip_cov, GENOME, NEW, COV_MIN
from family_def_feature_discovery import auc
from family_def_read_filters import dna_homology

HERE = os.path.dirname(os.path.abspath(__file__))
SAM = "/home/juanfra/miniforge3/bin/samtools"
CIG = re.compile(r"(\d+)([=XIDNMSH])")


def placement(f):
    de = next((float(t[5:]) for t in f[11:] if t.startswith("de:f:")), None)
    if de is None or de > DE_MAX:
        return None
    nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), 0)
    rl = next((int(t[5:]) for t in f[11:] if t.startswith("rl:i:")), 0)
    eq = x = qa = qt = 0
    pos = int(f[3]); cur = pos; bs = pos; blocks = []
    for n, op in CIG.findall(f[5]):
        n = int(n)
        if op == "=":
            eq += n; qa += n; qt += n; cur += n
        elif op == "X":
            x += n; qa += n; qt += n; cur += n
        elif op == "I":
            qa += n; qt += n
        elif op == "S":
            qt += n
        elif op in "DN":
            if cur > bs:
                blocks.append((bs, cur))
            cur += n; bs = cur
    if cur > bs:
        blocks.append((bs, cur))
    aln = eq + x
    if aln < 100 or qt == 0:
        return None
    return dict(de=de, nm_rate=nm / aln, rl=rl, aln=aln, qaln=qa / qt, blocks=blocks)


def scan(by, chroms):
    import subprocess
    mm = collections.defaultdict(dict)
    for c in chroms:
        for flag in ("-f", "0x100"), ("-F", "0x900"):
            p = subprocess.Popen([SAM, "view", flag[0], flag[1], NEW, c],
                                 stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
            for line in p.stdout:
                f = line.split("\t")
                if len(f) < 11:
                    continue
                pf = placement(f)
                if pf is None:
                    continue
                g = best_gene(by, f[2], int(f[3]) - 1, int(f[3]) - 1 + ref_span(f[5]))
                if g is None:
                    continue
                d = mm[f[0]]
                if g not in d or pf["de"] < d[g]["de"]:
                    d[g] = pf
            p.wait()
    return mm


def covered(blocks):
    iv = sorted(b for bl in blocks for b in bl)
    if not iv:
        return 0
    tot, cs, ce = 0, iv[0][0], iv[0][1]
    for s, e in iv[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            tot += ce - cs; cs, ce = s, e
    return tot + ce - cs


def logistic_cv(X, y, groups, l2=1.0, seed=0):
    """L2-logistic, grouped 5-fold CV (groups never split); returns pooled out-of-fold probs."""
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    X = np.hstack([X, np.ones((len(X), 1))])
    gids = list(dict.fromkeys(groups))
    np.random.default_rng(seed).shuffle(gids)
    folds = {g: i % 5 for i, g in enumerate(gids)}      # 5 grouped folds (shuffled)
    fold = np.array([folds[g] for g in groups])
    oof = np.zeros(len(y))
    for k in range(5):
        tr, te = fold != k, fold == k
        if te.sum() == 0 or y[tr].sum() in (0, tr.sum()):
            oof[te] = y[tr].mean() if tr.sum() else 0.5
            continue
        def loss(w):
            z = X[tr] @ w
            ll = np.sum(np.logaddexp(0, z) - y[tr] * z)
            return ll + l2 * np.sum(w[:-1] ** 2)
        w = minimize(loss, np.zeros(X.shape[1]), method="L-BFGS-B").x
        oof[te] = 1 / (1 + np.exp(-(X[te] @ w)))
    return oof


def main():
    chroms = sys.argv[1:] or ["NC_073240.2", "NC_073244.2", "NC_073248.2"]
    by, coord = collections.defaultdict(list), {}
    with open(GENES_BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4 and p[0] in set(chroms):
                by[p[0]].append((int(p[1]), int(p[2]), p[3])); coord[p[3]] = (p[0], int(p[1]), int(p[2]))
    for c in by:
        by[c].sort()
    mm = scan(by, set(chroms))
    ev = pair_evidence({q: {g: pf["de"] for g, pf in gd.items()} for q, gd in mm.items()})
    edges, _ = components(ev, DELTA, DE_MAX, MIN_READS)
    G = nx.Graph()
    for ga, gb, n in edges:
        G.add_edge(ga, gb)
    deg = dict(G.degree())
    comp_id = {}
    for i, cc in enumerate(nx.connected_components(G)):
        for g in cc:
            comp_id[g] = i
    Hd, _ = dna_homology()
    bam = pysam.AlignmentFile(NEW, "rb"); genome = pysam.FastaFile(GENOME)
    mcache = {}

    def model(g):
        if g not in mcache and g in coord:
            mcache[g] = build_model(bam, genome, *coord[g])
        return mcache.get(g, "")

    def dna(ga, gb):
        r = Hd.get((ga, gb) if ga < gb else (gb, ga))
        if r is None or r.get("id", 0) == 0:
            return 1
        return 0 if (r["id"] >= 0.90 and max(r["cov_a"], r["cov_b"]) >= 0.30) else None

    feat_names = ["recip_cov", "aln_ratio", "cov_frac_min", "qaln_min",
                  "jaccard_nbr", "adamic_adar", "rl_frac", "nm_tie", "de_mag"]
    rows, labels, groups = [], [], []
    for ga, gb, n in edges:
        lab = dna(ga, gb)
        if lab is None or ga not in coord or gb not in coord:
            continue
        recs = [(mm[q][ga], mm[q][gb]) for q in mm if ga in mm[q] and gb in mm[q]
                and abs(mm[q][ga]["de"] - mm[q][gb]["de"]) <= DELTA and max(mm[q][ga]["de"], mm[q][gb]["de"]) <= DE_MAX]
        if len(recs) < MIN_READS:
            continue
        eu_i = covered([pa["blocks"] for pa, pb in recs]); eu_j = covered([pb["blocks"] for pa, pb in recs])
        Na, Nb = set(G[ga]) - {gb}, set(G[gb]) - {ga}
        common = Na & Nb
        feat = [
            recip_cov(model(ga), model(gb)),
            float(np.median([min(pa["aln"], pb["aln"]) / max(pa["aln"], pb["aln"]) for pa, pb in recs])),
            min(covered([pa["blocks"] for pa, pb in recs]) / max(eu_i, 1),
                covered([pb["blocks"] for pa, pb in recs]) / max(eu_j, 1)),
            float(np.median([min(pa["qaln"], pb["qaln"]) for pa, pb in recs])),
            len(common) / max(len(Na | Nb), 1),
            sum(1.0 / math.log(deg[w]) for w in common if deg[w] > 1),
            float(np.median([(pa["rl"] / pa["aln"] + pb["rl"] / pb["aln"]) / 2 for pa, pb in recs])),
            float(np.median([abs(pa["nm_rate"] - pb["nm_rate"]) for pa, pb in recs])),
            float(np.median([(pa["de"] + pb["de"]) / 2 for pa, pb in recs])),
        ]
        rows.append(feat); labels.append(lab); groups.append(comp_id[ga])
    bam.close(); genome.close()

    X = np.array(rows); y = np.array(labels)
    n_par, n_brg = int((y == 0).sum()), int((y == 1).sum())
    print(f"=== JOINT multi-dim test: {len(rows)} ~R edges ({n_par} paralog / {n_brg} bridge), "
          f"{X.shape[1]} features, {len(set(groups))} family-groups ===")
    print(f"\n  single-feature AUC (for reference):")
    singles = []
    for i, fn in enumerate(feat_names):
        a = auc(X[:, i], y); singles.append(a)
        print(f"    {fn:14} {a:.3f}")
    best_single = max(singles)
    best_name = feat_names[int(np.argmax(singles))]
    # in-sample joint (overfitting check) vs CV at several regularizations
    Xs = (X - X.mean(0)) / (X.std(0) + 1e-9); Xa = np.hstack([Xs, np.ones((len(Xs), 1))])
    def fit(Xtr, ytr, l2):
        return minimize(lambda w: np.sum(np.logaddexp(0, Xtr @ w) - ytr * (Xtr @ w)) + l2 * np.sum(w[:-1] ** 2),
                        np.zeros(Xtr.shape[1]), method="L-BFGS-B").x
    insample = auc(1 / (1 + np.exp(-(Xa @ fit(Xa, y, 1.0)))), y)
    print(f"\n  best single feature: {best_name} AUC={best_single:.3f}")
    print(f"  JOINT logistic IN-SAMPLE AUC = {insample:.3f}  (vs best single {best_single:.3f})")
    print(f"  JOINT REPEATED grouped CV (10 random fold splits, mean +/- std):")
    best_cv = (0, None)
    for l2 in [1.0, 3.0, 10.0]:
        aucs = [auc(logistic_cv(X, y, groups, l2, seed=s), y) for s in range(10)]
        m, sd = np.mean(aucs), np.std(aucs)
        print(f"    L2={l2:>5}: CV-AUC = {m:.3f} +/- {sd:.3f}")
        if m > best_cv[0]:
            best_cv = (m, l2)
    # pooled OOF at the best L2 (averaged over splits) for the operating point
    oofs = np.array([logistic_cv(X, y, groups, best_cv[1], seed=s) for s in range(10)])
    oof = oofs.mean(0); cv_auc = best_cv[0]
    # operating point on the CV scores
    par_s = sorted(oof[y == 0]); brg_s = sorted(oof[y == 1])
    print(f"  CV score: paralog median={np.median(par_s):.2f}, bridge median={np.median(brg_s):.2f}")
    for q in [0.5, 0.7, 0.9]:
        thr = np.quantile(brg_s, 1 - q)
        bc = int((np.array(brg_s) >= thr).sum()); pc = int((np.array(par_s) >= thr).sum())
        print(f"    score>={thr:.2f}: flags {bc}/{n_brg} bridges ({100*bc/n_brg:.0f}%) at {pc}/{n_par} paralog cost ({100*pc/n_par:.0f}%)")
    json.dump(dict(n_par=n_par, n_brg=n_brg, n_features=X.shape[1], n_groups=len(set(groups)),
                   single_auc=dict(zip(feat_names, [round(a, 3) for a in singles])),
                   joint_cv_auc=round(cv_auc, 3), best_single=round(best_single, 3)),
              open(os.path.join(HERE, "family_def_joint.json"), "w"), indent=2)
    print("\n[+] wrote family_def_joint.json")


if __name__ == "__main__":
    main()
