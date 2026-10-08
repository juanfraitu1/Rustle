#!/usr/bin/env python3
"""Does Soto's copy-number cut follow duplicon boundaries? (docs/archive/2026-09/PREREG_soto_cn_vs_duplicon_2026-09-30.md, KEY=cnduplicon)

    python3 bench/soto_m2/soto_m2_duplicons.py --geneset elig.tsv --full-geneset full.tsv \
        --exons sd98_gene_exons.tsv --duplicons chm13.draft_v1.0_plus38Y_dupmasker_colors.bed

Clusters = sequence families (exon map-back edges, no copy-number gate) holding >= 2 Soto families with >= 2 clean members.
Per gene: exonic bases covered by each duplicon (overlapping duplicons all count). Similarity = weighted Jaccard. Per-cluster
statistic = mean within-Soto-family similarity - mean between-family similarity; pooled = mean over clusters; null = Soto
labels shuffled within each cluster (sizes kept), 10,000 shuffles, seed 20260930. Decision rule and halves: the prereg, §4.
"""
import argparse
import bisect
import csv
import os
import random
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "soto"))
import soto_replication as sr  # noqa: E402

SOTO = os.path.join(HERE, "..", "soto")
N_SHUFFLE, SEED = 10000, 20260930


def load_duplicons(path):
    col = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            col[f[0]].append((int(f[1]), int(f[2]), f[3]))
    for c in col:
        col[c].sort()
    maxlen = {c: max(b - a for a, b, _ in v) for c, v in col.items()}
    return col, maxlen


def compositions(exon_path, col, maxlen, genes):
    ex = defaultdict(list)
    with open(exon_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene_id"] in genes:
                ex[r["gene_id"]].append((r["chrom"], int(r["start"]), int(r["end"])))
    comp = {}
    for g in genes:
        hit = defaultdict(int)
        for c, s, e in ex.get(g, []):
            v = col.get(c, [])
            i = bisect.bisect_left(v, (s - maxlen.get(c, 0),))
            for a, b, d in v[i:]:
                if a >= e:
                    break
                o = min(b, e) - max(a, s)
                if o > 0:
                    hit[d] += o
        comp[g] = dict(hit)
    return comp


def wjac(x, y):
    if not x or not y:
        return 0.0
    keys = set(x) | set(y)
    mx = sum(max(x.get(k, 0), y.get(k, 0)) for k in keys)
    return sum(min(x.get(k, 0), y.get(k, 0)) for k in keys) / mx if mx else 0.0


def set_jac(x, y):
    if not x or not y:
        return 0.0
    return len(x & y) / len(x | y)


def delta(sim, labels):
    """mean within-family similarity - mean between-family similarity; None if either side has no pair."""
    n = len(labels)
    w = b = 0.0
    nw = nb = 0
    for i in range(n):
        for j in range(i + 1, n):
            if labels[i] == labels[j]:
                w += sim[i][j]
                nw += 1
            else:
                b += sim[i][j]
                nb += 1
    if not nw or not nb:
        return None
    return w / nw - b / nb


def run_test(clusters, simf, rng):
    """clusters: list of (genes, labels). Returns per-cluster (delta, p), pooled delta, pooled p."""
    mats = []
    for genes, labels in clusters:
        sim = [[simf(a, b) for b in genes] for a in genes]
        mats.append((sim, labels))
    obs = [delta(s, l) for s, l in mats]
    keep = [k for k, d in enumerate(obs) if d is not None]
    pooled = sum(obs[k] for k in keep) / len(keep)
    ge_pooled = 0
    ge_each = [0] * len(mats)
    for _ in range(N_SHUFFLE):
        tot = 0.0
        for k in keep:
            sim, labels = mats[k]
            perm = labels[:]
            rng.shuffle(perm)
            d = delta(sim, perm)
            tot += d
            if d >= obs[k] - 1e-12:
                ge_each[k] += 1
        if tot / len(keep) >= pooled - 1e-12:
            ge_pooled += 1
    per = [(obs[k], (1 + ge_each[k]) / (1 + N_SHUFFLE)) if obs[k] is not None else (None, None) for k in range(len(mats))]
    return per, pooled, (1 + ge_pooled) / (1 + N_SHUFFLE), len(keep)


def verdict(per, pooled_p):
    ds = [d for d, _ in per if d is not None]
    frac_pos = sum(d > 0 for d in ds) / len(ds) if ds else 0.0
    if pooled_p < 0.01 and frac_pos >= 2 / 3:
        return "FOLLOWS", frac_pos
    if pooled_p < 0.05:
        return "PARTIAL", frac_pos
    return "DOES NOT FOLLOW", frac_pos


def spearman(x, y):
    def ranks(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            for k in range(i, j + 1):
                r[o[k]] = (i + j) / 2
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    return sxy / (sxx * syy) ** 0.5 if sxx and syy else float("nan")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", default=os.path.join(SOTO, "soto_famCN_S1C.tsv"))
    ap.add_argument("--shared", default=os.path.join(SOTO, "shared_exons_5154_exon_mapback.tsv"))
    ap.add_argument("--split", default=os.path.join(SOTO, "soto_split_2026-09-29.tsv"))
    ap.add_argument("--geneset", required=True)
    ap.add_argument("--full-geneset", required=True)
    ap.add_argument("--exons", required=True)
    ap.add_argument("--duplicons", required=True)
    ap.add_argument("--out-tsv", help="per-cluster rows (descriptive output only; the test is unchanged)")
    a = ap.parse_args(argv)

    genes, _ = sr.load_geneset(a.geneset)
    full, _ = sr.load_geneset(a.full_geneset)
    edges = sr.read_edges(a.shared)
    clean, ambiguous = sr.load_truth(a.truth)
    rows = list(csv.DictReader(open(a.truth), delimiter="\t"))
    name = {r["Gene ID"]: r["Gene Name"] for r in rows}
    unit = {r["Gene ID"]: r["SD Unit"] for r in rows}
    cn = sr.load_famcn(a.truth, "Median famCN")
    half = sr.load_split(a.split)
    members = defaultdict(set)
    for g, f in clean.items():
        if not f.startswith("Unassigned"):
            members[f].add(g)
    soto = {f: m for f, m in members.items() if len(m) >= 2}
    g2f = {g: f for f, m in soto.items() for g in m}
    cover, kept, leaf_of = sr.pair_families(edges, genes, full, {}, gate=False)
    pred = sr.collapse_cover(cover, leaf_of, genes)
    universe = (full | genes) - ambiguous
    seqfam = defaultdict(set)
    for g in universe:
        if pred.get(g):
            seqfam[pred[g]].add(g)
    clusters = []
    for c in sorted(seqfam):
        fams = sorted({g2f[g] for g in seqfam[c] if g in g2f})
        if len(fams) >= 2:
            gs = sorted(g for f in fams for g in soto[f])
            clusters.append((c, gs, [g2f[g] for g in gs], fams))
    nfam = sum(len(x[3]) for x in clusters)
    assert (len(clusters), nfam) == (33, 87), (len(clusters), nfam)

    col, maxlen = load_duplicons(a.duplicons)
    allg = {g for x in clusters for g in x[1]}
    comp = compositions(a.exons, col, maxlen, allg)
    empty = sorted(g for g in allg if not comp[g])
    print(f"## KEY=cnduplicon result\n\n{len(clusters)} clusters, {nfam} Soto families, {len(allg)} genes; "
          f"{len(empty)} genes with no exonic duplicon\n")

    def assign_half(x):
        _, gs, labels, fams = x
        big = max(fams, key=lambda f: (len(soto[f]), f))
        hs = [half.get(g) for g in soto[big]]
        return "heldout" if hs.count("heldout") > hs.count("dev") else "dev"

    def report(title, simf, subset):
        rng = random.Random(SEED)
        cl = [(gs, labels) for (_, gs, labels, _) in subset]
        per, pooled, p, n_used = run_test(cl, simf, rng)
        v, frac = verdict(per, p)
        print(f"### {title}\n\nclusters used {n_used}; pooled delta {pooled:+.4f}; pooled p {p:.4f}; clusters with delta > 0: "
              f"{frac:.2f} -> **{v}**\n")
        return per, v

    simc = lambda x, y: wjac(comp[x], comp[y])
    dev = [x for x in clusters if assign_half(x) == "dev"]
    ho = [x for x in clusters if assign_half(x) == "heldout"]
    per_all, v_all = report("PRIMARY, all clusters (descriptive)", simc, clusters)
    _, v_dev = report(f"PRIMARY, DEV half ({len(dev)} clusters)", simc, dev)
    _, v_ho = report(f"PRIMARY, HELD-OUT half ({len(ho)} clusters)", simc, ho)
    final = v_ho if v_ho == v_dev else f"SPLIT (dev {v_dev}, held-out {v_ho})"
    print(f"**VERDICT (prereg §4): {final}**\n")

    if a.out_tsv:
        with open(a.out_tsv, "w") as fh:
            fh.write("half\tn_genes\tsoto_families\tdelta\tp\texample_genes\n")
            for x, (d, pv) in zip(clusters, per_all):
                c, gs, labels, fams = x
                fh.write(f"{assign_half(x)}\t{len(gs)}\t{','.join(fams)}\t{'' if d is None else f'{d:.4f}'}\t"
                         f"{'' if pv is None else f'{pv:.4f}'}\t{', '.join(sorted({name[g] for g in gs})[:4])}\n")
    print("| cluster genes | Soto families (sizes) | delta | p | genes, example |\n|---|---|---|---|---|")
    for (c, gs, labels, fams), (d, pv) in zip(clusters, per_all):
        sizes = ",".join(f"{f}({len(soto[f])})" for f in fams)
        ex = ", ".join(sorted({name[g] for g in gs})[:4])
        dd = "n/a" if d is None else f"{d:+.3f}"
        pp = "n/a" if pv is None else f"{pv:.3f}"
        print(f"| {len(gs)} | {sizes} | {dd} | {pp} | {ex} |")
    print()

    # secondary 1: S1C SD Unit labels as sets ('.' dropped)
    lab = {g: set(u.split(",")) - {".", ""} for g, u in unit.items()}
    sub = []
    for c, gs, labels, fams in clusters:
        keep = [k for k, g in enumerate(gs) if lab.get(g)]
        if len(keep) >= 3 and len({labels[k] for k in keep}) >= 2:
            sub.append((c, [gs[k] for k in keep], [labels[k] for k in keep], fams))
    if sub:
        report(f"SECONDARY, S1C SD Unit labels ({len(sub)} clusters with >= 3 labelled genes)",
               lambda x, y: set_jac(lab[x], lab[y]), sub)
    # secondary 2: without genes lacking an exonic duplicon
    sub2 = []
    for c, gs, labels, fams in clusters:
        keep = [k for k, g in enumerate(gs) if comp[g]]
        if len(keep) >= 3 and len({labels[k] for k in keep}) >= 2:
            sub2.append((c, [gs[k] for k in keep], [labels[k] for k in keep], fams))
    report(f"SECONDARY, genes with an exonic duplicon only ({len(sub2)} clusters)", simc, sub2)
    # secondary 3: copy-number gap vs duplicon dissimilarity
    xs_all, ys_all, xs_x, ys_x = [], [], [], []
    for c, gs, labels, fams in clusters:
        for i in range(len(gs)):
            for j in range(i + 1, len(gs)):
                gi, gj = gs[i], gs[j]
                if gi in cn and gj in cn:
                    dcn, dis = abs(cn[gi] - cn[gj]), 1 - simc(gi, gj)
                    xs_all.append(dcn)
                    ys_all.append(dis)
                    if labels[i] != labels[j]:
                        xs_x.append(dcn)
                        ys_x.append(dis)
    print(f"### SECONDARY, copy-number gap vs duplicon dissimilarity\n\nSpearman, all pairs in the clusters: "
          f"{spearman(xs_all, ys_all):+.3f} (n={len(xs_all)}); pairs across a Soto boundary: "
          f"{spearman(xs_x, ys_x):+.3f} (n={len(xs_x)})\n")


if __name__ == "__main__":
    main()
