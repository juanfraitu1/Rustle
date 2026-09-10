#!/usr/bin/env python3
"""PREREG_identity_weighted_density (b4b65ebd): sweep density_weighted = sum(identity)/C(n,2) over every
cluster in a mcl_families --dump-pairs run, sanity-check it against the reported density, correlate the
gap (density - density_weighted) against identity_gap.py's split verdict and the corroborated column.

  python3 bench/identity_weighted_density.py --clusters X.clusters.tsv --pairs X.pairs.tsv \
      [--min-size 4] [--max-size 60] [--tbc1d3-only]
"""
import argparse, csv, math, statistics as st, subprocess, sys, tempfile, os
from collections import defaultdict

IDENTITY_GAP = os.path.join(os.path.dirname(os.path.abspath(__file__)), "identity_gap.py")


def load_clusters(path):
    cl = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        c = cl.setdefault(r["cluster_id"], {"size": int(r["size"]), "density": float(r["density"]),
                                             "frac_in": float(r["frac_in"]),
                                             "corroborated": None if r["corroborated"] == "NA" else float(r["corroborated"]),
                                             "members": []})
        c["members"].append((r["chrom"], r["start"], r["end"]))
    return cl


def load_pairs(path):
    p = defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        p[r["cluster_id"]].append((r["a"], r["b"], float(r["identity"])))
    return p


def run_identity_gap(pairs):
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        for a, b, i in pairs:
            fh.write(f"{a}\t{b}\t{i}\n")
        path = fh.name
    try:
        out = subprocess.run(["python3", IDENTITY_GAP, path], capture_output=True, text=True, timeout=240).stdout
    except subprocess.TimeoutExpired:
        os.unlink(path)
        return None, None, "TIMEOUT"
    os.unlink(path)
    p_worst, split = None, None
    for line in out.splitlines():
        if line.startswith("p = ") and "worst null governs" in line:
            p_worst = float(line.split()[2])
            # "NO SPLIT" contains "SPLIT" as a substring -- check the negative form first.
            split = "NO SPLIT" not in line
    return p_worst, split, out


def mannwhitney_u_onesided(a, b):
    """P(a tends to be larger than b), one-sided Mann-Whitney via normal approximation. Returns (U, p)."""
    all_v = sorted([(x, 0) for x in a] + [(y, 1) for y in b])
    ranks = {}
    i = 0
    n = len(all_v)
    while i < n:
        j = i
        while j < n and all_v[j][0] == all_v[i][0]:
            j += 1
        r = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[k] = r
        i = j
    rank_a = sum(ranks[k] for k, (v, g) in enumerate(all_v) if g == 0)
    na, nb = len(a), len(b)
    U = rank_a - na * (na + 1) / 2.0
    mu = na * nb / 2.0
    sigma = math.sqrt(na * nb * (na + nb + 1) / 12.0)
    if sigma == 0:
        return U, float("nan")
    z = (U - mu) / sigma
    p_one = 1 - 0.5 * (1 + math.erf(z / math.sqrt(2)))
    return U, p_one


def spearman(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    rx = {i: r for r, i in enumerate(sorted(range(n), key=lambda k: xs[k]))}
    ry = {i: r for r, i in enumerate(sorted(range(n), key=lambda k: ys[k]))}
    d2 = sum((rx[i] - ry[i]) ** 2 for i in range(n))
    return 1 - 6 * d2 / (n * (n * n - 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--min-size", type=int, default=4)
    ap.add_argument("--max-size", type=int, default=60)
    ap.add_argument("--holdout-id", help="cluster_id to score separately as the held-out check (P5)")
    a = ap.parse_args()
    clusters = load_clusters(a.clusters)
    pairs = load_pairs(a.pairs)

    # P0: sanity gate
    mism = []
    for cid, c in clusters.items():
        n = c["size"]
        possible = n * (n - 1) / 2
        n_pairs = len(pairs.get(cid, []))
        recomputed = n_pairs / possible if possible else 0.0
        if abs(recomputed - c["density"]) > 0.0002:
            mism.append((cid, n, c["density"], recomputed, n_pairs, possible))
    print(f"P0 sanity gate: {len(clusters) - len(mism)}/{len(clusters)} clusters match (density recomputed from pairs.tsv)")
    for m in mism[:5]:
        print(f"   MISMATCH {m}")

    rows = []  # cid, size, density, density_weighted, gap, corroborated
    for cid, c in clusters.items():
        n = c["size"]
        if not (a.min_size <= n <= a.max_size):
            continue
        ps = pairs.get(cid, [])
        possible = n * (n - 1) / 2
        dens = len(ps) / possible if possible else 0.0
        dw = sum(i for _, _, i in ps) / possible if possible else 0.0
        rows.append(dict(cid=cid, size=n, density=dens, dw=dw, gap=dens - dw, corroborated=c["corroborated"], pairs=ps))
    print(f"\n{len(rows)} clusters in [{a.min_size},{a.max_size}] with >= 6 pairs (identity_gap.py's floor): "
          f"{sum(1 for r in rows if len(r['pairs']) >= 6)}")
    rows = [r for r in rows if len(r["pairs"]) >= 6]

    print(f"\nrunning identity_gap.py on all {len(rows)} clusters...")
    for r in rows:
        p, split, _ = run_identity_gap(r["pairs"])
        r["p_gap"] = p
        r["split"] = split

    scored = [r for r in rows if r["p_gap"] is not None]
    split_gaps = [r["gap"] for r in scored if r["split"]]
    nosplit_gaps = [r["gap"] for r in scored if not r["split"]]
    print(f"\nP1: {len(split_gaps)} SPLIT clusters, {len(nosplit_gaps)} NO-SPLIT clusters")
    if split_gaps and nosplit_gaps:
        print(f"   median gap: split={st.median(split_gaps):.4f}  no-split={st.median(nosplit_gaps):.4f}")
        U, p = mannwhitney_u_onesided(split_gaps, nosplit_gaps)
        print(f"   Mann-Whitney U={U:.1f}, one-sided p(split > no-split)={p:.4f}  (pass p<0.10 AND split median > no-split median)")

    print(f"\nP2: spearman(size, gap) = {spearman([r['size'] for r in rows], [r['gap'] for r in rows]):.3f}  (pass < 0.4, refuted >= 0.6)")

    withc = [r for r in rows if r["corroborated"] is not None]
    lowc = [r for r in withc if r["corroborated"] < 0.5]
    highc = [r for r in withc if r["corroborated"] >= 0.9]
    print(f"\nP3: {len(lowc)} clusters corroborated<0.5, {len(highc)} corroborated>=0.9")
    if lowc and highc:
        ratio_low = st.median([r["dw"] / r["density"] if r["density"] else 1.0 for r in lowc])
        ratio_high = st.median([r["dw"] / r["density"] if r["density"] else 1.0 for r in highc])
        print(f"   median dw/density: low-corrob={ratio_low:.3f}  high-corrob={ratio_high:.3f}  diff={ratio_high-ratio_low:.3f}  (pass >=0.05, refuted <0.02 or wrong sign)")

    print("\nP4: TBC1D3 cluster (search by presence of the known coordinate)")
    tbc = None
    for r in rows + [rr for rr in [] ]:
        pass
    # search across ALL clusters (any size) for the one containing the TBC1D3 gene coordinate
    for cid, c in clusters.items():
        if any(m[0] == "NC_073244.2" and m[1] == "39044723" for m in c["members"]) or \
           any(abs(int(m[1]) - 37113575) < 5000 and m[0] == "NC_073244.2" for m in c["members"]):
            tbc = cid
            break
    if tbc:
        n = clusters[tbc]["size"]; ps = pairs.get(tbc, [])
        possible = n * (n - 1) / 2
        dw = sum(i for _, _, i in ps) / possible if possible else 0.0
        print(f"   found cluster {tbc} (size {n}): density={clusters[tbc]['density']:.3f} density_weighted={dw:.3f}  "
              f"(pass: within 0.08 of artifact's 0.724 -> |{dw:.3f}-0.724|={abs(dw-0.724):.3f})")
    else:
        print("   TBC1D3 not found in this run's clusters (coordinates may differ from the artifact's chr17 human locus -- this IS gorilla, check separately)")

    if a.holdout_id:
        r = next((r for r in rows if r["cid"] == a.holdout_id), None)
        if r:
            print(f"\nP5 held-out {a.holdout_id}: size={r['size']} gap={r['gap']:.4f} predicted={'SPLIT' if r['gap']>st.median([x['gap'] for x in rows]) else 'NO SPLIT'} actual={'SPLIT' if r['split'] else 'NO SPLIT'} p={r['p_gap']}")

    print("\nfull table (size, density, density_weighted, gap, corroborated, identity_gap p, split):")
    for r in sorted(rows, key=lambda r: -r["gap"]):
        print(f"  {r['cid']:6s} n={r['size']:3d} dens={r['density']:.3f} dw={r['dw']:.3f} gap={r['gap']:.3f} "
              f"corrob={r['corroborated']}  p_gap={r['p_gap']}  split={r['split']}")


if __name__ == "__main__":
    main()
