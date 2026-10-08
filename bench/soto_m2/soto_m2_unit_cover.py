#!/usr/bin/env python3
"""Pre-registered KEY=unitcover (docs/archive/2026-09/PREREG_unit_cover_2026-09-30.md): families as a partition of duplicon units, genes as paths.

Owner of a duplicon = the Soto family whose single-family ("clean") genes have the most exonic bases on it (ties: lower family
number; leave-one-out when the predicted gene is itself clean). A gene's predicted families = the owners of the duplicons its exons
touch. Scored against Soto's multi-family genes (mean Jaccard per split half) with a permutation null and a location baseline.

    python3 bench/soto_m2/soto_m2_unit_cover.py --data families.json --duplicons dupmasker_colors.bed --out-tsv unit_cover.tsv
"""
import argparse
import bisect
import collections
import csv
import json
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
SOTO = os.path.join(HERE, "..", "soto")


def merge(ivs):
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def fam_num(f):
    return int(f.split("_")[1])


def jacc(s, p):
    u = s | p
    return len(s & p) / len(u) if u else 0.0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--duplicons", required=True)
    ap.add_argument("--split", default=os.path.join(SOTO, "soto_split_2026-09-29.tsv"))
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260930)
    a = ap.parse_args(argv)

    G = json.load(open(a.data))["genes"]
    half = {r["gene_id"]: r["half"] for r in csv.DictReader(open(a.split), delimiter="\t")}
    chroms = {g["c"] for g in G}
    segs = collections.defaultdict(list)
    for line in open(a.duplicons):
        f = line.split("\t")
        if f[0] in chroms:
            segs[f[0]].append((int(f[1]), int(f[2]), f[3]))
    starts, maxlen = {}, {}
    for c in segs:
        segs[c].sort()
        starts[c] = [s for s, _, _ in segs[c]]
        maxlen[c] = max(e - s for s, e, _ in segs[c])

    def comp(c, exons):
        out = collections.Counter()
        for x, y in exons:
            k = bisect.bisect_left(starts.get(c, []), x - maxlen.get(c, 0))
            while k < len(segs.get(c, [])) and segs[c][k][0] < y:
                s, e, d = segs[c][k]
                o = min(y, e) - max(x, s)
                if o > 0:
                    out[d] += o
                k += 1
        return out

    genes = []
    for g in G:
        S = {f for f in g["sf"] if not f.startswith("Unassigned")}
        genes.append(dict(id=g["id"], name=g["n"], c=g["c"], s=g["s"], e=g["e"], S=S, comp=comp(g["c"], merge(g["x"])),
                          half=half.get(g["id"])))
    clean = [g for g in genes if len(g["S"]) == 1]
    multi = [g for g in genes if len(g["S"]) >= 2]
    tot = collections.defaultdict(collections.Counter)    # duplicon -> family -> clean exonic bp
    for g in clean:
        F = next(iter(g["S"]))
        for d, bp in g["comp"].items():
            tot[d][F] += bp

    def owner(d, minus=None):
        cnt = tot.get(d)
        if not cnt:
            return None
        if minus:
            F, bp = minus
            cnt = cnt.copy()
            cnt[F] -= bp
            cnt = +cnt
            if not cnt:
                return None
        return min(cnt.items(), key=lambda kv: (-kv[1], fam_num(kv[0])))[0]

    own = {d: owner(d) for d in tot}
    own = {d: f for d, f in own.items() if f}

    def predict(g, owners):
        return {owners[d] for d in g["comp"] if d in owners}

    # location baseline: families of clean genes overlapping the gene's span, else the nearest clean gene on the chromosome
    by_chr = collections.defaultdict(list)
    for g in clean:
        by_chr[g["c"]].append(g)

    def p_loc(g):
        cand = [h for h in by_chr[g["c"]] if h is not g]
        ov = {next(iter(h["S"])) for h in cand if h["s"] < g["e"] and g["s"] < h["e"]}
        if ov or not cand:
            return ov
        near = min(cand, key=lambda h: max(h["s"] - g["e"], g["s"] - h["e"], 0))
        return set(near["S"])

    for g in multi:
        g["P"], g["Ploc"] = predict(g, own), p_loc(g)
        g["J"], g["Jloc"] = jacc(g["S"], g["P"]), jacc(g["S"], g["Ploc"])

    def mean_by_half(key_fn):
        out = {}
        for h in ("dev", "heldout"):
            xs = [key_fn(g) for g in multi if g["half"] == h]
            out[h] = (sum(xs) / len(xs) if xs else float("nan"), len(xs))
        return out

    obs = mean_by_half(lambda g: g["J"])
    loc = mean_by_half(lambda g: g["Jloc"])
    rng = random.Random(a.seed)
    ds, fs = list(own), [own[d] for d in own]
    ge = {h: 0 for h in ("dev", "heldout")}
    for _ in range(a.perms):
        rng.shuffle(fs)
        perm = dict(zip(ds, fs))
        pm = mean_by_half(lambda g: jacc(g["S"], predict(g, perm)))
        for h in ge:
            ge[h] += pm[h][0] >= obs[h][0]
    pval = {h: (1 + ge[h]) / (a.perms + 1) for h in ge}

    def verdict(h):
        if pval[h] >= 0.01:
            return "FAILS"
        return "HOLDS" if obs[h][0] > loc[h][0] else "PARTIAL"

    print(f"# KEY=unitcover\ngenes {len(genes)}: clean {len(clean)}, multi {len(multi)}; owned duplicons {len(own)} "
          f"(of {len(tot)} touched by clean genes); families owning >= 1 duplicon {len(set(own.values()))}\n")
    print("| half | multi genes | mean Jaccard (structure) | mean Jaccard (location baseline) | permutation p | verdict |")
    print("|---|---|---|---|---|---|")
    for h in ("dev", "heldout"):
        print(f"| {h} | {obs[h][1]} | {obs[h][0]:.3f} | {loc[h][0]:.3f} | {pval[h]:.4f} | {verdict(h)} |")
    vh, vd = verdict("heldout"), verdict("dev")
    print(f"\nVERDICT (held-out decides): {vh}" + (f"  [dev: {vd} -> SPLIT]" if vd != vh else f"  [dev agrees: {vd}]"))

    def prf(key):
        ex = sum(g["S"] == g[key] for g in multi) / len(multi)
        rec = sum(len(g["S"] & g[key]) / len(g["S"]) for g in multi) / len(multi)
        pr = [len(g["S"] & g[key]) / len(g[key]) for g in multi if g[key]]
        return ex, rec, (sum(pr) / len(pr) if pr else float("nan")), sum(1 for g in multi if not g[key])
    for key, lab in (("P", "structure"), ("Ploc", "location baseline")):
        ex, rec, pr, emp = prf(key)
        print(f"- {lab}, all {len(multi)} multi genes: exact set {ex:.3f}, recall {rec:.3f}, precision {pr:.3f}, empty {emp}")
    n_ok = n_multi = n_empty = 0
    for g in clean:
        F = next(iter(g["S"]))
        P = {owner(d, (F, g["comp"][d])) for d in g["comp"]} - {None}
        n_ok += P == {F}
        n_multi += len(P) >= 2
        n_empty += not P
    print(f"- clean genes (leave-one-out), {len(clean)}: predicted exactly their family {n_ok / len(clean):.3f}, "
          f"two or more families {n_multi / len(clean):.3f}, empty {n_empty / len(clean):.3f}")
    print("\nNPIP side (ID_149-ID_155):")
    npip = {f"ID_{i}" for i in range(149, 156)}
    for g in multi:
        if g["S"] & npip:
            srt = lambda s: ", ".join(sorted(s, key=fam_num)) or "-"
            print(f"- {g['name']}: Soto {srt(g['S'])} | structure {srt(g['P'])} | location {srt(g['Ploc'])}")
    with open(a.out_tsv, "w") as out:
        out.write("gene_id\tname\thalf\tsoto_families\tpredicted\tlocation_baseline\tjaccard\tjaccard_location\tn_units\n")
        for g in sorted(multi, key=lambda g: (g["c"], g["s"])):
            j = lambda s: ",".join(sorted(s, key=fam_num))
            out.write(f"{g['id']}\t{g['name']}\t{g['half']}\t{j(g['S'])}\t{j(g['P'])}\t{j(g['Ploc'])}\t{g['J']:.3f}\t{g['Jloc']:.3f}"
                      f"\t{len(g['comp'])}\n")


if __name__ == "__main__":
    main()
