#!/usr/bin/env python3
"""Pre-registered KEY=unitcovercn (docs/archive/2026-09/PREREG_unit_cover_cn_2026-09-30.md): finer units for the family cover.

Arms: A = KEY=unitcover (any-overlap usage, duplicon units; must reproduce soto_m2_unit_cover.py), B = sliver-free usage (a gene uses
a duplicon only if it is the dominant duplicon of one of its exons), C = duplicon x copy-number class units, D = B + C. Copy-number
classes: a duplicon's single-family genes with a famCN, single-linkage at Soto's |difference| < 2; a gene with famCN uses the classes
holding a member within < 2 of it, a gene without famCN uses every class. Each arm with S1C famCN and with our famCN.

    python3 bench/soto_m2/soto_m2_unit_cover_cn.py --data families.json --duplicons dupmasker_colors.bed
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


def argmax_family(cnt):
    cnt = +cnt
    return min(cnt.items(), key=lambda kv: (-kv[1], fam_num(kv[0])))[0] if cnt else None


def cn_classes(members):
    """members: [(cn, family, bp)] -> list of classes, each a list of members, cut where consecutive famCN differ by >= 2."""
    ms = sorted(members, key=lambda m: m[0])
    out = []
    for m in ms:
        if out and m[0] - out[-1][-1][0] < 2:
            out[-1].append(m)
        else:
            out.append([m])
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True)
    ap.add_argument("--duplicons", required=True)
    ap.add_argument("--split", default=os.path.join(SOTO, "soto_split_2026-09-29.tsv"))
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260930)
    ap.add_argument("--out-tsv")
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
        ex = merge(g["x"])
        cp = comp(g["c"], ex)
        dom = set()
        for x, y in ex:
            ce = comp(g["c"], [[x, y]])
            if ce:
                dom.add(min(ce.items(), key=lambda kv: (-kv[1], kv[0]))[0])
        genes.append(dict(id=g["id"], name=g["n"], c=g["c"], s=g["s"], e=g["e"], half=half.get(g["id"]),
                          S={f for f in g["sf"] if not f.startswith("Unassigned")}, comp=cp, dom=dom,
                          cn={"s1c": g.get("cs"), "ours": g.get("co")}))
    clean = [g for g in genes if len(g["S"]) == 1]
    multi = [g for g in genes if len(g["S"]) >= 2]

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
        g["Ploc"] = p_loc(g)

    def run(arm, cnsrc):
        sliver_free = arm in ("B", "D")
        use_cn = arm in ("C", "D")
        used = lambda g: [d for d in g["comp"] if (d in g["dom"] or not sliver_free)]
        # votes per duplicon from clean genes, in KEY=unitcover's order (clean genes in table order, duplicons in exon order)
        votes = collections.defaultdict(list)          # D -> [(gene, family, bp, cn)]
        for g in clean:
            F = next(iter(g["S"]))
            for d in used(g):
                votes[d].append((g, F, g["comp"][d], g["cn"][cnsrc]))

        def units_of_D(d, minus=None):
            """[(unit key, owner, [cn values])] for duplicon d, without gene `minus`."""
            vs = [v for v in votes.get(d, []) if v[0] is not minus]
            if not vs:
                return []
            withcn = [(v[3], v[1], v[2]) for v in vs if v[3] is not None]
            if not use_cn or not withcn:
                cnt = collections.Counter()
                for v in vs:
                    cnt[v[1]] += v[2]
                o = argmax_family(cnt)
                return [((d, "all"), o, None)] if o else []
            out = []
            for k, cls in enumerate(cn_classes(withcn)):
                cnt = collections.Counter()
                for cn, F, bp in cls:
                    cnt[F] += bp
                o = argmax_family(cnt)
                if o:
                    out.append(((d, k), o, [m[0] for m in cls]))
            return out

        def gene_units(g, minus=None):
            out = []
            c = g["cn"][cnsrc]
            for d in used(g):
                for key, o, cns in units_of_D(d, minus):
                    if cns is None or c is None or any(abs(c - x) < 2 for x in cns):
                        out.append((key, o))
            return out

        owner = {}
        for d in votes:
            for key, o, _ in units_of_D(d):
                owner[key] = o
        mu = {g["id"]: [k for k, _ in gene_units(g)] for g in multi}
        pred = {g["id"]: {owner[k] for k in mu[g["id"]] if k in owner} for g in multi}

        def means(pr):
            out = {}
            for h in ("dev", "heldout"):
                xs = [jacc(g["S"], pr[g["id"]]) for g in multi if g["half"] == h]
                out[h] = (sum(xs) / len(xs), len(xs))
            return out
        obs = means(pred)
        rng = random.Random(a.seed)
        keys, owns = list(owner), [owner[k] for k in owner]
        ge = {h: 0 for h in obs}
        for _ in range(a.perms):
            rng.shuffle(owns)
            perm = dict(zip(keys, owns))
            pm = means({gid: {perm[k] for k in ks if k in perm} for gid, ks in mu.items()})
            for h in ge:
                ge[h] += pm[h][0] >= obs[h][0]
        pval = {h: (1 + ge[h]) / (a.perms + 1) for h in ge}
        ex = sum(g["S"] == pred[g["id"]] for g in multi) / len(multi)
        rec = sum(len(g["S"] & pred[g["id"]]) / len(g["S"]) for g in multi) / len(multi)
        prs = [len(g["S"] & pred[g["id"]]) / len(pred[g["id"]]) for g in multi if pred[g["id"]]]
        ok = mul = emp = 0
        for g in clean:
            P = {o for _, o in gene_units(g, minus=g)}
            F = next(iter(g["S"]))
            ok += P == {F}
            mul += len(P) >= 2
            emp += not P
        return dict(obs=obs, p=pval, exact=ex, recall=rec, precision=sum(prs) / len(prs), clean_ok=ok / len(clean),
                    clean_multi=mul / len(clean), clean_empty=emp / len(clean), pred=pred, n_units=len(owner))

    locm = {}
    for h in ("dev", "heldout"):
        xs = [jacc(g["S"], g["Ploc"]) for g in multi if g["half"] == h]
        locm[h] = sum(xs) / len(xs)
    print(f"# KEY=unitcovercn\nmulti genes {len(multi)}, clean {len(clean)}; location baseline mean Jaccard dev {locm['dev']:.3f} / "
          f"held-out {locm['heldout']:.3f}\n")
    print("| arm | famCN | units | Jaccard dev / held-out | p dev / held-out | exact | recall | precision | clean: own family / "
          "2+ families / empty |")
    print("|---|---|---|---|---|---|---|---|---|")
    res = {}
    for arm in ("A", "B", "C", "D"):
        for src in (("s1c", "ours") if arm in ("C", "D") else ("s1c",)):
            r = run(arm, src)
            res[(arm, src)] = r
            print(f"| {arm} | {src if arm in ('C', 'D') else '-'} | {r['n_units']} | {r['obs']['dev'][0]:.3f} / {r['obs']['heldout'][0]:.3f} | "
                  f"{r['p']['dev']:.4f} / {r['p']['heldout']:.4f} | {r['exact']:.3f} | {r['recall']:.3f} | {r['precision']:.3f} | "
                  f"{r['clean_ok']:.3f} / {r['clean_multi']:.3f} / {r['clean_empty']:.3f} |")
    A, D = res[("A", "s1c")], res[("D", "s1c")]
    assert abs(A["obs"]["heldout"][0] - 0.473) < 0.0005 and abs(A["obs"]["dev"][0] - 0.439) < 0.0005, "arm A does not reproduce KEY=unitcover"
    up_j = D["obs"]["heldout"][0] > A["obs"]["heldout"][0]
    up_s = D["clean_multi"] < A["clean_multi"]
    sig = D["p"]["heldout"] < 0.01
    verdict = "REFINES" if sig and up_j and up_s else ("PARTIAL" if sig and (up_j or up_s) else "NO GAIN")
    print(f"\nVERDICT (arm D, S1C famCN, held-out): {verdict}  [Jaccard {D['obs']['heldout'][0]:.3f} vs A {A['obs']['heldout'][0]:.3f}; "
          f"clean false-multi {D['clean_multi']:.3f} vs A {A['clean_multi']:.3f}; p {D['p']['heldout']:.4f}]")
    npip = {f"ID_{i}" for i in range(149, 156)}
    srt = lambda s: ", ".join(sorted(s, key=fam_num)) or "-"
    print("\nNPIP side, arm D (S1C | ours famCN):")
    for g in multi:
        if g["S"] & npip:
            print(f"- {g['name']}: Soto {srt(g['S'])} | D/s1c {srt(D['pred'][g['id']])} | D/ours {srt(res[('D', 'ours')]['pred'][g['id']])}")
    if a.out_tsv:
        with open(a.out_tsv, "w") as out:
            out.write("gene_id\tname\thalf\tsoto\tA\tB\tC_s1c\tD_s1c\tC_ours\tD_ours\n")
            for g in sorted(multi, key=lambda g: (g["c"], g["s"])):
                cols = [srt(res[k]["pred"][g["id"]]).replace(", ", ",") for k in
                        (("A", "s1c"), ("B", "s1c"), ("C", "s1c"), ("D", "s1c"), ("C", "ours"), ("D", "ours"))]
                out.write("\t".join([g["id"], g["name"], str(g["half"]), srt(g["S"]).replace(", ", ",")] + cols) + "\n")


if __name__ == "__main__":
    main()
