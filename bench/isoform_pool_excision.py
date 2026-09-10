#!/usr/bin/env python3
"""P4 of PREREG isoform_pool: for each copy k that received isoform assignments, re-pool the reads under the
dump made WITHOUT copy k and count how many of those isoform groups abstain (not `assigned` to any copy).
  python3 bench/isoform_pool_excision.py --pooled isoforms_pooled.tsv --dir excise_all/pool --assign-dir excise_all/pool \
      --remap-dir excise_all --bam hsa16.bam --copies copies16.tsv
"""
import argparse, csv, os, subprocess, sys
from collections import Counter, defaultdict
sys.path.insert(0, os.path.dirname(__file__))
from isoform_pool import binom_tail, introns_of, LR, ALPHA


def pooled_verdicts(dump, assign, bam, remap):
    assign = {r["read_name"]: r for r in csv.DictReader(open(assign), delimiter="\t")}
    ev = {}
    for r in csv.DictReader(open(dump), delimiter="\t"):
        st = assign.get(r["read_name"])
        if st is None or st["origin_rejected"] != "0" or int(st["n_candidates"]) < 2:
            continue
        cands = [remap.get(c, c) for c in r["candidates"].split(",")]
        cols = [c.split(":") for c in r["columns"].split(",") if c] if r["columns"] else []
        pe = {}
        for i in range(len(cands)):
            for j in range(i + 1, len(cands)):
                n = ka = kb = 0
                for pos, o, al in cols:
                    x, y = al[i], al[j]
                    if x == "." or y == "." or x == y: continue
                    n += 1
                    if o == x: ka += 1
                    elif o == y: kb += 1
                if n: pe[(cands[i], cands[j])] = (n, ka, kb)
        ev[r["read_name"]] = (cands, pe)
    chain_of = {}
    out = subprocess.run(["samtools", "view", "-F", "2308", bam], capture_output=True, text=True).stdout
    for ln in out.splitlines():
        f = ln.split("\t", 6)
        if f[0] not in chain_of: chain_of[f[0]] = (f[2],) + introns_of(int(f[3]) - 1, f[5])
    groups = defaultdict(list)
    for name in ev:
        if name in chain_of and len(chain_of[name]) >= 2: groups[chain_of[name]].append(name)
    verdicts = {}
    for key, names in groups.items():
        C = sorted({c for n in names for c in ev[n][0]}, key=int)
        pool = defaultdict(lambda: [0, 0, 0])
        for n in names:
            for (A, B), (nn, ka, kb) in ev[n][1].items():
                p = pool[(A, B)]; p[0] += nn; p[1] += ka; p[2] += kb
        def stats(A, B):
            if (A, B) in pool: n, ka, kb = pool[(A, B)]
            elif (B, A) in pool: n, kb, ka = pool[(B, A)]
            else: return (0, 0, 0)
            return (n, ka, kb)
        def llr(A, B):
            n, ka, kb = stats(A, B); return LR * (ka - kb)
        def worst(A): return min((llr(A, B) for B in C if B != A), default=float("inf"))
        if not C: continue
        bk = max(C, key=lambda A: (worst(A), sum(llr(A, B) for B in C if B != A), -int(A)))
        thr = ALPHA / max(1, len(C) - 1); k0 = False; p_read = 0.0; margin = float("inf")
        for B in C:
            if B == bk: continue
            n, ka, kb = stats(bk, B)
            if n == 0: k0 = True
            p_read = max(p_read, binom_tail(n, ka)); margin = min(margin, llr(bk, B))
        verdicts[key] = ("tied" if k0 or len(C) < 2 else ("assigned" if p_read < thr and margin > 0 else "ambiguous"), bk, margin, len(names))
    return verdicts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pooled", required=True); ap.add_argument("--dir", required=True); ap.add_argument("--remap-dir", required=True)
    ap.add_argument("--bam", required=True)
    a = ap.parse_args()
    base = {}
    for r in csv.DictReader(open(a.pooled), delimiter="\t"):
        key = (r["chrom"],) + tuple(tuple(int(x) for x in s.split("-")) for s in r["chain"].split(";") if s)
        base[key] = r
    tot = abst = 0; detail = Counter()
    for k in sorted({r["isoform_copy"] for r in base.values() if r["isoform_status"] == "assigned"}, key=int):
        remap = {}
        for l in open(f"{a.remap_dir}/remap_no{k}.tsv"):
            new, old = l.split(); remap[new] = old
        v = pooled_verdicts(f"{a.dir}/no{k}.star_reads.tsv", f"{a.dir}/no{k}.assignments.tsv", a.bam, remap)
        mine = [key for key, r in base.items() if r["isoform_status"] == "assigned" and r["isoform_copy"] == k]
        for key in mine:
            tot += 1
            nv = v.get(key)
            if nv is None or nv[0] != "assigned":
                abst += 1; detail[f"no{k}: " + (nv[0] if nv else "group gone")] += 1
            else:
                detail[f"no{k}: assigned -> {nv[1]} (m={nv[2]:.0f})"] += 1
        print(f"copy {k}: {len(mine)} isoform assignments; abstain without the copy: {sum(1 for key in mine if v.get(key) is None or v[key][0] != 'assigned')}")
    print(f"P4 pooled excision: {abst}/{tot} = {100*abst/max(1,tot):.1f}% abstain  (pass >=90, refuted <80)")
    print("   detail:", sorted(detail.items()))


if __name__ == "__main__":
    main()
