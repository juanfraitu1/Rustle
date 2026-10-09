#!/usr/bin/env python3
"""Amendment 43: features of every flagged cluster of the binned run; prints DEVELOPMENT-set statistics only (the held-out set stays unseen). Miniforge python."""
import collections
import csv
import glob
import json
import os
import statistics as st
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import discover as D  # noqa: E402
import run_polish as RP  # noqa: E402
import seeds as SD  # noqa: E402
import specfeat as F  # noqa: E402

O = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"
R = f"{O}/net_run_binned"


def split_of():
    """bin contig of every net read -> 'dev' (odd chromosomes, unplaced) or 'held' (even, X, Y, chrM, unmapped)"""
    num = {ln.split("\t")[0]: ln.split("\t")[1] for ln in list(open("/mnt/linuxdisk/tmp/rna_allele/chrmap.tsv"))[1:]}
    side = {}
    for fn in glob.glob(f"{O}/netpos/*.tsv"):
        c = os.path.basename(fn)[:-4]
        n = num.get(c)
        s = "dev" if (n is None and c.startswith("NW_")) or (n is not None and n.isdigit() and int(n) % 2 == 1) else "held"
        for ln in open(fn):
            side[ln.split("\t", 1)[0]] = s
    return side


def main():
    ev = {v["k"]: v for v in json.load(open(f"{O}/eval_binned.json"))["verdicts"]}
    rows = {r["k"]: r for r in json.load(open(f"{R}/classes.json"))}
    side = split_of()
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{R}/cons.fa").items()}
    fa = f"{O}/spec_flags.fa"
    SD.write_fa(fa, {k: cons[k] for k in ev}, sorted(ev))
    paf = f"{O}/spec_flags.cs.paf"
    if not os.path.exists(paf + ".done"):
        subprocess.run(f"minimap2 -c --cs=long -x splice:hq -uf -t 4 {D.RP.GENOME_IDX} {fa} > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
        open(paf + ".done", "w").write("ok")
    best = {}
    for ln in open(paf):
        f = ln.rstrip("\n").split("\t")
        if "tp:A:P" in f[12:] and (f[0] not in best or int(f[9]) > int(best[f[0]][9])):
            best[f[0]] = f
    feats = {}
    for k, v in ev.items():
        key = rows[k]["key"]
        s = side.get(key, "held")                       # the centre read; unmapped reads have no bin contig -> held-out
        f = best.get(k)
        x = dict(k=k, side=s, verdict=v["verdict"], cls=v["cls"], reads=v["reads"], length=v["length"], R=rows[k]["R"], med=rows[k]["median_read_divergence"])
        if f:
            cs = next(t[5:] for t in f[12:] if t.startswith("cs:Z:"))
            x.update(F.cs_features(cs))
            x["qcov"] = (int(f[3]) - int(f[2])) / int(f[1])
            x["head"], x["tail"] = int(f[2]), int(f[1]) - int(f[3])
        feats[k] = x
    json.dump(feats, open(f"{O}/spec_features.json", "w"))
    dev = [x for x in feats.values() if x["side"] == "dev" and "match" in x]
    print(f"flags: {len(feats)}; development {sum(x['side'] == 'dev' for x in feats.values())} (aligned {len(dev)}), held-out {sum(x['side'] == 'held' for x in feats.values())}")
    for verdict in ("TRUE", "WRONG"):
        sel = [x for x in dev if x["verdict"] == verdict]
        aln = lambda x: x["match"] + x["sub"] + x["ins_bp"] + x["del_bp"]
        m = lambda fn: round(st.median(fn(x) for x in sel), 4)
        print(f"DEV {verdict:5s} n {len(sel):4d} | per 1 kb aligned: subs {m(lambda x: 1000 * x['sub'] / aln(x))}  indel events {m(lambda x: 1000 * (x['ins'] + x['dele']) / aln(x))}  "
              f"indel bp {m(lambda x: 1000 * (x['ins_bp'] + x['del_bp']) / aln(x))} | share of indels in homopolymers {m(lambda x: x['hp_indel'] / max(1, x['ins'] + x['dele']))} | "
              f"a>g share of subs {m(lambda x: x['ag'] / max(1, x['sub']))} | query covered {m(lambda x: x['qcov'])} head {m(lambda x: x['head'])} tail {m(lambda x: x['tail'])} | reads {m(lambda x: x['reads'])} "
              f"| read-to-consensus divergence {m(lambda x: x['med'] or 0)} | R {m(lambda x: x['R'] or 0)}")


if __name__ == "__main__":
    main()
