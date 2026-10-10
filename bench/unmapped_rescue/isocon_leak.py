#!/usr/bin/env python3
"""Amendment 47: IsoCon on leaked-transcript and real-flag loci. Miniforge python.
    isocon_leak.py select     -> W/isocon_leak/loci.json + one reads FASTA per locus
    isocon_leak.py score      -> per-locus outcome after iso_leak.sh ran IsoCon"""
import json
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import locus_support as L  # noqa: E402
import run_polish as RP  # noqa: E402
import seeds as SD  # noqa: E402
import spec_heldout as SH  # noqa: E402

O = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"
W = f"{O}/isocon_leak"
RCT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
D = 0.00958


def select():
    import collections
    import csv
    import pysam
    os.makedirs(f"{W}/reads", exist_ok=True)
    lab = json.load(open(f"{O}/dna_labels.json"))
    rows = {r["k"]: r for r in json.load(open(f"{O}/net_run_binned/classes.json"))}
    S = {"own": SH.load(f"{O}/spec_flags.cs.paf"), "human": SH.load(f"{O}/spec_dev.hsa.paf"), "chimp": SH.load(f"{O}/spec_dev.ptr.paf"),
         "orangutan": SH.load(f"{O}/spec_dev.ppy.paf"), "siamang": SH.load(f"{O}/spec_dev.ssy.paf")}
    best = lambda k: max(S, key=lambda s: S[s].get(k, 0))
    bins = json.load(open(f"{O}/net_run_binned/bins.json"))
    bin_of = {n: i for i, b in enumerate(bins) for n in b}
    cl = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{O}/net_run_binned/clusters.tsv"), delimiter="\t"):
        cl[r["cluster"]].append(r["read"])
    dev = {k: x for k, x in lab.items() if x["side"] == "dev" and rows[k]["R_hit"]}
    leak = sorted(k for k, x in dev.items() if x["verdict"] == "WRONG" and x["lab_all"] == "RNA-ONLY" and best(k) != "own")
    true = sorted(k for k, x in dev.items() if x["verdict"] == "TRUE")
    rnd = random.Random(47)
    rnd.shuffle(leak)
    rnd.shuffle(true)
    pick, used = [], set()
    for group, ks, n in (("leak", leak, 40), ("true", true, 20)):
        got = 0
        for k in ks:
            b = bin_of.get(rows[k]["key"])
            if b in used or got >= n:
                continue
            used.add(b)
            pick.append(dict(k=k, group=group, best=best(k)))
            got += 1
    net = SD.read_fa(f"{O}/net_all.fa")
    bam = pysam.AlignmentFile(L.BAM, "rb")
    for p in pick:
        k = p["k"]
        mine = cl[rows[k]["key"]]
        c, s, e = rows[k]["R_hit"]
        others, seqs = [], {}
        for a in bam.fetch(c, s, e):
            if a.is_secondary or a.is_supplementary or a.query_sequence is None or a.query_name in mine:
                continue
            q = a.query_sequence
            seqs[a.query_name] = q.translate(RCT)[::-1] if a.is_reverse else q
            others.append(a.query_name)
        keep = L.sample_every(others, max(0, 300 - len(mine)))
        for n in mine:
            seqs[n] = net[n]
        SD.write_fa(f"{W}/reads/{k}.fa", seqs, list(mine) + keep)
        p.update(flag_reads=len(mine), other_reads=len(keep))
    json.dump(pick, open(f"{W}/loci.json", "w"), indent=1)
    print(f"selected {sum(p['group'] == 'leak' for p in pick)} leak + {sum(p['group'] == 'true' for p in pick)} true loci; reads per locus median "
          f"{sorted(p['flag_reads'] + p['other_reads'] for p in pick)[len(pick) // 2]}")


def score():
    pick = json.load(open(f"{W}/loci.json"))
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{O}/net_run_binned/cons.fa").items()}
    pv = {}
    for side in ("dev", "held"):
        pv.update(SD.read_fa(f"{O}/privers.{side}.fa"))
    res = {"leak": [0, 0, 0], "true": [0, 0, 0]}
    rows = []
    for p in pick:
        k = p["k"]
        fc = f"{W}/iso/{k}/final_candidates.fa"
        if not os.path.exists(fc) or os.path.getsize(fc) == 0:
            res[p["group"]][2] += 1
            rows.append(dict(k=k, group=p["group"], ran=False))
            continue
        cands = SD.read_fa(fc)
        SD.write_fa(f"{W}/t_{k}.fa", {"flag": cons[k], "ref": pv[k]}, ["flag", "ref"])
        lines = subprocess.run(f"minimap2 -c -x map-hifi -N 5 -t 2 {W}/t_{k}.fa {fc}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout.splitlines()
        os.remove(f"{W}/t_{k}.fa")
        best = {}
        for ln in lines:
            f = ln.split("\t")
            ident = int(f[9]) / int(f[10])
            cov = int(f[10]) / min(int(f[1]), int(f[6]))
            if cov < 0.5:
                continue
            key = (f[0], f[5])
            best[key] = max(best.get(key, 0), ident)
        hit = any(best.get((c, "flag"), 0) >= 0.99 and best.get((c, "flag"), 0) - best.get((c, "ref"), 0) > D for c in cands)
        res[p["group"]][0 if hit else 1] += 1
        rows.append(dict(k=k, group=p["group"], ran=True, candidates=len(cands), reports_flag=hit, flag_reads=p["flag_reads"], other_reads=p["other_reads"]))
    json.dump(rows, open(f"{W}/scores.json", "w"), indent=1)
    for g, (h, m, nr) in res.items():
        print(f"{g}: IsoCon reports the flagged transcript at {h} of {h + m} loci ({h / max(1, h + m):.0%}); IsoCon did not finish / no output at {nr}")


if __name__ == "__main__":
    {"select": select, "score": score}[sys.argv[1]]()
