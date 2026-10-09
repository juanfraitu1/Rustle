#!/usr/bin/env python3
"""Amendment 37 S2 / S3: reads within delta of the consensus, by abPOA mode, for the clusters of the control runs. Miniforge python.

    consensus_modes_probe.py <class[,class...]> [sample N seed]   writes W/modes_probe_<classes>.json"""
import collections
import csv
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import consensus_support as CS  # noqa: E402
import partition as PT  # noqa: E402
import run_polish as RP  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
D = 0.00958
MODES = ("g", "l", "e")


def support(cons, reads, tmp):
    """-> (median divergence of the reads on the consensus, number of reads within delta)"""
    import run_augment as RA
    import subprocess
    os.makedirs(tmp, exist_ok=True)
    open(f"{tmp}/r.fa", "w").write("".join(f">{n}\n{s}\n" for n, s in reads.items()))
    open(f"{tmp}/c.fa", "w").write(f">c\n{cons}\n")
    sam = subprocess.run(f"minimap2 -ax splice:hq -uf --eqx -t 2 {tmp}/c.fa {tmp}/r.fa", shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout.splitlines()
    prim = RA.primaries(sam)
    des = sorted(prim[n]["de"] for n in reads if n in prim)
    med = des[len(des) // 2] if des else None
    return med, sum(x <= D for x in des)


def main(classes, sample=None, seed=1):
    runs = [a + s for a in ("ggo_testis", "a119b", "ptr", "ppy") for s in ("_control", "_control_s6", "_control_s7")]
    items = []
    for run in runs:
        d = f"{W}/discover_{run}"
        for r in json.load(open(f"{d}/classes.a35c.json")):
            if r["cls"] in classes:
                items.append((run, r["k"], r["reads"], r["cls"]))
    if sample:
        items = random.Random(seed).sample(sorted(items), sample)
    cache = {}
    out = []
    for run, k, n, cls in items:
        d = f"{W}/discover_{run}"
        if run not in cache:
            cl = collections.defaultdict(list)
            for r in csv.DictReader(open(f"{d}/clusters.tsv"), delimiter="\t"):
                cl["cl" + r["cluster"]].append(r["read"])
            cache[run] = (cl, SD.read_fa(open(f"{d}/reads_path.txt").read().strip()),
                          {x.split("|")[0]: s for x, s in RP.read_cons(f"{d}/cons.fa").items()})
        cl, reads, cons = cache[run]
        names = sorted(cl[k])[:100]
        sub = {x: reads[x] for x in names}
        row = dict(run=run, k=k, cls=cls, n=len(names))
        for m in MODES:
            c = cons[k] if m == "g" else PT.abpoa_consensus([sub[x] for x in names], mode=m)
            med, within = support(c, sub, f"/tmp/modes_{run}")
            row[m] = dict(med=med, within=within, cons=c)
        out.append(row)
    tag = "_".join(classes) + (f"_s{sample}_{seed}" if sample else "")
    json.dump(out, open(f"{W}/modes_probe_{tag}.json", "w"))
    tot = sum(r["n"] for r in out)
    print(f"{','.join(classes)}: {len(out)} clusters, {tot} reads")
    for m in MODES:
        w = sum(r[m]["within"] for r in out)
        sup = sum(r[m]["med"] is not None and r[m]["med"] <= D for r in out)
        print(f"  mode {m}: reads within delta {w} ({w / tot:.3f}); clusters supported (median <= delta) {sup} of {len(out)}")


if __name__ == "__main__":
    classes = sys.argv[1].split(",")
    main(classes, int(sys.argv[2]) if len(sys.argv) > 2 else None, int(sys.argv[3]) if len(sys.argv) > 3 else 1)
