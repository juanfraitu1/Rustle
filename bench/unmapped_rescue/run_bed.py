#!/usr/bin/env python3
"""Cluster the pool of a bed and score it (docs/PREREG_unmapped_rescue_2026-10-08.md sections 1-2, 4).

    run_bed.py <bed> [--delta D] [--tag T]    # bed = bedA | bedH | bedM; resumable: exit 75 = run again (one fresh minimap2 round per call)
Writes W/<bed>/<tag>/clusters.tsv (read, cluster, size) and metrics.json.
"""
import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import graph as G  # noqa: E402
import score as S  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
DELTA = 0.00958


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bed")
    ap.add_argument("--delta", type=float, default=DELTA)
    ap.add_argument("--tag", default="registered")
    ap.add_argument("--seeds", type=int, default=2000)
    a = ap.parse_args()
    d = f"{W}/{a.bed}"
    out = f"{d}/{a.tag}"
    os.makedirs(f"{out}/rounds", exist_ok=True)
    seqs = SD.read_fa(f"{d}/pool.fa")
    lens = {n: len(s) for n, s in seqs.items()}
    lab = {r["read"]: r for r in csv.DictReader(open(f"{d}/labels.tsv"), delimiter="\t")}
    truth = {n: ("bg" if lab[n]["role"] == "bg" else (lab[n]["family"] or None)) for n in seqs}
    dread = {n for n in seqs if lab[n]["role"] == "D"}
    map_fn = SD.minimap_map_fn(seqs, f"{out}/rounds", a.delta, 0.5)
    try:
        comp = SD.run_rounds(lens, map_fn, n_seeds=a.seeds, min_size=3, max_rounds=6)
    except SD.Pause as e:
        print(f"round {e} finished; run again")
        sys.exit(75)
    cl = G.clusters(comp, 3)
    m = S.cluster_metrics(cl, truth, dread)
    maj = {c: S.majority(rs, truth) for c, rs in cl.items()}
    m["background_in_family_clusters"] = sum(1 for c, rs in cl.items() if maj[c] for r in rs if truth.get(r) == "bg")
    m["delta"] = a.delta
    with open(f"{out}/clusters.tsv", "w") as o:
        o.write("read\tcluster\tsize\tmajority\n")
        for c, rs in cl.items():
            for r in rs:
                o.write(f"{r}\t{c}\t{len(rs)}\t{maj[c] or ''}\n")
    json.dump(m, open(f"{out}/metrics.json", "w"), indent=1)
    print({k: v for k, v in m.items() if k != "clusters_per_family"})
    print("clusters per family: min/median/max", min(m["clusters_per_family"].values()), sorted(m["clusters_per_family"].values())[len(m["clusters_per_family"]) // 2], max(m["clusters_per_family"].values()))


if __name__ == "__main__":
    main()
