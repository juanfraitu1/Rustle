#!/usr/bin/env python3
"""Amendment 23/24 on a real bed: star clustering (every gap counted) inside the components of <= 60 reads of an existing clustering.
    run_chain.py <bedA|bedH> <src_tag> <dst_tag>      e.g. bedH chain chain2   (src = the --proper clustering); then consensus.py and eval_clusters.py"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chain as CH  # noqa: E402
import score as S  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
DELTA = 0.00958


def main(bed, src, dst):
    d = f"{W}/{bed}"
    seqs = SD.read_fa(f"{d}/pool.fa")
    cl = {}
    for r in csv.DictReader(open(f"{d}/{src}/clusters.tsv"), delimiter="\t"):
        cl.setdefault(r["cluster"], []).append(r["read"])
    lab = {r["read"]: r for r in csv.DictReader(open(f"{d}/labels.tsv"), delimiter="\t")}
    truth = {n: ("bg" if lab[n]["role"] == "bg" else (lab[n]["family"] or None)) for n in seqs}
    dread = {n for n in seqs if lab[n]["role"] == "D"}
    new = CH.refine(cl, seqs, CH.minimap_allvsall(f"{d}/{dst}_tmp"), DELTA)
    os.makedirs(f"{d}/{dst}", exist_ok=True)
    maj = {c: S.majority(rs, truth) for c, rs in new.items()}
    with open(f"{d}/{dst}/clusters.tsv", "w") as o:
        o.write("read\tcluster\tsize\tmajority\n")
        for c, rs in new.items():
            for r in rs:
                o.write(f"{r}\t{c}\t{len(rs)}\t{maj[c] or ''}\n")
    m = S.cluster_metrics(new, truth, dread)
    m["background_in_family_clusters"] = sum(1 for c, rs in new.items() if maj[c] for r in rs if truth.get(r) == "bg")
    m["delta"] = DELTA
    json.dump(m, open(f"{d}/{dst}/metrics.json", "w"), indent=1)
    print(f"{bed}: {len(cl)} clusters -> {len(new)}; D reads clustered {m['d_clustered']} of {m['n_d']} ({m['coverage']:.1%}), purity {m['purity']:.4f}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
