#!/usr/bin/env python3
"""Amendment 42: the pipeline on the whole net with locus-binned clustering instead of the seed rounds; blind to the haplotype assemblies. Miniforge python.

    run_binned.py        one call = as much as fits in the budget; exit 75 = run again; prints DONE at the end. Output OUT/net_run_binned/."""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import binclust as B  # noqa: E402
import run_net as RN  # noqa: E402
import seeds as SD  # noqa: E402

OUT = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"
R = f"{OUT}/net_run_binned"


def main():
    os.makedirs(R, exist_ok=True)
    seqs = SD.read_fa(RN.NET)
    lens = {n: len(s) for n, s in seqs.items()}
    # 1. bins: each read's own primary alignment; the unmapped reads are one more bin
    if not os.path.exists(f"{R}/bins.json"):
        recs = []
        for fn in sorted(glob.glob(f"{OUT}/netpos/*.tsv")):
            for ln in open(fn):
                n, c, s, e, st = ln.rstrip("\n").split("\t")
                recs.append((n, c, int(s), int(e), st))
        placed = {r[0] for r in recs}
        bins = B.bin_reads(recs) + [sorted(n for n in seqs if n not in placed)]
        bins = [b for b in bins if len(b) >= 3]
        json.dump(bins, open(f"{R}/bins.json", "w"))
        print(f"bins: {len(bins)} with >= 3 reads, {sum(len(b) for b in bins)} reads; unmapped bin {len(seqs) - len(placed)} reads; largest {sorted((len(b) for b in bins), reverse=True)[:5]}")
    bins = json.load(open(f"{R}/bins.json"))
    # 2. cluster every bin (resumable, smallest bins first)
    done = RN.jl_load(f"{R}/bins.jsonl")
    if len(done) < len(bins):
        cf, af = B.compat_fn(seqs, f"{R}/tmp_ava", threads=4), B.assign_fn(seqs, f"{R}/tmp_assign", threads=4)
        with open(f"{R}/bins.jsonl", "a") as o:
            for i in sorted(range(len(bins)), key=lambda j: (len(bins[j]), j)):
                if str(i) in done:
                    continue
                if RN.left() < 30:
                    RN.pause(f"binned clustering ({len(done)} of {len(bins)} bins)")
                sp = f"{R}/state_{i}.json"               # a large bin is checkpointed pass by pass
                state = json.load(open(sp)) if os.path.exists(sp) else {}
                cl = B.cluster_bin(bins[i], lens, cf, af, cap=int(os.environ.get("BIN_CAP", "500")), state=state, stop=lambda: RN.left() < 150)
                if cl is None:
                    json.dump(state, open(sp, "w"))
                    RN.pause(f"bin {i} ({len(bins[i])} reads) after pass {state['passes']}; {len(done)} of {len(bins)} bins")
                o.write(json.dumps([str(i), cl]) + "\n")
                o.flush()
                done[str(i)] = cl
                if os.path.exists(sp):
                    os.remove(sp)
        print("binned clustering done")
    clusters = {}
    for i, cl in done.items():
        for k, rs in cl.items():
            clusters[k] = rs
    if not os.path.exists(f"{R}/clusters.tsv"):
        with open(f"{R}/clusters.tsv", "w") as o:
            o.write("read\tcluster\tsize\n")
            for k, rs in clusters.items():
                for r in rs:
                    o.write(f"{r}\t{k}\t{len(rs)}\n")
        print(f"{len(clusters)} clusters, {sum(len(v) for v in clusters.values())} reads")
    RN.finish(R, clusters, seqs)


if __name__ == "__main__":
    main()
