#!/usr/bin/env python3
"""Amendment 45: IsoCon-style read correction inside the Amendment 42 bins, then the same clustering, consensus from corrected reads, support test with the
original reads. Resumable; exit 75 = run again; DONE at the end. Output OUT/net_run_corr/."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import binclust as B  # noqa: E402
import readcorr as RC  # noqa: E402
import run_net as RN  # noqa: E402
import seeds as SD  # noqa: E402

OUT = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"
R = f"{OUT}/net_run_corr"
CHUNK = 300


def main():
    os.makedirs(f"{R}/tmp", exist_ok=True)
    orig = SD.read_fa(RN.NET)
    bins = json.load(open(f"{OUT}/net_run_binned/bins.json"))
    order = sorted(range(len(bins)), key=lambda j: (len(bins[j]), j))
    big = set(order[-132:])                                   # Amendment 42b: the 132 largest bins ran at cap 200
    # 1. correction, bin by bin (smallest first)
    dp = f"{R}/corr_done.jsonl"
    done = {json.loads(l) for l in open(dp)} if os.path.exists(dp) else set()
    if len(done) < len(bins):
        with open(f"{R}/corrected.fa", "a") as fa, open(dp, "a") as o:
            for i in order:
                if i in done:
                    continue
                if RN.left() < 40:
                    RN.pause(f"correction ({len(done)} of {len(bins)} bins)")
                names = sorted(bins[i], key=lambda r: (-len(orig[r]), r))
                cp_ = f"{R}/corr_chunks.jsonl"                 # a large bin is checkpointed chunk by chunk
                cdone_ = {tuple(json.loads(l)) for l in open(cp_)} if os.path.exists(cp_) else set()
                for c in range(0, len(names), CHUNK):
                    if (i, c) in cdone_:
                        continue
                    if RN.left() < 40:
                        RN.pause(f"correction: bin {i} ({len(bins[i])} reads) at chunk {c // CHUNK}; {len(done)} of {len(bins)} bins")
                    part = {n: orig[n] for n in names[c:c + CHUNK]}
                    for n, s in RC.correct_set(part, f"{R}/tmp").items():
                        fa.write(f">{n}\n{s}\n")
                    fa.flush()
                    with open(cp_, "a") as cpo:
                        cpo.write(json.dumps([i, c]) + "\n")
                o.write(json.dumps(i) + "\n")
                o.flush()
                done.add(i)
        print("correction done")
    corr = dict(orig)
    corr.update(SD.read_fa(f"{R}/corrected.fa"))              # the last copy of a read wins (a resumed bin is written again)
    lens = {n: len(s) for n, s in corr.items()}
    # 2. clustering on the corrected reads (as Amendment 42 ran: cap 500, 200 for the 132 largest bins)
    cp = f"{R}/bins.jsonl"
    cdone = RN.jl_load(cp)
    if len(cdone) < len(bins):
        cf, af = B.compat_fn(corr, f"{R}/tmp_ava", threads=4), B.assign_fn(corr, f"{R}/tmp_assign", threads=4)
        with open(cp, "a") as o:
            for i in order:
                if str(i) in cdone:
                    continue
                if RN.left() < 30:
                    RN.pause(f"clustering ({len(cdone)} of {len(bins)} bins)")
                sp = f"{R}/state_{i}.json"
                state = json.load(open(sp)) if os.path.exists(sp) else {}
                cl = B.cluster_bin(bins[i], lens, cf, af, cap=200 if i in big else 500, state=state, stop=lambda: RN.left() < 150)
                if cl is None:
                    json.dump(state, open(sp, "w"))
                    RN.pause(f"bin {i} ({len(bins[i])} reads) after pass {state['passes']}")
                o.write(json.dumps([str(i), cl]) + "\n")
                o.flush()
                cdone[str(i)] = cl
                if os.path.exists(sp):
                    os.remove(sp)
        print("clustering done")
    clusters = {k: rs for cl in cdone.values() for k, rs in cl.items()}
    if not os.path.exists(f"{R}/clusters.tsv"):
        with open(f"{R}/clusters.tsv", "w") as o:
            o.write("read\tcluster\tsize\n")
            for k, rs in clusters.items():
                for r in rs:
                    o.write(f"{r}\t{k}\t{len(rs)}\n")
        print(f"{len(clusters)} clusters, {sum(len(v) for v in clusters.values())} reads")
    RN.finish(R, clusters, corr, support_seqs=orig)


if __name__ == "__main__":
    main()
