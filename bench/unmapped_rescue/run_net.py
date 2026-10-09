#!/usr/bin/env python3
"""Amendment 41: the adopted pipeline on the whole net of KB3781 (526,772 reads), blind to the haplotype assemblies, resumable at every stage. Miniforge python.

    run_net.py            one call = as much as fits in TIME_BUDGET; exit 75 = run again; prints DONE at the end
Stages: rounds (seed-round clustering, proper edge) -> refine (star step, per component) -> consensus (abPOA local mode) -> align (consensus on the primary) -> classify
(identity x coverage with terminal-exon rescue and the library gate; support test). Output OUT/net_run/{clusters.tsv, cons.fa, classes.json}."""
import collections
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import chain as CH  # noqa: E402
import consensus_support as CS  # noqa: E402
import discover as D  # noqa: E402
import flagmetric as FM  # noqa: E402
import graph as GR  # noqa: E402
import partition as PT  # noqa: E402
import seeds as SD  # noqa: E402

NET = "/mnt/linuxdisk/tmp/o3_rescue/mattruth/net_all.fa"
R = "/mnt/linuxdisk/tmp/o3_rescue/mattruth/net_run"
BUDGET = int(os.environ.get("RUN_NET_BUDGET", "500"))
T0 = time.time()


def left():
    return BUDGET - (time.time() - T0)


def pause(stage):
    print(f"{stage}: paused after {time.time() - T0:.0f} s; run again")
    sys.exit(75)


def jl_load(p):
    out = {}
    if os.path.exists(p):
        for ln in open(p):
            if ln.strip():
                k, v = json.loads(ln)
                out[k] = v
    return out


def main():
    os.makedirs(f"{R}/rounds", exist_ok=True)
    seqs = SD.read_fa(NET)
    lens = {n: len(s) for n, s in seqs.items()}
    # 1. rounds (the query FASTA of a finished round is not needed again: removed to save disk)
    for i in range(1, 7):
        if os.path.exists(f"{R}/rounds/round{i}.paf.done") and os.path.exists(f"{R}/rounds/query{i}.fa"):
            os.remove(f"{R}/rounds/query{i}.fa")
    if not os.path.exists(f"{R}/components.json"):
        try:
            comp = SD.run_rounds(lens, SD.minimap_map_fn(seqs, f"{R}/rounds", D.DELTA, 0.5, proper=True), n_seeds=2000, min_size=3, max_rounds=6)
        except SD.Pause as e:
            pause(f"round {e}")
        cl = GR.clusters(comp, 3)
        json.dump({str(k): v for k, v in cl.items()}, open(f"{R}/components.json", "w"))
        print(f"rounds done: {len(cl)} components of >= 3 reads, {sum(len(v) for v in cl.values())} reads")
    comps = json.load(open(f"{R}/components.json"))
    # 2. refine, per component (resumable)
    ref_p = f"{R}/refine.jsonl"
    done = jl_load(ref_p)
    if len(done) < len(comps):
        allvsall = CH.minimap_allvsall(f"{R}/chain_tmp")
        with open(ref_p, "a") as o:
            for cid in sorted(comps, key=lambda c: (len(comps[c]), c)):
                if cid in done:
                    continue
                if left() < 30:
                    pause(f"refine ({len(done)} of {len(comps)} components)")
                new = CH.refine({cid: comps[cid]}, seqs, allvsall, D.DELTA)
                o.write(json.dumps([cid, {str(k): v for k, v in new.items()}]) + "\n")
                done[cid] = new
        print("refine done")
    clusters = {}
    for cid, new in done.items():
        for k, rs in new.items():
            clusters[k] = rs
    if not os.path.exists(f"{R}/clusters.tsv"):
        with open(f"{R}/clusters.tsv", "w") as o:
            o.write("read\tcluster\tsize\n")
            for k, rs in clusters.items():
                for r in rs:
                    o.write(f"{r}\t{k}\t{len(rs)}\n")
        print(f"{len(clusters)} clusters, {sum(len(v) for v in clusters.values())} reads")
    finish(R, clusters, seqs)


def finish(R, clusters, seqs):
    """stages 3-5 on a clustering {key: [reads]}: local-mode consensus, alignment on the primary, classification (resumable; R = run directory)"""
    # 3. consensus, local mode (resumable)
    cons_p = f"{R}/cons.jsonl"
    cons = jl_load(cons_p)
    if len(cons) < len(clusters):
        with open(cons_p, "a") as o:
            for k in sorted(clusters):
                if k in cons:
                    continue
                if left() < 20:
                    pause(f"consensus ({len(cons)} of {len(clusters)})")
                c = PT.abpoa_consensus([seqs[r] for r in clusters[k]], mode=D.CONS_MODE)
                o.write(json.dumps([k, c]) + "\n")
                cons[k] = c
        print("consensus done")
    keys = sorted(clusters)
    name = {k: f"cl{i}" for i, k in enumerate(keys)}
    if not os.path.exists(f"{R}/cons.fa"):
        SD.write_fa(f"{R}/cons.fa", {name[k]: cons[k] for k in keys}, [name[k] for k in keys])
    # 4. align the consensus sequences on the primary (one call)
    if not os.path.exists(f"{R}/cons.R.paf.done"):
        if left() < 240:
            pause("align (needs a fresh call)")
        subprocess.run(f"minimap2 -c -x splice:hq -uf -N 5 -t 4 {D.RP.GENOME_IDX} {R}/cons.fa > {R}/cons.R.paf", shell=True, check=True, stderr=subprocess.DEVNULL)
        open(f"{R}/cons.R.paf.done", "w").write("ok")
        print("align done")
    # 5. classify (resumable): score on R with rescue + library gate, support test
    import pysam
    gate = json.load(open(f"{D.W}/o3hap_mat/gate.json"))["gate"]
    recs = D.best_records(f"{R}/cons.R.paf")
    fa = pysam.FastaFile(D.PRIMARY_FA)
    cls_p = f"{R}/classes.jsonl"
    rows = jl_load(cls_p)
    if len(rows) < len(keys):
        with open(cls_p, "a") as o:
            for k in keys:
                n = name[k]
                if n in rows:
                    continue
                if left() < 10:
                    pause(f"classify ({len(rows)} of {len(keys)})")
                sc = D.idcov_with_rescue(cons[k], recs[n], fa, D.ident, gate, f"{R}/tmp_R") if n in recs else None
                med = CS.median_read_divergence(cons[k], {r: seqs[r] for r in sorted(clusters[k])[:100]}, f"{R}/tmp_support")
                cls = "UNSUPPORTED" if med is None or med > D.DELTA else FM.elsewhere_class(sc, [])
                hit = recs.get(n)
                row = dict(k=n, key=k, reads=len(clusters[k]), length=len(cons[k]), R=sc, median_read_divergence=med, cls=cls,
                           R_hit=(hit[5], int(hit[7]), int(hit[8])) if hit else None)
                o.write(json.dumps([n, row]) + "\n")
                rows[n] = row
    out = [rows[name[k]] for k in keys]
    json.dump(out, open(f"{R}/classes.json", "w"))
    by = collections.Counter(r["cls"] for r in out)
    rd = collections.Counter()
    for r in out:
        rd[r["cls"]] += r["reads"]
    print("classes (clusters / reads):", {c: (by[c], rd[c]) for c in by})
    print("DONE")


if __name__ == "__main__":
    main()
