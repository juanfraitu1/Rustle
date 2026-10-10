#!/usr/bin/env python3
"""Amendment 44: IsoCon-style read support of a flagged candidate among all reads of its locus, against the primary's own sequence there.

Pure (tested in test_locus_support.py): assign, share, minor, sample_every. Driver: main() over every flagged cluster of the binned run (resumable)."""
import json
import math
import os
import subprocess
import sys
import time


def assign(paf_lines, cons_name, ref_name):
    """-> (reads whose primary record is on the consensus, reads whose primary record is on the reference)"""
    best = {}
    for ln in paf_lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 12 or "tp:A:P" not in f[12:] or f[5] not in (cons_name, ref_name):
            continue
        if f[0] not in best or int(f[9]) > best[f[0]][1]:
            best[f[0]] = (f[5], int(f[9]))
    n_c = sum(t == cons_name for t, _ in best.values())
    return n_c, sum(t == ref_name for t, _ in best.values())


def share(n_cons, n_ref):
    n = n_cons + n_ref
    return n_cons / n if n else None


def minor(k, n, alpha=0.01):
    """one-sided binomial: are k of n reads too few for an allele carrying half the locus (P(X <= k | n, 1/2) < alpha)?"""
    if n == 0:
        return False
    p = sum(math.comb(n, i) for i in range(k + 1)) / 2 ** n
    return p < alpha


def sample_every(names, cap):
    if cap <= 0:
        return []
    if len(names) <= cap:
        return list(names)
    step = len(names) / cap
    return [names[int(i * step)] for i in range(cap)]


O = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"
BAM = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
RCT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def main(budget=520, cap=300):
    import pysam
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import run_polish as RP
    import seeds as SD
    t0 = time.time()
    rows = {r["k"]: r for r in json.load(open(f"{O}/net_run_binned/classes.json"))}
    ev = json.load(open(f"{O}/eval_binned.json"))["verdicts"]
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{O}/net_run_binned/cons.fa").items()}
    pv = {}
    for side in ("dev", "held"):
        if os.path.exists(f"{O}/privers.{side}.fa"):
            pv.update(SD.read_fa(f"{O}/privers.{side}.fa"))
    outp = f"{O}/locus_support.jsonl"
    done = {json.loads(l)["k"] for l in open(outp)} if os.path.exists(outp) else set()
    bam = pysam.AlignmentFile(BAM, "rb")
    tmp = f"{O}/tmp_locus"
    os.makedirs(tmp, exist_ok=True)
    with open(outp, "a") as o:
        for v in ev:
            k = v["k"]
            if k in done:
                continue
            if time.time() - t0 > budget:
                print(f"paused: {len(done)} of {len(ev)}")
                sys.exit(75)
            hit = rows[k]["R_hit"]
            rec = dict(k=k, n_cons=None, n_ref=None, locus_reads=0)
            if hit and k in pv and pv[k]:
                c, s, e = hit
                names, seqs = [], {}
                for a in bam.fetch(c, s, e):
                    if a.is_secondary or a.is_supplementary or a.is_unmapped or a.query_sequence is None:
                        continue
                    q = a.query_sequence
                    seqs[a.query_name] = q.translate(RCT)[::-1] if a.is_reverse else q
                    names.append(a.query_name)
                pick = sample_every(names, cap)
                SD.write_fa(f"{tmp}/r.fa", seqs, pick)
                SD.write_fa(f"{tmp}/t.fa", {"cons": cons[k], "ref": pv[k]}, ["cons", "ref"])
                lines = subprocess.run(f"minimap2 -c -x splice:hq -uf --secondary=no -t 4 {tmp}/t.fa {tmp}/r.fa", shell=True, stdout=subprocess.PIPE,
                                       stderr=subprocess.DEVNULL, text=True).stdout.splitlines()
                n_c, n_r = assign(lines, "cons", "ref")
                rec.update(n_cons=n_c, n_ref=n_r, locus_reads=len(names), sampled=len(pick))
            o.write(json.dumps(rec) + "\n")
            o.flush()
            done.add(k)
    print("DONE")


if __name__ == "__main__":
    main()
