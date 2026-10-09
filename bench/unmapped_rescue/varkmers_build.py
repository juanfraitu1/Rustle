#!/usr/bin/env python3
"""Amendment 43b: variant k-mers of every flagged consensus of the binned run -> query set for kmerhit. Miniforge python.
Writes /mnt/linuxdisk/tmp/dnaverify2/{qset.u64, var.npz}."""
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import discover as D  # noqa: E402
import run_polish as RP  # noqa: E402
import varkmers as V  # noqa: E402
from o3_maternal import common as C  # noqa: E402

O = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"
OUT = "/mnt/linuxdisk/tmp/dnaverify2"
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
FLANK = 1000


def best(paf):
    b = {}
    for ln in open(paf):
        f = ln.rstrip("\n").split("\t")
        if "tp:A:P" in f[12:] and (f[0] not in b or int(f[9]) > int(b[f[0]][9])):
            b[f[0]] = f
    return b


def blocks(ts, cigar):
    """target intervals of the aligned blocks between introns"""
    out, pos, start = [], ts, ts
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "M=XD":
            pos += n
        elif op == "N":
            out.append((start, pos))
            pos += n
            start = pos
    out.append((start, pos))
    return out


def main():
    import pysam
    os.makedirs(OUT, exist_ok=True)
    ev = json.load(open(f"{O}/eval_binned.json"))["verdicts"]
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{O}/net_run_binned/cons.fa").items()}
    al = C.alias()
    recs = {"pri": best(f"{O}/spec_flags.cs.paf"), "mat": best(f"{O}/flag_binned_cons.mat.paf"), "pat": best(f"{O}/flag_binned_cons.pat.paf")}
    fas = {"pri": (pysam.FastaFile(D.PRIMARY_FA), lambda n: n), "mat": (pysam.FastaFile(D.HAP_FA.format("mat")), lambda n: C.accession(n, al)),
           "pat": (pysam.FastaFile(D.HAP_FA.format("pat")), lambda n: C.accession(n, al))}
    var_all, var_pri, keys = {}, {}, []
    for v in ev:
        k = v["k"]
        pos, codes = V.kmer_codes(cons[k])
        jn = V.query_junctions(next(t[5:] for t in recs["pri"][k][12:] if t.startswith("cg:Z:"))) if k in recs["pri"] else []
        codes = codes[V.not_across(pos, jn)]
        ref = {}
        for w in ("pri", "mat", "pat"):
            f = recs[w].get(k)
            if not f:
                ref[w] = np.zeros(0, dtype=np.uint64)
                continue
            fa, nm = fas[w]
            c = nm(f[5])
            L = fa.get_reference_length(c)
            cg = next(t[5:] for t in f[12:] if t.startswith("cg:Z:"))
            parts = [V.kmer_codes(fa.fetch(c, max(0, a - FLANK), min(L, b + FLANK)).upper())[1] for a, b in blocks(int(f[7]), cg)]
            ref[w] = np.unique(np.concatenate(parts)) if parts else np.zeros(0, dtype=np.uint64)
        u = np.unique(codes)
        var_pri[k] = np.setdiff1d(u, ref["pri"])
        var_all[k] = np.setdiff1d(var_pri[k], np.union1d(ref["mat"], ref["pat"]))
        keys.append(k)
    q = np.unique(np.concatenate([var_pri[k] for k in keys]))       # var_all is a subset of var_pri
    q.astype("<u8").tofile(f"{OUT}/qset.u64")
    np.savez(f"{OUT}/var.npz", keys=np.array(keys), **{f"all_{i}": np.searchsorted(q, var_all[k]) for i, k in enumerate(keys)},
             **{f"pri_{i}": np.searchsorted(q, var_pri[k]) for i, k in enumerate(keys)})
    print(f"{len(keys)} consensus sequences; query set {len(q)} k-mers; median var_pri {int(np.median([len(var_pri[k]) for k in keys]))}, var_all {int(np.median([len(var_all[k]) for k in keys]))}")


if __name__ == "__main__":
    main()
