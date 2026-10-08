#!/usr/bin/env python3
"""abPOA consensus of a cluster (docs/PREREG_unmapped_rescue_2026-10-08.md section 1, step 2). Needs pyabpoa: run with /home/juanfra/miniforge3/bin/python.

    consensus.py <pool.fa> <clusters.tsv> <out.fa>    # one consensus per cluster, <= 100 reads, longest first
"""
import csv
import sys


def pick_reads(reads, lens, k=100):
    """at most k reads, longest first, ties by name"""
    return sorted(reads, key=lambda r: (-lens[r], r))[:k]


def main(pool, clusters, out):
    import pyabpoa
    seqs, cur = {}, None
    for ln in open(pool):
        if ln[0] == ">":
            cur = ln[1:].strip().split()[0]
            seqs[cur] = []
        else:
            seqs[cur].append(ln.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    lens = {k: len(v) for k, v in seqs.items()}
    cl = {}
    for r in csv.DictReader(open(clusters), delimiter="\t"):
        cl.setdefault(r["cluster"], []).append(r["read"])
    aln = pyabpoa.msa_aligner(aln_mode="g", is_aa=False, cons_algrm="HB")
    with open(out, "w") as o:
        for cid, reads in sorted(cl.items(), key=lambda kv: -len(kv[1])):
            pick = pick_reads(reads, lens, 100)
            res = aln.msa([seqs[r] for r in pick], out_cons=True, out_msa=False)
            if res.cons_seq:
                o.write(f">cl{cid}|n={len(reads)}\n{res.cons_seq[0]}\n")
    print("consensus written for", len(cl), "clusters")


if __name__ == "__main__":
    main(*sys.argv[1:4])
