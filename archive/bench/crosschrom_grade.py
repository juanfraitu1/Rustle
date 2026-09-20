#!/usr/bin/env python3
"""Grade the cross-chromosome confirmed pairs by POA core-block IDENTITY (the discriminator that
separates real homology from chance end-to-end alignment, which core_cov alone cannot). Recomputes
poa_pair_stats on each cross-chrom pair (parallel) and adds core_ident; characterizes the regime.

Run with /home/juanfra/miniforge3/bin/python.
"""
import csv
import multiprocessing as mp
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import poa_family_definition as P

FA = "/tmp/gene_reps_gw.fa"
IN = os.path.join(os.path.dirname(__file__), "crosschrom_pairs.tsv")
OUT = os.path.join(os.path.dirname(__file__), "crosschrom_graded.tsv")

_S = None
_A = None
def _init():
    global _S, _A
    _S = P.load_fasta(FA)
    _A = P.make_aligner()
def _work(t):
    a, b = t
    sa, sb = _S.get(a), _S.get(b)
    if not sa or not sb:
        return None
    st = P.poa_pair_stats(_A, sa, sb)
    return (a, b, round(st["core_recip"], 3), round(st["core_ident"], 3), round(st["ident"], 3))


def main():
    rows = list(csv.DictReader(open(IN), delimiter="\t"))
    pairs = [(r["gene_a"], r["gene_b"]) for r in rows]
    meta = {(r["gene_a"], r["gene_b"]): r for r in rows}
    with mp.Pool(5, initializer=_init) as pool:
        res = [x for x in pool.map(_work, pairs, chunksize=64) if x]

    with open(OUT, "w") as fh:
        fh.write("gene_a\tchrom_a\tgene_b\tchrom_b\tcore_cov\tcore_ident\taln_ident\tjaccard\n")
        for a, b, cc, ci, ai in res:
            r = meta[(a, b)]
            fh.write(f"{a}\t{r['chrom_a']}\t{b}\t{r['chrom_b']}\t{cc}\t{ci}\t{ai}\t{r['jaccard']}\n")

    ci = [x[3] for x in res]
    print(f"graded {len(res)} cross-chrom pairs -> {OUT}")
    # distribution of core_ident
    import numpy as np
    arr = np.array(ci)
    for lo in (0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3):
        print(f"  core_ident >= {lo}: {(arr>=lo).sum()}")
    print(f"  core_ident < 0.3 (chance baseline ~0.25): {(arr<0.3).sum()}")
    # anchor: where does RABL2 land?
    print("=== RABL2 anchor pairs (core_ident) ===")
    for a, b, cc, cid, ai in res:
        if "RABL2" in a or "RABL2" in b:
            print(f"    {a} <-> {b}: core_cov={cc} core_ident={cid} aln_ident={ai}")
    # the chance suspects (high core, low jaccard) — their core_ident?
    print("=== chance suspects (housekeeping<->LOC, was core=1.0 jacc<0.1) ===")
    for a, b, cc, cid, ai in res:
        if a in ("ATP5MF", "ATP5PD", "EEF1A1", "GNG5B") or b in ("ATP5MF", "ATP5PD", "EEF1A1", "GNG5B"):
            print(f"    {a} <-> {b}: core_cov={cc} core_ident={cid} aln_ident={ai}")


if __name__ == "__main__":
    main()
