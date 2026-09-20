#!/usr/bin/env python3
"""Complete annotation-based RNA-level family definition: grade ALL candidate gene pairs (same- AND
cross-chromosome) by POA contiguous-core coverage + core-block identity, so families = connected
components over the full homology graph (not just the cross-chrom subset graded earlier).

Reads /tmp/cc_candidates.tsv (minimizer-LSH candidate pairs over all 22,983 annotated gene reps).
Run with /home/juanfra/miniforge3/bin/python (Bio.Align). Parallel + deterministic.
"""
import csv
import multiprocessing as mp
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import poa_family_definition as P

FA = "/tmp/gene_reps_gw.fa"
META = "/tmp/gene_reps_gw.meta.tsv"
CAND = "/tmp/cc_candidates.tsv"
OUT = os.path.join(os.path.dirname(__file__), "family_pairs_graded.tsv")
T_CORE = 0.13
LEN_CAP = 20000

_S = None
_A = None
def _init():
    global _S, _A
    _S = P.load_fasta(FA); _A = P.make_aligner()
def _work(t):
    a, b = t
    sa, sb = _S.get(a), _S.get(b)
    if not sa or not sb or min(len(sa), len(sb)) > LEN_CAP:
        return None
    st = P.poa_pair_stats(_A, sa, sb)
    if st["core_recip"] >= T_CORE:
        return (a, b, round(st["core_recip"], 3), round(st["core_ident"], 3))
    return None


def main():
    g2chrom = {}
    with open(META) as fh:
        fh.readline()
        for line in fh:
            c = line.rstrip("\n").split("\t"); g2chrom[c[0]] = c[1]
    pairs = []
    with open(CAND) as fh:
        fh.readline()
        for line in fh:
            a, b, s, j = line.rstrip("\n").split("\t"); pairs.append((a, b))
    print(f"grading {len(pairs)} candidate pairs (all same+cross chrom)...")
    with mp.Pool(5, initializer=_init) as pool:
        res = [r for r in pool.map(_work, pairs, chunksize=64) if r]
    with open(OUT, "w") as fh:
        fh.write("gene_a\tchrom_a\tgene_b\tchrom_b\tcore_cov\tcore_ident\tsame_chrom\n")
        for a, b, cc, ci in sorted(res, key=lambda r: -r[3]):
            ca, cb = g2chrom.get(a, "?"), g2chrom.get(b, "?")
            fh.write(f"{a}\t{ca}\t{b}\t{cb}\t{cc}\t{ci}\t{int(ca == cb)}\n")
    same = sum(1 for a, b, cc, ci in res if g2chrom.get(a) == g2chrom.get(b))
    print(f"confirmed homologous pairs (core_cov>={T_CORE}): {len(res)} "
          f"(same-chrom={same}, cross-chrom={len(res)-same}) -> {OUT}")


if __name__ == "__main__":
    main()
