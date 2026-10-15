#!/usr/bin/env python3
"""Own-read evidence of a candidate copy (spec section 1, step 4).

A read is the candidate's own when its PRIMARY alignment (`-F 2308`: not unmapped, secondary or supplementary) has an
aligned base (M / = / X) on the candidate's projected exons: figure 8's read-support rule (C2). Its exon blocks are its
alignment split at `N` (a deletion stays inside an exon), 0-based half-open.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "figures"))
import _liftoff as LF  # noqa: E402
from templates import merge  # noqa: E402

EXCLUDE = 0x4 | 0x100 | 0x800  # 2308
M_OPS, D_OP, N_OP = (0, 7, 8), 2, 3
FETCH_SPAN = 200_000


def blocks(reference_start, cigartuples):
    """-> (aligned blocks: M/=/X, merged; exon blocks: split at N, D inside)."""
    pos = cur = reference_start
    aligned, exons = [], []
    for op, n in cigartuples:
        if op in M_OPS:
            aligned.append((pos, pos + n))
            pos += n
        elif op == D_OP:
            pos += n
        elif op == N_OP:
            if pos > cur:
                exons.append((cur, pos))
            pos += n
            cur = pos
    if pos > cur:
        exons.append((cur, pos))
    return merge(aligned), exons


def own_reads(bam, contig, exons):
    """The candidate's own reads in `bam` (pysam.AlignmentFile) -> [(read name, exon blocks)], sorted by name.
    A span over FETCH_SPAN bp is fetched exon by exon (a read with an aligned base on the exons overlaps one of them,
    so the reads are the same) instead of walking the reads of a long projected intron."""
    out = {}
    windows = [(exons[0][0], exons[-1][1])] if exons[-1][1] - exons[0][0] <= FETCH_SPAN else exons
    for s, e in windows:
        for r in bam.fetch(contig, s, e):
            if r.query_name in out or r.flag & EXCLUDE or r.cigartuples is None:
                continue
            al, ex = blocks(r.reference_start, r.cigartuples)
            if LF.iv_inter(al, exons) > 0:
                out[r.query_name] = ex
    return sorted(out.items())
