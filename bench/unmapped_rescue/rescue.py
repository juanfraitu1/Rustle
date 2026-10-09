#!/usr/bin/env python3
"""Terminal-exon rescue for the discovery coverage (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 29).

A spliced aligner can leave a small terminal exon (46 or 90 bases) of a consensus unaligned although it lies in the same locus. end_segments finds the unaligned ends of
the best record; the pieces are mapped (short-read preset) to a window around the record; rescued turns the good pieces into consensus intervals; combine adds them to the
main alignment without counting a base twice."""
import subprocess


def union_length(intervals):
    tot, cur = 0, None
    for a, b in sorted(intervals):
        if cur is None or a > cur[1]:
            if cur is not None:
                tot += cur[1] - cur[0]
            cur = [a, b]
        else:
            cur[1] = max(cur[1], b)
    return tot + (cur[1] - cur[0] if cur else 0)


def end_segments(qs, qe, qlen, min_len=20):
    """unaligned ends of the best record, (start, end) in consensus coordinates, each at least min_len long"""
    out = []
    if qs >= min_len:
        out.append((0, qs))
    if qlen - qe >= min_len:
        out.append((qe, qlen))
    return out


def rescued(segments, paf_lines, min_ident=0.98):
    """paf_lines: records of the segment sequences (named seg0, seg1, ...) mapped to the window. -> (consensus intervals, matches, block) of the records at >= min_ident"""
    iv, m, b = [], 0, 0
    for ln in paf_lines:
        f = ln.rstrip("\n").split("\t")
        if not f[0].startswith("seg"):
            continue
        s0 = segments[int(f[0][3:])][0]
        matches, blk = int(f[9]), int(f[10])
        if blk and matches / blk >= min_ident:
            iv.append((s0 + int(f[2]), s0 + int(f[3])))
            m += matches
            b += blk
    return iv, m, b


def combine(qlen, main, pieces):
    """main = (qs, qe, matches, block); pieces = (intervals, matches, block) -> (aligned bases, identity over everything aligned)"""
    qs, qe, mm, mb = main
    iv, pm, pb = pieces
    total = union_length([(qs, qe)] + list(iv))
    return total, (mm + pm) / max(1, mb + pb)


def map_pieces(window_seq, segments, cons, workdir):
    """minimap2 -x sr of the unaligned segments against the window sequence"""
    import os
    os.makedirs(workdir, exist_ok=True)
    wf, qf = f"{workdir}/window.fa", f"{workdir}/segs.fa"
    open(wf, "w").write(f">w\n{window_seq}\n")
    with open(qf, "w") as o:
        for i, (a, b) in enumerate(segments):
            o.write(f">seg{i}\n{cons[a:b]}\n")
    return subprocess.run(f"minimap2 -c -x sr -N 5 -t 2 {wf} {qf}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout.splitlines()
