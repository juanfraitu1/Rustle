#!/usr/bin/env python3
"""Amendment 38: trim a consensus to the columns covered by >= 2 of its reads. Needs minimap2; miniforge python.

support_window(length, intervals)      pure: (start, end) of the columns with coverage >= 2, or (0, length) when there are none
ref_span(pos0, cigar)                  pure: reference interval of an alignment
trim_consensus(cons, reads, tmp)       aligns the reads to the consensus (splice:hq -uf, as the support test) and trims"""
import os
import re
import subprocess

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
MIN_COVER = 2


def ref_span(pos0, cigar):
    """[start, end) on the reference of an alignment starting at 0-based pos0 (M, =, X, D and N consume the reference)"""
    n = sum(int(k) for k, op in CIG.findall(cigar) if op in "M=XDN")
    return pos0, pos0 + n


def support_window(length, intervals, min_cover=MIN_COVER):
    """first and last column (end exclusive) covered by >= min_cover intervals; (0, length) when no column is, so the consensus is kept"""
    diff = [0] * (length + 1)
    for a, b in intervals:
        a, b = max(0, a), min(length, b)
        if a < b:
            diff[a] += 1
            diff[b] -= 1
    cov, first, last = 0, None, None
    for i in range(length):
        cov += diff[i]
        if cov >= min_cover:
            if first is None:
                first = i
            last = i
    return (0, length) if first is None else (first, last + 1)


def trim_consensus(cons, reads, tmp, min_cover=MIN_COVER):
    """reads: list of sequences. -> the consensus cut to its support window"""
    os.makedirs(tmp, exist_ok=True)
    open(f"{tmp}/r.fa", "w").write("".join(f">r{i}\n{s}\n" for i, s in enumerate(reads)))
    open(f"{tmp}/c.fa", "w").write(f">c\n{cons}\n")
    sam = subprocess.run(f"minimap2 -ax splice:hq -uf -t 2 {tmp}/c.fa {tmp}/r.fa", shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True).stdout.splitlines()
    iv = []
    for ln in sam:
        if ln[0] == "@":
            continue
        f = ln.split("\t")
        if int(f[1]) & 2308 or f[2] == "*":
            continue
        iv.append(ref_span(int(f[3]) - 1, f[5]))
    a, b = support_window(len(cons), iv, min_cover)
    return cons[a:b]
