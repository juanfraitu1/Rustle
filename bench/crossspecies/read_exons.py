#!/usr/bin/env python3
"""SAM on stdin (already -F 2308 filtered) -> BED of per-read aligned blocks (exons).

N in an RNA CIGAR is an intron spliced OUT and BREAKS the block. D does not.
SAM POS is 1-based; BED is 0-based half-open.
"""
import sys

REF_CONSUMING_IN_BLOCK = {"M", "=", "X", "D"}

for line in sys.stdin:
    if line.startswith("@"):
        continue
    f = line.split("\t", 6)
    if len(f) < 6:
        continue
    rname, pos, cigar = f[2], int(f[3]) - 1, f[5]
    if cigar == "*":
        continue
    ref = pos
    bs = ref
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        n = int(num or 0)
        num = ""
        if ch in REF_CONSUMING_IN_BLOCK:
            ref += n
        elif ch == "N":
            if ref > bs:
                sys.stdout.write(f"{rname}\t{bs}\t{ref}\t1\n")
            ref += n
            bs = ref
        # S/H/I/P consume no reference
    if ref > bs:
        sys.stdout.write(f"{rname}\t{bs}\t{ref}\t1\n")
