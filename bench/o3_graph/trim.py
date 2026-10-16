#!/usr/bin/env python3
"""Library adapters off read ends before mapping (2026-10-16: the public testis library keeps its Tn5 mosaic-end / SMART
adapters at both ends of 99% of reads; Iso-Seq's lima/refine steps would remove them, and a read-to-path O2 cannot work
around them).

    trim.py fastq IN.fq OUT.fq ADAPTERS.json      ADAPTERS.json: {"adapter5": ..., "adapter3": ...} (calibrate.py adapters)

trim_read: a suffix of the 5' adapter (a partial adapter is its 3' part) aligned by edit distance to the read's first bases,
the read anchored at its start; the read is cut after the best end when the alignment covers >= MIN_ALN read bases with edits
<= MAX_ERR of them. The 3' adapter mirrors it (a partial adapter is its 5' part) at the read's end.
"""
import json
import sys

MAX_ERR, SLACK, MIN_ALN = 0.2, 20, 12


def _cut(read, adapter):
    """The read length the adapter (consumed to its END, free at its start) explains at the read's start, or 0."""
    a, r = len(adapter), min(len(read), len(adapter) + SLACK)
    prev = [0] * (a + 1)                                    # no read base yet: any adapter suffix may start here
    best = (10 ** 9, 0)
    for i in range(1, r + 1):
        cur = [i] + [0] * a
        b = read[i - 1]
        for j in range(1, a + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (b != adapter[j - 1]))
        if i >= MIN_ALN and cur[a] <= MAX_ERR * i:
            best = min(best, (cur[a], i))
        prev = cur
    return best[1]


def trim_read(read, adapter5, adapter3):
    c5 = _cut(read, adapter5)
    rest = read[c5:]
    c3 = _cut(rest[::-1], adapter3[::-1])
    return rest[:len(rest) - c3] if c3 else rest


def trim_fastq(src, dst, adapter5, adapter3):
    n = cut = 0
    with open(src) as fi, open(dst, "w") as fo:
        for h in fi:
            s = next(fi).rstrip("\n")
            next(fi)
            next(fi)
            t = trim_read(s, adapter5, adapter3)
            n += 1
            cut += t != s
            fo.write(f"{h}{t}\n+\n{'~' * len(t)}\n")
    return n, cut


def main(argv):
    if argv[0] == "fastq":
        a = json.load(open(argv[3]))
        n, cut = trim_fastq(argv[1], argv[2], a["adapter5"], a["adapter3"])
        print(f"[trim] {argv[2]}: {cut} of {n} reads trimmed", flush=True)
        return 0
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
