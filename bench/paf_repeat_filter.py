#!/usr/bin/env python3
"""Non-repeat aligned bases of gene-span PAF records (exploration after Addendum AJ).

usage: paf_repeat_filter.py annotate RMSK.out.gz IN.paf OUT.nr.tsv
       paf_repeat_filter.py filter IN.paf OUT.nr.tsv MIN_BP OUT.paf
annotate: for every record (query/target names `chrom:start-end`, 1-based spans; `cg:Z` CIGAR), walk the aligned blocks
  (M/=/X runs); a block's non-repeat bases = min(target-side bases outside interspersed repeats, query-side bases outside
  interspersed repeats). Interspersed = RepeatMasker classes LINE, SINE, LTR, Retroposon, DNA, RC, Unknown (the set
  `mcl_families --rmsk` uses). Writes line number and non-repeat aligned bases per record.
filter: keep records with non-repeat aligned bases >= MIN_BP.
"""
import bisect
import collections
import gzip
import re
import sys

CLASSES = {"LINE", "SINE", "LTR", "Retroposon", "DNA", "RC", "Unknown"}


def load_rmsk(path, contigs):
    iv = collections.defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split()
            if len(f) < 11 or not f[0].isdigit() or f[4] not in contigs:
                continue
            if f[10].split("/")[0] in CLASSES:
                iv[f[4]].append((int(f[5]) - 1, int(f[6])))
    merged = {}
    for c, v in iv.items():
        m = []
        for s, e in sorted(v):
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        # prefix sums of masked length for O(log n) masked-bases-in-interval queries
        starts = [s for s, _ in m]
        ends = [e for _, e in m]
        cum = [0]
        for s, e in m:
            cum.append(cum[-1] + e - s)
        merged[c] = (starts, ends, cum)
    return merged


def masked_in(idx, c, s, e):
    if c not in idx or e <= s:
        return 0
    starts, ends, cum = idx[c]
    i = bisect.bisect_right(ends, s)  # first interval ending after s
    j = bisect.bisect_left(starts, e)  # intervals [i, j) start before e
    if i >= j:
        return 0
    tot = cum[j] - cum[i]
    tot -= max(0, s - starts[i])
    tot -= max(0, ends[j - 1] - e)
    return max(0, tot)


def span(name):
    c, r = name.rsplit(":", 1)
    a, b = r.split("-")
    return c, int(a) - 1, int(b)


def cmd_annotate(rmsk, paf, out):
    contigs = set()
    for line in open(paf):
        f = line.split("\t", 6)
        contigs.add(f[0].rsplit(":", 1)[0])
        contigs.add(f[5].rsplit(":", 1)[0])
    idx = load_rmsk(rmsk, contigs)
    op = re.compile(r"(\d+)([MIDNSHP=X])")
    with open(out, "w") as fh:
        for ln, line in enumerate(open(paf)):
            f = line.rstrip("\n").split("\t")
            cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), None)
            qc, q0, _ = span(f[0])
            tc, t0, _ = span(f[5])
            qs, qe, ts = int(f[2]), int(f[3]), int(f[7])
            rev = f[4] == "-"
            nr = 0
            if cg is not None:
                qoff, toff = 0, 0  # consumed query bases (from qs forward, or from qe backward if rev), target bases
                for n, o in op.findall(cg):
                    n = int(n)
                    if o in "M=X":
                        tg = (t0 + ts + toff, t0 + ts + toff + n)
                        qg = (q0 + qe - qoff - n, q0 + qe - qoff) if rev else (q0 + qs + qoff, q0 + qs + qoff + n)
                        nr += min(n - masked_in(idx, tc, *tg), n - masked_in(idx, qc, *qg))
                        qoff += n
                        toff += n
                    elif o == "I":
                        qoff += n
                    elif o in "DN":
                        toff += n
            fh.write(f"{ln}\t{nr}\n")


def cmd_filter(paf, nr_tsv, min_bp, out):
    keep = {int(a) for a, b in (l.split() for l in open(nr_tsv)) if int(b) >= min_bp}
    n = 0
    with open(out, "w") as fh:
        for ln, line in enumerate(open(paf)):
            if ln in keep:
                fh.write(line)
                n += 1
    print(f"{out}: kept {n} records (non-repeat aligned >= {min_bp})")


if __name__ == "__main__":
    if sys.argv[1] == "annotate":
        cmd_annotate(*sys.argv[2:5])
    else:
        cmd_filter(sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5])
