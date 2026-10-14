#!/usr/bin/env python3
"""CHM13 v1.0 -> v2.0 coordinates through the T2T v1.0 -> v1.1 chain (v2.0 = v1.1 + the HG002 chrY).

    v1_to_v2.py lift CHAIN IN OUT [--pairs] [--one-based] [--report R.json] [--verify V1_FA V2_FA [--n 1000] [--seed 7]]

An interval lifts only if both of its ends fall inside aligned blocks of its chromosome's chain. An end inside a chain
gap (an indel, the rDNA arrays) or on a chromosome without a chain (GRCh38 chrY in the plus38Y files) makes it
unliftable: the row is dropped and counted, never shifted. An interval whose ends sit in different blocks spans an
indel; it lifts with the indel and is counted: `rows_crossing_gap` (kept rows with such a side) and
`sides_crossing_gap` (those sides). `--pairs` lifts columns 1-3 and 4-6 (SEDEF pairs; the row is dropped if
either side fails); `--one-based` reads and writes GFF-style closed starts.
"""
import argparse
import bisect
import collections
import json
import random
import subprocess
import sys


def read_chain(lines):
    """{chrom: [(v1_start, v2_start, size)]} — one + strand chain per chromosome, same name on both sides (asserted)."""
    out, cur, t, q = {}, None, 0, 0
    for line in lines:
        f = line.split()
        if not f:
            continue
        if f[0] == "chain":
            assert f[2] == f[7] and f[4] == f[9] == "+", line
            assert f[2] not in out, f"more than one chain for {f[2]}"
            cur, t, q = out.setdefault(f[2], []), int(f[5]), int(f[10])
            continue
        size = int(f[0])
        cur.append((t, q, size))
        if len(f) == 3:
            t += size + int(f[1])
            q += size + int(f[2])
    return out


class Lifter:
    def __init__(self, chain):
        self.blocks = chain
        self.starts = {c: [b[0] for b in bs] for c, bs in chain.items()}

    def _block(self, c, x):
        i = bisect.bisect_right(self.starts[c], x) - 1
        b = self.blocks[c]
        return i if i >= 0 and x < b[i][0] + b[i][2] else None

    def lift(self, chrom, s, e):
        """0-based half-open [s, e) on v1.0 -> (s2, e2, crosses_gap) on v2.0, or None."""
        if chrom not in self.blocks or e <= s:
            return None
        i, j = self._block(chrom, s), self._block(chrom, e - 1)
        if i is None or j is None:
            return None
        b = self.blocks[chrom]
        s2 = s - b[i][0] + b[i][1]
        e2 = e - 1 - b[j][0] + b[j][1] + 1
        return (s2, e2, i != j) if e2 > s2 else None


def lift_rows(lifter, lines, pairs=False, one_based=False):
    """-> (lifted lines, counts, [(v1 region, v2 region)] per lifted interval) for TSV/BED lines."""
    off = 1 if one_based else 0
    out, c, regions = [], collections.Counter(), []
    for line in lines:
        if line.startswith(("#", "track", "browser")) or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        new, got, cross = list(f), [], []
        for ci, si, ei in [(0, 1, 2)] + ([(3, 4, 5)] if pairs else []):
            r = lifter.lift(f[ci], int(f[si]) - off, int(f[ei]))
            if r is None:
                c[f"unliftable|{f[ci]}"] += 1
                got = None
                break
            new[si], new[ei] = str(r[0] + off), str(r[1])
            cross.append(r[2])
            got.append(((f[ci], int(f[si]) - off, int(f[ei])), (f[ci], r[0], r[1])))
        if got is None:
            c["rows_dropped"] += 1
            continue
        c["rows_crossing_gap"] += any(cross)
        c["sides_crossing_gap"] += sum(cross)
        if not pairs and len(f) >= 8 and f[6] == f[1] and f[7] == f[2]:  # BED thickStart/thickEnd follow the span
            new[6], new[7] = new[1], new[2]
        out.append("\t".join(new) + "\n")
        c["rows_lifted"] += 1
        regions.extend(got)
    return out, dict(sorted(c.items())), regions


def fetch(fa, regions):
    """samtools faidx -> {region string: upper-case sequence}; regions are 0-based half-open (chrom, s, e)."""
    names = [f"{c}:{s + 1}-{e}" for c, s, e in regions]
    seqs, cur = {}, None
    for i in range(0, len(names), 500):
        res = subprocess.run(["samtools", "faidx", fa] + names[i:i + 500], capture_output=True, text=True, check=True)
        for line in res.stdout.splitlines():
            if line.startswith(">"):
                cur = line[1:]
                seqs[cur] = []
            else:
                seqs[cur].append(line.upper())
    return {k: "".join(v) for k, v in seqs.items()}


def verify(v1_fa, v2_fa, pairs_of_regions, n, seed):
    """Sequence identity of n sampled lifted intervals (those not spanning an indel)."""
    rng = random.Random(seed)
    sample = [p for p in pairs_of_regions if (p[0][2] - p[0][1]) == (p[1][2] - p[1][1])]
    sample = rng.sample(sample, min(n, len(sample)))
    a, b = fetch(v1_fa, [p[0] for p in sample]), fetch(v2_fa, [p[1] for p in sample])
    same = sum(a[f"{p[0][0]}:{p[0][1] + 1}-{p[0][2]}"] == b[f"{p[1][0]}:{p[1][1] + 1}-{p[1][2]}"] for p in sample)
    return {"checked": len(sample), "identical": same, "fraction": same / len(sample) if sample else None}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("lift")
    p.add_argument("chain"); p.add_argument("src"); p.add_argument("dst"); p.add_argument("--report")
    p.add_argument("--pairs", action="store_true"); p.add_argument("--one-based", action="store_true")
    p.add_argument("--verify", nargs=2, metavar=("V1_FA", "V2_FA"))
    p.add_argument("--n", type=int, default=1000); p.add_argument("--seed", type=int, default=7)
    a = ap.parse_args(argv)
    lifter = Lifter(read_chain(open(a.chain)))
    out, counts, regions = lift_rows(lifter, open(a.src), a.pairs, a.one_based)
    with open(a.dst, "w") as fh:
        fh.writelines(out)
    rep = {"src": a.src, "dst": a.dst, "counts": counts}
    if a.verify:
        rep["verify"] = verify(a.verify[0], a.verify[1], regions, a.n, a.seed)
    if a.report:
        with open(a.report, "w") as fh:
            json.dump(rep, fh, indent=1)
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    sys.exit(main())
