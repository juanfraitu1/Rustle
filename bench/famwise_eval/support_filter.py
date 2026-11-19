#!/usr/bin/env python3
"""Question 3 of docs/GOAL_multicopy_families_2026-11-19.md, step 0: emulate "secondary alignments as support only" on an
assembly made with secondaries in the seeding pool: keep a multi-exon transcript only when >= --min-primary PRIMARY
records carry exactly its intron chain (contig, strand, chain); single-exon transcripts are kept when a primary record
overlaps them on the same strand. Secondaries may then have lifted a primary-backed chain over the read floor, but no
chain without a primary survives. An emulation, not the assembler mode: secondaries still shaped boundaries.

  support_filter.py IN.gtf SEL.bam OUT.gtf [--min-primary 1]"""
import argparse, collections, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools"))
from gtf_names import chain, opn

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gtf"); ap.add_argument("bam"); ap.add_argument("out"); ap.add_argument("--min-primary", type=int, default=1)
    a = ap.parse_args(argv)
    prim, spans = collections.Counter(), collections.defaultdict(list)
    for rs in R.load_records(a.bam).values():
        for r in rs:
            if r.primary:
                if r.chain:
                    prim[(r.chrom, r.strand, r.chain)] += 1
                spans[(r.chrom, r.strand)].append((r.start, r.end))
    for v in spans.values():
        v.sort()
    lines, ex, at = list(opn(a.gtf)), collections.defaultdict(list), {}
    for l in lines:
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            ex[t].append((int(f[3]), int(f[4]))); at[t] = (f[0], f[6])
    import bisect
    keep = set()
    for t, e in ex.items():
        e.sort(); c, st = at[t]
        if len(e) > 1:
            if prim[(c, st, chain(e))] >= a.min_primary:
                keep.add(t)
        else:
            v = spans.get((c, st), [])
            i = bisect.bisect_left(v, (e[0][1],))
            if any(s < e[0][1] and e[0][0] - 1 < x for s, x in v[max(0, i - 2000):i]):
                keep.add(t)
    with open(a.out, "w") as fo:
        for l in lines:
            m = re.search(r'transcript_id "([^"]+)"', l)
            if not m or m.group(1) in keep:
                fo.write(l)
    print(f"support_filter: kept {len(keep)} of {len(ex)} transcripts (multi-exon needs >= {a.min_primary} primary with the exact chain)")

if __name__ == "__main__":
    main()
