#!/usr/bin/env python3
"""What is a corrected primary worth? (user, 2026-11-19: tied multimappers pick their primary arbitrarily; for the reads that
can be told apart, does correcting it gain anything?) The ORACLE version: the truth decides which reads are corrected, so the
gain measured is the ceiling of any correction rule.

  correct_primary.py --truth READS.tsv --tie-split T.tsv --bam SEL.bam --out PREFIX [--near-only]

PREFIX.bam: every record of SEL.bam; for a molecule the truth decides (kind `one`) whose primary record is not on its origin
locus while another record is, the best-ms record on the origin becomes the primary and the old primary a secondary (flags
only; nothing is removed). --near-only: only molecules whose two best places are within 10 ms (the near ties).
PREFIX.moved.tsv: the moved molecules (read, near-tied, from, to). Seed the result with `--seed-pool primary`."""
import argparse, collections, os, sys
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True); ap.add_argument("--tie-split", required=True); ap.add_argument("--bam", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--near-only", action="store_true")
    a = ap.parse_args(argv)
    T, recs = R.load_truth(a.truth, a.tie_split), R.load_records(a.bam)
    move, n = {}, collections.Counter()
    with open(a.out + ".moved.tsv", "w") as fh:
        fh.write("read\tnear_tied\tfrom\tto\n")
        for q, rs in recs.items():
            t = T.get(q)
            if t is None or t.kind != "one":
                continue
            near = R.near_tied(rs, "ms")
            n["decided"] += 1; n["decided_near"] += near
            p = next((i for i, r in enumerate(rs) if r.primary), None)
            if p is None:
                continue
            o = next(iter(t.origin))
            if R.locus_index(rs[p], t) == o:
                continue
            n["primary_off"] += 1; n["primary_off_near"] += near
            home = [i for i, r in enumerate(rs) if R.locus_index(r, t) == o]
            if not home or (a.near_only and not near):
                continue
            to = max(home, key=lambda i: rs[i].ms)
            move[q] = (p, to)
            n["moved"] += 1; n["moved_near"] += near
            fh.write(f"{q}\t{int(near)}\t{rs[p].chrom}:{rs[p].start}\t{rs[to].chrom}:{rs[to].start}\n")
    tmp = os.environ.get("TMPDIR", os.path.dirname(os.path.abspath(a.out)))
    with pysam.AlignmentFile(a.bam, check_sq=False) as fh:
        w = pysam.AlignmentFile(f"{tmp}/cp.{os.getpid()}.bam", "wb", template=fh)
        seen = collections.Counter()
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.is_supplementary:
                w.write(r); continue
            i = seen[r.query_name]; seen[r.query_name] += 1
            m = move.get(r.query_name)
            if m:
                if i == m[1]:
                    r.flag &= ~0x100
                elif i == m[0]:
                    r.flag |= 0x100
            w.write(r)
        w.close()
    pysam.sort("-o", a.out + ".bam", w.filename.decode()); pysam.index(a.out + ".bam"); os.remove(w.filename)
    print(f"decided {n['decided']} (near-tied {n['decided_near']}); primary off its origin {n['primary_off']} (near {n['primary_off_near']}); "
          f"moved {n['moved']} (near {n['moved_near']})")

if __name__ == "__main__":
    main()
