"""PREREG r1794 (docs/PREREG_famwise_origin_2026-11-18.md), the re-toss (descriptive): what another minimap2 --seed changes.

  retoss.py records A.bam B.bam              records keyed (read, contig, start, strand, CIGAR, ms): in both, only in A,
                                             only in B; reads whose primary record differs
  retoss.py claims NAME=GTF[,GTF...] [...]   per arm: intron chains in any of its GTFs against chains in every one"""
import os, sys
import pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from origin_score import load_claims

def keyed(bam):
    """-> (set of record keys, {read: key of its primary record})"""
    recs, prim = set(), {}
    with pysam.AlignmentFile(bam, check_sq=False) as fh:
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.is_supplementary:
                continue
            k = (r.query_name, r.reference_name, r.reference_start, r.is_reverse, r.cigarstring, r.get_tag("ms"))
            recs.add(k)
            if not r.is_secondary:
                prim[r.query_name] = k
    return recs, prim

def main():
    if sys.argv[1] == "records":
        (a, pa), (b, pb) = keyed(sys.argv[2]), keyed(sys.argv[3])
        rows = [("records_in_both", len(a & b)), ("only_A", len(a - b)), ("only_B", len(b - a)),
                ("reads_primary_differs", sum(pa.get(q) != pb.get(q) for q in set(pa) | set(pb))), ("reads", len(set(pa) | set(pb)))]
    else:
        rows = [("arm", "gtfs", "chains_in_any", "chains_in_every", "not_in_every", "share")]
        for x in sys.argv[2:]:
            n, gs = x.split("=", 1)
            cs = [set(load_claims(g)) for g in gs.split(",")]
            every, anyc = set.intersection(*cs), set.union(*cs)
            rows.append((n, len(cs), len(anyc), len(every), len(anyc - every), f"{len(anyc - every) / max(1, len(anyc)):.4f}"))
    print("\n".join("\t".join(map(str, r)) for r in rows))

if __name__ == "__main__":
    main()
