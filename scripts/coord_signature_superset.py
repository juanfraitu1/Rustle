#!/usr/bin/env python3
"""Assert intron-chain coord-signatures of GTF A are a superset of GTF B.

Usage: coord_signature_superset.py SUPERSET.gtf SUBSET.gtf
Exit 0 if every (chrom,strand,intron-chain) in SUBSET is present in SUPERSET.
Exit 1 (and print the missing chains) otherwise. Mirrors
pipeline.rs intron_chain_from_exons (donor=exon_end, acceptor=next_start).
"""
import sys
from collections import defaultdict


def chains(path):
    tx = defaultdict(list)  # tid -> [(start,end)]
    meta = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            chrom, strand = f[0], f[6]
            start, end = int(f[3]) - 1, int(f[4])  # GTF 1-based incl -> 0-based half-open
            tid = None
            for kv in f[8].split(";"):
                kv = kv.strip()
                if kv.startswith("transcript_id"):
                    tid = kv.split('"')[1]
                    break
            if tid is None:
                continue
            tx[tid].append((start, end))
            meta[tid] = (chrom, strand)
    sigs = set()
    for tid, exons in tx.items():
        exons.sort()
        chrom, strand = meta[tid]
        introns = tuple(
            (exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)
        )
        if introns:  # single-exon transcripts have no chain; skip (matches union idiom)
            sigs.add((chrom, strand, introns))
    return sigs


def main():
    sup, sub = sys.argv[1], sys.argv[2]
    s_sup, s_sub = chains(sup), chains(sub)
    missing = s_sub - s_sup
    if missing:
        print(f"FAIL: {len(missing)} baseline chain(s) absent from VG output", file=sys.stderr)
        for m in sorted(missing)[:20]:
            print(f"  MISSING {m[0]} {m[1]} {len(m[2])} introns", file=sys.stderr)
        sys.exit(1)
    print(f"OK: VG superset baseline ({len(s_sub)} multi-exon baseline chains all present)")
    sys.exit(0)


if __name__ == "__main__":
    main()
