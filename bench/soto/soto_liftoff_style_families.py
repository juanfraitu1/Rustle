#!/usr/bin/env python3
"""Adapt LiftOff's own primitive -- align a gene's OWN sequence against the whole genome, keep hits
above an identity/coverage threshold -- into a FAMILY-DEFINITION mechanism instead of a copy-count
mechanism (docs/o1_ledger.md §6iv follow-up).

LiftOff's `-copies` asks "does gene G have undiscovered extra copies of itself in the genome?" and reports
them as satellites of G. This script asks a different question with the SAME alignment primitive: "does
gene G's own sequence land, well enough, on a DIFFERENT ALREADY-ANNOTATED gene H's own locus?" -- if so,
that is direct, real evidence G and H are related copies of one family, independent of any SD98/SEDEF
duplication call. Family = connected components of this graph (reusing soto_cluster_from_shared.py's
eligible/non-eligible-bridging + single-eligible-seed-island rules unchanged, since this is just a
DIFFERENT edge source feeding the SAME, already-validated clustering step).

WHY THIS IS A GENUINELY DIFFERENT SIGNAL FROM SEDEF+shared-exon (not just a rebrand): it aligns the
gene's OWN sequence directly, not a wide SD98-called duplication block -- avoiding exactly the unit-
averaging dilution that cost ID_328 (§6in: a 99.99%-identical gene body diluted to ~91.4% by averaging
over unnecessary flanking sequence).

Inputs (built earlier this session, reused as-is):
  --gene-bed   gene_v2_coords.bed (chrom,start,end,gene_id,score,strand) -- each SD98 gene's own v2.0
               locus (LiftOff's own lifted coordinates from cat_v4.bed, §6iv)
  --bam        gene_selfalign.bam -- every gene's own sequence aligned against the whole v2.0 genome
               (minimap2 asm20, -N 50 --secondary=yes -p 0.5, one call, no per-gene re-indexing)
"""
import argparse, csv, sys
from collections import defaultdict

import pysam


def cigar_identity(cigartuples):
    """Fraction of aligned (M/=/X, op codes 0/7/8) bases that are exact matches (=, op 7), pysam CIGAR
    ops. Falls back to None if no eqx info (plain M only, can't tell matches from mismatches)."""
    eq = x = 0
    has_eqx = False
    for op, ln in cigartuples:
        if op == 7:
            eq += ln
            has_eqx = True
        elif op == 8:
            x += ln
            has_eqx = True
    if not has_eqx or (eq + x) == 0:
        return None
    return eq / (eq + x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene-bed", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--min-identity", type=float, default=0.98)
    ap.add_argument("--min-cov", type=float, default=0.90,
                     help="hit must cover >=this fraction of the OTHER gene's own length to count as "
                          "landing on it (not just clipping its edge)")
    ap.add_argument("--out-shared", required=True)
    a = ap.parse_args()

    genes = {}  # gene_id -> (chrom, start, end)
    with open(a.gene_bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            genes[f[3]] = (f[0], int(f[1]), int(f[2]))

    by_chrom = defaultdict(list)
    for gid, (c, s, e) in genes.items():
        by_chrom[c].append((s, e, gid))
    for c in by_chrom:
        by_chrom[c].sort()

    shared = defaultdict(set)
    n_records = n_selfhit = n_other_hit = n_passed = 0
    bam = pysam.AlignmentFile(a.bam, "rb")
    for rec in bam.fetch(until_eof=True):
        if rec.is_unmapped:
            continue
        n_records += 1
        # bedtools getfasta -name emits "gene_id::chrom:start-end" as the FASTA record name --
        # strip the coordinate suffix before any gene_id lookup/comparison.
        query_gid = rec.query_name.split("::")[0]
        chrom = rec.reference_name
        hs, he = rec.reference_start, rec.reference_end
        own_c, own_s, own_e = genes.get(query_gid, (None, None, None))
        ident = cigar_identity(rec.cigartuples)
        if ident is None or ident < a.min_identity:
            continue
        for s2, e2, gid2 in by_chrom.get(chrom, ()):
            if e2 <= hs:
                continue
            if s2 >= he:
                break
            if gid2 == query_gid:
                n_selfhit += 1
                continue
            ov = min(he, e2) - max(hs, s2)
            cov_other = ov / (e2 - s2) if (e2 - s2) > 0 else 0
            n_other_hit += 1
            if cov_other >= a.min_cov:
                n_passed += 1
                shared[query_gid].add(gid2)
                shared[gid2].add(query_gid)

    print(f"[records] {n_records} mapped alignment records, {n_selfhit} self-hits (skipped), "
          f"{n_other_hit} hits on a different gene's locus at >=identity {a.min_identity}, "
          f"{n_passed} passed >=cov {a.min_cov}", file=sys.stderr)

    edge_pairs = sorted({tuple(sorted((g, p))) for g, partners in shared.items() for p in partners})
    with open(a.out_shared, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["gene_a", "gene_b"])
        w.writerows(edge_pairs)
    print(f"[done] {len(edge_pairs)} unique edges -> {a.out_shared}", file=sys.stderr)


if __name__ == "__main__":
    main()
