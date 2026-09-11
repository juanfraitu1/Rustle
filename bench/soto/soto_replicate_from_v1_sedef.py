#!/usr/bin/env python3
"""Replicate Soto's SD98/shared-exon step NATIVELY on CHM13 v1.0 -- no liftover anywhere, using the
official, PUBLISHED SEDEF-derived SD calls from Vollger et al. 2021 ("Segmental duplications and their
variation in a complete human genome", Zenodo 10.5281/zenodo.4726156,
`chm13.draft_v1.0_plus38Y.SDs.bed`) in place of a fresh SEDEF run on v2.0
(soto_replicate_from_sedef.py). Since both this file's coordinates AND cat_v4.bed's gene annotation are
v1.0-native, this sidesteps every acrocentric v2.0->v1.0 liftover issue investigated this session
(docs/o1_ledger.md §6il/§6in/§6iq) by construction -- there is nothing to lift.

DEVIATION FROM soto_replicate_from_sedef.py, DISCLOSED. This file has NO CIGAR string (44 columns, only
aggregate match/mismatch/indel COUNTS per pair, not positions) -- exon coordinates are projected between
the two sides of a pair via LINEAR-FRACTION interpolation across each side's own outer alignment bounds
(start/end), not a CIGAR walk. This is a real, different-in-kind approximation, not equivalent precision
by another route: indels within the aligned block are not individually accounted for, so a projected
position can drift from the true, CIGAR-exact position by however much indel content sits between the
block's start and the fraction point. For this file's own >=90%-identity floor, drift is typically small
relative to block length, but this is disclosed, not claimed away.

INPUT SCHEMA (44 cols, header row present; verified via the header row itself, not assumed):
field(1-idx) 1 chrom1, 2 start1, 3 end1, 6 strand1 (always '+', checked: 82578/82578), 10 chrom2,
11 start2, 12 end2, 14 strand2 (+/-), 24 fracMatch (identity). `chm13.draft_v1.0_plus38Y.SDs.lowid.bed`
is a SEPARATE run (0 shared exact (chrom1,start1,end1,chrom2,start2,end2) tuples with SDs.bed, checked),
not a strict identity-filtered subset/superset -- pass both via repeated --sedef if using the
lower-identity supplement; don't assume one contains the other.

Usage: soto_replicate_from_v1_sedef.py --sedef chm13.draft_v1.0_plus38Y.SDs.bed [--sedef ...more]
                                        --geneset soto_2334_geneset.tsv --cat-bed cat_v4.bed
                                        --out-shared shared_exons_v1native.tsv
"""
import argparse, csv, os, sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from soto_replicate_clustering import load_exons  # noqa: E402


def project_linear(gs, ge, from_lo, from_hi, to_lo, to_hi, forward):
    """Map [gs,ge] (within [from_lo,from_hi)) onto [to_lo,to_hi) by linear-fraction interpolation.
    `forward` = (strand2 == '+') works identically for EITHER projection direction (side1->side2 or
    side2->side1): the fraction is always measured from `from_lo`, and applied from `to_lo` when the
    pair aligns in the same orientation, or from `to_hi` (backwards) when it doesn't -- verified by hand
    on a synthetic reverse-strand round-trip (project forward then back recovers the original interval
    exactly) before trusting it on real data.
    """
    span_from = from_hi - from_lo
    span_to = to_hi - to_lo
    if span_from <= 0 or span_to <= 0:
        return None
    frac_s = (gs - from_lo) / span_from
    frac_e = (ge - from_lo) / span_from
    if forward:
        return (to_lo + frac_s * span_to, to_lo + frac_e * span_to)
    return (to_hi - frac_e * span_to, to_hi - frac_s * span_to)


def find_shared(own_chrom, own_lo, own_hi, other_chrom, other_lo, other_hi, forward,
                 per_chrom, min_cov, shared):
    """For genes on the OTHER side whose exon is covered >=min_cov by [other_lo,other_hi) (as a fraction
    of the OTHER gene's OWN exon length -- same convention as soto_replicate_from_sedef.py), project the
    covered portion onto the OWN side and link to whichever OWN-side gene's exon overlaps it AT ALL (no
    coverage minimum on the own side either, matching that same convention).
    """
    hits = per_chrom.get(other_chrom)
    if not hits:
        return 0
    n_proj = 0
    for s, e, g_other in hits:
        if e <= other_lo:
            continue
        if s >= other_hi:
            break
        cs, ce = max(s, other_lo), min(e, other_hi)
        if (ce - cs) / (e - s) < min_cov:
            continue
        proj = project_linear(cs, ce, other_lo, other_hi, own_lo, own_hi, forward)
        if proj is None:
            continue
        p_lo, p_hi = sorted(proj)
        n_proj += 1
        for s2, e2, g_own in per_chrom.get(own_chrom, ()):
            if e2 <= p_lo:
                continue
            if s2 >= p_hi:
                break
            if g_own != g_other:
                shared[g_own].add(g_other)
                shared[g_other].add(g_own)
    return n_proj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sedef", action="append", required=True,
                     help="native v1.0 SD BED (44 cols, header row); repeat for multiple files "
                          "(e.g. SDs.bed + SDs.lowid.bed)")
    ap.add_argument("--min-identity", type=float, default=0.98, help="SD98 floor (field 24, 1-indexed)")
    ap.add_argument("--geneset", required=True, help="gene_id, biotype TSV")
    ap.add_argument("--cat-bed", required=True, help="CAT v4 BED (v1.0)")
    ap.add_argument("--min-cov", type=float, default=0.99, help="bedtools -f equivalent")
    ap.add_argument("--out-shared", required=True, help="TSV: gene_a<TAB>gene_b (one row per edge, deduped)")
    ap.add_argument("--limit", type=int, default=0, help="stop after N qualifying rows total (0 = no limit)")
    a = ap.parse_args()

    genes = set()
    with open(a.geneset) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes.add(r["gene_id"])
    exons, _meta = load_exons(a.cat_bed, genes)
    per_chrom = defaultdict(list)
    for g, evs in exons.items():
        for c, s, e in evs:
            per_chrom[c].append((s, e, g))
    for c in per_chrom:
        per_chrom[c].sort()
    print(f"[genes] {len(genes)} SD98 genes, {sum(len(v) for v in exons.values())} exons", file=sys.stderr)

    n_rows = n_ident = n_proj_total = 0
    shared = defaultdict(set)

    for sedef_path in a.sedef:
        with open(sedef_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                n_rows += 1
                f = line.rstrip("\n").split("\t")
                if len(f) < 24:
                    continue
                chrom1, start1, end1 = f[0], int(f[1]), int(f[2])
                chrom2, start2, end2 = f[9], int(f[10]), int(f[11])
                strand2 = f[13]
                try:
                    identity = float(f[23])
                except ValueError:
                    continue
                if identity < a.min_identity:
                    continue
                n_ident += 1
                forward = strand2 == "+"

                n_proj_total += find_shared(chrom1, start1, end1, chrom2, start2, end2, forward,
                                             per_chrom, a.min_cov, shared)
                n_proj_total += find_shared(chrom2, start2, end2, chrom1, start1, end1, forward,
                                             per_chrom, a.min_cov, shared)

                if a.limit and n_ident >= a.limit:
                    break

    print(f"[rows] {n_rows} total, {n_ident} >= identity {a.min_identity}", file=sys.stderr)
    print(f"[link] {n_proj_total} exon projections -> {len(shared)} genes with >=1 shared exon",
          file=sys.stderr)

    edge_pairs = sorted({tuple(sorted((g, p))) for g, partners in shared.items() for p in partners})
    with open(a.out_shared, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["gene_a", "gene_b"])
        w.writerows(edge_pairs)
    print(f"[done] {len(edge_pairs)} unique shared-exon edges -> {a.out_shared}", file=sys.stderr)


if __name__ == "__main__":
    main()
