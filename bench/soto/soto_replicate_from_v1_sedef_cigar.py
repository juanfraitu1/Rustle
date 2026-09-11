#!/usr/bin/env python3
"""Replace soto_replicate_from_v1_sedef.py's linear-fraction approximation with a REAL, exact CIGAR for
every Vollger et al. 2021 SD call, by directly realigning each pair's own declared interval (mappy,
minimap2's Python binding) -- so the shared-exon projection can reuse soto_replicate_from_sedef.py's
already-validated, CIGAR-exact `cigar_ops`/`build_blocks`/`find_shared_exons` VERBATIM, unmodified, instead
of duplicating that logic. This keeps the two real advantages of the v1.0-native approach (the official,
published SD calls; zero liftover, since both sides are already v1.0) while dropping the ONE disclosed
weakness (§6ir: no CIGAR in the source file, approximation costs precision) that made the first v1.0-native
attempt (ARI 0.6316) worse than the v2.0-based replication (ARI 0.6959).

MECHANISM. For each qualifying (>=min-identity) row, extract side1's and side2's declared interval from
the v1.0 genome and align side1 (query) against a small on-the-fly mappy index built from side2 (target)
-- NOT reusing Vollger's own reported bounds as if they were an exact alignment; the actual aligned
sub-region (hit.q_st/q_en, hit.r_st/r_en) is used, which can be narrower than the declared interval if the
true homology doesn't extend to its edges.

CIGAR CONVENTION, ADAPTED, DISCLOSED. mappy/minimap2 CIGARs use the STANDARD SAM/PAF convention: D consumes
the TARGET (side2) only, I consumes the QUERY (side1) only. `build_blocks` (imported unchanged from
soto_replicate_from_sedef.py) expects the REVERSE (D=side1-only, I=side2-only, matching that project's own
native SEDEF file's own convention). Swapping D<->I in the CIGAR string before parsing (`flip_indels`)
reconciles the two, letting the exact same, already-verified block-walking function run on a real,
freshly-derived CIGAR -- verified against a hand round-trip (project through then back recovers the
original interval, both strands) before trusting it on real data.

Usage: soto_replicate_from_v1_sedef_cigar.py --sedef chm13.draft_v1.0_plus38Y.SDs.bed
                                              --genome t2t-chm13-v1.0.fa.gz
                                              --geneset soto_2334_geneset.tsv --cat-bed cat_v4.bed
                                              --out-shared shared_exons_v1native_cigar.tsv
"""
import argparse, csv, os, re, sys
from collections import defaultdict

import mappy
import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from soto_replicate_from_sedef import cigar_ops, build_blocks, find_shared_exons  # noqa: E402
from soto_replicate_clustering import load_exons  # noqa: E402


def flip_indels(cigar_str):
    """mappy CIGARs: D consumes TARGET only, I consumes QUERY only (standard SAM/PAF). build_blocks
    expects the reverse (D=query/side1-only, I=target/side2-only). Swap D<->I to reconcile."""
    return re.sub(r"[DI]", lambda m: "I" if m.group() == "D" else "D", cigar_str)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sedef", required=True, help="chm13.draft_v1.0_plus38Y.SDs.bed (44 cols, header row)")
    ap.add_argument("--min-identity", type=float, default=0.98, help="SD98 floor (field 24, 1-indexed)")
    ap.add_argument("--genome", required=True, help="v1.0 genome FASTA (indexed .fai), for re-alignment")
    ap.add_argument("--geneset", required=True, help="gene_id, biotype TSV")
    ap.add_argument("--cat-bed", required=True, help="CAT v4 BED (v1.0)")
    ap.add_argument("--min-cov", type=float, default=0.99, help="bedtools -f equivalent")
    ap.add_argument("--out-shared", required=True)
    ap.add_argument("--limit", type=int, default=0, help="stop after N qualifying rows (0 = no limit)")
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

    fa = pysam.FastaFile(a.genome)
    shared = defaultdict(set)
    n_rows = n_ident = n_aligned = n_cigar_bad = n_proj_total = 0

    with open(a.sedef) as fh:
        next(fh)  # header row
        for line in fh:
            n_rows += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 24:
                continue
            chrom1, start1, end1 = f[0], int(f[1]), int(f[2])
            chrom2, start2, end2 = f[9], int(f[10]), int(f[11])
            try:
                identity = float(f[23])
            except ValueError:
                continue
            if identity < a.min_identity:
                continue
            # chrY in this assembly is GRCh38's (CHM13 itself has none) -- a separate reference this
            # script doesn't load. Confirmed earlier (docs/o1_ledger.md §6ip) that zero of our 2,334
            # genes are on chrY, so skipping these rows costs nothing and matches the already-disclosed,
            # unfixable chrY gap rather than crashing on a missing sequence.
            if chrom1 == "chrY" or chrom2 == "chrY":
                continue
            n_ident += 1

            seq1 = fa.fetch(chrom1, start1, end1)
            seq2 = fa.fetch(chrom2, start2, end2)
            aligner = mappy.Aligner(seq=seq2, preset="asm20")
            hits = sorted(aligner.map(seq1), key=lambda h: -(h.q_en - h.q_st))
            if not hits:
                continue
            hit = hits[0]
            n_aligned += 1

            own_lo = start1 + hit.q_st
            own_hi = start1 + hit.q_en
            other_lo = start2 + hit.r_st
            other_hi = start2 + hit.r_en
            other_strand = "+" if hit.strand == 1 else "-"

            ops = list(cigar_ops(flip_indels(hit.cigar_str)))
            blocks, p1_end, bw_end = build_blocks(own_lo, ops)
            if p1_end != own_hi or bw_end != (other_hi - other_lo):
                n_cigar_bad += 1
                continue

            n_proj_total += find_shared_exons(
                blocks, other_lo, other_hi, other_strand,
                own_chrom=chrom1, own_lo_v1=own_lo, own_hi_v1=own_hi, own_offset=0,
                own_axis_is_bwalk=False,
                other_chrom=chrom2, other_lo_v1=other_lo, other_hi_v1=other_hi, other_offset=0,
                other_axis_is_bwalk=True,
                per_chrom=per_chrom, min_cov=a.min_cov, shared=shared,
            )
            n_proj_total += find_shared_exons(
                blocks, other_lo, other_hi, other_strand,
                own_chrom=chrom2, own_lo_v1=other_lo, own_hi_v1=other_hi, own_offset=0,
                own_axis_is_bwalk=True,
                other_chrom=chrom1, other_lo_v1=own_lo, other_hi_v1=own_hi, other_offset=0,
                other_axis_is_bwalk=False,
                per_chrom=per_chrom, min_cov=a.min_cov, shared=shared,
            )

            if a.limit and n_ident >= a.limit:
                break
            if n_ident % 500 == 0:
                print(f"  {n_ident} rows realigned...", file=sys.stderr)

    print(f"[rows] {n_rows} total, {n_ident} >= identity {a.min_identity}, "
          f"{n_aligned} realigned OK, {n_cigar_bad} CIGAR-length mismatches rejected", file=sys.stderr)
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
