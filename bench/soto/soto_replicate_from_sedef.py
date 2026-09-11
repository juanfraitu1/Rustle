#!/usr/bin/env python3
"""Redo Soto's SD98/shared-exon/famCN replication using a fresh, unmerged, native CHM13 v2.0 SEDEF
output, in place of the merged UCSC SD98 track that caused the merge-artifact bug documented in
tile_sd98_regions.py: bedtools-merging ~11k SD units into 817 blocks (mean 119.7 kb, max 4.25 Mb) made
each block's own perfect self-alignment the minimap2 primary, so `-p 0.5` discarded every true paralog
hit and made 69 entire Soto families structurally invisible. Tiling (20 kb / 10 kb windows) was a proxy
fix; this script uses the real thing: "Soto's unmerged SD98 unit BED", per that file's own closing note.

DEVIATION FROM SOTO'S LITERAL RECIPE, DISCLOSED. Their step 3 is "extract SD98 region FASTA, map back
to the genome with minimap2". This script skips that second alignment pass entirely and instead uses
SEDEF's OWN already-computed pairwise CIGAR (one row = one real duplication call between two specific
regions) to project exons from one side to the other -- the same "project through the CIGAR" principle
soto_replicate_clustering.py already uses (its own "ONE DECISION THAT MATTERS"), just fed by SEDEF's
alignment instead of a second minimap2 run. This avoids re-introducing ANY risk of the exact `-p 0.5`
bug above, since minimap2 is not invoked here at all. Report this as a deliberate, disclosed choice.

INPUT FORMAT (verified empirically against this file, not assumed from generic PAF/SAM convention):
34-column native SEDEF output, identical schema to the gorilla side's GGO_sedef_final.bed.
  field(1-idx)  1    2      3    4      5      6    7  8       9        ...  21        23           33
  content       chr1 start1 end1 chrom2 start2 end2 sc strand1 strand2  ...  identity  divergence   CIGAR
strand1 is always '+' (checked: 36143/36143). CIGAR semantics confirmed by direct arithmetic against
this file's own field 11/12 (M+D total == end1-start1; M+I total == end2-start2): **D consumes side1,
I consumes side2** -- the reverse of typical minimap2/PAF convention (D=target-only, I=query-only) --
so this file's CIGAR must NOT be walked with the usual SAM assumption.

COORDINATES. Input is CHM13 v2.0 (chrY present; Soto's own v1.0 SD track has none). Lifted to v1.0
using the SAME per-chromosome constant-offset table already built and validated for famCN
(famcn_from_wssd.py, fitted from soto_parCN_S1E.tsv's dual v1.0/v2.0 coordinates) -- reused verbatim,
not re-derived. A pair is dropped (and counted) if EITHER side fails to lift (straddles a regime switch,
or its chromosome has no anchors) -- never guessed, matching famcn_from_wssd.py's own stated policy.
"""
import argparse, csv, os, sys
from collections import defaultdict

BENCH = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BENCH)
from famcn_from_wssd import build_liftover, lift  # noqa: E402
from soto_replicate_clustering import load_exons  # noqa: E402


def cigar_ops(cg):
    n = 0
    for ch in cg:
        if ch.isdigit():
            n = n * 10 + (ord(ch) - 48)
        else:
            yield n, ch
            n = 0


def build_blocks(start1, ops):
    """Walk a SEDEF CIGAR (D consumes side1 only, I consumes side2 only, M/=/X consumes both).

    Returns (blocks, p1_end, bwalk_end). `blocks` is a list of (a_from, a_to, bwalk_from, bwalk_to)
    for every matched run: `a` is side1's plain genomic position (side1 is always '+'). `bwalk` is a
    strand-agnostic walk counter for side2 -- 0 at the alignment's own start, increasing monotonically
    with the CIGAR regardless of strand2 -- converted to a genomic position only via
    bwalk_to_genomic/genomic_to_bwalk, so the block list itself never needs to know strand2.
    p1_end/bwalk_end let the caller verify the CIGAR reconciles with the row's own interval lengths
    (p1_end must equal end1, bwalk_end must equal end2-start2) before trusting the row.
    """
    p1 = start1
    bw = 0
    blocks = []
    for n, op in ops:
        if op in "M=X":
            blocks.append((p1, p1 + n, bw, bw + n))
            p1 += n
            bw += n
        elif op == "D":
            p1 += n
        elif op == "I":
            bw += n
        # S/H/P not expected in an internal SEDEF self-alignment CIGAR; ignored if present.
    return blocks, p1, bw


def bwalk_to_genomic(bw, start2, end2, strand2):
    return start2 + bw if strand2 == "+" else end2 - bw


def genomic_to_bwalk(pos, start2, end2, strand2):
    return pos - start2 if strand2 == "+" else end2 - pos


def project_a_to_bwalk(blocks, a_lo, a_hi):
    """Project a-axis interval [a_lo,a_hi) through the matched blocks to a bwalk interval, or None."""
    lo = hi = None
    for a_from, a_to, bw_from, _bw_to in blocks:
        if a_to <= a_lo:
            continue
        if a_from >= a_hi:
            break
        if lo is None:
            lo = bw_from + max(a_lo, a_from) - a_from
        hi = bw_from + min(a_hi, a_to) - a_from
    return (lo, hi) if lo is not None and hi is not None and hi > lo else None


def project_bwalk_to_a(blocks, bw_lo, bw_hi):
    """Project a bwalk interval [bw_lo,bw_hi) through the matched blocks to an a-axis interval, or None."""
    lo = hi = None
    for a_from, _a_to, bw_from, bw_to in blocks:
        if bw_to <= bw_lo:
            continue
        if bw_from >= bw_hi:
            break
        if lo is None:
            lo = a_from + max(bw_lo, bw_from) - bw_from
        hi = a_from + min(bw_hi, bw_to) - bw_from
    return (lo, hi) if lo is not None and hi is not None and hi > lo else None


def find_shared_exons(blocks, side2_start_v2, side2_end_v2, side2_strand,
                       own_chrom, own_lo_v1, own_hi_v1, own_offset, own_axis_is_bwalk,
                       other_chrom, other_lo_v1, other_hi_v1, other_offset, other_axis_is_bwalk,
                       per_chrom, min_cov, shared):
    """For genes on the OTHER side (v1.0 space, restricted to [other_lo_v1,other_hi_v1)) whose exon
    is covered >=min_cov by the whole aligned span, project the exon through `blocks` onto the OWN
    side and link to whichever OWN-side gene's exon it lands in. Exactly one of
    own_axis_is_bwalk/other_axis_is_bwalk must be True: `blocks`' "a" axis is always side1 (plain),
    its "bwalk" axis is always side2 (strand-aware) -- call once with own=side1/other=side2 and once
    with the roles swapped. Returns the number of successful exon projections (for reporting only).
    """
    hits = per_chrom.get(other_chrom)
    if not hits:
        return 0
    n_proj = 0
    for s, e, g_other in hits:
        if e <= other_lo_v1:
            continue
        if s >= other_hi_v1:
            break
        if (min(e, other_hi_v1) - max(s, other_lo_v1)) / (e - s) < min_cov:
            continue
        v2_s, v2_e = max(s, other_lo_v1) - other_offset, min(e, other_hi_v1) - other_offset
        if other_axis_is_bwalk:
            bw_lo, bw_hi = sorted((
                genomic_to_bwalk(v2_s, side2_start_v2, side2_end_v2, side2_strand),
                genomic_to_bwalk(v2_e, side2_start_v2, side2_end_v2, side2_strand),
            ))
            proj = project_bwalk_to_a(blocks, bw_lo, bw_hi)
            if proj is None:
                continue
            p_v2_lo, p_v2_hi = proj
        else:
            proj = project_a_to_bwalk(blocks, v2_s, v2_e)
            if proj is None:
                continue
            g_lo = bwalk_to_genomic(proj[0], side2_start_v2, side2_end_v2, side2_strand)
            g_hi = bwalk_to_genomic(proj[1], side2_start_v2, side2_end_v2, side2_strand)
            p_v2_lo, p_v2_hi = min(g_lo, g_hi), max(g_lo, g_hi)
        p_lo, p_hi = p_v2_lo + own_offset, p_v2_hi + own_offset
        # NOTE: unlike soto_replicate_clustering.py's map-back scenario (where a region's own trivial
        # self-alignment to its OWN originating location must be excluded), there is no analogous
        # "self-mapping" artifact here: side1 and side2 are two DISTINCT loci SEDEF itself reported as
        # a real duplication pair, so a projection from side2 landing inside side1's own span is the
        # EXPECTED, correct outcome (that is literally where side1's gene lives), not a self-hit. The
        # only real self-link risk -- a gene ending up linked to itself -- is guarded below by
        # `g_own != g_other`. (An earlier version of this function incorrectly rejected every
        # same-chromosome pair here, since a projection landing inside its own OWN span is always
        # true; caught via a same-chromosome case, ID_2/CHM13_G0020704 x CHM13_G0020810 on chr16,
        # that should have linked and did not until this was removed.)
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
    ap.add_argument("--sedef", required=True, help="native CHM13 v2.0 SEDEF output (34 columns)")
    ap.add_argument("--min-identity", type=float, default=0.98, help="SD98 floor (field 21, 1-indexed)")
    ap.add_argument("--s1e", required=True, help="soto_parCN_S1E.tsv (builds the v2.0->v1.0 liftover)")
    ap.add_argument("--geneset", required=True, help="sd98_geneset_v1.tsv (gene_id, name, biotype, in_soto)")
    ap.add_argument("--cat-bed", required=True, help="CAT v4 BED (v1.0)")
    ap.add_argument("--min-cov", type=float, default=0.99, help="bedtools -f equivalent")
    ap.add_argument("--out-shared", required=True, help="TSV: gene_a<TAB>gene_b (one row per edge, deduped)")
    ap.add_argument("--limit", type=int, default=0, help="stop after N qualifying rows (0 = no limit; for smoke tests)")
    ap.add_argument("--extra-anchors",
                     help="OPT-IN: TSV (chrom, v2_pos, offset columns) of extra, individually-validated "
                          "liftover anchors -- e.g. from directly aligning a specific gene's own sequence "
                          "against both genome versions (docs/o1_ledger.md §6il) -- merged into the S1E "
                          "anchor table before fitting regimes. Omit for the original behaviour "
                          "(byte-identical to before this flag existed).")
    a = ap.parse_args()

    extra_anchors = None
    if a.extra_anchors:
        extra_anchors = []
        with open(a.extra_anchors) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                extra_anchors.append((r["chrom"], int(r["v2_pos"]), int(r["offset"])))
    table, spans = build_liftover(a.s1e, extra_anchors=extra_anchors)
    print(f"[liftover] {len(table)} chromosomes with anchors", file=sys.stderr)

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

    n_rows = n_ident = n_lift_ok = n_cigar_bad = 0
    n_proj_total = 0
    shared = defaultdict(set)

    with open(a.sedef) as fh:
        for line in fh:
            n_rows += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 33:
                continue
            chrom1, start1, end1 = f[0], int(f[1]), int(f[2])
            chrom2, start2, end2 = f[3], int(f[4]), int(f[5])
            strand2 = f[9]
            try:
                identity = float(f[20])
            except ValueError:
                continue
            if identity < a.min_identity:
                continue
            n_ident += 1

            lift1 = lift(table, spans, chrom1, start1, end1)
            lift2 = lift(table, spans, chrom2, start2, end2)
            if lift1 is None or lift2 is None:
                continue
            v1_start1, v1_end1 = lift1
            v1_start2, v1_end2 = lift2
            offset1 = v1_start1 - start1
            offset2 = v1_start2 - start2
            n_lift_ok += 1

            blocks, p1_end, bw_end = build_blocks(start1, cigar_ops(f[32]))
            if p1_end != end1 or bw_end != (end2 - start2):
                n_cigar_bad += 1
                continue

            n_proj_total += find_shared_exons(
                blocks, start2, end2, strand2,
                own_chrom=chrom1, own_lo_v1=v1_start1, own_hi_v1=v1_end1, own_offset=offset1,
                own_axis_is_bwalk=False,
                other_chrom=chrom2, other_lo_v1=v1_start2, other_hi_v1=v1_end2, other_offset=offset2,
                other_axis_is_bwalk=True,
                per_chrom=per_chrom, min_cov=a.min_cov, shared=shared,
            )
            n_proj_total += find_shared_exons(
                blocks, start2, end2, strand2,
                own_chrom=chrom2, own_lo_v1=v1_start2, own_hi_v1=v1_end2, own_offset=offset2,
                own_axis_is_bwalk=True,
                other_chrom=chrom1, other_lo_v1=v1_start1, other_hi_v1=v1_end1, other_offset=offset1,
                other_axis_is_bwalk=False,
                per_chrom=per_chrom, min_cov=a.min_cov, shared=shared,
            )

            if a.limit and n_ident >= a.limit:
                break

    print(f"[rows] {n_rows} total, {n_ident} >= identity {a.min_identity}, "
          f"{n_lift_ok} both-sides lifted (dropped {n_ident - n_lift_ok}), "
          f"{n_cigar_bad} CIGAR-length mismatches rejected", file=sys.stderr)
    print(f"[link] {n_proj_total} exon projections -> {len(shared)} genes with >=1 shared exon",
          file=sys.stderr)

    # sorted, not dict/set iteration order: this project's own convention (rustlib.py, e.g.) treats
    # hash-order-dependent output as a reproducibility defect, not a cosmetic one -- the underlying
    # edge set was already confirmed identical run-to-run, but row ORDER previously varied with
    # Python's per-run string-hash randomization, which a byte-identity check would wrongly flag.
    edge_pairs = sorted({tuple(sorted((g, p))) for g, partners in shared.items() for p in partners})
    with open(a.out_shared, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["gene_a", "gene_b"])
        w.writerows(edge_pairs)
    print(f"[done] {len(edge_pairs)} unique shared-exon edges -> {a.out_shared}", file=sys.stderr)


if __name__ == "__main__":
    main()
