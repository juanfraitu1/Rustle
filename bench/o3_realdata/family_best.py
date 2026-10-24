#!/usr/bin/env python3
"""An independent re-implementation of `copy_assign --psv-family-best` (prereg docs/archive/2026-10/merged/PREREG_o3_family_best_2026-10-15.md): which
hosted (family, read) pairs have a best alignment in their family, and where the others' best alignments lie.

    family_best.py PREFIX BAM COPIES.tsv ALL_COPIES.tsv OUT.tsv

Hosted pairs: PREFIX.assignments.tsv rows whose catalog_copy_idx names a copy of the family in COPIES.tsv (partner rows,
member_status == partner, are no copy). Per read: its primary and secondary records (supplementary and unmapped skipped),
the records with its highest AS (every tie kept; no AS on any record -> the primary); kept when one of them has a block
(the CIGAR's exon blocks, as `crate::bam::exons_from_cigar`) inside a copy span (start-end) of the family.
OUT.tsv: family, read, host, kept, best (the best loci) and best_families (the ALL_COPIES families those loci touch, or
`outside_catalog`).
"""
import collections
import csv
import sys

import pysam


def blocks(ref_start, cigartuples):
    """`cigar_exons` of src/bam.rs: exon blocks split at N (D inside a block, I never), 0-based half-open."""
    exons, n_ref, exstart, started, intron, ins = [], 0, ref_start, False, False, False
    for op, n in cigartuples:
        if op in (0, 7, 8):
            started, n_ref, intron, ins = True, n_ref + n, False, False
        elif op == 2:
            n_ref, ins = n_ref + n, False
        elif op == 1:
            ins = True
        elif op == 3:
            if not started:
                continue
            if not ins or not intron:
                if ref_start + n_ref > exstart:
                    exons.append((exstart, ref_start + n_ref))
            n_ref += n
            exstart, intron = ref_start + n_ref, True
        elif op in (4, 5):
            ins = False
    if not intron and ref_start + n_ref > exstart:
        exons.append((exstart, ref_start + n_ref))
    return exons


def best_loci(recs):
    """recs: [(is_primary, AS or None, chrom, blocks)] -> the (chrom, blocks) at the highest AS, or the primary's when no record
    carries AS."""
    scored = [r for r in recs if r[1] is not None]
    if not scored:
        return [(c, b) for p, _, c, b in recs if p][:1]
    top = max(r[1] for r in scored)
    return [(c, b) for _, s, c, b in scored if s == top]


def in_family(loci, spans):
    return any(sc == c and any(be > s and bs < e for bs, be in bl) for c, bl in loci for sc, s, e in spans)


def spans_of(copies_p, skip_partners=True):
    out = collections.defaultdict(list)
    with open(copies_p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if skip_partners and r.get("member_status") == "partner":
                continue
            out[r["family_id"]].append((int(r["copy_idx"]), r["chrom"], int(r["start"]), int(r["end"])))
    return out


def main(argv):
    prefix, bam_p, copies_p, all_p, out_p = argv
    fam = spans_of(copies_p)
    allc = spans_of(all_p, skip_partners=False)
    hosted = []
    with open(prefix + ".assignments.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            idx = r["catalog_copy_idx"]
            if idx.isdigit() and any(c[0] == int(idx) for c in fam.get(r["family_id"], ())):
                hosted.append((r["family_id"], r["read_name"], idx))
    need = {h[1] for h in hosted}
    recs = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_p) as bam:
        for rec in bam.fetch(until_eof=True):
            if rec.is_unmapped or rec.is_supplementary or rec.query_name not in need or rec.cigartuples is None:
                continue
            recs[rec.query_name].append((not rec.is_secondary, rec.get_tag("AS") if rec.has_tag("AS") else None,
                                         rec.reference_name, blocks(rec.reference_start, rec.cigartuples)))
    by_chrom = collections.defaultdict(list)
    for f, cs in allc.items():
        for _, c, s, e in cs:
            by_chrom[c].append((s, e, f))
    rows, kept = [], 0
    for f, read, host in hosted:
        loci = best_loci(recs.get(read, []))
        k = in_family(loci, [(c, s, e) for _, c, s, e in fam[f]])
        kept += k
        touched = sorted({g for c, bl in loci for s, e, g in by_chrom.get(c, ()) if any(be > s and bs < e for bs, be in bl)})
        rows.append(dict(family=f, read=read, host=host, kept=int(k),
                         best=";".join(f"{c}:{bl[0][0]}-{bl[-1][1]}" for c, bl in loci if bl),
                         best_families=",".join(touched) if touched else "outside_catalog"))
    with open(out_p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["family", "read", "host", "kept", "best", "best_families"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"{kept} of {len(hosted)} hosted (family, read) pairs have a best alignment in the family -> {out_p}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
