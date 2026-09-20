#!/usr/bin/env python3
"""filter_families_by_psv.py — classify catalog families by PSV-based homology evidence.

For each family in family_rna_refine.tsv, extract copy transcript sequences,
align them with mafft, call PSV columns, and record evidence metrics.
Output is a TSV that can be used to define/reject multi-copy families.
"""
import argparse
import csv
import os
import sys
import tempfile

# reuse PSV logic from psv_split_stringtie
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import psv_split_stringtie as psv


OUT_TSV = os.path.join(psv.HERE, "family_psv_evidence.tsv")


def classify_family(family_id):
    members = psv.load_family_members(family_id)
    if not members:
        return None

    n_members = len(members)
    chroms = sorted(set(m[1] for m in members))
    n_chroms = len(chroms)
    genes = sorted(set(m[4] for m in members if m[4] != "NA"))
    annotated_gene = members[0][4] if members[0][4] != "NA" else "NA"

    # genomic span per chromosome
    by_chrom = {}
    for _, c, s, e, _ in members:
        by_chrom[c] = (min(by_chrom.get(c, (s, e))[0], s), max(by_chrom.get(c, (s, e))[1], e))
    max_span = max(e - s for s, e in by_chrom.values()) if by_chrom else 0

    with tempfile.TemporaryDirectory(dir=psv.SCRATCH, prefix="psv_evidence_") as tmp:
        try:
            names, members, refseq, seq_info = psv.extract_copy_sequences(members, tmp, use_spliced=True)
        except Exception as exc:
            return {
                "family_id": family_id,
                "n_members": n_members,
                "n_chroms": n_chroms,
                "chroms": ",".join(chroms),
                "genes": ",".join(genes),
                "status": f"extract_error:{type(exc).__name__}",
                "n_psv": 0,
                "n_aligned_cols": 0,
                "ref_len": 0,
                "max_span": max_span,
                "mean_member_len": 0,
            }

        mean_member_len = sum(len(psv.fetch_spliced_sequence if False else refseq) for _ in names)  # placeholder
        # compute mean length from seq_info
        lens = []
        for nm in names:
            info = seq_info[nm]
            if info["type"] == "spliced":
                lens.append(sum(e - s + 1 for s, e in info["exons"]))
            else:
                lens.append(info["end"] - info["start"])
        mean_member_len = sum(lens) / len(lens) if lens else 0

        try:
            psv.align_copies_mafft(tmp)
            psvs, n_aligned = psv.call_psvs_mafft(tmp, names, refseq)
        except Exception as exc:
            return {
                "family_id": family_id,
                "n_members": n_members,
                "n_chroms": n_chroms,
                "chroms": ",".join(chroms),
                "genes": ",".join(genes),
                "status": f"align_error:{type(exc).__name__}",
                "n_psv": 0,
                "n_aligned_cols": 0,
                "ref_len": len(refseq),
                "max_span": max_span,
                "mean_member_len": mean_member_len,
            }

        # simple evidence-based status
        if n_psv := len(psvs):
            if n_chroms == 1 and max_span < 100_000:
                status = "tandem_like_psv"
            elif n_chroms == 1:
                status = "same_chrom_segmental_psv"
            else:
                status = "inter_chromosomal_psv"
        else:
            if n_aligned > 0:
                status = "aligned_no_psv"
            else:
                status = "no_alignment"

        return {
            "family_id": family_id,
            "n_members": n_members,
            "n_chroms": n_chroms,
            "chroms": ",".join(chroms),
            "genes": ",".join(genes),
            "status": status,
            "n_psv": n_psv,
            "n_aligned_cols": n_aligned,
            "ref_len": len(refseq),
            "max_span": max_span,
            "mean_member_len": round(mean_member_len, 1),
        }


def main():
    parser = argparse.ArgumentParser(description="Classify families by PSV evidence")
    parser.add_argument("--family", type=int, help="single family_id (default: all)")
    args = parser.parse_args()

    if args.family is not None:
        fams = [args.family]
    else:
        fams = sorted(set(int(r["family_id"]) for r in csv.DictReader(open(psv.REFINE_TSV), delimiter="\t")))

    rows = []
    for fid in fams:
        row = classify_family(fid)
        if row is not None:
            rows.append(row)
            print(f"fam {fid:>3}: {row['status']:>22} n_psv={row['n_psv']:>4} aligned={row['n_aligned_cols']:>5} "
                  f"members={row['n_members']} chroms={row['n_chroms']} genes={row['genes']}", file=sys.stderr)

    fieldnames = ["family_id", "n_members", "n_chroms", "chroms", "genes",
                  "status", "n_psv", "n_aligned_cols", "ref_len", "max_span", "mean_member_len"]
    with open(OUT_TSV, "w", newline="") as o:
        w = csv.DictWriter(o, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"[+] wrote {OUT_TSV} ({len(rows)} families)", file=sys.stderr)


if __name__ == "__main__":
    main()
