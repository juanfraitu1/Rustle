#!/usr/bin/env python3
"""Build a BED panel of loci where the query assembly has many transcripts vs reference.

Reads gffcompare ``*.loci`` (4+ columns: XLOC, region, ref_ids, query_ids).
Ranks loci by excess query transcripts (n_query - n_ref_transcripts) or by n_query.

Example:
  ./scripts/build_overprediction_panel_from_gffcompare.py \\
    GGO_19_after_termfix_cmp.loci --top 25 -o bench/ggo19_overprediction_panel

Writes:
  {out}.bed   — BED6 (0-based) for samtools -L / extract_mini_bam_from_bed.sh
  {out}.tsv — ranked table with counts and XLOC id

For fast rustle iteration, extract **one** BED row at a time (``head -1 …``); a merged
multi-row panel mini-BAM can still be expensive to assemble. ``n_query`` in the TSV is
from the full gffcompare run — check ``samtools view -c`` on your mini-BAM for local depth.

When you re-run ``gffcompare`` on a **regional** query GTF (mini-BAM output), use
``gffcompare -r <ref.gtf> -R -Q`` so sensitivity/precision are evaluated on the **same
locus overlap** as the reference (Sn correction and precision correction; ``-Q`` ignores
query transcripts that do not overlap any reference). See ``BLOCKER_PATCH_GUIDE.md``.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


REGION_RE = re.compile(
    r"^(?P<chrom>[^\[]+)\[(?P<strand>[+-])\](?P<start>\d+)-(?P<end>\d+)$"
)


def parse_region(field: str) -> tuple[str, str, int, int]:
    m = REGION_RE.match(field.strip())
    if not m:
        raise ValueError(f"unrecognized region field: {field!r}")
    chrom = m.group("chrom")
    strand = m.group("strand")
    start_1 = int(m.group("start"))
    end_1 = int(m.group("end"))
    return chrom, strand, start_1, end_1


def count_ids(column: str) -> int:
    if not column or column.strip() == "-":
        return 0
    parts = [p for p in column.split(",") if p.strip()]
    return len(parts)


def gff_to_bed6(chrom: str, strand: str, start_1: int, end_1: int, name: str, score: int) -> str:
    """GFF 1-based inclusive → BED 0-based, half-open."""
    bed_start = start_1 - 1
    bed_end = end_1
    return f"{chrom}\t{bed_start}\t{bed_end}\t{name}\t{score}\t{strand}\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("loci", type=Path, help="gffcompare *.loci file")
    ap.add_argument(
        "-o",
        "--output-prefix",
        type=Path,
        required=True,
        help="Write {prefix}.bed and {prefix}.tsv",
    )
    ap.add_argument(
        "--top",
        type=int,
        default=25,
        help="Keep top N loci by excess query transcripts (default: 25)",
    )
    ap.add_argument(
        "--min-query",
        type=int,
        default=8,
        help="Minimum query transcript count to include (default: 8)",
    )
    args = ap.parse_args()

    rows: list[dict] = []
    with args.loci.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            xloc, region, ref_col, q_col = parts[0], parts[1], parts[2], parts[3]
            try:
                chrom, strand, s1, e1 = parse_region(region)
            except ValueError:
                continue
            n_ref = count_ids(ref_col)
            n_q = count_ids(q_col)
            excess = n_q - max(1, n_ref)
            rows.append(
                {
                    "xloc": xloc,
                    "chrom": chrom,
                    "strand": strand,
                    "start_1": s1,
                    "end_1": e1,
                    "n_ref": n_ref,
                    "n_query": n_q,
                    "excess": excess,
                }
            )

    rows.sort(key=lambda r: (r["excess"], r["n_query"]), reverse=True)
    filtered = [r for r in rows if r["n_query"] >= args.min_query]
    top = filtered[: args.top]

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    bed_path = args.output_prefix.with_suffix(".bed")
    tsv_path = args.output_prefix.with_suffix(".tsv")

    with bed_path.open("w") as bed:
        for i, r in enumerate(top):
            name = f"{r['xloc']}|q{r['n_query']}|r{r['n_ref']}"
            bed.write(gff_to_bed6(r["chrom"], r["strand"], r["start_1"], r["end_1"], name, i + 1))

    with tsv_path.open("w", newline="") as tsv:
        w = csv.DictWriter(
            tsv,
            fieldnames=[
                "rank",
                "xloc",
                "chrom",
                "strand",
                "start_1based",
                "end_1based",
                "n_ref_ids",
                "n_query_ids",
                "excess_query",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for i, r in enumerate(top, start=1):
            w.writerow(
                {
                    "rank": i,
                    "xloc": r["xloc"],
                    "chrom": r["chrom"],
                    "strand": r["strand"],
                    "start_1based": r["start_1"],
                    "end_1based": r["end_1"],
                    "n_ref_ids": r["n_ref"],
                    "n_query_ids": r["n_query"],
                    "excess_query": r["excess"],
                }
            )

    print(f"Wrote {bed_path} ({len(top)} regions)")
    print(f"Wrote {tsv_path}")


if __name__ == "__main__":
    main()
