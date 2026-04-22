#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PanelLocus:
    label: str
    category: str
    chrom: str
    start_1based: int
    end_1based: int
    strand: str
    gene_hint: str
    note: str

    @property
    def samtools_region(self) -> str:
        return f"{self.chrom}:{self.start_1based}-{self.end_1based}"


DEFAULT_PANEL = [
    PanelLocus(
        label="good_strg88",
        category="good",
        chrom="NC_073243.2",
        start_1based=20429749,
        end_1based=20533661,
        strand="+",
        gene_hint="STRG.88",
        note="calibration locus with strong current parity",
    ),
    PanelLocus(
        label="medium_strg171",
        category="medium",
        chrom="NC_073243.2",
        start_1based=23628385,
        end_1based=23783928,
        strand="-",
        gene_hint="STRG.171",
        note="mixed seed-flow and late-family divergence",
    ),
    PanelLocus(
        label="bad_47560501",
        category="bad",
        chrom="NC_073243.2",
        start_1based=47560501,
        end_1based=47786347,
        strand="-",
        gene_hint="",
        note="top whole-run failure-rank locus dominated by zero-flux plus late collapse",
    ),
]


def parse_attrs(attr_field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in attr_field.split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, value = part.split(" ", 1)
        out[key] = value.strip().strip('"')
    return out


def overlaps(feature_chrom: str, start_1based: int, end_1based: int, locus: PanelLocus) -> bool:
    return (
        feature_chrom == locus.chrom
        and not (end_1based < locus.start_1based or start_1based > locus.end_1based)
    )


def write_panel_manifest(path: Path, loci: list[PanelLocus]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "label",
                "category",
                "chrom",
                "start_1based",
                "end_1based",
                "strand",
                "gene_hint",
                "samtools_region",
                "note",
            ]
        )
        for locus in loci:
            writer.writerow(
                [
                    locus.label,
                    locus.category,
                    locus.chrom,
                    locus.start_1based,
                    locus.end_1based,
                    locus.strand,
                    locus.gene_hint,
                    locus.samtools_region,
                    locus.note,
                ]
            )


def build_mini_gtf(reference_gtf: Path, out_gtf: Path, loci: list[PanelLocus]) -> None:
    kept_gene_ids: set[str] = set()
    transcript_rows: list[str] = []
    exon_rows: list[str] = []
    other_rows: list[str] = []

    with reference_gtf.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0]
            feature = fields[2]
            start_1based = int(fields[3])
            end_1based = int(fields[4])
            attrs = parse_attrs(fields[8])
            gene_id = attrs.get("gene_id", "")

            matched = any(overlaps(chrom, start_1based, end_1based, locus) for locus in loci)
            if feature == "transcript" and matched:
                if gene_id:
                    kept_gene_ids.add(gene_id)
                transcript_rows.append(line)
            elif gene_id and gene_id in kept_gene_ids:
                if feature == "exon":
                    exon_rows.append(line)
                else:
                    other_rows.append(line)

    with out_gtf.open("w") as handle:
        for row in transcript_rows:
            handle.write(row)
        for row in exon_rows:
            handle.write(row)
        for row in other_rows:
            handle.write(row)


def run_checked(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def unlink_if_exists(path: Path) -> None:
    try:
        path.unlink()
    except FileNotFoundError:
        pass


def build_mini_bam(input_bam: Path, out_prefix: Path, loci: list[PanelLocus]) -> tuple[Path, Path]:
    unsorted_bam = out_prefix.with_suffix(".unsorted.bam")
    sorted_bam = out_prefix.with_suffix(".bam")
    bai = sorted_bam.with_suffix(".bam.bai")

    cmd = ["samtools", "view", "-b", "-o", str(unsorted_bam), str(input_bam)]
    cmd.extend(locus.samtools_region for locus in loci)
    run_checked(cmd)
    run_checked(["samtools", "sort", "-o", str(sorted_bam), str(unsorted_bam)])
    run_checked(["samtools", "index", str(sorted_bam)])
    unlink_if_exists(unsorted_bam)
    return sorted_bam, bai


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-bam", type=Path, default=Path("GGO_19.bam"))
    ap.add_argument("--reference-gtf", type=Path, default=Path("GGO_19.gtf"))
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = out_dir / "panel_loci.tsv"
    mini_gtf = out_dir / "panel_reference.gtf"
    mini_prefix = out_dir / "panel_reads"

    write_panel_manifest(manifest, DEFAULT_PANEL)
    build_mini_gtf(args.reference_gtf, mini_gtf, DEFAULT_PANEL)
    build_mini_bam(args.input_bam, mini_prefix, DEFAULT_PANEL)

    print(f"panel_manifest\t{manifest}")
    print(f"panel_gtf\t{mini_gtf}")
    print(f"panel_bam\t{mini_prefix.with_suffix('.bam')}")
    print(f"panel_bai\t{mini_prefix.with_suffix('.bam.bai')}")


if __name__ == "__main__":
    main()
