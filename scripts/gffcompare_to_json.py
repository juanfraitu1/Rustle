#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path


def parse_attrs(attr_field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in attr_field.split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, value = part.split(" ", 1)
        out[key] = value.strip().strip('"')
    return out


def load_reference_transcript_ids(reference_gtf: Path) -> list[str]:
    ids: list[str] = []
    seen: set[str] = set()
    with reference_gtf.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature != "transcript":
                continue
            attrs = parse_attrs(fields[8])
            tid = attrs.get("transcript_id")
            if not tid or tid in seen:
                continue
            seen.add(tid)
            ids.append(tid)
    return ids


@dataclass(frozen=True)
class StatsSummary:
    reference_mrnas: int | None = None
    query_mrnas: int | None = None
    loci: int | None = None
    multi_exon_reference_mrnas: int | None = None
    transcript_level_sensitivity: float | None = None
    transcript_level_precision: float | None = None
    matching_transcripts: int | None = None
    missed_exons_num: int | None = None
    missed_exons_den: int | None = None
    missed_introns_num: int | None = None
    missed_introns_den: int | None = None


def parse_gffcompare_stats(stats_path: Path) -> StatsSummary:
    ref_re = re.compile(r"^#\s*Reference mRNAs\s*:\s*([0-9]+)\s+in\s+([0-9]+)\s+loci\s+\(([^)]*)\)")
    qry_re = re.compile(r"^#\s*Query mRNAs\s*:\s*([0-9]+)\s+in\s+([0-9]+)\s+loci")
    txlvl_re = re.compile(r"^\s*Transcript level:\s*([0-9.]+)\s*\|\s*([0-9.]+)")
    match_re = re.compile(r"^\s*Matching transcripts:\s*([0-9]+)")
    missed_exons_re = re.compile(r"^\s*Missed exons:\s*([0-9]+)\s*/\s*([0-9]+)")
    missed_introns_re = re.compile(r"^\s*Missed introns:\s*([0-9]+)\s*/\s*([0-9]+)")

    reference_mrnas: int | None = None
    query_mrnas: int | None = None
    loci: int | None = None
    multi_exon_reference_mrnas: int | None = None
    transcript_level_sensitivity: float | None = None
    transcript_level_precision: float | None = None
    matching_transcripts: int | None = None
    missed_exons_num: int | None = None
    missed_exons_den: int | None = None
    missed_introns_num: int | None = None
    missed_introns_den: int | None = None

    with stats_path.open() as handle:
        for raw in handle:
            line = raw.rstrip("\n")

            m = ref_re.match(line)
            if m:
                reference_mrnas = int(m.group(1))
                loci = int(m.group(2))
                trailer = m.group(3)
                # e.g. "31 multi-exon"
                mm = re.search(r"([0-9]+)\s+multi-exon", trailer)
                if mm:
                    multi_exon_reference_mrnas = int(mm.group(1))
                continue

            m = qry_re.match(line)
            if m:
                query_mrnas = int(m.group(1))
                # loci may be present in both ref/qry; keep ref loci if already parsed.
                if loci is None:
                    loci = int(m.group(2))
                continue

            m = txlvl_re.match(line)
            if m:
                transcript_level_sensitivity = float(m.group(1))
                transcript_level_precision = float(m.group(2))
                continue

            m = match_re.match(line)
            if m:
                matching_transcripts = int(m.group(1))
                continue

            m = missed_exons_re.match(line)
            if m:
                missed_exons_num = int(m.group(1))
                missed_exons_den = int(m.group(2))
                continue

            m = missed_introns_re.match(line)
            if m:
                missed_introns_num = int(m.group(1))
                missed_introns_den = int(m.group(2))
                continue

    return StatsSummary(
        reference_mrnas=reference_mrnas,
        query_mrnas=query_mrnas,
        loci=loci,
        multi_exon_reference_mrnas=multi_exon_reference_mrnas,
        transcript_level_sensitivity=transcript_level_sensitivity,
        transcript_level_precision=transcript_level_precision,
        matching_transcripts=matching_transcripts,
        missed_exons_num=missed_exons_num,
        missed_exons_den=missed_exons_den,
        missed_introns_num=missed_introns_num,
        missed_introns_den=missed_introns_den,
    )


def load_tmap_exact_ref_ids(tmap_path: Path) -> set[str]:
    out: set[str] = set()
    with tmap_path.open() as handle:
        header = handle.readline()
        if not header:
            return out
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            ref_id = fields[1]
            class_code = fields[2]
            if ref_id != "-" and class_code == "=":
                out.add(ref_id)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Convert gffcompare outputs to a stable JSON summary.")
    ap.add_argument("--reference-gtf", type=Path, required=True)
    ap.add_argument("--stats", type=Path, required=True, help="gffcompare .stats file")
    ap.add_argument("--tmap", type=Path, required=True, help="gffcompare .tmap file")
    ap.add_argument("--out-json", type=Path, required=True)
    args = ap.parse_args()

    stats = parse_gffcompare_stats(args.stats)
    ref_ids = load_reference_transcript_ids(args.reference_gtf)
    exact_ref = load_tmap_exact_ref_ids(args.tmap)
    missing_exact = [rid for rid in ref_ids if rid not in exact_ref]

    payload = {
        "inputs": {
            "reference_gtf": str(args.reference_gtf),
            "stats": str(args.stats),
            "tmap": str(args.tmap),
        },
        "gffcompare": {
            "reference_mrnas": stats.reference_mrnas,
            "query_mrnas": stats.query_mrnas,
            "loci": stats.loci,
            "multi_exon_reference_mrnas": stats.multi_exon_reference_mrnas,
            "transcript_level_sensitivity": stats.transcript_level_sensitivity,
            "transcript_level_precision": stats.transcript_level_precision,
            "matching_transcripts": stats.matching_transcripts,
            "missed_exons": {
                "num": stats.missed_exons_num,
                "den": stats.missed_exons_den,
            },
            "missed_introns": {
                "num": stats.missed_introns_num,
                "den": stats.missed_introns_den,
            },
        },
        "reference": {
            "transcripts_total": len(ref_ids),
            "exact_matched_total": len(exact_ref),
            "missing_exact_total": len(missing_exact),
            "missing_exact_transcript_ids": missing_exact,
        },
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

