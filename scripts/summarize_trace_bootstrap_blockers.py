#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-tsv", required=True, type=Path)
    ap.add_argument("--emitted-tsv", required=True, type=Path)
    ap.add_argument("--out-prefix", required=True, type=Path)
    return ap.parse_args()


def split_chain(chain: str) -> list[str]:
    if not chain:
        return []
    return [part for part in chain.split(",") if part]


def classify_blocker(ref_chain: list[str], emitted_chain: list[str]) -> str:
    ref_set = set(ref_chain)
    emitted_set = set(emitted_chain)
    if emitted_set < ref_set:
        return "emitted_subchain_of_ref"
    if ref_set < emitted_set:
        return "emitted_superset_of_ref"
    if ref_set == emitted_set:
        return "same_chain_boundary_shift"
    if ref_set & emitted_set:
        return "internal_junction_variant"
    return "no_shared_junctions"


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    emitted_by_id: dict[str, dict[str, str]] = {}
    with args.emitted_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            emitted_by_id[row["isoform_id"]] = row

    blocker_rows: list[dict[str, object]] = []
    with args.reference_tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["status"] not in {"partial_chain_overlap", "same_chain_boundary_mismatch", "unmatched"}:
                continue
            emitted = emitted_by_id.get(row["best_emitted"], {})
            ref_chain = split_chain(row["junction_chain"])
            emitted_chain = split_chain(emitted.get("junction_chain", ""))
            shared = len(set(ref_chain) & set(emitted_chain))
            blocker_rows.append(
                {
                    "gene_id": row["gene_id"],
                    "transcript_id": row["transcript_id"],
                    "ref_status": row["status"],
                    "best_emitted": row["best_emitted"],
                    "blocker_class": classify_blocker(ref_chain, emitted_chain),
                    "ref_num_exons": row["num_exons"],
                    "emitted_num_exons": emitted.get("num_exons", ""),
                    "emitted_read_count": emitted.get("read_count", ""),
                    "shared_junctions": shared,
                    "ref_only_junctions": len(set(ref_chain) - set(emitted_chain)),
                    "emitted_only_junctions": len(set(emitted_chain) - set(ref_chain)),
                    "ref_start": row["start"],
                    "ref_end": row["end"],
                    "emitted_start": emitted.get("start", ""),
                    "emitted_end": emitted.get("end", ""),
                    "ref_junction_chain": row["junction_chain"],
                    "emitted_junction_chain": emitted.get("junction_chain", ""),
                }
            )

    blocker_rows.sort(key=lambda row: (row["gene_id"], row["transcript_id"]))
    out_tsv = args.out_prefix.with_suffix(".tsv")
    out_md = args.out_prefix.with_suffix(".md")
    write_tsv(out_tsv, blocker_rows)

    lines = ["# Trace Bootstrap Blockers", ""]
    current_gene = None
    for row in blocker_rows:
        gene = str(row["gene_id"])
        if gene != current_gene:
            if current_gene is not None:
                lines.append("")
            lines.append(f"## {gene}")
            lines.append("")
            current_gene = gene
        lines.append(
            "- `{}`: `{}` via `{}` against `{}` (shared junctions `{}`, emitted reads `{}`)".format(
                row["transcript_id"],
                row["blocker_class"],
                row["ref_status"],
                row["best_emitted"] or ".",
                row["shared_junctions"],
                row["emitted_read_count"] or ".",
            )
        )
    out_md.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
