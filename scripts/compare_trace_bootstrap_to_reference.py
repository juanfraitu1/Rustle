#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class RefTx:
    gene_id: str
    transcript_id: str
    chrom: str
    strand: str
    start: int
    end: int
    exons: tuple[tuple[int, int], ...]
    junction_chain: str


@dataclass(frozen=True)
class EmittedTx:
    isoform_id: str
    gene_id: str
    chrom: str
    strand: str
    start: int
    end: int
    num_exons: int
    read_count: int
    junction_chain: str


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--emitted-tsv", required=True, nargs="+", type=Path)
    ap.add_argument("--out-prefix", required=True, type=Path)
    return ap.parse_args()


def parse_attrs(attr_field: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for chunk in attr_field.strip().split(";"):
        chunk = chunk.strip()
        if not chunk or " " not in chunk:
            continue
        key, value = chunk.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def parse_reference(gtf_path: Path) -> dict[str, list[RefTx]]:
    tx_meta: dict[str, dict[str, object]] = {}
    tx_exons: dict[str, list[tuple[int, int]]] = defaultdict(list)

    with gtf_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs_raw = fields
            attrs = parse_attrs(attrs_raw)
            tx_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id")
            if not tx_id or not gene_id:
                continue
            start_i = int(start)
            end_i = int(end)
            if feature == "transcript":
                tx_meta[tx_id] = {
                    "gene_id": gene_id,
                    "chrom": chrom,
                    "strand": strand,
                    "start": start_i,
                    "end": end_i,
                }
            elif feature == "exon":
                tx_exons[tx_id].append((start_i, end_i))

    out: dict[str, list[RefTx]] = defaultdict(list)
    for tx_id, meta in tx_meta.items():
        exons = tuple(sorted(tx_exons.get(tx_id, [])))
        if not exons:
            continue
        junctions = []
        for left, right in zip(exons, exons[1:]):
            junctions.append(f"{left[1]}-{right[0]}")
        ref = RefTx(
            gene_id=str(meta["gene_id"]),
            transcript_id=tx_id,
            chrom=str(meta["chrom"]),
            strand=str(meta["strand"]),
            start=int(meta["start"]),
            end=int(meta["end"]),
            exons=exons,
            junction_chain=",".join(junctions),
        )
        out[ref.gene_id].append(ref)
    return out


def parse_emitted(tsv_path: Path) -> tuple[str, list[EmittedTx]]:
    emitted: list[EmittedTx] = []
    gene_id = ""
    with tsv_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene_id = row["gene_id"]
            emitted.append(
                EmittedTx(
                    isoform_id=row["isoform_id"],
                    gene_id=row["gene_id"],
                    chrom=row["chrom"],
                    strand=row["strand"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    num_exons=int(row["num_exons"]),
                    read_count=int(row["read_count"]),
                    junction_chain=row["junction_chain"],
                )
            )
    return gene_id, emitted


def classify_emitted(emitted: EmittedTx, refs: list[RefTx]) -> tuple[str, str]:
    for ref in refs:
        if (
            emitted.chrom == ref.chrom
            and emitted.strand == ref.strand
            and emitted.start == ref.start
            and emitted.end == ref.end
            and emitted.junction_chain == ref.junction_chain
        ):
            return "exact_ref_match", ref.transcript_id

    same_chain = [ref for ref in refs if emitted.junction_chain == ref.junction_chain]
    if same_chain:
        best = min(
            same_chain,
            key=lambda ref: (abs(ref.start - emitted.start) + abs(ref.end - emitted.end), ref.transcript_id),
        )
        return "same_chain_boundary_mismatch", best.transcript_id

    partial = []
    emitted_juncs = set(filter(None, emitted.junction_chain.split(",")))
    for ref in refs:
        ref_juncs = set(filter(None, ref.junction_chain.split(",")))
        if emitted_juncs and ref_juncs and emitted_juncs & ref_juncs:
            partial.append(ref)
    if partial:
        return "partial_chain_overlap", partial[0].transcript_id

    return "unmatched", ""


def best_emitted_for_ref(ref: RefTx, emitted: list[EmittedTx]) -> tuple[str, str]:
    exact = [
        tx for tx in emitted
        if tx.chrom == ref.chrom
        and tx.strand == ref.strand
        and tx.start == ref.start
        and tx.end == ref.end
        and tx.junction_chain == ref.junction_chain
    ]
    if exact:
        exact.sort(key=lambda tx: (-tx.read_count, tx.isoform_id))
        return "exact_ref_match", exact[0].isoform_id

    same_chain = [tx for tx in emitted if tx.junction_chain == ref.junction_chain]
    if same_chain:
        same_chain.sort(key=lambda tx: (abs(tx.start - ref.start) + abs(tx.end - ref.end), -tx.read_count, tx.isoform_id))
        return "same_chain_boundary_mismatch", same_chain[0].isoform_id

    ref_juncs = set(filter(None, ref.junction_chain.split(",")))
    partial = []
    for tx in emitted:
        tx_juncs = set(filter(None, tx.junction_chain.split(",")))
        if ref_juncs and tx_juncs and ref_juncs & tx_juncs:
            partial.append(tx)
    if partial:
        partial.sort(key=lambda tx: (-tx.read_count, tx.isoform_id))
        return "partial_chain_overlap", partial[0].isoform_id

    return "unmatched", ""


def main() -> None:
    args = parse_args()
    refs_by_gene = parse_reference(args.reference)
    out_prefix = args.out_prefix
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    emitted_rows = []
    ref_rows = []

    for emitted_path in args.emitted_tsv:
        gene_id, emitted = parse_emitted(emitted_path)
        refs = refs_by_gene.get(gene_id, [])
        emitted_counts = Counter()
        ref_counts = Counter()

        for tx in emitted:
            status, best_ref = classify_emitted(tx, refs)
            emitted_counts[status] += 1
            emitted_rows.append({
                "gene_id": gene_id,
                "emitted_tsv": str(emitted_path),
                "isoform_id": tx.isoform_id,
                "status": status,
                "best_ref": best_ref,
                "start": tx.start,
                "end": tx.end,
                "num_exons": tx.num_exons,
                "read_count": tx.read_count,
                "junction_chain": tx.junction_chain,
            })

        for ref in refs:
            status, best_emitted = best_emitted_for_ref(ref, emitted)
            ref_counts[status] += 1
            ref_rows.append({
                "gene_id": gene_id,
                "transcript_id": ref.transcript_id,
                "status": status,
                "best_emitted": best_emitted,
                "start": ref.start,
                "end": ref.end,
                "num_exons": len(ref.exons),
                "junction_chain": ref.junction_chain,
            })

        summary_rows.append({
            "gene_id": gene_id,
            "emitted_tsv": str(emitted_path),
            "n_emitted": len(emitted),
            "n_reference": len(refs),
            "emitted_exact_ref_match": emitted_counts["exact_ref_match"],
            "emitted_same_chain_boundary_mismatch": emitted_counts["same_chain_boundary_mismatch"],
            "emitted_partial_chain_overlap": emitted_counts["partial_chain_overlap"],
            "emitted_unmatched": emitted_counts["unmatched"],
            "ref_exact_ref_match": ref_counts["exact_ref_match"],
            "ref_same_chain_boundary_mismatch": ref_counts["same_chain_boundary_mismatch"],
            "ref_partial_chain_overlap": ref_counts["partial_chain_overlap"],
            "ref_unmatched": ref_counts["unmatched"],
        })

    summary_tsv = out_prefix.with_suffix(".summary.tsv")
    emitted_tsv = out_prefix.with_suffix(".emitted.tsv")
    ref_tsv = out_prefix.with_suffix(".reference.tsv")
    summary_md = out_prefix.with_suffix(".summary.md")

    write_tsv(summary_tsv, summary_rows)
    write_tsv(emitted_tsv, emitted_rows)
    write_tsv(ref_tsv, ref_rows)

    with summary_md.open("w") as fh:
        fh.write("# Trace Bootstrap vs Reference\n\n")
        for row in summary_rows:
            fh.write(f"## {row['gene_id']}\n\n")
            fh.write(f"- emitted isoforms: `{row['n_emitted']}`\n")
            fh.write(f"- reference isoforms: `{row['n_reference']}`\n")
            fh.write(f"- emitted exact ref matches: `{row['emitted_exact_ref_match']}`\n")
            fh.write(f"- emitted same-chain boundary mismatches: `{row['emitted_same_chain_boundary_mismatch']}`\n")
            fh.write(f"- emitted partial-chain overlaps: `{row['emitted_partial_chain_overlap']}`\n")
            fh.write(f"- emitted unmatched: `{row['emitted_unmatched']}`\n")
            fh.write(f"- reference exact matches recovered: `{row['ref_exact_ref_match']}`\n")
            fh.write(f"- reference same-chain boundary mismatches: `{row['ref_same_chain_boundary_mismatch']}`\n")
            fh.write(f"- reference partial-chain overlaps: `{row['ref_partial_chain_overlap']}`\n")
            fh.write(f"- reference unmatched: `{row['ref_unmatched']}`\n\n")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
