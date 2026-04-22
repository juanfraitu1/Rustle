#!/usr/bin/env python3
"""Compare reference and query GTF transcripts at the junction-chain level."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple


ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


@dataclass
class Tx:
    tx_id: str
    gene_id: str
    chrom: str
    strand: str
    start: int
    end: int
    exons: List[Tuple[int, int]] = field(default_factory=list)
    attrs: Dict[str, str] = field(default_factory=dict)

    def finalize(self) -> None:
        self.exons.sort()

    @property
    def junctions(self) -> List[Tuple[int, int]]:
        return [
            (self.exons[i][1], self.exons[i + 1][0])
            for i in range(len(self.exons) - 1)
        ]


def parse_attrs(raw: str) -> Dict[str, str]:
    return {m.group(1): m.group(2) for m in ATTR_RE.finditer(raw)}


def load_gtf(path: str) -> List[Tx]:
    txs: Dict[str, Tx] = {}
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attr_raw = parts
            attrs = parse_attrs(attr_raw)
            tx_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id", "")
            if not tx_id:
                continue
            if feature == "transcript":
                txs[tx_id] = Tx(
                    tx_id=tx_id,
                    gene_id=gene_id,
                    chrom=chrom,
                    strand=strand,
                    start=int(start),
                    end=int(end),
                    attrs=attrs,
                )
            elif feature == "exon":
                tx = txs.get(tx_id)
                if tx is None:
                    tx = Tx(
                        tx_id=tx_id,
                        gene_id=gene_id,
                        chrom=chrom,
                        strand=strand,
                        start=int(start),
                        end=int(end),
                        attrs=attrs,
                    )
                    txs[tx_id] = tx
                tx.exons.append((int(start), int(end)))
                tx.start = min(tx.start, int(start))
                tx.end = max(tx.end, int(end))
    out = list(txs.values())
    for tx in out:
        tx.finalize()
    return out


def same_junction(a: Tuple[int, int], b: Tuple[int, int], tol: int) -> bool:
    return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol


def count_prefix(ref: Sequence[Tuple[int, int]], qry: Sequence[Tuple[int, int]], tol: int) -> int:
    n = min(len(ref), len(qry))
    i = 0
    while i < n and same_junction(ref[i], qry[i], tol):
        i += 1
    return i


def count_suffix(ref: Sequence[Tuple[int, int]], qry: Sequence[Tuple[int, int]], tol: int) -> int:
    n = min(len(ref), len(qry))
    i = 0
    while i < n and same_junction(ref[-(i + 1)], qry[-(i + 1)], tol):
        i += 1
    return i


def count_total_exact(ref: Sequence[Tuple[int, int]], qry: Sequence[Tuple[int, int]], tol: int) -> int:
    return sum(
        1 for a, b in zip(ref, qry) if same_junction(a, b, tol)
    )


def junc_str(j: Optional[Tuple[int, int]]) -> str:
    if j is None:
        return ""
    return f"{j[0]}-{j[1]}"


def classify(ref: Tx, qry: Tx, tol: int) -> str:
    ref_j = ref.junctions
    qry_j = qry.junctions
    prefix = count_prefix(ref_j, qry_j, tol)
    if prefix == len(ref_j) == len(qry_j):
        return "exact_chain"
    if prefix == min(len(ref_j), len(qry_j)):
        if len(ref_j) < len(qry_j):
            return "query_extension"
        if len(ref_j) > len(qry_j):
            return "query_truncation"
    if prefix >= max(1, min(len(ref_j), len(qry_j)) - 2):
        return "late_tail_divergence"
    if prefix >= 1:
        return "internal_divergence"
    return "no_prefix_match"


def best_query_for_ref(ref: Tx, queries: Sequence[Tx], tol: int, flag_filter: Optional[str]) -> Tuple[Optional[Tx], Dict[str, object]]:
    best_tx: Optional[Tx] = None
    best_key: Optional[Tuple[int, int, int, int, int]] = None
    best_info: Dict[str, object] = {}
    for qry in queries:
        flags = qry.attrs.get("flags", "")
        if flag_filter and flag_filter not in flags:
            continue
        prefix = count_prefix(ref.junctions, qry.junctions, tol)
        suffix = count_suffix(ref.junctions, qry.junctions, tol)
        total_exact = count_total_exact(ref.junctions, qry.junctions, tol)
        min_len = min(len(ref.junctions), len(qry.junctions))
        key = (
            prefix,
            total_exact,
            suffix,
            -abs(len(ref.junctions) - len(qry.junctions)),
            int(qry.attrs.get("FL_count", "0") or "0"),
        )
        if best_key is None or key > best_key:
            best_key = key
            best_tx = qry
            best_info = {
                "prefix": prefix,
                "suffix": suffix,
                "total_exact": total_exact,
                "min_len": min_len,
            }
    return best_tx, best_info


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("reference_gtf")
    ap.add_argument("query_gtf")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol", type=int, default=10)
    ap.add_argument("--flag-filter", default=None)
    args = ap.parse_args()

    refs = load_gtf(args.reference_gtf)
    queries = load_gtf(args.query_gtf)

    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "ref_id",
            "best_query_id",
            "best_query_flags",
            "best_query_fl_count",
            "ref_junctions",
            "query_junctions",
            "common_prefix_junctions",
            "common_suffix_junctions",
            "total_exact_junctions",
            "first_divergence_index_1based",
            "ref_first_divergence",
            "query_first_divergence",
            "classification",
            "ref_start",
            "ref_end",
            "query_start",
            "query_end",
        ])
        for ref in refs:
            best_qry, info = best_query_for_ref(ref, queries, args.tol, args.flag_filter)
            if best_qry is None:
                writer.writerow([
                    ref.tx_id, "", "", "", len(ref.junctions), "", 0, 0, 0, "", "", "", "no_candidate",
                    ref.start, ref.end, "", "",
                ])
                continue
            prefix = int(info["prefix"])
            ref_j = ref.junctions
            qry_j = best_qry.junctions
            first_idx: Optional[int] = None
            ref_div: Optional[Tuple[int, int]] = None
            qry_div: Optional[Tuple[int, int]] = None
            if prefix < min(len(ref_j), len(qry_j)):
                first_idx = prefix + 1
                ref_div = ref_j[prefix]
                qry_div = qry_j[prefix]
            elif len(ref_j) != len(qry_j):
                first_idx = prefix + 1
                ref_div = ref_j[prefix] if prefix < len(ref_j) else None
                qry_div = qry_j[prefix] if prefix < len(qry_j) else None
            writer.writerow([
                ref.tx_id,
                best_qry.tx_id,
                best_qry.attrs.get("flags", ""),
                best_qry.attrs.get("FL_count", ""),
                len(ref_j),
                len(qry_j),
                prefix,
                info["suffix"],
                info["total_exact"],
                first_idx or "",
                junc_str(ref_div),
                junc_str(qry_div),
                classify(ref, best_qry, args.tol),
                ref.start,
                ref.end,
                best_qry.start,
                best_qry.end,
            ])


if __name__ == "__main__":
    main()
