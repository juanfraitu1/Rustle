#!/usr/bin/env python3
"""Summarize exact-chain boundary candidates from locus path tables and emitted GTFs."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Tuple

import pysam


ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


@dataclass
class Tx:
    tx_id: str
    chrom: str
    strand: str
    start: int
    end: int
    attrs: Dict[str, str] = field(default_factory=dict)
    exons: List[Tuple[int, int]] = field(default_factory=list)

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
            if not tx_id:
                continue
            if feature == "transcript":
                txs[tx_id] = Tx(
                    tx_id=tx_id,
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


def normalized_chain(chain: Iterable[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    return tuple(sorted((min(a, b), max(a, b)) for a, b in chain))


def parse_path_chain(raw: str) -> Tuple[Tuple[int, int], ...]:
    pairs = []
    for chunk in raw.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        a, b = chunk.split("-")
        pairs.append((int(a), int(b)))
    return normalized_chain(pairs)


def parse_path_bounds(nodes_raw: str) -> Tuple[int, int]:
    values = [int(v) for v in nodes_raw.split(",") if v]
    if not values:
        return (0, 0)
    return (min(values[0], values[-1]), max(values[0], values[-1]))


def path_row_chain(row: dict[str, str]) -> Tuple[Tuple[int, int], ...]:
    raw = row.get("junction_chain") or row.get("derived_junction_chain") or row.get("path_junction_chain") or ""
    return parse_path_chain(raw) if raw else tuple()


def path_row_bounds(row: dict[str, str]) -> Tuple[int, int]:
    if row.get("derived_start") and row.get("derived_end"):
        return (int(row["derived_start"]), int(row["derived_end"]))
    raw = row.get("nodes") or row.get("path_nodes") or ""
    return parse_path_bounds(raw)


def path_row_exact_support(row: dict[str, str]) -> int:
    return int(row.get("count") or row.get("exact_count") or row.get("total_count") or 0)


def path_row_total_support(row: dict[str, str]) -> int:
    return int(row.get("total_count") or row.get("count") or row.get("exact_count") or 0)


def iter_read_junctions_1bp(rec) -> List[Tuple[int, int]]:
    pos = rec.reference_start + 1
    out = []
    for op, length in rec.cigartuples or []:
        if op in (0, 7, 8):
            pos += length
        elif op == 3:
            out.append((pos - 1, pos + length))
            pos += length
        elif op in (2, 6):
            pos += length
    return out


def exact_chain_boundary_counts(
    bam_path: str,
    chrom: str,
    start: int,
    end: int,
    chain: Tuple[Tuple[int, int], ...],
    tol: int,
) -> Counter[Tuple[int, int]]:
    bam = pysam.AlignmentFile(bam_path, "rb")
    counts: Counter[Tuple[int, int]] = Counter()
    for rec in bam.fetch(chrom, max(0, start - 1), end):
        if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
            continue
        juncs = iter_read_junctions_1bp(rec)
        if len(juncs) != len(chain):
            continue
        ok = True
        for (a, b), (c, d) in zip(juncs, chain):
            if abs(a - c) > tol or abs(b - d) > tol:
                ok = False
                break
        if ok:
            counts[(rec.reference_start + 1, rec.reference_end)] += 1
    return counts


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-gtf", required=True)
    ap.add_argument("--query-gtf", required=True)
    ap.add_argument("--path-table", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--tol", type=int, default=10)
    args = ap.parse_args()

    refs = load_gtf(args.reference_gtf)
    queries = load_gtf(args.query_gtf)

    queries_by_chain: Dict[Tuple[Tuple[int, int], ...], List[Tx]] = defaultdict(list)
    for tx in queries:
        queries_by_chain[normalized_chain(tx.junctions)].append(tx)

    path_rows_by_chain: Dict[Tuple[Tuple[int, int], ...], List[dict]] = defaultdict(list)
    with open(args.path_table, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chain = path_row_chain(row)
            if not chain:
                continue
            path_rows_by_chain[chain].append(row)

    summary_path = f"{args.out_prefix}.summary.tsv"
    candidate_path = f"{args.out_prefix}.candidates.tsv"
    with open(summary_path, "w", newline="", encoding="utf-8") as summary_handle, open(
        candidate_path, "w", newline="", encoding="utf-8"
    ) as candidate_handle:
        summary_writer = csv.writer(summary_handle, delimiter="\t")
        candidate_writer = csv.writer(candidate_handle, delimiter="\t")
        summary_writer.writerow(
            [
                "ref_id",
                "ref_start",
                "ref_end",
                "junctions",
                "raw_exact_boundary_count",
                "raw_has_ref_boundary",
                "path_exact_candidate_count",
                "path_has_ref_boundary",
                "best_path_support",
                "best_path_boundary",
                "emitted_same_chain_count",
                "emitted_has_ref_boundary",
                "best_emitted_fl_count",
                "best_emitted_boundary",
            ]
        )
        candidate_writer.writerow(
            [
                "ref_id",
                "source",
                "support",
                "start",
                "end",
                "is_ref_boundary",
                "kind",
                "tx_id",
            ]
        )

        for ref in refs:
            chain = normalized_chain(ref.junctions)
            raw_counts = exact_chain_boundary_counts(
                args.bam, ref.chrom, ref.start, ref.end, chain, args.tol
            )
            path_rows = path_rows_by_chain.get(chain, [])
            emitted = queries_by_chain.get(chain, [])

            for (start, end), count in raw_counts.most_common():
                candidate_writer.writerow(
                    [ref.tx_id, "bam_exact_chain", count, start, end, int((start, end) == (ref.start, ref.end)), "raw", ""]
                )
            for row in sorted(
                path_rows,
                key=lambda r: (path_row_exact_support(r), path_row_total_support(r)),
                reverse=True,
            ):
                start, end = path_row_bounds(row)
                exact_support = path_row_exact_support(row)
                total_support = path_row_total_support(row)
                if row.get("source"):
                    candidate_writer.writerow(
                        [
                            ref.tx_id,
                            row["source"],
                            exact_support,
                            start,
                            end,
                            int((start, end) == (ref.start, ref.end)),
                            "path",
                            "",
                        ]
                    )
                else:
                    candidate_writer.writerow(
                        [
                            ref.tx_id,
                            "exact_full",
                            exact_support,
                            start,
                            end,
                            int((start, end) == (ref.start, ref.end)),
                            "path",
                            "",
                        ]
                    )
                    candidate_writer.writerow(
                        [
                            ref.tx_id,
                            "total_support",
                            total_support,
                            start,
                            end,
                            int((start, end) == (ref.start, ref.end)),
                            "path",
                            "",
                        ]
                    )
            for tx in sorted(
                emitted,
                key=lambda t: int(t.attrs.get("FL_count", "0") or "0"),
                reverse=True,
            ):
                candidate_writer.writerow(
                    [
                        ref.tx_id,
                        tx.attrs.get("flags", ""),
                        int(tx.attrs.get("FL_count", "0") or "0"),
                        tx.start,
                        tx.end,
                        int((tx.start, tx.end) == (ref.start, ref.end)),
                        "emitted",
                        tx.tx_id,
                    ]
                )

            best_path = None
            if path_rows:
                best_row = max(path_rows, key=lambda r: (path_row_exact_support(r), path_row_total_support(r)))
                start, end = path_row_bounds(best_row)
                best_path = (path_row_exact_support(best_row), f"{start}-{end}")

            best_emitted = None
            if emitted:
                best_tx = max(emitted, key=lambda t: int(t.attrs.get("FL_count", "0") or "0"))
                best_emitted = (
                    int(best_tx.attrs.get("FL_count", "0") or "0"),
                    f"{best_tx.start}-{best_tx.end}",
                )

            summary_writer.writerow(
                [
                    ref.tx_id,
                    ref.start,
                    ref.end,
                    len(chain),
                    sum(raw_counts.values()),
                    int((ref.start, ref.end) in raw_counts),
                    len(path_rows),
                    int(any(path_row_bounds(row) == (ref.start, ref.end) for row in path_rows)),
                    best_path[0] if best_path else 0,
                    best_path[1] if best_path else "",
                    len(emitted),
                    int(any((tx.start, tx.end) == (ref.start, ref.end) for tx in emitted)),
                    best_emitted[0] if best_emitted else 0,
                    best_emitted[1] if best_emitted else "",
                ]
            )


if __name__ == "__main__":
    main()
