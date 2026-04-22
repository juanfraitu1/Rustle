#!/usr/bin/env python3
"""Audit where terminal boundary drift first appears for a locus."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

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


def parse_nodes(path: str) -> Dict[int, Tuple[int, str]]:
    out: Dict[int, Tuple[int, str]] = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for idx, row in enumerate(reader):
            node_idx = int(row["node_idx"]) if "node_idx" in row and row["node_idx"] else idx
            node_type = row.get("node_types") or row.get("node_type") or ""
            out[node_idx] = (int(row["position"]), node_type)
    return out


def parse_path_table(path: str) -> Dict[Tuple[Tuple[int, int], ...], List[dict]]:
    by_chain: Dict[Tuple[Tuple[int, int], ...], List[dict]] = defaultdict(list)
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chain = tuple(
                sorted(
                    tuple(sorted(map(int, chunk.split("-"))))
                    for chunk in row["junction_chain"].split(",")
                    if chunk
                )
            )
            by_chain[chain].append(row)
    return by_chain


def parse_node_list(raw: str) -> List[int]:
    return [int(x) for x in raw.split(",") if x]


def best_path_candidate(rows: List[dict]) -> Optional[dict]:
    if not rows:
        return None
    ranked = sorted(
        rows,
        key=lambda r: (
            int(r["count"]),
            float(r["weight"]),
            len(parse_node_list(r["nodes"])),
        ),
        reverse=True,
    )
    return ranked[0]


def nearest_node(target: int, nodes: Dict[int, Tuple[int, str]], wanted: str) -> Tuple[Optional[int], Optional[int]]:
    candidates = [(pos, idx) for idx, (pos, kinds) in nodes.items() if wanted in kinds.split(",")]
    if not candidates:
        return (None, None)
    pos, _idx = min(candidates, key=lambda item: (abs(item[0] - target), item[0]))
    return (pos, abs(pos - target))


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
    bam.close()
    return counts


def start_terminal_type(strand: str) -> str:
    return "TSS" if strand == "+" else "TES"


def end_terminal_type(strand: str) -> str:
    return "TES" if strand == "+" else "TSS"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-gtf", required=True)
    ap.add_argument("--query-gtf", required=True)
    ap.add_argument("--denovo-nodes", required=True)
    ap.add_argument("--transcript-nodes", required=True)
    ap.add_argument("--path-table", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--tol", type=int, default=10)
    args = ap.parse_args()

    refs = load_gtf(args.reference_gtf)
    queries = load_gtf(args.query_gtf)
    denovo_nodes = parse_nodes(args.denovo_nodes)
    transcript_nodes = parse_nodes(args.transcript_nodes)
    paths_by_chain = parse_path_table(args.path_table)

    queries_by_chain: Dict[Tuple[Tuple[int, int], ...], List[Tx]] = defaultdict(list)
    for tx in queries:
        queries_by_chain[normalized_chain(tx.junctions)].append(tx)

    out_path = Path(f"{args.out_prefix}.tsv")
    with out_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "ref_id",
                "ref_start",
                "ref_end",
                "junctions",
                "raw_exact_boundary_count",
                "raw_top_exact_boundary",
                "raw_has_ref_boundary",
                "denovo_has_ref_tss_node",
                "denovo_nearest_tss",
                "denovo_nearest_tss_delta",
                "denovo_has_ref_tes_node",
                "denovo_nearest_tes",
                "denovo_nearest_tes_delta",
                "transcript_has_ref_tss_node",
                "transcript_nearest_tss",
                "transcript_nearest_tss_delta",
                "transcript_has_ref_tes_node",
                "transcript_nearest_tes",
                "transcript_nearest_tes_delta",
                "best_path_source",
                "best_path_support",
                "best_path_boundary",
                "best_emitted_boundary",
                "best_emitted_fl_count",
            ]
        )

        for ref in refs:
            chain = normalized_chain(ref.junctions)
            raw_counts = exact_chain_boundary_counts(
                args.bam, ref.chrom, ref.start, ref.end, chain, args.tol
            )
            raw_top = raw_counts.most_common(1)
            raw_top_boundary = ""
            raw_exact_boundary_count = 0
            if raw_top:
                raw_top_boundary = f"{raw_top[0][0][0]}-{raw_top[0][0][1]}"
                raw_exact_boundary_count = raw_top[0][1]

            start_type = start_terminal_type(ref.strand)
            end_type = end_terminal_type(ref.strand)
            denovo_tss_pos, denovo_tss_delta = nearest_node(ref.start, denovo_nodes, start_type)
            denovo_tes_pos, denovo_tes_delta = nearest_node(ref.end, denovo_nodes, end_type)
            transcript_tss_pos, transcript_tss_delta = nearest_node(ref.start, transcript_nodes, start_type)
            transcript_tes_pos, transcript_tes_delta = nearest_node(ref.end, transcript_nodes, end_type)

            path_rows = paths_by_chain.get(chain, [])
            best_path = best_path_candidate(path_rows)
            best_path_boundary = ""
            best_path_source = ""
            best_path_support = 0
            if best_path:
                nodes = parse_node_list(best_path["nodes"])
                if nodes:
                    first = nodes[0]
                    last = nodes[-1]
                    a, b = sorted((first, last))
                    best_path_boundary = f"{a}-{b}"
                best_path_source = best_path["source"]
                best_path_support = int(float(best_path["count"]))

            emitted_same_chain = sorted(
                queries_by_chain.get(chain, []),
                key=lambda tx: int(tx.attrs.get("FL_count", "0") or "0"),
                reverse=True,
            )
            best_emitted_boundary = ""
            best_emitted_fl = 0
            if emitted_same_chain:
                tx = emitted_same_chain[0]
                best_emitted_boundary = f"{tx.start}-{tx.end}"
                best_emitted_fl = int(tx.attrs.get("FL_count", "0") or "0")

            writer.writerow(
                [
                    ref.tx_id,
                    ref.start,
                    ref.end,
                    len(chain),
                    raw_exact_boundary_count,
                    raw_top_boundary,
                    int((ref.start, ref.end) in raw_counts),
                    int(any(pos == ref.start and start_type in kinds.split(",") for pos, kinds in denovo_nodes.values())),
                    denovo_tss_pos or "",
                    denovo_tss_delta if denovo_tss_delta is not None else "",
                    int(any(pos == ref.end and end_type in kinds.split(",") for pos, kinds in denovo_nodes.values())),
                    denovo_tes_pos or "",
                    denovo_tes_delta if denovo_tes_delta is not None else "",
                    int(any(pos == ref.start and start_type in kinds.split(",") for pos, kinds in transcript_nodes.values())),
                    transcript_tss_pos or "",
                    transcript_tss_delta if transcript_tss_delta is not None else "",
                    int(any(pos == ref.end and end_type in kinds.split(",") for pos, kinds in transcript_nodes.values())),
                    transcript_tes_pos or "",
                    transcript_tes_delta if transcript_tes_delta is not None else "",
                    best_path_source,
                    best_path_support,
                    best_path_boundary,
                    best_emitted_boundary,
                    best_emitted_fl,
                ]
            )


if __name__ == "__main__":
    main()
