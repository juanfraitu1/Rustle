#!/usr/bin/env python3
"""Normalize old/current assembly trace formats into one TSV schema.

Supported inputs:
- legacy StringTie/rlink-style logs such as `trace_GGO_19_deep.log`
- old Rust port `TRACE_*` stderr logs
- current `RUSTLE_BNODE_TRACE_TSV` outputs
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Iterator, Optional


OUTPUT_FIELDS = [
    "source_kind",
    "line_no",
    "stage",
    "event",
    "decision",
    "reason",
    "entity_id",
    "parent_id",
    "family_key",
    "bundle_chrom",
    "bundle_strand",
    "span_start",
    "span_end",
    "boundary_start",
    "boundary_end",
    "target_start",
    "target_end",
    "support",
    "mode_total",
    "coverage",
    "abundance",
    "guide",
    "longread",
    "longstart",
    "longend",
    "hardstart",
    "hardend",
    "source_state",
    "sink_state",
    "node_count",
    "exon_count",
    "node_coords",
    "path_coords",
    "exons",
    "raw",
]

KEYVAL_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=")
SPAN_RE = re.compile(r"(\d+)-(\d+)")
COORD_PAIR_RE = re.compile(r"(\d+)-(\d+)")


@dataclass
class Row:
    source_kind: str = ""
    line_no: str = ""
    stage: str = ""
    event: str = ""
    decision: str = ""
    reason: str = ""
    entity_id: str = ""
    parent_id: str = ""
    family_key: str = ""
    bundle_chrom: str = ""
    bundle_strand: str = ""
    span_start: str = ""
    span_end: str = ""
    boundary_start: str = ""
    boundary_end: str = ""
    target_start: str = ""
    target_end: str = ""
    support: str = ""
    mode_total: str = ""
    coverage: str = ""
    abundance: str = ""
    guide: str = ""
    longread: str = ""
    longstart: str = ""
    longend: str = ""
    hardstart: str = ""
    hardend: str = ""
    source_state: str = ""
    sink_state: str = ""
    node_count: str = ""
    exon_count: str = ""
    node_coords: str = ""
    path_coords: str = ""
    exons: str = ""
    raw: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("-o", "--output", type=Path)
    parser.add_argument(
        "--format",
        choices=["tsv", "jsonl"],
        default="tsv",
        help="Output format. `tsv` matches older tooling; `jsonl` is better for 1:1 diffing.",
    )
    parser.add_argument(
        "--mode",
        choices=["auto", "legacy_log", "rust_trace", "current_tsv"],
        default="auto",
    )
    parser.add_argument(
        "--locus",
        help="Only keep rows overlapping start-end.",
    )
    parser.add_argument(
        "--event-regex",
        help="Only keep rows whose event or stage matches this regex.",
    )
    return parser.parse_args()


def detect_mode(path: Path) -> str:
    with path.open() as handle:
        for line in handle:
            text = line.rstrip("\n")
            if not text:
                continue
            if text.startswith("event\tbundle_chrom\tbundle_strand\tfamily_key"):
                return "current_tsv"
            if text.startswith("[TRACE_"):
                return "rust_trace"
            return "legacy_log"
    return "legacy_log"


def parse_keyvals_loose(text: str) -> tuple[str, dict[str, str]]:
    matches = list(KEYVAL_RE.finditer(text))
    if not matches:
        return text.strip(), {}
    prefix = text[: matches[0].start()].strip()
    out: dict[str, str] = {}
    for idx, match in enumerate(matches):
        key = match.group(1)
        start = match.end()
        end = matches[idx + 1].start() if idx + 1 < len(matches) else len(text)
        value = text[start:end].strip()
        out[key] = value
    return prefix, out


def normalize_span_text(text: str) -> tuple[str, str]:
    match = SPAN_RE.search(text)
    if not match:
        return "", ""
    return match.group(1), match.group(2)


def count_csvish(text: str) -> str:
    stripped = text.strip()
    if not stripped:
        return ""
    return str(len([part for part in stripped.split(",") if part]))


def span_from_coords(text: str) -> tuple[str, str]:
    coords = [(int(a), int(b)) for a, b in COORD_PAIR_RE.findall(text)]
    if not coords:
        return "", ""
    lo = min(a for a, _ in coords)
    hi = max(b for _, b in coords)
    return str(lo), str(hi)


def safe_pair(start: str, end: str) -> tuple[Optional[int], Optional[int]]:
    try:
        return int(start), int(end)
    except (TypeError, ValueError):
        return None, None


def row_span(row: Row) -> tuple[Optional[int], Optional[int]]:
    for start_key, end_key in (
        ("span_start", "span_end"),
        ("boundary_start", "boundary_end"),
        ("target_start", "target_end"),
        ("longstart", "longend"),
    ):
        start = getattr(row, start_key)
        end = getattr(row, end_key)
        if start and end:
            pair = safe_pair(start, end)
            if pair != (None, None):
                return pair
    for text in (row.node_coords, row.path_coords, row.exons):
        start, end = span_from_coords(text)
        if start and end:
            return int(start), int(end)
    return None, None


def locus_overlaps(row: Row, locus: tuple[int, int]) -> bool:
    start, end = row_span(row)
    if start is None or end is None:
        return False
    lo, hi = locus
    return start <= hi and end >= lo


def parse_legacy_line(line: str, line_no: int) -> Optional[Row]:
    text = line.rstrip("\n")
    if not text.startswith("--- "):
        return None
    body = text[4:]
    if ": " not in body:
        return None
    stage, rest = body.split(": ", 1)
    parts = rest.split(None, 1)
    event = parts[0]
    payload = parts[1] if len(parts) > 1 else ""
    payload = payload.replace("path:", "path_coords=")
    payload = payload.replace("nodes= ", "nodes=")
    prefix, fields = parse_keyvals_loose(payload)

    row = Row(source_kind="legacy_log", line_no=str(line_no), stage=stage, event=event, raw=text)
    row.reason = prefix

    if event == "PRED_FATE":
        n_match = re.search(r"\bn=(\d+)\s+(\d+-\d+)", payload)
        if n_match:
            row.entity_id = n_match.group(1)
            row.span_start, row.span_end = normalize_span_text(n_match.group(2))
        row.coverage = fields.get("cov", "")
        row.bundle_strand = fields.get("strand", "")
        row.exon_count = fields.get("exons", "")
        row.guide = fields.get("guide", "")
        row.decision = fields.get("fate", "")
        return row

    row.entity_id = fields.get("t") or fields.get("n") or fields.get("f") or ""
    row.parent_id = fields.get("k", "")
    row.coverage = fields.get("cov", "") or fields.get("flux", "")
    row.abundance = fields.get("abund", "")
    row.support = fields.get("keeptrf", "")
    row.guide = fields.get("guide", "")
    row.longread = fields.get("longread", "")
    row.longstart = fields.get("longstart", "")
    row.longend = fields.get("longend", "")
    row.hardstart = fields.get("hardstart", "")
    row.hardend = fields.get("hardend", "")
    row.source_state = fields.get("source", "")
    row.sink_state = fields.get("sink", "")
    row.node_count = fields.get("nodes", "") or fields.get("path_len", "")
    row.node_coords = fields.get("node_coords", "") or fields.get("node", "")
    row.path_coords = fields.get("path_coords", "")
    row.decision = fields.get("result", "") or fields.get("used_direct", "")
    if event == "SRCSINK_DECISION" and row.node_coords:
        row.span_start, row.span_end = span_from_coords(row.node_coords)
    elif row.path_coords:
        row.span_start, row.span_end = span_from_coords(row.path_coords)
    elif row.node_coords:
        row.span_start, row.span_end = span_from_coords(row.node_coords)
    elif row.longstart and row.longend:
        row.span_start, row.span_end = row.longstart, row.longend
    return row


def parse_rust_trace_line(line: str, line_no: int) -> Optional[Row]:
    text = line.rstrip("\n")
    if not text.startswith("[TRACE_"):
        return None
    close = text.find("]")
    if close == -1:
        return None
    event = text[1:close]
    payload = text[close + 1 :].strip()
    prefix, fields = parse_keyvals_loose(payload)
    row = Row(source_kind="rust_trace", line_no=str(line_no), stage="path_extract", event=event, raw=text)
    row.reason = prefix
    row.entity_id = fields.get("t") or fields.get("idx") or fields.get("tf") or ""
    row.parent_id = fields.get("k") or fields.get("out_idx") or ""
    row.coverage = fields.get("cov", "") or fields.get("tx_cov", "") or fields.get("keep_cov", "")
    row.abundance = fields.get("abund", "") or fields.get("abundsum", "")
    row.support = fields.get("keeptrf", "")
    row.guide = fields.get("guide", "")
    row.longstart = fields.get("longstart", "")
    row.longend = fields.get("longend", "")
    row.node_count = fields.get("nodes", "") or fields.get("tf_nodes", "")
    row.exon_count = fields.get("exons", "")
    row.node_coords = fields.get("coords", "") or fields.get("nodes", "")
    row.path_coords = fields.get("path_coords", "")
    row.exons = fields.get("exons", "")
    row.decision = fields.get("mode", "") or fields.get("reject", "")
    row.reason = "|".join(filter(None, [row.reason, fields.get("bundle", "")]))
    if "span" in fields:
        row.span_start, row.span_end = normalize_span_text(fields["span"])
    elif row.node_coords:
        row.span_start, row.span_end = span_from_coords(row.node_coords)
    elif row.exons:
        row.span_start, row.span_end = span_from_coords(row.exons)
    elif row.longstart and row.longend:
        row.span_start, row.span_end = row.longstart, row.longend
    return row


def parse_current_tsv(path: Path) -> Iterator[Row]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for idx, rec in enumerate(reader, 2):
            emit_nodes = rec.get("emit_path_nodes", "")
            emit_chain = rec.get("emit_chain", "")
            row = Row(
                source_kind="current_bnode",
                line_no=str(idx),
                stage=rec.get("event", ""),
                event=rec.get("event", ""),
                decision=rec.get("decision", ""),
                reason=rec.get("reason", ""),
                entity_id=rec.get("path_idx", "") or rec.get("selected_path_idx", ""),
                parent_id=rec.get("selected_path_idx", ""),
                family_key=rec.get("family_key", ""),
                bundle_chrom=rec.get("bundle_chrom", ""),
                bundle_strand=rec.get("bundle_strand", ""),
                span_start=rec.get("boundary_start", ""),
                span_end=rec.get("boundary_end", ""),
                boundary_start=rec.get("boundary_start", ""),
                boundary_end=rec.get("boundary_end", ""),
                target_start=rec.get("target_start", ""),
                target_end=rec.get("target_end", ""),
                support=rec.get("support", ""),
                mode_total=rec.get("mode_total_reads", ""),
                node_count=count_csvish(emit_nodes),
                node_coords=emit_nodes,
                path_coords=rec.get("boundary_chain", ""),
                exons=emit_chain,
                raw="\t".join(rec.get(field, "") for field in reader.fieldnames or []),
            )
            yield row


def parse_text_lines(path: Path, mode: str) -> Iterator[Row]:
    parser = parse_legacy_line if mode == "legacy_log" else parse_rust_trace_line
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            row = parser(line, line_no)
            if row is not None:
                yield row


def normalize_rows(path: Path, mode: str) -> Iterator[Row]:
    if mode == "current_tsv":
        yield from parse_current_tsv(path)
        return
    yield from parse_text_lines(path, mode)


def main() -> int:
    args = parse_args()
    mode = detect_mode(args.input) if args.mode == "auto" else args.mode
    locus = None
    if args.locus:
        start, end = normalize_span_text(args.locus)
        if not start or not end:
            raise SystemExit(f"invalid --locus: {args.locus}")
        locus = (int(start), int(end))
    event_re = re.compile(args.event_regex) if args.event_regex else None

    out_handle = args.output.open("w", newline="") if args.output else sys.stdout
    close_out = args.output is not None
    try:
        writer = None
        for row in normalize_rows(args.input, mode):
            if locus and not locus_overlaps(row, locus):
                continue
            if event_re and not (
                event_re.search(row.event or "") or event_re.search(row.stage or "")
            ):
                continue
            rec = asdict(row)
            if args.format == "tsv":
                if writer is None:
                    writer = csv.DictWriter(out_handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
                    writer.writeheader()
                writer.writerow(rec)
            else:
                out_handle.write(json.dumps(rec, sort_keys=True) + "\n")
    finally:
        if close_out:
            out_handle.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
