#!/usr/bin/env python3
"""Join rustle snapshot JSONL with normalized trace JSONL for quick 1:1 investigation.

This is intentionally lightweight: it does not attempt deep semantic equivalence.
It helps answer "what did rustle build (graph/transfrags/transcripts) at this locus"
alongside "what did StringTie log as decisions near this locus".

Inputs:
  - rustle snapshot JSONL from `--snapshot-jsonl`
  - normalized trace JSONL from `scripts/trace_decision_extract.py --format jsonl`

Example:
  python3 scripts/compare_snapshot_trace_json.py \
    --locus 17159009-17178368 \
    --snap /tmp/rustle_mini_panel/rustle_snapshots.jsonl \
    --trace /tmp/rustle_mini_panel/stringtie_trace_norm.jsonl \
    -o /tmp/rustle_mini_panel/join_report.json
"""

from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path
from typing import Any, Iterable, Iterator, Optional, Tuple


def parse_span(text: str) -> Optional[Tuple[int, int]]:
    if not text:
        return None
    if "-" not in text:
        return None
    a, b = text.split("-", 1)
    try:
        start = int(a.strip())
        end = int(b.strip())
    except ValueError:
        return None
    if start > end:
        start, end = end, start
    return start, end


def overlaps(a: Tuple[int, int], b: Tuple[int, int]) -> bool:
    return a[0] <= b[1] and a[1] >= b[0]


def read_jsonl(path: Path) -> Iterator[dict[str, Any]]:
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            yield json.loads(line)


def row_span_from_trace_row(row: dict[str, Any]) -> Optional[Tuple[int, int]]:
    for a, b in (
        ("span_start", "span_end"),
        ("boundary_start", "boundary_end"),
        ("target_start", "target_end"),
        ("longstart", "longend"),
    ):
        sa = row.get(a) or ""
        sb = row.get(b) or ""
        try:
            if sa and sb:
                return int(sa), int(sb)
        except ValueError:
            pass
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--locus", required=True, help="start-end (1-based inclusive intent)")
    ap.add_argument("--snap", type=Path, required=True, help="rustle snapshot jsonl")
    ap.add_argument("--trace", type=Path, required=True, help="normalized trace jsonl")
    ap.add_argument("-o", "--output", type=Path, help="output JSON report (default: stdout)")
    args = ap.parse_args()

    locus = parse_span(args.locus)
    if locus is None:
        raise SystemExit(f"invalid --locus: {args.locus}")

    # Rustle snapshots: keep any bundle overlapping locus.
    snap_rows = []
    for rec in read_jsonl(args.snap):
        b = rec.get("bundle") or {}
        try:
            bspan = (int(b.get("start")), int(b.get("end")))
        except Exception:
            continue
        if overlaps(locus, bspan):
            snap_rows.append(rec)

    # StringTie trace rows: keep those with a usable span overlapping locus.
    trace_rows = []
    for rec in read_jsonl(args.trace):
        span = row_span_from_trace_row(rec)
        if span is None:
            continue
        if overlaps(locus, span):
            trace_rows.append(rec)

    # Summaries for scanning.
    snap_stage_counts = collections.Counter((r.get("stage") or "") for r in snap_rows)
    trace_stage_event = collections.Counter(
        (r.get("stage") or "", r.get("event") or "", r.get("decision") or "")
        for r in trace_rows
    )

    out = {
        "locus": {"start": locus[0], "end": locus[1]},
        "rustle_snapshot": {
            "rows": snap_rows,
            "stage_counts": dict(snap_stage_counts),
        },
        "stringtie_trace": {
            "rows": trace_rows,
            "stage_event_decision_counts": [
                {"stage": k[0], "event": k[1], "decision": k[2], "count": v}
                for k, v in trace_stage_event.most_common(200)
            ],
        },
    }

    text = json.dumps(out, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(text + "\n")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

