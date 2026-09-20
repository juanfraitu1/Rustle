#!/usr/bin/env python3
"""Summarize bench needy-locus runs: gffcompare -R -Q stats + trace summary lines.

Looks for files:
  {prefix}_run_STRG_*.gtf  -> paired with
  {prefix}_run_STRG_*_cmp.stats
  {prefix}_run_STRG_*.trace.txt
  {prefix}_run_STRG_*.needy_bundle_summary.tsv (auto from rustle needy-locus mode; optional)

Example (repo root):
  python3 scripts/summarize_needy_locus_runs.py --prefix bench/ggo19_needy_top5 \\
    --out-tsv bench/ggo19_needy_top5_batch_summary.tsv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def parse_stats(path: Path) -> dict[str, str]:
    text = path.read_text()
    out: dict[str, str] = {}
    m = re.search(r"#\s+Query mRNAs\s*:\s*(\d+)", text)
    if m:
        out["query_mrnas"] = m.group(1)
    m = re.search(r"#\s+Reference mRNAs\s*:\s*(\d+)", text)
    if m:
        out["ref_mrnas"] = m.group(1)
    m = re.search(r"Transcript level:\s*([\d.]+\s*\|\s*[\d.]+)", text)
    if m:
        out["transcript_sn_pr"] = m.group(1).strip()
    m = re.search(r"Matching transcripts:\s*(\d+)", text)
    if m:
        out["matching_transcripts"] = m.group(1)
    return out


def parse_trace_summary(path: Path) -> str:
    lines = path.read_text().splitlines()
    data_rows: list[str] = []
    summary_bits: list[str] = []
    for line in lines:
        if line.startswith("#"):
            s = line[1:].lstrip()
            if re.match(r"^\s+\w[\w_]*:.+", s) and "blocked /" not in line:
                summary_bits.append(s.strip())
            continue
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) >= 2 and not parts[0].startswith("#"):
            data_rows.append(f"{parts[0]}:{parts[1]}")
    if data_rows:
        tail = "; ".join(data_rows[:6])
        if len(data_rows) > 6:
            tail += f" … (+{len(data_rows) - 6} refs)"
        return tail
    # No per-ref rows: look for "0 blocked" in strand lines
    for line in lines:
        if "blocked /" in line and "0 blocked" in line:
            return "no_ref_blockers"
    if summary_bits:
        return "; ".join(summary_bits[:5])
    return "no_ref_blockers"


def parse_needy_bundle_summary(path: Path) -> str:
    """Count key stages in parity-stage TSV (bundle_context + seed_flow rows)."""
    if not path.is_file():
        return ""
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        return f"empty_or_header_only lines={len(lines)}"
    bundle_ctx = 0
    seed_flow = 0
    for line in lines[1:]:
        stage = line.split("\t", 1)[0]
        if stage == "bundle_context":
            bundle_ctx += 1
        elif stage == "seed_flow":
            seed_flow += 1
    return f"bundle_context_rows={bundle_ctx} seed_flow_rows={seed_flow}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--prefix",
        type=Path,
        required=True,
        help="Stem like bench/ggo19_needy_top5 (expects {prefix}_run_STRG_* files)",
    )
    ap.add_argument(
        "--out-tsv",
        type=Path,
        default=None,
        help="Write TSV table (default: stdout only)",
    )
    args = ap.parse_args()
    stem = args.prefix.name
    parent = args.prefix.parent
    pattern = f"{stem}_run_STRG_*.trace.txt"
    traces = sorted(parent.glob(pattern))
    rows: list[list[str]] = []
    for tr in traces:
        tag = tr.name.replace(f"{stem}_run_", "").replace(".trace.txt", "")
        stats = parent / f"{stem}_run_{tag}_cmp.stats"
        if not stats.is_file():
            continue
        st = parse_stats(stats)
        ts = parse_trace_summary(tr)
        bundle_tsv = tr.with_name(tr.name.replace(".trace.txt", ".needy_bundle_summary.tsv"))
        bundle_note = parse_needy_bundle_summary(bundle_tsv)
        rows.append(
            [
                tag,
                st.get("ref_mrnas", ""),
                st.get("query_mrnas", ""),
                st.get("matching_transcripts", ""),
                st.get("transcript_sn_pr", ""),
                ts,
                bundle_note,
            ]
        )

    header = [
        "locus_run",
        "ref_mrnas",
        "query_mrnas",
        "matching_transcripts",
        "transcript_sn_pr",
        "trace_notes",
        "needy_bundle_summary",
    ]
    lines_out = ["\t".join(header)]
    for r in rows:
        lines_out.append("\t".join(r))
    text = "\n".join(lines_out) + "\n"
    if args.out_tsv:
        args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
        args.out_tsv.write_text(text)
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
