#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit terminal junction variants across native bundle trace, native GTF, "
            "rlink GTF, and legacy print_predcluster final-fate trace."
        )
    )
    parser.add_argument("--native-trace", required=True)
    parser.add_argument("--native-gtf", required=True)
    parser.add_argument("--rlink-gtf", required=True)
    parser.add_argument("--trace-log", required=True)
    parser.add_argument("--terminal-donor", type=int, required=True)
    parser.add_argument("--output-stem", required=True)
    return parser.parse_args()


def parse_attrs(attr_text: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for part in attr_text.strip().split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, value = part.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def last_junction(chain: str) -> str:
    if not chain:
        return ""
    parts = [part for part in chain.split(",") if part]
    return parts[-1] if parts else ""


def chain_has_terminal_donor(chain: str, donor: int) -> bool:
    lj = last_junction(chain)
    return lj.startswith(f"{donor}-")


def load_native_trace(path: str, donor: int) -> List[dict]:
    by_chain: Dict[str, dict] = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chain = row.get("emit_chain", "")
            if not chain_has_terminal_donor(chain, donor):
                continue
            entry = by_chain.setdefault(
                chain,
                {
                    "last_junction": last_junction(chain),
                    "family_input_rows": 0,
                    "family_input_support_sum": 0,
                    "mode_member_rows": 0,
                    "mode_member_support_sum": 0,
                    "selected_rows": 0,
                    "selected_support_sum": 0,
                    "selected_mode_total_reads": 0,
                    "selected_boundary_chain": "",
                    "selected_reason": "",
                },
            )
            event = row.get("event", "")
            decision = row.get("decision", "")
            support = int(row.get("support") or 0)
            mode_total_reads = int(row.get("mode_total_reads") or 0)
            if event == "family_input":
                entry["family_input_rows"] += 1
                entry["family_input_support_sum"] += support
            elif event == "mode_eval":
                if decision == "selected":
                    entry["selected_rows"] += 1
                    entry["selected_support_sum"] += support
                    entry["selected_mode_total_reads"] = max(
                        entry["selected_mode_total_reads"], mode_total_reads
                    )
                    entry["selected_boundary_chain"] = row.get("boundary_chain", "")
                    entry["selected_reason"] = row.get("reason", "")
                elif decision == "mode_member":
                    entry["mode_member_rows"] += 1
                    entry["mode_member_support_sum"] += support
    rows = list(by_chain.values())
    rows.sort(key=lambda row: row["last_junction"])
    return rows


def load_gtf(path: str, donor: int) -> List[dict]:
    rows: List[dict] = []
    with open(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attrs = parse_attrs(fields[8])
            chain = attrs.get("junction_chain", "")
            if not chain_has_terminal_donor(chain, donor):
                continue
            rows.append(
                {
                    "transcript_id": attrs.get("transcript_id", ""),
                    "gene_id": attrs.get("gene_id", ""),
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "num_exons": int(attrs.get("num_exons", "0") or 0),
                    "fl_count": int(attrs.get("FL_count", "0") or 0),
                    "min_junction_support": int(
                        attrs.get("min_junction_support", "0") or 0
                    ),
                    "flags": attrs.get("flags", ""),
                    "junction_chain": chain,
                    "last_junction": last_junction(chain),
                }
            )
    rows.sort(key=lambda row: (row["last_junction"], row["start"], row["end"], row["transcript_id"]))
    return rows


def load_trace_final_fates(path: str) -> Dict[Tuple[int, int], List[dict]]:
    pattern = re.compile(
        r"print_predcluster: FINAL_FATE n=(?P<n>\d+) "
        r"(?P<start>\d+)-(?P<end>\d+) cov=(?P<cov>[0-9.]+) "
        r"strand=(?P<strand>[+-]) exons=(?P<exons>\d+) .* fate=(?P<fate>[A-Z]+)"
    )
    by_span: Dict[Tuple[int, int], List[dict]] = defaultdict(list)
    with open(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            match = pattern.search(line)
            if not match:
                continue
            start = int(match.group("start"))
            end = int(match.group("end"))
            by_span[(start, end)].append(
                {
                    "line": lineno,
                    "pred_n": int(match.group("n")),
                    "cov": float(match.group("cov")),
                    "strand": match.group("strand"),
                    "exons": int(match.group("exons")),
                    "fate": match.group("fate"),
                }
            )
    return by_span


def write_tsv(path: Path, fieldnames: List[str], rows: List[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()

    native_trace_rows = load_native_trace(args.native_trace, args.terminal_donor)
    native_gtf_rows = load_gtf(args.native_gtf, args.terminal_donor)
    rlink_gtf_rows = load_gtf(args.rlink_gtf, args.terminal_donor)
    trace_fates = load_trace_final_fates(args.trace_log)

    trace_path = Path(args.native_trace)
    native_path = Path(args.native_gtf)
    rlink_path = Path(args.rlink_gtf)

    native_trace_out = Path(f"{args.output_stem}.native_trace.tsv")
    native_gtf_out = Path(f"{args.output_stem}.native_selected.tsv")
    rlink_gtf_out = Path(f"{args.output_stem}.rlink_selected.tsv")
    summary_out = Path(f"{args.output_stem}.summary.tsv")

    write_tsv(
        native_trace_out,
        [
            "last_junction",
            "family_input_rows",
            "family_input_support_sum",
            "mode_member_rows",
            "mode_member_support_sum",
            "selected_rows",
            "selected_support_sum",
            "selected_mode_total_reads",
            "selected_boundary_chain",
            "selected_reason",
        ],
        native_trace_rows,
    )

    native_selected_rows: List[dict] = []
    for row in native_gtf_rows:
        row = dict(row)
        matches = trace_fates.get((row["start"], row["end"]), [])
        row["trace_final_fates"] = ",".join(match["fate"] for match in matches)
        row["trace_final_lines"] = ",".join(str(match["line"]) for match in matches)
        native_selected_rows.append(row)
    write_tsv(
        native_gtf_out,
        [
            "transcript_id",
            "gene_id",
            "start",
            "end",
            "num_exons",
            "fl_count",
            "min_junction_support",
            "flags",
            "junction_chain",
            "last_junction",
            "trace_final_fates",
            "trace_final_lines",
        ],
        native_selected_rows,
    )

    rlink_selected_rows: List[dict] = []
    for row in rlink_gtf_rows:
        row = dict(row)
        matches = trace_fates.get((row["start"], row["end"]), [])
        row["trace_final_fates"] = ",".join(match["fate"] for match in matches)
        row["trace_final_lines"] = ",".join(str(match["line"]) for match in matches)
        rlink_selected_rows.append(row)
    write_tsv(
        rlink_gtf_out,
        [
            "transcript_id",
            "gene_id",
            "start",
            "end",
            "num_exons",
            "fl_count",
            "min_junction_support",
            "flags",
            "junction_chain",
            "last_junction",
            "trace_final_fates",
            "trace_final_lines",
        ],
        rlink_selected_rows,
    )

    summary_rows: List[dict] = []
    all_last = sorted(
        {
            row["last_junction"] for row in native_trace_rows
        }
        | {row["last_junction"] for row in native_gtf_rows}
        | {row["last_junction"] for row in rlink_gtf_rows}
    )
    native_trace_by_last = {row["last_junction"]: row for row in native_trace_rows}
    native_gtf_by_last = defaultdict(list)
    for row in native_selected_rows:
        native_gtf_by_last[row["last_junction"]].append(row)
    rlink_gtf_by_last = defaultdict(list)
    for row in rlink_selected_rows:
        rlink_gtf_by_last[row["last_junction"]].append(row)

    for last in all_last:
        nt = native_trace_by_last.get(last, {})
        native_selected_rows_by_last = native_gtf_by_last.get(last, [])
        rlink_selected_rows_by_last = rlink_gtf_by_last.get(last, [])
        summary_rows.append(
            {
                "last_junction": last,
                "native_family_input_rows": nt.get("family_input_rows", 0),
                "native_mode_member_rows": nt.get("mode_member_rows", 0),
                "native_selected_rows": nt.get("selected_rows", 0),
                "native_selected_ids": ",".join(
                    row["transcript_id"] for row in native_selected_rows_by_last
                ),
                "native_selected_spans": ",".join(
                    f"{row['start']}-{row['end']}" for row in native_selected_rows_by_last
                ),
                "rlink_selected_ids": ",".join(
                    row["transcript_id"] for row in rlink_selected_rows_by_last
                ),
                "rlink_selected_spans": ",".join(
                    f"{row['start']}-{row['end']}" for row in rlink_selected_rows_by_last
                ),
                "rlink_trace_final_fates": ",".join(
                    sorted(
                        {
                            fate
                            for row in rlink_selected_rows_by_last
                            for fate in row["trace_final_fates"].split(",")
                            if fate
                        }
                    )
                ),
            }
        )

    write_tsv(
        summary_out,
        [
            "last_junction",
            "native_family_input_rows",
            "native_mode_member_rows",
            "native_selected_rows",
            "native_selected_ids",
            "native_selected_spans",
            "rlink_selected_ids",
            "rlink_selected_spans",
            "rlink_trace_final_fates",
        ],
        summary_rows,
    )

    print(f"wrote {native_trace_out}")
    print(f"wrote {native_gtf_out}")
    print(f"wrote {rlink_gtf_out}")
    print(f"wrote {summary_out}")
    print(
        "summary: "
        f"{trace_path.name} native_trace_rows={len(native_trace_rows)} "
        f"{native_path.name} native_selected={len(native_gtf_rows)} "
        f"{rlink_path.name} rlink_selected={len(rlink_gtf_rows)}"
    )


if __name__ == "__main__":
    main()
