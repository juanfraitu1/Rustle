#!/usr/bin/env python3
"""Explain one bundle using staged backward-process layer tables.

This script is intentionally analysis-first. It does not try to rebuild the old
graph objects exactly. Instead it lifts the bundle summary, reference GTF,
legacy deep-log context, and normalized old Rust traces into one staged package:

- `rawgraph.nodes.tsv`: node-level source/sink and hard-boundary state
- `rawgraph.transfrags.tsv`: raw transfrags from `TRACE_RAWTRF`
- `checktrf.seeds.tsv`: deferred/rescue seeds from `TRACE_CHECKTRF_ADD`
- `checktrf.edges.tsv`: graph-like links from seeds to keep/output states
- `checktrf.outputs.tsv`: output-side summaries from `TRACE_CHECKTRF_*`
- `transcripts.tsv`: final transcript table copied from the provided old output
- `transcript_support.tsv`: final transcripts linked back to checktrf evidence
- `reference.tsv`: reference transcripts overlapping the focused current bundle
- `bundle_context.tsv`: current-bundle and old-log bundle context
- `stages.tsv`: raw stage counts
- `layers.tsv`: typed backward-process layers derived from the staged tables
- `explain.md`: human-readable interpretation
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


BUNDLE_BEGIN_RE = re.compile(
    r">bundle\s+(?P<chrom>\S+):(?P<start>\d+)-(?P<end>\d+)\s+"
    r"\[(?P<alignments>\d+)\s+alignments\s+\((?P<distinct>\d+)\s+distinct\),\s+"
    r"(?P<junctions>\d+)\s+junctions,\s+(?P<guides>\d+)\s+guides\]\s+begins processing"
)
BUNDLE_DONE_RE = re.compile(
    r"\^bundle\s+(?P<chrom>\S+):(?P<start>\d+)-(?P<end>\d+)\s+done\s+"
    r"\((?P<transcripts>\d+)\s+processed potential transcripts\)"
)
KEYVAL_RE = re.compile(r"([A-Za-z_][A-Za-z0-9_]*)=")
BOOL_PAIR_RE = re.compile(r"(?P<a>true|false|0|1)/(?P<b>true|false|0|1)")
COORD_RE = re.compile(r"(\d+)-(\d+)")
READ_RE = re.compile(r"read\[(?P<read>\d+)\]")
GATE_RE = re.compile(
    r"gate\[(?P<gate_idx>\d+)\]\s+t=(?P<seed_id>\d+)\s+abund=(?P<abundance>[-0-9.]+)\s+guide=(?P<guide>\d+)\s+pass=(?P<passed>\d+)"
)
JDEC_RE = re.compile(
    r"jdec\[(?P<idx>\d+)\]\s+(?P<junction>\d+-\d+)\s+strand=(?P<strand>-?\d+)\s+nm=(?P<nm>[-0-9.]+)\s+mm=(?P<mm>[-0-9.]+)\s+"
    r"nreads=(?P<nreads>[-0-9.]+)\s+nreads_good=(?P<nreads_good>[-0-9.]+)\s+left=(?P<left>[-0-9.]+)\s+right=(?P<right>[-0-9.]+)\s+"
    r"guide=(?P<guide>\d+)\s+reason=(?P<reason>\S+)"
)
GOOD_JUNC_RE = re.compile(
    r"good_junc:\s+ACCEPTED\s+(?P<junction>\d+-\d+):(?P<strand>-?\d+)\s+nreads=(?P<nreads>[-0-9.]+)\s+left=(?P<left>[-0-9.]+)\s+right=(?P<right>[-0-9.]+)"
)
JUNC_COLOR_BREAK_RE = re.compile(
    r"JUNC_COLOR_BREAK\s+read\[(?P<read_idx>\d+)\]\s+i=(?P<junc_order>\d+)\s+"
    r"junc=(?P<junction>\d+-\d+)\s+strand_now=(?P<strand_now>-?\d+)\s+mm=(?P<mm>[-0-9.]+)\s+"
    r"changeleft=(?P<changeleft>-?\d+)\s+changeright=(?P<changeright>-?\d+)\s+"
    r"nreads=(?P<nreads>[-0-9.]+)\s+nreads_good=(?P<nreads_good>[-0-9.]+)\s+good_junc=(?P<good_junc>-?\d+)"
)
PRED_LINE_RE = re.compile(
    r"pred\[(?P<pred_idx>\d+)\]:(?P<start>\d+)-(?P<end>\d+)\s+"
    r"\(cov=(?P<cov>[-0-9.]+),\s+readcov=(?P<readcov>[-0-9.]+),\s+"
    r"strand=(?P<strand>\S+)\s+falseflag=(?P<falseflag>\d+)\)"
)
PRED_FATE_RE = re.compile(
    r"PRED_FATE\s+n=(?P<pred_idx>\d+)\s+(?P<start>\d+)-(?P<end>\d+)\s+"
    r"cov=(?P<cov>[-0-9.]+)\s+strand=(?P<strand>\S+)\s+exons=(?P<exons>\d+)\s+"
    r"guide=(?P<guide>\S+)\s+fate=(?P<fate>\S+)"
)
PRED_FATE_SUMMARY_RE = re.compile(
    r"PRED_FATE_SUMMARY\s+npred=(?P<npred>\d+)\s+kept=(?P<kept>\d+)\s+killed=(?P<killed>\d+)"
)
FILTER_KILLED_RE = re.compile(
    r"FILTER\s+pred\[(?P<killed_idx>\d+)\]\s+(?P<killed_start>\d+)-(?P<killed_end>\d+)\s+"
    r"cov=(?P<killed_cov>[-0-9.]+)\s+strand=(?P<strand>\S+)\s+exons=(?P<killed_exons>\d+)\s+"
    r"KILLED_BY\s+pred\[(?P<killer_idx>\d+)\]\s+(?P<killer_start>\d+)-(?P<killer_end>\d+)\s+"
    r"cov=(?P<killer_cov>[-0-9.]+)\s+reason=(?P<filter_reason>\S+)"
)
SOURCE_PREFIX_RE = re.compile(r"^(?:---\s*)?(?P<source>[A-Za-z_]+(?:\[[0-9]+\])?):")
SPAN_KV_RE = re.compile(r"\b(?:bundle|bdata)=(?P<start>\d+)-(?P<end>\d+)")
ROUTE_BUNDLE_RE = re.compile(r"\bgrid=\d+\((?P<start>\d+)-(?P<end>\d+)\)")
NEW_BUNDLE_RE = re.compile(r"NEW_BUNDLE\s+bno=(?P<bno>\d+)\s+bid=(?P<bid>\d+)\s+(?P<start>\d+)-(?P<end>\d+)\s+cov=(?P<cov>[-0-9.]+)")
NEW_BNODE_RE = re.compile(r"NEW_BNODE\s+bid=(?P<bid>\d+)\s+(?P<start>\d+)-(?P<end>\d+)")
NEW_NODE_RE = re.compile(r"NEW_NODE\s+s=(?P<s>\d+)\s+g=(?P<g>\d+)\s+nodeid=(?P<nodeid>\d+)\s+(?P<start>\d+)-(?P<end>\d+)")
NODE_FINAL_END_RE = re.compile(r"NODE_FINAL_END\s+s=(?P<s>\d+)\s+g=(?P<g>\d+)\s+nodeid=(?P<nodeid>\d+)\s+(?P<start>\d+)-(?P<end>\d+)")
SRCSINK_RE = re.compile(
    r"SRCSINK_DECISION\s+node=(?P<nodeid>\d+)\((?P<start>\d+)-(?P<end>\d+)\)\s+"
    r"source=(?P<source_state>\S+)\s+addsource=(?P<addsource>[-0-9.]+)\s+hardstart=(?P<hardstart>\d+)\s+hassource=(?P<hassource>\d+)\s+"
    r"sink=(?P<sink_state>\S+)\s+addsink=(?P<addsink>[-0-9.]+)\s+hardend=(?P<hardend>\d+)\s+hassink=(?P<hassink>\d+)"
)
ROUTE_RE = re.compile(
    r"ROUTE_(?P<route_kind>STRANDED|UNSTRANDED)\s+pass=(?P<pass>\d+)\s+sno=(?P<sno>\d+)\s+"
    r"grid=(?P<grid_idx>\d+)\((?P<start>\d+)-(?P<end>\d+)\)"
)
JUNC_DELETE_RE = re.compile(
    r"(?P<event>E?JUNC_DELETE)\s+i=(?P<idx>\d+)\s+(?P<start>\d+)-(?P<end>\d+):(?P<strand>-?\d+)\s+nreads_good=(?P<nreads_good>[-0-9.]+)"
)
JUNC_DEMOTE_RE = re.compile(
    r"JUNC_DEMOTE\s+i=(?P<from_idx>\d+)\((?P<from_start>\d+)-(?P<from_end>\d+)\)\s+->\s+j=(?P<to_idx>\d+)\((?P<to_start>\d+)-(?P<to_end>\d+)\)\s+leftsup=(?P<left_from>[-0-9.]+)->(?P<left_to>[-0-9.]+)"
)
FLOW_STEP_RE = re.compile(
    r"---\s+(?P<phase>fwd_to_sink|back_to_source):\s+step\s+i=(?P<i>\d+)\s+->\s+maxc=(?P<maxc>-?\d+)\s+\((?P<start>\d+)-(?P<end>\d+)\)\s+maxcov=(?P<maxcov>[-0-9.]+)"
)
PARITY_SEED_PROC_RE = re.compile(
    r"PARITY_SEED_PROC\s+bundle=(?P<bundle_start>\d+)-(?P<bundle_end>\d+)\s+"
    r"idx=(?P<seed_idx>\d+)\s+span=(?P<span_start>\d+)-(?P<span_end>\d+)\s+"
    r"abund=(?P<abundance>[-0-9.]+)\s+nodes=(?P<nodes>\d+)\s+"
    r"back=(?P<back>\S+)\s+fwd=(?P<fwd>\S+)\s+flux=(?P<flux>[-0-9.]+)\s+"
    r"path=(?P<path>\[[^\]]*\])\s+outcome=(?P<outcome>\S+)"
)
PARITY_BUNDLE_RE = re.compile(
    r"PARITY_BUNDLE\s+sno=(?P<sno>\d+)\s+b=(?P<bundle_idx>\d+)\s+"
    r"start=(?P<bundle_start>\d+)\s+end=(?P<bundle_end>\d+)\s+"
    r"cov=(?P<cov>[-0-9.]+)\s+len=(?P<length>\d+)\s+bnodes=(?P<bnodes>\d+)"
)
PARITY_GRAPH_EVOLUTION_RE = re.compile(
    r"PARITY_GRAPH_EVOLUTION\s+s=(?P<sno>\d+)\s+g=(?P<graph_idx>\d+)\s+"
    r"op=(?P<event>\S+)\s+node=(?P<node_id>\d+)\s+start=(?P<span_start>\d+)\s+"
    r"end=(?P<span_end>\d+)\s+bid=(?P<bid>\d+)\s+parents=(?P<parents>\d+)"
)
PARITY_GRAPH_POST_RE = re.compile(
    r"PARITY_GRAPH_POST\s+sno=(?P<sno>\d+)\s+b=(?P<bundle_idx>\d+)\s+"
    r"nodes=(?P<nodes>\d+)\s+transfrags=(?P<transfrags>\d+)\s+"
    r"active_transfrags=(?P<active_transfrags>\d+)\s+trflong=(?P<trflong>\d+)"
)
PARITY_BOUNDARY_CHECK_RE = re.compile(
    r"PARITY_BOUNDARY_CHECK\s+nodes=(?P<nodes>\d+)\s+first=(?P<first>\d+)\s+last=(?P<last>\d+)\s+"
    r"rstart=(?P<span_start>\d+)\s+rend=(?P<span_end>\d+)\s+"
    r"is_source=(?P<is_source>\d+)\s+is_sink=(?P<is_sink>\d+)\s+"
    r"hardstart=(?P<hardstart>\d+)\s+hardend=(?P<hardend>\d+)"
)
PARITY_CHECKTRF_START_RE = re.compile(
    r"PARITY_CHECKTRF_START\s+keeptrf=(?P<keeptrf>\d+)\s+checktrf=(?P<checktrf>\d+)"
)
PARITY_CHK_RE = re.compile(
    r"PARITY_CHK\s+t=(?P<seed_idx>\d+)\s+abund=(?P<abundance>[-0-9.]+)\s+outcome=(?P<outcome>\S+)"
)
PARITY_CHECKTRF_OUTCOME_RE = re.compile(
    r"PARITY_CHECKTRF_OUTCOME\s+t=(?P<seed_idx>\d+)\s+outcome=(?P<outcome>\S+)"
)
TIEBREAK_BACK_RE = re.compile(
    r"TIEBREAK_BACK\s+parent=(?P<candidate>\d+)\s+prev_maxp=(?P<prev>\d+)\s+"
    r"parentcov=(?P<candidate_cov>[-0-9.]+)\s+maxcov=(?P<prev_cov>[-0-9.]+)\s+"
    r"winner=(?P<winner>\d+)\s+reason=(?P<reason>\S+)"
)
TIEBREAK_FWD_RE = re.compile(
    r"TIEBREAK_FWD\s+child=(?P<candidate>\d+)\s+prev_maxc=(?P<prev>\d+)\s+"
    r"childcov=(?P<candidate_cov>[-0-9.]+)\s+maxcov=(?P<prev_cov>[-0-9.]+)\s+"
    r"winner=(?P<winner>\d+)\s+reason=(?P<reason>\S+)"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle-summary", required=True, type=Path)
    parser.add_argument("--bundle-id", required=True)
    parser.add_argument("--bundle-log", required=True, type=Path)
    parser.add_argument("--reference-gtf", required=True, type=Path)
    parser.add_argument("--legacy-trace-log", type=Path)
    parser.add_argument("--parity-trace-log", type=Path)
    parser.add_argument("--legacy-trace-norm", required=True, type=Path)
    parser.add_argument("--rust-trace-norm", required=True, type=Path)
    parser.add_argument("--old-transcripts", required=True, type=Path)
    parser.add_argument("--output-stem", required=True)
    return parser.parse_args()


def parse_int(text: str) -> int:
    try:
        return int(text)
    except (TypeError, ValueError):
        return 0


def parse_float(text: str) -> float:
    try:
        return float(text)
    except (TypeError, ValueError):
        return 0.0


def parse_attrs(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in field.split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, value = part.split(" ", 1)
        out[key] = value.strip().strip('"')
    return out


def extract_field(text: str, key: str) -> str:
    match = re.search(rf"\b{re.escape(key)}=([^\s]+)", text)
    return match.group(1) if match else ""


def extract_segment(text: str, key: str) -> str:
    match = re.search(rf"\b{re.escape(key)}=\s*(.+?)(?=\s+[A-Za-z_][A-Za-z0-9_]*=|$)", text)
    return match.group(1).strip() if match else ""


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
        out[key] = text[start:end].strip()
    return prefix, out


def extract_trace_decision(raw: str) -> tuple[str, str]:
    match = re.search(r"^\[TRACE_CHECKTRF\]\s+t=(\d+)\s+([A-Z_]+)\b", raw)
    if not match:
        return "", ""
    return match.group(1), match.group(2)


def span_from_coords(text: str) -> tuple[int, int]:
    coords = [(int(a), int(b)) for a, b in COORD_RE.findall(text)]
    if not coords:
        return 0, 0
    return min(a for a, _ in coords), max(b for _, b in coords)


def spans_overlap(start: int, end: int, focus_start: int, focus_end: int) -> bool:
    return start > 0 and end > 0 and not (end < focus_start or start > focus_end)


def bool_text(value: str) -> str:
    lowered = (value or "").strip().lower()
    if lowered in {"1", "true", "yes"}:
        return "1"
    if lowered in {"0", "false", "no"}:
        return "0"
    return ""


def line_window_from_norm(path: Path, focus_start: int = 0, focus_end: int = 0, pad: int = 20000) -> tuple[int, int]:
    line_nos: list[int] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if not row.get("line_no"):
                continue
            span_start = parse_int(row.get("span_start", ""))
            span_end = parse_int(row.get("span_end", ""))
            if focus_start and focus_end and span_start and span_end:
                if span_end < focus_start or span_start > focus_end:
                    continue
            line_nos.append(parse_int(row["line_no"]))
    if not line_nos:
        return 1, 0
    line_nos.sort()
    cluster = [line_nos[0]]
    max_gap = 50000
    for line_no in line_nos[1:]:
        if line_no - cluster[-1] > max_gap:
            break
        cluster.append(line_no)
    return max(1, cluster[0] - pad), cluster[-1] + pad


def read_bundle_summary(path: Path, bundle_id: str) -> dict[str, str]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("bundle_id") == bundle_id:
                return row
    raise SystemExit(f"bundle_id {bundle_id} not found in {path}")


def load_bundle_log(path: Path) -> list[dict[str, str]]:
    bundles: list[dict[str, str]] = []
    pending: dict[tuple[str, str, str], dict[str, str]] = {}
    with path.open(errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.rstrip("\n")
            begin = BUNDLE_BEGIN_RE.search(text)
            if begin:
                row = {
                    "log_line": str(line_no),
                    "chrom": begin.group("chrom"),
                    "start": begin.group("start"),
                    "end": begin.group("end"),
                    "alignments": begin.group("alignments"),
                    "distinct": begin.group("distinct"),
                    "junctions": begin.group("junctions"),
                    "guides": begin.group("guides"),
                    "processed_transcripts": "",
                }
                bundles.append(row)
                pending[(row["chrom"], row["start"], row["end"])] = row
                continue
            done = BUNDLE_DONE_RE.search(text)
            if done:
                key = (done.group("chrom"), done.group("start"), done.group("end"))
                row = pending.get(key)
                if row is not None:
                    row["processed_transcripts"] = done.group("transcripts")
    return bundles


def choose_expected_bundle(current_bundle: dict[str, str], bundles: list[dict[str, str]]) -> dict[str, str]:
    chrom = current_bundle["chrom"]
    current_start = parse_int(current_bundle["start"])
    current_end = parse_int(current_bundle["end"])
    best: dict[str, str] | None = None
    best_score = (-1, -1)
    for row in bundles:
        if row["chrom"] != chrom:
            continue
        start = parse_int(row["start"])
        end = parse_int(row["end"])
        overlap = min(current_end, end) - max(current_start, start) + 1
        if overlap <= 0:
            continue
        contains = int(start <= current_start and end >= current_end)
        score = (contains, overlap)
        if score > best_score:
            best = row
            best_score = score
    if best is None:
        raise SystemExit("could not map current bundle to any bundle in bundle log")
    return best


def load_reference_overlaps(path: Path, chrom: str, start: int, end: int) -> list[dict[str, str]]:
    transcripts: list[dict[str, str]] = []
    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            tx_chrom = parts[0]
            tx_start = int(parts[3])
            tx_end = int(parts[4])
            if tx_chrom != chrom or tx_start > end or tx_end < start:
                continue
            attrs = parse_attrs(parts[8])
            transcripts.append(
                {
                    "gene_id": attrs.get("gene_id", ""),
                    "transcript_id": attrs.get("transcript_id", ""),
                    "chrom": tx_chrom,
                    "strand": parts[6],
                    "start": str(tx_start),
                    "end": str(tx_end),
                    "num_exons": attrs.get("num_exons", ""),
                    "within_current_bundle": "1" if tx_start >= start and tx_end <= end else "0",
                }
            )
    transcripts.sort(key=lambda row: (parse_int(row["start"]), parse_int(row["end"]), row["transcript_id"]))
    return transcripts


def load_legacy_nodes(path: Path) -> list[dict[str, str]]:
    nodes: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("event") != "SRCSINK_DECISION":
                continue
            node_id = row.get("entity_id", "") or row.get("line_no", "")
            start = row.get("span_start", "")
            end = row.get("span_end", "")
            nodes[node_id] = {
                "node_id": node_id,
                "start": start,
                "end": end,
                "hardstart": bool_text(row.get("hardstart", "")),
                "hardend": bool_text(row.get("hardend", "")),
                "source_state": row.get("source_state", ""),
                "sink_state": row.get("sink_state", ""),
                "raw": row.get("raw", ""),
            }
    return sorted(nodes.values(), key=lambda row: (parse_int(row["start"]), parse_int(row["end"]), row["node_id"]))


def split_hard_pair(raw: str) -> tuple[str, str]:
    match = re.search(r"\bhard=(true|false|0|1)/(true|false|0|1)\b", raw)
    if not match:
        return "", ""
    return bool_text(match.group(1)), bool_text(match.group(2))


def load_raw_transfrags(path: Path) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row.get("event") != "TRACE_RAWTRF":
                continue
            hardstart, hardend = split_hard_pair(row.get("raw", ""))
            out.append(
                {
                    "tf_idx": row.get("entity_id", ""),
                    "start": row.get("span_start", ""),
                    "end": row.get("span_end", ""),
                    "n_nodes": row.get("node_count", ""),
                    "abundance": row.get("abundance", ""),
                    "read_count": "",
                    "srabund": "",
                    "longstart": row.get("longstart", ""),
                    "longend": row.get("longend", ""),
                    "guide": bool_text(row.get("guide", "")),
                    "longread": bool_text(row.get("longread", "") or "1"),
                    "real": bool_text(re.search(r"\breal=(true|false|0|1)\b", row.get("raw", "") or "") and re.search(r"\breal=(true|false|0|1)\b", row.get("raw", "")).group(1) or ""),
                    "weak": "",
                    "coverage_weak": "",
                    "trflong_seed": "",
                    "usepath": "",
                    "group_size": "",
                    "junction_chain": "",
                    "node_ids": "",
                    "hardstart": hardstart,
                    "hardend": hardend,
                    "raw": row.get("raw", ""),
                }
            )
    return out


def parse_trace_records(path: Path) -> tuple[dict[str, dict[str, str]], list[dict[str, str]], dict[str, dict[str, str]]]:
    seeds: dict[str, dict[str, str]] = {}
    edges: list[dict[str, str]] = []
    outputs: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            event = row.get("event", "")
            raw = row.get("raw", "")
            prefix, fields = parse_keyvals_loose(raw)
            _ = prefix
            if event == "TRACE_CHECKTRF_ADD":
                seed_id = fields.get("t", row.get("entity_id", ""))
                seeds[seed_id] = {
                    "seed_id": seed_id,
                    "reason": fields.get("reason", ""),
                    "saved_abund": fields.get("saved_abund", ""),
                    "abundance": fields.get("abund", row.get("abundance", "")),
                    "guide": bool_text(fields.get("guide", row.get("guide", ""))),
                    "longstart": fields.get("longstart", row.get("longstart", "")),
                    "longend": fields.get("longend", row.get("longend", "")),
                    "back": fields.get("back", ""),
                    "fwd": fields.get("fwd", ""),
                    "flux": fields.get("flux", ""),
                    "cov": fields.get("cov", row.get("coverage", "")),
                    "tocheck": fields.get("tocheck", ""),
                    "tf_node_count": fields.get("tf_nodes", row.get("node_count", "")),
                    "tf_coords": fields.get("tf_coords", ""),
                    "path_node_count": fields.get("path_nodes", ""),
                    "path_coords": fields.get("path_coords", ""),
                    "raw": raw,
                }
                if fields.get("path_coords") or fields.get("path_nodes", "") not in {"", "0", "-"}:
                    edges.append(
                        {
                            "source_kind": "seed",
                            "source_id": seed_id,
                            "target_kind": "path",
                            "target_id": fields.get("k", "") or seed_id,
                            "edge_kind": "candidate_path",
                            "reason": fields.get("reason", ""),
                            "coverage": fields.get("cov", ""),
                            "out_idx": "",
                            "raw": raw,
                        }
                    )
                continue

            if not event.startswith("TRACE_CHECKTRF"):
                continue

            decision_seed_id, decision_label = extract_trace_decision(raw)
            seed_id = decision_seed_id or fields.get("t", row.get("entity_id", ""))
            keep_id = fields.get("k", row.get("parent_id", ""))
            out_idx = fields.get("out_idx", "")
            out_key = out_idx or f"seed:{seed_id}"

            if event in {"TRACE_CHECKTRF_KEEP", "TRACE_CHECKTRF_OUT_KEEP"}:
                edges.append(
                    {
                        "source_kind": "seed",
                        "source_id": seed_id,
                        "target_kind": "keep",
                        "target_id": keep_id,
                        "edge_kind": "matched_keep",
                        "reason": row.get("reason", ""),
                        "coverage": fields.get("keep_cov", ""),
                        "out_idx": out_idx,
                        "raw": raw,
                    }
                )
                if out_idx:
                    edges.append(
                        {
                            "source_kind": "keep",
                            "source_id": keep_id,
                            "target_kind": "tx",
                            "target_id": out_idx,
                            "edge_kind": "emits_tx",
                            "reason": row.get("reason", ""),
                            "coverage": fields.get("tx_cov", ""),
                            "out_idx": out_idx,
                            "raw": raw,
                        }
                    )

            if event in {"TRACE_CHECKTRF_OUT", "TRACE_CHECKTRF_OUT_KEEP"}:
                out_key = out_idx or f"seed:{seed_id}:out"
                start, end = span_from_coords(fields.get("exons", ""))
                rec = outputs.setdefault(
                    out_key,
                        {
                            "output_key": out_key,
                            "out_idx": out_idx,
                            "seed_ids": "",
                            "keep_ids": "",
                            "mode": "keep",
                            "coverage": fields.get("cov", ""),
                            "tx_cov": fields.get("tx_cov", ""),
                        "guide": bool_text(fields.get("guide", "")),
                        "start": str(start) if start else "",
                        "end": str(end) if end else "",
                        "exon_count": str(len([p for p in fields.get("exons", "").split(",") if p])) if fields.get("exons") else "",
                        "exons": fields.get("exons", ""),
                        "raw": raw,
                    },
                )
                seed_ids = set(filter(None, rec["seed_ids"].split(",")))
                keep_ids = set(filter(None, rec["keep_ids"].split(",")))
                if seed_id:
                    seed_ids.add(seed_id)
                if keep_id:
                    keep_ids.add(keep_id)
                rec["seed_ids"] = ",".join(sorted(seed_ids, key=parse_int))
                rec["keep_ids"] = ",".join(sorted(keep_ids, key=parse_int))
                rec["mode"] = rec["mode"] or fields.get("mode", "") or "keep"
                rec["coverage"] = rec["coverage"] or fields.get("cov", "")
                rec["tx_cov"] = rec["tx_cov"] or fields.get("tx_cov", "")

            if event == "TRACE_CHECKTRF":
                decision = decision_label or row.get("reason", "") or fields.get("mode", "")
                if decision == "RESCUED_INDEP":
                    out_key = out_idx or f"seed:{seed_id}:rescued_indep"
                    edges.append(
                        {
                            "source_kind": "seed",
                            "source_id": seed_id,
                            "target_kind": "tx",
                            "target_id": out_key,
                            "edge_kind": "rescued_independent",
                            "reason": decision or "RESCUED_INDEP",
                            "coverage": fields.get("cov", row.get("coverage", "")),
                            "out_idx": out_idx,
                            "raw": raw,
                        }
                    )
                    out = outputs.setdefault(
                        out_key,
                        {
                            "output_key": out_key,
                            "out_idx": out_idx,
                            "seed_ids": seed_id,
                            "keep_ids": "",
                            "mode": "rescued_independent",
                            "coverage": fields.get("cov", row.get("coverage", "")),
                            "tx_cov": "",
                            "guide": bool_text(fields.get("guide", "")),
                            "start": fields.get("first", "").split("-")[0] if fields.get("first") else "",
                            "end": fields.get("last", "").split("-")[-1] if fields.get("last") else "",
                            "exon_count": fields.get("exons", ""),
                            "exons": "",
                            "raw": raw,
                        },
                    )
                    if seed_id and seed_id not in set(filter(None, out["seed_ids"].split(","))):
                        out["seed_ids"] = ",".join(filter(None, sorted(set(filter(None, out["seed_ids"].split(","))) | {seed_id}, key=parse_int)))
    return seeds, edges, outputs


def parse_path_update_line(line_no: int, text: str) -> dict[str, str] | None:
    if "PATH_update_abund" not in text:
        return None
    phase_match = re.search(r"PATH_update_abund\[(\d+)\]", text)
    phase = phase_match.group(1) if phase_match else ""
    event = "HEADER"
    if "PATH_update_abund:" in text:
        event = text.split("PATH_update_abund:", 1)[1].strip().split()[0]
    read_match = READ_RE.search(text)
    row = {
        "line_no": str(line_no),
        "phase": phase,
        "event": event,
        "read_idx": read_match.group("read") if read_match else "",
        "s": extract_field(text, "s"),
        "g": extract_field(text, "g"),
        "gno": extract_field(text, "gno"),
        "trf_idx": extract_field(text, "trf_idx"),
        "abundance": extract_field(text, "abund"),
        "total": extract_field(text, "total"),
        "rstart": extract_field(text, "rstart"),
        "rend": extract_field(text, "rend"),
        "longstart": extract_field(text, "longstart"),
        "longend": extract_field(text, "longend"),
        "nnodes": extract_field(text, "nnodes"),
        "node_ids": extract_segment(text, "nodes"),
        "node": extract_field(text, "node"),
        "node0_start": extract_field(text, "node0_start"),
        "nodeLast_end": extract_field(text, "nodeLast_end"),
        "node0_hardstart": extract_field(text, "node0_hardstart"),
        "nodeLast_hardend": extract_field(text, "nodeLast_hardend"),
        "is_source": extract_field(text, "is_source"),
        "is_sink": extract_field(text, "is_sink"),
        "last_i": extract_field(text, "last_i"),
        "is_lr": extract_field(text, "is_lr"),
        "is_sr": extract_field(text, "is_sr"),
        "raw": text,
    }
    if not row["read_idx"] and event == "HEADER":
        return None
    return row


def parse_path_read_pattern_line(line_no: int, text: str) -> dict[str, str] | None:
    if "PATH_read_pattern:" not in text:
        return None
    read_match = READ_RE.search(text)
    if not read_match:
        return None
    payload = text.split("PATH_read_pattern:", 1)[1].strip()
    parts = payload.split()
    event = parts[2] if len(parts) >= 3 else ""
    graph_idx = extract_field(text, "g") or extract_field(text, "ngraph")
    node_field = extract_field(text, "gnode") or extract_field(text, "bnode_gnode")
    node_id = ""
    node_start = ""
    node_end = ""
    node_match = re.match(r"(?P<node_id>\d+)\((?P<start>\d+)-(?P<end>\d+)\)", node_field)
    if node_match:
        node_id = node_match.group("node_id")
        node_start = node_match.group("start")
        node_end = node_match.group("end")
    return {
        "line_no": str(line_no),
        "event": event,
        "read_idx": read_match.group("read"),
        "s": extract_field(text, "s"),
        "graph_idx": graph_idx,
        "node_field": node_field,
        "node_id": node_id,
        "node_start": node_start,
        "node_end": node_end,
        "bp": extract_field(text, "bp"),
        "total_nodes": extract_field(text, "total_nodes"),
        "k": extract_field(text, "k"),
        "ncoord": extract_field(text, "ncoord"),
        "raw": text,
    }


def parse_legacy_lineage(
    legacy_trace_log: Path | None,
    legacy_trace_norm: Path,
    focus_start: int,
    focus_end: int,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    if legacy_trace_log is None:
        return [], [], [], [], [], []

    start_line, end_line = line_window_from_norm(legacy_trace_norm, focus_start, focus_end)
    if end_line <= 0:
        return [], [], [], [], [], []

    read_pattern_rows: list[dict[str, str]] = []
    read_events: list[dict[str, str]] = []
    read_to_rawtrf: dict[tuple[str, str, str, str], dict[str, str]] = {}
    latest_trf_by_read_graph: dict[tuple[str, str, str], str] = {}
    seed_rows: dict[str, dict[str, str]] = {}
    gate_rows: list[dict[str, str]] = []
    junction_rows: list[dict[str, str]] = []
    current_gate_header: dict[str, str] | None = None

    with legacy_trace_log.open(errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if line_no < start_line:
                continue
            if line_no > end_line:
                break
            text = line.rstrip("\n")

            path_pattern_row = parse_path_read_pattern_line(line_no, text)
            if path_pattern_row is not None:
                read_pattern_rows.append(path_pattern_row)
                continue

            path_row = parse_path_update_line(line_no, text)
            if path_row is not None:
                read_events.append(path_row)
                read_key3 = (path_row["read_idx"], path_row["s"], path_row["g"])
                trf_idx = path_row["trf_idx"]
                if trf_idx:
                    latest_trf_by_read_graph[read_key3] = trf_idx
                else:
                    trf_idx = latest_trf_by_read_graph.get(read_key3, "")
                if trf_idx:
                    key = (path_row["read_idx"], path_row["s"], path_row["g"], trf_idx)
                    agg = read_to_rawtrf.setdefault(
                        key,
                        {
                            "read_idx": path_row["read_idx"],
                            "s": path_row["s"],
                            "g": path_row["g"],
                            "trf_idx": trf_idx,
                            "first_line": path_row["line_no"],
                            "last_line": path_row["line_no"],
                            "event_count": "0",
                            "events": "",
                            "abundance": "",
                            "total": "",
                            "rstart": "",
                            "rend": "",
                            "longstart": "",
                            "longend": "",
                            "nnodes": "",
                            "node_ids": "",
                            "is_lr": "",
                            "is_sr": "",
                        },
                    )
                    agg["last_line"] = path_row["line_no"]
                    events = set(filter(None, agg["events"].split(",")))
                    events.add(path_row["event"])
                    agg["events"] = ",".join(sorted(events))
                    agg["event_count"] = str(parse_int(agg["event_count"]) + 1)
                    for field in ["abundance", "total", "rstart", "rend", "longstart", "longend", "nnodes", "node_ids", "is_lr", "is_sr"]:
                        if path_row.get(field):
                            agg[field] = path_row[field]
                continue

            if text.startswith("--- parse_trflong: SEED "):
                prefix, fields = parse_keyvals_loose(text.split("SEED ", 1)[1])
                _ = prefix
                seed_id = fields.get("t", "")
                seed_rows[seed_id] = {
                    "seed_id": seed_id,
                    "line_seed": str(line_no),
                    "abundance": fields.get("abund", ""),
                    "nodes": fields.get("nodes", ""),
                    "guide": bool_text(fields.get("guide", "")),
                    "longread": bool_text(fields.get("longread", "")),
                    "longstart": fields.get("longstart", ""),
                    "longend": fields.get("longend", ""),
                    "node_coords": fields.get("node_coords", ""),
                    "back_result": "",
                    "back_path_len": "",
                    "back_path": "",
                    "fwd_result": "",
                    "fwd_path_len": "",
                    "fwd_path": "",
                    "used_direct": "",
                    "branch": "",
                    "reason": "",
                    "back": "",
                    "fwd": "",
                    "flux": "",
                    "cov": "",
                    "len": "",
                    "checktrf_gate_count": "",
                    "checktrf_gate_keeptrf": "",
                    "checktrf_gate_readthr": "",
                    "gate_pass": "",
                    "raw_seed": text,
                    "raw_decision": "",
                }
                current_gate_header = None
                continue

            if text.startswith("--- parse_trflong: back_to_source "):
                prefix, fields = parse_keyvals_loose(text.split("back_to_source ", 1)[1])
                _ = prefix
                seed = seed_rows.setdefault(fields.get("t", ""), {"seed_id": fields.get("t", "")})
                seed["back_result"] = fields.get("result", "")
                seed["back_path_len"] = fields.get("path_len", "")
                seed["back_path"] = text.split("path:", 1)[1].strip() if "path:" in text else ""
                continue

            if text.startswith("--- parse_trflong: fwd_to_sink "):
                prefix, fields = parse_keyvals_loose(text.split("fwd_to_sink ", 1)[1])
                _ = prefix
                seed = seed_rows.setdefault(fields.get("t", ""), {"seed_id": fields.get("t", "")})
                seed["fwd_result"] = fields.get("result", "")
                seed["fwd_path_len"] = fields.get("path_len", "")
                seed["fwd_path"] = text.split("path:", 1)[1].strip() if "path:" in text else ""
                continue

            if text.startswith("--- parse_trflong: SEED_DECISION "):
                prefix, fields = parse_keyvals_loose(text.split("SEED_DECISION ", 1)[1])
                _ = prefix
                seed = seed_rows.setdefault(fields.get("t", ""), {"seed_id": fields.get("t", "")})
                seed["used_direct"] = fields.get("used_direct", "")
                seed["branch"] = fields.get("branch", "")
                seed["reason"] = fields.get("reason", "")
                seed["back"] = fields.get("back", "")
                seed["fwd"] = fields.get("fwd", "")
                seed["flux"] = fields.get("flux", "")
                seed["cov"] = fields.get("cov", "")
                seed["len"] = fields.get("len", "")
                seed["raw_decision"] = text
                continue

            if text.startswith("--- get_trf_long: CHECKTRF_GATE "):
                prefix, fields = parse_keyvals_loose(text.split("CHECKTRF_GATE ", 1)[1])
                _ = prefix
                current_gate_header = {
                    "line_gate": str(line_no),
                    "count": fields.get("count", ""),
                    "keeptrf": fields.get("keeptrf", ""),
                    "readthr": fields.get("readthr", ""),
                }
                continue

            gate_match = GATE_RE.search(text.strip())
            if gate_match:
                seed_id = gate_match.group("seed_id")
                seed = seed_rows.setdefault(seed_id, {"seed_id": seed_id})
                if current_gate_header:
                    seed["checktrf_gate_count"] = current_gate_header["count"]
                    seed["checktrf_gate_keeptrf"] = current_gate_header["keeptrf"]
                    seed["checktrf_gate_readthr"] = current_gate_header["readthr"]
                seed["gate_pass"] = gate_match.group("passed")
                gate_rows.append(
                    {
                        "line_no": str(line_no),
                        "line_gate": current_gate_header["line_gate"] if current_gate_header else "",
                        "gate_idx": gate_match.group("gate_idx"),
                        "seed_id": seed_id,
                        "abundance": gate_match.group("abundance"),
                        "guide": gate_match.group("guide"),
                        "passed": gate_match.group("passed"),
                        "gate_count": current_gate_header["count"] if current_gate_header else "",
                        "keeptrf": current_gate_header["keeptrf"] if current_gate_header else "",
                        "readthr": current_gate_header["readthr"] if current_gate_header else "",
                        "raw": text,
                    }
                )
                continue

            jdec_match = JDEC_RE.search(text.strip())
            if jdec_match:
                junction_rows.append(
                    {
                        "line_no": str(line_no),
                        "event": "jdec",
                        "junction_idx": jdec_match.group("idx"),
                        "junction": jdec_match.group("junction"),
                        "strand": jdec_match.group("strand"),
                        "nm": jdec_match.group("nm"),
                        "mm": jdec_match.group("mm"),
                        "nreads": jdec_match.group("nreads"),
                        "nreads_good": jdec_match.group("nreads_good"),
                        "left": jdec_match.group("left"),
                        "right": jdec_match.group("right"),
                        "guide": jdec_match.group("guide"),
                        "reason": jdec_match.group("reason"),
                        "raw": text,
                    }
                )
                continue

            good_match = GOOD_JUNC_RE.search(text.strip())
            if good_match:
                junction_rows.append(
                    {
                        "line_no": str(line_no),
                        "event": "good_junc_accepted",
                        "junction_idx": "",
                        "junction": good_match.group("junction"),
                        "strand": good_match.group("strand"),
                        "nm": "",
                        "mm": "",
                        "nreads": good_match.group("nreads"),
                        "nreads_good": "",
                        "left": good_match.group("left"),
                        "right": good_match.group("right"),
                        "guide": "",
                        "reason": "ACCEPTED",
                        "raw": text,
                    }
                )

    kept_read_keys: set[tuple[str, str, str, str]] = set()
    read_to_rawtrf_rows: list[dict[str, str]] = []
    for row in sorted(
        read_to_rawtrf.values(),
        key=lambda rec: (parse_int(rec["read_idx"]), parse_int(rec["s"]), parse_int(rec["g"]), parse_int(rec["trf_idx"])),
    ):
        keep = False
        if spans_overlap(parse_int(row["rstart"]), parse_int(row["rend"]), focus_start, focus_end):
            keep = True
        elif spans_overlap(parse_int(row["longstart"]), parse_int(row["longend"]), focus_start, focus_end):
            keep = True
        else:
            span_start, span_end = span_from_coords(row.get("node_ids", ""))
            if spans_overlap(span_start, span_end, focus_start, focus_end):
                keep = True
        if keep:
            kept_read_keys.add((row["read_idx"], row["s"], row["g"], row["trf_idx"]))
            read_to_rawtrf_rows.append(row)

    filtered_read_events: list[dict[str, str]] = []
    for row in read_events:
        key = (row["read_idx"], row["s"], row["g"], row.get("trf_idx", "") or latest_trf_by_read_graph.get((row["read_idx"], row["s"], row["g"]), ""))
        if key in kept_read_keys:
            filtered_read_events.append(row)
            continue
        row_span_start, row_span_end = span_from_coords(" ".join(filter(None, [row.get("node", ""), row.get("node_ids", "")])))
        if spans_overlap(row_span_start, row_span_end, focus_start, focus_end):
            filtered_read_events.append(row)
            continue
        if spans_overlap(parse_int(row.get("rstart", "")), parse_int(row.get("rend", "")), focus_start, focus_end):
            filtered_read_events.append(row)
            continue
        if spans_overlap(parse_int(row.get("longstart", "")), parse_int(row.get("longend", "")), focus_start, focus_end):
            filtered_read_events.append(row)

    kept_seed_ids: set[str] = set()
    seed_rows_list: list[dict[str, str]] = []
    for row in sorted(seed_rows.values(), key=lambda rec: parse_int(rec["seed_id"])):
        keep = False
        if spans_overlap(parse_int(row.get("longstart", "")), parse_int(row.get("longend", "")), focus_start, focus_end):
            keep = True
        else:
            span_start, span_end = span_from_coords(row.get("node_coords", ""))
            if spans_overlap(span_start, span_end, focus_start, focus_end):
                keep = True
        if keep:
            kept_seed_ids.add(row["seed_id"])
            seed_rows_list.append(row)

    gate_rows = [row for row in gate_rows if row["seed_id"] in kept_seed_ids]
    gate_rows.sort(key=lambda row: (parse_int(row["line_gate"]), parse_int(row["gate_idx"])))

    filtered_junction_rows: list[dict[str, str]] = []
    for row in junction_rows:
        j_start, j_end = span_from_coords(row.get("junction", ""))
        if spans_overlap(j_start, j_end, focus_start, focus_end):
            filtered_junction_rows.append(row)
    filtered_junction_rows.sort(key=lambda row: parse_int(row["line_no"]))

    filtered_read_pattern_rows: list[dict[str, str]] = []
    for row in read_pattern_rows:
        row_span_start, row_span_end = span_from_coords(" ".join(filter(None, [row.get("node_field", ""), row.get("raw", "")])))
        if spans_overlap(row_span_start, row_span_end, focus_start, focus_end):
            filtered_read_pattern_rows.append(row)

    return filtered_read_pattern_rows, filtered_read_events, read_to_rawtrf_rows, seed_rows_list, gate_rows, filtered_junction_rows


def parse_deep_bundle_layers(
    legacy_trace_log: Path | None,
    legacy_trace_norm: Path,
    focus_start: int,
    focus_end: int,
) -> tuple[
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    if legacy_trace_log is None:
        return [], [], [], [], [], [], [], []

    start_line, end_line = line_window_from_norm(legacy_trace_norm, focus_start, focus_end)
    if end_line <= 0:
        return [], [], [], [], [], [], [], []

    graph_rows: list[dict[str, str]] = []
    mutation_rows: list[dict[str, str]] = []
    prediction_rows: list[dict[str, str]] = []
    prediction_edge_rows: list[dict[str, str]] = []
    junction_trace_rows: list[dict[str, str]] = []
    route_rows: list[dict[str, str]] = []
    junction_mutation_rows: list[dict[str, str]] = []
    flow_trace_rows: list[dict[str, str]] = []

    with legacy_trace_log.open(errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if line_no < start_line:
                continue
            if line_no > end_line:
                break
            text = line.rstrip("\n")
            stripped = text.strip()
            source_match = SOURCE_PREFIX_RE.match(stripped)
            source = source_match.group("source") if source_match else ""
            event = ""
            if ":" in stripped:
                after_colon = stripped.split(":", 1)[1].strip()
                if after_colon:
                    event = after_colon.split()[0]
            span_start, span_end = span_from_coords(text)
            bundle_match = SPAN_KV_RE.search(text)
            route_match = ROUTE_BUNDLE_RE.search(text)
            pred_match = PRED_LINE_RE.search(text)
            pred_fate_match = PRED_FATE_RE.search(text)
            pred_summary_match = PRED_FATE_SUMMARY_RE.search(text)
            filter_killed_match = FILTER_KILLED_RE.search(text)
            new_bundle_match = NEW_BUNDLE_RE.search(text)
            new_bnode_match = NEW_BNODE_RE.search(text)
            new_node_match = NEW_NODE_RE.search(text)
            node_final_match = NODE_FINAL_END_RE.search(text)
            srcsink_match = SRCSINK_RE.search(text)
            jdec_match = JDEC_RE.search(stripped)
            good_match = GOOD_JUNC_RE.search(stripped)
            junc_break_match = JUNC_COLOR_BREAK_RE.search(text)
            route_trace_match = ROUTE_RE.search(text)
            junc_delete_match = JUNC_DELETE_RE.search(text)
            junc_demote_match = JUNC_DEMOTE_RE.search(text)
            flow_step_match = FLOW_STEP_RE.search(text)

            row_base = {
                "line_no": str(line_no),
                "source": source,
                "event": event,
                "span_start": str(span_start) if span_start else "",
                "span_end": str(span_end) if span_end else "",
                "bundle_start": bundle_match.group("start") if bundle_match else (route_match.group("start") if route_match else ""),
                "bundle_end": bundle_match.group("end") if bundle_match else (route_match.group("end") if route_match else ""),
                "pred_idx": "",
                "cov": "",
                "readcov": "",
                "exons": "",
                "strand": "",
                "guide": "",
                "fate": "",
                "node_id": "",
                "bno": "",
                "bid": "",
                "grid": extract_field(text, "grid"),
                "color": extract_field(text, "color"),
                "reason": extract_field(text, "reason"),
                "raw": text,
            }

            if pred_match:
                prediction_rows.append(
                    {
                        **row_base,
                        "event": "PRED",
                        "pred_idx": pred_match.group("pred_idx"),
                        "span_start": pred_match.group("start"),
                        "span_end": pred_match.group("end"),
                        "cov": pred_match.group("cov"),
                        "readcov": pred_match.group("readcov"),
                        "strand": pred_match.group("strand"),
                    }
                )
                continue

            if filter_killed_match:
                prediction_edge_rows.append(
                    {
                        **row_base,
                        "event": "FILTER_KILLED_BY",
                        "pred_idx": filter_killed_match.group("killed_idx"),
                        "span_start": filter_killed_match.group("killed_start"),
                        "span_end": filter_killed_match.group("killed_end"),
                        "cov": filter_killed_match.group("killed_cov"),
                        "readcov": filter_killed_match.group("killer_cov"),
                        "exons": filter_killed_match.group("killed_exons"),
                        "strand": filter_killed_match.group("strand"),
                        "guide": filter_killed_match.group("killer_idx"),
                        "fate": f"killer={filter_killed_match.group('killer_idx')}",
                        "reason": filter_killed_match.group("filter_reason"),
                        "node_id": "",
                        "bno": "",
                        "bid": "",
                        "grid": "",
                        "color": "",
                        "raw": text,
                    }
                )

            if pred_fate_match:
                prediction_rows.append(
                    {
                        **row_base,
                        "event": "PRED_FATE",
                        "pred_idx": pred_fate_match.group("pred_idx"),
                        "span_start": pred_fate_match.group("start"),
                        "span_end": pred_fate_match.group("end"),
                        "cov": pred_fate_match.group("cov"),
                        "exons": pred_fate_match.group("exons"),
                        "strand": pred_fate_match.group("strand"),
                        "guide": pred_fate_match.group("guide"),
                        "fate": pred_fate_match.group("fate"),
                    }
                )
                continue

            if pred_summary_match:
                prediction_rows.append(
                    {
                        **row_base,
                        "event": "PRED_FATE_SUMMARY",
                        "pred_idx": pred_summary_match.group("npred"),
                        "cov": pred_summary_match.group("kept"),
                        "readcov": pred_summary_match.group("killed"),
                    }
                )
                continue

            if jdec_match:
                junction_trace_rows.append(
                    {
                        **row_base,
                        "event": "jdec",
                        "span_start": jdec_match.group("junction").split("-", 1)[0],
                        "span_end": jdec_match.group("junction").split("-", 1)[1],
                        "pred_idx": jdec_match.group("idx"),
                        "cov": jdec_match.group("nreads"),
                        "readcov": jdec_match.group("nreads_good"),
                        "strand": jdec_match.group("strand"),
                        "guide": jdec_match.group("guide"),
                        "fate": jdec_match.group("reason"),
                        "node_id": "",
                        "reason": f"mm={jdec_match.group('mm')},left={jdec_match.group('left')},right={jdec_match.group('right')}",
                        "raw": text,
                    }
                )
                continue

            if good_match:
                junction_trace_rows.append(
                    {
                        **row_base,
                        "event": "good_junc_accepted",
                        "span_start": good_match.group("junction").split("-", 1)[0],
                        "span_end": good_match.group("junction").split("-", 1)[1],
                        "cov": good_match.group("nreads"),
                        "strand": good_match.group("strand"),
                        "fate": "ACCEPTED",
                        "reason": f"left={good_match.group('left')},right={good_match.group('right')}",
                        "raw": text,
                    }
                )
                continue

            if junc_break_match:
                junction_trace_rows.append(
                    {
                        **row_base,
                        "event": "JUNC_COLOR_BREAK",
                        "span_start": junc_break_match.group("junction").split("-", 1)[0],
                        "span_end": junc_break_match.group("junction").split("-", 1)[1],
                        "pred_idx": junc_break_match.group("read_idx"),
                        "cov": junc_break_match.group("nreads"),
                        "readcov": junc_break_match.group("nreads_good"),
                        "strand": junc_break_match.group("strand_now"),
                        "guide": junc_break_match.group("good_junc"),
                        "fate": f"order={junc_break_match.group('junc_order')}",
                        "reason": (
                            f"mm={junc_break_match.group('mm')},"
                            f"changeleft={junc_break_match.group('changeleft')},"
                            f"changeright={junc_break_match.group('changeright')}"
                        ),
                        "raw": text,
                    }
                )
                continue

            if route_trace_match:
                route_rows.append(
                    {
                        **row_base,
                        "event": f"ROUTE_{route_trace_match.group('route_kind')}",
                        "span_start": route_trace_match.group("start"),
                        "span_end": route_trace_match.group("end"),
                        "pred_idx": route_trace_match.group("pass"),
                        "node_id": route_trace_match.group("sno"),
                        "bno": extract_field(text, "bno"),
                        "grid": route_trace_match.group("grid_idx"),
                        "color": extract_field(text, "color"),
                        "guide": extract_field(text, "route"),
                        "cov": extract_field(text, "neg_prop"),
                        "readcov": extract_field(text, "lastnodeid"),
                        "reason": extract_field(text, "reason"),
                        "raw": text,
                    }
                )
                continue

            if junc_delete_match:
                junction_mutation_rows.append(
                    {
                        **row_base,
                        "event": junc_delete_match.group("event"),
                        "span_start": junc_delete_match.group("start"),
                        "span_end": junc_delete_match.group("end"),
                        "pred_idx": junc_delete_match.group("idx"),
                        "strand": junc_delete_match.group("strand"),
                        "cov": junc_delete_match.group("nreads_good"),
                        "raw": text,
                    }
                )
                continue

            if junc_demote_match:
                junction_mutation_rows.append(
                    {
                        **row_base,
                        "event": "JUNC_DEMOTE",
                        "span_start": junc_demote_match.group("from_start"),
                        "span_end": junc_demote_match.group("from_end"),
                        "pred_idx": junc_demote_match.group("from_idx"),
                        "strand": "",
                        "cov": junc_demote_match.group("left_from"),
                        "readcov": junc_demote_match.group("left_to"),
                        "guide": junc_demote_match.group("to_idx"),
                        "fate": f"{junc_demote_match.group('to_start')}-{junc_demote_match.group('to_end')}",
                        "raw": text,
                    }
                )
                continue

            if flow_step_match:
                flow_trace_rows.append(
                    {
                        **row_base,
                        "event": f"{flow_step_match.group('phase')}_step",
                        "span_start": flow_step_match.group("start"),
                        "span_end": flow_step_match.group("end"),
                        "pred_idx": flow_step_match.group("i"),
                        "node_id": flow_step_match.group("maxc"),
                        "cov": flow_step_match.group("maxcov"),
                        "raw": text,
                    }
                )
                continue

            if (
                "HARD_BOUNDARY_CHECK" in text
                or "PATHPAT_OR" in text
                or "TIEBREAK" in text
            ):
                flow_trace_rows.append(
                    {
                        **row_base,
                        "event": event or "FLOW_TRACE",
                        "raw": text,
                    }
                )
                continue

            if source in {"print_predcluster", "infer_transcripts"}:
                prediction_rows.append(row_base)
                continue

            if source in {"build_graphs", "create_bundle", "create_graph", "merge_read_to_group"}:
                row = dict(row_base)
                if new_bundle_match:
                    row["bno"] = new_bundle_match.group("bno")
                    row["bid"] = new_bundle_match.group("bid")
                    row["span_start"] = new_bundle_match.group("start")
                    row["span_end"] = new_bundle_match.group("end")
                    row["cov"] = new_bundle_match.group("cov")
                elif new_bnode_match:
                    row["bid"] = new_bnode_match.group("bid")
                    row["span_start"] = new_bnode_match.group("start")
                    row["span_end"] = new_bnode_match.group("end")
                elif new_node_match:
                    row["node_id"] = new_node_match.group("nodeid")
                    row["span_start"] = new_node_match.group("start")
                    row["span_end"] = new_node_match.group("end")
                elif node_final_match:
                    row["node_id"] = node_final_match.group("nodeid")
                    row["span_start"] = node_final_match.group("start")
                    row["span_end"] = node_final_match.group("end")
                graph_rows.append(row)
                continue

            if source == "process_transfrags" or "LONGTRIM_" in text or "eliminate_trf" in text:
                row = dict(row_base)
                if srcsink_match:
                    row["node_id"] = srcsink_match.group("nodeid")
                    row["span_start"] = srcsink_match.group("start")
                    row["span_end"] = srcsink_match.group("end")
                    row["guide"] = ""
                    row["reason"] = f"source={srcsink_match.group('source_state')},sink={srcsink_match.group('sink_state')}"
                    row["cov"] = srcsink_match.group("addsource")
                    row["readcov"] = srcsink_match.group("addsink")
                    row["fate"] = f"hardstart={srcsink_match.group('hardstart')},hardend={srcsink_match.group('hardend')}"
                mutation_rows.append(row)
                continue

    return (
        graph_rows,
        mutation_rows,
        prediction_rows,
        prediction_edge_rows,
        junction_trace_rows,
        route_rows,
        junction_mutation_rows,
        flow_trace_rows,
    )


def parse_parity_checktrf_layers(
    parity_trace_log: Path | None,
    focus_start: int,
    focus_end: int,
) -> tuple[
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    if parity_trace_log is None:
        return [], [], [], [], [], []

    seed_rows: list[dict[str, str]] = []
    batch_rows: list[dict[str, str]] = []
    event_rows: list[dict[str, str]] = []
    tiebreak_rows: list[dict[str, str]] = []
    graph_rows: list[dict[str, str]] = []
    boundary_rows: list[dict[str, str]] = []

    recent_bundle_start = 0
    recent_bundle_end = 0
    recent_bundle_match = False
    parity_bundle_spans: dict[tuple[str, str], tuple[int, int]] = {}
    active_batch_id = ""
    active_batch_match = False
    active_batch_kind = ""
    batch_counter = 0

    def open_batch(line_no: int, raw: str, batch_kind: str, keeptrf: str = "", checktrf: str = "") -> None:
        nonlocal batch_counter, active_batch_id, active_batch_match, active_batch_kind
        batch_counter += 1
        active_batch_id = str(batch_counter)
        active_batch_match = recent_bundle_match
        active_batch_kind = batch_kind
        if active_batch_match:
            batch_rows.append(
                {
                    "line_no": str(line_no),
                    "batch_id": active_batch_id,
                    "batch_kind": batch_kind,
                    "bundle_start": str(recent_bundle_start) if recent_bundle_start else "",
                    "bundle_end": str(recent_bundle_end) if recent_bundle_end else "",
                    "keeptrf": keeptrf,
                    "checktrf": checktrf,
                    "raw": raw,
                }
            )

    def close_synthetic_batch() -> None:
        nonlocal active_batch_id, active_batch_match, active_batch_kind
        if active_batch_kind == "synthetic":
            active_batch_id = ""
            active_batch_match = False
            active_batch_kind = ""

    with parity_trace_log.open(errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.rstrip("\n")

            parity_bundle_match = PARITY_BUNDLE_RE.search(text)
            if parity_bundle_match:
                sno = parity_bundle_match.group("sno")
                bundle_idx = parity_bundle_match.group("bundle_idx")
                bundle_start = parse_int(parity_bundle_match.group("bundle_start"))
                bundle_end = parse_int(parity_bundle_match.group("bundle_end"))
                parity_bundle_spans[(sno, bundle_idx)] = (bundle_start, bundle_end)
                continue

            graph_evolution_match = PARITY_GRAPH_EVOLUTION_RE.search(text)
            if graph_evolution_match:
                sno = graph_evolution_match.group("sno")
                graph_idx = graph_evolution_match.group("graph_idx")
                bundle_start, bundle_end = parity_bundle_spans.get((sno, graph_idx), (0, 0))
                if spans_overlap(bundle_start, bundle_end, focus_start, focus_end):
                    graph_rows.append(
                        {
                            "line_no": str(line_no),
                            "row_kind": "PARITY_GRAPH_EVOLUTION",
                            "sno": sno,
                            "bundle_idx": graph_idx,
                            "graph_idx": graph_idx,
                            "bundle_start": str(bundle_start),
                            "bundle_end": str(bundle_end),
                            "event": graph_evolution_match.group("event"),
                            "node_id": graph_evolution_match.group("node_id"),
                            "span_start": graph_evolution_match.group("span_start"),
                            "span_end": graph_evolution_match.group("span_end"),
                            "bid": graph_evolution_match.group("bid"),
                            "parents": graph_evolution_match.group("parents"),
                            "nodes": "",
                            "transfrags": "",
                            "active_transfrags": "",
                            "trflong": "",
                            "raw": text,
                        }
                    )
                continue

            graph_post_match = PARITY_GRAPH_POST_RE.search(text)
            if graph_post_match:
                sno = graph_post_match.group("sno")
                bundle_idx = graph_post_match.group("bundle_idx")
                bundle_start, bundle_end = parity_bundle_spans.get((sno, bundle_idx), (0, 0))
                if spans_overlap(bundle_start, bundle_end, focus_start, focus_end):
                    graph_rows.append(
                        {
                            "line_no": str(line_no),
                            "row_kind": "PARITY_GRAPH_POST",
                            "sno": sno,
                            "bundle_idx": bundle_idx,
                            "graph_idx": bundle_idx,
                            "bundle_start": str(bundle_start),
                            "bundle_end": str(bundle_end),
                            "event": "POST",
                            "node_id": "",
                            "span_start": "",
                            "span_end": "",
                            "bid": "",
                            "parents": "",
                            "nodes": graph_post_match.group("nodes"),
                            "transfrags": graph_post_match.group("transfrags"),
                            "active_transfrags": graph_post_match.group("active_transfrags"),
                            "trflong": graph_post_match.group("trflong"),
                            "raw": text,
                        }
                    )
                continue

            boundary_check_match = PARITY_BOUNDARY_CHECK_RE.search(text)
            if boundary_check_match:
                span_start = parse_int(boundary_check_match.group("span_start"))
                span_end = parse_int(boundary_check_match.group("span_end"))
                if spans_overlap(span_start, span_end, focus_start, focus_end):
                    boundary_rows.append(
                        {
                            "line_no": str(line_no),
                            "bundle_start": str(span_start),
                            "bundle_end": str(span_end),
                            "span_start": boundary_check_match.group("span_start"),
                            "span_end": boundary_check_match.group("span_end"),
                            "nodes": boundary_check_match.group("nodes"),
                            "first": boundary_check_match.group("first"),
                            "last": boundary_check_match.group("last"),
                            "is_source": boundary_check_match.group("is_source"),
                            "is_sink": boundary_check_match.group("is_sink"),
                            "hardstart": boundary_check_match.group("hardstart"),
                            "hardend": boundary_check_match.group("hardend"),
                            "raw": text,
                        }
                    )
                continue

            seed_match = PARITY_SEED_PROC_RE.search(text)
            if seed_match:
                recent_bundle_start = parse_int(seed_match.group("bundle_start"))
                recent_bundle_end = parse_int(seed_match.group("bundle_end"))
                recent_bundle_match = spans_overlap(recent_bundle_start, recent_bundle_end, focus_start, focus_end)
                if recent_bundle_match:
                    seed_rows.append(
                        {
                            "line_no": str(line_no),
                            "bundle_start": seed_match.group("bundle_start"),
                            "bundle_end": seed_match.group("bundle_end"),
                            "seed_idx": seed_match.group("seed_idx"),
                            "span_start": seed_match.group("span_start"),
                            "span_end": seed_match.group("span_end"),
                            "abundance": seed_match.group("abundance"),
                            "nodes": seed_match.group("nodes"),
                            "back": seed_match.group("back"),
                            "fwd": seed_match.group("fwd"),
                            "flux": seed_match.group("flux"),
                            "path": seed_match.group("path"),
                            "outcome": seed_match.group("outcome"),
                            "raw": text,
                        }
                    )
                close_synthetic_batch()
                continue

            batch_match = PARITY_CHECKTRF_START_RE.search(text)
            if batch_match:
                open_batch(
                    line_no,
                    text,
                    batch_kind="explicit",
                    keeptrf=batch_match.group("keeptrf"),
                    checktrf=batch_match.group("checktrf"),
                )
                continue

            back_tiebreak_match = TIEBREAK_BACK_RE.search(text)
            if back_tiebreak_match and recent_bundle_match:
                tiebreak_rows.append(
                    {
                        "line_no": str(line_no),
                        "source": "parity_trace",
                        "event": "TIEBREAK_BACK",
                        "span_start": "",
                        "span_end": "",
                        "bundle_start": str(recent_bundle_start),
                        "bundle_end": str(recent_bundle_end),
                        "pred_idx": back_tiebreak_match.group("candidate"),
                        "cov": back_tiebreak_match.group("candidate_cov"),
                        "readcov": back_tiebreak_match.group("prev_cov"),
                        "exons": "",
                        "strand": "",
                        "guide": "",
                        "fate": f"prev={back_tiebreak_match.group('prev')}",
                        "node_id": back_tiebreak_match.group("winner"),
                        "bno": "",
                        "bid": "",
                        "grid": "",
                        "color": "",
                        "reason": back_tiebreak_match.group("reason"),
                        "raw": text,
                    }
                )
                close_synthetic_batch()
                continue

            fwd_tiebreak_match = TIEBREAK_FWD_RE.search(text)
            if fwd_tiebreak_match and recent_bundle_match:
                tiebreak_rows.append(
                    {
                        "line_no": str(line_no),
                        "source": "parity_trace",
                        "event": "TIEBREAK_FWD",
                        "span_start": "",
                        "span_end": "",
                        "bundle_start": str(recent_bundle_start),
                        "bundle_end": str(recent_bundle_end),
                        "pred_idx": fwd_tiebreak_match.group("candidate"),
                        "cov": fwd_tiebreak_match.group("candidate_cov"),
                        "readcov": fwd_tiebreak_match.group("prev_cov"),
                        "exons": "",
                        "strand": "",
                        "guide": "",
                        "fate": f"prev={fwd_tiebreak_match.group('prev')}",
                        "node_id": fwd_tiebreak_match.group("winner"),
                        "bno": "",
                        "bid": "",
                        "grid": "",
                        "color": "",
                        "reason": fwd_tiebreak_match.group("reason"),
                        "raw": text,
                    }
                )
                close_synthetic_batch()
                continue

            chk_match = PARITY_CHK_RE.search(text)
            if chk_match:
                if not active_batch_id:
                    open_batch(
                        line_no,
                        "[synthetic] PARITY_CHK run",
                        batch_kind="synthetic",
                    )
                if active_batch_match:
                    event_rows.append(
                        {
                            "line_no": str(line_no),
                            "batch_id": active_batch_id,
                            "event": "PARITY_CHK",
                            "seed_idx": chk_match.group("seed_idx"),
                            "abundance": chk_match.group("abundance"),
                            "outcome": chk_match.group("outcome"),
                            "matched": extract_field(text, "matched"),
                            "abundsum": extract_field(text, "abundsum"),
                            "cov": extract_field(text, "cov"),
                            "nexons": extract_field(text, "nexons"),
                            "exons": extract_field(text, "exons"),
                            "matches": "",
                            "nodes": "",
                            "raw": text,
                        }
                    )
                continue

            outcome_match = PARITY_CHECKTRF_OUTCOME_RE.search(text)
            if outcome_match:
                if not active_batch_id:
                    open_batch(
                        line_no,
                        "[synthetic] PARITY_CHECKTRF_OUTCOME run",
                        batch_kind="synthetic",
                    )
                if active_batch_match:
                    event_rows.append(
                        {
                            "line_no": str(line_no),
                            "batch_id": active_batch_id,
                            "event": "PARITY_CHECKTRF_OUTCOME",
                            "seed_idx": outcome_match.group("seed_idx"),
                            "abundance": "",
                            "outcome": outcome_match.group("outcome"),
                            "matched": "",
                            "abundsum": "",
                            "cov": "",
                            "nexons": "",
                            "exons": "",
                            "matches": extract_field(text, "matches"),
                            "nodes": extract_field(text, "nodes"),
                            "raw": text,
                        }
                    )
                continue

            close_synthetic_batch()

    return seed_rows, batch_rows, event_rows, tiebreak_rows, graph_rows, boundary_rows


def summarize_read_graph_paths(read_events: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: dict[tuple[str, str, str], dict[str, str]] = {}
    for event in read_events:
        key = (event.get("read_idx", ""), event.get("s", ""), event.get("g", ""))
        rec = rows.setdefault(
            key,
            {
                "read_idx": event.get("read_idx", ""),
                "s": event.get("s", ""),
                "g": event.get("g", ""),
                "gno": "",
                "first_line": event.get("line_no", ""),
                "last_line": event.get("line_no", ""),
                "event_count": "0",
                "events": "",
                "trf_idx": "",
                "rstart": "",
                "rend": "",
                "node0_start": "",
                "nodeLast_end": "",
                "node0_hardstart": "",
                "nodeLast_hardend": "",
                "is_source": "",
                "is_sink": "",
                "nnodes": "",
                "node_ids": "",
                "longstart": "",
                "longend": "",
                "source_check_entered": "0",
                "sink_check_entered": "0",
                "skip_single_node": "0",
                "new_trf_seen": "0",
                "final_nodes_seen": "0",
            },
        )
        rec["last_line"] = event.get("line_no", "")
        rec["event_count"] = str(parse_int(rec["event_count"]) + 1)
        events = set(filter(None, rec["events"].split(",")))
        events.add(event.get("event", ""))
        rec["events"] = ",".join(sorted(events))
        for field in [
            "gno",
            "trf_idx",
            "rstart",
            "rend",
            "node0_start",
            "nodeLast_end",
            "node0_hardstart",
            "nodeLast_hardend",
            "is_source",
            "is_sink",
            "nnodes",
            "node_ids",
            "longstart",
            "longend",
        ]:
            if event.get(field):
                rec[field] = event[field]
        event_name = event.get("event", "")
        if event_name == "SOURCE_CHECK_ENTER":
            rec["source_check_entered"] = "1"
        elif event_name == "SINK_CHECK_ENTER":
            rec["sink_check_entered"] = "1"
        elif event_name == "SKIP_SINGLE_NODE":
            rec["skip_single_node"] = "1"
        elif event_name == "NEW_TRF":
            rec["new_trf_seen"] = "1"
        elif event_name == "FINAL_NODES":
            rec["final_nodes_seen"] = "1"
    return sorted(rows.values(), key=lambda row: (parse_int(row["read_idx"]), parse_int(row["s"]), parse_int(row["g"])))


def summarize_read_pattern_paths(read_pattern_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    rows: dict[tuple[str, str, str], dict[str, str]] = {}
    for event in read_pattern_rows:
        key = (event.get("read_idx", ""), event.get("s", ""), event.get("graph_idx", ""))
        rec = rows.setdefault(
            key,
            {
                "read_idx": event.get("read_idx", ""),
                "s": event.get("s", ""),
                "graph_idx": event.get("graph_idx", ""),
                "first_line": event.get("line_no", ""),
                "last_line": event.get("line_no", ""),
                "event_count": "0",
                "events": "",
                "add_node_count": "0",
                "intersect_count": "0",
                "exon_skip_count": "0",
                "total_nodes_max": "",
                "node_ids": "",
                "node_spans": "",
                "skip_bnodes": "",
                "max_bp": "",
            },
        )
        rec["last_line"] = event.get("line_no", "")
        rec["event_count"] = str(parse_int(rec["event_count"]) + 1)
        events = set(filter(None, rec["events"].split(",")))
        events.add(event.get("event", ""))
        rec["events"] = ",".join(sorted(events))
        if event.get("event") == "ADD_NODE":
            rec["add_node_count"] = str(parse_int(rec["add_node_count"]) + 1)
            if event.get("node_id"):
                node_ids = set(filter(None, rec["node_ids"].split(",")))
                node_ids.add(event["node_id"])
                rec["node_ids"] = ",".join(sorted(node_ids, key=parse_int))
            if event.get("node_start") and event.get("node_end"):
                spans = set(filter(None, rec["node_spans"].split(";")))
                spans.add(f"{event['node_start']}-{event['node_end']}")
                rec["node_spans"] = ";".join(
                    sorted(spans, key=lambda s: (parse_int(s.split("-")[0]), parse_int(s.split("-")[-1])))
                )
        elif event.get("event") == "INTERSECT":
            rec["intersect_count"] = str(parse_int(rec["intersect_count"]) + 1)
        elif event.get("event") == "EXON_SKIP":
            rec["exon_skip_count"] = str(parse_int(rec["exon_skip_count"]) + 1)
            if event.get("node_field"):
                skip_bnodes = set(filter(None, rec["skip_bnodes"].split(";")))
                skip_bnodes.add(event["node_field"])
                rec["skip_bnodes"] = ";".join(sorted(skip_bnodes))
        total_nodes = parse_int(event.get("total_nodes", ""))
        if total_nodes and total_nodes > parse_int(rec["total_nodes_max"]):
            rec["total_nodes_max"] = str(total_nodes)
        bp = parse_int(event.get("bp", ""))
        if bp and bp > parse_int(rec["max_bp"]):
            rec["max_bp"] = str(bp)
    return sorted(rows.values(), key=lambda row: (parse_int(row["read_idx"]), parse_int(row["s"]), parse_int(row["graph_idx"])))


def summarize_rawtrf_provenance(
    read_to_rawtrf: list[dict[str, str]],
    read_graph_paths: list[dict[str, str]],
) -> list[dict[str, str]]:
    path_by_key = {
        (row["read_idx"], row["s"], row["g"]): row
        for row in read_graph_paths
    }
    rows: dict[str, dict[str, str]] = {}
    for link in read_to_rawtrf:
        trf_idx = link.get("trf_idx", "")
        if not trf_idx:
            continue
        rec = rows.setdefault(
            trf_idx,
            {
                "trf_idx": trf_idx,
                "read_count": "0",
                "graph_count": "0",
                "graphs": "",
                "first_line": link.get("first_line", ""),
                "last_line": link.get("last_line", ""),
                "events": "",
                "abundance_sum": "0",
                "total_max": "",
                "min_rstart": "",
                "max_rend": "",
                "min_longstart_nonzero": "",
                "max_longend_nonzero": "",
                "any_source": "0",
                "any_sink": "0",
                "any_hardstart": "0",
                "any_hardend": "0",
                "source_check_reads": "0",
                "sink_check_reads": "0",
                "skip_single_node_reads": "0",
                "node_ids_examples": "",
            },
        )
        rec["read_count"] = str(parse_int(rec["read_count"]) + 1)
        rec["last_line"] = link.get("last_line", "")
        events = set(filter(None, rec["events"].split(",")))
        events.update(filter(None, link.get("events", "").split(",")))
        rec["events"] = ",".join(sorted(events))
        rec["abundance_sum"] = f"{parse_float(rec['abundance_sum']) + parse_float(link.get('abundance', '')):.4f}"
        total = parse_float(link.get("total", ""))
        if total and total > parse_float(rec.get("total_max", "") or "0"):
            rec["total_max"] = link.get("total", "")

        for field, mode in [("rstart", "min"), ("rend", "max"), ("longstart", "min_nonzero"), ("longend", "max_nonzero")]:
            value = parse_int(link.get(field, ""))
            target = {
                "rstart": "min_rstart",
                "rend": "max_rend",
                "longstart": "min_longstart_nonzero",
                "longend": "max_longend_nonzero",
            }[field]
            current = parse_int(rec.get(target, ""))
            if mode == "min":
                if value and (current == 0 or value < current):
                    rec[target] = str(value)
            elif mode == "max":
                if value and value > current:
                    rec[target] = str(value)
            elif mode == "min_nonzero":
                if value and (current == 0 or value < current):
                    rec[target] = str(value)
            elif mode == "max_nonzero":
                if value and value > current:
                    rec[target] = str(value)

        path = path_by_key.get((link["read_idx"], link["s"], link["g"]))
        if path:
            graph_key = f"{path.get('s','')}:{path.get('g','')}"
            graphs = set(filter(None, rec["graphs"].split(",")))
            graphs.add(graph_key)
            rec["graphs"] = ",".join(sorted(graphs, key=lambda x: tuple(parse_int(v) for v in x.split(":"))))
            rec["graph_count"] = str(len(graphs))
            if bool_text(path.get("is_source", "")) == "1":
                rec["any_source"] = "1"
            if bool_text(path.get("is_sink", "")) == "1":
                rec["any_sink"] = "1"
            if bool_text(path.get("node0_hardstart", "")) == "1":
                rec["any_hardstart"] = "1"
            if bool_text(path.get("nodeLast_hardend", "")) == "1":
                rec["any_hardend"] = "1"
            if path.get("source_check_entered") == "1":
                rec["source_check_reads"] = str(parse_int(rec["source_check_reads"]) + 1)
            if path.get("sink_check_entered") == "1":
                rec["sink_check_reads"] = str(parse_int(rec["sink_check_reads"]) + 1)
            if path.get("skip_single_node") == "1":
                rec["skip_single_node_reads"] = str(parse_int(rec["skip_single_node_reads"]) + 1)
            examples = set(filter(None, rec["node_ids_examples"].split(" || ")))
            if path.get("node_ids") and len(examples) < 3:
                examples.add(path["node_ids"])
            rec["node_ids_examples"] = " || ".join(sorted(examples))
    return sorted(rows.values(), key=lambda row: parse_int(row["trf_idx"]))


def load_transcripts(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_transcript_support(transcripts: list[dict[str, str]], outputs: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    by_out_idx = {row["out_idx"]: row for row in outputs.values() if row.get("out_idx")}
    support_rows: list[dict[str, str]] = []
    for tx in transcripts:
        support = by_out_idx.get(tx.get("tx_idx", ""))
        support_rows.append(
            {
                "tx_idx": tx.get("tx_idx", ""),
                "start": tx.get("start", ""),
                "end": tx.get("end", ""),
                "num_exons": tx.get("num_exons", ""),
                "coverage": tx.get("coverage", ""),
                "longcov": tx.get("longcov", ""),
                "junction_chain": tx.get("junction_chain", ""),
                "checktrf_mode": support.get("mode", "") if support else "",
                "checktrf_seed_ids": support.get("seed_ids", "") if support else "",
                "checktrf_keep_ids": support.get("keep_ids", "") if support else "",
                "checktrf_cov": support.get("coverage", "") if support else "",
                "checktrf_tx_cov": support.get("tx_cov", "") if support else "",
            }
        )
    return support_rows


def scan_longtrim(path: Path, expected_bundle: dict[str, str]) -> dict[str, str]:
    target = f"bundle={expected_bundle['start']}-{expected_bundle['end']}"
    counts = Counter()
    examples: dict[str, str] = {}
    with path.open(errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            if target not in line:
                continue
            if "LONGTRIM_BNODE" in line:
                counts["LONGTRIM_BNODE"] += 1
                examples.setdefault("LONGTRIM_BNODE", f"{line_no}:{line.rstrip()}")
            if "LONGTRIM_BOUND" in line:
                counts["LONGTRIM_BOUND"] += 1
                examples.setdefault("LONGTRIM_BOUND", f"{line_no}:{line.rstrip()}")
    return {
        "longtrim_bnode": str(counts["LONGTRIM_BNODE"]),
        "longtrim_bound": str(counts["LONGTRIM_BOUND"]),
        "longtrim_bnode_example": examples.get("LONGTRIM_BNODE", ""),
        "longtrim_bound_example": examples.get("LONGTRIM_BOUND", ""),
    }


def write_tsv(path: Path, rows: Iterable[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def build_layer_rows(
    raw_nodes: list[dict[str, str]],
    raw_transfrags: list[dict[str, str]],
    read_pattern_rows: list[dict[str, str]],
    read_events: list[dict[str, str]],
    read_graph_paths: list[dict[str, str]],
    read_to_rawtrf: list[dict[str, str]],
    rawtrf_provenance: list[dict[str, str]],
    seed_lineage: list[dict[str, str]],
    gate_rows: list[dict[str, str]],
    junction_rows: list[dict[str, str]],
    graph_build_rows: list[dict[str, str]],
    graph_mutation_rows: list[dict[str, str]],
    prediction_rows: list[dict[str, str]],
    prediction_edge_rows: list[dict[str, str]],
    junction_trace_rows: list[dict[str, str]],
    route_rows: list[dict[str, str]],
    junction_mutation_rows: list[dict[str, str]],
    flow_trace_rows: list[dict[str, str]],
    seeds: dict[str, dict[str, str]],
    edges: list[dict[str, str]],
    outputs: dict[str, dict[str, str]],
    transcripts: list[dict[str, str]],
    transcript_support: list[dict[str, str]],
    parity_checktrf_seeds: list[dict[str, str]],
    parity_checktrf_batches: list[dict[str, str]],
    parity_checktrf_events: list[dict[str, str]],
    parity_graph_rows: list[dict[str, str]],
    parity_boundary_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    return [
        {
            "layer_id": "L0_bundle_graph_context",
            "count": str(len(graph_build_rows) + len(route_rows)),
            "primary_tables": "lineage.graph_build.tsv,lineage.routes.tsv",
            "detail": "Bundle/group/bnode construction and routing context before path extraction.",
        },
        {
            "layer_id": "L1_junction_policy",
            "count": str(len(junction_rows) + len(junction_trace_rows) + len(junction_mutation_rows)),
            "primary_tables": "lineage.junction_decisions.tsv,lineage.junction_trace.tsv,lineage.junction_mutation.tsv",
            "detail": "Junction acceptance, rejection, and mutation history.",
        },
        {
            "layer_id": "L2_path_lineage",
            "count": str(len(read_pattern_rows) + len(read_events) + len(read_graph_paths) + len(read_to_rawtrf)),
            "primary_tables": "lineage.read_patterns.tsv,lineage.read_events.tsv,lineage.read_graph_paths.tsv,lineage.read_to_rawtrf.tsv",
            "detail": "Read-pattern and read-to-graph-path lineage before raw-transfrag consolidation.",
        },
        {
            "layer_id": "L3_raw_transfrag_state",
            "count": str(len(raw_nodes) + len(raw_transfrags) + len(rawtrf_provenance) + len(graph_mutation_rows) + len(parity_graph_rows)),
            "primary_tables": "rawgraph.nodes.tsv,rawgraph.transfrags.tsv,lineage.rawtrf_provenance.tsv,lineage.graph_mutation.tsv,lineage.parity_graph.tsv",
            "detail": "Raw transfrags, node source/sink state, boundary provenance, and graph mutation history from deep and parity traces.",
        },
        {
            "layer_id": "L4_seed_flow",
            "count": str(len(seed_lineage) + len(flow_trace_rows) + len(parity_boundary_rows)),
            "primary_tables": "lineage.seed_decisions.tsv,lineage.flow_trace.tsv,lineage.parity_boundary_checks.tsv",
            "detail": "parse_trflong seed expansion, direct-vs-deferred choice, and low-level flow/boundary/tiebreak traces from both deep and parity logs.",
        },
        {
            "layer_id": "L5_checktrf_rescue",
            "count": str(len(seeds) + len(gate_rows) + len(edges) + len(outputs) + len(parity_checktrf_seeds) + len(parity_checktrf_batches) + len(parity_checktrf_events)),
            "primary_tables": "checktrf.seeds.tsv,checktrf.edges.tsv,checktrf.outputs.tsv,lineage.checktrf_gate.tsv,lineage.parity_checktrf_*.tsv",
            "detail": "Deferred rescue gate, seed/output mapping, and parity checktrf batches/events, including synthetic batches when the newer parity log omits explicit batch markers.",
        },
        {
            "layer_id": "L6_prediction_fate",
            "count": str(len(prediction_rows) + len(prediction_edge_rows)),
            "primary_tables": "lineage.predictions.tsv,lineage.prediction_edges.tsv",
            "detail": "Prediction filtering and killer-vs-killed fate edges before final emission.",
        },
        {
            "layer_id": "L7_transcript_emission",
            "count": str(len(transcripts) + len(transcript_support)),
            "primary_tables": "transcripts.tsv,transcript_support.tsv",
            "detail": "Final old-path transcript emission linked back to checktrf support when available.",
        },
    ]


def write_report(
    path: Path,
    current_bundle: dict[str, str],
    expected_bundle: dict[str, str],
    references: list[dict[str, str]],
    layer_rows: list[dict[str, str]],
    raw_nodes: list[dict[str, str]],
    raw_transfrags: list[dict[str, str]],
    seeds: dict[str, dict[str, str]],
    edges: list[dict[str, str]],
    outputs: dict[str, dict[str, str]],
    transcripts: list[dict[str, str]],
    transcript_support: list[dict[str, str]],
    longtrim: dict[str, str],
    read_pattern_rows: list[dict[str, str]],
    read_events: list[dict[str, str]],
    read_to_rawtrf: list[dict[str, str]],
    read_graph_paths: list[dict[str, str]],
    rawtrf_provenance: list[dict[str, str]],
    seed_lineage: list[dict[str, str]],
    gate_rows: list[dict[str, str]],
    junction_rows: list[dict[str, str]],
    graph_build_rows: list[dict[str, str]],
    graph_mutation_rows: list[dict[str, str]],
    prediction_rows: list[dict[str, str]],
    route_rows: list[dict[str, str]],
    junction_mutation_rows: list[dict[str, str]],
    flow_trace_rows: list[dict[str, str]],
    parity_checktrf_seeds: list[dict[str, str]],
    parity_checktrf_batches: list[dict[str, str]],
    parity_checktrf_events: list[dict[str, str]],
    parity_graph_rows: list[dict[str, str]],
    parity_boundary_rows: list[dict[str, str]],
) -> None:
    kept_tx = sum(1 for row in transcript_support if row["checktrf_mode"])
    output_modes = Counter(row.get("mode", "") for row in outputs.values())
    boundary_hard = sum(
        1
        for row in raw_nodes
        if row.get("hardstart") == "1" or row.get("hardend") == "1"
    )
    reasons = Counter(seed.get("reason", "") for seed in seeds.values())
    rescued_indep = sum(1 for edge in edges if edge["edge_kind"] == "rescued_independent")
    with path.open("w") as handle:
        handle.write(f"# Bundle {current_bundle['bundle_id']} Stage Explanation\n\n")
        handle.write("## Bundle context\n\n")
        handle.write(
            f"- Focused analysis window: `{current_bundle['chrom']}:{current_bundle['start']}-{current_bundle['end']}`"
            f" strand `{current_bundle.get('strand', '')}` with `{current_bundle.get('reads', '')}` reads.\n"
        )
        handle.write(
            f"- Canonical old StringTie bundle from `GGO_19.log`: `{expected_bundle['chrom']}:{expected_bundle['start']}-{expected_bundle['end']}`"
            f" at log line `{expected_bundle['log_line']}` with `{expected_bundle['alignments']}` alignments,"
            f" `{expected_bundle['distinct']}` distinct reads, `{expected_bundle['junctions']}` junctions,"
            f" `{expected_bundle['processed_transcripts']}` processed transcripts.\n"
        )
        handle.write(
            f"- All staged tables below now use the canonical old bundle span for overlap and trace slicing;"
            f" the focused window is retained only as a narrower local view.\n\n"
        )

        handle.write("## Reference context\n\n")
        handle.write(
            f"- Overlapping reference transcripts from `GGO_19.gtf`: `{len(references)}`"
            f" (`{sum(1 for row in references if row['within_current_bundle'] == '1')}` fully inside the canonical bundle).\n"
        )
        gene_count = len({row["gene_id"] for row in references if row["gene_id"]})
        handle.write(f"- Distinct overlapping reference genes: `{gene_count}`.\n\n")

        handle.write("## Layer model\n\n")
        handle.write(
            "- Primary model: typed backward-process layers. The older `G_raw -> G_checktrf -> G_tx` view is kept below only as a coarse projection.\n"
        )
        for row in layer_rows:
            handle.write(
                f"- `{row['layer_id']}`: `{row['count']}` rows across `{row['primary_tables']}`."
                f" {row['detail']}\n"
            )
        handle.write("\n## Three-Graph Projection\n\n")
        handle.write(
            f"- `G_raw`: `{len(raw_transfrags)}` raw transfrags from `TRACE_RAWTRF`,"
            f" with `{len(raw_nodes)}` node-level `SRCSINK_DECISION` rows and `{boundary_hard}` hard-boundary nodes.\n"
        )
        if longtrim:
            handle.write(
                f"- `LONGTRIM` in the old deep trace for the expected bundle: "
                f"`BNODE={longtrim['longtrim_bnode']}`, `BOUND={longtrim['longtrim_bound']}`.\n"
            )
        handle.write(
            f"- `G_checktrf`: `{len(seeds)}` deferred seeds,"
            f" `{len(edges)}` graph edges,"
            f" `{len(outputs)}` output-side events.\n"
        )
        if reasons:
            top_reasons = ", ".join(f"`{reason}`={count}" for reason, count in reasons.most_common(6))
            handle.write(f"- Dominant checktrf seed reasons: {top_reasons}.\n")
        handle.write(
            f"- `G_tx`: `{len(transcripts)}` final transcripts,"
            f" `{kept_tx}` of them linked directly to explicit `checktrf` output rows,"
            f" `{output_modes.get('keep', 0)}` keep-style outputs,"
            f" `{output_modes.get('rescued_independent', 0)}` rescued-independent outputs,"
            f" `{rescued_indep}` independent rescue edges.\n\n"
        )

        if read_pattern_rows or read_events or seed_lineage or junction_rows or graph_build_rows or graph_mutation_rows or prediction_rows or route_rows or junction_mutation_rows or flow_trace_rows or parity_checktrf_seeds or parity_checktrf_batches or parity_checktrf_events or parity_graph_rows or parity_boundary_rows:
            handle.write("## Deep-Trace Lineage\n\n")
            handle.write(
                f"- Bundle-scoped raw deep-trace lineage rows: `{len(read_pattern_rows)}` `PATH_read_pattern` events,"
                f" `{len(read_events)}` `PATH_update_abund` events,"
                f" `{len(read_graph_paths)}` read→graph path summaries,"
                f" `{len(read_to_rawtrf)}` aggregated read→raw-transfrag links,"
                f" `{len(rawtrf_provenance)}` raw-transfrag provenance summaries,"
                f" `{len(seed_lineage)}` seed-decision rows,"
                f" `{len(gate_rows)}` gate rows,"
                f" `{len(junction_rows)}` junction-decision rows,"
                f" `{len(graph_build_rows)}` graph-construction rows,"
                f" `{len(graph_mutation_rows)}` graph-mutation rows,"
                f" `{len(prediction_rows)}` prediction rows,"
                f" `{len(route_rows)}` routing rows,"
                f" `{len(junction_mutation_rows)}` junction-mutation rows,"
                f" `{len(flow_trace_rows)}` low-level flow rows,"
                f" `{len(parity_graph_rows)}` parity graph rows,"
                f" `{len(parity_boundary_rows)}` parity boundary-check rows,"
                f" `{len(parity_checktrf_seeds)}` parity seed rows,"
                f" `{len(parity_checktrf_batches)}` parity checktrf batches,"
                f" `{len(parity_checktrf_events)}` parity checktrf events.\n"
            )
            handle.write(
                "- This adds the missing pointer chain from read-level path updates into graph-local path states and raw transfrags, then into `parse_trflong` seed decisions, then into `checktrf` and final outputs.\n\n"
            )

        handle.write("## Interpretation\n\n")
        handle.write(
            "- The staged data fit a typed layer model better than a strict three-graph model: the intermediate evidence is not just one middle graph, but separate junction, path, raw-transfrag, seed-flow, rescue, and prediction layers.\n"
        )
        handle.write(
            "- Keep the three-graph view as a summary projection only: it is still useful for high-level blocker placement, but it is too coarse for faithful backward replay and late-stage loss attribution.\n"
        )
        handle.write(
            "- The bundle/log context shows the clean calibration bundles are not failing because the old run lacked a proper covering bundle.\n"
        )
        handle.write(
            "- The checktrf tables expose exactly the fields that matter for parity work: `hardstart`, `hardend`, `longstart`, `longend`, rescue reason, and output mapping.\n"
        )
        if read_pattern_rows or read_events or seed_lineage or junction_rows or graph_build_rows or graph_mutation_rows or prediction_rows:
            handle.write(
                "- The added lineage tables make the deeper old-process provenance explicit instead of inferred: read-level `PATH_read_pattern`, read-level `PATH_update_abund`, read→graph path state, raw-transfrag boundary provenance (`longstart`, `longend`, source/sink and hard-boundary evidence), seed expansion/decisions, gate membership, and junction acceptance are now part of the same package.\n"
            )
        if graph_build_rows or graph_mutation_rows or prediction_rows or parity_graph_rows:
            handle.write(
                "- The package now also exposes broader graph-state layers that were previously unused: deep-trace graph construction (`build_graphs`, `create_bundle`, `create_graph`, `merge_read_to_group`), deep-trace graph mutation / trimming (`process_transfrags`, `LONGTRIM_*`, `SRCSINK_DECISION`), parity-side graph evolution (`PARITY_GRAPH_EVOLUTION`, `PARITY_GRAPH_POST`), and prediction fate (`pred[i]`, `print_predcluster`, `PRED_FATE`).\n"
            )
        if route_rows or junction_mutation_rows or flow_trace_rows or parity_boundary_rows:
            handle.write(
                "- Additional replay-oriented layers are now staged explicitly: `ROUTE_*`, `JUNC_DELETE` / `EJUNC_DELETE` / `JUNC_DEMOTE`, low-level `fwd_to_sink` / `back_to_source` step traces, parity-side `PARITY_BOUNDARY_CHECK`, and both deep-trace boundary/tiebreak rows plus parity-trace `TIEBREAK_BACK` / `TIEBREAK_FWD` rows when present.\n"
            )
        if parity_checktrf_seeds or parity_checktrf_batches or parity_checktrf_events:
            handle.write(
                "- `trace_full.log` parity rows are now staged separately, so `checktrf` no longer has to be inferred only from gate side-effects when an explicit parity batch exists for the canonical bundle.\n"
            )
        handle.write(
            "- This package is intended to replace repeated one-off comparisons: once a blocker locus is mapped into a bundle, the same staged tables can explain whether the issue is already present in `G_raw`, introduced by `G_checktrf`, or only visible in `G_tx`.\n"
        )


def main() -> int:
    args = parse_args()
    output_stem = Path(args.output_stem)

    current_bundle = read_bundle_summary(args.bundle_summary, args.bundle_id)
    bundles = load_bundle_log(args.bundle_log)
    expected_bundle = choose_expected_bundle(current_bundle, bundles)

    chrom = expected_bundle["chrom"]
    start = parse_int(expected_bundle["start"])
    end = parse_int(expected_bundle["end"])
    focused_start = parse_int(current_bundle["start"])
    focused_end = parse_int(current_bundle["end"])

    references = load_reference_overlaps(args.reference_gtf, chrom, start, end)
    raw_nodes = load_legacy_nodes(args.legacy_trace_norm)
    raw_transfrags = load_raw_transfrags(args.rust_trace_norm)
    seeds, edges, outputs = parse_trace_records(args.rust_trace_norm)
    transcripts = load_transcripts(args.old_transcripts)
    transcript_support = build_transcript_support(transcripts, outputs)
    longtrim = scan_longtrim(args.legacy_trace_log, expected_bundle) if args.legacy_trace_log else {}
    read_pattern_rows, read_events, read_to_rawtrf, seed_lineage, gate_rows, junction_rows = parse_legacy_lineage(
        args.legacy_trace_log,
        args.legacy_trace_norm,
        start,
        end,
    )
    graph_build_rows, graph_mutation_rows, prediction_rows, prediction_edge_rows, junction_trace_rows, route_rows, junction_mutation_rows, flow_trace_rows = parse_deep_bundle_layers(
        args.legacy_trace_log,
        args.legacy_trace_norm,
        start,
        end,
    )
    (
        parity_checktrf_seeds,
        parity_checktrf_batches,
        parity_checktrf_events,
        parity_tiebreak_rows,
        parity_graph_rows,
        parity_boundary_rows,
    ) = parse_parity_checktrf_layers(
        args.parity_trace_log,
        start,
        end,
    )
    flow_trace_rows.extend(parity_tiebreak_rows)
    read_pattern_paths = summarize_read_pattern_paths(read_pattern_rows)
    read_graph_paths = summarize_read_graph_paths(read_events)
    rawtrf_provenance = summarize_rawtrf_provenance(read_to_rawtrf, read_graph_paths)
    layer_rows = build_layer_rows(
        raw_nodes,
        raw_transfrags,
        read_pattern_rows,
        read_events,
        read_graph_paths,
        read_to_rawtrf,
        rawtrf_provenance,
        seed_lineage,
        gate_rows,
        junction_rows,
        graph_build_rows,
        graph_mutation_rows,
        prediction_rows,
        prediction_edge_rows,
        junction_trace_rows,
        route_rows,
        junction_mutation_rows,
        flow_trace_rows,
        seeds,
        edges,
        outputs,
        transcripts,
        transcript_support,
        parity_checktrf_seeds,
        parity_checktrf_batches,
        parity_checktrf_events,
        parity_graph_rows,
        parity_boundary_rows,
    )

    bundle_context_rows = [
        {
            "scope": "focused_window",
            "bundle_id": current_bundle.get("bundle_id", ""),
            "chrom": current_bundle.get("chrom", ""),
            "strand": current_bundle.get("strand", ""),
            "start": current_bundle.get("start", ""),
            "end": current_bundle.get("end", ""),
            "reads": current_bundle.get("reads", ""),
            "alignments": "",
            "distinct": "",
            "junctions": "",
            "processed_transcripts": "",
            "log_line": "",
        },
        {
            "scope": "canonical_bundle",
            "bundle_id": current_bundle.get("bundle_id", ""),
            "chrom": expected_bundle.get("chrom", ""),
            "strand": current_bundle.get("strand", ""),
            "start": expected_bundle.get("start", ""),
            "end": expected_bundle.get("end", ""),
            "reads": "",
            "alignments": expected_bundle.get("alignments", ""),
            "distinct": expected_bundle.get("distinct", ""),
            "junctions": expected_bundle.get("junctions", ""),
            "processed_transcripts": expected_bundle.get("processed_transcripts", ""),
            "log_line": expected_bundle.get("log_line", ""),
        },
    ]

    stage_rows = [
        {"stage": "reference_overlap", "count": str(len(references)), "detail": "transcripts overlapping canonical old-log bundle"},
        {"stage": "rawgraph_nodes", "count": str(len(raw_nodes)), "detail": "legacy SRCSINK_DECISION rows"},
        {"stage": "rawgraph_transfrags", "count": str(len(raw_transfrags)), "detail": "TRACE_RAWTRF rows"},
        {"stage": "checktrf_seeds", "count": str(len(seeds)), "detail": "TRACE_CHECKTRF_ADD rows"},
        {"stage": "checktrf_edges", "count": str(len(edges)), "detail": "derived seed/keep/output edges"},
        {"stage": "checktrf_outputs", "count": str(len(outputs)), "detail": "derived output-side rows"},
        {"stage": "transcripts", "count": str(len(transcripts)), "detail": "final old-path transcripts"},
    ]
    if longtrim:
        stage_rows.extend(
            [
                {"stage": "longtrim_bnode", "count": longtrim["longtrim_bnode"], "detail": "LONGTRIM_BNODE rows in old deep log"},
                {"stage": "longtrim_bound", "count": longtrim["longtrim_bound"], "detail": "LONGTRIM_BOUND rows in old deep log"},
            ]
        )
    if (
        read_pattern_rows
        or read_pattern_paths
        or read_events
        or read_to_rawtrf
        or seed_lineage
        or gate_rows
        or junction_rows
        or graph_build_rows
        or graph_mutation_rows
        or prediction_rows
        or prediction_edge_rows
        or junction_trace_rows
        or route_rows
        or junction_mutation_rows
        or flow_trace_rows
        or parity_graph_rows
        or parity_boundary_rows
        or parity_checktrf_seeds
        or parity_checktrf_batches
        or parity_checktrf_events
    ):
        stage_rows.extend(
            [
                {"stage": "lineage_read_pattern_rows", "count": str(len(read_pattern_rows)), "detail": "bundle-scoped PATH_read_pattern rows from old deep log"},
                {"stage": "lineage_read_pattern_paths", "count": str(len(read_pattern_paths)), "detail": "aggregated read-pattern path rows"},
                {"stage": "lineage_read_events", "count": str(len(read_events)), "detail": "bundle-scoped PATH_update_abund rows from old deep log"},
                {"stage": "lineage_read_graph_paths", "count": str(len(read_graph_paths)), "detail": "aggregated read to graph-path state rows"},
                {"stage": "lineage_read_to_rawtrf", "count": str(len(read_to_rawtrf)), "detail": "aggregated read to raw-transfrag links"},
                {"stage": "lineage_rawtrf_provenance", "count": str(len(rawtrf_provenance)), "detail": "aggregated raw-transfrag boundary provenance"},
                {"stage": "lineage_seed_decisions", "count": str(len(seed_lineage)), "detail": "bundle-scoped parse_trflong seed rows with decisions"},
                {"stage": "lineage_gate_rows", "count": str(len(gate_rows)), "detail": "bundle-scoped CHECKTRF gate rows"},
                {"stage": "lineage_junction_rows", "count": str(len(junction_rows)), "detail": "bundle-scoped junction decision rows"},
                {"stage": "lineage_graph_build_rows", "count": str(len(graph_build_rows)), "detail": "bundle-scoped graph construction rows from old deep log"},
                {"stage": "lineage_graph_mutation_rows", "count": str(len(graph_mutation_rows)), "detail": "bundle-scoped graph mutation rows from old deep log"},
                {"stage": "lineage_prediction_rows", "count": str(len(prediction_rows)), "detail": "bundle-scoped prediction / fate rows from old deep log"},
                {"stage": "lineage_prediction_edge_rows", "count": str(len(prediction_edge_rows)), "detail": "structured killer-vs-killed prediction edges from old deep log"},
                {"stage": "lineage_junction_trace_rows", "count": str(len(junction_trace_rows)), "detail": "deep-trace junction build / decision rows from old deep log"},
                {"stage": "lineage_route_rows", "count": str(len(route_rows)), "detail": "build_graphs route rows from old deep log"},
                {"stage": "lineage_junction_mutation_rows", "count": str(len(junction_mutation_rows)), "detail": "JUNC_DELETE / EJUNC_DELETE / JUNC_DEMOTE rows from old deep log"},
                {"stage": "lineage_flow_trace_rows", "count": str(len(flow_trace_rows)), "detail": "low-level fwd/back step and boundary/tiebreak rows from old deep log plus parity-trace tiebreak rows"},
                {"stage": "lineage_parity_graph_rows", "count": str(len(parity_graph_rows)), "detail": "PARITY_GRAPH_EVOLUTION / PARITY_GRAPH_POST rows from trace_full.log for the canonical bundle"},
                {"stage": "lineage_parity_boundary_check_rows", "count": str(len(parity_boundary_rows)), "detail": "PARITY_BOUNDARY_CHECK rows from trace_full.log overlapping the canonical bundle"},
                {"stage": "lineage_parity_checktrf_seed_rows", "count": str(len(parity_checktrf_seeds)), "detail": "PARITY_SEED_PROC rows from trace_full.log for the canonical bundle"},
                {"stage": "lineage_parity_checktrf_batch_rows", "count": str(len(parity_checktrf_batches)), "detail": "PARITY_CHECKTRF_START rows or synthetic PARITY_CHK-run batches from trace_full.log for the canonical bundle"},
                {"stage": "lineage_parity_checktrf_event_rows", "count": str(len(parity_checktrf_events)), "detail": "PARITY_CHK / PARITY_CHECKTRF_OUTCOME rows from trace_full.log for the canonical bundle"},
            ]
        )

    write_tsv(
        output_stem.with_suffix(".bundle_context.tsv"),
        bundle_context_rows,
        [
            "scope",
            "bundle_id",
            "chrom",
            "strand",
            "start",
            "end",
            "reads",
            "alignments",
            "distinct",
            "junctions",
            "processed_transcripts",
            "log_line",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".reference.tsv"),
        references,
        ["gene_id", "transcript_id", "chrom", "strand", "start", "end", "num_exons", "within_current_bundle"],
    )
    write_tsv(
        output_stem.with_suffix(".layers.tsv"),
        layer_rows,
        ["layer_id", "count", "primary_tables", "detail"],
    )
    write_tsv(
        output_stem.with_suffix(".rawgraph.nodes.tsv"),
        raw_nodes,
        ["node_id", "start", "end", "hardstart", "hardend", "source_state", "sink_state", "raw"],
    )
    write_tsv(
        output_stem.with_suffix(".rawgraph.transfrags.tsv"),
        raw_transfrags,
        [
            "tf_idx",
            "start",
            "end",
            "n_nodes",
            "abundance",
            "read_count",
            "srabund",
            "longstart",
            "longend",
            "guide",
            "longread",
            "real",
            "weak",
            "coverage_weak",
            "trflong_seed",
            "usepath",
            "group_size",
            "junction_chain",
            "node_ids",
            "hardstart",
            "hardend",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".checktrf.seeds.tsv"),
        seeds.values(),
        [
            "seed_id",
            "reason",
            "saved_abund",
            "abundance",
            "guide",
            "longstart",
            "longend",
            "back",
            "fwd",
            "flux",
            "cov",
            "tocheck",
            "tf_node_count",
            "tf_coords",
            "path_node_count",
            "path_coords",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".checktrf.edges.tsv"),
        edges,
        ["source_kind", "source_id", "target_kind", "target_id", "edge_kind", "reason", "coverage", "out_idx", "raw"],
    )
    write_tsv(
        output_stem.with_suffix(".checktrf.outputs.tsv"),
        outputs.values(),
        ["output_key", "out_idx", "seed_ids", "keep_ids", "mode", "coverage", "tx_cov", "guide", "start", "end", "exon_count", "exons", "raw"],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.read_patterns.tsv"),
        read_pattern_rows,
        [
            "line_no",
            "event",
            "read_idx",
            "s",
            "graph_idx",
            "node_field",
            "node_id",
            "node_start",
            "node_end",
            "bp",
            "total_nodes",
            "k",
            "ncoord",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.read_pattern_paths.tsv"),
        read_pattern_paths,
        [
            "read_idx",
            "s",
            "graph_idx",
            "first_line",
            "last_line",
            "event_count",
            "events",
            "add_node_count",
            "intersect_count",
            "exon_skip_count",
            "total_nodes_max",
            "node_ids",
            "node_spans",
            "skip_bnodes",
            "max_bp",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.read_events.tsv"),
        read_events,
        [
            "line_no",
            "phase",
            "event",
            "read_idx",
            "s",
            "g",
            "gno",
            "trf_idx",
            "abundance",
            "total",
            "rstart",
            "rend",
            "longstart",
            "longend",
            "nnodes",
            "node_ids",
            "node",
            "node0_start",
            "nodeLast_end",
            "node0_hardstart",
            "nodeLast_hardend",
            "is_lr",
            "is_sr",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.read_graph_paths.tsv"),
        read_graph_paths,
        [
            "read_idx",
            "s",
            "g",
            "gno",
            "first_line",
            "last_line",
            "event_count",
            "events",
            "trf_idx",
            "rstart",
            "rend",
            "node0_start",
            "nodeLast_end",
            "node0_hardstart",
            "nodeLast_hardend",
            "is_source",
            "is_sink",
            "nnodes",
            "node_ids",
            "longstart",
            "longend",
            "source_check_entered",
            "sink_check_entered",
            "skip_single_node",
            "new_trf_seen",
            "final_nodes_seen",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.read_to_rawtrf.tsv"),
        read_to_rawtrf,
        [
            "read_idx",
            "s",
            "g",
            "trf_idx",
            "first_line",
            "last_line",
            "event_count",
            "events",
            "abundance",
            "total",
            "rstart",
            "rend",
            "longstart",
            "longend",
            "nnodes",
            "node_ids",
            "is_lr",
            "is_sr",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.rawtrf_provenance.tsv"),
        rawtrf_provenance,
        [
            "trf_idx",
            "read_count",
            "graph_count",
            "graphs",
            "first_line",
            "last_line",
            "events",
            "abundance_sum",
            "total_max",
            "min_rstart",
            "max_rend",
            "min_longstart_nonzero",
            "max_longend_nonzero",
            "any_source",
            "any_sink",
            "any_hardstart",
            "any_hardend",
            "source_check_reads",
            "sink_check_reads",
            "skip_single_node_reads",
            "node_ids_examples",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.seed_decisions.tsv"),
        seed_lineage,
        [
            "seed_id",
            "line_seed",
            "abundance",
            "nodes",
            "guide",
            "longread",
            "longstart",
            "longend",
            "node_coords",
            "back_result",
            "back_path_len",
            "back_path",
            "fwd_result",
            "fwd_path_len",
            "fwd_path",
            "used_direct",
            "branch",
            "reason",
            "back",
            "fwd",
            "flux",
            "cov",
            "len",
            "checktrf_gate_count",
            "checktrf_gate_keeptrf",
            "checktrf_gate_readthr",
            "gate_pass",
            "raw_seed",
            "raw_decision",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.checktrf_gate.tsv"),
        gate_rows,
        ["line_no", "line_gate", "gate_idx", "seed_id", "abundance", "guide", "passed", "gate_count", "keeptrf", "readthr", "raw"],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.junction_decisions.tsv"),
        junction_rows,
        ["line_no", "event", "junction_idx", "junction", "strand", "nm", "mm", "nreads", "nreads_good", "left", "right", "guide", "reason", "raw"],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.graph_build.tsv"),
        graph_build_rows,
        [
            "line_no",
            "source",
            "event",
            "span_start",
            "span_end",
            "bundle_start",
            "bundle_end",
            "pred_idx",
            "cov",
            "readcov",
            "exons",
            "strand",
            "guide",
            "fate",
            "node_id",
            "bno",
            "bid",
            "grid",
            "color",
            "reason",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.graph_mutation.tsv"),
        graph_mutation_rows,
        [
            "line_no",
            "source",
            "event",
            "span_start",
            "span_end",
            "bundle_start",
            "bundle_end",
            "pred_idx",
            "cov",
            "readcov",
            "exons",
            "strand",
            "guide",
            "fate",
            "node_id",
            "bno",
            "bid",
            "grid",
            "color",
            "reason",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.predictions.tsv"),
        prediction_rows,
        [
            "line_no",
            "source",
            "event",
            "span_start",
            "span_end",
            "bundle_start",
            "bundle_end",
            "pred_idx",
            "cov",
            "readcov",
            "exons",
            "strand",
            "guide",
            "fate",
            "node_id",
            "bno",
            "bid",
            "grid",
            "color",
            "reason",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.prediction_edges.tsv"),
        prediction_edge_rows,
        [
            "line_no",
            "source",
            "event",
            "span_start",
            "span_end",
            "bundle_start",
            "bundle_end",
            "pred_idx",
            "cov",
            "readcov",
            "exons",
            "strand",
            "guide",
            "fate",
            "node_id",
            "bno",
            "bid",
            "grid",
            "color",
            "reason",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.junction_trace.tsv"),
        junction_trace_rows,
        [
            "line_no",
            "source",
            "event",
            "span_start",
            "span_end",
            "bundle_start",
            "bundle_end",
            "pred_idx",
            "cov",
            "readcov",
            "exons",
            "strand",
            "guide",
            "fate",
            "node_id",
            "bno",
            "bid",
            "grid",
            "color",
            "reason",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.routes.tsv"),
        route_rows,
        [
            "line_no",
            "source",
            "event",
            "span_start",
            "span_end",
            "bundle_start",
            "bundle_end",
            "pred_idx",
            "cov",
            "readcov",
            "exons",
            "strand",
            "guide",
            "fate",
            "node_id",
            "bno",
            "bid",
            "grid",
            "color",
            "reason",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.junction_mutation.tsv"),
        junction_mutation_rows,
        [
            "line_no",
            "source",
            "event",
            "span_start",
            "span_end",
            "bundle_start",
            "bundle_end",
            "pred_idx",
            "cov",
            "readcov",
            "exons",
            "strand",
            "guide",
            "fate",
            "node_id",
            "bno",
            "bid",
            "grid",
            "color",
            "reason",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.flow_trace.tsv"),
        flow_trace_rows,
        [
            "line_no",
            "source",
            "event",
            "span_start",
            "span_end",
            "bundle_start",
            "bundle_end",
            "pred_idx",
            "cov",
            "readcov",
            "exons",
            "strand",
            "guide",
            "fate",
            "node_id",
            "bno",
            "bid",
            "grid",
            "color",
            "reason",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.parity_graph.tsv"),
        parity_graph_rows,
        [
            "line_no",
            "row_kind",
            "sno",
            "bundle_idx",
            "graph_idx",
            "bundle_start",
            "bundle_end",
            "event",
            "node_id",
            "span_start",
            "span_end",
            "bid",
            "parents",
            "nodes",
            "transfrags",
            "active_transfrags",
            "trflong",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.parity_boundary_checks.tsv"),
        parity_boundary_rows,
        [
            "line_no",
            "bundle_start",
            "bundle_end",
            "span_start",
            "span_end",
            "nodes",
            "first",
            "last",
            "is_source",
            "is_sink",
            "hardstart",
            "hardend",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.parity_checktrf_seeds.tsv"),
        parity_checktrf_seeds,
        [
            "line_no",
            "bundle_start",
            "bundle_end",
            "seed_idx",
            "span_start",
            "span_end",
            "abundance",
            "nodes",
            "back",
            "fwd",
            "flux",
            "path",
            "outcome",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.parity_checktrf_batches.tsv"),
        parity_checktrf_batches,
        [
            "line_no",
            "batch_id",
            "batch_kind",
            "bundle_start",
            "bundle_end",
            "keeptrf",
            "checktrf",
            "raw",
        ],
    )
    write_tsv(
        output_stem.with_suffix(".lineage.parity_checktrf_events.tsv"),
        parity_checktrf_events,
        [
            "line_no",
            "batch_id",
            "event",
            "seed_idx",
            "abundance",
            "outcome",
            "matched",
            "abundsum",
            "cov",
            "nexons",
            "exons",
            "matches",
            "nodes",
            "raw",
        ],
    )
    shutil.copyfile(args.old_transcripts, output_stem.with_suffix(".transcripts.tsv"))
    write_tsv(
        output_stem.with_suffix(".transcript_support.tsv"),
        transcript_support,
        [
            "tx_idx",
            "start",
            "end",
            "num_exons",
            "coverage",
            "longcov",
            "junction_chain",
            "checktrf_mode",
            "checktrf_seed_ids",
            "checktrf_keep_ids",
            "checktrf_cov",
            "checktrf_tx_cov",
        ],
    )
    write_tsv(output_stem.with_suffix(".stages.tsv"), stage_rows, ["stage", "count", "detail"])
    write_report(
        output_stem.with_suffix(".explain.md"),
        current_bundle,
        expected_bundle,
        references,
        layer_rows,
        raw_nodes,
        raw_transfrags,
        seeds,
        edges,
        outputs,
        transcripts,
        transcript_support,
        longtrim,
        read_pattern_rows,
        read_events,
        read_to_rawtrf,
        read_graph_paths,
        rawtrf_provenance,
        seed_lineage,
        gate_rows,
        junction_rows,
        graph_build_rows,
        graph_mutation_rows,
        prediction_rows,
        route_rows,
        junction_mutation_rows,
        flow_trace_rows,
        parity_checktrf_seeds,
        parity_checktrf_batches,
        parity_checktrf_events,
        parity_graph_rows,
        parity_boundary_rows,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
