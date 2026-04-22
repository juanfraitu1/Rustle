#!/usr/bin/env python3
"""Build a locus-scoped backward-vs-current oracle table for annotated transcripts.

This is the locus analogue of the bundle backward-oracle tooling used for bundle 904.
It compares:

- legacy deep-trace junction behavior (`jdec`, `good_junc`, bad/delete markers)
- legacy final `PRED_FATE` spans from a normalized backward trace TSV
- annotated reference transcripts for one locus
- current best query matches from the divergence table + current query GTF
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict


TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')
JUNCTION_CHAIN_RE = re.compile(r'junction_chain "([^"]*)"')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-log", required=True)
    parser.add_argument("--backward-trace-tsv", required=True)
    parser.add_argument("--reference-gtf", required=True)
    parser.add_argument("--current-query-gtf", required=True)
    parser.add_argument("--divergence-tsv", required=True)
    parser.add_argument("--output-stem", required=True)
    return parser.parse_args()


@dataclass
class TranscriptRecord:
    transcript_id: str
    start: int
    end: int
    exon_count: int
    junctions: list[tuple[int, int]]

    @property
    def span(self) -> str:
        return f"{self.start}-{self.end}"


def parse_gtf(path: Path) -> dict[str, TranscriptRecord]:
    exon_map: DefaultDict[str, list[tuple[int, int]]] = defaultdict(list)
    tx_spans: dict[str, tuple[int, int]] = {}
    tx_chains: dict[str, list[tuple[int, int]]] = {}

    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = line.rstrip("\n").split("\t")
            match = TRANSCRIPT_ID_RE.search(attrs)
            if not match:
                continue
            tid = match.group(1)
            if feature == "transcript":
                tx_spans[tid] = (int(start), int(end))
                chain_match = JUNCTION_CHAIN_RE.search(attrs)
                if chain_match and chain_match.group(1):
                    junctions = []
                    for token in chain_match.group(1).split(","):
                        donor, acceptor = token.split("-")
                        junctions.append((int(donor), int(acceptor)))
                    tx_chains[tid] = junctions
            elif feature == "exon":
                exon_map[tid].append((int(start), int(end)))

    out: dict[str, TranscriptRecord] = {}
    for tid, exons in exon_map.items():
        ordered = sorted(exons)
        start, end = tx_spans.get(tid, (ordered[0][0], ordered[-1][1]))
        junctions = tx_chains.get(
            tid,
            [(ordered[i][1], ordered[i + 1][0]) for i in range(len(ordered) - 1)],
        )
        out[tid] = TranscriptRecord(
            transcript_id=tid,
            start=start,
            end=end,
            exon_count=len(ordered),
            junctions=junctions,
        )
    return out


def load_divergence(path: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            out[row["ref_id"]] = row
    return out


def load_backward_pred_fate(path: Path) -> tuple[dict[str, list[dict[str, str]]], dict[str, list[dict[str, str]]]]:
    by_span: DefaultDict[str, list[dict[str, str]]] = defaultdict(list)
    by_span_exons: DefaultDict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["event"] != "PRED_FATE":
                continue
            span = f'{row["span_start"]}-{row["span_end"]}'
            by_span[span].append(row)
            key = f'{span}|{row["exon_count"]}'
            by_span_exons[key].append(row)
    return dict(by_span), dict(by_span_exons)


def collect_trace_metrics(path: Path, target_junctions: set[str]) -> dict[str, dict[str, int | str]]:
    metrics: dict[str, dict[str, int | str]] = {
        junction: {
            "good_junc_accept": 0,
            "jdec_good": 0,
            "bad_junc": 0,
            "junc_delete": 0,
            "ejunc_delete": 0,
            "example_accept": "",
            "example_bad": "",
        }
        for junction in target_junctions
    }
    with path.open(errors="replace") as handle:
        for line_no, line in enumerate(handle, 1):
            text = line.rstrip("\n")
            for junction in target_junctions:
                if junction not in text:
                    continue
                entry = metrics[junction]
                if "good_junc: ACCEPTED" in text:
                    entry["good_junc_accept"] += 1
                    if not entry["example_accept"]:
                        entry["example_accept"] = f"{line_no}:{text}"
                if "reason=GOOD" in text and "jdec[" in text:
                    entry["jdec_good"] += 1
                if "BAD_JUNC" in text:
                    entry["bad_junc"] += 1
                    if not entry["example_bad"]:
                        entry["example_bad"] = f"{line_no}:{text}"
                if "EJUNC_DELETE" in text:
                    entry["ejunc_delete"] += 1
                elif "JUNC_DELETE" in text:
                    entry["junc_delete"] += 1
    return metrics


def format_junctions(junctions: list[tuple[int, int]]) -> str:
    return ",".join(f"{a}-{b}" for a, b in junctions)


def format_metric_summary(junctions: list[tuple[int, int]], metrics: dict[str, dict[str, int | str]]) -> str:
    parts = []
    for donor, acceptor in junctions:
        key = f"{donor}-{acceptor}"
        metric = metrics.get(key)
        if not metric:
            continue
        parts.append(
            f"{key}:accept={metric['good_junc_accept']},jdec={metric['jdec_good']},bad={metric['bad_junc']},del={metric['junc_delete']},edel={metric['ejunc_delete']}"
        )
    return ";".join(parts)


def verdict_for_row(
    backward_exact_span_rows: list[dict[str, str]],
    backward_exact_span_exon_rows: list[dict[str, str]],
    current_classification: str,
    ref_unique_metrics: list[tuple[str, dict[str, int | str]]],
    query_unique_metrics: list[tuple[str, dict[str, int | str]]],
) -> str:
    if backward_exact_span_exon_rows and current_classification == "exact_chain":
        return "boundary_only_current_miss"
    if backward_exact_span_rows and current_classification == "internal_divergence":
        ref_good = sum(int(metric["good_junc_accept"]) for _, metric in ref_unique_metrics)
        query_good = sum(int(metric["good_junc_accept"]) for _, metric in query_unique_metrics)
        if ref_good > 0 and query_good == 0:
            return "canonical_ref_family_preferred_backward"
        if ref_good >= max(4, query_good * 2):
            return "canonical_ref_family_preferred_backward"
        return "backward_keeps_family_but_not_clear_junction_preference"
    if backward_exact_span_rows:
        return "backward_keeps_exact_span"
    return "no_exact_backward_span_match"


def main() -> int:
    args = parse_args()
    reference = parse_gtf(Path(args.reference_gtf))
    current_query = parse_gtf(Path(args.current_query_gtf))
    divergence = load_divergence(Path(args.divergence_tsv))
    backward_by_span, backward_by_span_exons = load_backward_pred_fate(Path(args.backward_trace_tsv))

    target_junctions: set[str] = set()
    per_ref_pairs: dict[str, tuple[list[tuple[int, int]], list[tuple[int, int]]]] = {}
    for ref_id, ref_tx in reference.items():
        div = divergence.get(ref_id)
        if not div:
            per_ref_pairs[ref_id] = ([], [])
            continue
        query_tx = current_query.get(div["best_query_id"])
        if not query_tx:
            per_ref_pairs[ref_id] = ([], [])
            continue
        ref_unique = [junction for junction in ref_tx.junctions if junction not in query_tx.junctions]
        query_unique = [junction for junction in query_tx.junctions if junction not in ref_tx.junctions]
        per_ref_pairs[ref_id] = (ref_unique, query_unique)
        target_junctions.update(f"{a}-{b}" for a, b in ref_unique)
        target_junctions.update(f"{a}-{b}" for a, b in query_unique)

    trace_metrics = collect_trace_metrics(Path(args.trace_log), target_junctions)

    ref_rows: list[dict[str, str]] = []
    junction_rows: list[dict[str, str]] = []
    for junction in sorted(target_junctions, key=lambda item: tuple(int(x) for x in item.split("-"))):
        metric = trace_metrics[junction]
        junction_rows.append(
            {
                "junction": junction,
                "good_junc_accept": str(metric["good_junc_accept"]),
                "jdec_good": str(metric["jdec_good"]),
                "bad_junc": str(metric["bad_junc"]),
                "junc_delete": str(metric["junc_delete"]),
                "ejunc_delete": str(metric["ejunc_delete"]),
                "example_accept": str(metric["example_accept"]),
                "example_bad": str(metric["example_bad"]),
            }
        )

    for ref_id, ref_tx in sorted(reference.items()):
        div = divergence.get(ref_id, {})
        current_id = div.get("best_query_id", "")
        current_tx = current_query.get(current_id)
        ref_unique, query_unique = per_ref_pairs.get(ref_id, ([], []))
        ref_unique_metrics = [(f"{a}-{b}", trace_metrics[f"{a}-{b}"]) for a, b in ref_unique]
        query_unique_metrics = [(f"{a}-{b}", trace_metrics[f"{a}-{b}"]) for a, b in query_unique]

        span_rows = backward_by_span.get(ref_tx.span, [])
        span_exon_rows = backward_by_span_exons.get(f"{ref_tx.span}|{ref_tx.exon_count}", [])
        verdict = verdict_for_row(
            span_rows,
            span_exon_rows,
            div.get("classification", ""),
            ref_unique_metrics,
            query_unique_metrics,
        )

        ref_rows.append(
            {
                "ref_id": ref_id,
                "ref_span": ref_tx.span,
                "ref_exons": str(ref_tx.exon_count),
                "ref_junctions": str(len(ref_tx.junctions)),
                "backward_exact_span_kept": "1" if span_rows else "0",
                "backward_exact_span_exon_kept": "1" if span_exon_rows else "0",
                "backward_kept_covs": ",".join(row["coverage"] for row in span_rows),
                "backward_kept_fates": ",".join(row["decision"] for row in span_rows),
                "current_best_query_id": current_id,
                "current_best_query_span": current_tx.span if current_tx else "",
                "current_best_query_exons": str(current_tx.exon_count) if current_tx else "",
                "current_classification": div.get("classification", ""),
                "first_divergence_ref": div.get("ref_first_divergence", ""),
                "first_divergence_query": div.get("query_first_divergence", ""),
                "ref_unique_junctions": format_junctions(ref_unique),
                "query_unique_junctions": format_junctions(query_unique),
                "ref_unique_trace_summary": format_metric_summary(ref_unique, trace_metrics),
                "query_unique_trace_summary": format_metric_summary(query_unique, trace_metrics),
                "oracle_verdict": verdict,
            }
        )

    stem = Path(args.output_stem)
    refs_path = stem.with_suffix(".refs.tsv")
    junctions_path = stem.with_suffix(".junctions.tsv")

    with refs_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "ref_id",
                "ref_span",
                "ref_exons",
                "ref_junctions",
                "backward_exact_span_kept",
                "backward_exact_span_exon_kept",
                "backward_kept_covs",
                "backward_kept_fates",
                "current_best_query_id",
                "current_best_query_span",
                "current_best_query_exons",
                "current_classification",
                "first_divergence_ref",
                "first_divergence_query",
                "ref_unique_junctions",
                "query_unique_junctions",
                "ref_unique_trace_summary",
                "query_unique_trace_summary",
                "oracle_verdict",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(ref_rows)

    with junctions_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "junction",
                "good_junc_accept",
                "jdec_good",
                "bad_junc",
                "junc_delete",
                "ejunc_delete",
                "example_accept",
                "example_bad",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(junction_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
