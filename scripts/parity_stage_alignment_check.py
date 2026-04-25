#!/usr/bin/env python3
"""
Join StringTie PARITY_STAGE_TSV output with Rustle --debug-stage-tsv and report
per-stage alignment. Optionally verify Rustle parity_shadow internal diff rates.

StringTie coordinates are treated as 1-based closed intervals matching GTF;
Rustle bundle coordinates are matched as (chrom, st_start - 1, st_end - 1).

Usage:
  python3 parity_stage_alignment_check.py \\
    --stringtie stringtie_parity.tsv \\
    --rustle rustle.stage.tsv \\
    [--parity-shadow rustle.parity_shadow.tsv] \\
    [--strict] [--json report.json]

Exit code: 0 if all checks pass under tolerances, else 1 (--strict) or 0 with warnings.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from dataclasses import dataclass, asdict
from typing import Any


def st_to_rust_key(chrom: str, start: int, end: int) -> tuple[str, int, int]:
    return (chrom, start - 1, end - 1)


def load_stringtie_by_stage(path: str) -> dict[str, dict[tuple[str, int, int], dict[str, str]]]:
    by_stage: dict[str, dict[tuple[str, int, int], dict[str, str]]] = defaultdict(dict)
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            stage = row["stage"]
            chrom = row["chrom"]
            start = int(row["start"])
            end = int(row["end"])
            k = st_to_rust_key(chrom, start, end)
            by_stage[stage][k] = row
    return by_stage


def load_rustle_by_stage(path: str) -> dict[str, list[dict[str, str]]]:
    by_stage: dict[str, list[dict[str, str]]] = defaultdict(list)
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            by_stage[row["stage_id"]].append(row)
    return by_stage


def i(row: dict[str, str], key: str) -> int:
    return int(row[key])


def median(xs: list[float]) -> float:
    if not xs:
        return 0.0
    s = sorted(xs)
    m = len(s) // 2
    if len(s) % 2:
        return float(s[m])
    return (s[m - 1] + s[m]) / 2.0


@dataclass
class CheckResult:
    name: str
    ok: bool
    bundles_checked: int
    mismatches: int
    detail: dict[str, Any]


def check_reads(
    st: dict[str, Any], ru: dict[str, Any], max_abs_delta: int
) -> CheckResult:
    """bundle_input.n_reads vs bundle_context.n_reads"""
    st_bundles = st["bundle_input"]
    ru_rows = {(r["chrom"], int(r["start"]), int(r["end"])): r for r in ru["bundle_context"]}
    mismatches = 0
    n = 0
    deltas = []
    for k, srow in st_bundles.items():
        rrow = ru_rows.get(k)
        if not rrow:
            continue
        n += 1
        a, b = i(srow, "n_reads"), i(rrow, "n_reads")
        d = abs(a - b)
        deltas.append(float(d))
        if d > max_abs_delta:
            mismatches += 1
    return CheckResult(
        name="reads_bundle_input_vs_bundle_context",
        ok=mismatches == 0,
        bundles_checked=n,
        mismatches=mismatches,
        detail={
            "max_abs_delta_allowed": max_abs_delta,
            "median_abs_delta": median(deltas),
            "max_abs_delta": max(deltas) if deltas else 0,
        },
    )


def check_junctions_kept(
    st: dict[str, Any],
    ru: dict[str, Any],
    max_abs_delta: int,
) -> CheckResult:
    """post_junction_filter.n_junctions_kept vs sum(graph_evolution.n_junctions) per bundle."""
    st_pf = st["post_junction_filter"]
    evo = ru["graph_evolution"]
    sum_j: dict[tuple[str, int, int], int] = defaultdict(int)
    nrows: dict[tuple[str, int, int], int] = defaultdict(int)
    for r in evo:
        k = (r["chrom"], int(r["start"]), int(r["end"]))
        sum_j[k] += i(r, "n_junctions")
        nrows[k] += 1
    mismatches = 0
    deltas = []
    n = 0
    for k, srow in st_pf.items():
        if k not in sum_j:
            continue
        n += 1
        a = i(srow, "n_junctions_kept")
        b = sum_j[k]
        d = abs(a - b)
        deltas.append(float(d))
        if d > max_abs_delta:
            mismatches += 1
    single = sum(1 for k in st_pf if k in sum_j and nrows[k] == 1)
    single_deltas = []
    for k, srow in st_pf.items():
        if k not in sum_j or nrows[k] != 1:
            continue
        single_deltas.append(
            abs(i(srow, "n_junctions_kept") - sum_j[k])
        )
    return CheckResult(
        name="junctions_kept_vs_graph_evolution_sum",
        ok=mismatches == 0,
        bundles_checked=n,
        mismatches=mismatches,
        detail={
            "tolerance_max_abs_delta": max_abs_delta,
            "median_abs_delta": median(deltas),
            "max_abs_delta": max(deltas) if deltas else 0,
            "bundles_single_graph_evolution_row": single,
            "single_row_max_abs_delta": max(single_deltas) if single_deltas else 0,
        },
    )


def check_partition_bnodes(st: dict[str, Any], ru: dict[str, Any]) -> CheckResult:
    """post_bundle_partition.n_bnodes vs bundle_partition.bundle_bnodes."""
    st_pp = st["post_bundle_partition"]
    ru_rows = {(r["chrom"], int(r["start"]), int(r["end"])): r for r in ru["bundle_partition"]}
    mismatches = 0
    d_bn = []
    n = 0
    for k, srow in st_pp.items():
        rrow = ru_rows.get(k)
        if not rrow:
            continue
        n += 1
        st_bn, ru_bn = i(srow, "n_bnodes"), i(rrow, "bundle_bnodes")
        d_bn.append(abs(st_bn - ru_bn))
        if st_bn != ru_bn:
            mismatches += 1
    return CheckResult(
        name="partition_bnodes_ST_vs_RU",
        ok=mismatches == 0,
        bundles_checked=n,
        mismatches=mismatches,
        detail={
            "median_abs_delta_bnodes": median([float(x) for x in d_bn]),
            "max_abs_delta_bnodes": max(d_bn) if d_bn else 0,
        },
    )


def check_partition_subbundles_vs_colors(st: dict[str, Any], ru: dict[str, Any]) -> CheckResult:
    """post_bundle_partition.n_subbundles vs bundle_partition.bundle_color_count (weaker signal)."""
    st_pp = st["post_bundle_partition"]
    ru_rows = {(r["chrom"], int(r["start"]), int(r["end"])): r for r in ru["bundle_partition"]}
    mismatches = 0
    d_sub = []
    n = 0
    for k, srow in st_pp.items():
        rrow = ru_rows.get(k)
        if not rrow:
            continue
        n += 1
        st_sub, ru_col = i(srow, "n_subbundles"), i(rrow, "bundle_color_count")
        d_sub.append(abs(st_sub - ru_col))
        if st_sub != ru_col:
            mismatches += 1
    return CheckResult(
        name="partition_subbundles_ST_vs_RU_colors",
        ok=mismatches == 0,
        bundles_checked=n,
        mismatches=mismatches,
        detail={
            "median_abs_delta": median([float(x) for x in d_sub]),
            "max_abs_delta": max(d_sub) if d_sub else 0,
        },
    )


def check_graph_topology(
    st: dict[str, Any],
    ru: dict[str, Any],
) -> CheckResult:
    """Optional: ST post_graph_create vs one Rustle graph_evolution row (single-graph bundles only).

    Node/edge/transfrag counts can still differ because the two implementations use different
    graph accounting; this check is best-effort for regression spotting, not strict parity.
    """
    st_pg = st["post_graph_create"]
    evo_by_bundle: dict[tuple[str, int, int], list[dict[str, str]]] = defaultdict(list)
    for r in ru["graph_evolution"]:
        k = (r["chrom"], int(r["start"]), int(r["end"]))
        evo_by_bundle[k].append(r)
    mismatches = 0
    n = 0
    for k, srow in st_pg.items():
        if i(srow, "n_graphs") != 1:
            continue
        rows = evo_by_bundle.get(k)
        if not rows or len(rows) != 1:
            continue
        n += 1
        rrow = rows[0]
        sn, se, stf = i(srow, "n_nodes"), i(srow, "n_edges"), i(srow, "n_transfrags")
        rn, re, rtf = i(rrow, "n_nodes"), i(rrow, "n_edges"), i(rrow, "n_transfrags")
        if sn != rn or se != re or stf != rtf:
            mismatches += 1
    return CheckResult(
        name="graph_topology_single_component_advisory",
        ok=mismatches == 0,
        bundles_checked=n,
        mismatches=mismatches,
        detail={
            "note": "Advisory only — graph metrics are not guaranteed 1:1 across assemblers"
        },
    )


def check_parity_shadow(path: str, max_diff_rate: dict[str, float]) -> list[CheckResult]:
    """Per-layer rows: fraction with diff==1 must be <= max_diff_rate[layer]."""
    totals: dict[str, int] = defaultdict(int)
    diffs: dict[str, int] = defaultdict(int)
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            layer = row["layer"]
            totals[layer] += 1
            if row.get("diff") == "1":
                diffs[layer] += 1
    results = []
    for layer, max_rate in max_diff_rate.items():
        t = totals.get(layer, 0)
        if t == 0:
            continue
        rate = diffs[layer] / t
        ok = rate <= max_rate
        results.append(
            CheckResult(
                name=f"parity_shadow:{layer}",
                ok=ok,
                bundles_checked=t,
                mismatches=diffs[layer],
                detail={"diff_rate": rate, "max_allowed": max_rate},
            )
        )
    return results


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--stringtie", required=True, help="StringTie PARITY_STAGE_TSV")
    ap.add_argument("--rustle", required=True, help="Rustle --debug-stage-tsv")
    ap.add_argument("--parity-shadow", help="Rustle parity_shadow TSV (optional)")
    ap.add_argument(
        "--junction-tolerance",
        type=int,
        default=50,
        help="Max |ST n_junctions_kept - sum(rustle graph n_junctions)| per bundle (default 50; multi-component bundles sum several graph rows)",
    )
    ap.add_argument(
        "--read-max-delta",
        type=int,
        default=0,
        help="Max allowed |Δ reads| per bundle between StringTie bundle_input and Rustle bundle_context (default 0)",
    )
    ap.add_argument(
        "--check-graph-topology",
        action="store_true",
        help="Also compare post_graph_create vs graph_evolution for single-graph bundles (often noisy)",
    )
    ap.add_argument(
        "--strict",
        action="store_true",
        help="Exit 1 if any check fails",
    )
    ap.add_argument("--json", help="Write JSON report path")
    args = ap.parse_args()

    st = load_stringtie_by_stage(args.stringtie)
    ru = load_rustle_by_stage(args.rustle)

    checks: list[CheckResult] = []
    checks.append(check_reads(st, ru, args.read_max_delta))
    checks.append(check_junctions_kept(st, ru, args.junction_tolerance))
    checks.append(check_partition_bnodes(st, ru))
    checks.append(check_partition_subbundles_vs_colors(st, ru))
    if args.check_graph_topology:
        checks.append(check_graph_topology(st, ru))

    if args.parity_shadow:
        checks.extend(
            check_parity_shadow(
                args.parity_shadow,
                {
                    "layer2_graph_shadow": 0.01,
                    "layer1_bundle_partition": 0.55,
                    "layer3_transfrag_shadow": 0.25,
                    "layer4_transcript_shadow": 0.25,
                },
            )
        )

    report = {
        "checks": [asdict(c) for c in checks],
        "all_ok": all(c.ok for c in checks),
    }

    for c in checks:
        status = "OK" if c.ok else "FAIL"
        print(f"[{status}] {c.name} checked={c.bundles_checked} mismatches={c.mismatches} {c.detail}")

    if args.json:
        with open(args.json, "w") as f:
            json.dump(report, f, indent=2)

    if args.strict and not report["all_ok"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
