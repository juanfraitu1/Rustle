#!/usr/bin/env python3
"""Compare junction-set TSV rows (StringTie vs Rustle or two Rustle snapshots).

Expected schema (see `parity_junction_dump`):
  kind, source, stage, chrom, start, end, strand, junctions

`junctions` is `donor-acceptor` pairs sorted lexicographically, joined with `|`.

Join key: (chrom, int(start), int(end), strand). By default the same `stage` name is
required in both files; use `--ref-stage` / `--cand-stage` when comparing StringTie
(`PARITY_JUNCTION_TSV`, see `stringtie/rlink.cpp`) to Rustle (`RUSTLE_PARITY_JUNCTION_TSV`),
e.g. `aggregate_graph_ready_after_read_pass` vs `cgroup_feed`, or `read_union_post_count_good`
vs `read_union`.

Usage:
  python3 compare_junction_parity.py ref.tsv cand.tsv
  python3 compare_junction_parity.py st.tsv ru.tsv --ref-stage graph_input --cand-stage cgroup_feed
  python3 compare_junction_parity.py a.tsv b.tsv --slop 3   # ±bp on donor and acceptor
  python3 compare_junction_parity.py a.tsv b.tsv --strict   # exit 1 on any mismatch
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


def load_junction_tsv(
    path: Path,
    *,
    stage: str | None,
    source: str | None,
) -> dict[tuple[str, int, int, str], str]:
    """(chrom, start, end, strand) -> junctions string."""
    out: dict[tuple[str, int, int, str], str] = {}
    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("kind") != "junction_set_v1":
                continue
            if stage is not None and row.get("stage") != stage:
                continue
            if source is not None and row.get("source") != source:
                continue
            k = (
                row["chrom"],
                int(row["start"]),
                int(row["end"]),
                row["strand"],
            )
            out[k] = row.get("junctions") or ""
    return out


def parse_junction_pairs(junctions: str) -> list[tuple[int, int]]:
    s = junctions.strip()
    if not s:
        return []
    out: list[tuple[int, int]] = []
    for part in s.split("|"):
        part = part.strip()
        if not part:
            continue
        d, _, a = part.partition("-")
        out.append((int(d), int(a)))
    return out


def max_bipartite_matching(
    left_size: int, right_size: int, adj: list[list[int]]
) -> int:
    """Maximum cardinality matching, left i to right j."""
    match_r = [-1] * right_size

    def aug(u: int, seen: list[bool]) -> bool:
        for v in adj[u]:
            if v < 0 or v >= right_size:
                continue
            if seen[v]:
                continue
            seen[v] = True
            if match_r[v] < 0 or aug(match_r[v], seen):
                match_r[v] = u
                return True
        return False

    n = 0
    for u in range(left_size):
        if aug(u, [False] * right_size):
            n += 1
    return n


def slop_equivalent(a: str, b: str, slop: int) -> bool:
    if slop <= 0:
        return a == b
    pa = parse_junction_pairs(a)
    pb = parse_junction_pairs(b)
    if len(pa) != len(pb):
        return False
    if not pa:
        return True
    n, m = len(pa), len(pb)
    adj: list[list[int]] = [[] for _ in range(n)]
    for i, (da, aa) in enumerate(pa):
        for j, (db, ab) in enumerate(pb):
            if abs(da - db) <= slop and abs(aa - ab) <= slop:
                adj[i].append(j)
    return max_bipartite_matching(n, m, adj) == n == m


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("reference_tsv", type=Path, help="e.g. StringTie junction_set TSV")
    ap.add_argument("candidate_tsv", type=Path, help="e.g. Rustle RUSTLE_PARITY_JUNCTION_TSV")
    ap.add_argument(
        "--stage",
        type=str,
        default=None,
        help="require this stage in both files (ignored if --ref-stage/--cand-stage set)",
    )
    ap.add_argument(
        "--ref-stage",
        type=str,
        default=None,
        help="stage column value for reference file",
    )
    ap.add_argument(
        "--cand-stage",
        type=str,
        default=None,
        help="stage column value for candidate file",
    )
    ap.add_argument("--ref-source", type=str, default=None, help="filter reference rows by source")
    ap.add_argument("--cand-source", type=str, default=None, help="filter candidate rows by source")
    ap.add_argument(
        "--slop",
        type=int,
        default=0,
        help="treat junction sets as equal if same count and one-to-one ±slop bp on donor and acceptor",
    )
    ap.add_argument("--strict", action="store_true", help="exit 1 if any joined key mismatches")
    ap.add_argument(
        "--show",
        type=int,
        default=12,
        help="max mismatch examples to print",
    )
    args = ap.parse_args()

    if args.ref_stage is None and args.cand_stage is None:
        stage_ref = stage_cand = args.stage or "graph_input"
    else:
        stage_ref = args.ref_stage or args.stage or "graph_input"
        stage_cand = args.cand_stage or args.stage or "graph_input"

    ref = load_junction_tsv(
        args.reference_tsv, stage=stage_ref, source=args.ref_source
    )
    cand = load_junction_tsv(
        args.candidate_tsv, stage=stage_cand, source=args.cand_source
    )
    keys = sorted(set(ref) & set(cand))

    def same(ja: str, jb: str) -> bool:
        if args.slop > 0:
            return slop_equivalent(ja, jb, args.slop)
        return ja == jb

    matched = sum(1 for k in keys if same(ref[k], cand[k]))
    missing_ref = sorted(set(cand) - set(ref))
    missing_cand = sorted(set(ref) - set(cand))
    mism = [k for k in keys if not same(ref[k], cand[k])]

    mode = f"slop={args.slop}" if args.slop else "exact string"
    print(
        f"Reference stage={stage_ref!r} ({args.reference_tsv}); "
        f"candidate stage={stage_cand!r} ({args.candidate_tsv})"
    )
    print(f"Compare mode: {mode}")
    print(f"Joined bundles: {len(keys)}")
    print(
        f"Junction-set match: {matched} / {len(keys)} "
        f"({100.0 * matched / len(keys) if keys else 0:.2f}%)"
    )
    print(f"Only in reference file: {len(missing_cand)}")
    print(f"Only in candidate file: {len(missing_ref)}")
    if mism and args.show:
        print(f"First mismatches (up to {args.show}):")
        for k in mism[: args.show]:
            print(f"  {k}")
            ra, ca = ref[k], cand[k]
            print(f"    REF ({len(parse_junction_pairs(ra))} pairs): {ra[:160]}{'...' if len(ra) > 160 else ''}")
            print(f"    CAND ({len(parse_junction_pairs(ca))} pairs): {ca[:160]}{'...' if len(ca) > 160 else ''}")

    if args.strict and (mism or missing_ref or missing_cand):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
