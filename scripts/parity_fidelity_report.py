#!/usr/bin/env python3
"""Consolidated StringTie vs Rustle parity dashboard (partition + junction + rubric).

Reads paired instrumentation TSVs from the same BAM run and reports:
  * Partition geometry (`PARITY_PARTITION_TSV` vs `RUSTLE_PARITY_PARTITION_TSV`)
  * Optional junction sets (`PARITY_JUNCTION_TSV` vs `RUSTLE_PARITY_JUNCTION_TSV`) for
    several stage alignments (exact match and Jaccard on splice-pair sets)
  * A short rubric of what is / is not covered toward “full” transcript-level parity

Nothing is executed except this script — point it at your existing TSVs.

Usage:
  python3 parity_fidelity_report.py \\
    --partition-st st_part.tsv --partition-ru ru_part.tsv \\
    --junction-st st_junc.tsv --junction-ru ru_junc.tsv

  python3 parity_fidelity_report.py --partition-st st.tsv --partition-ru ru.tsv --json out.json
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


# --- partition (same logic as compare_partition_geometry.py) -----------------


def load_partition_rows(path: Path) -> dict[tuple[str, int, int], str]:
    out: dict[tuple[str, int, int], str] = {}
    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("kind") != "partition_geometry_v1":
                continue
            k = (row["chrom"], int(row["start"]), int(row["end"]))
            out[k] = row.get("signature") or ""
    return out


def parse_partition_chains(sig: str) -> list[list[tuple[int, int]]]:
    """Signature: chains joined by '|', segments by '+', each segment 'start-end'."""
    sig = (sig or "").strip()
    if not sig:
        return []
    chains: list[list[tuple[int, int]]] = []
    for ch in sig.split("|"):
        ch = ch.strip()
        if not ch:
            continue
        segs: list[tuple[int, int]] = []
        for seg in ch.split("+"):
            seg = seg.strip()
            if not seg or "-" not in seg:
                continue
            a, _, b = seg.partition("-")
            segs.append((int(a), int(b)))
        if segs:
            chains.append(segs)
    return chains


def partition_mismatch_class(st_sig: str, ru_sig: str) -> str:
    stc = parse_partition_chains(st_sig)
    ruc = parse_partition_chains(ru_sig)
    if not stc and ruc:
        return "st_empty_ru_nonempty"
    if stc and not ruc:
        return "st_nonempty_ru_empty"
    if not stc and not ruc:
        return "both_empty"
    if len(stc) != len(ruc):
        return f"chain_count_diff(st={len(stc)},ru={len(ruc)})"
    st_segs = sum(len(c) for c in stc)
    ru_segs = sum(len(c) for c in ruc)
    if st_segs != ru_segs:
        return f"segment_count_diff(st={st_segs},ru={ru_segs})"
    return "same_counts_topology_diff"


@dataclass
class PartitionReport:
    joined_bundles: int
    union_bundles: int
    exact_signature_matches: int
    exact_match_pct: float
    only_in_stringtie: int
    only_in_rustle: int
    mismatch_histogram: dict[str, int]
    ru_empty_signature_on_mismatch: int


def analyze_partition(st: Path, ru: Path) -> PartitionReport:
    d_st = load_partition_rows(st)
    d_ru = load_partition_rows(ru)
    keys_st = set(d_st)
    keys_ru = set(d_ru)
    joined = keys_st & keys_ru
    union = keys_st | keys_ru
    mism_hist: dict[str, int] = {}
    ru_empty_mism = 0
    exact = 0
    for k in joined:
        if d_st[k] == d_ru[k]:
            exact += 1
        else:
            cls = partition_mismatch_class(d_st[k], d_ru[k])
            mism_hist[cls] = mism_hist.get(cls, 0) + 1
            if not (d_ru[k] or "").strip():
                ru_empty_mism += 1
    n_joined = len(joined)
    return PartitionReport(
        joined_bundles=n_joined,
        union_bundles=len(union),
        exact_signature_matches=exact,
        exact_match_pct=(100.0 * exact / n_joined) if n_joined else 0.0,
        only_in_stringtie=len(keys_st - keys_ru),
        only_in_rustle=len(keys_ru - keys_st),
        mismatch_histogram=mism_hist,
        ru_empty_signature_on_mismatch=ru_empty_mism,
    )


# --- junction (extends compare_junction_parity metrics) ----------------------


def load_junction_rows(
    path: Path, *, stage: str | None, source: str | None = None
) -> dict[tuple[str, int, int, str], str]:
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
            key = (
                row["chrom"],
                int(row["start"]),
                int(row["end"]),
                row.get("strand") or ".",
            )
            out[key] = row.get("junctions") or ""
    return out


def parse_junction_pairs(junctions: str) -> list[tuple[int, int]]:
    s = (junctions or "").strip()
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


def pair_set(junctions: str) -> set[tuple[int, int]]:
    return set(parse_junction_pairs(junctions))


def jaccard_sets(a: set[tuple[int, int]], b: set[tuple[int, int]]) -> float:
    if not a and not b:
        return 1.0
    u = len(a | b)
    if u == 0:
        return 1.0
    return len(a & b) / u


def max_bipartite_matching(left_size: int, right_size: int, adj: list[list[int]]) -> int:
    match_r = [-1] * right_size

    def aug(u: int, seen: list[bool]) -> bool:
        for v in adj[u]:
            if v < 0 or v >= right_size or seen[v]:
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


def _slop_equivalent(a: str, b: str, slop: int) -> bool:
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


@dataclass
class JunctionStageReport:
    ref_stage: str
    cand_stage: str
    slop: int
    joined_bundles: int
    exact_junction_string_matches: int
    exact_match_pct: float
    mean_jaccard_on_joined: float
    only_in_stringtie: int
    only_in_rustle: int


def analyze_junction_stage_pair_noslop_import(
    st_path: Path,
    ru_path: Path,
    *,
    ref_stage: str,
    cand_stage: str,
    slop: int,
) -> JunctionStageReport:
    """Avoid importing compare_junction_parity (package path issues when run as script)."""
    ref = load_junction_rows(st_path, stage=ref_stage)
    cand = load_junction_rows(ru_path, stage=cand_stage)
    keys = set(ref) & set(cand)
    exact = 0
    jaccs: list[float] = []
    for k in keys:
        ra, ca = ref[k], cand[k]
        if slop <= 0:
            match = ra == ca
        else:
            match = _slop_equivalent(ra, ca, slop)
        if match:
            exact += 1
        jaccs.append(jaccard_sets(pair_set(ra), pair_set(ca)))
    n = len(keys)
    mean_j = sum(jaccs) / len(jaccs) if jaccs else 0.0
    return JunctionStageReport(
        ref_stage=ref_stage,
        cand_stage=cand_stage,
        slop=slop,
        joined_bundles=n,
        exact_junction_string_matches=exact,
        exact_match_pct=(100.0 * exact / n) if n else 0.0,
        mean_jaccard_on_joined=mean_j,
        only_in_stringtie=len(set(ref) - set(cand)),
        only_in_rustle=len(set(cand) - set(ref)),
    )


DEFAULT_JUNCTION_STAGE_PAIRS: list[tuple[str, str]] = [
    ("aggregate_graph_ready_after_read_pass", "cgroup_feed"),
    ("aggregate_graph_ready_post_bundle_partition", "cgroup_feed"),
    ("read_union_post_count_good", "read_union"),
    ("aggregate_stranded_post_count_good", "graph_input"),
    ("read_union_pre_count_good", "read_union"),
]


def print_rubric() -> None:
    print(
        "\n--- Rubric: what “full parity” means ---\n"
        "  1) Partition geometry exact match  → same bundlenode chains before splice graphs.\n"
        "  2) Junction-set match (right stage) → same intron graph fed to cgroup / graph build.\n"
        "  3) Not covered here: read exon arrays, cgroup/merge_close_groups, 3-strand filters,\n"
        "     graph edges, path extraction, transfrag assembly, GTF line-by-line identity.\n"
        "  Full transcript parity requires all of the above; this report bounds (1) and (2).\n"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--partition-st", type=Path, default=None)
    ap.add_argument("--partition-ru", type=Path, default=None)
    ap.add_argument("--junction-st", type=Path, default=None)
    ap.add_argument("--junction-ru", type=Path, default=None)
    ap.add_argument(
        "--junction-slop",
        type=int,
        nargs="+",
        default=None,
        help="slop values for junction comparisons (default: 0 2 if --junction-* given)",
    )
    ap.add_argument("--json", type=Path, default=None, help="write machine-readable summary")
    args = ap.parse_args()

    if not (
        (args.partition_st and args.partition_ru)
        or (args.junction_st and args.junction_ru)
    ):
        print(
            "Provide at least one pair: (--partition-st + --partition-ru) "
            "or (--junction-st + --junction-ru).",
            file=sys.stderr,
        )
        return 2

    report: dict = {
        "partition": None,
        "junction_stage_rows": [],
        "summary": {},
        "rubric": "see stdout",
    }
    pr: PartitionReport | None = None

    if args.partition_st and args.partition_ru:
        pr = analyze_partition(args.partition_st, args.partition_ru)
        report["partition"] = asdict(pr)
        print("=== Partition geometry (partition_geometry_v1) ===")
        print(f"  StringTie: {args.partition_st}")
        print(f"  Rustle:    {args.partition_ru}")
        print(f"  Joined bundle keys:     {pr.joined_bundles}")
        print(f"  Union bundle keys:      {pr.union_bundles}")
        print(
            f"  Exact signature match:  {pr.exact_signature_matches} / {pr.joined_bundles} "
            f"({pr.exact_match_pct:.2f}%)"
        )
        print(f"  Keys only in StringTie: {pr.only_in_stringtie}")
        print(f"  Keys only in Rustle:    {pr.only_in_rustle}")
        if pr.joined_bundles:
            cov = 100.0 * pr.joined_bundles / pr.union_bundles if pr.union_bundles else 0.0
            print(f"  Key coverage (joined/union): {cov:.2f}%")
        if pr.mismatch_histogram:
            print("  Mismatch taxonomy (joined keys where signatures differ):")
            for k, v in sorted(pr.mismatch_histogram.items(), key=lambda x: -x[1]):
                print(f"    {k}: {v}")
            if pr.ru_empty_signature_on_mismatch:
                print(
                    f"  (Rustle empty signature on those mismatches: {pr.ru_empty_signature_on_mismatch})"
                )
        print()

    junction_slops = args.junction_slop if args.junction_slop is not None else [0, 2]
    best_junction_sl0: tuple[float, str, str, float] | None = None  # exact%, st_s, ru_s, jacc%

    if args.junction_st and args.junction_ru:
        print("=== Junction sets (junction_set_v1) ===")
        print(f"  StringTie: {args.junction_st}")
        print(f"  Rustle:    {args.junction_ru}")
        for sl in junction_slops:
            print(f"\n  --- slop={sl} (exact% uses slop-tolerant match when sl>0) ---")
            for ref_s, cand_s in DEFAULT_JUNCTION_STAGE_PAIRS:
                jr = analyze_junction_stage_pair_noslop_import(
                    args.junction_st,
                    args.junction_ru,
                    ref_stage=ref_s,
                    cand_stage=cand_s,
                    slop=sl,
                )
                if jr.joined_bundles == 0:
                    continue
                report["junction_stage_rows"].append(asdict(jr))
                print(
                    f"    ST {ref_s!r:50} → RU {cand_s!r:22}  "
                    f"exact {jr.exact_match_pct:6.2f}%  "
                    f"mean Jaccard {jr.mean_jaccard_on_joined*100:5.1f}%  "
                    f"joined={jr.joined_bundles}  onlyST={jr.only_in_stringtie} onlyRU={jr.only_in_rustle}"
                )
                if sl == 0 and (
                    best_junction_sl0 is None or jr.exact_match_pct > best_junction_sl0[0]
                ):
                    best_junction_sl0 = (
                        jr.exact_match_pct,
                        ref_s,
                        cand_s,
                        jr.mean_jaccard_on_joined * 100.0,
                    )
        print()

    # Bounded “how close” summary (does not claim GTF identity).
    summ: dict = {}
    if pr is not None:
        summ["partition_exact_pct_on_joined"] = round(pr.exact_match_pct, 3)
        summ["partition_joined_bundles"] = pr.joined_bundles
    if best_junction_sl0 is not None:
        summ["best_junction_exact_pct_slop0"] = round(best_junction_sl0[0], 3)
        summ["best_junction_stage_pair_slop0"] = [best_junction_sl0[1], best_junction_sl0[2]]
        summ["mean_jaccard_pct_at_best_pair_slop0"] = round(best_junction_sl0[3], 3)
    if pr is not None and best_junction_sl0 is not None:
        summ["bounded_blend_pct"] = round(
            (pr.exact_match_pct + best_junction_sl0[0]) / 2.0, 3
        )
        summ["bounded_blend_note"] = (
            "average of partition exact% and best default stage-pair junction exact% (slop=0); "
            "not a guarantee of GTF or path parity"
        )
    report["summary"] = summ
    if summ:
        print("=== Summary (instrumentation bounds; not full GTF parity) ===")
        for k, v in summ.items():
            print(f"  {k}: {v}")
        print(
            "  For transcript-level closeness, also run gffcompare (or your preferred GTF diff) "
            "on paired GTFs; see scripts/compare_gtf_attrs_exact_matches.py if applicable.\n"
        )

    print_rubric()

    if args.json:
        with args.json.open("w") as jf:
            json.dump(report, jf, indent=2)
        print(f"Wrote JSON summary to {args.json}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
