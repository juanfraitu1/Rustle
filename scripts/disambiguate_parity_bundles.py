#!/usr/bin/env python3
"""Join partition mismatches to junction sets and Rustle trace rows for layer attribution.

For each bundle where StringTie vs Rustle partition signatures differ, adds columns that
separate:

* **Junction table** — same bundle key in `PARITY_JUNCTION_TSV` / `RUSTLE_PARITY_JUNCTION_TSV`
  at chosen stages (exact string, Jaccard, optional ±slop bipartite match).
* **Rustle graph evolution** — `graph_topology_v1` rows: `post_graph_build` vs `pre_read_map`
  fingerprints differ when coverage / terminal synthesis / pruning changed adjacency before
  read mapping.
* **Rustle read mapping** — `read_map_summary_v1` (`read_mapped_tfs` per `run_tag`).
* **Rustle cgroup** — `cgroup_summary_v1` (sub-bundle and bundlenode counts).

StringTie does not emit graph/transfrag trace rows yet; those columns are Rustle-only anchors.

Example:

  python3 scripts/disambiguate_parity_bundles.py \\
    --partition-st st_part.tsv --partition-ru ru_part.tsv \\
    --junction-st st_junc.tsv --junction-ru ru_junc.tsv \\
    --rustle-trace ru_trace.tsv \\
    --st-junc-stage aggregate_graph_ready_after_read_pass \\
    --ru-junc-stage cgroup_feed \\
    -o disambig.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


# --- partition (aligned with parity_fidelity_report / compare_partition_geometry) ---


def load_partition_rows(path: Path) -> dict[tuple[str, int, int], dict[str, str]]:
    out: dict[tuple[str, int, int], dict[str, str]] = {}
    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("kind") != "partition_geometry_v1":
                continue
            k = (row["chrom"], int(row["start"]), int(row["end"]))
            out[k] = dict(row)
    return out


def parse_partition_chains(sig: str) -> list[list[tuple[int, int]]]:
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


def _parse_flat_segments(signature: str) -> list[tuple[int, int]]:
    if not signature or not signature.strip():
        return []
    segs: list[tuple[int, int]] = []
    for chain in signature.split("|"):
        chain = chain.strip()
        if not chain:
            continue
        for part in chain.split("+"):
            part = part.strip()
            if not part or "-" not in part:
                continue
            a, b = part.split("-", 1)
            segs.append((int(a), int(b)))
    return segs


def segment_multiset_fuzzy_match(sig_a: str, sig_b: str, slop: int) -> bool:
    if slop <= 0:
        return sig_a == sig_b
    sa = sorted(_parse_flat_segments(sig_a))
    sb = sorted(_parse_flat_segments(sig_b))
    if len(sa) != len(sb):
        return False
    for (x0, x1), (y0, y1) in zip(sa, sb):
        if abs(x0 - y0) > slop or abs(x1 - y1) > slop:
            return False
    return True


# --- junction ---


def load_junction_stage(path: Path, stage: str) -> dict[tuple[str, int, int, str], str]:
    out: dict[tuple[str, int, int, str], str] = {}
    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("kind") != "junction_set_v1":
                continue
            if row.get("stage") != stage:
                continue
            k = (
                row["chrom"],
                int(row["start"]),
                int(row["end"]),
                row.get("strand") or ".",
            )
            out[k] = row.get("junctions") or ""
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


def junction_slop_equivalent(a: str, b: str, slop: int) -> bool:
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


# --- Rustle trace ---


def parse_payload(payload: str) -> dict[str, str]:
    d: dict[str, str] = {}
    for part in (payload or "").split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        k, _, v = part.partition("=")
        d[k.strip()] = v.strip()
    return d


def load_trace(path: Path | None) -> dict[str, dict[tuple[str, int, int], list[dict[str, str]]]]:
    """kind -> bundle_key -> list of row dicts (chrom/start/end/strand/stage/payload_parsed)."""
    if path is None or not path.is_file():
        return {}
    by_kind: dict[str, dict[tuple[str, int, int], list[dict[str, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            kind = row.get("kind") or ""
            if kind not in (
                "cgroup_summary_v1",
                "graph_topology_v1",
                "read_map_summary_v1",
            ):
                continue
            k = (row["chrom"], int(row["start"]), int(row["end"]))
            payload = parse_payload(row.get("payload") or "")
            rec = {
                "strand": row.get("strand") or ".",
                "stage": row.get("stage") or "",
                **payload,
            }
            by_kind[kind][k].append(rec)
    return dict(by_kind)


def load_trace_events(path: Path | None) -> dict[tuple[str, int, int], list[dict[str, str]]]:
    if path is None or not path.is_file():
        return {}
    out: dict[tuple[str, int, int], list[dict[str, str]]] = defaultdict(list)
    with path.open(newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row.get("kind") != "trace_event_v1":
                continue
            k = (row["chrom"], int(row["start"]), int(row["end"]))
            payload = parse_payload(row.get("payload") or "")
            rec = {"stage": row.get("stage") or "", **payload}
            out[k].append(rec)
    return dict(out)


def summarize_trace_events(rows: list[dict[str, str]]) -> str:
    if not rows:
        return ""
    counts: dict[str, int] = defaultdict(int)
    for rec in rows:
        stage = rec.get("stage", "")
        event = rec.get("event", "")
        key = f"{stage}:{event}" if event else stage
        counts[key] += 1
    parts = [f"{k}={counts[k]}" for k in sorted(counts)]
    return ",".join(parts)


def _fps_signature(rows: list[dict[str, str]], stage: str) -> str:
    """Sorted 'run_tag:edge_fp' for stable multiset comparison."""
    xs = []
    for rec in rows:
        if rec.get("stage") != stage:
            continue
        rt = rec.get("run_tag", "")
        fp = rec.get("edge_fp", "")
        xs.append(f"{rt}:{fp}")
    xs.sort()
    return "|".join(xs)


def _aggregate_read_map_tfs(rows: list[dict[str, str]]) -> tuple[int, str]:
    total = 0
    parts = []
    for rec in rows:
        if rec.get("stage") != "post_read_map":
            continue
        try:
            n = int(rec.get("read_mapped_tfs", "0"))
        except ValueError:
            n = 0
        total += n
        parts.append(f"{rec.get('run_tag','')}={n}")
    parts.sort()
    return total, ",".join(parts)


def _first_cgroup(rows: list[dict[str, str]]) -> dict[str, str]:
    for rec in rows:
        if rec.get("stage") == "post_cgroup":
            return {k: rec.get(k, "") for k in ("sub_bundles", "bnodes", "color_roots", "assigned_reads")}
    return {}


def attribution_note(
    *,
    part_class: str,
    seg_fuzzy: bool,
    junc_exact: bool,
    junc_slop_ok: bool,
    graph_sig_same: bool | None,
) -> str:
    """Single-line hint for humans (not ground truth)."""
    bits: list[str] = []
    if not junc_exact:
        if junc_slop_ok:
            bits.append("junction_strings_differ_slop_bipartite_match_same_cardinality")
        else:
            bits.append("junction_intron_sets_differ_beyond_slop_at_paired_stages")
    elif graph_sig_same is False:
        bits.append("RU_graph_edge_fp_changed_post_graph_build_to_pre_read_map")
    elif seg_fuzzy and junc_exact:
        bits.append("partition_multiset_fuzzy_ok_junction_exact_topology_or_chain_order")
    elif part_class.startswith("chain_count_diff") and junc_exact:
        bits.append("subbundle_CBundle_chain_count_diff_junction_table_matches")
    elif part_class == "same_counts_topology_diff" and junc_exact:
        bits.append("same_segment_counts_topology_only_diff_junction_table_matches")
    if not bits:
        bits.append("inspect_partition_class_ru_trace_and_junction_columns")
    return ";".join(bits)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--partition-st", type=Path, required=True)
    ap.add_argument("--partition-ru", type=Path, required=True)
    ap.add_argument("--junction-st", type=Path, default=None)
    ap.add_argument("--junction-ru", type=Path, default=None)
    ap.add_argument("--rustle-trace", type=Path, default=None)
    ap.add_argument("--stringtie-trace", type=Path, default=None)
    ap.add_argument(
        "--st-junc-stage",
        default="aggregate_graph_ready_after_read_pass",
        help="StringTie PARITY_JUNCTION_TSV stage column",
    )
    ap.add_argument(
        "--ru-junc-stage",
        default="cgroup_feed",
        help="Rustle junction dump stage column",
    )
    ap.add_argument("--junction-slop", type=int, default=2, help="slop bp for junction_slop_ok column")
    ap.add_argument("--segment-slop", type=int, default=3, help="partition multiset segment slop")
    ap.add_argument(
        "--all-bundles",
        action="store_true",
        help="emit every joined bundle key (not only partition mismatches)",
    )
    ap.add_argument("-o", "--output", type=Path, required=True)
    args = ap.parse_args()

    st_part = load_partition_rows(args.partition_st)
    ru_part = load_partition_rows(args.partition_ru)
    joined = sorted(set(st_part) & set(ru_part))

    if (args.junction_st is None) ^ (args.junction_ru is None):
        ap.error("use both --junction-st and --junction-ru, or omit both")

    st_j: dict[tuple[str, int, int, str], str] = {}
    ru_j: dict[tuple[str, int, int, str], str] = {}
    if args.junction_st is not None:
        st_j = load_junction_stage(args.junction_st, args.st_junc_stage)
        ru_j = load_junction_stage(args.junction_ru, args.ru_junc_stage)

    trace = load_trace(args.rustle_trace)
    cgroup_by = trace.get("cgroup_summary_v1", {})
    graph_by = trace.get("graph_topology_v1", {})
    readmap_by = trace.get("read_map_summary_v1", {})
    ru_events_by = load_trace_events(args.rustle_trace)
    st_events_by = load_trace_events(args.stringtie_trace)

    fieldnames = [
        "chrom",
        "start",
        "end",
        "strand",
        "partition_mismatch_class",
        "partition_exact",
        "partition_segment_fuzzy_ok",
        "st_junc_count",
        "ru_junc_count",
        "junction_exact_string",
        "junction_jaccard",
        f"junction_slop{args.junction_slop}_bipartite_ok",
        "ru_cgroup_sub_bundles",
        "ru_cgroup_bnodes",
        "ru_cgroup_color_roots",
        "ru_cgroup_assigned_reads",
        "ru_graph_rows_post_build",
        "ru_graph_rows_pre_read_map",
        "ru_graph_fp_sig_post_build",
        "ru_graph_fp_sig_pre_read_map",
        "ru_graph_topology_changed_pre_read_map",
        "ru_read_mapped_tfs_sum",
        "ru_read_mapped_tfs_by_run_tag",
        "st_trace_events",
        "ru_trace_events",
        "attribution_hint",
    ]

    slop_j = max(0, int(args.junction_slop))
    seg_slop = max(0, int(args.segment_slop))

    with args.output.open("w", newline="") as outf:
        w = csv.DictWriter(outf, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        for k in joined:
            st_sig = st_part[k].get("signature") or ""
            ru_sig = ru_part[k].get("signature") or ""
            exact = st_sig == ru_sig
            if not args.all_bundles and exact:
                continue

            strand = ru_part[k].get("strand") or st_part[k].get("strand") or "."
            jk = (*k, strand)
            st_js = st_j.get(jk, "")
            ru_js = ru_j.get(jk, "")
            st_set = pair_set(st_js)
            ru_set = pair_set(ru_js)
            junc_exact = st_js == ru_js
            j_slop_ok = junction_slop_equivalent(st_js, ru_js, slop_j) if slop_j > 0 else junc_exact
            seg_fuzzy = segment_multiset_fuzzy_match(st_sig, ru_sig, seg_slop)

            g_rows = graph_by.get(k, [])
            sig_post = _fps_signature(g_rows, "post_graph_build")
            sig_pre = _fps_signature(g_rows, "pre_read_map")
            n_post = sum(1 for r in g_rows if r.get("stage") == "post_graph_build")
            n_pre = sum(1 for r in g_rows if r.get("stage") == "pre_read_map")
            topo_changed: bool | None
            if not g_rows:
                topo_changed = None
            else:
                topo_changed = sig_post != sig_pre

            cg = _first_cgroup(cgroup_by.get(k, []))
            rmt_sum, rmt_detail = _aggregate_read_map_tfs(readmap_by.get(k, []))
            st_evt = summarize_trace_events(st_events_by.get(k, []))
            ru_evt = summarize_trace_events(ru_events_by.get(k, []))

            pclass = partition_mismatch_class(st_sig, ru_sig)
            hint = attribution_note(
                part_class=pclass,
                seg_fuzzy=seg_fuzzy,
                junc_exact=junc_exact,
                junc_slop_ok=j_slop_ok,
                graph_sig_same=(not topo_changed) if topo_changed is not None else None,
            )

            w.writerow(
                {
                    "chrom": k[0],
                    "start": k[1],
                    "end": k[2],
                    "strand": strand,
                    "partition_mismatch_class": pclass,
                    "partition_exact": int(exact),
                    "partition_segment_fuzzy_ok": int(seg_fuzzy),
                    "st_junc_count": len(st_set),
                    "ru_junc_count": len(ru_set),
                    "junction_exact_string": int(junc_exact),
                    "junction_jaccard": f"{jaccard_sets(st_set, ru_set):.6f}",
                    f"junction_slop{args.junction_slop}_bipartite_ok": int(j_slop_ok),
                    "ru_cgroup_sub_bundles": cg.get("sub_bundles", ""),
                    "ru_cgroup_bnodes": cg.get("bnodes", ""),
                    "ru_cgroup_color_roots": cg.get("color_roots", ""),
                    "ru_cgroup_assigned_reads": cg.get("assigned_reads", ""),
                    "ru_graph_rows_post_build": n_post,
                    "ru_graph_rows_pre_read_map": n_pre,
                    "ru_graph_fp_sig_post_build": sig_post,
                    "ru_graph_fp_sig_pre_read_map": sig_pre,
                    "ru_graph_topology_changed_pre_read_map": (
                        "" if topo_changed is None else int(topo_changed)
                    ),
                    "ru_read_mapped_tfs_sum": rmt_sum,
                    "ru_read_mapped_tfs_by_run_tag": rmt_detail,
                    "st_trace_events": st_evt,
                    "ru_trace_events": ru_evt,
                    "attribution_hint": hint,
                }
            )

    print(f"Wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
