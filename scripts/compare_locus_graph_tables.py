#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def parse_chain(text: str):
    if not text:
        return []
    out = []
    for item in text.split(","):
        if not item:
            continue
        a, b = item.split("-")
        out.append((int(a), int(b)))
    return out


def parse_nodes(text: str):
    if not text:
        return []
    return [int(x) for x in text.split(",") if x]


def matching_junctions(ref_chain, query_chain, tol):
    n = min(len(ref_chain), len(query_chain))
    matched = 0
    for i in range(n):
        if (
            abs(ref_chain[i][0] - query_chain[i][0]) <= tol
            and abs(ref_chain[i][1] - query_chain[i][1]) <= tol
        ):
            matched += 1
        else:
            break
    return matched


def classify_match(ref_nodes, ref_chain, query_nodes, query_chain, tol):
    if ref_nodes == query_nodes:
        return "exact_path"
    if ref_chain == query_chain:
        return "exact_chain"
    if len(ref_chain) == len(query_chain) and all(
        abs(rd - qd) <= tol and abs(ra - qa) <= tol
        for (rd, ra), (qd, qa) in zip(ref_chain, query_chain)
    ):
        return "fuzzy_chain"
    mj = matching_junctions(ref_chain, query_chain, tol)
    if mj > 0:
        return f"prefix_{mj}"
    return "none"


def first_divergence(ref_chain, query_chain, tol):
    n = min(len(ref_chain), len(query_chain))
    for i in range(n):
        rd, ra = ref_chain[i]
        qd, qa = query_chain[i]
        if abs(rd - qd) > tol or abs(ra - qa) > tol:
            return f"junction_{i + 1}:{rd}-{ra}!={qd}-{qa}"
    if len(ref_chain) != len(query_chain):
        return f"chain_len:{len(ref_chain)}!={len(query_chain)}"
    return "boundary_only"


def source_rank(source: str):
    order = {
        "exact_full": 0,
        "stitched_full": 1,
        "graph_completed": 2,
        "total_support": 3,
    }
    return order.get(source, 9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--transcript-prefix", required=True)
    ap.add_argument("--denovo-prefix", required=True)
    ap.add_argument("--output-prefix", required=True)
    ap.add_argument("--junction-tolerance", type=int, default=0)
    args = ap.parse_args()

    transcript_paths = []
    with open(f"{args.transcript_prefix}.paths.tsv") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            transcript_paths.append(
                {
                    "transcript_id": row["transcript_id"],
                    "nodes": parse_nodes(row["nodes"]),
                    "chain": parse_chain(row["junction_chain"]),
                }
            )

    denovo_paths = []
    with open(f"{args.denovo_prefix}.paths.tsv") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            denovo_paths.append(
                {
                    "source": row["source"],
                    "count": int(row["count"]),
                    "weight": float(row["weight"]),
                    "nodes": parse_nodes(row["nodes"]),
                    "chain": parse_chain(row["junction_chain"]),
                }
            )

    ref_out = Path(f"{args.output_prefix}.refs.tsv")
    sum_out = Path(f"{args.output_prefix}.summary.tsv")

    rows = []
    for ref in transcript_paths:
        best = None
        best_key = None
        for query in denovo_paths:
            match_kind = classify_match(
                ref["nodes"],
                ref["chain"],
                query["nodes"],
                query["chain"],
                args.junction_tolerance,
            )
            mj = matching_junctions(ref["chain"], query["chain"], args.junction_tolerance)
            terminal_delta = abs(ref["nodes"][0] - query["nodes"][0]) + abs(ref["nodes"][-1] - query["nodes"][-1])
            rank = (
                1 if match_kind == "exact_path" else 0,
                1 if match_kind == "exact_chain" else 0,
                1 if match_kind == "fuzzy_chain" else 0,
                mj,
                -terminal_delta,
                -query["count"],
                -query["weight"],
                -len(query["nodes"]),
                -source_rank(query["source"]),
            )
            if best is None or rank > best_key:
                best = query
                best_key = rank

        match_kind = classify_match(
            ref["nodes"],
            ref["chain"],
            best["nodes"],
            best["chain"],
            args.junction_tolerance,
        )
        rows.append(
            {
                "transcript_id": ref["transcript_id"],
                "ref_nodes": ",".join(map(str, ref["nodes"])),
                "ref_junction_chain": ",".join(f"{a}-{b}" for a, b in ref["chain"]),
                "best_source": best["source"],
                "best_count": best["count"],
                "best_weight": f"{best['weight']:.3f}",
                "best_nodes": ",".join(map(str, best["nodes"])),
                "best_junction_chain": ",".join(f"{a}-{b}" for a, b in best["chain"]),
                "match_kind": match_kind,
                "matching_junctions": matching_junctions(ref["chain"], best["chain"], args.junction_tolerance),
                "ref_junction_count": len(ref["chain"]),
                "first_divergence": first_divergence(ref["chain"], best["chain"], args.junction_tolerance)
                if match_kind not in ("exact_path", "exact_chain")
                else ("boundary_only" if match_kind == "exact_chain" else "none"),
                "ref_start": ref["nodes"][0],
                "ref_end": ref["nodes"][-1],
                "best_start": best["nodes"][0],
                "best_end": best["nodes"][-1],
            }
        )

    with ref_out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "transcript_id",
                "match_kind",
                "matching_junctions",
                "ref_junction_count",
                "first_divergence",
                "best_source",
                "best_count",
                "best_weight",
                "ref_start",
                "ref_end",
                "best_start",
                "best_end",
                "ref_nodes",
                "best_nodes",
                "ref_junction_chain",
                "best_junction_chain",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    summary_counts = {}
    for row in rows:
        summary_counts[row["match_kind"]] = summary_counts.get(row["match_kind"], 0) + 1

    with sum_out.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["reference_transcripts", len(rows)])
        for key in sorted(summary_counts):
            writer.writerow([f"match_kind:{key}", summary_counts[key]])


if __name__ == "__main__":
    main()
