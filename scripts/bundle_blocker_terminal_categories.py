#!/usr/bin/env python3

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Classify exact backward-supported bundle blockers by terminal super-node behavior."
        )
    )
    parser.add_argument("--bnodes", required=True)
    parser.add_argument("--exact-neighbors-summary", required=True)
    parser.add_argument("--transfrags", required=True)
    parser.add_argument("--output-stem", required=True)
    return parser.parse_args()


def load_bnodes(path: str) -> List[Dict[str, int]]:
    rows: List[Dict[str, int]] = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "id": int(row["bnode_id"]),
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "support": int(row["support_reads"]),
                }
            )
    return rows


def build_super_nodes(bnodes: List[Dict[str, int]]) -> List[Dict[str, object]]:
    if not bnodes:
        return []

    supports = sorted(row["support"] for row in bnodes)
    q25 = supports[(len(supports) - 1) // 4]
    low_support_max = max(2, min(8, q25))

    super_nodes: List[Dict[str, object]] = []
    run_start = 0
    for idx in range(1, len(bnodes) + 1):
        if idx == len(bnodes):
            should_break = True
        else:
            prev = bnodes[idx - 1]
            nxt = bnodes[idx]
            weak_boundary = (
                prev["support"] <= low_support_max
                or nxt["support"] <= low_support_max
                or prev["support"] * 4 <= nxt["support"]
                or nxt["support"] * 4 <= prev["support"]
            )
            should_break = nxt["start"] != prev["end"] or not weak_boundary

        if should_break:
            members = bnodes[run_start:idx]
            super_nodes.append(
                {
                    "id": len(super_nodes),
                    "start": members[0]["start"],
                    "end": members[-1]["end"],
                    "max_support": max(member["support"] for member in members),
                    "raw_count": len(members),
                }
            )
            run_start = idx

    return super_nodes


def parse_chain(chain_text: str) -> List[Tuple[int, int]]:
    if not chain_text:
        return []
    out: List[Tuple[int, int]] = []
    for part in chain_text.split(","):
        donor, acceptor = part.split("-", 1)
        out.append((int(donor), int(acceptor)))
    return out


def last_acceptor(chain_text: str) -> Optional[int]:
    chain = parse_chain(chain_text)
    if not chain:
        return None
    return chain[-1][1]


def find_super_node(super_nodes: List[Dict[str, object]], pos: Optional[int]) -> Optional[Dict[str, object]]:
    if pos is None:
        return None
    for node in super_nodes:
        if int(node["start"]) <= pos <= int(node["end"]):
            return node
    return None


def load_terminal_support(path: str) -> Dict[str, Dict[int, int]]:
    support: Dict[str, Dict[int, int]] = defaultdict(lambda: defaultdict(int))
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["mode"] != "corrected":
                continue
            super_path = row["super_path"].split(",")
            if not super_path:
                continue
            terminal_super = int(super_path[-1])
            family_key = row["family_key"]
            support[family_key][terminal_super] += int(row["read_count"])
    return support


def classify(
    blocker_super_id: Optional[int],
    gtf_super_id: Optional[int],
    blocker_support: int,
    gtf_support: int,
) -> str:
    if blocker_super_id is None or gtf_super_id is None:
        return "unmapped_terminal"
    if blocker_super_id == gtf_super_id:
        return "within_supernode_boundary"
    if gtf_support == 0:
        return "unsupported_gtf_terminal_supernode"
    if blocker_support == 0:
        return "unsupported_blocker_terminal_supernode"
    return "supported_terminal_supernode_choice"


def main() -> None:
    args = parse_args()

    bnodes = load_bnodes(args.bnodes)
    super_nodes = build_super_nodes(bnodes)
    terminal_support = load_terminal_support(args.transfrags)

    output_stem = Path(args.output_stem)
    summary_path = output_stem.with_suffix(".summary.tsv")
    categories_path = output_stem.with_suffix(".categories.tsv")

    category_counts: Counter[str] = Counter()
    rows: List[Dict[str, object]] = []

    with open(args.exact_neighbors_summary, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            blocker_acceptor = last_acceptor(row["normalized_chain"])
            gtf_acceptor = last_acceptor(row["best_gtf_chain"])
            blocker_super = find_super_node(super_nodes, blocker_acceptor)
            gtf_super = find_super_node(super_nodes, gtf_acceptor)

            family_key = row["candidate_family_keys"]
            family_support = terminal_support.get(family_key, {})
            blocker_super_id = int(blocker_super["id"]) if blocker_super else None
            gtf_super_id = int(gtf_super["id"]) if gtf_super else None
            blocker_terminal_support = family_support.get(blocker_super_id, 0) if blocker_super else 0
            gtf_terminal_support = family_support.get(gtf_super_id, 0) if gtf_super else 0

            category = classify(
                blocker_super_id,
                gtf_super_id,
                blocker_terminal_support,
                gtf_terminal_support,
            )
            category_counts[category] += 1

            rows.append(
                {
                    "rank": row["rank"],
                    "backward_transcripts": row["backward_transcripts"],
                    "category": category,
                    "corrected_family_read_sum": row["corrected_family_read_sum"],
                    "blocker_last_acceptor": blocker_acceptor or "",
                    "blocker_super_id": blocker_super_id if blocker_super_id is not None else "",
                    "blocker_super_start": blocker_super["start"] if blocker_super else "",
                    "blocker_super_end": blocker_super["end"] if blocker_super else "",
                    "blocker_terminal_support": blocker_terminal_support,
                    "gtf_tx_id": row["best_gtf_tx_id"],
                    "gtf_last_acceptor": gtf_acceptor or "",
                    "gtf_super_id": gtf_super_id if gtf_super_id is not None else "",
                    "gtf_super_start": gtf_super["start"] if gtf_super else "",
                    "gtf_super_end": gtf_super["end"] if gtf_super else "",
                    "gtf_terminal_support": gtf_terminal_support,
                    "super_id_delta": (
                        abs(blocker_super_id - gtf_super_id)
                        if blocker_super_id is not None and gtf_super_id is not None
                        else ""
                    ),
                    "relationship": row["relationship"],
                    "blocker_only_junctions": row["blocker_only_junctions"],
                    "gtf_only_junctions": row["gtf_only_junctions"],
                    "candidate_family_keys": family_key,
                }
            )

    with open(summary_path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["exact_blockers", len(rows)])
        for category, count in sorted(category_counts.items()):
            writer.writerow([category, count])

    fieldnames = [
        "rank",
        "backward_transcripts",
        "category",
        "corrected_family_read_sum",
        "blocker_last_acceptor",
        "blocker_super_id",
        "blocker_super_start",
        "blocker_super_end",
        "blocker_terminal_support",
        "gtf_tx_id",
        "gtf_last_acceptor",
        "gtf_super_id",
        "gtf_super_start",
        "gtf_super_end",
        "gtf_terminal_support",
        "super_id_delta",
        "relationship",
        "blocker_only_junctions",
        "gtf_only_junctions",
        "candidate_family_keys",
    ]
    with open(categories_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"exact blockers={len(rows)}")
    for category, count in sorted(category_counts.items()):
        print(f"{category}={count}")
    print(f"wrote {summary_path}")
    print(f"wrote {categories_path}")


if __name__ == "__main__":
    main()
