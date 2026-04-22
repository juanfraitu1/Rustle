#!/usr/bin/env python3

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple


Junction = Tuple[int, int]


@dataclass
class Blocker:
    rank: int
    backward_transcripts: str
    backward_genes: str
    corrected_family_read_sum: int
    normalized_chain: str
    category: str
    old_tx_ids: str
    candidate_family_keys: str
    candidate_super_paths: str


@dataclass
class GtfTranscript:
    transcript_id: str
    gene_id: str
    fl_count: int
    start: int
    end: int
    chain: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find the nearest scoped-GTF neighbors for exact backward-supported blocker chains."
    )
    parser.add_argument("--blockers-tsv", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output-stem", required=True)
    parser.add_argument("--top-k", type=int, default=3)
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


def parse_chain(chain_text: str) -> List[Junction]:
    if not chain_text:
        return []
    out: List[Junction] = []
    for part in chain_text.split(","):
        donor, acceptor = part.split("-", 1)
        out.append((int(donor), int(acceptor)))
    return out


def format_chain(chain: Sequence[Junction]) -> str:
    return ",".join(f"{donor}-{acceptor}" for donor, acceptor in chain)


def load_exact_blockers(path: str) -> List[Blocker]:
    blockers: List[Blocker] = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["backward_match_state"] != "exact_match":
                continue
            blockers.append(
                Blocker(
                    rank=int(row["rank"]),
                    backward_transcripts=row["backward_transcripts"],
                    backward_genes=row["backward_genes"],
                    corrected_family_read_sum=int(row["corrected_family_read_sum"]),
                    normalized_chain=row["normalized_chain"],
                    category=row["category"],
                    old_tx_ids=row["old_tx_ids"],
                    candidate_family_keys=row["candidate_family_keys"],
                    candidate_super_paths=row["candidate_super_paths"],
                )
            )
    blockers.sort(key=lambda blocker: blocker.rank)
    return blockers


def load_gtf_transcripts(path: str) -> List[GtfTranscript]:
    transcripts: List[GtfTranscript] = []
    with open(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attrs = parse_attrs(fields[8])
            transcripts.append(
                GtfTranscript(
                    transcript_id=attrs.get("transcript_id", ""),
                    gene_id=attrs.get("gene_id", ""),
                    fl_count=int(attrs.get("FL_count", "0")),
                    start=int(fields[3]),
                    end=int(fields[4]),
                    chain=attrs.get("junction_chain", ""),
                )
            )
    return transcripts


def prefix_shared(a: Sequence[Junction], b: Sequence[Junction]) -> int:
    count = 0
    for left, right in zip(a, b):
        if left != right:
            break
        count += 1
    return count


def suffix_shared(a: Sequence[Junction], b: Sequence[Junction], prefix: int) -> int:
    count = 0
    max_suffix = min(len(a), len(b)) - prefix
    for offset in range(1, max_suffix + 1):
        if a[-offset] != b[-offset]:
            break
        count += 1
    return count


def classify_relationship(
    blocker_chain: Sequence[Junction],
    gtf_chain: Sequence[Junction],
    prefix: int,
    suffix: int,
    blocker_only: Sequence[Junction],
    gtf_only: Sequence[Junction],
    max_shift: int,
) -> str:
    if blocker_chain == gtf_chain:
        return "exact"

    overlap_len = min(len(blocker_chain), len(gtf_chain))
    unmatched_middle = overlap_len - prefix - suffix
    if len(blocker_only) == 1 and len(gtf_only) == 1 and max_shift <= 10:
        return "microshift"
    if prefix == overlap_len and len(blocker_chain) > len(gtf_chain):
        return "blocker_terminal_extension"
    if prefix == overlap_len and len(gtf_chain) > len(blocker_chain):
        return "gtf_terminal_extension"
    if prefix + suffix == overlap_len and blocker_only and gtf_only:
        return "terminal_mode_substitution"
    if unmatched_middle == 0 and blocker_only and not gtf_only:
        return "blocker_terminal_only"
    if unmatched_middle == 0 and gtf_only and not blocker_only:
        return "gtf_terminal_only"
    return "internal_substitution"


def compare_chains(
    blocker_chain_text: str, gtf_chain_text: str
) -> Dict[str, object]:
    blocker_chain = parse_chain(blocker_chain_text)
    gtf_chain = parse_chain(gtf_chain_text)
    blocker_set = set(blocker_chain)
    gtf_set = set(gtf_chain)
    shared = sorted(blocker_set & gtf_set)
    blocker_only = sorted(blocker_set - gtf_set)
    gtf_only = sorted(gtf_set - blocker_set)
    prefix = prefix_shared(blocker_chain, gtf_chain)
    suffix = suffix_shared(blocker_chain, gtf_chain, prefix)

    sum_shift = 0
    max_shift = 0
    overlap_len = min(len(blocker_chain), len(gtf_chain))
    for idx in range(overlap_len):
        if blocker_chain[idx] == gtf_chain[idx]:
            continue
        donor_shift = abs(blocker_chain[idx][0] - gtf_chain[idx][0])
        acceptor_shift = abs(blocker_chain[idx][1] - gtf_chain[idx][1])
        sum_shift += donor_shift + acceptor_shift
        max_shift = max(max_shift, donor_shift, acceptor_shift)

    return {
        "shared_count": len(shared),
        "blocker_only_count": len(blocker_only),
        "gtf_only_count": len(gtf_only),
        "shared_junctions": format_chain(shared),
        "blocker_only_junctions": format_chain(blocker_only),
        "gtf_only_junctions": format_chain(gtf_only),
        "prefix_shared": prefix,
        "suffix_shared": suffix,
        "sum_coord_shift": sum_shift,
        "max_coord_shift": max_shift,
        "relationship": classify_relationship(
            blocker_chain,
            gtf_chain,
            prefix,
            suffix,
            blocker_only,
            gtf_only,
            max_shift,
        ),
    }


def main() -> None:
    args = parse_args()

    blockers = load_exact_blockers(args.blockers_tsv)
    transcripts = [tx for tx in load_gtf_transcripts(args.gtf) if tx.chain]

    output_stem = Path(args.output_stem)
    summary_path = output_stem.with_suffix(".summary.tsv")
    neighbors_path = output_stem.with_suffix(".neighbors.tsv")

    summary_rows: List[Dict[str, object]] = []
    neighbor_rows: List[Dict[str, object]] = []

    for blocker in blockers:
        comparisons: List[Dict[str, object]] = []
        blocker_chain = parse_chain(blocker.normalized_chain)
        blocker_start = blocker_chain[0][0] if blocker_chain else 0
        blocker_end = blocker_chain[-1][1] if blocker_chain else 0

        for tx in transcripts:
            comparison = compare_chains(blocker.normalized_chain, tx.chain)
            comparison.update(
                {
                    "gtf_tx_id": tx.transcript_id,
                    "gtf_gene_id": tx.gene_id,
                    "gtf_fl_count": tx.fl_count,
                    "gtf_start": tx.start,
                    "gtf_end": tx.end,
                    "gtf_chain": tx.chain,
                    "span_shift": abs(blocker_start - tx.start) + abs(blocker_end - tx.end),
                }
            )
            comparisons.append(comparison)

        comparisons.sort(
            key=lambda row: (
                -int(row["shared_count"]),
                int(row["blocker_only_count"]) + int(row["gtf_only_count"]),
                int(row["sum_coord_shift"]),
                int(row["span_shift"]),
                -int(row["gtf_fl_count"]),
                row["gtf_tx_id"],
            )
        )

        best = comparisons[0]
        summary_rows.append(
            {
                "rank": blocker.rank,
                "backward_transcripts": blocker.backward_transcripts,
                "backward_genes": blocker.backward_genes,
                "corrected_family_read_sum": blocker.corrected_family_read_sum,
                "category": blocker.category,
                "old_tx_ids": blocker.old_tx_ids,
                "normalized_chain": blocker.normalized_chain,
                "best_gtf_tx_id": best["gtf_tx_id"],
                "best_gtf_fl_count": best["gtf_fl_count"],
                "best_gtf_chain": best["gtf_chain"],
                "shared_count": best["shared_count"],
                "blocker_only_count": best["blocker_only_count"],
                "gtf_only_count": best["gtf_only_count"],
                "prefix_shared": best["prefix_shared"],
                "suffix_shared": best["suffix_shared"],
                "max_coord_shift": best["max_coord_shift"],
                "sum_coord_shift": best["sum_coord_shift"],
                "span_shift": best["span_shift"],
                "relationship": best["relationship"],
                "blocker_only_junctions": best["blocker_only_junctions"],
                "gtf_only_junctions": best["gtf_only_junctions"],
                "candidate_family_keys": blocker.candidate_family_keys,
                "candidate_super_paths": blocker.candidate_super_paths,
            }
        )

        for neighbor_rank, neighbor in enumerate(comparisons[: args.top_k], start=1):
            neighbor_rows.append(
                {
                    "rank": blocker.rank,
                    "neighbor_rank": neighbor_rank,
                    "backward_transcripts": blocker.backward_transcripts,
                    "corrected_family_read_sum": blocker.corrected_family_read_sum,
                    "normalized_chain": blocker.normalized_chain,
                    "gtf_tx_id": neighbor["gtf_tx_id"],
                    "gtf_fl_count": neighbor["gtf_fl_count"],
                    "gtf_chain": neighbor["gtf_chain"],
                    "shared_count": neighbor["shared_count"],
                    "blocker_only_count": neighbor["blocker_only_count"],
                    "gtf_only_count": neighbor["gtf_only_count"],
                    "prefix_shared": neighbor["prefix_shared"],
                    "suffix_shared": neighbor["suffix_shared"],
                    "max_coord_shift": neighbor["max_coord_shift"],
                    "sum_coord_shift": neighbor["sum_coord_shift"],
                    "span_shift": neighbor["span_shift"],
                    "relationship": neighbor["relationship"],
                    "blocker_only_junctions": neighbor["blocker_only_junctions"],
                    "gtf_only_junctions": neighbor["gtf_only_junctions"],
                }
            )

    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "backward_transcripts",
                "backward_genes",
                "corrected_family_read_sum",
                "category",
                "old_tx_ids",
                "normalized_chain",
                "best_gtf_tx_id",
                "best_gtf_fl_count",
                "best_gtf_chain",
                "shared_count",
                "blocker_only_count",
                "gtf_only_count",
                "prefix_shared",
                "suffix_shared",
                "max_coord_shift",
                "sum_coord_shift",
                "span_shift",
                "relationship",
                "blocker_only_junctions",
                "gtf_only_junctions",
                "candidate_family_keys",
                "candidate_super_paths",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    with neighbors_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "rank",
                "neighbor_rank",
                "backward_transcripts",
                "corrected_family_read_sum",
                "normalized_chain",
                "gtf_tx_id",
                "gtf_fl_count",
                "gtf_chain",
                "shared_count",
                "blocker_only_count",
                "gtf_only_count",
                "prefix_shared",
                "suffix_shared",
                "max_coord_shift",
                "sum_coord_shift",
                "span_shift",
                "relationship",
                "blocker_only_junctions",
                "gtf_only_junctions",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(neighbor_rows)

    print(f"exact backward-supported blockers={len(summary_rows)}")
    print(f"wrote {summary_path}")
    print(f"wrote {neighbors_path}")


if __name__ == "__main__":
    main()
