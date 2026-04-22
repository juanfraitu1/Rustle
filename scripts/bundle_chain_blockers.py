#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Rank bundle chain discrepancies against backward StringTie transcript support."
    )
    parser.add_argument("--delta-tsv", required=True)
    parser.add_argument("--stringtie-gtf", required=True)
    parser.add_argument("--supported-gfa-dir", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
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


def load_overlapping_stringtie_transcripts(
    gtf_path: str, chrom: str, start: int, end: int
) -> Tuple[Dict[str, dict], Dict[str, List[str]]]:
    transcript_meta: Dict[str, dict] = {}
    transcript_exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    with open(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] != chrom:
                continue
            feature = fields[2]
            feature_start = int(fields[3])
            feature_end = int(fields[4])
            attrs = parse_attrs(fields[8])
            tx_id = attrs.get("transcript_id")
            if not tx_id:
                continue

            if feature == "transcript":
                if feature_start > end or feature_end < start:
                    continue
                transcript_meta[tx_id] = {
                    "gene_id": attrs.get("gene_id", ""),
                    "start": feature_start,
                    "end": feature_end,
                }
            elif feature == "exon":
                transcript_exons[tx_id].append((feature_start, feature_end))

    chain_to_txs: Dict[str, List[str]] = defaultdict(list)
    tx_info: Dict[str, dict] = {}
    for tx_id, meta in transcript_meta.items():
        exons = sorted(transcript_exons.get(tx_id, []))
        if len(exons) < 2:
            continue
        chain = ",".join(f"{exons[i][1]}-{exons[i+1][0]}" for i in range(len(exons) - 1))
        tx_info[tx_id] = {"gene_id": meta["gene_id"], "chain": chain}
        chain_to_txs[chain].append(tx_id)
    return tx_info, chain_to_txs


def load_gfa_support(gfa_dir: str, genes: List[str]) -> Dict[str, dict]:
    support: Dict[str, dict] = {}
    for gene_id in genes:
        gfa_path = Path(gfa_dir) / f"{gene_id}.supported.gfa"
        if not gfa_path.exists():
            continue
        with gfa_path.open() as handle:
            for line in handle:
                if not line.startswith("P\t"):
                    continue
                fields = line.rstrip("\n").split("\t")
                tx_id = fields[1]
                fsm = 0
                ism = 0
                for tag in fields[4:]:
                    if tag.startswith("FSM:i:"):
                        fsm = int(tag.split(":")[-1])
                    elif tag.startswith("ISM:i:"):
                        ism = int(tag.split(":")[-1])
                support[tx_id] = {"fsm": fsm, "ism": ism}
    return support


def write_summary(path: Path, rows: List[dict]) -> None:
    matched = [row for row in rows if row["backward_match_state"] == "exact_match"]
    unmatched = [row for row in rows if row["backward_match_state"] == "no_exact_match"]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["same_bridge_blockers", len(rows)])
        writer.writerow(["exact_backward_chain_matches", len(matched)])
        writer.writerow(["no_exact_backward_chain_match", len(unmatched)])
        writer.writerow(["matched_with_fsm_gt_0", sum(1 for row in matched if int(row["backward_total_fsm"]) > 0)])
        writer.writerow(["matched_with_ism_gt_0", sum(1 for row in matched if int(row["backward_total_ism"]) > 0)])


def main() -> None:
    args = parse_args()

    tx_info, chain_to_txs = load_overlapping_stringtie_transcripts(
        args.stringtie_gtf, args.chrom, args.start, args.end
    )
    gfa_support = load_gfa_support(
        args.supported_gfa_dir,
        sorted({info["gene_id"] for info in tx_info.values()}),
    )

    rows: List[dict] = []
    with open(args.delta_tsv, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["corrected_pair_state"] != "present":
                continue

            chain = row["normalized_chain"]
            matching_txs = chain_to_txs.get(chain, [])
            matching_genes = sorted({tx_info[tx]["gene_id"] for tx in matching_txs})
            total_fsm = sum(gfa_support.get(tx, {}).get("fsm", 0) for tx in matching_txs)
            total_ism = sum(gfa_support.get(tx, {}).get("ism", 0) for tx in matching_txs)
            best_fsm = max((gfa_support.get(tx, {}).get("fsm", 0) for tx in matching_txs), default=0)
            best_ism = max((gfa_support.get(tx, {}).get("ism", 0) for tx in matching_txs), default=0)

            rows.append(
                {
                    "category": row["category"],
                    "normalized_chain": chain,
                    "junction_count": row["junction_count"],
                    "corrected_family_read_sum": row["corrected_family_read_sum"],
                    "raw_family_read_sum": row["raw_family_read_sum"],
                    "corrected_pair_read_sum": row["corrected_pair_read_sum"],
                    "corrected_pair_transfrags": row["corrected_pair_transfrags"],
                    "old_tx_ids": row["old_tx_ids"],
                    "gtf_tx_ids": row["gtf_tx_ids"],
                    "candidate_family_keys": row["candidate_family_keys"],
                    "candidate_super_paths": row["candidate_super_paths"],
                    "backward_match_state": "exact_match" if matching_txs else "no_exact_match",
                    "backward_genes": ",".join(matching_genes),
                    "backward_transcripts": ",".join(matching_txs),
                    "backward_total_fsm": total_fsm,
                    "backward_total_ism": total_ism,
                    "backward_best_fsm": best_fsm,
                    "backward_best_ism": best_ism,
                }
            )

    rows.sort(
        key=lambda row: (
            -int(row["corrected_family_read_sum"]),
            -int(row["raw_family_read_sum"]),
            0 if row["backward_match_state"] == "exact_match" else 1,
            -int(row["backward_total_fsm"]),
            -int(row["backward_total_ism"]),
            row["category"],
            row["normalized_chain"],
        )
    )

    for idx, row in enumerate(rows, start=1):
        row["rank"] = idx

    output_stem = Path(args.output_stem)
    summary_path = output_stem.with_suffix(".summary.tsv")
    blockers_path = output_stem.with_suffix(".blockers.tsv")

    write_summary(summary_path, rows)

    fields = [
        "rank",
        "category",
        "normalized_chain",
        "junction_count",
        "corrected_family_read_sum",
        "raw_family_read_sum",
        "corrected_pair_read_sum",
        "corrected_pair_transfrags",
        "old_tx_ids",
        "gtf_tx_ids",
        "candidate_family_keys",
        "candidate_super_paths",
        "backward_match_state",
        "backward_genes",
        "backward_transcripts",
        "backward_total_fsm",
        "backward_total_ism",
        "backward_best_fsm",
        "backward_best_ism",
    ]

    with blockers_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"same-bridge blockers={len(rows)} exact_backward_matches="
        f"{sum(1 for row in rows if row['backward_match_state'] == 'exact_match')}"
    )
    print(f"wrote {summary_path}")
    print(f"wrote {blockers_path}")


if __name__ == "__main__":
    main()
