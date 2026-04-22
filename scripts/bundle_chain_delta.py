#!/usr/bin/env python3

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


@dataclass(frozen=True)
class BNode:
    id: int
    start: int
    end: int
    support_reads: int


@dataclass(frozen=True)
class SuperBNode:
    id: int
    start: int
    end: int
    support_reads: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare bundle-local old-path and scoped-GTF chains against bridge transfrags."
    )
    parser.add_argument("--bundle-id", type=int, default=904)
    parser.add_argument("--bnodes", required=True)
    parser.add_argument("--oldpath", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--transfrags", required=True)
    parser.add_argument("--families", required=True)
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


def parse_chain(chain: str) -> List[Tuple[int, int]]:
    if not chain:
        return []
    out = []
    for pair in chain.split(","):
        donor, acceptor = pair.split("-", 1)
        out.append((int(donor), int(acceptor)))
    return out


def format_chain(chain: Sequence[Tuple[int, int]]) -> str:
    return ",".join(f"{donor}-{acceptor}" for donor, acceptor in chain)


def parse_exons(exons_text: str) -> List[Tuple[int, int]]:
    out = []
    for exon in exons_text.split(","):
        start, end = exon.split("-", 1)
        out.append((int(start), int(end)))
    return out


def normalize_old_chain(chain: str) -> str:
    normalized = [(donor, acceptor + 1) for donor, acceptor in parse_chain(chain)]
    return format_chain(normalized)


def normalize_old_exons(exons_text: str) -> List[Tuple[int, int]]:
    return [(start + 1, end) for start, end in parse_exons(exons_text)]


def adaptive_low_support_max_from_supports(supports: Sequence[int]) -> int:
    if not supports:
        return 2
    ordered = sorted(supports)
    q25 = ordered[(len(ordered) - 1) // 4]
    return max(2, min(8, q25))


def build_super_bnodes(bnodes: Sequence[BNode]) -> Tuple[List[SuperBNode], List[int]]:
    if not bnodes:
        return [], []

    low_support_max = adaptive_low_support_max_from_supports(
        [node.support_reads for node in bnodes]
    )
    raw_to_super = [0] * len(bnodes)
    super_bnodes: List[SuperBNode] = []
    run_start = 0

    for idx in range(1, len(bnodes) + 1):
        if idx == len(bnodes):
            should_break = True
        else:
            prev = bnodes[idx - 1]
            nxt = bnodes[idx]
            weak_boundary = (
                prev.support_reads <= low_support_max
                or nxt.support_reads <= low_support_max
                or prev.support_reads * 4 <= nxt.support_reads
                or nxt.support_reads * 4 <= prev.support_reads
            )
            should_break = nxt.start != prev.end or not weak_boundary

        if should_break:
            members = bnodes[run_start:idx]
            super_id = len(super_bnodes)
            for raw_idx in range(run_start, idx):
                raw_to_super[raw_idx] = super_id
            super_bnodes.append(
                SuperBNode(
                    id=super_id,
                    start=members[0].start,
                    end=members[-1].end,
                    support_reads=max(node.support_reads for node in members),
                )
            )
            run_start = idx

    return super_bnodes, raw_to_super


def remap_path_to_super(path: Sequence[int], raw_to_super: Sequence[int]) -> List[int]:
    super_path: List[int] = []
    for raw_idx in path:
        super_idx = raw_to_super[raw_idx]
        if not super_path or super_path[-1] != super_idx:
            super_path.append(super_idx)
    return super_path


def trim_terminal_low_support(
    path: Sequence[int], node_supports: Sequence[int], terminal_low_support_max: int
) -> Sequence[int]:
    if not path:
        return path
    left = 0
    right = len(path)
    while right - left > 2 and node_supports[path[left]] <= terminal_low_support_max:
        left += 1
    while right - left > 2 and node_supports[path[right - 1]] <= terminal_low_support_max:
        right -= 1
    return path[left:right]


def internal_chain_family_key(
    path: Sequence[int],
    super_bnodes: Sequence[SuperBNode],
    node_supports: Sequence[int],
    terminal_low_support_max: int,
) -> List[int]:
    trimmed = list(trim_terminal_low_support(path, node_supports, terminal_low_support_max))
    if not trimmed:
        return []
    if len(trimmed) <= 2:
        key: List[int] = []
        for idx in trimmed:
            key.extend([super_bnodes[idx].start, super_bnodes[idx].end])
        return key

    key: List[int] = []
    for left_idx, right_idx in zip(trimmed, trimmed[1:]):
        left = super_bnodes[left_idx]
        right = super_bnodes[right_idx]
        if left.end < right.start:
            key.extend([left.end, right.start])

    if not key:
        for idx in trimmed[1:-1]:
            key.extend([super_bnodes[idx].start, super_bnodes[idx].end])

    return key


def segment_covers_interval(segment: Tuple[int, int], interval: Tuple[int, int]) -> bool:
    return segment[0] <= interval[0] and segment[1] >= interval[1] and interval[0] < interval[1]


def path_for_segments(segments: Sequence[Tuple[int, int]], bnodes: Sequence[BNode]) -> List[int]:
    path: List[int] = []
    for node in bnodes:
        interval = (node.start, node.end)
        if any(segment_covers_interval(segment, interval) for segment in segments):
            path.append(node.id)
    return path


def list_to_key(values: Iterable[int]) -> str:
    return ",".join(str(v) for v in values)


def family_key_from_chain(chain: str) -> str:
    flat: List[int] = []
    for donor, acceptor in parse_chain(chain):
        flat.extend([donor, acceptor])
    return list_to_key(flat)


def load_bnodes(path: str, bundle_id: int) -> List[BNode]:
    bnodes: List[BNode] = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if int(row["bundle_id"]) != bundle_id:
                continue
            bnodes.append(
                BNode(
                    id=int(row["bnode_id"]),
                    start=int(row["start"]),
                    end=int(row["end"]),
                    support_reads=int(row["support_reads"]),
                )
            )
    bnodes.sort(key=lambda node: node.id)
    return bnodes


def load_oldpath_entries(
    path: str,
    bnodes: Sequence[BNode],
    super_bnodes: Sequence[SuperBNode],
    raw_to_super: Sequence[int],
    super_supports: Sequence[int],
    super_low_support_max: int,
) -> List[dict]:
    entries: List[dict] = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            raw_chain = row["junction_chain"]
            normalized_chain = normalize_old_chain(raw_chain) if raw_chain else ""
            segments = normalize_old_exons(row["exons"])
            raw_path = path_for_segments(segments, bnodes)
            super_path = remap_path_to_super(raw_path, raw_to_super)
            family_key = list_to_key(
                internal_chain_family_key(
                    super_path,
                    super_bnodes,
                    super_supports,
                    super_low_support_max,
                )
            )
            entries.append(
                {
                    "id": f"old_{row['tx_idx']}",
                    "coverage": float(row["coverage"]),
                    "normalized_chain": normalized_chain,
                    "segments": segments,
                    "raw_path": raw_path,
                    "super_path": super_path,
                    "family_key": family_key,
                }
            )
    return entries


def load_gtf_entries(
    path: str,
    bnodes: Sequence[BNode],
    super_bnodes: Sequence[SuperBNode],
    raw_to_super: Sequence[int],
    super_supports: Sequence[int],
    super_low_support_max: int,
) -> List[dict]:
    transcript_meta: Dict[str, dict] = {}
    transcript_exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)

    with open(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            attrs = parse_attrs(fields[8])
            tx_id = attrs.get("transcript_id")
            if not tx_id:
                continue
            if feature == "transcript":
                transcript_meta[tx_id] = {
                    "id": tx_id,
                    "fl_count": float(attrs.get("FL_count", "0")),
                    "normalized_chain": attrs.get("junction_chain", ""),
                    "flags": attrs.get("flags", ""),
                }
            elif feature == "exon":
                transcript_exons[tx_id].append((start, end))

    entries: List[dict] = []
    for tx_id, meta in transcript_meta.items():
        chain = meta["normalized_chain"]
        if not chain:
            continue
        segments = sorted(transcript_exons.get(tx_id, []))
        raw_path = path_for_segments(segments, bnodes)
        super_path = remap_path_to_super(raw_path, raw_to_super)
        family_key = list_to_key(
            internal_chain_family_key(
                super_path,
                super_bnodes,
                super_supports,
                super_low_support_max,
            )
        )
        entries.append(
            {
                "id": tx_id,
                "coverage": meta["fl_count"],
                "normalized_chain": chain,
                "segments": segments,
                "raw_path": raw_path,
                "super_path": super_path,
                "family_key": family_key,
                "flags": meta["flags"],
            }
        )
    return entries


def aggregate_entries(entries: Sequence[dict]) -> Dict[str, dict]:
    grouped: Dict[str, dict] = {}
    for entry in entries:
        chain = entry["normalized_chain"]
        if not chain:
            continue
        bucket = grouped.setdefault(
            chain,
            {
                "count": 0,
                "ids": [],
                "max_coverage": 0.0,
                "sum_coverage": 0.0,
                "family_keys": set(),
                "super_paths": set(),
                "pairs": set(),
            },
        )
        bucket["count"] += 1
        bucket["ids"].append(entry["id"])
        bucket["max_coverage"] = max(bucket["max_coverage"], entry["coverage"])
        bucket["sum_coverage"] += entry["coverage"]
        family_key = entry["family_key"]
        super_path = list_to_key(entry["super_path"])
        if family_key:
            bucket["family_keys"].add(family_key)
        if super_path:
            bucket["super_paths"].add(super_path)
        if family_key or super_path:
            bucket["pairs"].add((family_key, super_path))
    return grouped


def load_transfrag_support(path: str, bundle_id: int) -> Dict[Tuple[str, str, str], dict]:
    support: Dict[Tuple[str, str, str], dict] = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if int(row["bundle_id"]) != bundle_id:
                continue
            key = (row["mode"], row["family_key"], row["super_path"])
            bucket = support.setdefault(
                key,
                {"transfrag_count": 0, "read_sum": 0, "transfrag_ids": []},
            )
            bucket["transfrag_count"] += 1
            bucket["read_sum"] += int(row["read_count"])
            bucket["transfrag_ids"].append(row["transfrag_id"])
    return support


def load_family_support(path: str, bundle_id: int) -> Dict[Tuple[str, str], dict]:
    support: Dict[Tuple[str, str], dict] = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if int(row["bundle_id"]) != bundle_id:
                continue
            support[(row["mode"], row["family_key"])] = {
                "read_count": int(row["read_count"]),
                "super_path_count": int(row["super_path_count"]),
                "boundary_start_min": int(row["boundary_start_min"]),
                "boundary_start_max": int(row["boundary_start_max"]),
                "boundary_end_min": int(row["boundary_end_min"]),
                "boundary_end_max": int(row["boundary_end_max"]),
            }
    return support


def aggregate_pair_support(
    pairs: Iterable[Tuple[str, str]],
    transfrag_support: Dict[Tuple[str, str, str], dict],
    mode: str,
) -> dict:
    out = {"transfrag_count": 0, "read_sum": 0, "pair_hits": [], "transfrag_ids": []}
    for family_key, super_path in sorted(set(pairs)):
        bucket = transfrag_support.get((mode, family_key, super_path))
        if not bucket:
            continue
        out["transfrag_count"] += bucket["transfrag_count"]
        out["read_sum"] += bucket["read_sum"]
        out["pair_hits"].append(f"{family_key}|{super_path}")
        out["transfrag_ids"].extend(bucket["transfrag_ids"])
    return out


def aggregate_family_support(
    family_keys: Iterable[str], family_support: Dict[Tuple[str, str], dict], mode: str
) -> dict:
    out = {
        "family_count": 0,
        "read_sum": 0,
        "super_path_sum": 0,
        "family_hits": [],
        "boundary_spans": [],
    }
    for family_key in sorted(set(family_keys)):
        if not family_key:
            continue
        bucket = family_support.get((mode, family_key))
        if not bucket:
            continue
        out["family_count"] += 1
        out["read_sum"] += bucket["read_count"]
        out["super_path_sum"] += bucket["super_path_count"]
        out["family_hits"].append(family_key)
        out["boundary_spans"].append(
            f"{bucket['boundary_start_min']}-{bucket['boundary_start_max']}|"
            f"{bucket['boundary_end_min']}-{bucket['boundary_end_max']}"
        )
    return out


def presence_state(raw_count: int, corrected_count: int) -> str:
    if raw_count and corrected_count:
        return "both"
    if raw_count:
        return "raw_only"
    if corrected_count:
        return "corrected_only"
    return "neither"


def write_summary(
    path: Path,
    old_grouped: Dict[str, dict],
    gtf_grouped: Dict[str, dict],
    overlap: int,
    old_only: Sequence[str],
    gtf_only: Sequence[str],
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["old_unique_multi_exon_chains", len(old_grouped)])
        writer.writerow(["gtf_unique_multi_exon_chains", len(gtf_grouped)])
        writer.writerow(["normalized_exact_overlap", overlap])
        writer.writerow(["old_only_chains", len(old_only)])
        writer.writerow(["gtf_only_chains", len(gtf_only)])


def write_deltas(
    path: Path,
    old_grouped: Dict[str, dict],
    gtf_grouped: Dict[str, dict],
    old_only: Sequence[str],
    gtf_only: Sequence[str],
    transfrag_support: Dict[Tuple[str, str, str], dict],
    family_support: Dict[Tuple[str, str], dict],
) -> None:
    fields = [
        "category",
        "normalized_chain",
        "junction_count",
        "old_tx_count",
        "old_tx_ids",
        "old_max_coverage",
        "old_sum_coverage",
        "old_family_keys",
        "old_super_paths",
        "gtf_tx_count",
        "gtf_tx_ids",
        "gtf_max_fl_count",
        "gtf_sum_fl_count",
        "gtf_family_keys",
        "gtf_super_paths",
        "candidate_family_keys",
        "candidate_super_paths",
        "raw_pair_state",
        "raw_pair_transfrags",
        "raw_pair_read_sum",
        "raw_pair_hits",
        "corrected_pair_state",
        "corrected_pair_transfrags",
        "corrected_pair_read_sum",
        "corrected_pair_hits",
        "family_state",
        "raw_family_matches",
        "raw_family_read_sum",
        "raw_family_super_path_sum",
        "corrected_family_matches",
        "corrected_family_read_sum",
        "corrected_family_super_path_sum",
        "raw_boundary_spans",
        "corrected_boundary_spans",
    ]

    rows = []
    for category, chains in (("old_only", old_only), ("gtf_only", gtf_only)):
        for chain in chains:
            old_bucket = old_grouped.get(chain, {})
            gtf_bucket = gtf_grouped.get(chain, {})
            candidate_pairs = old_bucket.get("pairs") or gtf_bucket.get("pairs") or set()
            candidate_family_keys = old_bucket.get("family_keys") or gtf_bucket.get("family_keys") or set()
            candidate_super_paths = old_bucket.get("super_paths") or gtf_bucket.get("super_paths") or set()

            raw_pair = aggregate_pair_support(candidate_pairs, transfrag_support, "raw")
            corrected_pair = aggregate_pair_support(candidate_pairs, transfrag_support, "corrected")
            raw_family = aggregate_family_support(candidate_family_keys, family_support, "raw")
            corrected_family = aggregate_family_support(candidate_family_keys, family_support, "corrected")

            rows.append(
                {
                    "category": category,
                    "normalized_chain": chain,
                    "junction_count": len(parse_chain(chain)),
                    "old_tx_count": old_bucket.get("count", 0),
                    "old_tx_ids": ",".join(old_bucket.get("ids", [])),
                    "old_max_coverage": f"{old_bucket.get('max_coverage', 0.0):.4f}",
                    "old_sum_coverage": f"{old_bucket.get('sum_coverage', 0.0):.4f}",
                    "old_family_keys": ";".join(sorted(old_bucket.get("family_keys", set()))),
                    "old_super_paths": ";".join(sorted(old_bucket.get("super_paths", set()))),
                    "gtf_tx_count": gtf_bucket.get("count", 0),
                    "gtf_tx_ids": ",".join(gtf_bucket.get("ids", [])),
                    "gtf_max_fl_count": f"{gtf_bucket.get('max_coverage', 0.0):.4f}",
                    "gtf_sum_fl_count": f"{gtf_bucket.get('sum_coverage', 0.0):.4f}",
                    "gtf_family_keys": ";".join(sorted(gtf_bucket.get("family_keys", set()))),
                    "gtf_super_paths": ";".join(sorted(gtf_bucket.get("super_paths", set()))),
                    "candidate_family_keys": ";".join(sorted(candidate_family_keys)),
                    "candidate_super_paths": ";".join(sorted(candidate_super_paths)),
                    "raw_pair_state": "present" if raw_pair["transfrag_count"] else "absent",
                    "raw_pair_transfrags": raw_pair["transfrag_count"],
                    "raw_pair_read_sum": raw_pair["read_sum"],
                    "raw_pair_hits": ";".join(raw_pair["pair_hits"]),
                    "corrected_pair_state": "present" if corrected_pair["transfrag_count"] else "absent",
                    "corrected_pair_transfrags": corrected_pair["transfrag_count"],
                    "corrected_pair_read_sum": corrected_pair["read_sum"],
                    "corrected_pair_hits": ";".join(corrected_pair["pair_hits"]),
                    "family_state": presence_state(raw_family["family_count"], corrected_family["family_count"]),
                    "raw_family_matches": raw_family["family_count"],
                    "raw_family_read_sum": raw_family["read_sum"],
                    "raw_family_super_path_sum": raw_family["super_path_sum"],
                    "corrected_family_matches": corrected_family["family_count"],
                    "corrected_family_read_sum": corrected_family["read_sum"],
                    "corrected_family_super_path_sum": corrected_family["super_path_sum"],
                    "raw_boundary_spans": ";".join(raw_family["boundary_spans"]),
                    "corrected_boundary_spans": ";".join(corrected_family["boundary_spans"]),
                }
            )

    rows.sort(
        key=lambda row: (
            row["category"],
            -int(row["corrected_family_read_sum"]),
            -int(row["raw_family_read_sum"]),
            row["normalized_chain"],
        )
    )

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()

    bnodes = load_bnodes(args.bnodes, args.bundle_id)
    super_bnodes, raw_to_super = build_super_bnodes(bnodes)
    super_supports = [node.support_reads for node in super_bnodes]
    super_low_support_max = adaptive_low_support_max_from_supports(super_supports)

    old_entries = load_oldpath_entries(
        args.oldpath,
        bnodes,
        super_bnodes,
        raw_to_super,
        super_supports,
        super_low_support_max,
    )
    gtf_entries = load_gtf_entries(
        args.gtf,
        bnodes,
        super_bnodes,
        raw_to_super,
        super_supports,
        super_low_support_max,
    )

    old_grouped = aggregate_entries(old_entries)
    gtf_grouped = aggregate_entries(gtf_entries)
    old_chains = set(old_grouped)
    gtf_chains = set(gtf_grouped)

    overlap = len(old_chains & gtf_chains)
    old_only = sorted(old_chains - gtf_chains)
    gtf_only = sorted(gtf_chains - old_chains)

    transfrag_support = load_transfrag_support(args.transfrags, args.bundle_id)
    family_support = load_family_support(args.families, args.bundle_id)

    output_stem = Path(args.output_stem)
    write_summary(output_stem.with_suffix(".summary.tsv"), old_grouped, gtf_grouped, overlap, old_only, gtf_only)
    write_deltas(
        output_stem.with_suffix(".deltas.tsv"),
        old_grouped,
        gtf_grouped,
        old_only,
        gtf_only,
        transfrag_support,
        family_support,
    )

    print(
        f"bundle {args.bundle_id}: old={len(old_grouped)} gtf={len(gtf_grouped)} "
        f"overlap={overlap} old_only={len(old_only)} gtf_only={len(gtf_only)}"
    )
    print(f"wrote {output_stem.with_suffix('.summary.tsv')}")
    print(f"wrote {output_stem.with_suffix('.deltas.tsv')}")


if __name__ == "__main__":
    main()
