#!/usr/bin/env python3
import argparse
import csv
from collections import Counter

import pysam


def parse_chain(text: str):
    out = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        left, right = item.split("-")
        out.append((int(left), int(right)))
    return out


def parse_targets(text: str):
    out = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        name, chain = item.split(":", 1)
        left, right = chain.split("-")
        out.append((name, int(left), int(right)))
    return out


def read_introns(read, chrom):
    if read.is_unmapped or read.reference_name != chrom:
        return []
    ref = read.reference_start + 1
    introns = []
    for op, length in read.cigartuples or []:
        if op == 3 or (op == 2 and 5 <= length <= 1_000_000):
            donor = ref - 1
            ref += length
            acceptor = ref
            introns.append((donor, acceptor))
        elif op in (0, 2, 7, 8):
            ref += length
    return introns


def prefix_match_count(introns, prefix, tol):
    matched = 0
    for observed, expected in zip(introns, prefix):
        if (
            abs(observed[0] - expected[0]) <= tol
            and abs(observed[1] - expected[1]) <= tol
        ):
            matched += 1
        else:
            break
    return matched


def load_edge_support(path):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    out = {}
    for row in rows:
        if "from_pos" in row:
            left = int(row["from_pos"])
            right = int(row["to_pos"])
            support = row.get("support", "0")
        elif "left_coord" in row:
            left = int(row["left_coord"])
            right = int(row["right_coord"])
            support = row.get("support", "0")
        else:
            left = int(row["from"])
            right = int(row["to"])
            support = row.get("support", row.get("isoform_count", "0"))
        out[(left, right)] = support
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--region-start", type=int, required=True)
    ap.add_argument("--region-end", type=int, required=True)
    ap.add_argument("--prefix", required=True, help="comma-separated donor-acceptor list")
    ap.add_argument(
        "--targets",
        required=True,
        help="comma-separated name:donor-acceptor list",
    )
    ap.add_argument(
        "--prefix-threshold",
        type=int,
        action="append",
        default=[],
        help="Prefix match count threshold to report",
    )
    ap.add_argument(
        "--prefix-tolerance",
        type=int,
        default=0,
        help="Coordinate tolerance for prefix matching",
    )
    ap.add_argument(
        "--edge-table",
        action="append",
        default=[],
        help="label=path to edge TSV",
    )
    ap.add_argument("--output", required=True)
    ap.add_argument("--notes-output", required=True)
    args = ap.parse_args()

    prefix = parse_chain(args.prefix)
    targets = parse_targets(args.targets)
    thresholds = sorted(set(args.prefix_threshold))
    total = Counter()
    by_threshold = {thr: Counter() for thr in thresholds}
    direct_after = Counter()

    bam = pysam.AlignmentFile(args.bam, "rb")
    for read in bam.fetch(args.chrom, args.region_start - 1, args.region_end):
        introns = read_introns(read, args.chrom)
        if not introns:
            continue
        observed = set(introns)
        prefix_count = prefix_match_count(introns, prefix, args.prefix_tolerance)
        for name, left, right in targets:
            if (left, right) in observed:
                total[name] += 1
                for thr in thresholds:
                    if prefix_count >= thr:
                        by_threshold[thr][name] += 1
        for thr in thresholds:
            if prefix_count >= thr and (
                len(introns) == prefix_count
                or (len(introns) > prefix_count and introns[prefix_count][0] > prefix[-1][1])
            ):
                direct_after[thr] += 1
    bam.close()

    edge_tables = {}
    for item in args.edge_table:
        label, path = item.split("=", 1)
        edge_tables[label] = load_edge_support(path)

    headers = ["feature", "left", "right", "bam_total"]
    headers.extend(f"bam_prefix_ge_{thr}" for thr in thresholds)
    headers.extend(edge_tables.keys())

    with open(args.output, "w") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        for name, left, right in targets:
            row = [name, left, right, total.get(name, 0)]
            row.extend(by_threshold[thr].get(name, 0) for thr in thresholds)
            row.extend(edge_tables[label].get((left, right), 0) for label in edge_tables)
            writer.writerow(row)

    with open(args.notes_output, "w") as fh:
        for thr in thresholds:
            fh.write(f"direct_after_prefix_ge_{thr}\t{direct_after[thr]}\n")


if __name__ == "__main__":
    main()
