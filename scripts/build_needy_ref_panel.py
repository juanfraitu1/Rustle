#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

from analyze_gffcompare_gap import parse_gff_reference


def locus_id(transcript_id: str) -> str:
    if transcript_id.startswith("STRG."):
        parts = transcript_id.split(".")
        if len(parts) >= 2:
            return ".".join(parts[:2])
    return transcript_id


def load_subopt_classes(path: Path) -> dict[str, str]:
    by_tx: dict[str, str] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tid = row["transcript_id"]
            by_tx[tid] = row["dominant_class"]
    return by_tx


def join_items(values: Iterable[str]) -> str:
    items = [v for v in values if v]
    return ",".join(items)


def parse_transcript_id(attr_field: str) -> str | None:
    match = re.search(r'transcript_id "([^"]+)"', attr_field)
    if match:
        return match.group(1)
    return None


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build a ranked needy-reference locus panel from missing-ref analysis outputs."
    )
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--ref-buckets-tsv", type=Path, required=True)
    ap.add_argument("--subopt-tsv", type=Path, required=True)
    ap.add_argument("--out-prefix", type=Path, required=True)
    ap.add_argument("--top-loci", type=int, default=15)
    ap.add_argument("--top-refs-per-locus", type=int, default=6)
    ap.add_argument("--locus-slack-bp", type=int, default=2000)
    args = ap.parse_args()

    ref_by_id = parse_gff_reference(args.reference)
    dominant_class_by_tx = load_subopt_classes(args.subopt_tsv)

    locus_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    bucket_counts: Counter[str] = Counter()
    class_counts: Counter[str] = Counter()

    with args.ref_buckets_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            tid = row["transcript_id"]
            lid = locus_id(tid)
            row = dict(row)
            row["locus_id"] = lid
            row["dominant_class"] = row["dominant_class"] or dominant_class_by_tx.get(tid, "")
            locus_rows[lid].append(row)
            bucket_counts[row["bucket"]] += 1
            if row["dominant_class"]:
                class_counts[row["dominant_class"]] += 1

    panel_rows: list[dict[str, object]] = []
    for lid, rows in locus_rows.items():
        refs = [ref_by_id[row["transcript_id"]] for row in rows if row["transcript_id"] in ref_by_id]
        if not refs:
            continue
        chrom = refs[0].chrom
        strand = refs[0].strand
        tx_start = min(ref.start for ref in refs)
        tx_end = max(ref.end for ref in refs)
        panel_start = max(1, tx_start - args.locus_slack_bp)
        panel_end = tx_end + args.locus_slack_bp
        absent = sum(1 for row in rows if row["recovery_category"] == "absent")
        subopt = sum(1 for row in rows if row["recovery_category"] == "subopt_only")
        class_counter = Counter(
            row["dominant_class"] for row in rows if row.get("dominant_class")
        )
        bucket_counter = Counter(row["bucket"] for row in rows)
        priority = (
            absent * 5
            + subopt * 4
            + class_counter.get("c", 0) * 3
            + class_counter.get("j", 0) * 2
            + class_counter.get("k", 0) * 2
            + sum(1 for row in rows if row["bucket"] == "absent_full_junction_overlap_query")
        )
        refs_sorted = sorted(
            rows,
            key=lambda row: (
                row["recovery_category"] != "absent",
                row.get("dominant_class", "") not in {"c", "j", "k"},
                row["transcript_id"],
            ),
        )
        focus_refs = [row["transcript_id"] for row in refs_sorted[: args.top_refs_per_locus]]
        panel_rows.append(
            {
                "locus_id": lid,
                "chrom": chrom,
                "strand": strand,
                "panel_start": panel_start,
                "panel_end": panel_end,
                "n_needy_refs": len(rows),
                "absent_refs": absent,
                "subopt_refs": subopt,
                "priority_score": priority,
                "classes": join_items(f"{k}:{v}" for k, v in sorted(class_counter.items())),
                "buckets": join_items(f"{k}:{v}" for k, v in sorted(bucket_counter.items())),
                "focus_refs": join_items(focus_refs),
                "debug_bundle": f"{chrom}:{panel_start}-{panel_end}",
                "trace_locus": f"{panel_start}-{panel_end}",
            }
        )

    panel_rows.sort(
        key=lambda row: (
            -int(row["priority_score"]),
            -int(row["n_needy_refs"]),
            str(row["chrom"]),
            int(row["panel_start"]),
        )
    )

    out_tsv = args.out_prefix.with_suffix(".panel.tsv")
    out_md = args.out_prefix.with_suffix(".panel.md")
    out_cmds = args.out_prefix.with_suffix(".commands.sh")
    out_ref_dir = args.out_prefix.parent / f"{args.out_prefix.name}_refs"
    out_ref_dir.mkdir(parents=True, exist_ok=True)

    with out_tsv.open("w", newline="") as handle:
        fieldnames = [
            "locus_id",
            "chrom",
            "strand",
            "panel_start",
            "panel_end",
            "n_needy_refs",
            "absent_refs",
            "subopt_refs",
            "priority_score",
            "classes",
            "buckets",
            "focus_refs",
            "debug_bundle",
            "trace_locus",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in panel_rows[: args.top_loci]:
            writer.writerow(row)

    lines = [
        "# Needy Ref Panel",
        "",
        f"- source ref buckets: `{args.ref_buckets_tsv}`",
        f"- source subopt classes: `{args.subopt_tsv}`",
        f"- total needy loci: `{len(panel_rows)}`",
        f"- top loci written: `{min(args.top_loci, len(panel_rows))}`",
        "",
        "## Global Counts",
        "",
    ]
    for bucket, count in bucket_counts.most_common():
        lines.append(f"- bucket `{bucket}`: `{count}`")
    lines.extend(["", "## Dominant Classes", ""])
    for cls, count in class_counts.most_common():
        lines.append(f"- class `{cls}`: `{count}`")
    subset_lines_by_locus: dict[str, list[str]] = defaultdict(list)
    wanted_tids: set[str] = set()
    for row in panel_rows[: args.top_loci]:
        wanted_tids.update(str(row["focus_refs"]).split(","))
    with args.reference.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n")
            parts = fields.split("\t")
            if len(parts) < 9:
                continue
            tid = parse_transcript_id(parts[8])
            if tid and tid in wanted_tids:
                subset_lines_by_locus[locus_id(tid)].append(fields)

    lines.extend(["", "## Top Loci", ""])
    for idx, row in enumerate(panel_rows[: args.top_loci], start=1):
        subset_path = out_ref_dir / f"{row['locus_id'].replace('.', '_')}.gtf"
        subset_path.write_text("\n".join(subset_lines_by_locus.get(str(row["locus_id"]), [])) + "\n")
        lines.append(
            f"{idx}. `{row['locus_id']}` `{row['chrom']}:{row['panel_start']}-{row['panel_end']}({row['strand']})` "
            f"priority `{row['priority_score']}` needy `{row['n_needy_refs']}` absent `{row['absent_refs']}` "
            f"subopt `{row['subopt_refs']}` classes `{row['classes'] or '-'}` focus `{row['focus_refs']}` subset `{subset_path}`"
        )
    out_md.write_text("\n".join(lines) + "\n")

    cmd_lines = [
        "#!/usr/bin/env bash",
        "# Run from repo root. Compare panel output to the locus reference slice with -R -Q.",
        "set -euo pipefail",
        "",
    ]
    for row in panel_rows[: args.top_loci]:
        subset_path = out_ref_dir / f"{row['locus_id'].replace('.', '_')}.gtf"
        locus_tag = str(row["locus_id"]).replace(".", "_")
        out_base = args.out_prefix.parent / f"{args.out_prefix.name}_run_{locus_tag}"
        cmd_lines.append(f"# {row['locus_id']} refs={row['focus_refs']}")
        cmd_lines.append(
            "RUSTLE_TRACE_LOCUS={trace} \\\n".format(trace=row["trace_locus"])
            + "  target/release/rustle GGO_19.bam -L --chrom {chrom} \\\n".format(
                chrom=row["chrom"]
            )
            + "  -o {out}.gtf \\\n".format(out=out_base)
            + "  --debug-bundle \"{bundle}\" --only-debug-bundle \\\n".format(
                bundle=row["debug_bundle"]
            )
            + "  --trace-reference {tr} \\\n".format(tr=subset_path)
            + "  --trace-output {out}.trace.txt \\\n".format(out=out_base)
            + "  --parity-stage-tsv {out}.parity.tsv \\\n".format(out=out_base)
            + "  --snapshot-jsonl {out}.snap.jsonl\n".format(out=out_base)
        )
        cmd_lines.append(
            "gffcompare -r {ref} -R -Q -o {out}_cmp {out}.gtf\n".format(
                ref=subset_path, out=out_base
            )
        )
        cmd_lines.append("")
    out_cmds.write_text("\n".join(cmd_lines))


if __name__ == "__main__":
    main()
