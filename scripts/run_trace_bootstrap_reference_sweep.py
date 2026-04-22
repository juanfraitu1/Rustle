#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class GeneSpec:
    gene_id: str
    transcript_count: int


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--binary", type=Path, default=Path("target/debug/isoseq-assembler"))
    ap.add_argument("--binary-extra-arg", action="append", default=[])
    ap.add_argument("--input-bam", type=Path, default=Path("GGO_19.bam"))
    ap.add_argument("--reference", type=Path, default=Path("GGO_19.gtf"))
    ap.add_argument("--deep-trace-log", type=Path, default=Path("trace_GGO_19_deep.log"))
    ap.add_argument("--full-trace-log", type=Path, default=Path("trace_full.log"))
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--limit", type=int)
    ap.add_argument("--min-isoforms", type=int, default=1)
    ap.add_argument("--gene-id", action="append", default=[])
    ap.add_argument("--skip-existing", action="store_true")
    ap.add_argument("--locus-slack-bp", type=int, default=10000)
    ap.add_argument(
        "--emitted-tsv-suffix",
        default="backward_trace_emitted_isoforms.tsv",
        help="Sibling suffix under each output stem used for reference comparison",
    )
    ap.add_argument(
        "--compare-suffix",
        default="_compare",
        help="Suffix appended to each output stem for comparison outputs",
    )
    ap.add_argument(
        "--blocker-suffix",
        default="_blockers",
        help="Suffix appended to each output stem for blocker outputs",
    )
    ap.add_argument("--label", default="Trace Bootstrap Sweep")
    ap.add_argument("--compare-script", type=Path, default=Path("scripts/compare_trace_bootstrap_to_reference.py"))
    ap.add_argument("--blocker-script", type=Path, default=Path("scripts/summarize_trace_bootstrap_blockers.py"))
    ap.add_argument("--skip-blockers", action="store_true")
    return ap.parse_args()


def parse_attrs(attr_field: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for chunk in attr_field.strip().split(";"):
        chunk = chunk.strip()
        if not chunk or " " not in chunk:
            continue
        key, value = chunk.split(" ", 1)
        attrs[key] = value.strip().strip('"')
    return attrs


def parse_reference_genes(reference_path: Path) -> list[GeneSpec]:
    counts: Counter[str] = Counter()
    with reference_path.open() as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = parse_attrs(fields[8])
            gene_id = attrs.get("gene_id")
            if gene_id:
                counts[gene_id] += 1
    return [
        GeneSpec(gene_id=gene_id, transcript_count=count)
        for gene_id, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    ]


def run_cmd(cmd: list[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log_fh:
        proc = subprocess.run(cmd, stdout=log_fh, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed ({proc.returncode}): {' '.join(cmd)}")


def load_single_row(path: Path) -> dict[str, str]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            return row
    raise RuntimeError(f"no rows found in {path}")


def run_gene(args: argparse.Namespace, spec: GeneSpec) -> dict[str, object]:
    safe_id = spec.gene_id.replace("/", "_").replace(".", "_")
    stem = args.out_dir / safe_id
    output_gtf = Path(f"{stem}.gtf")
    emitted_tsv = args.out_dir / f"{safe_id}.{args.emitted_tsv_suffix}"
    compare_prefix = args.out_dir / f"{safe_id}{args.compare_suffix}"
    compare_summary = compare_prefix.with_suffix(".summary.tsv")
    compare_reference = compare_prefix.with_suffix(".reference.tsv")
    compare_emitted = compare_prefix.with_suffix(".emitted.tsv")
    blocker_prefix = args.out_dir / f"{safe_id}{args.blocker_suffix}"
    blocker_tsv = blocker_prefix.with_suffix(".tsv")
    run_log = args.out_dir / f"{safe_id}.run.log"
    compare_log = args.out_dir / f"{safe_id}.compare.log"
    blocker_log = args.out_dir / f"{safe_id}.blockers.log"

    if not (args.skip_existing and compare_summary.exists() and compare_reference.exists() and compare_emitted.exists()):
        cmd = [
            str(args.binary),
            "--input", str(args.input_bam),
            "--output", str(output_gtf),
            "--reference", str(args.reference),
            "--isoseq-backward-trace",
            "--isoseq-backward-gene-id", spec.gene_id,
            "--isoseq-backward-locus-slack-bp", str(args.locus_slack_bp),
            "--isoseq-backward-deep-trace-log", str(args.deep_trace_log),
            "--isoseq-backward-full-trace-log", str(args.full_trace_log),
        ]
        cmd.extend(args.binary_extra_arg)
        run_cmd(cmd, run_log)

        if not emitted_tsv.exists():
            raise RuntimeError(f"missing emitted TSV for {spec.gene_id}: {emitted_tsv}")

        compare_cmd = [
            sys.executable,
            str(args.compare_script),
            "--reference", str(args.reference),
            "--emitted-tsv", str(emitted_tsv),
            "--out-prefix", str(compare_prefix),
        ]
        run_cmd(compare_cmd, compare_log)

        if not args.skip_blockers:
            blocker_cmd = [
                sys.executable,
                str(args.blocker_script),
                "--reference-tsv", str(compare_reference),
                "--emitted-tsv", str(compare_emitted),
                "--out-prefix", str(blocker_prefix),
            ]
            run_cmd(blocker_cmd, blocker_log)

    row = load_single_row(compare_summary)
    return {
        "gene_id": spec.gene_id,
        "transcript_count": spec.transcript_count,
        "n_emitted": int(row["n_emitted"]),
        "n_reference": int(row["n_reference"]),
        "ref_exact_ref_match": int(row["ref_exact_ref_match"]),
        "ref_same_chain_boundary_mismatch": int(row["ref_same_chain_boundary_mismatch"]),
        "ref_partial_chain_overlap": int(row["ref_partial_chain_overlap"]),
        "ref_unmatched": int(row["ref_unmatched"]),
        "summary_tsv": str(compare_summary),
        "reference_tsv": str(compare_reference),
        "blocker_tsv": str(blocker_tsv),
        "run_log": str(run_log),
    }


def choose_genes(all_genes: list[GeneSpec], wanted: Iterable[str], min_isoforms: int, limit: int | None) -> list[GeneSpec]:
    wanted_set = {gene_id for gene_id in wanted if gene_id}
    chosen = [spec for spec in all_genes if spec.transcript_count >= min_isoforms and (not wanted_set or spec.gene_id in wanted_set)]
    if limit is not None:
        chosen = chosen[:limit]
    return chosen


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_md(path: Path, completed: list[dict[str, object]], failed: list[dict[str, object]], label: str) -> None:
    ref_total = sum(int(row["n_reference"]) for row in completed)
    exact_total = sum(int(row["ref_exact_ref_match"]) for row in completed)
    same_total = sum(int(row["ref_same_chain_boundary_mismatch"]) for row in completed)
    partial_total = sum(int(row["ref_partial_chain_overlap"]) for row in completed)
    unmatched_total = sum(int(row["ref_unmatched"]) for row in completed)
    emitted_total = sum(int(row["n_emitted"]) for row in completed)

    lines = [
        f"# {label}",
        "",
        f"- loci completed: `{len(completed)}`",
        f"- loci failed: `{len(failed)}`",
        f"- emitted isoforms total: `{emitted_total}`",
        f"- reference isoforms total: `{ref_total}`",
        f"- reference exact matches recovered: `{exact_total}`",
        f"- reference same-chain boundary mismatches: `{same_total}`",
        f"- reference partial-chain overlaps: `{partial_total}`",
        f"- reference unmatched: `{unmatched_total}`",
        "",
    ]

    residual = [
        row for row in completed
        if int(row["ref_same_chain_boundary_mismatch"]) > 0
        or int(row["ref_partial_chain_overlap"]) > 0
        or int(row["ref_unmatched"]) > 0
    ]
    if residual:
        lines.extend(["## Residual Loci", ""])
        for row in sorted(
            residual,
            key=lambda item: (
                -(int(item["ref_unmatched"]) + int(item["ref_partial_chain_overlap"]) + int(item["ref_same_chain_boundary_mismatch"])),
                -int(item["n_reference"]),
                str(item["gene_id"]),
            ),
        ):
            lines.append(
                "- `{}`: exact `{}/{}`, same-chain `{}`, partial `{}`, unmatched `{}`, emitted `{}`".format(
                    row["gene_id"],
                    row["ref_exact_ref_match"],
                    row["n_reference"],
                    row["ref_same_chain_boundary_mismatch"],
                    row["ref_partial_chain_overlap"],
                    row["ref_unmatched"],
                    row["n_emitted"],
                )
            )
        lines.append("")

    if failed:
        lines.extend(["## Failed Loci", ""])
        for row in failed:
            lines.append(f"- `{row['gene_id']}`: `{row['error']}`")
        lines.append("")

    path.write_text("\n".join(lines) + "\n")


def concatenate_blockers(completed: list[dict[str, object]], out_prefix: Path) -> None:
    rows: list[dict[str, str]] = []
    for row in completed:
        blocker_path = Path(str(row["blocker_tsv"]))
        if not blocker_path.exists() or blocker_path.stat().st_size == 0:
            continue
        with blocker_path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows.extend(reader)
    rows.sort(key=lambda item: (item["gene_id"], item["transcript_id"]))
    write_tsv(out_prefix.with_suffix(".tsv"), rows)

    lines = ["# Trace Bootstrap Sweep Blockers", ""]
    current_gene = None
    for row in rows:
        gene_id = row["gene_id"]
        if gene_id != current_gene:
            if current_gene is not None:
                lines.append("")
            lines.append(f"## {gene_id}")
            lines.append("")
            current_gene = gene_id
        lines.append(
            "- `{}`: `{}` via `{}` against `{}`".format(
                row["transcript_id"],
                row["blocker_class"],
                row["ref_status"],
                row["best_emitted"] or ".",
            )
        )
    out_prefix.with_suffix(".md").write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    all_genes = parse_reference_genes(args.reference)
    selected = choose_genes(all_genes, args.gene_id, args.min_isoforms, args.limit)
    if not selected:
        raise SystemExit("no genes selected")

    completed: list[dict[str, object]] = []
    failed: list[dict[str, object]] = []

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {pool.submit(run_gene, args, spec): spec for spec in selected}
        for future in as_completed(futures):
            spec = futures[future]
            try:
                completed.append(future.result())
                print(
                    f"[ok] {spec.gene_id} ref={completed[-1]['n_reference']} exact={completed[-1]['ref_exact_ref_match']} "
                    f"same={completed[-1]['ref_same_chain_boundary_mismatch']} "
                    f"partial={completed[-1]['ref_partial_chain_overlap']} unmatched={completed[-1]['ref_unmatched']}",
                    flush=True,
                )
            except Exception as exc:  # noqa: BLE001
                failed.append({"gene_id": spec.gene_id, "transcript_count": spec.transcript_count, "error": str(exc)})
                print(f"[fail] {spec.gene_id}: {exc}", file=sys.stderr, flush=True)

    completed.sort(key=lambda row: (-int(row["transcript_count"]), str(row["gene_id"])))
    failed.sort(key=lambda row: (-int(row["transcript_count"]), str(row["gene_id"])))

    completed_path = args.out_dir / "sweep.summary.tsv"
    failed_path = args.out_dir / "sweep.failed.tsv"
    write_tsv(completed_path, completed)
    write_tsv(failed_path, failed)
    write_md(args.out_dir / "sweep.summary.md", completed, failed, args.label)
    concatenate_blockers(completed, args.out_dir / "sweep.blockers")


if __name__ == "__main__":
    main()
