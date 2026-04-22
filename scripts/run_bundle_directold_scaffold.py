#!/usr/bin/env python3

import argparse
import csv
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the focused rlink/StringTie-CPP parity scaffold for one targeted bundle: "
            "generate the rlink bundle output, ensure the current BNODE "
            "comparison artifacts exist, run the default soft-filter+augment parity "
            "probe, and emit chain-delta/blocker summaries."
        )
    )
    parser.add_argument("--bundle-id", type=int, required=True)
    parser.add_argument("--input-bam", required=True)
    parser.add_argument(
        "--assembler",
        default="target/debug/isoseq-assembler",
        help="Path to the built isoseq-assembler binary.",
    )
    parser.add_argument(
        "--bundle-prefix",
        help=(
            "Base prefix for bundle-local artifacts. Defaults to "
            "GGO_19_bundle<BUNDLE_ID>."
        ),
    )
    parser.add_argument(
        "--stringtie-gtf",
        default="GGO_19_stringtie.gtf",
        help="Backward StringTie GTF used by bundle_chain_blockers.py.",
    )
    parser.add_argument(
        "--supported-gfa-dir",
        default="chr19_stringtie_supported_gfas",
        help="Directory of backward supported GFAs used by bundle_chain_blockers.py.",
    )
    parser.add_argument(
        "--bundle-gap",
        type=int,
        default=50,
        help="Gap threshold for strand-split bundle creation.",
    )
    parser.add_argument(
        "--dedup-tolerance",
        type=int,
        default=5,
        help="Deduplication tolerance used by the focused modes.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-run steps even when the expected outputs already exist.",
    )
    return parser.parse_args()


def run(cmd: List[str]) -> None:
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def ensure_outputs(outputs: Iterable[Path], cmd: List[str], force: bool) -> None:
    outputs = list(outputs)
    if not force and all(path.exists() for path in outputs):
        print(f"= reuse {outputs[0].parent / outputs[0].name}")
        return
    run(cmd)


def bundle_meta_from_transfrag(path: Path, bundle_id: int) -> Tuple[str, int, int]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if int(row["bundle_id"]) == bundle_id:
                return row["chrom"], int(row["start"]), int(row["end"])
    raise SystemExit(f"Bundle {bundle_id} not found in {path}")


def main() -> None:
    args = parse_args()

    prefix = args.bundle_prefix or f"GGO_19_bundle{args.bundle_id}"
    assembler = args.assembler
    bundle_id = str(args.bundle_id)

    bnode_diag_stem = f"{prefix}_bnode_diag"
    transfrag_stem = f"{prefix}_transfrag"
    familytrace_stem = f"{prefix}_familytrace"
    bnode_scoped_gtf = f"{prefix}_bnode_scoped.gtf"
    parityprobe_stem = f"{prefix}_bnode_softaugment"
    parityprobe_gtf = f"{parityprobe_stem}.gtf"
    rlink_main_stem = f"{prefix}_rlink_main"
    rlink_final_gtf = f"{rlink_main_stem}.gtf"
    rlink_transcripts_tsv = f"{rlink_main_stem}.rlink.transcripts.tsv"
    delta_stem = f"{prefix}_rlink_chain_delta"
    blockers_stem = f"{prefix}_rlink_chain_blockers"
    parity_delta_stem = f"{prefix}_softaugment_chain_delta"
    parity_blockers_stem = f"{prefix}_softaugment_chain_blockers"

    ensure_outputs(
        [
            Path(f"{bnode_diag_stem}.bnodes.tsv"),
            Path(f"{bnode_diag_stem}.bundles.tsv"),
        ],
        [
            assembler,
            "--input",
            args.input_bam,
            "--output",
            f"{bnode_diag_stem}.gtf",
            "--bundle-bnode-diagnostic",
            "--bundle-bnode-target-id",
            bundle_id,
            "--bundle-bnode-gap",
            str(args.bundle_gap),
            "--bundle-bnode-output-stem",
            bnode_diag_stem,
        ],
        args.force,
    )

    ensure_outputs(
        [
            Path(f"{transfrag_stem}.transfrags.tsv"),
            Path(f"{transfrag_stem}.bundles.tsv"),
        ],
        [
            assembler,
            "--input",
            args.input_bam,
            "--output",
            f"{transfrag_stem}.gtf",
            "--bundle-transfrag-diagnostic",
            "--bundle-bnode-target-id",
            bundle_id,
            "--bundle-bnode-gap",
            str(args.bundle_gap),
            "--bundle-bnode-output-stem",
            transfrag_stem,
        ],
        args.force,
    )

    ensure_outputs(
        [
            Path(f"{familytrace_stem}.family_members.tsv"),
            Path(f"{familytrace_stem}.bundles.tsv"),
        ],
        [
            assembler,
            "--input",
            args.input_bam,
            "--output",
            f"{familytrace_stem}.gtf",
            "--bundle-bnode-family-trace",
            "--bundle-bnode-target-id",
            bundle_id,
            "--bundle-bnode-gap",
            str(args.bundle_gap),
            "--bundle-bnode-output-stem",
            familytrace_stem,
        ],
        args.force,
    )

    ensure_outputs(
        [Path(bnode_scoped_gtf)],
        [
            assembler,
            "--input",
            args.input_bam,
            "--output",
            bnode_scoped_gtf,
            "--bundle-bnode-assembly",
            "--bundle-bnode-target-id",
            bundle_id,
            "--bundle-bnode-gap",
            str(args.bundle_gap),
            "--dedup",
            "--dedup-tolerance",
            str(args.dedup_tolerance),
        ],
        args.force,
    )

    ensure_outputs(
        [
            Path(rlink_final_gtf),
            Path(f"{rlink_main_stem}.rlink.summary.tsv"),
            Path(f"{rlink_main_stem}.rlink.nodes.tsv"),
            Path(rlink_transcripts_tsv),
            Path(f"{rlink_main_stem}.rlink.gtf"),
        ],
        [
            assembler,
            "--input",
            args.input_bam,
            "--output",
            rlink_final_gtf,
            "--bundle-rlink-assembly",
            "--bundle-bnode-target-id",
            bundle_id,
            "--bundle-bnode-gap",
            str(args.bundle_gap),
            "--bundle-bnode-output-stem",
            rlink_main_stem,
            "--dedup",
            "--dedup-tolerance",
            str(args.dedup_tolerance),
        ],
        args.force,
    )

    ensure_outputs(
        [Path(parityprobe_gtf)],
        [
            assembler,
            "--input",
            args.input_bam,
            "--output",
            parityprobe_gtf,
            "--bundle-bnode-assembly",
            "--bundle-bnode-target-id",
            bundle_id,
            "--bundle-bnode-gap",
            str(args.bundle_gap),
            "--bundle-bnode-output-stem",
            parityprobe_stem,
            "--bundle-rlink-guide-soft-filter",
            "--bundle-rlink-guide-augment",
        ],
        args.force,
    )

    ensure_outputs(
        [
            Path(f"{delta_stem}.summary.tsv"),
            Path(f"{delta_stem}.deltas.tsv"),
        ],
        [
            sys.executable,
            "scripts/bundle_chain_delta.py",
            "--bundle-id",
            bundle_id,
            "--bnodes",
            f"{bnode_diag_stem}.bnodes.tsv",
            "--oldpath",
            rlink_transcripts_tsv,
            "--gtf",
            bnode_scoped_gtf,
            "--transfrags",
            f"{transfrag_stem}.transfrags.tsv",
            "--families",
            f"{familytrace_stem}.family_members.tsv",
            "--output-stem",
            delta_stem,
        ],
        args.force,
    )

    chrom, start, end = bundle_meta_from_transfrag(
        Path(f"{transfrag_stem}.bundles.tsv"), args.bundle_id
    )

    ensure_outputs(
        [
            Path(f"{blockers_stem}.summary.tsv"),
            Path(f"{blockers_stem}.blockers.tsv"),
        ],
        [
            sys.executable,
            "scripts/bundle_chain_blockers.py",
            "--delta-tsv",
            f"{delta_stem}.deltas.tsv",
            "--stringtie-gtf",
            args.stringtie_gtf,
            "--supported-gfa-dir",
            args.supported_gfa_dir,
            "--chrom",
            chrom,
            "--start",
            str(start),
            "--end",
            str(end),
            "--output-stem",
            blockers_stem,
        ],
        args.force,
    )

    ensure_outputs(
        [
            Path(f"{parity_delta_stem}.summary.tsv"),
            Path(f"{parity_delta_stem}.deltas.tsv"),
        ],
        [
            sys.executable,
            "scripts/bundle_chain_delta.py",
            "--bundle-id",
            bundle_id,
            "--bnodes",
            f"{bnode_diag_stem}.bnodes.tsv",
            "--oldpath",
            rlink_transcripts_tsv,
            "--gtf",
            parityprobe_gtf,
            "--transfrags",
            f"{transfrag_stem}.transfrags.tsv",
            "--families",
            f"{familytrace_stem}.family_members.tsv",
            "--output-stem",
            parity_delta_stem,
        ],
        args.force,
    )

    ensure_outputs(
        [
            Path(f"{parity_blockers_stem}.summary.tsv"),
            Path(f"{parity_blockers_stem}.blockers.tsv"),
        ],
        [
            sys.executable,
            "scripts/bundle_chain_blockers.py",
            "--delta-tsv",
            f"{parity_delta_stem}.deltas.tsv",
            "--stringtie-gtf",
            args.stringtie_gtf,
            "--supported-gfa-dir",
            args.supported_gfa_dir,
            "--chrom",
            chrom,
            "--start",
            str(start),
            "--end",
            str(end),
            "--output-stem",
            parity_blockers_stem,
        ],
        args.force,
    )

    print()
    print("Scaffold outputs:")
    print(f"  final rlink:           {rlink_final_gtf}")
    print(f"  raw rlink:             {rlink_main_stem}.rlink.gtf")
    print(f"  baseline chain delta:  {delta_stem}.summary.tsv")
    print(f"  baseline blockers:     {blockers_stem}.summary.tsv")
    print(f"  parity probe gtf:      {parityprobe_gtf}")
    print(f"  parity probe delta:    {parity_delta_stem}.summary.tsv")
    print(f"  parity probe blockers: {parity_blockers_stem}.summary.tsv")


if __name__ == "__main__":
    main()
