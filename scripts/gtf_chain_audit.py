#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Tx:
    gene_id: str
    tx_id: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]]  # 1-based inclusive

    def introns(self) -> list[tuple[int, int]]:
        ex = sorted(self.exons)
        out: list[tuple[int, int]] = []
        for (s1, e1), (s2, e2) in zip(ex, ex[1:]):
            # GTF is 1-based inclusive; intron is (e1, s2) boundary pair.
            out.append((e1, s2))
        return out

    def span(self) -> tuple[int, int]:
        ex = sorted(self.exons)
        return ex[0][0], ex[-1][1]


def parse_attrs(attr_field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in attr_field.split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, value = part.split(" ", 1)
        out[key] = value.strip().strip('"')
    return out


def load_gtf(path: Path) -> list[Tx]:
    tx: dict[str, Tx] = {}
    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, s, e, _, strand, _, attrs_raw = fields
            if feature != "exon":
                continue
            start = int(s)
            end = int(e)
            attrs = parse_attrs(attrs_raw)
            gene_id = attrs.get("gene_id")
            tx_id = attrs.get("transcript_id")
            if not gene_id or not tx_id:
                continue
            rec = tx.get(tx_id)
            if rec is None:
                rec = Tx(gene_id=gene_id, tx_id=tx_id, chrom=chrom, strand=strand, exons=[])
                tx[tx_id] = rec
            rec.exons.append((start, end))
    return list(tx.values())


def overlaps(tx: Tx, chrom: str, start: int, end: int) -> bool:
    if tx.chrom != chrom:
        return False
    ts, te = tx.span()
    return not (te < start or ts > end)


def main() -> int:
    ap = argparse.ArgumentParser(description="Audit intron-chain signatures from a GTF.")
    ap.add_argument("--gtf", type=Path, required=True)
    ap.add_argument("--gene-id", default="", help="Filter: exact gene_id")
    ap.add_argument(
        "--region",
        default="",
        help="Filter: REGION like 'chr:start-end' (1-based inclusive).",
    )
    ap.add_argument("--max-exons-print", type=int, default=12)
    args = ap.parse_args()

    region_chrom = ""
    region_start = 0
    region_end = 0
    if args.region:
        chrom, rest = args.region.split(":", 1)
        a, b = rest.split("-", 1)
        region_chrom = chrom
        region_start = int(a.replace(",", ""))
        region_end = int(b.replace(",", ""))

    txx = load_gtf(args.gtf)
    txx.sort(key=lambda t: (t.chrom, t.span()[0], t.span()[1], t.gene_id, t.tx_id))

    for t in txx:
        if args.gene_id and t.gene_id != args.gene_id:
            continue
        if args.region and not overlaps(t, region_chrom, region_start, region_end):
            continue
        ex = sorted(t.exons)
        intr = t.introns()
        ts, te = t.span()
        ex_preview = ex[: args.max_exons_print]
        ex_str = ",".join(f"{a}-{b}" for a, b in ex_preview)
        if len(ex) > len(ex_preview):
            ex_str += f",...(+{len(ex)-len(ex_preview)} more)"
        intr_str = ",".join(f"{a}-{b}" for a, b in intr[: args.max_exons_print])
        if len(intr) > args.max_exons_print:
            intr_str += f",...(+{len(intr)-args.max_exons_print} more)"
        print(
            "\t".join(
                [
                    t.gene_id,
                    t.tx_id,
                    t.chrom,
                    t.strand,
                    f"{ts}-{te}",
                    f"exons={len(ex)}",
                    f"introns={len(intr)}",
                    f"exon_preview={ex_str}",
                    f"intron_preview={intr_str}",
                ]
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

