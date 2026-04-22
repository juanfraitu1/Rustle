#!/usr/bin/env python3
"""Compare GTF attributes on exact-matched transcripts (gffcompare class_code '=').

This answers: for transcripts that rustle got structurally correct, are the
reported attributes (cov/longcov/TPM/FPKM) on the same scale as the reference?

Usage:
  python3 scripts/compare_gtf_attrs_exact_matches.py \
    --reference-gtf GGO_19.gtf \
    --query-gtf out_rustle_GGO_19_20260410.gtf \
    --tmap cmp_rustle_GGO_19_20260410.out_rustle_GGO_19_20260410.gtf.tmap \
    --out-tsv /tmp/gtf_attr_cmp.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, Optional, Tuple


ATTR_RE = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]+)"')


def parse_attrs(attr_field: str) -> Dict[str, str]:
    return {k: v for k, v in ATTR_RE.findall(attr_field)}


def parse_float(x: Optional[str]) -> Optional[float]:
    if x is None:
        return None
    try:
        return float(x)
    except ValueError:
        return None


def iter_transcripts(gtf: Path) -> Iterator[Tuple[str, Dict[str, str]]]:
    with gtf.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            if parts[2] != "transcript":
                continue
            attrs = parse_attrs(parts[8])
            tid = attrs.get("transcript_id")
            if not tid:
                continue
            yield tid, attrs


def load_transcript_attrs(gtf: Path) -> Dict[str, Dict[str, str]]:
    return dict(iter_transcripts(gtf))


def iter_exact_tmap_pairs(tmap: Path) -> Iterator[Tuple[str, str]]:
    with tmap.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cols = {name: i for i, name in enumerate(header)}
        # gffcompare tmap headers: ref_id, class_code, qry_id are typical
        ref_i = cols.get("ref_id")
        qry_i = cols.get("qry_id")
        cc_i = cols.get("class_code")
        if ref_i is None or qry_i is None or cc_i is None:
            raise SystemExit(f"tmap missing expected columns: have {header}")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(ref_i, qry_i, cc_i):
                continue
            if parts[cc_i] != "=":
                continue
            ref_id = parts[ref_i]
            qry_id = parts[qry_i]
            if ref_id and qry_id and ref_id != "-" and qry_id != "-":
                yield ref_id, qry_id


def safe_ratio(a: Optional[float], b: Optional[float]) -> Optional[float]:
    if a is None or b is None:
        return None
    if b == 0:
        return None
    return a / b


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-gtf", type=Path, required=True)
    ap.add_argument("--query-gtf", type=Path, required=True)
    ap.add_argument("--tmap", type=Path, required=True)
    ap.add_argument("--out-tsv", type=Path, required=True)
    args = ap.parse_args()

    ref = load_transcript_attrs(args.reference_gtf)
    qry = load_transcript_attrs(args.query_gtf)

    rows = []
    for ref_id, qry_id in iter_exact_tmap_pairs(args.tmap):
        ra = ref.get(ref_id, {})
        qa = qry.get(qry_id, {})
        rec = {
            "ref_id": ref_id,
            "qry_id": qry_id,
            "ref_cov": ra.get("cov", ""),
            "qry_cov": qa.get("cov", ""),
            "ref_longcov": ra.get("longcov", ""),
            "qry_longcov": qa.get("longcov", ""),
            "ref_tpm": ra.get("TPM", ""),
            "qry_tpm": qa.get("TPM", ""),
            "ref_fpkm": ra.get("FPKM", ""),
            "qry_fpkm": qa.get("FPKM", ""),
        }

        rc = parse_float(ra.get("cov"))
        qc = parse_float(qa.get("cov"))
        rlc = parse_float(ra.get("longcov"))
        qlc = parse_float(qa.get("longcov"))
        rt = parse_float(ra.get("TPM"))
        qt = parse_float(qa.get("TPM"))
        rf = parse_float(ra.get("FPKM"))
        qf = parse_float(qa.get("FPKM"))

        rec["cov_ratio_qry_over_ref"] = safe_ratio(qc, rc) or ""
        rec["longcov_ratio_qry_over_ref"] = safe_ratio(qlc, rlc) or ""
        rec["tpm_ratio_qry_over_ref"] = safe_ratio(qt, rt) or ""
        rec["fpkm_ratio_qry_over_ref"] = safe_ratio(qf, rf) or ""
        rows.append(rec)

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0].keys()) if rows else [])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # Print compact summary to stdout.
    def collect(key: str) -> list[float]:
        out = []
        for r in rows:
            v = r.get(key)
            if v == "" or v is None:
                continue
            try:
                out.append(float(v))
            except ValueError:
                pass
        return out

    def quantiles(xs: list[float]) -> dict[str, float]:
        if not xs:
            return {}
        ys = sorted(xs)
        def q(p: float) -> float:
            i = int(round(p * (len(ys) - 1)))
            return ys[max(0, min(len(ys) - 1, i))]
        return {"n": float(len(ys)), "q10": q(0.10), "q50": q(0.50), "q90": q(0.90)}

    for key in (
        "cov_ratio_qry_over_ref",
        "longcov_ratio_qry_over_ref",
        "tpm_ratio_qry_over_ref",
        "fpkm_ratio_qry_over_ref",
    ):
        qs = quantiles(collect(key))
        if qs:
            print(key, qs)
    print("wrote", args.out_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

