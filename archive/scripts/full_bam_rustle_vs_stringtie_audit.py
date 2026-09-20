#!/usr/bin/env python3
"""
Full-BAM systemic audit: Rustle vs StringTie (same BAM), using gffcompare.

Why this over “needy loci” alone
--------------------------------
Regional mini-BAMs stress-test one gene but distort precision (query GTF scope). A **full BAM**
run measures **global** Sn/Pr and, via `.tmap`, **which gffcompare class_codes dominate** when
Rustle disagrees with StringTie — that is the right scale for **systemic** blockers (bundling,
junction graph, predcluster, boundary snap), not one-off locus quirks.

Recommended pipeline
--------------------
  # 1) Same BAM, de novo (no -G on either tool for apples-to-apples assembly parity)
  stringtie -L -p "$(nproc)" -o bench/full/st.gtf in.bam
  ./target/release/rustle -o bench/full/rs.gtf in.bam -L -p "$(nproc)"

  # 2) StringTie as reference: “how much of StringTie’s transcript set does Rustle recover?”
  gffcompare -r bench/full/st.gtf -o bench/full/cmp_rs bench/full/rs.gtf

  # 3) Optional reverse: Rustle as reference (extra novel calls from StringTie’s perspective)
  gffcompare -r bench/full/rs.gtf -o bench/full/cmp_st bench/full/st.gtf

  # 4) This script (systemic class table + optional trace merge)
  python3 scripts/full_bam_rustle_vs_stringtie_audit.py \\
    --stats bench/full/cmp_rs.stats \\
    --tmap bench/full/cmp_rs.rs.gtf.tmap \\
    --out-tsv bench/full/cmp_rs.systemic_audit.tsv

Trace + rlink.cpp
-----------------
  - Run Rustle with `--trace-reference <ref.gtf>` when you need **NotExtracted** / **Filter**
    reasons for the **same** reference IDs StringTie emitted (use `st.gtf` or a gene subset).
  - Cross-check hotspots against **rlink.cpp** in this repo and the `rlink.cpp:LINE` comments
    in `src/rustle/*.rs` — those are the intended StringTie parity anchors.

class_code (gffcompare): see `gffcompare -h` / manual; `=` exact, `c` contained, `j` novel splice
match, `u` novel unmapped, etc.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


def parse_stats(path: Path) -> dict[str, Any]:
    out: dict[str, Any] = {}
    if not path.is_file():
        return out
    ref_re = re.compile(r"^#\s*Reference mRNAs\s*:\s*(\d+)\s+in\s+(\d+)\s+loci")
    qry_re = re.compile(r"^#\s*Query mRNAs\s*:\s*(\d+)\s+in\s+(\d+)\s+loci")
    tx_re = re.compile(r"^\s*Transcript level:\s*([0-9.]+)\s*\|\s*([0-9.]+)")
    match_re = re.compile(r"^\s*Matching transcripts:\s*(\d+)")
    with path.open() as f:
        for line in f:
            m = ref_re.match(line)
            if m:
                out["reference_mrnas"] = int(m.group(1))
                out["reference_loci"] = int(m.group(2))
                continue
            m = qry_re.match(line)
            if m:
                out["query_mrnas"] = int(m.group(1))
                out["query_loci"] = int(m.group(2))
                continue
            m = tx_re.match(line)
            if m:
                out["transcript_sensitivity_pct"] = float(m.group(1))
                out["transcript_precision_pct"] = float(m.group(2))
                continue
            m = match_re.match(line)
            if m:
                out["matching_transcripts"] = int(m.group(1))
                continue
    return out


def load_tmap_rows(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            rows.append({k: (v or "") for k, v in row.items()})
    return rows


def best_class_for_ref(classes: list[str]) -> str:
    """One status per StringTie ref transcript: '=' wins; else pick a representative miss."""
    if any(c == "=" for c in classes):
        return "="
    if not classes:
        return "?"
    # Prefer informative structural codes for reporting (not perfect ordering).
    order = "ucnjomfpyxesrik"
    ranked = sorted(set(classes), key=lambda c: (0 if c in order else 1, order.index(c) if c in order else 99))
    return ranked[0]


def aggregate_by_ref(rows: list[dict[str, str]]) -> tuple[Counter[str], Counter[str], dict[str, str]]:
    """Returns (raw class counts, per-ref resolved counts, ref_id -> resolved class)."""
    raw = Counter()
    by_ref: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        rid = row.get("ref_id", "").strip()
        cc = row.get("class_code", "").strip()
        if not rid or rid == "-":
            continue
        raw[cc] += 1
        by_ref[rid].append(cc)
    resolved: dict[str, str] = {rid: best_class_for_ref(cls) for rid, cls in by_ref.items()}
    per_ref = Counter(resolved.values())
    return raw, per_ref, resolved


def load_trace_blockers(path: Path | None) -> dict[str, str]:
    """Best-effort: ref_id -> blocker string from rustle trace (# lines or tab)."""
    if path is None or not path.is_file():
        return {}
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 2 and not parts[0].startswith("#"):
            out[parts[0]] = parts[1]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Systemic Rustle vs StringTie audit from gffcompare outputs.")
    ap.add_argument("--stats", type=Path, required=True, help="gffcompare .stats (reference = StringTie)")
    ap.add_argument("--tmap", type=Path, required=True, help="gffcompare .tmap")
    ap.add_argument("--out-tsv", type=Path, required=True, help="Summary TSV (systemic table)")
    ap.add_argument(
        "--trace",
        type=Path,
        default=None,
        help="Optional rustle --trace-reference report to join first column as trace hint",
    )
    args = ap.parse_args()

    st = parse_stats(args.stats)
    rows = load_tmap_rows(args.tmap)
    raw_cnt, per_ref_cnt, resolved = aggregate_by_ref(rows)

    n_ref = len(resolved)
    n_eq = sum(1 for c in resolved.values() if c == "=")
    n_miss = n_ref - n_eq

    trace_map = load_trace_blockers(args.trace)
    # Top systemic classes among **misses** (per-ref, not '=')
    miss_classes = Counter(c for c in resolved.values() if c != "=")

    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w") as f:
        f.write("key\tvalue\n")
        for k, v in sorted(st.items()):
            f.write(f"gffcompare_stats.{k}\t{v}\n")
        f.write("reference_transcripts_in_tmap\t{}\n".format(n_ref))
        f.write("reference_transcripts_exact_match\t{}\n".format(n_eq))
        f.write("reference_transcripts_no_exact_match\t{}\n".format(n_miss))
        f.write("\n# per_reference_class\tcount\n")
        for cls, c in per_ref_cnt.most_common():
            f.write(f"per_ref_class.{cls}\t{c}\n")
        f.write("\n# miss_only_class\tcount\n")
        for cls, c in miss_classes.most_common():
            f.write(f"miss_per_ref.{cls}\t{c}\n")
        f.write("\n# raw_tmap_row_class\tcount\n")
        for cls, c in raw_cnt.most_common():
            f.write(f"raw_row.{cls}\t{c}\n")

    # Optional: sample of worst refs with trace
    sample_path = args.out_tsv.with_suffix(".miss_refs.tsv")
    with sample_path.open("w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["ref_id", "resolved_class", "trace_blocker_if_any"])
        for rid in sorted(resolved.keys()):
            rc = resolved[rid]
            if rc == "=":
                continue
            w.writerow([rid, rc, trace_map.get(rid, "")])

    print("=== Full-BAM systemic audit (StringTie = reference, Rustle = query) ===")
    print(f"  reference transcripts (tmap): {n_ref}")
    print(f"  exact (=) per ref: {n_eq}  miss: {n_miss}")
    if st:
        print(
            f"  gffcompare transcript level: Sn {st.get('transcript_sensitivity_pct', '?')}% | "
            f"Pr {st.get('transcript_precision_pct', '?')}%"
        )
    print("  Top per-ref classes (all):")
    for cls, c in per_ref_cnt.most_common(12):
        print(f"    {cls}\t{c}")
    print("  Top miss-only per-ref classes:")
    for cls, c in miss_classes.most_common(12):
        print(f"    {cls}\t{c}")
    print(f"  Wrote {args.out_tsv}")
    print(f"  Wrote {sample_path} (non-exact refs + optional trace)")
    print()
    print("Next: map dominant miss classes to code areas —")
    print("  j / boundary shifts  -> junction_correction_window (-E), longtrim, guide snap")
    print("  c / o / containment  -> predcluster pairwise, isofrac, dedup_subset")
    print("  u / novel              -> graph connectivity, emit_junction_paths, bundle boundaries")
    print("  Trace: rustle --trace-reference st.gtf --trace-output trace.txt  (same BAM)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
