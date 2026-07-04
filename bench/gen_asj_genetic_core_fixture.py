#!/usr/bin/env python3
"""Golden-fixture generator for the Rust port of bench/asj_genetic_core.py (O3 DELIVERABLE:
the reproducible 54-call genetic core).

Runs the REAL bench.asj_genetic_core logic (load + the row-aligned join + the 3-filter funnel +
the exact csv.DictWriter output) over the shipped bench/asj_calls_verified.tsv +
asj_calls_strandbias.tsv, and captures a SELF-CONTAINED fixture:
  * the two INPUT tables (only the columns the funnel/output need, as rows) so the Rust #[test]
    reconstructs them without the shipped files,
  * the EXACT 14-col output rows (values in csv.DictWriter order),
  * the exact \\r\\n output BYTES (csv.DictWriter reproduction; byte-identical to
    asj_genetic_core.tsv), so the Rust test proves byte fidelity without the shipped output,
  * the funnel counts (n_all/n_transversion/n_nonloc/n_core) + sorted core-gene list, and
  * the retired-flagship rows (PSMD2 & DAXX first high_confidence: transversion=1 but EXCLUDED).

This is a PURE-TSV join -- no BAM, no floats recomputed -- so no pysam/scipy is needed.
Run with /home/juanfra/miniforge3/bin/python (PYTHONHASHSEED=0).
"""
import csv
import io
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import asj_genetic_core as gc  # noqa: E402  (the REAL module under port)

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "asj_genetic_core_fixture.json")

# The 7-tuple row-alignment key (bench/asj_genetic_core.py:36).
KEY = ("gene", "chrom", "anchor", "alleleX", "alleleY", "donor", "acceptor")
# Columns the funnel/output actually read from each input table.
VER_COLS = list(KEY) + ["dPSI", "q", "high_confidence"]
SB_COLS = list(KEY) + ["sor", "sor_pass"]


def main():
    # ---- load EXACTLY like asj_genetic_core.main() (the REAL loader) ----
    ver = gc.load("asj_calls_verified.tsv")   # has high_confidence (transversion)
    sb = gc.load("asj_calls_strandbias.tsv")   # has sor_pass / sor

    # ---- reproduce main()'s two structural asserts ----
    assert len(ver) == len(sb), "call tables differ in length"
    key = lambda r: tuple(r[k] for k in KEY)
    assert all(key(a) == key(b) for a, b in zip(ver, sb)), "call tables are not row-aligned"

    # ---- the 3-filter funnel (verbatim from main()) ----
    cols = ["gene", "chrom", "anchor", "alleleX", "alleleY", "donor", "acceptor", "dPSI", "q",
            "high_confidence", "non_loc_copy_confound_clean", "sor", "sor_pass", "in_genetic_core"]
    rows = []
    for v, s in zip(ver, sb):
        transversion = v["high_confidence"] == "1"
        non_loc = not v["gene"].startswith("LOC")
        sor_pass = s["sor_pass"] == "1"
        in_core = transversion and non_loc and sor_pass
        rows.append({**v, "non_loc_copy_confound_clean": int(non_loc),
                     "sor_pass": int(sor_pass), "sor": s["sor"], "in_genetic_core": int(in_core)})

    n_all = len(rows)
    n_trans = sum(r["high_confidence"] == "1" for r in rows)
    n_nonloc = sum(r["high_confidence"] == "1" and r["non_loc_copy_confound_clean"] for r in rows)
    core = [r for r in rows if r["in_genetic_core"]]
    genes = sorted({r["gene"] for r in core})

    # ---- the exact csv.DictWriter output BYTES (reproduces main()'s writer verbatim) ----
    buf = io.StringIO(newline="")
    w = csv.DictWriter(buf, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow(r)
    out_bytes = buf.getvalue()   # header + rows, each \r\n-terminated (trailing \r\n too)
    header = "\t".join(cols)
    assert out_bytes.startswith(header + "\r\n"), "unexpected header/terminator"

    # cross-check against the shipped file if present (main() writes it)
    shipped_path = os.path.join(HERE, "asj_genetic_core.tsv")
    if os.path.exists(shipped_path):
        with open(shipped_path, "rb") as fh:
            shipped = fh.read()
        assert out_bytes.encode() == shipped, "reproduced bytes != shipped asj_genetic_core.tsv"

    # ---- captured 14-col output rows (str(value) in DictWriter order) ----
    out_rows = [[str(r[c]) for c in cols] for r in rows]

    # ---- retired-flagship rows (main()'s next(... high_confidence=='1' ...) per gene) ----
    flagships = []
    for g in ("PSMD2", "DAXX"):
        hit = next((r for r in rows if r["gene"] == g and r["high_confidence"] == "1"), None)
        assert hit is not None, f"expected a high_confidence {g} row"
        flagships.append({
            "gene": g,
            "high_confidence": hit["high_confidence"],
            "sor": hit["sor"],
            "sor_pass": int(hit["sor_pass"]),
            "in_genetic_core": int(hit["in_genetic_core"]),
        })

    fixture = {
        "source": "bench/asj_genetic_core.py",
        "in_verified_basename": "asj_calls_verified.tsv",
        "in_strandbias_basename": "asj_calls_strandbias.tsv",
        "out_tsv_basename": "asj_genetic_core.tsv",
        "cols": cols,
        "header": header,
        "line_terminator": "\r\n",
        "n_all": n_all,
        "n_transversion": n_trans,
        "n_nonloc": n_nonloc,
        "n_core": len(core),
        "core_genes": genes,
        "flagships": flagships,
        # self-contained INPUT tables (only the needed columns, as rows)
        "ver_rows": [{c: v[c] for c in VER_COLS} for v in ver],
        "sb_rows": [{c: s[c] for c in SB_COLS} for s in sb],
        # exact output rows (14 cols) + full \r\n bytes
        "out_rows": out_rows,
        "out_tsv_bytes": out_bytes,
    }
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=1)
        fh.write("\n")

    print(f"[wrote {OUT}]")
    print(f"  funnel: {n_all} -> {n_trans} transversion -> {n_nonloc} non-LOC -> {len(core)} core "
          f"in {len(genes)} genes")
    print(f"  in_genetic_core sum = {sum(int(r['in_genetic_core']) for r in rows)}")
    print(f"  core genes: {', '.join(genes)}")
    for f in flagships:
        print(f"  [flagship] {f['gene']}: hc={f['high_confidence']} sor={f['sor']} "
              f"sor_pass={f['sor_pass']} in_core={f['in_genetic_core']}")
    print(f"  out bytes = {len(out_bytes)} ({out_bytes.count(chr(13)+chr(10))} CRLF lines)")


if __name__ == "__main__":
    main()
