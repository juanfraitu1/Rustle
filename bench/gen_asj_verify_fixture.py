#!/usr/bin/env python3
"""Golden-fixture generator for the Rust port of bench/asj_verify.py (O3 confound control).

Imports the REAL bench.asj_verify (frac_mq0_at) + pysam, replays main()'s per-anchor +
per-row logic over the shipped bench/asj_calls.tsv against GGO.bam, and emits raw f64 bits
AND the exact Python str(float) for every frac_mq0 so the Rust #[test] can assert:
  * frac_mq0 BIT-EXACT (per unique anchor) + its written str() reproduced,
  * dist (int) exact, high_confidence exact, and the full verified-TSV line text exact.

Per-anchor n/z (and the 600-cap flag) are recomputed here for COVERAGE metadata only; the
parity value itself comes from the REAL asj_verify.frac_mq0_at. Run with
/home/juanfra/miniforge3/bin/python (PYTHONHASHSEED=0; pysam).
"""
import json
import os
import struct
import sys
import time

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import asj_verify as av  # noqa: E402

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "asj_verify_fixture.json")


def bits(x):
    """Little-endian IEEE-754 double bit pattern as lowercase hex (16 chars)."""
    return struct.pack("<d", x).hex()


def counts_at(bam, chrom, pos):
    """COVERAGE-ONLY replica of frac_mq0_at's loop returning (n, z, capped). The parity
    frac value comes from the real av.frac_mq0_at; this only reports cap coverage."""
    n = z = 0
    capped = False
    for aln in bam.fetch(chrom, pos, pos + 1):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        n += 1
        if aln.mapping_quality == 0:
            z += 1
        if n >= 600:
            capped = True
            break
    return n, z, capped


def main():
    t0 = time.time()
    bam = pysam.AlignmentFile(av.BAM, "rb")

    # ---- read calls EXACTLY like asj_verify.main() (no blank-line skip) ----
    with open(av.CALLS) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(hdr, line.rstrip("\n").split("\t"))) for line in fh]

    # ---- per unique anchor: frac_mq0 via the REAL function (first-occurrence order) ----
    order = []
    anchors = {}
    for r in rows:
        key = (r["chrom"], int(r["anchor"]))
        if key not in anchors:
            frac = av.frac_mq0_at(bam, key[0], key[1])  # the REAL parity value
            n, z, capped = counts_at(bam, key[0], key[1])  # coverage metadata only
            anchors[key] = frac
            order.append({
                "chrom": key[0], "anchor": key[1],
                "frac_mq0_bits": bits(frac), "frac_mq0_str": str(frac),
                "n": n, "z": z, "capped": bool(capped),
            })

    # ---- per row: dist, frac_mq0, high_confidence, exact written line ----
    cols = ["gene", "chrom", "anchor", "alleleX", "alleleY", "donor", "acceptor",
            "psiX", "psiY", "dPSI", "q", "anchor_type", "dist", "frac_mq0"]
    calls = []
    n_hc = 0
    for i, r in enumerate(rows):
        anc = int(r["anchor"]); d = int(r["donor"]); a = int(r["acceptor"])
        r["dist"] = min(abs(anc - d), abs(anc - a))
        r["frac_mq0"] = anchors[(r["chrom"], anc)]
        h = int(r["anchor_type"] == "transversion" and r["frac_mq0"] < av.MQ0_THR)
        n_hc += h
        line = "\t".join(str(r[c]) for c in cols) + f"\t{h}"
        calls.append({
            "idx": i, "gene": r["gene"], "chrom": r["chrom"],
            "anchor": anc, "donor": d, "acceptor": a,
            "anchor_type": r["anchor_type"],
            "dist": r["dist"],
            "frac_mq0_bits": bits(r["frac_mq0"]), "frac_mq0_str": str(r["frac_mq0"]),
            "high_confidence": h,
            "verified_line": line,
        })

    verified_header = "\t".join(cols) + "\thigh_confidence"
    cap_hits = sum(1 for a in order if a["capped"])
    max_n = max((a["n"] for a in order), default=0)

    fixture = {
        "source": "bench/asj_verify.py",
        "bam_path": av.BAM,
        "in_tsv_basename": os.path.basename(av.CALLS),
        "out_tsv_basename": os.path.basename(av.OUT),
        "mq0_thr": av.MQ0_THR,
        "prox": av.PROX,
        "n_calls": len(calls),
        "n_unique_anchors": len(order),
        "cap_600_hit_count": cap_hits,
        "max_reads_n": max_n,
        "n_high_confidence": n_hc,
        "verified_header": verified_header,
        "anchors": order,
        "calls": calls,
    }
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=1)
        fh.write("\n")

    dt = time.time() - t0
    print(f"[wrote {OUT}]  {len(calls)} calls / {len(order)} unique anchors in {dt:.1f}s")
    n_tv = sum(1 for c in calls if c["anchor_type"] == "transversion")
    n_zero = sum(1 for a in order if a["frac_mq0_bits"] == bits(0.0))
    print(f"  transversion={n_tv} non-transversion={len(calls) - n_tv} high_confidence={n_hc}")
    print(f"  anchors: zero-frac={n_zero} pos-frac={len(order) - n_zero} "
          f"max_reads_n={max_n} cap_600_hits={cap_hits}")
    if cap_hits == 0:
        print(f"  NOTE: no anchor exceeds 600 reads (max_n={max_n}); the 600-cap is exercised "
              f"structurally but never triggers on this BAM.")


if __name__ == "__main__":
    main()
