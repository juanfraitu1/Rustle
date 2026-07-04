#!/usr/bin/env python3
"""Golden-fixture generator for the Rust port of bench/asj_strand_bias.py (O3 SOR filter).

Runs the REAL bench.asj_strand_bias functions (strand_table + sor + fisher_strand_p) over
the shipped bench/asj_calls.tsv against GGO_mm.bam, replicating main()'s per-row logic
EXACTLY (donor/acceptor normalise, chrom-not-in-refs -> NaN, zero-used tally). Emits raw
f64 bits for sor/fisher so the Rust #[test] can assert byte-parity of the SOR gate and
tolerance-parity of the (non-gate) fisher cross-check.

Also appends a few SYNTHETIC edge rows -- a fake chrom (NaN branch) and a donor>acceptor
row (normalise/swap branch) -- computed with the SAME real functions, because the real
475 calls happen to have neither a missing chrom nor a swap.

Run with /home/juanfra/miniforge3/bin/python (PYTHONHASHSEED=0; pysam + scipy).
"""
import json
import math
import os
import struct
import sys
import time

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import asj_strand_bias as sb  # noqa: E402
from scipy.stats import fisher_exact  # noqa: E402

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "asj_strand_bias_fixture.json")


def bits(x):
    """Little-endian IEEE-754 double bit pattern as lowercase hex (16 chars)."""
    return struct.pack("<d", x).hex()


# ---- Python replica of the Rust allele_specific_junctions::fisher_exact_2x2 --------------
# (Lanczos lgamma + the "pk <= p_obs*(1+1e-7)" two-sided sum) to MEASURE the divergence of
# the reused Rust routine from scipy, so the Rust test's tolerance is empirically justified.
_LG_C = [
    0.9999999999998099, 676.5203681218851, -1259.1392167224028, 771.3234287776531,
    -176.6150291621406, 12.507343278686905, -0.1385710952657201, 9.984369578019572e-6,
    1.5056327351493116e-7,
]


def _lgamma_lanczos(x):
    if x < 0.5:
        return math.log(math.pi) - math.log(math.sin(math.pi * x)) - _lgamma_lanczos(1.0 - x)
    x -= 1.0
    t = x + 7.0 + 0.5
    a = _LG_C[0]
    for i in range(1, 9):
        a += _LG_C[i] / (x + i)
    return 0.5 * math.log(2.0 * math.pi) + (x + 0.5) * math.log(t) - t + math.log(a)


def _lchoose(n, k):
    if k > n:
        return float("-inf")
    return _lgamma_lanczos(n + 1.0) - _lgamma_lanczos(k + 1.0) - _lgamma_lanczos(n - k + 1.0)


def fisher_rust_replica(a, b, c, d):
    r1, r2, col1, n = a + b, c + d, a + c, a + b + c + d
    if n == 0:
        return 1.0
    lcn = _lchoose(n, col1)

    def pmf(k):
        return math.exp(_lchoose(r1, k) + _lchoose(r2, col1 - k) - lcn)

    p_obs = pmf(a)
    lo = max(0, col1 - r2)
    hi = min(col1, r1)
    total = 0.0
    for k in range(lo, hi + 1):
        pk = pmf(k)
        if pk <= p_obs * (1.0 + 1e-7):
            total += pk
    return min(total, 1.0)


def process(bam, refs, chrom, donor, acceptor):
    """Replicate main()'s per-row computation (normalise + branches). Returns a dict."""
    if donor > acceptor:
        donor, acceptor = acceptor, donor
    if chrom not in refs:
        rf = rr = af = ar = 0
        s = float("nan")
        fp = float("nan")
    else:
        rf, rr, af, ar = sb.strand_table(bam, chrom, donor, acceptor)
        s = sb.sor(rf, rr, af, ar)
        fp = sb.fisher_strand_p(rf, rr, af, ar)
    sor_pass = (not math.isnan(s)) and s < sb.SOR_KEEP
    return {
        "chrom": chrom, "donor": donor, "acceptor": acceptor,
        "chrom_in_refs": chrom in refs,
        "R_fwd": rf, "R_rev": rr, "A_fwd": af, "A_rev": ar,
        "sor_bits": bits(s), "sor_is_nan": math.isnan(s),
        "fisher_bits": bits(fp), "fisher_is_nan": math.isnan(fp),
        "sor_pass": bool(sor_pass),
    }


def main():
    t0 = time.time()
    bam = pysam.AlignmentFile(sb.DEFAULT_BAM, "rb")
    refs = set(bam.references)

    with open(sb.IN_TSV) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(hdr, line.rstrip("\n").split("\t")))
                for line in fh if line.strip()]

    calls = []
    max_fisher_div = 0.0
    n_keep_3 = n_keep_2 = n_zero_used = n_high = 0
    for i, r in enumerate(rows):
        chrom = r["chrom"]
        donor = int(r["donor"])
        acceptor = int(r["acceptor"])
        rec = process(bam, refs, chrom, donor, acceptor)
        rec["idx"] = i
        rec["gene"] = r["gene"]
        rec["donor_raw"] = donor
        rec["acceptor_raw"] = acceptor
        calls.append(rec)
        # tallies + fisher divergence (finite tables only)
        s = struct.unpack("<d", bytes.fromhex(rec["sor_bits"]))[0]
        if not math.isnan(s):
            if s < sb.SOR_KEEP:
                n_keep_3 += 1
            if s < sb.SOR_STRICT:
                n_keep_2 += 1
            if s >= sb.SOR_KEEP:
                n_high += 1
        if (rec["R_fwd"] + rec["R_rev"]) == 0 and rec["chrom_in_refs"]:
            n_zero_used += 1
        if not rec["fisher_is_nan"]:
            fp = struct.unpack("<d", bytes.fromhex(rec["fisher_bits"]))[0]
            fr = fisher_rust_replica(rec["R_fwd"], rec["R_rev"], rec["A_fwd"], rec["A_rev"])
            max_fisher_div = max(max_fisher_div, abs(fp - fr))

    # ---- synthetic edge rows (real functions; the 475 real calls have no NaN/swap) ----
    real_chrom = rows[0]["chrom"]
    real_d = int(rows[0]["donor"])
    real_a = int(rows[0]["acceptor"])
    synth = []
    # (1) chrom NOT in refs -> NaN sor/fisher, zero counts, sor_pass=False.
    s1 = process(bam, refs, "NC_FAKE_CONTIG.9", real_d, real_a)
    s1["note"] = "chrom_not_in_refs"
    s1["donor_raw"] = real_d
    s1["acceptor_raw"] = real_a
    synth.append(s1)
    # (2) donor>acceptor (swapped) -> must normalise to the SAME result as row 0.
    s2 = process(bam, refs, real_chrom, real_a, real_d)  # pass swapped
    s2["note"] = "swapped_donor_acceptor"
    s2["donor_raw"] = real_a
    s2["acceptor_raw"] = real_d
    synth.append(s2)

    fixture = {
        "source": "bench/asj_strand_bias.py",
        "bam_basename": os.path.basename(sb.DEFAULT_BAM),
        "bam_path": sb.DEFAULT_BAM,
        "in_tsv_basename": os.path.basename(sb.IN_TSV),
        "sor_keep": sb.SOR_KEEP,
        "sor_strict": sb.SOR_STRICT,
        "n_calls": len(calls),
        "max_fisher_div_vs_rust_replica": max_fisher_div,
        "stats": {
            "n_keep_3": n_keep_3, "n_keep_2": n_keep_2,
            "n_zero_used": n_zero_used, "n_high_sor": n_high,
        },
        "calls": calls,
        "synthetic": synth,
    }
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=1)
        fh.write("\n")

    dt = time.time() - t0
    print(f"[wrote {OUT}]  {len(calls)} real + {len(synth)} synthetic calls in {dt:.1f}s")
    print(f"  stats: {fixture['stats']}")
    print(f"  max |scipy - rust_fisher_replica| = {max_fisher_div:.3e}")
    # sanity: how many distinct sor values, min/max
    sors = [struct.unpack("<d", bytes.fromhex(c["sor_bits"]))[0]
            for c in calls if not c["sor_is_nan"]]
    print(f"  sor finite range: min={min(sors):.4f} max={max(sors):.4f} "
          f"(n_keep3={n_keep_3} n_keep2={n_keep_2} n_high={n_high} n_zero_used={n_zero_used})")


if __name__ == "__main__":
    main()
