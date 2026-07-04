#!/usr/bin/env python
"""BYTE-PARITY golden-fixture generator for the Rust port of the canonical
(k, w)-minimizers function (`vg_family/minimizers.rs`) -- the FOUNDATION + byte-
parity crux of the O1 over-merge-gate migration (RETIREMENT_AND_MIGRATION.md
Tier-2 #7, full-port path). Every repeat-bridge gate rests on this ONE function.

It imports the SHIPPED Python module verbatim -- `vg_repeat_catalog.minimizers`
(+ its K/W/_B constants) -- and runs the REAL function over a DIVERSE, exhaustive
sequence set, emitting to
`src/rustle/vg_family/testdata/minimizers_fixture.json`:

  { "k": 15, "w": 10, "B": {...}, "diversity": {class: count},
    "cases": [ {"name","class","seq","out":[[canon,masked_int],...]}, ... ] }

Edge classes covered (>= 60 sequences):
  empty, below_k, exactly_k, homopolymer (all-same-base ties),
  random_acgt, with_n (single N / N-runs / N at start & end -> window resets),
  mixed_case (lowercase-majority masking both ways at the k/2 boundary),
  revcomp (a seq + its reverse-complement -- canon MULTISET asserted equal here),
  tandem (repeated minimizer selections -> run-collapse),
  real_exon (a few soft-masked slices pulled from GGO.fasta if available, else skip).

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/gen_minimizers_fixture.py
Writes: src/rustle/vg_family/testdata/minimizers_fixture.json
Deterministic (fixed random seed; sorted selection; PYTHONHASHSEED=0).
"""
import os
import sys

# Deterministic hashing/selection regardless of the caller's env.
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import json
import random

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import vg_repeat_catalog as V  # noqa: E402  (the SHIPPED minimizers + K/W/_B spec)

OUT = os.path.join(HERE, "..", "src", "rustle", "vg_family", "testdata",
                   "minimizers_fixture.json")

K, W = V.K, V.W
RNG = random.Random(0xC0FFEE)  # fixed seed -> deterministic sequences

_COMP = {"A": "T", "C": "G", "G": "C", "T": "A",
         "a": "t", "c": "g", "g": "c", "t": "a", "N": "N", "n": "n"}


def revcomp(s):
    """Case-preserving reverse complement (unknown chars pass through)."""
    return "".join(_COMP.get(ch, ch) for ch in reversed(s))


def rand_seq(n, alphabet="ACGT"):
    return "".join(RNG.choice(alphabet) for _ in range(n))


def canon_multiset(out):
    """Multiset (sorted list) of the canon VALUES in a minimizers() output. This is
    reverse-complement invariant (canon = min(fwd, revcomp) is strand-agnostic; the
    per-window min VALUE is RC-mirror-invariant and run-collapse commutes with
    reversal), which we assert on every revcomp pair as a self-check."""
    return sorted(v for v, _ in out)


CASES = []          # list of dict(name, class, seq)
_seen_names = set()


def add(name, cls, seq):
    assert name not in _seen_names, f"duplicate case name {name}"
    _seen_names.add(name)
    CASES.append(dict(name=name, cls=cls, seq=seq))


# ---------------------------------------------------------------- empty / below_k
add("empty", "empty", "")
add("below_k_1", "below_k", "A")
add("below_k_k_minus_1", "below_k", "A" * (K - 1))          # 14 -> []
add("below_k_random14", "below_k", rand_seq(K - 1))
add("below_k_with_n", "below_k", "ACGT" * 3 + "NN")          # 14, still < k

# ---------------------------------------------------------------- exactly_k
add("exactly_k_polyA", "exactly_k", "A" * K)                 # one k-mer, one run
add("exactly_k_polyT", "exactly_k", "T" * K)
add("exactly_k_random_a", "exactly_k", rand_seq(K))
add("exactly_k_random_b", "exactly_k", rand_seq(K))
add("exactly_k_mixed_case", "exactly_k", "ACGTacgtACGTacg")  # majority lower? 7 lower of 15

# ---------------------------------------------------------------- homopolymer (ties)
for base in ("A", "C", "G", "T"):
    add(f"homopolymer_{base}_60", "homopolymer", base * 60)
add("homopolymer_A_200", "homopolymer", "A" * 200)
add("homopolymer_lower_g_60", "homopolymer", "g" * 60)       # masked homopolymer
# near-homopolymer: one differing base creates a distinguishable minimizer
add("near_homopolymer_A_mid", "homopolymer", "A" * 30 + "C" + "A" * 30)
add("near_homopolymer_dimer", "homopolymer", "AC" * 40)      # AC repeat = periodic ties

# ---------------------------------------------------------------- random_acgt
for i in range(12):
    add(f"random_acgt_{i}", "random_acgt", rand_seq(RNG.randint(K, 400)))
add("random_acgt_exactly_kw", "random_acgt", rand_seq(K + W - 1))  # first full window
add("random_acgt_long", "random_acgt", rand_seq(1000))

# ---------------------------------------------------------------- with_n (resets)
base = rand_seq(120)
add("with_n_single_mid", "with_n", base[:60] + "N" + base[60:])
add("with_n_run_mid", "with_n", base[:50] + "NNNNN" + base[50:])
add("with_n_start", "with_n", "NNNNN" + base)
add("with_n_end", "with_n", base + "NNNNN")
add("with_n_start_and_end", "with_n", "NNN" + base + "NNN")
add("with_n_many_runs", "with_n",
    rand_seq(40) + "NN" + rand_seq(40) + "N" + rand_seq(40) + "NNN" + rand_seq(40))
add("with_n_only", "with_n", "N" * 50)                        # no valid k-mer -> []
add("with_n_lowercase_n", "with_n", rand_seq(40) + "n" + rand_seq(40))  # 'n' not in _B
add("with_n_dense", "with_n", ("ACGTN" * 30))                 # never reaches valid>=k
# ambiguity codes (not in _B) also reset -- R, Y, plus a stray non-ACGTN byte
add("with_ambiguity_codes", "with_n", rand_seq(40) + "R" + rand_seq(40) + "Y" + rand_seq(40))
add("with_gap_dash", "with_n", rand_seq(40) + "-" + rand_seq(40))

# ---------------------------------------------------------------- mixed_case (masking)
# The mask fires iff lc_sum*2 > k, i.e. >= 8 lowercase of the 15-base window (k/2 boundary).
add("mixed_all_lower", "mixed_case", ("acgt" * 30))
add("mixed_all_upper", "mixed_case", ("ACGT" * 30))
# exactly at the boundary: build a 30-mer with a moving 8-lower / 7-lower window.
add("mixed_boundary_8lower", "mixed_case", "aaaaaaaa" + "TTTTTTT" + "ACGTACGTACGT")   # 8 lower then 7 upper
add("mixed_boundary_7lower", "mixed_case", "aaaaaaa" + "TTTTTTTT" + "ACGTACGTACGT")   # 7 lower then 8 upper
add("mixed_alternating", "mixed_case", ("aCgTaCgT" * 15))     # ~half lower -> boundary sweep
add("mixed_lower_block_then_upper", "mixed_case", rand_seq(60).lower() + rand_seq(60))
add("mixed_upper_block_then_lower", "mixed_case", rand_seq(60) + rand_seq(60).lower())
add("mixed_soft_island", "mixed_case", rand_seq(40) + rand_seq(30).lower() + rand_seq(40))

# ---------------------------------------------------------------- tandem (run-collapse)
add("tandem_dimer_AC", "tandem", "AC" * 100)
add("tandem_trimer_ACG", "tandem", "ACG" * 80)
add("tandem_5mer", "tandem", "ACGTA" * 50)
add("tandem_ATTT_micro", "tandem", "ATTT" * 60)
add("tandem_unit_then_break", "tandem", "ACGTACG" * 20 + "TTTTT" + "ACGTACG" * 20)
add("tandem_kmer_repeat", "tandem", rand_seq(K) * 15)         # exact k-mer tandem
add("tandem_with_n_break", "tandem", "ACGTAC" * 20 + "NNN" + "ACGTAC" * 20)
add("tandem_period16", "tandem", rand_seq(16) * 20)           # period > k

# ---------------------------------------------------------------- revcomp pairs
REVCOMP_PAIRS = []
_rc_specs = [
    ("rc_random_a", rand_seq(200)),
    ("rc_random_b", rand_seq(137)),
    ("rc_random_c", rand_seq(88)),
    ("rc_with_n", rand_seq(50) + "NNN" + rand_seq(50)),
    ("rc_mixed_case", rand_seq(60) + rand_seq(60).lower()),
    ("rc_tandem", "ACGTG" * 40),
    ("rc_exactly_kw", rand_seq(K + W - 1)),
    ("rc_palindrome_ish", rand_seq(80) + revcomp(rand_seq(80))),
]
for nm, s in _rc_specs:
    add(nm + "_fwd", "revcomp", s)
    add(nm + "_rev", "revcomp", revcomp(s))
    REVCOMP_PAIRS.append((nm + "_fwd", nm + "_rev"))

# ---------------------------------------------------------------- real_exon (optional)
def load_real_exons():
    """Pull a few short soft-masked slices from GGO.fasta if available (else skip).
    Uses pysam (already a dep of vg_repeat_catalog). Deterministic coordinates."""
    fa_path = getattr(V, "GENOME", None)
    if not fa_path or not os.path.exists(fa_path):
        return []
    try:
        import pysam
        fa = pysam.FastaFile(fa_path)
    except Exception:
        return []
    out = []
    try:
        # first few contigs, deterministic offsets/lengths; keep native soft-mask/N.
        specs = [
            (0, 100000, 300),
            (0, 500000, 220),
            (0, 1200000, 411),
            (1, 200000, 260) if len(fa.references) > 1 else None,
            (2, 300000, 333) if len(fa.references) > 2 else None,
            (0, 2500000, 180),
        ]
        for j, spec in enumerate(specs):
            if spec is None:
                continue
            ci, start, ln = spec
            if ci >= len(fa.references):
                continue
            chrom = fa.references[ci]
            clen = fa.get_reference_length(chrom)
            if start + ln > clen:
                continue
            s = fa.fetch(chrom, start, start + ln)
            if s and len(s) >= K:
                out.append((f"real_exon_{j}", s))
    finally:
        fa.close()
    return out


for nm, s in load_real_exons():
    add(nm, "real_exon", s)


# ---------------------------------------------------------------- run the REAL fn
def main():
    diversity = {}
    cases_out = []
    for c in CASES:
        out = V.minimizers(c["seq"])  # THE SHIPPED FUNCTION (k=K, w=W defaults)
        diversity[c["cls"]] = diversity.get(c["cls"], 0) + 1
        cases_out.append(dict(
            name=c["name"], cls=c["cls"], seq=c["seq"],
            out=[[int(v), int(m)] for (v, m) in out],
        ))

    # ---- revcomp canon-multiset self-check (must hold for canonical minimizers) ----
    by_name = {c["name"]: V.minimizers(c["seq"]) for c in CASES}
    for fwd, rev in REVCOMP_PAIRS:
        mf, mr = canon_multiset(by_name[fwd]), canon_multiset(by_name[rev])
        assert mf == mr, (
            f"REVCOMP canon-multiset mismatch for {fwd}/{rev}:\n  fwd={mf}\n  rev={mr}")

    fixture = dict(
        k=K, w=W,
        B={str(kk): vv for kk, vv in sorted(V._B.items())},
        diversity=diversity,
        n_cases=len(cases_out),
        cases=[dict(name=c["name"], **{"class": c["cls"]}, seq=c["seq"], out=c["out"])
               for c in cases_out],
    )

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(fixture, fh, indent=1, sort_keys=True)

    total_runs = sum(len(c["out"]) for c in cases_out)
    print(f"[fixture] wrote {os.path.normpath(OUT)}")
    print(f"          {len(cases_out)} sequences, {total_runs} emitted runs total")
    print(f"          classes: " +
          ", ".join(f"{k}={v}" for k, v in sorted(diversity.items())))
    print(f"          revcomp pairs self-checked: {len(REVCOMP_PAIRS)} (canon multiset equal)")


if __name__ == "__main__":
    main()
