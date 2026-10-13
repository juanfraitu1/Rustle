# PREREG — a log-odds, semi-global, anchored-band engine for O2's copy-vs-copy PSV alignment (2026-10-12)

**Why.** `PREREG_psv_significant_blocks_2026-10-12.md` (B2 FAILED, mask opt-in) showed that `banded_msa_pair`'s cost
model (mismatch 1, gap open 32, extend 1; global) makes a ~30-base misalignment cheaper than one gap: inside stretches
an independent aligner supports, 9,703 of 16,183 of O2's apparent copy differences (60%) sit at the wrong position, and
14% of all PSV pairs come from non-homologous sequence forced into the global path. The user also asked that this part
be VERY FAST.

## Engine (fixed before any held-out data is touched; may be tuned on DEV only, each change logged as an amendment)

`semiglobal_msa_pair(a, b)` — same 2-row gapped-MSA output as today's engines, so `discover_psvs` is unchanged downstream:
- **Scores** (maximised): match +2, mismatch or non-ACGT -3; **two-piece affine gaps** cost min(5 + 2k, 24 + k) for a gap
  of k (minimap2-style: a long gap — an exon one copy lacks — costs ~1 / base, less than the -1.75 / base expectation of
  aligning unrelated sequence).
- **Free end gaps** in both sequences (overlap alignment): an unshared 5' / 3' extent is left unaligned, not forced.
- **Anchored band**: unique shared 11-mers, longest colinear chain; row i between anchors (i1, j1), (i2, j2) is
  computed for j in [i + min(d1, d2) - W, i + max(d1, d2) + W] (d = j - i), W = 100; before the first / after the last
  anchor the band follows that anchor's diagonal. No anchor -> no alignment (no shared 11-mer: no detectable homology).
- **Guard**: `significant_block_mask` (Karlin-Altschul, alpha 1e-3) on the result.
- Selected by `RUSTLE_PSV_ENGINE=semiglobal` (opt-in) until this prereg passes; default unchanged.

## Harness (no BAM)

`utilities psv-discover --families COPIES.tsv --copies-fa COPIES.fa --out PREFIX [--engine global|semiglobal]`
runs `discover_psvs_oriented` per family and writes `PREFIX.psv_cols.tsv` / `.psv_copies.tsv` (the dump's formats),
`PREFIX.psv_pos.tsv` (every copy's spliced offset and base per column) and `PREFIX.time.tsv` (wall per family).
Validated first: `--engine global` reproduces the 10-10 `copy_assign --dump-psv` columns on the dev subset (a column
the RNA-editing pre-pass removed in copy_assign may be extra in the harness; nothing else may differ).

## Yardstick (`bench/psv_sun/engine_eval.py`)

Per family, copy 0 vs each copy X in copy 0's frame, the colinear chain of `minimap2 -cx asm20` pieces (independent of
both engines). A **difference pair** = (X, column) where X's allele differs from copy 0's.
- In a chained piece: **V** if the engine's X offset equals minimap2's mapping of copy 0's offset, else **M**.
- Outside every piece: Smith-Waterman (`swscore`, BLASTN scores) of the surrounding unaligned stretches vs 20 shuffles:
  **H** (beats all; divergent homolog) or **N** (does not; unrelated).

## Substrates

- DEV (may be looked at freely): the 56 KB3781 families of `PREREG_psv_sun_pileup_2026-10-12.md` (`run/S.txt`).
- **HELD-OUT 1**: the 269 other KB3781 families with 2-20 copies (`rustle_figures/runs/gorilla_KB3781/gorilla_KB3781.fam.copies.tsv`).
- **HELD-OUT 2** (another species and sample): the 2,101 human A119b families with 2-20 copies
  (`rustle_figures/runs/human_A119b/human_A119b.fam.copies.tsv`).

## Bars (each held-out substrate separately; semiglobal vs global)

- **P1 precision**: V / (V + M) >= 0.95 for semiglobal.
- **P2 no loss**: V(semiglobal) >= 0.99 x V(global).
- **P3 junk**: N(semiglobal) <= 0.10 x N(global).
- **P4 speed**: total harness wall (same threads) semiglobal <= global; and no family slower than 2 s.
- **P5**: library tests pass; on DEV `copy_assign` default assignments stay byte-identical with the engine switched.
- Reported, no bar: H counts, total columns, per-family time distribution.
All of P1-P4 on BOTH held-out substrates -> the engine becomes the default. Otherwise it stays opt-in and the failing
bar is reported.

## Dev amendments (2026-10-12, DEV only; written before any held-out family is aligned or scored)

- **E1 engine (dev test failure):** free end gaps let the path start on the first row OR column only, so two unshared
  prefixes could not both be skipped and 60 unrelated bases were aligned (unit test). The engine is now LOCAL at its ends
  (Smith-Waterman: start / end anywhere) and the long-gap piece is a CONSTANT, cost min(5 + 2k, 30) (like an intron in a
  spliced alignment): a missing exon of any length costs 30, and unrelated stretches on both sides can be skipped by two
  gaps. Scores, anchors, band and the significant-block guard unchanged. 6 engine tests + 1,029 library tests pass.
  Anchor limit documented: a PERIODIC 20% difference pattern shares no exact 11-mer and gets no alignment.
- **Harness validated:** `psv-discover --engine global` reproduces the dev `copy_assign --dump-psv` tables exactly
  (sorted: 33,016 / 209 lines) in 6 s instead of a 5-min BAM run.
- **A1 yardstick (dev P1 = 0.736):** 1,567 of the 1,725 M pairs sat > 20 bp from minimap2's position, in copies with
  INTERNAL repeats: MCL8's copy 0 is a ~42-bp-period tandem array (323 off-diagonal self-alignments), MCL37 an 84-bp
  shift, MCL32's copy X a duplicated segment matching copy 0's offset equally well at two places (local identity 0.84
  vs 0.84). Positional homology is undefined there for ANY aligner. New class **R**: a pair inside a minimap2 piece whose
  copy-0 offset lies in an off-diagonal self-alignment of copy 0, or whose minimap2 X position lies in one of X
  (`minimap2 -cx asm20 -X -N 50 -p 0.1`; engine-independent). R is reported and excluded from V and M. Bars unchanged.
- Dev with A1: global V 3,192 / M 7,543 (0.297), R 5,448, H 27,138, N 15,058; semiglobal V 4,239 / M 29 (**0.993**),
  R 2,259, H 1,501, N 10; harness 10.16 s vs 1.90 s summed, slowest family 2.27 s vs 0.72 s.
