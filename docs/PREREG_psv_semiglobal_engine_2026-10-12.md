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

## Held-out result 1 (2026-10-12) and Amendment E2 + fresh held-out 3 (written before any chimp family is aligned)

**HELD-OUT 1 (gorilla, 269 families) — P1-P4 PASS:** global V 16,445 / M 22,848 (precision 0.419), R 7,090, H 63,784,
N 86,288; semiglobal V **18,317** / M **31** (**0.998**), R 2,187, H 2,939, N **119**; harness 34.53 s vs **5.15 s**
summed, 0 semiglobal families > 2 s (`eng/ho1_eval/engine_eval.json`).

**HELD-OUT 2 (human, 2,101 families) — P4 FAILS for semiglobal:** 10 families > 2 s (MCL1887 34.1 s: a 265-kb copy vs an
11-kb one; MCL1904 24.5 s). The GLOBAL engine could not be run as comparator: it aborted allocating 135 GB after 158 s
(whole set) and did not finish 2,099 families (without MCL1887 / MCL1904) within 590 s on 2 threads — semiglobal did all
2,101 in 115.7 s wall. **P1 on human (semiglobal alone, 7 chunks): V 377,664 / M 1,128 = 0.997 PASS** (R 49,458, H 89,694,
N 4,505; `eng/ho2_eval/engine_eval_semi_total.json`). P2 / P3 are not computable on human (no comparator output). **Verdict as registered: not all bars pass on both held-out substrates -> the engine stays
OPT-IN.**

**Amendment E2 (band; motivated by held-out 2, so judged only on held-out 3):** between anchors the diagonal envelope
spans every column across a huge one-sided insertion (250k rows x 11k columns). Alignment paths are monotone, so each
row's band is intersected with the columns between the flanking anchors (+- W): before the first anchor columns
<= its end, after the last >= its start. New unit test: a 30-kb one-sided insertion — banded score equals the full DP
and < 10 M cells.

**HELD-OUT 3 (fresh; another species, never aligned or scored):** the 393 chimp PTR families with 2-20 copies
(`rustle_figures/runs/chimp_PTR/chimp_PTR.fam.copies.tsv`). Bars P1-P4 exactly as above (P4: total semiglobal <= global
and no semiglobal family > 2 s); yardstick with A1. If the global comparator cannot finish a family, that family is
excluded from P2/P3 and listed, and counts against global in P4. **All of P1-P4 pass on held-out 3 (with held-out 1
already passed) -> the engine becomes the default; otherwise it stays opt-in.**
