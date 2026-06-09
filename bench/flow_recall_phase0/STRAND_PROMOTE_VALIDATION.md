# RUSTLE_STRAND_PURE_MINORITY — default-promotion validation

**Date:** 2026-06-09. Validates the strand-aware bundling flag as a default candidate after the
flow/filter recall levers were falsified (see PHASE0_RESULTS.md). Baseline `-L`, per-chrom
`c.bam`, gffcompare vs RefSeq annotation. Default vs `+RUSTLE_STRAND_PURE_MINORITY`.

## Targeted: recovers the strand_swamp st_only misses

`strand_rescue_check.py` on the 9 strand_swamp cases (antisense-dominated convergent cut-offs):
default `-L` recovers **2/9** target chains; `+flag` recovers **8/9** → **+6**. (Earlier
"ruled out strand-bundling" note was from NC_073247.2, which has NO strand_swamp cases — wrong
test chromosome.)

## Genome-wide safety: 8 chromosomes, +20 FSM, ZERO losses

| chrom | def FSM | +flag FSM | net | losses | IC Pr def→flag |
|---|---:|---:|---:|---:|---|
| NC_086017.1 | 1422 | 1428 | +6 | 0 | 33.7→33.8 |
| NC_073244.2 | 1430 | 1435 | +5 | 0 | 36.9→37.0 |
| NC_073233.2 | 1300 | 1303 | +3 | 0 | 38.6→38.7 |
| NC_073236.2 | 850 | 853 | +3 | 0 | 29.7→29.7 |
| NC_073242.2 | 1013 | 1015 | +2 | 0 | 31.5→31.5 |
| NC_073227.2 (control) | 1023 | 1024 | +1 | 0 | 34.7→34.6 |
| NC_073240.2 | 787 | 787 | 0 | 0 | — |
| NC_073231.2 (control) | 909 | 909 | 0 | 0 | — |
| **total** | | | **+20** | **0** | flat-to-+ |

- **Purely additive:** 0 FSM losses on ANY chromosome, including the two controls (no
  strand_swamp). Worst case is neutral, never negative — the non-reshape property.
- **Precision neutral-to-positive:** intron-chain Pr flat or +0.1; only +2–5 predicted tx per
  chrom. Control NC_073227 Pr −0.1 (negligible).
- **Safe:** ~16s/run, peak RSS ≤1.7GB, no OOM.

Extrapolates to roughly +60 FSM genome-wide (26 chroms), consistent with the
`cutoff_bundle_investigation` "+30 cut-offs / 57 rescuable."

## GENOME-WIDE validation (all 26 chroms) — net +48 FSM / 2 losses

`/tmp/strand_safety/allchrom.sh`, `-L` ± flag, gffcompare vs annotation:

- **TOTAL: default 24,374 → flag 24,422 = net +48 FSM** (50 gains − 2 losses, ~25:1).
- 24/26 chroms strictly net ≥ 0 with 0 losses; biggest gains NC_073234 +7, NC_073224 +6,
  NC_086017 +6, NC_073244 +5. Peak RSS ≤2.1GB, no OOM.
- **2 losses, both convergent-collision behavior** (the cut-off memory's predicted "5 collisions"):
  - NC_073237.2: lost TRMT10B (−) — its annotated isoform displaced by rescuing the OVERLAPPING
    convergent partner EXOSC3 (+) (25.18M overlap). A strand SWAP; chrom still net +1 (unrelated
    clean gain PIERCE1). Inherent tradeoff of strand-aware assembly at a fully-overlapping pair.
  - NC_073247.2: lost GLUD2 (+) at 132.69M — collision at its own locus (rescued partner didn't
    FSM-match). Non-colocated with that chrom's gain (MAP3K15, 29.8M). Pure small collision cost.

**Verdict: default-worthy.** +48 net, losses understood and minimal (~25:1), OOM-safe,
precision neutral-to-positive. Proceed to flip with an opt-out env var.

## Promotion options
1. Flip default (opt-OUT env var) — 3 gate sites: `cross_strand_predcluster.rs:435`,
   `pipeline.rs:12067`, `pipeline.rs:17513`. Validate all 26 chroms first OR trust the 8-chrom +
   prior cut-off validation.
2. Keep opt-in, document the validated win.

Artifacts: `strand_rescue_check.py`; sweep harness `/tmp/strand_safety/multichrom.sh`.
