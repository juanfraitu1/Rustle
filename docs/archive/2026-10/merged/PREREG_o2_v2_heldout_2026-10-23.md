# Pre-registration: O2 v2 (forced matches, weighed on every tied copy) on fresh erased-copy simulations (P8), 2026-10-23

**Status:** committed before the families are picked or run. Goal milestone 3.

**Method.** `o3port.py` as of 4b6ba85f, default options. Its O2 output: a tied read goes to its haplotype's copy only when that
haplotype was weighed for every copy the read ties to and the edge is forced, or to the certified new copy it belongs to; otherwise
ABSTAIN.

**Panel.** `o3sim_one.py pick` (its rule: the simplest >= 3-copy family whose erased copy is >= 97% identical over >= 80% to >= 2
survivors and has >= half of its reads tied on >= 2 survivors).
- Excluded: GWFAM248 and every family used in P1-P7 or the tie scan's development families (O3SIM_EXCLUDE).
- Picked repeatedly, excluding each pick, until 3 families are found.
- Per family, the one-copy sweep's 40 configurations (expression, divergence, same vs other individual, controls), at the code's
  defaults.

**Truth.** Each read's origin. **Baseline:** minimap2's primary.

## Decision rules (pooled over every run)
- **O2-1:** wrong <= 2% of placed tied reads.
- **O2-2:** fixes > breaks against minimap2.
- **Verdict.** Pass iff both hold. Run once.

## Result (2026-10-23): cannot be run
No fresh gorilla family meets the pick rule once GWFAM248, P1-P7's families and the 53 development families are excluded. In the
remaining >= 3-copy families, no erased copy is both >= 97% identical over >= 80% to two survivors and hidden as tied multimappers.
So the fibroblast data hold no fresh held-out for O2 v2. Nothing was simulated or scored.

## After the user chose "another held-out substrate" (2026-10-23, development only; nothing fresh run)
- **Fresh pool.** 10 fresh >= 3-copy fibroblast families remain (GWFAM34, 36, 264, 332, 135, 419, 441, 433, 2, 1). The tie scan
  gives 0 AS-tied reads in every one with their real copies, so real fresh ties do not exist.
- **Synthetic duplicate substrate** (`bench/o3_graph/o2dup.py`): a fresh family plus X' = copy X with d substitutions on its own
  contig; cases control / missing; optional conversion (`O2DUP_CONV=1`: X's B allele carries X''s base at one substitution).
- **Development (GWFAM415, GWFAM248), d = 0-10, with and without conversion: O2 v2 places 0 tied reads.**
  - d = 0: the tie is symmetric (X and X' identical in the reference); RNA cannot decide; abstaining is right.
  - d >= 2 without conversion: minimap2 is already right (1 / 396 wrong); nothing to fix.
  - with conversion (minimap2 52-100 wrong): X' is homozygous, so its second slot is free and X's B haplotype fits it (Theorem A);
    no edge is forced, so v2 abstains, correctly.
- **So on this substrate O2-2 (fixes > breaks) cannot pass by construction** unless the substrate is tuned (heterozygous X', fewer
  free slots) until v2 places reads. That would be fitting the test to the method, so no P9 was pre-registered.
