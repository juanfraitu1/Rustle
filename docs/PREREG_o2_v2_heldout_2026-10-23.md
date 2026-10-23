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
