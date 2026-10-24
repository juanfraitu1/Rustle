# Pre-registration: O2 `--semi` (two-site lead) on fresh gorilla fibroblast families (P11), 2026-10-23

**Status:** committed before any P11 family is simulated or run. Goal milestone 3, under the user's rule "the totally tied ones are
hopeless" (keep semi-tied reads and MAPQ-0 secondaries; drop only exact ties).

**Method.** `o3port.py --semi` as of 319be331 (register 1389). A tied read (>= 2 copies within 10 AS):
- certified new copy or forced haplotype edge: as O2 v2;
- exact tie (same AS and edit distance on its two best copies): HOPELESS (abstain);
- semi-tied: its best AS / fewest-edit copy, only when the runner-up needs >= 2 more edits (two distinguishing sites); else abstain.

**Development** (GWFAM248, GWFAM415; d 0 / 1 / 2): 584 placed, 2 wrong (0.34%); minimap2 MAPQ >= 1 on the same reads 40.6% wrong.

**Panel.** The first 10 families (by name) of P10's panel: GWFAM117, 118, 119, 120, 132, 133, 139, 154, 171, 172. Two-copy gorilla
fibroblast families, used for O3 (P10) but never for O2.

**Substrate.** `bench/o3_graph/o2dup.py` with `O2DUP_CONV=1`, `O2DUP_PORT_ARGS=--semi`, 100 reads per allele, seed 1, d in {1, 2, 3}:
the synthetic duplicate X' of the family's most heterozygous copy, a conversion of X's B allele at X''s first substitution; cases
control and missing (the sample lacks X'). 60 runs. A run that fails (no isoform, time limit) is reported and skipped.

**Truth.** Each read's origin. **Baselines:** minimap2's primary; minimap2's primary where MAPQ >= 1.

## Decision rules (pooled over every completed run)
- **S1:** O2 wrong <= 2% of its placed tied reads.
- **S2:** 0 tied reads placed on X' in the missing case.
- **S3:** O2's wrong rate < minimap2 MAPQ >= 1's wrong rate on the same tied reads.
- **Power:** >= 100 tied reads placed by O2, else underpowered (no pass).
- **Verdict.** Pass iff S1, S2, S3 and power. Run once; nothing changed after looking.

**Reported, not decided:** fixes and breaks against minimap2's primary; hopeless and abstain counts; per d.
