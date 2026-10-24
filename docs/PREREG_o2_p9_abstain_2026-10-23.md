# Pre-registration: O2 v2 on fresh gorilla fibroblast families with a synthetic near-identical duplicate (P9), 2026-10-23

**Status:** committed before any fresh family is run with `o2dup.py`. Goal milestone 3.

**Criterion 3, reframed by the user (2026-10-23, option i).** O2 places AS-tied reads with <= 2% wrong and abstains on ties that
RNA cannot decide; it never places reads on a copy the sample lacks. "More fixes than breaks" is no longer required. In a copy set
with free allele slots (homozygous copies), no tie is decidable from RNA (Theorem A), so abstaining there is the right answer.

**Method.** `o3port.py` as of ee6fb354, default options (O2 v2).

**Substrate** (`bench/o3_graph/o2dup.py`, `O2DUP_CONV=1`, 100 reads per allele, seed 1; development: GWFAM248, GWFAM415):
- copy X (most allele differences; synthetic B allele with 6 substitutions when < 2) duplicated as X' with d substitutions on its
  own contig, homozygous;
- the sample's X B allele carries X''s base at X''s first substitution (a conversion) when d > 0;
- case control (the sample has X'), case missing (the sample lacks X'; X' stays in the reference and the table).

**Panel.** Every fresh >= 3-copy fibroblast family with <= 13 copies (outside the 53 development families, refabsent's 34,
GWFAM248 and P7's four), sorted by name: GWFAM135, GWFAM264, GWFAM332, GWFAM34, GWFAM36, GWFAM419, GWFAM433, GWFAM441.
GWFAM2 (27 copies) and GWFAM1 (36) are left out for runtime (P7's 27-copy run exceeded the time limit). A family whose run fails
(no isoform, time limit) is reported and skipped. d in {0, 2}.

**Truth.** Each read's origin.

## Decision rules (pooled over every run that completes)
- **Q1:** wrong <= 2% of placed tied reads (with nothing placed, Q1 holds: 0 wrong).
- **Q2:** 0 tied reads placed on X' in the missing case (a copy the sample lacks).
- **Verdict.** Pass iff Q1 and Q2. Run once; nothing changed after looking.

**Reported, not decided:** tied reads, placed, abstained, fixes and breaks against minimap2, minimap2's wrong count, runtime.
