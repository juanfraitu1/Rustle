# PREREG lr-kallisto PLACEMENT ON A NON-SELF SAMPLE, OR6737 (register 1650, 2026-11-12)

Why: r1647's k-mer rule was precise on KB3781 and HG002, where reads and reference are one individual; reach fell
four-fold on HG002's non-reference haplotype (16 molecules). User: run the other gorilla (OR6737 testis on Jim's
assembly).

Limit, stated first: OR6737 has NO read-level truth. This row cannot measure the wrong rate. It measures reach and
agreement only.

Set: `panel.py` on a rebuilt extended `as_table` of `GGO_mm.bam` (amendment before any result: the workspace table lacks the position columns, the first job found 0 target reads and stopped) + `~/runs/gorilla_OR6737/or.fam.copies.tsv` (top 80, 2 Mb cost
guard, first 40), `copy_assign --tie-ms 10` (binary `target_r17`); molecules = contested rows. Reference =
`multi_iso_ref.py ... other` on `or.*` + `or.gtf`. Arm = r1647's kallisto line and rule, unchanged.
Baseline, like for like: the same kallisto rule on ALL KB3781 contested molecules of r1514's runs (r1647's 0.60 was on
molecules filtered to those with truth).

Measures:
1. Reach = assigned / contested molecules, beside the KB3781 baseline. Bar: OR6737 reach >= half of KB3781's.
2. Disagreement with O2 where both assign (pick's span does not overlap O2's copy). Bar: <= 5% (0 of 259 and 1 of 346
   in r1647).
3. Description only: share of assigned molecules whose pick's span does not hold the primary.
Verdict: both bars = the rule transfers in reach and agrees with the PSV certificate on a second individual; accuracy
stays unmeasured. Any bar failing = recorded as such.
Truth-bearing alternative, not run here: HG002 reads on CHM13 with the r1644 haplotype indexes as truth.

Outputs: `/anvme/workspace/v138da13-rustle/r1650_kallisto_or/`.

## RESULT r1650 (2026-11-12, job 4467054, 228 s; `r1650_kallisto_or/result.txt`)
OR6737 panel: 40 families, 41 regions, 13.5 Mb; `copy_assign --tie-ms 10` 6,496 contested rows = 5,971 molecules, O2
assigns 729. kallisto pseudoaligned 4,497 of 5,971 (0.753).

| | OR6737 (not the reference individual) | KB3781 baseline (all 37,611 contested molecules) |
|---|---|---|
| kallisto assigns | 2,259 (reach 0.378) | 9,479 (0.252) |
| abstains: class spans several loci | 2,215 | 23,273 |
| both kallisto and O2 assign | 212 | 667 |
| of those, different locus | 3 (0.014) | 19 (0.028) |
| kallisto-only molecules | 2,047 | 8,812 |
| pick does not overlap the primary record (description) | 306 of 2,259 (13.5%) | 320 of 9,479 (3.4%) |

Bar 1 PASS (0.378 >= half of 0.252). Bar 2 PASS (0.014 <= 0.05). **Both bars pass**: on a second individual the rule
keeps its reach and agrees with the PSV certificate where both speak.
What this does not show: the wrong rate. 13.5% of OR6737 picks move the read off its primary's locus, against 3.4% on
KB3781; without truth those moves cannot be called right or wrong, and they are exactly where allele differences
between OR6737 and the reference would mislead an exact-k-mer rule. The KB3781 reach is lower because its contested
set is dominated by exact-tie walls (23,273 multi-locus classes), not because the rule transfers better.
Next, if pursued: the truth-bearing non-self design (HG002 reads on CHM13, r1644 haplotype indexes as truth).
