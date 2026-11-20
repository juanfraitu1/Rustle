# Row 1837: isoform fraction and sub-chain ratio re-set on real read counts (pre-registered 2026-11-20)

User (2026-11-20), after row 1836 (`docs/PREREG_dup_count_2026-11-20.md`): re-set them on real counts. Under
`RUSTLE_DUP_COUNT=1` (floor on distinct records, `reads` = every record) precision rises everywhere but correct chains fall 1-8%,
because the isoform fraction (0.02 of the locus's best transcript) and the sub-chain ratio (0.3 of the containing chain) were set
on de-duplicated counts. Workspace `r1837_real_count_grid`; binaries of `rustle_target_b`. Registered before any arm ran.

## Step 1, the grid (development: the seven datasets of rows 1833 / 1836)
Arms under `RUSTLE_DUP_COUNT=1`: `--polish-isoform-fraction` 0.02, 0.01, 0.005, 0.002 x `--polish-ism-ratio` 0.3, 0.15 (eight
arms), against `base` (today's defaults, collapsed counts). Stage `assemble`, `RUSTLE_PTC=0`. Unit: exact intron chains,
`PREFIX.gtf`; `PREFIX.final.gtf` and the SIRV contigs shown; gained and lost chains separately.

**Choice, fixed now.** An arm *qualifies* on a dataset when its precision is higher than base's and it has at least 98% of
base's correct chains (net). The choice is among the arms that qualify on A119b and on at least five of the other six
datasets: the one with the largest fraction, then the largest ratio (the least change). None = no candidate; the table goes
to the user.
Why 98% and not row 1793's 95%: the aim is a default that is better on precision at nearly no cost in recall, on A119b too.

## Step 2, check of the chosen arm
Unseen under any counting variant: bonobo testis (`PPA_mm`), KB3781 gorilla fibroblast. Reused, said: mouse adrenal, mouse
cortex, chimpanzee testis (check data of row 1836 at the old values). Run only after the choice is written here.
**Bar:** the chosen arm qualifies (precision higher, >= 98% of base's correct chains) on all five. **S1:** SIRV precision >= 0.80.
**Decision.** Both = candidate default (`RUSTLE_DUP_COUNT=1` with the chosen values); the user rules. Otherwise reported.

## Limits, stated in advance
- Lower fractions give back correct chains and wrong chains together; the grid may show no arm that keeps the precision gain.
- Single-exon transcripts are not re-set here: their floor still compares de-duplicated unspliced counts with multi-exon counts
  that include copies (row 1836's side effect). Single-exon models per arm are reported for the chosen arm only.
- The chain score filter and the 5' truncation rule keep their fit on de-duplicated counts.
- The simulated sets carry the simulator's own duplicate structure (row 1833's diagnosis).

## Result (2026-11-20): no arm qualifies by the rule (best: fraction 0.01, ratio 0.3, one dataset short); step 2 not run

Jobs 4484265-4484271, all nine arms exit 0 everywhere; `r1837_real_count_grid/score.NAME.txt`. `PREFIX.gtf`, each arm against
`base`: precision difference / share of base's correct chains kept; **bold** = qualifies (precision higher, >= 0.98 kept).

| arm (fraction, ratio) | A119b | mouse ES | fly sim | mouse sim | zebrafish sim | mouse heart | OR6737 | SIRV precision |
|---|---|---|---|---|---|---|---|---|
| 0.02, 0.3 (row 1836) | +0.014 / 0.919 | +0.043 / 0.962 | +0.014 / 0.966 | +0.024 / 0.926 | +0.023 / 0.980 | **+0.025 / 0.987** | +0.031 / 0.953 | 0.877 |
| 0.02, 0.15 | -0.003 / 0.927 | +0.030 / 0.978 | -0.007 / 0.970 | +0.003 / 0.941 | -0.009 / 0.986 | **+0.004 / 0.998** | +0.003 / 0.955 | 0.864 |
| **0.01, 0.3** | **+0.002 / 0.992** | -0.007 / 0.995 | **+0.013 / 0.988** | +0.017 / 0.966 | **+0.021 / 0.990** | **+0.006 / 0.996** | **+0.006 / 0.998** | 0.877 |
| 0.01, 0.15 | -0.013 / 1.001 | -0.017 / 1.012 | -0.008 / 0.992 | -0.004 / 0.982 | -0.012 / 0.997 | -0.013 / 1.007 | -0.018 / 1.001 | 0.864 |
| 0.005, 0.3 | -0.007 / 1.037 | -0.048 / 1.016 | **+0.011 / 1.003** | **+0.011 / 0.991** | **+0.019 / 0.998** | -0.006 / 1.001 | -0.013 / 1.021 | 0.820 |
| 0.005, 0.15 | -0.021 / 1.046 | -0.055 / 1.033 | -0.010 / 1.007 | -0.011 / 1.007 | -0.013 / 1.004 | -0.025 / 1.011 | -0.034 / 1.025 | 0.810 |
| 0.002, 0.3 | -0.014 / 1.063 | -0.081 / 1.031 | **+0.008 / 1.014** | **+0.000 / 1.009** | **+0.015 / 1.002** | -0.018 / 1.003 | -0.029 / 1.033 | 0.658 |
| 0.002, 0.15 | -0.027 / 1.072 | -0.086 / 1.048 | -0.013 / 1.018 | -0.021 / 1.024 | -0.018 / 1.009 | -0.036 / 1.014 | -0.048 / 1.037 | 0.654 |

**Choice by the rule: none.** Fraction 0.01 with ratio 0.3 qualifies on A119b and on four of the other six (needed: five);
it misses on mouse ES (precision 0.331 -> 0.324) and mouse simulated (0.966 of the correct chains). Step 2 was not run:
bonobo testis and KB3781 stay unseen under any counting variant.

**What the grid shows (my interpretation).**
- The sub-chain ratio should stay 0.3: every 0.15 arm is less precise than its 0.3 partner by 0.010-0.021 for 0.3-1.5% more
  correct chains.
- The isoform fraction is the dial. Under real counts 0.02 is a strict setting (precision +0.014 to +0.043, 2-8% of the correct
  chains lost), 0.01 is about neutral against today's default (A119b: 297 correct chains fewer, 5,376 wrong chains fewer; the
  other datasets within +-0.02 precision and >= 0.966 kept), 0.005 and lower give more correct chains than the default
  (A119b +3.7%) at lower precision on the real libraries.
- SIRV keeps its gain down to fraction 0.01 (0.877 against base 0.589) and loses it below (0.82, 0.66): the shifted-junction
  chains there sit between 0.5% and 1% of the true isoform's reads.
- So real counts do not give a default that is better on both axes everywhere; they give a cleaner dial. At 0.01 the SIRV
  problem is fixed at almost no cost elsewhere, with mouse ES and mouse simulated the two datasets slightly worse.
Not done: a fraction between 0.01 and 0.02 (it would be chosen after seeing this table); the single-exon side.
