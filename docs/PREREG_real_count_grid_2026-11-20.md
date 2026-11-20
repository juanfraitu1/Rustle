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
