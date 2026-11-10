# PREREG r1604: tie groups with near ties on AS, and on AS or ms (2026-11-10)

Follows r1603 (`docs/PREREG_tie_groups_2026-11-10.md`, FAIL: reach 39/72, purity 0.858, TAF11L 11/31). Post hoc there: exact AS
equality keeps 1-2 molecules per TAF11L unit. User (2026-11-10): ties are read on AS and on ms (minimap2's ranking score), so
the tie set may need to be bigger. Written and committed before any arm runs.

## Arms (only the tied-record rule changes; loci, groups, models, panel and M as in r1603)
- **as98**: a record is tied when AS >= 0.98 x the molecule's best AS (the r1521 panel's T rule).
- **asms98**: as98, or ms >= 0.98 x the molecule's best ms. A superset of as98.
- **exact** (r1603) is rebuilt from the saved records as the control; it must reproduce 39 / 0.858 / 11.
A molecule is tied when >= 2 of its records are tied under the arm's rule.

## Run
The molecules' records are saved once (`samtools view -N names.txt -b`, r1603's `names.txt`), so each arm is a build from that
file, not a BAM pass.

## Bars (per arm, as r1603)
- B1 reach >= 58 of 72; B2 purity >= 0.90; B3 TAF11L models <= 5 kb reach >= 20 of 31; B4 the whole job <= 1 h wall.
- An arm PASSES on all four. If both pass, asms98 is kept only if its purity is not below as98's; otherwise as98.
- Any miss in an arm = that arm FAILS as registered, with the numbers.
- Descriptive: tied records, loci, groups, single-locus groups, loci with no CAT gene, per arm.

## Outcome
(empty until the run)
