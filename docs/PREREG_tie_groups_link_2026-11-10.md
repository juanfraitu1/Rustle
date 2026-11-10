# PREREG r1605: tie groups — a scorer that sees the whole family, and a link rule against chaining (2026-11-10)

Follows r1604 (`docs/PREREG_tie_groups_near_2026-11-10.md`): near ties reach 72/72 and TAF11L 31/31 but purity is 0.80.
Post hoc there: the scorer's family table lists expressed members only and has no pseudogenes; groups chain (one of 831 loci).
User (2026-11-10): do both fixes. Written and committed before any arm runs. Tied rule = as98 throughout (asms98 added
88k records and no reach in r1604).

## 1. Scorer (`tie_groups.py score ... full`)
- Family of a gene = its component in the full primates-level truth (`famtruth.truth_pairs` over every CAT gene), not the
  expressed panel rows. A group's family = the most common such family among the genes its loci overlap.
- A pseudogene (biotype contains `pseudogene`) counts as that family when the alphabetic prefix of its symbol equals a
  family member's prefix, or equals it after dropping a final `P` (TSPY13P ~ TSPY1, AMYP1 ~ AMY1A). No other biotype gets it.
- Judged groups, as before: those with a locus reaching a panel member. A group is pure when every locus that overlaps a CAT
  gene overlaps one in the family under these two rules. lncRNA, TEC and other genes outside the family stay foreign.
- **Grouped reach**: a member counts only when a locus reaching it is in a group of >= 2 loci (an ambiguous copy alone in
  its group is not a family-level product).

## 2. Link rule (`tie_groups.py build PREFIX as98 FRACTION`)
Two loci are linked when they share >= 2 tied molecules **and** that count is >= FRACTION x the molecule count of the smaller
locus. Arms: FRACTION 0 (control = r1604 as98), 0.25, 0.5.

## Bars (per arm; M = the 72 indistinguishable members)
- B1g grouped reach >= 58 of 72.
- B2f purity (new scorer) >= 0.90.
- B3 TAF11L members reached by a model <= 5 kb >= 20 of 31 (as r1603).
- B5 largest group <= 60 loci.
- B4 the job <= 1 h wall.
- An arm PASSES on all five. If 0.25 and 0.5 both pass, the one with the higher purity is kept; a tie goes to 0.25.
- The control under the old scorer must reproduce r1604 (72 / 0.803 / 31).
- Not addressed, reported only: groups holding different genes of one segmental duplication (NPIPB with BOLA2 / EIF3C).

## Outcome
(empty until the run)
