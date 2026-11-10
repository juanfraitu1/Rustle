# PREREG r1602 (step 0, descriptive, no bars): does any read prefer each hard member over every other placement? (2026-11-10)

Follows r1601 (`docs/PREREG_hard_members_2026-11-10.md`). User goal (2026-11-10): find the isoforms and loci of multi-copy
families that StringTie and FLAIR leave as ambiguous. Before proposing any emitter: which of the 150 hard members carry a
read-level signal at all. Every read-level rule tried on this panel failed (register 1516-1536); a per-copy model where no read
prefers the copy is a guess.

## Measure (per hard member, the 150 of r1601, grouped by r1601 level on `a1.gtf`: L3+ 35 / L1-2 33 / missed 82)
Molecules = reads with a mapped, non-supplementary record (primary or secondary) on the member's strand overlapping its CAT exon
blocks (the r1521 panel's rule). Each such record's AS is compared with the molecule's genome-wide `best_as` and `second_as`
from the as_table (`~/runs/human_A119b/a1.molecules.tsv`, columns 2-3):
- decisive: AS here == best_as > second_as (or no second record);
- tied: AS here == best_as == second_as;
- weaker: AS here < best_as.
Member class: **distinguishable** >= 2 decisive molecules; **weak** exactly 1; **indistinguishable** 0.
Reported: the class counts per r1601 group and by read class (T-only / MAPQ0 / ordinary), with the family names.

## Already seen (2026-11-10, before this file), TAF11L on chr5:17.457-17.600 Mb
- CAT/RefSeq members are 597 bp single-exon (CDS-only) blocks at a 3,434 bp period; reads there are 1.5-3.5 kb, unspliced
  (1,777 of 1,862 records), 44 primaries and 1,818 secondaries.
- Read ends cluster at unit offset 2,200-2,400 (731 of 1,462) and starts at 2,400-2,800: one ~3 kb transcript per unit,
  longer than the CAT block, so r1601's single-exon L3 rule cannot pass an honest per-unit model there.
- `a1.gtf` has two single-exon models over the array: 113 kb (`+`) and 30 kb (`+`, over the six `-` strand members; primaries
  there 11 reverse / 6 forward).

## Use
No bars. The partition decides what to propose next: an emitter only for distinguishable members; for indistinguishable
groups, an equivalence-class product (one structure, the set of compatible loci, a family-level count).

## Outcome
(empty until the run)
