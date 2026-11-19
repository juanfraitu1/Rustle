# Secondary alignments as shadow loci: real missed loci, or shadows? (2026-11-19, PSU, register 1825)

The advisor's idea: even if a gene is not expressed, its paralogs or pseudogenes should be findable through the secondary
alignments of the expressed copy's reads, and in some cases this could help seeding. The user's question: are there real loci
covered only by multi-mappers (real, missed), or are they all shadows of other loci that should stay undiscovered?
Development only: HG002 (all four quarters pooled, spent; reuse per the user's ruling) and KB3781, own-haplotype truth.
Script `bench/famwise_eval/shadow_census.py`; workspace `goal_multicopy_dev/shadow` (linked from `r1825_shadow_loci`).

## What was known before today
- Loci the ASSEMBLER builds from secondaries alone are mostly shadows: KB3781 strict truth 0 real / 155 shadow (row 1671);
  the shadow flag failed on HG002, its real ones being rDNA scaffold copies (rows 1672, 1673). Wider seeding bands add chains
  mostly at shadow loci (row 1767).
- As STRUCTURE (the pseudogene idea proper): adding the undecidable reads' secondaries gives the annotated structure to 1 copy
  primaries miss (row 1820, C1 fails); borrowed structure at paralogs adds nothing judged (row 1822, B1 fails).
- Family members found only when secondaries seed: 11 of 15 have no read from them (`MULTICOPY_FAMILIES_2026-11-19.md` step 7).

## New: a census from the reads, not from an assembly
Secondary-only place = records within 10 ms of their molecule's best, merged on contig + strand + overlap, >= 2 molecules,
and no primary record of any molecule there. Classed by the truth: REAL_DECIDED (>= 1 molecule whose single origin is the
place; its primary went elsewhere), REAL_SHARED (a molecule fitting this place and another copy equally), SHADOW (judged
molecules, every origin elsewhere), UNJUDGED.
| sample | places | REAL_DECIDED | REAL_SHARED | SHADOW | UNJUDGED |
| HG002 on CHM13 (170,727 molecules) | 51 | 6 | 11 (median 4 molecules) | 3 | 31 |
| HG002, exact ms ties only | 19 | 2 | 6 | 2 | 9 |
| KB3781 on its own reference | 12 | 0 | 4 (2-4) | 2 | 6 |
The six real places on HG002, molecules truly from the place / molecules at the place: CBWD4P (unprocessed pseudogene) 32 /
782; NUDT4P2 3 / 117; POTEB 2 / 3; MSTRG.3185 2 / 86; SRGAP2D 1 / 32; LRRC37A4P 1 / 20. Neither Rustle (primaries) nor
StringTie has a multi-exon transcript at any of the 51.

## Answer
- **Yes, a few real loci are covered only by multi-mappers**, mostly pseudogenes (CBWD4P, NUDT4P2, POTEB, LRRC37A4P) whose own
  reads all lost the primary to a paralog. On another individual's reference (HG002 on CHM13) there are 6 in a whole sample;
  on the individual's own reference, none.
- **But at each of them the reads are overwhelmingly another copy's** (CBWD4P: 96% of 782). Seeding the place from its
  secondaries would build a real locus from the wrong reads: its structure and depth would be the paralog's. So they should
  not be seeded by default, and a shadow place should not be reported as expression.
- **The rest are shared (a copy set, not separable) or shadows.** For the pseudogene idea as structure (find a paralog's
  structure through secondaries) rows 1820 and 1822 already say: little beyond projection, and it must stay labelled.
- **Where it could still help:** the few real places are exactly where a primary was wrong and a copy-specific allele could
  separate the real minority (CBWD4P's 32 molecules). A certificate like copy_phasing's, restricted to secondary-only places,
  could seed only those molecules. The ceiling is tiny (the 6 places; row 1824: correcting every wrong primary adds 0.2%).

## Caveats
HG002 is a cell line on another individual's reference; the selection holds reads with a primary on a truth-family member;
UNJUDGED places (31 of 51) may hide more; no testis sample has such a truth.
