# Landing loci with near-best support, P2 control (2026-10-23)

**Question (user):** are there loci that only secondary alignments support, and are they copies the sample has ("phantom")?

**Data.** The port's landing-loci scan on the P2 control (53 real gorilla KB3781 families, unmasked reference; run
`/mnt/linuxdisk/tmp/rna_allele/p31/cand.loci.tsv`): loci outside a family's copies where >= 6 of its reads land (primary or
secondary). Annotation: `winloci_data/GGO_genomic.gff` genes / pseudogenes; other families' copies from `intervals.tsv`.

## Counts
| | loci |
|---|---|
| all landing loci | 746 in 40 families |
| only weaker (secondary, > 10 AS below the read's best) records | **722** |
| >= 6 reads within 10 AS of their best (near-best) | 13 |

**The 722 secondary-only loci** are reference loci the sample's reads resemble but do not come from: silent paralogs and pseudogenes,
as in the one-copy simulation (a silent pseudogene collected 516 secondaries, none near-best). They mark family members for O1; they
are not copies the sample expresses, and not O3 copies (a copy the reference lacks puts its reads on the nearest copy, primary or
tied, which the certified / putative tiers read).

## The 13 near-best loci
| class | loci | what they are |
|---|---|---|
| **another family's copy** | 7 (GWFAM169 x3, GWFAM401 x2, GWFAM402, GWFAM247) | the locus is a copy of GWFAM401 / 407 / 169 / 269: the reads are shared across two O1 families (an O1 boundary, the 16q NPIP-region families and TRIM5 / TRIM6 / TRIM34), not a new copy |
| **annotated gene, no family, strong** | 2: GWFAM123 at NC_073240.2:60.13-60.24 Mb (296 near-best of 516; MYO1E / FAM81A region incl. pseudogene LOC115930967); GWFAM99 at NC_073247.2:103,873,272-103,874,657 (151 of 323; PCDH11X with pseudogene LOC101144990) | expressed loci missing from the copy table: candidate unlisted copies (likely the annotated pseudogenes being transcribed) |
| **annotated gene, no family, weak** | 4: GWFAM169 / 401 at NC_073242.2:16.19-16.31 Mb (19 and 10), GWFAM268 (8 of 778), GWFAM28 (14 of 972) | a handful of near-best reads beside hundreds of secondaries: mostly paralog spill-over, a few reads possibly from the locus |
| **unannotated** | 0 | |

**Answer.** Secondary-only loci exist in large numbers (722 / 746) and are silent reference paralogs / pseudogenes, not phantom copies.
Of the 13 with near-best support, none is unannotated: 7 are other families' copies (O1 boundary), 2 are strong candidates for
expressed copies missing from the copy table (both overlap annotated pseudogenes), 4 are weak. Nothing here is a copy absent from
the reference. Not followed up: whether the GWFAM123 / GWFAM99 loci belong in those families' copy tables (O1 completeness).
