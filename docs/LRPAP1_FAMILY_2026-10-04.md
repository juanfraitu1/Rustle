# The LRPAP1 family in gorilla: eight copies, one annotated name — what O1, O2 and O3 say about it (2026-10-04, for the advisor)

LRPAP1 (LDL receptor related protein associated protein 1; RAP) is single-copy in human (CHM13 chr4:3,503,365-3,532,077; 1,717 A119b reads),
chimpanzee and orangutan (Liftoff self-lift: `in_place` only). **In gorilla (mGorGor1 primary assembly) its 20-kb, intron-containing body is
present at eight loci** on six chromosomes: the Liftoff self-lift of the gorilla annotation (`figures/_liftoff.py`, Fig. 8) projects LRPAP1 at
0.957-0.990 sequence identity with 0.978-1.000 coverage to seven other sites, five of which RefSeq annotates as "alpha-2-macroglobulin
receptor-associated protein-like" LOC genes (two protein-coding, three pseudogenes) and two of which carry no annotation at all. These are
segmental duplicates (exon-intron structure preserved), not retrocopies, and the youngest (copyB, LOC134757218: 1.4-1.5% from LRPAP1) are
nearly identical to each other (0.998). One copy sits on the Y chromosome.

Work dir `/mnt/linuxdisk/tmp/lrpap1/` (copies table `lrpap1.copies.tsv`, Liftoff-projected models `lrpap1.truth.gtf`, de novo family runs
`kb3781_6contigs.*` / `or6737_6contigs.*`, O2 runs `o2_fibro.*` / `o2_testis.*`, haplotype alignments `copies8.{mat,pat}.paf`, read support
`support_{fibro,testis}.*`). Page: LRPAP1 Family (artifact). Nothing here was pre-registered: it is a description of one family with the
shipped tools; every number is reproducible from the files named.

## The eight copies

| copy | gorilla chr | span (primary) | identity to LRPAP1 | annotation (RefSeq mGorGor1) | reads fibroblast / testis | de novo locus fibroblast | testis | dominant expressed chain (testis) | chain-support reads (C) f / t | TSS-anchored (D') f / t | expressed start |
|---|---|---|---|---|---|---|---|---|---|---|---|
| LRPAP1 | chr3 (hsa4) | NC_073227.2:12,090,719-12,110,719 (-) | 1.000 | RefSeq gene LRPAP1 (protein-coding) | 1,488 / 1,977 | DN_NC_073227.2_12086114_8 (8 exons, 206 reads) | DN_NC_073227.2_12088225_8 (8 exons, 329 reads) | 1799 reads, 7 introns | 1471 / 1932 | 1468 / 1943 | 12110552 |
| LOC134756753 | chr12 (hsa2a) | NC_073236.2:22,553,726-22,573,741 (-) | 0.990 | RefSeq LOC134756753, protein-coding, "alpha-2-macroglobulin receptor-associated protein-like" | 836 / 1,546 | DN_NC_073236.2_22551320_8 (8 exons, 142 reads) | DN_NC_073236.2_22549336_8 (8 exons, 257 reads) | 1440 reads, 7 introns | 823 / 1524 | 821 / 1528 | 22573573 |
| copyB | chr12 (hsa2a) | NC_073236.2:24,785,460-24,805,469 (+) | 0.985 | not annotated | 19 / 22 | DN_NC_073236.2_24785601_8 (8 exons, 15 reads) | DN_NC_073236.2_24785618_8 (8 exons, 9 reads) | 10 reads, 7 introns | 15 / 19 | 18 / 19 | 24785624 |
| copyC | chr12 (hsa2a) | NC_073236.2:30,205,371-30,225,550 (+) | 0.988 | not annotated | 57 / 121 | DN_NC_073236.2_30205515_7 (7 exons, 23 reads) | DN_NC_073236.2_30205531_7 (7 exons, 38 reads) | 45 reads, 6 introns | 51 / 113 | 50 / 119 | 30205536 |
| LOC129526389 | chr14 (hsa13) | NC_073238.2:23,392,677-23,412,705 (+) | 0.990 | RefSeq LOC129526389, protein-coding, "…-like" | 105 / 84 | DN_NC_073238.2_23392840_7 (7 exons, 16 reads) | DN_NC_073238.2_23392840_7 (7 exons, 28 reads) | 69 reads, 6 introns | 97 / 81 | 97 / 81 | 23392848 |
| LOC134757218 | chr16 (hsa15) | NC_073240.2:15,822,038-15,842,048 (-) | 0.986 | RefSeq LOC134757218, pseudogene | 28 / 1 | DN_NC_073240.2_15821047_8 (8 exons, 15 reads) | DN_NC_073240.2_15820894_8 (8 exons, 12 reads) | 1 reads, 6 introns | 18 / 0 | 24 / 0 | 15841883 |
| LOC129523574 | chr22 (hsa21) | NC_073246.2:11,292,022-11,312,055 (-) | 0.957 | RefSeq LOC129523574, pseudogene | 0 / 45 | no locus (no reads) | DN_NC_073246.2_11287585_8 (8 exons, 16 reads) | 44 reads, 6 introns | 0 / 45 | 0 / 45 | 11311888 |
| LOC129530227 | chrY | NC_073248.2:45,277,394-45,297,429 (+) | 0.980 | RefSeq LOC129530227, pseudogene (chrY) | 1 / 18 | no locus (no reads) | DN_NC_073248.2_45277558_8 (8 exons, 6 reads) | 7 reads, 7 introns | 0 / 14 | 0 / 14 | 45277560 |


(chain-support reads = Amendment C of `docs/PREREG_spliced_copy_support_2026-10-04.md`: reads whose junction chain is an expressed chain of
the copy or a 5' piece of it, uniquely placed; TSS-anchored = D', reads starting within 150 bp of the chain's modal start carrying its first
three introns — the gorilla libraries have no cap signal, so the capped-start rule E does not apply. Fibroblast = KB3781, the assembly's own
animal; testis = OR6737.)

## O1 — the family, de novo

- **Fibroblast assembly (KB3781), the families stage on the six contigs** (`tools/rustle_pipeline.sh families`, shipped defaults, bridge regroup
  off because the genome-wide assembly predates f1v2): **family MCL2 = exactly the six expressed copies** (LRPAP1, LOC134756753, copyB, copyC,
  LOC129526389, LOC134757218), size 6, density 1.000, no foreign member. The two copies without fibroblast reads (LOC129523574 0 reads,
  LOC129530227 1) have no locus — nothing to cluster.
- **Testis assembly (OR6737): family MCL9 = all eight copies**, density 0.893, no foreign member.
- Every de novo locus is represented by a 7- or 8-exon transcript (`DN_…_8`): the full LRPAP1 structure, unlike the NPIP fragments of
  `docs/SPLICED_COPY_SUPPORT_2026-10-04.md`. The dominant expressed chain at LRPAP1 is the RefSeq model XM_031007070.3's intron chain exactly
  (1,724 of 1,977 testis reads; 7 introns), and 98-99% of the reads at every expressed copy are that copy's dominant chain or a 5' piece of it.
- Guided mode (the annotation) can only know six of the eight (LRPAP1 + five LOCs); the two unannotated chr12 copies (copyB 24.79 Mb, copyC
  30.21 Mb) are found by the reads alone — the semi-guided case of `project_two_modes_scope` in one family. The legacy GWFAM catalog (378
  families) does not contain LRPAP1.

## O2 — copy assignment

- `copy_assign --families` on MCL2 (fibroblast) and MCL9 (testis): **the AS-tied gate admits 12 molecules in the fibroblast run and 8 in the
  testis run — none of them a read whose primary alignment lies in a copy** (all are secondary-only visitors; 9 contested, all `ambiguous`; 8
  origin-rejected). At 1.0-4.3% divergence between copies the aligner places essentially every read uniquely (LRPAP1's tie rate in the
  Fig. 3 table: 41 of 2,027 reads), and the family carries 4,257 (fibroblast) / 1,507 (testis) PSV columns. There is nothing for O2 to
  arbitrate here; this is the regime where the copies are distinguishable by sequence and O2 correctly stays out of the way. Expression per
  copy is therefore the primary read count: LRPAP1 and LOC134756753 carry 90% of the reads in both tissues; LOC129526389 and copyC are
  expressed at 5-8%; copyB, LOC134757218 at ~1%; the chr22 and chrY pseudogene copies only in testis (45 and 18 reads).

## O3 — are all copies in the reference?

- Aligning the eight primary-assembly copies to the animal's own haplotype assemblies (`gorilla_haps/{mat,pat}.fa`, minimap2 asm20, hits
  >= 95% identity over >= 80% of the 20-kb copy): **maternal haplotype 7 loci (chr3, chr12 x3, chr14, chr16, chr22), paternal 8 (the same plus
  chrY)**. The primary assembly is the paternal haplotype at seven of the eight (identity 1.000 at the same coordinates) and the maternal one on
  chr22 (NC_073246.2 is a maternal contig). Every primary copy has a counterpart on the other haplotype at 0.988-0.999 identity; no haplotype
  carries an LRPAP1 locus the primary lacks. **Diploid copy number 15 (7 + 8); no reference-absent copy** — the honest O3 statement is "none
  detected, and the matched haplotypes agree". (The WGS trio dosage table does not cover this family: it is not a GWFAM family.)

## What the family is good for

- A clean positive example for O1 de novo: a gorilla-specific expansion to eight full-length copies, recovered completely from reads in the
  testis library and completely among the expressed copies in fibroblasts, with two copies the annotation does not have.
- A clean negative example for O2: divergence 1-4% ⇒ no AS ties ⇒ the aligner is the assignment; O2 reports it as such instead of inventing
  decisions. Contrast NPIP (59% tied).
- A clean negative control for O3: haplotype-resolved copy numbers agree with the reference.
- Open: the Y-linked copy (LOC129530227, 2% diverged, testis-expressed at 18 reads) and the chr22 copy (4.3%, oldest) are the pseudogene-annotated
  members that are nevertheless transcribed in testis; whether their chains are full length (8 exons at both loci in the testis assembly) says
  they are not degraded transcripts. Dating the duplications from the identities (1.0-4.3%) is a separate exercise.
