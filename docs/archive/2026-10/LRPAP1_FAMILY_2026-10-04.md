# The LRPAP1 family in gorilla: eight copies, one annotated name — what O1, O2 and O3 say about it (2026-10-04, for the advisor)

LRPAP1 (LDL receptor related protein associated protein 1; RAP) is single-copy in human (CHM13 chr4:3,503,365-3,532,077; 1,717 A119b reads),
chimpanzee and orangutan (Liftoff self-lift: `in_place` only). **In gorilla (mGorGor1 primary assembly) its 20-kb, intron-containing body is
present at eight loci** on six chromosomes: the Liftoff self-lift of the gorilla annotation (`figures/_liftoff.py`, Fig. 8) projects LRPAP1 at
0.957-0.990 sequence identity with 0.978-1.000 coverage to seven other sites, all seven of which RefSeq annotates as "alpha-2-macroglobulin
receptor-associated protein-like" LOC genes (four protein-coding, three pseudogenes; see the correction below — the first version of this page called two of them unannotated). These are
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


(chain-support reads = Amendment C of `docs/archive/2026-10/PREREG_spliced_copy_support_2026-10-04.md`: reads whose junction chain is an expressed chain of
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
  `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md`. The dominant expressed chain at LRPAP1 is the RefSeq model XM_031007070.3's intron chain exactly
  (1,724 of 1,977 testis reads; 7 introns), and 98-99% of the reads at every expressed copy are that copy's dominant chain or a 5' piece of it.
- Guided mode (the annotation) knows all eleven loci by name (correction below); the de novo mode finds the same eleven from reads alone. The legacy GWFAM catalog (378
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

## How many copies are found in gorilla, by each definition

| how the copies are counted | copies in gorilla |
|---|---|
| RefSeq annotation (mGorGor1): named loci | **11** — LRPAP1 + 10 "alpha-2-macroglobulin receptor-associated protein-like" LOCs (6 protein-coding, 4 pseudogene); the first version of this table said 6: the copies' RefSeq name was missed (row 1246) |
| Sequence: Liftoff self-lift of LRPAP1's body (>= 95.7% identity, >= 98% coverage) | **8** full-length; a genome-wide minimap2 search adds **3 partial copies** (the 5' three-quarters of the body) = **11** |
| Haplotype assemblies of the same animal (>= 95% identity over >= 80% of the copy) | maternal **7**, paternal **8** (the Y copy); diploid 15 |
| Expressed: >= 3 uniquely placed reads with one identical >= 2-intron chain (Amendment C) | full-length: fibroblast **6**, testis **7**, either **8**; partial: **3 / 3 / 3** → **11 of 11** in at least one library |
| A de novo locus built at the copy | fibroblast **9** (6 + 3), testis **11** (8 + 3) |
| A member of a de novo LRPAP1 family | fibroblast **9** in two families (MCL2 = the 6 full-length, MCL11 = the 3 fragments), testis **11** in two (MCL9 = 8, MCL44 = 3) |
| FOUND, strict: the locus representative is an expressed chain of the copy or a 5' piece of it (C) | fibroblast **9 of 9** expressed, testis **10 of 10**, in at least one library **11 of 11** |
| FOUND, TSS-anchored: the representative starts at the copy's modal start with its first three introns (D') | fibroblast **8**, testis **9** — LRPAP1's (fibroblast) and LOC134756753's (testis) representatives carry 6 of the 7 introns and miss the first; all three fragments pass |
| In the legacy GWFAM catalog | **0** (the family is not in it) |
| The advisor's reference (SD-based, Supplementary Table VIII.38) | **10** = 1 ancestral + 9; our 11 minus the chr16 5'-fragment at 17.18 Mb |

(`found_fibro.*`, `found_testis.*` in the work dir: `bench/copy_support.py` with each library's de novo loci as the `--loci` set.)


## Correction and the paper's count (2026-10-04 16:10): eleven LRPAP1-like loci, not eight — and all of them annotated

The advisor's reference reports **10 copies** in gorilla (1 ancestral on chr3 + 9: four on chr12, two on chr14, one each on chr16, chr22 and
Y; copies 2 (chr12) and 6 (chr14) "solitary", without the flanking DOK7/HGFAC). What this dossier missed, and why:

1. **Three partial copies.** A genome-wide search with the 20-kb LRPAP1 body (minimap2 asm20 against the primary assembly) finds **11 loci**:
   the eight full-length ones above (coverage 1.0, identity 0.969-0.977 to LRPAP1 at asm20 scoring) **and three copies of the 5' three-quarters
   of the body** (body positions 4,971-20,001 of 20,001 = 15 kb carrying exons 1-5, lacking the last four exons; identity 0.959-0.969): chr12
   NC_073236.2:23,071,592-23,086,720 (RefSeq LOC134756368), chr14 NC_073238.2:25,252,147-25,267,314 (LOC115932954), chr16
   NC_073240.2:17,177,201-17,192,472 (LOC115932756) — all three protein-coding RefSeq models of 5-6 exons, all three "solitary" (no HGFAC-like /
   DOK7 neighbour; every full-length copy has an HGFAC-like gene 60-90 kb away). The Liftoff self-lift that this dossier started from reports a
   copy only when the whole gene lifts (coverage ~0.98+), so the 5'-fragment copies never entered the list.
2. **The two "unannotated" copies were annotated.** RefSeq names every full-length copy: chr12 24.79 Mb = LOC115933156 (pseudogene), chr12 30.21
   Mb = LOC129523503 (protein-coding). The first grep of this dossier matched the description "LDL receptor related protein associated protein
   1" and missed RefSeq's name for the copies, "alpha-2-macroglobulin receptor-associated protein-like"; register row 1245 is corrected by row
   1246. RefSeq's own count: LRPAP1 + 10 "…-like" LOCs = 11.
3. **Mapping to the paper's ten:** chr3 LRPAP1 = ancestral; chr12 x4 = our 22.55 (LOC134756753), 23.07 (partial, "solitary" = their copy 2), 24.79
   (LOC115933156), 30.21 (LOC129523503); chr14 x2 = 23.39 (LOC129526389) and 25.25 (partial, "solitary" = their copy 6); chr16 = 15.82
   (LOC134757218); chr22 = LOC129523574; Y = LOC129530227. **Our eleventh, chr16 NC_073240.2:17.18 Mb (LOC115932756, a third solitary 5'
   fragment, 40 / 218 reads), is not in the paper's list**; their protein-identity tiers (98-99% for copies 1, 4, 7; 84-91% for 5, 8, 9; 65-79% for
   2, 3, 6) are consistent with the 5' fragments being the low-identity, solitary members.
4. **The partial copies are expressed and found.** Reads fibroblast / testis: 83 / 158, 168 / 272, 40 / 217; dominant chains of 4 introns
   (FSM of their RefSeq models at chr12 and chr14, a novel splice site at chr16); 94-98% of their reads are chain support; de novo loci of 5
   exons at all three in both assemblies; FOUND under C and D' in both libraries (6 of 6).
5. **Our O1 puts them in a second family.** In both assemblies the three fragments form their own de novo family (MCL11 fibroblast, MCL44
   testis; size 3) beside the full-length family (MCL2 size 6, MCL9 size 8): MCL partitions the eleven into "full-length" and "5'-fragment"
   clusters — the fragments are more alike (one breakpoint, identity 0.995-0.999 among themselves) than they are to the full copies (0.96-0.97).
   Against the paper's one family of ten this is a split, the same partition-vs-cover question as NPIP's subfamilies (register 1236): at the
   family level the count is 11 in two clusters; at the paper's level one family of 10 (11 with the chr16 fragment).
6. Haplotypes of the fragments: paternal at the same coordinates (identity 1.000); maternal counterparts at 0.9988 (chr14) and 0.9991 (chr16),
   the chr12 fragment's maternal counterpart not resolved by this crude test (its best hit is the chr14 one at 0.9956).

### All eleven copies

| copy | kind | gorilla chr | span | identity | RefSeq | reads f / t | de novo locus f / t | chain-support f / t | TSS-anchored f / t |
|---|---|---|---|---|---|---|---|---|---|
| LRPAP1 | full-length (20 kb, 8 exons) | chr3 (hsa4) | NC_073227.2:12,090,719-12,110,719 (-) | 1.000 | RefSeq gene LRPAP1 (protein-coding) | 1,488 / 1,977 | DN_NC_073227.2_12086114_8 (8 exons, 206 reads) / DN_NC_073227.2_12088225_8 (8 exons, 329 reads) | 1471 / 1932 | 1468 / 1943 |
| LOC134756753 | full-length (20 kb, 8 exons) | chr12 (hsa2a) | NC_073236.2:22,553,726-22,573,741 (-) | 0.990 | RefSeq LOC134756753, protein-coding, "alpha-2-macroglobulin receptor-associated protein-like" | 836 / 1,546 | DN_NC_073236.2_22551320_8 (8 exons, 142 reads) / DN_NC_073236.2_22549336_8 (8 exons, 257 reads) | 823 / 1524 | 821 / 1528 |
| LOC115933156 | full-length (20 kb, 8 exons) | chr12 (hsa2a) | NC_073236.2:24,785,460-24,805,469 (+) | 0.985 | RefSeq LOC115933156, pseudogene, "…-like" (not reported by the Liftoff pairs: a pseudogene) | 19 / 22 | DN_NC_073236.2_24785601_8 (8 exons, 15 reads) / DN_NC_073236.2_24785618_8 (8 exons, 9 reads) | 15 / 19 | 18 / 19 |
| LOC129523503 | full-length (20 kb, 8 exons) | chr12 (hsa2a) | NC_073236.2:30,205,371-30,225,550 (+) | 0.988 | RefSeq LOC129523503, protein-coding, "…-like" | 57 / 121 | DN_NC_073236.2_30205515_7 (7 exons, 23 reads) / DN_NC_073236.2_30205531_7 (7 exons, 38 reads) | 51 / 113 | 50 / 119 |
| LOC129526389 | full-length (20 kb, 8 exons) | chr14 (hsa13) | NC_073238.2:23,392,677-23,412,705 (+) | 0.990 | RefSeq LOC129526389, protein-coding, "…-like" | 105 / 84 | DN_NC_073238.2_23392840_7 (7 exons, 16 reads) / DN_NC_073238.2_23392840_7 (7 exons, 28 reads) | 97 / 81 | 97 / 81 |
| LOC134757218 | full-length (20 kb, 8 exons) | chr16 (hsa15) | NC_073240.2:15,822,038-15,842,048 (-) | 0.986 | RefSeq LOC134757218, pseudogene | 28 / 1 | DN_NC_073240.2_15821047_8 (8 exons, 15 reads) / DN_NC_073240.2_15820894_8 (8 exons, 12 reads) | 18 / 0 | 24 / 0 |
| LOC129523574 | full-length (20 kb, 8 exons) | chr22 (hsa21) | NC_073246.2:11,292,022-11,312,055 (-) | 0.957 | RefSeq LOC129523574, pseudogene | 0 / 45 | no locus (no reads) / DN_NC_073246.2_11287585_8 (8 exons, 16 reads) | 0 / 45 | 0 / 45 |
| LOC129530227 | full-length (20 kb, 8 exons) | chrY | NC_073248.2:45,277,394-45,297,429 (+) | 0.980 | RefSeq LOC129530227, pseudogene (chrY) | 1 / 18 | no locus (no reads) / DN_NC_073248.2_45277558_8 (8 exons, 6 reads) | 0 / 14 | 0 / 14 |
| LOC134756368 | partial (5′ 15 kb, 5 exons) | chr12 (hsa2a) | NC_073236.2:23,070,423-23,086,604 (-) | 0.969 | RefSeq LOC134756368, protein-coding, "…-like"; 5′ three-quarters of the body (15 kb, exons 1–5), no HGFAC/DOK7 neighbour ("solitary") | 83 / 158 | DN_NC_073236.2_23068850_5 (5 exons) / DN_NC_073236.2_23070390_5 (5 exons) | 79 / 152 | 78 / 154 |
| LOC115932954 | partial (5′ 15 kb, 5 exons) | chr14 (hsa13) | NC_073238.2:25,250,978-25,267,164 (-) | 0.966 | RefSeq LOC115932954, protein-coding, "…-like"; 5′ three-quarters of the body (15 kb, exons 1–5), no HGFAC/DOK7 neighbour ("solitary") | 168 / 272 | DN_NC_073238.2_25249406_5 (5 exons) / DN_NC_073238.2_25250945_5 (5 exons) | 160 / 263 | 159 / 262 |
| LOC115932756 | partial (5′ 15 kb, 5 exons) | chr16 (hsa15) | NC_073240.2:17,176,045-17,192,323 (-) | 0.959 | RefSeq LOC115932756, protein-coding, "…-like"; 5′ three-quarters of the body (15 kb, exons 1–5), no HGFAC/DOK7 neighbour ("solitary") | 40 / 217 | DN_NC_073240.2_17176025_5 (5 exons) / DN_NC_073240.2_17174470_5 (5 exons) | 37 / 203 | 35 / 206 |


## What the family is good for

- A clean positive example for O1 de novo: a gorilla-specific expansion to eight full-length copies, recovered completely from reads in the
  testis library and completely among the expressed copies in fibroblasts, with two copies the annotation does not have.
- A clean negative example for O2: divergence 1-4% ⇒ no AS ties ⇒ the aligner is the assignment; O2 reports it as such instead of inventing
  decisions. Contrast NPIP (59% tied).
- A clean negative control for O3: haplotype-resolved copy numbers agree with the reference.
- Open: the Y-linked copy (LOC129530227, 2% diverged, testis-expressed at 18 reads) and the chr22 copy (4.3%, oldest) are the pseudogene-annotated
  members that are nevertheless transcribed in testis; whether their chains are full length (8 exons at both loci in the testis assembly) says
  they are not degraded transcripts. Dating the duplications from the identities (1.0-4.3%) is a separate exercise.
