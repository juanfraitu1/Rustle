# CHM13 CAT/Liftoff vs NCBI RefSeq (native CHM13): why they differ (2026-11-02, descriptive, no decision rule)

Script `bench/annot_compare/cat_vs_refseq.py` (inputs in its usage line; T2T-unique = UCSC hs1 `hgUnique.hg38`, SD = CHM13 SD90).
Location = >= 50% of the gene span in T2T-unique, else in SD90, else other. "Our copies" counts A119b family copies
(`a1.fam.copies.tsv`), which exist only for multi-copy families, so 0 is expected at single-copy loci.

Reading: (1) genes missing from CAT are almost all a CATALOGUE difference (RefSeq LOC / Gnomon models with no Ensembl gene), not
projection failure: 204 of 214 missing coding genes; projection-related coding losses are ~10 (SD: NPIPB3, FAM86B1, OR2G6; GSTT1).
(2) CAT merges genes that SHARE EXONS: the largest merges are alternative-promoter clusters with shared constant exons (PCDHG 22,
PCDHA 15, UGT1A 9 source genes in one CAT gene), which RefSeq keeps as separate genes. (3) Longer CAT genes are mostly long introns
from extra isoforms ("stretched intron" heuristic, 711) more than merges (35); RefSeq is the longer one in 499 pairs, so length
differences go both ways.

## RefSeq genes with no CAT gene at their exons (same strand)
| cause | location | coding | non-coding |
|---|---|---|---|
| catalogue (no Ensembl counterpart) | T2T-unique | 45 | 1098 |
| catalogue (no Ensembl counterpart) | SD90 | 36 | 381 |
| catalogue (no Ensembl counterpart) | other | 123 | 6439 |
| projected elsewhere | T2T-unique | 0 | 2 |
| projected elsewhere | SD90 | 4 | 12 |
| projected elsewhere | other | 1 | 16 |
| not projected | T2T-unique | 1 | 4 |
| not projected | SD90 | 0 | 8 |
| not projected | other | 4 | 43 |
RefSeq genes 44904 (coding 20090); missing in CAT 8217 (coding 214)
  e.g. catalogue (no Ensembl counterpart) / SD90: LOC124903856 LOC124903857 LOC124904581 LOC124904409 LOC124905947 LOC124906326 LOC124900997 LOC124906507
  e.g. catalogue (no Ensembl counterpart) / other: LOC124905687 LOC124903828 LOC124903840 LOC124905703 LOC107984936 LOC107985043 LOC124903813 LOC124905421
  e.g. catalogue (no Ensembl counterpart) / T2T-unique: LOC124900386 LOC124905693 LOC124905556 LOC128966598 LOC124905938 LOC124905939 LOC124906401 LOC124906629
  e.g. not projected / SD90: 
  e.g. not projected / other: PYDC5 RNF228 LY6S PANO1
  e.g. projected elsewhere / other: SPRY3
  e.g. projected elsewhere / SD90: OR2G6 FAM86B1 NPIPB3 QTGAL
  e.g. not projected / T2T-unique: GSTT1
  e.g. projected elsewhere / T2T-unique: 

## CAT genes with no RefSeq gene at their exons
| biotype | T2T-unique | SD90 | other |
|---|---|---|---|
| protein_coding | 28 | 23 | 30 |
| lncRNA | 181 | 671 | 7157 |
| other | 1074 | 3718 | 14698 |
CAT genes 64213; missing in RefSeq 27580

## Same gene, CAT span >= 1.5x and >= 10 kb longer than RefSeq (coding RefSeq genes)
| merged | 35 | MICOS10(145k vs 32k) PLPP3(150k vs 84k) TRIM5(276k vs 95k) PRB1(43k vs 3k) MIP(19k vs 7k) CLN5(28k vs 13k) TMEM63C(142k vs 77k) NOX5(132k vs 48k) |
| stretched intron | 711 | KLHL21(23k vs 12k) NBPF1(235k vs 71k) CROCC(320k vs 59k) NBL1(61k vs 15k) CDC42(75k vs 48k) ZBTB8OS(51k vs 31k) AK2(73k vs 28k) SMIM12(147k vs 9k) |
| isoform extent | 271 | ISG15(13k vs 1k) FBXO2(18k vs 6k) CELA2B(25k vs 15k) HMGCL(36k vs 23k) TSSK3(12k vs 1k) TMEM269(31k vs 15k) TMEM53(39k vs 21k) ZSWIM5(289k vs 190k) |
CAT longer: 1017; RefSeq longer (reverse): 499

## CAT collapsed genes: how many genes RefSeq and our assembly have at the same exons
paralog merges: 136 CAT genes standing for 370 source genes; RefSeq has >= 2 of the merged genes at 24 (18%); our assembly has >= 2 copies at 2 (1%)
all collapsed: 1667 CAT genes standing for 4210 source genes; RefSeq has >= 2 of the merged genes at 360 (22%); our assembly has >= 2 copies at 61 (4%)
| CAT gene | locus | source genes merged | RefSeq genes (named among merged) | our copies |
|---|---|---|---|---|
| PCDHGB4 | chr5:141913341-142038959 | 22 | 22 (22) | 0 |
| PCDHAC1 | chr5:141451597-141537634 | 15 | 16 (15) | 0 |
| U2 | chr17:44116164-44116355 | 11 | 1 (0) | 0 |
| UGT1A1 | chr2:234248900-234261929 | 9 | 9 (9) | 0 |
| ZNF17 | chr19:60508195-60518972 | 6 | 1 (1) | 0 |
| Y_RNA | chr13:5626961-5627062 | 5 | 1 (0) | 0 |
| KRTAP10-7 | chr21:42961054-42962632 | 4 | 1 (1) | 0 |
| FP236383.1 | chr13:5774166-5774256 | 4 | 0 (0) | 0 |
| MIR6724-2 | chr13:9245847-9245939 | 4 | 0 (0) | 0 |
| OR4F21 | chr1:111902-112896 | 3 | 1 (0) | 0 |
| CORT | chr1:9993647-9995614 | 3 | 2 (2) | 0 |
| EFNA3 | chr1:154218210-154226913 | 3 | 2 (1) | 0 |
| NOS1AP | chr1:161414003-161714449 | 3 | 1 (1) | 0 |
| SNORA70 | chr1:201789106-201789224 | 3 | 1 (0) | 0 |
| FXYD6 | chr11:117853383-117893879 | 3 | 2 (2) | 0 |
| HOXC5 | chr12:53998417-54000923 | 3 | 1 (1) | 0 |
| ERCC5 | chr13:102064026-102094203 | 3 | 2 (2) | 0 |
| SYNJ2BP | chr14:64572380-64623200 | 3 | 2 (2) | 0 |
| BDKRB2 | chr14:90435927-90475415 | 3 | 1 (1) | 0 |
| MIR3179-2 | chr16:16305715-16305799 | 3 | 2 (2) | 0 |
| MIR3670-2 | chr16:16311926-16311991 | 3 | 2 (2) | 0 |
| CMTM1 | chr16:72360687-72373432 | 3 | 2 (2) | 0 |
| DDX19B | chr16:76100651-76146326 | 3 | 1 (1) | 0 |
| RNU6-355P | chr16:96232696-96232800 | 3 | 0 (0) | 0 |
| ZNF878 | chr19:12172011-12181178 | 3 | 1 (1) | 0 |

## Merge audit of our own loci (A119b, vs coding RefSeq genes; `bench/annot_compare/merge_audit.py`)
| GTF | loci | spanning >= 2 coding RefSeq genes | fused (genes share no exon) | shared-exon cluster | mixed |
|---|---|---|---|---|---|
| a1.gtf (assembly) | 76,549 | 1,285 | 691 | 509 | 85 |
| a1.families.gtf (families input) | 76,161 | 1,112 | 543 | 502 | 67 |
684 of the 691 fused loci contain a transcript that itself joins >= 2 genes (read-supported junctions, not grouping alone). Largest
fused loci are tandem arrays (chr4:9.13-9.39 Mb 41 LOC genes; chr5:17.46-17.57 Mb 31). Truth side: famtruth gives a CAT collapsed
gene ONE truth gene with all its ENSG ids, so pairs inside PCDHA/PCDHG/UGT1A are never scored; our assembly makes the same merge,
so the score cannot see it. Per-locus tables: /anvme/workspace/v138da13-rustle/annot/merge_audit.a1{,fam}.tsv.
