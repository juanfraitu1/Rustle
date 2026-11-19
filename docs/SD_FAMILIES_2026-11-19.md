# Are segmental duplications where recent multi-copy gene families live? (2026-11-19, PSU, register 1827)

User's question: look for ways to correlate segmental duplications (SDs) with multi-copy gene families; can we affirm that SDs are
the most likely places to harbor recent multi-copy gene families? Descriptive. Script `bench/denovo_eval/sd_families.py`; workspace
`goal_multicopy_dev/sd` (linked from `r1827_sd_families`).

## Data
- Genes: 20,012 protein-coding CAT genes on CHM13 with an Ensembl ID that Compara knows (famtruth's judgeable universe).
- Age of a gene: the youngest duplication level among its Compara within-species paralogue pairs (`$DENOVO_PARALOGS`; Homo
  sapiens < Homininae < Hominidae < Hominoidea < Catarrhini < Simiiformes < ... < Opisthokonta); genes with no pair = no paralogue.
- SDs: SEDEF on CHM13 (`chm13v2.0.fa.sedef.bed`, 88,756 pairs, both copies; identity = fracMatch). A gene is in an SD when >= 50%
  of its exonic bp lie in SD copies; its SD identity = the highest among the SDs covering it.

## 1. Recent family members sit in SDs; the older the family, the less
| youngest duplication | genes | in an SD | in an SD >= 95% identity |
| human-specific | 675 | **97.0%** | 89.9% |
| Homininae / Hominidae | 248 | 92.3% | 45.2% |
| Hominoidea / Catarrhini | 320 | 81.6% | 20.6% |
| Simiiformes / Haplorrhini / Primates | 244 | 66.8% | 4.1% |
| Euarchontoglires .. Mammalia | 821 | 21.2% | 5.4% |
| Amniota .. Opisthokonta | 13,329 | 6.7% | 2.6% |
| no paralogue (baseline) | 4,375 | 9.1% | 3.7% |
Inverse: of the 1,353 genes in an SD of >= 95% identity, 785 (58%) have a recent (<= Catarrhini) paralogue.

## 2. The two copies of a recent pair are the two copies of one SD
A paralogue pair is joined when ONE SEDEF pair has one copy on each gene:
| pair's duplication level | pairs | joined by one SD | median SD identity |
| human-specific | 3,052 | **99.5%** | 0.996 |
| Homininae / Hominidae | 1,066 | 98.7% | 0.990 |
| Hominoidea / Catarrhini | 578 | 87.5% | 0.951 |
| Simiiformes .. Primates | 1,048 | 55.4% | 0.916 |
| Euarchontoglires .. Mammalia | 3,040 | 9.1% | 0.885 |
| older | 129,310 | 0.8% | 0.878 |
The SD's identity falls with the age of the duplication: the SD is the duplication event itself, eroding with time.

## 3. Not circular: recent copies were duplicated as genomic DNA
SEDEF calls any >= 1 kb, >= 90% identical copy an SD, whatever made it, so recent paralogues would be "in SDs" even if they were
retrocopies. The mechanism separates the two: a segmental (DNA-based) duplication copies introns and flanks; a retrocopy (RNA-based)
copies exons only. Joined pairs whose SD carries >= 50% of the gene's intronic bp (or >= 1 kb of flank for single-exon genes):
| level | joined pairs | genomic (DNA-based) | exonic-only (retrocopy-like) |
| human-specific | 3,036 | **99.0%** | 31 |
| Homininae / Hominidae | 1,052 | 97.1% | 31 |
| Hominoidea / Catarrhini | 506 | 85.2% | 75 |
| Simiiformes .. Primates | 581 | 63.9% | 210 |
| Euarchontoglires .. Mammalia | 278 | 40.3% | 166 |
| older | 984 | 38.6% | 604 |

## 4. The exceptions
98 of 1,243 recent-family genes (<= Catarrhini) are less than half covered by SDs (44 of them at the oldest recent level,
Catarrhini; 69 partly covered): small genes in tandem clusters below SEDEF's 1 kb / 90% thresholds (olfactory receptors 9, UGT 5,
KRTAP 8, LCE 3, SPRR 2) and retrocopies (CSNK2A3, MAP1LC3, UBE2, EBLN).

## Answer
Yes, with numbers: SDs are where recent multi-copy gene families live. 97% of genes with a human-specific paralogue and 92% with a
great-ape paralogue lie in SDs (baseline 9% for genes with no paralogue, ~10x enrichment); 99% of human-specific paralogue pairs are
the two copies of one SD, and 99% of those SDs copied the genomic locus (introns, flanks), not the mRNA. The link weakens
monotonically with the age of the duplication (SD identity 0.996 -> 0.878; pairs joined 99.5% -> 0.8%) as SDs erode below detection,
and is weakest for small tandem-cluster genes and retrocopies. For the thesis: SD-embedded families are the young subset of
multi-copy gene families (vault `Topics/Terminology.md`), and this measures how complete that subset is (about 92% of recent ones).

## Caveats
CHM13 only (one genome; SEDEF's own thresholds: >= 1 kb, >= 90%); Compara levels come from Ensembl paralogue calls mapped to CAT
genes through Ensembl IDs (genes without an Ensembl ID are out); "in an SD" is a 50% exonic-cover threshold; the pair test needs one
SEDEF pair to touch both genes, so pairs split across several SEDEF records count as not joined (a lower bound).
