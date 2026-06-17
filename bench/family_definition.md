# Annotation-based RNA-level multi-copy gene family definition (completeness)

**Definition:** a multi-copy gene family = a maximal connected set of ANNOTATED genes whose representative transcripts pairwise share a POA contiguous exon core ≥ 0.13 at core identity ≥ 0.7. Roster = the 22,983 annotated gene models; grouping = the validated POA criterion; families = connected components. Annotated-only (extensible: the same criterion runs against de-novo/unannotated loci later to add unannotated copies).

- homology edges: 16,354 ; **families (components, size≥2): 1,337** ; genes in a family: 5,060
- component size distribution (top): [250, 114, 72, 47, 45, 45, 34, 32, 29, 29, 27, 26]

## Completeness — are KNOWN multi-copy families recovered?
- **curated textbook families: 8/11 recovered** (members land in one component) ; missed: CRYBG, DEFB*, SIGLEC*
- similarity-built universe families: **46/53 recovered** ; missed (sample): ASDURF, CASP8, CDPF1, CREB1, GCA, GPR39, LOC129529456

### curated families & their recovery
| family | members (annotated) | in one component? |
|---|---|---|
| APOBEC3 | 6 (APOBEC3B, APOBEC3C, APOBEC3D, APOBEC3F, APOBEC3G, APOBEC3H) | YES (5/6) |
| CRYBG | 3 (CRYBG1, CRYBG2, CRYBG3) | no |
| DAZ | 2 (DAZ1, DAZL) | YES (2/2) |
| DEFB* | 20 (DEFB1, DEFB108B, DEFB110, DEFB112, DEFB113, DEFB115…) | no |
| GGT | 6 (GGT1, GGT5, GGT6, GGT7, GGTLC1, GGTLC2) | YES (3/6) |
| MAGEA* | 5 (MAGEA1, MAGEA10, MAGEA12, MAGEA4, MAGEA9) | YES (3/5) |
| PRAMEF* | 4 (PRAMEF12, PRAMEF17, PRAMEF19, PRAMEF20) | YES (2/4) |
| RABL2 | 2 (RABL2A, RABL2B) | YES (2/2) |
| RFPL | 4 (RFPL1, RFPL2, RFPL3, RFPL4A) | YES (3/4) |
| SIGLEC* | 6 (SIGLEC1, SIGLEC10, SIGLEC12, SIGLEC15, SIGLEC5, SIGLECL1) | no |
| TAS2R* | 18 (TAS2R1, TAS2R10, TAS2R13, TAS2R14, TAS2R16, TAS2R20…) | YES (2/18) |

## Over-merge control (the precision side of 'complete')
- mega-components (size≥25, likely domain-hub chains, NOT single families): 14 (largest 250)
- these are the transitive-closure artifacts; a tighter clustering (mutual-core / community detection) would split them. Size-2/3 components are clean copy sets.

## Honest scope
- ANNOTATED genes only (this definition); unannotated/de-novo copies are added later by running the SAME POA criterion against de-novo loci (the cross-chrom discovery already does this).
- RNA-structural operational family (shared contiguous exon core), NOT the formal gene-tree/protein family — a parallel DNA tier would supply that (and the ground truth).
- universe is itself similarity-built (recovery there is partly internal consistency); the curated textbook set is the independent completeness check.
