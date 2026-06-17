# Two-tier overlay: RNA layer on DNA/protein-defined multi-copy families

**DNA tier** = protein clusters (mmseqs2 easy-cluster, ≥30% id / 50% cov, on translated CDS) →
formal multi-copy families, INCLUDING ancient/diverged families the RNA-similarity definition
missed. **RNA overlay** = per copy, is it transcribed (real GGO IsoSeq), and how well.

- DNA multi-copy families (protein clusters ≥2): **3,587** (14,545 gene copies)
- RNA-tier families (POA, for comparison): 1,337

## What actually transcribes (genome vs transcriptome)
- of 14,545 copies in DNA-defined families: **transcribed (≥5 reads): 10,490 (72%)** ; well-expressed (≥40): 6,824 (47%) ; **silent: 4,055 (28%)**
- i.e. the DNA tier enumerates every copy; the RNA layer shows which are live — many copies are
  present in the genome but transcriptionally silent.

## Ancient-family gain (the whole point of the DNA tier)
- curated families recovered by the DNA tier that the RNA tier MISSED: **DEFB*, SIGLEC***
| family | DNA tier | RNA tier | DNA copies | transcribed |
|---|---|---|---|---|
| APOBEC3 | YES | YES | 4 | 4 |
| CRYBG (ANCIENT) | no | no | 0 | 0 |
| DAZ | no | YES | 1 | 1 |
| DEFB* (ANCIENT) | YES | no | 4 | 0 |
| MAGEA* | YES | YES | 5 | 5 |
| PRAMEF* | YES | YES | 4 | 0 |
| RABL2 | YES | YES | 2 | 2 |
| RFPL | YES | YES | 4 | 3 |
| SIGLEC* (ANCIENT) | YES | no | 3 | 2 |
| TAS2R* (ANCIENT) | YES | YES | 18 | 7 |

## Per-family 3-number summary (sample: largest + curated)
| DNA family | copies | transcribed | well-expressed | example members |
|---|---|---|---|---|
| DFAM0 | 501 | 83 | 4 | LOC101123691, LOC101123789, LOC101123793, LOC101124039, LOC101124044… |
| DFAM1 | 229 | 219 | 180 | KRBOX5, LOC101123988, LOC101124084, LOC101124732, LOC101124778… |
| DFAM2 | 136 | 134 | 108 | LOC101126065, LOC101126415, LOC101127631, LOC101128844, LOC101129578… |
| DFAM3 | 93 | 82 | 50 | CDC42, IFT27, LOC101128843, LOC101133567, LOC101142457… |
| DFAM4 | 74 | 19 | 8 | BFSP2, DES, GFAP, GHAA, INA… |
| DFAM5 | 64 | 31 | 2 | ACR, CELA1, CFD, CTRL, F10… |
| DFAM6 | 47 | 12 | 1 | LOC101124648, LOC101136188, LOC101144932, LOC101146296, LOC101147007… |
| DFAM7 | 47 | 0 | 0 | LHB, LOC101123748, LOC101124600, LOC101154318, LOC109024055… |
| DFAM8 | 46 | 42 | 31 | CDK1, CDK10, CDK14, CDK15, CDK16… |
| DFAM9 | 45 | 26 | 13 | ACKR2, ACKR4, AGTR1, APLNR, CCR10… |
| DFAM10 | 44 | 0 | 0 | LOC101123968, LOC101128508, LOC115932853, LOC115932855, LOC129523578… |
| DFAM11 | 43 | 17 | 0 | LOC101126783, LOC101127524, LOC129524344, LOC129524346, LOC129524347… |

## Honest scope
- DNA tier = protein clustering of annotated CODING genes → catches ancient coding families +
  currently-silent coding copies. Non-coding/pseudogene + UNANNOTATED copies need a genomic
  self-alignment pass (next extension).
- 'transcribed' uses the REAL IsoSeq (not the ideal-coverage synthetic) — so silent = silent in
  this sample; ideal-coverage GGO would test detectability, not biology.
- copy-resolvability (which transcribed copies are distinguishable per-read) follows from the
  identifiability theorem: dispersed copies resolve by locus; co-located need PSVs (rare here).
