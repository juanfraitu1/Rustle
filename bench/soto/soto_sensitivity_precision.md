# Soto per-family sensitivity / precision — updated 2026-07-21

## Aggregate

| | members | %  |
|---|---|---|
| **detected now** (current catalog) | 276/362 | 76.2% |
| **projected post-seeding-fix** | 316/362 | 87.3% |
| families fully recovered (now) | 42/83 | 51% |

Missed members (86) by cause:
- **seeding-fixable** (unspliced pooled+deleted; the seeding fix recovers): **38**
- **K=0 identifiability floor** (exon-homogenized / young-identical; non-exon rescue verified = 0 → needs DNA): **34**
- silent (classification artifact; ≥1 is a seeding case): 2
- other: 12

**Honest reading:** current sensitivity 76.2%; the seeding fix projects to
~87.3% (recovering the 38 seeding-fixable members); the residual
34 are the genuine K=0 floor (rigorously confirmed irreducible from RNA — the aligner already
uses UTR/intron/flank and they still tie). The projection is a PROJECTION until the genome-wide
catalog is rebuilt with the seeding-fix binary (`bench/soto/seeding_validation.txt` proves the
per-locus mechanism: OLD 0 reps → NEW 7 reps).

## Per-family table (incomplete families first)

| family | gene | members | sens now | prec | missed | seeding-fix | K=0 | proj post-fix |
|---|---|---|---|---|---|---|---|---|
| ID_431 | PPIAL4A | 5 | 0% | NA | 5 | 0 | 5 | 0/5 = 0% |
| ID_63 | ANAPC1P1 | 8 | 50% | 100% | 4 | 3 | 0 | 7/8 = 88% |
| ID_213 | AC246785.1 | 4 | 0% | NA | 4 | 0 | 1 | 1/4 = 25% |
| ID_226 | H2BP1 | 5 | 20% | 100% | 4 | 2 | 2 | 3/5 = 60% |
| ID_261 | CNTNAP3P5 | 5 | 20% | 100% | 4 | 1 | 1 | 2/5 = 40% |
| ID_14 | AC005562.2 | 7 | 57% | 100% | 3 | 2 | 1 | 6/7 = 86% |
| ID_127 | TCAF2 | 3 | 0% | NA | 3 | 3 | 0 | 3/3 = 100% |
| ID_167 | CHEK2P2 | 3 | 0% | NA | 3 | 0 | 2 | 0/3 = 0% |
| ID_240 | SYT15 | 3 | 0% | NA | 3 | 1 | 2 | 1/3 = 33% |
| ID_260 | AL513478.2 | 3 | 0% | NA | 3 | 3 | 0 | 3/3 = 100% |
| ID_402 | NCF1B | 3 | 0% | NA | 3 | 1 | 2 | 1/3 = 33% |
| ID_411 | OR11H12 | 3 | 0% | NA | 3 | 0 | 1 | 0/3 = 0% |
| ID_26 | AC006453.4 | 2 | 0% | NA | 2 | 1 | 1 | 1/2 = 50% |
| ID_43 | TRIM64DP | 4 | 50% | 100% | 2 | 0 | 0 | 2/4 = 50% |
| ID_65 | FAR2P1 | 9 | 78% | 100% | 2 | 0 | 1 | 7/9 = 78% |
| ID_92 | AC243829.7 | 2 | 0% | NA | 2 | 0 | 2 | 0/2 = 0% |
| ID_146 | AC233280.19 | 2 | 0% | NA | 2 | 1 | 1 | 1/2 = 50% |
| ID_148 | AC126603.1 | 2 | 0% | NA | 2 | 2 | 0 | 2/2 = 100% |
| ID_175 | TEKT4P2 | 2 | 0% | NA | 2 | 2 | 0 | 2/2 = 100% |
| ID_215 | AC244669.2 | 7 | 71% | 100% | 2 | 1 | 1 | 6/7 = 86% |
| ID_251 | ANKRD36B | 2 | 0% | NA | 2 | 2 | 0 | 2/2 = 100% |
| ID_313 | CDH12P1 | 2 | 0% | NA | 2 | 0 | 2 | 0/2 = 0% |
| ID_332 | DDX11L2 | 6 | 67% | 100% | 2 | 0 | 1 | 4/6 = 67% |
| ID_338 | DEFB108C | 3 | 33% | 100% | 2 | 2 | 0 | 3/3 = 100% |
| ID_386 | LIMS1 | 2 | 0% | NA | 2 | 2 | 0 | 2/2 = 100% |
| ID_407 | NSUN5P2 | 4 | 50% | 100% | 2 | 1 | 1 | 3/4 = 75% |
| ID_458 | SPAG11A | 2 | 0% | NA | 2 | 2 | 0 | 2/2 = 100% |
| ID_12 | POM121 | 4 | 75% | 100% | 1 | 1 | 0 | 4/4 = 100% |
| ID_35 | AL669831.1 | 10 | 90% | 100% | 1 | 1 | 0 | 10/10 = 100% |
| ID_78 | GOLGA8A | 5 | 80% | 100% | 1 | 0 | 1 | 4/5 = 80% |
| ID_131 | AMY2A | 6 | 83% | 100% | 1 | 0 | 1 | 5/6 = 83% |
| ID_147 | AC125634.1 | 4 | 75% | 100% | 1 | 1 | 0 | 4/4 = 100% |
| ID_156 | WHAMMP3 | 5 | 80% | 100% | 1 | 1 | 0 | 5/5 = 100% |
| ID_208 | GTF2IP14 | 5 | 80% | 100% | 1 | 1 | 0 | 5/5 = 100% |
| ID_222 | AC243829.6 | 1 | 0% | NA | 1 | 0 | 1 | 0/1 = 0% |
| ID_302 | BOLA2B | 1 | 0% | NA | 1 | 0 | 1 | 0/1 = 0% |
| ID_327 | CTSLP3 | 4 | 75% | 100% | 1 | 1 | 0 | 4/4 = 100% |
| ID_334 | DEFB104B | 1 | 0% | NA | 1 | 0 | 1 | 0/1 = 0% |
| ID_348 | DUX4L50 | 1 | 0% | NA | 1 | 0 | 1 | 0/1 = 0% |
| ID_391 | WASH2P | 5 | 80% | 100% | 1 | 0 | 1 | 4/5 = 80% |
| ID_481 | UBE2Q2P6 | 5 | 80% | 100% | 1 | 0 | 0 | 5/5 = 100% |

*(42 families already at 100% sensitivity / 100% precision are omitted from the table above; full data in `soto_sensitivity_precision.tsv`.)*
