# Pre-registration: are O1's families real duplications? Copy pairs against independent segmental-duplication annotations (2026-10-21)

**User (2026-10-21):** after register 1367 (the multi-mapping loci are duplications, but that says nothing about O1's families):
"yes, please proceed with what you suggest". The suggestion: compare O1's family copies with SEDEF and the CHM13 SD track, and
count how many copy pairs are annotated duplication pairs. The concern behind it: the advisor knows NPIP and the other mainstays,
and might think all other "families" are artefacts.

Written before any overlap is computed.

## Catalogs (O1 outputs, no SD input)

- **human_A119b:** our pipeline's families stage on CHM13. 2,150 families, 9,907 copies, 7,774 with `max_family_identity` >= 0.90.
- **gorilla_KB3781:** the same stage on the GGO T2T. 326 families, 918 copies, 773 with identity >= 0.90.
- **gorilla fibro_cat:** the KB3781 catalog used by the O2 / O3 tests (`rna_allele/fibro_cat.copies.tsv`). 666 families, 3,981
  copies, no identity column.

The families are built from RNA-supported loci and their sequence similarity. The default stage reads no SD annotation (`sd_depth`
and `core_hull` are NA, and `families.key` has no arguments), so the test is not circular.

## SD annotations (independent of our code)

- **human:** the published CHM13 SD track at >= 90% identity (`CHM13_SDs.SD90.bed`; both sides of every row). The
  satellite-filtered version (`SD90.nosat`) is reported as sensitivity.
- **gorilla:** SEDEF on the GGO T2T (`winloci_data/GGO_sedef_final.bed`, rows with fracMatch, column 21, >= 0.90; both sides).
  Register 1278 chose SEDEF as gorilla's SD source.

## Test

- **SD-paired copy:** a copy is SD-paired in its family when one SD row has one side overlapping one of the copy's exons and the
  other side overlapping an exon of another copy of the same family (>= 1 bp each).
- **Bar** (the pipeline catalogs): "O1's families are real duplications" when, among copies with `max_family_identity` >= 0.90,
  >= 90% are SD-paired in their family, in each catalog.
- **Specificity null:** the same test after shuffling family labels among copies, keeping family sizes (seed 1). Expected near 0.
- **Reported (not bars):**
  - all copies, by identity bin (>= 0.98, 0.95-0.98, 0.90-0.95, < 0.90);
  - families in which every >= 0.90 copy is SD-paired;
  - fibro_cat over all its copies;
  - the mainstays by gene or description: NPIP, TBC1D3, GOLGA6 / 8, SMN, LRPAP1, TSPY, RBMY, DAZ, HSFY, BPY2, CDY;
  - the copies that fail, by class: chrY, chrX, de novo copy (no annotated gene), and other.

## What it can and cannot show

**It shows** whether the copies O1 groups are duplications that an independent DNA method also pairs. That is the evidence that
the families are not artefacts of our grouping.

**It cannot show:**
- that every copy is expressed (O2);
- that families below 90% identity are real (SD90 cannot pair them; they are reported, not barred);
- anything about samples other than the two references.

## Rules

- Descriptive with one bar; nothing tuned.
- Light jobs only.
- Code: `bench/o1_sd/sd_pairs.py`, with tests.
- The metric-auditor checks the result before it is reported. Never pool human and gorilla. Push only when asked.

## Results (tool ba7ab571; run once; `/mnt/linuxdisk/tmp/o1_sd/sd_pairs.json`): **bar FAILS in both pipeline catalogs**

| catalog | copies >= 0.90 identity | SD-paired in their family | shuffled-label null | overlap no SD at all | families with every >= 0.90 copy paired |
|---|---|---|---|---|---|
| human_A119b (CHM13 SD90) | 7,774 | **5,605 (72.1%)** | 1,092 (14.0%) | 1,932 (24.9%) | 1,186 / 1,755 |
| gorilla_KB3781 (GGO SEDEF >= 0.90) | 773 | **637 (82.4%)** | 31 (4.0%) | 93 (12.0%) | 214 / 277 |
| gorilla fibro_cat, all copies (no identity column) | 3,981 | 3,557 (89.3%) | 763 (19.2%) | – | 543 / 666 (every copy) |

- **By identity bin, human:**
  - >= 0.98: 74.5% (null 16.7%);
  - 0.95-0.98: 80.0%;
  - 0.90-0.95: 51.2%;
  - < 0.90: 17.1%.
- **By identity bin, gorilla:** >= 0.98 82.2%; 0.95-0.98 88.1%; 0.90-0.95 76.6%; < 0.90 50.0%.
- **Satellite-filtered human SDs:** 52.5% of >= 0.90 copies are paired, so 27% of the human pairings rest on satellite SDs.
- **Correction to the design:** every copy in the two pipeline catalogs is a de novo copy (built from reads), so the registered
  failure class "de novo copy" says nothing.
- **What the unpaired copies are** (descriptive):
  - human: 2,169 unpaired; 1,932 overlap no SD; 1,191 multi-exon and 741 single-exon; 989 lie in an annotated gene, 943 do not;
  - gorilla: 136 unpaired, 93 with no SD.
  
  Their nature (retrocopies, short shared segments, repeat-driven groupings) is not tested here.
- **Mainstays** (descriptive; matched by overlap with annotated genes, after fixing a chromosome-naming mismatch in the first pass):

  | catalog | copies SD-paired in their family |
  |---|---|
  | human_A119b | NPIP 46 / 46 (12 families); GOLGA8 24 / 24; GOLGA6 24 / 25; TBC1D3 22 / 24; SMN1 2 / 2; SMN2 4 / 4; DAZ 6 / 6; TSPY 7 / 7; HSFY 8 / 9; RBMY 5 / 6; CDY 4 / 5; BPY2 3 / 3 (all: 155 / 161, 96%) |
  | gorilla fibro_cat (curated families) | NPIP 25 / 25; LRPAP1 (SM22) 22 / 22; TBC1D3 11 / 14 |
  | gorilla_KB3781 | only 2 mainstay copies matched (NPIP 1, LRPAP1 1), both paired |

**Reading.**
- The bar ("O1's families are real duplications": >= 90% of high-identity copies SD-paired) fails: 72.1% human, 82.4% gorilla.
- The pairing is specific: shuffled families pair 4-14%.
- The advisor's mainstays are 96% (human) and 58 / 61 (gorilla curated) SD-paired: those families are duplications by an
  independent DNA method.
- But a quarter of the human high-identity copies (an eighth in gorilla) lie in no annotated SD. Those groupings are not
  confirmed as duplications by this test, and are the place to look for artefacts.

Audit pending.
