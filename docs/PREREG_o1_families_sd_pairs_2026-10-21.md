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
