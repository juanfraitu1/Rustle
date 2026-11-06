# PREREG: annotated-only assign (row 1551, 2026-11-06 21:30, user-proposed)

## Question
The assign cost is concentrated in de novo families with no annotated gene: of 1,565 A119b families, 1,353 have at
least one copy overlapping a GENCODE/CAT gene (DENOVO_CAT_GFF, gene features); the chr13 13-16.5 Mb hotspot
holds 104 families of which 12 are annotated, and the nine expensive panel families are all unannotated. User's
proposal: assign only annotated families until the pipeline is fast, then revisit. Output-changing by design.

## Rule (opt-in `RUSTLE_ANNOTATED_ONLY=GFF`; unset = unchanged)
tools/annotated_families.py writes the families with no copy overlapping any `gene` feature of GFF to
PREFIX.assign.unannotated.txt; the assign stage adds them to the existing `--skip-families` list (the depth-skip
path: their regions are not swept, their copies stay in the table). Genome-wide: 212 families skipped.

## Bars (gw16: singlenode, gw7 recipe, bin13, RUSTLE_DISPATCH_DEFER=1, staged BAM, wall limit 2 h)
- assign wall <= 2,700 s (half of gw12's 6,228; the rest of the genome took ~2,100 s of wall in gw12).
- Kept families: every row of gw_off's six tables whose family is annotated appears byte-identical in gw16's
  tables, and gw16 has no row for a skipped family (the check is on the family-filtered tables). Rows that
  differ are reported by family; the pre-registered expectation is 0 differing rows, since --skip-families is
  the shipped depth-skip path, but reads shared with a skipped family can change (documented in the pipeline).
- RSS < 250 GB.
## Prediction
wall 1,800-2,600 s (PASS); chr13 detect < 300 s; 0 differing rows among kept families; the 212 skipped families
are listed in PREFIX.assign.unannotated.txt. Not a default: adoption is a catalogue decision the user makes.
