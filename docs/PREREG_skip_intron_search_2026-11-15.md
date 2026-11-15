# PREREG row 1754 (2026-11-15): find a skipped exon of <= 15 bp by searching the intron for the read's extra bases

Row 1751 (register 1753) failed on reach: an exon dropped in every read has no junction in the locus to realign to
(reach 0.26-0.52 of row 1746's exons). This row searches the intron itself. Written and committed before the search was
run on any sample (one small gorilla contig was run for mechanics only, output not read). Rows 1746 and 1751 results on
all four samples were already seen, so no sample is unseen here; nothing in the rule is tuned on them.

## Rule (`bench/famcmp/skip_search.py`)

- A primary read's intron (d, a) with canonical ends and L = 3..15 extra read bases in the junction window (row 1751's
  window, 25 bp each side, clamped to the read's flanking exons).
- Sites: every L bp stretch of the intron with the intron's own splice signals around it (AG before, GT after, on the
  intron's strand) and >= 50 bp of intron on each side.
- The search is used only when sites x 4^-L <= 0.01 (expected chance matches per read). Short exons in long introns are
  therefore not searched.
- A site is taken only when it is the single one whose sequence is in the read, with <= 1 edit left in the window and
  none of the L extra bases unexplained.

## Step 0 (read level, four samples)

- **Gate (decoy):** the same search on the reversed intron. Reads placed by the decoy / reads placed by the search
  <= 0.05 on each sample, else INVALID.
- **Precision:** found units whose outer junctions (d, a) flank an annotated exon of <= 15 bp: the found exon is the
  annotated one in >= 0.95 on each sample with >= 20 such units.
- **Reach:** row 1746's exons (>= half the skip reads carry the exon; 26 / 27 / 34 / 27) reached by this search or
  by row 1751's rule: >= 0.75 on >= 3 of the 4 samples. The search alone and the lengths of the exons not reached are
  reported.
- Reported, no bar: found units by annotation; "neither" = unannotated short exons.

## Step 1 (assembly; read only if step 0 passes all three)

As row 1751: panel = reads within 100 kb of a fired unit, arms `base` and `fix` (`skip_fix.py fix` on the units of both
rules), default assemble, same binary. S1 skip chains at exon-only annotated fired units fall by >= 50%; S3 panel
precision of `fix` >= `base` - 0.0005. S2 is changed from row 1751 on what its losses turned out to be (all but 2 of
15 lost correct chains were annotated skip forms at units where both forms are annotated, their reads moved to the exon
form): lost correct chains that are not such skip forms <= new correct chains. PASS = S1-S3 on >= 3 of 4 samples. A PASS
makes the two rules together a candidate for an opt-in stage (user decision). No default changes here.

## Limits

Exons of 3-5 bp are almost never searchable (the chance bound), so reach is expected to miss them. Only an exact (d, a)
skip intron with the bases inserted next to it is reached. The rewrite's read test (2-edit margin) is looser than the
search, so it also moves reads of the same intron that carry the exon with an error.
