# PREREG row 1563: row 1562's rule against repeat-class truth, plus the assign wall it buys

Written after row 1562's numbers (recall 0.171 against a truth that was 71% centromere-TRANSITION families with
no repeat content) and before any number of THIS row.

## Rule (unchanged from 1562): family flagged iff median within-copy 15-mer duplicated fraction >= 0.5.
## Truth: censat majority-in restricted to repeat classes hsat1A, hsat1B, hsat2, hsat3, hor, bsat, mon, rDNA,
censat (ct excluded). Computed by bench/satseq.py with a class-filtered id list.
## Bars
- precision >= 0.9, recall >= 0.8 on that truth (the 1562 post-hoc figures, 94/103 and 90/110 on the three
  tandem classes, are NOT this row's numbers: bsat and rDNA count against recall here).
- Assign wall with `RUSTLE_EXCLUDE_SELFREPEAT=1` (the flag list fed to the same skip hatch) <= 1,200 s on
  singlenode, non-repeat assigned reads within 1% of gw_off's count on the same families.
## Prediction
recall 0.75-0.85 (bsat/rDNA copies are 0/23 at dup15 >= 0.5: the bar is at risk from those two classes),
precision ~0.9; wall 850-1,100 s (the hsat1A hotspot is 77/94 caught; the missed 17 may hold a long family).
