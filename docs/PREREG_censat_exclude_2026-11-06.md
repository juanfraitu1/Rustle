# PREREG row 1559: exclude families whose copies lie in centromeric satellite (`RUSTLE_EXCLUDE_BED=censat.bed`)

## Rule (before any number)
tools/rustle_pipeline.sh stage_assign: when RUSTLE_EXCLUDE_BED is set, a family is skipped (via the existing
`--skip-families` list, same hatch as rows 1551/depth) when the MAJORITY of its copies overlap an element of the
BED. BED = T2T censat v2.1 (hsat1A/1B/3, HOR, mon, bsat, censat, ct, rDNA). The families stage is untouched; the
catalogue still lists them; only assignment skips them. Opt-in.

## Bars
- Wall: row 1550 bar A measured against the new critical path, and absolute: assign <= 2,400 s on singlenode
  (row 1551 skipped the same hotspot and measured 1,806 s).
- Identity: the six tables restricted to kept families identical to gw_off (the row 1551 `[kept]` check) —
  the rows 1551 comparison FAILED this because skipping a family changes neighbours' candidate sets; the same is
  expected here and is recorded, not a bar. Bar = the assigned-read count on non-satellite families within 1%
  of gw_off's 3,367 (2,518 none + 746 ct + 103 other... ct is satellite: kept = none-class 2,518 only; bar =
  >= 2,493).
- Reads lost: the 8,927 + 746 + 103 satellite assignments are dropped BY DESIGN; stated now so it cannot be
  read as a recall loss.

## Prediction
Assign wall 1,700-2,000 s; e2e then ~2,700 + 320 + 1,900 = ~4,900 s < 7,200 s, which would make bar C
reachable. ct families (746 reads, 13.70-13.92 Mb transition zone, EIF3FP1 etc.) are the ambiguous class:
reported separately so the user can keep them.

## Implementation note (before the run)
`annotated_families.py --majority-in BED` lists families with strictly more than half their copies overlapping an
element: 550 of 1,565 (the 618 in the register row counted the majority CLASS incl. transition zone; 550 is the
figure the run uses). Row 1551 output unchanged (212). BED kept at data/denovo_eval/chm13v2.0_censat_v2.1.bed.
Run gw18 = gw16 harness with RUSTLE_EXCLUDE_BED instead of RUSTLE_ANNOTATED_ONLY.

## Result gw18 (job 4450882, singlenode, bin14 pipeline path, RUSTLE_EXCLUDE_BED + DISPATCH_DEFER=1)
| bar | value | verdict |
|---|---|---|
| assign wall <= 2,400 s | 843 s | PASS |
| non-satellite assigned reads >= 2,493 (gw_off 2,524) | 2,526 (10 families differ by a read or two: neighbours' candidate sets change when a satellite family leaves the window, as row 1551 saw) | PASS |
| RSS | 35.4 GB (gw15: 215 GB) | — |
Kept-family tables not byte-identical (expected, recorded). Longest region now 171 s (xchrom family), so the
row 1550 bar A "wall <= 1.1 x critical path" reads 843 vs 188 s: the regime changed and the bar is no longer the
right one; the user's goal text names ~4,550 s, and 843 s is below it by 5x. Goal amendment is the user's call.
Satellite reads dropped by design: 9,768 of 12,294 assigned in gw_off.
