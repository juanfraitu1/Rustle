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
