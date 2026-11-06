# PREREG row 1558: release the read-budget estimate once the fetch count is known (`RUSTLE_BUDGET_SHRINK=1`)

## Evidence (gw15, row 1550 Result)
chr13's budget permit is its BAM-index ESTIMATE, 3.3M records, held for the whole region (4,282 s). The fetch
returns 1.42M records after 241 s. The 1.9M-record difference kept MCL3/5/8/10 (500-580k each, 460-495 s each)
and the chr1/2/7/19 hulls out of the 4M pool until chr13 finished; that tail is the 2,356 s by which gw15 missed
bar A.

## Rule (before any number)
After `reads_in_region` returns, the region's `PermitGuard` shrinks to `min(held, bam_reads.len())` and returns
the difference to the pool (`PermitGuard::shrink_to`). The budget still bounds what is resident: the record count
actually fetched is the quantity the budget was estimating. Opt-in `RUSTLE_BUDGET_SHRINK=1`; unset = gw15 code
path. Paired with `RUSTLE_DISPATCH_DEFER=1`.

## Bars (row 1550's, unchanged)
A. wall <= 1.1 x (chr13 reads_in_region + detect_and_assign) of the same run.
B. six tables identical to gw_off; RSS < 250 GB.

## Prediction
chr13 releases 1.9M permits at t ~ 250 s; the xchrom families and hulls (sum ~ 2,300 s of work over 60+ free
threads) complete inside chr13's 4,000 s. Wall 4,300-4,600 s, A PASS. RSS rises above gw15's 215 GB because more
regions are resident alongside chr13; the risk is the 250 GB bar, stated now. Identity by construction (schedule
only), to be confirmed by the six comparisons.

## Pre-check
pre15: the 14-region panel with budget 50k, DEFER=1 + SHRINK=1 vs unset: identical tables; `[timing] budget
shrink` lines present.
