# PREREG row 1545: per-family read-star pilot with early exit (2026-11-06 08:10, before any build or run)

## Question
The eight chr13 hotspot families (plus MCL0; the `hot9` panel) align ~30k molecules each for 0-35 unique
assignments (baseline `~/asg_prof/h9out/auto.*`: 104 Assigned of 212,691 rows; per family 2, 0, 31, 2, 18, 35,
6, 2, 8). Can a per-family pilot skip the remainder where the pilot assigns nothing, at a bounded loss?

## Rule (opt-in, `RUSTLE_STAR_PILOT=N`, unset = current behaviour, byte-identical)
When a family has more than N representative molecules: run read-star and the per-read tail on a deterministic
stride sample of N molecules (every k-th rep, k = ceil(reps/N)). If that pilot yields zero `Assigned` reads,
the remaining molecules get empty read-star observations and go through the same per-read tail (they land on
the existing no-candidate path, n_candidates 0, status ambiguous). Otherwise the remainder is aligned as now.
Output order is reps order either way. N = 6,000, fixed now, not tuned after the panel.

## Bars (hot9 panel, singlenode, `hot9.sbatch` recipe, wall limit 2 h)
- Identity: the unset arm reproduces `h9out/auto.*` on all five tables (`cmp`). Required before reading the pilot arm.
- Loss: Assigned reads lost vs baseline <= 10 of 104 (10%), counted per family; no family with >= 10 baseline
  assignments may skip.
- Time: pilot-arm wall <= 1,000 s vs 1,664 s baseline (0.6x).
## Predictions
Families with 0-2 assignments skip (MCL104, MCL0, MCL106, MCL11; expected ~0.67 skip probability each at 2);
MCL1 (8) and MCL1167 (6) skip with probability ~0.2-0.3; MCL105/1165/1166 run in full. Expected loss 5-9
Assigned reads, wall ~900-1,100 s. If loss > 10: FAIL as registered, the gate needs a second signal (not a
retune of N). Genome-wide arm only if the panel passes; adoption as a default is the user's call (it changes
outputs: quant and famcn counts move with the skipped reads).

## Result (2026-11-06 09:20): loss bar PASS, time bar FAIL as registered; lever closed

hot10 (job 4447570, singlenode). Unset arm: 1,654 s, five tables identical to h9out/auto (identity precondition
holds; the code path with the flag unset is unchanged). Pilot arm (N = 6,000): 1,685 s, RSS 76 GB, Assigned
97 of 104 (lost 7: MCL11 2, MCL1167 6; MCL1166 +1), no family with >= 10 assignments skipped. Loss <= 10 PASS.
Time FAIL: only the three cheapest families skipped (4,505 of 13,514; 5,330 of 15,989; 5,603 of 22,409
molecules; 14-63 s each in full). The six expensive families (20k-30k molecules, 104-1,376 s each) all
assigned something in the pilot and ran in full, so the wall did not move. Prediction (MCL0/MCL106/MCL11 skip,
wall ~1,000 s) was wrong: with 2 true assignments in 28k reads the stride pilot found them in 2 of 3 cases.
Conclusion: the hotspot cost sits in families that do yield assignments; a zero-yield early exit cannot
remove it. The lever stays opt-in and unadopted. A yield-rate gate (skip when the pilot's rate predicts < k
assignments in the remainder) would trade those assignments for time and is a different rule, not a retune.
Timing label fix: the pilot arm printed "{reps} of {aligned}"; order corrected after the run (cosmetic).
