# PREREG row 1546: unit-form read-star on the hotspot panel (2026-11-06 09:35, before the run)

## Question
`--read-star-unit` (the §6fc spliced-unit form, default before 2026-09-05) was measured 10x faster than the
genomic form and was replaced by user decision because the genomic form gave 0 wrong anchors on the paired 35
(docs/archive/2026-11/MEMORY_DIGEST.md). r1544 measured 469 vs 39 s on MCL1's reads. What is the panel-level
time and assignment trade on the nine hotspot families? This row quantifies a standing user decision; it
cannot adopt anything.

## Arms (hot9 panel, singlenode, hot10 recipe, bin10, wall limit 1 h)
- baseline: h9out/auto (genomic two-form, 1,654 s, 104 Assigned)
- unit: `--read-star-unit`, all else equal
## Bars
- Time: wall <= 660 s (0.4x).
- Assignments: Assigned per family vs baseline; count how many unit-form Assigned reads name a DIFFERENT copy
  than the baseline for the same read (disagreements), and how many baseline Assigned reads are lost. Report
  both; no pass threshold is set because the 09-05 criterion (anchors against truth) is not measurable here.
## Prediction
Wall 300-600 s; Assigned count within +-30% of 104; disagreements > 0 (the reason the genomic form was adopted).

## Result (2026-11-06 10:05): time bar FAIL (1,053 s); the forms assign almost disjoint read sets; closed

hot11 (job 4447631): `--read-star-unit` wall 1,053 s vs 1,654 baseline, RSS 36 GB. chr13 detect fell to
191 s (read-star sum 476 s vs ~3,500 s), but the region's BAM fetch took 779 s, so the wall bar (660 s) FAILS;
the bar ignored the fetch floor. Assigned rows: baseline 104, unit 42; 101 of the baseline's 104 (read,
family) keys are not assigned by the unit form and 39 of its 42 are new; 1 shared key disagrees on the copy.
Per family: MCL105 31 -> 3, MCL1166 35 -> 8, MCL1167 6 -> 0, MCL1165 18 -> 19, MCL1 8 -> 10. The unit form is
not a faster version of the same answer; it is a different answer, which is what the 2026-09-05 decision
already said. Not adopted; the decision stands.

Side finding (identity-preserving): the chr13 fetch (`reads_in_region`, 1.42M reads) varies 269-868 s across
runs of the same code (gw7 868, gw11 269, e2e 283, this panel 779), and whole-run fetch sums are 8,152-11,308 s.
The BAM (96 GB) is on the atuin NFS home. This explains Result 7's 724 s gap (gw7 vs e2e assign) better than
"run-to-run timing of the budget": gw7's chr13 fetch alone was 585 s slower. Next row: stage the BAM to
node-local storage before assign.
