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
