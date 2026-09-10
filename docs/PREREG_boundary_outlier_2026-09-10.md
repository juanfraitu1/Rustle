# PREREG — `--min-boundary-fraction`: flag a locus's most-distal boundary when it sits alone, far past the rest (2026-09-10, before the code change)

**Why this design, not the junction-crossing test discussed last night.** A true per-junction crossing-vs-
terminating test needs raw per-read CIGAR/position data at assembly time (`pass1_skeletons`) — a materially
bigger change. Re-examining the copy-4 case shows a cheaper, transcript-level test captures the SAME
phenomenon: the readthrough transcript's own exons up to 15386124 are ordinary (shared with siblings' EARLIER
stopping points), but no OTHER transcript in its locus reaches anywhere NEAR its final boundary (15401113) —
the next-most-distal sibling stops at 15382725, an **18.4 kb gap**. This is checkable directly from
`TranscriptRec` (n_reads, start, end, gene_tid), already grouped by `collapse_loci_groups`, no BAM re-analysis.

**Rule (behind `--min-boundary-fraction <f>` + `--min-boundary-gap <bp>`, defaults `0.0`/`1000` = off,
byte-identical GTF).** Per `gene_tid` group (≥ 2 members), independently on each side:
- RIGHT: bucket members by `end` (round to 50 bp); `reads_far` = Σ n_reads in the farthest bucket;
  `gap` = farthest bucket's position − second-farthest bucket's position.
- LEFT: same on `start`, nearest/farthest reversed.
A transcript in the farthest bucket is `boundary_low_confidence` iff `gap > min_boundary_gap` AND
`reads_far / group_total < min_boundary_fraction`. **Both conditions required** — depth alone
(`--min-isoform-fraction`, §6hr) over-flags ordinary heterogeneity; the gap requirement is what makes this a
different, narrower test. Combines with `--min-isoform-fraction`'s `low_confidence` by OR; both feed the
existing `--gtf-copy-set` exclusion (unmodified downstream). Emits `boundary_gap_left/right "N"` and folds into
`low_confidence` (adds `low_confidence_reason "boundary"|"depth"|"both"` when either fires).

| # | prediction | refuted by |
|---|---|---|
| P0 | flags at 0.0 / gap default: GTF byte-identical to `ours_final3.gtf` (`888794bf`); assignments untouched | any byte differs |
| P1 | at floor 0.10 / gap 1000: `DN_chr16_15368428_10` (copy 4) is flagged `boundary_low_confidence "true"` (gap ≈ 18.4 kb, reads_far/total = 2/41 ≈ 0.049) | not flagged |
| P2 | **≤ 5 %** of human MCL0 family transcripts newly flagged (this test is deliberately narrower than the depth-only one — row 808 refuted at 34–75 %; the gap requirement should keep it rare) | > 15 % |
| P3 | re-running `bench/isoform_copy_lift.py`'s phantom table on a GTF built with this flag active (excluding `boundary_low_confidence` transcripts from `--gtf-copy-set` evidence) shows the copy-4-class phantom/width inflation gone, **without** a comparable rise in flagged-but-legitimate isoforms (spot check the flagged set by hand, ≤ 2 false positives among the flagged) | ≥ 3 clearly-legitimate isoforms flagged |
| P4 | `bench/width_deficit.py`'s containment sweep, re-run with this flag's GTF as the `--gtf` input at containment floor **0.5** (the loose floor that originally let the readthrough through), now shows median \|relative deficit\| **within 2×** of the 0.8-floor result from §6hq (i.e. the flag substitutes for the containment floor) | > 5× worse |
| P5 | suite passes; O2 assignment outputs untouched on both species (this only touches the `--gtf` emitter) | any assignment-side number changes |
Human MCL0 primary; gorilla MCL1 secondary, report only.
