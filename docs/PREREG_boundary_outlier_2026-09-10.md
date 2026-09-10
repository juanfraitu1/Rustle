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

## Outcome (2026-09-10) — `bakeoff/human/ours_bnd10*`, `bakeoff/mcl1_bnd10`; default `--min-boundary-gap` raised 1000→5000 mid-measurement
| # | verdict |
|---|---|
| P0 | ✓ GTF byte-identical `888794bf`, assignments `8a057f68`, on both the gap=1000 and the revised gap=5000 default builds; gorilla MCL1 assignments `c87c7f71` unchanged |
| P1 | ✓ `DN_chr16_15368428_10` flagged `boundary_gap_right "18400"; low_confidence_reason "boundary"` — matches the hand-computed 18.4 kb gap and 2/41 = 4.9 % share exactly |
| P2 | ✓✓ **9/676 = 1.33 %** at gap 1000 (well under the 5 % ceiling); **7/676 = 1.03 %** at the revised gap 5000; gorilla MCL1 **2/403 = 0.5 %** |
| P3 | ⚠ **needed one round of tuning, honestly reported**: at the pre-registered gap floor (1000 bp), 2 of 9 flags were ordinary smooth-tail heterogeneity misread as isolated outliers — copies 9 and 25, gaps 1400/1150 bp, sitting in a locus with a CONTINUUM of ~30 alternative termini where the last two points happened to be > 1000 bp apart by chance (typical local spacing there is ~900 bp). **Raising `--min-boundary-gap` to 5000 removes both**; the remaining 7 (human) + 2 (gorilla) are all clean, isolated single-or-few-transcript outliers 5.5–139 kb from a well-supported majority cluster, spot-checked by hand (full boundary distribution per locus), 0 remaining concerns |
| | ⚠ **one genuinely uncertain case, reported rather than hidden**: copy 0 (the catalog's largest, 49.4 kb) has 15 independently-supported transcripts (2–9 reads each, summing far more than the outliers) terminating in a tight ~40 bp window at 11977724–11977765, and only 2 low-read transcripts (2, 4 reads) reaching close to the RefSeq-annotated end (12012756, off by ~15 bp). The tight 40 bp clustering across 15 independent transcripts is itself a strong signal — this may indicate the annotation's own boundary over-extends past what the RNA evidence supports for this LOC-named/predicted gene, not (or not only) a read-through artifact. Flagged either way, correctly conservative; not resolved further here |
| P4 | ✓ (narrowly, and reported as a narrower-purpose mechanism, not a substitute): excluding `low_confidence` transcripts at containment 0.5 (loose) reduces median \|relative deficit\| 5.7 % → 4.8 % and outliers-over-15 % 3 → 2, vs. tightening the containment floor to 0.8 (§6hq) which reaches 3.7 %/1 outlier on its own. Within the pre-registered 2× bound, but this flag targets the GTF EMITTER's own artifacts (independent of any copy catalog) and is complementary to, not a replacement for, a copy-aware containment floor in a downstream scorer |
| P5 | ✓ suite 871 / 0 / 11; assignments unchanged on both species |

### Reading
Unlike `--min-isoform-fraction` (§6hr, refuted — 75 % over-flagged), the boundary test works as intended:
**requiring BOTH a clear positional gap (≥ 5 kb) AND low relative depth catches exactly the readthrough class
(1–1.3 % of transcripts) without touching ordinary alternative-TSS/TES heterogeneity.** This is the practical
answer to the advisor's StringTie comparison: not a full per-junction flow-decomposition model, but a much
cheaper transcript-level test that reproduces its key behaviour — an isolated, low-support extension far past
a locus's well-supported boundary gets demoted, and everything else is untouched. Default: `min_isoform_fraction`
stays 0.0 (off, refuted for its purpose); `min_boundary_fraction` stays 0.0 (off) pending the user's decision,
with `min_boundary_gap` defaulting to the measured value of 5000.
