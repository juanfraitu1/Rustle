# PREREG — indel PSV columns in read-star (`--indel-psv`, 2026-09-09, before the code change and before any run)

**Why (user, step 2 hold):** "PSVs can also be indels." Today they are not: `read_star_columns` makes a column
only where ≥ 2 candidates carry an aligned BASE that differs (`I` in the read-vs-copy CIGAR leaves the read
position `None`, `D` is invisible). Between reference copies, indel events are 0.13 (human MCL0) / 0.05
(gorilla MCL1) per substitution column over all pairs, **0.2–0.5 among the most similar pairs** (genomic unit
spans, asm20 `--eqx`, non-self pairs ≥ 5 kb aligned); no pair in either family is separated by indels ONLY
(0 pairs with X = 0). The best_by_psv PREREG (f2cb3c19) had already found 3/15 mis-tied molecules that "differ
only by an indel — a policy case". ⚠ Distinct from `--origin-drop-indels`, which touches only the origin
certificate (read-vs-reference indels the individual carries, absent from every copy, non-discriminating).

## Rule (behind `--indel-psv`, default OFF; OFF = byte-identical)
1. In `read_star_stream`, per (molecule, candidate) hit, record every `I` (read bases the copy lacks) and `D`
   (copy bases the read lacks) event of length **≥ 3 bp** as `(read_pos, kind, len)` in forward read
   coordinates (`I`: the run's first read base; `D`: the boundary before the next read base). `N` never counts.
2. Events of all candidates are clustered along the read: sorted by position, a new cluster starts when an
   event begins more than **20 bp** after the current cluster's end. One column per cluster, at its first
   position, appended after the substitution columns.
3. Alleles: candidate has ≥ 1 event in the cluster → `'1'`; covers the column position (aligned base at `pos`
   or `pos − 1`) with no event → `'0'`; otherwise `None`. The read's observation is always `'0'` (the read is
   its own coordinate system: a copy either needs a gap against it here or it does not). Emitted only if both
   `'1'` and `'0'` occur among the candidates (decisive by construction).
4. Two-form star (§6fp): a column is emitted only among candidates taken in the SAME form (genomic or unit);
   a molecule whose candidates mix forms gets indel columns within each form-set separately. (A retained
   intron is an `I` against the spliced unit and nothing against the locus — a fake column otherwise.)
5. Error model: the same per-column rate as substitution columns (`error_rate` 0.003, mismatch = e/3 in the
   LLR), i.e. one decisive indel column is worth one substitution column. No new weights.
6. Everything downstream is untouched: `BubbleGraph`, `copy_pair_significance`, `best_by_psv`, the origin
   certificate, the AS-tied gate, `--dump-star` (indel columns appear as `'0'`/`'1'` alleles).

## Predictions (contested set = `origin_rejected == 0 ∧ n_candidates ≥ 2`, reproduces the stderr line:
human `ours_odi` 1,143 = 230 assigned / 531 tied / 382 ambiguous; gorilla MCL1 33 = 4/28/1; MCL7 11 = 0/10/1.
⚠ 134 human / 5 MCL1 / 5 MCL7 of the `tied` are `tie_outside_catalog` and are forced Tied by rule — the
convertible tied pool is **397 / 23 / 5**.)
| # | prediction | refuted by |
|---|---|---|
| P0 | flag OFF byte-identical: human `4048494b`, gorilla MCL1 `d6605062`, MCL7 `17b92323`; full escapes `91081887` / `ff0b8f16` | any md5 differs |
| P1 | ≤ 30 % of the 1,143 human contested molecules gain ≥ 1 indel column (indels are rarer than substitutions and must fall inside the read footprint) | > 50 % |
| P2 | the 230 human assignments are stable: ≥ 99 % keep status AND copy | > 2 % change copy |
| P3 | tied → assigned conversions among the 397 convertible: **between 8 (2 %) and 80 (20 %)** | 0, or > 120 |
| P4 | excision control on the copy receiving the MOST new (tied→assigned) reads: catalog without it, re-run with the flag — ≥ 95 % of those reads abstain | < 90 % |
| P5 | gorilla MCL1: 1–8 of the 23 convertible tied become assigned, the 4 base assignments unchanged; MCL7: ≤ 2 | MCL1 0 or > 12; any base assignment changes copy |
| P6 | indel-column lengths: ≥ 70 % of decisive indel columns are ≥ 10 bp (SDs differ by Alu-scale and larger events; a spike at 3–5 bp = homopolymer/alignment noise) | < 50 % |
Runs: human with `--origin-drop-indels --indel-psv` (arm A, the comparison to `ours_odi`) and without
`--origin-drop-indels` (arm B vs `ours_psv`, `1f8157bd`) — the two flags are independent and the default of
the first is still held. Human and gorilla never pooled.

## Outcome (2026-09-09) — `bakeoff/human/ours_indelA` (+`_dump`), `ours_indelB`, `bakeoff/mcl1_indelA`, `mcl7_indelA`; scorer `bench/o2_indel_psv_score.py`
| # | verdict |
|---|---|
| P0 | ✓ flag OFF byte-identical on the 22 shared columns (`4048494b` / `d6605062` / `17b92323`; the PREREG md5s predate step 1's `aligner_disagreement` column, the 23rd); full escapes `91081887` / `ff0b8f16` ✓ |
| P1 | ✓ 229/1,143 = 20.0 % of the contested molecules gain ≥ 1 indel column (median 1, max 4); family-wide 1,902 molecules / 7,217 columns |
| P2 | ⛔ **220/230 = 95.7 % keep status and copy** (< 99); 7 `assigned`→`tied`, 1 →`ambiguous`, 2 change copy (0.9 %, not refuted) |
| P3 | ⛔ **0 of 397 convertible `tied` become `assigned`**; 9 become `ambiguous`; 1 `ambiguous`→`assigned` |
| P4 | — moot (no conversions) |
| P5 | ⛔ gorilla MCL1 0 of 23 (predicted 1–8), the 4 base assignments unchanged; MCL7 0 (✓) |
| P6 | ⛔ **3,292/7,217 = 45.6 %** of indel columns have an event ≥ 10 bp (predicted ≥ 70 %) |
Arm B (no `--origin-drop-indels`, vs `ours_psv`): 14 → 15 assigned, 0 conversions, 7 non-contested rows change
status (the bk moved, so the certificate's target moved).

### Why (from `--dump-star`, `scratchpad/indel_dump_analysis.py`)
* **The tied pool has no indel to use.** 462 of the tied molecules that gained indel columns have every column
  CONSISTENT with their substitution-best copy — and identical for its K = 0 twin. Twins (e.g. 6/7/8) do not
  differ by an indel inside a 1.5–3 kb read footprint: the reference-copy measurement (0.2–0.5 indel events per
  substitution column among the most similar pairs, X median 5 over ~20 kb) predicts ≈ 0.2 events per read.
  P3 = 0 was the expected value; the prediction was wrong, not the channel.
* **Terminal events are artifacts.** 319 columns lie within 20 bp of a read end (306 at the START); **33.9 %**
  of them contradict the substitution-best copy vs **9.9 %** of the 5,904 internal ones. Five of the seven
  lost copy-2 assignments were tipped by ONE such column (read pos 8–16: copies 0–4, 9–12, 19–22 "need a
  gap", copies 5–8 do not) — the aligner's end-gap placement, not a PSV.
* **The bk score is not pairwise (the §6ha weakness, again).** For those seven, the pairwise LLR still favours
  copy 2 over copy 6 by 34.5 (5 columns net, `margin` = −34.5 with bk = 6), but `psv_score` counts every
  column a candidate carries: copies 6/7/8 carry the read's 57 bp segment that reference copy 2 lacks, so
  they collect matches at substitution columns inside it where copy 2 is `None`; the scores were already
  near-equal and one artifact column flipped bk to 6, whose twin 7 then forces `Tied` (K = 0).
* The 9.9 % internal disagreement is the individual's polymorphic indels (a read from copy 2 carrying the
  6/7/8 allele of a 57 bp indel, §6hc) plus alignment representation (a divergent tract as an indel against
  one copy, as substitutions against another).

### Reading
As pre-registered, indel PSV columns add nothing to O2 (0 conversions in three families) and cost 7 human
assignments through two diagnosable defects: terminal-event artifacts (fixable inside the flag: drop events
within ~30 bp of the read end or of the hit's own query bounds) and the non-pairwise `bk` score (a shared-
machinery change: choose `bk` by pairwise duels). Neither is applied here. **Flag stays OFF.**
