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
