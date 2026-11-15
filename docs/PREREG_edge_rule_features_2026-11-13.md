# Row 1674 — an empirical multi-feature rule for family edges the default refuses (2026-11-13)

User goal: "find empirical rules to improve edge recall in the multi-copy gene families definition".
Families layer (O1). Anything found stays opt-in (`fam.copies.*` feeds assign; a default change is the user's call).
Not for `advisor_share`. No read is assigned to a copy.

**Why this and not another single lever.** Closed as single levers (register): aligner sensitivity (1662 0b), union of
records (1662 0c), protein edges (1589, 1662 0d/0f), family join (1662), shared reads (1512, 1664), neighbourhood
Jaccard (931/943), triangle/k-core, re-clustering (2026-11-07), SD-added edges (1581), wider truth (1663). The register
has no row that combines several de novo features in one rule ("decision tree" / "logistic": 0 hits).

**Recount on today's default (job 4470676, `r1674_edgerule/em.*.txt`), missed truth pairs by the stage that drops them:**

| sample | default P / R | no record | low SHARE | no exon base | other gates | passes (clustering) |
|---|---|---|---|---|---|---|
| A119b | 0.907 / 0.809 | 218 | 147 | 12 | 4 | 34 |
| HG002 | 0.727 / 0.444 | 77 | 86 | 30 | 3 | 34 |
| OR6737 (protein-hit truth) | 0.359 / 0.759 as scored | 116 | 107 | 5 | 1 | 15 |
| KB3781 (protein-hit truth) | 0.219 / 0.450 as scored | 22 | 90 | 9 | 1 | 10 |

**Ceiling.** The rule can only reach pairs with an alignment record that the gates refuse: A119b 163, HG002 119,
OR6737 113, KB3781 100 missed true pairs, against 2,473 / 1,734 / 1,081 / 1,358 false pairs at low SHARE alone and
~200,000 at "no exon base". The no-record class (218 / 77 / 116 / 22) is out of reach.

**Population (frozen).** Gene pairs as `edgemiss.py` counts them (expressed coding genes of the truth universe,
different symbols) with >= 1 PAF record between their loci, not in one default family, not at stage `edge`.
Gorilla: protein-hit truth (`DENOVO_GGO_ENS_MAP`), pairs whose genes map to the same human gene set aside (undecidable,
r1655 amendment 1).

**Features (frozen, `bench/denovo_eval/edgefeat.py`; loci + PAF + optional SEDEF BED only, no annotation, no reads):**
share, share_union, identity, aln_bp, n_records, exon_bp, exon_frac (exon-to-exon bases / shorter exon length),
span_cov_min / max (aligned share of each locus span), nonexon_frac (aligned bases outside exons), exons_min / max,
exon_len_ratio, span_ratio, same_chrom, log_dist, sd_held (one SEDEF pair has a mate over each locus), rank_max (worse of
the two partner ranks by pooled matches), degree_min / max (number of partner loci). Shared molecules are left out
(row 1664; tissue-dependent).

**Model class (frozen, `bench/denovo_eval/edgerules.py`).** A conjunction of <= 3 threshold conditions; thresholds =
deciles of the training true pairs; selected = most training true pairs with training precision >= 0.80 and >= 20 true
(ties: fewer conditions, then precision). One rule per training sample, no second pick.

**Step 0 (edge level, human, spent samples used as each other's test).** Fit on A119b, score on HG002; fit on HG002,
score on A119b. **GO** if a selected rule gives, on the sample it was not fitted on, >= 20 new true pairs at precision
>= 0.75. If both go, the one with more held-out true pairs is THE rule. Otherwise **STOP**: report and ask the user; no
refit, no looser bar.

**Step 1 (gorilla transfer, edge level, frozen rule, no refit).** OR6737 + KB3781 pooled: new true pairs >= 20 and
precision >= 0.60. Void if the pooled reachable true pairs in the population are < 50.

**Step 2 (family level, base track where an inflation curve exists).** The rule's locus pairs fed as `--extra-edges`
to `mcl_families` with the shipped flags, on A119b and HG002 (`port1593` loci). Bars, each sample: recall >= baseline
+0.01; precision >= the base inflation curve (`r1662_edgemiss/{A,H}.I*`) at the arm's recall; true pairs lost <= 10% of
true pairs gained. Then the same route on OR6737 and KB3781: pooled true gained >= 10, pairs added precise >= 0.60
(same-human-gene pairs set aside), lost <= 10% of gained.

**Found** = step 0 GO and steps 1 and 2 pass. Then: opt-in flag proposal to the user, nothing wired by default.

## Result: STOP at step 0 (2026-11-13, job 4470684, `r1674_edgerule/feat.{A,H,OR,KB}.tsv`)

Populations (record, unlinked, gate-refused): A119b 202,653 pairs / 159 true; HG002 167,424 / 118; OR6737 8,152 / 111;
KB3781 226,839 / 99.

- Fit on A119b: no conjunction reaches the training floor (precision >= 0.80 with >= 20 true). No rule to test.
- Fit on HG002: selected `span_cov_min>=0.2727 & log_dist>=5.7322 & log_dist<=7.4034` (train 25 true / 5 false);
  on A119b 15 true / 49 false, precision 0.234. Bar: >= 20 true at >= 0.75. FAIL.

STOP as registered. Steps 1 and 2 not run; nothing wired.

**Descriptive, after the verdict (not a rule; each sample scored on itself, the most favourable case).** Best
precision any conjunction of <= 3 thresholds reaches:

| sample | with >= 10 true | with >= 20 true | with >= 40 true |
|---|---|---|---|
| A119b | 12 / 0 | 20 / 14 (0.59) | 41 / 87 (0.32) |
| HG002 | 12 / 0 | 20 / 12 (0.63) | 47 / 82 (0.36) |
| OR6737 | 11 / 0 | 20 / 9 (0.69) | 55 / 118 (0.32) |
| KB3781 | 11 / 1 | 20 / 13 (0.61) | 41 / 83 (0.33) |

Even in-sample no sample reaches 20 true pairs at 0.75, and the conditions chosen differ in every sample (distance on
HG002, exon count on A119b, record count on gorilla). The features do shift between true and false pairs (true pairs:
same chromosome, SD-held, partner rank ~10 vs ~280, exon bases shared 0.2-0.4 of the shorter member), but a few
hundred false pairs share each of those properties, against ~100-160 true pairs.
Half of the true pairs in this population have SHARE 0: they share exon bases (median 0.2-0.4 of the shorter member) but
no exon has a 50% reciprocal partner.
