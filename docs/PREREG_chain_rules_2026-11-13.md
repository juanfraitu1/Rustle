# PREREG r1677: empirical precision / recall rules for assembly from the tool lessons (2026-11-13)

User goal (2026-11-13): "find some empirical rules to precision and recall in assembly based on what we have learned from
stringtie, flair, isoquant, isocon". Written after the one-chromosome smoke test of the feature script (row counts only),
before any feature was tabulated against the label.

This is r1619's question again. r1619 screened features computable from our GTF alone: +0.025 over the read-count curve on
gorilla, -0.006 held out on A119b (FAIL). That screen is not rerun. New here: features from the reads, taken from what each
tool does, and the recall half (1-read chains) scored on the same scale.

## Data
Multi-exon exact intron chains vs RefSeq. Candidates = chains of the 1-read-floor run (`RUSTLE_MIN_CHAIN_READS=1`) plus any
default chain it lacks; `in_default` marks the default set. Seen samples: gorilla OR6737 (`r1628/orrun/or.gtf`,
`r1616/floor1/or.gtf`), human A119b (`r1617/{default,floor1}/a.gtf`). Frozen test: gorilla KB3781 (never used for a chain
rule). Features: `bench/famcmp/chain_rules_feat.py` (one pass over primary spliced records), workspace `r1677_chainrules`.

| tool lesson | feature |
|---|---|
| StringTie: minor-isoform fraction, retained intron, 5' fragments (nascent, "lacks an end intron" r1636) | `locus_frac`, `retention`, `cont5` / `cont3` = reads whose chain strictly contains the chain and extends it 5' / 3' |
| StringTie / r1618: partial reads as support | `partial` = reads whose chain is a strict contiguous run of the chain |
| FLAIR: every read chain a candidate, 3-read cut, subset filter | `full` (full-chain reads), `cont5` + `cont3` |
| IsoQuant: intron support floor, junction correction, polyA | `min_isup` / `min_isup_rel` over all reads, `wobble`, `pas`, `a_rich`, `end3_conc` |
| IsoCon: read sequence consensus | `de` (mean read divergence), `mq0_share` (tied reads) |
| bambu: chains per gene, relative support | `locus_chains`, `locus_frac` |

## Part 1: rules table (descriptive; the deliverable whatever Part 2 says)
One row per feature class (fixed thresholds below, set now): share correct yes vs no within read strata 2 / 3-5 / 6-20 / 21+
on default chains, both samples. A row is called a **stable rule** when the yes/no ratio of share correct is on the same side
of 1 in every stratum of both samples with >= 50 chains on each side, and the stratum-pooled ratio is <= 0.67 or >= 1.5 in
both. Classes: `cont5 >= full`; `cont3 >= full`; `partial = 0`; `partial >= 5 x full`; `min_isup_rel < 0.05`; `wobble`;
`retention`; `locus_frac < 0.05`; `locus_chains >= 10`; `end3_conc < 0.5`; `pas = 0`; `a_rich >= 12`; `mq0_share >= 0.5`;
`de` above the sample median; `n_introns = 1`. Closed levers (floors, snap, subset merge, nascent, end mode, IsoCon
consensus, tool agreement) enter as rows with their recorded figures, no new arm.

## Part 2: one score (step 0 with a STOP)
Logistic regression on the features, each log-transformed where a count and standardised **within the sample** (no absolute
depth). Two fits: default chains only (precision), and the superset (swap).
- **Step 0, cross-sample:** fit on OR6737 and score A119b, and the reverse. Measure A = precision of the score curve minus
  the read-count curve at 0.90 of the default correct chains. **GO when A >= +0.02 in both directions** (r1619's held-out
  bar; r1619 had -0.006). Also reported: 5-fold by chromosome within each sample, and the r1619-feature subset.
- **Swap (recall half), cross-sample:** rank the superset by the score fitted on the other sample's superset; keep as many
  chains as the default has (a sample-internal cut). Measure B = correct among them / default correct - 1.
  **GO when B >= +0.02 in both directions.** New and lost chains are reported separately.
- Neither GO: STOP, ship Part 1, nothing ported.
- **Frozen test (only on a GO):** fit on OR6737 + A119b pooled, coefficients frozen, KB3781 scored once (needs one 1-read
  assemble for the swap). PASS = the same bar (A >= +0.02, or B >= +0.02) for the part that had its GO. A PASS proposes an
  opt-in score column / filter; defaults stay the user's call.

## Limits stated now
RefSeq is incomplete and uneven between the samples (gorilla XM models; A119b's default 2-read chains are 0.066 correct), so
anything tied to expression looks good; hence strata and curve-against-curve. Both fit samples are testis. Multi-exon only.

## Step 0 outcome (2026-11-13, job 4470697, `r1677_chainrules/{rules,simple}.txt`)
- **A: GO.** Score fitted on the other sample vs the read-count curve at 0.90 of the default correct chains: A119b 0.257 vs
  0.217 (+0.039), OR6737 0.447 vs 0.399 (+0.048). r1619's features alone, same protocol: +0.006 / +0.017. 5-fold by
  chromosome within a sample: +0.064 / +0.071.
- **B: no GO, STOP.** Equal-count swap on the 1-read superset: A119b +1.8% correct (new 38,378 chains at 0.040, lost at
  0.022), OR6737 +0.55% (0.097 vs 0.086). Not tested further; the 1-read class stays closed.
- Part 1 stable classes: `cont5 >= full` (pooled ratio 0.13 / 0.10), `cont3 >= full` (0.44 / 0.20), `retention` (0.56 / 0.48).

## Amendment 1 (2026-11-13, written before any KB3781 feature table existed)
Seen-sample look (post hoc, `simple.py`): one class carries most of the gain without a fitted number. Dropping default chains
with `cont5 >= full` (at least as many reads run through the whole chain and continue 5' of it as end with it) keeps 0.990 /
0.973 of the correct chains and lifts precision 0.352 -> 0.382 (OR6737) and 0.176 -> 0.217 (A119b), +0.026 / +0.033 over the
read-count curve at the same correct count.
**Frozen test on KB3781** (`r1622/kb/kb.gtf`, default chains, scored once):
- T1, score: pooled OR6737 + A119b fit, coefficients frozen. PASS: A >= +0.02 at 0.90 of the default correct chains.
- T2, plain rule: drop `cont5 >= full`. PASS: keeps >= 0.95 of the default correct chains AND precision >= the read-count
  curve at its correct count + 0.02.
- Descriptive, no bar: T2 on `PREFIX.final.gtf` of the three samples (does the r1498 sub-chain drop already remove the class?).
A T2 PASS proposes the rule as an opt-in filter / tag on the emitted GTF; T1 PASS alone proposes a score column.

## Outcome on KB3781 (2026-11-13, job 4470748, `r1677_chainrules/{rules.kb,simple.kb,table.kb}.txt`): FAIL as registered
- **T1 FAIL.** Frozen pooled score vs the read-count curve at 0.90 of the correct chains: 0.424 vs 0.418 (+0.006; bar +0.02).
  +0.024 at 0.75, +0.048 at 0.50: the score only helps at the strict end, as F did in r1619.
- **T2 FAIL on the gain, PASS on the keep.** Drop `cont5 >= full`: 3,959 of 72,102 chains dropped, 218 correct among them
  (0.055); keeps 0.992 of the correct chains; precision 0.358 -> 0.376, read-count curve there 0.363 (+0.013; bar +0.02).
- **On `final.gtf` (descriptive):** the r1498 sub-chain drop already removes about half the class on the gorillas and a third on A119b. Left over: KB3781 1,904
  chains at 0.102 (+0.004 over the curve), OR6737 2,980 at 0.080 (+0.009), A119b 30,436 at 0.031 (+0.019; precision 0.191 ->
  0.219 keeping 0.976 of the correct chains).
- Nothing ported by this verdict, defaults unchanged (the user's later opt-in port is the last section). The two fit samples are testis, KB3781 is fibroblast with a class half the size
  (5.5% of chains vs 8.8% and 21%); the class is as wrong there as elsewhere, there is just less of it to remove.

## The rules (Part 1; share correct yes / no pooled over read strata as a ratio; OR6737 / A119b / KB3781)
Hold in every read stratum of all three samples:
1. **5' truncation.** A chain with at least as many reads running through it and continuing 5' as reads ending with it is
   wrong at any depth: ratio 0.13 / 0.10 / 0.18; at 21+ reads 0.069 vs 0.761, 0.033 vs 0.598, 0.077 vs 0.703.
   (StringTie's missing-end-intron and nascent classes, FLAIR's subset filter.)
2. **3' truncation** (`cont3 >= full`): ratio 0.44 / 0.20 / 0.42. Weaker: real alternative 3' ends are in the class.
3. **Retained intron** (an exon covers an intron with >= 5x the chain's reads): 0.56 / 0.48 / 0.47.
   Rules 2 and 3 hold in direction but buy nothing over a read cut as filters (+0.004 / +0.014 and 0.000 / -0.002 on the
   seen samples); only rule 1 does.
4. **Read count** orders correctness in every tool, but the level is the sample's: 2-read chains 0.17 / 0.066 / 0.15, so no
   absolute floor transfers; r1616-r1618 floors and the 1-read swap here all admit chains at or below the weakest kept class.
5. **Depth without agreement** (r1633 tables, two samples): chains with 21+ full-chain reads that neither StringTie nor
   FLAIR outputs are 0.115 / 0.051 correct against 0.795 / 0.665 when both output them. Whether rules 1-3 account for that class was not measured (the feature tables carry no tool flags).
Sample-dependent, not rules: weakest intron below 5% of the locus top and locus fraction below 5% (0.58 / 0.84 / 0.58 and
0.51 / 0.76 / 0.53: gorillas only), junction wobble (0.70 / 0.99 / 0.77), no polyA signal and A-rich 3' end (A119b only,
0.56 and 0.49; gorillas 0.9-1.0), tied reads (A119b only). No signal anywhere: read divergence (0.99 / 0.98 / 1.02, the
IsoCon lesson), 3' end concentration above 2 reads, partial-read support.
Recall: no rule. 69% of the reference chains other tools find and ours misses have one full-chain read (r1614), and one
score over the 1-read superset trades in chains at 0.040 / 0.097 for chains at 0.022 / 0.086.

## Post hoc check: chains no primary read spells (2026-11-13, `r1677_chainrules/full0.txt`; after the KB3781 outcome)
`cont5 >= full` is true when both are 0, so default chains with no full-chain primary record sat in the truncation classes
and in `pas = 0` / `end3_conc < 0.5` by construction. They are few and are their own weak class: 720 / 9,501 / 454 chains
at 0.147 / 0.019 / 0.088 correct (OR6737 / A119b / KB3781). On chains with `full >= 1` only, rule 1 is stronger (ratio 0.08 /
0.10 / 0.16; every stratum of every sample 0.02-0.07 correct), rule 2 is 0.40 / 0.24 / 0.43, and the polyA rows do not move
(no signal 0.96 / 0.61 / 0.97; A-rich 0.88 / 0.47 / 0.90: still A119b only). The `end3_conc` signal at 2 reads was this
artefact. Rule 1 on `full >= 1` as a drop rule: +0.025 / +0.024 / +0.013 over the read-count curve, keeping 0.994 / 0.977 /
0.993 of the correct chains; the registered T2 verdict does not change.

## By tool: what each one's filter is worth on our chains (recorded figures; no new arm)
| tool | its rule | on our chains |
|---|---|---|
| StringTie | minor-isoform fraction | 2% cut costs 947 correct chains for +0.05 precision on OR6737 (r1614); locus fraction < 5% is a class on the gorillas only (here) |
| StringTie | retained-intron filter | class is 0.47-0.56 as often correct, no gain over a read cut as a filter (here) |
| StringTie | 5' fragments: nascent tag, transfrag missing an end intron | nascent tag 211 / 2,030 chains, 0.142 / 0.073 correct, FAIL held out (r1635); the read-level 5' class here is 10-20x larger and 0.02-0.07 correct |
| StringTie | partial reads extend and support chains | extension +630 / +1,371 correct at new-chain precision 0.173 / 0.064, FAIL by 0.002 (r1618); partial support carries nothing as a feature (here) |
| StringTie | flow | port reaches Jaccard 0.81 gorilla / 0.59 A119b and loses on both axes against the default (r1658) |
| FLAIR | every read chain a candidate + realignment + 3-read cut | 7.9k new correct for 7.4k lost at precision 0.062 vs 0.176 on A119b (r1640) |
| FLAIR | subset reads merged into the longer chain | 15,464 new chains at 0.049, FAIL (r1637); the same containment used as a drop rule is rule 1 |
| IsoQuant | intron support floor, junction snap | floor costs a third of our correct chains; snap 10 bp loses 1,598 correct (r1618); wobble is a class on the gorillas only |
| IsoQuant | polyA-aware ends | its 3' ends beat ours (0.66 / 0.58 vs 0.53 / 0.50 within 10 bp); no read-end statistic transfers (r1660); polyA signal marks correct chains on A119b only |
| IsoCon | read sequence consensus | worse in both classes (r1478), 9.6% recovery with 58% artifacts (r1515); read divergence carries nothing (ratio 0.98-1.02) |
| bambu | trained score on read-class features, strict cut | ahead of ours only at its strict cuts on A119b (r1641); our score likewise helps only at the strict end (+0.05 at half the correct chains, +0.006 at 0.90 held out) |
| any | absolute read floor | 1-read chains 0.04-0.10 correct; no score swaps them in at a gain (B: +1.8% / +0.55%) |

## Port (2026-11-13, user: "lets make the 5' truncation drop opt in")
Driver stage `trunc5` (`tools/rustle_pipeline.sh trunc5 --bam B --fasta G --out PREFIX`, `bench/famcmp/trunc5.py`) writes
`PREFIX.trunc5.gtf` = `PREFIX.gtf` minus multi-exon chains with at least one full-chain primary read and at least as many
primary reads continuing 5' of the chain. Differs from the registered T2 rule in one point: chains no primary read spells
are kept (the post hoc `full >= 1` form). Named stage only, never part of `all`; `PREFIX.gtf` untouched.

## Row 1692 (2026-11-13, registered before the run; user: start "HG002 trunc5"): the shipped `trunc5` stage on a fourth sample
Sample: HG002 Kinnex on CHM13 (`HG002_kinnex.t2t.bam`), default `assemble` through the driver, then the driver's `trunc5`
stage; truth = CHM13 RefSeq exact intron chains, as for A119b. HG002 was not used to find or tune the rule (found on
OR6737 + A119b, frozen test KB3781).
Bar, the KB3781 bar verbatim: precision of `PREFIX.trunc5.gtf` at least +0.02 above the read-count curve of `PREFIX.gtf` at
the same number of correct chains. Reported with it: chains dropped and the share of them correct, correct chains kept.
Parity: the stage's dropped count must equal `cont5 >= full and full >= 1` on the feature table, or the row is INVALID.
Neither outcome changes the stage: it stays an opt-in named stage, never part of `all`.
Outcome: **FAIL as registered.** HG002 Kinnex on CHM13 (job 4470864; default assemble 1,140 s, 153,033 multi-exon chains, 35,745 correct, precision 0.234). The stage drops 16,889 chains, 546 of them correct (0.032); keeps 0.985 of the correct chains; precision 0.259 against 0.244 for the read-count curve at the same correct chains: +0.015, bar +0.02. Parity holds (16,889 = the feature-table rule). Four samples now: +0.025 (OR6737), +0.024 (A119b), +0.013 (KB3781), +0.015 (HG002) on the `full >= 1` form: the class is real everywhere (3-6% of dropped chains correct) but the gain over simply raising the read floor is under the bar on both samples the rule was not found on. Stage stays opt-in.
