# Row 1793: which assembler rules are overfit to the primate samples? (pre-registered 2026-11-18)

Rows 1779 / 1780 failed two of three transfer bars: on mouse, fly and SIRV reads ours is the most sensitive tool but less precise
than StringTie and IsoQuant, and the chain score filter loses too many correct chains. User (2026-11-18): establish which rules
are overfit; the four datasets of rows 1779 / 1780 may be used for this diagnosis, and it may lead to proposed default changes
that the user rules on. Registered before any arm ran.

**Consequence, accepted by the user.** The four datasets were "never mined for a rule". After this row they are diagnostic data:
a default changed because of it must be checked on data not used here (for example LRGASP human WTC11, or a newly simulated
species) before it is quoted as transferring.

**Arms.** `tools/rustle_pipeline.sh assemble` at today's defaults, one change per arm (`r1793_rule_ablation/arms.sbatch`; a rule
without a switch is changed by `sed` on a snapshot of the driver, nothing in the repository changes):

| arm | change | rule judged |
|---|---|---|
| base | none | |
| floor1 | `RUSTLE_MIN_CHAIN_READS=1` | 2-read floor |
| frac0 | `--polish-isoform-fraction 0` | isoform fraction 0.02 |
| ism0 | `--polish-ism-ratio 0` | sub-chain collapse at ratio 0.7 |
| ret_off | `RUSTLE_POLISH_RETAINED_RATIO=1e9` | retained-intron filter at 10 |
| gate_off | `--assembly-junctions majority` | strict junction gate |
| ext0 | `RUSTLE_EXTEND_TRANSFRAGS=0` | transfrag extension |
| tss_off | `RUSTLE_POLISH_TSS=off` | read-proven start rescue |
| seed_primary | `RUSTLE_SEED_POOL=primary` | seed pool `good` |
| sub_off | `RUSTLE_FINAL_SUBCHAIN=off` | sub-chain drop of `final.gtf` |
| cs_off | `RUSTLE_CHAINSCORE=0` | chain score filter of `final.gtf` |
| cs_r1677 | `RUSTLE_CHAINSCORE_MODEL=r1677` | model r1769 against r1677 |
| on_jmaj | `RUSTLE_JUNCTION_MAJORITY=1` | opt-in, failed on primates (r1723) |
| on_psnap | `RUSTLE_POLISH_JUNCTION_SNAP=reads` | opt-in, in-sample only |
| on_rsnap | `RUSTLE_READ_JUNCTION_SNAP=10` | opt-in, failed on primates (r1618) |

The three `on_` arms are there because 3,075 of ours' 5,940 wrong chains on simulated mouse have a new intron within 10 bp of a
true one (row 1779), the error these rules address.

**Data.** New: fly simulated, mouse simulated, mouse ES (mouse contigs against GENCODE vM27; SIRV contigs against the SIRV
annotation). Tuning samples for comparison: OR6737 (RefSeq, `apetools/GGO.refseq.gtf`) and A119b (RefSeq on CHM13).

**Unit and table.** Exact intron chains (`curve.load`). Per arm and dataset, for `PREFIX.gtf` and `PREFIX.final.gtf`: chains,
correct, precision, and against base the correct chains gained and lost and the wrong chains gained and lost, counted separately
(`r1793_rule_ablation/score.py`).

**Reading, fixed now.** For an "off" arm the rule's effect is base against the arm (for the three `final.gtf` arms on
`PREFIX.final.gtf`, otherwise on `PREFIX.gtf`).
- A rule **pays** on a dataset when base has the higher precision and the rule costs at most 5% of the arm's correct chains
  (for floor1, frac0 and cs_off, whose purpose is a large cut: at most 10%).
- **OVERFIT** = pays on both tuning samples and fails to pay on at least two of the three large new datasets.
- **GENERAL** = pays on both tuning samples and on at least two of the three large new datasets. Anything else: MIXED, with the table.
- An `on_` arm is a **candidate** when, on at least two of the three large new datasets, it raises precision and loses at most
  2% of base's correct chains; its effect on the two tuning samples is shown next to it. A candidate is not a default.
SIRV is shown, not counted (60 chains).

**Chain score, feature by feature (no assembly run).** For model r1769, per feature the mean contribution to the score
(weights times design columns) over correct chains, on each dataset; a feature is named as a transfer problem when its mean
contribution on the correct chains of both simulated sets is lower than on both tuning samples by more than the largest
difference between the two tuning samples. Descriptive; it names where to look.

**Limits.** RefSeq is incomplete, so "wrong" on the tuning samples includes real chains, while the simulated truth is complete:
a rule can look worse on primates than it is. Rules interact; one-at-a-time arms do not show that. Two tuning samples stand for
twelve. No proposal from this row becomes a default without the user's ruling and a check on unused data.

## Result (2026-11-18): two rules OVERFIT, four GENERAL, one opt-in CANDIDATE

Jobs 4482349-4482353, all 75 arms exit 0; feature job 4482359. Tables: workspace `r1793_rule_ablation/`
(`score.{or,a1,dme,mmu,es,sirv}.txt`, `verdict.txt`, `feat.*.tsv`). Cell = precision the rule adds / share of the arm's correct
chains it costs (for `on_` arms: precision gained / share of base's correct chains lost).

| rule (arm) | OR6737 | A119b | fly sim | mouse sim | mouse ES | verdict |
|---|---|---|---|---|---|---|
| chain score filter (cs_off) | +0.015 / 1.0% | +0.064 / 4.2% | +0.129 / **15.3%** | +0.078 / **12.9%** | +0.084 / 7.8% | **OVERFIT** |
| ISM ratio 0.7 (ism0) | +0.128 / 1.4% | +0.060 / 3.1% | +0.043 / 1.4% | +0.054 / **9.8%** | +0.089 / **9.2%** | **OVERFIT** |
| isoform fraction 0.02 (frac0) | +0.052 / 3.5% | +0.028 / 7.1% | +0.039 / 2.2% | +0.095 / 3.0% | +0.110 / 4.4% | GENERAL |
| junction gate strict (gate_off) | +0.015 / 0.0% | +0.012 / 0.2% | +0.058 / 0.4% | +0.069 / 0.6% | +0.026 / 0.0% | GENERAL |
| retained-intron ratio (ret_off) | +0.015 / 0.4% | +0.016 / 1.2% | +0.000 / 0.3% | +0.004 / 2.4% | +0.030 / 1.6% | GENERAL |
| sub-chain drop (sub_off) | +0.017 / 0.8% | +0.011 / 1.1% | +0.002 / 0.2% | +0.002 / 2.1% | +0.003 / 0.2% | GENERAL |
| read floor 2 (floor1) | +0.125 / **10.3%** | +0.051 / 6.1% | +0.025 / 5.5% | +0.035 / 6.4% | +0.077 / 5.5% | MIXED |
| transfrag extension (ext0) | **-0.005** | **-0.007** | **-0.004** | **-0.025** | **-0.005** | MIXED: lowers precision on all five |
| seed pool (seed_primary) | **-0.003** | **-0.008** | **-0.084** | **-0.005** | **-0.010** | MIXED: lowers precision on all five |
| chain score r1769 over r1677 (cs_r1677) | **-0.018** | +0.041 / 2.6% | +0.034 / **10.5%** | +0.023 / 2.3% | +0.059 / -0.5% | MIXED |
| TSS rescue (tss_off) | no change | **-0.003** | no change | no change | no change | MIXED: acts on A119b only |
| polish junction snap (on_psnap) | +0.003 / 0.0% | +0.000 / 0.0% | +0.022 / 0.1% | +0.074 / 0.1% | +0.005 / 0.0% | **CANDIDATE** |
| read junction snap (on_rsnap) | +0.003 / 6.4% | -0.002 / 6.1% | +0.049 / 2.7% | +0.074 / 4.6% | +0.039 / 4.5% | not a candidate |
| junction majority (on_jmaj) | same counts as gate_off on all five | | | | | not a candidate |

Gained and lost chains separately, the cases that carry the reading (correct +/- ; wrong +/- against base):
- Transfrag extension, the rule against no extension, `PREFIX.gtf`: mouse sim correct +18 / -16, wrong +1,002 / -8; fly +2 / -0,
  +39 / -0; ES +30 / -14, +1,154 / -193; OR6737 +656 / -26, +3,132 / -971; A119b +1,666 / -295, +24,228 / -8,815. With the
  filter on top (`final.gtf`) extension costs correct chains on the simulated sets (mouse 22,892 with, 24,091 without; precision
  0.916 / 0.929) and gains them on the tuning samples (OR6737 25,944 / 24,305).
- Seed pool, against primary-only seeds: fly correct +5, wrong +781; mouse +44 / -3, +219 / -3; ES +55, +1,916; OR6737 +164,
  +1,050; A119b +286, +12,758.
- Polish junction snap on: mouse sim correct +2 / -33, wrong +11 / -2,740; fly +1 / -6, +2 / -223; ES +3 / -9, +60 / -894;
  OR6737 +10 / -15, +19 / -682; A119b +4 / -12, +48 / -508. Interaction: on `final.gtf` it loses 751 correct chains on mouse sim
  (3.3%) and 207 on fly (3.5%), though precision rises to 0.953 / 0.980.
- SIRV (shown, not counted; 60 true chains): base 47 correct / 0.398; ext0 49 / 0.505; on_psnap 48 / 0.471; ism0 51 / 0.227.

**Chain score features.** Share of correct chains scoring below the cut 0.055: OR6737 0.011, A119b 0.045, ES 0.082, mouse sim
0.140, fly 0.156. Features that meet the transfer-problem rule (mean contribution on correct chains; tuning OR6737 / A119b
against simulated fly / mouse): `r_cont5` 1.13 / 1.42 vs 0.76 / 0.70; `l_reads` 0.55 / 0.70 vs 0.25 / 0.28; `l_partial`
0.53 / 0.56 vs 0.13 / 0.30; `r_cont3` 0.19 / 0.33 vs -0.03 / -0.05; `l_cont3`, `l_length` (below 0.04 each). On real mouse ES
these features sit with the tuning samples (`r_cont5` 1.30, `l_reads` 0.71, `l_partial` 0.51, `r_cont3` 0.29). So the filter
rewards a chain for carrying truncated reads; real libraries of either clade have them, the simulated reads have fewer.

**Reading.** The chain score filter's loss is mostly a real-versus-simulated difference, not primate against mouse: on real mouse
it pays (7.8%), on simulated reads it does not. The ISM ratio is the one rule that fails on real non-primate data. Extension and
the seed pool are recall rules; by the precision reading they pay nowhere, the tuning samples included, so they are not overfit,
they are a trade. Nothing here changes a default. Proposals for the user's ruling, each needing a check on unused data:
(1) polish junction snap on; (2) ISM ratio re-set on non-primate data; (3) the truncated-read features of the chain score.
