# Pre-registration: a per-read power analysis of O2's abstention (KEY=o2_power)

**Written 2026-09-29 BEFORE any O2 status was tabulated by `n_decisive`, by copy count, or against simulated truth.**
Analysis only: no change to O2's rule, its defaults, `src/` or `bench/`; nothing is committed. Code and outputs live in
`/mnt/linuxdisk/tmp/rustle_figures_dev/o2_power/`. Origin of the question: `docs/VAULT_METHODS_AUDIT_2026-09-28.md`
item 3 ("Coverage-power calculation (longcallR: power above 80% only over about 50× at SOR = 2): could justify the O2
abstain rate. A framing point, not a lever."). Human and gorilla numbers are never pooled.

## 0. What the longcallR number is (quoted, so it is not stretched)

Huang & Li 2026 (Nat Methods, doi 10.1038/s41592-026-03045-6), Methods, allele-specific junctions (PDF p. 8):
"R = |Ra1|×|Rp2| / (|Rp1|×|Ra2|), SOR = ln(R + 1/R). ... To simplify the power analysis, we assume that the number of
reads from both haplotypes is identical. ... Significant ASJs are identified using Fisher's exact test (P value <0.01).
Our results show that for an SOR of 2.0, the statistical power exceeds 80% when the total read coverage is greater than
50×." SOR = 2 means R + 1/R = e² ⇒ **odds ratio R ≈ 7.25**. It is a LOCUS-level test (a 2×2 table of haplotype × junction
presence over already-phased reads). It says nothing about whether one read can be phased. O2 decides one READ at a
time, so the analogue of "coverage" for O2's abstention is not locus depth.

## 1. The rule being modelled (read from the code, `src/rustle/vg_family/copy_assign_pipeline.rs` 2565-2930, defaults)

Per molecule that passed the AS-tied gate, with k family candidates (partners removed) aligned by read-star:
columns = read positions where ≥ 2 candidates carry an aligned base and they do not all agree; observation = the read's
base. `bk` = the maximin of the pairwise LLRs (λ = ln(3(1−e)/e) per column, +λ if the read shows bk's allele, −λ if
the competitor's, 0 otherwise), ties broken by the column score. For each competitor c: d_c = columns where bk and c
both carry a base and differ; K_c = those where the read shows bk's allele; p_c = P(Binom(d_c, ε) ≥ K_c), ε = e/3;
llr_c as above. With **e = 0.003, α = 1e-3** (the shipped `--error-rate` / `--alpha`):

* **Tied** if k = 1 without the sole-candidate rule, or **some competitor has d_c = 0** (twin over the footprint).
* **Assigned** iff max_c p_c < α/(k−1) (Bonferroni over the competitors) and min_c llr_c > 0; else **Ambiguous**.
* Then: origin certificate (edits X + unaligned vs Binomial(n, e), normal tail < α) may demote to Ambiguous
  (`origin_rejected`); a sole candidate that passes it is Assigned; an Assigned read with a tie partner outside every
  family unit is demoted to Tied (`tie_outside_catalog`).
* `n_decisive` in `assignments.tsv` = m := min_c d_c (0 for every Tied-by-twin row; reset to 0 for sole candidates).

## 2. The model (derived before looking)

Read from its true copy T; per-base substitution error e_t, uniform over the 3 other bases; columns independent.
Then, per competitor, K_c ~ Binom(d_c, 1−e_t) and the read shows c's allele with prob. e_t/3 per column.

* **m = 0 ⇒ abstain with probability 1** (Tied). No depth changes this: it is identifiability, not power.
* **m = 1 ⇒ abstain with probability 1 for every k ≥ 2 and every e_t.** The best attainable p-value from one column
  is ε = e/3 = 0.001, and the gate needs p < α/(k−1) ≤ 0.001 strictly. At the shipped defaults e/3 = α exactly
  (IEEE: `0.003/3.0 == 0.001`), so a single column can never clear the gate. This is a knife-edge of two defaults,
  not a sampling property.
* **m ≥ 2:** K*(d, k) = min{K : P(Binom(d, ε) ≥ K) < α/(k−1)}; K*(2..9, 26) = 2, K*(2..45, 2) = 2. So an all-agreeing
  read always passes, and P(abstain | m = 2) ≈ 1 − (1−e_t)² ≈ 2e_t per competitor at distance 2 (0.2–0.6%);
  P(abstain | m ≥ 3) = O(e_t²). P(wrong) ≤ Σ_c C(d_c, K*)(e_t/3)^{K*} ≈ 1e-6.
* **Origin certificate under the same model:** X ~ Binom(L, e_t) with no unaligned bases; P(reject) for L = 2 kb is
  ≈ 4e-3 at e_t = e and ≈ 1e-8 at e_t = 0.001.
* **Per-read "power curve" in identity terms:** if the columns separating T from its nearest competitor fall along the
  read as a Poisson process of rate δ (1 − identity), m ~ Poisson(λ = L·δ) and P(assign) = 1 − (1 + λ)e^{−λ};
  **80% at λ ≈ 3** (e.g. a 2-kb read needs the nearest copy ≤ 99.85% identical).

Error rate used (sources): e = 0.003 is O2's own null (`copy_assign --error-rate` default, "HiFi ~0.003").
The truth simulation injects e_t = 0.001 substitutions + 0.0003 insertions + 0.0003 deletions (`bench/sim.py`
`simulate_reads(err=0.001, indel=0.0003)`). The real libraries' mismatch rate on MAPQ-60 primary reads (`=`/`X`
CIGAR, 20,000 reads, measured today as an input, `hifi_err.py`): human A119b chr16 X/(=+X) = 0.00138 (median read
0.00092), gorilla GGO NC_073242.2 0.00198 (median 0.00090). These include true sample-vs-reference differences, so
they are upper bounds on sequencing error. Predictions use e_t ∈ {0.001, 0.003}.

## 3. Data (fixed now)

Population = one row per molecule and family. Real: `primary_local = 1` and `contested = 1` (the "contested" set, memory
`project_o2_as_tied_gate`). Sim: every MAPQ-0 simulated read scored by `bench/score.py reads`, row of its TRUE family
(the OWN reading).

| id | substrate | table | truth |
|---|---|---|---|
| S-H | human A119b chr16 catalog, 1,398 copies, 28,453 reads | `rustle_figures/o2sim/human/o2.assignments.tsv` + `score/per_read_human_o2.tsv` | yes |
| S-G | gorilla NC_073244.2 catalog | `rustle_figures/o2sim/gorilla/o2.*` (30 rows) + `gw22/o2sim/gsd_o2.*` (SD catalog, 165 rows) | yes (read name `family\|copy\|i`) |
| R-H | human A119b NPIP (MCL0, 26 copies) | `bakeoff/human/ours_final2.assignments.tsv` (md5 8a057f68) + `ours_final2_dump.star_reads.tsv` | no |
| R-G | gorilla GGO NPIP (MCL1_073242), MCL7_073242, MCL58_073242 | `bakeoff/mcl1_final2` (c87c7f71), `mcl7_final` (eb34c7af), `mcl58/final2` (1934cbf4) | no |

Strata: k = `n_candidates` ≤ 1 (sole / orphan) reported apart; for k ≥ 2, m = `n_decisive` ∈ {0, 1, 2, 3–4, ≥ 5};
k bins {2, 3–5, 6–10, ≥ 11}. Outcomes: assigned · ambiguous-PSV (ambiguous, not origin-rejected) · ambiguous-origin ·
tied-outside (tied with `tie_outside_catalog` = 1) · tied. Abstain = anything not assigned.

**Seen before writing (disclosed):** the headline totals in memory (R-H contested 1,118 = 262 assigned / 512 tied /
344 ambiguous; R-G MCL1 33 = 4/28/1; MCL7 11 = 0/10/1; MCL58 40 = 2/38/0); `ours_final2.params.tsv` (origin_rejected
1,257 rows over all rows); S-H `score_reads_human_o2.txt` OWN by divergence bin (ALL 163 correct / 0 wrong / 1,084
abstain; coverage <0.5% 0.053, 0.5–1% 0.772, 1–2% 0.804, 2–5% 0.306, ≥5% 0.278). No row was tabulated by m or k.

## 4. Predictions and falsifiers

* **P0 (instrument).** A Python port of the §1 PSV decision, fed each R-H dump row's columns and the reported copy as bk,
  reproduces the table's status on ≥ 99% of dump rows with k ≥ 2 that are not origin-rejected and not tied-outside;
  and the reported copy is a maximin maximizer on ≥ 99%. *Falsifier:* < 99% ⇒ stop, the model describes a rule I misread.
* **P1 (rule).** 0 assigned rows with k ≥ 2 and m ≤ 1, in every table. *Falsifier:* any.
* **P2 (sim, per-read model).** S-H, k ≥ 2, m ≥ 2: ambiguous-PSV ≤ 2% overall, and in every m stratum with n ≥ 30 the
  observed count ≤ the 95% upper binomial bound of the model's expectation at e_t = 0.003. Wrong OWN assignments:
  expected < 0.01; *falsifier* ≥ 2.
* **P3 (sim, origin).** S-H ambiguous-origin at m ≥ 2 ≤ 1% (model: < 0.5%). *Falsifier:* > 2% ⇒ a non-error mechanism
  (locus extent, unaligned tails) that I then decompose, not a power effect.
* **P4 (sim, identity → m).** λ = read length × (1 − closest_identity of the source copy, `sim.sibling.tsv`).
  Observed P(m ≤ 1 or tied) within 0.15 of (1+λ)e^{−λ} in the λ bins [0, 0.5) and [0.5, 2); in λ ≥ 5 the observed
  P(m ≤ 1) exceeds the Poisson value by > 0.15 (already implied by the disclosed ≥ 2% coverage of 0.28–0.31: reads
  of divergent copies still tie). If the low-λ bins miss too, the advisor text must be stated in m, never in identity.
* **P5 (real, human NPIP R-H).** (a) ≥ 60% of contested abstentions have m ≤ 1 (identifiability or a single column).
  *Falsifier:* < 60% ⇒ the "not a power problem" framing fails for NPIP. (b) At m ≥ 2 the ambiguous-PSV fraction is
  ≥ 5× the model's expectation at e_t = 0.003 (real reads disagree with every candidate at PSV columns far more than
  HiFi error: sample-vs-reference variation, conversion, alignment). *Falsifier:* within the P2 band ⇒ real data behave
  like sequencing error only. (c) The per-column disagreement rate measured on m ≥ 5 dump rows (read base ≠ bk's
  allele at bk-vs-c columns), plugged into the model, under-predicts the m ∈ {2, 3–4} ambiguous-PSV fraction by ≥ 2×
  (disagreements cluster within molecules). Direction only; reported either way.
* **P6 (locus depth).** My Fisher-exact simulation of longcallR's design (two haplotypes, equal reads, P < 0.01,
  R = 7.25, symmetric usage p = √R/(1+√R) vs 1−p) crosses 80% power at a total coverage in [35, 70]×. *Falsifier:*
  outside ⇒ their unstated usage fractions differ; I report the dependence on base usage instead of one number.
  The O2 analogue (the same 2×2 with copy in place of haplotype) needs N80 **assignable** reads, i.e. a depth
  D80 = N80 / q where q is a copy's assignable fraction (unique placements + O2-assigned) — reported per R-H copy;
  a copy whose reads are all m = 0 has q = 0 and no finite D80.

## 5. Hostile self-review (before)

1. **The AS-tied gate conditions on the read's bases.** A read from T that differs from its tie partner at PSV columns
   rarely ties in AS unless its own bases split between the two copies; so reads that reach O2 with m ≥ 2 are enriched
   for exactly the reads that contradict themselves, and the iid-error model is the wrong null for them. P5(b) is the
   test of that, stated in advance.
2. **m is reported relative to bk, not to T.** For correct assignments bk = T; for abstentions bk may be a twin. m is the
   rule's own quantity, so stratifying by it is the rule's view of its evidence, not truth.
3. **Stratifying by m is partly tautological:** m = 0 is Tied by construction (P1 is a rule check, not a finding). The
   content is P(m) itself (what fraction of the abstention is identifiability) and the m ≥ 2 behaviour.
4. **Copy count k and m are confounded** (large families have more near twins). Report both margins, no regression.
5. **The sim has no allelic variation**; its m ≥ 2 behaviour is the ideal case. The real tables have no truth; their
   m ≥ 2 abstention cannot be called wrong or right, only attributed to a mechanism.
6. **The longcallR analogy is loose by design:** its 50× is phased reads at a heterozygous locus with a known effect
   size; O2's per-read abstention is not a test at all. P6 converts only the locus-level half.
7. **Sole candidates (k = 1)** are assigned without any PSV: kept apart, never mixed into the m strata.

## Amendment 1 (2026-09-29, after P0 and a population count, BEFORE any tabulation by m or k)

The prereg (sha1 6cfcae42) is bound. P0 was run (it passed; see Outcome). While fixing the population I found that §3's
real-data population is wrong: the memory's "contested" set (R-H 1,118 = 262 / 512 / 344) is **not** the TSV column
`contested` (primary MAPQ < 60). It is the binary's own stderr set **CONTESTED = AS-tied molecules that are neither
origin-rejected nor single-candidate** (R-H: 7,657 − 1,257 origin-rejected − 5,282 single-candidate = 1,118; MCL1:
636 − 603 − 0 = 33), i.e. rows with `n_candidates` ≥ 2 and `origin_rejected` = 0, all rows. Seen while establishing
this (no m or k involved): the R-H row counts by (`primary_local`, `contested`, status) — e.g. (1, 1): 153 assigned /
3,063 tied / 729 ambiguous — and the four `*.err` decompositions.

**Change:** the PRIMARY real-data population is the binary's CONTESTED set (the one the quoted abstain rate is over).
The §3 literal population (`primary_local` = 1 ∧ `contested` = 1) is reported as SECONDARY, where origin-rejected and
single-candidate rows appear as their own abstention categories. P5(a)'s 60% bar applies to the primary population;
the secondary one is reported without a bar. Predictions, strata and falsifiers are otherwise unchanged. For the
gorilla simulations the scorer's own per-read file (`o2sim/gorilla/score/per_read_gorilla_o2.tsv`) is used, like S-H;
`gw22/o2sim/gsd_o2` (no scorer file) uses the true family's rows with `contested` = 1.

## Outcome (2026-09-29)

Instruments (`/mnt/linuxdisk/tmp/rustle_figures_dev/o2_power/`, sha1 prefix): `o2_rule.py` 632e3531 (port of the PSV
decision), `o2_model.py` cff5e326, `load.py` 222e3ab6, `p0_instrument.py` fe044a11, `o2_strata.py` 977e6ac9,
`p5c_dump.py` 14594d7a, `p4_lambda.py` 19b1bb88, `p6_locus_power.py` 1d0fb012, `hifi_err.py` f8682fb2. Outputs
`model_table.txt`, `strata.txt`, `p5c.txt`, `p4.txt`, `p6.txt`. Prereg sha1 6cfcae42 (Amendment 1: 53657d25).

### The per-read decision function (model_table.txt; exact, checked by Monte Carlo on the full rule)

| m (columns vs nearest competitor) | P(abstain), e_t = 0.001 | P(abstain), e_t = 0.003 | P(wrong) |
|---|---|---|---|
| 0 | 1 (tie) | 1 | 0 |
| 1 | **1 for every k** (one column gives p = e/3 = 0.001 = α, never < α/(k−1)) | 1 | 0 |
| 2 | 0.002 (nearest only) … 0.049 (all 25 competitors at 2) | 0.006 … 0.140 | ≤ 3e-5 |
| ≥ 3 | ≤ 1e-4 | ≤ 7e-4 | ≤ 3e-4 (union bound, k = 26) |

MC of the full rule (maximin choice included) agrees with the exact values to ≤ 0.001; 0–8 wrong per 20,000 even at
e_t = 0.03–0.05. **O2's per-read decision is a step function of m: abstain at m ≤ 1, assign at m ≥ 2 unless the read's
own bases conflict.** HiFi-level error moves abstention by < 1% at m ≥ 3.

### Verdicts

| # | verdict | numbers |
|---|---|---|
| P0 | ✅ | port reproduces 952/952 R-H statuses (k ≥ 2, not origin-rejected, not tied-outside); m identical 952/952; reported copy in the maximin set 952/952 |
| P1 | ✅ | 0 assigned rows at k ≥ 2 and m ≤ 1 in all 7 tables (S-H m = 1: 30/30 ambiguous; R-H m = 1: 63/63) |
| P2 | ✅ | S-H m ≥ 2 (n = 250): PSV-ambiguous 1 (0.4%); m = 2: 1 of 32 vs expected 0.06–0.24 (95% ub 1); m ≥ 5: 0 of 193. Wrong 0 of 163 assigned (153 at m ≥ 2 + 10 sole candidates) |
| P3 | ✅ | S-H origin-rejected at m ≥ 2: 0 of 250 |
| P4 | ✅ as stated / ⚠ per bin | λ < 0.5: obs P(m ≤ 1) 0.999 vs 1.000 (n = 759); [0.5, 2): 0.719 vs 0.582 (+0.138); λ ≥ 5 pooled 0.230 vs 0.002 (+0.228) — but [5, 20) alone is +0.111 and [2, 5) +0.216; copies with no sibling hit (n = 110): 0.80 |
| P5a | ⛔ **falsified (human)** | R-H: 467 of 856 abstentions at m ≤ 1 = **54.6% < 60%**. Gorilla: MCL1 23/29 (79%), MCL7 7/11 (64%), MCL58 38/38 (100%) |
| P5b | ✅ | R-H m ≥ 2: 281 PSV-ambiguous of 651 vs expected 1.6–17.8 at e_t = 0.003 (95% ub ≤ 26) — **16–173×**; m = 2: 257/271 |
| P5c | ✅ (direction) | per-column disagreement on m ≥ 5 rows r5 = 635/42,013 = 1.51% (5–15× HiFi error); plugged into the model it expects 8.1 PSV-ambiguous at m = 2 (observed 257, 32×) and 0.03 at m 3–4 (observed 19) |
| P6 | ✅ | Fisher design reproduced: N80 = 60 at symmetric usage (0.729 / 0.271), power 0.688 at 50×; 62–128 over base usage 0.05–0.4 |

### What the human NPIP abstentions are (CONTESTED 1,118; 856 abstain; `strata.txt`, `p5c.txt`)

| mechanism | rows | share | does depth help? |
|---|---|---|---|
| m = 0, twin over the footprint (tied) | 376 + 28 tied-outside = 404 | 47.2% | no — identifiability |
| m = 1, one column (knife-edge of e/3 = α) | 63 | 7.4% | no (per read); only a per-base error below 0.003 would |
| m ≥ 2, **read splits evenly between two copies** (margin exactly 0) | 281 | 32.8% | no — contradicting evidence |
| m ≥ 2, PSV-assigned but a tie partner outside every family unit (§6gz) | 108 | 12.6% | no — an unscored competitor |
| m ≥ 2, sequencing error (model) | ≤ 18 expected, e_t = 0.003, worst config | ≤ 2% | this is the only "power" share |

**Every m ≥ 2 PSV-ambiguous row in all seven tables has margin exactly 0 (283/283):** the read shows bk's allele and
the competitor's allele at equally many columns. 257 are (d, K, J) = (2, 1, 1); **244 of the 281 come from one copy pair,
copies 6 / 7 (chr16:16.39 Mb / 18.33 Mb)**, whose reads carry copy-X's base at one PSV and copy-Y's at another ~837 bp
downstream (constant gap ± 5 bp; signatures C…G 107, C…G 99 in the two orientations). Hundreds of molecules share it, so
it is the sample's own sequence, matching neither reference copy (allelic variation at a PSV, conversion, or a reference
error — not distinguished here). This is hostile-review item 1 realised: an exact AS tie selects reads whose bases split
between copies, and O2 abstaining on them is the correct answer, not a lack of power. The truth simulation, which has no
such haplotypes, shows 1 balanced read in 250.

Secondary population (`primary_local` ∧ MAPQ < 60, R-H n = 3,945): 2,651 are single-candidate rows, 2,615 of them
tied-outside — the family has one copy there and the tie partner is outside the catalog; plus 360 origin-rejected rows.
Neither is a power effect.

### Locus depth (the longcallR connection, `p6.txt`)

The same 2×2 with copy in place of haplotype needs N80 = 60 **assignable** reads, so total depth D80 = 60 / q:

| population | q = assigned / contested | D80 |
|---|---|---|
| human NPIP, all contested | 262/1,118 = 0.234 | 256 |
| … m ≤ 1 | 0/467 | none (no finite depth) |
| … m = 2 | 3/271 = 0.011 | 5,420 |
| … m 3–4 | 6/49 = 0.122 | 490 |
| … m ≥ 5 | 253/331 = 0.764 | 78 |
| gorilla NPIP (MCL1), contested | 4/33 = 0.121 | 495 |
| gorilla MCL58 / MCL7 | 2/40 / 0/11 | 1,200 / none |

Per human NPIP copy (coarse: span-overlap primaries, abstentions attributed by the chosen copy): 11 of 26 copies have
(MAPQ-60 + O2-assigned)/primaries < 0.1, i.e. D80 > 600 reads by the conservative count; 8 copies have ≥ 0.69.

### Hostile self-review (after)

1. **P5a failed on one copy pair.** 244 of the 281 balanced reads are copies 6/7 in one human sample; without that pair
   the m ≤ 1 share would be 467/612 = 76%. The failure is honest as stated, but it is one locus, and "identifiability"
   vs "contradiction" is a split of the same non-power abstention.
2. **The balanced reads' origin is not established.** Allelic variation at a PSV, gene conversion in A119b, and a CHM13
   assembly error all give the same reads. No DNA of A119b was used.
3. **The sim agrees with the model because it was built with the model's assumptions** (iid substitutions, no
   haplotypes). It validates the rule's arithmetic and the 0-wrong claim, not the real-data behaviour.
4. **m is measured against bk,** which for a balanced read is one of the two copies arbitrarily; m is symmetric for the
   pair, so the stratum does not move.
5. **The D80 numbers use q from the AS-tied set only**; loci also carry clear-best reads the aligner places, so the
   family-level D80 is the hard-read worst case, and the per-copy table is coarse (other genes inside copy spans; bk
   attribution overshoots copy 6).
6. **longcallR's 50× is not exactly reproduced** (60 at symmetric usage; their usage split is unstated). The conversion
   D80 = N80 / q does not depend on which N80 is used.
7. **P4's identity form is only a lower bound on abstention** above λ ≈ 2 (PSVs cluster; reads fall in shared stretches),
   so the advisor sentence is stated in m.

### Advisor paragraph

Copy assignment decides one read at a time from the positions inside that read where the candidate copies differ, so
its abstentions are not a sequencing-depth problem and more reads at a locus would not reduce them. Under the shipped
rule (per-base error 0.003, family-wide false-assignment level 0.001 split over the competitors), a read is never
assigned from zero or one differing position — one position can at best reach a p-value of 0.001, which is exactly
the level — and is assigned from two or more unless its own bases disagree with its best copy; with HiFi-like errors
(measured mismatch ≤ 0.0014 in human, 0.002 in gorilla) that disagreement happens for under 1% of reads, and a wrong
assignment is rarer still (model bound ≤ 3 in 10,000; none among 163 simulated assignments). On a simulation with known origins (human chr16, 1,247 reads the aligner placed at mapping quality 0) this
held exactly: 163 assigned, 0 wrong, 89% of abstentions at zero or one differing position. On real human NPIP reads, 55% of
the 856 abstentions have zero or one differing position, 33% are reads whose bases support two copies equally (244 of
them from one copy pair whose reads match neither reference copy), and 13% have an equally good placement outside the
family; sequencing error accounts for about 2% at most. The longcallR figure (80% power above ~50× at their effect size SOR = 2,
an odds ratio of about 7.25; we reproduce 60×) applies only to locus-level questions built from assigned reads: a copy-specific splicing test
at that effect size would need about 60 assignable reads, i.e. about 256 multi-mapping NPIP reads in human and 495 in
gorilla at the observed assigned fractions, and no depth at all suffices for copies whose reads never span a differing position.
