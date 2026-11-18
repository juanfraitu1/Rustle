# Row 1769: does a chain score fitted across species beat the frozen two-sample score? (pre-registered 2026-11-17)

User goal (2026-11-17): with several samples and species at hand, look for more empirical rules that raise precision.
Registered before any feature table of the nine ape samples was read. Nothing here changes a default; a PASS makes a candidate.

**Why this and not a new single rule.** The single read-derived rules are closed (rows 1616-1619, 1677, 1723, 1725-1732, 1766;
rule E of row 1499 for template switching and A-runs). What is new is the number of labelled samples: the frozen chain score
(row 1677, default on `PREFIX.final.gtf` since row 1761) was fitted on OR6737 + A119b only. Rules must not need a second sample
at analysis time (user 2026-11-16): a fixed weight file applied to one sample complies.

**Samples (11).** Default-pool assemblies of today (`artifact_refresh_2026-11-16`: `g.def`, `h.def`, `apes/STEM/def`): A119b,
OR6737, KB3781, and testis + reference-individual cell line of PTR, PPA, PPY, PAB. Truth = exact multi-exon RefSeq chain of the
sample's own reference. Features = `bench/famcmp/chain_rules_feat.py` (the 23 of the frozen score), standardised within the
sample as the shipped filter does. HG002 is not included (its default assembly was not finished).

**Part A, descriptive, no bar.** The frozen score on each sample (correct kept, precision of kept); per feature, the precision
of the chains the frozen filter keeps, by feature quantile, per sample.

**Part B, registered: leave-one-species-out refit.** For each of the six species, fit the same logistic model (same 23
features, same fit routine) on the pooled default chains of the other five species' samples, and score the held-out species'
samples once. Comparison at matched yield: cut the refit score so that it keeps exactly as many correct chains as the frozen
filter keeps on that sample; compare the precision of the kept chains.
- PASS: precision gain >= +0.005 on at least 9 of 11 samples and no sample below -0.002.
- Otherwise FAIL: more samples do not improve the linear score; closed.
Disclosed bias: the frozen score is in-sample on OR6737 and A119b, which favours it there.

**Part C, registered: additive non-linear score with three added genome-only features.** Same protocol and the same bar
against the frozen score, and reported against Part B. Each of the 23 features is replaced by its within-sample decile
indicators (an additive model, fitted by the same routine with the same ridge), plus: `ptc` and ORF length / transcript
length (`tools/gtf_ptc.py`), count of GC-AG or AT-AC introns in the chain, and log of the shortest intron.
- PASS / FAIL as in Part B. If B fails and C passes, C is the candidate; if both pass, the one with the larger median gain.

**New versus lost** (standing rule): for a passing part, report per sample the chains it keeps that the frozen filter drops
and the reverse, each with its correct count, not the net.

**Limits.** RefSeq is predicted on the apes and incomplete (row 1750: unannotated orangutan introns are real at 0.90), so
"wrong" includes real unannotated chains; on A119b the result is also given for curated (NM / NR) chains alone. Every sample
is spent for earlier rules; the standard here is transfer to a species absent from the fit.

## Result (2026-11-17, `r1769_multisample_rules/loso.txt`, jobs 4478260 features, loso rerun after a script fix)

Precision of kept chains at the frozen filter's correct count; B and C fitted without the sample's species.
| sample | default chains / correct | frozen kept / precision | B precision (gain) | C precision (gain) | C keeps, frozen drops / correct | frozen keeps, C drops / correct |
|---|---|---|---|---|---|---|
| A119b | 236,301 / 40,086 | 197,147 / 0.1995 | 0.1948 (-0.0047) | 0.2207 (+0.0212) | 7,519 / 355 | 26,476 / 355 |
| OR6737 | 76,139 / 26,425 | 63,262 / 0.3985 | 0.4062 (+0.0077) | 0.4498 (+0.0513) | 3,491 / 716 | 10,700 / 716 |
| KB3781 | 74,288 / 26,316 | 63,223 / 0.3878 | 0.4262 (+0.0384) | 0.4817 (+0.0939) | 2,674 / 1,105 | 14,996 / 1,105 |
| PTR testis | 57,039 / 21,904 | 47,243 / 0.4468 | 0.4500 (+0.0032) | 0.4771 (+0.0303) | 3,116 / 507 | 6,120 / 507 |
| PTR AG18354 | 92,517 / 31,155 | 79,452 / 0.3704 | 0.4049 (+0.0346) | 0.4513 (+0.0809) | 3,710 / 1,255 | 17,958 / 1,255 |
| PPA testis | 58,953 / 18,545 | 50,342 / 0.3495 | 0.3699 (+0.0204) | 0.4198 (+0.0703) | 2,262 / 546 | 10,697 / 546 |
| PPA PR00251 | 102,121 / 22,152 | 83,221 / 0.2611 | 0.2594 (-0.0017) | 0.3104 (+0.0494) | 4,208 / 263 | 17,449 / 263 |
| PPY testis | 113,109 / 21,347 | 92,410 / 0.2255 | 0.2281 (+0.0026) | 0.2615 (+0.0360) | 3,988 / 283 | 16,700 / 283 |
| PPY AG05252 | 80,604 / 24,366 | 66,736 / 0.3524 | 0.3686 (+0.0162) | 0.4184 (+0.0660) | 2,972 / 521 | 13,499 / 521 |
| PAB testis | 75,338 / 16,931 | 59,967 / 0.2759 | 0.2657 (-0.0102) | 0.2889 (+0.0129) | 4,258 / 258 | 6,945 / 258 |
| PAB AG06213 | 73,064 / 25,253 | 61,431 / 0.3871 | 0.4183 (+0.0312) | 0.4730 (+0.0860) | 2,474 / 892 | 13,637 / 892 |

- **Part B: FAIL** (gain >= +0.005 on 6 of 11, minimum -0.0102, median +0.0077). More species do not improve the linear score.
- **Part C: PASS** (11 of 11, minimum +0.0129, median +0.0513). Candidate, not a default.
- Part A (same file): among frozen-kept chains, `ptc` chains match RefSeq 0.04-0.13 against 0.26-0.51 for the rest; a chain with
  no full-chain read 0.05-0.20 against 0.22-0.46; retention-flagged 0.09-0.19; `mq0_share` and `l_cont5` carry nothing.

## Ablation (post hoc, not registered; `abl.txt`, job 4478414)
C0 = decile model of the 23 features only; C1 = C0 + non-GT-AG count + shortest intron; C2 = C1 + `ptc` + ORF share (= C).
| | gain >= +0.005 | minimum | median |
|---|---|---|---|
| C0 | 9 of 11 | +0.0031 | +0.0164 |
| C1 | 10 of 11 | +0.0038 | +0.0180 |
| C2 | 11 of 11 | +0.0129 | +0.0513 |
About a third of C's gain is the non-linear form alone; two thirds come from the two ORF features. The ORF features move
correct chains from non-coding to coding RefSeq records at the same total: correct predicted non-coding (XR) chains kept,
frozen -> C0 -> C2: OR6737 872 -> 835 -> 635, KB3781 885 -> 796 -> 577, PTR AG18354 1,821 -> 1,731 -> 1,275, PAB AG06213
776 -> 729 -> 528; A119b curated non-coding (NR) 3,983 -> 3,993 -> 3,892 and XR 1,591 -> 1,596 -> 1,544, while curated coding
(NM) 26,196 -> 26,181 -> 26,244. C2 therefore conflicts in spirit with the user's decision of 2026-11-16 that `ptc` is an
attribute and not a filter, and its gain against a predicted, ORF-built annotation is partly circular. C1 has neither problem.
Gains are larger on the cell lines (C0 +0.045 to +0.050 on KB3781, PTR AG18354, PAB AG06213) than on testis.

## Part D, registered 2026-11-17 before any HG002 feature table of today's assembly was made: confirmation on HG002
HG002 (Kinnex, CHM13, default-pool assembly of 16 Nov, `apes/HG002_kinnex.t2t/def/x.gtf`) took no part in any fit above. Fit C1
and C2 once on all 11 samples, freeze the weights (`r1769_multisample_rules/model_C1.json`, `model_C2.json`), score HG002 once.
At the frozen filter's correct count on HG002:
- C1 CONFIRMED if precision gain >= +0.005; C2 likewise, reported with the correct chains kept per RefSeq class.
- A confirmed C1 is the candidate to offer the user as an alternative weight file for the chain score filter; C2 is offered only
  with its non-coding cost stated. Neither becomes a default without the user.
Disclosed: HG002 was the held-out sample of row 1761 for the frozen score, so the frozen score is not in-sample there either.

## Part D result (2026-11-17, `confirm.txt`, job 4478422): both CONFIRMED on HG002
HG002 default chains 164,142, correct 36,727, precision 0.2238; frozen filter keeps 137,826 with 35,632 correct, 0.2585.
| model (fitted on all 11) | kept at 35,632 correct | precision | gain | keeps, frozen drops / correct | frozen keeps, model drops / correct | correct kept NM / XM / NR / XR (frozen 25,096 / 7,022 / 3,047 / 445) |
|---|---|---|---|---|---|---|
| C1 (no ORF features) | 120,003 | 0.2969 | +0.0384 | 3,881 / 471 | 21,704 / 471 | 25,222 / 6,965 / 2,995 / 431 |
| C2 (with `ptc`, ORF share) | 109,952 | 0.3241 | +0.0655 | 4,558 / 570 | 32,432 / 570 | 25,364 / 7,062 / 2,788 / 395 |
Curves, precision at 0.97 / 0.95 / 0.9 / 0.8 of the correct default chains: frozen 0.259 / 0.274 / 0.315 / 0.405; C1 0.298 / 0.330 /
0.381 / 0.463; C2 0.324 / 0.362 / 0.423 / 0.506.
C1 is the candidate: +0.038 precision on an unseen sample at the same correct count, slightly more curated coding chains, 52
fewer curated non-coding. C2 costs 259 curated non-coding chains (8.5% of them) for a further +0.027. Weights:
`r1769_multisample_rules/model_C1.json`, `model_C2.json` (design = `designC9` in `confirm.py`; not ported, not a default).

## Port (2026-11-17, user: "yes lets port the no-orf one"): opt-in weight file, default unchanged
`RUSTLE_CHAINSCORE_MODEL=r1769` (or a path to a model file) selects `tools/data/chainscore_r1769.json` for the `chainscore` stage
and for the filter on `PREFIX.final.gtf` inside `assemble`; unset = `r1677`, as before. Code: `design_decile` in
`bench/famcmp/chain_rules.py`, `intron_extras` and the `design` switch in `bench/famcmp/chainscore.py`, `chainscore_model` in
the driver; driver tests extended (25 pass).

**Cut.** TAU 0.074 = the largest cut (grid 0.001) at which the 11 fit samples together keep at least as many correct default
chains as the r1677 filter (263,590). A fixed cut does not match the correct count per sample (`port.txt`):
| sample | correct kept r1677 (share) | r1769 (share) | precision r1677 -> r1769 |
|---|---|---|---|
| HG002 (unseen) | 35,632 (0.970) | 35,484 (0.966) | 0.2585 -> 0.3048 |
| A119b | 39,328 (0.981) | 37,347 (0.932) | 0.1995 -> 0.2653 |
| OR6737 | 25,211 (0.954) | 25,947 (0.982) | 0.3985 -> 0.3853 |
| KB3781 | 24,518 (0.932) | 25,576 (0.972) | 0.3878 -> 0.3999 |
| PTR testis | 21,106 (0.964) | 21,417 (0.978) | 0.4468 -> 0.4401 |
| PTR AG18354 | 29,426 (0.945) | 30,236 (0.971) | 0.3704 -> 0.3916 |
| PPA testis | 17,593 (0.949) | 18,011 (0.971) | 0.3495 -> 0.3526 |
| PPA PR00251 | 21,725 (0.981) | 20,724 (0.936) | 0.2611 -> 0.3502 |
| PPY testis | 20,843 (0.976) | 20,029 (0.938) | 0.2255 -> 0.2696 |
| PPY AG05252 | 23,516 (0.965) | 23,661 (0.971) | 0.3524 -> 0.3703 |
| PAB testis | 16,546 (0.977) | 16,137 (0.953) | 0.2759 -> 0.3093 |
| PAB AG06213 | 23,778 (0.942) | 24,692 (0.978) | 0.3871 -> 0.3910 |
Rows other than HG002 are in-sample. The cut is harder on noisy samples and gentler on clean ones; three samples keep under
0.95 of their correct chains, the floor the r1677 cut was set by. The cut is one number in the model file and is the user's call.

**Parity** (`parity.sbatch`, job 4478469): the driver stage with `RUSTLE_CHAINSCORE_MODEL=r1769` on HG002 drops 47,739 chains,
the same 47,739 as the analysis script; the default model through the changed code drops 26,316, as before the change.

## Default since 2026-11-18 (user: "lets make 0.055 default")
- `tools/data/chainscore_r1769.json` is the default model of the chain score filter, cut 0.055 (was opt-in at 0.074).
  `RUSTLE_CHAINSCORE_MODEL=r1677` gives the earlier default (cut 0.061); `RUSTLE_CHAINSCORE=0` still switches the filter off.
- Why 0.055: loosest cut on a 0.005 grid at which every one of the 12 samples keeps at least 0.95 of its correct default
  chains (lowest A119b 0.955; at 0.074 A119b 0.932, PPA PR00251 0.936, PPY testis 0.938).
- Kept share of correct chains, precision of kept, r1677 at 0.061 -> r1769 at 0.055: A119b 0.981, 0.199 -> 0.955, 0.242;
  HG002 0.970, 0.259 -> 0.978, 0.280; OR6737 0.954, 0.399 -> 0.989, 0.373; KB3781 0.932, 0.388 -> 0.983, 0.385;
  PTR testis 0.964, 0.447 -> 0.988, 0.426; PTR AG18354 0.945, 0.370 -> 0.982, 0.376; PPA testis 0.949, 0.349 -> 0.983, 0.341;
  PPA PR00251 0.981, 0.261 -> 0.956, 0.319; PPY testis 0.976, 0.226 -> 0.960, 0.251; PPY AG05252 0.965, 0.352 -> 0.983, 0.351;
  PAB testis 0.977, 0.276 -> 0.972, 0.289; PAB AG06213 0.942, 0.387 -> 0.987, 0.377. Eleven of these samples were used to fit
  the model; HG002 was not.
- Parity (job 4480508, HG002): the driver with no environment variable drops 36,101 of 164,142 chains, identical to the
  analysis list; `RUSTLE_CHAINSCORE_MODEL=r1677` drops 26,316 as before. `cargo test --test pipeline_driver_dispatch chain`: 2 passed.
- Not rerun with the new default: the published comparison pages and the advisor package (their `final.gtf` is the r1677 filter).
