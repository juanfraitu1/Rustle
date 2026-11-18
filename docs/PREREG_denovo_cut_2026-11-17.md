# Row 1774: can the chain score's cut be set per sample without an annotation? (pre-registered 2026-11-17)

User (2026-11-17): the advisor is against using an annotation for guidance; "is there any way to establish a threshold de novo?"
Background: one fixed cut lands at a different operating point per sample (r1677 at 0.061 keeps 0.932-0.981 of the correct
default chains; r1769 at 0.074 keeps 0.932-0.982). Registered before any estimate was computed. Nothing here changes a default.

**Target.** Per sample, the loosest cut that still keeps 95% of the correct default chains (what the r1677 cut was set by on
its fit samples). RefSeq is used to grade the outcome only.

**Arm A, decoys.** Decoy = a default chain with at least one intron that is not GT-AG (annotation-free; not a feature of the
r1677 score). With N(t) chains and D(t) decoys kept at cut t, and two shipped constants pi = P(decoy | wrong), rho = P(decoy | correct):
`D = pi * W + rho * C`, `N = W + C` gives the estimated correct count `C(t) = (pi * N - D) / (pi - rho)`. The cut is the largest
t (grid 0.001) with `C(t) >= 0.95 * C(0)`. pi and rho are the pooled values of the samples of the OTHER species
(leave-one-species-out), i.e. constants that would ship with the model file; the sample itself contributes no labels.
**Arm B, mixture.** Two-component Gaussian mixture on logit(score) of the sample's default chains (EM, started at the lower and
upper score quartile); `C(t)` = summed posterior of the higher component above t; same rule for the cut.
Both arms on the r1677 score. Arm A on the r1769 score is reported for information only (its score uses the non-GT-AG count).

**Outcome**, per sample: the share of truly correct default chains kept at the chosen cut (target 0.95), against the fixed cut.
- PASS (per arm): share within 0.95 +- 0.015 on at least 10 of 12 samples AND largest deviation from 0.95 smaller than the fixed
  cut's largest deviation. A PASS makes a candidate (opt-in per-sample cut), not a default.
- Otherwise FAIL.
Also reported: estimated against true precision at the chosen cut; pi and rho per sample (the stability the method rests on).

**Limits.** All 12 samples are spent for the scores themselves. GC-AG and AT-AC introns are real, so rho > 0; decoys are assumed
equally likely at every score within a class, which the score's features (intron support) can violate. RefSeq on the apes is
predicted and incomplete, so "wrong" includes real unannotated chains.

## Result (2026-11-17, job 4478581, `r1774_denovo_cut/cut.txt`): both arms FAIL
Share of correct default chains kept, target 0.95, r1677 score:
- fixed cut 0.061: within 0.015 on 5 of 12, largest deviation 0.031 (0.932-0.981).
- Arm A (decoys): unusable. Within 0.015 on 0 of 12; cuts at 0.000-0.012 on nine samples (keeps everything), 0.395 on PPA PR00251
  (keeps 0.48) and no valid cut on A119b and HG002. Reason: a non-GT-AG intron is carried by 0.066-0.091 of the CORRECT chains
  (GC-AG and AT-AC are real, and a chain has many introns) against 0.086-0.161 of the wrong ones, so pi - rho is 0.02-0.08 and
  differs between samples (0.020 PAB testis, 0.078 HG002); the estimate divides by it.
- Arm B (mixture): within 0.015 on 5 of 12, largest deviation 0.049 (keeps 0.953-0.999), worse than the fixed cut; its precision
  estimates are far too high (0.47-0.96 against a true 0.17-0.44): the score's two groups are not Gaussian on the logit scale.
r1769 score: decoy arm 0 of 12 (information only), mixture 0 of 12 (keeps 0.98-1.00).
Closed for this decoy definition and this mixture. Not tested: a decoy restricted to introns that are none of GT-AG, GC-AG,
AT-AC (would bring rho near 0); that would be a new, post hoc arm.

**Post hoc arm (user: "ok lets run it"; job 4478586, `motif.py`): not possible.** A decoy defined as an intron that is none of
GT-AG, GC-AG, AT-AC does not exist in the output: 0 such default chains on all 12 samples (Rustle emits only those three motif
pairs). Every "non-GT-AG" chain of arm A is a GC-AG or AT-AC chain, which explains rho of 7-9%.
