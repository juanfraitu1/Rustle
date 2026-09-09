# PREREG — a whole-distribution test of subfamily structure (2026-09-09, before any run)

## Why
The largest-gap statistic (`bench/identity_gap.py`) uses one number out of n pairs. It has power on deep gaps
(18 gorilla families at p < 0.01, gaps 0.022–0.089) but NPIP's gap is 0.0014 (p = 0.072) even though the
A/B partition is visibly present in the identity distribution (every A–A ≥ 0.9759, every A–B ≤ 0.9744).

## Two statistics, both using the WHOLE distribution
1. **Silverman's critical-bandwidth test** (Silverman 1981): find the smallest Gaussian-KDE bandwidth at which
   the within-family identity distribution has ONE mode; bootstrap from that KDE (with Silverman's variance
   correction) and ask how often a unimodal sample needs a bandwidth at least as large. Non-parametric,
   tests unimodality directly, no threshold chosen.
2. **1- vs 2-component Gaussian mixture, ΔBIC** (EM): ΔBIC = BIC₁ − BIC₂; ≥ 10 is Kass–Raftery "very strong".
   Parametric, cheap, a second opinion.

Both are added to `bench/modality.py`; `identity_gap.py` and `gw_subfamily_scan.py` report them alongside
the gap p. **The gap test is not removed** — three statistics, all reported.

## Sensitivity set (already known, gorilla)
The 18 gorilla families at gap-p < 0.01 are treated as POSITIVES; the families at gap-p_worst > 0.80 (the top
quartile) as NEGATIVES. This is circular only if the new statistics are tuned on them — they are not; they have
no free parameters beyond the bootstrap count.

## Predictions
| # | prediction | refuted by |
|---|---|---|
| P1 | Silverman reports multimodal (p < 0.05) on **≥ 15 of the 18** gap-positives | < 12/18 |
| P2 | Silverman reports unimodal (p ≥ 0.05) on **≥ 80 %** of the gap-negatives | < 65 % |
| P3 | ΔBIC ≥ 10 on ≥ 15/18 positives; ΔBIC < 10 on ≥ 80 % negatives | as above |
| **P4** | **Human NPIP**: Silverman **p < 0.05** — the A/B separation is visible in the whole distribution even though its gap is shallow | p ≥ 0.05 ⟹ the structure is real in the partition but NOT statistically separable from one mode, and that is the honest statement |
| P5 | **Human TBC1D3**: unimodal on both statistics | any multimodal call |
| P6 | **Gorilla NPIP** (gw MCL12): unimodal on both, consistent with the gap result | multimodal |
| P7 | Silverman and ΔBIC agree on ≥ 85 % of the 244 gorilla families | < 70 % |

## ⛔ What must be reported regardless
If P4 fails, the NPIP subcluster claim stands ONLY as "the components at an identity floor reproduce the
published partition" (§6gw) — with no statistical certificate. That is a demonstration, not a detector.

## Rules
Human and gorilla never pooled. Gene symbols read out only. No parameter is tuned on any family.

---
## Amendment 1 (2026-09-09, written AFTER P4 failed and BEFORE the genome-wide partition run)
⚠⚠ **The partition permutation test below is POST-HOC.** It was designed after seeing that Silverman gave
human NPIP p = 0.125 despite a perfect A/B partition (every A–A ≥ 0.9759, every A–B ≤ 0.9744). It is
reported as such: its human-NPIP result (p = 0.006) is a demonstration that the 1-D tests measure the wrong
thing, NOT a pre-registered confirmation. Its genome-wide calibration is pre-registered here, now:
| # | prediction | refuted by |
|---|---|---|
| P8 | partition p < 0.05 on **≥ 15/18** gap-positives | < 12 |
| P9 | partition p ≥ 0.05 (or no valid partition) on **≥ 80 %** of the 25 gap-negatives | < 65 % |
| P10 | gorilla NPIP (MCL12): partition p < 0.05, matching its Silverman 0.045 | p ≥ 0.10 |
**Statistic**: max over threshold-induced component partitions of (mean within − mean between identity),
components ≥ 3 members. **Null**: permute identities across edges (identical marginal). No parameter tuned on
any family; `MIN_COMP = 3` was fixed before any run.
