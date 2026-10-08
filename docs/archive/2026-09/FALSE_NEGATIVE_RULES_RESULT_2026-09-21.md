# False negatives: where they are, and why the obvious rules do not transfer

Run 2026-09-21 against `docs/PREREG_false_negative_rules_2026-09-21.md` (committed `5723b4ff` before
any rule was scored). Held-out chr2/chr8/chr10, Soto S1C truth, same scorer as every other arm.

## The decomposition — the durable result of this run

Every within-family truth pair traced through the pipeline (263 pairs, **88 false negatives**):

| stage | pairs | of all | **of the 88 FNs** |
|---|---|---|---|
| recovered (same cluster) | 175 | 66.5% | — |
| ⭐ **edge kept, but MCL SPLIT them** | **40** | 15.2% | **45.5%** |
| no alignment at all | 23 | 8.7% | 26.1% |
| identity/cov gate dropped it | 20 | 7.6% | 22.7% |
| exon conjunct dropped it | 5 | 1.9% | 5.7% |

⭐ **Two facts worth carrying forward.** First, the largest single bucket is a **grouping** loss: for
45.5% of false negatives the edge survives every filter and MCL then separates the members
(ANAPC1 ~ ANAPC1P1, ~ANAPC1P4, ~ANAPC1P5 are all edged and all split). Second, **48.8% are
edge-construction losses** — no alignment at all, or killed by the identity/coverage gate — which is
exactly §6o8's priority #1. Only **5.7%** are the exon conjunct, the filter that has absorbed most of
this session's attention.

## The rules, and the test-set trap they fell into

Three one-line post-clustering merges over pairs of MCL clusters. Scored on the held-out set, best
threshold per rule:

| rule | F | sens | prec | FN pairs recovered | FP pairs added | exchange |
|---|---|---|---|---|---|---|
| baseline MCL | 0.7016 | 0.808 | 0.722 | — | — | — |
| R1 edge-count ≥ 3 | 0.7067 | 0.828 | 0.730 | 6 | 24 | 0.25 : 1 |
| R2 edge-fraction ≥ 0.75 | 0.7165 | 0.842 | 0.730 | 14 | 24 | 0.58 : 1 |
| **R3 neighbourhood Jaccard ≥ 0.2** | **0.7222** | **0.852** | 0.715 | **14** | **2** | **7 : 1** |

R3 looked like the first clean recall win of the session: **+0.0206 F, sensitivity +0.044, 14 false
negatives recovered for 2 false positives.** It also validated §6u3 operationally — J_N used as a
post-clustering merge, abstaining below component size 5.

⛔ **Then the threshold was selected honestly, on chr16, as every other arm this session has done:**

| threshold | chr16 (development) F | held-out F |
|---|---|---|
| baseline | 0.6698 | 0.7016 |
| r3 = 0.2 | **0.5626** | 0.7222 |
| r3 = 0.3 | 0.6484 | 0.7114 |
| **r3 = 0.5 (best on development)** | **0.6698** | **0.7016** |

**The development-selected threshold is a no-op: +0.0000 on held-out.** The value that helps the held-out
set (0.2) is the *worst* on chr16, costing it 0.107 F and dropping precision 0.691 → 0.541. **No threshold
helps both.** The +0.0206 was entirely test-set selection.

⭐ The mechanism is visible: chr16 is NPIP/PKD1P/SMG1 territory where families genuinely interrelate, so
merging on shared neighbours over-merges there, while chr2/chr8/chr10 families are more separable. This
is r344's verdict in a new place — *"no operating point does both"*.

## What this means for false-negative work

- ⛔ **Post-clustering merges are not the answer.** All three rules buy recall by selling precision, and
  the only favourable exchange rate does not transfer.
- ⭐ **The decomposition says where to look instead**: 48.8% of false negatives never get an edge at all
  (no alignment, or the identity/coverage gate), which is edge construction — §6o8's stated priority #1,
  and the one axis this session has not touched.
- ⚠ **The exon conjunct is not the problem.** It accounts for 5 of 88 false negatives (5.7%), which is
  consistent with §6t7's finding that its rejections are the right call.

## Limits

- 263 truth pairs, 88 false negatives, three chromosomes; individual buckets are 5–40 pairs.
- One development chromosome (chr16) for threshold selection, and it is an unusually family-dense one —
  which is precisely why it disagrees with the held-out set, but a second development chromosome would
  make the transfer test stronger.
