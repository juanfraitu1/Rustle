# An RNA-level family definition scored against RNA-DERIVED truth

**Result: connected components of the L3 edge graph at `w_98 >= 0.995` reach bipartite R = P = F = 0.833
with pairwise sensitivity 1.000 against the Dishuck Iso-Seq groups — an RNA-level truth derived
independently of our alignment gate.** The shipped cut (0.98) gives F 0.222 on the same truth.

## Why the earlier numbers were so low — the truth, not the definition

§6o8 measured a pairwise recall CEILING of 0.052 (human) / 0.229 (gorilla) against the guided truth. §6o9
found the cause: that truth is an MCL over annotated GENE SPANS (genomic, intron-containing), and **only
6.7% of its within-family pairs align at all as spliced RNA** (5.0% pass the edge gate, against a 5.2%
ceiling). Two alternative explanations were refuted first — the node→truth mapping (ceiling still 0.047 at
99% overlap) and node fragmentation (edge rate flat 20-33% across every completeness stratum).

**Scoring an RNA-level definition against a DNA-level truth demands pairs that do not exist as RNA.**

## Against RNA-derived truth (Dishuck NPIP), two independent views

| L3 cut | iso: R | P | **F** | sens | prec | | lit: R | P | **F** | sens | prec |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **0.980 (SHIPPED)** | 0.222 | 0.222 | **0.222** | 1.000 | 0.157 | | 0.667 | 0.692 | **0.679** | 0.910 | 0.557 |
| 0.985 | 0.667 | 0.667 | **0.667** | 1.000 | 0.500 | | 0.630 | 1.000 | **0.773** | 0.462 | 1.000 |
| 0.990 | 0.667 | 0.667 | **0.667** | 1.000 | 0.500 | | 0.593 | 1.000 | **0.744** | 0.387 | 1.000 |
| **0.995** | **0.833** | **0.833** | **0.833** | **1.000** | 0.667 | | 0.519 | 1.000 | **0.683** | 0.276 | 1.000 |
| 0.999 | 0.722 | 1.000 | 0.839 | 0.500 | 1.000 | | 0.185 | 1.000 | 0.312 | 0.070 | 1.000 |

`iso` = Dishuck Iso-Seq groups (5 multi-copy groups, 18 nodes); `lit` = Dishuck subfamilies (2 groups,
27 nodes). Mean F across the two views: shipped 0.98 → **0.451**; 0.985 → 0.720; 0.990 → 0.706;
**0.995 → 0.758**.

**Both views agree the shipped 0.98 is too low.** They disagree on the exact optimum (iso 0.995,
lit 0.985), so the defensible statement is **0.985-0.995**, not a single value.

## The headline definition

> **Family = connected component of the L3 edge graph (t1 AND f_ex >= 0.30 AND gap-excluded identity
> w_98 >= 0.995).** On the Iso-Seq truth: **bipartite R 0.833, P 0.833, F 0.833; pairwise sensitivity
> 1.000, precision 0.667.**

Recall is perfect at the pairwise level — every true Iso-Seq pair is recovered — and precision is the
remaining axis.

## Caveats that must travel with this

- **NPIP is a development family** and the cut was chosen on these truths. The trend is monotone and large
  (F 0.222 → 0.833), not a knife-edge, but it is not held-out.
- **Small**: 18 nodes / 5 groups (iso), 27 nodes / 2 groups (lit).
- The two views disagree on the optimum, which is why the recommendation is a range.
- Earlier in this run an RNA truth was derived from alignability under the same gate that builds the edges;
  that is **circular** and was discarded. The Dishuck truths are not — they come from published Iso-Seq
  grouping, independent of anything computed here.

**Next: a held-out RNA-derived truth** (TBC1D3 Guitart groups, or Dishuck on a second family) to confirm
0.985-0.995 without choosing it on the reported substrate.
