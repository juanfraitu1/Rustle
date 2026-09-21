# Pre-registration — empirical rules against FALSE NEGATIVES

**Written 2026-09-21 before any rule is scored.** User: *"we have been checking if the definition has
false positives and determining empirical rules to remove them, however we have not been focusing on
false negatives, can we find some easy empirical rules that can help diminish them too?"*

The observation is correct: §6t3–§6u3 were almost entirely precision work (the exon conjunct, containment
rejection, operator precision). This is the recall side.

## Where the false negatives actually are — measured first

Every within-family truth pair on held-out chr2/chr8/chr10 (263 pairs), traced through the pipeline:

| stage | pairs | of all | **of the 88 FNs** |
|---|---|---|---|
| recovered (same cluster) | 175 | 66.5% | — |
| ⭐ **edge kept, but MCL SPLIT them** | **40** | 15.2% | **45.5%** |
| no alignment at all | 23 | 8.7% | 26.1% |
| identity/cov gate dropped it | 20 | 7.6% | 22.7% |
| exon conjunct dropped it | 5 | 1.9% | 5.7% |

⭐ **The largest single bucket is a GROUPING loss, not an edge loss.** For 40 of 88 false negatives the
edge survives every filter and MCL then separates the two members — e.g. ANAPC1 ~ ANAPC1P1, ~ANAPC1P4,
~ANAPC1P5 are all edged and all split. Those 40 are recoverable without touching edge construction,
which is the cheapest place to look.

## The rules to test — post-clustering merges, each one line

For every pair of MCL clusters (A,B) in the same component, merge them if:

| rule | statistic |
|---|---|
| **R1 — edge count** | ≥ *m* shipped edges run between A and B |
| **R2 — edge fraction** | edges(A,B) / min(\|A\|,\|B\|) ≥ *f* |
| **R3 — neighbourhood Jaccard** | mean J_N over the A–B edges ≥ *j*, using §6u3's metric |

Swept: m ∈ {2,3,5,10}, f ∈ {0.25,0.50,0.75,1.0}, j ∈ {0.10,0.20,0.30,0.50}.

⚠ **R3 must abstain on small components.** §6u3 measured J_N at **AUC 0.500 — exact chance — at component
size 2**, and 60% of truth families are pairs. R3 is therefore only applied where both clusters sit in a
component of ≥ 5 nodes, fixed now, and that restriction is part of the rule, not a tuning knob.

## What is reported, and the bar

Each rule reports **FN recovered**, **FP introduced**, and the **pooled F** through the same scorer, on
the same held-out chromosomes and Soto truth as every other arm.

Baseline: shipped MCL through the port, **F 0.7123** (sens 0.838 / prec 0.715).

| outcome | verdict |
|---|---|
| F ≥ 0.7123 + 0.02 | ⭐⭐ **A REAL RECALL RULE** |
| F ≥ 0.7123 and sensitivity up by ≥ 0.02 | ⭐ **RECALL GAINED FREE** |
| F within ±0.02 | ⚠ **A TRADE, not a gain** — report the exchange rate |
| F < 0.7123 − 0.02 | ⛔ **NO** |

⚠ Declared now: a rule that raises sensitivity while lowering F is **not** a false-negative fix, it is a
precision sale, and I will report it as such with the exchange rate (FN recovered per FP introduced). The
bar is F, because §6t3 measured that precision is the axis this definition lives on.
