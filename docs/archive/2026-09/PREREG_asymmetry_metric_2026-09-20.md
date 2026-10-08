# Pre-registration — a pair metric that tolerates length asymmetry without saturating

**Written 2026-09-20 before any new metric is scored.** User: *"find a better jaccard-like metric that
does not penalize length asymmetry and see if it helps with family definition."*

## The two horns, both measured

| | metric | failure |
|---|---|---|
| **r359** | Jaccard `m/(la+lb−m)` | **penalises** asymmetry — EEF1A1 retrocopies score J ≈ 0.07 at core identity 0.92 |
| **r913** | containment `m/min(la,lb)` | **saturates** on asymmetry — a short gene inside a long one scores ~1.0, long genes become hubs (864 chr8 pairs ≥ 1.0; largest component 171 vs 23; precision 0.164 vs 0.639) |

A usable metric needs asymmetry *tolerance* without saturation. `cov_longer = m/max(la,lb)`, which O1
already gates on at 0.30, is on the Jaccard horn.

## The candidates — all parameter-free except E3, all on the same best PAF record per pair

| arm | metric | why |
|---|---|---|
| **B** | `m / (la + lb − m)` — Jaccard | the incumbent, already scored **0.6044** |
| **E1** | **`m / sqrt(la · lb)`** — **Ochiai / cosine** | the **geometric mean of the two containments**; equals containment when `la = lb`, and unlike containment cannot reach 1.0 unless the alignment covers both genes. This is the principled answer to the two horns and the one I expect to win |
| E2 | `2m / (la + lb)` — Dice | symmetric like Jaccard but milder; included to separate "milder penalty" from "geometric shape" |
| E3 | containment `m/min(la,lb)`, **guarded** by `min(la,lb)/max(la,lb) ≥ 0.25` | r913's metric with an explicit anti-hub guard; the only arm with a parameter, included because it is the obvious engineering fix and should be beaten on its own terms |

Families = **connected components** of pairs scoring ≥ t — identical machinery to §6t3 arm B, so the
metric is the only thing that changes. Swept over **t ∈ {0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70}**.

## Comparators, already measured on this exact substrate, truth and scorer

| arm | pooled held-out F |
|---|---|
| A — shipped: MCL + exon conjunct | **0.7016** |
| A− — MCL alone | 0.6394 |
| **B — Jaccard, components** | **0.6044** |
| D — chained, containment | 0.2132 |

## The bar — committed now

Let **F_E** be a candidate's pooled held-out F with **t selected on chr16**, never on the held-out set.

| outcome | verdict |
|---|---|
| F_E ≥ 0.6394 | ⭐⭐ **MAJOR** — a one-line metric matches MCL-alone; report for adoption |
| 0.6244 ≤ F_E < 0.6394 | ⭐ **IMPROVED** — better scalar, worth having |
| 0.5844 ≤ F_E < 0.6244 | ⚠ **NEUTRAL** — the horn is real but this does not fix it |
| F_E < 0.5844 | ⛔ **WORSE** |

Secondary, committed now: **the metric must also raise the score of a known asymmetric true pair**. I
will report each metric's value on the r359 case class — pairs whose length ratio `min/max ≤ 0.5` that
are in the same Soto family — as *recall on asymmetric true pairs at the selected t*. A metric that wins
overall while losing there has not fixed the stated problem and will be reported as such.

⚠ I will not re-tune γ or E3's 0.25 guard, drop a chromosome, or change the truth after seeing scores.
If more than one arm clears the bar I will report all of them and prefer the parameter-free one.
