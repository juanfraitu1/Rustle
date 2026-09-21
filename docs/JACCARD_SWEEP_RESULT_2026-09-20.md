# "A simple Jaccard sweep should work" — tested head-to-head

Run 2026-09-20 against `docs/PREREG_jaccard_sweep_2026-09-20.md` (committed `0a38c1d3` before either arm
was scored). Same substrate, same external truth (Soto S1C, ≥ 3 members), same scorer
(`bench/heldout_family_score.py`) for every arm. Held-out chromosomes: chr2, chr8, chr10.

## He is substantially right

| arm | pooled held-out **F** | sens | prec | exact |
|---|---|---|---|---|
| **A — shipped**: MCL + `--min-shared-exon-frac 0.60` | **0.7016** | 0.808 | **0.722** | 3 |
| A− — MCL alone, exon conjunct OFF | 0.6394 | 0.814 | 0.586 | 3 |
| **B — his rule**: connected components of `J ≥ t`, nothing else | **0.6044** | 0.713 | 0.639 | 2 |

**A one-line rule reaches 86% of the shipped rule's F.** Arm B has no MCL, no exon conjunct, no core
step and no read corroboration — it is `J = nmatch / (len_a + len_b − nmatch) ≥ t`, union-find, done.
The pre-registered verdict is **⚠ CLOSE** (−0.097, against a −0.10 boundary).

⚠ **Fairness matters here and it cuts both ways.** Arm B's best threshold on the held-out data is
t = 0.30 (F 0.6109) — but picking it there fits B to its own test set. Selecting t on chr16 exactly as
arm A's 0.60 was selected on chr5/7/21 gives **t = 0.50 → 0.6044 held out**, so B's threshold transfers.
Conversely, at the t = 0.05 that r772 fixed for *edge recovery*, B collapses to **F 0.246** — that
threshold was tuned for a different objective and is not a fair grouping threshold.

Register **r772** already supports him upstream: the objection *"Jaccard cannot replace pairwise
alignment for O1's edges"* was **refuted on measurement** — J ≥ 0.05 recovers 94.1% of O1 edges, J ≥ 0.01
recovers 99.8%, because O1's own `cov_longer ≥ 0.30` makes both criteria global-coverage criteria.

## But the decomposition says *which* complication to drop

| what it buys, pooled held-out | ΔF |
|---|---|
| the clustering **operator** (MCL over connected components) | **+0.0350** |
| the **exon conjunct** (`--min-shared-exon-frac 0.60`) | **+0.0622** — 64% of the gap |

⭐ **The conjunct is a pure precision device and it is nearly free**: precision **0.586 → 0.722
(+0.136)** for a sensitivity cost of **−0.007**. It is not machinery — it is a one-line biological
requirement that the pair's best record cover ≥ 60% of the *smaller gene's exonic length* with
exon-to-exon evidence, which is what stops a co-duplicated neighbour riding one shared base of flanking
sequence into the family.

## What to say to him

1. **He is right about the operator.** MCL buys +0.035 over union-find on a threshold. On this evidence
   the clustering algorithm is close to interchangeable, and §6o8/§6o9 independently found grouping is
   **saturated** — recall is bounded by the truth's divergence, not by the operator. Defending MCL as
   sophistication is not supported.
2. **He is wrong that the sweep alone suffices** — but the missing ingredient is not complexity. It is
   one conjunct worth twice what the operator is worth, and it buys precision, the axis a family
   definition lives or dies on.
3. ⭐ **The synthesis is simpler than what ships**: Jaccard threshold + the exon conjunct + connected
   components would keep the precision win and drop MCL entirely. That is a smaller rule than the
   current one and is the obvious thing to test next.

## The synthesis I proposed was tested, and it FAILED

The decomposition above suggested an obvious move: keep the cheap conjunct, drop the expensive operator.
Tested as arm C (`bench/jaccard_plus_exon.py`, J ≥ 0.50 AND f_ex ≥ 0.60, connected components):

| arm | F | sens | prec |
|---|---|---|---|
| A — shipped: MCL + conjunct | **0.7016** | 0.808 | 0.722 |
| B — Jaccard alone | 0.6044 | 0.713 | 0.639 |
| **C — Jaccard + conjunct, no MCL** | **0.5869** | **0.606** | 0.656 |

⛔ **C is worse than B.** Adding the conjunct to connected components *costs* 0.018 F, where adding it to
MCL *gained* 0.062. Sensitivity is what breaks: 0.713 → 0.606.

⭐ **The reason vindicates MCL in a way the decomposition alone did not show.** A hard pairwise conjunct
deletes edges. MCL's flow can route around a deleted edge — two members still reach each other through
a third — whereas connected components cannot: the deleted edge is the only path, and the family splits.
**So MCL is not +0.035 of operator sophistication; it is what makes the conjunct affordable.** The two
are not separable components whose gains add, and the decomposition table above must not be read as if
they were.

⚠ This is not an artifact of my approximation. `jaccard_plus_exon.py` computes f_ex as an **upper bound**
on the true shared-exon fraction (it projects the aligned interval onto each gene's exons independently
rather than requiring a base be exonic on both sides), so arm C sees a **more permissive** conjunct than
`mcl_families` applies — and sensitivity still collapsed.

## Limits

- 27 Soto families over three chromosomes; one assembly, one annotation.
- The conjunct's contribution is measured by turning it off inside `mcl_families` (`f_ex = 0.0`), so it
  is measured **in the presence of MCL**. Whether it contributes the same +0.062 on top of plain
  connected components is the untested step, and is exactly the synthesis arm above.
- Arm B is scored only at thresholds on the committed grid; a finer sweep could move it by a little.
