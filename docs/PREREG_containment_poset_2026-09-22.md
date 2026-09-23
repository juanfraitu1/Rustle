# Pre-registration — is a POSET a better family object than a partition of an undirected graph?

**Written 2026-09-22, §6y0, before any family is emitted.** User: *"is there any data structure that would
fit the spliced reads"* → *"yes that sounds good"* (test the antichain-as-family definition).

## What r1020/r1021 established, and what is still only a hypothesis

**Measured (r1020), chr16 de novo:** the containment relation is a **strict partial order** — antisymmetry
**0 violations**, transitivity **92.9%**, longest chain **38**, only 59.3% of container pairs nesting (a
**DAG, not a forest**). It predicts same-referee-family at **AUC 0.943** against a length-ratio null of
**0.396**, holding in every length band, so it is not a size artefact. **183 of 1,192 elements lie below
≥2 incomparable containers** — multi-membership a partition cannot express.

**Hypothesis (r1021), untested:** §6s9 (a partition cannot hold a fusion), §6u8/§6u9 (Soto's truth is a
cover every scorer drops) and §6t7/r1002 (containment is asymmetric) are ONE limitation — *the output is a
partition of an undirected graph while the data is a poset*. **That a poset SCORES better is exactly what
this file tests. Nothing above licenses it.**

⚠**This is NOT §6u8's refuted cover test.** That one PATCHED a finished MCL partition (node v joins
cluster C iff v has ≥ k edges into C) and was refuted on two truths. Here the cover is not a patch: it
**falls out of the structure**, because an element below two incomparable maximal elements is in two
down-sets by construction.

## The edge set is HELD FIXED — only the structure varies

Both arms consume the **shipped `--dump-graph` output**, which has already passed every shipped conjunct
(identity ≥0.70, `cov_longer` ≥0.30, alen ≥300, `--min-exonic-bp 1 --min-shared-exon-frac 0.60`).
Direction for each surviving edge is then read off the PAF. **No edge is added or removed by either arm**,
so any difference is attributable to the structure and to nothing else.

⚠**Comparator is `mcl_port` MCL on that same graph, never the shipped Rust F** (register 917: `mcl_port`
is not bit-identical to the Rust MCL — 245 vs 168 clusters on chr2). The shipped numbers are context only.

## Candidate family objects — fixed now, all standard poset objects, no new thresholds

| # | object | why it is a candidate |
|---|---|---|
| **P1** | **principal down-set ⇓m of each maximal element m** | the natural cover: an element below two incomparable maximals joins both |
| **P2** | **maximal antichains** (Mirsky/Dilworth) | a family as a set of mutually non-containing copies — the "same level" reading |
| **P3** | weakly connected components of the comparability graph | the NULL. If P3 ≈ P1 the poset added nothing beyond connectivity |
| **M** | `mcl_port` MCL on the identical graph | the comparator |

Each run both **raw** and under **transitive closure** (the measured relation is 92.9% transitive; closing
it is a provable operation and the 7.1% are plausibly gate noise — reported separately, never pooled).

`C`, the containment floor, is swept **on development only** over {0.50, 0.60, 0.70, 0.80, 0.90}.

## Truth and metrics

**Truth = the COVER**, built by `bench/cover_and_jn_definition.py`'s existing machinery: Soto S1C's
`No. Assigned Families` column (149/2,334 gene IDs are multi-family) and the protein-family referee.
⚠**Scoring a cover prediction against a partition-restricted truth is rigged** (§6u8's own warning).

- **PRIMARY: pairwise precision / recall / F** — the only metric that is well defined for a cover: a pair
  `(a,b)` is predicted together iff they co-occur in **any** family.
- **SECONDARY: one-to-one bipartite F** — reported for comparability with every prior number, and
  explicitly flagged as **penalising a cover by construction**.
- **Guards:** largest family size (r913's hub signature) and number of families. A down-set definition
  with a 38-long chain could emit one family containing everything; that is a failure, not a win.

## Substrates

Development: **chr16** (de novo and guided). Held out: **chr2 / chr8 / chr10**, no re-tuning, verdict taken
there. ⚠Never pooled with gorilla.

## The bar — committed now

| outcome | verdict |
|---|---|
| a poset object beats **M** on held-out **pairwise F** by ≥0.02, largest family ≤2× M's, and beats **P3** | ⭐ **THE STRUCTURE IS THE ANSWER** — pursue a poset definition |
| beats M by 0.00–0.02, or beats M but not P3 | ⚠ **PARTIAL** — the gain is connectivity, not order |
| no poset object beats M on held-out, **or** the winner's largest family exceeds 2× M's | ⛔ **NO** — the relation is real (r1020) but the object is wrong |

**Predicted, before looking — ⚠ PARTIAL, and I expect BOTH naive objects to fail in opposite directions.**
With a longest chain of 38, **P1's down-sets should over-merge** (a top element drags in most of a
component); **P2's antichains should under-merge**, because two genuine copies at different truncation
levels are comparable and would be split apart. My honest expectation is that r1020's relation is right and
these two off-the-shelf objects are the wrong cut of it — which, if it happens, is a result about *which*
poset object a gene family is, not a refutation of r1021.

I will not change the objects, the edge set, the truth, the metrics or the bar after seeing any number.

---

# OUTCOME (2026-09-22) — ⚠ **PARTIAL. The gain is CONNECTIVITY, not ORDER — and development did not transfer.**

Edge set held fixed (shipped `--dump-graph`), comparator `mcl_port` on the identical graph.

## Development — chr16 de novo, where it looks decisive

| object | Soto cover F | referee F | fams | largest |
|---|---|---|---|---|
| M (`mcl_port`) | 0.516 | 0.464 | 275 | 37 |
| **P1 down-sets (C=0.70)** | **0.714** | **0.524** | 289 | 37 |
| P2 antichains | 0.074 | — | 80 | 16 |
| P3 components (NULL) | 0.137→0.475 | 0.407 | 261 | 64 |

On development P1 beats M by **+0.198** and the null by **+0.239**, at the *same* largest family (37).
C is insensitive again (0.697–0.714 across 0.50–0.90). ⭐**Transitive closure is a no-op to three decimals
everywhere** — because a down-set traversal IS a reachability computation, so P1 already closes
transitively. A useful internal consistency check, not a finding.

## Held out — chr2 / chr8 / chr10 at the development-selected C = 0.70

| cell | M | **P1** | P3 (null) |
|---|---|---|---|
| chr2 · Soto cover | 0.502 | **0.546** | 0.520 |
| chr8 · Soto cover | 0.296 | 0.337 | **0.432** |
| chr10 · Soto cover | **0.558** | 0.522 | 0.519 |
| chr2 · referee | **0.417** | 0.385 | 0.407 |
| chr8 · referee | **0.934** | **0.934** | 0.896 |
| chr10 · referee | 0.206 | 0.236 | **0.246** |

⚠**P1 beats M in 3 of 6 cells, loses 2, ties 1 — and it does NOT beat P3**, which wins three cells
outright. Per the bar's row *"beats M but not P3"*, the verdict is ⚠ **PARTIAL: the gain is connectivity,
not order.**

⛔⛔**The development margin over the null collapsed on held-out**: P1 − P3 was **+0.239** on chr16 and is
**+0.026 / −0.095 / +0.003** on the three held-out chromosomes. **This is exactly what holding a substrate
back is for**, and it is the second time this session a chr16-selected margin failed to transfer (r933 was
the first).

## What survives, and what does not

⭐**Survives — r1020's measurement.** Containment IS a strict partial order (antisymmetry 0, transitivity
92.9%, AUC 0.943, not a size artefact). Nothing here touches that.

⛔**Does not survive — r1021's inference.** "The data is a poset, therefore a poset object is the right
family" is **not supported by scoring**. Plain **weakly connected components of the containment relation
(P3) do as well as the ordered object (P1)**, so the DIRECTION — the very thing r1021 argued an undirected
graph could not express — **is not what carries the signal.** What carries it is *which pairs are related
at all*, which an undirected graph represents perfectly well.

⛔**P2 antichains fail outright** (F 0.025–0.184 everywhere), exactly as predicted: two genuine copies at
different truncation levels are comparable, so an antichain splits them apart.

## Prediction scorecard

I predicted ⚠ PARTIAL with **both** naive objects failing — P1 over-merging and P2 under-merging.
**Half right.** P2 under-merged as predicted. **P1 did NOT over-merge** (largest family 37, identical to
M's, and 17–57 held-out); it failed for a different reason — it is not distinguishable from its own null.
⭐**Predicting the right verdict for the wrong mechanism is not a successful prediction**, and the null arm
(P3) is the only reason the distinction was visible at all. It was worth including.
