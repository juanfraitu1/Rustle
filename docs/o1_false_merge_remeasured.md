# The false-merge rate, re-measured under the 2026-08-20 defaults

**Status 2026-08-20. Re-run of the frozen 150-window negative panel with the current binary.**
Harness: `/home/juanfra/winloci_scratch/o1neg/` (`run2.sh`, `score.py`, seed 101, panel unchanged).
July outputs preserved as `out.jul17` / `dump.jul17` and re-scored for the comparison.

## 0. ⚠ A provenance correction

`ONE_METHOD.md` was updated on 2026-08-19 to say the rates "were measured on the shipped 494-family
catalog". **That is wrong for the false-merge rate.** It is measured on **HUMAN CHM13 v2.0**, over 150
gene-tight single-locus windows drawn from 1,630 eligible, with **A119b** IsoSeq — not on the GGO
catalog at all. It is also a **specificity and a LOWER bound**, not a precision: the panel has no
positive stratum, so there is no prevalence and no precision to compute.

## 1. The headline is unchanged

| | July (pre-guard) | 2026-08-20 defaults |
|---|---:|---:|
| **windows emitting ≥1 family** | **2/150 = 0.0133** [0.0037, 0.0473] | **2/150 = 0.0133** [0.0037, 0.0473] |
| the two windows | W063 `ZNF492`, W106 `ANKHD1` | **the same two** |
| co-membership assertions | 4 | 4 |
| self-overlap (false by coordinates) | 1 | 1 |

**Everything shipped on 2026-08-19/20 — the transcript-orientation guard as default, one path with
refine rejected, substrate typing, λ recomputation, the streamed PAF — leaves the false-merge
specificity exactly where it was.** The definition's measured error rate is stable under all of it.

## 2. But the spurious-edge burden collapsed

Same panel, same scorer, `dump/*.edges.tsv`:

| | July | now | change |
|---|---:|---:|---|
| E_r edges emitted on the negative panel | **28** | **3** | **−89%** |
| edges between OVERLAPPING spans | **26** | **1** | **−96%** |
| self-identity CERTIFIED false | **7** | **1** | −86% |

This is the orientation guard doing exactly what the coordinate analysis predicted: overlapping
antisense pairs are ~46× enriched among minus-only edges, and on this human negative panel they were
**25 of the 28** edges. The guard removes them at the edge layer.

## 3. Why the family rate did not move — and why that is the right outcome

Those 26 spurious edges were already being stopped downstream: only 2 of 28 ever reached `copies.tsv`.
The guard cleaned a layer that was not the binding constraint on the emitted catalog.

And the two survivors were **never orientation artifacts**, so no version of this guard could have
fixed them:

* **W063 `ZNF492`** — self-identity: the aligned block IS the 1,204 bp genomic intersection of the two
  node spans, identity **exactly 1.000000**. That is census **pathology (a)**, one locus emitted as two.
* **W106 `ANKHD1`** — the 206 bp linking node is **100% soft-masked** interspersed repeat.

Both were predicted not to move: pathology (a) is a coordinate signature, not a strand one.

## 4. ⚠ A rate that rose while the burden fell

The **certified-false edge rate** went **7/28 = 0.2500 → 1/3 = 0.3333** while the absolute count fell
**7 → 1**. The denominator shrank faster than the numerator. Quoting the rate alone would say the edge
layer got *worse*; quoting the count alone would miss that what remains is proportionally more
concentrated. **Report both, or neither.**

## 5. What this does and does not establish

* ✅ The false-merge specificity is **stable at 1.33%** across a substantial change of defaults — the
  strongest evidence yet that it is a property of the definition and not of an invocation.
* ✅ The orientation guard's precision benefit is now **measured on an independent human negative
  panel**, not only on the GGO FP arm it was derived from: **28 → 3 edges**.
* ❌ It is **still a specificity and a lower bound**. No positive stratum, no prevalence, no precision.
* ❌ It does **not** transfer to gorilla. This panel is human CHM13/A119b. The GGO 627-family catalog
  has no comparable negative panel.
