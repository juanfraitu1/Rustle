# The corroboration layer — external certificates for `E_r`

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog. The current default emits **627 families / 2,019 copies**. See
> [`NUMBERS.md`](NUMBERS.md) and [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

**Status 2026-08-20. All candidates are T8 — offline, not through the shipped binary — except the
transcript-orientation guard, which is the default.** Companion to
[`o1_false_positive_rules.md`](o1_false_positive_rules.md) (the veto side) and
[`o1_ledger.md`](o1_ledger.md) (the one-page verdict list).

## 1. Why a layer, and not a better rule

Two measurements force the architecture.

**Every statistic computed from the two nodes is length in disguise.** The FP class *is* the
short-node class, so an internal statistic half-works and then collapses when length-normalised —
which is what happened to coverage on all four substrate × denominator cells, to block count
(0.7548 → 0.6410 per kb), to junction crossing (100× exon-count bias) and to read tiling.

**And the discriminating clause does not order the classes at all.** True GFPT1×GFPT2 scores
**0.5353**, false ATP1A1×ATP4A scores **0.5689** — a true pair below a false one, so no threshold on
it separates them.

⟹ **The object is one high-sensitivity RELATION plus low-sensitivity, high-precision CERTIFICATES that
mark which edges independent evidence corroborates.** Every attempt to *replace* the relation with a
certificate failed, because they are not the same kind of object: a relation must reach everything, a
certificate must be right.

## 2. The certificates

Unit = **pair**, GGO, FP arm 14 / TP arm 150. ⚠ The FP arm is **not held out** — this line of work was
derived on it — so the TP column is the load-bearing one and every 0/14 carries **Wilson95
[0.0000, 0.2153]**: bounded, not proven zero.

| certificate | references | FP admitted | TP covered | portable? |
|---|---|---:|---:|---|
| exon-sum coverage ≥ 0.50 *(the relation, for contrast)* | nothing external | 14/14 | **150/150** | — |
| genomic-span coverage ≥ 0.50 | the genome | **0/14** | 30/150 = 0.2000 | ✅ genome only |
| SD containment ≥ 0.90 | an SD catalog | **0/14** | 24/150 = 0.1600 | ❌ species-specific |
| **flank homology** | the genome | **0/14** | 30/150 = 0.2000 | ✅ **genome only** |
| **UNION of the three** | — | **0/14** | **41/150 = 0.2733** | |

### 2a. Flank homology — the portable one

Take **5 kb outside each locus, exclude the gene bodies entirely**, and ask whether the extra-genic
sequence aligns (**≥ 300 bp at ≥ 0.90 identity**). A segmental duplication copies a *neighbourhood*, so
a real pair shares flanks; a mobile element inserted into two unrelated genes does not.

| min bp | min id | FP admitted | TP kept |
|---:|---:|---:|---:|
| 300 | 0.90 | **0/14** | 30/150 |
| 1000 | 0.80 | 1/14 | 36/150 |
| 2000 | 0.80 | **0/14** | 27/150 |

It reaches SD containment's specificity **using only the genome**, which removes the portability
objection that stopped SD containment from being folded in. ⚠ **A negative is weak** — a duplication
whose boundaries fall inside the gene, or an old one with diverged flanks, looks identical. **The
positive is the informative direction.**

### 2b. SD containment — what it licenses

**Sayable:** for **136/627 = 21.7%** of families every member locus lies entirely inside one segmental
duplication at ≥ 0.90 — *the family and the duplicated segment are the same event.*
**Not sayable:** which copy is ancestral. With one genome every duplication edge is symmetric, so this
is `DNA_DUPLICATION`, never `DERIVED_FROM`. ⚠ Nor is later **gene conversion** excluded — containment
dates the **segment**, not the family.

### 2c. The vetoes, for completeness

| veto | references | FP rejected | cost |
|---|---|---:|---|
| **transcript orientation** — SHIPPED DEFAULT | strand | 6/14 | 4/9,032 edges; **2.05%** genome-wide |
| genome-anchored repeat multiplicity ≥ 50 | the genome | 10/12 | 0/135 arms; **3.67%** genome-wide |

⚠ The repeat veto is a **veto, never an admission criterion**: TP median `gmult` is 2, so a criterion
strict enough to exclude repeats excludes most real paralogues.

## 3. Limits, stated plainly

* **They overlap heavily.** 30 + 24 + 30 = 84 if disjoint; **41** in fact. All three are genomic
  co-duplication evidence, not three independent lines — do not present them as convergent proof.
* **Coverage is a quarter.** 27.3% of true pairs are corroborated; the rest are **uncorroborated, not
  rejected.** The layer adds confidence where it fires and says nothing where it does not.
* **T8.** None of the certificates has run through the shipped binary. The orientation *veto* has.
* **0/14 is bounded, not zero**, and the arm is not held out.

## 4. Reproduce

```bash
python3 bench/o1_flank_extension.py          # flank homology
python3 bench/o1_sd_anchor.py 0.90           # SD containment, identity floor as argv
python3 bench/o1_sd_certificate.py           # duplication-mechanism classes on a catalog
python3 bench/o1_substrate_denominator.py    # genomic-span coverage
```
