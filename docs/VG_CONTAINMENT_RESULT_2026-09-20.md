# Would a variation graph solve the containment problem? — No, and the number says why

Run 2026-09-20 against `docs/PREREG_vg_containment_2026-09-20.md` (committed `de0f6949` before any
multiplicity was computed). Tool: `bench/vg_multiplicity_containment.py`. Population: the **76 pairs
the shipped rule rejects at containment ≥ 0.90** on held-out chr2/chr8/chr10, both endpoints
Soto-labelled — **19 TRUE, 57 FALSE**.

## The test, and why it needs no VG

A variation graph's distinctive contribution over a pairwise alignment is **multi-way**: how many other
sequences traverse a segment. That is node multiplicity, and the existing all-vs-all PAF already
determines it. **If multiplicity does not separate the classes, the main extra information a VG carries
does not separate them, and building one cannot help.**

## Result: the signal is real, in the right direction, and far too weak

| | TRUE (n=19) | FALSE (n=57) |
|---|---|---|
| median multiplicity | **4.0** | **6.0** |

| mult ≤ k | TRUE kept | FALSE kept | precision | recall |
|---|---|---|---|---|
| 0 | 2 | 2 | 0.500 | 0.105 |
| 3 | 6 | 9 | 0.400 | 0.316 |
| **5** | **11** | **23** | **0.324** | 0.579 |
| 10 | 19 | 40 | 0.322 | 1.000 |
| any | 19 | 57 | 0.250 | 1.000 |

**Best precision retaining ≥ 10 of the 19 TRUE pairs: 0.324**, against a baseline of 0.250 and a
pre-registered bar of **0.60 for "VG HELPS"** and 0.40 for PARTIAL. ⛔ **Pre-registered verdict: NO.**
The one point above 0.40 (mult ≤ 3) keeps only 6 TRUE pairs.

⭐ **AUC = 0.681.** Register **r384** measured minimizer multiplicity as a *general* edge separator at
**AUC 0.686** — an independent population, an independent question, and the same answer to three decimal
places. Multiplicity is a ~0.68-AUC signal wherever it is measured. That is not an edge-admission
criterion; it is a weak prior.

## What this does and does not settle

⛔ **Settled: the multi-way "how many share this segment" channel does not solve containment.** It is the
channel people mean when they say a VG would disambiguate repeats, it is exactly what a VG makes
cheap to read, and at AUC 0.68 it cannot carry an admission decision that pairwise identity
(0.987 vs 0.989) already failed to carry.

⚠ **Not settled: path topology.** A VG carries more than multiplicity — it carries the *order and
structure* of shared segments. This test says nothing about that channel, and it is the same
**structure** direction §6t7 arrived at independently (matching intron chains, exon-to-exon evidence).
If a VG is to help here, that is where the argument has to be made, and it would have to beat the exon
conjunct, which is already a structural test and already rejects these pairs.

⚠ Also unchanged: **r472** — a VG cannot *define* families (it presupposes its members), and nothing
here proposes that. This was a graph-derived signal for edge admission, not a definition.

## Limits

- 76 pairs, three chromosomes, one truth (Soto). "FALSE" means "different Soto family"; a few could be
  real relationships Soto omits, which would raise precision without changing the AUC much.
- Multiplicity is computed as distinct other genes overlapping ≥ 50% of the long gene's aligned
  interval. A different overlap floor moves the counts; the AUC is a ranking statistic and is stable to
  that choice.
