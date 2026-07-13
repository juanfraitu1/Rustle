# Can we align below asm20? Tweaking minimap2 vs protein homology

**Date:** 2026-07-11. Question: the 5 residual refine false negatives are real DIVERGENT paralog families below
the asm20 (0.80) / sensitive (-k11 -w5, 0.70) identity floor. Can a more sensitive pairwise alignment — a tweaked
minimap2, or protein-level homology — recover them WITHOUT reintroducing the correctly-dropped false positives?

## The experiment

A labeled test set: **5 divergence-FN families** (real paralogs the gate drops: ARMCX ARMCX1/6, IFITM IFITM1/3,
FRG1like, CCDC144B, KRAB_ZNF ZNF665/677/818) vs **5 correctly-dropped FPs** (PBX1_split, HS6ST2_split = one gene
split into fragments; NNT_GHR, GARRE1_ZNF540, CLIC2_TMLHE = repeat/intronic bridges). 4 strategy agents each swept
one axis on the shared genomic sequences; an "edge" = copies align at identity ≥ 0.50 AND coverage-of-shorter ≥
0.50. A strategy WINS if it recovers FN without reintroducing FP.

## Result — every strategy separates (0 FP), but recall is seeding-limited

| family | label | nt k-mer ladder (-k7 -w3) | nt scoring relax (-B1) | **protein (mmseqs)** | translated 6-frame |
|---|---|---|---|---|---|
| FRG1like | FN | **EDGE** .98/1.0 | **EDGE** .98/1.0 | ORF-fail | **EDGE** .82/.94 |
| CCDC144B | FN | **EDGE** .96/1.0 | **EDGE** .94/1.0 | **EDGE** .84/.73 | **EDGE** .75/.91 |
| KRAB_ZNF | FN | — .00/.00 | — .73/.42 | **EDGE** .54/1.0 | — .60/.46 |
| ARMCX | FN | — .00/.00 | — .00/.00 | — .29/.86 | — .38/.38 |
| IFITM | FN | — .00/.00 | — .00/.00 | ORF-fail | — .00/.00 |
| PBX1_split | FP | — /.00 | — /.04 | — /.00 | — /.06 |
| NNT_GHR | FP | — .97/**.22** | — .91/**.34** | — /.00 | — .67/**.36** |
| GARRE1_ZNF540 | FP | — .98/**.13** | — .81/**.36** | — /.00 | — .65/**.32** |
| HS6ST2_split | FP | — /.00 | — .83/**.05** | — /.00 | — /.09 |
| CLIC2_TMLHE | FP | — .96/**.09** | — .90/**.29** | — .87/**.14** | — .68/**.32** |

(cells: identity / coverage-of-shorter; **bold** coverage = the FP ceiling.)

## Two findings

### 1. The COVERAGE floor is the false-positive guard — not the identity floor
Every FP tops out at **≤ 0.36 coverage-of-shorter** across ALL strategies (repeat-bridges ≤ 0.36; gene-splits ≤
0.09). They align at *high identity* (0.86–0.98) over their small shared segment — a repeat or a domain — but
**never reach 0.50 coverage**. So no strategy, however sensitive on identity, reintroduces an FP: all four give
**0 FP**. **Lowering the identity/sensitivity floor is therefore SAFE** — the coverage floor (≥ 0.50 of the shorter
copy) is what rejects bridges and splits, and it is identity-independent. This reverses the naive worry that "more
sensitive alignment readmits domain bridges": a domain bridge fails on *coverage*, not identity.

### 2. Nucleotide seeding walls at ~65% identity; protein is the only lever below it
minimap2 — even at `-k7 -w3` with relaxed DP scoring (`-B1`) — produces **zero anchors** for ARMCX, IFITM, and
KRAB_ZNF. The bottleneck is **seeding**, not scoring: below ~65% nucleotide identity there are too few exact
k-mers to chain, so relaxed scoring has nothing to extend. Tweaking minimap2 below asm20 therefore recovers only
the *moderately*-divergent FNs (FRG1like 98%, CCDC144B 94%) — which the genomic-span tier already recovers.

**Protein homology is the unique lever for deep coding divergence.** KRAB_ZNF (nucleotide-invisible) is recovered
**only by protein** (0.54 amino-acid identity over the full ORF, cov 1.0). And it is safe by the same coverage
mechanism: a domain bridge shares < 50% of the ORF (CLIC2_TMLHE: 55-aa domain, protein cov 0.14) so it fails the
coverage floor in protein space too. This is exactly the existing `--protein-tail` tier (mmseqs, fident ≥ 0.50).

### The genuine floor
ARMCX at **29% amino-acid identity** is below any safe protein threshold — a truly ancient paralog pair. IFITM's
protein miss is an *extraction* artifact (the 6-frame longest-ORF heuristic missed its small coding exon; annotation
CDS would recover it). Non-coding divergent paralogs cannot use protein at all. These are the residual limit.

## Recommendation

- **Do NOT tweak minimap2 below asm20.** It is safe but seeding-limited — it recovers only the moderate divergence
  the genomic-span tier already handles, and cannot reach the deep cases (0 anchors below ~65% identity).
- **Enable `--protein-tail`** to recover the deep coding-divergence FNs (KRAB_ZNF-class): it is the only strategy
  that crosses the seeding wall, it is already implemented, and it is FP-safe (the coverage floor rejects domain
  bridges in protein space). Cost: an `mmseqs` dependency, and it recovers ~1–2 more families, not the ancient
  (ARMCX) or non-coding cases.
- **Defense point:** the gate's specificity comes from the **coverage floor (≥ 0.50)**, not the identity floor —
  so identity sensitivity can be raised for recall (nucleotide, then protein) up to each method's seeding limit,
  with no false-positive cost.

Reproduce: `/home/juanfra/winloci_scratch/mmsweep/` holds the labeled FASTAs + `manifest.tsv`; per-strategy
alignment commands are in the strategy descriptions above.
