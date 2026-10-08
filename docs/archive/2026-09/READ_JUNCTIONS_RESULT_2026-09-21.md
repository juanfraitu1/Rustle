# Read-derived junctions for the junction family rule — coverage 5×, F worse

Run 2026-09-21 against `docs/PREREG_read_junctions_2026-09-21.md` (committed `66f33f34` before any read
junction was extracted). Substrate: **A119b.t2t.bam** on held-out chr2/chr8/chr10 (5.12M / 2.21M / 2.41M
mapped reads), junctions from CIGAR `N` at `-F 2308`, ≥ 3 supporting reads. Everything downstream is
§6u1's rule unchanged.

## The census, which is what the hypothesis was about

**At the GENE level, reads add almost nothing.**

| | pooled over the three chromosomes |
|---|---|
| genes | 7,020 |
| spliced by annotation | 5,652 (**80.5%**) |
| spliced by reads | 4,203 (59.9%) |
| **spliced by reads but NOT by annotation** | **35** |
| both-spliced fraction of shipped homologous edges | **48.6% → 48.7%, +1 edge** |

⭐ **The intronless genes are genuinely intronless.** 35 of the 1,368 gain a junction. Reads see *fewer*
spliced genes than the annotation does (59.9% vs 80.5%) because expression and coverage are the limit,
not annotation quality. §6u1's ceiling is biological — processed pseudogenes are retrotransposed from
mRNA and have no introns by construction — and no better annotation will lift it.

**At the JUNCTION level, reads add a great deal.**

| | pooled |
|---|---|
| annotated junctions | 36,552 |
| read-confirmed (± 10 bp) | 28,733 (**78.6%**) |
| **read-only, in genes that were already spliced** | **54,794 (150% of the annotated count)** |

So the two censuses disagree, and that is the interesting part: reads do not create new spliced *genes*,
they add 2.5× as many *junctions* to the genes that were already spliced.

## And the extra junctions make the definition worse

Union of annotated and read junctions, same rule, same k sweep:

| k | F (reads ∪ annot) | coverage | **§6u1 annotated only** |
|---|---|---|---|
| 1 | 0.2769 | **33.3%** | **0.4513** / 6.5% |
| 2 | 0.2342 | 25.1% | 0.3915 / 4.5% |
| 3 | **0.3080** | 20.2% | 0.3287 / 3.5% |
| 5 | 0.2414 | 13.6% | 0.2187 / 1.6% |

⛔ **Pre-registered verdict: WORSE.** Coverage rises ~5× (6.5% → 33.3%), exactly as the junction census
predicted, but precision collapses — **0.541 → 0.307 at k = 1** — and F falls with it. The 150% of extra
junctions are alternative splice forms, low-abundance isoforms and alignment noise, and **each one is
another chance for two unrelated loci to match by accident.** More junction evidence is not better
junction evidence.

## Where this leaves the junction definition

- §6u1's annotated-junction arm (F 0.4513 at 6.5% coverage) remains the best junction-based definition,
  and it is far below the shipped rule (**0.7123 at 95.4%**).
- The ceiling is now known to be **two independent walls**: half of family relationships have a member
  with **no introns at all** (biological, confirmed here), and adding more junction evidence to the rest
  **costs precision faster than it buys recall**.
- ⚠ One thing this does *not* say: that read junctions are useless. They confirm **78.6%** of annotated
  junctions, which is a strong validation signal for the annotation, and that is a different use from
  family definition.

## A correction to something I suspected

I expected the containment problem (§6t7–§6u0) and this junction ceiling to be the same phenomenon —
processed pseudogenes. **They are not.** The rejected containment pairs are *more* likely to be
both-spliced (65.8%) than the accepted edges (59.6%), so the containment problem is not the intronless
population. Two separate walls, not one.

## Limits

- One tissue (testis) and one library; a different tissue would move the read census, though not the
  intronless finding.
- 3-read floor fixed in advance; a higher floor would trim the noisy extra junctions, but choosing it
  after seeing these numbers is exactly the fit this pre-registration forbids.
