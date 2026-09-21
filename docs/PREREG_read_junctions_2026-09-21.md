# Pre-registration — READ-DERIVED junctions for the junction family rule

**Written 2026-09-21 before any read junction is extracted.** Follows §6u1, which named this as the
follow-up: *"I used annotated junctions. Read-derived junctions would change the census entirely."*

## The specific hypothesis

§6u1 built a drift-tolerant junction family rule and found its ceiling is not drift but **missing
introns**: of 3,785 shipped homologous edges, only **48.6% have both genes spliced at all** by the
annotation, and 19.5% of genes carry zero annotated junctions. Pseudogene fragments — the §6t7–§6u0
population — are annotated as single-exon.

> **H:** many of those genes are **transcribed and spliced**, and the annotation simply does not record
> it. Read-derived junctions would then raise the both-spliced fraction, and with it the junction rule's
> coverage and F.

If instead the intronless genes are genuinely intronless — processed pseudogenes, which are
retrotransposed from mRNA and *have* no introns by construction — reads will add nothing, and §6u1's
ceiling is biological rather than annotational. **Both outcomes are informative and the prediction is
made before looking.**

## Method — frozen

- Substrate: **A119b.t2t.bam** (human testis, the lab's deep library), held-out chr2/chr8/chr10, which
  carry 5.12M / 2.21M / 2.41M mapped reads.
- Junctions from the CIGAR **`N`** operator (`N` in an RNA CIGAR = intron spliced out), reads filtered
  **`-F 2308`** (primary, mapped, non-supplementary) per the standing invariant.
- A junction is kept at **≥ 3 supporting reads** — the same floor `RUSTLE_GATE_MIN_READS` uses.
- A junction is assigned to a gene if both endpoints fall inside that gene's span.
- Everything downstream is §6u1's rule unchanged: project through the PAF CIGAR, match at ±10 bp,
  edge iff ≥ k shared, k swept 1–5, connected components, same scorer, same Soto truth.

## What is reported

1. **The census first**: genes with ≥1 read junction vs ≥1 annotated junction, and the both-spliced
   fraction of shipped homologous edges under each — this is the quantity the hypothesis is about.
2. Then F / sensitivity / precision / node coverage, exactly as §6u1 reported them.

## The bar — committed now

Comparators: §6u1's annotated-junction arm **F 0.4513 at 6.5% coverage**; the shipped rule **F 0.7123 at
95.4%**.

| outcome | verdict |
|---|---|
| F ≥ 0.7123 **and** coverage ≥ 80% | ⭐⭐ **A REAL ALTERNATIVE DEFINITION** |
| F ≥ 0.5513 (annotated + 0.10) **and** coverage ≥ 13% (2× annotated) | ⭐ **READS RESCUE IT MATERIALLY** |
| within ±0.10 of 0.4513 | ⚠ **NEUTRAL** — the ceiling is biological, not annotational |
| F < 0.3513 | ⛔ **WORSE** |

⚠ Declared now: if the census shows the intronless genes stay intronless under reads, that is the
**answer to the question**, not a failed run — it would establish that §6u1's ceiling is a property of
processed pseudogenes and cannot be lifted by better annotation. I will report the census even if the
downstream F is unchanged, and I will not lower the 3-read floor after seeing the numbers.
