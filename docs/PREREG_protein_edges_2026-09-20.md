# Pre-registration — PROTEIN EDGES: can they close §6o8's no-edge gap?

**Written 2026-09-20 before any protein alignment is run.** User goal: *"fix dominant gap for families."*

## The gap, as the ledger states it

**§6o8**: over the guided truth, only **27.6% (human) / 34.5% (gorilla)** of truth-family nodes carry ANY
edge; **>50% of truth families have no edge on any member**; pairwise recall is capped at **0.052 / 0.229**
whatever the grouping rule. Its priority list is explicit: **(1) edge construction, (2) node completeness,
(3) grouping — already saturated.**

**§6o9** then showed the RNA ceiling is a truth-level mismatch: only **6.7%** of guided-truth family pairs
align at all as spliced RNA and **5.0%** pass the shipped gate, against a 5.2% ceiling — *"edge construction
is already at its limit"* **for nucleotide alignment**.

## The hypothesis

The limit in §6o9 is a limit of **nucleotide** alignment. Protein-level similarity survives synonymous
divergence that nucleotide identity loses, so a protein edge may exist where a nucleotide edge cannot.

**H:** Adding a protein-level edge materially reduces the fraction of truth-family pairs with no edge.

## What is already refuted, and why this is not that

- **r826** — protein rescue of *nucleotide-unreachable old retrogenes* — **not supported, 3/7 rescued**,
  with three distinct causes (metric, saturation, real protein divergence). ⚠ That was **7 retrogenes from
  Addendum Z**, a specific hard set, and it asked whether protein rescues *those*. It did **not** measure
  the no-edge rate over a general truth-family population, which is what §6o8's ceiling is about. This test
  is therefore a different question on a different population, and r826's 3/7 is the prior to beat.
- **r1072** — cDNA/protein homology **replacing de novo** construction. Not proposed. Here protein edges
  are an ADDITION to the existing edge set, and the de novo path is untouched.
- **§6ko** already builds protein families reproducibly (`bench/protein_families.py`; fresh hold-out
  chr1/2/3 at 0.935 / 0.999 / **F 0.977**), so the protein machinery itself is not what is under test.

## Method — frozen

- **Substrate:** chr2, chr8, chr10 (the zero/near-zero-exposure held-out set from §6s8) plus chr16 (the
  development chromosome) as the comparator.
- **Truth:** Soto et al. 2025 published families (S1C `Family ID`), ≥ 2 members on the chromosome — the
  same external truth used in §6s8, unchanged.
- **Nucleotide edge (the incumbent):** the shipped gene-body gate already computed for §6s8 —
  `minimap2 -x asm20 -c --eqx -P`, identity ≥ 0.7, cov_longer ≥ 0.3, ≥ 300 bp.
- **Protein edge (the addition):** one protein per gene = longest CDS, translated; all-vs-all
  `blastp -evalue 1e-5`; edge iff non-overlapping HSPs cover **≥ 0.30 of the longer protein** — §6ko's
  own rule, copied, not re-tuned.
- **Measure, per chromosome and pooled:** the fraction of within-truth-family PAIRS with (a) a nucleotide
  edge, (b) a protein edge, (c) either. And the §6o8 statistic: fraction of truth families with **no edge
  on any member**, under nucleotide alone vs nucleotide ∪ protein.

## The bar — committed now

| pooled no-edge FAMILY fraction, nucleotide ∪ protein vs nucleotide alone | verdict |
|---|---|
| absolute drop ≥ 15 points | ⭐ **MATERIAL** — protein edges are a real fix for the dominant gap |
| 5–15 points | ⚠ **PARTIAL** — worth having, does not close the gap |
| < 5 points | ⛔ **NO** — report it and stop; do not tune the coverage floor to manufacture a win |

⚠ **Declared before the run:** a protein edge between two genes is evidence they are paralogues, **not**
evidence their spliced RNAs are mutually alignable. If this passes, the honest claim is *"the family
relation is recoverable at the protein level where nucleotide alignment fails"* — **not** that the RNA
edge graph improved. Any O1 rule change would be a separate, separately pre-registered question.
⚠ I will not drop a chromosome, re-tune the 0.30 coverage floor, or switch truth after seeing the numbers.
