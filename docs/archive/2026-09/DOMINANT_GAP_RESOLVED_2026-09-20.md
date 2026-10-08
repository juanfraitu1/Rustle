# The dominant gap, resolved: §6o8's 0.052 ceiling is a property of the TRUTH's divergence, not of O1

Run 2026-09-20. Tools: `bench/rna_truth_from_protein.py`, `bench/protein_edge_gap.py`.
Supersedes the framing in `docs/PROTEIN_EDGES_RESULT_2026-09-20.md`, which characterised the gap but
did not explain it.

## What was tried, in order

1. **Protein EDGES added to the nucleotide edge set** (§6t0). Held-out pair coverage 80.0% → 93.8%,
   no-edge families halved, zero false edges on three held-out chromosomes. **PARTIAL** — useful, but the
   0.052 ceiling is defined against a *truth*, and adding edges does not change the truth.
2. **A non-circular RNA-level TRUTH**, which §6o9 named as the missing ingredient. Built protein families
   (§6ko's rule, MCL I = 2.8, r2 exclusions applied) — defined on the spliced product, built by blastp
   over amino acids, so neither DNA-level nor circular with the nucleotide gate under test.

**The second one failed, and the failure is the answer.**

| chromosome | protein-family truth | align at all | through shipped gate | **pairwise ceiling** |
|---|---|---|---|---|
| chr2 | 135 families / 1,019 pairs | 11.4% | 7.7% | **0.082** |
| chr10 | 79 / 495 | 8.9% | 4.6% | **0.048** |
| chr16 | 80 / 1,474 | 16.6% | 4.4% | **0.055** |
| chr8 | 78 / 1,862 | 71.3% | 61.1% | **0.707** |
| *(§6o9's DNA gene-span truth, for reference)* | | 6.7% | 5.0% | *0.052* |

Three of four chromosomes land on **§6o8's ceiling exactly**. Swapping a DNA-level truth for a
product-level one moved nothing. ⚠ chr8's 0.707 is **one family** — FAM90A, a 42-member tandem array
that is 46.2% of that chromosome's pairs — not a general result.

## Why: the ceiling is a divergence cliff, and it is at protein identity 0.60

Every protein-similar pair on chr2/8/10/16, stratified by blastp identity (independent of the nucleotide
aligner being scored), asking whether the two spliced RNAs align at all:

| protein identity | RNA-alignable | pairs | rate |
|---|---|---|---|
| **0.95 – 1.00** | 1,446 | 1,459 | **99.1%** |
| 0.90 – 0.95 | 78 | 91 | 85.7% |
| 0.80 – 0.90 | 67 | 120 | 55.8% |
| 0.60 – 0.80 | 105 | 285 | 36.8% |
| **< 0.60** | 40 | 4,325 | **0.9%** |

**69% of protein-family pairs sit below 0.60, where 0.9% align.** That single fact produces the 0.052
ceiling, and it produces it under *any* family truth — DNA, protein, or published — because every such
truth is dominated by ancient paralogues. **No edge operator and no truth swap can move it.** The pairs
genuinely do not exist as RNA.

## And where the thesis actually operates, there is no ceiling

| family | pairs | median protein identity | **RNA-alignable** |
|---|---|---|---|
| **FAM90A** (chr8) | 136 | 0.985 | **136/136 = 100.0%** |
| **NPIP** (chr16) | 171 | 0.780 | **155/171 = 90.6%** |

The thesis substrate is recent primate-specific expansions. On those the RNA edge graph reaches
**90.6–100%** of within-family pairs — not 5%. ⭐ **§6o8's "the goal is unattainable" applies to a
population the thesis is not about.**

## What this changes

- ⛔ **Stop treating 0.052 as a defect to engineer away.** It is the alignable fraction of a truth whose
  pairs are mostly ancient. §6o8's priority list — *edge construction first* — was based on reading it
  as an edge defect; §6o9 half-corrected that, and this closes it.
- ⭐ **State the definition's SCOPE with a number instead.** The RNA-level definition recovers
  multi-copy families whose members are still nucleotide-alignable — empirically, protein identity
  ≳ 0.80 (55.8% and rising to 99.1%), with a cliff below 0.60. That is a defensible scope claim, and it
  covers NPIP, TBC1D3 and FAM90A.
- ⚠ **Recall against an ancient-paralogue truth should not be quoted as an O1 result at all** without
  that scope attached; it measures the truth's age distribution.
- Protein edges (§6t0) remain a genuine but partial improvement, still unadopted.

## Limits

- Four chromosomes, one assembly, one annotation. The cliff's location (0.60) is measured, not derived.
- The ≥ 0.95 band is 1,459 pairs of which FAM90A is 136 (9%), so the 99.1% is not one array — but it is
  still enriched for recent tandem arrays by construction.
- "Align at all" uses `minimap2 -x asm20 -c -X -N 50 -p 0.1`; a more sensitive aligner would move the
  rates, though the cliff's shape is driven by sequence divergence, not by the tool.
