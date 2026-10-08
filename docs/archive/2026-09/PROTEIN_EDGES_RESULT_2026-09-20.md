# Protein edges vs the §6o8 no-edge gap — result

Run 2026-09-20 against `docs/PREREG_protein_edges_2026-09-20.md` (committed `97a8eac9` before any
alignment). Tool: `bench/protein_edge_gap.py`. Truth: Soto et al. 2025 published families, unchanged
from §6s8. Protein rule copied verbatim from §6ko (longest CDS, blastp e ≤ 1e-5, non-overlapping HSPs
covering ≥ 0.30 of the longer protein) — not re-tuned.

## First, the correction: §6o8's gap is not on this layer

§6o8 reports **> 50% of truth families with no edge on any member** and a pairwise-recall ceiling of
0.052. That measurement is over the **RNA node graph** scored against a **DNA gene-span truth**. Measured
instead over the **DNA gene-body graph** against an **external published truth**, the same statistic is:

| | no-edge families, nucleotide only |
|---|---|
| §6o8 (RNA nodes vs guided gene-span truth) | **> 50%** |
| here (gene bodies vs Soto families), held out chr2/chr8/chr10 | **4/27 = 14.8%** |
| here, chr16 (development) | 2/8 = 25.0% |

⚠ **So this test did not measure the gap as §6o8 states it.** The >50% figure is specific to the RNA
layer, and §6o9 had already diagnosed why: the DNA-level truth *"demands pairs that do not exist as
RNA"* — only 6.7% of guided-truth family pairs align at all as spliced RNA. The dominant gap is a
**layer mismatch**, not missing edges, and it is not fixable by adding edges of any kind.

## What the test did establish: protein edges materially extend edge construction

| arm | families | pairs | nucleotide | **nucleotide ∪ protein** | gain |
|---|---|---|---|---|---|
| **held out** (chr2, chr8, chr10) | 27 | 65 | 52 (80.0%) | **61 (93.8%)** | **+13.8 pts** |
| development (chr16) | 8 | 77 | 71 (92.2%) | 75 (97.4%) | +5.2 pts |

| no-edge families | nucleotide | nucleotide ∪ protein |
|---|---|---|
| **held out pooled** | 4/27 (14.8%) | **2/27 (7.4%)** — halved |
| chr16 | 2/8 (25.0%) | 2/8 (25.0%) — no change |

Pre-registered bar on the held-out no-edge drop (≥ 15 MATERIAL / 5–15 PARTIAL / < 5 NO): **7.4 points →
⚠ PARTIAL.** Worth having; does not by itself close a gap.

## And they are precise — the risk the pre-registration did not cover

Protein edges among Soto-labelled genes, split by whether both endpoints share a family:

| chromosome | within-family | **cross-family (false)** | precision |
|---|---|---|---|
| chr2 | 26 | **0** | **1.000** |
| chr8 | 18 | **0** | **1.000** |
| chr10 | 17 | **0** | **1.000** |
| chr16 (development) | 75 | 12 | 0.862 |
| **pooled** | 136 | 12 | **0.919** |

⭐ **Zero false edges on all three held-out chromosomes.** chr16's 12 cross-family edges are the one
weak cell, and chr16 is NPIP/PKD1P/SMG1 territory where the families genuinely interrelate — the same
boundary the fusion work (§6s9) is about.

⚠ This precision is measured **only over genes Soto labels**. Edges touching unlabelled genes are not
scored, so this is not a genome-wide false-merge rate; §6bt's gene-tight false-merge protocol would be
the way to get one.

## Where the dominant gap actually is, and what would move it

§6o9's own conclusion stands and this run reinforces it: **a non-circular RNA-level truth is the
missing ingredient.** Scoring an RNA-level definition against a DNA-level truth is bounded at ~5%
by construction; scoring it against an alignability-derived truth is circular. §6p0 solved this for
exactly one family (NPIP, against the Dishuck Iso-Seq groups). **n = 1.**

The honest next step is not another edge operator but more RNA-derived truths — the candidates §6o9
already named (Soto restricted to expressed loci, protein-level families, Dishuck-style Iso-Seq groups)
— so that the RNA-level definition can be scored on more than one family without circularity.

⚠ Nothing here is adopted. Adding protein edges to the shipped O1 rule is a separate question that would
need its own pre-registration, a genome-wide false-merge measurement, and a decision about whether an
edge that exists only in protein space belongs in a definition the thesis calls *topological*.
