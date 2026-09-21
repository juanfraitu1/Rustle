# Pre-registration — a family definition based on SHARED SPLICE JUNCTIONS

**Written 2026-09-20 before any junction edge is built.** User: *"find an option for defining a spliced
family based on the splice junctions."*

## Why now, and what the two prior attempts actually refuted

§6u0 closed the containment question with: *"what would still be new is not another pair statistic but
intron-chain identity — do the two genes share the same junctions, not merely exonic bases."* Two prior
attempts exist and **neither tested a drift-tolerant rule**:

- **r344** — *"splice-junction / intron architecture can DEFINE families with precision"*: ⛔ **strict
  concordance** breaks the 389-member bridge 389→32 but leaves APOBEC3/RFPL with **0 members**, because
  *"paralog intron lengths and counts DRIFT after duplication; no operating point does both."* The
  refutation is specifically of **strict** concordance — matching whole chains, lengths included.
- **r85** — a **junction-crossing predicate** in the definition: rejected 12.80% of shipped edges with a
  **100× monotone bias by exon count** (0.3555 at 2 exons → 0.0036 at > 6). That is a different
  predicate (does an edge cross a junction), not a chain comparison.

## The rule — drift-tolerant by construction

A junction is a donor/acceptor pair. Copies sit at different loci, so junctions cannot be compared by
coordinate; they are compared **through the pairwise alignment**:

```
for a pair (A,B) with a PAF record (CIGAR available, -c was used):
    project each junction of A into B's frame through the alignment
    a junction MATCHES if BOTH its endpoints land within ±10 bp of a junction of B
    shared(A,B) = number of matched junctions
EDGE iff shared(A,B) >= k
```

This addresses r344's exact failure mode. **Intron length drift is tolerated** because junctions are
matched by projected position, not by intron length — an intron that grew still has its flanking
junctions projecting correctly if the exons align. **Intron count drift is tolerated** because the rule
requires *k shared*, never *all shared*. The ±10 bp tolerance is the only parameter and is fixed now,
at the same order as the `junction_canonical_tolerance` the assembler already uses.

Swept over **k ∈ {1, 2, 3, 4, 5}**. Grouping: **connected components** and **MCL I=2.8**, both reported,
since §6t9 measured that the operator matters (F spread 0.41–0.73).

## The coverage limit, stated before the run

On the held-out chromosomes, **1,368 of 7,020 genes (19.5%) have no junction at all**. A junction rule
**structurally cannot place any of them**. §6n5's lesson — every triangle operator scored well on F while
dissolving all 2-member groups — is the precedent: a rule that wins on F while covering a fraction of the
nodes has not defined families. **Node coverage is therefore reported beside F, never instead of it.**

## Substrate, truth, comparators

Held-out chr2/chr8/chr10, Soto S1C families (≥ 3 members), `bench/heldout_family_score.py` — identical to
every other arm today. Comparators on this exact substrate:

| arm | F | coverage |
|---|---|---|
| MCL on the shipped graph (shipped rule, via the port) | **0.7123** | 95.4% |
| label propagation, shipped graph | 0.7295 | 100% |
| connected components, shipped graph | 0.5811 | 100% |

## The bar — committed now

| outcome | verdict |
|---|---|
| F ≥ 0.7123 **and** coverage ≥ 80% of the spliced genes | ⭐⭐ **A REAL ALTERNATIVE DEFINITION** |
| F ≥ 0.7123 but coverage < 80% | ⚠ **BOUGHT WITH COVERAGE** — §6n5's trap; report, do not adopt |
| 0.60 ≤ F < 0.7123 | ⚠ **VIABLE BUT WORSE** — a junction definition exists and costs F |
| F < 0.60 | ⛔ **NO** — junctions do not define families at this resolution |

Secondary, committed now: the rule's behaviour **on the 19 TRUE / 57 FALSE containment pairs** of
§6t7–§6u0, reported as precision/recall, since that is the population every other channel failed on.

⚠ I will not change the ±10 bp tolerance, the k grid, or the truth after seeing scores, and I will not
report F without the coverage column.
