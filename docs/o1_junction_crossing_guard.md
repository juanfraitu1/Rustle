# A threshold-free edge predicate that works at n = 2

**Status 2026-08-19. Offline (T8), nothing through the shipped binary, no default moved.**
Companion to [`o1_genome_anchored_repeat_gate.md`](o1_genome_anchored_repeat_gate.md) and
[`o1_error_case_census.md`](o1_error_case_census.md).

## 1. The problem this addresses

At **n = 2** — **348 of 494 = 70.45%** of the gorilla catalog — the γ-quasi-clique machinery is inert
and the entire definition reduces to **one coverage number**. That number is provably non-separating:
the accepted true pair GFPT1×GFPT2 scores **0.5353** while the rejected false pair ATP1A1×ATP4A scores
**0.5689**. No threshold on it orders the classes correctly.

The census prescribed the shape of any fix: **change the denominator or the substrate, not the
threshold.** This is the first candidate that does.

## 2. The predicate

> **The homology must cross a splice junction.** Reject an edge iff the passing alignment lies
> entirely within a **single exon on both sides**.

`max_exon_frac(side) = max over exons of (alignment ∩ exon) / alignment length`; reject iff it is
**1.0 on both sides**.

**Why it escapes the named hole: it has no length denominator at all.** The scale-free defect exists
because a ~1 kb repeat is ≥ 0.50 of any node under 2 kb. A repeat confined to one exon spans one exon
**at every node size**. The statistic is structurally immune to the defect rather than tuned against it.

**It is threshold-free** — "crosses a junction" is discrete. No fitted number enters.
**It is pair-local**, depending only on the two nodes and their exon structures ⟹ **P1 is untouched**.
**It is an edge predicate**, so it works identically at n = 2, where no graph structure exists.

## 3. Result on the frozen arms (unit = pair, GGO)

| | FP rejected / 14 | TP lost / 150 |
|---|---:|---:|
| coverage-of-longer @ 0.20 | 1 | damages many |
| transcript-orientation guard | 6 | 4 / 9,032 edges |
| genome-anchored repeat veto @ M=50 | 10 (of 12 scored) | **0** |
| **junction-crossing** | **12** | **9 = 0.0600** |
| **union of the last two** | **13** | 9 = 0.0600 |

The single FP the union misses is **FP00058**, the LAGE3 processed pseudogene against its own parent —
a **truth-label failure, not a false merge**, so it *should* survive. Effectively **13/13 real FPs.**

The two guards are complementary, not redundant: junction-only 3, gmult-only 1, both 9.
The genome-anchored veto abstains on 2 FP and 15 TP pairs (no shared 21-mer at identity 0.69–0.80).

## 4. The cost is real — 9 true pairs, and they are not an artifact

⚠ **First reported here as "all 9 are single-exon genes, so the cost is zero". That was WRONG**, and the
error is worth recording: the TP arm's `a_nex` column is **exons TOUCHED**, not total exons (that is
`a_tot_ex`), and touched-count is 1 **by construction** whenever `max_exon_frac = 1.0`. The check was
circular. With the correct column there are **0 single-exon nodes in either arm** — consistent with the
shipped catalog having **0/1415** single-exon nodes.

The 9 lost pairs are genuine: multi-exon nodes (2–7 exons) at identity 0.7117–0.7617 and coverage
0.5132–0.6586 whose homology genuinely sits inside one exon on both sides. That is a real biological
class — paralogues sharing a single domain-bearing exon — and rejecting them is a **recall cost, not
noise**. No scoping removes it.

**T13 verified:** both arms compute `max_exon_frac` by the identical formula
(`max(exon overlap) / alignment length`), from two different generators.

## 5. What is NOT established

1. **T8.** Offline. Nothing through the shipped binary; `E_r` is unchanged.
2. ⚠ **The FP side is not held out.** This discriminator was *found* on these same 14 FPs (the
   "≥99% of the alignment inside one exon on both sides, FP 13/14 vs TP 9/150" characterisation), so
   **12/14 is a description of a known set, not a rate**. **The load-bearing number is the TP side,
   9/150**, because the TP arm was never selected for exon structure.
3. ⚠ **Partly entailed by the truth predicate.** True families are annotated multi-exon protein-coding
   genes, so "shares more than one exon" is somewhat implied by being a real family. This does not
   dissolve the result — the rule rejects 12 pairs that the coverage clause *accepted* — but it means
   the FP rate needs an independent FP set.
4. **Exon structure is read-derived**, so "exon" is a model feature. Not more circular than the
   incumbent (the reps are exon-sums), but not independent evidence either.
5. **Gate vs flag undecided**, and the choice is a real one:
   - **gmult as a gate, junction as a flag** — 10/14 rejected at **zero** TP cost, 3 more flagged.
     The conservative ship.
   - **union as a gate** — 13/14 at 6.00% TP cost. A deliberate recall-for-precision trade.
   - **both as flags** — 13/14 surfaced, flag precision 13/22 ≈ 0.59.

## 6. What this changes about O1's defensibility

At n = 2 the definition was one non-separating number. It can become **one coverage number plus two
pair-local structural predicates**, one of them threshold-free and immune to the scale-free defect by
construction. That does not repair the coverage clause — the named hole and its 8.30% exposure ceiling
stand — but it means the *edge* carries structural evidence rather than a single scalar, and it does so
at exactly the family size where the graph machinery offers nothing.

## 7. Reproduction

Columns are already in the frozen arms: `fp14_detail.tsv` (`a_maxexonfrac`/`b_maxexonfrac`,
`a_nexon` = total exons) and `tp150_detail.tsv` (`a_mx`/`b_mx`, `a_tot_ex` = total exons,
⚠ `a_nex` = exons **touched**). Generators: `o1_antifp/analyze.py` and `o1_antifp/tpstats.py`.
