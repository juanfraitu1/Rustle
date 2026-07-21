# Identifiability-Respecting Locus Merge — Design

**Date:** 2026-07-21
**Status:** approved design, pre-implementation
**Motivation.** The advisor's standing objection: *the method must not merge things that should
be kept separate.* The Soto floor decomposition (`bench/soto/soto_floor_decomposition.tsv`) shows
**28 "distinguishable-but-MERGED" members** — copies with unique-mapper reads (`unique_own_reads`
of 40, 24, 20, 11…) that were nonetheless absorbed into a co-located locus. This is a genuine
over-merge, and its cause is precise.

**Root cause (verified in code).** `distinct_locus_reps` (`denovo_pipeline.rs:2881`) collapses two
co-located same-strand copies **unconditionally** (`collapse = true`), consulting only coordinates
and strand — it never looks at reads. So a copy with 40 unique mappers merges into its overlapping
neighbour purely because their spans overlap.

**The fix, in one line.** Merge two co-located copies **iff no read distinguishes them** — i.e.
merge iff χ(H) says they are one colour. This applies the thesis's own identifiability theorem
(MCC = χ(H)) to the one stage that currently ignores it, and it is coordinate-threshold-free.

---

## 1. Scope

**In scope:** the **28 "distinguishable-but-MERGED"** cases — co-located copies wrongly collapsed
despite read evidence separating them.

**Out of scope:**
- The **10 "distinguishable-but-unseeded"** cases (no locus was created at all — a *detection*
  gap, not a merge). Different fix, deferred.
- The **36 true-K=0** members (genuinely indistinguishable from RNA — they *must* still merge;
  this design must not touch them).
- The 5-ε unification (deliverable C) and the assembler carve (separate efforts).

**Supersedes B2.3.** Deliverable B deferred the `distinct_locus_reps` (any-overlap) vs
`distinct_loci` (reciprocal-50%) choice to C. This design resolves it: the answer is neither
overlap threshold — it is read-separability. `consolidation_divergences.md` will be updated to
record that B2.3 is resolved here.

---

## 2. The mechanism

### 2.1 The merge criterion

Replace the unconditional `collapse = true` for same-strand co-located copies with:

> **collapse(i, j) = co-located AND same-strand AND NOT reads_distinguish(i, j)**

where `reads_distinguish(i, j)` is true when **any** of the following holds, each at the
pipeline's *existing* read-support floors (no new constant is introduced):
- **unique-mapper support:** copy i has ≥ `min_reads` (the conflict floor, currently 3) reads that
  map uniquely (MAPQ > 0) to i and not j, or vice-versa — read from the existing
  `family_mapq0_support` signal; **or**
- **read-supported PSV:** there is a column where ≥ 2 reads (the existing `PSV_MIN_ALLELE_READS`)
  carry an allele one copy has and the other does not; **or**
- **copy-specific junction:** a junction pinned to one copy and not the other (the existing
  `copy_junction_support` signal).

If none hold, the two copies are indistinguishable from the reads (true K=0) and collapse — the
honest, correct merge. The antisense-minority artifact rule (`ANTISENSE_MINORITY_DENOM`) for
opposite-strand pairs is **preserved unchanged**.

### 2.2 Why this is χ(H), not a new heuristic

`reads_distinguish` is exactly the edge predicate of the read-conflict graph restricted to a
co-located pair: two copies that no read separates fall in the same conflict colour (merge);
copies some read separates are different colours (keep separate). The family copy count becomes
`max(coordinate-distinct loci, χ(H))` over the block — the theorem the thesis already argues,
now enforced at merge time. There is **no arbitrary overlap cutoff**.

---

## 3. Components & data flow

| Unit | Responsibility | Input | Output |
|---|---|---|---|
| `reads_distinguish(i, j, ev)` | the pairwise separability oracle | two copies + read evidence | bool |
| `distinct_locus_reps(copies, ev)` | merge loop, now evidence-aware | copies + evidence | copy set |
| evidence threading | carry `placements` / mapq0-support / PSV/junction to the merge | already computed at `denovo_pipeline.rs:1426–1441` | `ev` handed to the two call sites |

**Data flow:** reads → `build_read_placements` + `conflict_edges` + `family_mapq0_support`
(all existing, `:1426–1441`) → packaged as the evidence `ev` → `distinct_locus_reps` consults
`reads_distinguish` per co-located pair → copy set.

**Touch points:** `distinct_locus_reps` signature (adds `ev`); its two call sites
(`denovo_pipeline.rs:2265`, `:2522`); `refine_families_exon_sum`'s signature if `:2522`'s caller
must forward `ev`. The evidence source functions are **unchanged** — this only *consumes* them.

**Feasibility note.** `family_mapq0_support` already returns per-copy unique-vs-MAPQ-0 support;
PSV/junction evidence is on the assembled family. Where the merge runs *before* per-family PSV
assignment, the unique-mapper leg (from `placements`/`family_mapq0_support`) is the primary
signal and is sufficient for the 28 cases (they are flagged by `unique_own_reads`); the PSV/
junction legs are additive where that evidence is in scope. The plan resolves the exact evidence
available at each of the two call sites and wires the legs that are reachable there.

---

## 4. Validation (result-changing on purpose — NOT byte-identical)

This deliberately changes the catalog. It is gated by a **Soto before/after**, not byte-identity:

1. **Recovery up:** Soto member sensitivity rises from **276/362 (76.2%)** — the merged copies
   with genuine unique support become their own recovered members. Report the new number.
2. **Precision held:** recovery precision stays at **100%** (no false-positive copies); no
   single-copy locus false-splits. This is the guardrail against over-correcting into over-splitting.
3. **K=0 untouched:** the **36 true-K=0** members still merge (they have no distinguishing
   reads) — spot-check named K=0 families before/after; their graphs must be unchanged.
4. **Floor decomposition re-run:** `python3 bench/soto/soto_floor_decompose.py` — the
   "distinguishable-but-MERGED" bucket shrinks; the shrinkage lands in *recovered*, not in a new
   false-split category.
5. **`cargo test` green**, plus a new unit test for `reads_distinguish` (a pair with unique reads
   → kept separate; a pair with no distinguishing reads → merged).

The B byte-identity gate (`byte_identity_gate.sh`) will now FAIL by design on the affected
families — that is the intended improvement. Re-baseline the gate to the improved output **only
after** the Soto before/after confirms recovery↑ / precision-held, and record the delta.

---

## 5. Risks

- **Over-splitting** (the opposite error the fix could introduce): a single copy with a few noisy
  unique-ish reads gets split. Mitigation: the unique-mapper leg uses the existing `min_reads`
  floor (3), not 1; the Soto precision check (must stay 100%) is the backstop.
- **Evidence not in scope at a call site.** If `:2522` runs before PSV assignment, only the
  unique-mapper leg is available there; the plan wires what's reachable and documents any leg
  deferred. The unique-mapper leg alone covers the 28 flagged cases.
- **Antisense rule regression.** The opposite-strand minority-artifact collapse must be preserved
  exactly; the change is scoped to the same-strand branch.

---

## 6. Success criteria

- `distinct_locus_reps` merges a co-located same-strand pair only when `reads_distinguish` is
  false; unit-tested both ways.
- Soto member sensitivity **> 76.2%**, precision **= 100%**, the 36 K=0 members still merged.
- `soto_floor_decompose.py` shows the merged bucket shrunk into recovered members.
- `cargo test` green; the change is scoped to the merge decision + evidence threading (no
  unrelated refactor).
- `consolidation_divergences.md` updated: B2.3 resolved by read-separability.
