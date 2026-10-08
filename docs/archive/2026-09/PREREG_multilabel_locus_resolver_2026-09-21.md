# PREREG — multi-label locus-to-gene resolution for scoring

**Written 2026-09-21 before running any comparison.** Register 962/963 found that a max-overlap,
one-name-per-locus resolver used throughout this session's own scoring (`score_family_def.py`,
`cover_and_jn_definition.py`) undercounts recall: when a de novo locus spans multiple genes (readthrough),
only the single largest-overlap gene gets credit, so a real, correctly-assembled gene riding in the same
locus as a bigger neighbour scores as a false negative even though its sequence is very likely present.
This fixes the RESOLVER, not the family definition or the assembler — it is an evaluation-methodology
change, scoped to how a locus's identity is read for scoring.

## The fix

Replace winner-take-all with **multi-label, coverage-of-the-GENE**:

    for a locus with exonic union E_locus, and a candidate gene g with its OWN exonic union E_g
    (unioned across every annotated transcript of g):
        credit(locus, g)  iff  |E_locus ∩ E_g| / |E_g| >= FLOOR

A locus can now credit **every** gene it substantially covers, not just the biggest. The floor is on
**EXONIC overlap divided by the gene's own exonic length** — not genomic span, which is exactly what let
a 300 kb readthrough locus outvote a 3 kb gene it barely touches. `FLOOR = 0.50`, chosen before looking,
matching the existing project convention (§6u7/§6v1's per-copy scorer already uses 50% coverage as its
correctness floor, applied in the complementary direction).

⚠ This does not touch the TRUTH side (Soto/protein-referee labels are unchanged) — only the PREDICTION
side's locus-to-name mapping changes, so it cannot condition the denominator on the prediction
(register 770's trap).

## What could go wrong, stated before looking

Multi-labeling a locus is a PRECISION risk, not just a recall gain: if a locus's multi-label set spans
two genes that are NOT the same truth family, the cluster containing that locus now looks like it covers
BOTH families, which can inflate the bipartite match's apparent sensitivity to one family at the cost of
appearing to cover a family it does not truly represent. This must be measured, not assumed away.

## Bars — stated before looking

1. **Node coverage (real chr16, `RUSTLE_JUNCTION_MAJORITY` default-on)**: the 5 identified naming-collision
   genes (`NPIPB4`, `RRN3P1`, `RRN3P2`, `SLC7A5P2`, `SLX1A-SULT1A3`) must now show as covered. If any of
   the 5 is still uncovered, the floor or the exon-union construction has a bug, not a biological limit.
2. **Precision must not collapse.** Report sensitivity, precision and F together, never sensitivity alone.
   A drop of more than 0.03 in pooled precision (matching §6u8's guard) without a compensating F gain
   means the multi-label credit is loose enough to be admitting spurious matches, and the floor should be
   raised, not the result accepted as-is.
3. **CDR2 is a stated NEGATIVE CONTROL, not a target.** §6v8 already showed CDR2's problem is that the
   homology signal belongs to a different gene's exons entirely (`RRN3P3`, not in Soto's truth) — no
   locus-naming fix can create sequence overlap with CDR2's own exons where none exists. CDR2 remaining
   uncovered here is evidence the fix is not overreaching, not evidence the fix failed.
4. Report pooled family F on BOTH truths (Soto cover, protein referee), old resolver vs new, side by side.

---

## OUTCOME (appended 2026-09-21, after scoring)

**REFUTED at every floor tested (0.5, 0.7, 0.9), on both truths, on both arms.** Bar 1 (the 5
naming-collision genes must show as covered) partially passed after two implementation bugs were found
and fixed along the way (see below) — 2 of 5 (`NPIPB4`, `RRN3P1`) were recovered; `RRN3P2` was found to be
a DIFFERENT failure mode (its locus is real and well-matched but has ZERO surviving homology edges at
all, unrelated to naming); `SLC7A5P2` exposed a genuine two-pass-parsing bug (fixed); `SLX1A-SULT1A3` not
re-examined individually. Bar 2 (precision must not collapse) FAILED decisively:

| arm | truth | old prec | floor=0.5 | floor=0.7 | floor=0.9 |
|---|---|---|---|---|---|
| GUIDED | Soto | 0.4294 | 0.3281 (-0.101) | 0.3449 (-0.085) | 0.3449 (-0.085) |
| GUIDED | referee | 0.3203 | 0.2856 (-0.035) | 0.2738 (-0.047) | 0.2738 (-0.047) |
| REAL de novo | Soto | 0.3313 | 0.2399 (-0.091) | 0.2518 (-0.080) | 0.2981 (-0.033) |
| REAL de novo | referee | 0.2183 | 0.2056 (-0.013) | 0.1873 (-0.031) | 0.1676 (-0.051) |

Pooled F is WORSE than the existing max-overlap resolver in **all 4 comparisons at every floor**:
0.3761/0.2041/0.2186/0.1391 (best floor per cell) vs 0.4174/0.2248/0.3184/0.1817 old. Raising the floor
to 0.9 (near-total gene coverage required before crediting) does not rescue precision back to baseline in
any cell.

**Why it fails, and what that means.** Multi-labeling doesn't just add correct credit — it inherits the
locus's own ambiguity. A readthrough-fused locus genuinely contains sequence from multiple, often
UNRELATED genes (a real paralog riding alongside an uncharacterized neighbour with no family
relationship). Multi-labeling credits BOTH to the same cluster, which the bipartite/pairwise scoring then
reads as evidence the cluster represents both truth families at once — a false pairing that costs more
in precision than the naming-collision recall it recovers. This holds even at floor=0.9, meaning the
problem isn't a badly-chosen threshold; it's that per-locus labeling (any scheme) cannot correctly score
an object that is genuinely a mixture. **The fix is not at the scoring layer.** The underlying node still
needs to be split at construction time before it can be scored correctly — the standing, still-open
readthrough-aware node-split target from §6v1/§6v4, not a resolver change.

**Two real implementation bugs found and fixed along the way, unrelated to the refutation's substance:**
1. `gene_exon_unions` used a single-pass parser assuming parent-before-child GFF order. RefSeq's own file
   violates this for at least one record (`SLC7A5P2`'s sole exon sits one line above its own gene and
   transcript records, all three sharing identical coordinates) — fixed with a two-pass parse.
2. The locus-reconstruction script initially copied `heldout_family_score.py`'s `start + 1` coordinate
   convention unverified. `mcl_families` (`src/bin/mcl_families.rs:893-901`) writes cluster-member
   coordinates VERBATIM, matching this session's own 1-based GTF-derived locus GFF3 exactly — the +1
   shifted every member 1 bp past its own locus's true start, silently excluding a locus's OWN raw record
   from its self-containment search and collapsing every score to near-zero. Fixed by removing the
   erroneous offset; verified directly against a raw GTF record (`DN_chr16_32747532_2`, span
   32747533-32751315, matching its clusters.tsv row byte for byte).

**Conclusion: register 962/963's naming-collision diagnosis stands — it is real and independently
verified — but the natural scoring-side fix for it is refuted.** The node itself needs to change, not
how it's read.
