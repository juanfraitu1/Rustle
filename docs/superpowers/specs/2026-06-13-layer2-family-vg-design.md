# Layer-2 Family Variation Graph — Design

**Date:** 2026-06-13
**Status:** Design (awaiting review → implementation plan). No code yet.

## Goal

Re-architect Rustle's VG (paralog-family) mode into two strictly separated layers so that
VG mode **always ⊇ baseline** by construction, and copy/isoform recovery becomes a *purely
additive* layer on top of an untouched baseline:

- **Layer 1 — baseline, untouched.** Every locus gets exactly the baseline per-locus splice graph
  and transcripts it would get with VG off. Coordinates and bundle membership are sacred. (Already
  shipped — commit `664919c`, Phase 1.)
- **Layer 2 — family variation graph over a *subset*.** Take Layer-1's per-locus splice graphs,
  identify the subset that are genuinely multi-copy paralog families, merge each family's homologous
  graphs into one variation graph (shared exons/junctions → shared nodes; copy-specific bases/junctions
  → private nodes), and use secondary/supplementary alignments — fed from a **side-index**, never from
  `bundle.reads` — as extra evidence to (a) recover starved/reference-collapsed copies and (b) surface
  extra isoforms. Layer-2 transcripts are unioned with Layer-1 by intron chain: only genuinely novel
  copies/isoforms emit; a baseline transcript is never dropped or altered.

This is the thesis object (`project_thesis_framing_family_vg`): one variation graph per family,
copies = paths, ambiguous secondaries = shared evidence, recovery = constrained flow-decomposition on
the family graph. Layer 2 is where that machinery lives; Layer 1 guarantees we can never do harm.

---

## Background / why this re-architecture

### The bug Phase 1 fixed

Un-guided VG used to ingest secondary/supplementary reads into `bundle.reads`
(`detect_bundles_from_bam_with_snp`, `bundle.rs:~1422`, the old `vg_include_secondary =
config.vg_mode || ...`). Those cross-mapped secondaries inflated per-gene-graph transfrags
(measured: gene 645 had 17 transfrags vs baseline's 10), and flow-extraction then yielded **0**
transcripts for a whole ~200 kb / 5-gene region on chr19 (`NC_073243.2:111.86–112.07 Mb`). VG
*dropped* 17 baseline transcripts — a direct VG ⊄ baseline violation
(`project_vg_drops_baseline_region_rootcause`).

**Phase 1 (commit `664919c`) removed `config.vg_mode ||` from that one ingest gate.** VG no longer
puts secondaries into `bundle.reads` by default: VG-default coord-signature == baseline (chr19
2013 == 2013 tx, dropped region 0 → 17), suite 302/0, `RUSTLE_PRECISE=1` byte-identical to
`4705ab1`. Copy recovery (e.g. RABL2 1 → 7 isoforms) is **paused** — it relied on
secondaries-in-bundle, which is exactly what we removed.

### The caveat Phase 1 created (and that Layer 2 must answer)

`build_multimap_index` (`vg.rs:70-93`) reads `bundle.reads`. With secondaries no longer in
`bundle.reads`, family-discovery signal collapses (measured: 2125 → 313 multimap reads). So Layer 2
**cannot** read secondaries from bundles — it must read them from a separate side-index. This is the
central plumbing decision of this design.

### The user's invariants (verbatim intent)

- *"The idea of the vg mode is to enhance baseline and only make a family of things that have
  evidence of being them; it should not drop bundles from baseline."*
- *"Original bundles should not be touched; secondary alignments should be used for amending graphs
  instead."*
- *"It should not touch bundles at all… bundles that are dropped if all of the reads that compose them
  are secondary, those could be extra in vg mode if we have enough proof that they are real."*

Layer 2 obeys all three by being additive-by-construction and never re-entering the bundling path.

---

## Architecture

```
                 BAM (primary + secondary/supplementary)
                          │
          ┌───────────────┴────────────────┐
          │                                │
   primary records only           secondary/supplementary
          │                                │
          ▼                                ▼
  ┌─────────────────┐            ┌─────────────────────┐
  │ LAYER 1         │            │ secondary side-index │   (NOT bundle.reads)
  │ baseline bundling│           │  - per-locus overlap │
  │ + splice graphs  │           │  - cross-map links   │
  │ + transcripts    │           │  - all-secondary regs│
  └────────┬─────────┘           └──────────┬──────────┘
           │  Layer-1 graphs G_i            │
           │  (clean, baseline coords)      │
           ▼                                ▼
        ┌──────────────────────────────────────────┐
        │ LAYER 2                                   │
        │ 1. family discovery (subset of G_i)        │
        │ 2. merge family graphs → variation graph   │
        │ 3. amend with side-index secondaries       │
        │ 4. constrained flow-decomp → copy/iso paths│
        │ 5. union-by-chain vs Layer 1 (additive)    │
        └────────────────────┬─────────────────────┘
                             ▼
              Layer-1 transcripts  ∪  novel Layer-2 transcripts
```

The key structural commitment: **Layer 2 reads Layer-1's already-built per-locus splice graphs
(clean, baseline coordinates), not bundles and not raw reads.** It never re-bundles, never re-splits,
never shifts a coordinate. The only new reads it sees are secondaries, and only through the
side-index.

---

## Components

### C1. Secondary side-index (the plumbing keystone)

A structure built once per chromosome from the secondary/supplementary records that Layer 1 drops.
It is the *only* place secondaries live; `bundle.reads` stays primary-only (baseline).

**What it stores, per secondary/supplementary alignment:**
- read id / qname (to link a secondary back to the primary it shadows)
- mapped span + intron chain (junctions) on this placement
- the Layer-1 locus (graph index) it overlaps, if any
- alignment score / NM (for PSV / decisive-evidence later)

**What it must support (the two consumers):**
1. **Family discovery** — replace the lost `build_multimap_index` signal. A cross-mapped read = a
   read whose primary lands in locus A and whose secondary lands in locus B → an A–B family link.
   This restores the 2125-read signal *without* touching bundles.
2. **Graph amendment** — for a given Layer-1 locus, return the secondaries overlapping it and their
   candidate junctions/sub-paths.

**Where it is built:** at the point Layer 1 currently *discards* long sec/supp (the same gate Phase 1
narrowed, `bundle.rs:~1422`). Instead of only dropping them, route them into the side-index. Phase 1
already proved dropping them from bundles is correct; this just stops *throwing them away*.

**Design choice — global vs per-bundle.** Cross-map links are inherently cross-locus (siblings can be
>16 Mb apart, as in the chr19 region), so the index must be at least chromosome-global for discovery.
Per-locus *views* are derived from it for amendment.

> ⚠ Open decision for review: a chromosome-global side-index over all secondaries can be large on
> repeat-rich chromosomes (chrY). Mitigation options: (a) keep only secondaries whose primary or
> placement falls in/near a candidate family locus; (b) cap per-locus secondary count with a logged
> drop (no silent truncation — `project` memory rule). Recommend (a): prune the index to family-candidate
> loci after discovery's first pass.

### C2. Family discovery (subset selection)

Choose which Layer-1 loci form genuine multi-copy paralog families. Two signals, **both** required
(advisor Canzar is "really suspicious of methods with no similarity threshold for merging graphs" —
`reference_advisor_canzar`):

1. **Sequence/exon similarity** — exon-set k-mer similarity between the two loci's Layer-1 graphs
   exceeds a threshold. Reuses the exon-restricted discriminative-minimizer machinery
   (`project_minimizer_exon_restriction`): compute on **exons only** (spliced reads carry no introns),
   canonical minimizer hashing. Gated by a CLI threshold (working name `--family-exon-similarity`,
   default conservative).
2. **Cross-map linkage** — at least N reads in the side-index link the two loci (a primary in A with
   a secondary in B). This is the empirical "these loci actually share reads" check.

A pair passing both becomes a family edge; connected components over family edges = families. A locus
in no family stays pure Layer 1 and never enters Layer 2.

> The similarity threshold is the principled-merge knob the advisor wants. It must be explicit,
> defaulted conservatively, and reported (how many families formed, at what similarity).

### C3. Family graph merge (the variation graph)

For each family, merge its members' **Layer-1 splice graphs** into one variation graph:
- homologous exons/junctions (aligned across copies by the similarity map) → **shared nodes**
- copy-specific bases/junctions → **private nodes**
- each copy is recoverable as a path; the reference-collapsed copy is the one whose path is currently
  thin.

This reuses the existing `FamilyGraph` machinery (`vg.rs` / `vg_family/`), **but sourced from
Layer-1 graphs instead of from bundles.** That is the one substantive change to the merge: its input
is `G_i` (clean per-locus graphs), not a re-bundled read set. No coordinate is invented — shared
nodes are an *alignment* between existing Layer-1 coordinates, and each emitted path carries its own
copy's real coordinates.

### C4. Secondary amendment + constrained flow-decomposition

With the family variation graph built:
1. **Amend** — fold the side-index secondaries overlapping the family in as candidate edges/transfrags
   (their junctions become candidate edges; their spans become flow on shared/private nodes). This is
   how an under-covered copy "borrows strength": a read ambiguous between copies adds flow to the
   *shared* backbone, while a copy-specific (PSV-bearing) read pins a *private* path.
2. **Decompose** — recover copies (paths) and isoforms (sub-paths) by flow-decomposition of the
   amended graph under allele-linkage constraints (PSV co-occurrence forbids inconsistent path
   combinations). This is the thesis "constrained flow-decomposition," and it reuses the existing EM /
   decisive-evidence / PSV machinery (`enumerate_diagnostic_sites` `vg.rs:4102`, `classify_family`
   `vg.rs:2198`, `compute_copy_ownership` `vg.rs:1857`, the restricted EM `em_family_qualifies`).

All-secondary regions (C5) are handled here too: a candidate region that exists *only* in secondaries
becomes a candidate new copy path, admitted only on decisive proof.

### C5. All-secondary regions → candidate new copies

A region every read of which is secondary is dropped by Layer 1 (baseline drops these as cross-map
artifacts — correct). Layer 2 reconsiders such a region **only** when the side-index shows decisive
evidence it is a real, separately-existing copy (PSV / decisive-ownership gate, the same proof bar as
C4). This is the user's *"bundles dropped if all reads are secondary… could be extra in vg mode if we
have enough proof they are real."* Gated, logged, off by default until validated.

### C6. Union-by-chain (additivity enforcement)

The final emission is `Layer-1 transcripts ∪ novel Layer-2 transcripts`, deduplicated by **intron
chain (coord signature)**:
- every Layer-1 transcript emits unchanged (identity preserved by signature match);
- a Layer-2 transcript emits **only if its chain is not already in Layer 1** (a genuinely novel copy
  or isoform);
- a Layer-2 path that re-derives a Layer-1 chain is silently folded into the Layer-1 transcript (no
  duplicate, no coordinate change).

This is the structural guarantee that Layer 2 can only ever *add*. Reuses the existing union-by-chain
idiom (`build_fresh_baseline_subbundle` / `RescueClass::UnionBaseline`, coord-signature diff harness).

---

## Properties (why this is safe by construction)

- **VG ⊇ baseline, always.** Layer 1 is byte-identical to baseline (Phase 1, proven). Layer 2 only
  unions in novel chains. There is no path by which a baseline transcript is dropped or moved.
- **Bundles never touched.** Secondaries live in the side-index, never `bundle.reads`. No
  re-bundling, no re-splitting, no mega-bundle merge (the original bug class is structurally gone).
- **Scoped.** Only loci passing C2 (similarity **and** cross-map linkage) enter Layer 2. Everything
  else is pure baseline.
- **Principled merge.** The similarity threshold (C2) is explicit and conservative — directly answers
  the advisor's "no similarity threshold" critique.
- **Built from clean graphs.** Layer 2 consumes Layer-1's graphs (baseline coords), not contaminated
  reads — the root cause of the region-drop cannot recur.

---

## What changes vs. today

| Concern | Today (pre-Layer-2) | Layer 2 |
|---|---|---|
| Secondary storage | (Phase 1) dropped from bundles | side-index (C1) |
| Family discovery | `build_multimap_index` ← `bundle.reads` | side-index cross-map links (C1/C2) |
| Subset selection | implicit (any multimap) | similarity **+** linkage gate (C2) |
| Family graph input | bundles (could contaminate) | Layer-1 splice graphs (C3) |
| Copy recovery (RABL2 1→7) | secondaries-in-bundle (removed) | amend + flow-decomp (C4) |
| New copies from all-sec regions | n/a | proof-gated candidates (C5) |
| Emission | VG could drop baseline | union-by-chain, additive only (C6) |

---

## Open decisions for review

1. **Side-index scope/size** (C1 ⚠) — chromosome-global vs pruned-to-family-candidates. Recommend
   prune after discovery's first pass; cap per-locus with a logged drop.
2. **Similarity threshold default** (C2) — what value of `--family-exon-similarity`, and is exon-only
   k-mer similarity the right metric vs. graph-node overlap? (Lean exon-only minimizers — already
   validated.)
3. **All-secondary new-copy gate** (C5) — exact decisive-proof bar; ship default-off until validated.
4. **Flagging** — is Layer 2 always-on under `--vg`, or behind an opt-in during development?
   Recommend: Layer-2 plumbing default-off behind a flag until C4/C5 validated genome-wide, with
   Layer 1 (baseline-identical) as the `--vg` default meanwhile.
5. **EM/decisive reuse** — confirm the restricted EM and PSV machinery operate unchanged on the merged
   family graph (they were written against bundle-sourced family graphs).

---

## Testing strategy

- **Invariant test (the floor):** for chr19 (and a chrY family chrom), VG-mode coord-signature ⊇
  baseline coord-signature, *always* — assert no baseline chain is ever missing. This is the
  regression that would have caught the original bug; it must be a standing test.
- **Phase-1 preservation:** with Layer 2 disabled, VG remains byte-identical to baseline and
  `RUSTLE_PRECISE=1` byte-identical to `4705ab1`.
- **Recovery test:** RABL2 family recovers its starved copy's isoforms (1 → 7) via the side-index +
  amendment path, *without* any secondary in `bundle.reads`. This proves the re-expression of the
  paused capability.
- **Additivity test:** every Layer-2 emission is either a Layer-1 chain (folded) or a chain absent
  from Layer 1 (novel) — never a modified Layer-1 chain.
- **Determinism:** `RAYON_NUM_THREADS=1 -p1` reproducible; side-index iteration order sorted/tie-broken.
- **Whole-genome:** per-chrom serial only (whole-genome `-L` OOMs ~18 GB).

---

## Thesis tie-in

Layer 2 *is* the thesis object (`project_thesis_framing_family_vg`):
- **The flip:** ambiguous secondaries are not noise to resolve (Canzar 2016 facility-location) but
  *shared evidence* about a family's common backbone — they add flow to shared nodes (C4).
- **The hard problem in his taste:** one variation graph per family; copies = paths, isoforms =
  sub-paths; recovery = minimum constrained flow-decomposition under allele-linkage constraints (C3/C4).
- **Borrow strength:** the under-covered copy shares 90%+ of its graph with paralogs; pool shared-node
  evidence, require only copy-specific (PSV) reads to pin divergences (C4) — the pooling-vs-separation
  tradeoff is the research question.
- **Principled merge:** the similarity threshold (C2) is the explicit, defensible knob the advisor
  demands — not arbitrary 1/k.

This design keeps the engineering invariant (VG ⊇ baseline) and the thesis goal (exploit ambiguity to
co-assemble a family) in the *same* architecture: Layer 1 guarantees safety, Layer 2 carries the
novelty.
