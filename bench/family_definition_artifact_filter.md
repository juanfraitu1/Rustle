# Refining the multi-copy-family definition — exclude ARTIFACTS, not non-coding (the VG does it)

**Goal correction (user):** the target is to drop **artifacts** (repeats, low-complexity, TEs,
chimeras, over-merges), NOT to require coding. A real **lncRNA multi-copy family is a real family** —
requiring an ORF would wrongly delete it. So the filter must be coding-agnostic + intrinsic.

## Two things that DON'T work (measured, `refine_family_definition.py`)
- **ORF is a bad filter.** On de-novo loci labelled by biotype: **76% of lncRNA have an ORF ≥100 aa**
  (spurious ORFs in long transcripts) and 6% of *coding* loci don't. ORF *keeps* lncRNA and *drops*
  coding — it conflates non-coding with artifact.
- **k-mer (4-mer) complexity SATURATES.** Median 0.988 across all biotypes (only 256 possible 4-mers →
  any transcript >300 bp hits the ceiling); the repeat artifacts also score 0.96. Not discriminative.

## The fix: artifacts are TOPOLOGICAL — the variation graph exposes them (no ORF, no annotation)
The same VG object that *defines* the family and *assigns* the reads also *filters* the artifacts —
the thesis's graph does triple duty, with no ad-hoc thresholds or coding assumptions.

1. **Repeat / tandem / TE → a CYCLE in the variation graph.** A real family is an *acyclic* bubble-chain
   (backbone + localised PSV bubbles); a repeat makes the graph loop back on itself. Intrinsic signal =
   **long-k-mer self-recurrence** within a transcript. Measured (`15-mer recurrence`, 2,000 de-novo tx):
   median **0.000** (real = unique k-mers, acyclic), tail to **0.57** — flags internal-repeat
   transcripts (e.g. `DN_NC_073224.2_130224868_27`, 4 kb, 57% recurrence) that 4-mer complexity (0.96)
   called "complex." **Sensitive exactly where the naive measure saturates.**
2. **Over-merge / repeat-bridge → a SPARSE, non-clique component.** A real family is a *clique* in the
   ~B backbone graph (every pair aligns reciprocally); an artifact bridge is a sparse connector through
   an articulation point (the comp-238 lncRNA blob). Intrinsic signal = graph density / articulation
   points / community detection (partly shipped: 238→57).
3. **Chimera → an INCOHERENT backbone.** A chimeric transcript's pieces map to unrelated loci → no
   single colinear spine in the VG. Intrinsic signal = backbone coherence (one dominant colinear path).

## Refined definition (coding-agnostic, intrinsic, topological)
> A multi-copy gene family = **~R ∩ ~B** AND the family's variation graph is **(a) acyclic** (no
> internal-repeat cycle — `k-mer self-recurrence < ~0.10`), **(b) a clique** (dense reciprocal ~B, not a
> sparse articulation bridge), and **(c) copy-count-bounded** (2 ≤ copies ≤ ~12; dispersed TEs make many
> copies). No ORF requirement → real lncRNA families are KEPT (they are acyclic, clean cliques);
> repeats/TEs/chimeras/over-merges are DROPPED on graph structure alone.

This keeps the definition purely topological (interest I) and unifies it with copy-assignment
(interest II): one variation graph, three jobs — define, assign, and filter — no coding model, no
arbitrary sequence thresholds.

> **Scope of the "clique" / "triple-duty" claim (L2 — now partly closed).** The clique density bar is
> enforced as a **DROP on BOTH** pipelines: the DNA cDNA-homology manifest (`make_dna_family_manifest.py`,
> `overmerge_sparse` at `n≥4 ∧ density<0.30`) AND the de-novo split path
> (`src/rustle/vg_family/family_split.rs` + `denovo_family_split.py`) — whose `web` density bar is now
> **aligned to 0.30** (was 0.15), and `Web` families are **excluded from copy-assignment**
> (`denovo_pipeline.rs`). Measured effect: the de-novo split goes 695 family / 3 web → **691 / 7**; the 4
> newly-dropped are the large multi-chromosome domain-sharing over-merges (DSFAM0 = a 164-member ZNF
> across 19 chromosomes, plus ZNF/APOL families) — the exact blobs that embarrassed "family = clique."
> `web_min_size` stays 10 on the de-novo side (vs the manifest's 4) **on purpose**: a SMALL sparse group
> can be a real divergent family (e.g. a 7-copy single-chrom MAGEB at density 0.24), so only
> *large-and-sparse* is dropped. Still partial: only the density/cliqueness leg is wired on the de-novo
> side — the k-mer-self-recurrence (cycle) leg remains cDNA-manifest-only (L5).

## Wired into the pipeline (`family_def_artifact_filter.py` → `make_dna_family_manifest.py`)

SHIPPED as a default-on filter with an opt-out (mirrors the retrocopy filter). Per family it computes
**cliqueness (density), copy-count, and k-mer self-recurrence** (genomic, capped 20 kb) and drops:
`repeat_cyclic` (recur ≥ 0.15), `overmerge_sparse` (n ≥ 4 and density < 0.30), `te_highcopy` (a high
backstop, n > 80, dense). **DENSITY — not copy-count — is the discriminator**: a real family is a dense
clique *at any size*, so a 41-copy clique (density 1.0) is KEPT and only sparse bridges are dropped.

Result on real data: **1,216 → 1,167 real families; 49 dropped** (39 repeat_cyclic + 10
overmerge_sparse — incl. the comp-238 lncRNA-repeat blob, density 0.117). Large real cliques (n =
40–47) correctly KEPT. Every family is written to the sidecar
`dna_family_artifact_class.tsv` (`family_id  class  n  density  recur  kept`).

**Opt-out / visualise:** `FAMILY_DEF_NO_ARTIFACT_FILTER=1` keeps *all* families (artifacts included) but
still tags each with its class — so you can pull a `repeat_cyclic` / `overmerge_sparse` / lncRNA family
from the sidecar and **show its graph** (a repeat family's variation graph has *cycles*; an over-merge
is a *sparse bridge*, not a clique) next to a clean `real` family. Thresholds are env-tunable
(`FAMILY_DEF_MIN_DENSITY`, `FAMILY_DEF_REPEAT_RECUR`, `FAMILY_DEF_MAX_COPIES`). Coding-agnostic
throughout — real **lncRNA** families pass (acyclic clean cliques) and are KEPT.

Artifacts: `family_def_artifact_filter.py` (filter + self-test) · `make_dna_family_manifest.py` (wired) ·
`dna_family_artifact_class.tsv` (sidecar) · `refine_family_definition.py` (the ORF/4-mer negative).
