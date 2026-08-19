# Deferred O1 implementation: broad families and recent-copy subfamilies

**Status:** specified, not implemented (2026-08-16).

**Motivating evidence:** [`bench/o1_golga2_subfamily_audit.md`](../bench/o1_golga2_subfamily_audit.md).

This hierarchy is the first implementation slice of the richer block-aware provenance network
specified in [`o1_duplication_provenance_model.md`](o1_duplication_provenance_model.md). The hierarchy
can separate recent copies from broad homologs; the provenance model additionally represents mosaic
duplicons and ancestry without conflating a shared block with gene-family membership.

## Problem

One flat partition cannot express both biologically valid statements:

1. GOLGA2 is homologous to, and is the documented ancestral source closest to, GOLGA6/8.
2. GOLGA2 is not one of the recent high-identity GOLGA6/8 duplication copies.

Calling the chr9 node an unrelated false positive makes the audit artificially pure. Removing it
with a global RNA identity floor of 0.80 also damages named MAGEA and NBPF members. Conversely,
keeping it without qualification makes a recent-copy subfamily look impure.

## Proposed output contract

Preserve the current RNA O1 family as the broad homology object and emit a second, nested DNA-typed
subfamily view:

```text
broad_family_id       RNA forward E_r + current gamma quasi-clique
recent_subfamily_id   connected/cohesive block of high-identity genomic-span edges
copy_relation         RECENT_COPY | BROAD_ONLY | DNA_UNRESOLVED
```

- `RECENT_COPY`: the locus belongs to a DNA-supported subfamily containing at least two distinct
  physical loci.
- `BROAD_ONLY`: the RNA edge places the locus in the broad family, but it has no qualifying edge
  into a multi-locus recent-copy block.
- `DNA_UNRESOLVED`: the genomic-span substrate is missing, truncated, or otherwise inadequate;
  absence of a DNA edge is not treated as negative evidence.

For GOLGA, the intended result is:

```text
GOLGA6/8 chr15 loci  broad_family=GOLGA-wide  recent_subfamily=GOLGA6/8  RECENT_COPY
GOLGA2 chr9          broad_family=GOLGA-wide  recent_subfamily=.         BROAD_ONLY
```

The labels above describe relationships derived from sequence. Gene symbols are validation labels
only and must not enter the production decision.

## Proposed algorithm

1. Construct RNA nodes and broad families exactly as now. Do not change the shipping O1 partition.
2. Within each emitted broad family, extract the transcript-normalised complete genomic span for
   every locus—the same substrate already used by `--joint-dna-rna`.
3. Build `E_recent` with the DNA detector's existing recent-duplication operating point:
   identity ≥0.90, coverage of shorter span ≥0.50, either alignment orientation.
4. Partition `E_recent` with the same deterministic gamma-quasi-clique engine used by O1.
5. Emit only subfamilies containing at least two distinct physical loci. Label remaining broad
   members `BROAD_ONLY`, unless substrate-quality checks require `DNA_UNRESOLVED`.
6. Keep the RNA/DNA edge class (`RNA_DNA`, `RNA_ONLY`, `DNA_ONLY`) as evidence. It must not silently
   become an intersection rule for broad-family membership.

The 0.90 value is not a new fitted GOLGA threshold: it is the existing DNA recent-segmental-
duplication floor. It defines a different object from the sensitive RNA homology family and must
be named as such.

## Interface sketch

Add an opt-in flag provisionally named:

```text
--emit-recent-subfamilies
```

Proposed outputs:

```text
<out>.subfamilies.tsv
<out>.subfamily_edges.tsv
<out>.subfamily_rule.tsv
```

Suggested `subfamilies.tsv` columns:

```text
broad_family_id  copy_idx  recent_subfamily_id  copy_relation
dna_substrate_status  n_recent_edges  max_recent_identity  max_recent_coverage
```

The flag-off contract is strict: existing `families.tsv`, `copies.tsv`, `copies.fa`, family ids,
and partitions must remain byte-identical.

## Required safeguards

- RNA remains transcript-oriented and forward-only when that guard is requested; DNA remains
  orientation-agnostic so inverted segmental duplications are valid.
- Do not use chromosome identity as a rule. Genuine family copies can be cross-chromosomal.
- Do not use gene-name prefixes, annotations, Soto family ids, or the GOLGA2 coordinate as rules.
- Do not promote the RNA 0.80 counterfactual to a default: it reduced the HSA regional catalog
  from 165 to 151 family copies, lost two named MAGEA members, and split NBPF further.
- A missing DNA edge cannot distinguish divergence from incomplete node boundaries. Implement and
  certify the `DNA_UNRESOLVED` state before interpreting `BROAD_ONLY` as biological divergence.
- Keep segmental-duplication membership and multi-copy-gene-family membership distinct. A DNA
  subfamily is nested evidence about duplication recency, not proof of expression or gene status.

## Acceptance tests

### GOLGA positive discriminator

- Broad RNA family still contains the freshly emitted GOLGA2 node.
- Recent-copy output contains the intact 18-copy chr15 GOLGA6/8 block.
- GOLGA2 is `BROAD_ONLY`, not deleted and not assigned to that recent-copy block.
- SWI5 is not treated as the node identity merely because of its small boundary overlap.

### Cross-family safety panel

- Flag-off catalog is byte-identical.
- Broad-family partitions and the 72/75 named-node fresh-emission result are unchanged.
- MAGEA and NBPF broad-family recall does not inherit the losses caused by global RNA identity 0.80.
- RABL2, AMY1, PCDHB, TBC1D3, RANBP2/RGPD, RBMY1, and the remaining expanded panel receive explicit
  `RECENT_COPY`, `BROAD_ONLY`, or `DNA_UNRESOLVED` labels.
- The nine re-emitted NBPF repeat-bridge contaminants remain outside the broad NBPF target family;
  the subfamily layer must not reattach them through DNA repeats.

### Structural certificates

- Every emitted recent subfamily reports node count, edge count, density, lambda, identity/coverage
  floors, substrate, and orientation policy.
- Every `RECENT_COPY` assignment has a path in the archived `E_recent` graph.
- Every `BROAD_ONLY` assignment states whether usable DNA substrate was actually present.

## Implementation locations

- Reuse the RNA/DNA substrate construction in
  `src/rustle/vg_family/denovo_pipeline.rs::write_joint_rna_dna_certificate`.
- Reuse `family_split::gamma_quasi_clique_partition`; do not create a second clustering algorithm.
- Add CLI plumbing and output emission in `src/bin/gw_family_catalog.rs` after the broad catalog is
  emitted, mirroring the current report-only joint certificate.
- Extend the fresh O1 panel to score broad-family recall and recent-subfamily precision separately.

## Thesis claim boundary

Until implemented and validated, O1 should claim that Rustle emits a broad RNA homology family and
that a nested recent-copy view is specified but deferred. The current audit may show GOLGA2 as a
related outgroup for interpretation, but must not imply that production Rustle already emits the
proposed `recent_subfamily_id` or `copy_relation` fields.
