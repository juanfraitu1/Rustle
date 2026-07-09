# FamilyGraph models exon combinations a linear reference cannot

**Claim under test** (the thing the advisor credits the VG for, not O2 copy-assignment):
`FamilyGraph` (`src/rustle/vg_family/family_graph.rs`) represents a gene family as one graph
where nodes are exon-equivalence classes and each copy is a path through it. A shared exon
is ONE node touched by every copy that carries it; a copy-specific exon is a singleton-
contributor bubble (`ExonClass.copy_specific == true`). A linear reference has no such
shared node — it needs a separate transcript per exon combination, and every exon that
recurs across combinations gets physically duplicated once per transcript that uses it.

This is a demonstration/verification, not a new production feature.

## Mechanism, demonstrated

Test: `family_graph_models_exon_combinations` in
`src/rustle/vg_family/family_graph.rs` (`#[cfg(test)]`), built via the real builder
(`build_family_graph_from_layer1_graphs`) using the same input-construction pattern as
the existing `merge_from_layer1_graphs_shares_homologous_exon_and_has_edges` test
(synthetic Layer-1 splice graphs + a deterministic-hash genome via
`tests_support::make_layer1_graph` / `make_two_copy_genome`).

Toy family, 3 copies, deliberately combinatorial exon usage:

- copy A = exons [1,2,3]
- copy B = exons [1,2,4] (shares 1,2 with A; 4 is private)
- copy C = exons [1,5,3] (shares 1 with A/B, 3 with A; 5 is private)

Exon "1" is identical across all three copies; "2" is shared by A,B only; "3" is shared
by A,C only; "4" and "5" are copy-specific.

Run:

```
PATH=/home/juanfra/miniforge3/bin:$PATH cargo test --lib family_graph_models_exon_combinations -- --nocapture
```

Printed output:

```
=== FamilyGraph exon-combination demonstration ===
Copy A = exons [1,2,3]   Copy B = exons [1,2,4]   Copy C = exons [1,5,3]
--- Nodes (exon-equivalence classes) ---
  exon1 (shared A,B,C)   span=(100, 160) copy_specific=false contributors=3
  exon2 (shared A,B)     span=(300, 360) copy_specific=false contributors=2
  exon5 (private C)      span=(400, 460) copy_specific=true  contributors=1
  exon3 (shared A,C)     span=(500, 560) copy_specific=false contributors=2
  exon4 (private B)      span=(700, 760) copy_specific=true  contributors=1
n_nodes = 5   (shared=3, copy-specific bubbles=2)
n_edges (junctions) = 5
--- Reconstructed per-copy paths (= exon combinations, genomic-ascending) ---
  copy A path: [NodeIdx(1), NodeIdx(0), NodeIdx(2)]
  copy B path: [NodeIdx(1), NodeIdx(0), NodeIdx(3)]
  copy C path: [NodeIdx(1), NodeIdx(4), NodeIdx(2)]
n_distinct_paths = 3  (3 copies -> 3 exon combinations, all in ONE graph)
--- Linear-reference contrast ---
  Linear reference would need: 3 separate transcripts, 9 total exon-instances (shared exons 1/2/3 duplicated across the copies that carry them)
  FamilyGraph instead needs:   5 exon-class nodes + 5 junction edges representing all 3 combinations in ONE graph
====================================================
```

The test asserts the structure directly, not just the printout: exon "1" has
`per_copy_spans.len() == 3` and `copy_specific == false`; exons "2" and "3" have
`per_copy_spans.len() == 2`; exons "4" and "5" have `per_copy_spans.len() == 1` and
`copy_specific == true`; `n_nodes (5) < summed per-copy exon count (9)`; and the three
`recover_paralog_path` results are pairwise distinct (3 copies -> 3 distinct exon
combinations, all reconstructable from the one graph's `per_copy_spans`).

**Reading the numbers**: 3 copies x 3 exons each = 9 exon-instances a linear reference
would carry (one full transcript per combination, shared exons re-written into every
transcript that uses them). The FamilyGraph needs only 5 nodes — the compression comes
entirely from the 3 shared exons collapsing from 3+2+2=7 instances down to 3 nodes,
leaving 2 more nodes for the true singletons (4,5). 5 nodes + 5 edges represent all 3
combinations simultaneously; a linear model has no node that is "the same" across
transcripts, so it cannot express that copies A and C are identical at exon 3 without
either merging A and C into one transcript (losing B's distinctness) or writing exon 3's
sequence out twice.

## Real-data corroboration: GSTM

`bench/bundle_isoform_probe.tsv` (method: `bench/BUNDLE_ISOFORM_PROBE.md`,
`bench/bundle_isoform_probe.py`) compares, per known multi-copy family region, the
number of distinct **exact intron/exon chains** recovered from reads against
StringTie's per-locus bundle assembly:

| family | exact chains | StringTie chains |
|---|---:|---:|
| GSTM | 349 | 273 |

GSTM (glutathione S-transferase mu, a tandem paralog cluster) carries **349 distinct
exon/intron chains** across its copies versus 273 chains from StringTie's linear,
per-locus assembly — the same underlying phenomenon the toy test isolates: multiple
paralogs recombine a shared exon pool into different combinations, and a per-locus
linear assembler compresses/misses some of that combinatorial diversity relative to
following the exact chain each read supports.

**Honesty note** (carried over from `BUNDLE_ISOFORM_PROBE.md` itself, which already
flags this): part of the 349-vs-273 gap is exact-chain fragmentation of the same
underlying StringTie isoforms, not 76 net-new biological combinations — the source
document's own interpretation (#2) says so explicitly. The number is cited here as
real-data evidence that multi-copy family regions carry more distinct exon
combinations than a linear per-locus assembler reports, not as a literal count of
`FamilyGraph` paths (the probe does not build a `FamilyGraph`; it compares exact-chain
transcript counts to StringTie's GTF). The toy test above is what establishes the
mechanism precisely; GSTM is corroboration that the same shape of problem (more
combinations than a linear assembler natively expresses) shows up on real data.

## Info-sharing framing

Because a shared exon is one node, a read that lands on exon "1" is evidence usable by
every copy compatible with that node — copy A, B, and C simultaneously — rather than an
ambiguous multi-mapping the linear reference has to arbitrate by picking one
"best" transcript (or spreading 1/k weight, which the pipeline deliberately avoids
everywhere else). The multimapping read is *shared evidence across the paths it's
compatible with*, not noise to resolve away. Disambiguation (which path a given read's
own junction/PSV evidence supports) is the separate O2 copy-assignment problem this
demonstration does not touch — here we only show that the O1 structure the VG builds is
the one a linear reference cannot represent without either over-collapsing distinct
copies or duplicating shared content per combination.

## Verification

- `cargo test --lib family_graph_models_exon_combinations -- --nocapture`: passes, prints
  the table above.
- `cargo test --lib`: full suite green (see commit for the run this was verified against).
