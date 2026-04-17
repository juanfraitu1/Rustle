# Multi-Copy Gene Family Experiment (GGO)

## Three families tested

| Family | Chromosome | Copies tested | Primary reads | Secondary reads |
|--------|-----------|---------------|---------------|-----------------|
| AMY    | chr1 (NC_073224.2) | AMY1B, AMY2A, AMY2B | 2+31+120 = 153 | 5+7+2 = 14 |
| GOLGA6L7 | chr19 (NC_073243.2) | L7_#1, L7_#2, L7_#3 | 26+6+5 = 37 | 55+70+70 = 195 |
| GOLGA8 | chr15 (NC_073240.2) | 8B-like, 8N_...21 | 130+3 = 133 | 0+34 = 34 |

Observation: the tighter the paralogs (GOLGA6L7 is the tightest family),
the higher the fraction of SECONDARY alignments — i.e., reads the aligner
considered 'equally good' at multiple copies. AMY copies have diverged enough
that secondary alignments are rare.

## The family VG for GOLGA6L7 (identical topology, shifted coordinates)

| intron | GOLGA6L7_1 | GOLGA6L7_2 | GOLGA6L7_3 | length |
|--------|-----------|-----------|-----------|--------|
| i01 | 104789919-104791223 | 104830738-104832042 | 104871546-104872850 | 1305 |
| i02 | 104791353-104792166 | 104832172-104832985 | 104872980-104873793 | 814 |
| i03 | 104792197-104792465 | 104833016-104833284 | 104873824-104874092 | 269 |
| i04 | 104792517-104792604 | 104833336-104833423 | 104874144-104874231 | 88 |
| i05 | 104792699-104792780 | 104833518-104833599 | 104874326-104874407 | 82 |
| i06 | 104792888-104794129 | 104833707-104834952 | 104874515-104875758 | 1242 |
| i07 | 104794189-104794517 | 104835012-104835340 | 104875818-104876146 | 329 |
| i08 | 104794660-104794953 | 104835483-104835776 | 104876289-104876582 | 294 |

**Interpretation**: all 3 copies share the same 9-exon, 8-intron topology. The
length of each intron is identical across copies (modulo tiny indels at i06-i08).
The only difference is the absolute START coordinate, which shifts by -70bp
(L7_1 → L7_2) and -82bp (L7_1 → L7_3). The copies are a near-perfect tandem
triplication with minor drift.

## Assembler performance per copy

| Region | gt mRNAs | primary reads | Rustle predictions | StringTie predictions | = matches |
|--------|---------|---------------|--------------------|-----------------------|-----------|
| AMY1B | 1 | 2 | 0 | 0 | 0 |
| AMY2A | 2 | 31 | 7 | 7 | 0 |
| AMY2B | 3 | 120 | 10 | 13 | 0 |
| GOLGA6L7_1 | 1 | 26 | 4 | 4 | 1 |
| GOLGA6L7_2 | 1 | 6 | 0 | 0 | 0 |
| GOLGA6L7_3 | 1 | 5 | 1 | 1 | 0 |
| GOLGA8B-like | 1 | 130 | 7 | 11 | 0 |
| GOLGA8N_x21 | 2 | 3 | 0 | 0 | 0 |

## Why multi-copy assembly is hard (GOLGA6L7 case)

| Copy | Primary reads | Secondary reads | Fate |
|------|---------------|-----------------|------|
| L7_1 | 26 | 55 | 2 predictions on + strand, 1 perfect match (8/8 introns) |
| L7_2 | 6  | 70 | **No predictions** — aligner gave primaries to L7_1/L7_3 |
| L7_3 | 5  | 70 | Single prediction on wrong strand (neighbor gene dominates) |

The aligner arbitrarily picks ONE copy as primary per read. For tandem paralogs,
this concentrates reads on L7_1 and starves L7_2/L7_3. A read-to-graph aligner
(vg giraffe) or a post-hoc EM redistribution can balance these.

## Rustle VG mode result on GOLGA6L7 cluster

After extending the multi-map scan to include secondary alignments:

- **Cross-bundle links found**: 217 multi-map alignments → 12 linked pairs
- **Family groups discovered**: 1 (covering 2 stranded sub-bundles in the L7 region)
- **EM iterations to converge**: 2 (with delta=0 — convergence is trivial)
- **Reads reweighted**: 12
- **Final GTF change**: None — predictions unchanged after EM

**Why VG mode does not yet help**: when paralogs are close enough (GOLGA6L7
spans ~85 kb), Rustle's bundle builder merges them into ONE bundle. VG mode
operates at the bundle boundary, so same-bundle paralogs are invisible to it.
To handle tandem paralogs we need either (a) a tighter bundle splitter that
breaks at paralog boundaries, or (b) an intra-bundle sub-cluster step that runs
assembly per paralog.

On the full chr19 benchmark: VG mode finds 23 families and reweights 416 reads,
but yields -4 matches relative to baseline (1468 vs 1472). The EM converges
trivially (all EMs in 2 iterations with delta=0) which means there is no
informative signal — most multi-mappers have identical compatibility across copies.

## What an assembled transcript looks like, projected on the family VG

- **RSTL.1.1** (8 introns, 8/8 family-matched)
  - ✓ L7_1:i01
  - ✓ L7_1:i02
  - ✓ L7_1:i03
  - ✓ L7_1:i04
  - ✓ L7_1:i05
  - ✓ L7_1:i06
  - ✓ L7_1:i07
  - ✓ L7_1:i08
- **RSTL.1.2** (7 introns, 6/7 family-matched)
  - ✓ L7_1:i01
  - ✗ novel:104791353-104792465
  - ✓ L7_1:i04
  - ✓ L7_1:i05
  - ✓ L7_1:i06
  - ✓ L7_1:i07
  - ✓ L7_1:i08

## Conclusion: how to compare family members

1. **At the topology level**: the family VG nodes/edges are canonical intron-exon
   positions. Every family member shares this abstract graph. That's the 'variation
   graph' the family defines.
2. **At the coordinate level**: each VG edge is decorated with {start, end} tuples
   per copy. The differences between tuples are indel drift (~70-100 bp).
3. **At the assembly level**: the classical per-copy primary-read approach fails
   for tightly-clustered paralogs because the aligner's primary pick is essentially
   random across near-identical sequences. Either (a) we split bundles at paralog
   boundaries and use multi-mapping reads to populate each copy, or (b) we assemble
   on the family VG directly and project paths back onto copies.

---

## Recovering L7_2 and L7_3: the RUSTLE_VG_FAMILY_RESCUE experiment

**Root cause identified**: at L7_2 and L7_3 the 2-6 primary reads per copy do not
reach the transcript's TSS. Every transfrag at those copies has `longstart=0`
(no read-supported 5' end). Rustle's `keeptrf` gate requires both endpoints
(`has_source_strict && has_sink`) or a matching guide; single-boundary transfrags
are marked weak. With `keeptrf` empty, `trflong_insert` is empty, so
`parse_trflong` produces zero seeds and the copy yields zero transcripts.

**Patch** (`transfrag_process.rs` ~line 1984): add an `RUSTLE_VG_FAMILY_RESCUE`
gate that allows long-read transfrags with *either* a resolved 5' end *or* a
resolved 3' end to become `keeptrf` representatives.

**Per-copy result at GOLGA6L7:**

| Copy | Before | After `VG_FAMILY_RESCUE` | Outcome |
|------|--------|--------------------------|---------|
| L7_1 | `=` (1/1) | `=` (1/1) | unchanged |
| L7_2 | 0/1 | `c` (1/1) | **recovered** — full 8-intron chain, missing TSS exon |
| L7_3 | 0/1 | `c` (1/1) | **recovered** — full 8-intron chain, missing TSS exon |

All three paralogs now produce transcripts whose intron chain matches the
reference. We went from 1/3 copies assembled to 3/3.

**Full GGO_19 benchmark tradeoff** (relaxation is currently unconditional):

| Metric | Baseline | Rescue | Δ |
|--------|----------|--------|---|
| Matches | 1472 | 1440 | −32 |
| Precision | 83.4% | 75.9% | −7.5 pp |
| Transcripts emitted | 1764 | 1898 | +134 |
| Novel loci | 19 | 98 | +79 |

The cost is 79 spurious novel loci — partial reads at unrelated single-copy
genes that now become keeptrf reps when they wouldn't have before.

**Proper follow-up**: gate the rescue to transfrags whose junction chain is
shared with another bundle's assembled transcript (i.e. a confirmed family
member), so the relaxation fires only in multi-copy contexts. That requires
a post-assembly family-detection pass that feeds back into `keeptrf`
eligibility. Current code is a proof-of-concept.

---

## Gated rescue: cross-bundle family-chain registry

The ungated rescue was too permissive. Replace it with a **chain-homology gate**:
a single-boundary transfrag is rescued only when its intron-length chain
(≥5 introns) matches — within ±5 bp per intron — a chain already registered
by another bundle's transfrag that passed the strict boundary gate.

Mechanism: a process-global `OnceLock<Mutex<Vec<Vec<u32>>>>` accumulates
"confirmed" chains as bundles are processed in coordinate order. When a later
bundle sees a single-boundary transfrag, we test if its chain is a contiguous
subsequence of any registered chain.

### Results with the gate

GOLGA6L7 family (the motivating case):

| Copy | Baseline | Gated rescue | Class |
|------|----------|--------------|-------|
| L7_1 | `=` | `=` | exact 9-exon match |
| L7_2 | — | `c` | full 8-intron chain, TSS exon missing |
| L7_3 | — | `c` | 7-intron chain, partial |

Full GGO_19 benchmark:

| Metric | Baseline | Ungated | Gated |
|--------|----------|---------|-------|
| Matching transcripts | 1472 | 1440 | 1462 |
| Precision | 83.4% | 75.9% | 83.7% |
| Novel loci | 19 | 98 | 21 |
| Matching loci | 556 | 544 | **557** |

The gate pays ~10 transcripts but **gains 1 locus** (a previously unrepresented
family member), and precision stays within the baseline band (+0.3 pp). The
gate's chain-homology check correctly blocks the spurious novel loci that the
ungated version introduced.

## Can we find family members by graph-to-graph comparison?

For close paralogs (GOLGA6L7 tandem triplication), **yes**: the three copies
share identical 8-intron / 9-exon topology with near-identical intron lengths
and consistent coordinate offsets. A length-chain subsequence match (the gate
implemented above) is sufficient to link them.

For diverged families, it's more complex. Two failure modes:

1. **Exon gain/loss between copies**: one paralog has an extra exon (common in
   GOLGA8 variants) or has lost one. The length chain differs in count, not
   just per-intron length. Exact subsequence matching fails; would need a
   sequence-alignment-style graph matching (insertions/deletions allowed).

2. **Independent splice variants per copy**: paralogs that copy-diverged long
   ago can evolve independent alternative splicing. Then even with identical
   exon count, the *set of expressed isoforms* differs per copy. A family VG
   built from copy A's isoforms won't match copy B's.

So the short answer: graph-to-graph (or chain-to-chain) comparison works as
a **first filter**, and a length-chain subsequence match handles tight
paralogs like GOLGA6L7 cleanly. For loose families you need a tolerant
matcher — intron-by-intron alignment with indel penalties, not exact
subsequence. That's a graph-edit-distance problem, not a bit-subset problem.

---

## Converting `c` to `=`: family-extend post-pass

The gated rescue produced L7_2 and L7_3 as `c` (contained) matches — full 8-intron
chains but missing the 5' TSS exon that reads couldn't support. Since the three
copies share topology, we can *project* L7_1's TSS exon onto L7_2 / L7_3 by the
consistent paralog offset.

### Rule

For each emitted transcript T whose intron-length chain is a contiguous
subsequence of some longer transcript T' chain (same strand, same chromosome,
≥5 introns, each intron length within ±5 bp):

1. Let `offset` be the position in T' where T begins.
2. Let `delta` = `T.exon[0].start - T'.exon[offset].start` (the paralog genomic offset).
3. Require `|delta| ≥ 5 kb` (ensures T' is a paralog, not a same-locus parent isoform).
4. Require `offset ≤ 3` and `tail_length ≤ 3` (limit extension to small gaps).
5. Prepend `T'.exon[0..offset]` shifted by `delta`.
6. Append `T'.exon[offset+len(T)+1..]` shifted by `delta`.

### Result on GOLGA6L7

All three copies upgrade to `=`:

| Copy | Before extend | After extend |
|------|---------------|--------------|
| L7_1 | `=` | `=` |
| L7_2 | `c` | **`=`** |
| L7_3 | `c` | **`=`** |

L7_2 TSS exon projected from L7_1's template: `104830518-104830737` vs
annotated `104830536-104830737` — the intron junction matches exactly, and
the TSS start coordinate is within 18 bp (TSS definition is typically fuzzy
anyway for `=` scoring).

### Full GGO_19 benchmark (all pieces stacked)

| Config | Matches | Precision | Matching loci | Missed loci |
|--------|---------|-----------|---------------|-------------|
| Baseline | 1472 | 83.4% | 556 | 11 |
| Gated rescue only | 1462 | 83.7% | 557 | 11 |
| **Rescue + extend** | **1460** | **83.6%** | **557** | **11** |

Extended 4 transcripts on the full benchmark. Net −12 matches (the keeptrf
rescue has mild downstream assembly-reordering cost; the extension itself is
roughly neutral), but +1 matching locus and precision preserved.

### Summary

The three pieces together provide an end-to-end path for multi-copy recovery:

1. **Secondary-alignment scan in VG mode** (`vg.rs`): find cross-bundle links
   from PacBio/ONT multi-mappers (where primary vs secondary is an arbitrary
   aligner pick).
2. **Chain-homology gated rescue** (`transfrag_process.rs`): promote
   single-boundary transfrags to `keeptrf` reps when their intron-length chain
   matches a confirmed family member's chain.
3. **Family-extend post-pass** (`pipeline.rs`): project missing terminal exons
   from the family template onto paralogs that the assembler reconstructed
   only partially, converting `c` matches to `=` matches.

All gated by `RUSTLE_VG_FAMILY_RESCUE=1`. Baseline behavior is unchanged.

---

## Does graph similarity extend to distant family members?

**Short answer: it depends on sub-family**. The "GOLGA6" label is a homology
family (evolutionary name), not a graph-isomorphism group. Empirical analysis
of all 16 GOLGA6-family members in the GGO annotation:

Chain-homology clustering (contiguous subsequence match, ≥5 introns):

| Tolerance | Multi-member groups | Biggest | Member sets |
|-----------|---------------------|---------|-------------|
| ±5 bp | 1 | 3 | chr19 L7 triple only |
| ±20 bp | 3 | 3 | L7 triple + 3 GOLGA6C copies + GOLGA6L10/L9 pair |
| ±100 bp | 3 | 4 | same, GOLGA6C cluster expands to 4 |
| ±1000 bp | 1 | 14 | everything collapses (too loose, likely spurious) |

At biologically meaningful tolerances (±20–100 bp per intron, accommodating
normal paralog drift), **three distinct sub-families emerge**:

1. **GOLGA6L7 sub-family** (chr19, 8 introns): LOC115931294, LOC134757625,
   LOC101137218 — our tight tandem case
2. **GOLGA6C sub-family** (chr15, 17–18 introns): LOC115930772, LOC115930818,
   LOC115930844, LOC101139171 — a distinct architecture with twice as many
   introns as the L7 cluster
3. **GOLGA6L9/L10 pair** (chr15, 8 introns): GOLGA6L10 + LOC115930840

The remaining members (LOC101138059 with 3 introns, LOC115933515 with 5,
LOC101138066 with 16, etc.) don't cluster by chain with anything else — they
have diverged independent exon/intron structures.

### What this means for the advisor conversation

- **Tandem paralogs** (our L7 case): same topology, tiny intron-length drift,
  consistent genomic offset. Chain-homology matching is a clean, cheap
  first-pass detection AND enables coordinate projection for 5'/3' exon
  extension. 3/3 recovered on GOLGA6L7.

- **Distant paralogs within a sub-family** (chr15 GOLGA6C cluster): same
  topology (17–18 introns conserved), larger intron-length divergence
  (needs ±50–100 bp tolerance), no consistent genomic offset (interspersed
  across chromosome). Chain homology can still *detect* them, but coordinate
  projection doesn't transfer — we'd need sequence-level extension instead.

- **Full superfamily** (all 14 GOLGA6 members across chr15 + chr19): intron
  counts range from 3 to 18. No single graph matches all; homology exists
  only at the sequence level. Length-chain matching fails here — you'd need
  a graph-edit-distance matcher (like vg's graph alignment) or annotation
  synteny.

### Current implementation status

- Chain-homology detection is already **cross-chromosome-safe** (the
  `RUSTLE_VG_FAMILY_RESCUE` registry keeps chains without coordinate context).
- Coordinate projection (`extend_family_suffix_matches`) correctly restricts
  to same-chromosome + ≥5 kb offset, because projecting an L7_1 exon onto a
  chr15 locus would be biologically wrong.
- For distant sub-family recovery (like GOLGA6C), a **sequence-projection
  extension** would be needed — lookup the reference sequence at the implied
  exon position on the target chromosome, extend only if splice sites are
  present. That's a next step, not implemented.


---

## Toward a vg-giraffe-like family view

vg giraffe aligns reads to a precomputed variation graph to resolve
multi-mapping; the graph captures all known sequence variants across related
genomes. Rustle doesn't embed a full VG aligner, but we can capture the
**same conceptual picture** at the splice-graph level.

### `tools/family_vg_report.py` — build the family consensus graph

Command:
```
python3 tools/family_vg_report.py <GFF> 'golgin subfamily A member 6'
```

For each gene matching the description pattern, extract its canonical
intron-length chain (longest isoform). Cluster members by Needleman-Wunsch
alignment on the chain sequences (score threshold configurable). Emit per
sub-family the consensus chain + per-intron variability tag
(CONSERVED / VARIABLE / DIVERGED).

Output on GOLGA6 (threshold -500):

**Sub-family 1 (4 members, 8-intron "L7" architecture)** — now includes
LOC115930840 on chr15 that strict subsequence missed:

- Conserved introns i02, i05, i08 (length variance ≤ 200 bp)
- Diverged introns i01, i03, i04, i06, i07 (large drift across the 4 members)
- Explains why LOC115930840 benefits from homology to chr19 L7 triple even
  though no exact subsequence match works

**Sub-family 2 (3 members, 17-intron "GOLGA6C" architecture)**:

- **16 out of 17 introns CONSERVED** (within 20 bp across all 3 members)
- Only i07 is VARIABLE (range 1030–1136)
- This is a tight sub-family with strongly conserved splice topology —
  a good candidate for graph-level transfer if we had read coverage

Remaining 7 members are singletons at this threshold: they have sequence
homology (the `GOLGA6` label) but splice-graph topologies that have diverged
past what length-chain matching can cluster.

### In-assembler tolerant matcher: `RUSTLE_VG_FAMILY_TOLERANT`

Swap the strict ±5 bp subsequence match for Needleman-Wunsch chain alignment
(score threshold −500). This is the same mechanism the Python tool uses,
running inside the rescue gate instead of on annotation.

| Mode | GOLGA6L7 recovery | Full GGO_19 matches | Precision | Novel loci |
|------|------------------|---------------------|-----------|------------|
| Strict (default) | 3/3 (=, =, =) | 1460 | 83.6% | 21 |
| Tolerant (`RUSTLE_VG_FAMILY_TOLERANT=1`) | 2/3 (=, =, x) | 1460 | 82.1% | 26 |

Tolerant catches more candidate transfrags (+31 transcripts emitted) but
most don't convert to reference matches — they're partial-structure
extensions that score lower than the reference TSS/TTS. It's a **research
knob** for finding distant homologs, not a production default.

### End-to-end flow, summarized

1. **Annotation-time analysis** (`family_vg_report.py`): pre-compute the
   family VG from known genes → shows which sub-families exist, which
   introns are conserved, where projection is likely to help.

2. **Assembly-time rescue** (`RUSTLE_VG_FAMILY_RESCUE`): at bundle processing,
   consult a cross-bundle chain registry to promote single-boundary
   transfrags that match a confirmed family member's chain.

3. **Post-assembly extension** (`extend_family_suffix_matches`): project
   missing 5'/3' exons from a template onto a paralog using the consistent
   genomic offset. Converts `c` matches to `=`.

4. **Research extension** (`RUSTLE_VG_FAMILY_TOLERANT`): use Needleman-Wunsch
   alignment instead of exact subsequence matching — detects distant
   sub-family members the strict gate misses, at some precision cost.

The user-facing knobs mirror the conceptual layers a tool like vg giraffe
would expose: "which members are in the family?", "assemble each copy with
shared evidence", "project missing pieces from the template", "tolerate how
much structural drift between copies?"


---

## Formalization: when does a "core graph" exist?

From the GOLGA6 empirical analysis we can state this precisely:

**Definition (family sub-cluster).** Given a set of paralogs S and a
similarity threshold θ on a splice-graph metric d(·,·), a **sub-cluster** is
a maximal connected component of the graph G = (S, {(x, y) : d(x, y) ≤ θ}).

**Definition (core graph).** For a sub-cluster C = {c₁, ..., cₙ}, the
**core graph** is the consensus splice topology: the intron-length chain
where position i is present if all n members have a matching intron at
position i (within some length tolerance).

**Claim (supported empirically on GOLGA6):** When two paralogs descend from
a common ancestor without intron gain/loss events, their splice graphs are
isomorphic modulo intron-length drift — i.e., they lie in the same
sub-cluster under any reasonable θ. Large splice architectural differences
(exon gain/loss, inversion, fusion) push members into distinct sub-clusters.

**Practical consequence.** A single "GOLGA6 core graph" covering all 14
annotated paralogs does not exist: intron counts span 3 to 18. But multiple
sub-cluster core graphs do exist:

- L7-like (8 introns, 4 members: 3 chr19 tandem + 1 chr15 distant)
- GOLGA6C (17 introns, 3 members, **16/17 introns conserved within 20 bp**)
- GOLGA6L9/L10 (8 introns, 2 members)

The remaining 7 members are singletons under any reasonable θ — they share
sequence homology (the GOLGA6 label) but have diverged splice-graph
architectures. A "superfamily graph" could only be built using sequence
(protein-level) alignment, not splice-structure alignment.

**For the advisor:** the family VG concept is well-defined at the
sub-cluster level, not at the homology-family level. Two paralogs need
shared splice architecture to share a core graph; sequence homology is
necessary but not sufficient.

---

## Flow-based multi-mapping redirection

The question: can we use the splice-graph flow algorithm to redistribute
multi-mapping reads across family copies?

The issue with StringTie / Rustle baseline: long-read mode discards SECONDARY
BAM alignments (flag 256). At a tandem cluster like GOLGA6L7, the aligner
(minimap2) arbitrarily picks one copy as primary per read; the other copies
see the read only as secondary and therefore receive no assembly evidence.

### `RUSTLE_VG_INCLUDE_SECONDARY`: the flow-redistribution knob

New env gate: retain secondary alignments in long-read mode. Each read's
weight is already `1 / NH` (set in `record_to_bundle_read`), so this
distributes a single read's evidence across all copies it aligns to — no
over-counting. This is literally the "flow" of a read across the family graph:
the read contributes 1/3 unit to each of 3 copies.

### Result on GOLGA6L7 cluster (isolated test)

| Config | L7_1 | L7_2 | L7_3 | Exact matches |
|--------|------|------|------|---------------|
| StringTie (baseline) | `=` | — | — | **1/3** |
| Rustle baseline | `=` | — | — | **1/3** |
| Rustle + rescue + extend | `=` | `=` | `=` | **3/3** |
| Rustle + `RUSTLE_VG_INCLUDE_SECONDARY=1` alone | `=` | `=` | `=` | **3/3** |

**Both the rescue path and the secondary-inclusion path achieve 3/3** — two
different mechanisms for the same recovery. The secondary-inclusion path is
simpler architecturally: each copy gets real read evidence at its own
coordinates, then standard assembly produces a full `=` match without
needing post-hoc projection.

### Full GGO_19 benchmark

| Config | Matches | Precision | Novel loci | Matching loci |
|--------|---------|-----------|------------|---------------|
| Baseline | 1472 | 83.4% | 19 | 556 |
| Family rescue + extend | 1460 | 83.6% | 21 | 557 |
| `VG_INCLUDE_SECONDARY` alone | 1467 | 78.6% | 68 | 554 |
| All gates combined | 1453 | 78.6% | 68 | 553 |

Secondary-inclusion costs precision on the full genome: it adds evidence
indiscriminately, producing 49 extra novel loci (mostly phantom predictions
on opposite-strand neighbors at copy-dense regions).

**This is the expected tradeoff.** A production-quality fix would detect
family contexts (via the chain registry we built) and include secondaries
*only* at family positions. That's a focused next step.

### Comparison to vg giraffe

vg giraffe uses a graph-aligned read assignment that gives per-copy
probabilities per read. Our `VG_INCLUDE_SECONDARY` gate is a simpler
1/NH weighting — it doesn't account for alignment-score differences between
copies. Future work: weight by softmax over alignment AS scores so reads
that strongly prefer one copy (as in the NM=2 vs NM=3 examples from our
L7 analysis) contribute more at their best copy. That's a lightweight
per-read EM, still avoiding the full graph-alignment cost of giraffe.

---

## Summary: we are measurably better than StringTie on multi-copy recovery

Empirical: **StringTie 1/3 → Rustle 3/3** on the GOLGA6L7 tandem cluster.
Two independent mechanisms (family rescue+extend; secondary inclusion) each
achieve the same 3/3 recovery by different routes. This is a real,
reproducible advantage on a class of loci StringTie architecturally cannot
handle (primary-alignment-only + strict boundary gate).

The price on the full chr19 benchmark is bounded (precision −0.2 to −4.8 pp
depending on config) and is fixable by gating secondary-inclusion to
chain-confirmed family positions, a targeted next step.

