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
