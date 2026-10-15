# De novo analogue of Liftoff's `-copies` (2026-10-14)

Pre-registration, written before any number below was computed. Code: `bench/template_copies/copies.py` at `f0833180`
(9 tests). Branch `denovo-copies` (on top of `template-copies`).

## Question

The user asked (2026-10-14): "in order to prove that the pipeline is also very similar to liftoff, should we include
something similar to the -copies functionality of liftoff?" and then "yes lets do that".

Liftoff's `-copies` lifts every annotated gene and keeps its extra placements at sequence identity ≥ `-sc`. The de novo
analogue lifts our assembled loci instead of the annotation: each base locus's representative transcript is aligned to
the genome, and every alignment that passes Liftoff's own thresholds is a projected copy. Two questions:

1. For the genes we assemble, does the projection find the extra copies Liftoff finds?
2. What does it find that Liftoff does not?

## Known before running (recorded now)

- **Truth.** Figure 8's Liftoff tables (`/mnt/linuxdisk/tmp/rustle_figures/liftoff/{human,gorilla}/liftoff_loci.tsv`;
  Liftoff v1.6.3, `-copies -sc 0.95 -a 0.5 -s 0.5 -overlap 0.1`, sharded, with figure 8's cross-shard rules M1-M3).
  Extra copies that are not short (exon union ≥ 200 bp): human 1,888, gorilla 1,129, all at sequence identity ≥ 0.95.
- **Templates.** The template-copies templates (one representative per base locus: the transcript with the most reads),
  built 2026-10-14: A119b_shipped 94,685 (projected: 16 PAF batches, minimap2 `-x splice -c -N 50 -p 0.1`, CHM13 v2.0
  splice index), testis 13,048 and KB3781 17,672 (built, not projected).
- **Seen already on A119b_shipped**, in the template-copies study (a different question): the number of alignments at
  identity ≥ 0.80 and aligned ≥ 0.50 (561,980), the candidate counts, the template-guided loci and their Liftoff scores
  (register 1285), including that 63 of 106 read-supported template-reachable copies were reached by an assembled
  locus. **Not seen**: any alignment at Liftoff's thresholds, any coverage of a Liftoff copy by an alignment, any site.
- **A design difference, expected and not measured.** Liftoff aligns each gene's genomic sequence (exons and introns);
  the projection aligns the spliced transcript. An intronless retrocopy of a multi-exon gene should therefore be easy
  for the projection and hard for Liftoff.

## Method (fixed)

**Projected copy.** An alignment of a template that passes Liftoff's thresholds, computed as Liftoff computes them
(`find_best_mapping.convert_all_children_coords`) over the template's bases from the PAF CIGAR:
- coverage = aligned template bases (M) / template length ≥ 0.5 (`-a`);
- identity = matching bases (PAF column 10) / (template length + target-only bases, CIGAR D) ≥ 0.95 (`-sc` of the
  figure-8 run).

Both are normalised by the whole template, as Liftoff normalises by the whole gene, so a projected copy matches ≥ 95% of
its template's bases.

**Recall.**
- A Liftoff extra copy (not short) is **reachable** when its source gene's placement (`in_place` or `moved`, not short,
  same `source_id`) is covered ≥ 50% by a base locus. The **own templates** are the templates of every such base locus:
  the analogue of lifting that gene.
- **Found (own)**: a projected copy of an own template covers ≥ 50% of the copy's exon union (figure 8's `MATCH_COV`).
  **Found (any)**: a projected copy of any template does.
- Copies whose source is not placed, or is placed but not assembled, are counted and never scored.
- **Misses**, in this order: an alignment of an own template covers the copy but fails a threshold (`below_threshold`,
  with its best identity and coverage); an own template has minimap2's cap of 51 alignments (`capped`); neither
  (`no_alignment`).

**Sites** (what the projection finds beyond its own loci).
- Off-locus projected copies: the template's own base locus covers < 50% of the alignment's exon union.
- They merge into **sites**: same contig, either strand, exon unions overlapping by ≥ 50% of the smaller (transitively).
- Each site's best alignment is the one with the highest identity, then coverage, then template id. Its exons class the
  site, in this order:
  1. a Liftoff extra copy covers ≥ 50% (`liftoff_copy`);
  2. a Liftoff gene placement (`in_place` / `moved`) does (`liftoff_gene`);
  3. a base locus does (`assembled`);
  4. none of them does (`none`).
- **Retro-like**: the best alignment is unspliced (one block) and its template has ≥ 2 exons.
- **Unexplained**: class `none` and not retro-like.

**Guard S0.** At least 0.90 of a sample's templates must have an alignment covering ≥ 50% of their own exons; otherwise
the index and the base GTF are not the same assembly, and that sample has no verdict until the index is fixed.

## Bars (per sample; species never pooled)

- **C1 (primary): recall (own) ≥ 0.80**, with ≥ 30 reachable copies. Fewer than 30: not powered, no verdict for that
  sample.
- **C2 (secondary): unexplained sites ≤ 0.10 of sites.**

Reported, not barred:
- recall (any);
- reachable share (reachable / all extra copies);
- recall by the source gene's number of Liftoff copies (1, 2-5, ≥ 6);
- the misses and their best identity and coverage;
- the templates at the alignment cap;
- each site class with its retro-like count;
- the templates with at least one site.

## Samples

| sample | role | templates (run dir under `/mnt/linuxdisk/tmp/o1_template_copies/`) | base GTF | index | truth |
|---|---|---|---|---|---|
| A119b | development | `A119b_shipped` (projected) | `runs/human_A119b/human_A119b.gtf` | `npip_ladder/idx/target.splice.mmi` | human |
| testis | held out (same genome, another sample) | `testis` | `runs/human_testis/human_testis.gtf` | the same | human |
| KB3781 | held out (gorilla) | `KB3781` | `runs/gorilla_KB3781/gorilla_KB3781.gtf` | `winloci_data/GGO.splice.mmi` | gorilla |

The held-out templates are projected with `run_tg.py project` (same code and settings as A119b) after this file is
committed.

## Readings (fixed now)

1. **C1 passes in A119b and in both held-out samples** → claimed: "for the genes we assemble, the de novo projection
   finds Liftoff's extra copies", quoted with the reachable share (the claim says nothing about genes we do not
   assemble). C2 adds "and what it finds beyond Liftoff is mostly annotated genes or retrocopy-like", only in the samples
   where C2 passes.
2. **C1 passes in A119b but fails in a held-out sample** → not claimed (sample-specific); all numbers reported.
3. **C1 fails in A119b** → not claimed; the miss decomposition says why. The held-out samples still run (they are cheap)
   and are reported.

## Rulings (made while designing)

- **R1.** The templates are the template-copies templates as built (no rebuild): one representative per base locus.
- **R2.** The projection settings are template-copies' (`-x splice -c -N 50 -p 0.1`; `-p 0.1` is the guided finder's
  setting). A119b's projection is reused as it is.
- **R3.** No read evidence is used anywhere: Liftoff's copies need none, and the copy does not need to be assembled.
- **R4.** The guided analogue (annotation records as templates) is G1, figure 8's separately registered test; it is not
  part of this one.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.
