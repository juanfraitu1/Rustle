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

## Amendment 1 (2026-10-14, speed only, before any result)

The first development run (A119b) was stopped by its 10-minute limit before it wrote anything. Profiling the stages
showed: inputs 2 s, own templates 1 s, the PAF pass 12 s, recall 1 s; the site merge did not finish. Seen while
profiling, and nothing else: 489,897 projected copies (self alignments included), 85,729 kept alignments of 24,801 own
templates, and projected-copy spans of 5,050 bp (median), 106,416 bp (99th percentile) and 1,640,982 bp (maximum). No
recall, miss, site or class number was computed or seen.

The cause is the merge's sweep over spans: a copy with long introns stays open across its whole span, so the sweep
compares nearly every pair. Diagnosing it showed where the pairs come from. Batch 2 holds 6,000 templates from one
tandem array (chr13 13.9-15.5 Mb, assembled as ~5,000 overlapping loci), and 5,035 of them reach minimap2's cap of 51
alignments, all inside the array. That gives 252,592 off-locus copies in batch 2, and 251,540 of them chain into one
cluster when any shared exonic base links two copies.

The merge (`copies.site_groups`) now sweeps exon blocks and keeps the open blocks per component. A block is compared
only with other components' open blocks. The exonic bases a pair shares so far are summed, and the pair joins as soon
as the sum reaches the bar: the sum only grows, so this is the rule itself. The rule is unchanged (exon unions
overlapping by >= 50% of the smaller, transitively), and so are the sites:
- `test_site_merge_uses_exons_not_spans` pins the two cases that differ between spans and exons;
- `test_site_groups_equal_brute_force_on_random_pileups` compares 60 random pileups against every pair;
- on real off-locus copies the new merge gave the same partition as the span sweep in batches 0-1 (2,166 copies), 5
  (5,594) and 8 (20,581);
- batch 2 alone now takes 225 s at 0.8 GB.

Also seen while checking, and nothing else: the number of sites in those batches (1,048, 816, 3,266 and 564). No recall,
miss, class or C2 number was computed or seen.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

## Results (2026-10-14)

Run: code `3576bb9f` (Amendment 1). A119b's projection was reused (16 batches). The held-out templates were projected
after this file was committed, with the same code and settings (3 batches each): testis in 392 s at 16.2 GB, and KB3781
in 339 s at 19.3 GB (`winloci_data/GGO.splice.mmi`). Scoring took 253 s at 1.3 GB for A119b and under 10 s for each
held-out sample. Guard S0 held: 0.9956, 0.9999 and 1.0000 of templates align to their own exons.

| | A119b (development) | testis (held out) | KB3781 (held out) |
|---|---|---|---|
| templates / alignments / projected copies | 94,685 / 636,800 / 489,897 | 13,048 / 54,244 / 24,807 | 17,672 / 45,933 / 22,740 |
| templates at the 51-alignment cap | 8,682 | 347 | 189 |
| Liftoff extra copies (not short) | 1,888 | 1,888 | 1,129 |
| reachable (source placement assembled) | 420 (0.222) | 62 (0.033) | 52 (0.046) |
| **C1: recall (own)** | **268 / 420 = 0.638, FAIL** | **49 / 62 = 0.790, FAIL** | **44 / 52 = 0.846, pass** |
| recall (any) | 365 / 420 = 0.869 | 58 / 62 = 0.935 | 46 / 52 = 0.885 |
| misses: below threshold / capped / no alignment | 44 / 46 / 62 | 7 / 5 / 1 | 5 / 0 / 3 |
| recall (own) by the source's Liftoff copies: 1 / 2-5 / ≥ 6 | 71/104, 109/157, 88/159 | 9/10, 21/28, 19/24 | 9/12, 3/8, 32/32 |
| sites | 11,463 | 1,861 | 1,600 |
| `liftoff_copy` / `liftoff_gene` / `assembled` / `none` (retro-like) | 375 (1) / 2,643 (47) / 3,504 (62) / 4,941 (266) | 193 (0) / 1,061 (64) / 226 (1) / 381 (43) | 44 (1) / 762 (30) / 212 (1) / 582 (7) |
| **C2: unexplained sites** | **4,675 = 0.408, FAIL** | **338 = 0.182, FAIL** | **575 = 0.359, FAIL** |

The below-threshold misses had a median best identity and coverage of 0.898 and 0.930 in A119b, 0.498 and 0.515 in
testis, and 0.778 and 0.786 in KB3781.

### Verdict (2026-10-14)

**Reading 3: not claimed.** C1 fails in development (0.638 < 0.80). It also fails in testis (0.790) and passes only in
KB3781 (0.846). C2 fails in all three samples. The expected retrocopy advantage did not show: retro-like sites are rare
(266 of A119b's 4,941 `none` sites).

### Post hoc (descriptive, not pre-registered; `bench/template_copies/copies_posthoc.py`)

- **Recall (any) is not all projection.** It also counts a copy that is itself assembled, found by the self alignment of
  the template at the copy: 20 of A119b's 365, none in testis, 1 of KB3781's 46. Found by projection proper: 345, 58,
  45.
- **The misses follow the template's shape.** Call an own template gene-like when it covers ≥ 80% of the source
  placement's exon union and that union covers ≥ 80% of it.
  - Copies whose best own template is gene-like are found 15 of 16 times (A119b), 18 of 18 (testis) and 40 of 40
    (KB3781).
  - Otherwise: 253 of 404, 31 of 44 and 4 of 12.
  - Only 16 of A119b's 420 reachable copies have a gene-like own template. The rule is strict: a Liftoff placement's
    exons are the union of every isoform of the gene, while a template is one transcript.
- **Where the unexplained sites' templates sit.** Each `none` site's best template, placed against Liftoff's gene
  placements:

  | template | A119b (4,941 sites) | testis (381) | KB3781 (582) |
  |---|---|---|---|
  | ≥ 50% inside one gene's exons | 258 | 43 | 157 |
  | shares exonic bases with a gene, less than that | 1,089 | 97 | 175 |
  | inside a gene's span only (intronic) | 996 | 60 | 18 |
  | intergenic | 2,598 | 181 | 232 |

  - Only the intergenic templates (53%, 48% and 40% of the sites) are clearly loci the annotation lacks, which Liftoff
    cannot copy. The rest overlap or sit inside annotated genes. C2's unexplained share is a mix of the two; it is not
    a difference in input alone.
  - The ten templates with the most `none` sites make 443 (9%), 178 (47%) and 256 (44%) of them.
  - In A119b, 2,735 of the 4,941 come from single-exon templates. The contigs holding the most are chrY 492, chr13 459,
    chr1 393, chr9 361, chr15 320 and chr21 286.
  - 36.6% of A119b's templates are single-exon, against 3.4% (testis) and 3.0% (KB3781), so C2 is not comparable across
    samples.
- **The gene-like counts are per sample and concentrated.** KB3781's 40 gene-like copies come from 9 genes, 32 of them
  from two (LOC134756906 19, LOC134757412 13). The two human samples share one Liftoff table, so the same copies can
  appear in both (LINC01784 5 in each). The three samples are not pooled.

### Final-review findings (verified, post hoc, 2026-10-14)

A fresh reviewer re-derived every recall number and stress-tested the site merge against all pairs, forcing the
pruning step, with no failure. The code stands; the text above was corrected. Each count below was recomputed with
`copies_posthoc.py` unless it says "reviewer".
- **Recall (own) includes copies found inside the template's own locus** (as registered, so the verdict is unchanged).
  In A119b, 30 of the 268 are found only by an alignment inside the own template's locus: the reviewer found the
  template itself spans the copy in 17 of them and a secondary hit inside a fused locus in 13. By projection off that
  locus, recall (own) is 238 / 420 = 0.567. Testis and KB3781 have none of these.
- **The held-out pass rests on two genes.** KB3781's 0.846 is per copy, and two genes give 32 of its 52 reachable
  copies, all found. Without the top gene it is 25 / 33 = 0.758; without both, 12 / 20 = 0.60. Testis misses the bar by
  one copy (50 / 62 would pass).
- **Unconditional recall** (found (own) / every Liftoff extra copy): 268 / 1,888 = 0.142, 49 / 1,888 = 0.026 and
  44 / 1,129 = 0.039.
- **The coverage threshold never binds** (reviewer). Identity ≤ coverage by construction (matches ≤ aligned bases, and
  the denominator ≥ template length), and Liftoff applies only `-sc` to extra copies (`liftover_types.py`). No number
  changes.
- **The `no_alignment` label is loose** (reviewer). In A119b, 13 of those 62 misses have a passing own alignment that
  covers 13-49% of the copy, and 4 more have a failing partial one.
- **Two design differences from Liftoff** (reviewer). Liftoff chains several alignment records into one placement; the
  analogue scores single PAF records. A genome-only stretch that minimap2 calls an intron costs nothing here but would
  lower Liftoff's identity.
- **Templates with at least one site** (a pre-registered reported item, omitted above): 13,771, 1,216 and 928.
- **Erratum to Amendment 1.** During profiling, `recall()` ran (1 s) but nothing it returned was printed. "No recall
  ... number was computed or seen" should read "printed or seen".
