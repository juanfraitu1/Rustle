# Figure 3 — Transcripts assembled in multi-mapping loci

**Building loci also from secondary alignments within 2% of the best alignment score rebuilds annotated transcripts
that no primary alignment reaches.**

> **Status (2026-09-25).** The build is now genome-wide for both samples (every contig the sample's annotation
> covers). Gorilla already was; the gorilla numbers below stay. The human numbers below are from **chr20–22** and are
> replaced by the genome-wide rebuild (`make.py data fig3`: one pass over the 96 GB human BAM, 20–40 min in bounded
> calls); `python3 figures/fig_secondary.py summary` then prints every number this caption quotes.

> **Mode (apples to apples).** Annotation-free (de novo) on every side: Rustle assemble (reads + genome), StringTie
> `-L` without `-G`, FLAIR collapse without annotation (`flair correct` skipped), IsoSeq collapse (checked in
> `benchmark_collapse/run_stringtie.sbatch` and `run_flair.sbatch`; the figure's first line says the same). **Guided
> comparison: not available (guided StringTie/FLAIR GTFs not supplied).** When the user registers annotation-guided
> runs (`samples.tsv` columns `stringtie_guided_gtf`, `flair_guided_gtf`), the same build scores them into a separate
> table and figure (`fig3_guided_bins`, drawn as fig3g_guided; benchmark samples only, where the reads per transcript are counted), never in a panel with the annotation-free methods. Rustle has no annotation-guided
> transcript assembly, so that figure has no Rustle row, and the guided tools are there scored against the annotation
> they were given (docs/archive/2026-09/PREREG_guided_transcript_comparison_2026-09-25.md; GLOSSARY *Modes*).

**Claim.** The reference is each genome's RefSeq annotation, and a match is an exact intron chain (gffcompare `=`).
In gorilla OR6737 testis IsoSeq, genome-wide (88,387 multi-exon transcripts with at least one read), Rustle builds
loci from the primary alignments plus the secondary alignments that score within 2% of the read's best alignment
score (AS) anywhere in the genome. With them it rebuilds 50 of the 839 transcripts that no primary alignment reaches
(6.0%; 95% interval 2.6–10.9%, resampling read-sharing groups or genes), in 23 independent units (21 read-sharing
groups and 2 genes). No other method matches any of the 839: Rustle with primary alignments only by construction,
and StringTie, FLAIR and IsoSeq collapse, run on the same alignments, match none either. Without NC_073244.2 (gorilla
chr20), the contig the seeding default was decided on (and Fig. 6d's contig), the result reads 33 of 794 (4.2%) in 21
units, and again no other method matches any.

Where at least one primary alignment reaches the transcript, the methods are not separated. In the top tie-fraction
bin (more than 90% of the reads have a second alignment within 2% of their best; 1,262 transcripts, 235 read-sharing
groups), Rustle reproduces 291 chains in 79 groups, IsoSeq collapse 248 in 83, Rustle with primary alignments only
201 in 71, StringTie 194 in 73 and FLAIR 162 in 59. Rustle reproduces the most chains, IsoSeq collapse reaches the
most groups, and their 95% intervals overlap (16.8–30.2% vs 14.1–25.7%). Set against any baseline (StringTie, FLAIR or
IsoSeq collapse), the one-sided matches are 56 Rustle-only (31 groups) vs 62 baseline-only (36 groups) in the top
bin and 12 (10) vs 57 (33) in the (0.5, 0.9] bin; only the no-primary category is one-sided (50 in 23 units vs 0).

Against Rustle with primary alignments only the gain is paired and one-sided: 160 high-tie transcripts gained in 66
read-sharing groups on 23 contigs, 3 lost. The largest group, the CGB-like tandem array LOC129528600–LOC129528625 on
NC_073244.2, gives 20 gains rebuilt from the same 4 reads; without it the top bin reads 287 chains (78 groups) for
Rustle against 244 (82) for IsoSeq collapse. Without NC_073244.2 the gains are 137 in 63 groups, and the top bin reads
280 for Rustle against 237 for IsoSeq collapse (n = 1,227). Where no read has a second alignment within 2% (tie
fraction 0) the two Rustle configurations match 24,408 and 24,407 chains. Adding the secondary alignments lowers
intron-chain precision (Fig. 1). In gorilla, the 1,074 extra multi-exon transcripts add 162 matching chains (15.1%);
in human chr20–22 (current table), the 1,017 extra add 4 (0.4%). Human chr20–22 is too small to separate the
methods (152 top-bin transcripts in 40 read-sharing groups), and no human method, Rustle included, matches any of
the 69 transcripts that no primary alignment reaches there.

## Panels

Definitions used in all panels:
- *Candidate alignment.* A read's primary alignment, or a secondary alignment whose AS is at least 98% of the read's
  best AS anywhere in the genome. This is exactly the set Rustle builds loci from. Supplementary and unmapped
  records never count.
- *Counting.* A read counts at a transcript when one of its candidate alignments has aligned bases (a gapless
  M/=/X run) that overlap one of the transcript's exons. Intron `N` and deletion `D` spans do not count, so a read
  that splices over the transcript does not count.
- *Tie fraction.* The share of a transcript's reads that have a second alignment somewhere in the genome scoring at
  least 98% of their best AS (from the whole-BAM table of each read's best and second-best AS). A read with a single
  alignment never has one.
- *Reached by a primary alignment.* At least one of the transcript's reads is counted there through its **primary**
  alignment.
- *Read-sharing group.* Transcripts with tie fraction > 0.5 that share any read with a second alignment within 2% are
  one group (union-find over the whole genome). "Groups" in the panels means read-sharing groups.
- *Independent unit.* A transcript's unit is its read-sharing group when it has one, and otherwise its gene. Panels
  a–d all count these units.

**a** (gorilla OR6737, testis; mGorGor1 GCF_029281585.2; genome-wide) **and b** (human A119b, T2T-CHM13 v2.0; current
table chr20–22, genome-wide after the rebuild). The fraction of each category's multi-exon RefSeq transcripts with at
least one read whose intron chain each method reproduces exactly. Transcripts with at least one read: gorilla 88,387
of 95,833; human 9,711 of 9,846. Both panels share one x range.

*Facets.* The first five facets are tie-fraction bins among transcripts reached by at least one primary alignment:

| tie fraction | gorilla n (units) | human n (units), chr20–22 |
|---|---|---|
| 0 | 81,893 (22,763 genes) | 9,024 (2,028 genes) |
| (0, 0.1] | 1,978 (489 genes) | 223 (72 genes) |
| (0.1, 0.5] | 1,450 (470 genes) | 165 (31 genes) |
| (0.5, 0.9] | 965 (150 groups) | 78 (14 groups) |
| > 0.9 | 1,262 (235 groups) | 152 (40 groups) |

The last facet ("all") holds every transcript that no primary alignment reaches, whatever its tie fraction:
- gorilla: 839 transcripts in 218 units (read-sharing groups, or genes where a transcript has none). 780 have tie
  fraction 1, 14 lie in (0.5, 1) and 45 are at ≤ 0.5.
- human (chr20–22): 69 transcripts in 5 read-sharing groups, all with tie fraction 1.

Each transcript is in exactly one facet.

*Methods* are the rows:
- Rustle (filled circle): the pipeline default.
- Rustle, primary alignments only (open circle): the same assembler run with `--no-seed-secondaries`.
- StringTie 3.0.1, FLAIR 3.0.1 and IsoSeq collapse 26.2.0, run by the lab on the same alignments.

*Intervals.* Bars are 95% intervals from resampling independent units: a percentile bootstrap with 2,000 resamples,
the same resamples for every method. Where no transcript or every transcript is matched, the interval is Wilson's
with n = the number of units. In the facets with tie fraction > 0.5, the printed values are transcripts matched
(units with a match).

*No-primary facet.* Only Rustle with primary alignments only, 0 by construction, has no interval. Every other 0
carries its Wilson interval: gorilla 0–1.7% over 218 units (drawn, but hidden under the markers), human 0–43% over 5
read-sharing groups (visible for Rustle, StringTie, FLAIR and IsoSeq collapse).

*Other bins.* In gorilla, IsoSeq collapse reproduces the most chains in each of the three bins up to 0.5. In the
tie-fraction-0 bin: IsoSeq collapse 27,938, Rustle 24,408, FLAIR 23,128 and StringTie 22,169 of 81,893. Fig. 1 covers
this difference.

**c** (gorilla) **and d** (human). The categories of a and b with tie fraction > 0.5. Each bar counts transcripts
matched by one side only:
- Right of zero: matched by Rustle but by none of the three baselines (upper bar of each pair, "Any baseline"), or
  matched by Rustle but not by Rustle with primary alignments only (lower bar).
- Left of zero: the reverse. For "Any baseline": matched by at least one of StringTie, FLAIR and IsoSeq collapse,
  and not by Rustle.
- "N in G groups": the read-sharing groups behind the count; in the no-primary block, "N in G units" (read-sharing
  groups, or genes where a transcript has none), the same 218 units as panel a.

Titles: gorilla 3,021 transcripts with tie fraction > 0.5 in 373 read-sharing groups; human (chr20–22) 299 in 45.

| category | vs any baseline: Rustle only / baseline only | vs Rustle, primary alignments only: Rustle only / other only |
|---|---|---|
| gorilla (0.5, 0.9] | 12 in 10 groups / 57 in 33 | 19 in 14 / 0 |
| gorilla > 0.9 | 56 in 31 / 62 in 36 | 93 in 42 / 3 in 3 |
| gorilla, no primary alignment | 50 in 23 units / 0 | 50 in 23 / 0 |
| human (0.5, 0.9], chr20–22 | 1 in 1 / 3 in 2 | 1 in 1 / 0 |
| human > 0.9, chr20–22 | 2 in 2 / 2 in 1 | 3 in 3 / 0 |
| human, no primary alignment, chr20–22 | 0 / 0 | 0 / 0 |

In gorilla, the 23 units behind the 50 no-primary gains are 21 read-sharing groups and 2 transcripts with tie
fraction ≤ 0.5.

Below the panels is the precision cost of the added secondary alignments. It is read at plot time from
`fig1_gffcompare.tsv`: gffcompare intron-chain precision (matching chains over multi-exon transcripts), Rustle vs
Rustle with primary alignments only, and what the extra multi-exon transcripts add.
- Gorilla genome-wide: 35.3% vs 35.6% (25,872 of 73,345 vs 25,710 of 72,271); the 1,074 extra transcripts add 162
  matching chains (15.1%).
- Human chr20–22 (current table): 17.3% vs 18.7% (2,263 of 13,093 vs 2,259 of 12,076); the 1,017 extra add 4 (0.4%).

**e** One example locus, chosen by a fixed rule (`pick_example`). Among the transcripts no primary alignment reaches
and Rustle rebuilds exactly, it takes those with tie fraction > 0.5 (the typical case: 794 of the 839) whose
annotated introns are all at least 20 bp, then the one with the most reads. Ties go to the shortest span, then the
id. The intron condition drops RefSeq models with a 1–2 bp "intron", an indel correction of a model "modified
relative to the genomic sequence", which shows no splicing; no candidate's shortest intron lies between 3 and 81 bp,
so the cut does not choose among real introns.

The rule picks gorilla LOC115933275 / XM_063702354.1 (NC_073224.2:223,067,331–223,078,481, minus strand; RefSeq
product "neuroblastoma breakpoint family member 3, transcript variant X2"), 9 exons.
- It has 140 reads, all with tie fraction 1. None has its primary alignment on its exons; every one has another
  alignment within 2% of its best AS, and the aligner put the primary elsewhere. Its read-sharing group holds 18
  transcripts (tie fraction > 0.5) of 11 genes that share 418 reads.
- The primary-alignment depth track is 0 across the window.
- Rustle builds four transcripts there (three drawn). One matches the example exactly; one of the faded ones has the
  intron chain of the gene's other RefSeq transcript XM_055385705.2, which is also a no-primary gain. No other method
  has a transcript in the window.

Tracks:
- Top two tracks: read depth per 50 bp, counting each alignment once per bin in which it has aligned bases. The
  first track counts primary alignments, which every method receives. The second counts the secondary alignments
  with AS ≥ 98% of the read's best ("Secondary ≥ 98% of best"), which only Rustle's default configuration uses.
- Below: the gene's RefSeq transcripts (the example is dark and starred), then each method's transcripts that
  overlap the window. At most three are drawn per method, exact matches first. `=` marks an exact match to the
  example's intron chain; other chains are faded.

## Methods essentials

- **Reads and alignments.**
  - Gorilla: OR6737 testis IsoSeq on mGorGor1 (`GGO_mm.bam`); reference annotation RefSeq GCF_029281585.2
    (`GGO.GCF_029281585.2_RefSeq.gtf.gz`).
  - Human: A119b IsoSeq on T2T-CHM13 v2.0 (`A119b.t2t.bam`); reference annotation CHM13 v2.0 RefSeq
    (`A119b.chm13v2.0_RefSeq.gtf.gz`). Genome-wide = the 24 annotated contigs (the annotation has no chrM record);
    current table chr20–22.
  - Both BAMs: minimap2 2.31 `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`.
- **Methods.**
  - Rustle: `tools/rustle_pipeline.sh assemble`, run genome-wide, default configuration, and again with
    `--no-seed-secondaries`.
  - StringTie 3.0.1, FLAIR 3.0.1 and IsoSeq collapse 26.2.0: the lab's GTFs.
  - Every method and the annotation are restricted to the same contigs.
- **Matching.** Each method is scored with gffcompare 0.12.10 against the restricted annotation. A reference
  transcript is matched when a transcript of the method has class `=` to it or to a reference transcript with the
  identical intron chain. gffcompare reports one `ref_id` per transcript, so matches are propagated to annotation
  duplicates, for every method alike.
- **Counting.**
  - `as_table` scans the full BAM once (genome-wide for both samples) and writes each read's best and second-best AS
    (the pipeline's `assemble` stage writes it; `${work}/runs/<sample>/<sample>.molecules.tsv`).
  - `python3 figures/fig_secondary.py count` (pysam) streams the evaluation contigs one at a time (cached per contig;
    `--budget-s` bounds one call): gorilla 10,649,751 mapped non-supplementary records on every annotated contig,
    human (current table) 4,282,433 on chr20–22.
  - It discards secondary alignments below 98% of the read's best AS: gorilla 6,156,058, human (chr20–22) 2,810,890.
  - It keeps multi-exon reference transcripts with at least one read. No method matches any transcript without a
    read (gorilla 7,446, human chr20–22 135).
- **Commands.**
  - `python3 figures/make.py data fig3` builds the four tables (heavy).
  - `python3 figures/make.py plot fig3` renders the figure.
  - `python3 figures/fig_secondary.py summary` prints every number in this caption from the tables, including the
    lines without NC_073244.2.
- **Bin edges.** The edges are fixed in `TIE_BINS` and were set before any per-method match rate was examined.
  Gorilla transcripts with at least one read: 81,926 (92.7%) at exactly 0, 3,440 in (0, 0.5], 1,369 in (0.5, 1) and
  1,652 at exactly 1.
  - Bin 0 isolates the spike at 0.
  - 0.5 separates "most reads have one best place" from "most reads have a second alignment within 2%".
  - > 0.9 captures the point mass at 1.
  - The primary-support split was added afterwards and is not tuned: it is the definition of "reachable from
    primary alignments".

| Table | Produced by |
|---|---|
| `fig3_ref_tie` | `fig_secondary.build` → `species_rows`: per-transcript counts, tie fraction, read-sharing group, matches |
| `fig3_bins`, `fig3_gain` | `fig_secondary.summarize` over `fig3_ref_tie` (categories, bootstrap over units, one-sided matches vs any baseline and vs Rustle with primary alignments only) |
| `fig3_example` | `fig_secondary.pick_example` + `example_rows` (region fetch of the BAM and the method GTFs; gorilla) |
| note under c–d | `fig_secondary.precision_cost`, read at plot time from `fig1_gffcompare.tsv` |

## Caveats

- **Transcripts are not independent.** Near-identical copies are rebuilt from the same reads, so quote read-sharing
  groups alongside transcripts. The CGB-like array's 40 transcripts with tie fraction > 0.5 share 4 reads in total.
- **"By construction" applies only to Rustle with primary alignments only.** StringTie, FLAIR and IsoSeq collapse
  received the same alignments, secondary alignments included. StringTie reads secondary alignments (verified in
  `bench/COPY_ASSIGNMENT_AND_GATE.md`), so their 0 of 839 is observed, and panel a gives it its Wilson interval
  (0–1.7% over 218 units).
- **Fig. 6d's seeding contrast is not independent of this figure.** Its whole sensitivity gain (the RFPL4A-like
  array, protein-homology family PF82) and 210 of its 275 precision pairs (the CGB-like array) come from the two
  NC_073244.2 arrays whose transcripts are 23 of this figure's 160 seeded-only gains, in 3 of its 66 read-sharing
  groups. Figure 6a–c (human chr16) shares nothing with this figure.
  - The CGB-like array, LOC129528600–LOC129528625 (20 genes; NC_073244.2:60,735,364–60,864,845), gives 20 gains here
    in one read-sharing group. In Fig. 6d its 20 genes, with LOC101154318 upstream, form the 21-locus cluster behind
    the 210 precision pairs.
  - The RFPL4A-like array ("ret finger protein-like 4A") gives 3 gains here in 2 read-sharing groups: LOC129528712 and
    LOC129528713 in the group spanning NC_073244.2:67,785,174–67,816,085 (IsoSeq collapse also matches both), and
    LOC101151087 in the one spanning 67,821,739–67,845,225. In Fig. 6d, LOC129528713 is the one gene of the scored set that
    the added secondary alignments bring into the PF82 cluster.
- **An exact match at every copy is not a claim about which copy is expressed.** Rustle places the transcript of
  reads with several near-equal alignments at each copy the aligner reports. Choosing between copies is the
  copy-assignment stage.
- **Union asymmetry in c–d.** "Any baseline" is the union of three methods set against Rustle alone, which favours
  the left side. In the tie-fraction-0 bin that union matches transcripts Rustle does not: 9,363 vs 216 the other way
  in gorilla, and 968 vs 26 in human (chr20–22).
- **Human chr20–22 is uninformative here.** Every human 95% interval in the facets with tie fraction > 0.5 overlaps
  every other method's. The gain from the added secondary alignments is 4 transcripts in 4 read-sharing groups. The
  genome-wide rebuild replaces this caveat.
- **Read pool.** Secondary alignments below 98% of the read's best AS are candidate alignments for no method, so
  they do not count as aligned to a transcript. The counts describe the pool Rustle actually builds loci from.
- **`=` compares intron chains only**, so transcript ends are free.
