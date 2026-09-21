# Node construction: the de novo defect is OVER-merge, not fragmentation

**§6u6, 2026-09-21.** Substrate: A119b IsoSeq vs CHM13, chr16 de novo (`--assemble-only`),
2,550 de novo loci, 2,080 annotated chr16 genes. Pre-registration:
`docs/PREREG_read_bridged_node_merge_2026-09-21.md` (md5 `10c85736`).
⚠ HUMAN substrate — do not pool with the gorilla numbers.

## Headline

The standing picture was that de novo node construction loses to fragmentation — *"median 2 de novo loci
per NPIP gene"*, and the §6p1 5 kb same-strand merge exists to consolidate. **On a fixed-universe,
one-to-one, collision-penalised scorer that is measurably sensitive to fragmentation, fragmentation
costs nothing and over-merge costs everything.**

| per-copy correctness by how many de novo loci overlap the gene | genes | correct | rate |
|---|---|---|---|
| 1 locus | 627 | 269 | **0.429** |
| 2 | 337 | 194 | 0.576 |
| 3 | 157 | 102 | 0.650 |
| 4 | 89 | 57 | 0.640 |
| >= 5 | 119 | 83 | **0.697** |

The rate RISES with fragmentation. The instrument is not blind to the axis — it points the other way.
Baseline failure decomposition over the 1,329-gene universe: **collision 490 = 36.9%**, node spills
outside its gene 134 = 10.1%, correct 705 = 53.0%.

## The pre-registered rule is REFUTED

RBM: merge two same-strand loci when >= N primary reads (`-F 2308`, MAPQ >= Q) each cover >= 25 bp of
exonic sequence unique to each side. Bar was **>= +2.0pp** per-copy correctness on chr16.

| arm | merges | nodes | correct | rate | delta | collisions | false-merge rate | junk |
|---|---|---|---|---|---|---|---|---|
| baseline | 0 | 2550 | 705 | 0.5305 | — | 214 | 26.6% | 405 |
| Q>=60 N>=1 | 1191 | 1359 | 359 | 0.2701 | **-26.03pp** | 283 | 33.6% | 248 |
| Q>=60 N>=3 | 753 | 1797 | 500 | 0.3762 | -15.43pp | 267 | 28.8% | 334 |
| Q>=60 N>=10 | 437 | 2113 | 602 | 0.4530 | -7.75pp | 251 | 26.9% | 381 |
| Q>=60 N>=25 | 256 | 2294 | 653 | 0.4913 | -3.91pp | 233 | 26.7% | 398 |

Every arm negative, monotone in merge aggressiveness. The **post-hoc** variant restricted to
exon-overlapping pairs (not in the prereg) is also negative on every arm — it only approaches baseline
as it stops merging: N>=3 **-1.96pp**, N>=10 -1.05pp, N>=25 -0.60pp.

chr19 (the held-out substrate) was **not run**: the rule failed its development bar, so there is nothing
to confirm. It remains unexamined and available for the next rule.

## Why it fails — the bridges are readthrough, not splits

| Q>=60 bridged pairs | pairs | same home gene | exon-overlapping | median start-distance |
|---|---|---|---|---|
| N>=1 | 2,033 | 1,003 = **49.3%** | 27.3% | 30,283 bp |
| N>=10 | 641 | 427 = 66.6% | 44.8% | 17,655 bp |
| N>=50 | 217 | 153 = 70.5% | 47.9% | 13,806 bp |

At N>=1 more than half the bridged pairs join **different genes** a median 30 kb apart. Read count does
not separate the populations — 29.5% still cross gene boundaries at N>=50. This is register 690's
*"read-through molecules (>= 3) chain neighbours"*, measured at scale.

## What the over-merge actually is

Of 369 colliding annotated-gene pairs: **172 = 46.6% are genuinely disjoint genes fused into one node**
(median intergenic gap 6,110 bp), 180 = 48.8% have overlapping spans with disjoint exons, 17 = 4.6% are
true annotation nesting. And the fusions are not thin:

| reads linking the two disjoint genes inside one node | median | q25 | q75 | <= 2 | >= 3 |
|---|---|---|---|---|---|
| any MAPQ | 46 | 10 | 98 | 17 | 155 |
| MAPQ >= 60 | 38 | 7 | 89 | 26 | 146 |
| spliced | 16 | 9 | 82 | 17 | 155 |

**A >= 3 MAPQ-60-read support floor leaves 84.9% of the fusions standing** — the fusions are
heavily-supported readthrough transcription, so no support threshold cuts them.

Nor does the curated label rescue them: only **4.7%** (6 spanned by a curated readthrough gene + 2 with a
readthrough-labelled member) of the 172 are annotated as readthrough at all; **83.1% have no annotated
record spanning both**. Examples: `ERI2`+`THUMPD1` (38,375 bp gap), `ERI2`+`LOC124903661` (57,342 bp),
`ADAT1`+`TMEM231` (40,623 bp).

## Two construction invariants established on the way

1. **Junction-disjointness between de novo loci is BY CONSTRUCTION** — 0 of 13,685 chr16 junctions are
   used by more than one locus. Any "these two loci share no junction" statistic is a tautology, not
   evidence. This retires a whole shape of proposal.
2. **The split is mid-exon, not at a boundary.** Of 728 intra-gene same-strand spliced-spliced
   exon-overlapping pairs, 685 = 94.1% share NO exact exon — they partially overlap one. And 53.7% of
   adjacent same-gene locus pairs overlap in span (median gap **-1,575 bp**), so the §6p1 5 kb DISTANCE
   merge was never addressing the dominant class.

## Consequences

- The §6p1 5 kb same-strand merge should stay off, and register 879's *"should not be turned on to chase
  this number"* now has a mechanism: consolidation attacks a defect that does not cost the metric and
  amplifies one that does.
- Register 487 (orphan locus-stitching, empty class) re-measured and re-confirmed: **12 pairs = 0.9%**.
- The open lever is a node SPLIT at read-supported long-range chains, which needs its own
  pre-registration and must answer register 846 ([[project_node_cut_rule]]: cutting a chimera at its
  parent boundary DOUBLES rather than separates, and short pieces become hubs).

## Reproduce

```
python3 bench/read_bridged_merge.py --gtf dn16.gtf --bam chr16.bam \
    --gff chr16.genes.gff --chrom chr16 --mapq-tiers 60 --sweep 1,3,5,10,25
```
⚠ The max-overlap tie-break must be deterministic (sorted by node id); dict iteration order moved 5
genes between runs and made the headline irreproducible (705 vs 710).
