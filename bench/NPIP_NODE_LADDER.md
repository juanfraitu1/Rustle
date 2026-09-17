# NPIP node ladder: annotated DNA nodes → de novo shared-definition nodes (human CHM13, 2026-09-17)

This is a measurement on a development family (NPIP, 27 truth copies). It is not validation. Nothing is committed and no pipeline code was changed. Declarations were written before any metric: `/mnt/linuxdisk/home/juanfraitu/npip_ladder/DECLARATIONS.txt`, first written 00:28:21, md5 3249e522. Addendum 1 was appended at 00:58:12, before any arm was scored (md5 ab17e692). Everything above the addendum is unchanged.

## Read first

Terms: **F** is the bipartite F with NPIP-only precision. **F strict** counts every node in the matched family, including nodes that are not truth copies. Scoring is FAMILY level with triangle leaders.

| arm | nodes | families | R | P NPIP-only | **F** | P strict | F strict | components F / F strict | matched family = 27 copies + NPIP-exon fragments + other loci | NPIP pairs with a direct edge: exon-only / body-only / both (of 351) | copies with a full-length node (≥ 0.8 of copy exon bp) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A0 ANN | 240 | 33 | 1.000 | 1.000 | **1.000** | 0.600 | 0.750 | 1.000 / 0.415 | 45 = 27 + 7 + 11 | 7 / 161 / 168 | 27/27 |
| A1 ANN-EXPR | 167 | 23 | 1.000 | 1.000 | **1.000** | 0.675 | 0.806 | 1.000 / 0.540 | 40 = 27 + 7 + 6 | 7 / 161 / 168 | 27/27 |
| A2 ANN-ALLTX | 240 | 33 | 1.000 | 1.000 | **1.000** | 0.600 | 0.750 | 1.000 / 0.415 | 45 = 27 + 7 + 11 | 9 / 141 / 188 | 27/27 |
| A3a READ-EXTENT | 158 | 21 | 1.000 | 1.000 | **1.000** | 0.692 | 0.818 | 1.000 / 0.524 | 39 = 27 + 8 + 4 | 7 / 161 / 168 | 26/27 |
| A3b READ-REP | 167 | 20 | 1.000 | 1.000 | **1.000** | 0.659 | 0.794 | 1.000 / 0.524 | 41 = 27 + 7 + 7 | 19 / 143 / 186 | 27/27 |
| A4 DN-SD | 560 | 68 | 1.000 | 1.000 | **1.000** | 0.290 | 0.450 | 1.000 / 0.196 | 93 = 27 + 23 + 43 | 9 / 68 / 187 | 5/27 |

Drop in F from one step to the next (a positive number is a drop):

| step | drop in F | drop in F strict |
|---|---|---|
| A0→A1 (expression) | +0.000 | −0.056 (rises) |
| A1→A3a (read-derived exons/body) | +0.000 | −0.012 (rises) |
| A1→A3b (read-derived rep transcript) | +0.000 | +0.012 |
| **A1→A4 (full de novo nodes)** | +0.000 | **+0.356** |
| A0→A2 (transcript set) | +0.000 | +0.000 |

**Decision-rule outcome:** F is 1.000 in all six arms, so all four declared steps drop by exactly 0.000. That is a four-way tie, and the pre-registered rule does not pick a first fix target on this substrate. The only step where F strict (reported alongside) drops substantially is A1→A4 (0.806 → 0.450).

**A2 vs A0:** using the full transcript set changes nothing at family level. The triangle families and the components are identical as node sets. Only the edges move: 31 more exon edges (7 more pairs, 953 → 960). Among NPIP pairs, 20 change from body-only to both, and 2 pairs with no direct edge in A0 gain an exon-only edge.

What the numbers show (every item is a count from the tables):
- **Edges and grouping are not the NPIP problem on this substrate.** Every arm puts all 27 copies into one family, under both triangle leaders and connected components. So no copy is cut off from the family and there is no family split. Not every copy pair has a direct edge, though. Pairs with one: 336/351 in A0, A1 and A3a; 338 in A2; 348 in A3b; 264 in A4. The other pairs are joined through shared neighbours.
- **The de novo gap is node quality.** Only 5/27 copies get a full-length node in A4, against 27/27 in A1. A4 maps 4 copies to an opposite-strand node (NPIPB2, LOC128966608, NPIPB6, NPIPB12). The A4 family also carries 23 non-mapped nodes that overlap NPIP copy exons (fragments, versus 7 in A1). It carries 43 other loci (versus 6 in A1):
  - 23 have no RefSeq exon but overlap an NPIP copy's gene span: 21 lie fully inside it and 2 only partly; 16 are antisense and 7 on the same strand;
  - 11 are intergenic;
  - 3 sit inside another gene's span;
  - 6 overlap other RefSeq exons (for example SMG1P6, BOLA2-SMG1P6, LOC124907845, PKD1P5-LOC105376752).
- The annotated arms also carry non-NPIP neighbours (PKD1P1/2/3, PKD1P readthroughs, CLN3, MIR6770/MIR6511 in A0). Those neighbours are why P strict is 0.600–0.692 even with annotated nodes.
- Read-derived extents (A3a) and read-derived representative transcripts (A3b) do not reproduce the A4 loss. A3a loses 1 full-length copy; A3b adds exon edges (19 exon-only NPIP pairs, versus 7). Both keep F strict at 0.79–0.82. The loss comes with the de novo node set itself: fragmentation, split loci, and extra unannotated read loci.

## Per-copy, A1 (ANN-EXPR) and A4 (DN-SD) side by side

Columns: node (0-based half-open, strand, exon count) · covered = M1 node exon bp over the copy's exons / copy exon bp · #ov = number of nodes whose exons overlap the copy's exons · fam = predicted family · in = sits in the matched FAMILY group · why = reason when not counted. Reads = primary MAPQ ≥ 1 same-strand reads on the copy's exons.

| copy | Iso-Seq group | reads | A1 node | A1 covered | A1 #ov | A1 fam | A1 in | A1 why | A4 node | A4 covered | A4 #ov | A4 fam | A4 in | A4 why |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NPIPB2 | B2 | 422 | chr16:11963310-12012756 - 12ex | 1702/1702 (1.00) | 1 | GWFAM0 | yes | - | chr16:11963703-11968998 + 2ex | 416/1702 (0.24) | 2 | GWFAM0 | yes | - |
| NPIPA2 | A2/3 | 350 | chr16:14740933-14763869 + 16ex | 4517/4517 (1.00) | 1 | GWFAM0 | yes | - | chr16:14749588-14763868 + 7ex | 2827/4517 (0.63) | 3 | GWFAM0 | yes | - |
| NPIPA1 | A1 | 741 | chr16:14938509-14953106 + 8ex | 1084/1084 (1.00) | 2 | GWFAM0 | yes | - | chr16:14938814-14953123 + 7ex | 1021/1084 (0.94) | 1 | GWFAM0 | yes | - |
| PKD1P6-NPIPP1 | A4 | 456 | chr16:15105352-15141806 - 30ex | 5391/5391 (1.00) | 2 | GWFAM0 | yes | - | chr16:15115422-15132313 - 10ex | 1898/5391 (0.35) | 5 | GWFAM0 | yes | - |
| NPIPA5 | A5 | 146 | chr16:15368419-15386444 - 11ex | 3741/3741 (1.00) | 1 | GWFAM0 | yes | - | chr16:15368407-15382704 - 7ex | 2915/3741 (0.78) | 2 | GWFAM0 | yes | - |
| NPIPA6 | A6-9 | 127 | chr16:16340280-16359022 + 10ex | 1570/1570 (1.00) | 3 | GWFAM0 | yes | - | chr16:16344720-16359014 + 8ex | 1088/1570 (0.69) | 2 | GWFAM0 | yes | - |
| NPIPA7 | A6-9 | 81 | chr16:16391363-16406190 + 8ex | 1264/1264 (1.00) | 1 | GWFAM0 | yes | - | chr16:16391852-16406187 + 7ex | 1093/1264 (0.86) | 1 | GWFAM0 | yes | - |
| NPIPA8 | A6-9 | 38 | chr16:18325161-18343979 - 10ex | 1588/1588 (1.00) | 2 | GWFAM0 | yes | - | chr16:18325177-18339582 - 7ex | 1080/1588 (0.68) | 2 | GWFAM0 | yes | - |
| NPIPA9 | A6-9 | 849 | chr16:18372371-18391116 - 10ex | 2869/2869 (1.00) | 2 | GWFAM0 | yes | - | chr16:18372354-18376104 - 2ex | 704/2869 (0.25) | 1 | GWFAM0 | yes | - |
| NPIPB3 | B3-5 | 1298 | chr16:21337399-21360419 - 8ex | 3674/3674 (1.00) | 1 | GWFAM0 | yes | - | chr16:21337422-21344055 - 3ex | 2989/3674 (0.81) | 2 | GWFAM0 | yes | - |
| LOC128966608 | B3-5 | 485 | chr16:21680685-21703416 - 7ex | 3557/3557 (1.00) | 2 | GWFAM0 | yes | - | chr16:21692009-21693757 + 1ex | 129/3557 (0.04) | 2 | GWFAM0 | yes | - |
| NPIPB4 | B3-5 | 648 | chr16:22359821-22382362 + 12ex | 3634/3634 (1.00) | 2 | GWFAM0 | yes | - | chr16:22382178-22387825 + 2ex | 184/3634 (0.05) | 4 | GWFAM0 | yes | - |
| NPIPB5 | B3-5 | 280 | chr16:22781606-22814310 + 9ex | 9237/9237 (1.00) | 1 | GWFAM0 | yes | - | chr16:22788631-22792526 + 1ex | 3031/9237 (0.33) | 4 | GWFAM0 | yes | - |
| NPIPB6 | B6-9 | 659 | chr16:28623032-28645151 - 14ex | 7025/7025 (1.00) | 1 | GWFAM0 | yes | - | chr16:28625753-28628690 + 1ex | 2635/7025 (0.38) | 5 | GWFAM0 | yes | - |
| NPIPB7 | B6-9 | 124 | chr16:28737008-28778611 - 11ex | 2952/2952 (1.00) | 2 | GWFAM0 | yes | - | chr16:28736995-28771999 - 8ex | 1274/2952 (0.43) | 2 | GWFAM0 | yes | - |
| NPIPB8 | B6-9 | 44 | chr16:28928824-28939457 + 7ex | 1293/1293 (1.00) | 1 | GWFAM0 | yes | - | chr16:28935124-28939424 + 3ex | 781/1293 (0.60) | 1 | GWFAM0 | yes | - |
| NPIPB9 | B6-9 | 614 | chr16:29016080-29053452 + 13ex | 4193/4193 (1.00) | 2 | GWFAM0 | yes | - | chr16:29048645-29050801 + 1ex | 520/4193 (0.12) | 3 | GWFAM0 | yes | - |
| NPIPB10P | B10 | 66 | chr16:29318967-29333452 + 9ex | 880/880 (1.00) | 1 | GWFAM0 | yes | - | chr16:29328121-29333923 + 1ex | 441/880 (0.50) | 1 | GWFAM0 | yes | - |
| NPIPB11 | B11 | 155 | chr16:29663350-29688504 - 8ex | 6117/6117 (1.00) | 1 | GWFAM0 | yes | - | chr16:29663364-29679826 - 7ex | 3236/6117 (0.53) | 2 | GWFAM0 | yes | - |
| NPIPB12 | B12/13 | 27 | chr16:29765317-29788442 - 10ex | 3737/3737 (1.00) | 2 | GWFAM0 | yes | - | chr16:29765341-29767376 + 1ex | 1783/3737 (0.48) | 3 | GWFAM0 | yes | - |
| LOC124907834 | B12/13 | 487 | chr16:30507433-30529793 - 7ex | 3244/3244 (1.00) | 1 | GWFAM0 | yes | - | chr16:30507437-30523997 - 7ex | 3120/3244 (0.96) | 1 | GWFAM0 | yes | - |
| NPIPB13 | B12/13 | 151 | chr16:30609396-30634670 - 9ex | 6144/6144 (1.00) | 2 | GWFAM0 | yes | - | chr16:30609406-30625966 - 7ex | 3108/6144 (0.51) | 2 | GWFAM0 | yes | - |
| NPIPB14P | B14 | 1242 | chr16:75785805-75805631 - 9ex | 2418/2418 (1.00) | 2 | GWFAM0 | yes | - | chr16:75785725-75790348 - 2ex | 1574/2418 (0.65) | 2 | GWFAM0 | yes | - |
| NPIPB15 | B15 | 168 | chr16:80194083-80209885 + 9ex | 4340/4340 (1.00) | 2 | GWFAM0 | yes | - | chr16:80195301-80209898 + 7ex | 1605/4340 (0.37) | 1 | GWFAM0 | yes | - |
| LOC124907808 | B15 | 22 | chr16:80308276-80324260 + 10ex | 4981/4981 (1.00) | 1 | GWFAM0 | yes | - | chr16:80319739-80324246 + 2ex | 2777/4981 (0.56) | 2 | GWFAM0 | yes | - |
| LOC124907807 | B15 | 77 | chr16:80421481-80438591 + 9ex | 4571/4571 (1.00) | 1 | GWFAM0 | yes | - | chr16:80424045-80438604 + 7ex | 1805/4571 (0.39) | 1 | GWFAM0 | yes | - |
| NPIPB1P | B1 | 407 | chr18:11781940-11796470 - 8ex | 981/981 (1.00) | 1 | GWFAM0 | yes | - | chr18:11781592-11796212 - 7ex | 810/981 (0.83) | 1 | GWFAM0 | yes | - |

All 27 copies have ≥ 22 reads, so A1 drops no copy. No copy fits any "why" category in A1 or A4: every copy is in the matched family.

## SUBFAMILY levels

No arm splits NPIP, so the matched group is always the single NPIP family. Only the strict columns change between arms.

| grouping | level | arm | R | P NPIP-only | F | P strict | F strict | pair sens | pair prec strict | matching |
|---|---|---|---|---|---|---|---|---|---|---|
| triangle | SUBFAMILY-1 | A0 | 0.704 | 0.704 | 0.704 | 0.422 | 0.528 | 1.000 | 0.201 | NPIPB→GWFAM0 (J 0.704) |
| triangle | SUBFAMILY-1 | A1 | 0.704 | 0.704 | 0.704 | 0.475 | 0.567 | 1.000 | 0.255 | NPIPB→GWFAM0 (J 0.704) |
| triangle | SUBFAMILY-1 | A2 | 0.704 | 0.704 | 0.704 | 0.422 | 0.528 | 1.000 | 0.201 | NPIPB→GWFAM0 (J 0.704) |
| triangle | SUBFAMILY-1 | A3a | 0.704 | 0.704 | 0.704 | 0.487 | 0.576 | 1.000 | 0.269 | NPIPB→GWFAM0 (J 0.704) |
| triangle | SUBFAMILY-1 | A3b | 0.704 | 0.704 | 0.704 | 0.463 | 0.559 | 1.000 | 0.243 | NPIPB→GWFAM0 (J 0.704) |
| triangle | SUBFAMILY-1 | A4 | 0.704 | 0.704 | 0.704 | 0.204 | 0.317 | 1.000 | 0.047 | NPIPB→GWFAM0 (J 0.704) |
| triangle | SUBFAMILY-2 | A0 | 0.148 | 0.148 | 0.148 | 0.089 | 0.111 | 1.000 | 0.024 | A6-9→GWFAM0 (J 0.148) |
| triangle | SUBFAMILY-2 | A1 | 0.148 | 0.148 | 0.148 | 0.100 | 0.119 | 1.000 | 0.031 | A6-9→GWFAM0 (J 0.148) |
| triangle | SUBFAMILY-2 | A2 | 0.148 | 0.148 | 0.148 | 0.089 | 0.111 | 1.000 | 0.024 | A6-9→GWFAM0 (J 0.148) |
| triangle | SUBFAMILY-2 | A3a | 0.148 | 0.148 | 0.148 | 0.103 | 0.121 | 1.000 | 0.032 | A6-9→GWFAM0 (J 0.148) |
| triangle | SUBFAMILY-2 | A3b | 0.148 | 0.148 | 0.148 | 0.098 | 0.118 | 1.000 | 0.029 | A6-9→GWFAM0 (J 0.148) |
| triangle | SUBFAMILY-2 | A4 | 0.148 | 0.148 | 0.148 | 0.043 | 0.067 | 1.000 | 0.006 | A6-9→GWFAM0 (J 0.148) |
| components | FAMILY | A0 / A1 / A2 / A3a / A3b / A4 | 1.000 (all) | 1.000 (all) | 1.000 (all) | 0.262 / 0.370 / 0.262 / 0.355 / 0.355 / 0.108 | 0.415 / 0.540 / 0.415 / 0.524 / 0.524 / 0.196 | 1.000 (all) | 0.067 / 0.134 / 0.067 / 0.123 / 0.123 / 0.011 | NPIP→GWFAM0 |
| components | SUBFAMILY-1 | A0 / A1 / A2 / A3a / A3b / A4 | 0.704 (all) | 0.704 (all) | 0.704 (all) | 0.184 / 0.260 / 0.184 / 0.250 / 0.250 / 0.076 | 0.292 / 0.380 / 0.292 / 0.369 / 0.369 / 0.138 | 1.000 (all) | 0.038 / 0.076 / 0.038 / 0.070 / 0.070 / 0.006 | NPIPB→GWFAM0 |
| components | SUBFAMILY-2 | A0 / A1 / A2 / A3a / A3b / A4 | 0.148 (all) | 0.148 (all) | 0.148 (all) | 0.039 / 0.055 / 0.039 / 0.053 / 0.053 / 0.016 | 0.062 / 0.080 / 0.062 / 0.078 / 0.078 / 0.029 | 1.000 (all) | 0.005 / 0.009 / 0.005 / 0.008 / 0.008 / 0.001 | A6-9→GWFAM0 |

FAMILY pairwise with triangle leaders: sens 1.000 in every arm; prec strict 0.355 / 0.450 / 0.355 / 0.474 / 0.428 / 0.082 (A0 / A1 / A2 / A3a / A3b / A4).

Edge totals (exon edges / gene-body edges / distinct pairs / loci in triangle families): A0 571/860/953/177 · A1 461/711/762/119 · A2 602/860/960/177 · A3a 530/736/784/112 · A3b 564/711/819/118 · A4 1349/1885/1919/364.

## Rust vs mirror (A4)

| check | Rust (`RUSTLE_SHARED_DEFINITION=1`) | mirror | agree |
|---|---|---|---|
| reps → gene-level loci + read-locus nodes = nodes | 534 → 522 + 38 = 560 (stderr) | 534 → 522 + 38 = 560 | yes |
| node set: tx query keys (chrom\|strand\|rep_exons) and sequence md5 | 560 in captured `tx.fa` | 560 | 560/560 keys, 0 md5 mismatches |
| node set: body query keys (chrom\|start\|end) and sequence md5 | 560 in captured `body.fa` | 560 | 560/560, 0 mismatches |
| edges (exon / body / pairs) | 1349 / 1885 / 1919 (stderr) | 1349 / 1885 / 1919 on Rust's PAF; identical on the union PAF | yes |
| families (sets of chrom, start, end, strand, exons) | 68 families, 364 loci | 68 on Rust's PAF, 68 on the union PAF | **exact, both** (symmetric difference 0) |

- **minimap2 differs slightly between runs.** Rust's internal run and the union run give the same record sets for 558/560 tx queries and 560/560 body queries, once MAPQ and tags other than `cg` are ignored. Two tx queries differ by 4 secondary records in total, which changes no edge. This is Addendum 1.
- **The Rust run used a wrapper as `RUSTLE_MINIMAP2`.** The wrapper swapped `build()`'s `target.fa` for prebuilt `.mmi` indexes of the same sorted-contig genome. Indexing inline would have exceeded the 10-minute call cap: the splice index alone took 403 s and 19.2 GB.

Differences between Rust `shared_definition.rs` and the prototype `bench/denovo_shared_def.py` (Rust wins; none of them affects this run):
1. `transcript_hits` and gene-body `emit`: Rust guards `bl > 0` / `qlen > 0`; the prototype divides unguarded.
2. Target: Rust maps against every contig of `--fasta`, in sorted-name order (overridable with `RUSTLE_SD_TARGET_ORDER`). The prototype maps against `SHAREDEF_CONTIGS` only (`ggo3.fa`) and restricts nodes and reads to those contigs.
3. Reads: Rust scans the whole BAM in file order. The prototype uses `pysam.fetch` per contig, in coordinate order. The order changes group order only, not membership.
4. MAPQ 255: noodles returns None and Rust treats it as 0, so the read fails MAPQ ≥ 1. pysam keeps 255, so the prototype passes it.
5. The prototype has `--min-mapq` and `--alpha` (a relative link floor in `split_linked`). Rust hard-codes MAPQ ≥ 1 and alpha 0.
6. The prototype's `families` default grouping is `components`, with `leaders`, `bridges` and `triangle` as options. Rust ships only `triangle_leaders`.
7. The prototype names queries by md5[:16] of the key and aligns in 4 Mb batches against a `.mmi`. Rust uses raw key strings in a single run.
8. `ExonIndex.hits`: the prototype uses a bisect window of maxlen, Rust a prefix scan. The results are identical.

## Readings chosen (full text in DECLARATIONS.txt)

- **Regions:** the 19 lit windows and 15 copyregions windows, merged into 18 intervals (bookended intervals join). All 27 truth copies' exons lie fully inside them (S3 check: 0 exons outside).
- **Reads:** primary, MAPQ ≥ 1 (255 counts as 0), strand from `ts:A` as in `read_blocks()`. A node's reads are same-strand reads with a block overlapping a node exon by ≥ 1 bp. In total, 52,917 of 504,090 records.
- **A0 nodes:** one node per gene or pseudogene feature with ≥ 1 exon descendant overlapping the regions. Exons are the union of all descendants, clipped to the gene span.
  - Truth copies use the frozen exons. NPIPB14P gets the PDXDC2P-NPIPB14P readthrough exons clipped to its span.
  - rep = MANE Select, else RefSeq Select, else most exonic bp; ties go to the first in the file.
  - Node row order puts the 27 copies first, then coordinate order. This settles the exact M1 tie between NPIPB14P and the readthrough, and it is the last tie-break key for leaders.
- **A1:** A0 nodes with n_reads ≥ 3.
- **A2:** every candidate transcript is its own splice query.
- **A3a:** depth-2 read exons, intersected with the gene span; nodes under 100 bp are dropped (9 dropped, all miRNA or LOC loci; no truth copy).
- **A3b:** the most-supported exact intron chain.
  - Its model runs from the lower median of the first-block starts to the lower median of the last-block ends.
  - Ties go to more exonic bp, then leftmost.
  - This is a stand-in for the pipeline's representative transcript.
- **A4:** mirror nodes and the union PAF are scored. Rust families are identical, so no second scoring was needed.
- **Edge procedure:** as `shared_definition.rs`. The target is all 25 CHM13 contigs in sorted-name order. One minimap2 run per preset on the union of all arms' queries, deduplicated by sequence md5:
  - tx: 1,477 sequences, 4.0 Mb, 2 chunks;
  - body: 862 sequences, 12.9 Mb, 3 chunks.
- **Components:** connected components of the same edge pairs, size ≥ 2.
- **Scorer:** the frozen M1 loop, copied. Before any arm was scored, it reproduced the frozen `group_level.tsv` rows for DN0_lit and DN1rB_lit exactly (assert).
- **"why" categories:** unexpressed, missing, fragment (< 0.8), merged, split. None applied here.
- **NPIP pair edges:** pairs of M1 nodes that are distinct and directly joined. No pair shares a node.

## Provenance

- Repo `/mnt/c/Users/jfris/Desktop/Rustle`, branch `dna-from-genome`, HEAD `111f6727` (src clean). `gw_family_catalog` was rebuilt with `CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin gw_family_catalog` (exit 0, 1 m 18 s). It is copied to `npip_ladder/rust/gw_family_catalog.bin`, sha256 `9d4482840ef65f6978260001dbd129d06cfb532ec1d57f274c84dd1dc0ea154b`. minimap2 2.30-r1287.
- Work directory `/mnt/linuxdisk/home/juanfraitu/npip_ladder/`:
  - `regions.bed`, `reads.bam` (+ `.bai`);
  - `idx/target.{fa,splice.mmi,asm20.mmi}`;
  - `rust/{rust_default.*,dump/,rust_sd.*,mm2_wrapper.sh,capture/}`;
  - `union/{tx,body}_union.{fa,paf}` with chunks, `build.pkl`, `reads.pkl`;
  - `results.json`, `report_tables.md`;
  - `scripts/{ladder.py,make_report.py}`;
  - `logs/`.
- Inputs: `winloci_data/A119b.t2t.bam`, `winloci_data/chm13v2.0.fa`, `o1_falsemerge/lit/refseq_c16_17_18.gff`, `o1_falsemerge/lit/windows.bed`, `bakeoff/human/copyregions.bed`, `npip_per_member/truth_exons.tsv`, `docs/lit_subclusters_npip_dishuck_check.tsv`, and the frozen scorer `npip_per_member/wholechr16/scripts/npip_per_member.py`.

| step | command | wall time |
|---|---|---|
| reads | `samtools view -b -M -L regions.bed -@ 4 A119b.t2t.bam` + `samtools index` | 22 s |
| target | `samtools faidx -n 0 chm13v2.0.fa <25 contigs, LC_ALL=C sorted>` | 58 s |
| splice index | `minimap2 -x splice -t 4 -d target.splice.mmi target.fa` | 403 s (19.2 GB) |
| asm20 index | `minimap2 -x asm20 -t 4 -d target.asm20.mmi target.fa` | 196 s (11.9 GB) |
| Rust rep dump | `RUSTLE_ER_EDGE_DUMP=dump/x gw_family_catalog --bam reads.bam --fasta chm13v2.0.fa --out rust_default --homology-primary --threads 4` | 3 m 35 s |
| Rust shared definition | `RUSTLE_SHARED_DEFINITION=1 RUSTLE_MINIMAP2=mm2_wrapper.sh TMPDIR=npip_ladder/tmp gw_family_catalog … --out rust_sd …` | 6 m 00 s (minimap2 105 s + 49 s) |
| build arms | `ladder.py verify-scorer` ; `ladder.py build` | 3 s ; 11 s |
| union mapping | `minimap2 -c -N 50 -p 0.1 -x splice -uf -t 4 target.splice.mmi tx_chunk{0,1}.fa` ; `minimap2 -c -N 50 -p 0.1 -x asm20 -t 4 target.asm20.mmi body_chunk{0,1,2}.fa` | 55 + 62 s ; 48 + 54 + 44 s |
| score | `ladder.py score` ; `make_report.py` | 3 s |

## Caveats

- **Development family:** NPIP is where the definition was developed. These numbers are not validation, and F = 1.000 here says nothing about held-out families.
- **One haplotype and one library:** CHM13 reference, a single testis Iso-Seq library (A119b). All 27 copies have ≥ 22 reads here, so A1 cannot test unexpressed copies.
- **A3b is a stand-in** for the pipeline's representative-transcript rule, not the rule itself.
- **This substrate hides the loss at NPIP-only level.** The windows are copy-centred, so every locus in them is near an NPIP copy. A saturated F (1.000 everywhere) can only show losses through P strict and node quality. The 43 non-NPIP A4 loci show up only in the strict columns.
- **Not comparable with earlier catalogs:** the frozen DN0 and DN1r+B results (F 0.773 and 0.826 on the lit windows) used the E_r de novo catalog and a substrate missing 5 copies. They are not a rung of this ladder.
- **Strand in mapping:** M1 ignores strand. 4 A4 copies map to opposite-strand nodes and still count.

## Verification (independent recompute)

An independent verifier ran this check on 2026-09-17. It did not read or import the builder's `scripts/`. Its own code is in `npip_ladder/verify/`:
- `v_nodes.py` rebuilds the nodes from the primary inputs.
- `v_edges.py` rebuilds edges and families from `shared_definition.rs`, read directly.
- `v_score.py` rescores every arm.
- `v_paf_diff.py` and `v_samtools.py` run the tool-level checks.

Nothing was committed.

**1. Declarations before metrics.**
- The part of `DECLARATIONS.txt` above ADDENDUM 1 still has the original md5, `3249e522…` (00:28:21), so nothing above the addendum changed.
- The full file (md5 `ab17e692…`) was last modified at 00:58:12. That is before `results.json` and `logs/score.log` (00:58:30) and `report_tables.md` (00:59:35).
- The union PAFs (00:57:39) and `build.pkl` (00:52:43) are not metrics.
- Addendum 1 is disclosed above.
- `ladder.py` has an mtime of 00:58:26, later than `build.pkl`. It does not matter, because the rebuilt node tables below match `build.pkl` exactly.

**2. Nodes, edges and families: every arm matches exactly.**
- *Nodes.* Every arm was rebuilt from the GFF, `reads.bam` and the Rust dump, not from `build.pkl`. All six arms match the builder's tables row for row: order, chrom, strand, exons, rep_exons and n_reads, plus the A2 transcript sets. The counts are A0 240, A1 167, A2 240, A3a 158, A3b 167 and A4 560 (534 reps → 522 + 38). There are 52,917 reads out of 504,090 records. The same 9 A3a nodes are dropped, and the regions are rebuilt to the same 18 intervals.
- *Union FASTA.* Every union FASTA name equals the md5 of its sequence. Every arm's query sequences, taken from the genome, are present in it, with none missing and none extra.
- *Edges and families.* `transcript_hits`, `gene_body_chains`, `edges`, `triangle_leaders` and connected components were rewritten from the Rust source. On the union PAFs they give exactly the builder's counts:

  | arm | exon edges | body edges | pairs | triangle families (loci) | components (loci) |
  |---|---|---|---|---|---|
  | A0 | 571 | 860 | 953 | 33 (177) | 21 (184) |
  | A1 | 461 | 711 | 762 | 23 (119) | 14 (121) |
  | A2 | 602 | 860 | 960 | 33 (177) | 21 (184) |
  | A3a | 530 | 736 | 784 | 21 (112) | 12 (116) |
  | A3b | 564 | 711 | 819 | 20 (118) | 13 (122) |
  | A4 | 1349 | 1885 | 1919 | 68 (364) | 31 (394) |

  The triangle and component families match `results.json` in membership, id and order for all six arms.

**3. Scoring: every number matches.**
- The copied scorer reproduces the frozen `group_level.tsv` FAMILY, SUBFAMILY-1 and SUBFAMILY-2 rows for DN0_lit and DN1rB_lit.
- Rescored from the verifier's own families, all 36 `core` rows (6 arms × 2 groupings × 3 levels) are character-identical to `results.json`.
- Every number in the tables above agrees: the Read-first table, drops, SUBFAMILY and components rows, pairwise precision and edge totals.
- The per-copy table (27 rows, A1 and A4) matches cell for cell, including read counts.
- So do full-length counts, NPIP pair kinds (no pair shares a node), family compositions, the A4 breakdown (23 / 11 / 3 / 6 with the same example genes), the 4 opposite-strand copies, the A0 neighbour names, A2 = A0 families, and the caveat numbers F 0.773 / 0.826 with 5 missing copies.
- Three sentences were imprecise and are corrected in place:
  1. "no missing edge between NPIP copies" was false. 3 to 87 copy pairs have no direct edge, depending on the arm.
  2. A2 vs A0 also gives 2 NPIP pairs a first, exon-only, edge.
  3. Not all 23 unannotated A4 loci lie inside an NPIP copy span: 21 are fully inside and 2 only partly.

**4. Arm faithfulness.**
- *A0/A2.* 12 genes were checked against the raw GFF with awk, independently of the Python parser: NPIPB14P, PDXDC2P-NPIPB14P, NPIPA7, NPIPA8, NPIPB2 (6 transcripts), NPIPB5 (4 transcripts), NPIPA1, CLN3 (the 253 bp `gene-CLN3` record; GeneID 1201 has a second record), PKD1P1, SMG1P6, MIR6770-3 (primary transcript plus 2 miRNAs) and LOC124907845. In every case, node exons = the union of exon descendants clipped to the span, and rep = the RefSeq Select transcript, else the one with most exonic bp. NPIPB14P has no transcript of its own, so its rep is its node exons, taken from the readthrough clipped to its span. The GFF has no MANE Select tags.
- *A1.* Reads for all 27 copies were recounted from `samtools view -F 2308 -q 1 -M -L <copy exons>` text output. The count keeps only reads on the copy's strand (flag 0x10 flipped when `ts:A:-`) whose CIGAR blocks overlap the copy's exons. All 27 equal the report, the minimum is 22, and every copy has ≥ 3 reads. The plain `samtools view -c` span counts, which ignore strand, are higher (for example NPIPB12 has 245 vs 27), as expected. The BAM has 0 primary MAPQ-255 records.
- *A3a.* Depth-2 exons were recomputed with `samtools depth -J` on each copy's same-strand reads, for NPIPB2 (11 exons), NPIPA7 (8) and NPIPB14P (5). They equal the builder's A3a exons exactly.
- *A3b.* The chain choice was recomputed from samtools text for NPIPA1, NPIPB3 and NPIPB14P, and the chosen model equals the builder's. For all three, the most-supported "chain" is the unspliced one (74, 398 and 39 reads), so the representative is a single block. **Across the ladder, 13/27 copy nodes and 59/167 A3b nodes get a single-exon representative.** A3b is therefore a weak stand-in for a spliced representative transcript. Its "does not reproduce the A4 loss" result should be read with that in mind.
- *A4.* The verifier's own A4 nodes match Rust's captured `tx.fa` and `body.fa`: 560/560 keys each and 0 sequence mismatches. On Rust's captured PAFs they give 1349 / 1885 / 1919 edges, with an edge-pair symmetric difference of 0 against the union PAF. The 68 families (364 loci) equal Rust's `rust_sd.copies.tsv` as sets, both on the union PAF and on Rust's PAFs. Rust's `GWFAM` ids follow a different order from the mirror's (−size, first index) order. GWFAM0, the 93-node NPIP family, has the same id in both, but only 4 of the 68 ids coincide. The sets are the same, and no reported number depends on the ids.
- *Addendum 1.* When records are compared per query key (MAPQ and all tags except `cg` ignored), Rust's run and the union run differ on exactly 2 tx queries, by 4 records in total, and on 0 body queries. This confirms Addendum 1.
- *Prototype differences.* The 8 listed Rust-vs-prototype differences match `bench/denovo_shared_def.py` and `bench/guided_pipeline.py`: unguarded `nm/bl`, `SHAREDEF_CONTIGS` + `fetch` per contig, `--min-mapq`/`--alpha`, `components` as the default in `families --arms`, md5[:16] names, `BATCH_BP = 4_000_000`, and a bisect window of `maxlen`. With no MAPQ-255 reads and no `=`/`X` in the PAF CIGARs, none of them affects this run.

**5. Substrate and provenance.**
- `regions.bed` equals the merge of `windows.bed` (19) and `copyregions.bed` (15), which gives 18 intervals. All 27 copies' exons lie fully inside them.
- The BAM `@PG` line is `samtools view -b -M -L regions.bed -@ 4 -o reads.bam …/A119b.t2t.bam` (samtools 1.22.1), which matches the table apart from `-o`.
- `idx/target.fa` holds the 25 CHM13 contigs in sorted-name order, with lengths equal to the genome `.fai`.
- The binary's sha256 is `9d448284…`, and minimap2 is 2.30-r1287.
- Index times and memory (403 s / 19.2 GB, 196 s / 11.9 GB) and mapping times (Rust 105 + 49 s; union 55 + 62 s, 48 + 54 + 44 s) match the minimap2 logs.

**Verdict:** the numbers reproduce. Three sentences were corrected, and the A3b single-exon issue is added. The decision-rule outcome stands as reported: F is 1.000 in every arm, so the four steps tie and the rule picks no target; F strict drops only at A1→A4. NPIP is the development family, so none of this is validation.
