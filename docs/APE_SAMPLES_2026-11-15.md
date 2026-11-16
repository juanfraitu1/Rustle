# Four new ape samples: chimp, bonobo, two orangutans (2026-11-15)

One arc. Registered tests and their results are in `docs/PREREG_ptr_ppa_frozen_rules_2026-11-15.md` (main text +
addenda 1-4; register rows 1747, 1748, 1750, 1751, 1752). Paths are in `docs/DATA.md`. Nothing here changed a default.

## Data
| | PTR chimp | PPA bonobo | PPY Bornean orangutan | PAB Sumatran orangutan |
|---|---|---|---|---|
| assembly (RefSeq NC_ names) | mPanTro3-v2.1, GCF_028858775.2 | mPanPan1-v2.1, GCF_029289425.2 | mPonPyg2-v2.1, GCF_028885625.2 | mPonAbe1-v2.1, GCF_028885655.2 |
| RefSeq gene + pseudogene / CAT transcripts | 41,564 / 157,663 | 38,541 / 156,944 | 37,509 / 155,236 | 41,590 / 155,425 |
| reads (ccs, lima, refine, same lab as OR6737) | 1,509,097 | 1,575,695 | 4,200,083 | 3,112,611 |
| primary mapped / secondary records | 99.97% / 2.44M | 99.99% / 2.13M | 99.45% / 12.6M | 97.86% / 10.1M |
| mean read length / reads with a splice | 2,454 bp / 0.78 | 2,763 / 0.85 | 1,755 / 0.70 | 1,120 / 0.63 |
Fetch: `genomes/fetch_primate.sh TAG ASSEMBLY GENOMEARK_ID` (md5, contig and length checks). Alignment:
`sbatch -J TAG genomes/map_reads.sbatch UNALIGNED.bam`, minimap2 2.30 with the flags of every other sample. Truth GTF:
`bench/famcmp/refseq_gff2gtf.py` (reproduces the gorilla truth file model for model). Tissue: PTR and PPA are testis (user, 2026-11-15). PPY and PAB come from the same library set, but their tissue
has NOT been confirmed: read every "testis" statement below as PTR / PPA plus OR6737 and A119b.

## Default assembly against RefSeq (gffcompare -r, PREFIX.gtf)
| sensitivity / precision (%) | PTR | PPA | PPY | PAB | GGO OR6737 | HSA A119b |
|---|---|---|---|---|---|---|
| intron | 47.4 / 85.6 | 51.7 / 80.5 | 56.7 / 69.9 | 42.6 / 75.9 | 59.7 / 82.6 | 61.2 / 62.7 |
| intron chain | 16.1 / 39.0 | 19.7 / 32.3 | 22.3 / 19.3 | 16.1 / 22.9 | 27.1 / 35.3 | 22.1 / 17.6 |
| locus | 26.9 / 64.3 | 26.8 / 61.3 | 30.9 / 41.6 | 25.6 / 43.1 | 31.2 / 70.5 | 34.9 / 31.4 |
| matching intron chains | 21,049 | 18,143 | 20,246 | 16,116 | 25,873 | 38,727 |
| novel loci | 10.0% | 14.7% | 33.6% | 24.4% | 8.9% | 47.4% |
Sensitivity is not comparable across species (reference sizes differ: PTR 142k mRNAs, PPA 103k). Outputs:
`ptr_ppa_rules/gffcmp/`, `pongo_asm/gffcmp/`.

## Frozen rules on the new samples (RefSeq label)
| | PTR | PPA | PPY | PAB |
|---|---|---|---|---|
| `cont5 >= full` / `cont3 >= full` / retained intron (pooled ratio, all PASS) | 0.10 / 0.22 / 0.32 | 0.24 / 0.37 / 0.37 | 0.09 / 0.20 / 0.35 | 0.08 / 0.22 / 0.42 |
| `trunc5` gain over the read-count curve (bar +0.02) | +0.034 | +0.013 FAIL | +0.027 | +0.043 |
| frozen r1677 logistic score at 0.90 (bar +0.02) | +0.085 | +0.041 | +0.072 | +0.097 |
| `shiftdrop` indel rule, four bars | FAIL one (dropped in RefSeq 0.079) | PASS | PASS | PASS |
The three chain classes hold on eight samples. The frozen score passes on the four new samples (PTR / PPA testis,
PPY / PAB tissue unconfirmed) and failed on the fibroblast (KB3781 +0.006): not shown across tissues. `trunc5` clears its bar on five of eight samples. All four new
samples are spent for every existing rule.

## Why the orangutans score lower (0.19 / 0.23 chain precision)
- Not unannotated loci (6-7% of transcripts), not a thin RefSeq (adding CAT lifts every sample alike), not depth:
  the exact share is lower in every read stratum (21+ reads: 0.57 / 0.54 against 0.78 / 0.72 / 0.74 for PTR / PPA / GGO).
- Truncation classes are two to three times more common (shorter reads), but outside them the 21+ stratum is still
  0.63 / 0.62 against 0.81 / 0.74.
- Decomposition of 21+ read chains against RefSeq or CAT (`pongo_asm/scripts/gap.py`), PPY / PAB vs PTR / PPA: annotated
  chain missing 5' introns 0.087 / 0.142 vs 0.052 / 0.036; one new intron sharing a site 0.081 / 0.070 vs 0.019 /
  0.043 (mostly skipped exons or an exon end > 50 bp away); no annotated intron 0.036 / 0.064 vs 0.008 / 0.016. New
  introns are canonical (non GT-AG / GC-AG <= 0.2%).
- **The unannotated introns are real** (row 1750): of PPY introns lifted to PAB, >= 2 PAB reads splice exactly there
  for 0.989 of annotated-chain introns, 0.914 of the one-new-intron class and 0.898 of the no-annotated-intron class
  (21+ reads; 0.917 / 0.693 / 0.655 at 6-20 reads); only 19-30% of them are annotated in PAB. That part of the gap is
  missing annotation.
- **The 5' shortened part is unexplained** (row 1752): the end-position test was not informative, because annotated
  chains have a concentrated 5' end (10 bp) in only 40-44% of cases.
- Not looked at: where the ~3 secondary alignments per orangutan read land.

## Sister-species reads as truth (rows 1750, 1751)
Maps `genomes/xmap/{PPY_to_PAB,PAB_to_PPY,PTR_to_PPA,PPA_to_PTR}.paf` (minimap2 asm5), `xmap/lift.sh` (paftools
liftover; partial and double lifts dropped). Annotated introns landing exactly on an annotated sister intron:
0.844 / 0.822 / 0.871 / 0.785.
- **Valid for intron structure.** Introns (`xspecies_introns.py`) and chain confirmation used for the 1-read class.
- **Invalid for transcript ends.** Chain-level confirmation (`xspecies_chains.py`: whole chain lifts, >= 1 sister read
  has exactly it) confirms the `cont5` truncation class at 0.77-0.95 of the RefSeq-exact rate: a read cut anywhere in
  the same exon has the same intron chain. Never score `trunc5` or an end rule with it. RefSeq-or-sister precision
  (0.715 / 0.671 / 0.499 / 0.673) is an upper bound.
- **1-read class is noise:** confirmed 0.28-0.35 x as often as the default's 2-read chains on all four samples. The
  2-read floor stays.
- `shiftdrop`'s dropped chains on the orangutans are sister-confirmed 0.30 / 0.50 (all chains 0.47 / 0.67): unresolved
  whether part of the class is real or the indel recurs at the same sequence.

## Open
1. 5' shortened class: needs a negative control (a class known to be cut molecules) and a start measure for broad
   start regions.
2. `shiftdrop` dropped chains that the sister confirms.
3. The orangutan secondary-alignment excess.
4. Any new rule: mine on the four older samples, score once per new sample; the new samples are unseen only for
   rules that do not exist yet.

## Cell-line Iso-Seq of the T2T individuals themselves (2026-11-15 / 16)
Counterparts of gorilla KB3781: reads from the animal whose genome is the reference. Source: BioProject PRJNA1016395 (Yoo 2025 data
statement, vault `Fulltext/Yoo_2025_Complete_Ape_Genomes_fulltext.md`). No registered test; nothing here changed a default.
| | PPY AG05252 | PAB AG06213 | PPA PR00251 | PTR AG18354 |
|---|---|---|---|---|
| cell type (marker reads COL1A1 / COL1A2 vs PTPRC / CD79A) | fibroblast (91,586 / 154,125 vs 8 / 10) | fibroblast (134,257 / 140,906 vs 31 / 14) | fibroblast (59,310 / 102,740 vs 15 / 8) | lymphoblastoid (1,282 / 2,447 vs 4,941 / 6,321) |
| UW runs used (Kinnex, Iso-Seq) | SRR27438214, SRR27438215 | SRR27438208, SRR27438209 | SRR27438210, SRR27438211 | SRR27438206, SRR27438207 |
| reads / primary mapped | 12,543,631 / 99.99% | 11,752,713 / 99.99% | 14,757,139 / 99.98% | 12,847,301 / 99.86% |
| BAM under `reads/TAG/` | `PPY_AG05252_fib_mm.bam` | `PAB_AG06213_fib_mm.bam` | `PPA_PR00251_fib_mm.bam` | `PTR_AG18354_lcl_mm.bam` |
| default assembly: chains / exact vs RefSeq / precision | 77,528 / 23,491 / 0.303 | 70,941 / 24,673 / 0.348 | 97,744 / 20,799 / 0.213 | 89,021 / 30,351 / 0.341 |
| mean read length / introns per read / unspliced reads | 2,163 / 7.3 / 0.18 | 2,413 / 9.1 / 0.08 | 1,881 / 5.5 / 0.39 | 2,307 / 8.9 / 0.11 |
| chains that are a truncated piece of a RefSeq chain | 0.18 | 0.12 | 0.25 | 0.12 |
Method (workspace, outside the repo): `genomes/fetch_fibro.sh TAG STEM RUN...` (ENA FASTQ, md5 per run, joined into one file),
`sbatch -J TAG genomes/map_reads.sbatch X.fastq.gz STEM_mm` (project flags; takes FASTQ and an output stem, refuses to overwrite),
`sbatch -J TAG cellline_asm/job.sbatch STEM` (default `assemble`, writes the `as_table` sidecar, `curve.py` vs RefSeq). The human and gorilla
BAMs are linked into `reads/HSA`, `reads/GGO` (`GGO_KB3781_fib_mm.bam` = `GCA_029281585.2_flnc_mm.bam`).
- **Run sets.** KB3781's BAM holds four runs (UW SRR27438212 / 13, PSU SRR27178662 / 63). The UW runs are oriented FLNC. The PSU runs
  (`mXxx_IsoSeq`, `_IsoSeq_long`) are deposited as raw CCS with primers and polyA in mixed orientation: `-uf` needs them refined first
  (KB3781's were, by the group that built that BAM). Only the UW runs were used for the four new samples.
- **Orientation check.** minimap2's `ts` tag does not show read orientation; count GT-AG vs CT-AC intron motifs relative to the read on a
  sample aligned without `-uf` (oriented: ~2,000 sense vs ~40 antisense; unoriented: about equal).
- **Bonobo is fragmentary, in the RNA.** Lower precision at every read support (0.52 at >= 21 reads vs 0.66-0.70), 14% of all reads
  unspliced, intronic and on the gene's sense strand (others 1-5%), sense : antisense 3.6 : 1, so RNA fragments and not DNA. RefSeq size and
  the share of chains with a novel intron (0.46-0.50) are the same across the four. The PSU bonobo runs (SRR27178652 / 51; FLNC by
  `cutadapt --revcomp --discard-untrimmed --poly-a -g "GCAATGAAGTCGCAGGGTTGGG...GTACTCTGCGTTGATACCACTGCTT"`, 2,644,957 of 3,424,967 kept,
  the rest 5'-primer-on-both-ends artifacts) show the same pattern (unspliced 0.32, 13% intronic sense). Merged with UW
  (`PPA_PR00251_fib_all_mm.bam`, 17.40M reads): 107,007 chains, 22,084 exact, 0.206: +1,285 correct for +9,263 chains. Keep the UW-only
  BAM as the bonobo sample. PSU alone scores 0.330 at 2.64M reads, a depth effect. Not measured: absorbing truncated chains.
- **Not run on these samples:** families, assign, any rule.

## Other tools on every ape sample, de novo and guided (2026-11-16, register 1763, report only)

User request: StringTie, FLAIR, IsoQuant and IsoCon on the ape samples that lacked them, guided and unguided. Workspace `apetools/`
(= `r1763_apetools`), one directory per sample; `job.sbatch` (`sbatch -J TAG job.sbatch STEM [REGION]`) runs every arm on one node and
ends with the tables below. Recipes are those of `benchmark_collapse` (unguided), r1622 and r1646 (guided): StringTie 3.0.3 `-L [-G]`,
FLAIR 3.0.1 `bam2bed` + `collapse` (guided: `correct` + `collapse -f`), IsoQuant 4.1.0 `-d pacbio_ccs [--genedb]`, ours
`assemble --guide` with today's defaults. One guide per sample for every guided arm: RefSeq with 5/16 of the transcripts hidden
(`guide_proto.py split`, `GUIDE_HIDE=tx`). All nine jobs ended with every arm non-empty.

De novo, exact intron chains against the full RefSeq, correct chains / precision. Raw operating points, not matched counts.

| Sample | ours.final.gtf | ours.gtf | StringTie | FLAIR | IsoQuant |
|---|---|---|---|---|---|
| PTR testis | 20,993 / 0.462 | 21,904 / 0.384 | 21,532 / 0.350 | 19,730 / 0.273 | 15,761 / 0.441 |
| PPA testis | 17,501 / 0.359 | 18,545 / 0.315 | 17,288 / 0.281 | 15,934 / 0.229 | 13,109 / 0.369 |
| PPY testis | 20,639 / 0.239 | 21,347 / 0.189 | 20,561 / 0.167 | 19,881 / 0.109 | 14,094 / 0.245 |
| PAB testis | 16,436 / 0.294 | 16,930 / 0.225 | 16,981 / 0.234 | 16,513 / 0.149 | 12,866 / 0.309 |
| PPY fib AG05252 | 23,405 / 0.368 | 24,366 / 0.302 | 20,787 / 0.314 | 23,246 / 0.101 | 14,803 / 0.351 |
| PAB fib AG06213 | 23,661 / 0.398 | 25,253 / 0.346 | 21,510 / 0.388 | 23,993 / 0.131 | 15,667 / 0.426 |
| PPA fib PR00251 | 21,592 / 0.280 | 22,151 / 0.217 | 18,780 / 0.195 | 21,171 / 0.063 | 13,139 / 0.217 |
| PTR lcl AG18354 | 29,200 / 0.382 | 31,155 / 0.337 | 25,433 / 0.366 | 30,658 / 0.130 | 16,579 / 0.404 |
| GGO fib KB3781 | 24,383 / 0.397 | 26,315 / 0.354 | 22,087 / 0.380 | 24,734 / 0.125 | 15,785 / 0.427 |

`ours.gtf` has the most correct chains on every sample but PAB testis (StringTie +51). IsoQuant is more precise than `ours.final.gtf`
on five of nine samples at 55-75% of our correct chains. FLAIR de novo emits very large sets on the deep samples (337,026 chains on
PPA fib).

Guided, hidden chains found / precision of the non-guide chains (`guide_proto.py score`; whole-set guided precision is inflated by the
guide and is not the comparison).

| Sample | ours guided | StringTie -G | FLAIR | IsoQuant |
|---|---|---|---|---|
| PTR testis | 3,979 / 0.446 | 3,554 / 0.295 | 4,473 / 0.180 | 4,223 / 0.179 |
| PPA testis | 3,516 / 0.318 | 3,231 / 0.225 | 3,198 / 0.155 | 3,486 / 0.139 |
| PPY testis | 4,047 / 0.174 | 3,543 / 0.145 | 3,464 / 0.061 | 3,898 / 0.085 |
| PAB testis | 3,111 / 0.191 | 2,887 / 0.165 | 2,667 / 0.090 | 3,644 / 0.116 |
| PPY fib | 5,432 / 0.290 | 4,136 / 0.258 | 4,403 / 0.076 | 4,010 / 0.133 |
| PAB fib | 5,715 / 0.302 | 4,296 / 0.306 | 4,571 / 0.102 | 4,297 / 0.176 |
| PPA fib | 4,527 / 0.210 | 3,577 / 0.152 | 3,847 / 0.050 | 3,535 / 0.071 |
| PTR lcl | 7,031 / 0.313 | 4,658 / 0.285 | 7,329 / 0.101 | 4,392 / 0.158 |
| GGO KB3781 | 6,206 / 0.301 | 4,555 / 0.304 | 4,874 / 0.095 | 4,334 / 0.177 |

Ours finds the most hidden chains on six of nine samples (FLAIR on PTR testis and PTR lcl, IsoQuant on PAB testis); StringTie finds the
most guide chains everywhere. Our de novo assembly finds more hidden chains than our guided one (e.g. 5,787 vs 3,516 on PPA testis):
the guided mode's 3-read floor on non-guide chains.

Reference-free arm. IsoCon has no guided mode and runs on read pools of a target panel (row 1478), so it was run only where a panel
exists: KB3781 with the OR6737 gorilla panel (`apetools/GGO_KB3781_isocon`, `isocon_kb.sbatch`; the panel was chosen for OR6737 testis
expression, 53 of 75 annotated copies are expressed here). Copies modelled / junction precision / spliced transcripts: ours 50 / 0.74 /
233, StringTie 49 / 0.70 / 287, FLAIR 47 / 0.63 / 365, IsoCon 44 / 0.79 / 663, IsoQuant 42 / 0.71 / 149. For the four testis samples
the row 1759 replacement `isoseq cluster2` ran on the same FLNC reads (`apetools/TAG_cluster2/cluster2.gtf`, `cluster2.sbatch`):
PTR 62,932 chains / 17,417 correct / 0.277; PPA 59,301 / 15,339 / 0.259; PPY 144,851 / 17,737 / 0.122; PAB 89,523 / 13,773 / 0.154.
The cell lines have no reference-free arm (SRA FASTQ, `cluster2` untried on it); IsoCon on chimp, bonobo and orangutan needs a
per-species panel (`isocon_targets.py select`).

Not done: matched-count curves on the nine samples above. The A119b and OR6737 tool tables were redone with today's defaults on
2026-11-16 (next section); their StringTie arm is still 3.0.1. Guided IsoQuant spends about 2.5 h converting the guide GTF (single thread).

Follow-up, register 1764 (`docs/PREREG_isocon_junctions_2026-11-16.md`): IsoCon's higher junction precision on the gorilla panel is its
operating point. With our junctions cut by read support to IsoCon's count per copy, ours is ahead on both samples (KB3781 0.926 vs
0.893; OR6737 0.894 vs 0.878, RefSeq intron truth); alignment slips are 2-4 junctions per sample in either tool. Nothing to port.

## A119b and OR6737 with today's defaults, IsoQuant added (2026-11-16, register 1765)

Rerun of the two base samples after the transfrag-extension and chain-score defaults, for the assembler artifact
(claude.ai/artifact/BKZMLTF4aqfVPD3yYjcYfd, rebuilt by `build.py human.json` in workspace `artifact_refresh_2026-11-16`; runs
`{g,h}.{def,pri}`, scorer `score.py` = `curve.load` over any number of GTFs, which reproduces the 10 Nov figures on the old files).
`--seed-pool primary`, exact intron chains against the full RefSeq: distinct multi-exon chains / correct / precision. StringTie 3.0.1
and FLAIR 3.0.1 are the `benchmark_collapse` files; IsoQuant 4.1.0 is row 1646's `-d pacbio_ccs --polya_trimmed all` arm.

| Sample | ours.gtf | ours.final.gtf | StringTie | FLAIR | IsoQuant |
|---|---|---|---|---|---|
| A119b | 224,494 / 39,823 / 0.177 | 175,447 / 38,543 / 0.220 | 211,865 / 31,781 / 0.150 | 626,983 / 39,230 / 0.063 | 84,539 / 18,903 / 0.224 |
| OR6737 | 75,018 / 26,263 / 0.350 | 60,442 / 24,930 / 0.413 | 67,491 / 23,248 / 0.345 | 138,125 / 24,293 / 0.176 | 37,089 / 16,677 / 0.450 |

Default pool (tied secondary alignments in): `ours.gtf` 236,301 / 40,086 / 0.170 and 76,139 / 26,424 / 0.347; `ours.final.gtf`
182,771 / 38,873 / 0.213 and 61,337 / 25,049 / 0.408. Against the 10 Nov runs the primary-pool `ours.gtf` gains 1,385 correct chains
on A119b (0.185 -> 0.177) and 629 on OR6737 (0.355 -> 0.350). The chain-score cut was fixed on these two samples, so their
`final.gtf` figures are in-sample.

IsoQuant's precision is its operating point (register 1765, `docs/PREREG_isoquant_junctions_2026-11-16.md`,
`bench/denovo_eval/isoquant_junc.py`, workspace `r1765_isoquant_junc`), at both levels:
- Chains: ours cut by reads to IsoQuant's correct count, 0.540 vs 0.450 (OR6737, >= 6 reads) and 0.414 vs 0.224 (A119b, >= 12
  reads); `final.gtf` 0.567 and 0.504.
- Junctions, RefSeq genes with a spliced model in both tools (gffcompare tmap), ours cut by read support to IsoQuant's count per
  gene, RefSeq intron truth: 0.925 vs 0.910 (OR6737, 14,037 genes) and 0.818 vs 0.731 (A119b, 20,444 genes).
- IsoQuant keeps almost no slips (8 / 14 against our 192 / 311: a wrong junction within 10 bp of an annotated one at <= 0.2 of its
  reads), under 2% of our wrong junctions at the cut. It keeps junctions that are neither GT-AG nor GC-AG / AT-AC (803 / 5,751;
  ours 0). Nothing to port.

Not done: the loci counts of the artifact's tabs 2 and 3 (434 / 1,928 new loci) were not recounted on the new runs; our slips were
not checked against `shiftdrop`; no IsoCon panel for A119b; KB3781 is outside row 1765 (no primary-pool run, other IsoQuant arm).
A `tools/` snapshot made for a job needs `bench` linked beside it, or the chain score step warns and leaves `final.gtf` unfiltered.

