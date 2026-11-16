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
