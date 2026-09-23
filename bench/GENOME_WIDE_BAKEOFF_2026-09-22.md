# Genome-wide `--assemble-only` bakeoff: ours vs StringTie vs FLAIR

**§6w4, 2026-09-22.** User goal: *"run the pipeline as assembly only, no O2 yet, for all the genome to
compare fully with stringtie and flair."*

## Substrate

| | |
|---|---|
| reads | `A119b.t2t.bam` — **68,026,217 records, 25 contigs**, 96 GB |
| genome | `chm13v2.0.fa` |
| reference | `chm13v2.0_RefSeq_full.gff.gz` → **188,903 mRNAs in 46,481 loci** (`gff_to_gtf` per contig, concatenated) |

⚠ **A119b, not `human_testis.t2t.bam`** — no testis BAM is on this disk, and A119b is the library the
lab's tool comparison used. Register 867 stands: the two human libraries are **not comparable** (A119b is
~6× deeper), so nothing here may be quoted beside a testis-panel number.

⭐ **StringTie 3.0.3 for every contig** — this run does NOT inherit register 870's defect, where the old
six-chromosome panel mixed 3.0.1 (chr20) with 3.0.3 (the rest).

## Ours — `--assemble-only`, shipped polish, no O2

```sh
copy_assign --assemble-only --assembly-polish full --polish-isoform-fraction 0.02 \
  --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 \
  --bam A119b.t2t.bam --fasta chm13v2.0.fa --region "$c:1-$len" --out $c
```
All 25 contigs, **zero errors**, and every log ends `0 families, 0 read assignments` — O2 confirmed off.
⚠ `RUSTLE_JUNCTION_MAJORITY` is now **default ON** (flipped 2026-09-21, register 960), so these numbers
are on the NEW default and are not comparable to pre-09-21 assembly numbers.

| metric | ours |
|---|---|
| transcripts | **264,996** in 69,899 loci (231,111 multi-exon) |
| base level (sens/prec) | 48.6 / 17.1 |
| exon level | 53.7 / 44.5 |
| intron level | 61.2 / 62.1 |
| **intron chain level** | **22.2 / 16.8** |
| transcript level | 20.7 / 14.7 |
| locus level | 34.6 / 22.4 |
| **matching intron chains** | **38,859** |
| matching transcripts | 39,043 |
| matching loci | 16,104 |
| missed loci | 20,650/46,481 (44.4%) |
| novel loci | 43,633/69,899 (62.4%) |

## Head to head — ours vs StringTie 3.0.3, genome-wide, same BAM, same reference

| metric (sens / prec) | **ours** | StringTie 3.0.3 |
|---|---|---|
| transcripts | 264,996 in 69,899 loci | 248,460 in 77,813 loci |
| base level | 48.6 / 17.1 | **59.1** / 14.5 |
| exon level | **53.7 / 44.5** | 53.0 / 41.1 |
| intron level | 61.2 / **62.1** | **63.3** / 59.0 |
| **intron chain level** | **22.2 / 16.8** | 18.8 / 15.8 |
| **transcript level** | **20.7 / 14.7** | 17.5 / 13.3 |
| locus level | **34.6 / 22.4** | 34.0 / 19.9 |
| **matching intron chains** | **38,859** | 32,871 |
| matching transcripts | **39,043** | 33,111 |
| matching loci | **16,104** | 15,795 |
| missed loci | 44.4% | **39.9%** |
| novel exons | **23.9%** | 36.0% |
| peak RSS | ~12.6 GB (chr1, per contig) | 3.7 GB (whole genome) |

⭐**Ours wins every chain-level metric** — intron chain and transcript, on BOTH sensitivity and precision
— and recovers **5,988 more matching intron chains (+18.2%)** from 6.7% more transcripts.
⚠**StringTie wins raw coverage**: base sensitivity 59.1 vs 48.6, intron sensitivity 63.3 vs 61.2, and it
misses fewer loci (39.9% vs 44.4%). It pays for that in precision everywhere and in **36.0% novel exons
vs our 23.9%**.
⚠StringTie is also far leaner: **3.7 GB peak for the whole genome in one process** against our ~12.6 GB
for chr1 alone (register 976).

## StringTie / FLAIR

- **StringTie**: `stringtie -L -p 4` genome-wide, single pass — **DONE**, 2h53m wall, 3.7 GB peak RSS.
- **FLAIR**: **NOT RUN genome-wide — deferred to a cluster (user's call, 2026-09-22).** Measured here so
  the cluster job can be sized:

  | | |
  |---|---|
  | stock recipe (`flair align` + `collapse`), chrY 41,464 primary reads | **11m35s** wall / 44m33s CPU |
  | `collapse` only, from a BAM-derived bed12 | **5m17s** |
  | throughput | **~131 reads/s** |
  | extrapolated genome-wide | **~25-57 h** on this box |

  ⚠ **Monolithic genome-wide FLAIR is impossible on this disk regardless of time**: at ~4.5 KB/read a
  whole-genome FASTQ is ~200-300 GB against 260 GB free, before FLAIR writes its own alignments.

  ⭐**`flair align` can be skipped.** `bamToBed -bed12` on our BAM reproduces `flair align`'s own bed
  **exactly in columns 1-4, 6-8 and 10-12** — every coordinate, read name and splice block, in 0.98 s
  instead of minutes. It differs only in col5 and col9, and the col5 difference is **real, not
  cosmetic**: col5 is MAPQ, and FLAIR's stock recipe aligns each contig against a **single-contig**
  reference, which inflates it (**31,823 reads at MAPQ 60 vs our genome-wide 29,888**) because no
  competing locus elsewhere in the genome is visible. Cross-check on chrY: stock **2,410** transcripts
  vs this route **2,436** — 1.1% apart.

  Four contigs completed before the run was stopped, kept in `/mnt/linuxdisk/tmp/gw22/flair/`:
  **chr20 20,894 · chr22 17,180 · chr21 16,457 · chrY 2,410** transcripts. Not scored here — a 4-contig
  subset is not comparable to the genome-wide ours/StringTie numbers above, and must not be quoted as if
  it were.

## Operational finding — register 976

`tools/genome_wide_sweep.sh` documents "3.5-4.5 GB per slot, 4 concurrent fit". **That is a property of
the shallower library it was measured on, not of `--assemble-only`.** On A119b: **chr1 12.6 GB, chr13
10.3 GB**, and two concurrent contigs drove the 25 GB box to 0 available with 8 GB of swap. The
documented `--jobs 4` would have OOM-killed it. Empirical scaling here is **~1.8 GB per million
records**; the run was repacked into 8 batches under a ~14 GB concurrent cap at `--jobs 3`, after which
memory never exceeded ~5 GB used and batches landed in **150-180 s each**.

## Reproduce

```sh
/mnt/linuxdisk/tmp/gw22/run_ours_rest.sh      # memory-aware batched sweep, then merge
/mnt/linuxdisk/tmp/gw22/run_stringtie.sh      # waits for ours to clear memory
/mnt/linuxdisk/tmp/gw22/flair_fast.sh CHR ...  # per contig; genome-wide needs a cluster
gffcompare -r ref_genome.gtf -o cmp_ours ours_genome.gtf
```

## Addendum (§6z5/§6z6, later the same day) — the lab's own tool runs, both species

The lab's genome-wide StringTie 3.0.1 / FLAIR 3.0.1 / isoseq collapse outputs (`../benchmark_collapse/`,
`../isoseq_upload/`, same BAMs) were scored against the annotation with the pre-registered 15-cell rule in
`docs/PREREG_external_tool_bakeoff_2026-09-22.md`: **gorilla 10/15, human chr20/21/22 12/15**; the human
whole-genome FLAIR/isoseq rows are scored on the cluster (`benchmark_collapse/score_vs_annotation.sh`, our
GTFs shipped there). Ours vs the lab's StringTie 3.0.1 genome-wide: 22.2/16.8 vs 18.1/15.0 intron chain,
38,859 vs 31,792 matching chains. What the tools match and we do not is single-read isoforms (register
1057), plus one real defect — `copy_assign`'s coordinate de-duplication of primary reads (register 1058,
`docs/PREREG_primary_dedupe_2026-09-22.md`), fixed behind `--keep-coordinate-duplicates` and NOT made the
default because the shipped polish was fitted on the de-duplicated counts.
