# Pre-registration: %LRC (LRGASP long-read coverage) as a reporting metric — 2026-09-29

Written 2026-09-29 before any %LRC number of ours exists. Source: `docs/archive/2026-09/VAULT_METHODS_AUDIT_2026-09-28.md`, "Not in
the code", item 1. **A reporting metric only: no filter, no default change, no decision rule, no ranking of tools.**
The audit already says it moves none of O1/O2/O3; it only tightens assembly-level statements ("the models are read
evidence end to end").

## Source definition (quoted)

LRGASP, Pardo-Palacios et al. 2024, *Nat. Methods* 21:1349-1363, Box 1 (Challenge 1 metrics):

> %LRC | Fraction of the transcript model sequence length mapped by one or more long reads.

Results text: "We tested this by measuring long-read coverage (LRC) of the transcript predictions from our read
alignments. FLAMES, Iso_IB, IsoTools, LyRic and Mandalorion showed nearly complete LRC for their transcript models
(>98% coverage for all transcripts). In contrast, FLAIR, Spectra, TALON, IsoQuant, Bambu and StringTie2 had lower
coverage rates (~90, 90, 85, 75, 60 and 45%, respectively) (Extended Data Fig. 2), suggesting that they may use
different alignment strategies or additional information (for example, reference annotation or short reads) to
finalize transcript models." Extended Data Fig. 2: "Percentage of transcript models with different ranges of sequence
coverage by long reads", in three classes **> 98 %, 98-75 %, < 75 %**.

What the sources leave open (checked 2026-09-29): the paper gives no per-model procedure. The published figure code
(`LRGASP/Challenge1_Figures_Code`, `Functions_Supplementary_Figures_Challenge1_v4.R`, function `LRC()`) only reads
precomputed per-model fractions from `coverage_files/*.txt` and bins them (its third class is coded `x < 0.70` while
labelled "<75%"; we use the printed labels). The LRGASP challenge-1 evaluation script
(`sqanti3_lrgasp.challenge1.py`) does not compute it, and SQANTI3 5.5.4 (the local install) has no %LRC column. The
operational choices below are therefore ours, fixed here before looking.

## Operational definition

- **Model**: one GTF transcript (`transcript_id` on one contig); its `exon` records, merged where they overlap or touch.
  L = the model's exonic length (the "transcript model sequence length").
- **Reads**: the sample's registry BAM (`figures/samples.tsv`), the same BAM every arm was assembled from. **Primary
  alignments only**: records with none of the flags 0x4, 0x100, 0x800 (`samtools view -F 2308`, the repo invariant).
  No MAPQ filter and no other flag filter.
- **Multi-mapped reads**: a read counts once, at its primary placement (MAPQ-0 primaries included); its secondary and
  supplementary records are not read. Consequence, stated in advance: in duplicated regions the primary of a tied
  read sits on one copy chosen by the aligner, so a copy's model can have %LRC < 1 although the reads fit it equally
  well. %LRC does not resolve that ambiguity (O2 does), and the default Rustle arm, whose loci are seeded with
  secondary alignments within 2 % of the best score, is expected to lose %LRC on exactly those models;
  `rustle_primary` (`--no-seed-secondaries`) is the like-for-like arm for this metric. A supplementary piece of a
  chimeric read does not count either (fused / readthrough models can lose %LRC for that reason).
- **Covered base**: a reference base under an aligned read base, CIGAR `M`, `=` or `X` (pysam `get_blocks`, the same
  bases `samtools depth` counts without `-J`). Introns (`N`), deletions (`D`) and clips are not covered. **Strand is
  ignored** (a read covers a base whatever its orientation relative to the model).
- **Exon level, not span level**: coverage is by aligned blocks, so a read spliced across a model's retained intron
  does not cover that intron. A span-level variant (first to last aligned base of a read) is not computed.
- **%LRC(model) = covered exonic bases / L**, reported as a fraction in [0, 1].
- **Summaries per sample x arm** (all models, and the multi-exon models separately): number of models, mean, median,
  and the share of models in each LRGASP class: **> 0.98**, **0.75 to 0.98** (both ends included), **< 0.75**; plus
  the shares with %LRC exactly 1 and exactly 0.
- **Per SQANTI3 structural category** (FSM / ISM / NIC / NNC / Other, the fig. 2 fold): the same summary, only where a
  fig. 2 SQANTI3 classification is cached (today: the development scope, gorilla chr20 / chr22 / chrY and human
  chr20-22), scoring the exact per-contig GTF SQANTI3 classified and joining by isoform id; a classification whose
  input GTF changed since SQANTI3 ran is skipped and counted.

## Scope

- Samples: the six registry samples, each against its own BAM, on its evaluation contigs
  (`assembly.evaluation_contigs`: every contig its annotation covers; human leaves out chrM). **Never pooled across
  samples or species**; every row carries `species` and `sample`.
- Arms: Rustle's two genome-wide GTFs from the run cache (`rustle`, `rustle_primary`) on every sample; the lab's
  annotation-free StringTie (`-L`, no `-G`), FLAIR (collapse without annotation) and IsoSeq collapse GTFs on
  human_A119b and gorilla_OR6737 (the only samples that have them). Those three are trusted baselines, described next
  to Rustle, not competitors.
- Code: `figures/_lrc.py` (the BAM pass is cached per contig and resumable inside a per-call budget, so every call is
  a light job), unit test `figures/test_lrc.py` on a hand-built BAM and GTF with the expected values computed by hand.

## Expectations written before looking (descriptive, not decision rules)

- E1. `rustle_primary`: >= 98 % of models with %LRC > 0.98 on every sample (its models are built from the same
  primary alignments; a shortfall would point at exon edges the assembler places beyond the aligned bases, e.g. end
  extension or small-gap bridging).
- E2. `rustle` (default): below `rustle_primary`, with the gap concentrated in models supported by secondary
  alignments. The gap is a property of the primary-only read set, not a defect of the models.
- E3. StringTie, FLAIR and IsoSeq collapse were run without annotation on the same BAMs, so each is expected to be
  high as well (>= 90 % of models > 0.98); LRGASP's low StringTie2 / FLAIR values came from runs that could use the
  annotation or short reads, which ours did not.

## What it is used for

A column in assembly-level reporting (context for figures 1-2 and `bench/SQANTI3_POLISH.md`): the share of models
covered end to end by primary read alignments, per sample and arm. It changes no default, filters nothing, ranks no
tool and feeds no decision. A secondary-inclusive or span-level variant, a MAPQ filter, or any use as a filter needs an
amendment here first.

## Outcome (2026-09-29, same day, after the code and test above; nothing in the definition was changed after looking)

**Checks.** `python3 figures/test_lrc.py`: 7/7 pass (hand-built BAM + GTF, expected values by hand in the test's
docstring); a mutant that counts secondary and supplementary records fails 2 of them. The union of human_A119b chr21
(24,512,279 bases) equals `samtools depth -G 2308 -r chr21 | awk '$3>0' | wc -l` exactly. Every classified SQANTI3
isoform joined (0 stale). Cost: the BAM pass runs at ~60 k records/s (gorilla_OR6737 ~3 min, human_A119b ~15 min, in
calls of <= 110 s); scoring one arm takes <= 60 s (human IsoSeq collapse, 3.2 M models, 0.97 GB peak). Caches:
`${work}/lrc/` (a symlink to `/mnt/linuxdisk/tmp/rustle_figures_dev/lrc/cache`); tables:
`/mnt/linuxdisk/tmp/rustle_figures_dev/lrc/tables/lrc_summary.tsv`, `lrc_sqanti.tsv`.

**Genome-wide, all models, each sample on its evaluation contigs** (never pooled; `>0.98`, `<0.75`, `=0` = share of
models):

| species | sample | arm | models | mean %LRC | > 0.98 | < 0.75 | = 0 |
|---|---|---|---:|---:|---:|---:|---:|
| human | human_A119b | rustle | 252,178 | 0.9908 | 97.80% | 1.16% | 0.40% |
| human | human_A119b | rustle_primary | 238,367 | 0.99993 | 99.998% | 0 | 0 |
| human | human_A119b | StringTie | 250,491 | 0.99985 | 99.992% | 0 | 0 |
| human | human_A119b | FLAIR | 1,096,234 | 0.99993 | 99.999% | 0 | 0 |
| human | human_A119b | IsoSeq collapse | 3,195,506 | 0.99991 | 99.9995% | 0 | 0 |
| human | human_testis | rustle | 25,237 | 0.9839 | 97.01% | 1.79% | 1.08% |
| human | human_testis | rustle_primary | 24,233 | 0.99994 | 99.988% | 0 | 0 |
| gorilla | gorilla_OR6737 | rustle | 74,045 | 0.9939 | 98.94% | 0.72% | 0.34% |
| gorilla | gorilla_OR6737 | rustle_primary | 72,690 | 0.99994 | 100% | 0 | 0 |
| gorilla | gorilla_OR6737 | StringTie | 68,249 | 0.99986 | 99.991% | 0 | 0 |
| gorilla | gorilla_OR6737 | FLAIR | 153,402 | 0.99994 | 99.999% | 0 | 0 |
| gorilla | gorilla_OR6737 | IsoSeq collapse | 551,342 | 0.99985 | 99.999% | 0 | 0 |
| gorilla | gorilla_KB3781 | rustle | 72,622 | 0.9958 | 99.40% | 0.46% | 0.33% |
| gorilla | gorilla_KB3781 | rustle_primary | 71,855 | 0.99998 | 99.997% | 0 | 0 |
| chimpanzee | chimp_PTR | rustle | 54,821 | 0.9953 | 99.07% | 0.56% | 0.23% |
| chimpanzee | chimp_PTR | rustle_primary | 53,909 | 0.99995 | 99.996% | 0 | 0 |
| orangutan | orangutan_PPY | rustle | 111,629 | 0.9821 | 96.98% | 2.21% | 1.02% |
| orangutan | orangutan_PPY | rustle_primary | 106,467 | 0.99991 | 99.991% | 0 | 0 |

Multi-exon models only: the same picture (rustle 98.6% / 97.8% / 99.3% / 99.6% / 99.4% / 98.2% > 0.98 in the table's
sample order; every other arm >= 99.98%). The share at exactly 1 (79.6-97.7%, lowest for StringTie) differs from the
> 0.98 share only by slivers under 2% of a model (bases under read deletions, which the definition does not count,
and model ends); not interpreted.

**Per SQANTI3 category (fig. 2 development scope, `lrc_sqanti`).** rustle_primary and all three baselines: >= 99.88%
of models > 0.98 in every category, both species. The default rustle arm's shortfall sits in FSM and Other: human
chr20-22 FSM 90.3% > 0.98 (7.7% < 0.75) and Other 89.3% (8.0% < 0.75) vs ISM / NIC / NNC >= 99.3%; gorilla chr20 /
chr22 / chrY FSM 97.6% (1.6% < 0.75) and Other 91.1% (7.6% < 0.75) vs ISM / NIC / NNC >= 99.0%.

**Expectations.** E1 holds (rustle_primary >= 99.98% > 0.98 on all six samples). E2 holds: the default arm is below
rustle_primary on every sample (by 0.6 to 3.0 points of the > 0.98 share), and the gap is its secondary-seeded models:
of rustle's models below 0.75, 519/536 (96.8%, gorilla_OR6737) and 2,894/2,922 (99.0%, human_A119b) have an intron
chain (mono-exon: span) absent from rustle_primary, and every zero-coverage model (253 and 999) is such a chain
(diagnostic run 2026-09-29 on the two benchmark samples; script in the session scratchpad, not in the repo). E3 holds
(StringTie, FLAIR, IsoSeq collapse >= 99.99% > 0.98 in both species).

**Reading.** For assemblers that use only these reads and no annotation, %LRC on primary alignments is saturated
(>= 99.98% of models > 0.98), so it does not separate rustle_primary from the three baselines; LRGASP's spread came
from runs that could use the annotation or short reads. Its only signal here is the default arm: 0.6-3.0% of its
models are covered < 98% by PRIMARY alignments, and these are almost all the models built from secondary alignments
within 2% of the best score (mostly FSM copies of annotated genes and "Other"), i.e. the multi-copy models the default
seeding exists to build. That is the stated cost of primary-only counting, not missing read evidence. Caveat on
strength: coverage is strand-blind and the libraries carry broad background (primary aligned bases cover 1.79 Gb of
CHM13 for human_A119b, 0.70 Gb orangutan, 0.43 Gb gorilla_OR6737, 0.43 Gb chimpanzee, 0.36 Gb gorilla_KB3781,
0.05 Gb human_testis), so a high %LRC is weak evidence for a model's structure. It shows that the model's bases
carry reads, not that the model is correct. Nothing changes: no filter, no default, no ranking.
