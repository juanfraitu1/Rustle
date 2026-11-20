# Row 1831: are the single-exon rules and the end layers overfit to primates? (pre-registered 2026-11-20)

User (2026-11-20): reduce assembly rules that only benefit primates; take the levers no ablation covered. Rows 1793 and 1811
switched every chain-level lever; they scored exact intron chains, which cannot see single-exon transcripts or transcript ends.
The six levers below were all set on A119b (chr20 or held-out chr11 + chr7) and never tested elsewhere. Diagnostic: nothing here
changes a default. Workspace `r1831_mono_end_rules` (`arms.sbatch`, `score.py`). Registered before any arm ran.

## Arms (one switch against `base` = today's defaults, stage `assemble`, `RUSTLE_PTC=0`)
| arm | change | lever | read on |
|---|---|---|---|
| shadow_off | driver without `--polish-mono-shadow` | shadow rule: drop a single-exon transcript under a spliced gene | `x.gtf` |
| quant0 | `--polish-mono-quantile 0` (default 0.82) | read floor of single-exon transcripts = the 0.82 quantile of multi-exon support | `x.gtf` |
| primed_off | `RUSTLE_POLISH_MONO_PRIMED=0` | drop internally primed single-exon transcripts | `x.gtf` |
| monomodel_off | `RUSTLE_MONO_MODEL=off` | single-exon models rebuilt from read-end clusters (row 1489) | `x.final.gtf` |
| tssend_off | `RUSTLE_POLISH_TSS_END=off` | 5' end = median read start (row 1487) | `x.final.gtf` |
| tesext_off | `RUSTLE_POLISH_TES_EXTEND=off` | 3' UTR extension (row 1488) | `x.final.gtf` |

## Datasets
Non-primate (N): fly simulated, mouse simulated, zebrafish simulated (complete truths of rows 1779 / 1800), mouse ES real (mouse
contigs, GENCODE vM27), mouse heart real (GENCODE vM27). Primate (P): OR6737 (`apetools/GGO.refseq.gtf`), A119b (RefSeq on
CHM13, `ptr_ppa_rules/gffcmp/hsa.ref.gtf`). Zebrafish runs on the re-mapped BAM of row 1830 (the original is on an
unavailable storage target).

## Measures (`score.py`)
- **Single-exon:** a model MATCHES when it overlaps a truth single-exon transcript by >= 50% of both lengths, strand equal unless
  the model has none. Model precision = matching models / models; truth single-exon transcripts FOUND.
- **Ends:** among transcripts whose intron chain is correct, the share whose 5' (3') end lies within 50 bp of the nearest truth
  transcript of that chain.
- **Bases:** strand-aware exonic base sensitivity and precision.
- Multi-exon chains and correct chains per arm, to show whether a lever touches them.

## Reading, fixed now (lever effect = base against the arm without the lever)
- **Single-exon filters** (shadow, quantile, primed) and **mono-model:** the lever *pays* on a dataset when single-exon model
  precision is higher with it and it costs at most 5% of the truth single-exon transcripts found without it (row 1793's reading).
- **tssend:** *helps* when the 5' share within 50 bp is at least 1 point higher with it.
- **tesext:** *helps* when base sensitivity is at least 0.5 points higher with it and base precision not more than 0.5 points
  lower (the benefit row 1488 registered). The 3' share within 50 bp is reported beside it.
- **Verdict per lever:** GENERAL = pays / helps on both P and on at least 4 of 5 N. OVERFIT = on both P, and fails on at least 3
  of 5 N. NO EFFECT = the measured quantities differ by less than 0.1% on every dataset. Otherwise MIXED, with the table.
- **Decision.** A table for the user. OVERFIT or MIXED levers are candidates for toning down or switching off; a changed default
  is the user's call and, per the user's ruling on reuse, is checked on one more dataset not used here (mouse adrenal or cortex).

## Limits, stated in advance
- Simulated reads have no intronic, antisense or internally primed unspliced signal, the very thing the three single-exon filters
  remove. On the simulated sets these filters can show a cost but hardly a benefit; mouse ES and heart carry the benefit side.
- Annotated single-exon genes are few and annotations incomplete: an unmatched single-exon model on a real library can be real.
- Annotated transcript ends are weak truth on real libraries (RefSeq 3' ends on A119b; rows 1487 / 1488 used held-out
  chromosomes and CAGE / polyA sites). The simulated truths are exact at both ends. For P, the registered results of rows
  1487-1489 stand beside this table.
- One switch at a time; `final.gtf` also passes the chain score filter and the 5' truncation rule.
- The scorer was tested on the `base`-like fly files of row 1830 (default output, no arm) before registration.

## Result (2026-11-20): no lever is GENERAL; mono-model OVERFIT by the reading; the 3' extension acts on primates only

Jobs 4484140-4484146, every arm exit 0; tables `r1831_mono_end_rules/score.NAME.txt`. Repo `b7de7072`, binaries of 2026-11-20.
The `base` chain sets equal row 1817's on every dataset (checked chain by chain on mouse heart).

Datasets in the order fly, mouse simulated, zebrafish, mouse ES, mouse heart | OR6737, A119b. "with -> without" = base -> arm.

**Single-exon filters, `x.gtf`** (model precision with / without the lever; truth single-exon transcripts found with / without):
| lever | fly | mouse sim | zebrafish | mouse ES | heart | OR6737 | A119b | verdict |
|---|---|---|---|---|---|---|---|---|
| shadow, precision | 0.921 / 0.905 | 0.966 / 0.967 | 0.861 / 0.799 | 0.121 / 0.076 | 0.557 / 0.178 | 0.216 / 0.043 | 0.021 / 0.007 | |
| shadow, found | 343 / 409 | 85 / 117 | 144 / 165 | 275 / 375 | 59 / 82 | 74 / 142 | 234 / 341 | MIXED: pays nowhere (cost 13-48%) |
| quantile 0.82, precision | 0.921 / 0.877 | 0.966 / 0.958 | 0.861 / 0.614 | 0.121 / 0.032 | 0.557 / 0.144 | 0.216 / 0.057 | 0.021 / 0.011 | |
| quantile 0.82, found | 343 / 1,800 | 85 / 850 | 144 / 927 | 275 / 532 | 59 / 221 | 74 / 211 | 234 / 396 | MIXED: pays nowhere (cost 41-90%) |
| primed, precision | 0.921 / 0.926 | 0.966 / 0.961 | 0.861 / 0.864 | 0.121 / 0.095 | 0.557 / 0.576 | 0.216 / 0.210 | 0.021 / 0.011 | |
| primed, found | 343 / 423 | 85 / 126 | 144 / 155 | 275 / 287 | 59 / 120 | 74 / 85 | 234 / 241 | MIXED: pays on mouse ES and A119b only |

- All three raise single-exon precision on the real libraries (the primed rule not on heart) and all three fail the 5% cost bar
  nearly everywhere. None is primate-specific: the pattern is the same on mouse and on the apes.
- The quantile floor is the expensive one. On the simulated sets it removes 81-90% of the true single-exon transcripts that are
  found without it, where precision without it is still 0.88-0.96 on fly and mouse (0.61 on zebrafish). Fly truth has 3,656
  single-exon transcripts: 343 found by default, 1,800 without the floor.
- The primed rule removes true single-exon transcripts on simulated reads (7-33% of those found) and half of them on heart
  (59 / 120), at no precision gain there.
- None of the three changes a multi-exon chain of `x.gtf`.

**mono-model, `x.final.gtf`** (precision with / without; found with / without):
| fly | mouse sim | zebrafish | mouse ES | heart | OR6737 | A119b | verdict |
|---|---|---|---|---|---|---|---|
| 0.989 / 0.921; 342 / 343 | 0.969 / 0.966; 301 / 85 | 1.000 / 0.861; 80 / 144 | 0.103 / 0.121; 270 / 275 | 0.351 / 0.557; 128 / 59 | 0.222 / 0.216; 134 / 74 | 0.169 / 0.021; 315 / 234 | **OVERFIT** (pays on both P, fails zebrafish, ES, heart) |
The failures differ in kind: zebrafish loses 44% of the found truth at higher precision; ES and heart lose precision (heart finds
twice as many). The large gain is A119b's (0.021 -> 0.169).

**tssend, `x.final.gtf`** (share of correct-chain transcripts with the 5' end within 50 bp, with / without):
| fly | mouse sim | zebrafish | mouse ES | heart | OR6737 | A119b | verdict |
|---|---|---|---|---|---|---|---|
| 0.979 / 0.994 | 0.982 / 0.994 | 0.980 / 0.995 | 0.674 / 0.500 | 0.634 / 0.655 | 0.368 / 0.602 | 0.713 / 0.461 | MIXED: helps mouse ES (+17 points) and A119b (+25); worse on the other five, OR6737 by 23 points |
With the layer `final.gtf` also holds 3-49 fewer correct chains (it changes the features the chain score filter reads) and base
sensitivity is 0.4-1.5 points lower on all seven datasets.

**tesext, `x.final.gtf`** (base sensitivity with / without; base precision with / without; 3' share within 50 bp with / without):
| fly | mouse sim | zebrafish | mouse ES | heart | OR6737 | A119b | verdict |
|---|---|---|---|---|---|---|---|
| 0.6836 / 0.6835 | 0.7392 / 0.7389 | 0.7579 / 0.7578 | 0.3562 / 0.3552 | 0.2553 / 0.2549 | 0.4443 / 0.4388; 0.692 / 0.696; 0.608 / 0.601 | 0.4727 / 0.4534; 0.410 / 0.429; 0.596 / 0.557 | MIXED: helps OR6737; A119b +1.9 points sensitivity for -1.9 precision (outside the bound) |
On the five non-primate datasets the layer moves base sensitivity by at most 0.1 point: its effect exists on the two ape testis
libraries only (as transfrag extension, row 1811). It also removes single-exon models there (A119b 14,121 -> 9,212, model
precision 0.111 -> 0.169).

**Reading for the user (my interpretation, beyond the registered verdicts).**
1. The registered question "overfit to primates?" gets a yes for one lever (mono-model) and "acts on primates only" for a second
   (tesext). The other four are not primate-specific.
2. The larger finding is single-exon recall: the three `x.gtf` filters together keep 9% of fly's true single-exon transcripts.
   On real libraries they buy precision that stays low anyway (A119b 0.02, mouse ES 0.12).
3. tssend is right where it was tuned (A119b) and on mouse ES, and wrong on OR6737 and slightly on four others. The row cannot
   say why; A119b and mouse ES are the two libraries where the earlier cap-start work found usable 5' signal, and OR6737 is one
   where it did not (not tested here).
Limits as stated before the run; in addition, several models can match one truth transcript, so "matching models" exceeds
"found" in `final.gtf`, and precision there counts models.

Candidates, the user's call (none made): quantile floor lower than 0.82 (a grid, checked on mouse adrenal or cortex); mono-model
and tesext as opt-in layers; tssend conditional on a 5' signal test.
