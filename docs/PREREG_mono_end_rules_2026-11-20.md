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
