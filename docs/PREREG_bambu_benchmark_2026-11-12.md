# PREREG r1641: bambu as a fourth tool, de novo and guided kept apart (2026-11-12)

Advisor request (via user, 2026-11-12): compare the pipeline with bambu (github.com/GoekeLab/bambu); test tools for their
de novo ability; guided may be shown too, but as a separate block so every transcript is attributable to one mode.
Written after a one-chromosome pilot (gorilla NC_073242.2: runtime and file format only, no score read) and before any
whole-sample run.

## Runs (bambu 3.12.1, bioconda, env `/anvme/workspace/v138da13-rustle/bambu_env`; `bench/famcmp/bambu_run.R`)
Samples: gorilla OR6737 (`GGO_mm.bam`) and human A119b (`A119b.t2t.bam`), the BAMs the r1619 baselines used.
- **De novo:** `annotations = NULL`, all other options default. Two discovery calls on one set of read classes: bambu's
  own default NDR (the default point) and `NDR = 1` (every candidate, with its NDR, for the curve).
- **Guided:** `annotations` = the r1622 transcript-split guide (`t.guide.gtf` gorilla, `ht.guide.gtf` human; hidden
  transcripts are the truth for novel calls), same two calls (bambu's recommended NDR, and NDR = 1).
- The two modes never share read classes (separate `rcOutDir`): bambu's read-class step uses the annotation.
- Ours guided on the same splits: `utilities guide` (G1d, the shipped rule) on the r1614 / r1617 default de novo GTFs.
- Baselines already on disk and checked annotation-free: StringTie `-L` (no `-G`), FLAIR collapse (no `-f`).

## Measures (multi-exon exact intron chains vs RefSeq; descriptive, no adoption decision hangs on it)
- De novo block (`bench/famcmp/mode_table.py`): chains / correct / precision at each tool's default point; curves
  ours by reads, StringTie by cov, bambu by NDR. Read curve against curve at matched correct count, never a tuned point
  against a default point.
- Guided block (`guide_proto.py score`, GUIDE_HIDE=tx): hidden chains found and precision outside the guide, for ours
  guided, StringTie -G, bambu guided (default NDR and NDR = 1). Bambu writes every guide transcript whether expressed or
  not, so its "guide chains found" column is not comparable and is not reported.
- A transcript in a de novo file never saw an annotation; a transcript in a guided file is labelled by its file.

## Stated limits
- Without annotations bambu cannot train on the sample and scores with its shipped pretrained model; the default point
  is then the least reliable figure, the NDR curve the fairer one.
- bambu's defaults drop candidates with < 2 reads and transcripts contained in another (`remove.subsetTx`); left as is.
- Wall limits 4 h (gorilla) / 20 h (A119b); a run that hits its limit is reported as not finished, not rerun tuned.

## Outcome, gorilla OR6737 (2026-11-12, job 4466673, `r1641_bambu/g.*`; bambu 1,462 s de novo, 1,517 s guided, 72-core node)
First attempt failed on a missing R package (BiocFileCache, needed for `rcOutDir`); added to the env, rerun unchanged.

De novo block (no tool saw an annotation):
| tool, point | chains | correct | precision |
|---|---|---|---|
| ours, default | 73,348 | 25,794 | 0.352 |
| StringTie, default | 67,491 | 23,248 | 0.344 |
| FLAIR | 138,125 | 24,293 | 0.176 |
| bambu, every candidate (NDR 1) | 31,947 | 15,381 | 0.481 |
| bambu, default (NDR 0.5 chosen by bambu) | 16,202 | 10,982 | 0.678 |
| ours, reads >= 6 / >= 11 / >= 21 | 31,542 / 20,362 / 12,788 | 16,965 / 12,996 / 9,429 | 0.538 / 0.638 / 0.737 |
| StringTie, cov >= 3 / >= 6 | 31,578 / 18,303 | 15,644 / 11,121 | 0.495 / 0.608 |
- bambu's ceiling is 15,381 correct chains; ours and StringTie reach 25,794 and 23,248.
- At matched correct count (linear interpolation): at 15,381 correct bambu 0.481, StringTie about 0.50, ours about 0.58;
  at bambu's default 10,982 correct bambu 0.678, ours about 0.69, StringTie about 0.61. Below about 9,400 correct our
  curve has no point (last cut is 21 reads); bambu continues to 0.777 at 8,271 and 0.881 at 3,657.
- 15,143 of bambu's 15,381 correct chains are also ours; bambu-only chains 3,097, 57 correct.

Guided block (transcript split, 29,964 hidden chains; every tool saw the same guide):
| tool | hidden found | chains outside the guide | precision outside the guide |
|---|---|---|---|
| ours guided (G1d) | 5,524 | 15,531 | 0.356 |
| StringTie -G | 4,538 | 15,322 | 0.296 |
| bambu guided, NDR 1 | 4,093 | 14,219 | 0.288 |
| bambu guided, its recommended NDR | 2,420 | 4,246 | 0.570 |
- At about the same number of chains outside the guide, ours finds more hidden chains than StringTie -G and bambu. Bambu's
  recommended point is a stricter trade (fewer, more precise); our guided output was not cut to that size, so that point
  has no matched comparison.
A119b: job 4466674 running.

## Outcome, human A119b (2026-11-12, job 4466674, `r1641_bambu/a.*`; bambu 6,575 s de novo, 6,616 s guided)
De novo block:
| tool, point | chains | correct | precision |
|---|---|---|---|
| ours, default | 219,517 | 38,715 | 0.176 |
| StringTie, default | 211,865 | 31,781 | 0.150 |
| FLAIR | 626,983 | 39,230 | 0.063 |
| bambu, every candidate (NDR 1) | 41,745 | 12,793 | 0.306 |
| bambu, default (NDR 0.5) | 9,599 | 6,680 | 0.696 |
- 12,482 of bambu's 12,793 correct chains are also ours; bambu-only chains 5,147, 108 correct.

Guided block (transcript split, 54,491 hidden chains):
| tool | hidden found | chains outside the guide | precision outside the guide |
|---|---|---|---|
| ours guided (G1d) | 7,963 | 57,687 | 0.138 |
| StringTie -G | 5,334 | 42,238 | 0.126 |
| bambu guided, NDR 1 | 3,630 | 29,630 | 0.123 |
| bambu guided, its recommended NDR | 1,108 | 1,690 | 0.656 |

## De novo precision at bambu's correct counts (both samples; ours cut by reads, StringTie by cov, first cut reaching the count; post hoc, replaces the interpolated figures above)
| sample | correct chains (bambu NDR cut) | bambu | ours | StringTie |
|---|---|---|---|---|
| gorilla | 15,381 (1) | 0.481 | 0.563 | 0.495 |
| gorilla | 10,982 (0.5, default) | 0.678 | 0.686 | 0.608 |
| gorilla | 8,271 (0.3) | 0.777 | 0.766 | 0.681 |
| gorilla | 3,657 (0.1) | 0.881 | 0.870 | 0.826 |
| A119b | 12,793 (1) | 0.306 | 0.490 | 0.333 |
| A119b | 9,688 (0.7) | 0.491 | 0.554 | 0.378 |
| A119b | 6,680 (0.5, default) | 0.696 | 0.625 | 0.438 |
| A119b | 4,185 (0.3) | 0.861 | 0.683 | 0.496 |
| A119b | 1,169 (0.1) | 0.953 | 0.745 | 0.679 |
**Reading:** bambu de novo reaches a third (A119b) to 0.6 (gorilla) of our correct chains and cannot go further. Inside its
range: level with ours on gorilla; on human it is clearly more precise at its default cut and stricter (0.696 vs 0.625,
0.861 vs 0.683), ours is more precise when bambu is opened up. Its score ranks chains better than our read count does on
human. Guided: ours finds the most hidden chains on both samples; bambu's recommended cut is small and precise, unmatched
by any cut of ours that was run.

## Corrections and disclosures (2026-11-12, after both outcomes)
- Non-default bambu setting besides those registered: `lowMemory = TRUE` (reads processed per chromosome) in both modes.
- "Every candidate (NDR 1)" opens only the NDR filter. bambu's other discovery defaults still apply (min.readCount 2,
  remove.subsetTx, min.readFractionByGene 0.05), so bambu's top correct count is its ceiling at default discovery
  options, not what the method could reach with them relaxed (no relaxed arm run).
- Guided recommended NDR from the logs: 0.367 gorilla, 0.282 A119b. The de novo logs print no threshold; 0.5 is inferred from the default call matching the NDR <= 0.5 cut exactly on both samples. bambu calibrates it against the guide, in which the hidden
  transcripts are absent, so it counts them as novel when estimating annotation completeness.
- Multi-exon bambu chains without a strand (cannot match): 4 gorilla, 32 A119b at NDR 1. Negligible.
- The interpolated gorilla figures (0.58 / 0.69) are superseded by the matched-count table (0.563 / 0.686).

## Follow-up (2026-11-12)
IsoQuant and guided FLAIR were added as register row 1646: `docs/PREREG_isoquant_benchmark_2026-11-12.md` has the
five-tool tables for both modes, precision at IsoQuant's correct counts, and the reading on whether bambu "makes up"
transcripts (its non-RefSeq chains are the most corroborated by other tools, 0.83-0.84; it under-calls instead).
