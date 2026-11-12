# PREREG r1646: IsoQuant as a fifth tool, de novo and guided kept apart (2026-11-12)

User / advisor (2026-11-12): the advisor suspects bambu invents or over-fits to known transcripts; add IsoQuant in both
modes on gorilla OR6737 and A119b. Same design as r1641 (`docs/PREREG_bambu_benchmark_2026-11-12.md`): same BAMs, same
RefSeq truth, same r1622 transcript-split guide, same scorers; bambu stays in the tables. Written after a one-chromosome
pilot (gorilla NC_073242.2: runtime, file format, and how many chains the polyA setting moves; no score against RefSeq read).

## Runs (IsoQuant 4.1.0, bioconda, env `/anvme/workspace/v138da13-rustle/isoquant_env`)
`isoquant -d pacbio_ccs --bam BAM --reference FA [--genedb GUIDE]`, everything else default, four runs per sample:
- de novo (no `--genedb`) and guided (`--genedb` = the guide), each with
  - **`--polya_trimmed all` (primary arm):** the reads are polyA-trimmed FLNC (0.1% tails detected) and IsoQuant asks for this
    setting on such data; and
  - **shipped default (second arm):** what a user gets without reading the warning. On the pilot the two differ
    (1,636 vs 1,722 multi-exon chains, 1,501 shared).
- Guide file: transcript ids prefixed `T.` (IsoQuant refuses a GTF where a transcript id equals its gene id); coordinates unchanged.
- Scored file: `PREFIX.transcript_models.gtf` (in guided mode = expressed known + novel; not the extended annotation).

## Measures (multi-exon exact intron chains vs RefSeq; descriptive)
- De novo block: default point of every tool; IsoQuant cut by its own transcript counts; precision of every tool at
  matched correct counts (first cut reaching the count), as in the r1641 table.
- Guided block: hidden chains found, chains outside the guide, precision outside the guide.
- The advisor's two questions, read from the same tables: (a) invented transcripts = chains only one tool reports and
  their precision; (b) known-biased = hidden chains found in guided mode relative to StringTie -G and ours.
- Wall limits 4 h gorilla / 20 h A119b; a run that hits its limit is reported as not finished.

## Outcome, gorilla OR6737 (2026-11-12, job 4466747, `r1646_isoquant/g.*`; IsoQuant 402 s de novo, 7,556 s guided of which ~1 h 50 is the guide database build)
De novo block (pt = `--polya_trimmed all`, primary):
| tool, point | chains | correct | precision |
|---|---|---|---|
| ours, default | 73,348 | 25,794 | 0.352 |
| StringTie, default | 67,491 | 23,248 | 0.344 |
| FLAIR | 138,125 | 24,293 | 0.176 |
| bambu, NDR 1 / default | 31,947 / 16,202 | 15,381 / 10,982 | 0.481 / 0.678 |
| IsoQuant pt | 37,089 | 16,677 | 0.450 |
| IsoQuant default | 40,182 | 17,409 | 0.433 |
| IsoQuant pt, count >= 3 / 6 / 11 / 21 | 29,856 / 20,313 / 14,423 / 9,830 | 15,244 / 12,234 / 9,704 / 7,273 | 0.511 / 0.602 / 0.673 / 0.740 |
- Chains only one tool reports (and correct): ours 9,443 (805), StringTie 16,994 (1,844), FLAIR 82,045 (2,015), bambu 2,816 (44),
  IsoQuant pt 615 (25), IsoQuant default 428 (13).

Guided block (29,964 hidden chains):
| tool | guide chains found | hidden found | chains outside the guide | precision outside the guide |
|---|---|---|---|---|
| ours guided (G1d) | 24,959 | 5,524 | 15,531 | 0.356 |
| StringTie -G | 32,166 | 4,538 | 15,322 | 0.296 |
| bambu guided, NDR 1 / recommended | not comparable | 4,093 / 2,420 | 14,219 / 4,246 | 0.288 / 0.570 |
| IsoQuant guided pt | 20,676 | 4,761 | 24,516 | 0.194 |
| IsoQuant guided default | 23,969 | 4,768 | 26,677 | 0.179 |

## Amendment 1 (2026-11-12, user: "can we include flair in the comparisons?"; written before the run)
FLAIR is already in the de novo block (collapse without `-f`, r1619 baseline). Missing is FLAIR guided on the splits used
here (the r1622 FLAIR guided run used the gene split). Added: `flair correct -f GUIDE` + `flair collapse -f GUIDE
--trust_ends` on the existing alignments (BAM -> BED12 by `flair_bam2bed.py`, no realignment), the r1622 recipe, with
`t.guide.gtf` (gorilla) and `ht.guide.gtf` (A119b). FLAIR 3.0.1, the locally patched install (`ssPrep.py` skips reads whose
junction conversion fails). Scored in the guided block with the same columns. Outputs `r1646_flairguided/`; wall limits
6 h gorilla / 24 h A119b, a run that hits its limit is reported as not finished.

## Outcome, Amendment 1, gorilla OR6737 (job 4467064, `r1646_flairguided/g.fl/`, `g.score.txt`; correct 49 s + collapse 1,165 s)
| tool (guided) | guide chains found | hidden found | chains outside the guide | precision outside the guide |
|---|---|---|---|---|
| FLAIR guided | 20,528 | 4,646 | 40,159 | 0.116 |
FLAIR guided writes 88,157 isoforms (60,687 multi-exon chains) against 138,125 chains de novo. Among guided tools it finds
about as many hidden chains as IsoQuant (4,761) and StringTie -G (4,538), with the most chains outside the guide and the
lowest precision there.

## Outcome, human A119b (2026-11-12, job 4466748, `r1646_isoquant/a.*`; IsoQuant 2,479 s de novo, 2,824 s de novo pt, 19,724 s guided, 20,710 s guided pt, of which about 4.5 h is the guide database build)
De novo block:
| tool, point | chains | correct | precision |
|---|---|---|---|
| ours, default | 219,517 | 38,715 | 0.176 |
| StringTie, default | 211,865 | 31,781 | 0.150 |
| FLAIR | 626,983 | 39,230 | 0.063 |
| bambu, NDR 1 / default | 41,745 / 9,599 | 12,793 / 6,680 | 0.306 / 0.696 |
| IsoQuant pt | 84,539 | 18,903 | 0.224 |
| IsoQuant default | 109,919 | 20,180 | 0.184 |
| IsoQuant pt, count >= 3 / 6 / 11 / 21 | 63,114 / 42,339 / 30,660 / 21,191 | 17,670 / 15,123 / 12,811 / 10,337 | 0.280 / 0.357 / 0.418 / 0.488 |
- Chains only one tool reports (and correct): ours 60,038 (1,953), StringTie 78,868 (2,319), FLAIR 462,779 (6,245), bambu 3,867 (73),
  IsoQuant pt 2,610 (83), IsoQuant default 4,444 (50). 17,945 of IsoQuant pt's 18,903 correct chains are also ours.

Guided block (54,491 hidden chains; FLAIR rows from job 4467065, `r1646_flairguided/a.fl/`, `a.score.txt`; bed 1,215 s + correct 153 s + collapse 4,801 s):
| tool | guide chains found | hidden found | chains outside the guide | precision outside the guide |
|---|---|---|---|---|
| ours guided (G1d) | 40,084 | 7,963 | 57,687 | 0.138 |
| StringTie -G | 59,668 | 5,334 | 42,238 | 0.126 |
| bambu guided, NDR 1 / recommended | not comparable | 3,630 / 1,108 | 29,630 / 1,690 | 0.123 / 0.656 |
| IsoQuant guided pt | 26,098 | 5,285 | 69,083 | 0.077 |
| IsoQuant guided default | 35,523 | 5,295 | 92,465 | 0.057 |
| FLAIR guided | 33,873 | 8,787 | 138,907 | 0.063 |

For reference only, de novo outputs scored on the same split (not guided results): ours 11,833 hidden / 192,635 outside the guide / 0.061; FLAIR 12,142 / 599,894 / 0.020.

## De novo precision at IsoQuant pt's correct counts (both samples; post hoc; `bench/famcmp/matched_iq.py`: each other curve at its strictest cut reaching the count, correct count at that cut in brackets; ours by reads, StringTie by cov, bambu by NDR)
| sample | correct chains (IsoQuant count cut) | IsoQuant pt | ours | StringTie | bambu |
|---|---|---|---|---|---|
| gorilla | 16,677 (all) | 0.450 | 0.538 (16,965) | 0.407 (18,842) | beyond its ceiling |
| gorilla | 15,244 (>= 3) | 0.511 | 0.563 (15,929) | 0.495 (15,644) | 0.495 |
| gorilla | 12,234 (>= 6) | 0.602 | 0.652 (12,452) | 0.582 | 0.634 |
| gorilla | 9,704 (>= 11) | 0.673 | 0.724 (9,934) | 0.630 (10,281) | 0.729 |
| gorilla | 7,273 (>= 21) | 0.740 | 0.785 (7,441) | 0.716 | 0.803 |
| A119b | 18,903 (all) | 0.224 | 0.395 (19,428) | 0.230 (21,409) | beyond its ceiling |
| A119b | 17,670 (>= 3) | 0.280 | 0.414 (18,019) | 0.268 (18,730) | beyond its ceiling |
| A119b | 15,123 (>= 6) | 0.357 | 0.454 (15,286) | 0.309 (15,246) | beyond its ceiling |
| A119b | 12,811 (>= 11) | 0.418 | 0.490 (13,027) | 0.333 (13,013) | 0.306 at 12,793 (its ceiling) |
| A119b | 10,337 (>= 21) | 0.488 | 0.543 (10,345) | 0.370 | 0.455 |

**Reading (descriptive; no adoption decision):**
- De novo: IsoQuant reaches further than bambu (16.7k / 18.9k correct against 15.4k / 12.8k) but only 0.65 / 0.49 of our
  correct chains. At every IsoQuant cut on both samples our read-count curve is more precise at an equal or larger
  correct count (+0.05 to +0.09 gorilla, +0.06 to +0.17 A119b). IsoQuant is slightly above StringTie on gorilla (+0.02 to
  +0.04) and above it on A119b except at its loosest point. On gorilla at the strict end (>= 11 reads and tighter) bambu
  is above IsoQuant and level with ours.
- Guided: no tool recovers hidden transcripts at high precision except bambu at its recommended cut, which finds the
  fewest. Ours guided has the best precision among the open settings on both samples and the most hidden chains on
  gorilla; on A119b FLAIR guided finds more (8,787 vs 7,963) with 2.4 times the chains outside the guide at half the
  precision. IsoQuant guided finds about as many hidden chains as StringTie -G with more chains outside the guide.
- On "makes up transcripts" (post hoc): of the de novo chains not in RefSeq, another tool reports the same chain for
  bambu 0.83 / 0.84 (gorilla / A119b) and IsoQuant pt 0.92 / 0.83, against ours 0.81 / 0.66, StringTie 0.65 / 0.56,
  FLAIR 0.30 / 0.22. Bambu's and IsoQuant's non-RefSeq chains are the most corroborated and the fewest; their private
  chains (0.6k-4.4k, 1-4% in RefSeq) are far fewer than any other tool's. Read support of private chains was not checked.
- Guided bambu does carry annotation without read evidence, by design: it writes all 65,846 / 120,728 guide chains while
  24,885 / 39,960 have a full-chain primary alignment. Those are labelled guide transcripts, not novel calls.
- On "leaves out new isoforms": at its recommended cut bambu guided finds 2,420 / 1,108 hidden chains against ours
  5,524 / 7,963; opened to NDR 1 it finds 4,093 / 3,630. De novo it stops at 0.60 / 0.33 of our correct chains.
- Not run: IsoQuant/FLAIR curves in the guided block, relaxed bambu discovery options, any second truth besides RefSeq.
