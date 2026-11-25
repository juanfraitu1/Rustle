# PREREG: the 10-AS filter as the baseline, and the classes that follow from it (row 1954, 2026-11-25)

The advisor (relayed by the user, 2026-11-25) asks that we start from this baseline: the correct copy of a read is the
one of its primary alignment when the primary is 10 AS above every other alignment of the read. This is the Eichler
lab's 10-point filter (Guitart 2024), adapted to a single reference; the paper's own procedure compares primaries across
haplotype assemblies and has not been run here. The user: see how many near-tied reads have a primary that is 10 AS
from all the others, and classify from there.

Row 1786 profiled the filter by AS lead against the haplotype truth. New here: the cross with `ms` (the score minimap2
ranks by, on which our near tie is defined), over whole samples, and what each SECONDARY alignment is worth.

## Definitions (per read, from the BAM's score table)
- **AS lead** = AS of the primary minus the highest AS of any other alignment. **Rule fires** when AS lead >= 10.
- **Near tie on `ms`** (the pipeline's definition): another alignment within 10 `ms` of the read's best, or the primary is
  not the best on `ms`.
- Cells: rule fires / abstains x near tie / clear on `ms`. Reads with a single alignment are counted apart.

## Parts
1. **Census, whole samples** (score tables of OR6737, KB3781, A119b, HG002; `bench/o2_truth/as10_census.py`): the four cells,
   the AS lead in bins (< 0, 0, 1-9, 10-19, 20-49, >= 50), and for the cell "rule fires AND near tie on `ms`" whether the
   primary is unspliced and its rival spliced (the intron effect).
2. **Truth, real reads** (family reads of HG002 on CHM13 and of KB3781 on its own assembly; the repaired truth tables of
   row 1785, weighted to the sample; `as10_truth.py`): per cell, how often the primary is at the read's origin; where the
   origin is when it is not the primary; reads whose origin the truth cannot give, by reason.
3. **Truth, planted reads** (row 1952's sets, seeds 1967-1969): the same cells with the copy of origin known for every read,
   counting a copy identical in sequence to the true one apart.
4. **What a secondary alignment is worth** (parts 2 and 3): the share of secondary alignments that sit at the read's
   origin, by their AS distance below the primary (above it, 0, 1-9, 10-49, >= 50) and by whether they are within 10 `ms`.

## Claims
- **B1, the baseline holds where it speaks:** on HG002 and on KB3781, among reads where the rule fires and `ms` is clear,
  the primary is at the origin in at least 99% of decidable reads.
- **B2, it does not cover our near ties:** on each of the four samples the rule fires on fewer than 30% of the reads that
  are near-tied on `ms`.
- **B3, where it fires on a near tie it is less reliable:** on HG002, among reads where the rule fires and `ms` is a near
  tie, the primary is at the origin less often than in B1's cell, by at least 3 points.
- **B4, far secondaries carry no origin:** on HG002 and KB3781 fewer than 1% of the secondary alignments that are 10 or
  more AS below the primary sit at the read's origin.

## Predictions (mine)
B1 holds (row 1786: under 1% wrong at leads above 10). B2 holds, with 10-25%. B3: uncertain; I expect it to hold on HG002
through the spliced-against-intronless class (AS lead 24 per intron at an exact `ms` tie). B4 holds.

## Result (job 4504961)

### 1. Census, whole samples (reads with several alignments)
| | OR6737 | KB3781 | A119b | HG002 |
|---|---|---|---|---|
| reads with several alignments | 1,345,946 | 4,177,965 | 5,795,139 | 15,837,678 |
| rule fires, `ms` clear | 94.9% | 88.7% | 91.6% | 91.4% |
| rule fires, `ms` near tie | 0.14% | 0.12% | 0.84% | 0.11% |
| rule abstains, `ms` clear | 3.9% | 10.0% | 2.9% | 7.3% |
| rule abstains, `ms` near tie | 1.1% | 1.2% | 4.6% | 1.1% |
| **near-tied on `ms`** | 16,245 | 53,712 | 316,155 | 198,767 |
| **of them, the rule fires (primary >= 10 AS above every other)** | 1,903 (11.7%) | 5,133 (9.6%) | 48,743 (15.4%) | 17,267 (8.7%) |

- The near-tied reads by AS lead (OR6737 / KB3781 / A119b / HG002): primary BELOW a rival on AS 15% / 18% / 25% / 26%;
  exact AS tie 33% / 47% / 29% / 40%; lead 1-9 40% / 26% / 31% / 25%; lead 10-19 9% / 4% / 13% / 7%; lead 20 or more 3% / 5% / 3% / 1%.
- Where the rule fires on a near tie, the pair is spliced against spliced in 33-65% and unspliced against unspliced in
  27-42%; the intron case (unspliced primary, spliced rival) is 4-9% (29% on KB3781). Most of these leads are 10-19.
- The cell "rule abstains, `ms` clear" is almost all reads whose primary has a LOWER AS than a rival (OR6737 48,429 of 51,844;
  HG002 1,032,470 of 1,163,865): the aligner ranked them on `ms`, and AS disagrees.

### 2. Against the truth, real family reads (weighted to the sample)
| HG002 on CHM13 | reads | primary at the origin | not | several copies of the individual fit equally | truth cannot judge | primary right, of the judged |
|---|---|---|---|---|---|---|
| rule fires, `ms` clear | 114,816 | 55,967 | 63 | 1,764 | 57,021 | 99.89% |
| rule fires, `ms` near tie | 997 | 181 | 7 | 415 | 395 | 96.3% |
| rule abstains, `ms` clear | 3,911 | 1,626 | 5 | 340 | 1,940 | 99.69% |
| rule abstains, `ms` near tie | 6,399 | 291 | 43 | 3,897 | 2,168 | 87.2% |

| KB3781 on its own assembly | reads | at the origin | not | several copies fit equally | cannot judge | primary right, of the judged |
|---|---|---|---|---|---|---|
| rule fires, `ms` clear | 92,326 | 63,275 | 3 | 330 | 28,717 | 99.995% |
| rule fires, `ms` near tie | 113 | 3 | 0 | 45 | 64 | (3 of 3) |
| rule abstains, `ms` clear | 3,186 | 2,358 | 0 | 121 | 707 | 100% |
| rule abstains, `ms` near tie | 1,193 | 12 | 5 | 751 | 425 | (12 of 17) |

### 3. Planted reads (six sets; the copy of origin known for every read)
| | reads per set | primary at the true copy | at a copy identical in sequence | at another copy | on no planted copy |
|---|---|---|---|---|---|
| rule fires, `ms` clear | 1,380-2,828 | 99.7-100% | 0 | 0-6 | 0-7 |
| rule fires, `ms` near tie | 296-581 | 95.3-100% | 0 | 0-16 | 0-1 |
| rule abstains, `ms` near tie | 2,339-4,505 | 47-50% | 28-31% | 17-18% | 3-5% |

When the primary is not at the true copy, the true copy is at a secondary with EXACTLY the primary's AS in about 97% of
cases (2,252 of 2,308 in the first set), at a secondary with a higher AS in 1-2%, and 1-9 AS below in under 1%.

### 4. What a secondary alignment is worth (share that sits at the read's origin)
| AS distance below the primary | HG002, read near-tied on `ms` | HG002, `ms` clear | planted, alignment within 10 `ms` |
|---|---|---|---|
| above the primary | 4.6% (8 of 174) | 0.1% | 6-20% |
| 0 (exact AS tie) | 56% (19 of 34) | 0 of 9 | 16-20% |
| 1-9 | 5.6% (14 of 252) | 0.4% | 0.01-0.07% |
| 10-49 | 1.5% (8 of 520) | 0.3% | 0 |
| 50 or more | 0.07% | 0.02% | - |
(KB3781: 5 secondaries at the origin among 142 in near-tied reads, all at distance 0 or 1-9; 3 of 272,000 elsewhere.)

### Claims
- **B1 HOLDS.** Where the rule fires and `ms` is clear the primary is at the origin in 99.89% (HG002) and 99.995% (KB3781)
  of judged reads. The advisor's baseline is right where it speaks, and it speaks for about 90% of multi-aligned reads.
- **B2 HOLDS.** The rule fires on 11.7%, 9.6%, 15.4% and 8.7% of the reads near-tied on `ms`. For roughly nine in ten of
  the reads this project is about, the baseline says nothing.
- **B3 HOLDS on HG002** (96.3% against 99.89%, 3.6 points; 188 judged reads), not testable on KB3781 (3 judged reads).
  The stronger fact in that cell is not the error rate: in 42% of its reads (415 of 997) the individual's own haplotypes
  say several copies fit equally, yet the rule names one.
- **B4 HOLDS.** Secondaries 10 or more AS below the primary sit at the origin in 0.03% (HG002, 74 of 239,375) and 0.001%
  (KB3781).

### The classes this gives
1. **Settled by the baseline** (rule fires, `ms` clear; about 90% of multi-aligned reads): the primary is right, the
   secondaries carry no origin (they are evidence of homology only).
2. **Settled by `ms`, missed by the baseline** (rule abstains, `ms` clear; 3-10%): mostly reads whose primary is BELOW a
   rival on AS. The primary is right in 99.7-100%. The 10-AS filter would set these aside, or move them to the best-AS
   locus, which row 1786 found wrong 17-40% of the time.
3. **Named by the baseline, tied on `ms`** (0.1-0.8%): the primary is right in 96%, but in four in ten of these reads the
   individual's own copies fit equally.
4. **Open on both scores** (1-5%): the primary is at the origin in 87% of the judged reads on HG002 and at the exact copy in
   about half of the planted ones; in 61% of the real ones several copies of the individual fit equally. When the primary
   is wrong the origin is, almost always, a secondary with exactly the primary's AS.
So the secondaries that are worth something are the EXACT ties (and, far behind, those within 10 `ms` and under 10 AS);
every other secondary sits at the origin in well under 1%. That matches row 1951 (another session): a seeding pool of
exact `ms` ties alone raises planted isoform precision by 0.07-0.10 and loses annotated copies on three of four real
samples, a trade and not adopted.

### Is a primary with a lower AS than a secondary a mistake? (user's question; description, OR6737 score table + minimap2 2.30 source)
No. minimap2 does not use AS to choose the primary. It sorts a read's hits by the score of the best-scoring segment
(`dp_max`; the `ms` tag), after a re-ranking step when the two best are within 0.9 of each other (`mm_update_dp_max`,
`align.c:1022`), ties by a per-hit hash (`mm_hit_sort`, `hit.c:203`); going down that order, a hit overlapping a primary
on the read becomes its secondary (`mm_set_parent`, `hit.c:125`). AS is the full alignment score, which in splice mode
pays about 24 for every intron.
On OR6737 50,803 reads have a primary below its best rival on AS: in 88% the primary is spliced and the rival unspliced,
and the primary's `ms` is higher in 97% (equal in 1%, lower in 2%, the re-ranking step). These reads are most of the
class "rule abstains, `ms` clear", where the primary is at the origin in 99.7-100% of judged reads. The spliced alignment
at the gene won on bases; the intronless copy only looks better because AS does not charge it for introns.

### Can a locus that rests on tied primaries be made up, and does a mapping-quality filter fix it? (user's questions; descriptions, workspace `coin_toss_2026-11-25`)
- **The tie-break shown directly** (7,980 planted reads, same command): mapped twice, 0 primaries move; with another
  `--seed` 2,284 move (28.6%); with the reads renamed 2,308 (28.9%); every moved read has mapping quality 0; with
  `--no-hash-name` renaming moves none. Reads of seven pairs of identical copies land 343 on the true copy, 344 on the twin.
- **Every mapped read has a primary.** OR6737 and KB3781: 0 molecules without one; unmapped reads 5,506 and 959; 16,866
  and 26,582 reads reach the cap of 50 secondaries, so their full set of ties is not in the file.
- **Loci that rest on tied reads alone** (OR6737, from the provenance table): of 20,142 loci with at least one primary
  alignment carrying a transcript, 19,083 have no near-tied read, 514 some, and 545 (2.7%; 774 transcripts, 4,678 reads;
  116 with five or more reads) have ONLY near-tied reads behind them. The copy exists; that it is expressed is not shown.
- **Mapping quality 0 is the tie-break flag:** on OR6737 98.2% of exact `ms` ties have mapping quality 0 and 98.3% of the
  5,783 mapping-quality-0 primaries are exact ties. Near ties that are not exact have mapping quality 1-29 (92%).
- **Removing the mapping-quality-0 primaries** (planted, seeds 1967-1969; 42-44% of primaries; job 4505172):

| | StringTie | StringTie without them | Rustle primaries | Rustle primaries without them |
|---|---|---|---|---|
| silent copies with a transcript (of 63) | 23, 27, 25 | 0, 1, 0 | 23, 25, 24 | 0, 0, 0 |
| expressed copies found, all expressed (of 133) | 133, 133, 133 | 83, 83, 83 | 129, 129, 129 | 79, 79, 79 |
| isoform recall, all expressed | 0.669, 0.669, 0.657 | 0.397, 0.402, 0.393 | 0.820, 0.824, 0.816 | 0.452, 0.452, 0.452 |

  The made-up loci go, and with them every copy that has an identical twin. Three treatments of a tied read: keep the
  primary as placed (StringTie, FLAIR, bambu: silent copies modelled), remove mapping quality 0 (IsoQuant in part: tied
  copies lost), or keep the reads and label what rests on ties (Rustle).
