# STEP 0: ARE A READ'S TIED PLACES A PROPERTY OF ITS LOCUS? (register 1943, 2026-11-25)

Why: user goal of 2026-11-25 ("explore how to simplify the pipeline further using variation graphs ... graphs might
greatly simplify both assembly and assignment"). Every earlier graph row asked whether a graph PLACES reads better
(rows 1532, 1533, 1544, 1552, 1649, 1777: no; path compatibility = the exact k-mer rule; exact ties stay ties). This step
asks a different thing: whether the graph can REPLACE the bookkeeping the pipeline builds around minimap2's secondary
alignments. No accuracy claim is tested here.

## What the pipeline does today for multi-copy reads
1. minimap2 `-N 50 -p 0.1`: every read is aligned to up to 50 further places (OR6737: 6.27 M of 10.67 M records are secondary).
2. A score table over the BAM (best `ms` / AS per read).
3. Two bands on those records: the seeding pool (AS >= 0.98 x best) and the near tie (`ms` within 10).
4. A tie table (read, locus, gap), the EM labels, the sequence test, the second (floor-1) assembly: all read 1-3.

## The graph form
A family graph is the alignment of the copies to each other. If a read sits on locus L, the places it ties with are the
places L's own sequence aligns to nearly as well, and the difference between the two alignments of the read is the
difference between the two copies under the read (the bubbles it crosses). Then: map reads once (primary only), assemble,
align each TRANSCRIPT to the genome once (68,847 sequences instead of 4.4 M reads x 50), and read every tie off that.
Steps 1-3 become one small table, and "near tie" needs no band on a read's scores: it is the copies' difference under
the read.

## Measurement (OR6737, development sample; everything cached; no default touched)
- Loci: `r1920_mono_from_primaries/real_or/p.gtf` (assembly from primary alignments only, so no secondary enters).
- Edges: spliced transcript sequences (gffread) aligned to GGO.fa, `minimap2 -ax splice:hq -uf --eqx -N 50 -p 0.1`.
- Tied reads: as `tools/tie_table.py` (a second record within 10 `ms` of the best). A read's transcript T = the `p.gtf`
  transcript its primary alignment shares most exon bases with (same strand). The read's span is carried through T's own
  alignment into T's coordinates and through each other alignment H of T back to the genome.
- An in-band secondary record of the read is RECOVERED when a projected H overlaps >= 50% of its aligned bases with the
  matching relative strand.
- Predicted gap of H for the read = cost of H's differences inside the read's span (5 per mismatch, min(6 + n, 24) per
  gap, + n for inserted read bases, + 1 per base of the span H does not cover). Predicted tie = gap <= 10.
- Untied control: every 50th read (crc32 of the name) with no record in the band and a primary on a `p.gtf` transcript.

## Stop rule (fixed before the run)
PROCEED to a registered row (build the tie table from the projection, run stage `o2` on it, compare labels, copies and
sequences with the cached run) only if all four hold:
- (a) reach: >= 0.70 of tied reads have their primary on a `p.gtf` transcript;
- (b) recall: >= 0.90 of the in-band secondary records of those reads are recovered;
- (c) band: >= 0.90 of the recovered records get a predicted gap <= 10;
- (d) false ties: untied control reads with a predicted tie, scaled x 50, <= 0.20 x the number of tied reads in (a).
Else STOP and say which part fails.
Not judged, reported: exact-tie agreement (predicted 0 against actual 0); recall with any transcript of the locus
instead of the best one; minimap2 wall and output size on 200,000 reads with `-N 50 -p 0.1` against `--secondary=no`.

Prediction (mine): (a) holds (~0.75); (b) holds narrowly (0.85-0.95); (c) fails or is close (0.80-0.90: the read's own
errors and alleles are not in the transcript); (d) uncertain, the weakest point: a transcript-wide alignment can call a
tie where the read's 10-`ms` band has none. Overall: STOP on (c) or (d), with (b) showing the places themselves are
derivable.

Outputs: `/anvme/workspace/v138da13-rustle/r1943_graph_ties_step0/`. Script `bench/o2_families/graph_ties_step0.py`.

## RESULT (2026-11-25, jobs 4504539 and 4504549; `r1943_graph_ties_step0/step0.registered.summary.txt`, `posthoc.summary.txt`)
| bar | value | verdict |
|---|---|---|
| (a) reach, tied reads with a usable transcript | 12,714 of 16,275 (0.781) | holds |
| (b) recall, in-band secondary records recovered | 15,393 of 16,112 (0.955) | holds |
| (c) band, recovered records with predicted gap <= 10 | 9,398 of 15,393 (0.611) | FAILS (bar 0.90) |
| (d) false ties, control x 50 against tied reads | 1,350 against 12,714 (0.11) | holds |

**STOP, on (c).** Prediction held (STOP on the band); the band is lower than I expected (0.80-0.90).
- The places are a property of the locus: 95.5% of the in-band secondary records lie where the read's own transcript aligns
  a second time; a second transcript of the locus adds 9 records. Aligning the transcripts took 67 s.
- The closeness is not. For 5,995 recovered records the copies differ by more than 10 under the read although the read's two
  alignments are within 10 (1,268 of them exact ties). Exact ties: 4,201 of 5,776 predicted exact.
- Post hoc (after the numbers, no bar): counting a copy difference FOR the copy when the read itself differs from its own
  place at that base raises the band to 0.778 and the false ties to 0.24 (over the bar). The rest is in how each alignment
  ends and clips, which only aligning the read at the copy gives.
- Untied reads: of 94,487 out-of-band secondary records matched to a projection, 7 are predicted in band.
- 3,561 tied reads (0.22) have their primary alignment on no transcript of the two-read assembly.
- Added after the rule, descriptive: of the 445 copies the multi-copy stage rescued in the cached run, 370 are a secondary
  place of a transcript of the primaries-only assembly, 72 are at no transcript alignment.
- Not judged, mapping cost on 197,787 reads, 32 threads, same index: with secondaries 188 s wall / 4,180 s CPU / 431 MB BAM
  (479,453 records); primary only 120 s / 2,324 s / 258 MB (199,365 records).

Reading: a graph of the copies can say WHERE the copies of a locus are without asking the aligner for secondaries of every
read; it cannot say which reads are near-tied there. Replacing the secondaries would need the reads re-aligned at the
projected places, i.e. a second alignment step in place of the first one's secondaries: not simpler. One sample, development.
- Size of the alternative (re-align reads only at the projected places): 15,197 of 68,847 transcripts (0.22) have a second place
  at `-p 0.1`, and they carry 958,848 of 3,071,258 assembled reads (0.31); tied reads sit on 1,174 transcripts, but which ones
  is known only after the alignment. A targeted second alignment would touch a third of the reads.
- Assembly side, from earlier rows (not re-run): 8% of isoforms are shared by several members (row 1822, arm F) and borrowed
  structure fails (rows 1822 B1, FAMASM), so one assembly per family graph has little to merge.

## SECOND SAMPLE, A119b on CHM13 (2026-11-25, job 4504560, same script and bars; `r1943_graph_ties_step0/a119b/`)
| bar | OR6737 | A119b |
|---|---|---|
| (a) reach | 0.781 | 240,491 of 316,163 (0.761) holds |
| (b) recall of in-band secondary records | 0.955 | 585,807 of 977,313 (0.599) FAILS |
| (c) predicted gap <= 10 | 0.611 | 0.548 FAILS (post hoc read-aware 0.766) |
| (d) false ties | 0.11 | 0.20, at the bar (post hoc read-aware 0.30) |

STOP on A119b too, and on recall as well: this sample has 3.1 in-band secondary records per tied read (OR6737 1.3), i.e.
large arrays, where one alignment of the transcript with at most 50 secondaries does not list every copy. So "the places
are a property of the locus" held on OR6737 only; in large families even the places are incomplete. Rescued copies at a
transcript's secondary place: 413 of 539. Transcripts to genome: 579,586 records in about 5 minutes; projection 36 minutes
(single process).
