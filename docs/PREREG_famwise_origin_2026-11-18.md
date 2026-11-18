# PREREG r1796 / r1794 / r1795: the family-wise layer against the origin of the reads (2026-11-18)

Written before any count of these rows. Workspaces `r1796_famwise_origin` (row a; first claimed as r1793, which another session's row took the same hour), `r1794_famwise_oracle` (row b and the
re-toss), `r1795_famwise_planted` (row c). Design: `docs/superpowers/specs/2026-11-18-family-wise-assembly-design.md`
(sections "Scores" and "Measurement"); plan `docs/superpowers/plans/2026-11-18-famwise-evaluation.md`. Inputs: the
catalogue-free truth tables of r1785 (`r1785_truth_repair/{kb,hg}/reads.tsv`) and r1787 (orangutan), the tie split of
r1791, the selected-reads BAMs of the truth runs (they carry `ms`), the default assemblies of 2026-11-16 and the StringTie,
FLAIR and IsoQuant GTFs on disk (trusted baselines), `bench/famsim`. The stage `famwise` is commit 59639df4 and is not
changed by these rows. The four bars were drafted in the design and accepted by the user on 2026-11-18 ("yes these are
good"), before any run below.

## Why
Advisor, as relayed on 2026-11-18: the main problem of transcriptomics is the true provenance of reads with nearly tied
primary and secondary alignments, and an assembly that works on families instead of loci should gain several-fold. The
user's results say near ties are few and their effect on an assembly is small. Three questions separate the two views:
a. Is the layer's statement about a molecule (the set of its tied records) closer to where the read came from than the
   primary alignment alone, and do its statuses sort transcripts by how well their copy is supported?
b. How much could ANY provenance method change an assembly, on real reads, if it were perfect?
c. Where, in terms of differences between copies and of their expression, is a several-fold change possible at all?

## Rule, fixed now
- **Samples.** KB3781 develops (every script is debugged on it and then frozen). The judging reads of rows a and b are a
  quarter of the HG002 reads nobody has looked at: the truth job's selection with `NT_SELECT_PCT=25 NT_SELECT_OFFSET=25`
  (hash residues 25 to 49; all reads seen so far are residues 0 to 24), opened once. Orangutan AG06213 is reported
  beside, no bar.
- **Scores.** Every band, stratum and count is computed on minimap2's `ms` and on AS, both printed; the bars are read on
  `ms`. The truth tables define their strata on AS; the scorer re-reads both scores of each record from the BAM.
- **Truth kinds.** `one` = status `decidable` (one origin copy); `several` = status `tie` with class SEVERAL COPIES in
  the tie split; `one_place` = tie of class ONE PLACE; `tie_unsplit` = tie without a split row; `outside`; `unjudged` =
  any other status or an origin beyond the 12 kept loci. Judged = `one` and `several`.
- **Statement of a molecule.** The primary record alone; or the primary plus every record whose score is within W of the
  molecule's best, W in 0 (the exact band), 5, 10, 20. A statement holding a record on a place the truth never judged
  has no verdict: counted, never scored. `equals` = the statement's loci are exactly the truth's origin set; `sound` =
  the origin set and the statement share a locus.
- **Near-tied on a score** = the molecule's two best places are within 10 on that score; records on one contig and
  strand with overlapping spans are one place (the truth's own merge), so a second alignment of the same place does
  not make a near tie.
- **Claims.** A claim is a multi-exon transcript of an arm's GTF, keyed (contig, strand, intron chain). Its supporters
  are judged molecules with a record within 10 `ms` of their best (or the primary) carrying exactly that key: one rule
  for every tool, no tool's own read map. A supporting record on a place the truth never judged does not count. TRUE: a supporter's origin set is exactly that record's locus; else SHARED: a
  supporter's origin set holds it; else PHANTOM; no supporter: UNJUDGED. One verdict per key, whatever arm emits it. Arms: the default run's `families.gtf`
  (`rustle`, split by the layer's status from the stage run on the same selected reads), StringTie, FLAIR, IsoQuant.
- **Three inputs of row b**, made from one alignment: P = primary alignments only; D = the default (all records, default
  seeding); O = each read with a truth entry (kinds `one`, `several`, `tie_unsplit`) kept only at its copy or copies of
  origin, the others at their primary alignment. Assemblers: Rustle (`assemble`) and StringTie (`-L`). Regions = merged
  truth spans (origin places and the primary's place) of every read that is tied or whose primary is off its origin and
  that has a record at its origin. R = RefSeq-exact intron chains in the regions, O over P. Rustle is read on
  `PREFIX.gtf` (the assembler's output; `PREFIX.final.gtf`, after the chain score filter, is printed beside, no bar).
- **Row c**, `bench/famsim`: one synthetic 6-exon template (2,000 bp), two copies each with one exon-skipping isoform;
  differences per transcript k in 0, 1, 2, 4, 8, 20 (exonic SNPs on copy B) x expression A / B 200 / 200, 200 / 20,
  200 / 0 x full-length reads or 3' fragments (36 cells), plus an array of 8 identical copies at 25 reads each, both
  read kinds (2 cells); 5 seeds; inputs P, D, O as above with the planted copy as origin. `recovered` = planted
  (isoform, copy) pairs with expression > 0 whose exact chain is emitted on the copy. A cell is several-fold when the
  median over the seeds of recovered O over recovered P is >= 2, or P recovers none and O some.
- Intervals are Wilson 95%.

## Readings, fixed now
- **B1 (row a, held-out HG002 quarter).** Molecules near-tied on `ms`, kinds `one` and `several` pooled, bands on `ms`.
  PASS when the exact-band statement equals the truth set more often than the primary alone with disjoint intervals,
  AND, among the families with at least 20 judged near-tied molecules, in more than half the exact band's share is
  strictly the higher one. Molecules without a verdict for either statement are left out of both shares and
  counted. No family table, or no family of that size: undecidable.
- **B2 (row a, same reads).** (i) PASS when the share of PHANTOM + SHARED among claims of status `tied_only` is higher
  than among `certified` ones with disjoint intervals; fewer than 30 judged `tied_only` claims = undecidable, said so.
  (ii) PASS when `certified` + `best_only` claims hold at least 90% of the `rustle` arm's TRUE claims. The PHANTOM share
  of certified claims is printed as a check, not a bar.
- **B3 (row b, regions).** The advisor's expectation holds if R >= 2 for either assembler; the user's if R <= 1.1 for
  both; in between the number is reported and neither is claimed. Read on the held-out HG002 truth; KB3781 beside.
  Also printed: chains and RefSeq-exact chains per input in the regions and genome-wide, and the chains of D that O
  lacks with how many of them are RefSeq-exact.
- **B4 (row c).** PASS when several-fold cells lie only at k <= 2 or in the array row, and at least one exists. No
  several-fold cell anywhere: the simulator shows nothing and no claim is made. A several-fold cell at k >= 4: FAIL
  (provenance matters for diverged copies too, in the advisor's favour).
- **Descriptive, no bar:** the within-W envelope with both bounds on the reads the truth calls `outside`; the range
  table (copies whose true count lies in [lo, hi]); the baselines' claim tables; molecules near-tied on one score only;
  the re-toss (the selected KB3781 reads realigned with `--seed` 11, the control that must reproduce the records, and
  12, 13, 14: claims not in every seed per arm, rows of the layer that differ).

## Prediction
B1 PASS (a replication: on the first HG002 quarter, AS proxy, the exact band equalled the truth for 0.90 of `one` and
0.82 of `several` molecules, the primary for 0.92 and 0). B2 (i) PASS or undecidable on KB3781, PASS on HG002; (ii)
PASS. B3: R <= 1.1 for both assemblers. B4 PASS, with the several-fold cells at k = 0 (silent copy or array).

## Seen already
The AS-proxy table of the design (three samples; not on `ms`, not on the held-out quarter). On KB3781, O differs from P
for 937 of 100,810 selected reads. At tied loci the default seeding band over primary-only gives x1.10 to x1.22
RefSeq-exact chains. The whole-sample output of the stage (design, "Stage output"), which holds no truth.

## Limits stated beforehand
- Annotation families only (no truth read touches CD99); near ties complete, clear-lead reads a sample (30,000 per run);
  at most 12 loci per read; the truth keeps placements at AS >= 0.8 x best and defines its strata on AS, the layer on `ms`.
- HG002 is a 25% read sample on another individual's reference. "Judged" requires the origin among the kept placements:
  reads the truth calls `outside` are possible failures of the layer's assumption and give the upper bound on them.
- The sampled truth over-weights near ties, in the advisor's favour for B3. O is an oracle only for the reads the truth
  judges; every other read stays at its primary alignment in all three inputs.
- Row c is one template and one exon-skipping isoform; its map says where a several-fold change is possible, not how
  often real families sit there (the real two-member pairs are drawn on it from r1773).
- Any change to a scorer after its KB3781 run is listed here as an amendment before the held-out reads are opened.

## Amendments
- 2026-11-18, before the held-out reads were opened. Development runs on reads already seen (KB3781, the first HG002
  quarter, orangutan; workspace `r1796_famwise_origin/{kb,hg_seen,pab}`) reproduced the known counts (KB3781, AS: 45
  one-origin and 786 tied near ties, primary equal for 40 of 45). One change to the scorer after its first KB3781 run:
  a per-claim table `claims_tx.tsv` (no number of another table changes) and a faster range loop. Scorer frozen at
  `records.py` md5 5775699ad0fd, `origin_score.py` md5 3e3340e43052.
- Seen on the already-seen HG002 quarter with the frozen scorer (so the held-out run is a replication of these, not a
  first look): B1 0.82 against 0.12 and 8 of 15 families; B2 (i) 0.90 against 0.02, (ii) 0.97.
- 2026-11-18, written while the held-out job (4482393) ran and before any of its output was read; no script changes.
  The design's preview table is not reproduced in two cells of its one-origin rows: on the seen HG002 quarter the
  registered exact-band statement equals the truth for 400 of 600 one-origin molecules on AS (359 of 440 on `ms`), not
  0.90, and the within-10 statement for none, not 0.15 (a molecule whose two best places are within 10 always states
  two places). The probe behind that table used another statement, which I could not recover; the tie rows, the
  near-tie counts and the primary's shares are reproduced. The Prediction's "0.90 of `one`" was therefore wrong before
  the run. Consequence for reading B1, said now: the pooled share is carried by `several` molecules, where a single
  primary alignment cannot equal a set of two or more copies by construction; on `one` molecules the primary alone is
  the better statement (396 against 359 of 440 on the seen quarter). Both kinds are reported separately beside the bar.
- 2026-11-18, row b, found by the independent reviewer after the KB3781 and seen-quarter oracle runs and before the
  held-out one printed anything: `origin_bam.py` wrote the regions from every read O restricts or moves, which is not
  the Rule above (reads that are tied or whose primary is off its origin). It took in 20,650 one-origin KB3781 reads
  already at home (O = P for them) and left out 692 ties whose records all lie on their origin copies; the regions held
  95% of all chains. Fixed to the Rule's sentence (test `test_regions_are_the_tied_or_moved_reads`, failing first); the
  regenerated files equal the reviewer's independent ones (109 and 130 intervals). The counts on the wide regions are
  kept as `regions.shipped_wide.bed` and reported beside. No assembly was rerun; the rule did not change.
- Same review, not changed: `records.locus_index` takes the first truth locus a record overlaps, the truth merges
  greedily on growing spans; replaying the truth's merge moves 3 of 3,379 molecules of B1 on the seen quarter and no
  claim verdict, and 64 of 24,641 (KB3781) / 112 of 21,797 (seen quarter) oracle reads keep another record set. The
  scorer stays frozen for the held-out reads (they were open by then); stated as a limit of the result.

## Result (2026-11-18)
Predictions: B1, B2 and B3 as predicted; B4 not (no several-fold cell exists, the prediction named cells at k = 0).

**Row a (r1796), held-out HG002 quarter.** 136,243 selected reads (hash residues 25 to 49, none shared with the seen
quarter), 42,572 truth rows, ties split into 2,875 SEVERAL COPIES and 1,223 ONE PLACE, none unsplit. Bands on `ms`.
- **B1 PASS.** Molecules near-tied on `ms`, kinds `one` + `several`: the exact band equals the truth set for 2,662 of
  3,259 (0.817, 0.803 to 0.830), the primary alone for 341 of 3,259 (0.105, 0.095 to 0.116); the exact band is the
  higher one in 9 of the 15 families with 20 or more such molecules; no molecule without a verdict. By kind: `several`
  (2,873) primary 0 by construction, exact band 2,359 (0.82), within 10 2,761 (0.96); `one` (386) primary 341 (0.88),
  exact band 303 (0.78). The true copy is in the statement for 341 (primary), 366 (exact band) and 384 (within 10) of
  the 386. On AS: 498 `one`, primary 451, exact band 326; 2,874 `several`, exact band 2,359.
- **B2 (i) PASS.** PHANTOM + SHARED among judged claims: `tied_only` 113 of 135 (0.837, 0.766 to 0.890), `certified` 18
  of 1,316 (0.014, 0.009 to 0.022). **(ii) PASS:** `certified` + `best_only` hold 1,311 of the arm's 1,362 TRUE claims
  (0.963). Check: 12 of 1,316 judged certified claims are PHANTOM.
- Beside, no bar. Seen HG002 quarter: B1 2,772 of 3,379 against 396; B2 (i) 114 of 127 against 27 of 1,329; (ii) 0.970.
  KB3781 and orangutan: B1 undecidable (13 and 4 judged near-tied molecules, no family of 20), B2 (i) undecidable (no
  judged `tied_only` claim), (ii) 0.984 and 0.986.
- Descriptive, held-out quarter. Judged claims that are not TRUE: the default Rustle run 222 of 1,584 (0.14), its
  `certified` + `best_only` part 48 of 1,359 (0.035), StringTie 70 of 1,053 (0.066), FLAIR 152 of 1,989 (0.076),
  IsoQuant 31 of 541 (0.057). Range: the true count lies in [lo, hi] for 1,375 of 1,384 nodes (exact band) and 1,426 of
  1,433 (wide). 227 molecules near-tied on `ms` are `outside` (the upper bound on failures of the layer's assumption).
  Near-tied on one score only: 957 on AS, 80 on `ms`.

**Row b (r1794).** Regions as registered (amendment above); RefSeq-exact chains in the regions, P -> O.
- **B3: R <= 1.1 for both assemblers on the held-out HG002 truth, the user's expectation holds.** 138 regions. Rustle
  141 -> 147 (R = 1.04); StringTie 113 -> 113 (R = 1.00). Rustle after the chain score filter 130 -> 125.
- Beside. KB3781 (109 regions): Rustle 42 -> 47 (1.12, five chains), after the filter 41 -> 41, StringTie 31 -> 31. Seen
  HG002 quarter (130 regions): Rustle 142 -> 146 (1.03), after the filter 127 -> 126, StringTie 110 -> 112 (1.02). On
  the wide regions first written (253 and 333 intervals): 1.00 / 0.99 / 1.00 in both samples.
- The default input D gives more RefSeq-exact chains in the regions than the oracle: 154 against 147 (held-out), 153
  against 146 (seen), 49 against 47 (KB3781). Of the 132 chains of D that O lacks on the held-out quarter 8 are
  RefSeq-exact.
- Limits: small counts; D carries secondary records only for reads the truth judges (the default on every record of
  the KB3781 selection has 1,569 chains, D 1,496); StringTie's counts are the same on P and D, and its O differs from P
  by the 62 reads whose primary is off their origin; the `locus_index` limit above.
- **Re-toss** (KB3781 selected reads; seed 11 reproduces all 732,401 records). Seeds 12 / 13 / 14: the primary changes
  for 68 / 72 / 78 of 100,810 reads. Chains not in every seed: Rustle 0 of 1,569, after the chain score filter 2 of
  1,120, StringTie 22 of 1,372. Layer rows that differ from seed 11: at most 2 of 1,568 per band.

**Row c (r1795).** 38 cells x 5 seeds (five seed-1 cells had been seen while the script was debugged).
- **B4: no several-fold cell, so no claim is made.** The registered quantity is saturated: primary-only already
  recovers every planted pair with reads in 32 of the 38 cells; the largest median O over P is 12 / 11 (the array with
  fragments). Where P misses a pair (the 10-fold weaker copy read by fragments; the array with fragments), the default D
  recovers it (medians 4 of 4 up to k = 8, 16 of 16) and the oracle does not (3, 12).
- Descriptive, not a registered reading (`map/status_by_copy.tsv`, `map/planted_map.png`). At a copy with no reads the
  default run emits chains up to 8 differences per transcript (full-length reads: 2.0 per run up to k = 4, 0.6 at
  k = 8; fragments: 5.0 up to k = 4, 3.0 at k = 8), primary-only at k = 0 and, with fragments, falling from 4.2 to 0.6
  by k = 4, the oracle never. The layer's status
  on the default run, cells with k >= 1: none of the 124 transcripts at the silent copy is `certified` or `best_only`;
  165 of the 175 at the expressed copy are. At k = 0 all are `tied_only`. Of the 115 real two-member pairs with an
  identity (gorilla and orangutan, r1773), 95 lie beyond 40 differences per 2 kb and 10 at two or fewer.

**Reading.** On real reads a perfect knowledge of the origin of the reads with primary or secondary alignments would not
change what is assembled by a factor near two: 1.04 and 1.00 on the held-out quarter. What the layer adds is the other
side: it names the set of copies a molecule fits (right for 0.82 of the molecules that fit several) and sorts the
transcripts, 0.84 of the `tied_only` ones being shared or false against 0.014 of the `certified` ones. Not established:
anything for families outside the annotation's, the size of the `outside` reads' effect, and B4.
