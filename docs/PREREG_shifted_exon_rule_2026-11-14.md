# PREREG row 1726 (2026-11-14): a graded drop rule for weaker shifted-site chains, derived on two samples, tested on HG002

User goal: "derive empirical rules to improve assembly focusing on shifted internal exons". Follows row 1725 (FAIL).
Written before any table of this row.

Seen (row 1725, OR6737 / A119b): the class S1 + S2 and its marginal strata. In RefSeq: reads ratio to the partner
< 0.1 0.094 / 0.118; frame-breaking shift 0.157 / 0.137 vs frame-keeping 0.386 / 0.242; <= 10 bp and ratio < 0.1
0.075 / 0.109. The oracle ceiling over the read-count curve is +0.022 / +0.007, so no rule on this class reaches
row 1677's +0.02 on A119b: **this row cannot produce an assembly-level precision rule and does not claim one.** Not
seen: any cross of frame x ratio x shift, the read-cut's annotated share at a matched count, anything on HG002.

Question: is there a rule on this class that removes chains which are annotated clearly less often than the chains a
plain read cut of the same size removes, without costing annotated shifted-site pairs? That is the smallest sense in
which a class rule "improves" the assembly beyond the read count.

## Rule family (GTF only; class S1 + S2 of `shift_class.py`, terminal-exon class S1t never touched)
Drop a weaker chain c with partner d when: frame in {breaking only, any} AND reads(c) / reads(d) < Q, Q in
{0.02, 0.05, 0.1, 0.2} AND shift in {<= 10 bp, <= 50 bp, any}. 24 cells. Chains are dropped whole; none is rewritten.

## Derivation (discovery; never pooled)
On OR6737, a cell qualifies when (a) it drops >= 300 chains, (b) dropped chains in RefSeq <= 0.5 x the RefSeq share
of the read cut at the same count (lowest `reads` first, file order within ties), (c) annotated pairs lost <= 0.02
(class chains in RefSeq whose partner is in RefSeq too; share dropped). Among qualifying cells take the one that
drops most chains. That cell must meet (a)-(c) on A119b as it stands ((a) at >= 600); if it does not, take the next
largest OR6737 cell that does; if none, the row ends FAIL with HG002 unopened. The choice is recorded before HG002
is opened.

## Test (held out: HG002 on CHM13, `r1692_trunc5_hg/hg.gtf`, RefSeq + CAT as on A119b)
The frozen cell, bars:
1. dropped chains in RefSeq <= 0.5 x the read cut's at the same count;
2. annotated pairs lost <= 0.02;
3. correct chains lost <= 0.5% of the assembly's;
4. >= 600 chains dropped.
Reported beside them: precision before / after and the read curve at the same count (expected gain small, see
above); share of dropped chains whose moved intron is annotated (RefSeq + CAT), as the second yardstick, no bar.
PASS = an empirical rule for this class that beats a read cut chain for chain; it stays a script
(`bench/famcmp/shift_rule.py`), not a driver stage and not a default, until the user decides.
HG002 was used before for frozen rules (rows 1692, 1717-1722).

## DERIVATION RESULT (2026-11-14, before HG002; `r1726_shiftrule/{or,a}.txt`)
OR6737 qualifying cells, largest first (dropped; in RefSeq vs read cut; pairs lost): breaking / 0.05 / any 1,528
(0.029 vs 0.156; 0.014); breaking / 0.05 / 50 1,095 (0.031 vs 0.153; 0.010); breaking / 0.1 / 10 1,023 (0.034 vs
0.154; 0.011); any / 0.05 / 10 1,013 (0.029 vs 0.155; 0.009); breaking / 0.05 / 10 674 (0.021 vs 0.168; 0.004).
Q = 0.02 drops nothing on either sample.
On A119b as they stand: the first four fail (b) (ratio to the read cut 0.81, 0.63, 1.01, 0.67; the first also (c)
0.033). The fifth passes: 749 dropped, 0.037 in RefSeq vs 0.093 for the read cut (0.40), pairs lost 0.006.
**Frozen cell: frame-breaking shift (not a multiple of 3) of <= 10 bp, weaker chain under 0.05 of its partner's
reads.** Precision OR6737 0.3517 -> 0.3547 (read cut 0.3534), A119b 0.1764 -> 0.1768 (0.1766): small, as stated.
It is the smallest qualifying cell, taken by the registered order, and the A119b margin (0.40 vs 0.5) is thin.

## RESULT HG002 (held out; `r1726_shiftrule/hg.txt`, class table `hg.class.txt`)
Assembly 153,033 chains, 0.234 in RefSeq; S1 25,880 chains, 0.238 in RefSeq, moved intron annotated 0.633; oracle
ceiling over the read curve +0.010.
| bar | HG002 |
|---|---|
| 1. dropped in RefSeq <= 0.5 x the read cut's | 0.060 vs 0.062 (0.98) FAIL |
| 2. annotated pairs lost <= 0.02 | 0.009 PASS |
| 3. correct chains lost <= 0.5% | 0.15% PASS |
| 4. >= 600 dropped | 860 PASS |
Precision 0.2336 -> 0.2346, read cut at the same count 0.2345. Moved intron annotated in 0.226 of the dropped chains.

**VERDICT row 1726: FAIL held out.** On HG002 the frozen rule removes chains that are annotated exactly as often as
those a plain read cut removes: it adds nothing over the read count. The advantage shrank from OR6737 (0.12 of the
read cut's share) to A119b (0.40) to HG002 (0.98). A119b and HG002 are both human, so this is not a species effect;
the cell was the one survivor of five chosen on the discovery samples. No other
cell was scored on HG002. Not a stage, not a default; `shift_rule.py` stays as the scoring script.

Where the two rows leave the class: on the human samples the oracle ceiling over a read cut is +0.007 / +0.010, the
weaker chains are in RefSeq as often as the average chain (0.188 vs 0.176; 0.238 vs 0.234), and annotated pairs are
emitted at the cassette control (row 1714). Neither side has headroom that a GTF-level rule reaches. Row 1725's snap arm was already a read-level
test (it moved 1,390 of 6,989 near pairs on OR6737, 795 of 11,947 on A119b); untested are read features beyond its
edit-distance comparison, which would first need a yardstick for unannotated sites.

Comparator check after review: the registered read cut breaks ties in table order (sorted by contig), so at 674-860
chains it is the 2-read chains of the first contig. Genome-wide, 2-read chains (no chain has fewer) are in RefSeq
0.170 (OR6737), 0.066 (A119b), 0.066 (HG002). Against that the rule's dropped chains are 0.021 (0.12 x), 0.037
(0.56 x), 0.060 (0.91 x). HG002 does not flip; on A119b the cell would have missed bar (b) (0.56 vs 0.5), so by the
genome-wide comparator no cell held on both discovery samples and HG002 would not have been opened. HG002 was
scored for this one cell only. Row 1724 (another session, fraction-rule exemption on the recall side) is pending.

## Row 1727 (pre-registered 2026-11-14, before any KB3781 table of this class): does the rule hold on gorilla libraries?
Seen: the frozen cell on OR6737 / A119b / HG002 against genome-wide 2-read chains: 0.12 x / 0.56 x / 0.91 x; the
OR6737 oracle ceiling +0.022 against +0.007 / +0.010 on human. Not seen: anything of this class on KB3781.
Reading to test: the rule works on the gorilla testis libraries and not on the human ones (as `capstarts` works on
human libraries only). KB3781 (`r1622/kb/kb.gtf`, gorilla RefSeq) was spent for row 1677's chain rules; this is a
different rule, frozen elsewhere, and the sample's first use for this class.
- Primary: the frozen cell (frame-breaking, <= 10 bp, < 0.05 of the partner's reads), unchanged.
- Secondary, named now: OR6737's largest qualifying cell (frame-breaking, < 0.05, any shift; OR6737 1,528 dropped,
  0.029 in RefSeq), which failed on A119b.
- Comparator: RefSeq share of all chains at the lowest read count, genome-wide (`shift_rule.py` now sorts by reads
  and takes the tied count's share for the remainder, so table order no longer matters).
- Bars, each cell on its own: dropped chains in RefSeq <= 0.5 x the comparator; annotated pairs lost <= 0.02;
  correct chains lost <= 0.5%; >= 300 dropped (KB3781 has about 70k chains, as OR6737).
- PASS = an empirical rule for gorilla libraries of this kind only, two samples of one species and tissue; it says
  nothing for human libraries, where rows 1725 / 1726 stand. It stays a script until the user decides.

## RESULT row 1727 (KB3781; `r1726_shiftrule/kb.txt`, `kb.class.txt`)
Assembly 72,102 chains, 0.358 in RefSeq; S1 13,511 chains, 0.292 in RefSeq; oracle ceiling +0.023. Comparator
(2-read chains genome-wide) 0.154.
| bar | frozen cell (<= 10 bp) | secondary (any shift) |
|---|---|---|
| dropped (>= 300) | 549 | 1,611 |
| dropped in RefSeq <= 0.5 x comparator | 0.020 (0.13 x) PASS | 0.038 (0.25 x) PASS |
| annotated pairs lost <= 0.02 | 0.002 PASS | 0.015 PASS |
| correct chains lost <= 0.5% | 0.04% PASS | 0.24% PASS |
| precision: before -> rule (read cut at equal count) | 0.3584 -> 0.3610 (0.3600) | 0.3584 -> 0.3658 (0.3631) |
Same cells re-scored with the genome-wide comparator: OR6737 0.021 vs 0.170 (0.12 x), A119b 0.037 vs 0.066
(0.57 x), HG002 0.060 vs 0.066 (0.92 x).

**VERDICT row 1727: PASS, both cells, for gorilla libraries of this kind.** Rule: an emitted chain that differs from
a stronger emitted chain only by one internal exon end moved by a number of bases that is not a multiple of 3, and
has under 0.05 of that chain's reads, is in RefSeq 0.02-0.04 of the time on OR6737 and KB3781, against 0.15-0.17
for 2-read chains. Dropping it costs 0.2-1.5% of annotated shifted-site pairs.
Limits:
- Human libraries: rows 1725 / 1726 stand (0.57 x and 0.92 x; FAIL). The rule is not general.
- Two samples of one species and one tissue (testis), judged by RefSeq only (no CAT). The species difference sits
  in the comparator, not in the rule: dropped chains are rare in RefSeq everywhere (0.021 / 0.020 gorilla, 0.037 /
  0.060 human); 2-read chains are 0.170 / 0.154 on the gorillas against 0.066 / 0.066 on the humans. Why gorilla
  2-read chains score higher is the open part.
- The gain is small: +0.003 / +0.007 precision, +0.001 / +0.003 over a read cut of equal size.
- The comparator is the lowest-read chains, as registered; dropped chains were not matched on their own read count.
- The secondary cell was chosen on OR6737 and failed on A119b; KB3781 is its one test.
- A script only (`bench/famcmp/shift_class.py` then `shift_rule.py`); no stage, no default, no GTF written.

Corrections and controls after review (post hoc, verdict unchanged):
- KB3781 was not unseen for this kind of class: row 1677's `wobble` class (an intron within 10 bp of one with >= 5 x
  its reads) had shown 0.77 x on KB3781 and "gorillas only". The gorilla-yes, human-no direction was known; row 1727
  is a confirmation with a sharper cell, not a blind test. The "Not seen" line above is wrong on that point.
- Does the shift do the work, or only the low read share? Default chains under 0.05 of their locus's top chain
  (`r1677_chainrules/feat.{or,kb}.tsv`, `locus_frac`) are in RefSeq 0.143 (OR6737, 15,335 chains) and 0.151 (KB3781,
  18,671). The cell's chains are 0.021 / 0.020. So the low share alone does not explain it.
- Does the frame do work? Frame-keeping chains of the same cell: 0.044 (339) / 0.045 (446), about twice the frame-
  breaking ones. Both are far below the low-share control; the frame adds a factor of two.
- The gain over a read cut (+0.001 / +0.003) is far below row 1677's +0.02: a class rule, not an assembly-level one.

## End-to-end check of the wider cell on the gorilla GTFs (2026-11-14, descriptive, no bar; job 4471964)
`shift_rule.py CLASS breaking 0.05 any IN.gtf OUT.gtf` -> `r1726_shiftrule/{or,kb}.rule.gtf` (1,528 / 1,611 chains
removed; precision 0.3517 -> 0.3585 and 0.3584 -> 0.3658). `alt_events.py` on the filtered files, row 1714's
internal shifted-site pairs with minor >= 5 reads, both emitted, default -> rule:
| stratum | OR6737 | KB3781 |
|---|---|---|
| int3 <= 10 bp | 0.981 -> 0.971 | 0.930 -> 0.920 |
| int3 > 10 bp | 0.966 -> 0.957 | 0.938 -> 0.916 |
| int5 <= 10 bp | 0.974 -> 0.972 | 0.956 -> 0.956 |
| int5 > 10 bp | 0.954 -> 0.946 | 0.934 -> 0.926 |
So the wider cell costs up to 1 point of annotated pairs on OR6737 and up to 2.2 points on KB3781 (15 of 637 pairs in
the worst stratum), more than the GTF-only "pairs lost 0.015" suggested and past row 1725's 0.01 guard on KB3781.
The <= 10 bp strata lose 0-1 point. The narrow frozen cell was not run end to end.

## Row 1728 (pre-registered 2026-11-14, user: "Keep looking on human"): read-level features on A119b, step 0
Discovery only, A119b; HG002 stays closed unless a cell qualifies here, and then needs its own registration.
Seen: rows 1725-1727 (human ceiling +0.007 / +0.010; GTF-level cells fail on human; the snap's edit-distance test
moved 795 of 11,947 near pairs). Not seen: any per-read feature at the junctions of this class.
Class: A119b S1 chains (26,591; moved intron K, partner's intron J sharing one boundary). Reads = primary records
carrying K or J exactly (at most 300 each, BAM order). `bench/famcmp/shift_reads.py`. Features of K:
- F1 `indel_K`: share of K's reads with an insertion or deletion within 10 bp of either end of the intron;
- F2 `indel_excess`: F1 minus the same share over J's reads;
- F3 `homopolymer`: a run of >= 5 identical bases within 5 bp of the moved boundary on the reference;
- F4 `intron_ratio`: reads carrying K / reads carrying J (intron level, not chain level).
Cells (10, fixed): F1 >= 0.5; F2 >= 0.25; F3 yes; F4 < 0.05; F1 >= 0.5 or F3 yes; each over all shifts and at
shift <= 10 bp. A chain with no K read found is left out of every cell and counted.
A cell qualifies when it holds >= 500 chains, <= 0.033 of them in RefSeq (0.5 x the genome-wide 2-read share
0.066), <= 0.15 with the moved intron annotated (class 0.61), and <= 0.02 of annotated pairs (chain and partner in
RefSeq) fall in it. **STOP** if no cell qualifies: then read-level slip features do not give a human rule either,
and the class is closed for human libraries with the ceiling as the reason. No cell is added after the table.

## RESULT row 1728 (A119b; job 4471970, `r1728_shiftreads/a.txt`, per-chain `a.reads.tsv`)
S1 chains 26,591, a K read found for 26,383; annotated pairs 4,417; class 0.188 in RefSeq, 0.606 intron-annotated.
| cell | chains | in RefSeq | moved intron annotated | annotated pairs in cell | qualifies |
|---|---|---|---|---|---|
| F1 indel_K >= 0.5, all shifts | 578 | 0.033 | 0.100 | 0.003 | yes |
| F1, <= 10 bp | 310 | 0.032 | 0.061 | 0.002 | no (size) |
| F2 indel excess >= 0.25, all / <= 10 bp | 709 / 343 | 0.056 / 0.058 | 0.158 / 0.093 | 0.008 / 0.004 | no |
| F3 homopolymer, all / <= 10 bp | 1,713 / 468 | 0.201 / 0.254 | 0.583 / 0.641 | 0.069 / 0.024 | no |
| F4 intron ratio < 0.05, all / <= 10 bp | 4,258 / 1,252 | 0.089 / 0.062 | 0.280 / 0.177 | 0.082 / 0.017 | no |
| F1 or F3, all / <= 10 bp | 2,240 / 753 | 0.162 / 0.171 | 0.469 / 0.418 | 0.072 / 0.026 | no |
**One cell qualifies: F1, all shifts** (at least half of the weaker junction's reads carry an insertion or deletion
within 10 bp of the intron). It sits on the RefSeq bar (19 of 578 = 0.0329 vs 0.033), but its moved intron is
annotated in 0.10 against 0.61 for the class, the first feature that separates the class on a human sample. A
homopolymer at the site is no signal (0.58, the class rate). As implemented F3 is a run of >= 5 inside a +-9 bp
window. The cell is 2.2% of the class and 0.26% of emitted chains: whatever HG002 says, the assembly-level effect
is under +0.001 precision.

## Row 1729 (pre-registered 2026-11-14, before any HG002 read-level table): the F1 cell, frozen, on HG002
HG002 on CHM13 (`r1692_trunc5_hg/hg.gtf` and its BAM), S1 class by `shift_class.py`, `shift_reads.py` unchanged.
Bars for the cell "F1 indel_K >= 0.5, all shifts" as printed by the script: >= 300 chains; in RefSeq <= 0.033 (0.5 x
HG002's genome-wide 2-read share 0.066); moved intron annotated <= 0.15; annotated pairs in the cell <= 0.02.
PASS = an empirical rule on two human samples: a weaker shifted-site chain whose junction reads mostly carry an
indel beside the intron is an alignment slip, not an isoform. It would be a clean-up of a few hundred chains, a
script, not a stage or default. The other nine cells print too; they are not read for a decision. HG002 was used
before for frozen rules (rows 1692, 1717-1722, 1726).

## RESULT row 1729 (HG002 held out; job 4471986, `r1728_shiftreads/hg.txt`)
S1 chains 25,880 (K read found 25,831), class 0.238 in RefSeq, 0.633 intron-annotated; annotated pairs 5,573.
| bar | A119b (derivation) | HG002 |
|---|---|---|
| chains >= 300 | 578 | 473 PASS |
| in RefSeq <= 0.033 | 0.033 (19) | 0.036 (17; the bar allows 15) FAIL |
| moved intron annotated <= 0.15 | 0.100 | 0.099 PASS |
| annotated pairs in the cell <= 0.02 | 0.003 | 0.002 PASS |
**VERDICT row 1729: FAIL as registered, on the RefSeq bar by two chains.** The verdict is not softened. What
replicates on the held-out sample: weaker shifted-site chains whose junction reads mostly carry an indel beside the
intron have the moved intron annotated in 0.10 on both human samples (class 0.61 / 0.63) and are in RefSeq 0.033 /
0.036 (class 0.19 / 0.24; 2-read chains 0.066). The RefSeq bar was set at half the 2-read share and A119b sat
exactly on it, so a miss of this size was likely from the start. Size: 473-578 chains per sample, 0.3% of emitted
chains; dropping them would move precision by under +0.001. No GTF was filtered, no stage, no default.
Not tested: the same feature on the gorilla samples; whether these chains' reads are better explained by the
partner's junction (the snap arm of row 1725 asks that with edit distance and is a separate test).

Checks after review (post hoc; verdicts unchanged, no new cut taken):
- Not a read-count effect. Moved intron annotated, cell vs rest of the class, by reads on the weaker junction:
  A119b 2 reads 0.134 vs 0.224 (112 chains in the cell), 3-4 reads 0.086 vs 0.302 (81), 5+ reads 0.094 vs 0.660
  (385); HG002 0.118 vs 0.225 (68), 0.115 vs 0.276 (78), 0.092 vs 0.672 (327). In RefSeq at 5+ reads: 0.021 vs
  0.199 and 0.037 vs 0.253. Two thirds of the cell are junctions with 5 or more reads.
- Row 1728's qualification rested on rounding: A119b's exact 2-read share is 0.0657, half of it 0.03285, and the cell
  is 19 of 578 = 0.03287. Against the unrounded bar the cell misses by 0.00002 and HG002 would not have been opened.
- Overlap with the shipped snap switch (row 1725 arm, A119b): 206 of the 578 chains are already absent from that
  arm; 372 are not. The feature is partly the existing switch and partly new.
- The row shows these chains' moved introns are rarely annotated; it does not show they are alignment errors.
- HG002 is spent for this cell; no further human sample is left for a re-registered bar.

## Row 1730 (pre-registered 2026-11-14, user: "Indel feature on gorilla"; before any gorilla read-level table)
The cell of rows 1728 / 1729, frozen (S1 chains; >= half of the reads on the moved junction carry an insertion or
deletion within 10 bp of the intron; all shifts), on OR6737 (`r1628/orrun/or.gtf`, `GGO_mm.bam`) and KB3781
(`r1622/kb/kb.gtf`, `GCA_029281585.2_flnc_mm.bam`), gorilla RefSeq as the only annotation. `shift_reads.py`
unchanged. Seen: the human results above; on gorilla, rows 1725-1727 (GTF-level only).
Bars, each sample on its own:
1. >= 150 chains (the human cell is 1.8-2.2% of the class; the gorilla classes are 11,838 / 13,511);
2. in RefSeq <= 0.5 x the sample's exact genome-wide 2-read share, computed unrounded from the class table (about
   0.085 / 0.077);
3. moved intron annotated <= 0.15 (gorilla class 0.47 / 0.50);
4. annotated pairs in the cell <= 0.02.
The script's printed "qualifies" column uses the human constants and is ignored. PASS on both = the feature holds on
gorilla libraries too. Descriptive, no bar: the union with row 1727's narrow cell (chains, in RefSeq), and how many
cell chains the narrow cell already holds. No filter, stage or default follows without the user.

## RESULT row 1730 (job 4471992, `r1728_shiftreads/{or,kb}.txt`)
| bar | OR6737 | KB3781 |
|---|---|---|
| 1. chains >= 150 | 655 PASS | 308 PASS |
| 2. in RefSeq <= half the exact 2-read share | 0.026 vs 0.0849 PASS | 0.023 vs 0.0768 PASS |
| 3. moved intron annotated <= 0.15 | 0.040 (class 0.473) PASS | 0.042 (class 0.495) PASS |
| 4. annotated pairs in the cell <= 0.02 | 0.005 PASS | 0.001 PASS |
**VERDICT row 1730: PASS on both gorilla samples.** Descriptive: 326 / 108 of the cell's chains are also in row
1727's narrow cell; the union of the two is 1,003 chains at 0.030 in RefSeq (OR6737) and 749 at 0.024 (KB3781).
Across the four samples the cell's chains are in RefSeq 0.026 / 0.023 / 0.033 / 0.036 and their moved intron is
annotated in 0.040 / 0.042 / 0.100 / 0.099 (class 0.47-0.63): the feature behaves the same everywhere; the human
rows failed on the comparator (human 2-read chains are themselves only 0.066 in RefSeq), not on the feature.
Still true: human FAIL as registered (rows 1728 rounding, 1729); the cell is 0.3-0.9% of emitted chains; nothing
shows these chains are alignment errors, only that they are rarely annotated; gorilla judged by RefSeq alone.
Post hoc, junctions with 5+ reads: cell vs rest in RefSeq 0.011 vs 0.333 (OR6737, 371 chains), 0.016 vs 0.328
(KB3781, 187); moved intron annotated 0.024 vs 0.578 and 0.043 vs 0.562. Not a read-count effect on gorilla either.

## Row 1732 (pre-registered 2026-11-14, before any filtered file of the union): the two rules applied to the assemblies
Filter = drop the transcripts of chains in row 1727's narrow cell (S1 / S2, frame-breaking shift <= 10 bp, < 0.05 of
the partner's reads) or in row 1730's indel cell. Gorilla samples OR6737 and KB3781, bars on each:
1. precision of the filtered GTF above that of a read cut of the same number of chains (genome-wide tie share);
2. correct chains lost <= 0.5%;
3. row 1714's four internal shifted-site strata (minor >= 5 reads), both isoforms emitted: no stratum falls by more
   than 0.01 (`alt_events.py` on the filtered file).
Human, descriptive only (rows 1728 / 1729 stand as FAIL): the indel cell alone on A119b and HG002, same three
figures, pairs on A119b only. PASS = a measured assembly improvement on gorilla; the filter stays a script
(`bench/famcmp/shift_filter.py`) and its output a separate file, until the user decides on a stage.

## RESULT row 1732 (job 4472004; `r1728_shiftreads/{or,kb,a,hg}.filt.gtf`, `.filt.txt`, `.events.txt`)
| | OR6737 (both rules) | KB3781 (both rules) | A119b (indel only) | HG002 (indel only) |
|---|---|---|---|---|
| chains dropped, in RefSeq | 1,003, 0.030 | 749, 0.024 | 578, 0.033 | 473, 0.036 |
| precision -> filtered (read cut) | 0.3517 -> 0.3561 (0.3542) | 0.3584 -> 0.3619 (0.3606) | 0.1764 -> 0.1767 (0.1767) | 0.2336 -> 0.2342 (0.2341) |
| 1. above the read cut | PASS | PASS | equal | +0.0001 |
| 2. correct chains lost <= 0.5% | 0.12% PASS | 0.07% PASS | 0.05% | 0.05% |
| 3. internal pair strata, largest fall <= 0.01 | 0.020 (int3 <= 10 bp, 101 -> 99 of 103) FAIL | 0.011 (int3 <= 10 bp, 174 -> 172 of 187) FAIL | 0.004 | not run |
| other three strata, largest fall | 0.003 | 0.001 | 0.001 | |
**VERDICT row 1732: FAIL as registered, on bar 3, by two annotated pairs in the smallest stratum on each gorilla
sample.** Bars 1 and 2 pass on both: the filtered gorilla assemblies are +0.004 precision (+0.002 / +0.001 over a
read cut of the same size) for 0.1% of correct chains. On human the indel rule alone equals a read cut. This is the
measured size of what the two rules do to an assembly: a small, real clean-up on gorilla with a cost of two
annotated acceptor-shift pairs per sample, and nothing on human. The filter is `bench/famcmp/shift_filter.py`; not
a stage, not a default. Families and assignment were not run on the filtered files.

Numbering: this row was pre-registered as 1731; another session wrote a different row 1731 to the register while
its job ran, so it is recorded as row 1732. The pre-registration text and bars are unchanged.

Diagnostic after row 1732 (post hoc, no rule taken from it): annotated internal pairs (minor >= 5) whose weaker chain
the filter removed: OR6737 6 (3 by the frame cell, 3 by the indel cell), KB3781 3 (2 frame, 1 indel). In the failing
stratum (acceptor shift <= 10 bp): OR6737 one of each, KB3781 two by the frame cell. Neither rule alone carries the
loss, so dropping one of them is not a clean repair. One removed chain is an annotated 3 bp shift with all 14 of
its junction reads carrying an indel: the indel cell does remove some real isoforms.

## `shiftdrop` stage and its families / assign measurement (pre-registered 2026-11-14 before the job's output; user: "Port opt-in stage for gorilla")
Port: `bench/famcmp/shiftdrop.py GTF BAM OUT_GTF [THREADS]` (annotation-free; both rules) and driver stage
`shiftdrop` -> `PREFIX.shiftdrop.gtf`, named only, never part of `all`; `PREFIX.gtf` is not changed and no other
stage reads the new file. Register number taken when the result is written (numbers are being used by another
session). Job 4472020, workspace `r1733_shiftdrop` (directory name only).
- Parity: the stage's output must be byte-identical to row 1732's scored `r1728_shiftreads/{or,kb}.filt.gtf`.
- Families / assign, each gorilla sample, two arms from today's binaries on the same inputs: base (the assembly's
  `families.gtf`) and drop (the same file without the transcripts the stage removed). "No harm" bars: pair precision
  and pair recall (gorilla protein-hit truth, `pairscore.py`, both copies tables) each within 0.01 of base; assigned
  rows within 2% of base. A gain needs +0.01 on precision or recall with the other not below base - 0.01.
- Whatever the outcome the stage stays opt-in and gorilla-labelled; this measurement only says whether using its
  output as the families input would matter.

## RESULT `shiftdrop` stage (register row 1734; job 4472020, workspace `r1733_shiftdrop`)
- Port: `bench/famcmp/shiftdrop.py`, driver stage `shiftdrop` in `tools/rustle_pipeline.sh` (usage line, stage
  function, dispatch), dispatch test `shiftdrop_needs_the_assembly_and_then_writes_the_shiftdrop_gtf` (passes).
- Parity: byte-identical to row 1732's scored filter on both samples. OR6737 1,003 of 73,348 chains dropped (frame
  348, indel 329, both 326), 48 s; KB3781 749 of 72,102 (441 / 200 / 108), 97 s, 36 threads.
- Families / assign, base vs drop: the families input loses 1,003 / 749 transcripts; the copies table (family,
  copy, position columns) and the assignments file are **identical** in both samples (OR6737 pair score 0.358 /
  0.754 both arms, 127 assigned rows of 10,036; KB3781 0.219 / 0.450, 228 of 58,265). The loci files differ in a few
  records only. "No harm" bars met with zero change; no gain.
**Reading:** the stage cleans the transcript catalogue (about 1% of chains, 2-3% of them annotated) and does not
reach families or assignment at all: weaker shifted-site chains are never what a locus or a copy rests on. It is an
opt-in for the isoform output of gorilla testis libraries, with the row 1732 cost (two annotated acceptor-shift
pairs per sample). On human libraries each rule alone failed its bars (rows 1726, 1728, 1729); the two together were not run there. Uncommitted.
To remove: delete `stage_shiftdrop`, its usage lines and dispatch entry in `tools/rustle_pipeline.sh`, the test, and
`bench/famcmp/shiftdrop.py`.

Corrections after review: the bar named both copies tables. Compared as whole sorted files: KB3781 `fam.copies.tsv`
and `famreps.fam.copies.tsv` identical between arms; OR6737 both differ in one row, one field (the locus hull start
of one copy, 81588228 -> 81588462); assignments identical in both samples. The pair scores printed by the job
(0.358 / 0.754, 0.219 / 0.450) are equal between arms but far from the recorded gorilla baselines (this job's
families run was not set up like those), so they are not gorilla family scores and only their equality is used.
The full dispatch test file passes (21 tests). The stage was ported on the user's go although row 1732 failed its
annotated-pair bar by two pairs per sample.

## Second gorilla annotation (pre-registered 2026-11-14 before any table with it; user: "can we download or lift another gorilla annotation??")
Source: T2T primates CAT consensus, `mGorGor1_v2.0_CatRefConsensus_v3.bb` (GenomeArk
`species/Gorilla_gorilla/mGorGor1/assembly_curated/gene/`), 156,649 transcripts on GenBank names CM0554xx.2; the
NCBI assembly report maps them to our NC_ names with identical lengths for all 25 nuclear chromosomes (relationship
"="), so coordinates are used as they are. Workspace `r1735_ggo_cat` (register number taken when the result is
written). The set uses RefSeq among its inputs: a wider yardstick, not an independent one; it makes the gorilla
figures comparable with the human ones (RefSeq + CAT).
Measurement: `shift_class.py` with the CAT introns as the extra annotation, `shift_reads.py` and `shift_rule.py`
unchanged, OR6737 and KB3781. Yardstick = moved intron annotated in RefSeq or CAT. Chain identity stays RefSeq.
Bars, each sample:
1. indel cell (row 1730): moved intron annotated <= 0.15;
2. frame cell (row 1727, narrow): moved intron annotated <= 0.15;
3. sanity, not a bar on the rules: >= 0.80 of RefSeq introns of coding genes are in CAT (else the name map or the
   conversion is wrong and nothing is read).
PASS = the wider annotation agrees that the chains `shiftdrop` removes are rarely annotated. FAIL on 1 or 2 = part of
what the stage removes on gorilla is annotated alternative sites that RefSeq lacks, and the stage's doc must say so.
Reported beside: the class share (RefSeq only it was 0.473 / 0.495).

## RESULT second gorilla annotation (register row 1738; workspace `r1735_ggo_cat`: `cat.gff3`, `{or,kb}.class.tsv`)
Sanity: RefSeq introns of coding transcripts (NM_ / XM_) in CAT 202,094 of 243,769 = 0.829 (bar 0.80, met; non-coding
0.345; all 0.733). CAT has 300,862 introns. Moved intron annotated in RefSeq or CAT:
| | OR6737 | KB3781 |
|---|---|---|
| class S1 (RefSeq only: 0.473 / 0.495) | 0.612 | 0.638 |
| 1. indel cell <= 0.15 | 655 chains, 0.072 PASS | 308, 0.097 PASS |
| 2. frame cell (narrow, S1 part) <= 0.15 | 674, 0.120 PASS | 549, 0.202 FAIL |
| union = what `shiftdrop` removes (S1 part) | 1,003, 0.123 | 749, 0.186 |
**VERDICT: FAIL as registered on bar 2 (KB3781); bar 1 PASS on both.** With the wider annotation the gorilla class
looks like the human one (0.61 / 0.64 against 0.61 / 0.63), and the indel rule's chains are annotated 0.07 / 0.10
(human 0.10 / 0.10): that rule is the same on all four samples. The frame rule is weaker than RefSeq alone showed:
12% / 20% of the chains it removes carry a moved intron that CAT annotates. So part of what `shiftdrop` removes on
gorilla is annotated alternative sites that RefSeq lacks; the stage's documentation now says so. CAT includes RefSeq
models among its inputs, so this is a wider yardstick, not an independent one. Chain identity was not re-scored
against CAT.

## `shiftdrop` cut to the indel rule (pre-registered 2026-11-14 before the job; user: "Indel rule only")
The stage now drops only S1 chains with an indel within 10 bp of the moved intron in >= half of the reads carrying
it. Four samples, each on its own: OR6737, KB3781, A119b, HG002 (assemblies as above).
Already seen: the cell on all four (in RefSeq 0.026 / 0.023 / 0.033 / 0.036; moved intron annotated with RefSeq + CAT
0.072 / 0.097 / 0.100 / 0.099); the indel-only filter on A119b (precision 0.1764 -> 0.1767, read cut 0.1767, pairs
largest fall 0.004) and HG002 (0.2336 -> 0.2342, read cut 0.2341). Not seen: annotated pairs end to end on OR6737,
KB3781 and HG002 for the indel rule alone; gorilla precision against the read cut for it alone.
Bars:
0. parity: the stage's output is byte-identical to `shift_filter.py ... indel` on each sample;
1. correct chains lost <= 0.2%;
2. row 1714's four internal shifted-site strata (minor >= 5 reads), both isoforms emitted: no stratum falls by more
   than 0.01;
3. dropped chains in RefSeq <= 0.05 and moved intron annotated (RefSeq + CAT) <= 0.15;
4. precision of the filtered GTF not below the equal-size read cut (unrounded).
PASS on all four = the stage is a safe clean-up of the isoform catalogue on every sample tested; its gain is known to
be small (under +0.003 precision). It stays opt-in and outside `all` whatever the outcome.

## RESULT `shiftdrop`, indel rule only (register row 1739; job 4472048, workspace `r1739_shiftdrop_indel`)
| bar | OR6737 | KB3781 | A119b | HG002 |
|---|---|---|---|---|
| chains dropped | 655 of 73,348 | 308 of 72,102 | 578 of 219,517 | 473 of 153,033 |
| 0. parity with `shift_filter.py ... indel` | PASS | PASS | PASS | PASS |
| 1. correct chains lost <= 0.2% | 0.07% PASS | 0.03% PASS | 0.05% PASS | 0.05% PASS |
| 2. internal pair strata, largest fall <= 0.01 | 0.0097 (101 -> 100 of 103) PASS | 0.001 PASS | 0.004 PASS | 0.006 PASS |
| 3. dropped in RefSeq <= 0.05; moved intron annotated <= 0.15 | 0.026; 0.072 PASS | 0.023; 0.097 PASS | 0.033; 0.100 PASS | 0.036; 0.099 PASS |
| 4. precision vs equal-size read cut (unrounded) | 0.35460 vs 0.35330 PASS | 0.35988 vs 0.35932 PASS | 0.176743 vs 0.176657 PASS | 0.234190 vs 0.234098 PASS |
| precision before -> after | 0.3517 -> 0.3546 | 0.3584 -> 0.3599 | 0.1764 -> 0.1767 | 0.2336 -> 0.2342 |
| stage wall, 72 threads | 42 s | 95 s | 172 s | under 589 s (includes A119b's filter step) |
**VERDICT: PASS on all four samples, all five bars.** The stage is a safe clean-up of the isoform catalogue on two
gorilla and two human samples: it removes 0.2-0.9% of chains, 2-4% of them in RefSeq, 7-10% with an annotated moved
intron, and costs at most one annotated shifted-site pair per stratum.
Limits, stated with the verdict: the gain is small (+0.0003 to +0.003 precision; +0.0001 to +0.0013 over a read cut
of the same size; on the human samples the margin over the read cut is 0.0001); bars 3 and the human figures were
seen before this registration, and OR6737 passes bar 2 by one pair (a second lost pair would fail it); rows 1728 /
1729 remain FAIL as registered against their stricter comparator bar; nothing shows the removed chains are
alignment errors; families and assignment were shown unchanged for the two-rule version on gorilla (row 1734), not
re-run for this subset. Opt-in, never part of `all`, uncommitted. Dispatch tests: the `shiftdrop` test and 18 others pass; two `all` tests fail
because the driver now runs `capstarts` inside `all` (another session's edit in progress, not touched here).

Scope of row 1739's PASS, added after review: it is a no-harm confirmation, not a held-out validation. Bar 3 and the
human side of bar 4 were set after their values were seen; new information was the annotated pairs on OR6737, KB3781
and HG002 and the gorilla precision against the read cut. HG002 failed this cell as registered in row 1729. The
A119b assembly is the one of rows 1714-1723, not today's default. "Families / assign unchanged" was measured for
the two-rule version on gorilla only; for the shipped stage the fact is that no other stage reads its output.
