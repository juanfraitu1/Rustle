# PREREG r1622: guided assembly, an optional guide GTF for `assemble` (2026-11-10)

User (2026-11-10): besides the de novo mode, offer a guide annotation that helps assembly, as StringTie (`-G`) and FLAIR do;
separate from de novo, off by default; most end users have a usable annotation. Written before any code or run.
Nothing guided exists in `assemble` today (`families` has a guided mode; `--gff` only tags copies).

## What the guide may do (three arms, each adds one thing)
- **G1 known-chain rescue.** A read whose full intron chain equals a guide transcript's chain supports that transcript. Such a
  chain is emitted with >= 1 read (de novo needs 2) and is exempt from the 2% locus-fraction cut and the fragment (ISM) cut.
  Aimed at the largest de novo miss class (r1614: 69% of the chains StringTie or FLAIR find and we miss have one full-chain read).
- **G2 = G1 + junction correction.** A read junction within 10 bp of a guide junction is moved onto it before chains are compared
  (r1614 wobble class). De novo junction snapping failed (r1527, r1618); the guide gives the target site.
- **G3 = G2 + partial reads.** A read whose chain is a contiguous run of exactly ONE guide chain at its locus supports that
  transcript (the StringTie-like class, 8% of the misses). A run compatible with several guide chains supports none.

## Invariants (checked by a test before any scoring)
- Off by default; without a guide every product is byte-identical to today.
- The guide only ADDS: every de novo transcript line is in the guided GTF unchanged (rescues are added after polish, so they
  cannot shift a locus's fractions and drop a de novo chain: the r1618 lesson). G2/G3 may change read counts of de novo
  chains; then the invariant is on the set of chains, and the count change is reported.
- A guide transcript with no supporting read is never emitted. Emitted guide transcripts carry `guide "ID"`; ends come from the
  reads, not the guide.
- Families and assign keep reading the de novo GTF until a separate non-regression test (r1490: assembly layers broke families).

## Design: the measure cannot be "matches the annotation we guided with"
Guide and truth come from one annotation (RefSeq; gorilla GCF_029281585.2, human CHM13 RefSeq), split by gene: genes whose
gene_id hashes (md5, last hex digit < 5, about 30%) are HIDDEN from the guide and are the truth for novel discovery; the
other 70% are the guide. A second guide adds DECOYS: for every 20th guide transcript with >= 2 introns a copy whose one internal
acceptor is moved 4-10 bp (seeded) to a position that is in no RefSeq intron. Multi-exon exact intron chains throughout.
Dev = gorilla OR6737; held out = human A119b (no guided result looked at there before the arms are frozen).
StringTie `-G` and FLAIR (correct + collapse with the GTF) are run with the SAME guides as trusted reference points.

## Bars (per arm; an arm is offered only if it passes on gorilla and then on A119b)
- B1 no invention: guide transcripts emitted with zero supporting reads = 0; decoys emitted <= 1% of the decoys.
- B2 nothing lost: every de novo chain is in the guided output (count of missing = 0).
- B3 guided genes: correct chains among guide genes >= de novo's + 10%. (By construction a rescued chain is "correct"; this bar
  only shows the guide is used. The informative bars are B1, B4, B5.)
- B4 hidden genes: correct chains among hidden genes >= de novo's, and precision of the non-guide chains (correct hidden /
  all emitted chains not in the guide) >= de novo's - 0.005.
- B5 against the reference points, reported with one bar: non-guide-chain precision >= StringTie `-G`'s on the same guide.
- B6 cost: assemble wall time <= 1.10x de novo.
Descriptive: supported guide chains recovered (a guide chain with >= 1 full-chain primary alignment) by us / StringTie / FLAIR;
rescues by read count; what G2 and G3 add over G1; decoy emission by StringTie and FLAIR.

## Interface (if an arm passes)
`tools/rustle_pipeline.sh assemble --guide ANNOT.gtf [--guide-mode rescue|correct|partial]` -> `copy_assign --guide-gtf`;
product `PREFIX.gtf` with the `guide` attribute; `PREFIX.assemble.guide` records the guide's path and md5.

## Order of work
1. Prototype G1 on the GTF + BAM in Python (floor-1 run already exists: `r1616/floor1`, `r1617/floor1`), score gorilla.
2. G2, G3 prototypes; freeze the arms that pass; A119b once.
3. Rust port of the passing arm with the invariants as tests; the driver flag.

## Amendment 1 (2026-11-10, before any score; job 4459077)
- Decoy as built (`guide_proto.py split`, seed 1622): the genomic right end of one intron is moved 4-10 bp either way (an
  acceptor on +, a donor on -), new intron in no RefSeq transcript. Gorilla: 3,078 decoys; 13,023 of 41,535 genes hidden.
- G1 prototype counts support on primary alignments only (the near-tied secondary alignments of the seeding pool are not
  added), so it is conservative. G1 = de novo chains + guide chains with >= 1 full-chain primary alignment.
- FLAIR is not installed on this cluster; the guided FLAIR reference point is not run in the prototype step. StringTie 3.0.3
  `-L -G` is (the de novo StringTie GTF in the table is the earlier 3.0.1 run).
- For G1, B2 and B4 hold by construction (it only adds guide chains); they become real tests at the Rust port and for G2/G3.

## Amendment 2 (2026-11-10, before any score; job 4459089)
FLAIR 3.0.1 installed (`pixi global`, the version of the de novo run). Guided reference point: BED12 from the same BAM through
FLAIR's own `dofiltering` (`flairbin/flair_bam2bed.py`, rebuilt here: primary alignments, quality 0, supplementary halves
removed), `flair correct -f GUIDE`, `flair collapse -f GUIDE --trust_ends`, reads = `samtools fastq -F 0x900` of that BAM.
Run with the guide and with the guide + decoys. Its rows join the score table when the job ends.

## Outcome, G1 prototype on gorilla (2026-11-10, job 4459077, `r1622/`; FLAIR rows pending job 4459089)
Guide chains 65,968, of them 25,016 with >= 1 full-chain primary alignment; hidden chains 29,842; decoys 3,078.
| arm | chains | guide chains found | hidden found | non-guide chains | non-guide precision | decoys |
|---|---|---|---|---|---|---|
| ours de novo | 73,348 | 17,741 | 8,053 | 55,607 | 0.145 | - |
| ours G1 | 80,697 | 25,090 | 8,053 | 55,607 | 0.145 | 1 of 3,078 has a full-chain read (the scorer's 0 excludes decoys by construction) |
| StringTie de novo (3.0.1) | 67,491 | 15,879 | 7,369 | 51,612 | 0.143 | - |
| StringTie -G (3.0.3) | 47,802 | 30,897 | 5,404 | 16,905 | 0.320 | 34 of 3,078 (1.1%) |
- B1 PASS (1 decoy, 0.03%; no unsupported guide chain). B2, B4 hold by construction. B3 PASS: +7,349 guide chains (+41%);
  6,104 of the rescues rest on one read. **B5 FAIL: non-guide precision 0.145 against StringTie -G's 0.320.** G1 FAILS as registered.
- StringTie -G finds 5,881 guide chains that have no full-chain primary alignment (it builds them from partial reads: the G3
  class) and gives up 1,965 hidden chains against its own de novo run while dropping two thirds of its non-guide chains.
- Post hoc, matched recall: our non-guide chains cut by reads reach StringTie -G's 5,404 hidden chains at about 0.26
  precision (reads >= 5: 5,748 at 0.245; >= 6: 5,321 at 0.267), still 0.055 under it. De novo the two tools were level
  (0.145 / 0.143), so the guide itself cleans StringTie's novel set: reads that fit a known transcript are absorbed by it.
  G1 only adds and has no such step. Not registered yet: an arm that also drops non-guide chains explained by a guide chain
  (it breaks the "guide only adds" invariant, so it is the user's call).

## Amendment 3 (2026-11-10, user: "lets do that"; written before the arms were computed)
Two arms that also REMOVE chains the guide explains. "Guide only adds" stays as the strict mode (G1); these are separate modes.
- **G1a absorb fragments:** G1, minus every non-guide chain that is a contiguous run of a guide chain on its strand (a 5' or 3'
  fragment of a known transcript).
- **G1b = G1a + retained introns:** also minus every non-guide chain obtained from a guide chain by leaving out one or more
  introns (its introns are an ordered subset of one guide chain's and it spans from that chain's first kept to last kept intron
  with the skipped introns inside its exons).
Bars as registered, with two changes for these arms: B2 becomes "every removed de novo chain is an explained chain" (count of
others = 0), and B4's recall half becomes hidden found >= 0.99 x de novo's. B5 (non-guide precision >= StringTie -G's 0.320 on
gorilla) is the bar that failed G1.
Limit of the design, stated now: genes are hidden whole, so a true novel isoform of a GUIDE gene that looks like a fragment or
a retained intron cannot be seen as a loss here. Reported: how many removed chains carry >= 10 reads.

## Outcome, G1a / G1b on gorilla (2026-11-10, computed from the G1 files)
| arm | chains | hidden found | non-guide chains | non-guide precision |
|---|---|---|---|---|
| ours de novo | 73,348 | 8,053 | 55,607 | 0.145 |
| G1 | 80,697 | 8,053 | 55,607 | 0.145 |
| G1a | 74,077 | 8,053 | 48,987 | 0.164 |
| G1b | 73,738 | 8,053 | 48,648 | 0.166 |
| StringTie -G | 47,802 | 5,404 | 16,905 | 0.320 |
- Removed: 6,620 fragments (0 hidden chains, 0 in full RefSeq; 1,963 carry >= 10 reads) and 339 retained-intron chains
  (0 hidden; 101 carry >= 10 reads). Amended B2 and B4 PASS (nothing removed that is not explained; hidden found unchanged).
- **B5 FAIL for both: 0.164 / 0.166 against 0.320.** Fragments and retained introns explain a small part of the gap.
- Post hoc, where the gap sits (non-guide chains by class; correct = hidden RefSeq chains):

| class | ours de novo | StringTie de novo | StringTie -G |
|---|---|---|---|
| guide gene, all introns known | 9,324 | 7,485 | 2,548 |
| guide gene, some intron new | 20,808 | 20,065 | 2,614 |
| hidden gene | 21,340 (8,047 correct) | 19,151 (7,365) | 9,411 (5,403) |
| no annotated intron | 4,135 | 4,911 | 2,332 |

  De novo the two tools look alike in every class. With a guide StringTie drops 87% of the chains that carry a new intron
  inside a guide gene (20,065 -> 2,614) and half of its hidden-gene chains, the correct ones included (7,365 -> 5,403). Inside
  hidden genes alone (no guide present) ours is 8,047 / 21,340 = 0.377 and StringTie -G 5,403 / 9,411 = 0.574: it runs at a
  stricter operating point with a guide, so B5 compares our default point with a stricter one. The curve comparison
  (ours cut by reads at StringTie -G's recall) is the fair form and was 0.26 vs 0.32 before absorption; not redone yet.
- Guided FLAIR rows: job 4459089 failed (`flair correct` 3.0.1 takes no `-g`); resubmitted as 4459155.

## Amendment 4 (2026-11-10, user: "ok lets do that"; written before any of it was computed)
**Arm G1c(K):** G1b, minus every non-guide chain that shares >= 1 intron with a guide chain on its strand, carries >= 1 intron
no guide chain has, and has `reads` < K. K in {3, 5, 10, 20}. Rescued guide chains are untouched.
**Transcript split (new, because the gene split cannot see the cost):** same hash rule on `transcript_id` instead of the gene
(`GUIDE_HIDE=tx`, files `r1622/t.*`). Hidden chains are then mostly unseen isoforms of guide genes: the class this floor and
the G1a/G1b removals can destroy. StringTie -G is rerun on this guide; decoys are not rescored.
**Bars, on the transcript split, both axes against StringTie -G on the same split:** G1c(K) PASSES if hidden found >= StringTie
-G's AND non-guide precision >= StringTie -G's. Frozen K = the smallest passing K; no K passes = FAIL. A frozen K goes to A119b
held out once, same two bars.
Reported without a bar: G1, G1a, G1b on the transcript split (how many hidden chains the fragment and retained-intron removals
cost there), and G1c(K) on the gene split.

## Outcome, Amendment 4 on gorilla (2026-11-10, job 4459165, `r1622/t.*`)
Transcript split: guide chains 65,846, hidden chains 29,964. Our non-guide chains 55,693: fragments 8,368 (277 hidden),
retained 391 (23 hidden), new-intron class 32,752 (3,490 hidden).
| arm | chains | guide found | hidden found | non-guide | non-guide precision |
|---|---|---|---|---|---|
| ours de novo | 73,348 | 17,655 | 8,139 | 55,693 | 0.146 |
| G1 | 80,652 | 24,959 | 8,139 | 55,693 | 0.146 |
| G1a | 72,284 | 24,959 | 7,862 | 47,325 | 0.166 |
| G1b | 71,893 | 24,959 | 7,839 | 46,934 | 0.167 |
| G1c(3) | 60,160 | 24,959 | 7,413 | 35,201 | 0.211 |
| G1c(5) | 51,512 | 24,959 | 6,875 | 26,553 | 0.259 |
| G1c(10) | 45,046 | 24,959 | 6,141 | 20,087 | 0.306 |
| G1c(20) | 41,945 | 24,959 | 5,615 | 16,986 | 0.331 |
| StringTie -G | 47,488 | 32,166 | 4,538 | 15,322 | 0.296 |
- **G1c(10) PASSES both bars** (hidden 6,141 >= 4,538; precision 0.306 >= 0.296); K = 3 and 5 fail precision. **Frozen K = 10.**
- Cost, now visible: the fragment and retained-intron removals lose 300 real hidden isoforms (3.7%), the floor at 10 another
  1,698; G1c(10) keeps 75% of de novo's hidden chains, StringTie -G 56% of ours.
- Still behind on guide chains found: 24,959 vs 32,166 (StringTie -G builds known chains from partial reads; arm G3).
- Gene split, no bar: G1c(3/5/10/20) non-guide precision 0.199 / 0.232 / 0.262 / 0.278 with hidden found 8,051-8,047 of 8,053
  (StringTie -G 5,404 at 0.320).
- FLAIR 3.0.1 `correct` crashed twice: it takes no `-g` (job 4459089), then `ssPrep.py` runs `min(None)` on a read it had
  just reported as skipped (job 4459155). One-line local patch (`if sizes is None: continue`; original kept as `ssPrep.py.orig`
  in the pixi env), so those reads are skipped as the tool's own message says. Rerun = the job submitted after this line.

## Guided FLAIR reference rows, gorilla gene split (2026-11-10, job 4459240, `r1622/fl_guide*`)
| arm | chains | guide chains found | hidden found | non-guide chains | non-guide precision | decoys |
|---|---|---|---|---|---|---|
| FLAIR de novo (earlier run) | 138,125 | 16,737 | 7,556 | 121,387 | 0.062 | - |
| FLAIR correct + collapse with the guide | 48,529 | 20,393 | 1 | 28,136 | 0.000 | 10 of 3,078 |
With a GTF and no short-read junction file `flair correct` keeps only reads whose splice sites it can snap to the guide
("NO NOVEL SPLICE SITES WILL BE DETECTED"), so this mode cannot find a gene outside the guide: 1 hidden chain. Its 28,136
non-guide chains are new combinations of guide splice sites. It is a reference row for known-chain recovery only (20,393 vs
our G1 25,090 and StringTie -G 30,897), not for novel chains.

## Held out: frozen G1c(10) on A119b, transcript split (2026-11-10, job 4459201, `r1622/ht.*`): FAIL on precision
Guide chains 120,728, hidden chains 54,491. Our non-guide chains 192,635: fragments 30,957 (821 hidden), retained 1,639
(64 hidden), new-intron class 88,219 (4,686 hidden).
| arm | chains | guide found | hidden found | non-guide | non-guide precision |
|---|---|---|---|---|---|
| ours de novo | 219,517 | 26,882 | 11,833 | 192,635 | 0.061 |
| G1 | 232,719 | 40,084 | 11,833 | 192,635 | 0.061 |
| G1b | 200,123 | 40,084 | 10,948 | 160,039 | 0.068 |
| **G1c(10), frozen** | 127,754 | 40,084 | **8,713** | 87,670 | **0.099** |
| G1c(20) | 118,840 | 40,084 | 7,924 | 78,756 | 0.101 |
| StringTie -G | 101,906 | 59,668 | 5,334 | 42,238 | 0.126 |
- Hidden found PASS (8,713 >= 5,334); **non-guide precision FAIL (0.099 < 0.126).** G1c(10) FAILS held out; no K in the
  registered set passes on A119b. The floor only touches the new-intron class; on A119b the rest of the non-guide set
  (hidden-gene-like and unannotated chains at 2 reads) is much larger than on gorilla.
- Post hoc, not registered (advisor-model note: compare at matched recall): G1c(10) with EVERY non-guide chain cut by reads.

| non-guide reads >= | gorilla chains / hidden / precision | A119b chains / hidden / precision |
|---|---|---|
| 3 | 15,531 / 5,524 / 0.356 | 57,687 / 7,963 / 0.138 |
| 5 | 12,354 / 4,867 / 0.394 | 40,168 / 7,012 / 0.175 |
| 10 | 9,807 / 4,084 / 0.416 | 28,210 / 5,786 / 0.205 |
| StringTie -G (default point) | 15,322 / 4,538 / 0.296 | 42,238 / 5,334 / 0.126 |

  A 3-read floor on all non-guide chains is above StringTie -G's default point on both axes in both samples, and a 10-read
  floor on A119b still finds more hidden chains at 0.205. Both samples have now been seen for this cut, and StringTie -G was
  not cut by its own coverage, so this is a candidate rule, not a result. Known-chain recovery is still behind on A119b
  (40,084 vs 59,668): arm G3.

## Amendment 5 (2026-11-10, user: "lets register it for KB3781"; written before any KB3781 file for this was made)
**Frozen rule G1d:** G1c(10), then every non-guide chain needs `reads` >= 3. (Guide chains: one full-chain primary alignment.)
**Sample:** gorilla KB3781 (fibroblast; `GCA_029281585.2_flnc_mm.bam`), not used for any guided arm. Limits stated now: same
species, same RefSeq annotation and the same transcript split (`r1622/t.guide.gtf`) as OR6737; the read data are new. Our
de novo GTF is a fresh `rustle_pipeline.sh assemble` at current defaults (`r1622/kb/kb.gtf`); StringTie 3.0.3 `-L -G`.
**Bars:**
- R1 (default point): G1d hidden found >= StringTie -G's AND non-guide precision >= StringTie -G's.
- R2 (curve against curve): StringTie -G's non-guide chains cut at cov >= {0, 1, 2, 3, 5, 10}; ours = G1d cut at reads >=
  {3, 5, 10, 20, 50}. For each StringTie point that keeps >= half of its default hidden found, take our point with the
  smallest floor whose precision >= that StringTie point's; ours must find >= as many hidden chains there. PASS = holds at
  every such StringTie point (no qualifying point of ours = fail at that point).
- PASS both: G1d is the proposed rule for the guided mode's novel chains (port still the user's call). R1 only: reported as
  "above the default point, not across the curve". Otherwise FAIL.

## Amendment 6 (2026-11-10): G3 re-registered as a guide graph per locus. Written before any G3 code or count.
User (2026-11-10): "ideate it as pre-made variation graphs for all loci and then keep them if enough reads support the
models". This replaces the original one-line G3 (unique-run rule) and does not need G2 first.
Wording: this is a **guide splice graph** (exons as nodes, guide introns as edges, each guide transcript a path). It is built
the way the thesis's family variation graph is built from copies, but it lives inside one locus and is not that graph.

**Model.** Per guide locus, the graph of all its guide transcripts. A primary alignment is COMPATIBLE with a path when its
intron chain is a contiguous run of the path's introns and its two ends lie inside the path's exons at those positions (no
intron of the path skipped or read through). A read may be compatible with several paths.

**Keep rule G3(K), frozen form:** a guide path with no full-chain read is kept when
1. every pair of adjacent introns of the path is carried together by >= K compatible reads (single-intron paths: the intron
   by >= K reads), so the path is tiled link by link, and
2. the path is not EXPLAINED AWAY: after the paths with full-chain reads are kept (G1), at least one of its adjacent-intron
   pairs (or, for a path sharing every pair with kept paths, one intron-to-end exon segment) is carried by >= K compatible
   reads that are compatible with NO already-kept path.
K in {1, 2, 3}. G3 only adds guide paths; the novel-chain rules (G1b, G1c, G1d) are untouched and reported unchanged.
Stated now: rule 1 checks neighbouring introns only, so a path that recombines two expressed isoforms at a shared exon can
pass it; rule 2 and the recombinant decoys below are there to measure that.

**Step 0 (count, with a stop):** of the guide chains StringTie -G emits on gorilla OR6737 without a full-chain primary
alignment (5,881 on the gene split), how many pass rule 1 at K = 1. STOP if < 50%: the missing chains are then not a tiling
problem and the model is wrong.

**Tests (gorilla OR6737 gene split to choose K; A119b held out once with the frozen K):**
- T1 splice-shift decoys (the existing 3,078): kept <= 0.5% (StringTie -G keeps 1.1%).
- T2 recombinant decoys (new): for guide genes with >= 2 multi-exon chains sharing an exon, one artificial chain = the 5'
  part of one transcript joined to the 3' part of another at a shared exon, in no RefSeq transcript; up to 3,000, seed 1622.
  Added to the guide for our run and for a StringTie -G run. Bar: our kept share <= StringTie -G's kept share.
- T3 masked recall: 2,000 guide chains with >= 3 full-chain primary alignments (md5 order), their full-chain reads ignored;
  recovered by G3(K) >= 0.50. Reported beside it: how many of them StringTie -G emits on a BAM without those reads.
- T4 reach: guide chains found (G1 + G3) >= 0.95 x StringTie -G's (30,897 on gorilla).
- T5 no invention: kept paths with zero compatible reads = 0 (by construction; checked).
Frozen K = the smallest K passing T1 and T2; then T3 and T4 are read at that K. PASS = all five on gorilla and on A119b.
FAIL on T4 only is reported as "safe but short of StringTie -G" and the arm may still be offered as opt-in; FAIL on T1 or T2
closes the arm at that K.
Not run yet; the KB3781 test of G1d (Amendment 5) comes first.

## Outcome, Amendment 5: G1d on KB3781 (2026-11-10, job 4459336, `r1622/kb*`): FAIL as registered (R1 by 0.0005), R2 PASS
De novo on KB3781: 72,102 chains, 8,158 hidden found, non-guide 54,416 at 0.150. Guide chains found 24,127 (G1) against
StringTie -G's 29,551.
| point | non-guide chains | hidden found | non-guide precision |
|---|---|---|---|
| **G1d (reads >= 3), frozen** | 20,240 | **6,137** | **0.3032** |
| G1d reads >= 5 | 17,195 | 5,605 | 0.326 |
| G1d reads >= 10 | 14,368 | 4,955 | 0.345 |
| G1d reads >= 20 | 8,155 | 3,839 | 0.471 |
| G1d reads >= 50 | 3,810 | 2,433 | 0.639 |
| **StringTie -G default (cov >= 0)** | 15,000 | **4,555** | **0.3037** |
| StringTie -G cov >= 2 | 14,603 | 4,536 | 0.311 |
| StringTie -G cov >= 3 | 13,864 | 4,487 | 0.324 |
| StringTie -G cov >= 5 | 13,425 | 4,404 | 0.328 |
| StringTie -G cov >= 10 | 8,716 | 3,453 | 0.396 |
- **R1 FAIL:** hidden found 6,137 >= 4,555 holds; non-guide precision 6,137 / 20,240 = 0.3032 against 4,555 / 15,000 = 0.3037.
  Short by 0.0005 (11 chains) with 35% more hidden chains. The bar is the bar: FAIL.
- **R2 PASS at all six StringTie points:** at its precision 0.304 / 0.304 / 0.311 / 0.324 our reads >= 5 point finds 5,605
  (against 4,555 / 4,554 / 4,536 / 4,487); at 0.328 reads >= 10 finds 4,955 (4,404); at 0.396 reads >= 20 finds 3,839 (3,453).
- As registered ("PASS both" or "R1 only" or FAIL) this is a FAIL, so G1d with a 3-read floor is not the proposed rule.
  What the curve says is that the floor's value is the open choice: 5 reads is above StringTie -G on both axes at its default
  point on KB3781 (5,605 at 0.326), and was on OR6737 (4,867 at 0.394 vs 4,538 at 0.296) and A119b (7,012 at 0.175 vs 5,334
  at 0.126). All three samples have now been seen for a 5-read floor; no unseen sample is left for it in the current set
  except HG002.

## User decision (2026-11-10, after the KB3781 outcome)
"I don't think it is as important to make sure it works in a sample unseen like HG002, 3 reads sounds reasonable."
**G1d with the 3-read floor is the guided mode's rule for chains outside the guide**, by the user's decision. The record
stands as it is: R1 failed on KB3781 by 0.0005 precision with 35% more hidden chains, R2 passed, and no further held-out test
is planned. Guided mode stays opt-in and separate from de novo; nothing is in the binary yet. Next: G3 step 0 (Amendment 6).

## Outcome, Amendment 6 step 0 (2026-11-10, job 4459360, `r1622/g.tile.tsv`): STOP as registered
Gorilla OR6737, gene split. Guide chains 65,968; 40,952 have no full-chain primary alignment; StringTie -G emits 6,274 of
those (chains counted without strand here; 5,881 with strand in the G1 outcome).
| K | no-full-read guide chains passing rule 1 | of StringTie -G's 6,274 | passing but not emitted by StringTie -G | splice-shift decoys passing |
|---|---|---|---|---|
| 1 | 612 | 499 = **0.080** (stop bar 0.50) | 113 | 1 of 3,078 |
| 2 | 246 | 214 = 0.034 | 32 | 0 |
| 3 | 149 | 133 = 0.021 | 16 | 0 |
- **STOP: link-by-link tiling reaches 8% of the chains StringTie -G adds, not 50%.** The model as registered is wrong for
  this class: these chains are not full paths waiting to be stitched from overlapping partial reads.
- Check of the count: all 25,017 chains with full-chain reads have every pair covered at least that many times.
- Post hoc, what the 6,274 look like: 499 fully tiled, 1,657 miss exactly one adjacent-intron pair (mostly the first or the
  last link, i.e. a terminal exon no read connects), 787 miss two, 2,368 miss five or more; 4,703 of them have >= 6
  introns. Of all no-full-read guide chains, 62% of the links are uncovered. StringTie -G therefore emits known transcripts
  on much weaker evidence than a tiled path (it does not require reads to connect neighbouring introns), which is also why
  it keeps 34 splice-shift decoys where tiling keeps 1.
- Consequence: matching StringTie -G's count of known chains would need the same weak evidence. G3 as a tiled guide graph
  adds 612 chains at K = 1 (+2.4% over G1's 25,090) with 1 decoy; not built further. T1-T5 not run.

## Port (2026-11-10, user: "lets port that with a guide flag, we can optimize it later")
- `src/guide.rs` (rule G1d, 6 unit tests), `utilities guide --gtf DENOVO.gtf --guide GUIDE.gtf[.gz] --bam BAM --out GUIDED.gtf`,
  and `rustle_pipeline.sh assemble --guide GUIDE.gtf` -> `PREFIX.guided.gtf`. `PREFIX.gtf` and every later stage stay de novo.
- Parity with the prototype on gorilla OR6737 (job 4459404, gene-split guide): **47,369 multi-exon chains in both, 0 only in
  one.** De novo models 73,911: kept 40,583, fragments 6,620, retained 339, new intron below 10 reads 17,927, outside the
  guide below 3 reads 8,442; guide transcripts added 7,349. 563 single-exon models pass through (never scored).
- Not optimized: the BAM pass is single-threaded, 1,469 s on the 11.7 GB gorilla BAM (about 3 h expected on A119b's 96 GB).
- Choices in the port: thresholds are constants; rescued guide transcripts carry the guide's coordinates, source `guide`,
  `reads` = full-chain primary alignments, `guide "rescued"`; models sharing a chain share the best-supported one's verdict.
  Not covered: a resumed `all` that skips assemble does not write the guided file (run `assemble --guide` or the utility).
