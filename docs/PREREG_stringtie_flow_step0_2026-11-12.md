# PREREG r1658: StringTie flow, step 0 = the retired Rust port (registered 2026-11-12, before building or running)

User choice 2026-11-12: "StringTie flow step 0". r1634 / r1639 closed the StringTie match short of its flow (best Jaccard 0.519).
Found while orienting: commit `2e046730` (2026-07-14) deleted a Rust port of StringTie's long-read assembler with its network
flow (`src/rustle/max_flow.rs` etc., ~106k lines); the handoff note of that tree records 1,710 of 1,839 StringTie transcripts
matched on gorilla chr19. So step 0 is that tree, not a new prototype: same question (does the flow close the gap), far cheaper.

## Procedure
- `git archive 2e046730^` into the workspace `r1658/src` (no checkout, no worktree; `~/Rustle` untouched), build there.
  Wall limit: one short round of mechanical build fixes; if it does not build, stop and report.
- Reads: gorilla OR6737 (`GGO_mm.bam`), contig NC_073242.2. Target: `r1636/plain.gtf` (StringTie 3.0.1 `-L`, reproduced there).
- Arms (fixed): A = `rustle -L` as built, defaults (StringTie-exact mode is its default). B = A + `--vg --vg-solver em-hmm
  --genome-fasta` (the handoff note's benchmark command). No arm may read a StringTie GTF (the old `RUSTLE_CHIMERA_FILTER_GTF`
  filter chain is excluded: it is not a rule the assembler could apply alone).
- Score: `bench/famcmp/tool_match.py` chain reader; Jaccard, reached, inside vs StringTie; correct and precision vs RefSeq;
  new and lost chains vs today's default counted separately.

## Bars and stop rule
- "Almost exactly" = Jaccard >= 0.80, reached >= 0.90, inside >= 0.90 (as r1634), on the contig. Baseline: default 0.414.
- PASS on the contig -> whole gorilla, then A119b once with the arm frozen; then ask the user how the mode should exist
  (the 106k lines are not restored into `~/Rustle` without that decision).
- FAIL -> do not debug the old port. Report, and ask before a new prototype on StringTie's own graph.

## Result (2026-11-12, `r1658/{a,b}.gtf`; built with crate download after `--offline` failed on memmap2; 6 s / 34 s on 16 threads)
StringTie on the contig: 3,089 chains, 966 correct. RefSeq 4,324.
| arm | chains | correct | precision | reached | inside | Jaccard | new vs default (correct) | lost vs default (correct) |
|---|---|---|---|---|---|---|---|---|
| default (today) | 3,182 | 1,063 | 0.334 | 0.594 | 0.577 | 0.414 | | |
| A (retired port, `-L`) | 3,260 | 1,000 | 0.307 | 0.916 | 0.868 | 0.805 | 1,336 (130) | 1,258 (193) |
| B (A + `--vg`) | same chain set as A | | | | | | | |
FAIL as registered: Jaccard and reached pass, inside 0.868 < 0.90. Stop rule applies: the old port is not debugged and the
whole-genome / A119b runs are not launched without the user's ruling. What it shows: the flow is the lever (0.414 -> 0.805);
the price on this contig is -63 correct chains and -0.027 precision (130 gained, 193 lost).

## Amendment D1 (2026-11-12; user chose "Debug toward 0.90 inside", overriding the do-not-debug stop rule; registered before running)
Diagnosis (StringTie's decision log, `r1636/run1.jsonl`): of the port's 429 extra chains, 169 are chains StringTie builds and
then removes in its final filter (isofrac 61, retained intron 58, print_predcluster flag-off 24, printResults delete 19,
included 7); 215 it never builds; 45 are contiguous parts of its chains. StringTie-only 258: 117 different combinations of
introns the port emits, 74 with an intron the port lacks, 67 parts of a port chain.
Arms: the port's existing switches that its source labels as StringTie-parity for those two filters, each alone and all
together, nothing new written: `RUSTLE_ISOFRAC_ST_COV`, `RUSTLE_ISOFRAC_PER_MAXINT`, `RUSTLE_RI_USE_RAW_FLOW`,
`RUSTLE_RI_ZERO_INTRON_LOW`, `RUSTLE_RI_J0_UNIQUE_INTRON_OFF`, plus `RUSTLE_NOVEL_JX_MAX_USAGE=3` (the handoff note's
precision filter that reads no GTF). Pick = highest Jaccard on the contig among arms meeting all three bars; if none meets
them, STOP (no source edits to the old port in this amendment). The contig is then spent: the picked arm is frozen and run
once on whole gorilla and once on A119b, same bars; those two runs are the test.

### D1 result (2026-11-12, `r1658/d1.*.gtf`): STOP, no arm meets the inside bar
| arm | chains | correct | precision | reached | inside | Jaccard |
|---|---|---|---|---|---|---|
| A (no switch) | 3,260 | 1,000 | 0.307 | 0.916 | 0.868 | 0.805 |
| ISOFRAC_ST_COV | 3,203 | 990 | 0.309 | 0.901 | 0.869 | 0.793 |
| RI_USE_RAW_FLOW | 3,291 | 1,007 | 0.306 | 0.919 | 0.862 | 0.801 |
| RI_ZERO_INTRON_LOW | 3,243 | 1,000 | 0.308 | 0.914 | 0.870 | 0.804 |
| all six | 3,212 | 997 | 0.310 | 0.900 | 0.866 | 0.790 |
ISOFRAC_PER_MAXINT, RI_J0_UNIQUE_INTRON_OFF and NOVEL_JX_MAX_USAGE=3 leave the output unchanged. Closing the last 0.03 needs
source edits to the old port's final filter (169 chains StringTie builds then removes); not done here.

### Descriptive: arm A unchanged, whole gorilla and A119b (not a test; bars reported, no rule picked from it)
Whole gorilla (2026-11-12, `r1658/gw.or.txt`, 388 s at 64 threads; StringTie 67,491 chains, 23,248 correct):
| arm | chains | correct | precision | reached | inside | Jaccard | new vs default (correct) | lost vs default (correct) |
|---|---|---|---|---|---|---|---|---|
| default | 73,348 | 25,794 | 0.352 | 0.639 | 0.588 | 0.442 | | |
| A (retired port) | 72,500 | 24,309 | 0.335 | 0.929 | 0.865 | 0.812 | 27,077 (3,011) | 27,925 (4,496) |
A119b: the port is killed for memory on the whole BAM at 64 and at 8 threads (jobs 4468186, 4468218). Per-contig run with a
150 GB address-space cap (job 4468232): chr1-4, chr7 and chr12 abort on allocation, the others finish; the table on the
finished contigs goes below when the job ends (a partial genome, biased away from the largest contigs).

A119b, 18 of 25 contigs (2026-11-12, job 4468232, `r1658/gw.a.txt`; chr1-4, chr7, chr12, chr14 aborted on memory; StringTie,
RefSeq and the default restricted to the same contigs; StringTie there 123,224 chains, 18,835 correct):
| arm | chains | correct | precision | reached | inside | Jaccard | new vs default (correct) | lost vs default (correct) |
|---|---|---|---|---|---|---|---|---|
| default | 130,207 | 22,942 | 0.176 | 0.409 | 0.387 | 0.248 | | |
| A (retired port) | 137,334 | 20,008 | 0.146 | 0.784 | 0.704 | 0.590 | 81,462 (2,760) | 74,335 (5,694) |
Reading (descriptive): the port moves A119b from 0.248 to 0.590 but is far from the bars there; its gorilla level (0.81) does
not transfer, consistent with a port developed against gorilla chr19. Price on these contigs: -2,934 correct chains net
(2,760 gained, 5,694 lost), precision 0.176 -> 0.146. StringTie match: not reached on either sample as registered.

## Amendment D2 (2026-11-12; user chose "Fund the StringTie port work"; registered before any source edit)
Context: StringTie's own output scores below today's default on both axes on both samples (23,248 at 0.344 vs 25,794 at
0.352 gorilla; 31,781 at 0.150 vs 38,715 at 0.176 A119b), so this work buys a StringTie-like mode, not a better price.
Pre-checks (no edits): A119b rescored without strand = 0.590 (strand is not the cause); every finished contig sits at
0.49-0.71; intron-set Jaccard port vs StringTie 0.863 on A119b against 0.963 on gorilla; half of the port-only chains
carry an intron StringTie uses nowhere (0.525 A119b, 0.551 gorilla), 0.40 of StringTie-only chains an intron the port lacks.
- Dev panel: gorilla NC_073242.2 and A119b chr20 (StringTie traced with its decision log on both; chr20 must first be
  shown to reproduce the benchmark GTF).
- Order: stages (junction set, path choice, final filter) by their ceiling on chr20 = Jaccard if that stage's
  disagreements vanished.
- A fix is kept only if inside rises and reached does not fall on both panel contigs. Edits only in workspace `r1658/src`
  (pristine copy kept, `diff -ru` exported); nothing enters `~/Rustle`.
- Budget: 4 rounds; stop early when a round gains < 0.01 Jaccard on chr20, or when the stage ceilings show 0.80 is out of
  reach on chr20 without rewriting a stage.
- Frozen test afterwards: whole gorilla, the other finished A119b contigs, and KB3781 against a fresh StringTie `-L` run
  (both registered samples are tuning data now). Bars unchanged (0.80 / 0.90 / 0.90). Memory on the large A119b contigs is
  taken up only if the ceilings say A119b is reachable.

### D2 stage ceilings, A119b chr20 (2026-11-12; traced StringTie reproduces the benchmark there, 5,133 / 5,133 chains)
Port 5,763 chains, StringTie 5,133, shared 4,156: Jaccard 0.617, reached 0.810, inside 0.721.
Port-only 1,607: FILTER 689 (StringTie builds then removes: retained intron 283, isofrac 267, predcluster flag-off 78,
printResults delete 44, other 17); PATH 536 (its introns, a combination it never builds); JUNCTION 382 (an intron it never
puts in a path). StringTie-only 977: 560 the port has the introns for, 417 with an intron the port emits nowhere.
Ceilings (Jaccard if a class of port-only chains vanished): filter 0.687, path 0.670, junction 0.654, all three 0.810.
Reading of the stop rule: 0.80 on chr20 needs all three stages and the StringTie-only side; no single stage gets there.
Work continues on the filter (largest ceiling) because the user funded several rounds; the bars are not expected to be met.

### D2 round 1 (registered before running): the filter, by the port's own relaxations
The port's isofrac filter differs from StringTie by documented relaxations tuned on gorilla chr19: accumulated coverage capped
at 2.5 x the dominant (`RUSTLE_ISOFRAC_MULTICOV_CAP`, 0 = no cap, StringTie's behaviour), rescue protection for
`alt_splice_rescue` models, chain de-duplication. Coverage itself agrees (port / StringTie 1.03 median on the killed chains).
Arms, each alone and R1 = all together: `RUSTLE_ISOFRAC_MULTICOV_CAP=0`, `RUSTLE_ALT_SPLICE_PROTECT_OFF=1`,
`RUSTLE_RI_ALT_SPLICE_PROTECT_OFF=1`, `RUSTLE_END_RI_STRICT_OFF=1`, `RUSTLE_RI_J0_MIN_LONGCOV=0`. Kept per the D2 rule
(inside up, reached not down, on both panel contigs).

### D2 round 1 result (2026-11-12, job 4468433, `r1658/r1/`): no arm kept, D2 STOPS
| arm | gorilla contig: reached / inside / Jaccard | A119b chr20: reached / inside / Jaccard |
|---|---|---|
| base | 0.916 / 0.868 / 0.805 | 0.810 / 0.721 / 0.617 |
| MULTICOV_CAP=0 | 0.915 / 0.872 / 0.807 | 0.806 / 0.747 / 0.633 |
| RI_J0_MIN_LONGCOV=0 | 0.917 / 0.864 / 0.802 | 0.811 / 0.708 / 0.608 |
| all five (R1) | 0.916 / 0.868 / 0.804 | 0.807 / 0.738 / 0.628 |
ALT_SPLICE_PROTECT_OFF, RI_ALT_SPLICE_PROTECT_OFF and END_RI_STRICT_OFF change nothing on either contig.
Removing the coverage cap is the only lever that moves A119b (+0.016 Jaccard, inside +0.026) but it lowers reached on both
contigs (-0.001, -0.004), so the keep rule rejects it; kept gain this round = 0 < 0.01, which is the registered stop.
What is left is not a switch: the retained-intron rule (283 chains on chr20), the path stage (536) and the junction stage
(382 + 417) each need the port's code brought into line with StringTie's, stage by stage, with the decision logs.
Frozen test not run (nothing to freeze). StringTie match: not reached; best measured = the unchanged port, gorilla 0.812,
A119b 0.590 (18 contigs), or with the cap removed about +0.016 on A119b at a small cost in reached.

### D2 round 2 (2026-11-12; user chose "Align the retained-intron rule"; one round, code read + decision log, no edit made)
The port's pairwise rule (`transcript_filter.rs`, retained-intron predicate) is StringTie's `retainedintron()` (rlink.cpp)
plus spares that all have switches; with the spares off the output does not change on either panel contig (round 1, D1).
So the rule is already aligned. The disagreement is in what it is fed. A119b chr20, 298 multi-exon chains StringTie kills
as retained intron and the port emits:
| cause | chains |
|---|---|
| victim / killer coverage ratio >= 0.1 in the port, < 0.1 in StringTie | 151 |
| killer chain absent from the port's output (path stage) | 78 |
| ratio < 0.1 in the port too (a low-intron mark or evaluation order differs) | 51 |
| port gives the victim more coverage than the killer | 17 |
| other | 1 |
**Correction (same day).** The first version of this section said the port inflates the coverage of minor isoforms
1.5x and that the cause sits in the flow stage. That is withdrawn: it compared the port's GTF `cov` with StringTie's GTF
`cov`, and StringTie's coverage values are on different scales at different stages (GTF cov / `pred_output` cov median
0.66; `pred_create` cov / `pred_output` cov median ~2165, i.e. a base sum against a per-base value). Against StringTie's
`pred_output` cov the port's GTF cov is 0.92 at the median; longcov port / StringTie is 1.00 (73% of
shared chains within 10%). ("Flat with depth", written here first, was also wrong: see round 3.) The table above also mixes scales in its first row (port ratios from GTF cov, StringTie
ratios from kill-time cov), so its split between "ratio" and "other" is not established; the 78 "killer absent" rows do
not depend on coverage and stand. Round 3 redoes the comparison on like-for-like values (both tools' `pred_filter_stage`
dumps, same stage).

### D2 round 3 (2026-11-12; like-for-like stage dumps, no edit made)
Both tools write a `pred_filter_stage` dump (every prediction with cov and longcov at named filter stages:
`RUSTLE_PARITY_LOG` / `STRINGTIE_PARITY_LOG`, steps `pred_filter_stage,pred_kill`). Run single-threaded on the two
panel contigs (workspace `r1658/r3`, scripts `r1658/stages.py`, `r1658/flips.py`); the GTFs reproduce the base arm
(3,260 / 2,831 shared; 5,763 / 4,156 shared) and StringTie's GTFs are identical to the references. Multi-exon chains:
| stage | gorilla contig port / StringTie / shared / Jaccard | A119b chr20 port / StringTie / shared / Jaccard |
|---|---|---|
| ENTER (before any filter) | 6,592 / 6,456 / 5,699 / 0.775 | 13,966 / 13,364 / 10,462 / 0.620 |
| after pairwise | 4,400 / 3,899 / 3,691 / 0.801 | 8,334 / 7,043 / 5,990 / 0.638 |
| after isofrac | 3,780 / 3,310 / 3,106 / 0.780 | 6,757 / 5,514 / 4,643 / 0.609 |
| FINAL of the filter stage | 3,431 / 3,089 / 2,849 / 0.776 | 6,271 / 5,133 / 4,232 / 0.590 |
| GTF (port's global stages after) | 3,260 / 3,089 / 2,831 / 0.805 | 5,763 / 5,133 / 4,156 / 0.617 |
1. The sample gap is already there before any filter (0.620 against 0.775): the two tools propose different chains on
   the deeper sample. It is not created by the filters.
2. The filters are not neutral, their disagreements offset. Port-only chains at FINAL: StringTie never proposed 299 of
   582 (gorilla) and 1,084 of 2,039 (chr20); it proposed and then removed the rest. StringTie-only chains: the port
   never proposed 99 of 240 and 553 of 901. About half of each side is path stage, half filter stage, as the step-0
   ceilings said: no single stage reaches the bar.
3. Coverage, like for like at ENTER on shared chains: longcov port / StringTie 1.00 at every depth (0.83 / 0.81 within
   10%). cov 0.99 for StringTie cov < 3 but 0.90 (gorilla) and 0.73 (chr20) at cov >= 100; only 0.63 / 0.39 of chains
   within 10%. So the port gives dominant isoforms less cov, which makes relative filters more lenient. StringTie's
   kill-time cov equals its ENTER cov (median ratio 1.000), so this comparison is on one scale.
4. What that costs, measured: StringTie's multi-exon retained-intron kills recomputed with the port's pre-pairwise cov
   for the same victim and killer chains. Gorilla 1,445 comparable kills: 63 flip to "keep", 24 distinct chains reach the
   port's FINAL. chr20 3,128 comparable: 351 flip, 147 reach FINAL. Removing all 147 would lift chr20 Jaccard from 0.617
   to about 0.630. Not comparable: victim chain never in the port (558 / 1,906), killer chain not in the port (9 / 221),
   kill record that did not join to an ENTER record (250 / 693). The isofrac kills have no killer in the log and were
   not recomputed.
Reading: the coverage skew is real but small (at most ~150 chains on chr20 through this filter). The round-2 cause is
replaced by: the chain sets differ before filtering, and filter disagreements on both sides offset. Round closed without
an edit; 3 of 4 rounds used; bar unchanged, FAIL stands (inside 0.868 gorilla contig, Jaccard 0.590 on A119b).


### D2 round 4, registered before the panel run (2026-11-12; user chose "Round 4: coverage skew")
Lead, from one locus only (A119b chr20:53.37 Mb, a chain with cov 275.6 in StringTie and 21.0 in the port, same longcov,
same per-node flux records): StringTie's checktrf pass matches a single-node long-read transfrag (1,453 reads lying
inside the 3' exon) to the kept predictions that contain the node and adds its abundance to them
(`trf_redistribute`, site `long_checktrf`); the port's checktrf loop only attempts the match for transfrags with more
than one node (`path_extract.rs`, `tf_nodes.len() > 1`), so that coverage never reaches the spliced predictions.
Arm: let single-node long-read transfrags take the same match-and-redistribute branch
(`RUSTLE_CHECKTRF_SINGLE_NODE_REDIST=1`, workspace tree `r1658/src` only; default unchanged).
Dev panel: the gorilla contig and A119b chr20, against the same StringTie outputs.
Read-outs: (a) port / StringTie cov on shared chains at ENTER for StringTie cov >= 100 (now 0.90 / 0.73);
(b) retained-intron kills that flip and reach the output (now 24 / 147); (c) Jaccard, inside, reached.
Keep rule: (a) moves toward 1 on both contigs and neither inside nor reached falls on either contig. If kept, the
frozen test is the whole gorilla sample and the 18 A119b contigs with the registered bars unchanged (0.80 / 0.90 /
0.90). If not kept, the round closes and the D2 budget is spent.

### D2 round 4 panel result (2026-11-12)
At the lead locus the switch reproduces StringTie (274.9 against 275.6; the 28-exon chain 695 + 140 against 835).
| read-out | gorilla contig base -> arm | A119b chr20 base -> arm |
|---|---|---|
| (a) cov port / StringTie at ENTER, StringTie cov >= 100 | 0.902 -> 0.905 | 0.732 -> 0.750 |
| shared chains within 10% at ENTER | 0.635 -> 0.647 | 0.392 -> 0.418 |
| (b) flipped retained-intron kills reaching the output | 24 -> 21 | 147 -> 134 |
| (c) chains / shared | 3,260 / 2,831 -> 3,256 / 2,835 | 5,763 / 4,156 -> 5,809 / 4,194 |
| reached / inside / Jaccard | 0.916 / 0.868 / 0.805 -> 0.918 / 0.871 / 0.808 | 0.810 / 0.721 / 0.617 -> 0.817 / 0.722 / 0.622 |
Keep rule met on its letter (a moves toward 1 on both; inside and reached do not fall), so the frozen test runs
(workspace `r1658/f4`). The size is small: single-node redistribution is a real difference from StringTie but explains
little of the coverage skew (0.73 -> 0.75 on chr20), so most of the skew has other sources that this round did not find.

### D2 round 4 frozen test (2026-11-12, workspace `r1658/f4`): FAIL, D2 budget spent
Bars 0.80 Jaccard / 0.90 reached / 0.90 inside, unchanged. Port with `RUSTLE_CHECKTRF_SINGLE_NODE_REDIST=1` against the port as it was:
| sample | arm | chains | correct | precision | reached | inside | Jaccard | new (correct) vs default | lost (correct) |
|---|---|---|---|---|---|---|---|---|---|
| gorilla OR6737 | port | 72,500 | 24,309 | 0.335 | 0.929 | 0.865 | 0.812 | 27,077 (3,011) | 27,925 (4,496) |
| gorilla OR6737 | port + switch | 72,641 | 24,314 | 0.335 | 0.931 | 0.865 | 0.814 | 27,218 (3,019) | 27,925 (4,499) |
| A119b, 18 contigs | port | 137,334 | 20,008 | 0.146 | 0.784 | 0.704 | 0.590 | 81,462 (2,760) | 74,335 (5,694) |
| A119b, 18 contigs | port + switch | 139,571 | 20,057 | 0.144 | 0.795 | 0.702 | 0.594 | 83,454 (2,797) | 74,090 (5,682) |
Gorilla fails on inside (0.865 < 0.90); A119b fails all three bars. The switch is a faithful piece of StringTie that moves
Jaccard by 0.002 to 0.004. Four rounds used, no single cause found: the chain sets differ before the filters and the
filters' disagreements offset. The StringTie match is closed as NOT REACHED with the retired port; the price of the
closest point is on record (gorilla: 3,019 correct chains gained against 4,499 lost, precision 0.352 -> 0.335).
The switch exists only in workspace `r1658/src`; nothing was restored into `~/Rustle`.
