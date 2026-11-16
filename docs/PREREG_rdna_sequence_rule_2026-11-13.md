# PREREG r1691: a sequence rule for rRNA / rDNA families, the part the soft mask cannot do (2026-11-13)

Follows row 1690 (user goal: soft masking to avoid rRNA and other low-complexity regions while still calling gene
families). The mask reaches satellites everywhere but rRNA only where it is RepeatMasker-made, and rDNA spacer transcripts
nowhere. Written before any alignment to the reference below.

## Rule
Reference = one human rDNA unit with its spacer (CHM13 chr21:3,201,678-3,246,465, 44,787 bp: the RNA45SN1 start to the
next unit's start) plus one 5S unit (chr1:227,800,448-227,802,688). A copy counts when >= 50% of its spliced sequence
(`PREFIX.fam.copies.fa`) is covered by `minimap2 -x asm20` alignments to the reference; a family is listed when strictly
more than half of its copies count (base copies table), as for the BED and soft-mask lists.
- S1: the reference as it is.
- S2: the same reference with its interspersed and simple repeats replaced by N (every RepeatMasker interval of the unit
  except class rRNA). This is where the soft mask is used: the spacer holds Alu and simple repeats that would attract
  unrelated transcripts.

## A119b (reference cut from the same genome: development only; run `r1662_edgemiss/asg_def`, labels of row 1690)
- Reach: >= 13 of the 15 rDNA-label families listed.
- Cost: no protein-coding truth pair lost (pairscore unchanged with the listed families removed) and <= 5 listed families
  with a none or ct label.
- Freeze the arm that passes with fewer none / ct listings (tie: S2). Neither passes: STOP.

## Held out (frozen arm; other individuals and species, the reference stays human)
- HG002 paternal: the row-1673 real secondary-only loci (`r1673_sechap/L.fa`, `L_content.tsv`) whose CHM13 placement lies in
  a censat rDNA array: >= 0.90 of them count under the rule (per locus sequence, no family).
- HG002 default families (`r1662_edgemiss/e2e_hg`), OR6737 (`e2e_or`), KB3781 (`kb_all`): no protein-coding truth pair lost
  with the listed families removed; the rDNA-label family of each gorilla run (row 1690: MCL4, MCL97) is listed.
- PASS all: built as an opt-in next to `RUSTLE_EXCLUDE_SOFTMASK` on the assign skip list, the reference shipped in the repo.

## A119b outcome and freeze (2026-11-13, `r1691_rdna/a.s{1,2}.*`): both arms pass; S2 frozen
| arm | rDNA-label families listed (of 15) | none / ct listed | pair precision / recall (full 0.902 / 0.808) | ignored non-coding pairs (full 50,704) |
|---|---|---|---|---|
| S1 unit as is | 15 | 3 (2 of unknown repeat class, 1 SINE) | 0.902 / 0.808 | 4,874 |
| S2 repeats of the unit masked | 13 | 0 | 0.902 / 0.808 | 4,875 |
S2 frozen (fewer none / ct listings). Its two misses (MCL66, MCL69) are simple-repeat spacer families; MCL66 falls under the
soft-mask rule of row 1690, MCL69 under neither.

## Held out, frozen S2 (2026-11-13, `r1691_rdna/{hg,or,kb}.*`, `hgL.s2.tsv`): FAIL as registered on the HG002 loci reach
| test | result | bar |
|---|---|---|
| HG002 row-1673 loci placed in a CHM13 rDNA array (92) counted | 74 (0.804): rRNA-class 30 of 30, Alu-rich spacer 43 of 61, simple repeat 1 of 1 | >= 0.90 FAIL |
| HG002 families: listed / pairs lost | 17 of 537, every one on unplaced scaffolds (six of them hold 150-250 copies) / 0 (0.727 / 0.444 unchanged) | PASS |
| OR6737: rDNA-label family MCL4 listed / pairs lost | yes, plus MCL9 (13 copies, no label) / 0 (0.359 / 0.759) | PASS |
| KB3781: rDNA-label family MCL97 listed / pairs lost | yes / 0 (0.219 / 0.450) | PASS |
The masking that made S2 clean on A119b is what misses the spacer loci: they are the unit's Alu sequence.

## Post hoc: S1 (the unit as it is) on the same held-out sets; not the frozen arm
HG002 array loci 92 of 92; HG002 families 20 listed (the 3 extra are scaffold families of 4, 114 and 108 copies), 0 pairs
lost; OR6737 3 listed, KB3781 3 listed (the extras are 2 and 4, 16 copies, no label), 0 pairs lost in both. With A119b
(15 of 15, 3 none-label families, 0 pairs lost) S1 meets every bar of this row, chosen after seeing the held-out result.

## Reading
- A human rDNA unit finds the rDNA families of another human and of gorilla with no protein-coding pair lost in four
  samples; the soft mask finds none of them outside CHM13.
- With the spacer's repeats masked it is cleaner but misses Alu-rich spacer transcripts; unmasked it takes them all and
  lists a few small unlabelled families (3 / 0 / 1 / 2).
- Not built into the driver by this verdict (the frozen arm failed; Amendment 1 below wires it opt-in). `tools/sequence_families.py` is the standalone tool (any reference
  FASTA); the references are `r1691_rdna/rdna_s{1,2}.fa`. Not tested outside primates, and no assign run with the list.

## Amendment 1 (2026-11-13, written before the run): opt-in switch and its confirming run
Built under the user's goal "finish applying pending tasks for completing the objectives"; the verdict above is unchanged
(FAIL as registered on the frozen masked reference). The rule is wired as opt-in `RUSTLE_EXCLUDE_FASTA=REF.fa` on assign's
skip list (joined with the BED and soft-mask lists), the two references are kept in `tools/data/` (`rdna_human.fa` = S1,
`rdna_human.masked.fa` = S2). It does not change the region-parallel defaults (those follow the hotspot switches only).
Confirming run, A119b, same inputs as the row 1690 run: arm M = `RUSTLE_EXCLUDE_SOFTMASK=0.9` (the existing `asg_mask`
run, 601 s), arm MR = the same plus `RUSTLE_EXCLUDE_FASTA=tools/data/rdna_human.fa`. Bars, all must hold for PASS:
- B1 the driver's list (`a1.assign.seqexcl.txt`) is identical to the standalone tool's list on the same copies FASTA;
- B2 on families skipped in neither arm, assigned rows differ by at most 1% in total;
- B3 wall of MR at most 1.10 x a same-job rerun of M (rerun, so both arms share node and cache state);
- B4 `check_run_tables.py` reports no failed check on MR.
The largest single process (ru_maxrss of the driver's children, not the job total) is recorded for both arms, no bar. A FAIL leaves the switch in as opt-in
only if B1 and B4 hold; otherwise it is removed.

### Amendment 1 outcome (2026-11-13)
A119b confirming run (job 4470811, same node, M first): PASS on all four bars. B1 the driver's list = the standalone tool's list (18 families, 10 not already on the soft-mask list). B2 families skipped in neither arm: 5,784 vs 5,784 assigned rows, 0 families differ; the 10 newly skipped families held 4 assigned rows and 10,338 of 296,370 rows. B3 wall 373 s vs 675 s for the soft-mask-only rerun (M ran first on a cold cache, so the speed-up is an upper bound; the bar of 1.10 x holds either way). B4 no failed table check. Largest single process of the driver (ru_maxrss of its children, not the job total) 47.8 GB vs 107.7 GB: the largest process was working on rDNA families. The 18 families are also the list row 1691 scored on A119b (0 protein-coding truth pairs lost); the 10 not already on the soft-mask list are MCL3, 4, 34, 40, 60, 69, 89, 90, 1066, 1193, and the 4 assigned rows were all in MCL3.
The switch stays opt-in; defaults unchanged. Not tested: any sample other than A119b through the driver, non-primate rDNA.

## Amendment 2 (2026-11-15, written before the run; user: "lets do the HG002 run for the rDNA exclusion")
Question behind it: may `RUSTLE_EXCLUDE_FASTA=tools/data/rdna_human.fa` become a default of `assign`? This run does not
decide that; it supplies the missing measurement on a second human sample. The register number is taken when the row is
written.

Scope, stated first: **this is not a held-out test of the rule.** The shipped reference is S1 (unmasked), and row 1691
already scored S1 on HG002 post hoc (92 of 92 array loci, 20 families listed, 0 protein-coding pairs lost). Row 1691 stays
FAIL as registered on the frozen masked reference. Not seen on HG002: any `assign` product with the switch, the assigned
rows the listed families hold, the effect on other families, wall and memory.

Run: HG002 Kinnex on its paternal assembly (`~/data/hg002/hg002.pat.fa`, `~/work/jfim/hg002_run/hg002.pat.bam`, SEDEF BED
`sedef_build/hg002pat/final.bed`); the rDNA loci are on that assembly's unplaced scaffolds. Assembly = the existing
`~/runs/hg002_kinnex/hg.gtf` (2026-10-29, not today's assembler default); `families` is run fresh through today's driver.
Binaries = a copy of today's `target/release` (sole-candidate promotion off). Arm R = default + `RUSTLE_EXCLUDE_FASTA`,
run first (cold cache counts against R); arm D = default, same families products, `assign` only, under a 3 h limit. No
soft mask, no BED in either arm. 72 threads, `RUSTLE_TIMING=1`.

Bars, all must hold for PASS:
- B1 the driver's list (`hg.assign.seqexcl.txt`) is identical to the standalone tool's list on the same copies FASTA;
- B2 on families skipped in neither arm, assigned rows differ by at most 1% in total;
- B3 the listed families hold at most 1% of arm D's assigned rows (families and rows reported);
- B4 wall of R's `assign` at most 1.10 x D's; largest single process recorded, no bar;
- B5 `check_run_tables.py` reports no failed check on R.
If D hits the limit or is killed: a cost result; B2-B4 are "not measured" and the verdict is INCOMPLETE, not PASS.
Reported without a bar: number of listed families against row 1691's 20 (the families table is rebuilt, so the count may
differ). No default is changed by this run whatever the outcome; the reference is one human unit, gorilla was not run
through the driver, nothing outside primates was tested.

Amendment 2, note before any result (2026-11-15 03:35): the first job (4472167) was cancelled with arm R's `assign`
unfinished after 88 minutes. Cause: the BAM was read from `~/work` at under 1 MB/s (one thread in disk wait, 0.3 cores
used); 1,085 regions were done, the pending ones were a chr14 region of 8.0 M records and two chr14 / chr22 families of
10 M records each (MCL125 / MCL126, not on the rDNA list). No assignments file existed; nothing was scored. Resubmitted
unchanged except: the BAM is first copied to the workspace `stage/` (same file), and arm R also runs under the 3 h
limit. Bars unchanged; an arm that hits the limit makes the verdict INCOMPLETE. Partial products kept in `r_slowio/`.

### Amendment 2 outcome (2026-11-15; register row 1741; job 4472222, workspace `r1741_rdna_hg`)
**VERDICT: INCOMPLETE as registered. Arm D (default, no exclusion) was killed by the node's out-of-memory handler.**
| bar | result |
|---|---|
| B1 driver list = standalone tool's list | PASS: identical, 20 of 537 families (row 1691's S1 count was 20) |
| B2 families skipped in neither arm, assigned rows within 1% | not measured (arm D wrote no assignments) |
| B3 listed families hold <= 1% of D's assigned rows | not measured |
| B4 wall of R <= 1.10 x D | not measured as a ratio: R finished in 1,631 s; D was killed at 1,128 s |
| B5 no failed table check on R | PASS (one informational line) |

- Arm R (default + `RUSTLE_EXCLUDE_FASTA`): `assign` exit 0, 1,631 s, largest single process 31.5 GB, 213,987 rows,
  2,898 assigned (2,607 molecules), 162,716 ambiguous, 48,373 tied. Slowest region 66 s.
- Arm D (default): `assign` exit 137 after 1,128 s, largest single process 239.8 GB, on a 72-core `singlenode` node;
  Slurm recorded one oom_kill event. 525 regions were done. The last family logged was building profiles for 249
  copies; the cross-chromosome families the log names include all 20 listed families (MCL0-MCL7 hold 108-249 copies
  each and are all on the list). Which family's process reached 239.8 GB was not isolated.
- Cost result, stated as measured: in this run, at the driver's thread defaults (4 family threads, 16 star threads,
  72 threads), today's default `assign` was OOM-killed on HG002's paternal assembly; with the rDNA exclusion it
  completed in 27 minutes at 31.5 GB. An earlier default `assign` on this sample did complete (2026-10-31, older
  binaries and families table, settings not recovered), so the default can finish under other settings. Whether the
  239.8 GB comes from the rDNA families themselves or from several large families running in parallel is not known. What the exclusion costs in assigned rows
  on HG002 is unknown, because the comparison arm does not exist. On A119b the listed families held 4 assigned rows.
- Not a held-out test of the rule (see the scope above); row 1691 stays FAIL as registered. One sample, one node type,
  one run per arm. The assembly is the 2026-10-29 `hg.gtf`, not today's assembler default. No soft mask or BED in
  either arm. Gorilla was not run through the driver; nothing outside primates was tested.
- Side finding: `~/work` served the BAM at under 1 MB/s for the large regions (first job, cancelled); the BAM is now
  also at workspace `stage/hg002.pat.bam` (131 GB). The three largest regions (8-10 M records, chr14 / chr22, not
  rDNA) finished in arm R once the BAM was on the fast disk.
- No default changed. Whether the exclusion becomes a default of `assign` is the user's decision.

Amendment 2, descriptive addition (written 2026-11-15 before looking): the default `assign` of 2026-10-31 on this sample
(`~/runs/hg002_kinnex/hg.assign.assignments.tsv`) did complete. The rDNA rule is run on that run's `hg.fam.copies.fa`
and the rows of the listed families are counted in its assignments. It scores no bar (older families table, older
defaults, sole-candidate promotion on, so assigned rows are overstated); the verdict stays INCOMPLETE.
Result of the descriptive addition: on the 2026-10-31 run the rule lists 14 families; they hold 45,589 of 247,796
assignment rows (18%) and **0 of 17,960 assigned rows**. On that run the exclusion would have removed no assignment.

## Amendment 3 (2026-11-15 evening, written before the run; user: "yes lets do that" to "run the HG002 check for the rDNA skip and make it default if it passes")
Correction of what I told the user: I called the gap "one registered rerun with the unmasked reference". The loci reach of the
unmasked reference on HG002 is already known post hoc (92 of 92), so rerunning it would decide nothing. The measurement that
is missing is Amendment 2's comparison arm (bars B2-B4), lost to the out-of-memory kill. This amendment reruns it; like
Amendment 2 it is **not a held-out test of the rule** (row 1691 stays FAIL as registered on the masked reference).

Run: as Amendment 2 (same families products `r1741_rdna_hg/r/hg.fam.*`, same pinned binaries `r1741_rdna_hg/bin`, 72 threads,
3 h limit per arm), but both arms on one `spr1tb` node (1 TB), arm R first again, in new directories `r3/` and `d3/`; the driver
is a snapshot of today's `tools/` taken at submission (`snap/`), so later driver edits cannot reach the job.

Bars unchanged: B1 list identical to the tool's; B2 families skipped in neither arm, assigned rows within 1% in total; B3 the
listed families hold at most 1% of arm D's assigned rows; B4 wall of R at most 1.10 x D; B5 no failed table check on R.

Decision (the user's, given beforehand): **all five PASS -> `RUSTLE_EXCLUDE_FASTA` defaults to `tools/data/rdna_human.fa` in
`assign` (empty value = off).** Any bar FAIL -> stays opt-in. Arm D killed or over the limit again -> INCOMPLETE; no switch is
made, and the fact that the default cannot finish HG002 on a 1 TB node is put to the user.

### Amendment 3 outcome (2026-11-16; register row 1758; job 4473464 on spr1tb, workspace `r1741_rdna_hg/{r3,d3}`)
**VERDICT: FAIL on B3; the switch is not made (the registered decision).** Both arms completed.
| bar | result |
|---|---|
| B1 driver list = tool's list | PASS: identical, 20 families |
| B2 families skipped in neither arm, assigned rows within 1% | PASS: 2,898 vs 2,898, no family differs (rows 213,987 in both) |
| B3 listed families hold <= 1% of D's assigned rows | **FAIL: 51 of 2,949 = 1.73%** (MCL0 50, MCL6 1: the 249- and 114-copy rDNA array families on unplaced scaffolds); they hold 45,921 of D's 259,908 rows |
| B4 wall R <= 1.10 x D | PASS: 1,811 s vs 3,557 s (0.51) |
| B5 no failed table check on R | PASS (one informational line) |
Largest process: R 30.4 GB, D 417.7 GB. So the default arm does finish HG002, on a 1 TB node only; Amendment 2's kill at 239.8 GB
was the node's limit. Not a held-out test.
