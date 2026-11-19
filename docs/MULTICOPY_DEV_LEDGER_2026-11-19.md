# Development ledger of the multi-copy goal (copy of PSU workspace goal_multicopy_dev/LEDGER.md, 2026-11-19)

Summary: `MULTICOPY_FAMILIES_2026-11-19.md`. Paths below are PSU workspace paths. No bar was read on any number here.


## Step 1, the measure (2026-11-19, job 57495211, `m1/`)
`bench/famwise_eval/copyset_score.py` (+ tests). Seen quarter: true units one 1,130 / set 105 (>= 2 molecules).
| arm | one-copy precision | set precision | one resolved | set structure / resolved |
| full (whole-sample Rustle GTF) | 1356/1613 0.84 (NARROW 187) | - | 690 | 34 / 0 |
| Rustle P (sel) | 0.92 | - | 778 | 16 / 0 |
| Rustle D (sel, good) | 0.81 (NARROW 245) | - | 784 | 39 / 0 |
| StringTie P / D | 0.92 / 0.92 | - | 529 / 527 | 11 / 0 |
| famwise (exact band) | 0.979 | 53/126 0.42 (WIDE 65) | 642 | 34 / 34 |
Reading: per-copy output's error is NARROW (one copy claimed where reads fit several); famwise sets fix it but half are
WIDE: joined by a single tied molecule (min-link 1) while most supporters are decidable (example: 81 of 81).

## Step 2, question 1 call (`q1/`)
v1 on fw.families.tsv placed/inside: no effect (the wrong sets are mostly one-"locus" family transcripts at distant places).
v2 `resolution_call.py`: per family transcript at >= 2 places, tied share = supporters whose records within W ms reach >= 2
of its places. Sweep W in {0, 10} x share {0.2, 0.5, 0.8}:
| call | one prec | set prec | all prec | resolved one + set |
| famwise | 0.979 | 0.421 | 0.930 | 642 + 34 = 676 |
| W0 0.5 | 0.890 | 0.895 | 0.890 | 689 + 23 = 712 |
| W10 0.5 | 0.964 | 0.591 | 0.941 | 674 + 33 = 707 |
| W10 0.8 | 0.932 | 0.690 | 0.920 | 683 + 31 = 714 |
Candidate (chosen on development): W10 share 0.5. To be judged only in the final row on residues 75-99.

## Step 3, question 3 step 0 (`q3/`, `chain/`): secondaries as seeds vs support only (emulated by `support_filter.py`)
Exact matching: D_supp1 (D kept if >= 1 primary carries the chain) 0.885 / 783 one units vs P 0.921 / 778, D 0.806 / 784.
Chain arm -> famwise -> call (W10, 0.5), exact matching: P_call 0.947 / 789, D_call 0.934 / 800, D_supp1_call 0.942 / 797.
Set units missing structure were mostly read FRAGMENTS (unit chain = a read's chain; 48 of 66 missing in D had the chain
emitted at no copy) -> `copyset_score.py --match contained` (a read chain matches a transcript holding it as a contiguous run).
Contained matching:
| arm | all prec | resolved one + set | set structure |
| P | 0.901 | 1039 + 0 | 63/105 |
| D | 0.749 | 1041 + 0 | 92/105 |
| P_call | 0.920 | 1037 + 38 = 1075 | 63 |
| D_supp1_call | 0.907 | 1040 + 48 = 1085 | 72 |
| D_call | 0.852 | 1026 + 64 = 1090 | 92 |
Sweep (`chain/sweep2`, window 0 / 2 / 5 x share 0.5 / 0.8): narrower windows raise set precision (D_supp1 w0 0.8: 15/17)
but lower overall precision (one-copy NARROW claims stay): W10 0.5 stays the best overall. Candidates for the final row:
P_call and D_supp1_call (W10, 0.5). Answer so far to Q3: secondaries as seeds give copy-set structure (92 vs 63 of 105)
at a precision cost the call recovers only partly; support-only keeps 72.
Remaining error: copy-set claims WIDE (supporters decidable) and one-copy NARROW claims at copies outside any family
transcript.

## Step 4, question 4 ceiling (`diag_q4`): what borrowing could add, at most
Units missed by D under contained matching: 102 of 1,235 (one 89, set 13); 72 of them rest on 2 reads. 73 have a junction
the assembler emitted nowhere at that place (rare isoforms below the floor), 27 have every junction emitted but not the
combination, 2 have a copy with no transcript. So borrowing structure from a paralog can add at most ~8% of units on this
measure, mostly 2-read isoforms; the junction-column arm is worth building only as a labelled class, not for the
headline. Not yet measured: how many of the 73 have the same junction column emitted at a paralog.

## Step 5, question 2: simulated tandem arrays (`arrays_sim/`; famsim, 10 TSPY3 units 20 kb apart, 3 decoys, verify 0 FAIL)
| scenario | primaries on own unit | spread TV | silent units with a chain P / D / D_supp1 | call W10 | call W0 |
| identical, all 10 expressed | 0.10 | 0.05 | - | 1 copy set per isoform (10 units) | same |
| identical, 3 expressed | 0.10 | 0.74 | 7/7 / 7/7 / 7/7 | copy set x2 | copy set x2 |
| identical, 1 expressed | 0.10 | 0.90 | 9/9 / 9/9 / - | copy set x2 | - |
| 0.1% exonic SNPs, 3 expr | 0.81 | 0.14 | 3/7 / 7/7 / 3/7 | copy set x2 (5 places) WRONG | separable x2 |
| 0.2% | 1.00 | 0 | 0/7 / 7/7 / 0/7 | 1 separable, 1 copy set WRONG | separable x2 |
| 0.5% | 1.00 | 0 | 0/7 / 0/7 / - | nothing to call | - |
Answers (development): identical units -> one transcript set over N units is the right output, the call makes it; the
per-unit spread is not recoverable from placements (minimap2 scatters uniformly). Secondaries as seeds put chains on
every silent unit up to 0.2% divergence even when every primary is right; primary-backed filter removes them. The
W10 call is too loose in near-identical arrays; W0 is right there, W10 was better overall on HG002 (step 3).
Candidates for the final row: {P, D_supp1} x {call W0, call W10}, all 0.5 share.

## Step 6, review fixes (2026-11-19)
- Scorer: a tie for the plurality verdict is now MIXED (was EQUAL, which inflated set precision); molecule-weighted share
  added. Rescored (`rescore/`, contained): P 0.896 / 1,039 (weighted 0.655); P_fw 0.913 / 1,056; P_W10 0.914 / 1,075
  (weighted one 0.835, set 0.962, set claims 17/30); P_W0 0.908 / 1,068; S_fw 0.897 / 1,065; S_W10 0.901 / 1,085;
  S_W0 0.870 / 1,072; StringTie P 0.872 / 937. The call adds ~+1.8% over the uncalled layer, +3.5% over P.
- Literal-source truth on the array simulations (`sim_truth.py`, `*/lit/`): identical units -> per-copy claims OFF
  (20 of 23), copy-set claims WIDE (2 of 2) in every scenario; no read-based output is right there; the copy set is the
  least wrong (it holds the source). The goal's "origin" now reads "the copies its reads fit equally".

## Step 7 (2026-11-19, job 57496495, development only)
- Q3 new / lost chains against P (contained, `q3nl/new_lost.tsv`): D (secondaries seed) new 653 = EQUAL 50, NARROW 335,
  OFF 103, MIXED 11, unjudged 154; S (primary-backed) new 161 = EQUAL 27, NARROW 102, OFF 16, MIXED 2, unjudged 14; both
  lose the same 25 (EQUAL 11). Support-only cuts wrong-copy chains ~6x and keeps over half the right new ones; most new
  chains are NARROW (copy-set structure the call can group). The real assembler mode (spliced secondaries as support only,
  mirroring RUSTLE_SEED_UNSPLICED=support) is NOT built: the filter is an emulation.
- Q4 paralog columns (`diag_q4b.py`): of 102 units D misses, 4 have a read with a record at another place whose junctions
  D emits there; borrowing junction columns could add ~0.3% of units on this measure. Q4 closed on this measure (dev).
- Q2 real testis arrays (A119b slice, P -> famwise -> call, `q2a/`): 72 family transcripts at >= 2 places; W0 36 / W10 43
  copy sets: MAGEA 8, AMY 6 (2-7 places), TTTY 5, RBMY 5, GAGE 2-3 (up to 8 places), CT47A 2, CT45A (W10); TSPY's 4 family
  transcripts (2-3 places) are separable. No truth (A119b's Y differs from CHM13's): descriptive.

## Step 8: the real support-only mode (2026-11-19)
`RUSTLE_SEED_SPLICED=support` (opt-in) in `src/family.rs` (Pass1Acc.primary_chains; a spliced chain seeds only when >= 1
primary record carries it, secondaries still add support), guard in `src/copy_assign_bin.rs` (refused off the streaming
path), driver comment, unit test `spliced_support_mode_seeds_only_primary_chains_and_counts_secondary_support` (pass).
Build: `goal_multicopy_dev/target/release` (not the shared binary). Default unchanged: SIRV final.gtf md5 5cb0e7e4 = the
Fritz value; seen-quarter `--seed-pool good` GTF byte-identical to the shared binary's.
Seen quarter (contained): S_real 0.850 / 1,040 (emulation 0.851 / 1,040); new vs P: EQUAL 28, NARROW 105, OFF 16 (emulation
27 / 102 / 16); S_real -> famwise -> call W10: 0.914 / 1,085 (one 0.914, set 33/69, weighted 0.916). Arrays: silent units with
a chain at 0.1% 7 -> 3, 0.2% 7 -> 0 (identical units 7/7, left to the copy-set call). The emulation was faithful.

## Step 9: family-member completeness (user 2026-11-19; StringTie as the reference tool; `members/`, `members_full/`)
Members = 1,243 genes of 365 truth families (famtruth recent paralogues). Found = a transcript maps to the member.
| sample | expressed (>= 2 primaries) | StringTie | Rustle P | Rustle D | missed by all = nothing / reads_there / seedable / claimable | targeted seeding restores |
| HG002 seen quarter (truth) | 338 | 290 | 268 | 272 | 936 = 890 / 42 / 4 / 0 | 0 (StringTie +1 via claimed: 3 reads) |
| HG002 all reads (slice, no truth) | 409 | 348 | 312 | 332 | 768 = 626 / 44 / 98* / - | Rustle +3 of 98 |
| A119b testis (slice, no truth) | 835 | 721 | 646 | 702 | 425 = 289 / 70 / 66* / - | Rustle +7 of 66 |
(* upper bound: a read whose best record lies outside the slice has all its in-slice records counted as tied.)
D finds members with < 2 own primaries (HG002 all 110, A119b 52); on the truth quarter the 15 D-only members: 4 have own
primaries, 11 have no read whose origin is the member (secondary-built shadows, as rows 1671 / 1672).
Readings: missing members are mostly silent (no reads reach them); no member is missing because its reads were placed
elsewhere (claimable 0 on the truth quarter); seeding tied secondaries at the members nobody finds restores 3-7; secondary
seeding's extra members are mostly shadows. The real gap: StringTie finds 24-75 more expressed members than Rustle
(seen quarter: 24 found by StringTie, missed by P and D), all with their own primary reads (the assembler / floor loses them;
the r1534 class).

## Step 10: the member gap (seen quarter, development; `diag_gap`, `gap/`)
The 24 expressed members StringTie finds and Rustle P / D miss: 19 have NO Rustle transcript on the member (2-36 own primaries,
often multi-exon chains of few reads: row 1536's K2 class "unique low-support chains"), 5 sit inside a neighbour's transcript
(OR2T35, TAS2R43 / 46, SIGLEC5, PILRA: fused / neighbour).
Existing opt-in levers on P (members found of 338 expressed; one-copy precision, contained):
| P 268, 0.896 | RUSTLE_READ_ISOFORM_K=2: 277, 0.884 | K=3: 274, 0.888 | RUSTLE_MONO_MODEL=cluster: 268 (no member change), 0.896 |
| K=2 + mono: 275, 0.884 | StringTie: 290, 0.872 |
K=2 closes 9 of the 22-member gap at a precision still above StringTie's. Not judged: HG002 has no unseen quarter left; a
judging row would need KB3781 or another own-haplotype sample (reuse per the user's ruling).

## Step 11: question 4 built as a labelled class (`bench/famwise_eval/borrow.py`, `q4/`, seen quarter, development)
Paralog transcripts projected onto copies where a molecule has a record within 10 ms; kept when a read chain at the copy is a
contiguous run of the projection within 5 bp and the assembly does not hold it within 5 bp; snapped onto the read's
junctions; `fw_class "borrowed_structure"`. First version (exact junctions) kept 0: projected junctions sit a few bases off,
and 5,012 of 8,033 "unheld" secondary chains were junction-shifted copies of assembled chains.
P + 453 borrowed: copy-set units with structure 63 -> 68, members found 268 -> 275 (expressed 268 -> 271); read as
expression the borrowed chains are 15/90 EQUAL (0.17; 363 unjudged), so the class must stay labelled; P + borrowed overall
precision 0.896 -> 0.866 if they were counted.
