# The consensus defect behind Amendment 14's false flags — read-only diagnosis (Task 5 of the A13 plan), 2026-10-03

Archived from the A13 plan workspace (`task-5-report.md`); referenced by prereg Amendment 15. Scratch paths under `scratchpad/a14_diag/` were session-local.

# Task 5 report: the consensus defect behind Amendment 14's false flags (read-only diagnosis)

**Status: DONE, no code changed.** I reproduced the stage's GWFAM37:c1 consensus **byte for byte** offline (8,201 bp, identical
to `scratchpad/a14_diag/g37_cons.fa`). All 8 GWFAM37 clusters match, and so does the stage's log line. Four more families also
match: all 33 clusters the stage wrote for the 5 families are byte-identical.

## 0. Cause statement

The duplicated sequence enters at the **first polish (T)**, inside `consensus_from_template`. The refinement re-polish and the
merge re-polish reuse the same template and rebuild it.

- **The template.** The medoid is a typical read of the cluster's majority structure: GWFAM37:c1's template is 3,007 bp against a
  median member of 2,908 bp. A minority of members (85 of 368) carry a segment the template lacks.
- **The misplaced gaps.** `MM2_MEMBERS` (minimap2 2.30 `splice:hq -uf`) places each minority read's long insertion at a column
  that depends on where the read starts. It pairs the insertion with a template "intron" (`~`) over the template's own bases
  between that column and the true insertion point. The inserted piece therefore re-carries those template bases.
- **The vote.** The long-insertion rule inserts the most frequent identical insertion at each column once ≥ 3 members carry it,
  whatever its share (`o3_candidates.rs:469`). R5's `~` rule never votes the skipped template bases away (`:454`). Nothing links
  neighbouring columns.

So one biological insertion that passes the ≥ 3 bar at several columns is inserted several times. In GWFAM37:c1 it passes at six
columns, 618-747, with 3-13 identical carriers out of 255-332 covering members (1-5%). Each copy also re-inserts 71-200 bp of the
template. In total the consensus gains 5,213 bp.

## 1. Reproduction: byte-identical

**What was ported, literally:** `port.py` covers:
- the library functions: `parse_paf`, `parse_cs`, `cluster_reads`, `structural_scores` / `structural_template` / `refined_template`
  (eligibility, d, ties, fallbacks), `consensus_from_template` (:433-484), `refine_cluster`, `minimizer_sketch` / `sketch_share`,
  `variant_is_real` and `sample_net`;
- the binary's phase 1: `hits_on_own_target`, `polish`, `seeded_clusters`, `refine`, `merge`, `apply_absorptions` and
  `family_clusters`.

The port follows the binary's flow, file layout and minimap2 command lines. Two differences, neither of which changes a result:
- minimap2 PAFs are cached by content in my scratch directory;
- Python `**` stands in for Rust's `powi` in `variant_is_real`.

**Inputs:**
- **Net:** the family's 1,000 used reads, from `a14/net_reads.tsv` (`used` = 1). My `sample_net` port picks the same 1,000 of 1,315.
- **Read sequences:** from `a14/cand_g4.nets.fa`, written sorted and one record per line, as `Net::write` does.
- **minimap2:** 2.30-r1287, the binary on PATH. The stage ran without `RUSTLE_MINIMAP2`.
- **Source:** the stage binary was built from 75827d5d (logic = 0f5824a7). The worktree differs from it in one doc comment only.

**Result for GWFAM37** (seven minimap2 calls: AVA, T, C, S, merge round 1, R, merge round 2):
- The log line is the same: 5 read clusters, 18 reads split off into 4 new clusters, 1 cluster absorbed in 2 rounds, 1
  LongestAligned fallback, 8 clusters.
- All 8 final clusters match the stage: same members, same consensus bytes.

**Other families checked the same way:** A14 GWFAM425 and GWFAM99, A13 GWFAM37 and GWFAM100. Every cluster written to
`cand.clusters.fa` / `cand.reads.tsv` is byte-identical. The unmatched clusters number exactly the in-reference clusters, which the
stage does not write: 1, 2 and 3. The log lines are identical.

## 2. The mechanism

### 2.1 The template GWFAM37:c1 was polished on

The cluster was read cluster 2 of 5, with 368 members. The medoid is SRR27178663.1064406, 3,007 bp:
- 107 aligned partners, sum d 496, **mean d 4.6**;
- the next best member has mean d 16.3;
- all 368 members are eligible (aligned partners range from 96 to 361; the bar is 50).

The longest members score far worse:

| member | length | partners | mean d |
|---|---|---|---|
| SRR27438213.371558 | 5,700 bp | 194 | 1,323 |
| SRR27438212.2325525 | 5,066 bp | 152 | 1,433 |
| SRR27178663.505141 | 5,011 bp | 124 | 1,445 |

Member lengths: median 2,908 bp, maximum 5,700 bp.

The medoid is therefore a majority-structure read: 283 of the 368 members have no long insertion against it at columns 600-760.
The long reads carry the extra segments, and every partner charges them the pair's indel bases. The template choice is
structurally "right". The defect lies in what the vote does with the minority's insertions.

### 2.2 Where the duplication enters

- **The T consensus is already the final consensus:** 8,201 bp and byte-identical.
- **Refinement:** it kept 361 of 368 members together with their template, and re-polished from the same T hits on the same
  template. The output was the same bytes.
- **Merge:** c1 took no part in it. The only absorption put an S cluster into the cluster that became c2.

The other cases follow the same pattern:
- **GWFAM99:c1:** the merge re-polish (R) chose the same medoid and gave the T consensus byte for byte.
- **GWFAM425:c1:** the T consensus is 11,627 bp. The refinement re-polish keeps the structure and length and changes only bases.
- **A13's GWFAM37:c0** (8,946 bp, S-derived): the same template read, defective at T.

### 2.3 The member alignments: one insertion, placed at many columns

Read SRR27178663.1714971 (3,591 bp) against the template shows the pattern (`@` = template column):

| alignment | cs | matches |
|---|---|---|
| splice:hq, the stage's `MM2_MEMBERS` | `:614 +803@618 ~197@618 :366 -19 :3 * :1803` | 2,786 |
| map-hifi | `:809 +606@813 :368 -19 :3 * :1803` | 2,983 |
| splice:hq, same read trimmed to start at its base 600 | `:9 +809@612 ~203@612 …` | — |

- The trimmed read shows that the insertion's column follows the read's start.
- The 803 inserted bases are T[618..818] (200 bp, exact; asm20 maps the piece's first 200 bp to T[618..818]) followed by 603 bp the
  template lacks.
- Dropping `-uf`, using `splice` instead of `splice:hq`, `-r 2000` and `-z 1000,1000` all give the same placement.

Across the cluster:
- The 85 carriers of this segment get their insertion in **13 ten-bp bins between columns 610 and 750**, at roughly the read's own
  start + 611.
- map-hifi puts all 85 in one bin, at 810.
- Over all long member insertions (the final 361 members, each re-aligned to the template alone): 169 in 36 bins with splice:hq
  against 163 in 5 bins with map-hifi.
- In the stage's own T alignments (368 members), 108 of the 174 long member insertions (62%) have the form `+seq ~n` where `seq`
  begins with the n skipped template bases.

The same form appears in the other defective clusters. The other 22 T clusters of these five families have 0-56 long member
insertions (≤ 22 of this form) and insert 0-3 pieces. Two of them are flagged unions anyway, small-scale versions of the defect:
GWFAM37 T0 → GWFAM37_0 (+389 bp) and GWFAM425 T2 → GWFAM425_1 (+1,571 bp).

| cluster | long member insertions | of this form |
|---|---|---|
| GWFAM37 T1 (c1) | 174 | 108 |
| GWFAM425 T1 | 543 | 221 |
| GWFAM99 T2 | 133 | 67 |
| A13 GWFAM37 T0 | 171 | 106 |

### 2.4 Every long insertion of the final GWFAM37:c1 consensus

| column | length | identical carriers | members with any ≥ 20 bp insertion here (distinct seqs) | covering | prefix = template from the column | piece also in the other pieces (asm20) |
|---|---|---|---|---|---|---|
| 502 | 26 | 20 | 20 (1) | 341 | 1 | no |
| 618 | 803 | 3 | 4 (3) | 332 | **200** (asm20: piece 0-200 → T 618-818) | 717 bp |
| 704 | 717 | 13 | 21 (9) | 288 | **114** | 717 |
| 715 | 706 | 4 | 4 (1) | 279 | **103** | 706 |
| 721 | 700 | 3 | 6 (4) | 273 | **97** | 700 |
| 731 | 690 | 9 | 15 (7) | 258 | **87** | 690 |
| 747 | 674 | 3 | 4 (2) | 255 | **71** | 674 |
| 1531 | 836 | 4 | 7 (4) | 323 | 41 (members show `+836 ~40`) | no |
| 1569 | 61 | 7 | 7 (1) | 317 | 3 | no |

**Reading the table:**
- Every piece from 618 to 747 equals T[column..818] exactly, followed by **the same 603-bp tail** (byte-identical in all six):
  803 − 200 = 717 − 114 = … = 674 − 71 = 603.
- These six pieces add 4,290 bp from one event that 23% of the members carry.
- asm20 reports only the 200-bp prefix against the template; the shorter prefixes are exact matches below asm20's reporting floor.
- The 603-bp tail is not in the template. It is in the consensus six times.

**Where the 8,201 bp comes from:** template 3,007 + 5,213 inserted − 19 bp (a majority deletion of a 19-bp poly-C). No short
insertion was applied.

**Matching the controller's observations:**
- The consensus self-alignment is [644-4460] = [1447-5134] (3,687 matches) and so on.
- The longest reads align to the consensus in two pieces.

### 2.5 The other two bad cases

**GWFAM425:c1 (11,627 bp)**
- **T cluster:** 387 members.
- **Medoid:** SRR27178662.520845, 3,530 bp, mean d 710.6. The next best has mean d 748.0. The longest member (7,558 bp) has mean d
  1,353. Median member length is 3,577 bp.
- **Inserted:** 8,098 bp in 9 pieces, of which two nested series:
  - columns 9 and 12: 900 / 897 bp; 20 / 4 carriers of 207 / 206 covering; template prefixes 180 / 177 bp;
  - columns 210, 214, 216 and 218: 1,490 / 1,486 / 1,484 / 1,482 bp; 26 / 8 / 3 / 3 carriers of 181-147 covering; prefixes
    144-136 bp.
- **Aligner placement:** 522 long member insertions in 30 bins with splice:hq; map-hifi puts 256 in 7 bins and none at 0-20 or
  210-220.

**GWFAM99:c1 (4,240 bp)**
- **T cluster:** 309 members.
- **Medoid:** SRR27438212.1149756, 1,896 bp, mean d 0.0. Several members tie at 0, and the tie goes to the longest. The longest
  member (5,427 bp) has mean d 1,664. Median member length is 1,890 bp.
- **Inserted:** 2,346 bp in 6 pieces:
  - a nested pair at columns 963 and 1025: 976 / 914 bp; 5 / 4 carriers of 272 / 249 covering; prefixes 158 / 96 bp; the two pieces
    share their tail;
  - a pair at adjacent columns 1120 and 1121: 111 / 107 bp, 15 / 13 carriers;
  - 218 bp at column 1040 and 20 bp at column 121.
- **Merge:** the merge re-polish (R) used the same template and reproduced the T consensus.

### 2.6 Why the refinement cannot catch it

361 of 368 members "fit" the 8,201-bp consensus: median `de` 0.0010, median `shorter_cov` 0.998. Yet each of them skips a median
5,213 bp of the consensus (`~` / `-` ≥ 20), and 348 of 361 skip ≥ 1 kb.

The reason is in `refine_cluster` (:508). `de` is gap-compressed, and `shorter_cov` measures the read. Extra sequence in the
consensus is invisible to both.

## 3. How general it is: all flagged unions of A14 and A13

### 3.1 The tests

| test | definition |
|---|---|
| (a) | the union is longer than the longest member read of its candidate's clusters |
| (b) | `minimap2 -x asm20 -X -c` of the union against itself gives an off-diagonal hit ≥ 200 bp |
| (c) | ≥ 50% of the union's inserted bases (runs ≥ 20 bp) occur elsewhere in the same union. Insertions come from the union's genome hit in `iso.base.paf`, the best by identity × coverage (the classify rule), read from the `cg` CIGAR. "Elsewhere" means ≥ 50% of the insertion's 21-mers appear outside its own span. |
| (c') | the same as (c), restricted to insertions ≥ 100 bp |
| (d) | an unaligned end (≥ 20 bp) of that hit is a copy of the union's own sequence |

**The duplication signature:** identity × coverage < 0.999, together with (b) or (c') (call it **T1**), or (c) or (d) (**T2**).

### 3.2 Results

| set | n | (a) | (b) | (c) | **signature T1+T2** | inserted bp (≥ 20-bp runs) | of which own sequence | of which the copy is adjacent |
|---|---|---|---|---|---|---|---|---|
| A14 class c | 54 | 11 | 20 | 38 | **41 (76%)** = 29 T1 + 12 T2 | 23,553 | **22,144 (94%)** | 19,623 (83%) |
| A14 class a / b | 2 | 0 | 0 | 0 | 0 | 0 | — | — |
| A13 D-derived | 30 | 5 | 5 | 4 | **4 (13%)** (GWFAM104_0, 244_0, 331_0, 425_0) | 7,381 | 7,381 | 7,039 |
| A13 survivor-derived | 46 | 7 | 17 | 32 | **35 (76%)** = 24 T1 + 11 T2 | 15,778 | 14,760 (94%) | 12,207 |
| A13 "elsewhere" | 6 | 1 | 3 | 0 | 0 | 28 | 0 | — |
| control: A14 linked cluster consensus sequences (d ≤ δ) | 159 | 58 | **3 (2%)** | n/a | n/a | — | — | — |

**Reading the table:**
- **(a) does not discriminate.** 58 of 159 sound linked consensus sequences are also longer than every member read. Minority-exon
  insertions are routine and harmless when they are genomic.
- **The (b) hits that are not defects.** Two D-derived unions and three "elsewhere" unions are ≥ 0.999 to the genome. Their (b)
  hits are genuine internal repeats.
- **Insertions dominate the 41 class-c unions with the signature.** Median per union: 1 mismatch, 0 deletion bases, 155 inserted
  bases. Insertions are 98.9% of their mismatch + insertion + deletion bases (X = NM − I − D).

### 3.3 The other 13 class-c unions are a different story

**10 mismatch-dominated (X or D ≥ 10):** GWFAM98_0, 161_0, 163_0, 163_1, 244_0, 244_1, 269_0, 269_1, 317_1 and 439_0.
- Median 7.5 reads (8 of the 10 have ≤ 14 reads).
- Median 30 mismatches; insertions are only 42% of their divergence bases.
- This looks like few-read or paralog-mixed consensus, not this mechanism.

**3 other cases:**
- **GWFAM335_0 and 335_1:** the union is exact in two pieces, 631 bp and 3,253 bp, about 204 kb apart on NC_073239.2. That gap is
  beyond splice:hq's 200 kb (`-G 200k`), so the genome hit is split. This is a genome-alignment limit, not a consensus defect.
- **GWFAM47_1:** 310 inserted bp that are not a copy of the union. Unexplained.

### 3.4 A correction to the brief

By NM − I − D, X = 0 and D = 0 hold for 17 of the 54 class-c unions, not all of them. In the insertion-type group the statement is
true in effect.

### 3.5 Why the control fails at 76% while D-derived unions look sound

In a no-deletion control, a surviving copy's cluster is flagged **only** when its consensus sits more than δ (about 1%) from the
genome. So the 54 false flags are, by construction, the defective consensus sequences, and this mechanism makes up 41 of them.

D-derived unions are flagged for the right reason whatever their quality. They therefore show the defect at about its base rate:
4 of 30. The A13 survivor-derived flags and the A14 flags show the selected rate: 76% each.

## 4. Candidate fixes, measured offline

### 4.1 Method

For each variant I re-ran the family's whole phase 1 through the port. Templates, refinement and merge all interact, so the
target is not just re-polished.

- **Which cluster is reported:** the one holding the target's reads (overlap shown in the table file).
- **Genome comparison:** every final consensus is aligned with `MM2_GENOME` to a genome slice and judged by the classify rule on
  the best identity × coverage hit.
- **The slice:** all of the family's copies ± 200 kb from the unmasked `_pri`. It reproduces the stage unions' `iso.base.paf`
  numbers exactly: 0.4854, 0.4944, 0.5911, 1.0000, 0.9996.
- **The variants:**
  - (i) a ≥ 20 bp insertion also needs 2 × count ≥ covering;
  - (ii) the template is the longest eligible member (ties to the smaller name; same fallbacks), in every `structural_template`
    call;
  - (iii) both.

### 4.2 Results

Each cell shows consensus length / identity × coverage against the genome, then the classify fate:

| case | stage | (i) majority | (ii) longest template | (iii) both |
|---|---|---|---|---|
| A14 GWFAM37:c1 (bad) | 8,201 / 0.4854 NewCopy | **2,988 / 0.9997 InRef** | 5,659 / 0.9996 InRef | 5,633 / 0.9996 InRef |
| A14 GWFAM425:c1 (bad) | 11,627 / 0.4944 NewCopy | **3,529 / 0.9994 InRef** | 10,725 / 0.8262 **NewCopy** | 7,558 / 0.9988 Linked |
| A14 GWFAM99:c1 (bad) | 4,240 / 0.5911 NewCopy | **1,894 / 0.9952 Linked** | 5,405 / 0.9993 InRef | 5,405 / 0.9993 InRef |
| A13 GWFAM37:c1 (D, sound) | 3,249 / 1.0000 InRef | **3,033 / 1.0000 InRef** | 5,777 / 0.8558 **NewCopy (broken)** | 4,871 / 0.9994 InRef |
| A13 GWFAM100:c1 (D, sound) | 4,630 / 0.9996 InRef | **4,502 / 0.9996 InRef** | 5,340 / 0.9996 InRef | 5,212 / 0.9996 InRef |

### 4.3 Whole-family effect

This counts every final cluster of the five families judged on the unmasked slice; NewCopy clusters with ≥ 6 reads are would-be
flags.

The three A14 families carry 6 stage flags (GWFAM37_0/1, 425_0/1, 99_0/1):

| variant | would-be flags | of which new |
|---|---|---|
| (i) | **0** | 0 |
| (ii) | 3 (GWFAM425's two remain) | 1: the reads of GWFAM99's in-reference 317-read T cluster |
| (iii) | 1 | 1: the same cluster |

A13 GWFAM37's survivor-derived 8,946-bp flag disappears under (i) and (iii), and stays under (ii).

**What goes wrong under (ii) and (iii) in GWFAM99:**
- The longest member makes a consensus that only 5 of the 317 members fit.
- The refinement splits off the other 312, and their longest member (3,776 bp) gets a misplaced series at columns 214-231.
- The result is 8,144 bp / 0.434 under (ii) and 4,864 bp / 0.7667 under (iii).
- Under (iii), one piece passes the majority test at column 231: 44 carriers ≥ 27 covering. The carriers' own `~` takes them out of
  the covering count.

### 4.4 The A13 controls against the masked genome

This is the genome the A13 stage judged. Both D-derived unions stay NewCopy (still detected) under every variant. The stage itself
measured d 0.2829 and 0.3680 on these; with (i) they are 0.2911 and 0.3501.

### 4.5 What the table says

**(i)** is the only variant that:
- clears all three bad cases;
- keeps both controls at ≥ 0.9996;
- adds no false flag in these five families.

Its cost: minority isoform exons are no longer added. R2's "exon union" becomes the majority isoform; the controls shrink by 216
and 128 bp. The rule is not airtight either: misplaced carriers do not count as covering, so a misplaced insertion can still reach
a majority (seen under (iii)).

**(ii)** is unsafe: it fails GWFAM425, breaks a sound control and creates a false flag.

### 4.6 An extra measurement, not requested

**(iv) Undoing the misplaced gap.** `(+seq, ~n)` becomes `(:n, +seq[n:])` when `seq` begins with the n skipped template bases.
It recovers most of the length but does not clear the cases:

| case | (iv) |
|---|---|
| GWFAM37:c1 | 5,102 bp / 0.9875 |
| GWFAM425:c1 | 5,961 bp / 0.9745 |
| GWFAM99:c1 | 2,350 bp / 0.9481 |
| A13 GWFAM37:c1 | 3,249 bp / 1.0000 |
| A13 GWFAM100:c1 | 5,156 bp / 0.9996 |

(iv) combined with (i) gives exactly (i). The per-column minority insertion is the necessary part of the defect.

### 4.7 What I did not run

I did not run a batch- or control-wide sweep of the fixes. It would spend the A14 control before a preregistered rule exists. It
is cheap if wanted: about 10-50 s per family as light jobs.

## 5. Files

Scratch directory: `/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/774b64db-68a6-4b08-af48-b829dee19664/scratchpad/a14_diag/t5/`

| file | contents |
|---|---|
| `port.py` | the port, plus the variant switches; all off = the stage |
| `run_family.py` | phase 1 of one family, compared with the stage |
| `analyze_cluster.py`, `step2_inserts.py` | template scores and the insertion tables |
| `inserts37.py`, `refinecheck.py`, `misplaced_rate.py`, `tvsfinal.py` | the evidence of §2.3-2.6 |
| `step3_signature.py` → `step3_unions.tsv`, `step3_summary.txt`, `step3_linked_control.tsv`; `step3_classes.py` → `step3_classes.tsv` | §3 |
| `step4_fixes.py` → `step4_fixes.{0..4}.tsv`, `step4_family_fates.{0..4}.tsv`; `step4_masked.py` | §4 |
| `work/<run>_<family>/` | each family's `net.fa`, minimap2 PAFs and `clusters.<variant>.fa` |
| `work/slice_*.fa`, `work/mslice_*.fa` | the genome slices |

## 6. Caveats

- **The misplaced placement is minimap2's behaviour.** It belongs to 2.30-r1287 `splice:hq`. I confirmed it on a single
  read / template pair but did not trace it in minimap2's code. The vote assumes one event lands at one column.
- **Genome-slice fates are a stand-in for the stage's classify.** They match the stage on the 5 unions checked. The family-level
  counts in §4.3 ignore the haplotype classes (C1') and the merging of components.
- **Test (c) uses 21-mer containment, not alignment.** That makes it a sensitive test of "same sequence elsewhere". It would also
  count a genuine tandem duplication polymorphism. In a no-deletion control those should match a haplotype, and class c does not.
