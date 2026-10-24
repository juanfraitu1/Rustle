# Round 22: O2 only on ambiguous reads, with sequence ties, on fresh multigene families (2026-10-19)

This pre-registration was written before any world of round 22's families was built; no read of these families has been
simulated. Code: branch `o3-fast`, frozen at the commit that adds this file.

## Why

**Round 21 failed by one read** (register 1331): 5 of 6 cells passed, with C1 passing.
- **The read:** NCOR1P4. Its own copy aligns all 340 bases across 3 introns (AS 262, NM 0). The higher-scoring copy leaves
  38 bases soft-clipped and has 1 intron (AS 278, NM 0).
- **Why it was misplaced:** minimap2's AS charges every intron and nothing for a soft clip, so the margin of 16 made the read
  minimap2's, and it went to the copy that explains less of the read.

**The user (2026-10-19), on the proposed fix:** "ok, lets do that". A read is ambiguous, and so O2's, when its AS is within
10 of the best (Eichler's rule), or when an alignment on another locus explains it with no more edits than the best-AS one,
counting every base it leaves unaligned (`veto_distance`).

**Development (all spent):**
- **Round 21's 36 worlds:** 6 of 6 cells now pass.
- **Round 19's:** 4 of 6, unchanged. The failures are the PABPC1L2B reads, which fit the other copy better by sequence
  (NM 1-2 against 7-8): the split's accepted cost, under the strict B2 the user chose.
- **HG002's untied reads with a donor copy:** 44,086 right and 594 wrong, against v14's 44,021 and 599.
- **HG002's tied reads:** 50 and 7, as v14.

## The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v17** = v16 (round 21) + `--seq-tie`. That is v14 + `--o2-scope tied --as-slack 9 --seq-tie`
(`o2graph.py run W intega`).
- **A read's loci:** its alignments in `E0.bam` within 9 AS of its best, plus any alignment whose `veto_distance` is at
  most the best-AS alignment's.
- **The decision rules are round 21's:**
  - one seeded copy: minimap2's placement, unless its alignment does not explain the read;
  - a tie: v14's call when it is one of the tied copies;
  - a locus off the family's copies: abstain;
  - the silent copy alone: abstain.

**Integrated flags (`$V17`):** round 21's `$V16` plus `--seq-tie`:
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --ref-copies --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
--o2-scope tied --as-slack 9 --seq-tie
```

**The simulation:** as rounds 19-21 (profiles v4).

## Families (held out; drawn by the fixed rule)

**The pool:** the adopted multigene rule with all components (`r19seq/families19cal_full2_all.json`, 439 families).
- **Spent:** every family touching a gene of an earlier simulated set (27 files: `ls families*.json | grep -v _all.json`,
  round 21's included), or of the HG002 worlds.
- That leaves 98 unspent families with 303 copies. File: `r22seq/families22seq_all.json`.

**The draw:** the rule of rounds 17-21 with seed 22.
- **26 families, 90 copies:** 38 coding and 52 non-coding (34 transcribed processed, 16 transcribed unprocessed,
  1 unprocessed, 1 IG_V).
- **Copies per family:** {2: 19, 3: 3, 6: 2, 8: 1, 23: 1}.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families22.json` (set `r22`, sha1 4aadf1d72ccc).

| family | copies | coding genes | non-coding members | chromosomes |
|---|---|---|---|---|
| S19ABCC6P2 | 2 | - | ABCC6P2 (tr. unproc.), ABCC6P1 (tr. unproc.) | chr16 |
| S19AC1143242 | 2 | - | AC114324.2 (tr. unproc.), AC114324.2 (tr. unproc.) | chr5 |
| S19AC1364281 | 2 | AC136428.1 | IGHV3-11 (IG_V_gene) | chr14, chr16 |
| S19AC2337026 | 23 | - | AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.), AC233702.6 (tr. proc.) | chr17 |
| S19ANKRD18B | 2 | ANKRD18B | ANKRD18CP (unproc.) | chr9 |
| S19ANKRD30BL | 6 | ANKRD30BL | FO681491.1 (tr. unproc.), CR383658.1 (tr. unproc.), AP005901.1 (tr. unproc.), ANKRD30BP2 (tr. unproc.), CR383658.1 (tr. unproc.) | chr2, chr13, chr14, chr18, chr21, chr22 |
| S19CDK2AP2P1 | 8 | - | CDK2AP2P1 (tr. proc.), PTGER4P2 (tr. proc.), CDK2AP2P1 (tr. proc.), CDK2AP2P1 (tr. proc.), PTGER4P2 (tr. proc.), PTGER4P2 (tr. proc.), PTGER4P2 (tr. proc.), PTGER4P2 (tr. proc.) | chr9, chr13, chr14, chr21 |
| S19CYP11B1 | 2 | CYP11B1, CYP11B2 | - | chr8 |
| S19DEFB105A | 3 | DEFB105A, DEFB105B | - | chr8 |
| S19FAM205A | 2 | FAM205A | FAM205BP (tr. unproc.) | chr9 |
| S19FOLH1 | 2 | FOLH1 | FOLH1B (tr. unproc.) | chr11 |
| S19GPAA1P1 | 2 | - | GPAA1P1 (tr. proc.), GPAA1P2 (tr. proc.) | chr2 |
| S19KRT17P1 | 3 | - | KRT17P1 (tr. unproc.), KRT17P2 (tr. unproc.), KRT17P6 (tr. unproc.) | chr17 |
| S19LDHAL6B | 2 | LDHAL6B | LDHAL6CP (tr. proc.) | chr12, chr15 |
| S19LILRA2 | 2 | LILRA2, LILRA1 | - | chr19 |
| S19MT1G | 2 | MT1G, MT1H | - | chr16 |
| S19PCDH11X | 2 | PCDH11X, PCDH11Y | - | chrX, chrY |
| S19PDXDC1 | 2 | PDXDC1 | PDXDC2P (tr. unproc.) | chr16 |
| S19PRIM2 | 2 | PRIM2 | AC244258.1 (tr. unproc.) | chr6 |
| S19SPIN2B | 2 | SPIN2B, SPIN2A | - | chrX |
| S19TYW1 | 2 | TYW1, TYW1B | - | chr7 |
| S19ZNF208 | 3 | ZNF208, ZNF676, ZNF728 | - | chr19 |
| S19ZNF506 | 6 | ZNF506, ZNF93, ZNF486, ZNF737, ZNF626, ZNF66 | - | chr19 |
| S19ZNF611 | 2 | ZNF611, ZNF600 | - | chr19 |
| S19ZNF799 | 2 | ZNF799, ZNF443 | - | chr19 |
| S19ZNF98 | 2 | ZNF98, ZNF492 | - | chr19 |

## Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v4.json` and
  `profiles/public_testis_v4.json`.
- **Seeds:** 1201, 1202 and 1203, with no seed selection. That makes 6 cells and 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families22.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout22
for S in 1201 1202 1203; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r22 $S --local $O/local $V17 --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done; done
```

## Bars

**Correctness:** as round 21 (fastset's own e2e step).
- **Table and scoring:** the scoped table `o2.intega.m2.tsv`, per-read detectability less each world's error floor.
- **Set apart:** reads tied with a locus on no copy of the world (`--off-family-ties 9 --seq-tie`) are the class
  `off_family_tie`, left out of both systems' sums.
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**Verdict.** A cell PASSES when all six of its worlds pass both bars, and round 22 PASSES only when all six cells pass (72 bar
checks).

**Cost (C1),** as round 21, on the mainstay set (seed 306, both profiles), into `fast/r22cost/v17/`.
- **Bar:** v17's integrated CPU seconds (O1 + O2 + O2_genome + O2_scope + O3) at most the current method's, in each profile.

**Reported, no bar:**
- v14 on the same worlds and scoring.
- The off-family ties' fates.
- The split of placements between minimap2 (`as_resolved`) and O2: O2's domain size, its placements and errors.
- Wall times and peak RSS.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

## What a PASS would and would not claim

**A correctness PASS** would claim that O2 confined to ambiguous reads, with ambiguity defined by AS within 10 or by
equal-or-better sequence on another locus, does at least as well as the current method in every world type. It would hold
on 26 held-out multigene families, most of them with processed pseudogene copies.

**It would not claim:**
- **Anything about reads tied with a locus off every copy.** These are scored apart.
- **Accuracy where a donor's copy carries its paralog's alleles.** Its reads fit the paralog better by sequence; register
  1330 and round 19's development.
- **Real-read accuracy beyond the HG002 development numbers above.**

## Amendment 1 (2026-10-19, 11:25, before any crashed run was rerun)

**The crash.** `/mnt/linuxdisk` filled at 11:17, while cell public_s1201 was writing its masked genome.
- The other four cells and both cost runs then failed at their first write. Their output files are empty.
- **Cell A119b_s1201 had already completed: PASS, all six worlds.**

**The fix.** The user approved deleting the masked genomes of the finished rounds (regenerable from each cell's seed) and
the day's HG002 temporary files: 32 files, 86 GB. The six read-only round directories, rounds 9-14, were left alone.

**The rerun.** The crashed runs are rerun once, unchanged, as the rules require: public_s1201 (its partial `masked/` and
`set.log` cleared), A119b_s1202, public_s1202, A119b_s1203, public_s1203, and both cost runs.

## Results (2026-10-19; `/mnt/linuxdisk/tmp/o3_graph/fast/heldout22/`)

**Correctness: PASS.** All six cells pass and all 72 bar checks hold: A119b_s1201 in the first run, the other five in the
Amendment-1 rerun, with no crash.

| system | judged correct | judged errors | off-family ties: correct / errors |
|---|---|---|---|
| integrated v17 | 78,421 | 3,255 | 3,355 / 94 |
| current | 53,763 | 12,600 | 1,652 / 380 |

**C1: incomplete.**
- A119b passes: 953 integrated CPU s against the current method's 1,086.
- The public cost run was stopped by the user at 12:35 (2026-10-19), before it finished, to refocus on the families where
  O2 can be shown to work. The round's cost bar is therefore UNDECIDED.
- Rounds 19-21 passed C1, and the scope step adds about 5 s.

**What it means.** With ambiguity defined by Eichler's within-10 rule plus sequence ties, O2 confined to the ambiguous reads
does at least as well as the current method in every world type, on 26 held-out families that are mostly processed
pseudogenes. Adoption as the default needs the cost bar's public profile.
