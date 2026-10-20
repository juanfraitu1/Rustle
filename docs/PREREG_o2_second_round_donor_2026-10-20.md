# Pre-registration: the second round on a held-out donor (HG00733, simulated from its own assembly) (2026-10-20)

**User (2026-10-20):**
- After the audit of register 1345: "lets try to implement all of that, I think that is paramount to prove that O2 works".
- For the held-out data: "Donor simulation first (Recommended)".

## What is new against register 1345

1. **One read-crediting rule for every arm** (`bench/o3_graph/o2second.py score-reads`, 7c72b036). Each read contributes the
   records its arm's read pool admits:
   - primaries only;
   - `good`: primaries + secondaries with AS >= 0.98 x the read's best, the shipped pool;
   - `home`: a homed read's records at home, otherwise `good`.

   A record is credited to the copy of the arm's transcript covering at least half of its aligned bases. No tool's own read
   counting is used.
2. **Only the reads O2 can move are scored:** the corrector's domain reads (tied within 10 AS on two or more copies of one
   family), each scored when its true copy is among its tied copies. Each read is:
   - **right:** credited to its true copy alone;
   - **wrong:** credited to another copy only;
   - **echo:** credited to its true copy and another;
   - **uncredited.**
3. **A feasibility gate.** ORACLE, every read homed at its true copy, must beat the others before O2 is judged. Here the truth
   is exact for every read, so ORACLE is real perfect placement.

## The held-out donor

**HG00733** (HPRC year 1, phased hifiasm assemblies, maternal and paternal). It has never been used in this project.

Its reads were simulated from its own sequence (`bench/o3_graph/donorsim.py`):
- **Orthologs.** Each CHM13 copy of the 281 qualifying non-PAR families (`famident/qualify_nonpar.json`; every member pair at
  [0.90, 1.00) identity), with ±20 kb of flank, was aligned to each haplotype (minimap2 -x asm20 -c). The ortholog is the
  longest primary with MAPQ >= 20 through which every isoform lifts. Copies whose orthologs overlap on a haplotype were
  dropped (122 pairs).
  - Result: 1,252 copy-haplotype orthologs; 663 of the 724 copies have one or two. Median identity 0.998 over the whole
    flanked region.
- **Reads.** Each expressed copy's CHM13 isoforms were lifted onto its orthologs. 150 reads per copy, times realworld's
  expression multipliers, were split over haplotypes and isoforms (main 2.0, others 0.5). Each read is the donor
  transcript's own sequence, carrying the `A119b_testis_v4_ref` profile's:
  - retained introns;
  - truncation;
  - end jitter;
  - HiFi errors.

  The donor's alleles come from its genome, not from the profile.
  - Result: 83,700 reads (41,437 maternal, 42,263 paternal); 61 copies have no ortholog and 8 are silent.
- **Truth.** A read's truth is the CHM13 copy whose ortholog produced it.
- **Alignment.** The reads were aligned to CHM13 v2.0 as HG002's were: `-ax splice:hq -uf --secondary=yes -N 100 -p 0.5`,
  giving 271,686 records.

**Seen before this file:**
- the ortholog and build counts above;
- the alignment's flagstat.

No corrector call, assembly or score on these reads has been looked at.

The rules were set on HG002 (development): the corrector (version 2, bfa15f12), the home table (b15b5948) and the scorer
(7c72b036, development 0.805 vs 0.681 on 15 HG002 families).

## Arms (Rustle arms genome-wide on the simulated BAM)

| arm | what | read pool for the score |
|---|---|---|
| **LOOP-O2** | pass 2, every read homed (`loop_home.py o2 --every-read`): the corrector's copy when it decided, else minimap2's when not tied within 10 AS; ties O2 left echo | home |
| **P** | primaries only (`--seed-pool primary`) | primary |
| **GOOD** | pass 1, the shipped assembly | good |
| **ST** | StringTie 3.0.3 `-L` on the primaries | primary |
| **FL** | FLAIR 3.0.0 `bam2bed` + `collapse --trust_ends --generate_map` on the primaries | primary |
| ORACLE | pass 2, every read homed at its true copy (`--every-read --truth`) | home |

The corrector runs on the simulated BAM with its own confident reads as the sample (`o2correct.py donor`).

## Bars

- **G (gate):** ORACLE's right rate is at least 0.05 above the best of P, GOOD, ST and FL.
- **B1:** LOOP-O2's right rate is above P's, GOOD's, ST's and FL's.
- **B2:** LOOP-O2's wrong + echo rate is no higher than P's.
- **B3 (power):** at least 100 scored reads.
- **B4 (the user's bar, "more than 10 families"):** LOOP-O2's right rate is above P's in more than 10 families with at least
  10 scored reads.

**Verdict:**
- **UNDERPOWERED:** B3 fails.
- **INFEASIBLE:** G fails. No verdict on O2.
- **PASS:** B1, B2 and B4 hold.
- **FAIL:** otherwise.

**A PASS would claim:** on a donor never used in development, with real donor-vs-CHM13 copy differences and exact truth,
rebuilding loci after O2 gives the reads minimap2 leaves tied to their true copy's locus. That happens more often than with
minimap2's primaries, the shipped assembly, StringTie or FLAIR, and in more than 10 near-identical families.

**It would not claim:**
- real reads: expression and error profile are simulated, though calibrated on A119b;
- copies the donor has and CHM13 lacks: no reads are simulated from donor loci without a CHM13 ortholog.

**Rules:**
- Run script: `/mnt/linuxdisk/tmp/o3_graph/donorsim/HG00733/run.sh`.
- A bug is fixed and recorded as an amendment before any rerun.

## Results (`donorsim/HG00733/run.log`, `sim/second/score_reads.json`)

**FAIL. The gate holds, so the test was feasible (ORACLE 0.970). B1 and B2 fail; B3 and B4 hold.**

| HG00733, 11,596 scored reads (281 families) | right | wrong | echo | uncredited |
|---|---|---|---|---|
| **LOOP-O2** | 0.477 | 0.079 | **0.428** | 0.016 |
| P (primaries only) | 0.619 | 0.350 | 0 | 0.031 |
| GOOD (pass 1) | 0.077 | 0.023 | 0.889 | 0.012 |
| ST (StringTie) | 0.619 | 0.346 | 0 | 0.036 |
| FL (FLAIR) | 0.630 | 0.349 | 0 | 0.021 |
| ORACLE | 0.970 | 0.002 | 0 | 0.028 |

**Bars:**
- **B1 fails:** 0.477 < P's 0.619.
- **B2 fails:** wrong + echo is 0.507, against P's 0.350.
- **B4 holds:** LOOP-O2 credits more reads right than P in 20 families with at least 10 reads, but fewer in 71.

**Why it fails: the pre-registered home rule.** The corrector decided 6,440 of the 11,970 domain reads. The other 5,530
(4,937 `tie`, 593 `undominated`) have no home, so they keep pass 1's echo at every tied copy, and an echo is not "right".
Primaries only gives each of those reads minimap2's single pick, which is right about 62% of the time.

**On the reads it moved, the corrector is right.** Of the 455 reads O2 took off minimap2's copy:

| arm | right | wrong |
|---|---|---|
| LOOP-O2 | 406 | 47 |
| P, ST, FL | about 41 | 381-413 |

**Descriptive: the corrector itself on this held-out donor** (`o2correct.py score`, printed by the pre-registered run script;
not a bar):

| HG00733 | reads | corrector right | minimap2 right | fixes / breaks |
|---|---|---|---|---|
| all | 11,596 | 0.6695 | 0.6377 | 408 / 40 |
| exact ties | 4,975 | 0.488 | 0.449 | 225 / 29 |
| near ties | 6,621 | 0.806 | 0.780 | 183 / 11 |

52 families have more fixes than breaks; 7 have more breaks than fixes.

**Not claimed:** that loci rebuilt after O2 credit tied reads to their true copy more often than StringTie, FLAIR or
minimap2's primaries. Under this file's rule they do not, because echoed ties are not right.

**Post hoc, not a result:** with minimap2's primary as the home of the ties O2 leaves, LOOP-O2 would be P plus O2's moves:
about 0.619 + (406 - 41) / 11,596 = 0.650 right. That variant needs its own pre-registration on another donor.
