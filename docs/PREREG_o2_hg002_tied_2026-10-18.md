# Pre-registration: does O2 put AS-tied reads on the right copy? (HG002, the donor's own genome as truth, 2026-10-18)

**User (2026-10-18):**
- "can we prove that we O2 and O3 work? ... can we prove we can effectively assign reads to copies and point out possible
  cases in which there is a missing copy that would explain the PSVs?"
- "minimap2 already assigns them if they have mapq>0 the problem is the multimappers that have equal AS so we don't know to
  which copy do they truly belong and the distinction between primary and secondary is meaningless (coin toss)"
- "Ok we need to prove it in that way, any dataset it just needs to assign reads to copies based on PSVs for tied AS reads"
- "lets ensure that O2 can only work for tied AS reads, everything else is already resolved by minimap, then O3 is just an
  extension in which if we notice the PSVs nothing in the genome explains them"

**The split this tests:**
- **Resolved:** a read whose best-AS alignments all lie on one copy belongs to minimap2.
- **O2:** an AS-tied read. O2 may choose only among the copies it ties on, or abstain.
- **O3:** reads whose PSVs no genome copy explains. An O3 call overrides any placement.

## Why HG002

- **Real reads with their donor's genome.** PacBio Kinnex FLNC of HG002 (37,200,651 reads; Revio, 2023; the BAM header
  names no tissue) and HG002's diploid T2T assembly (Q100 v1.1, maternal and paternal), with the Q100 chain to CHM13 v2.0.
- **The 10-17 proof was simulation only** and said so (`docs/PREREG_o2_proof_2026-10-17.md`: "never claim O2 on real
  reads from this").
- **Held out.** No O2 rule was developed or tuned on HG002 reads. They were used once before, on 2026-09-29, to test a DNA
  presence gate. That gate failed and is not in the pipeline (register 1180).

## The population: reads minimap2 can only toss a coin for

1. **Candidates.** Every FLNC read with a minimap2 `map-hifi` hit to a transcript of the adopted multigene pool, at >= 0.90
   identity over >= 200 bp: 633,939 reads.
   - The pool is `seqfamilies.py build --rule multigene` (register 1328): 439 families, 1,685 copies, 3,583 transcripts
     (`r19seq/families19cal_full2_all.json`, set `seq19`).
2. **CHM13.** The candidates mapped to the whole of CHM13 v2.0 with
   `minimap2 -a -x splice:hq -uf --secondary=yes -N 100 -p 0.5`.
3. **Tied.** A read is tied when its best AS is shared by alignments on two or more copies of one family. A copy is matched
   by exon overlap.
   - Reads tied across two families are counted, not scored.
   - These reads' primary alignment is minimap2's coin toss among the tied copies.
   - **Result:** 23,424 reads tied within one family, plus 543 across two.
     - **Size of the tie:** 23,266 tie on 2 copies, 141 on 3, 15 on 4, 1 on 5 and 1 on 7.
     - **Primary MAPQ:** 0 for 23,131; 1-29 for 82; 30 or more for 211.
     - **Families:** 72 have tied reads.
4. **Pseudoautosomal families are a separate stratum.** A family is pseudoautosomal when every copy lies in CHM13
   v2.0's PAR1 or PAR2 on chrX or chrY (`hg002truth.pseudoautosomal`): 16 of the 439 families, 12 of them with tied reads.
   - **Why apart:** in this male donor, the chrX and chrY copies are the maternal and paternal alleles of one recombining
     locus. Telling them apart is phasing, not copy assignment. Also, CHM13 v2.0's chrY is HG002's own chrY.
   - **Size:** they hold 9,410 of the 23,424 tied reads. SLC25A6 alone has 7,338.
   - **The verdict is decided on the other 60 families with tied reads (14,014 tied reads).** The pseudoautosomal
     stratum runs separately, scored the same way, and is reported with no bar.

## The truth: the donor's own haplotypes (`bench/o3_graph/hg002truth.py`, frozen with this document)

- **Mapping.** Each tied read is mapped to each HG002 haplotype with the same settings.
- **Lift.** Its best-scoring placements, over both haplotypes, are lifted to CHM13 through the chain. Each placement is put
  on the copy its lifted blocks overlap most.
- **Classes:**
  - `unique`: every best placement lands on one copy, and CHM13 tied the read on it. **Truth = that copy.** Two
    haplotypes landing on one copy are the copy's two alleles.
  - `unique_elsewhere`: as `unique`, but on a copy of the family that CHM13 did not tie the read on. **Truth = that copy.**
  - `own_tie`: the best placements land on two or more copies. The donor's copies tie too, so the read is undecidable.
  - `elsewhere`: the best placements lift onto CHM13, but onto no copy of the family. The read comes from another CHM13
    locus.
  - `no_counterpart`: no best placement lifts onto CHM13 at all. This is a donor copy that CHM13 lacks: O3's case.
  - `poor_match`: the donor's best placement has NM > max(2, 0.3% of its aligned bases), arm 4's floor
    (`bench/o2_truth/haplotype_truth.py`). Not scored.
  - `unmapped`: no alignment to either haplotype. Not scored.
- **The donor margin** of a `unique` read is reported but not used: the best AS minus the best AS of a placement off the
  truth copy.

**Counts (the truth step only; no assigner's output had been seen):**

| stratum | tied | `unique` | `unique_elsewhere` | `own_tie` | `elsewhere` | `no_counterpart` | `poor_match` |
|---|---|---|---|---|---|---|---|
| non-pseudoautosomal | 14,014 | 466 | 118 | 9,948 | 948 | 31 | 2,503 |
| pseudoautosomal | 9,410 | 624 | 0 | 7,151 | 0 | 0 | 1,635 |

- **No read was unmapped.**
- **`poor_match`:** NM quartiles 6, 9 and 13 (the floor is kept as arm 4 set it).
- **Donor margin of the `unique` reads:** quartiles 10, 15 and 20 AS points, about 2-4 mismatches. 44 reads have no
  placement off their copy at all.
- **So most AS ties are ties in the donor too.** 71% of the non-pseudoautosomal tied reads are `own_tie`, as register 1180
  found for PAR and EIF3C reads.

## The world and the run

- **Families:** every non-pseudoautosomal `seq19` family with at least one scored tied read, set `hg002`
  (`hg002truth.py world truth.tsv ... seq19 1 world --par exclude`): 54 families, 280 copies.
  > AC0064533 AC1275021 AGAP4 ANKRD36C BMS1P1 BOLA2B CBWD2 CCT6P3 CFC1B CSAG2 CTAGE15 DMRTC1 DUSP22 EIF3C F8A2 FAM156B
  > FCGR2A GOLGA2P10 GOLGA6L10 GOLGA8A GTF2H2C GTF2IRD2 H2AB2 H2AC18 H4C14 HSFX2 HSFY1 HSPA6 IKBKG LRRC37A2 MAGED4B
  > MSTRG3150 NBPF1 NCF1 NOMO1 NPIPB2 NUTM2B OR4G2P PDPR PMS2P11 RASA4B RBMY1A1 RBMY2BP RPL6 SERF1B SLC35E2B SMG1 SMN2
  > SPDYE1 SULT1A1 TBC1D3B ULK4P2 USP17L11 ZC3H11A
- **Reads:** every candidate whose CHM13 primary lies within 1 kb of a copy of those families, as sequenced
  (`realreads.py pool`): 83,146 reads.
  - **Truth reads in the world:** 582 of the 584.
  - **Scored tied reads in the world:** 10,839 of 11,511. A read whose coin-toss primary lies off every copy is not
    pooled. The scorer reports these reads as `not_in_world`.
- **The world's truth** is the tied reads with a donor copy (`unique` and `unique_elsewhere`). Every other read is in the
  world but not scored.
- **One INTACT world of real reads**, with no simulation. The seed is 1 and unused.
- **Pipeline:** v14 as in round 19 (register 1321, not changed since), with `--step-timeout 21600`. The timeout does not
  change outputs.
  ```
  P=/mnt/linuxdisk/tmp/o3_graph/hg002proof; O=/mnt/linuxdisk/tmp/o3_graph/fast/hg002
  V14="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome --branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue --short-path-align --ref-copies --sim-het-draw --heavy-current --step-timeout 21600 --het 0 --jobs 2 --threads 2"
  python3 bench/o3_graph/fastset.py run $O/nonpar $P/world/families.json hg002 1 --worlds INTACT --real-reads $P/world/real --local $O/local_nonpar $V14
  python3 bench/o3_graph/o2proof.py score $O/nonpar_o2proof.json $O/nonpar/INTACT --tied $P/tied.tsv
  python3 bench/o3_graph/hg002truth.py score $P/truth.tsv $P/tied.tsv $O/nonpar/INTACT $O/nonpar_score.json
  python3 bench/o3_graph/hg002truth.py bars $O/nonpar_o2proof.json $O/nonpar_score.json $O/nonpar_bars.json
  ```
  - **The pseudoautosomal stratum** runs the same way after the premise test, with `--par only`, into `$O/par`.
  - **Known risk.** v14's genome veto trusts a primary at MAPQ >= 30 even when a secondary ties its AS (round 20,
    register 1327). 211 of the 23,424 tied reads have such a primary. The scoped systems keep such a call when it is on a tied copy.
- **Systems under test, each judged on its own:**
  - **Integrated O2:** the O1 -> O2 table `graph/o2.integv.m2.tsv`.
  - **Current O2:** `copy_assign` (`/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`, built 2026-10-15 from a tree with the
    exact-tie fix 7ac557e7 and the gate fix d5d982a3).
  - **The user's scoped O2,** `integrated_tied` and `current_tied`. These are the same calls, but a tied read's call stands
    only on a copy it ties on, and is otherwise an abstention; an untied read takes minimap2's copy. On this population
    every read is tied, so only the first rule applies. It is computed from the same run, from `tied.tsv`.
- **Baselines, with no bar:**
  - `mm2`: always the primary's copy (the coin toss).
  - `mm2_q1`: the primary's copy when MAPQ >= 1, else abstain.

## Scoring

- **O2 level** (`o2proof.py score ... --tied tied.tsv`; the oracle is unchanged since b3b5df59, and `--tied` only adds
  the scoped systems): the world's truth reads.
  - Each read is classed by what its own sequence can tell on the world's `E0.bam`: `decidable`, `weak`, `undecidable`,
    `off_family`, `true_unaligned`.
  - Each assigner's outcome is `correct` (the donor's copy), `wrong` (another copy, another family, or a locus that is no
    copy) or `abstain`. O3 is ignored at this level.
- **By donor class** (`hg002truth.py score`): every scored tied read in the world, for each assigner, with O3.
  - An absent-copy call takes precedence over a placement.
  - `unique` / `unique_elsewhere`: right only on the donor's copy. An O3 call on it is an error, since its copy is in
    CHM13.
  - `own_tie`: a placement on a tied copy is a guess; anywhere else, or an O3 call, it is an error.
  - `no_counterpart`: any placement is an error; an O3 call is right.
  - `elsewhere`: a placement is an error; an O3 call is reported apart.

## The premise: minimap2 resolves the untied reads

This test is run after the O2 run, with the same truth code and settings.
- **Reads:** the candidates whose best-AS CHM13 alignments all lie on one copy (`hg002truth.py untied`).
- **Truth:** mapped to both haplotypes (`truth_rest.sh`).
- **P1:** minimap2's copy is the donor's (`unique`) on >= 0.99 of the untied reads with a donor copy (`unique` +
  `unique_elsewhere`), with a floor of 100 reads.
- **Reported:**
  - `own_tie`: the donor's copies tie where CHM13's do not, so minimap2's choice cannot be checked.
  - `no_counterpart` and `elsewhere`: O3's domain.

## Bars

**O2 (each system), on the world's truth reads (`unique` + `unique_elsewhere`), O2 placements:**
- **T1, precise:** wrong / placed <= 0.02.
- **T2, better than the coin toss:** the upper Wilson 95% bound of the system's wrong / placed is below the lower bound of
  `mm2`'s wrong / all.
- **T3, calibrated:** placements on `own_tie` reads <= 0.05 of them, from the by-class score.
- **T4, uses the evidence:** correct on `decidable` reads >= 0.90 of them.

**Power floors:**
- T1 and T2 need >= 100 placed reads.
- T3 needs >= 100 `own_tie` reads.
- T4 needs >= 100 `decidable` reads.
- A bar below its floor is UNDERPOWERED, not failed.

**Validity:** `mm2` must be wrong on >= 25% of the truth reads. Otherwise the population is not a coin toss, and T2 is
void.

**Systems judged:** `integrated`, `current`, `integrated_tied`, `current_tied` (`hg002truth.py bars`, frozen with this
document: `FLOOR, T1_MAX, T3_MAX, T4_MIN, T5_MAX, VALID_MIN = 100, 0.02, 0.05, 0.90, 0.05, 0.25`).

**O2 verdict (per system):**
- **PASS:** T1 and T2 are powered and hold, and T3 and T4 hold or are underpowered.
- **FAIL:** any powered bar fails.
- **UNDERPOWERED:** otherwise.

**O3 (each system), from the by-class score:**
- **T5, safe:** O3 calls on `unique`, `unique_elsewhere` and `own_tie` reads <= 0.05 of them.
- **Reported:** O3 calls on `no_counterpart` reads, the "missing copy that would explain the PSVs" (recall, no bar), and
  on `elsewhere` reads.

**Reported with no bar:**
- Each assigner's correct / unique-truth reads (recall): the share of the coin-toss reads a system rescues.
- The `weak`, `undecidable`, `off_family` and `true_unaligned` classes.
- The donor margin, MAPQ and family breakdowns.

**Rules:**
- Nothing is changed after the run.
- A bug found in the truth, the scorers or the pipeline is recorded as an amendment, with the reason, before any rerun.

## What a PASS would and would not claim

**A PASS would claim** that on real HG002 reads whose CHM13 alignment score ties between copies, the system's placements
agree with the donor's own genome at >= 98%, better than minimap2's coin toss. It abstains on reads the donor's genome
cannot resolve, and places the reads its sequence decides.

**It would not claim:**
- **Families HG002's library barely expresses.** The testis-only Y families give 5 tied reads (HSFY1 3, RBMY1A1 1,
  RBMY2BP 1).
- **Power where CHM13 has no evidence.** Register 1182: a tied read often carries the donor's own alleles, which CHM13's
  copies lack equally. A CHM13-based assigner can only abstain on it.
- **Gorilla, or other donors.**

## Results: the non-pseudoautosomal world (2026-10-19, run once as pre-registered; `/mnt/linuxdisk/tmp/o3_graph/fast/hg002/`)

The run took 48 min (map_union 747 s, current 726 s, O2 650 s, O3 543 s), peak RSS 14.4 GB, with no failure. The code was
unchanged from 7e865cba.

**Verdict: UNDERPOWERED for all four systems.**
- **No system placed 100 of the 584 truth reads,** so T1 and T2 are below their floor. Their point estimates are in the
  table below.
- **Validity holds:** `mm2` is wrong on 341 of 582 (58.6%).
- **T3 holds for every system** (powered).
- **T4 is underpowered:** only 8 of the truth reads are `decidable`.
- **T5 holds for the integrated system and FAILS for the current one.**

| system | placed | wrong (T1) | T2: upper / coin toss's lower | own_tie placed (T3) | O3 calls on present-copy reads (T5) |
|---|---|---|---|---|---|
| integrated | 57 | 7 (0.123) | 0.232 / 0.544 | 96 / 9,607 (0.010) | 40 / 10,189 (0.004) |
| integrated_tied | 57 | 7 (0.123) | 0.232 / 0.544 | 96 / 9,607 (0.010) | 40 / 10,189 (0.004) |
| current | 41 | 5 (0.122) | 0.255 / 0.544 | 40 / 9,607 (0.004) | **6,999 / 10,189 (0.687)** |
| current_tied | 39 | 3 (0.077) | 0.203 / 0.544 | 37 / 9,607 (0.004) | **6,999 / 10,189 (0.687)** |
| `mm2`, always | 582 | 341 (0.586) | | | |
| `mm2` at MAPQ >= 1 | 22 | 12 (0.545) | | | |

**Reported breakdowns (no bar):**
- **Most AS ties are ties in the donor too, and both O2s abstain on them.** 9,607 of the world's scored tied reads are
  `own_tie`; the systems place 0.4-1.0% of them.
- **CHM13 holds almost no evidence for the rest.** The oracle finds 560 of the 584 truth reads `undecidable`, 14
  `off_family`, 8 `decidable` and 2 `true_unaligned`. These reads carry the donor's own alleles, which both CHM13 copies
  lack (register 1182).
- **MAPQ means nothing on a tie.** A tied read's primary at MAPQ >= 1 is wrong 12 times in 22, as the user said.
- **Placements by donor margin** (AS points between the truth copy and the best placement off it):

  | donor margin | reads | integrated: right / wrong | current: right / wrong | `mm2`: right / wrong |
  |---|---|---|---|---|
  | no placement off the truth copy | 44 | 42 / 0 | 0 / 0 | 22 / 22 |
  | <= 5 (about one site) | 169 | 1 / 5 | 3 / 0 | 83 / 86 |
  | 6-10 | 192 | 2 / 0 | 4 / 1 | 66 / 126 |
  | > 10 | 179 | 5 / 2 | 29 / 4 | 70 / 107 |

- **Each system's right placements come mostly from one family.**
  - The integrated system's 50 are 42 IKBKG (the donor has no second placement for them), 3 NBPF1, 3 GTF2H2C, 1 RPL6
    and 1 ZC3H11A.
  - Where the donor has a rival copy, the integrated system is 8 right and 7 wrong, no better than the coin toss.
  - The current system's 36 are 26 NBPF1, then 2 each of SPDYE1, GTF2H2C and GOLGA8A, and 1 each of EIF3C, RPL6,
    ZC3H11A and AC1275021.
  - Its PSV tie-break is right 36 times in 41, and 36 in 39 when scoped.
- **The integrated system's 7 wrong placements:**
  - Five are reads the donor separates by 4-5 AS points, about one site: EIF3C four times, SULT1A1 once.
  - Two are AC1275021 reads whose CHM13 truth copy aligns 228-492 bases worse than its rival. The donor's copy has a
    structure CHM13's lacks.
- **O3 on donor copies CHM13 lacks.** The world holds 30 `no_counterpart` reads: NPIPB2 12, GTF2H2C 9, SPDYE1 5, USP17L11
  3 and BMS1P1 1.
  - The integrated O3 flags 10 of them, all NPIPB2: a missing NPIPB copy. It places 7 and abstains on 13.
  - The current O3 flags none of them and places 16. Instead it flags 6,907 EIF3C reads, which tie in the donor too:
    allelic variation called a missing copy (register 1181).
- **The round-8 end-to-end bars on the truth reads** (fastset's own): integrated 50 correct and 12 errors; current 33 and
  37.

**What it means:**
- **On real AS-tied reads, CHM13's PSVs rarely carry the answer, and a CHM13-based O2 can mostly only abstain.** Both do:
  the calibration bar T3 holds on real data.
- **Where they place, they beat the coin toss (88% against 42%) but are far from the 98% bar,** on too few reads to
  decide. The current system's PSV tie-break is the better of the two on reads with a real rival (36/41).
- **The answer lies in the donor's own copies.** The reads that decide them are the tied reads' own alleles, so the lever
  is to learn the sample's copies from its uniquely placed reads and assign against those. That is a design question; it
  was not tested here.

## Results: the premise (2026-10-19, `hg002proof/truth_rest.sh`, as registered)

- **Reads:** 599,678 of the 633,939 candidates are resolved by minimap2, meaning their best-AS CHM13 alignments all lie on
  one copy.
- **Truth classes:** `unique` 215,514, `unique_elsewhere` 6,342, `own_tie` 3,112, `elsewhere` 1,353, `no_counterpart`
  259,889, `poor_match` 113,468.

**P1: FAILS.** minimap2's copy is the donor's on 215,514 of 221,856 untied reads with a donor copy (0.971, against the
0.99 bar; powered).

**Reported breakdowns (no bar):**

| families | `unique` | `unique_elsewhere` | P1 |
|---|---|---|---|
| not pseudoautosomal, no chrM copy | 207,998 | 4,481 | 0.979 |
| pseudoautosomal | 7,515 | 1,861 | 0.802 |
| with a chrM copy | 1 | 0 | |

- **`no_counterpart` here is a truth artifact.** The Q100 chain has no chrM chain. So the donor's chrM placements of
  258,334 mitochondrial-transcript reads (families MTND6, MTRNR2L8 and MTCO2, whose other copies are nuclear) do not lift.
  Outside them, 1,555 reads are `no_counterpart`. The tied reads were barely affected (31 `no_counterpart`).
- **Where minimap2 misses.** Its misses concentrate in copy-number-variable segmental duplications.
  - The largest are P2RY8 1,149 (pseudoautosomal), GPAT2 1,045, SMN2 606, AKAP17A 473 (pseudoautosomal), THOC3 416,
    FCGR2A 400, WASHC1 344, CCDC144A 210, TBC1D3B 199 and PLCXD1 178 (pseudoautosomal).
  - Donor margins of the misses: quartiles 5, 10 and 192. Half are 1-2 sites, the rest structural.
- **What it means.** "Everything else is already resolved by minimap2" holds for about 98% of the untied reads in
  ordinary families, not 99%.
  - The rest are reads whose own copy carries donor alleles that make another CHM13 copy fit better. The user's split
    gives them to minimap2 and gets them wrong.
  - Only a check of the sample's own variation, the O3 extension "PSVs nothing in the genome explains", can notice them.

## Results: the pseudoautosomal stratum (2026-10-19, `fast/hg002/run_par.sh`; reported, no bar)

- **The world:** 12 families, 24 copies, 22,302 reads and 624 truth reads (all `unique`). The run took 6 min.
- **Neither system attempts phasing.** The integrated system places 1 truth read (right) and the current one none. Both
  place 0 of the 7,151 `own_tie` reads.
- **O3** calls 3 (integrated) and 7 (current) of the 7,775 present-copy reads.
- **The coin toss is wrong on 287 of 624.**
- **CHM13 cannot separate them.** The oracle finds 620 of the 624 truth reads `undecidable` on CHM13: their X and Y alleles
  differ by donor variants that neither CHM13 copy carries.
