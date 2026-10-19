# Round 21: O2 only on ambiguous reads, by Eichler's rule (within 10 AS), on fresh multigene families (2026-10-19)

This pre-registration was written before any world of round 21's families was built; no read of these families has been
simulated. Code: branch `o3-fast`, frozen at the commit that adds this file.

## Why

**The user (2026-10-19):**
- "go ahead and use only for tied AS reads"
- "only AS tied reads are ambiguous and need provenance; O2 on mapq>0 reads is tautological since they are already placed
  at their right placement so only multimappers with tied AS where the distinction between primary and secondary is
  meaningless need O2"
- "according to the Eichler paper, their definition of a well placed multimapper is that where the AS is within 10"

**Eichler's rule, as implemented here before** (register 935): a read is placed on its best-AS alignment iff no other
alignment scores within 10 AS units of it.
- Under splice:hq a substitution moves AS by 5 (+1 match, -4 mismatch), so 10 is two sequence differences.
- On HG002's untied reads (development, register 1330 data), minimap2's copy at a margin of 5-9 was wrong 46.5% of the time,
  and at a margin of 10 or more 1.6% of the time.

**Development (all spent):**
- **HG002 (real reads, the donor's own copies as truth):** this system and v14 agree.
  - On 49,057 untied reads with a donor copy: 44,064 right and 597 wrong, against v14's 44,021 and 599.
  - On the 582 tied reads: 49 right and 7 wrong, against 50 and 7.
  - Equal-AS-only scoping instead was wrong on 3.2%.
- **Round 19's 36 worlds, rescored with the scorer below:** 4 of 6 cells pass; v14 passes all 6.
  - Both failures are B2: 1-2 PABPC1L2B reads in each of three worlds. Their own copy carries the paralog's alleles, so
    minimap2 places them at a margin of 12 on the other copy, where v14's paths abstained.
  - The user kept B2 strict and accepted these reads as the split's known cost.
- **Two fixes found in development, both in the code:**
  - A catalogued copy O1 dropped as silent ties like a copy, not like a locus off the catalogue (register 1315's lesson).
  - Reads tied within 10 AS with a locus on no copy of the world are scored apart (user's decision): the simulation never
    expresses that locus.

## The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v16** = v14 (round 19, register 1324) + `--o2-scope tied --as-slack 9` (`o2graph.py run W intega`). Per read and
family, the alignments in `E0.bam` within 9 AS of the read's best are its loci:
- **No alignment:** abstain.
- **A locus off the family's copies:** abstain.
- **One seeded copy:** that copy, minimap2's placement. The exception is when its alignment does not explain the read (a
  novel insertion or clip, or divergence beyond DELTA plus the world's error floor); then abstain, leaving it to O3.
- **The silent copy alone:** abstain.
- **Two or more copies (a tie):** v14's call when it is one of them; otherwise abstain.

O3 is unchanged: it reads `o2.integ`, and its calls take precedence.

**Integrated flags (`$V16`):**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --ref-copies --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
--o2-scope tied --as-slack 9
```

**The simulation:** as rounds 19-20 (profiles v4).

## Families (held out; drawn by the fixed rule)

**The pool:** the adopted multigene rule (`seqfamilies.py build --rule multigene`, register 1328), with all components
(`r19seq/families19cal_full2_all.json`, 439 families).
- **Spent:** every family touching a gene of an earlier simulated set (26 files, `ls families*.json | grep -v _all.json`),
  or of the HG002 worlds (`hg002proof/world*/families.json`, seen in development).
- That leaves 135 unspent families with 393 copies. File: `r21seq/families21seq_all.json`.

**The draw:** the rule of rounds 17-20 with seed 21: names sorted, `random.Random(21).shuffle`, taken in order until the
copy total reaches 90, skipping a family that would pass 150.
- **37 families, 90 copies:** 45 coding and 45 non-coding (33 transcribed unprocessed, 8 unprocessed, 3 transcribed
  processed, 1 processed).
- **Copies per family:** {2: 28, 3: 4, 4: 4, 6: 1}.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families21.json` (set `r21`, sha1 0e002f85636e).

| family | copies | coding genes | non-coding members | chromosomes |
|---|---|---|---|---|
| S19AC0928653 | 6 | - | AC092865.3 (tr. unproc.), AC079781.4 (tr. unproc.), AC068587.8 (tr. unproc.), AC068587.8 (tr. unproc.), AC105233.3 (tr. unproc.), AC068587.6 (tr. unproc.) | chr7, chr8, chr12 |
| S19AC1335553 | 3 | - | AC133555.3 (tr. unproc.), AC106782.1 (tr. unproc.), AC106782.1 (tr. unproc.) | chr16 |
| S19AC2455951 | 4 | - | AC245595.1 (tr. unproc.), AC245595.1 (tr. unproc.), AC027612.2 (tr. unproc.), BAGE2 (tr. unproc.) | chr1, chr2, chr21 |
| S19ADAMTS7P1 | 2 | - | ADAMTS7P1 (tr. unproc.), ADAMTS7P4 (tr. unproc.) | chr15 |
| S19ANAPC1 | 2 | ANAPC1 | ANAPC1P4 (tr. unproc.) | chr2 |
| S19AQP12B | 2 | AQP12B, AQP12A | - | chr2 |
| S19BMS1P8 | 4 | - | BMS1P8 (tr. unproc.), BMS1P8 (tr. unproc.), AC026273.1 (tr. unproc.), AC026273.1 (tr. unproc.) | chr2, chr16 |
| S19CCT8L2 | 2 | CCT8L2 | CCT8L1P (proc.) | chr7, chr22 |
| S19CD177 | 2 | CD177 | CD177P1 (unproc.) | chr19 |
| S19CSF2RA | 2 | CSF2RA | - | chrX, chrY |
| S19CYP21A2 | 2 | CYP21A2 | CYP21A1P (tr. unproc.) | chr6 |
| S19CYP2C19 | 2 | CYP2C19, CYP2C9 | - | chr10 |
| S19CYP3A7 | 2 | CYP3A7, CYP3A4 | - | chr7 |
| S19DHX40 | 2 | DHX40 | DHX40P1 (tr. unproc.) | chr17 |
| S19EOLA1 | 2 | EOLA1, EOLA2 | - | chrX |
| S19MALL | 2 | MALL | AC140479.5 (tr. unproc.) | chr2 |
| S19MSTRG2518 | 4 | MSTRG.2518 | USP32P1 (tr. unproc.), USP32P3 (tr. unproc.), AC053481.2 (tr. unproc.) | chr17 |
| S19MSTRG676 | 2 | MSTRG.676 | PARGP1 (tr. unproc.) | chr10 |
| S19NCOR1P4 | 3 | - | NCOR1P4 (tr. proc.), NCOR1P4 (tr. proc.), NCOR1P4 (tr. proc.) | chr13, chr14, chr21 |
| S19NLGN4X | 2 | NLGN4X, NLGN4Y | - | chrX, chrY |
| S19NXF2 | 2 | NXF2, NXF2B | - | chrX |
| S19OR1E1 | 2 | OR1E1, OR1E2 | - | chr17 |
| S19OR4A15 | 2 | OR4A15 | OR4A21P (unproc.) | chr11 |
| S19OR4C46 | 2 | OR4C46 | OR4C10P (unproc.) | chr11 |
| S19OR52L1 | 2 | OR52L1 | OR52L2P (unproc.) | chr11 |
| S19OR5F1 | 2 | OR5F1 | OR5F2P (unproc.) | chr11 |
| S19PDE4DIP | 3 | PDE4DIP | AC241952.1 (tr. unproc.), AC245297.1 (tr. unproc.) | chr1 |
| S19PROS1 | 2 | PROS1 | PROS2P (unproc.) | chr3 |
| S19RLN2 | 2 | RLN2, RLN1 | - | chr9 |
| S19S100A7A | 2 | S100A7A, S100A7 | - | chr1 |
| S19SORD | 2 | SORD | SORD2P (tr. unproc.) | chr15 |
| S19TBC1D26 | 2 | TBC1D26, TBC1D28 | - | chr17 |
| S19TEKT4P2 | 2 | - | TEKT4P2 (tr. unproc.), TEKT4P2 (tr. unproc.) | chr13, chr15 |
| S19TMEM191C | 2 | TMEM191C | TMEM191A (tr. unproc.) | chr22 |
| S19TPTE2 | 3 | TPTE2, TPTE | TPTE2P4 (unproc.) | chrY, chr13, chr21 |
| S19ZNF257 | 4 | ZNF257, ZNF724, ZNF675, ZNF273 | - | chr7, chr19 |
| S19ZNF285 | 2 | ZNF285 | AC018688.1 (unproc.) | chr19, chr20 |

## Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v4.json` and
  `profiles/public_testis_v4.json`.
- **Seeds:** 1101, 1102 and 1103, with no seed selection. That makes 6 cells and 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families21.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout21
for S in 1101 1102 1103; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r21 $S --local $O/local $V16 --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done; done
```

## Bars

**Correctness.** As rounds 15-20: per-read detectability less each world's error floor, on the scoped table
(`--o2 o2.intega.m2.tsv`). Reads tied within 9 AS with a locus on no copy of the world (`--off-family-ties 9`) are the
class `off_family_tie`, left out of both systems' sums. fastset's own e2e step does this.
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**Verdict.** A cell PASSES when all six of its worlds pass both bars, and round 21 PASSES only when all six cells pass (72 bar
checks).

**Cost (C1), in CPU seconds,** as rounds 16-20, on the mainstay set (seed 306, both profiles), with nothing else running:
```
M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
for P in A119b public; do
  python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r21cost/v16/$P \
    /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V16 \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v4.json
done
```
- **Bar:** v16's integrated CPU seconds (O1 + O2 + O2_genome + O2_scope + O3, summed over the six worlds) at most the
  current method's, in each profile.

**Reported, no bar:**
- v14 on the same worlds and scoring (`o2.integv.m2.tsv`, from the same runs).
- The `off_family_tie` reads' fates.
- The errors of reads minimap2 resolved (`as_resolved`) apart from O2's.
- Wall times and peak RSS.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

## What a PASS would and would not claim

**A correctness PASS** would claim that O2 confined to ambiguous reads (within 10 AS, Eichler's rule), with minimap2's
placement for the rest, does at least as well as the current method in every world type. It would hold on 37 held-out
multigene families, half of them with non-coding members.

**It would not claim:**
- **Anything about reads tied with a locus off every copy.** These are scored apart.
- **Accuracy where a donor's copy carries its paralog's alleles at two or more sites.** minimap2 places those reads on the
  paralog. Development shows such reads exist, and real data does too (register 1330).
- **Real-read accuracy beyond the HG002 development numbers above.**
