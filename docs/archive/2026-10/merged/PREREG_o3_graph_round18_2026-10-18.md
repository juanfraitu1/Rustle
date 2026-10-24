# Round 18: integrated v14 (copies the catalogue lacks, found in the genome) on fresh name-free families, three seeds per profile, and the cost in CPU seconds (2026-10-18)

Pre-registration, written before any world of round 18's families was built: no read of these families has been simulated.
Code: branch `o3-fast`, frozen at the commit that adds this file.

**Why.**
- **Round 17 failed every cell on OUTSIDE's B2** (register 1320). When a pair's twin is uncatalogued, O1 seeds nothing for it
  (no read is uniquely placed), and O2 gave the twin's reads to the one remaining copy.
- **The user chose to make it O3's job** (2026-10-18). Before O2, the catalogue is completed from the genome: a locus that holds a
  catalogued copy by the family rule becomes a new copy.
- **On round 17's six cells, now spent:** 0 -> 5 cells pass, and OUTSIDE errors fall from 367-781 to 1-159 (the current method:
  109-191).
  - B1 holds in all 36 worlds, at 1-8 correct reads a world.
  - The sixth cell fails REMOVED by one read (2 errors against 1), as v13 did.

## The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v14** = v13 (register 1319) + O1 `--ref-copies` (24bd5cde, 581570cf).
- Each kept curated copy's representative transcript is spliced-aligned to the world's homology-closed local reference (the
  masked one in REMOVED).
- An off-catalogue locus that holds it by the family rule becomes a new copy: identity >= 0.90 over >= 0.8 of it, and its exon
  chain carried by the locus, a junction at an unaligned end needing no intron.
- These copies get no `--twins` annotation prior. A read-seeded new copy that the search confirms loses that prior too.

**Integrated flags (`$V14`):**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue
--short-path-align --ref-copies --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

## Families (held out; name-free; drawn by a rule fixed before the pool was rebuilt)

**The pool:** `seqfamilies.py build` (51f9b3a6), round 17's rule with the amended structure test (581570cf: a junction at an
unaligned end needs no intron).
- 4,508 of 4,738 sequence links pass structure, forming 256 families.
- **Spent** are all families touching a gene of any earlier simulated set, round 17's included. That leaves 98 unspent families,
  210 copies.

**The draw** (fixed in 51f9b3a6, before the rebuild): the unspent families sorted by name, `random.Random(18).shuffle`, taken in
that order until the copy total reaches 90, a family skipped if it would take the total past 150.
- **43 families, 91 copies;** copies per family: {2: 39, 3: 3, 4: 1}.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families18.json` (set `r18`, sha1 e689764bcf06), and `r18seq/families18seq_all.json`.

| family | copies | genes | chromosomes |
|---|---|---|---|
| S17ACSM2A | 2 | ACSM2A, ACSM2B | chr16 |
| S17AKAP17A | 2 | AKAP17A | chrX, chrY |
| S17AKR1B10 | 2 | AKR1B10, AKR1B15 | chr7 |
| S17AKR1C1 | 2 | AKR1C1, AKR1C2 | chr10 |
| S17AKR7A3 | 2 | AKR7A2, AKR7A3 | chr1 |
| S17ALG10 | 2 | ALG10, ALG10B | chr12 |
| S17ASMT | 2 | ASMT | chrX, chrY |
| S17ATAD3B | 2 | ATAD3A, ATAD3B | chr1 |
| S17CCDC74B | 2 | CCDC74A, CCDC74B | chr2 |
| S17CCL4L2 | 2 | CCL4L2 | chr17 |
| S17CD99 | 2 | CD99 | chrX, chrY |
| S17CELA3B | 2 | CELA3A, CELA3B | chr1 |
| S17CRLF2 | 2 | CRLF2 | chrX, chrY |
| S17CYP4A11 | 2 | CYP4A11, CYP4A22 | chr1 |
| S17EIF3C | 2 | EIF3C, EIF3CL | chr16 |
| S17ERVV1 | 2 | ERVV-1, ERVV-2 | chr19 |
| S17FAM86B2 | 2 | FAM86B1, FAM86B2 | chr8 |
| S17FFAR3 | 2 | FFAR3, GPR42 | chr19 |
| S17GTF2IRD2 | 2 | GTF2IRD2, GTF2IRD2B | chr7 |
| S17HTN3 | 2 | HTN1, HTN3 | chr4 |
| S17IFITM2 | 2 | IFITM2, IFITM3 | chr11 |
| S17IL9R | 2 | IL9R | chrX, chrY |
| S17MAGEB6B | 2 | MAGEB6, MAGEB6B | chrX |
| S17MAGED4B | 2 | MAGED4, MAGED4B | chrX |
| S17MTRNR2L6 | 2 | MTRNR2L6 | chr7, chr9 |
| S17OCM | 2 | OCM, OCM2 | chr7 |
| S17OR10H2 | 3 | OR10H1, OR10H2, OR10H5 | chr19 |
| S17OR4F21 | 2 | OR4F16, OR4F21 | chr1, chr5 |
| S17OR52I2 | 2 | OR52I1, OR52I2 | chr11 |
| S17OR8B2 | 2 | OR8B2, OR8B3 | chr11 |
| S17PLGLB1 | 3 | PLG, PLGLB1, PLGLB2 | chr2, chr6 |
| S17PRAMEF1 | 4 | PRAMEF1, PRAMEF13, PRAMEF14, PRAMEF2 | chr1 |
| S17PRAMEF18 | 2 | PRAMEF18, PRAMEF19 | chr1 |
| S17RAB6C | 2 | RAB6C, RAB6D | chr2 |
| S17RFPL1 | 3 | RFPL1, RFPL2, RFPL3 | chr22 |
| S17ROPN1 | 2 | ROPN1, ROPN1B | chr3 |
| S17SYT15 | 2 | FO681492.1, SYT15 | chr10 |
| S17TCP11X1 | 2 | TCP11X1, TCP11X2 | chrX |
| S17TRIM16 | 2 | TRIM16, TRIM16L | chr17 |
| S17TRIM43B | 2 | TRIM43, TRIM43B | chr2 |
| S17TUBB2A | 2 | TUBB2A, TUBB2B | chr6 |
| S17USP18 | 2 | USP18, USP41 | chr22 |
| S17ZBED1 | 2 | ZBED1 | chrX, chrY |

## Worlds and runs

- **Worlds:** six world types under `profiles/A119b_testis_v2.json` and `profiles/public_testis_v2.json`.
- **Seeds:** 801, 802 and 803, fixed with no seed selection: 6 cells, 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families18.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout18
for S in 801 802 803; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r18 $S --local $O/local $V14 --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done; done
```

## Bars

**Correctness (as rounds 15-17).** Per-read detectability less each world's error floor, on the veto table
(`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 18 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1), in CPU seconds, as rounds 16-17,** on the mainstay set (seed 306, both profiles), nothing else running:
```
M=/mnt/linuxdisk/tmp/o3_graph/fast/r14dev/mainstay_npip306/local
for P in A119b public; do
  python3 bench/o3_graph/fastset.py run /mnt/linuxdisk/tmp/o3_graph/fast/r18cost/v14/$P \
    /mnt/linuxdisk/tmp/o3_graph/families_mainstay_npip.json mainstay 306 --local $M $V14 \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done
```
- **Bar:** v14's integrated CPU seconds (O1 + O2 + O2_genome + O3, summed over the six worlds) at most the current method's, in
  each profile.
- **Reported, no bar:** wall times, peak RSS, the cells' CPU seconds, and the genome copies O1 found per world.

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

## What a PASS would and would not claim

**A correctness PASS** would claim that v14 does at least as well as the current method in every world type, on 43
held-out name-free families, under both profiles and three seeds.

**A cost PASS** would claim no more CPU than the current method on NPIP-sized families.

**Neither would claim:**
- anything on real reads;
- anything about arrays of FAM90A's or DUX4's size beyond the spent runs.

The catalogue still starts from the annotation's copies; the genome search adds the copies it lacks.

## Results (2026-10-18 10:57-11:43; nothing changed after the first run, no crash)

**Verdict: FAIL -- 4 of 6 cells pass; A119b 803 and public 803 fail. C1 PASSES** in both profiles.

| cell | failing worlds (integrated vs current: correct / errors) |
|---|---|
| A119b 803 | CROSS 4,052 / 0 vs 4,093 / 33; INTACT 3,726 / 0 vs 3,771 / 34 |
| public 803 | CROSS 3,703 / 0 vs 3,799 / 42; INTACT 3,391 / 0 vs 3,498 / 29; SILENT 2,140 / 0 vs 2,226 / 21 |

- **B2 holds in all 36 worlds.** Over the 36 worlds v14 makes 130,722 correct calls against the current method's 86,384, and
  3,782 errors against 10,431.
- **Round 17's failure is gone on fresh families:** OUTSIDE errors are 0-11 against the current method's 69-198, and REMOVED
  errors 1-17 against 26-85.
- **Every failing bar is B1, all in seed 803.** In those worlds v14 makes no errors and the current method makes 21-42.

**The mechanism, post hoc: one family, S17OR4F21 (OR4F21 on chr1, OR4F16 on chr5).**
- **O1's genome search adds five copies to it:** OR4F2P (chr11), OR4F7P (chr6), OR4F8P (chr19), OR4F28P (chr15) and
  FO082796.2 (chr3), 98.8-99.4% identical to OR4F16's transcript. All are unprocessed pseudogenes in CAT, and they are 5 of the 13
  genome copies O1 finds in each INTACT world.
- **The simulator expresses only the catalogue's protein-coding copies,** so these five never have a read.
- **Most reads the current method gets right and v14 does not are this family's:** 102 of 111 in A119b 803 INTACT, 159 of 172 in
  public 803 INTACT.
  - Tied with a pseudogene copy: 64 and 101.
  - Inside the fixed 2-edit margin of one (`thin`, mostly one substitution): 32 and 57.
- **The current method runs with `--no-as-tied-only`** (O3 needs it). That turns off the shipped tie-outside-catalogue rule
  (2026-09-09), so copy_assign compares only the catalogued copies.
  - Read s0000757's primary alignment is on OR4F2P, at the same AS as on OR4F21 (424).
  - copy_assign assigns it to OR4F21 at p 1e-18, against OR4F16 alone.
- **Seed 803 drew OR4F21 at 253-286 reads a world,** against 28-29 in seeds 801-802, where the same mechanism costs a few reads
  (B1 margins of +4 to +20 in INTACT).
- **In the real A119b testis library** (CHM13 v2.0 primary alignments), OR4F21 has 33 reads and OR4F16 26.
  - The five pseudogenes have 41: OR4F7P 18 (17 at MAPQ >= 1), OR4F8P 7, FO082796.2 7, OR4F28P 6, OR4F2P 3.
  - The simulation's silent pseudogenes are not what real testis shows.

**Two post-hoc recounts on these runs** (spent; not a verdict):
- **B1 on decidable reads only** (a catalogued read is undecidable when one of its alignments elsewhere in E0.bam scores at least
  its best at its true copy): 4 of the 5 failing bars pass. Public 803 INTACT stays short, 3,383 against 3,392, by the
  one-substitution reads.
- **Thin reads decided by the calibrated significance gate alone,** in place of the fixed 2-edit margin (offline replay, veto and
  structure check not replayed):
  - +100 to +566 correct a world, with no added error in INTACT or SILENT and 2-7 in CROSS.
  - OUTSIDE is not scored by the replay, which does not map new copies to their truth copies.
  - Its genome-step form, `--significance-margin`, was rejected in register 1319: one error per world where the current method
    makes none.

**C1, the cost on the mainstay set:** v14 used 1,014 CPU seconds (A119b) and 820 (public), against the current method's 1,044 and
945 (-3% and -13%).
- O1's genome search costs 141-147 CPU seconds of these.
- Reported: wall 565 s and 473 s against 310 s and 283 s, at the flags' threads; peak RSS 2.42 and 2.35 GB.
- Round 18's cells used 1,336 CPU seconds in all against 1,254. The set's one-time local-reference build peaks at 14.4 GB.
