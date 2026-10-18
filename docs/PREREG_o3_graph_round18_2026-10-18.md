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
