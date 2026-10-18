# Round 15: integrated v11 on fresh families, three seeds per profile, and the cost of full families (2026-10-17)

Pre-registration, written before any world of round 15's families was built: no read of these families has been simulated.
Code: branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.**
- **Round 14 passed (register 1312), but its review called the pass one fragile draw:** A119b passed only 2 of 4 further
  seeds. v10 then passed all 18 spent seed and diploid-donor runs.
- **The goal check found the cost claim unshown on full families.** On the NPIP-sized mainstay set, the integrated stages took
  1,074 s against the current method's 369 s (A119b).
- **Round 15 tests both claims, each with its own pre-registered bar.**

## The two systems

**Current:** copy_assign's O2 and `--psv-clusters` O3, unchanged since round 1.

**Integrated v11** = v10 + two changes. v10 is round 14's v9 plus O3 `--haplotype-check-genome` (f867ec7c). Both changes were
developed on spent runs only:
- **O3 `--branch-linear`** (44798db2): the branch step without the family graph.
  - Candidates are spliced-aligned (minimap2 splice:hq) to the family's loci as linear sequences.
  - A read base is novel only where no locus alignment explains it, which stands for the graph's path switches between copies.
  - vg msga is gone. It was 121 of O3's 158 s on the mainstay CROSS world, and it was not deterministic.
  - **On all 34 spent runs** (16 rounding runs, round 14's 10 seed runs, 8 diploid-donor sets): every verdict is v10's, 14
    runs have identical bars, and 19 gain 1-9 correct reads with equal or fewer errors.
- **`fastset --aln-cache`** (f2e59107, 7839138c): exact alignment caches shared by a set's worlds.
  - O2's path alignments are cached per read sequence.
  - O3's catalog-isoform-to-host alignments are cached by window, query names and sequences.
  - The O2 tables were byte-identical, cold and warm, on the mainstay worlds.

**Measured on spent runs and not adopted**, because each turned a passing development run into a failure:
- **O2 `--twin-abstain`** (de03bf5e): a read that neither of two twins beats by the margin is abstained. Round 10 public
  OUTSIDE fell to 751 correct, below the current method.
- **O2 `--o2-genome guard`** (de03bf5e): a kept path call is abstained when a locus off every seeded copy explains the read as
  well. Round 10 A119b and public fell below the current method's correct count in several worlds; the guard abstains on
  silent array units' ties.

Both cut the O2 proof's OUTSIDE errors on near-identical twins, from 138-171 to 0-2 with the two together. That failure stays
open: an identical, uncatalogued, expressed twin gives O1 no uniquely placed read to seed it.

**Integrated flags:**
```
--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability
--twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome
--branch-linear --aln-cache --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2
```

## Families (held out; drawn by a rule fixed before the pool was built)

**The pool** (`candidates15.py`, e1fd5bb5):
- Every CAT protein-coding gene-name stem with >= 3 genes, none of them matched by any earlier group: 1,409 stems. The list's
  sha1 is 808e8f6d.
- Built with round 12's rule, `families.py build --cap 1000 --link locus --min-cover 0.8` (543 s).
- 105 stems formed a family; 30 have >= 3 copies (eligible), 7 of them >= 6 copies.

**The draw:**
- **Rule** (fixed before the build, applied once): eligible names sorted; seed 15. When >= 4 eligible families have >= 6
  copies, draw 4 of those, then 4 of the rest.
- **Large families drawn:** IFNA, FAM90A, SCYGR, SSU72P, from FAM90A, IFNA, OR2L, OR4A, SCYGR, SSU72P and VN1R.
- **Others drawn:** H2AB, KRT, ANP, BOLA, from 23.

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| IFNA | 7 | 7 | chr9 |
| FAM90A | 57 | 65 | chr8, chr12 |
| SCYGR | 10 | 10 | chr2 |
| SSU72P | 7 | 7 | chr7, chr11 |
| H2AB | 3 | 3 | chrX |
| KRT | 4 | 4 | chr17 |
| ANP | 3 | 6 | chr4, chr12, chr15 |
| BOLA | 3 | 12 | chr16 |

- 94 copies. FAM90A is the family the user parked for O2/O3 on 09-17, as a good near-identical copy-assignment case.
- **Files:** `/mnt/linuxdisk/tmp/o3_graph/families15.json` (this set) and `r15/families15_all.json`, `r15/build.log`.

## Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v2.json` and
  `profiles/public_testis_v2.json`.
- **Seeds:** 501, 502 and 503, fixed with no seed selection. That makes 6 cells (seed x profile) and 36 worlds.

```
F=/mnt/linuxdisk/tmp/o3_graph/families15.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout15
for S in 501 502 503; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r15 $S --local $O/local $V --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done; done
```

## Bars

**Correctness.** Read outcomes are round 14's: per-read detectability less each world's error floor, scored on the veto table
(`--o2 o2.integv.m2.tsv`). Bars per world:
- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

**A cell PASSES when all six of its worlds pass both bars. Round 15 PASSES only when all six cells pass** (72 bar checks).

**Cost (C1).** Measured on the mainstay set: NPIP as one family, TBC1D3, GOLGA8, GOLGA6L, CT47A, XAGE and SPATA31; seed
306; both profiles. These families are spent, but cost is not tuned on them.
- **Setup:** `fastset.py run` with v11's flags, the existing mainstay local reference, nothing else running.
- **Bar:** v11's integrated wall time, summed over the six worlds (O1 + O2 + O2_genome + O3), is at most the current method's
  summed wall time, in each profile.
- **Reported, no bar:** v10's run on the same set, and every step's wall time and peak RSS (`timing.json`).

**Rules:**
- Nothing is changed after the first run.
- A bug is recorded as an amendment, with its reason, before any rerun.
- A run that crashes is rerun once, and the crash is recorded.

## What a PASS would and would not claim

**A correctness PASS** would claim that v11 does at least as well as the current method in every world type, on 8 held-out
families, under both profiles and three seeds.

**A cost PASS** would claim that, on NPIP-sized families, the integrated O1 -> O2 -> O3 costs no more wall time than the
current method on the same machine and settings.

**Neither would claim:**
- anything on real reads;
- that near-identical uncatalogued twins are handled (see above).

## Amendment 1 (2026-10-17 17:25, before any rerun)

**What happened:** /mnt/linuxdisk filled up during the batch (0 bytes free). Cell A119b 501 finished; every later run, cells
2-6 and the v10 cost runs, exited with code 120 within a second, at its first write. C1's two v11 runs had already finished.

**The space freed:** 37 per-run masked genome copies (`masked/masked.fa`, 3.1 GB each, 108 GB in all) of runs completed and
recorded on 10-15 and 10-16, none referenced by an active directory.
- Each copy is the CHM13 reference with a few intervals hard-masked. It can be regenerated exactly by `realworld.py mask` from
  the run's families, set and seed.
- The list of deleted copies is `docs/deleted_masked_2026-10-17.txt`.

**Under the rule above:** the crashed runs are rerun once, unchanged. Cell public 501's partial directory is removed first:
its mask step died writing.

## Results (2026-10-17; cells 2-6 rerun once under Amendment 1)

**Verdict: FAIL.** 3 of the 6 cells pass; the bar was all six. **C1 FAILS too.**

| cell | verdict | failing worlds (integrated errors vs current) | integrated / current wall |
|---|---|---|---|
| A119b 501 | PASS | | 855 / 417 s |
| public 501 | PASS | | 662 / 374 s |
| A119b 502 | FAIL | CROSS 3 vs 0, INTACT 2 vs 0 | 944 / 468 s |
| public 502 | PASS | | 771 / 411 s |
| A119b 503 | FAIL | CROSS 5 vs 0, INTACT 5 vs 0, SILENT 5 vs 0 | 1,117 / 462 s |
| public 503 | FAIL | CROSS 2 vs 0, INTACT 1 vs 0, SILENT 1 vs 0 | 934 / 411 s |

**B1 holds in all 36 worlds**: integrated correct exceeds current everywhere.

**Every failing bar is a few errors where the current method makes none, and every one of those errors is in FAM90A** (57
near-identical copies). By the O2 proof's oracle:
- **Most are undecidable reads.** The read's alignment cost ties at the true copy and a sibling, yet path O2 assigned it by the
  margin, to a copy with a richer path set:
  - A119b 503: five reads of copies 20, 21, 29 and 41, all assigned to FAM90A.45;
  - public 503: copy 29 -> 39, copy 12 -> 11.
- **One read is decidable** but was called against the oracle (A119b 502: FAM90A.45 -> FAM90A.2).
- **This is a calibration gap in large near-identical families.** Path structure asymmetries that union projection does not
  remove decide the margin, and the margin is fixed (2 edits) whatever the number of rivals. copy_assign's PSV certificate
  divides alpha by the number of rivals, and it abstains on these reads.
- **Elsewhere the integrated pipeline stays far ahead:** OUTSIDE 4-66 errors against 101-263; BRANCH 13-60 against 44-195.

**C1, the cost on the mainstay set:** v11's integrated stages took 416 s (A119b) and 418 s (public), against the current
method's 320 s and 333 s.
- **That fails C1** by 96 s and 85 s.
- **It is 2.4-2.9x faster than v9's integrated stages** (1,200 s and 987 s): O3 825 -> 149 s and 763 -> 98 s, without the vg
  graph.
- **O2 is now the larger cost** (228 / 279 s). The cache helps only worlds whose path set an earlier world made.
- **The thread counts were not equal:** the current method's minimap2 ran 8 threads under `--heavy-current`, and each integrated
  step had 2. Both runs' bars are identical to v9's (PASS).
- **On round 15's families,** where FAM90A's 57 copies make O2 large, the integrated stages took 1.7-2.4x the current method.
