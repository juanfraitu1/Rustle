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
