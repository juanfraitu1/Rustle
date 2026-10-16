# Integrated O1 -> O2 -> O3 on an eighth handpicked family set, with alleles (2026-10-16)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User question (2026-10-16):** "would fully integrating O1, O2 and O3 help? I think O1 might reduce the possible universe by
only seeding what is supported by RNA reads". The plan was approved: "yes, please do all of that".

**Why it might help (from rounds 1-6).** O2's errors come from three places:
- a copy present in the genome but missing from the catalog (OUTSIDE). Seeding from RNA support could add it: 48 of 56 such
  copies over seven earlier sets have most of their reads' primary alignments on their own locus;
- isoforms annotated only on a paralog (CROSS). A copy's observed structures could add them;
- copies absent from the reference, near-identical to a paralog, and alleles. Integration cannot fix these.

Shrinking the universe to RNA-supported copies could reduce ties. Every catalogued copy was expressed in the earlier worlds,
so this is measured here for the first time (SILENT).

## Methods under test

**O1 seeding** (`o1seed.py WORLD` -> `graph/seeded.json`), per family:
- **family reads:** reads with any alignment (primary or secondary) on a curated copy's locus (+-1 kb);
- **kept:** a curated copy is kept when >= 5 family reads have their primary alignment on it; otherwise it is DROPPED
  (silent). Every primary counts here, so an identical expressed copy, whose reads split at mapq 0, is never dropped;
- **new:** a locus off every family's catalog becomes a new copy when >= 5 family reads have their primary alignment there
  UNIQUELY (mapq >= 1). Adding a copy needs positive evidence; a mapq-0 primary is a coin flip;
- **structures:** each copy's curated isoforms, plus the intron chains seen in >= 3 near-uniquely placed (mapq >= 50) primary
  alignments that this copy overlaps most and that lie within its span (+-100 bp). A new copy has its observed chains only,
  from uniquely placed reads.

**Integrated O2** (`o2graph.py run WORLD integ 2`): padded union O2 (round 3: every structure threaded through every copy, path
ends padded with each copy's genome) on the seeded catalog, over the family reads.

**Integrated O3** (`detect.py`, arm **I**): as round 4's H (pooled, mosaic rule, branch calls, genome check), on the reads the
integrated O2 leaves `unexplained`.

Beside them: the current O2 and CUR; padded union O2 on the curated catalog; H and S23.

## Worlds

All with heterozygous alleles (`--het 0.001`), seed 81:
- INTACT;
- CROSS (`--projection liftover`);
- OUTSIDE;
- REMOVED (special copy hard-masked, own index);
- **SILENT**: a seeded third of each family's copies express nothing; all of them stay in the catalog.

**Seed (detectability pre-check, round-5 review M6).** `precheck.py` aligns each candidate special copy's main isoform to the
genome; its divergence is that of its best alignment away from its own locus. Seed 81 is the smallest seed >= 81 with >= 3
detectable families: SPRR2 0.067, TUBB8 0.086, WASH 0.0116. The others are within delta.

## Families (held out, picked before any run)

The candidates were taken in a fixed order: CTAGE, AQP7, PRSS, SPRR2, FAM74A, TUBB8, DDX11L, WASH, SIGLEC, IFITM, FAM138,
FAM157. CTAGE was a round-6 candidate past its first 8 and was never built into a world. The set is the first 8 with >= 3
copies:

| family | copies / isoforms |
|---|---|
| CTAGE | 4 / 4 |
| AQP7 (AQP7 + pseudogenes) | 4 / 7 |
| SPRR2 | 5 / 9 |
| FAM74A | 4 / 4 |
| TUBB8 (TUBB8 + pseudogenes) | 6 / 11 |
| DDX11L | 10 / 13 |
| WASH | 9 / 22 |
| FAM138 | 8 / 17 |

`/mnt/linuxdisk/tmp/o3_graph/families7.json` (the earlier sets identical to `families6.json`).

## Scoring

- A call is correct when it names the read's own copy. A new seeded copy is mapped to the truth copy whose locus it overlaps
  (scoring only).
- **Reads of an absent copy that the genome explains within delta** (undetectable by O3's definition: its reads' median `de`
  on the masked genome <= delta) cannot be told from the locus that explains them. Placing them there is counted apart
  (`unknowable_placed`), for every method alike, and reported; it is not an error. All other wrong calls are errors.

## Known before running (development)

Development: round 6's families with alleles (INTACT, REMOVED, OUTSIDE; plus new SILENT and CROSS worlds, seed 71), round 5's
families with alleles (`dev_het`), and round 4's families without alleles. Correct / wrong calls (unknowable placements in
brackets):

| set | world | current O2 | padded union O2 | **integrated O2** | ratio of correct calls to the current O2's | uncatalogued reads on their own copy | new copies (on a truth copy) |
|---|---|---|---|---|---|---|---|
| round 6 | SILENT | 1,380 / 0 | 1,509 / 0 | **1,569 / 0** | x1.137 | - | 0 |
| round 6 | CROSS | 2,131 / 0 | 2,318 / 0 | **2,315 / 0** | x1.086 | - | 0 |
| round 6 | INTACT | 1,983 / 0 | 2,104 / 0 | **2,104 / 0** | x1.061 | - | 0 |
| round 6 | REMOVED | 1,683 / 0 (157) | 1,790 / 0 (118) | **1,790 / 0 (165)** | x1.064 | - | 2 (0) |
| round 6 | OUTSIDE | 1,683 / 157 | 1,790 / 118 | **2,074 / 0** | x1.232 | 335 / 480 | 8 (8) |
| round 5 | INTACT | 2,356 / 9 | 2,542 / 0 | **2,542 / 0** | x1.079 | - | 0 |
| round 5 | REMOVED | 2,065 / 0 (226) | 2,256 / 0 (297) | **2,265 / 0 (323)** | x1.097 | - | 1 (0) |
| round 5 | OUTSIDE | 2,065 / 192 | 2,256 / 297 | **2,542 / 0** | x1.231 | 347 / 481 | 7 (7) |
| round 4 | INTACT | 1,855 / 0 | 1,850 / 0 | **1,850 / 0** | x0.997 | - | 0 |
| round 4 | CROSS | 2,325 / 0 | 2,252 / 0 | **2,313 / 0** | x0.995 | - | 0 |
| round 4 | OUTSIDE | 1,434 / 60 | 1,431 / 60 | **1,795 / 0** | x1.252 | 364 / 481 | 7 (7) |
| round 4 | REMOVED | 1,437 / 61 (60) | 1,431 / 0 (60) | **1,431 / 0 (60)** | x0.996 | - | 1 (0) |

O1 dropped all 12 silent copies and none of the 344 expressed copies over these worlds.

Found and fixed on development, before this commit:
- **New-locus support counts unique placements only (mapq > 0).** The first version seeded 10 genomic paralogs in round 6's
  CROSS from mapq-0 coin flips. They competed with the true copies: 2,253 correct against padded union's 2,318.
- **Each alignment counts for the one curated copy it overlaps most, and a structure only from alignments inside that copy's
  span.** CSAG's copies 1 and 2 are 91 bp apart, so their padded windows took each other's reads and structures, and 105
  reads tied.
- **A curated copy's observed structures need near-unique reads (mapq >= 50),** the user's suggestion against fragments of
  ambiguous placement.
- **A new locus still counts any unique placement (mapq >= 1).** With >= 50, real uncatalogued copies near-identical to a
  paralog (ANKRD20A.4, SPDYE.2, TRIM64.2, DEFA.1, FRG1.2, CSAG.3: all reads unique at mapq 1-49) were not seeded, and
  round 5's OUTSIDE recovered only 161 of 481 reads.
- **New copies in REMOVED (1-2 per world) are genomic paralogs that absorbed a removed copy's reads.** All of those copies
  are within delta of the genome, so the placements are counted as unknowable.

O3 (arm I) on these worlds was re-run with this final seeding after the commit; it is reported under Results. With the
previous seeding (unique placements for new loci, no ownership rule), arm I had:
- no false call in any of the 12 worlds;
- detections of 3 of 3 detectable removed copies in round 6's REMOVED and 6 of 6 in round 4's.

## Bars (held out: the 8 families above, seed 81, five worlds)

| id | bar |
|---|---|
| J1 | the integrated O2 has >= 0.98 x the current O2's correct calls, in every world (as rounds 4-5's A1 / C1) |
| J2 | the integrated O2 has <= the current O2's wrong calls, in every world |
| J3 | in OUTSIDE the integrated O2 assigns >= 50% of the uncatalogued in-genome copies' reads to their own (seeded) copy |
| J4 | O1 drops >= 80% of SILENT's silent copies and <= 10% of expressed copies (pooled over the five worlds) |
| J5 | in SILENT the integrated O2 has >= padded union O2's correct calls on the curated catalog |
| J6 | arm I has a false call in <= 1 of the 32 INTACT, CROSS, OUTSIDE and SILENT families |
| J7 | arm I detects every detectable REMOVED copy CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |

## Readings (fixed now)

1. **J1-J5 all pass:** claimed, "on simulated reads of handpicked real families with heterozygous alleles, an O1 that seeds
   the copy catalog from RNA support makes O2 as right (>= 98% of the correct calls) and no more wrong than the current O2 in
   every world: intact,
   paralog-annotated isoforms, an uncatalogued in-genome copy, a removed copy, silent copies. It recovers >= half of the
   uncatalogued copies' reads, drops >= 80% of silent copies and <= 10% of expressed ones, and shrinking the universe does not
   cost correct calls." Any failing bar: not claimed.
2. **J6-J7 pass:** claimed, "integrated O3 detects >= 75% of detectable removed copies with at most one false call in 32
   intact, paralog-annotation, uncatalogued-copy and silent-copy families". J7 not measurable: claimed without the
   removed-copy clause.
3. Each claim stands or falls as a whole.

Scope: simulated reads (curated CAT isoforms, HiFi errors, heterozygous sites at 0.1%, uniform expression) on the real CHM13
genome. Real data's other species, allelic imbalance and expression variation are not modelled.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

## Results (held out, 2026-10-16)

Run: code and prereg `1e513942`, the frozen copy_assign binary, seed 81, `--het 0.001`. Outputs in
`/mnt/linuxdisk/tmp/o3_graph/real/heldout7/` (`integ.json`, `o3_score.json`, `o3_summary.txt`); summaries from `summary.py integ`
and `summary.py o3`.

**O2** (correct / wrong; unknowable placements in brackets):

| world | current O2 | padded union O2 | **integrated O2** | ratio of correct calls to the current O2's |
|---|---|---|---|---|
| INTACT | 2,360 / 0 | 2,476 / 0 | **2,476 / 0** | x1.049 |
| CROSS | 2,697 / 0 | 2,860 / 15 | **2,955 / 11** | x1.096 |
| OUTSIDE | 2,081 / 68 | 2,193 / 193 | **2,533 / 180** | x1.217 |
| REMOVED | 2,081 / 45 (n/a) | 2,193 / 0 (193) | **2,193 / 0 (193)** | x1.054 |
| SILENT | 1,723 / 0 | 1,762 / 0 | **1,834 / 0** | x1.064 |

- **OUTSIDE.** O1 seeded 7 new copies, all on a truth copy, and 348 of the 481 uncatalogued reads went to their own copy. The
  integrated O2's 180 wrong calls are:
  - 60 + 60: two uncatalogued copies that were not seeded (CTAGE.2, DDX11L.9), whose reads went to a near-identical paralog;
  - 60: WASH.7's reads called to WASH.7 under the DDX11L family. The WASH and DDX11L families are adjacent subtelomeric genes,
    and the uncatalogued WASH.7 was seeded in both, since a new locus is off every family's catalog.
- **CROSS.** 11 paralog-structure reads were called to a wrong copy: DDX11L 8, WASH 3.
- **SILENT.** O1 dropped 10 of the 14 silent copies and 0 of 220 expressed copies over the five worlds. The 4 kept silent
  copies had >= 5 primary alignments: mapq-0 coin flips of identical expressed copies. Keeping them is the rule by design (an
  identical expressed copy must never be dropped).
- **REMOVED.** The removed copies' median `de`: TUBB8 0.087, SPRR2 0.076, WASH 0.0124 (detectable); AQP7 0.0077, FAM74A
  0.0093, FAM138 0.0043, CTAGE 0.0023, DDX11L 0.0022 (unknowable).

**O3** (8 families per world):

| | CUR | S23 | H | **I** |
|---|---|---|---|---|
| REMOVED detected (3 detectable) | 1 (TUBB8) | 2 | 2 | **2** (SPRR2, TUBB8) |
| false-call families, INTACT + CROSS + OUTSIDE + SILENT (of 32) | 8 | 0 | 0 | **0** |

WASH (median `de` 0.0124, near delta) is missed by every arm.

| bar | value | |
|---|---|---|
| J1 | lowest world ratio 1.049 | pass |
| J2 | OUTSIDE 180 > 68; CROSS 11 > 0 | **fail** |
| J3 | 348 / 481 = 72% | pass |
| J4 | silent copies dropped 10 / 14 = 71%; expressed 0 / 220 | **fail** |
| J5 | 1,834 vs 1,762 | pass |
| J6 | 0 of 32 | pass |
| J7 | 2 of 3 detectable (0.67), CUR's one included | **fail** |

### Verdict (2026-10-16)

**O2: not claimed (J2 and J4 fail).** On this fresh set the integrated O1 -> O2:
- assigned 5-22% more reads correctly than the current O2 in every world (J1), with no wrong call in INTACT, REMOVED and
  SILENT;
- recovered 72% of the uncatalogued copies' reads (J3);
- shrank the universe without losing correct calls (J5);
- dropped no expressed copy.

But:
- **wrong calls:** it made more than the current O2 where uncatalogued copies are near-identical to catalogued ones and
  were not seeded (CTAGE, DDX11L), and where an uncatalogued copy was seeded into two adjacent families (WASH / DDX11L);
- **silent copies:** identical silent copies keep coin-flip support, so only 71% were dropped.

**O3: not claimed (J7 fails, 2 of 3).** No false call over 32 families. WASH's removed copy, near delta, was missed by
every arm, as was round 3's case near the boundary.

**What the integration showed, across development and this held-out set:**
- **What it fixes:** RNA-supported seeding fixes O2's main failure mode where it can, namely uncatalogued copies with
  uniquely placed reads. On this set wrong calls went from 193 for padded union to 180, with 348 reads now right.
- **What it does not fix:** copies near-identical at mapq 0 (unseedable), and alleles.
- **A bug it exposed:** a seeded locus must belong to one family only. Fixing that is a method change for a new prereg.
