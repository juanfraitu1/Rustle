# O1 -> O2 -> O3 integrated: the overall picture (2026-10-16)

User (2026-10-16): "lets look at the overall picture first, probably it is better to look at O1, O2 and O3 integrated since
that is the end goal".

This is a post hoc summary across every simulated world built in rounds 1-7 (`docs/PREREG_o3_graph_*`, register 1297-1305).
No new claim is made here; the next step proposed below is a pre-registered end-to-end test.

## The integrated pipeline (frozen at round 7, `1e513942`)

1. **O1, RNA-supported seeding** (`o1seed.py`):
   - keep each curated copy that has primary-alignment support, and drop the silent ones;
   - add off-catalog loci with >= 5 uniquely placed reads;
   - give each copy the intron chains of near-unique reads (mapq >= 50) that it owns.
2. **O2, padded union** (`o2graph.py run W integ`): every family structure is threaded through every seeded copy, with path
   ends padded with the copy's own genome. Each read is assigned at margin 2, or abstains, or is `unexplained` beyond delta.
3. **O3, arm I** (`detect.py`):
   - clusters what O2 leaves unexplained, pooled over the family, with the mosaic rule;
   - adds branch calls from the family's vg graph;
   - keeps a call only when the genome does not explain its reads.

The comparison is the current method: copy_assign's O2 plus its `--psv-clusters` O3 on the curated catalog.

## End to end: every read's fate (`summary.py e2e`)

Each read is classed by the copy it came from:
- catalogued;
- uncatalogued but in the genome;
- absent and detectable;
- absent and within delta of the genome, so unknowable.

Its fate under each system is one of: assigned to its own copy (correct), to another copy (error), inside an absent-copy call
(correct only for an absent copy), or left unassigned. Pooled over 43 worlds, 8 family sets:

| world type | current: correct / error / unassigned | integrated: correct / error / unassigned | correct | errors |
|---|---|---|---|---|
| INTACT | 17,990 / 483 / 5,727 | 19,377 / 10 / 4,813 | x1.077 | -98% |
| CROSS (isoforms annotated on a paralog) | 18,470 / 508 / 6,954 | 20,002 / 22 / 5,908 | x1.083 | -96% |
| OUTSIDE (uncatalogued copy in the genome) | 15,084 / 2,762 / 6,354 | 19,550 / 304 / 4,346 | x1.296 | -89% |
| REMOVED (copy absent from the reference) | 16,116 / 623 / 6,128 | 18,065 / 26 / 4,457 | x1.121 | -96% |
| BRANCH (absent copy with a novel exon) | 11,418 / 1,955 / 5,064 | 15,008 / 83 / 3,346 | x1.314 | -96% |
| SILENT (a third of the copies unexpressed) | 3,085 / 46 / 956 | 3,403 / 0 / 684 | x1.103 | -100% |
| **all** | **82,163 / 6,377 / 31,183** | **95,405 / 445 / 23,554** | **x1.161** | **-93%** |

Unknowable placements (reads of absent copies within delta of the genome) were 1,333 for the current method and 1,652 for the
integrated pipeline.

- **Per world.** The integrated pipeline has at least as many correct fates and at most as many errors in 38 of the 43
  worlds. The five exceptions:
  - round 3's INTACT and CROSS: 1 and 2 fewer correct;
  - round 4's INTACT and CROSS: 5 and 12 fewer correct;
  - round 6's REMOVED: 227 more correct, 1 more error.
- **Not tuned on.** The 30 worlds never used to develop the integrated pipeline (the original development set, rounds 1, 2,
  3 and 5 without alleles, and round 7's held-out set) give x1.178 correct and 92% fewer errors. The 13 worlds it was
  developed on give x1.122 and 98% fewer.

## Where the errors come from

All 43 worlds, by the same rules as `summary.py e2e` (each error counted once):

| error source | current | integrated |
|---|---|---|
| a present copy's read inside an absent-copy call (O3) | 3,412 | 61 |
| a detectable absent copy's read placed on a present copy | 1,695 | 82 (all novel-exon, BRANCH) |
| an uncatalogued copy's read placed on another copy | 1,213 | 291 |
| a catalogued copy's read placed on another copy (O2) | 57 | 11 |
| **all** | **6,377** | **445** |

There were 109,493 reads from catalogued copies and 4,335 from uncatalogued ones.

## What each objective contributes

- **O3 makes the largest change.** More than half of the current method's errors (3,412 of 6,377) are absent-copy calls
  on present copies' reads.
  - Arm I made 117 true and 6 false calls over the 43 worlds.
  - All 6 false calls are one cluster: NPIPA.5's reads hosted on the sister family NPIPB in round 1's held-out set, just
    above delta (median `de` 0.0097, known since round 4).
  - The earlier arms made the same number: H made 5 false calls (the same NPIPB cluster); K made 5 (one WASH cluster in each
    round-7 world, 5 reads each), which arm I does not make.
  - Branch calls, the genome check, the mosaic rule and family pooling catch detectable absent copies' reads: those placed
    on a present copy fall from 1,695 to 82.
- **O1 fixes the catalog.** Uncatalogued copies' reads placed on another copy fall from 1,213 to 291. With seeding they go to
  their own copy: 72% of them on round 7's held-out set (348 of 481).
  - O1 dropped 22 of 26 silent copies and none of the 1,740 expressed ones.
  - Its observed structures supply isoforms annotated on a paralog (CROSS).
- **O2 was already precise on catalogued copies** (57 wrong of 109,493 reads); the graph O2 brings that to 11 and abstains
  less.

## What no step fixes

- **Near-identical copies whose reads tie at mapq 0.** They cannot be seeded, told apart, or dropped (CTAGE.2 and DDX11L.9 on
  round 7's set). Reads alone cannot resolve them.
- **Alleles versus copies below delta** (round 6). The same signature needs the individual's own haplotypes. HG002 has them;
  its Kinnex reads are being aligned on the cluster.
- **Removed copies near delta** (WASH at 0.0124, round 7), at the boundary of detectability.
- **A sister family's copy hosted at delta** (NPIPA.5's reads on NPIPB, median `de` 0.0097): its reads look like an absent
  copy. This is the only source of the integrated pipeline's false calls.
- **One bug, known:** an uncatalogued locus can be seeded into two adjacent families (WASH and DDX11L). A seeded locus must
  belong to one family. *Fixed in 1e3960fe; see the update below.*

## Why the per-objective bars failed although the whole works

The pre-registered O2 bars compared O2 alone with the current O2 alone:
- they did not count the current method's O3 calls on present copies, its largest error source;
- they counted a near-identical uncatalogued copy's reads placed on a paralog against O2, although no method can separate
  them.

End to end, those reads move from "error" (current: called absent, or placed on the wrong copy) to "correct" (seeded) or to
the same unresolvable placement.

## Next

1. **A pre-registered end-to-end test on a fresh family set.** Bars on read fates per world type: correct >= current,
   errors <= current. Include the one-family seeding fix and BRANCH as well as the round-7 worlds.
2. **Real data: HG002 Kinnex.** Truth for each read comes from its unique alignment to HG002's own maternal and paternal
   haplotypes (`/mnt/linuxdisk/tmp/rustle_figures_dev/o2_presence/hg002/hg002v1.1.fa`). Copies present in HG002 but absent
   from CHM13 give real O3 targets. It needs the Kinnex reads aligned to CHM13 (the pipeline) and to the HG002 assembly
   (the truth).

## Update after round 8 (2026-10-16)

**Two fixes since the tables above were written:**
- the one-family seeding rule (1e3960fe). Set 1's NPIPA/NPIPB locus had also been seeded twice, in all five of its worlds;
- the e2e scorer: a call on a seeded copy that matches no truth copy is now a placement, not an abstention (round 8 review I3).

**The same 43 worlds, re-scored:**
- integrated: 95,418 correct and 432 errors (was 95,405 and 445):
  - set 1's OUTSIDE: 3 more correct and 13 fewer errors;
  - set 1's CROSS: 10 more correct;
- REMOVED: 76 reads moved from abstain to unknowable (unknowable 1,728, abstain 4,381). No correct or error count changed;
- the integrated pipeline's error sources: uncatalogued 281 (was 291), absent-copy calls on present copies 58 (was 61),
  novel-exon 82, catalogued 11;
- arm I's false calls: 5, all of them the NPIPB cluster. The fix removed one.

**The worlds of one set share most of their reads:**
- INTACT, OUTSIDE and REMOVED simulate the same reads;
- CROSS and BRANCH contain all of INTACT's reads;
- SILENT's reads are a subset.

Pooled totals are therefore not independent counts; each world is the unit of comparison.

**Round 8, the pre-registered end-to-end test on a ninth set: FAIL by one read** (CROSS: 1 error against 0). B1 held in all six
worlds. See `docs/PREREG_o3_graph_round8_2026-10-16.md`.

## Update after rounds 9-11 (2026-10-16, evening): the integrated pipeline v4 passes on held-out full families

| round | pipeline | reads | verdict | what decided it |
|---|---|---|---|---|
| 9 (r1307) | v2 | calibrated v1 profiles | FAIL | 1-3.5% fewer correct catalogued reads (O2 `unexplained` gate on realistic reads); O1 dropping expressed copies |
| 10 (r1308) | v3 | calibrated v2 profiles | FAIL | one false O3 call (a host-window artefact in a 45-copy array); a new copy within delta of a catalogued one |
| 11 (r1309) | v4 | calibrated v2 profiles | **PASS, 12/12 worlds** | also under round 10's scoring, at copy level, and without H2AC/H2BC. The review found the INTACT/CROSS/SILENT margins are the current method's O3 artefact; with O3 removed from both, O2 is at parity under A119b (+3 to +5) |

**v4, end to end** (`fastset.py ... --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
--o2-genome fallback --per-read-detectability --twins --error-floor`):
1. **O1** (`o1seed.py --twins`):
   - **Seeding:** the RNA-supported catalog, keeping a curated copy with >= 5 primaries or one owned uniquely placed read.
   - **Undecidable copies:** a curated copy within delta of a kept one is never dropped.
   - **Twins:** a new locus within delta of a curated copy is its twin.
2. **O2** (`o2graph.py run W integ 2 --distance gapcomp --union-exact --error-floor`, then `integl 2 --error-floor`):
   - **Paths:** every isoform structure is threaded through every seeded copy. Structures still project when the aligner
     drops a terminal micro-exon, and stay apart unless their junctions match exactly.
   - **Distance:** gap-compressed (indel events, clipped ends aligned against the path).
   - **Decision:** at margin 2. A twin wins a read from its curated copy only by the margin.
   - **Genome fallback:** a read the paths leave unexplained is placed on its own genome alignments (mismatches + indel
     events). Novel sequence >= 30 bp refuses, and the best alignment off the seeded copies keeps the read unexplained.
3. **O3** (`detect.py --member-genome-check --member-attribution --error-floor --arms I`):
   - **Clustering:** O2's unexplained reads, pooled over the family, with PSV clustering and the mosaic rule.
   - **Branch calls:** novel segments on the family's vg graph.
   - **Attribution:** a call keeps only the reads its consensus explains better than every catalogued copy that competes
     for them.
   - **Genome check:** per read, less the error floor.
4. **Delta and the error floor.** Delta (0.00958) is an error-free divergence: p99 between KB3781's haplotype assemblies.
   Every read-level comparison with it subtracts the world's sequencing-error floor first (median `de` of MAPQ-60 primaries).

**What remains (round 11's review and its own report):**
- **O3's read hosting still comes from copy_assign's run.** `--o1-hosts` removes that dependence, but in development it lost
  O3 recall in arrays; a tiered version is being measured.
- **Errors left in round 11:**
  - uncatalogued reads placed on a paralog (OUTSIDE 44 and 37, against the current method's 367 and 297);
  - BRANCH reads lacking the novel exon, which count as detectable by construction (118 and 119, against 313 and 297).
- **Not yet tested:** real reads end to end. Rounds 9-11 are simulations calibrated on two real libraries.

## Update after rounds 12-13 (2026-10-17): v5 and v5' fail

| round | pipeline | families | verdict | what decided it |
|---|---|---|---|---|
| 12 (r1310) | v5 = v4 + O1-hosted O3 | 8, 104 copies (80-copy DUX4) | FAIL, both profiles | O2: retained introns absorbed by a paralog's path (gap compression), a read DNA-identical at two copies, 2 margin-2 coincidences in DUX4 |
| 13 (r1311) | v5' = v5 + branch genome-novel prefilter | 8, 89 copies (43-copy OR7E) | FAIL (A119b); public 6/6 | one false O3 call: short haplotype-2 reads of a present copy carrying its two donor alleles |

- **The integrated pipeline is now independent of the current method:** O3's hosts come from O1.
- **It is cheap on full families:** 216-223 s for the integrated steps of a six-world set of 89 copies, with O3's longest
  world 72 s at 0.17 GB.
- **What remains is two confounds a sequence-only rule cannot always settle:**
  - an allele of a present copy against a low-divergence missing copy (O3);
  - a read's structure against a paralog's annotation (O2).
- **Development since round 13** (`docs/O1O2O3_FAST_REALISTIC_2026-10-16.md` section 14):
  - O3 `--allele-discount`: donor alleles, learned from the reads O2 assigns, leave the genome check;
  - O2 `--o2-genome veto`: a path call that the read's own genome alignments contradict by the margin is overruled.
