# Pre-registration: does O2 work? (held-out near-identical families, 2026-10-17)

**User (2026-10-17):** "In a previous session we found that GH, USP17L, TSPY, RBMY, BPY2, CDY, OR4F and LCE1; six of these
eight have near-identical copies; can we use these to prove that O2 works?" Then, on the proposed design: "yes please".

**Why those eight cannot be the proof:** they have been development set 9 since round 8, so the integrated O2's rules were
tuned on them. Five of the eight test only abstention (copies at divergence 0 to 0.0005). They are used here as development
evidence. The claim is decided once, on the held-out families below.

## The question

The integrated O2 (the O1 -> O2 table `graph/o2.integv.m2.tsv`, pipeline v10) assigns each read of a family to one of the
family's copies or abstains. **"O2 works" means three things, each with a bar:**
- **precise:** it rarely assigns a read to the wrong copy;
- **calibrated:** it does not guess on reads whose sequence cannot tell their copy apart;
- **powerful:** it assigns the reads whose sequence does tell.

## The oracle (`bench/o3_graph/o2proof.py`, frozen with this document)

Every read of a catalogued copy (truth kind `present` or `silent`; absent and uncatalogued copies are O3's) is classed by
what its own sequence can tell.
- **Costs:**
  - The read's alignments in `E0.bam` (primary and secondary; supplementary pieces left out) are placed on the copy whose
    locus each overlaps most. Off every copy, an alignment goes to its 10-kb genome bin.
  - A locus costs its best alignment's NM plus the read bases that alignment leaves unaligned.
- **Margin** = the best family rival's cost minus the true copy's. Family rivals are the other copies of the read's family in
  the reference genome, catalogued or not.
- **Classes:**
  - `off_family`: a locus outside the family costs no more than the true copy. Which loci form the family is O1's
    question, so these reads are reported only.
  - `undecidable`: margin <= 0. A family rival explains the read at least as well.
  - `weak`: margin 1. Reported only.
  - `decidable`: margin >= 2 (O2's own margin), or no family rival aligned.
  - `true_unaligned`: the aligner never placed the read on its true copy. Reported only.
- **Outcomes:** `correct` (the read's own copy), `wrong` (another copy, another family, or a locus that is no copy),
  `abstain`.
- **Baselines** (reported, no bar):
  - `mm2`: always the primary alignment's copy.
  - `mm2_q1`: the primary at MAPQ >= 1, else abstain.
  - `current`: copy_assign's O2.

**Development evidence** (spent: the 16 v10 development runs, sets 8-10 and rounds 10-13, both profiles;
`/mnt/linuxdisk/tmp/o2proof/*.json`):
- **Wrong / assigned:** 0 in 13 runs. The others are 1/13,150 (set 8 A119b), 1/13,456 (set 8 public) and 16/12,909
  (round 12 A119b: the DUX4 / KRTAP2 coincidences no sequence rule separates).
- **Wrong on undecidable reads:** at most 11 of 14,203 (round 12 A119b). `mm2` errs on 23-89% of them.
- **Correct on decidable reads:** at least 0.991 in every run.
- **Set 9 alone (the eight families above):** 0 wrong of 5,023 / 4,611 assigned; decidable 4,099/4,137 and 3,730/3,776.
  776 / 788 reads that tie with a locus outside the family are placed on their true copy.

## Families (held out; picked before any build)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `o2proof`, c647219b). No family is in any
earlier list; VCY, HSFY and PRY were scored once before, for copy_assign's O2 only, by another simulator (2026-10-12).
- **Build rule:** round 12's, `families.py build families_o2p_all.json --cap 1000 --link locus --min-cover 0.8`, every
  group (111 s).
- **The set** is the first 8 families with >= 2 copies, in list order. HBG gave no family.

| family | copies | isoforms | genes | chromosome |
|---|---|---|---|---|
| HBA | 2 | 8 | HBA1, HBA2 | chr16 |
| HSPA1 | 2 | 2 | HSPA1A, HSPA1B | chr6 |
| TUBA3 | 3 | 5 | TUBA3C, TUBA3D, TUBA3E | chr2, chr13 |
| VCY | 2 | 2 | VCY, VCY1B | chrY |
| TPSB | 2 | 7 | TPSAB1, TPSB2 | chr16 |
| C4 | 2 | 8 | C4A, C4B | chr6 |
| HSFY | 2 | 8 | HSFY1, HSFY2 | chrY |
| SERPINB34 | 2 | 6 | SERPINB3, SERPINB4 | chr18 |

17 copies. Files: `/mnt/linuxdisk/tmp/o3_graph/families_o2p.json` (this set), `families_o2p_all.json`, `families_o2p.log`.

## Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v2.json` and
  `profiles/public_testis_v2.json`, with a diploid donor (`--sim-het-draw`).
- **Seeds:** 401, 402, 403, fixed in advance with no seed selection. That makes 6 cells (seed x profile) and 36 worlds.
- **Pipeline:** v10, round 14's flags plus `--haplotype-check-genome`.

```
F=/mnt/linuxdisk/tmp/o3_graph/families_o2p.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/o2proof
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2"
for S in 401 402 403; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F o2proof $S --local $O/local $V --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done; done
python3 bench/o3_graph/o2proof.py score $O/<cell>.json $O/<cell>/{INTACT,CROSS,OUTSIDE,REMOVED,BRANCH,SILENT}
```

## Bars (integrated O2), per profile pooled over its 3 seeds and 6 worlds, and in each of the 6 cells

- **B1, precise:** wrong / assigned <= 0.005.
- **B2, calibrated:** wrong on `undecidable` reads <= 0.01 of the `undecidable` reads.
- **B3, powerful:** correct on `decidable` reads >= 0.95 of the `decidable` reads.

**Power floors:**
- A pooled profile needs >= 500 `undecidable` and >= 500 `decidable` reads.
- A cell needs >= 100 of each.
- A bar below its floor is UNDERPOWERED for that profile or cell, not failed.

**Validity:**
- B2 is a test only if the set is hard: `mm2` must err on >= 25% of the pooled `undecidable` reads. Otherwise B2 is
  UNDERPOWERED (an easy set).
- B1 needs >= 500 assigned reads per pooled profile.

**Verdict:**
- **PASS:** B1, B2 and B3 hold in both pooled profiles and in every powered cell.
- **FAIL:** any powered bar fails.
- **UNDERPOWERED:** otherwise.

**Reported with no bar:** copy_assign's O2 on the same metrics, `mm2` and `mm2_q1`, the `weak`, `off_family` and
`true_unaligned` classes, and each world and family.

**Rules:**
- Nothing is changed after the first run.
- A bug found in the oracle or the pipeline is recorded as an amendment, with the reason, before any rerun.
- Development sets are never rescored to choose between versions of the oracle.

## What a PASS would and would not claim

**A PASS would claim** that on simulated reads of 8 held-out near-identical families, under two read profiles and three
seeds, the integrated O2 is precise, calibrated and powerful in the senses above.

**It would not claim:**
- **Accuracy on real reads.** For the Y families there is no matched donor: the complete Y assemblies come from cell lines,
  and these genes are expressed in testis only.
- **Anything about off-family ties.** In the simulation an unlisted paralog never expresses; in real data it can.
- **Anything about the families' scope.** That is O1's question.

## Results (2026-10-17, run once as pre-registered; `/mnt/linuxdisk/tmp/o3_graph/fast/o2proof/`)

The six cells ran in 11 min in total, with no failure. The oracle was unchanged from b3b5df59.

**Verdict: PASS.** B1, B2 and B3 hold in both pooled profiles and in all six cells, every one of them powered. The
validity condition holds too: `mm2` errs on 53% of the undecidable reads in both profiles.

| cell | assigned | wrong (B1) | undecidable: wrong (B2) | decidable: correct (B3) |
|---|---|---|---|---|
| A119b 401 | 3,938 | 0 | 0 / 1,236 | 3,640 / 3,649 (0.998) |
| A119b 402 | 2,317 | 0 | 0 / 645 | 2,052 / 2,083 (0.985) |
| A119b 403 | 5,237 | 0 | 0 / 910 | 5,151 / 5,170 (0.996) |
| public 401 | 3,851 | 0 | 0 / 1,432 | 3,531 / 3,540 (0.997) |
| public 402 | 2,027 | 0 | 0 / 849 | 1,724 / 1,762 (0.978) |
| public 403 | 5,261 | 0 | 0 / 1,018 | 5,198 / 5,233 (0.993) |

**Pooled per profile, every assigner (catalogued reads of the 6 worlds x 3 seeds):**

| assigner | profile | assigned | wrong | undecidable: wrong / correct | decidable correct | weak: correct / wrong |
|---|---|---|---|---|---|---|
| integrated O2 | A119b | 11,492 | 0 | 0 / 381 of 2,791 | 0.9946 of 10,902 | 268 / 0 of 2,032 |
| integrated O2 | public | 11,139 | 0 | 0 / 375 of 3,299 | 0.9922 of 10,535 | 311 / 0 of 2,013 |
| copy_assign O2 | A119b | 8,428 | 0 | 0 / 0 | 0.7729 | 2 / 0 |
| copy_assign O2 | public | 8,054 | 0 | 0 / 0 | 0.7642 | 3 / 0 |
| mm2 at MAPQ >= 1 | A119b | 12,955 | 21 | 21 / 0 | 1.0000 | 2,032 / 0 |
| mm2 at MAPQ >= 1 | public | 12,582 | 34 | 34 / 0 | 1.0000 | 2,013 / 0 |
| mm2, always | A119b | 15,729 | 1,483 | 1,479 / 1,312 | 1.0000 | 2,032 / 0 |
| mm2, always | public | 15,847 | 1,761 | 1,761 / 1,538 | 1.0000 | 2,013 / 0 |

**Reading:**
- **The integrated O2 made no wrong assignment in 22,631.** On the reads whose sequence cannot tell their copy apart,
  always-assigning errs on 53%.
- **It also placed 756 reads the edit-distance oracle calls undecidable, all correctly.** Its gap-compressed distance
  counts an indel run once and aligns clipped ends. That separates copies where NM plus clipping ties: homopolymer errors
  inflate NM equally on both copies.
- **copy_assign is as precise but much less powerful,** with 76-77% of decidable reads against 99.2-99.5%.
- **A plain MAPQ >= 1 filter is the strong baseline here.** Every decidable read has MAPQ >= 1, so it recovers all of
  them, and all weak ones. It pays for that with 55 wrong assignments, all on undecidable reads.
  - The integrated O2 trades recall on weak reads (13-15% against 100%) and 0.5-0.8% of decidable reads for zero errors.
  - The claim is precision and calibration at near-MAPQ power, not more power than MAPQ.
- **The `off_family` class is nearly empty here** (4 reads): these families' paralogs are all in the family.

**Not part of this test, but seen in the same runs (the round-8 end-to-end bars, O1 -> O2 -> O3):**
- **OUTSIDE fails B2 in all six cells:** integrated 138-171 errors, against 4-77 for the current method.
- **The cause is O1 missing the uncatalogued twin.**
  - O1 seeds the uncatalogued copy as a new locus in 6 of the 8 families, but not in C4 or HSFY, whose twins are
    near-identical.
  - Those copies' reads then have a single candidate in the catalog, and O2 places them on the catalogued twin. The
    current method abstains on them.
- **A119b REMOVED also fails B2 in two cells** (3-4 errors against 0-1).
- **These are O1's family-scope misses, which this test excludes by construction.** They are the open item: a read whose
  genome alignment ties with a locus O1 did not seed should be abstained on.
