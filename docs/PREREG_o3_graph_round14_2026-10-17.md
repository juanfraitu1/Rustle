# Round 14: the integrated pipeline v9 on a fifteenth set of full families, with a diploid donor (2026-10-17)

Pre-registration, written before any world of the fifteenth set was built: no read of these families has been simulated.
Code: branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.** Round 13 (register 1311) failed under A119b on one false O3 call. Its review traced the call to O3's mosaic distance,
which counted one long read's sequencing indels at positions only that read covers, on a present copy's haplotype-2 reads.
Round 12 (register 1310) failed on O2 errors: retained-intron reads that a paralog's path absorbed as one gap-compressed edit,
and two coincidences in an 80-copy array. v9 adds one rule against each mechanism that development could fix. The simulator's
donor is now diploid under both profiles, and development under that donor found the allele confound again, which v9's
haplotype check addresses.

**Question.** On a fresh set of full families, does v9 do at least as well as the current method in every world type, with both
calibrated v2 profiles: at least as many correct read outcomes and no more wrong ones?

## The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v9** = round 13's v5' + three rules:
- **O3 `--mosaic-shared`** (551c7685). The mosaic rule's nearest-copy distance is taken only over positions that >= 2 of the
  cluster's members cover. One read's sequencing errors where it alone covers are no evidence of the cluster's copy.
- **O3 `--haplotype-check`** (1c604e47 + cd93a327). A PSV call is dropped when the haplotype it is stays within delta of a copy:
  - the family's reads that O2 explains are aligned to the call's host window;
  - those carrying one of the call's unexplained alts (alts no copy's map carries) join it;
  - the extended cluster's nearest-copy distance, over positions >= 2 of them cover, must exceed delta.
  A present copy's haplotype-2 reads carry its donor alleles over the whole copy, whereas a missing copy's PSVs are carried
  by its own reads.
- **O2 `--o2-genome veto`** (fe6804cc; `o2graph.py run W integv 2`, replacing the fallback step `integl`). It is the fallback
  for reads the paths leave unexplained, and it also rechecks every read the paths assign, on the read's own genome
  alignments. A path call is overruled only when three things hold:
  - another seeded copy's locus, or a locus off them, beats the called copy's by the margin, within the gate; clipped ends and
    long insertions cost their length (`veto_distance`);
  - the read's primary alignment lies on that winner;
  - at MAPQ >= 30.
  A genome tie or a thin margin keeps the path call.
- Everything else as v5':
  - `o1seed.py --twins`;
  - `o2graph.py run W integ 2 --distance gapcomp --union-exact --error-floor`;
  - `detect.py --member-genome-check --member-attribution --error-floor --o1-hosts --branch-genome-novel --arms I`.
- **Measured and not adopted** (each turned a passing development run into a failure, set 10 A119b REMOVED):
  - `--allele-discount` v7f;
  - `--member-attribution-delta`;
  - `--o3-pool final` (it does not remove round 13's call, and it adds REMOVED and BRANCH errors).

**The simulator (both systems): `--sim-het-draw`** (551c7685). Each copy's het and hom site counts are drawn (Binomial over its
exonic bases) instead of rounded. Round 13's review found the public donor effectively haploid under rounding: 86 of 89 copies
had no het site. The calibrated rates are unchanged.

```
fastset.py run OUT FAMILIES SET SEED --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
    --o2-genome veto --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared
    --haplotype-check --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2 --sim-profile PROFILE
```

## Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout14`, e69d201e), all in no earlier
list:
- segmental-duplication pseudogene families: SMG1P, PKD1P, CHEK2P, RRN3P, CDC27P;
- few-copy gene families: AKR1C, HLA-DQ, PCDHB, ANKRD36, TPTE, NXF, SERPINA, MICAB;
- olfactory subfamilies: OR51, OR8, OR10, OR5;
- retrocopy families as spares: NPM1P, HMGB1P, GAPDHP, RPL21P, HNRNPA1P.

**Build rule:** round 12's, `families.py build families14_all.json --cap 1000 --link locus --min-cover 0.8`, every group (102 s).
The set is the first 8 families with >= 3 copies, in list order:

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| SMG1P | 5 | 14 | chr16 |
| PKD1P | 3 | 6 | chr16 |
| CHEK2P | 7 | 11 | chr2, chr10, chr15, chr22, chrY |
| CDC27P | 13 | 13 | chr2, chr13, chr14, chr15, chr20, chr21, chr22, chrY |
| PCDHB | 14 | 19 | chr5 |
| OR8 | 3 | 5 | chr11 |
| OR10 | 4 | 5 | chr11 |
| NPM1P | 25 | 25 | 13 chromosomes |

74 copies. RRN3P, AKR1C and OR5 gave 2 copies; HLA-DQ, ANKRD36, TPTE, NXF, SERPINA, MICAB and OR51 no family.
Files: `/mnt/linuxdisk/tmp/o3_graph/families14.json` (this set), `families14_all.json`, `families14.log`.

## Worlds

Six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v2.json` and
`profiles/public_testis_v2.json` (`--het 0`: the profile's calibrated heterozygosity): 12 worlds.

**The seed:** `precheck.py families14.json heldout14 201 211` (`precheck14.txt`; 96 s, 13.3 GB) gives the smallest seed in
201-210 with >= 3 detectable removed copies: **seed 201**.
- Detectable: OR8 0.1158, OR10 0.0834, NPM1P 0.0529, PCDHB 0.0324, CHEK2P 0.0161.
- Within delta: PKD1P 0.0019, SMG1P 0.0003, CDC27P 0.0.

**Commands** (one per profile; the first run builds the local reference):
```
F=/mnt/linuxdisk/tmp/o3_graph/families14.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout14
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout14 201 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis_v2.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout14 201 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis_v2.json
```

## Read outcomes, bars and reading

Round 12's, unchanged. The e2e scorer reads the veto's table (`--o2 o2.integv.m2.tsv`).
- **Read outcomes:** per-read detectability less the world's error floor. Round 10's rule and copy level are secondary.
- **Bars:** B1 (integrated correct >= current) and B2 (integrated errors <= current), per world, all 12 worlds.
- **PASS** (B1 and B2 hold in all 12 worlds). Claim: "the integrated O1 -> O2 -> O3 pipeline v9, which uses nothing from the
  current method, gives in every world type of a fresh set of full families at least as many correct read outcomes as the
  current method and no more wrong ones, with reads calibrated on two real libraries".
- **FAIL** otherwise. No partial claim.

## Known before running (development)

**Development data:** 16 spent runs (sets 8-10 and rounds 10-13, both profiles, two untrimmed-library runs; round 13 is
development since its review), plus 8 spent family sets re-simulated with the diploid donor (`--sim-het-draw`: rounds 9-11
and 13, both profiles). Each lever was measured alone against v5'. The adoption rule was fixed before any result was read: a
lever must remove round 13's false call or reduce errors, and must not turn a passing run into a failure
(`docs/O1O2O3_FAST_REALISTIC_2026-10-16.md` section 14, `fast/r14dev/`).

| | 16 spent runs passing | het-draw runs passing |
|---|---|---|
| v5' (round 13) | 13 | not run |
| v8 = v5' + `--mosaic-shared` + veto | 15 | round 13's set: A119b FAIL (FTH1P.6 haplotype), public PASS |
| **v9 = v8 + `--haplotype-check`** | **15** | **8 of 8** |

- **The one development failure is round 12 A119b.** Its 2 DUX4 margin-2 coincidences and 1 KRTAP2 read, DNA-identical at
  two copies, are placed by O2 in INTACT, CROSS and SILENT. The current method abstains on them, and no sequence rule
  separates them.
- **v9's O3 rules against v5' on the 16 runs** (O2 unchanged): 13 runs are identical and 3 change, all for the better:
  - round 13 A119b: its false call removed, REMOVED +2;
  - round 11 A119b: REMOVED +11;
  - set 10 A119b: REMOVED +45.
- **The veto:**
  - it fixes round 12's four CDRT15 retained-intron errors (round 12 public passes);
  - it cuts OUTSIDE errors in 8 runs;
  - it costs 1 and 2 CROSS reads in set 8.
- **The haplotype check took three versions:**
  - in genome coordinates it missed a 45-copy array (TAF11L.10's haplotype-2 reads, whose MAPQ-0 primaries lie on three
    units);
  - selecting carriers by every alt pulled in other copies' reads;
  - the adopted version selects by unexplained alts on the host window.
- **Not adopted:**
  - `--allele-discount` and `--member-attribution-delta` each turned set 10 A119b into a failure;
  - `--o3-pool final` kept round 13's call and added errors.
- **Cost:** the haplotype check aligns the family's explained reads to each host window that has a candidate call. On round
  13 A119b, O3's longest world went from 26 to 58 s.
- **Untested in development:** these families; segmental-duplication pseudogene families; a protocadherin cluster; the diploid
  donor on round 12's and sets 8's families.

## Reported beside the bars (not gating)

- fates per read class and error sources, tabulated, every world, both profiles;
- O1's drops (of how many silent copies), undecidable copies, new copies and twins;
- arm I's calls, true or false, every world;
- the veto's overrulings and the haplotype check's removals, right or wrong;
- round 10's scoring and copy-level scoring;
- leave-one-family-out margins, both profiles;
- the cost: `timing.json`, the integrated steps apart from the current method's, and O3's peak memory.

## After the run

- An independent review by a fresh subagent, given read-only copies.
- Register row 1312.
- Results appended below a line, without editing anything above it.
