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

---

# Results (2026-10-17)

Both runs completed as pre-registered, with no amendment. Outputs are in `fast/heldout14/{A119b,public}/`, frozen read-only
after scoring. Wall time was 305 s for A119b, which also built the local reference, and 232 s for public.

## Verdict: PASS under both profiles (12 of 12 worlds). Claimed.

Bars from `e2e.json` (primary scoring):

| world | A119b: correct (current / integrated) | errors | public: correct | errors |
|---|---|---|---|---|
| INTACT | 3,389 / 3,447 | 0 / 0 | 3,334 / 3,391 | 0 / 0 |
| CROSS | 3,449 / 3,513 | 0 / 0 | 3,387 / 3,451 | 0 / 0 |
| OUTSIDE | 3,308 / 3,588 | 232 / 39 | 3,288 / 3,574 | 220 / 32 |
| REMOVED | 3,502 / 3,584 | 0 / 0 | 3,472 / 3,558 | 4 / 0 |
| BRANCH | 3,393 / 3,691 | 248 / 125 | 3,339 / 3,610 | 270 / 143 |
| SILENT | 2,551 / 2,609 | 0 / 0 | 2,584 / 2,636 | 0 / 0 |

- **B1 holds in all 12 worlds,** by +52 to +298.
- **B2 holds in all 12.** In INTACT, CROSS, SILENT and REMOVED the integrated pipeline makes no error under either profile.
- **The donor is diploid under both profiles.** 51 of A119b's 74 expressed copies and 19 of public's have haplotype-2 reads
  (round 13 under rounding: public 3 of 89).

**Reported beside the bars:**
- **Secondary scorings:**
  - public passes under round 10's rule and at copy level;
  - **A119b fails both in REMOVED B2:** 1 vs 0 errors (round 10's rule) and 3 vs 0 (copy level). Every other cell holds.
- **Leave one family out: both profiles pass in all 16 rows.**
  - CHEK2P carries most of A119b's margins. Without it: REMOVED +0 (B1 holds with equality), INTACT +3, CROSS +6, SILENT +2.
  - Without it, public is at INTACT +6 and SILENT +2.
- **Where the integrated errors are:**
  - OUTSIDE: 39 and 32 reads of uncatalogued copies placed on a catalogued paralog;
  - BRANCH: 123-142 reads of branch copies that do not carry the novel exon, placed on the parent copy, plus 1-2 placed by
    the genome step.
- **Arm I:**
  - **REMOVED:** every detectable removed copy called under both profiles (CHEK2P.4, PCDHB.0, OR8.2, OR10.0, NPM1P.12). No
    call holds a present read in any world.
  - **BRANCH:** all 7 novel-exon copies called (SMG1P, PKD1P, CHEK2P, PCDHB, OR8, OR10, NPM1P). CDC27P has no branch copy.
- **O1:**
  - **SILENT:** 18 of 22 silent copies dropped; CHEK2P.2, PCDHB.1 and PKD1P.0 kept as undecidable.
  - **OUTSIDE:** new copies in 6 families (A119b) and 7 (public), one of them a twin (SMG1P).
- **Fates per read class** (correct / error / abstain / unknowable, current -> integrated):

| world | A119b catalogued | A119b other class | public catalogued | public other class |
|---|---|---|---|---|
| INTACT | 3389/0/514/0 -> 3447/0/456/0 | | 3334/0/505/0 -> 3391/0/448/0 | |
| CROSS | 3449/0/551/0 -> 3513/0/487/0 | | 3387/0/536/0 -> 3451/0/472/0 | |
| OUTSIDE | 3308/0/282/0 -> 3350/0/240/0 | uncatalogued 0/232/81/0 -> 238/39/36/0 | 3288/1/264/0 -> 3334/0/219/0 | uncatalogued 0/219/67/0 -> 240/32/14/0 |
| REMOVED | 3314/0/276/0 -> 3367/0/223/0 | absent detectable 180/0/48/0 -> 217/0/11/0; undetectable 8/0/33/44 -> 0/0/16/69 | 3294/0/259/0 -> 3342/0/211/0 | absent detectable 178/4/38/0 -> 216/0/4/0; undetectable 0/0/24/42 -> 0/0/0/66 |
| BRANCH | 3389/0/514/0 -> 3447/0/456/0 | absent detectable 4/248/139/0 -> 244/125/22/0 | 3334/0/505/0 -> 3391/0/448/0 | absent detectable 5/270/103/0 -> 219/143/16/0 |
| SILENT | 2551/0/464/0 -> 2609/0/406/0 | | 2584/0/467/0 -> 2636/0/415/0 | |

- **Cost:**
  - **The integrated steps, summed over six worlds:** 159 s (A119b) and 149 s (public).
    - O3: 98 and 100 s, longest world 44 s, peak 0.28 GB.
    - O2 18-21 s (+ the genome veto 18-23 s); O1 11-19 s.
  - **The current method:** 133 and 129 s, peaking at 2.5 GB per world.

## Independent review (2026-10-17, from 09:30 PDT) and qualifications

Review: `/mnt/linuxdisk/tmp/o3_graph/review14/REVIEW.md` (fresh subagent, read-only copies, its own scorer, end-to-end
reproductions).

**The PASS stands as a pre-registered result, but it is one fragile draw.**
- **What holds:**
  - the reviewer's scorer reproduces all three scorings, both profiles, cell for cell (0 differences in 268 cells);
  - the protocol is clean: frozen code, one launch with the pre-registered flags, no re-runs, no trace of these families in
    any earlier truth file;
  - both runs reproduce byte-identically end to end, also with the current method's outputs hidden.
- **Post-hoc seed probe** (frozen code, same families, seeds 202-205; it cannot change the verdict): public passes 4 of 4.
  **A119b passes 2 of 4:**
  - seed 202 fails REMOVED 6 vs 0: present NPM1P.14 reads with drawn homozygous alleles, called on the window of a new copy
    that O1 seeds at the NPM1 parent gene; the genome check passes them by +0.0005;
  - seed 205 fails INTACT and CROSS 4 vs 0: CDC27P.8's haplotype-2 reads called on a paralog's window; the haplotype check
    misses by 0.00024 (0.00982 vs delta).
  Both are round 13's mechanism, a present copy's donor-allele reads called on another locus's window. **The allele confound is
  not solved.**
- **What decided this round's cells:**
  - public CROSS rests on `--haplotype-check` alone: without it B2 fails 4 vs 0, so round 13's v5' would have failed. The
    removed call is 4 hap-1 reads of CHEK2P.5, a 246-bp copy, carrying a homozygous allele that `--sim-het-draw` drew; the
    genome check let it through by +0.0011;
  - the veto overruled 0 reads in 12 worlds (`o2.integv` = `o2.integl`), and `--mosaic-shared` changed no bar: round 14 is no
    evidence for either;
  - A119b REMOVED B2 (0 vs 0) is decided by one placed PCDHB.0 read 0.00008 under the detectability threshold. At floor 0.0031
    instead of 0.0032 the cell fails. Round 10's rule and copy level fail it (1 and 3 vs 0);
  - most of the B1 margin in INTACT, CROSS, SILENT and REMOVED is one copy's reads (CHEK2P.5, which copy_assign mostly calls
    ambiguous). Without CHEK2P the margins are +0 to +14.
- **Corrections to the results:**
  - BRANCH's integrated errors are not DNA-identical: 105 of 125 (A119b) and 127 of 143 (public) carry the branch copy's SNPs.
    OUTSIDE's are DNA-identical.
  - The pre-registered report of the veto's overrulings and the haplotype check's removals is missing. It is given above: 0
    overrulings; removals = the public CROSS call above, plus PCDHB.B PSV calls whose reads are also in PCDHB's branch call.
- **Undisclosed in the pre-registration:** the family set and the seed were fixed (01:33-01:40) 7.5 h before the method was
  frozen, during v8/v9 development. The reviewer found no use of them in any development run.

**The claim is literally true for this draw under the primary scoring and must be cited with these qualifiers:**
- one simulated draw; at 4 further seeds, A119b fails 2;
- the two secondary scorings fail A119b REMOVED;
- public CROSS depends on the haplotype check;
- the veto and the shared-position mosaic distance were not exercised;
- it is no evidence that the allele confound is solved.
