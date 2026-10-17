# Round 10: the integrated pipeline v3 under calibrated v2 reads, end to end, on an eleventh set of full families (2026-10-16)

Pre-registration, written before any world of the eleventh set was built: no read of these families has been simulated. Code:
branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.** The user's goal (2026-10-16): "integrate O1, O2 and O3 more so it's faster and cheaper to test full families; also
improve the simulations to account for more realistic ones and determine if PSVs can also include SUNs and probably some
influence of the cigar string, nm, de tags". Round 9 failed under both profiles. The development work since then removes
the round-9 failure mechanisms; see "Known before running" below and `docs/O1O2O3_FAST_REALISTIC_2026-10-16.md` section 11.
The round-9 failures were:
- O2 left realistic reads `unexplained` on catalogued copies (1-3.5% fewer correct reads);
- O1 dropped expressed copies;
- O3 made false calls on present copies' reads.

**Question.** On a fresh set of full (uncapped) families, with reads from both calibrated v2 profiles, does the integrated
pipeline v3 do at least as well as the current method in every world type? That means at least as many correct read outcomes
and no more wrong ones.

## The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v3** (every change since round 9 was chosen on development sets 8-10, all spent):
- **O1:** `o1seed.py`. A curated copy is kept with >= 5 primaries OR >= 1 owned primary at MAPQ >= 1 (e40367ea).
- **O2:** `o2graph.py run W integ 2 --distance gapcomp --union-exact`, then `o2graph.py run W integl 2`.
  - `--union-exact` (36a430f8): a paralog isoform projects onto a copy even when the spliced aligner drops a terminal exon
    < 30 bp. A projected structure stays apart from the copy's own unless its junctions match exactly.
  - `integl` (eb4f0c95), the genome fallback for the reads the paths leave `unexplained`:
    - every primary and secondary genome alignment is scored as mismatches + indel events;
    - an insertion or clip >= 30 bp refuses (it is O3's evidence);
    - the read's best alignment off the family's seeded copies keeps it `unexplained`;
    - the read is placed only at margin 2 within delta.
  - O3 still reads `o2.integ.m2.tsv`, so the fallback changes no O3 evidence.
  - The set is scored on `o2.integl.m2.tsv`.
- **O3:** `detect.py W --member-genome-check --member-attribution --arms I`.
  - `--member-attribution` (12f8f2af): an absent-copy PSV call keeps only the reads its consensus explains better than every
    single catalogued copy, by O2's margin of 2.
  - The mosaic rule's nearest copy compares off-site positions with what the cluster's reads show, not with the host's base
    (f080a30f).
  - Speedups with identical calls: lazy isoform maps, faster PSV clustering (9829bc06).

Both systems run per world from one command. The O2 and O3 speedups (`--o2-procs 2`, `--o3-procs 2`) leave outputs
byte-identical.
```
fastset.py run OUT FAMILIES SET SEED --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
    --o2-genome fallback --per-read-detectability --het 0 --jobs 2 --threads 2 --sim-profile PROFILE
```

## Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout10`, committed in 70ebf2b1):
- CYP3A, FAM47, FAM9, DAZ: round-9 candidates past its first 8, never built into a world;
- HNRNPCL, PRR20, ZNF705, TAF11L, F8A, FAM236, CENPVL, FOXD4L, OR2T, KRTAP10, H2BC, FAM66, CXADRP, GOLGA2P: none was in any
  earlier list.

**Build rule:** the locus link of the rebuilt simulation sets.
- A gene's canonical transcript, spliced onto the other gene's locus, covers >= 0.8 of it in one direction, at identity >= 0.90.
- Span <= 60 kb.
- **No copy cap** (`--cap 1000`).
- `families.py build families10_all.json --cap 1000 --link locus --min-cover 0.8`, every group.

The set is the first 8 families with >= 3 copies, in list order. CYP3A links 2 copies; FAM47 and FAM9 form no family.

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| DAZ | 3 | 7 | chrY |
| HNRNPCL | 3 | 3 | chr1 |
| PRR20 | 3 | 4 | chr13 |
| ZNF705 | 6 | 10 | chr8, chr11, chr12 |
| TAF11L | 45 | 45 | chr5 |
| F8A | 3 | 3 | chrX |
| FAM236 | 4 | 7 | chrX |
| CENPVL | 3 | 3 | chrX |

70 copies. Files: `/mnt/linuxdisk/tmp/o3_graph/families10.json` (this set only), `families10_all.json`, `families10.log`.

## Worlds

Six world types (INTACT, CROSS with liftover projection, OUTSIDE, REMOVED, BRANCH, SILENT), each under two read profiles:
`bench/o3_graph/profiles/A119b_testis_v2.json` and `bench/o3_graph/profiles/public_testis_v2.json`. That makes 12 worlds.
- **The v2 profiles (round-9 review):**
  - errors and donor alleles are taken from every MAPQ-60 read;
  - the expression spread is measured over read copies, less Poisson noise;
  - truncation is redrawn below 200 bases;
  - plus round 9's calibrated truncation, homopolymer-biased indels, het and hom donor alleles, and retained introns.
- With a profile, `--het 0` means the donor's alleles come from it.
- **The seed:** `precheck.py families10.json heldout10 121 131` (`precheck10.txt`) gives the smallest seed in 121-130 with
  >= 3 detectable removed copies: **seed 122** (42.6 s, 14.4 GB peak):
  - detectable: HNRNPCL 0.0453, TAF11L 0.0444, ZNF705 0.0244;
  - within delta: DAZ 0.0, PRR20 0.0011, F8A 0.0, FAM236 0.0022, CENPVL 0.0006.

**Commands** (one per profile; the first run builds the local reference):
```
F=/mnt/linuxdisk/tmp/o3_graph/families10.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout10
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome fallback --per-read-detectability --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout10 122 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis_v2.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout10 122 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis_v2.json
```
`fastset.py` writes `e2e.json` (per-read detectability, on `o2.integl.m2.tsv`) and `bars.txt` (`summary.py bars`) per profile.
Full families with a 45-copy array may need copy_assign under the heavy lock. If a run's current step exceeds the light
limit, the run is repeated with `--heavy-current`, whose output does not depend on it.

## Read outcomes

These are round 9's (`summary.py` `read_fate`, `cmd_e2e`), unchanged:
- per-read detectability: a REMOVED copy's read is detectable by its own primary `de` on the masked genome (> delta);
- a BRANCH read is detectable by construction;
- an O3 call on a read takes precedence over its O2 placement;
- a read called in several families is judged by its call in its own family.

## Bars, per world, all 12 worlds

- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

Abstain and unknowable count for neither (`summary.py bars`).

## Reading

- **PASS** (B1 and B2 hold in all 12 worlds). Claim: "with reads calibrated on two real libraries, the integrated O1 -> O2 ->
  O3 pipeline v3 gives, in every world type of a fresh set of full families, at least as many correct read outcomes as the
  current method and no more wrong ones".
- **FAIL** otherwise. No partial claim; the failing worlds and margins are reported.

## Known before running (development sets 8, 9 and 10, v2 profiles, per-read scoring)

Every development run below was simulated fresh with the v2 profiles and the new O1 rule:
- `fast/fbdev/` (O3 runs: `calls.tsv.prev`, `calls.noattr.tsv`, `calls.tsv`);
- `fast/fbdev_exact/` (`--union-exact`).

| configuration | set 8 A119b | set 8 public | set 9 A119b | set 9 public | set 10 A119b | set 10 public |
|---|---|---|---|---|---|---|
| no fallback (round 9's O2), old O3 | FAIL | FAIL | FAIL | PASS | FAIL | FAIL |
| fallback, old O3 | FAIL | FAIL | FAIL | PASS | FAIL | PASS |
| fallback, O3 consensus fix | FAIL | FAIL | FAIL | PASS | FAIL | PASS |
| fallback, + member attribution | PASS | FAIL | PASS | PASS | FAIL | PASS |
| **v3: fallback, attribution, `--union-exact`** | **PASS** | **PASS** | **PASS** | **PASS** | **FAIL** | **PASS** |

- **v3 fails only set 10 A119b, by 2 correct reads in INTACT** (1,136 vs 1,138; errors 0 vs 2). The lost reads are short
  reads whose sequencing edits exceed delta's share of their length: 2-5 edits on 184-484 bp, MAPQ 60, on their own copy.
  The current method assigns them.
- **The rejected alternatives,** measured on the same runs:
  - `combined` mode, the genome on every read with the elsewhere guard: errors drop (OUTSIDE 76 -> 0), but correct reads
    collapse in arrays (set 9 INTACT 935 -> 534).
  - `seeded` mode, the genome on explained reads without the guard: fixes CROSS but costs INTACT 13-47 reads.
  - `resolve` mode, the genome on O2's tied reads: adds 0-1 read per world.
- **The development summary of v3's margins** (correct, integrated minus current):
  - INTACT: +24, +97, +28, +45, -2, +6;
  - CROSS: +157, +230, +23, +46, +8, +13;
  - SILENT: +57, +176, +34, +41, +13, +11.
- **Errors:**
  - in INTACT, CROSS and SILENT: 0-1 per world, against 0-9 for the current method;
  - in REMOVED: 0-2 (current 0-12);
  - in OUTSIDE and BRANCH: far fewer than the current method.
- **A per-world bar is decided by single reads** wherever the current method makes no error, and by a handful of reads in
  INTACT, where both systems see the full catalog.

## Reported beside the bars (not gating)

- the fates per read class and the error sources;
- O1's drops and new copies;
- arm I's calls, true or false;
- copy-level scoring (`summary.py e2e` without `--per-read-detectability`);
- **the run's cost:** `timing.json`, every step's wall time and peak RSS;
- **an untrimmed-library run:** the public profile with adapters (`public_testis_v3.json`) on the same seed, trimmed by
  `fastset.py --trim-adapters profiles/public_testis_adapters.json`, same options, scored the same way.

## After the run

- An independent review by a fresh subagent, given read-only copies.
- Register row 1308.
- Results appended below a line, without editing anything above it.

---

# Results (2026-10-16; run with the code of 58bca61d, started 19:30)

## Verdict: FAIL under both profiles. Not claimed.

Bars from `fastset.py`'s `e2e.json` (per-read detectability, on `o2.integl.m2.tsv`):

| world | A119b: correct (current / integrated) | errors | B1, B2 | public: correct | errors | B1, B2 |
|---|---|---|---|---|---|---|
| INTACT | 645 / 736 | 0 / **8** | yes, **no** | 576 / 654 | 0 / 0 | yes, yes |
| CROSS | 689 / 778 | 0 / **8** | yes, **no** | 616 / 694 | 0 / 0 | yes, yes |
| OUTSIDE | 828 / 830 | 142 / 61 | yes, yes | 802 / **751** | 115 / 57 | **no**, yes |
| REMOVED | 879 / 959 | 8 / 8 | yes, yes | 843 / 919 | 6 / 0 | yes, yes |
| BRANCH | 682 / 1,043 | 154 / 65 | yes, yes | 576 / 946 | 139 / 48 | yes, yes |
| SILENT | 464 / 475 | 0 / **8** | yes, **no** | 424 / 439 | 0 / 0 | yes, yes |

B1 holds in 11 of 12 worlds, by +2 to +370 correct reads. The integrated pipeline makes fewer errors in OUTSIDE, BRANCH
and REMOVED. Copy-level scoring (`e2e_copylevel.json`) gives the same verdict; there it also fails REMOVED's B2 (A119b 9 vs
0, public 2 vs 1).

**The two failures, each one mechanism, found post hoc:**
1. **A119b, B2 in INTACT, CROSS and SILENT: one false O3 call, in TAF11L, the 45-copy array.**
   - Arm I calls 8 reads of TAF11L.33, a present copy, an absent copy, on host 44 (est_div 0.012-0.016, members' median
     genome `de` 0.0103).
   - The same call recurs in BRANCH, where it costs nothing.
   - **Cause:** the family's unexplained pool (281 reads from about 40 copies) clusters on one host into a few large mixed
     clusters. Their consensus carries several copies' private alleles, so it is `unexplained` beyond delta.
   - Member attribution then keeps the TAF11L.33 reads whose own alleles match that consensus at >= 2 positions where
     TAF11L.33's catalogued sequence differs.
   - Development's arrays (set 9's USP17L and TSPY) did not produce such a call.
2. **public, B1 in OUTSIDE (-51): O1's new copy within delta of a catalogued one.**
   - O1 seeds the uncatalogued FAM236.1 correctly, and the integrated pipeline gets 40 reads of uncatalogued copies right.
   - FAM236.1 differs from FAM236.0 by about 0.2% (the pre-check's 0.0022), so 166 catalogued FAM236 reads become `thin` or
     `tied` between the two.
   - The current method, which does not know FAM236.1, assigns those reads, and makes 114 errors on uncatalogued reads.

**Reported beside the bars:**
- **Arm I:**
  - REMOVED: all 3 detectable removed copies called under both profiles (HNRNPCL.1, TAF11L.41, ZNF705.0);
  - BRANCH: 5 of 5 novel-exon copies under both profiles;
  - false calls: the TAF11L.33 call above (A119b: INTACT, CROSS, BRANCH, SILENT); none under the public profile.
- **O1:** 5 silent copies dropped per SILENT world and no expressed copy; 3-4 new copies in OUTSIDE.
- **Cost:**
  - A119b 327 s for six worlds, including the one-time local reference build (44 s, 14.4 GB);
  - public 273 s, peak 1.9 GB (copy_assign);
  - per set: O3 170-185 s, current 131-144 s, O1 53-61 s, O2 33-40 s, the genome fallback 2-3 s.
- **The untrimmed-library run:** below.

**The untrimmed-library run** (secondary; the public profile with its Tn5/SMART adapters, `public_testis_v3.json`, trimmed by
`--trim-adapters`, same seed and options): FAIL. Correct reads match the adapter-free public run (INTACT 656 vs 654, BRANCH
947 vs 946), so trimming restores the reads; two known mechanisms add errors:
- **INTACT, CROSS, SILENT, 2 errors each:** O1 dropped PRR20.1, an expressed copy all of whose reads tie, and O2 placed two of
  its reads on PRR20.0.
- **REMOVED, 17 errors:** a second mixed-cluster call in TAF11L (host 43, 25 members from present copies), round 10's
  failure 1.
- **OUTSIDE B1:** 753 vs 799, failure 2.
- Its O3 step ran while `detect.py` gained the opt-in `--member-attribution-delta`. The run did not pass it, so arm I behaved as
  at 58bca61d.

## Reading

**The pre-registered claim fails under both profiles. No partial claim is made.** Set beside round 9, on held-out families:
- **The catalogued-read loss under realistic reads is gone.** B1 holds in INTACT, CROSS and SILENT under both profiles, by +11 to
  +91 reads; round 9 lost 1-3.5% there.
- **The remaining errors are O3 calls in a 45-copy near-identical array** (A119b), and one O1 new copy within delta of a
  catalogued copy that makes that copy's reads ambiguous (public OUTSIDE).
- Neither mechanism occurred on the development sets.
