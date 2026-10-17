# Round 11: the integrated pipeline v4 under calibrated v2 reads, end to end, on a twelfth set of full families (2026-10-16)

Pre-registration, written before any world of the twelfth set was built: no read of these families has been simulated. Code:
branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.** The user's goal is unchanged: integrate O1, O2 and O3 so full families are fast and cheap to test, make the
simulations more realistic, and determine what PSVs/SUNs and the CIGAR, NM and `de` tags contribute.

Round 10 (register 1308) failed. Its independent review corrected the first failure's mechanism to a host-window artefact. v4
removes each round-10 failure mechanism found on development data (`docs/O1O2O3_FAST_REALISTIC_2026-10-16.md` section 12).

**Question.** On a fresh set of full (uncapped) families, with reads from both calibrated v2 profiles, does the integrated
pipeline v4 do at least as well as the current method in every world type? That means at least as many correct read outcomes
and no more wrong ones.

## The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v4** = round 10's v3 plus three changes, all chosen on development data (sets 8-10 and round 10's spent worlds):
- **O1:** `o1seed.py W --twins`.
  - **The keep rule:** >= 5 primaries, or >= 1 owned MAPQ >= 1 primary.
  - **A new copy within delta of a kept curated copy is its twin** (genomic spans, minimap2 asm20 `de`).
  - **A curated copy the rule would drop is kept when it lies within delta of a kept one.** The reads cannot say which of the
    two is expressed.
- **O2:** `o2graph.py run W integ 2 --distance gapcomp --union-exact --error-floor`, then `run W integl 2 --error-floor`.
  - This is round 10's union paths and genome fallback.
  - **A twin takes a read from its curated copy only by the margin.**
  - **The `unexplained` gate is delta plus the world's error floor.**
- **O3:** `detect.py W --member-genome-check --member-attribution --error-floor --arms I`.
  - **Attribution keeps no read that no catalogued copy competes for.**
  - **The genome check subtracts the error floor.**
  - **The read hosts still come from copy_assign's run** (`cur/A.assignments.tsv`), as in rounds 7-10. O1 hosting
    (`--o1-hosts`) was measured and not adopted: it loses REMOVED and BRANCH recall in arrays.
- **The error floor** is the world's median `de` over MAPQ-60 primary alignments, i.e. its reads' sequencing error. Delta
  (0.00958) is the p99 exonic divergence between KB3781's haplotype assemblies (register 1199), an error-free divergence.
  Every read-level comparison with delta therefore subtracts the floor first.

Both systems run per world from one command:
```
fastset.py run OUT FAMILIES SET SEED --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
    --o2-genome fallback --per-read-detectability --twins --error-floor --het 0 --jobs 2 --threads 2 --sim-profile PROFILE
```
`--o2-procs 2` and `--o3-procs 2` (defaults) leave outputs byte-identical.

## Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout11`, committed in c034cc5d):
- FOXD4L, OR2T, KRTAP10, H2BC, FAM66, CXADRP, GOLGA2P: round-10 candidates past its first 8, never built into a world;
- H3C, H2AC, KRTAP2, KRTAP1, PRR23, ETD, SULT1C, UGT2A, CDRT15, NAIP, HLA-DRB, TCEAL, ARMCX, OR4C, OR4K, OR11H, SEC14L, DUX4:
  none was in any earlier list.

**Build rule:** round 10's (locus link >= 0.8 at identity >= 0.90, span <= 60 kb, no copy cap):
`families.py build families11_all.json --cap 1000 --link locus --min-cover 0.8`, every group. The set is the first 8 families
with >= 3 copies, in list order (GOLGA2P links 2 copies):

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| FOXD4L | 6 | 6 | chr2, chr9 |
| OR2T | 3 | 3 | chr1 |
| KRTAP10 | 9 | 10 | chr21 |
| H2BC | 6 | 6 | chr6 |
| FAM66 | 3 | 10 | chr8, chr12 |
| CXADRP | 3 | 3 | chr15, chr18, chr21 |
| H3C | 3 | 3 | chr1 |
| H2AC | 4 | 4 | chr6 |

37 copies. Files: `/mnt/linuxdisk/tmp/o3_graph/families11.json` (this set only), `families11_all.json`, `families11.log`.

## Worlds

Six world types (INTACT, CROSS with liftover projection, OUTSIDE, REMOVED, BRANCH, SILENT), each under round 10's two read
profiles, `profiles/A119b_testis_v2.json` and `profiles/public_testis_v2.json` (`--het 0`). That makes 12 worlds.
- **The seed:** `precheck.py families11.json heldout11 141 151` (`precheck11.txt`; 93.9 s, 13.2 GB) gives the smallest seed
  in 141-150 with >= 3 detectable removed copies: **seed 141**.
  - Detectable: H2AC 0.1027, H2BC 0.0943, KRTAP10 0.0544, FAM66 0.0448, FOXD4L 0.0153, CXADRP 0.0113.
  - Within delta: OR2T 0.0061, H3C 0.0.

**Commands** (one per profile; the first run builds the local reference):
```
F=/mnt/linuxdisk/tmp/o3_graph/families11.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout11
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome fallback --per-read-detectability --twins --error-floor --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout11 141 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis_v2.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout11 141 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis_v2.json
```
`fastset.py` writes `e2e.json` and `bars.txt` per profile. The scoring is per-read detectability with the error floor, on
`o2.integl.m2.tsv`. If a run's current step exceeds the light limit, the run is repeated with `--heavy-current`, whose output
does not depend on it.

## Read outcomes

**These are round 10's, with one change:** a REMOVED copy's read is detectable when its own primary `de`, minus its world's
error floor, exceeds delta (`summary.py e2e --error-floor`).
- **Why:**
  - Delta is an error-free haplotype divergence.
  - A read's `de` adds its own sequencing errors.
  - Without the correction, a read whose copy lies within delta of a present copy counts as detectable merely through its
    errors.
- **The change applies to both systems.** It moves REMOVED reads with `de` in (delta, delta + floor] from detectable to
  undetectable: placing them becomes `unknowable` instead of an error, and an O3 call on them stays correct.
- **The rule was decided on development data after round 10.**
- **Secondary, not gating:** round 10's rule (`e2e_r10rule.json`, same command without `--error-floor`) and copy-level
  scoring are both reported.

Unchanged:
- a BRANCH read is detectable by construction;
- an O3 call on a read takes precedence over its O2 placement;
- a read called in several families is judged by its call in its own family.

## Bars, per world, all 12 worlds

- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

Abstain and unknowable count for neither.

## Reading

- **PASS** (B1 and B2 hold in all 12 worlds under the primary scoring). Claim: "with reads calibrated on two real libraries,
  the integrated O1 -> O2 -> O3 pipeline v4 gives, in every world type of a fresh set of full families, at least as many
  correct read outcomes as the current method and no more wrong ones, with per-read detectability corrected for sequencing
  error".
- **FAIL** otherwise. No partial claim; the failing worlds and margins are reported.
- **If the primary passes and round 10's rule fails,** the claim is stated with that qualifier.

## Known before running (development, v2 profiles)

Nine runs: development sets 8, 9 and 10 under both profiles, plus round 10's two profiles and its untrimmed-library run
(`fast/r11dev6/`):
- **Primary scoring:** v4 PASSES all 9.
- **Round 10's rule:** v4 passes 5 of 9. The 4 failures are all REMOVED B2: set 8 public 9 vs 7, set 9 A119b 1 vs 0, set 10
  A119b 33 vs 10, set 10 public 13 vs 4. Placements of reads with `de` between delta and delta + floor.
- **v4's margins (correct, integrated minus current):**
  - INTACT: +92, +78, +54, +101, +37, +48, +8, +8;
  - CROSS: +90, +78, +191, +234, +36, +49, +17, +15;
  - SILENT: +12, +15, +69, +114, +32, +42, +11, +7;
  - REMOVED: +90, +76, +67, +87, +20, +12, +1, +55.
- **Errors in INTACT, CROSS and SILENT:** 0-1 per world, against 0-9 for the current method.
- **The thinnest margins:** set 10 A119b INTACT +8 and REMOVED +1, against 1 error vs 4.
- **Untested in development:**
  - a family set without arrays or near-identical copies;
  - histone clusters (H2BC, H3C, H2AC), whose copies are short and single-exon.

## Reported beside the bars (not gating)

- the fates per read class and the error sources, tabulated;
- O1's drops (of how many silent copies), its kept-undecidable copies, its new copies and twins;
- arm I's calls, true or false, in every world;
- round 10's scoring and copy-level scoring;
- **the run's cost:** `timing.json`, every step's wall time and peak RSS; the integrated pipeline's own steps (O1, O2,
  O2_genome, O3) apart from the current method's;
- **an untrimmed-library run:** `public_testis_v3.json` + `--trim-adapters profiles/public_testis_adapters.json`, same seed
  and options.

## After the run

- An independent review by a fresh subagent, given read-only copies.
- Register row 1309.
- Results appended below a line, without editing anything above it.
