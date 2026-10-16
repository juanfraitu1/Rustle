# Round 9: the integrated pipeline v2 under realistic reads, end to end, on a tenth set of full families (2026-10-16)

Pre-registration, written before any world of the tenth set was built: no read of these families has been simulated. Code:
branch `o3-fast` (`bench/o3_graph/`), frozen at the commit that adds this file.

**Why (the user's goal, 2026-10-16):** "integrate O1, O2 and O3 more so it's faster and cheaper to test full families; also
improve the simulations to account for more realistic ones and determine if PSVs can also include SUNs and probably some
influence of the cigar string, nm, de tags".

The development work behind this test is `docs/O1O2O3_FAST_REALISTIC_2026-10-16.md` (sets 8 and 9):
- **Calibrated realistic reads remove the integrated pipeline's advantage on catalogued copies.**
- **CIGAR-aware evidence restores it:**
  - O2 counts indel EVENTS, not indel bases, and aligns clipped ends against the path;
  - O3 drops, read by read, the members the genome explains.

**Question.** On a fresh set of full (uncapped) families, with reads from two calibrated real libraries, does the integrated
pipeline v2 do at least as well as the current method in every world type? That means at least as many correct read outcomes
and no more wrong ones.

## The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v2:**
- **O1:** `o1seed.py` (round 8's, unchanged);
- **O2:** `o2graph.py run W integ 2 --distance gapcomp`;
- **O3:** `detect.py W --member-genome-check`, arm I.

Both run per world from `fastset.py run ... --o2-distance gapcomp --member-genome-check`, with the set's local reference
(`localref.py`, validated in the development doc).

## Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout9`):
- PPIAL4, FCGR1, PGA, H4C, LCE2, LCE3, KRTAP4, AGAP;
- ELOA3, FAM25, CYP3A, FAM47, FAM9, KRTAP9, DAZ.

PPIAL4 and FCGR1 were round-8 candidates past its first 8, never built into a world; none of the others was in any earlier
list.

The build rule is rounds 1-8's (identity >= 0.90 over >= half the shorter canonical transcript, span <= 60 kb), with
**no copy cap** (`--cap 1000`: full families). The set is the first 8 with >= 3 copies (H4C and LCE3 link 2 copies each):

| family | copies | isoforms | chromosome |
|---|---|---|---|
| PPIAL4 | 7 | 7 | chr1 |
| FCGR1 | 3 | 9 | chr1 |
| PGA | 3 | 12 | chr11 |
| LCE2 | 3 | 3 | chr1 |
| KRTAP4 | 5 | 5 | chr17 |
| AGAP | 7 | 15 | chr10 |
| ELOA3 | 12 | 12 | chr18 |
| FAM25 | 3 | 4 | chr10 |

Files: `/mnt/linuxdisk/tmp/o3_graph/families9.json` (round 8's sets plus this one), `families9_all.json`, `families9.log`.

## Worlds

Six world types (INTACT, CROSS with liftover projection, OUTSIDE, REMOVED, BRANCH, SILENT), each under two read profiles:
`bench/o3_graph/profiles/A119b_testis.json` and `bench/o3_graph/profiles/public_testis.json`. That makes 12 worlds.
- **Profiles:** calibrated truncation, homopolymer-biased indels, the donor's het and hom alleles, retained introns, and
  per-copy expression spread. With a profile, `--het 0` means the donor's alleles come from it.
- **Seed 101.** `precheck.py families9.json heldout9 101 111` (`precheck9.txt`) gives the smallest seed in 101-110 with >= 3
  detectable removed copies:
  - detectable: LCE2 0.068, KRTAP4 0.039, FAM25 0.017, AGAP 0.014;
  - within delta: PPIAL4, FCGR1, PGA, ELOA3.

**Commands** (one per profile; the local reference is built by the first run):
```
F=/mnt/linuxdisk/tmp/o3_graph/families9.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout9
V="--o2-distance gapcomp --member-genome-check --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout9 101 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout9 101 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis.json
python3 bench/o3_graph/summary.py e2e $O/<profile>/e2e_pr.json <its six worlds> --per-read-detectability
python3 bench/o3_graph/summary.py bars $O/<profile>/e2e_pr.json
```

## Read outcomes

These are round 8's (`summary.py` `read_fate`, `cmd_e2e`), with one change, made on development data before this
registration: **a REMOVED copy's read is detectable by its own primary `de` on the masked genome (> delta)**
(`--per-read-detectability`).
- **Why:** on set 9 the one remaining error under both profiles was a read of a detectable copy (median `de` 0.018) that is
  near-identical to a present copy over its own truncated span. No method can tell that read apart.
- **The same rule applies to both systems.**
- **Secondary report:** round 8's copy-level classification (`e2e.json`) is reported too, not gating.

A read called in several families is judged by its call in its own family, as before.

## Bars, per world, all 12 worlds

- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

Abstain and unknowable count for neither (`summary.py bars`).

## Reading

- **PASS** (B1 and B2 hold in all 12 worlds). Claim: "with reads calibrated on two real libraries, the integrated O1 -> O2 ->
  O3 pipeline v2 gives, in every world type of a fresh set of full families, at least as many correct read outcomes as the
  current method and no more wrong ones".
- **FAIL** otherwise. No partial claim; the failing worlds and margins are reported.

## Known before running (development sets 8 and 9)

- **Integrated v1** (round 8's code) fails under both profiles on both sets: see the development doc, section 3.
- **Integrated v2, set 9:**
  - B1 holds in all 12 worlds;
  - copy-level scoring: B2 fails only in REMOVED, 1 vs 0, under both profiles;
  - per-read detectability: PASS under both profiles.
- **Integrated v2, set 8:**
  - public profile: PASS under both scorers;
  - A119b profile: FAIL under both. B1 fails in INTACT, 2,215 vs 2,234 correct. Under copy-level scoring B2 also fails in
    REMOVED, 5 vs 3; per-read scoring makes it 3 vs 6.
  - The INTACT shortfall is O2's `unexplained` reads: 37 in SPRR2, 12 in AQP7. The current method assigns them.
  - A MAPQ-60 genome fallback for them was measured and rejected: 17 wrong of 329 outside BRANCH, 159 wrong of 220 in BRANCH.
- **Development summary:** under per-read scoring v2 passes 3 of the 4 realistic set x profile runs (all but set 8, A119b);
  v1 passed none.
- **Integrated v2 under idealized reads** (copy-level scoring): PASS on sets 8 and 9. On set 9 it fixes round 8's single
  CROSS error, at a cost of 4-6 correct reads per world. On set 8 it gets fewer reads right than v1 (INTACT 2,407 vs 2,476)
  with fewer errors (CROSS 3 vs 11, REMOVED 0 vs 13).
- **What was chosen on these sets:** the options (`gapcomp`, `--member-genome-check`) and the per-read rule. This tenth set
  is their first test.
- **A per-world error bar is decided by single reads** wherever the current method makes no error.

## Reported beside the bars (not gating)

- the fates per read class, and the error sources;
- O1's drops and new copies;
- arm I's calls, true or false;
- **the run's cost:** `timing.json`, every step's wall time and peak RSS.

## After the run

- An independent review by a fresh subagent.
- Register row 1307.
- Results appended below a line, without editing anything above it.
