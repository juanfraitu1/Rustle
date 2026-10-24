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

## Amendment 1 (2026-10-16 15:28, before any result)

**The first run crashed in O3 under both profiles.**
- `vg msga` (v1.73, a deprecated subcommand) fails an assertion (`find_mems_deep: lcp_maxima.size() == mems.size()`) on
  KRTAP4's five loci in FASTA order, both by default and with `-a`.
- Any two of the loci build, and every other order tried builds.

**The fix (`detect.py`, `msga_orders`, tested):**
- retry `vg msga -a` with the loci reversed, then each rotation;
- if every order fails, the family gets no graph and no branch calls; its PSV clustering still runs.
- The reversed order builds KRTAP4's graph.

**Nothing else changes.**
- On development worlds, where msga never failed twice, `detect.py`'s output is byte-identical.
- No outcome of round 9 was computed before this amendment: both runs stopped at O3, and no e2e or bars were written.
- The two crashed run directories were deleted; the run is repeated from the start with the same commands.

## After the run

- An independent review by a fresh subagent.
- Register row 1307.
- Results appended below a line, without editing anything above it.

---

# Results (2026-10-16; run after amendment 1 with the code of 8235f733)

## Verdict: FAIL under both profiles. Not claimed.

Bars from `summary.py bars e2e_pr.json` (per-read detectability). Copy-level scoring gives the same verdict, cell for cell
except REMOVED's error counts (A119b 10 vs 2, public 5 vs 0).

| world | A119b profile: correct (current / integrated) | errors | B1, B2 | public profile: correct | errors | B1, B2 |
|---|---|---|---|---|---|---|
| INTACT | 1,149 / 1,109 | 2 / 0 | **no**, yes | 994 / 976 | 1 / **5** | **no**, **no** |
| CROSS | 1,264 / 1,227 | 2 / 1 | **no**, yes | 1,085 / 1,071 | 2 / **5** | **no**, **no** |
| OUTSIDE | 868 / 1,069 | 177 / 3 | yes, yes | 740 / 944 | 174 / 8 | yes, yes |
| REMOVED | 1,038 / 1,013 | 11 / 3 | **no**, yes | 917 / 964 | 4 / 0 | yes, yes |
| BRANCH | 1,181 / 1,301 | 219 / 55 | yes, yes | 994 / 1,108 | 180 / 67 | yes, yes |
| SILENT | 698 / 686 | 2 / 0 | **no**, yes | 552 / 542 | 1 / **2** | **no**, **no** |

**Where the integrated pipeline wins.** On a fresh set of full families with realistic reads, it is far ahead wherever the
catalog is incomplete:
- OUTSIDE: +201 and +204 correct, 3 and 8 errors against 177 and 174;
- BRANCH: +120 and +114 correct, 55 and 67 errors against 219 and 180;
- REMOVED (public): +47 correct, 0 errors against 4.

**Where it loses: catalogued copies.** It gets 1-3.5% fewer reads right in INTACT, CROSS and SILENT, and in REMOVED under
A119b.
- **Fewer correct reads:** the losses are reads O2 leaves `unexplained`, spread over the families. For INTACT under A119b:
  AGAP 15, KRTAP4 12, LCE2 11, FCGR1 9, PPIAL4 5. The current method assigns them.
- **The public profile's extra errors:**
  - INTACT's 5 are all PGA reads (pepsinogen A3/A4/A5, near-identical copies) called to the wrong copy by O2;
  - CROSS and SILENT were not broken down.
- This is the same mechanism development set 8 showed under A119b (INTACT -19).

**Cost.** 168 s at 0.76 GB peak (A119b) and 149 s at 0.74 GB (public) for six worlds of 43 copies each. The local reference
was built once, in the first run.

## Reading

**The pre-registered claim fails under both profiles.** No partial claim is made.

Three determinations hold on the held-out set as on development data:
- **Realistic reads change the comparison:** the advantage that idealized reads showed on catalogued copies does not survive.
- **The integrated pipeline's strength is the incomplete catalog:** copies the catalog lacks (OUTSIDE), copies absent from the
  reference with a novel exon (BRANCH), and their errors.
- **The remaining weakness is O2's divergence gate on realistic reads.** These are reads with retained introns, long indels or
  short spans that no annotated path explains. A genome-placement fallback for them was measured on development data and
  rejected for its errors.

## Independent review (2026-10-16) and corrections

**The verdict holds.**
- The reviewer's own scorer, written from this document's text, reproduces `e2e.json` and `e2e_pr.json` exactly for both
  profiles.
- The committed code rebuilds all 12 worlds byte for byte. Re-running O1, O2 (`gapcomp`) and O3 (`--member-genome-check`) on
  copies reproduces their outputs.
- The families share no gene with any earlier set. The seed came from the pre-check before any run. The final runs started
  five seconds after the amendment commit.

**Corrections to the results above:**

1. **The public profile's extra errors are an O1 step, not an O2 one.**
   - All 12 integrated errors under the public profile (INTACT 5, CROSS 5, SILENT 2) are reads of PGA.2 (PGA5).
   - O1 dropped PGA.2: it had 5 reads but only 4 primaries on its locus, 2 of them uniquely placed, so it fell under the 5-read
     support rule.
   - O1 dropped 2-4 expressed copies with 4-5 reads each in every public-profile world (AGAP.0, LCE2.1, PGA.2, PPIAL4.3), and
     none under A119b.
   - Those drops account for 10 of INTACT's 18-read deficit, 6 of CROSS's 14, and all 10 of SILENT's.
   - The verdict still fails through B1 in INTACT and CROSS, and through every A119b failure.
   - Cause: the public profile's expression spread (1.16) was inflated by zero-read copies and Poisson counting noise (item 3).
     It produced copies with 4-5 reads, and O1's fixed support rule drops them.
2. **The pre-registered side reports:**
   - **Arm I:** no false call in any of the 12 worlds.
   - **REMOVED:** under A119b arm I called 2 of the 4 detectable copies (LCE2.0, KRTAP4.0, plus PPIAL4.4). It missed AGAP.5 and
     FAM25.2, whose median `de` is 0.0149 and 0.0164, although O2 left 54 of 57 and 21 of 21 of their reads unexplained.
     Under the public profile it called all four.
   - **BRANCH:** 6 of 6 novel-exon copies under A119b, 5 of 6 under the public profile.
   - **O1:** dropped 8 of 13 silent copies under both profiles; the expressed-copy drops are in item 1.
   - **"REMOVED losses are reads O2 leaves unexplained" is not the whole story.** The current method caught FAM25.2, which
     arm I missed; that accounts for 15 of the 55 correct reads lost in A119b's REMOVED.
3. **The expression spread.** Re-estimated over copies with reads, less Poisson noise
   (`calibrate.py` `expression_log_sd`): A119b 0.86 (unchanged), public 1.03 (not 1.16).
4. **Cost.** The per-run peaks (0.74-0.76 GB) leave out the one-time local reference build. That build loads the whole-genome
   index, 38-75 s at about 14 GB, once per family set. A new family set therefore still peaks near 14 GB, against 19.4 GB
   before. Round 9's own build ran inside the deleted crashed run.
5. **"Realistic reads change the comparison" was not tested on this set:** there was no idealized-read run of set 10. A post
   hoc idealized run is reported in the development doc.
6. **Smaller points:**
   - The local reference was validated against whole-genome mapping on idealized reads only. minimap2's repeat filter differs
     on the local reference (99 vs 969 occurrences), so realistic reads are untested there.
   - minimap2 breaks MAPQ-0 ties partly on the read NAME, so mapping by sequence moves tie choices; the bars were identical.
   - Before amendment 1, the crashed runs' current, O1 and O2 outputs existed and were deleted. The order retry cannot be tuned
     to the bars.
   - Per-read detectability partly reflects sequencing noise.
   - A BRANCH read without the novel exon still counts as detectable.
