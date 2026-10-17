# Round 12: the integrated pipeline v5, independent of the current method, on a thirteenth set of full families (2026-10-16)

Pre-registration, written before any world of the thirteenth set was built: no read of these families has been simulated.
Code: branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.** The user's goal is to integrate O1, O2 and O3 so that full families are fast and cheap to test. Round 11 (register
1309) passed with v4, whose O3 still took its read hosts from the current method's run (copy_assign,
`cur/A.assignments.tsv`). v5 hosts O3's reads from O1, so the integrated pipeline O1 -> O2 -> O3 needs nothing from the
current method; the current method runs only as the comparison.

**Question.** On a fresh set of full families, including an 80-copy array, does v5 do at least as well as the current method
in every world type, with both calibrated v2 profiles? That means at least as many correct read outcomes and no more wrong
ones.

## The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v5** = round 11's v4, with O3 hosted by O1 (`--o1-hosts`; f45e0fe5, f4f2b248):
- **O1:** `o1seed.py W --twins` (round 11).
- **O2:** `o2graph.py run W integ 2 --distance gapcomp --union-exact --error-floor`, then `run W integl 2 --error-floor`
  (round 11).
- **O3:** `detect.py W --member-genome-check --member-attribution --error-floor --o1-hosts --arms I`.
  - **Hosts:** O1's seeded copies, new ones included and dropped ones excluded.
  - **Reads:** every O1 family read is hosted on the seeded copy its primary alignment overlaps most. Failing that, on the
    copy a secondary alignment overlaps most; failing that, on the copy whose +-1 kb window its alignments overlap.
  - **No file of the current method's run is read.**
- **Since round 11's review:**
  - O3 opens one FASTA handle per process (7fbc236d; review I1). Two-worker and one-worker runs are byte-identical.
  - `--mapped-sites` (981708b5) was measured and is not adopted. It changes no bar and leaves the pre-check false calls as
    they were: these are cross-family branch calls and identical-copy clusters, all removed by the genome check.
- **Development: v5 passes all 12 runs under the primary scoring** (`fast/r12dev_base/`).
  - **The runs:** sets 8, 9 and 10 under both profiles, round 10's three runs and round 11's three, all spent.
  - **Round 10's rule:** v5 passes 8 of 12. The 4 failures are REMOVED B2 in sets 8 public, 9 A119b, 10 A119b and 10
    public, as for v4.
  - **Against v4 (copy_assign hosting), 72 world bars are identical in 68.** The four that differ:
    - round 11 A119b REMOVED 1,844 -> 1,845 (v4's run hit the FASTA race) and BRANCH 2,040 -> 2,041;
    - round 10 A119b BRANCH 1,043 -> 1,042;
    - set 8 A119b BRANCH 2,536 -> 2,538.
  - The -25 BRANCH reads seen earlier in set 9 public (`fast/r11dev7/`) do not recur with one FASTA handle per process.
  - Before the window fallback (f4f2b248), hosting only primaries lost array reads and passed 6 of 9.

The command:
```
fastset.py run OUT FAMILIES SET SEED --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
    --o2-genome fallback --per-read-detectability --twins --error-floor --o1-hosts --het 0 --jobs 2 --threads 2
    --sim-profile PROFILE
```
On development set 10 (A119b) this one command passes. The integrated steps take 118 s for six worlds (O1 7 s, O2 27 s,
genome fallback 1 s, O3 83 s), and the current method 47 s (`fast/r12dev_fastset/`).

## Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout12`, committed in 35bf1d95):
- KRTAP2, KRTAP1, ETD, CDRT15, NAIP, HLA-DRB, OR11H, DUX4: round-11 candidates past its first 8, never built into a world;
- POLR2J, KRTAP19, KRTAP3, MRGPRX: in no earlier list.

**Build rule:** round 11's, `families.py build families12_all.json --cap 1000 --link locus --min-cover 0.8`, every group. The
set is the first 8 families with >= 3 copies, in list order:

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| KRTAP2 | 4 | 5 | chr17 |
| KRTAP1 | 3 | 3 | chr17 |
| ETD | 3 | 5 | chrX |
| CDRT15 | 4 | 6 | chr17 |
| NAIP | 4 | 8 | chr5 |
| HLA-DRB | 3 | 3 | chr6 |
| OR11H | 3 | 3 | chr14, chr22 |
| DUX4 | 80 | 85 | chr1, chr3, chr4, chr9, chr10, chr16, chr20, chrY |

104 copies. DUX4's copies sit in D4Z4 tandem arrays, so neighbouring copies fall inside one another's +-2 kb host windows
(round 10's artefact). Files: `/mnt/linuxdisk/tmp/o3_graph/families12.json` (this set only), `families12_all.json`,
`families12.log`.

## Worlds

Six world types (INTACT, CROSS with liftover projection, OUTSIDE, REMOVED, BRANCH, SILENT), each under round 11's two
profiles, `profiles/A119b_testis_v2.json` and `profiles/public_testis_v2.json` (`--het 0`). That makes 12 worlds.

**The seed:** `precheck.py families12.json heldout12 161 171` (`precheck12.txt`; 76.5 s, 13.9 GB) gives the smallest seed in
161-170 with >= 3 detectable removed copies: **seed 161**.
- Detectable: KRTAP1 0.0514, HLA-DRB 0.0498, DUX4 0.0275, KRTAP2 0.0129.
- Within delta: CDRT15 0.0091, OR11H 0.0061, ETD 0.0032, NAIP 0.0.

**Commands** (one per profile; the first run builds the local reference):
```
F=/mnt/linuxdisk/tmp/o3_graph/families12.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout12
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome fallback --per-read-detectability --twins --error-floor --o1-hosts --heavy-current --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout12 161 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis_v2.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout12 161 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis_v2.json
```
- **Copy_assign's memory:** the current method runs under the heavy lock (`--heavy-current`) from the start. An 80-copy
  array may exceed the light jobs' 2 GB, and its output does not depend on the lock.
- **An 8-minute step timeout:** if a step exceeds rlock's 590 s timeout, the run is repeated with that step's timeout raised.
  The outputs do not depend on it. Either repeat is an amendment, recorded before any result is read.

## Read outcomes

Round 11's, unchanged:
- per-read detectability less the world's error floor (`summary.py e2e --error-floor`);
- a BRANCH read is detectable by construction;
- an O3 call on a read takes precedence over its O2 placement;
- a read called in several families is judged by its call in its own family.

**Secondary:** round 10's rule without the floor, and copy-level scoring.

## Bars, per world, all 12 worlds

- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

Abstain and unknowable count for neither.

## Reading

- **PASS** (B1 and B2 hold in all 12 worlds under the primary scoring). Claim: "the integrated O1 -> O2 -> O3 pipeline v5, which
  uses nothing from the current method, gives in every world type of a fresh set of full families, including an 80-copy
  array, at least as many correct read outcomes as the current method and no more wrong ones, with reads calibrated on two
  real libraries".
- **FAIL** otherwise. No partial claim; the failing worlds and margins are reported.

## Reported beside the bars (not gating)

- fates per read class, tabulated, and the error sources;
- O1's drops (of how many silent copies), kept-undecidable copies, new copies and twins;
- arm I's calls, true or false, in every world;
- round 10's scoring and copy-level scoring;
- **margins with each family left out in turn**, so that no single family carries the verdict (round 11's histone
  clusters did, for the margins);
- **the run's cost:** `timing.json`, with the integrated pipeline's own steps apart from the current method's.

## After the run

- An independent review by a fresh subagent, given read-only copies.
- Register row 1310.
- Results appended below a line, without editing anything above it.

---

# Results (2026-10-16; run with the code of 00c241db, started 22:29)

## Amendment 1 (23:12, recorded after the A119b run's result was read)

**The public run crashed in its OUTSIDE world's O3 step.**
- The step outlasted rlock's 590 s limit (DUX4, 80 copies) and was killed before `/usr/bin/time` wrote its line.
- `fastset.py` then failed to parse the empty timing file and stopped (exit 1). No `e2e.json` or `bars.txt` was written for
  the public profile.

**As pre-registered for a step over the limit, the public run is repeated with the limit raised:**
- `--step-timeout 1800`, a new driver option; outputs do not depend on it;
- with the driver's parsing fixed (`parse_time`), so a killed step is reported, not crashed on;
- into `heldout12/public_r2`, with the crashed run's directory kept.

**This amendment comes after the A119b result was read,** which the pre-registration's amendment rule did not foresee.
A119b already FAILs (below), so the round's verdict is FAIL whatever the public repeat gives. The repeat completes the
record and is reported as such.

## Verdict: FAIL under both profiles. Not claimed.

Bars from `e2e.json` (primary scoring), A119b and the public repeat (`public_r2`, Amendment 1):

| world | A119b: correct (current / integrated) | errors | public: correct | errors |
|---|---|---|---|---|
| INTACT | 2,189 / 2,380 | 0 / **3** | 2,171 / 2,399 | 0 / 0 |
| CROSS | 2,250 / 2,494 | 0 / **5** | 2,237 / 2,508 | 0 / **2** |
| OUTSIDE | 1,683 / 2,289 | 435 / 113 | 1,622 / 2,335 | 500 / 134 |
| REMOVED | 1,799 / 2,149 | 16 / 5 | 1,765 / 2,175 | 32 / 9 |
| BRANCH | 2,189 / 2,597 | 208 / 136 | 2,171 / 2,578 | 217 / 161 |
| SILENT | 1,415 / 1,695 | 0 / **2** | 1,339 / 1,695 | 0 / 0 |

- **B1 holds in all 12 worlds,** by +191 to +713.
- **B2 fails in four worlds,** where the current method makes no error: A119b INTACT, CROSS and SILENT, and public CROSS.
- **Round 10's rule and copy level give the same verdict, failing the same cells.**
- **No single family carries the failure.** Rescored with each family left out in turn, A119b fails every time; without
  DUX4 its errors are INTACT 1, CROSS 3, SILENT 1 against 0.

**The errors are O2 assignments at margin 2 to a near-identical wrong copy, never O3:**
- **DUX4 (80 copies):** a read's sequencing edits coincide with another unit's differences.
  - A DUX4.9 read is 2 edits from DUX4.19 and 4 from its own copy.
  - A DUX4.2 read matches DUX4.22 exactly and its own copy at 2.
  - With 80 near-identical competitors, chance coincidences beat a fixed margin of 2.
- **KRTAP2 (INTACT):** a KRTAP2.0 read aligns perfectly to the genome but extends past KRTAP2.0's annotated end, and only
  KRTAP2.1's longer isoform covers it. Its path distances are 2 against 85: structure decides, not sequence.
- **CDRT15 (CROSS):** CDRT15.0 and CDRT15.2 reads express structures annotated only on CDRT15.1 and CDRT15.3. The union paths
  do not give them those structures, as in round 10's WASH case.
- **The current method abstains on all of them.**

**Reported beside the bars:**
- **Arm I:**
  - **REMOVED:** all 4 detectable removed copies called under both profiles (DUX4.1, HLA-DRB.1, KRTAP1.1, KRTAP2.1), and none
    of the 4 within delta.
  - **BRANCH:** novel-exon copies of CDRT15, ETD, HLA-DRB, KRTAP2 and NAIP called under both profiles.
  - **OUTSIDE:** 3 (A119b) and 2 (public) PSV calls in KRTAP2 on present or uncatalogued reads (55 and 48 reads), where B2 still
    holds.
- **O1:**
  - **SILENT:** 6 of 33 silent copies dropped; 2-4 kept as undecidable (DUX4, ETD, NAIP, OR11H);
  - **OUTSIDE:** 5 new copies (DUX4, CDRT15, ETD, HLA-DRB, KRTAP2), two of them twins.
- **Cost: v5 is independent of the current method but not cheaper on an 80-copy array.**
  - **The integrated steps, summed over six worlds:** 1,095 s (A119b) and 1,624 s (public).
    - O3: 791 and 1,355 s, with one world's O3 alone over 590 s; it hosts reads on all 80 seeded units, and each host aligns
      the family pool.
    - O2: 174-208 s. O1: 90 s.
  - **The current method:** 383 and 356 s, peaking at 6.8-7.0 GB per world.
  - **Wall:** 1,098 s and 1,257 s per profile.

## Independent review (2026-10-16/17, from 23:32) and corrections

Fresh subagent; full report `/mnt/linuxdisk/tmp/o3_graph/review12/REVIEW.md`.

**The FAIL stands.**
- The reviewer's own scorer reproduces every bar, cell for cell, under all 3 rules.
- 00c241db reproduces O1, O2 and the fallback byte-identically.
- **The protocol is clean apart from Amendment 1, which could not change anything.** The crashed public run and its repeat
  are byte-identical in 71 of 72 outputs; the 72nd was never written by the crash.
- **Amendment 1's timing:** it was committed at 23:09:51, not at its 23:12 heading. A119b's bars and its error reads had been
  read first, because one tool output delivered both.
- **v5 uses nothing from the current method** (confirmed).

**Corrections:**
1. **The failure's mechanism was wrong for 5 of its 7 distinct reads (Critical).** The 12 failing error instances are 7 reads,
   repeated across worlds that share their reads.
   - **DUX4 (2 reads):** margin-2 coincidences, as stated. The genome's own primary alignment also puts both reads on the
     wrong unit (MAPQ 28-29).
   - **KRTAP2 (1 read):** it does not run past an annotated end.
     - It is KRTAP2.0's isoform with its intron RETAINED, DNA-identical at KRTAP2.0 and KRTAP2.1 (two NM-0 alignments), so
       abstaining is the only right answer.
     - Its 2 vs 85: KRTAP2.0's spliced path clips 83 bases, and a clip over 50 costs its length. KRTAP2.1's long-exon
       isoform absorbs the read's tandem-repeat difference as one 15-bp insertion, which costs one event.
   - **CDRT15, A119b (2 reads):** retained introns again, not missing structures.
     - The union does give CDRT15.0 the read's structure. The paralog's path absorbs the 265-bp retained intron as ONE
       142-bp insertion (gap-compressed 9), while CDRT15.0's own paths clip it at full cost.
     - One read went to a sole candidate: an O1-observed chain on CDRT15.1, learned from CROSS reads.
     - **The error floor contributes:** both reads sit at 0.0100 per base, between delta and delta + floor. Under delta
       alone they would be `unexplained` and placed correctly by the genome fallback.
   - **CDRT15, public (2 reads):** structures absent from the union, as stated.
   - **The current method abstains on the DUX4 and KRTAP2 reads, but places all 4 CDRT15 reads correctly.**
   - Retained introns cause 3 of A119b's 5 error reads. A119b retains introns at 6.4x the public rate.
2. **Leave one family out (Important):** public fails only through CDRT15. Without HLA-DRB, BRANCH B2 fails under both
   profiles.
3. **DUX4's novel-exon copy is missed, unreported (Important).** O2 places 57 and 62 of its reads that CARRY the novel exon,
   because a 150-bp novel exon costs one gap-compressed edit.
4. **The cost mechanism was wrong (Important).**
   - GraphAligner aligned 4% of the DUX4 reads in public OUTSIDE (61-96% elsewhere). 3.44 Mbp of whole reads became "novel
     segments".
   - The branch arm's all-against-all minimap2 then took 553 s and 1.7 GB, producing a 3,974-read call that the genome check
     removed.
   - O3's time tracks this novel-segment volume, not the 80 hosts, and peaked at 6.15 GB under the light lock.
5. **Arm I's KRTAP2.1 "detection" in REMOVED is identical to its OUTSIDE false calls (Important).**
   - KRTAP2.1's 52-57-bp exons in 15-bp tandem repeats mis-align even when the copy is present (INTACT median `de` 0.0137).
   - The OUTSIDE call members are all uncatalogued KRTAP2.1 reads: 35 and 24 distinct reads, not 55 and 48.
6. **The zero O3 errors in INTACT, CROSS and SILENT rest on the genome check alone (Important).** Before it, arm I's calls
   hold 1,203-1,462 present-copy reads per world, mostly DUX4 branch clusters of reads GraphAligner could not align.
7. **Minor:**
   - **O1 seeded a false new DUX4 copy in REMOVED** from the absent copy's reads (harmless to the bars).
   - **KRTAP1.1, NAIP.0 and OR11H.2 are unseeded outside copies.**
   - **The two twins keep few of their own reads.**
   - **The current method peaked at 6.2-6.5 GB under public**, not 6.8-7.0.
   - **The fates per read class are tabulated only in the review**, for the third round running.

**What this changes:**
- **Gap compression conflates sequencing indels with structure.** A retained intron or a novel exon becomes one edit when
  minimap2 bridges it, while the same sequence costs its full length when minimap2 clips it. That asymmetry, not near-identity,
  decided 5 of the 7 error reads, and it hides novel exons from O3.
- **The branch arm has no defence against GraphAligner failing on a graph;** its all-against-all then sets O3's cost.
