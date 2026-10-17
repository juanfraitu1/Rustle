# Round 13: the integrated pipeline v5', independent of the current method, on a fourteenth set of full families (2026-10-17)

Pre-registration, written before any world of the fourteenth set was built: no read of these families has been simulated.
Code: branch `o3-fast` (`bench/o3_graph/`, `bench/o3_unexplained/psvclust.py`), frozen at the commit that adds this file.

**Why.** The user's goal: integrate O1, O2 and O3 so that full families are fast and cheap to test, with realistic
simulations. Round 12 (register 1310) failed under both profiles. Its review found two problems:
- **Most of its error reads are structural.** Retained introns are absorbed as one gap-compressed edit, and two copies are
  DNA-identical.
- **O3's cost was set by the branch arm,** whose all-against-all alignment ran away when GraphAligner failed on the 80-copy
  DUX4 graph.

v5' removes the cost problem. It does not remove round 12's O2 errors (see "Known before running").

**Question.** On a fresh set of full families, including a 43-copy pseudogene array and retrocopy families, does v5' do at least
as well as the current method in every world type, with both calibrated v2 profiles? That means at least as many correct read
outcomes and no more wrong ones.

## The two systems

**Current** (unchanged since round 1): copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`)
on the curated catalog, plus its `--psv-clusters` O3.

**Integrated v5'** = round 12's v5 + `--branch-genome-novel` (16ba2508):
- **O1:** `o1seed.py W --twins`.
- **O2:** `o2graph.py run W integ 2 --distance gapcomp --union-exact --error-floor`, then `run W integl 2 --error-floor`.
- **O3:** `detect.py W --member-genome-check --member-attribution --error-floor --o1-hosts --branch-genome-novel --arms I`.
  - It is hosted by O1 and reads nothing of the current method's run.
  - The branch arm clusters only reads whose own genome alignment leaves >= `MIN_NOVEL` bases unexplained. The genome check
    applies the same criterion to every call afterwards. A read the genome explains carries no novel exon, and this bounds
    the all-against-all when GraphAligner fails on a graph.
- **Measured and not adopted:**
  - **`--distance gapstruct`** (gap compression only below 30 bp). It still fails round 12, adds an error in set 8 A119b
    INTACT and SILENT, and costs set 10 A119b REMOVED 3 correct reads.
  - **`--mapped-sites`.**

```
fastset.py run OUT FAMILIES SET SEED --o2-distance gapcomp --member-genome-check --member-attribution --union-exact
    --o2-genome fallback --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --heavy-current
    --step-timeout 1800 --het 0 --jobs 2 --threads 2 --sim-profile PROFILE
```
`--heavy-current` and `--step-timeout 1800` are set from the start: round 12's 80-copy array needed both, and outputs depend on
neither.

## Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout13`, 5c1f3ab9):
- **MRGPRX:** a round-12 candidate with >= 3 copies past its first 8, never built into a world.
- **S100A7, KRTAP12, KRTAP13, KRTAP20, KRTAP21, KRTAP6, HLA-ABC, H1, OR2A, SEPTIN7P, CTSLP, UBE2Q2P, PGAM1P, OR52, OR7E:** in
  no earlier list.
- **Spares, appended before any run (bb60c046)**, because the first build gave only 5 families with >= 3 copies: YWHAZP,
  PTMAP, TPT1P, FTH1P, EEF1A1P, KRT8P (retrocopy families, smallest first).

**Build rule:** round 12's, `families.py build families13_all.json --cap 1000 --link locus --min-cover 0.8`, every group. The set
is the first 8 families with >= 3 copies, in list order:

| family | copies | isoforms | chromosomes |
|---|---|---|---|
| MRGPRX | 4 | 5 | chr11 |
| SEPTIN7P | 5 | 8 | chr7, chr9 |
| CTSLP | 6 | 8 | chr9, chr10 |
| UBE2Q2P | 4 | 8 | chr15 |
| OR7E | 43 | 44 | chr2, chr3, chr4, chr7, chr8, chr9, chr10, chr11, chr13, chr14 |
| YWHAZP | 5 | 5 | chr2, chr6, chr9, chr10 |
| PTMAP | 9 | 9 | chr3, chr5, chr7, chr9, chr12, chr13, chr20 |
| FTH1P | 13 | 13 | chr1, chr2, chr3, chr5, chr6, chr8, chr9, chr11, chr13, chrX |

89 copies. Files: `/mnt/linuxdisk/tmp/o3_graph/families13.json` (this set only), `families13_all.json`, `families13.log`.

## Worlds

Six world types (INTACT, CROSS with liftover projection, OUTSIDE, REMOVED, BRANCH, SILENT), each under the two v2 profiles,
`profiles/A119b_testis_v2.json` and `profiles/public_testis_v2.json` (`--het 0`). That makes 12 worlds.

**The seed:** `precheck.py families13.json heldout13 181 191` (`precheck13.txt`; 70.5 s, 14.5 GB) gives the smallest seed in
181-190 with >= 3 detectable removed copies: **seed 181**.
- Detectable: MRGPRX 0.0809, YWHAZP 0.0597, OR7E 0.051, SEPTIN7P 0.0436, FTH1P 0.0435, UBE2Q2P 0.0433, CTSLP 0.0102.
- Within delta: PTMAP 0.0061.

**Commands** (one per profile; the first run builds the local reference):
```
F=/mnt/linuxdisk/tmp/o3_graph/families13.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout13
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome fallback --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2"
python3 bench/o3_graph/fastset.py run $O/A119b $F heldout13 181 --local $O/local $V --sim-profile bench/o3_graph/profiles/A119b_testis_v2.json
python3 bench/o3_graph/fastset.py run $O/public $F heldout13 181 --local $O/local $V --sim-profile bench/o3_graph/profiles/public_testis_v2.json
```

## Read outcomes, bars and reading

These are round 12's, unchanged:
- **Read outcomes:** per-read detectability less the world's error floor. Round 10's rule and copy level are secondary.
- **Bars:** B1 (integrated correct >= current) and B2 (integrated errors <= current), per world, all 12 worlds.
- **PASS** (B1 and B2 hold in all 12 worlds). Claim: "the integrated O1 -> O2 -> O3 pipeline v5', which uses nothing from the
  current method, gives in every world type of a fresh set of full families at least as many correct read outcomes as the
  current method and no more wrong ones, with reads calibrated on two real libraries".
- **FAIL** otherwise. No partial claim.

## Known before running (development; v5' on 14 spent runs, `fast/r13dev2/`)

**v5' passes 12 of 14:** sets 8-10 under both profiles, and rounds 10 and 11's six runs.
- **Round 12's two runs still fail**, by their O2 errors, which v5' does not address:
  - A119b: INTACT 3, CROSS 5, SILENT 2 against 0;
  - public: CROSS 2 against 0.
- **Against v5, the branch prefilter:**
  - **BRANCH in round 12:** correct +78 and +75, errors 136 -> 79 and 161 -> 99. DUX4's novel-exon copy is no longer lost in
    giant mixed clusters.
  - **Elsewhere:**
    - REMOVED moves by at most 16 reads (round 11 A119b: 1,845 -> 1,829, still far above current's 1,545);
    - BRANCH by up to 6;
    - 6 of the 14 runs are identical.
  - **O3's longest world:** at most 78 s, against 699 s before (round 12, public OUTSIDE).
- **Round 12's open O2 cases, with no fix in v5':**
  - copies DNA-identical over a read with a retained intron;
  - retained introns absorbed by a paralog's path as one gap-compressed edit;
  - margin-2 coincidences in an 80-copy array, which the genome's own alignment also misplaces.
- **Untested in development:** pseudogene and retrocopy families of this kind, and a 43-copy array.

## Reported beside the bars (not gating)

- fates per read class and error sources, **tabulated, every world, both profiles**;
- O1's drops (of how many silent copies), undecidable copies, new copies and twins;
- arm I's calls, true or false, every world;
- round 10's scoring and copy-level scoring;
- leave-one-family-out margins, **both profiles**;
- the cost: `timing.json`, the integrated steps apart from the current method's, and O3's peak memory.

## After the run

- An independent review by a fresh subagent, given read-only copies.
- Register row 1311.
- Results appended below a line, without editing anything above it.
