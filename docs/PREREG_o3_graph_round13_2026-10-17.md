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

---

# Results (2026-10-17)

Both runs completed as pre-registered, with no amendment. Outputs are in `fast/heldout13/{A119b,public}/` (frozen read-only after
scoring), wall 417.8 s (A119b, which also built the local reference) and 325.8 s (public).

## Verdict: FAIL (A119b). Not claimed.

Bars from `e2e.json` (primary scoring):

| world | A119b: correct (current / integrated) | errors | public: correct | errors |
|---|---|---|---|---|
| INTACT | 4,709 / 4,737 | 0 / **3** | 4,622 / 4,671 | 0 / 0 |
| CROSS | 4,783 / 4,860 | 2 / **3** | 4,690 / 4,781 | 3 / 0 |
| OUTSIDE | 3,990 / 4,383 | 683 / 16 | 3,817 / 4,230 | 786 / 8 |
| REMOVED | 4,657 / 4,682 | 10 / 0 | 4,467 / 4,618 | 95 / 1 |
| BRANCH | 4,718 / 4,888 | 95 / 63 | 4,622 / 4,785 | 102 / 72 |
| SILENT | 3,237 / 3,257 | 0 / **3** | 3,214 / 3,246 | 0 / 0 |

- **public passes all 6 worlds.**
- **A119b:**
  - **B1 holds in all 6 worlds,** by +20 to +393.
  - **B2 fails in three:** INTACT 3 vs 0, CROSS 3 vs 2, SILENT 3 vs 0.
- **Round 10's rule and copy level give the same verdict** (A119b FAIL, public PASS).

**The failure is one false O3 call, repeated in four worlds (found after the run, post hoc).** Every A119b error is arm I's PSV
call in MRGPRX on three reads of the present copy MRGPRX.1 (MRGPRX4). O2 makes no wrong-copy assignment in any of the 12 worlds.
- **The same three reads in every world.** The simulator draws an unchanged copy's reads from the same seed, so INTACT, CROSS,
  OUTSIDE and SILENT carry the same three reads under different names. No call is made in REMOVED, where MRGPRX.0 (whose window
  hosts the call) is the removed copy. BRANCH has the three reads too and makes no call; why was not traced.
- **They carry MRGPRX.1's heterozygous alleles.**
  - With a profile, `--het 0` means the profile's calibrated heterozygosity (`realworld.py`: `het_rate = 2 * het_per_base`),
    not a haploid donor.
  - MRGPRX.1's haplotype 2 carries two alleles 75 bp apart: C>G at chr11:18,269,598 and C>T at 18,269,673. 45 haplotype-2
    reads carry both, and 43 haplotype-1 reads carry neither.
- **Short reads turn two alleles into a delta-sized divergence.**
  - Of the 45 haplotype-2 reads covering both sites, O2 places 42 on MRGPRX.1. It leaves 3 unexplained: the two shortest (190
    and 231 bp) and a 506-bp read with 6 sequencing indels besides the two alleles.
  - Their primary alignments (MAPQ 60, on MRGPRX.1) have `de` 0.0158, 0.0129 and 0.0157. NM is 3, 3 and 11. That puts them in
    the world's top 0.2% (the median MAPQ-60 `de` is 0.0029, p99 0.0114).
  - **Both O2 stages leave them unexplained,** above delta plus floor: the paths and the genome fallback. Arm I clusters O2's
    unexplained reads, and without the call the three reads would be abstentions.
  - **The call:** arm I pools the family's unexplained reads across hosts; the three form a cluster on MRGPRX.0's window
    (`est_div` 0.024).
  - **Both checks pass:**
    - attribution keeps all three reads (its internals were not traced);
    - genome check: median `de` 0.0157 less the floor 0.0029 = 0.0128 > delta, and each read alone also passes the per-read check
      (the smallest is 0.0129 - 0.0029 = 0.0100).
  - The two alleles alone add 2/190 and 2/231 to `de` (0.0105 and 0.0087).
- **Delta is calibrated on whole copies.** It is the p99 exonic divergence between two haplotype assemblies over whole copies
  (register 1199). Two heterozygous sites within 200 bp of a short read exceed it, so the per-read genome check cannot separate an
  allele from a copy on short reads. This is the allele-versus-copy confound of the phased toy (2026-10-12). Round 10 met it too:
  TAF11L.33's haplotype-2 reads carried two donor alleles (median `de` 0.0103), and the error floor removed them (v4). Here the
  reads are short enough that two alleles exceed delta even after the floor.
- **Leave one family out:** A119b passes only without MRGPRX. Public passes with any family left out.
  - **Without OR7E (the 43-copy array), A119b's B1 margins turn negative:** INTACT -5, REMOVED -11, SILENT -11.
  - Its correct-read margin over the current method in those worlds is carried by one family.

**Reported beside the bars:**
- **Arm I:**
  - **REMOVED:** all 7 detectable removed copies called under both profiles (MRGPRX.0, YWHAZP.3, OR7E.15, SEPTIN7P.0, FTH1P.8,
    UBE2Q2P.3, CTSLP.1); PTMAP.8, within delta, not called.
  - **BRANCH:** all 6 novel-exon copies called under both profiles. PTMAP and FTH1P are single-exon retrocopies, with no BRANCH copy.
  - **Calls on present reads:** none under public. Under A119b, only the MRGPRX call above, in INTACT, CROSS, OUTSIDE and SILENT.
- **O1:**
  - **SILENT:** 24 of 27 silent copies dropped under both profiles; public keeps UBE2Q2P.1 as undecidable.
  - **OUTSIDE:** new copies in CTSLP, MRGPRX, OR7E, SEPTIN7P, UBE2Q2P and YWHAZP, with no twin.
  - **REMOVED:** new copies in FTH1P, SEPTIN7P (2), UBE2Q2P and YWHAZP.
- **Fates per read class** (correct / error / abstain / unknowable, current -> integrated):

| world | A119b catalogued | A119b other class | public catalogued | public other class |
|---|---|---|---|---|
| INTACT | 4709/0/499/0 -> 4737/3/468/0 | | 4622/0/479/0 -> 4671/0/430/0 | |
| CROSS | 4783/2/565/0 -> 4860/3/487/0 | | 4690/3/534/0 -> 4781/0/446/0 | |
| OUTSIDE | 3990/5/483/0 -> 4017/3/458/0 | uncatalogued 0/678/52/0 -> 366/13/351/0 | 3817/77/413/0 -> 3880/0/427/0 | uncatalogued 0/709/85/0 -> 350/8/436/0 |
| REMOVED | 3990/5/483/0 -> 4022/0/456/0 | absent detectable 636/5/43/0 -> 652/0/32/0; undetectable 31/0/2/13 -> 8/0/4/34 | 3817/77/413/0 -> 3886/0/421/0 | absent detectable 649/18/98/0 -> 727/1/37/0; undetectable 1/0/9/19 -> 5/0/4/20 |
| BRANCH | 4702/7/499/0 -> 4739/0/469/0 | absent detectable 16/88/119/0 -> 149/63/11/0 | 4622/0/479/0 -> 4671/0/430/0 | absent detectable 0/102/90/0 -> 114/72/6/0 |
| SILENT | 3237/0/311/0 -> 3257/3/288/0 | | 3214/0/304/0 -> 3246/0/272/0 | |

- **Cost: the branch prefilter holds on a 43-copy array.**
  - **The integrated steps, summed over six worlds:** 223 s (A119b) and 216 s (public), against round 12's 1,095 and 1,624 s.
    - O3: 164 and 160 s, longest world 72 s and 67 s (REMOVED), peak 0.17 GB.
    - O2 28-31 s (+ genome fallback 3-4 s); O1 22-27 s.
  - **The current method:** 188 and 183 s, peaking at 3.2 GB per world.
  - **Shared steps (simulate, map, mask):** 114 and 90 s. The local reference build (55 s, 14.4 GB) is paid once per set.
