# Round 8: the integrated O1 -> O2 -> O3, end to end, on a ninth handpicked family set (2026-10-16)

Pre-registration, written before any world of the ninth set was built: no read of these families has been simulated. Code:
branch `o3-e2e` (`bench/o3_graph/`), frozen at the commit that adds this file.

**User request (2026-10-16):** "fix the bug so an added copy belongs to only one family. Then pre-register an end-to-end test
on a fresh 9th family set, with bars on read outcomes per world type (correct ≥ current, errors ≤ current)".

**Question.** On a fresh family set, does the integrated pipeline do at least as well as the current method in every world
type? That means at least as many reads with a correct outcome (on their own copy, or inside a call for their absent copy) and
no more reads with a wrong one.

## The two systems

**Current** (unchanged since round 1):
- copy_assign's O2 (`run.py cur`, binary `/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`) on the curated catalog;
- its `--psv-clusters` O3: a cluster with status `unexplained` and estimated divergence > delta (0.00958) is an absent-copy
  call.

**Integrated** (round 7 plus the one-family rule):
- **O1**, `o1seed.py`:
  - keep each curated copy on which >= 5 family reads have their primary alignment (+-1 kb), and drop the others (silent);
  - add off-catalog loci with >= 5 uniquely placed reads;
  - give each copy the intron chains of near-unique reads (mapq >= 50) that it owns;
  - **new (1e3960fe):** a new locus found in several families stays only in the family whose curated copies its reads align
    to most (aligned bases, unpadded copy spans). A tie seeds it in none.
- **O2**, `o2graph.py run W integ 2`: padded union O2 on the seeded catalog.
- **O3**, `detect.py` arm I: the pooled clusters of reads that O2 leaves unexplained, with the mosaic rule, branch calls and
  the genome check.

Both run per world from one driver, `bench/o3_graph/round8.sh`.

## Families (held out, picked before any run)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `heldout8`, 1d391d5d):
- NOTCH2NL, SRGAP2, GH (GH1/2, CSH1/2, CSHL1), USP17L, TSPY, RBMY, BPY2, CDY, OR4F, LCE1;
- KRTAP5, PNMA6, PPIAL4, FCGR1, NOMO, CROCC, RGPD, FAM86, GPR89, POM121.

None of them was in any earlier world. FAM86, GPR89 and POM121 were round-4 candidates past its first 8.

The build rule is the same as in rounds 1-7:
- copies are linked at identity >= 0.90 over >= half the shorter canonical transcript;
- genes span <= 60 kb;
- a seeded sample of 10 copies is taken when a family has more.

NOTCH2NL and SRGAP2 form no family under that rule. The set is the first 8 with >= 3 copies:

| family | genes | copies (of linked) / isoforms |
|---|---|---|
| GH | CSH1, CSH2, CSHL1, GH2 (chr17) | 4 / 16 |
| USP17L | USP17L3, L4, L23, L25-27 (chr4, chr8) | 10 (of 59) / 10 |
| TSPY | TSPY1-4, 8, 10 (chrY) | 10 (of 35) / 38 |
| RBMY | RBMY1A1, 1D, 1E, 1F (chrY) | 5 / 16 |
| BPY2 | BPY2, BPY2B (chrY) | 3 / 12 |
| CDY | CDY1, CDY2A, CDY2B (chrY) | 4 / 4 |
| OR4F | OR4F4, OR4F5 (chr3, chr11, chr15) | 3 (of 8) / 3 |
| LCE1 | LCE1D, E, F (chr1) | 3 (of 6) / 5 |

Files: `/mnt/linuxdisk/tmp/o3_graph/families8.json` (round 7's file plus this set), `families8_all.json` and `families8.log`.
The rebuild reproduces every earlier set exactly.

**This is the hardest set so far.** Six of the eight families have near-identical copies:
- the four Y-ampliconic families;
- OR4F;
- the USP17L tandem array.

Their reads tie at mapq 0, the failure mode that no step of either system fixes.

## Worlds

All with heterozygous alleles (`--het 0.001`, as in round 7), seed **93**:
- INTACT;
- CROSS (`--projection liftover`);
- OUTSIDE;
- REMOVED (the special copy hard-masked, with its own masked genome and index);
- BRANCH (an extra absent copy with a novel 150-bp exon);
- SILENT.

**Seed (detectability pre-check, as in round 7).** `precheck.py families8.json heldout8 91 101` (`precheck8.txt`) gives seed 93,
the smallest seed in 91-100 with >= 3 detectable removed copies:
- detectable: GH 0.0367, USP17L 0.0163, LCE1 0.0391;
- within delta: TSPY 0, RBMY 0.0005, BPY2 0, CDY 0.0004, OR4F 0.

**Commands:**
```
R=/mnt/linuxdisk/tmp/o3_graph/real/heldout8; F=/mnt/linuxdisk/tmp/o3_graph/families8.json
python3 bench/o3_graph/realworld.py mask $R/masked $F --set heldout8 --seed 93
minimap2 -x splice:hq -t 4 -d $R/masked/masked.splice.mmi $R/masked/masked.fa
for sc in INTACT CROSS OUTSIDE REMOVED BRANCH SILENT; do bench/o3_graph/round8.sh $R $F heldout8 93 $R/masked $sc; done
python3 bench/o3_graph/summary.py e2e $R/e2e.json $R/INTACT $R/CROSS $R/OUTSIDE $R/REMOVED $R/BRANCH $R/SILENT
python3 bench/o3_graph/summary.py bars $R/e2e.json
```

## Read outcomes (frozen: `summary.py` `read_fate`, `cmd_e2e`)

**Class.** Each read is classed by the copy it came from:
- catalogued;
- uncatalogued but in the genome (OUTSIDE);
- absent and detectable: BRANCH's extra copy, or a REMOVED copy whose reads' median `de` on the masked genome is > delta;
- absent and undetectable.

**Fate.** Under each system, a read's fate is one of:
- **correct:**
  - inside an absent-copy call, for an absent copy's read;
  - otherwise, assigned to its own copy (an uncatalogued copy can only be its own once O1 has seeded it).
- **error:**
  - inside an absent-copy call, for a present copy's read;
  - assigned to another copy;
  - a detectable absent copy's read assigned to a present copy.
- **unknowable:** an undetectable absent copy's read assigned to a present copy.
- **abstain:** no call.

An absent-copy call takes precedence over an assignment.

**A read called in several families** is judged by its call in its own family, as in round 7 and the overview. This favours
the current method: over the 43 earlier worlds, 115 reads had calls in two families, all of them from the current method.

## Bars, per world, all six worlds

- **B1:** integrated correct >= current correct.
- **B2:** integrated errors <= current errors.

Correct and errors are summed over every read class. Abstain and unknowable count for neither. `summary.py bars` computes both
bars (`round8_bars`, tested in `test_summary.py`).

## Reading

- **PASS** (B1 and B2 hold in all six worlds). Claim: "on a fresh family set, the integrated O1 -> O2 -> O3 pipeline gives,
  in every world type, at least as many correct read outcomes as the current method and no more wrong ones".
- **FAIL** otherwise. No partial claim. The failing worlds and their margins are reported.

## Known before running

- **The 43 earlier worlds** (rounds 1-7, post hoc, `docs/O1O2O3_INTEGRATED_OVERVIEW_2026-10-16.md`): both bars held in 38.
  - B1 failed by 1, 2, 5 and 12 reads in round 3's and round 4's INTACT and CROSS, where neither system made an error.
  - B2 failed by 1 error in round 6's REMOVED, alongside 227 more correct reads.
  - A strict bar can therefore fail by a handful of reads.
- **Round 7's held-out set** passes these bars post hoc, although its pre-registered per-objective O2 bars failed.
- **The one-family rule on development worlds.** These are the six worlds where a new locus had been seeded twice: set 1's
  five worlds (NPIPA/NPIPB) and round 7's OUTSIDE (WASH/DDX11L).
  - It was never worse.
  - Set 1's OUTSIDE errors fell from 30 to 17, and set 1's CROSS gained 10 correct reads.
  - Nothing else changed.
- **BRANCH with alleles** had never been run before. On round 7's families (development):
  - current: 2,406 correct and 134 errors;
  - integrated: 2,888 correct and 9 errors;
  - arm I made 16 calls, all true.
- **The integrated pipeline's false absent-copy calls** over the 43 worlds were all one cluster (NPIPA.5's reads hosted on
  NPIPB, at delta). The fix does not change it.
- **Its errors on uncatalogued copies** (291) come from near-identical copies whose reads tie at mapq 0. This set has six such
  families.

## Reported beside the bars (not gating)

- the fates per read class, and the error sources (as in the overview);
- **O1:** copies kept and dropped (silent versus expressed), new copies (on a truth copy or not), and loci resolved by the
  one-family rule;
- **O3:** arm I's calls, true or false;
- the same per family.

## After the run

- An independent review by a fresh subagent: the code path, the numbers, and this document's reading.
- Register row 1306.
- Results appended below a line, without editing anything above it.
