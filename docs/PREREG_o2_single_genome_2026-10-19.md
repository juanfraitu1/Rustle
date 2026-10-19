# Pre-registration: O2 in one genome -- the fewest-edit rule on CHM13-only worlds (round 23, 2026-10-19)

**User (2026-10-19):**
- "in O2 we only need gene copies, the multi-mapping reads they carry that cannot be resolved by minimap2 so they are given
  an arbitrary coin toss assignation of which one is the primary. In this case I just need to get a way for O2 to try to
  find what is the most likely real placement (or the real primary) or if there are five multimappers only two could be
  actually informative and the others an artifact, so that is why I would prefer that until we make O2 work we stay with
  only 1 genome, in this case CHM13".
- Then, after the proposal below: "yes please".

## The rule (code d63df788, `fastset.py ... --fewest-edits` = `o2graph.py run W intega 2 ... --fewest-edits`)

For a read tied within 10 AS (or by sequence) on two or more copies of its family:
- **O2's own call stands**, its structural checks included.
- **Where O2 abstained for want of its two-edit margin** (`thin`, `tied`, `locus_thin`, `locus_tied`), the copy with the
  fewest edits takes the read when it is the only one. Edits count every unaligned base (`veto_distance`).
- **Otherwise** the read keeps the set of its fewest-edit copies. The other placements explain it worse and are dropped.
- **Any other abstention** (a structural conflict, a call off the tied copies) keeps every tied copy.

Output: `o2.intege.m2.tsv` (call, reason, v17's call as `base_call`, `candidates`, `set`), beside an unchanged
`o2.intega.m2.tsv` (v17).

## One genome

The worlds simulate every read from its CHM13 copy with no donor alleles (`profiles/*_v4_ref.json`: v4 with
`het_per_base` = `hom_per_base` = 0). The reads keep the profile's errors, truncations, retained introns and adapters.

Only the single-genome world types are run:
- **INTACT:** every copy in the genome and the catalogue;
- **CROSS:** copies also express paralogs' isoform structures;
- **SILENT:** a third of the copies express nothing.

REMOVED, OUTSIDE and BRANCH involve a copy the reference or catalogue lacks, which is O3's question.

**Development (spent: round 22's families, A119b profile, seed 1201, the three worlds; `fast/ref_dev`):**
- Of 846 reads tied in their own family, the rule places **428 right, 0 wrong**, against v17's 109 / 0.
- Every set holds the true copy (846 / 846), and the mean set is 1.65 copies against 2.46 tied.
- minimap2's primary is wrong on 232 of the 846.

## Families (held out; drawn by the fixed rule)

**The pool:** round 22's unspent pool (98) less every family sharing a gene with round 22's draw, leaving 72 families with
213 copies (`r23seq/families23seq_all.json`). None of them was run by the sample-copy step on real reads.

**The draw:** the rule of rounds 17-22 with seed 23 (`seqfamilies.py draw`): **37 families, 90 copies**.
- **Copies per family:** {2: 32, 3: 3, 4: 1, 13: 1}.
- **Copy biotypes:** 43 protein-coding, 22 transcribed unprocessed, 15 transcribed processed, 5 unprocessed, 2 processed,
  2 lncRNA, 1 polymorphic.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families23.json` (set `r23`, sha1 0e46c2d3b8d8).

The families: AC1387761, AL5913791, AP0005524, AP0065871, CASP1, CKMT1B, CTRB2, CXorf51B, CYB5B, CYP4Z1, FAM153B, KRT16,
LRRC37B, MAGEA9B, MRPL36, MST1, MSTRG5154, NOS2P2, NPY4R, OR1D5, OR4C13, OR5M10, OR5T3, OTOAP1, PRDM7, PRSS40B, REXO1L1P (13
copies), RPL38, RPS27AP5, RPS4Y1, SETSIP, SMPD4, THOC3, TPTE2P2, TRIM77, ZNF33A, ZNF852 (`S19` prefix).

## Runs

- **Cells:** 2 profiles (A119b, public) x 3 seeds (1301, 1302, 1303), with no seed selection.
- **Worlds:** 3 per cell, so 18 in all.

```
F=/mnt/linuxdisk/tmp/o3_graph/families23.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout23
for S in 1301 1302 1303; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r23 $S --local $O/local $V17 --fewest-edits \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v4_ref.json --worlds INTACT,CROSS,SILENT
done; done
python3 bench/o3_graph/summary.py sets $O/sets.json $O/*_s13*/{INTACT,CROSS,SILENT} --off-family-ties 9 --seq-tie
```

`$V17` is round 22's flag set, unchanged.

## Bars (pooled over the 18 worlds; reads tied in their own family, off-family ties apart)

- **F1, precise:** the rule's wrong placements are at most 0.5% of its placements. Floor: 200 placements; below it,
  UNDERPOWERED.
- **F2, sets keep the truth:** the set holds the true copy for at least 99% of these reads. A placement counts as a set of
  one.
- **F3, more placed:** the rule places more of these reads than v17.
- **F4, no harm end to end:** fastset's verdict, now scored on the rule's table, is PASS in every cell (B1 and B2 in every
  world).

**Verdict:**
- **PASS:** F1 is powered, and F1 to F4 hold.
- **FAIL:** a powered bar fails.
- **UNDERPOWERED:** otherwise.

**Reported:**
- per cell and world: placements right / wrong, v17's, the sets' sizes against the tied copies, and minimap2's primary
  right / wrong;
- the reasons in `o2.intege`.

**Rules:** code frozen at d63df788. A bug is recorded as an amendment before any rerun.

## What a PASS would and would not claim

**A PASS would claim:** in one genome (reads from CHM13's own copies), O2 places most of the reads minimap2 leaves within
10 AS on their true copy at no more than 0.5% error. For the rest, it returns the smaller set of equally explained copies,
and that set holds the true one.

**It would not claim:**
- real samples, whose copies differ from CHM13's (O3: on HG002 a one-edit margin picked the wrong copy 46.5% of the time);
- copies the reference or catalogue lacks;
- that identical copies can be split. No rule can split them; the set says so.

## Results (2026-10-19, run as frozen; `fast/heldout23/run.log`, `sets.json`, `sets.txt`, each cell's `bars.txt`)

**PASS: F1 to F4 hold.**
- **F1:** 0 wrong of 2,574 placements (limit 0.5%, floor 200).
- **F2:** sets hold the true copy for 8,360 of 8,360 reads.
- **F3:** 2,574 placed, against v17's 959.
- **F4:** fastset's verdict is PASS in all 6 cells.

| cell | tied reads | rule: right / wrong | v17: right / wrong | sets holding the true copy | mean copies: tied -> set | minimap2 primary wrong |
|---|---|---|---|---|---|---|
| A119b s1301 | 1,398 | 421 / 0 | 153 / 0 | 1,398 | 2.01 -> 1.70 | 519 (37.1%) |
| A119b s1302 | 1,583 | 592 / 0 | 308 / 0 | 1,583 | 2.04 -> 1.63 | 598 (37.8%) |
| A119b s1303 | 1,477 | 404 / 0 | 95 / 0 | 1,477 | 2.01 -> 1.73 | 621 (42.0%) |
| public s1301 | 1,149 | 285 / 0 | 71 / 0 | 1,149 | 2.01 -> 1.75 | 402 (35.0%) |
| public s1302 | 1,429 | 526 / 0 | 235 / 0 | 1,429 | 2.00 -> 1.63 | 445 (31.1%) |
| public s1303 | 1,324 | 346 / 0 | 97 / 0 | 1,324 | 2.00 -> 1.74 | 480 (36.3%) |
| **pooled** | **8,360** | **2,574 / 0** | **959 / 0** | **8,360** | **2.01 -> 1.69** | **3,065 (36.7%)** |

**Reading:**
- On reads minimap2 leaves within 10 AS, its primary is wrong 36.7% of the time.
- The rule places 31% of them (2.7 times v17), all on their true copy.
- The other 69% keep a set of equally explained copies, and the set always holds the true copy. Most tied sets are pairs
  (mean 2.01), so a set of two means the two copies are identical over the read; no read placement can split them.

**Reasons on every tied row (all families, all reads):**
- placed: fewest_edits 1,847, O2's own assigned 962 and short_path 20 -- all 2,829 on the true copy, in the read's own
  family;
- kept as a set: fewest_edit_set 6,042;
- every tied copy kept: structure_only 58, unexplained 15.

**End to end:** the only integrated error is one read in two public_s1303 worlds (s0000683, FAM153B). minimap2 resolved it
by 10 AS or more onto the wrong copy (`as_resolved`), outside O2's domain, and v17 makes the same call.
