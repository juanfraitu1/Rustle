# PREREG — does O2 recover the right copy? Truth from the read's own unique anchor (2026-10-13)

**Why.** The user asked for examples that prove O2 works. Existing evidence:
- O2's calibration (it never assigns undecidable reads; YAG Arm A, chr16 simulation 157/157);
- haplotype phasing on the assembly's own animal (ASJ, T1 discordance .001 held-out).

Real-read ACCURACY has never been measured with a truth independent of O2's own evidence. Simulated reads come from
the reference copies, and real reads have had no truth. The user chose three truth sources, run in this order, each
pre-registered before it runs:
1. anchor trim (this document);
2. copy-specific junctions;
3. split-half on the reference animal itself.

Arms 2 and 3 are added as amendments before they run.

## Arm 1 — anchor trim (gorilla = dev; human = held-out)

**Truth construction** (`bench/o2_anchor/anchor_trim.py`, 7 tests; run before this prereg, counts only):
- **Read selection.** Per catalog family with >= 2 copies (gorilla_KB3781 default catalog, 326 families / 918 copies),
  take primary MAPQ-60 reads of the fibroblast BAM (`-F 2308 -q 60`). The reads come from the assembly's own animal.
  Each read's reference span must touch exactly one family; truth copy = the copy it overlaps most. At most 20 reads per
  copy, the lowest md5 of the name.
- **Trimming.** Each read, in its original orientation, is aligned to the genomic spans of its family's copies
  (`minimap2 -c -x splice:hq -uf`, secondaries kept). It is trimmed to the interval covered by both the truth copy's best
  hit and the best sibling's best hit. The trimmed read is kept iff that interval is >= 300 bp and >= 100 bp were
  removed. The truth is therefore the removed, copy-unique part of the read. O2 sees only the part shared with a sibling.
- **Feasibility (gorilla).** 10,395 reads picked, 4,048 trimmed, from 240 families.
- **Re-alignment.** The trimmed reads are re-aligned genome-wide with the original BAM's settings (`-ax splice:hq -uf
  --eqx -Y -N 50 -p 0.1 --secondary=yes`, GGO.splice.mmi).
- **Primary MAPQ after trimming.** 3,152 stay MAPQ 60 (the shared part still places them), 744 are MAPQ 1-59 and 152
  are MAPQ 0.

**O2 run.** `copy_assign` exactly as the pipeline's assign stage:
- inputs: `--regions`, `--families`, `--copies-fa`, plus `--skip-poa-diagnostic --region-threads 4 --tied-prefilter`
  and `--dump-star`;
- BAM: the original records in every catalog copy's ±50 kb window (what `copy_assign` reads) merged with all records
  of the trimmed reads (names `<read>|trim`);
- families: `--only-families` = every family a trimmed read's record touches, run in batches if needed (rows are per
  read x family, so batching does not change them).

The original full reads stay in the BAM; each trimmed read is an extra molecule.

**Scoring** (`bench/o2_anchor/score.py`; locus-level: a copy is right iff its catalog span overlaps the truth copy's
span):
- **T0 truth gate.** Among trimmed reads still MAPQ 60, the aligner's primary is on the truth copy for >= 0.99. If not,
  the truth construction is suspect: stop, nothing below is read.
- **A1 precision.** Over the contested trimmed reads (`contested` = 1 in the truth family's row), count every `assigned`
  row, in the truth family and in any other family (a foreign assigned row is wrong).
  - PASS iff wrong / assigned <= 0.01 with >= 30 assigned. UNDERPOWERED if < 30 assigned.
  - The 95% upper bound is reported.
- **A2 calibration.** Some contested trimmed reads are NOT distinguishable: the truth copy is not a candidate, or some
  rival differs from it at < 2 read-star columns the read covers (YAG Arm A's definition, O2's m >= 2 step). Of
  those, <= 0.05 assigned.
- **Reported, no bar:**
  - recall over distinguishable contested reads;
  - the aligner primary's accuracy on the same contested reads, and on those O2 assigns (O2 against a coin toss);
  - per family;
  - an examples file: every wrong call, and up to 20 correct calls with their columns.

**Held-out: human A119b** (testis, catalog human_A119b, CHM13; the donor is not the reference man). Same code and bars
(H0-H2), run only after the gorilla results are written below. A donor mismatch can lower recall there. Its effect on
precision is the question.

Outputs `/mnt/linuxdisk/tmp/o2_anchor/`; the merged BAM is on C: (`/mnt/c/Users/jfris/Desktop/data/o2_anchor/`).

## Results — gorilla (dev), 2026-10-13

Run `/mnt/linuxdisk/tmp/o2_anchor/gorilla/` (copy_assign 206 s, 4.8 GB; scorer committed before scoring, b/o2_anchor).

**T0 truth gate: 3,105 / 3,152 = 0.985, FAIL (bar 0.99).** Per the registration, A1 and A2 are not verdicts. The 47
trimmed reads whose shared part is still placed uniquely, but away from the truth copy, are molecules whose anchor and
shared part point to different loci:
- 28 sit on the very sibling they were trimmed against. These are mosaics: e.g. 4 MCL1 reads are anchored at copy 3,
  but their shared part matches copy 8. Possible causes: allelic gene conversion (the primary assembly holds one
  haplotype), or chimeric molecules.
- 16 sit within 100 kb of the truth copy, outside its catalog span: a nearby duplicate segment.
- 3 sit on no catalog copy.

The anchor truth therefore disagrees with the molecule's own shared sequence about 1.5% of the time. A 1% precision bar
cannot be judged against it.

**Population (descriptive, not verdicts).**
- Of 4,048 trimmed reads: 3,152 stay MAPQ 60, 789 are MAPQ < 60 but not exact AS ties (the aligner still has a best
  placement), and **107 become exact AS ties**, O2's subject.
- Of those 107:
  - 95 are not distinguishable (some rival is identical over the read's columns);
  - 12 are distinguishable;
  - the aligner's primary is right for 55 (0.51, a coin toss).
- Printed by the scorer but not read as verdicts (T0 FAIL): 1 molecule assigned (right), 0 of the 95 undecidable
  assigned.

## Amendment 1 (2026-10-13, after the gorilla dev run, before the human held-out run)

The gorilla run measured the truth's own disagreement: d = 1 - T0 rate (gorilla 0.015). Molecules whose shared part
points to another locus would be "wrong" for any method that follows the shared sequence, O2 included. For the held-out
human:
- **T0 is reported as the truth-noise rate d, not as a gate.** It stops the arm only if d > 0.05.
- **A1 PASS** iff wrong / assigned <= d + 0.01, with >= 30 assigned (else UNDERPOWERED).
- **Each wrong call is classified:** `sibling` when O2's copy is the sibling the read was trimmed against (the mosaic
  signature), else `other`. The `other` count is reported separately.
- **A2 is unchanged.**

The gorilla numbers above are dev. Only the human run tests this amended rule.

**Amendment 1b (feasibility, before the human run).** The human catalog is 2,150 families / 9,907 copies, the BAM is
96 GB. So the human arm takes the 400 families with the lowest md5(family_id) among families with >= 2 copies, and per
copy streams at most the first 5,000 primary MAPQ-60 records (by position) before the md5 sampling of 20
(`anchor_trim.py ... 20 400 5000`). Gorilla defaults (all families, no cap) are unchanged. The index is CHM13's splice
index; the BAM's minimap2 settings are identical to gorilla's.

**Amendment 1c (feasibility, before the human copy_assign run).** The families touched by all human trimmed reads
come to 779 families / 17.3 M records, too large. `copy_assign` writes rows only for AS-tied molecules (gorilla: rows for
107 of 4,048 trimmed reads, all of them the tied ones). So the human `--only-families` list is restricted to:
- the truth families of trimmed reads whose best secondary AS >= 0.98 x primary AS (a superset of the exact ties);
- plus every family those reads' records touch.

That is 365 reads, 129 families, 81 Mb of windows. The scored rows are the same as with all families.

## Results — human A119b (held-out)

Run `/mnt/linuxdisk/tmp/o2_anchor/human/` (8 family batches, 133-247 s each, <= 4.7 GB). Selection: 4,031 trimmed reads
from 298 families.

**Amendment 1d.** MCL3 (260 copies, an array family) is excluded. It exhausted 20 GB alone, and it holds no truth
read. 55 trimmed reads have a record on it, so its possible foreign rows are unobserved.

- **T0 truth noise:** d = 1 - 3,523 / 3,572 = **0.0137**, below 0.05, so the arm continues (Amendment 1).
- **Population:** 53 contested (exact AS ties); the aligner's primary is right for 23 (0.43).
- **A2 calibration: PASS.** 0 of the 42 undecidable molecules were assigned.
- **A1 (amended): UNDERPOWERED.** 4 assigned (bar: 30).
  - 1 right: MCL871. The read matches the truth copy at 17 columns and the sibling at 4; p 6e-48; the aligner's
    primary was wrong.
  - 3 wrong: three MCL1623 molecules, all of the `sibling` kind.
- Recall over distinguishable reads: 1 / 11.

**The 3 wrong calls are an O2 BUG, not evidence.** Each read matches both candidates at exactly 5 of 10 columns, an
exact tie. Yet `assigned` with margin 0.000 and p 2.5e-13.
- The certificate's rule is "Assigned iff every rival p < alpha / #rivals AND LLR > 0". Summing +lr / -lr column by
  column in the order `+++--+--+-` gives LLR = 1.78e-15 instead of 0. Floating-point residue passes "> 0".
- The p-value test alone cannot stop it: 5 alleles of either copy are too many errors for the other, both ways.
- In the shipped runs: gorilla_KB3781 has 7 of 140 contested assigned calls with margin 0.000; chimp_PTR has 2 of 46.
- `fixed_only_certificate` (popfix) sums the same way.

Fix candidate (needs its own prereg / register row): compute the LLR from integer counts, (k - m) x lr (+ the junction
term likewise), so an exact tie is exactly 0.

**Gorilla dev, for comparison (not verdicts, T0 gate failed):** 1 assigned, right: MCL22. The read matches truth copy 5
at all 28 columns, and the rivals at 8-25 with 3-20 mismatches. The aligner's primary was wrong. 0 of the 95
undecidable molecules were assigned.

**Arm 1 reading.**
- Every assignment that does not come from the rounding bug is right: 2 of 2, both where the aligner was wrong.
- Calibration holds: 0 of 137 undecidable molecules assigned across both species.
- Removing a read's unique anchor rarely leaves an exact tie that PSVs can resolve, so precision stays UNDERPOWERED.
  Ties: 107 in gorilla, 53 in human; distinguishable: 12 and 11.
- The test found a real certificate bug.
