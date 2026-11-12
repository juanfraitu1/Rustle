# PREREG O2 CROSS-FAMILY CANDIDATES FOR SOLE-CANDIDATE READS (register 1657, 2026-11-12)

Why: r1656 found O2's wrong near-tie assignments on HG002 are sole-candidate rows; the follow-up found the truth copy
is in another family on 58 of 86. User (11-12): do the cross-family row.

## Step 0 (spent r1653 sets, job 4467169, `r1657_xfam/{hgpat,hg13}.step0.txt`): the demote-only guard STOPS
Blind rule tried: demote a sole-candidate row when one of the molecule's tied placements (`is_tied_placement`,
genome-wide, non-supplementary) lies on a unit of another family. Decidable sole rows, right / wrong:
| | tied placement on another family's unit | own family only |
|---|---|---|
| HG002 own haplotype | 16 / 12 | 3 / 3 |
| HG002 on CHM13 | 13 / 71 | 0 / 0 |
The criterion fires on 112 of 118 rows, right and wrong alike: a demote-only guard is `--no-sole-candidate` (86 / 32
removed) minus six rows. STOP, not built. It also measures the mechanism r1656 inferred: nearly every sole candidate
has a tie partner on another family's unit, which the per-family table never scores against it.
Limit: best AS / ms taken genome-wide here, per loaded region in the binary.

## The registered arm: decide across families instead of demoting (existing flag, no code)
Arm U = `copy_assign --tie-ms 10 --union-certificate`: one certificate over the union of the copies of every family a
tied placement touches (register 787's rule in full). Differs from closed rows: r1481 scored U on the A119b default
gate by assigned COUNT (59 < 63), never by truth and never under `--tie-ms`; `RUSTLE_XFAM_RECONCILE=abstain` only
demotes molecules assigned in two families (most wrong sole rows have an `ambiguous` partner row).
Comparator F = base with sole-candidate rows counted as tied (= `--no-sole-candidate`; derived from BASE, no run).

Sets (fresh): HG002 on its own paternal haplotype and HG002 on CHM13, the next 120 panel families (panel order, regions
<= 2 Mb) used by neither the r1651/r1644 panels nor the r1653 draws. BASE = r1653's command; U = BASE +
`--union-certificate`. Truth: `haplotype_truth.py` on BASE, `strict_truth.py` slack 0 on each arm, decidable assigned
rows.

Bars:
0. Power: >= 20 decidable BASE sole rows pooled, else NO VERDICT. Validity: U exits 0 and writes union rows.
1. On each set U has at least as many right rows as F.
2. On each set U's wrong rate <= F's wrong rate + 0.005.
3. Pooled, U right - F right >= 10 (else U buys nothing over the flag).
4. Bars 1-2 hold pooled with the family holding most BASE sole rows removed.
5. Cost: U wall <= 2 x BASE wall on each set.
Readout: fate under U of BASE sole rows (same copy / other family / not assigned, right / wrong); U vs BASE.
Verdict: all pass = PASS: recommend `--union-certificate` with `--tie-ms` over the flag; default stays the user's call.
Any of 1-4 fails = FAIL: `--no-sole-candidate` is the answer for this class. Bar 5 alone failing = PASS with cost note.
Outputs: `/anvme/workspace/v138da13-rustle/r1657_xfam/`.

## RESULT r1657 (2026-11-12, job 4467177, 76 min; `r1657_xfam/result.txt`)
Decidable assigned rows, right / wrong (wrong rate):
| | BASE | F (sole rows tied) | U (`--union-certificate`) | bar 1 | bar 2 |
|---|---|---|---|---|---|
| HG002 own haplotype | 160 / 26 (0.140) | 140 / 16 (0.103) | 170 / 14 (0.076) | PASS | PASS |
| HG002 on CHM13 | 81 / 75 (0.481) | 81 / 66 (0.449) | 108 / 105 (0.493) | PASS | FAIL |
- Bar 0 PASS (39 decidable BASE sole rows; thin). Bar 3 PASS (U +57 right over F pooled).
- Bar 4 FAIL (minus MCL499: F 216 / 80, 0.270; U 267 / 110, 0.292). Bar 5 PASS (U wall 1.01 x and 1.03 x BASE).
**FAIL as registered** (bar 2 on CHM13, bar 4).
Readout. BASE sole rows under U: own haplotype 20 right all kept, 10 wrong -> 9 still wrong, 1 moved to the right
family; CHM13 9 wrong -> 8 not assigned, 1 still wrong. So U does what it was asked on the sole rows of CHM13 and keeps
every right one, but does not fix the own-haplotype wrong ones, and it ADDS assignments: on CHM13 MCL720 goes from no
decidable assigned row to 17 right / 31 wrong, which is the bar-2 failure.
Side observation, not a bar: on this third CHM13 draw the wrong rate is 0.449 even without sole rows (MCL598 0 / 34,
MCL1339 0 / 20): in thinner families the cross-individual error is not the sole-candidate class.
Conclusion: for the sole-candidate class `--no-sole-candidate` stays the recommendation; `--union-certificate` is
better on the same individual (0.076 vs 0.103, +30 right) and worse across individuals; stays opt-in.

## Follow-up: why MCL598, MCL1339, MCL720 go wrong (2026-11-12, post hoc; `bench/o2_contested/truth_audit.py`, `r1657_xfam/fix/`)
- MCL598 (0 / 34): TRUTH ARTIFACT. The reads are 420 bp, spliced over 3 introns, NM 0 on the chr12 copy (copy 0, the
  pick and the primary) and NM 0-1 at the chr12 locus on both HG002 haplotypes. Copy 0's own sequence differs from
  HG002 there by 5 edits in 1,418 bp, over `strict_truth`'s near-perfect gate (max(2, 0.003 x aligned) = 4.3), so it
  gets no locus; the 424 bp retrocopy (copy 3) reaches the same locus within its +2 slack and becomes the only truth.
- MCL720 under U (17 / 31): 26 of the 31 are the same artifact (the read sits NM 0 on the pick's own best haplotype
  hit; the pick fails the gate at NM 7 of 1,859; the 2-exon neighbour MCL991 overlaps the read's end and is named).
- MCL1339 (2 / 20): REAL. 12 reads are NM 0 on HG002's copy-0 locus; on CHM13 copy 1 is the primary at AS 1659 /
  ms 1756 / NM 4 against copy 0 at AS 1658 / ms 1754 / NM 3: HG002's allele is 3-4 edits from BOTH reference copies,
  and the pick (sole candidate) goes to the wrong one. 8 more name the nested family MCL1335 (not examined).
The gate was written for reads and catalog from one individual (r1510); across individuals a copy that differs from
the truth haplotype by > 0.3% loses its locus and a shorter sibling inherits it.
Audit of every HG002 set (wrong rows whose pick's own best haplotype hit is the read's best locus):
| set | right / wrong as scored | artifacts | corrected |
|---|---|---|---|
| r1651 CHM13 panel | 170 / 175 (0.507) | 81 (MCL189 58, MCL601 23) | 251 / 94 (0.272) |
| r1653 CHM13 held-out | 1,010 / 141 (0.123) | 27 (MCL338) | 1,037 / 114 (0.099) |
| r1657 CHM13 draw, BASE | 81 / 75 (0.481) | 35 | 116 / 40 (0.256) |
| r1653 own haplotype | 972 / 63 (0.061) | 0 | unchanged |
| r1657 own haplotype, BASE | 160 / 26 (0.140) | 9 (MCL101) | 169 / 17 (0.091) |
Verdicts re-read on the corrected dumps (post hoc, registered verdicts stand):
- r1656 still PASS: CHM13 sole rows 37 / 47 (0.560) vs others 1,000 / 67 (0.063); own haplotype unchanged.
- r1657 still FAIL, on the other set: CHM13 U 169 / 44 (0.207) vs F 116 / 31 (0.211) passes; own haplotype U 171 / 13
  (0.071) vs F 148 / 8 (0.051) fails bar 2. U is +76 right pooled at a similar wrong rate.
- r1651's "O2 170 right / 175 wrong" becomes 251 / 94: O2 across individuals is wrong on about a quarter of decidable
  rows, not half; sole candidates are 63 of the 94.
Limits: "pick on the read's best locus = right" is my correction rule, not a registered truth; KB3781 (same
individual) also shows 7 and 14 such rows, not looked into.
