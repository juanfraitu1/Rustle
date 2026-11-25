# PREREG r1949: assigning reads to copies, as two questions, at today's defaults (2026-11-25)

Why: user, 2026-11-25: "when possible to assign do we do better than minimap, and in case of unassignables, can we
correctly assign it to everything where it is needed?", and the advisor's "sharing information across the family" read as
isoforms several members produce. Row 1796 answered a form of both on 2026-11-18 (exact `ms` band against the primary,
`famwise` statuses). Since then the default labels use an EM with the alignment scores inside it (row 1921) and the EM's
per-read shares can be written. This row re-measures on the same reads with today's products.

**Not unseen data.** HG002 has no unopened reads left: quarter 25-49 was row 1796's judging set, 0-24 its development
set. Both are reported; neither is a held-out test. Workspace `r1949_assign_two_questions`.

## Inputs
- Truth, records and kinds exactly as row 1796 (`bench/famwise_eval/records.py`; `r1796_famwise_origin/{hg_held,hg_seen}`:
  `sel.bam`, truth table, tie split). `one` = one origin copy (assignable); `several` = fits several copies (unassignable).
- Today's products: HG002 on CHM13 at the defaults of 2026-11-23 (`r1931_hg002_sequences/hg/g.multicopy.gtf`), labels
  remade on the whole BAM with `em_filter.py --label --weighted 1 --ties --shares` (nothing else rerun).
- A share row's locus is matched to the truth's places by contig, strand and overlap of the locus span; a candidate the
  truth never judged (also `outside`) gives the statement no verdict (counted, not scored), as in row 1796.

## Statements about a molecule (molecules near-tied on `ms`, as row 1796)
primary alone; exact band (W = 0); within 10; **EM top** = the locus of highest share (all loci within 0.001 of it);
**EM set** = every locus with share >= 0.10. A molecule without share rows (no rival in the band at assembled loci) is
given its primary for both EM statements. `equals` = the statement's places are exactly the origin set.

## Readings, fixed now (quarter 25-49; quarter 0-24 printed beside)
- **Q1, assignable reads (`one`):** BETTER than minimap2 if EM top equals the origin more often than the primary with
  disjoint Wilson 95% intervals; WORSE if less often with disjoint intervals; else SAME.
- **Q2, unassignable reads (`several`):** HOLDS if the EM set equals the origin set for >= 0.90; also printed: exact band,
  within 10.
- **Q1 + Q2 by one rule:** over `one` + `several` pooled, each statement's share of `equals`; the best is named. No bar.
- **Q3, labels on transcripts:** claims as row 1796 (multi-exon transcripts of `g.o2.gtf`, verdict TRUE / SHARED / PHANTOM
  from judged supporters). HOLDS if PHANTOM + SHARED is higher among `homology_only` + `shared` than among `expressed`
  with disjoint intervals AND `expressed` holds >= 0.90 of the TRUE claims.
- Printed, no bar: molecules with no verdict per statement; isoforms carried at several members = intron-chain claims
  whose judged supporters include a `several` molecule, by label.

Predictions (mine): Q1 SAME or WORSE (an EM moves ambiguous reads to the deeper copy, row 1554). Q2 does not hold
(EM set 0.80-0.88; within 10 highest on `several`). Pooled: exact band best. Q3 holds.

## RESULT (2026-11-25, job 4504735; `r1949_assign_two_questions/{hg_held,hg_seen}.tsv`)
Checks: the remade labels are identical to the cached run's; the molecule counts reproduce row 1796 (3,259 with exact
band 2,662 and primary 341; seen quarter 3,379 / 2,772 / 396).

Share of molecules whose statement equals the origin set (quarter 25-49; quarter 0-24 in brackets):
| statement | `one` (386) | `several` (2,873) | pooled |
|---|---|---|---|
| primary alone | 0.883 (0.900) | 0 | 0.105 |
| exact band | 0.785 (0.816) | 0.821 (0.821) | 0.817 |
| within 10 | 0.003 (0) | 0.961 (0.959) | 0.848 |
| EM top | 0.806 (0.834) | 0 | 0.095 |
| EM set (share >= 0.10) | 0.767 (0.800) | 0.813 (0.810) | 0.807 |

- **Q1: WORSE than minimap2.** EM top 311 of 386 (0.806, 0.763 to 0.842) against the primary 341 (0.883, 0.848 to 0.912),
  intervals disjoint; the same on the other quarter (0.834 against 0.900). Predicted (same or worse).
- **Q2: does NOT hold for the EM set** (0.813, bar 0.90). Within 10 `ms` reaches 0.961; the exact band 0.821. Predicted.
- **Pooled: within 10 is highest (0.848)**, not the exact band as I predicted; it is carried by `several` being 88% of the
  molecules and is wrong for every `one` molecule (a near-tied read with one origin always has a rival inside 10).
- **No single statement serves both kinds.** With the kind known, primary for `one` and within 10 for `several` would be
  right for 3,102 of 3,259 (0.952). The read's own scores do not give the kind: 21% of `one` molecules have an exact tie,
  14% of `several` molecules have a true copy 1-10 below the best.
- **Q3: HOLDS, weakly.** `expressed` claims 250 of 1,549 shared or phantom (0.161, 0.144 to 0.181), `homology_only` 9 of 24
  (0.375, 0.212 to 0.573), intervals disjoint; `expressed` holds 1,299 of 1,314 TRUE claims. But only 24 judged claims carry
  `homology_only` (15 of them TRUE) and no judged claim carries `shared`; the label is per locus. Row 1796's per-transcript
  statuses separated far better on the same reads (certified 0.014 shared or phantom, tied_only 0.84).
- Isoforms carried by a molecule that fits several members: 190 of 1,549 `expressed` claims (220 of 1,579), 6 of 24 `homology_only`.

Reading. Assignable reads: nothing here improves on the primary alignment; the EM, which helps per copy, harms per read.
Unassignable reads: the tie table's band (within 10) already gives every needed copy for 96%; the exact band and the EM set
miss about one in five. What is missing is a way to tell the two kinds apart, and a per-transcript status in the default
output (the `famwise` statuses do this; the `o2` labels do not).
