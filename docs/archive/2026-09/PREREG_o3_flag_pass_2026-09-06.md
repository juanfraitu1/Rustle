# PREREG — the O3 flag pass: every family's rejected reads through the missing-copy detector and the admission prototype (2026-09-06)

Objective O3 = detect and FLAG reference-absent / unannotated copies. Instruments already measured: the
consistent-site detector (§6ff, excision: 7–100× the controls at a removed 98.5–98.7 % copy) and the admission
prototype (§6fc, row 706: safe, helps only the admitted loci's own reads). This pass runs both over every family
of a sweep and writes a flag table. `bench/o3_flag_pass.py`.

## Definitions (fixed before the run)
- Unit read: a primary record (`-F 2308`) with an aligned block inside a unit of the family (the scorer's rule).
- Rejected read: unit read with `origin_rejected = 1` and `n_candidates ≥ 1`; grouped by its best candidate Y
  (`catalog_copy_idx`). Orphan: `n_candidates = 0`.
- Detector (verbatim §6ff): the rejected reads of (family, Y) aligned to Y's locus padded by the longest read
  (`minimap2 -x splice -c --eqx -N 1`); a consistent site = a position mismatched by ≥ 3 reads and by ≥ 50 % of
  the reads covering it; covered kb = positions covered by ≥ 3 reads; rate = sites / covered kb. Control = Y's own
  assigned MAPQ-60 unit reads, same alignment, same statistic. At most 500 reads per group (first 500 by name).
- Flag rule (no new constant): under H0 the rejected reads' sites arise at the control rate; P(Poisson(rate_ctl ×
  kb_rej) ≥ sites_rej), with rate_ctl floored at 1 site over the control's covered kb when the control has 0
  sites; flag `missing_copy` iff p < α / n_pairs with α = 0.001 (O2's alpha) and n_pairs = the pairs tested in
  the sweep. The §6ff 5/kb line is reported beside it, never used.
- Class of a rejected group: `divergent` if the median rejected read's mismatch bases exceed its unaligned
  bases against Y, else `structural` (an exon or extent Y lacks).
- Admission (verbatim §6fc/§6fb): rejected reads whose primary lies outside every unit of the family, plus the
  orphans, clustered by primary placement (≥ 3 reads within 5 kb); each cluster is a candidate locus, annotated
  with the RefSeq gene there (if any) and with the unit of another catalog family covering it (if any):
  `unannotated` / `annotated_no_unit` / `other_family`.

## Predictions
- P1 (reproduction, the gate): on `adj/excise/{npip_x7,npip_x11,npip_x2,npip_x13,znf_x3}` the pass reports the
  §6ff rates 47.7 / 2.9 / 17.7 / 14.0 / 117.6 sites per kb at the excised copy's best Y (± 0.1) and flags x7,
  x2, x13 and znf_x3 as `missing_copy`; whether x11 (2.9/kb) is flagged under the Poisson rule is reported.
- P2 (specificity): on `sweep_v13` (76 families, intact) the control rate is below 1 site/kb for ≥ 90 % of the
  (family, Y) pairs (heterozygosity + error), and the excised-copy positives' rates exceed every intact NPIP
  pair's rate except pairs that the pass classes `divergent` with ≥ 20 reads.
- P3 (what the intact catalog flags): the flagged pairs are ≤ 25 % of the pairs tested on the 3 contigs; the
  `structural` class is not flagged by the detector (its rejections carry unaligned bases, not sites); the
  admission candidates include the §6fb annotation-gap loci still unrepresented after L1 (ZNF600-region 419
  reads, PDXDC1/NTAN1 340 reads) unless L1's dropped members absorbed them (reported either way).
- P4 (genome-wide): the pass over `sweep_gw_v2` (754 families) completes; the flagged pairs and candidate loci
  are tabulated with their classes; no number is chosen after looking.
Fail rules: P1 rate mismatch > 0.1/kb ⟹ the pass does not implement the detector, fix before anything. P2's
control rate ≥ 1/kb in > 10 % of pairs ⟹ the control is not a null (report which pairs). P3 > 25 % flagged ⟹
the flag is not a rare event on an intact catalog and the rule is reported as such, not tightened after the
fact.
