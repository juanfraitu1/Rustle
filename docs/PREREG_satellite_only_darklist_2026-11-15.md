# PREREG row 1750 (2026-11-15): a satellite-only darklist, for families and for assign

**Why:** row 1749: the T2T censat BED is not a satellite list (its `censat` rows hold the TSPY / RBMY composite arrays,
its `ct` rows 211 Mb of transition sequence with NBPF, POTE, GOLGA), so as a darklist it removed 753 true pairs, and
the families assign skips with it (`RUSTLE_EXCLUDE_BED`, row 1559) hold 682 of 1,765 true pairs.

**List under test:** `annot_dl/censat.satellite_only.bed` = the censat v2.1 rows whose class is `hsat1A`, `hsat1B`,
`hsat2`, `hsat3`, `hor`, `dhor`, `mon`, `bsat`, `gsat` or `rDNA` (1,005 rows, 231.6 Mb); every `ct` and `censat` row is left out.
Fixed before any run; no row-level tuning afterwards.

**Sample:** A119b (seen; a safety and cost test). Workspace `r1750_satonly`.

**Bars:**
- A. Families switch (`RUSTLE_FAMILIES_EXCLUDE_BED` = the list; base arm = row 1749's `base`, same driver command):
  truth pairs lost <= 2, recall >= base - 0.002, precision >= base - 0.002 (`pairscore.py` defaults,
  `--gtf a1.families.gtf` in both arms).
- B. List content: the families of the base arm that the list names (`annotated_families.py --majority-in`), scored
  alone, hold <= 2 true pairs.
- C. Hotspot: of the base-arm families with most copies in chr13:13,000,000-16,000,000, >= 0.95 are named by the list.
- D. Assign cost: `assign` through the driver on the row 1690 base (`r1662_edgemiss/asg_def`, binaries `audit/g2/r76/bin`,
  as `r1690_softmask/asg_bed`: 346 s with the full BED) with `RUSTLE_EXCLUDE_BED` = the list finishes in <= 692 s
  (2 x); killed at 40 min = FAIL. And on families skipped by neither list, assigned rows within 1% of `asg_bed`.

**Verdict:** PASS = A, B, C and D. Each bar is reported on its own. Replacing the censat BED by this list in the
adopted assign route, or making the families switch a default, is the user's decision.

**Reported, no bar:** stage walls, loci dropped, families named, assigned rows in the families only the full BED
skips (what the full BED was costing), "more accurate" wording as in row 1749 (>= +0.010 on either axis).
