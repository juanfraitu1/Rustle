# PREREG row 1749 (2026-11-15): a satellite / rDNA darklist BEFORE the families stage

**Question (user, 2026-11-15):** assign already skips satellite / rDNA families (`RUSTLE_EXCLUDE_BED`, row 1559). If the
same loci are left out before families are built, do the remaining families get simpler or more accurate?

**Change under test:** opt-in `RUSTLE_FAMILIES_EXCLUDE_BED=BED` on the `families` stage (`tools/gtf_exclude_bed.py`):
a locus (gene_id) with >= 0.5 of its exon bases inside the BED is dropped from the families input
(`PREFIX.fam.input.gtf`). Unset = the same command as before.

**Sample and arms:** A119b (seen sample; this is a safety and cost test, not a held-out gain test), driver `families`
stage at its defaults (exon-share edge rule, SD damping with the CHM13 SEDEF sidecar, famreps), input
`~/runs/human_A119b/a1.families.gtf`, workspace `r1749_fam_darklist`:
- `base`: switch unset.
- `excl`: `RUSTLE_FAMILIES_EXCLUDE_BED=~/data/denovo_eval/chm13v2.0_censat_v2.1.bed` (the BED assign uses; it holds
  `ct` transition rows as well as satellites and rDNA).

**Score:** `pairscore.py` defaults (pair precision / recall, default level) on each arm's `fam.copies.tsv`, with
`--gtf a1.families.gtf` in BOTH arms so the truth universe is the same; hit dumps for the lost-pair count.

**Bars, fixed before the run:**
1. Safe: truth pairs hit by `base` and not by `excl` <= 2, and recall(excl) >= recall(base) - 0.002.
2. Precision(excl) >= precision(base) - 0.002.
3. "More accurate": precision or recall of `excl` >= base + 0.010. Otherwise the statement is "no accuracy effect".
4. "Simpler / cheaper": families stage wall of `excl` <= 0.90 x base (same node, arms run one after the other).

**Verdict:** PASS = bars 1 and 2. Bars 3 and 4 only decide the wording. The switch stays opt-in whatever the result;
making it a default is the user's decision and would need a second sample.

**Also reported, no bar:** loci dropped, families and copies per arm, families of `base` with no dropped locus whose
member set changed in `excl` (how much the satellite loci shaped gene families), copies of `excl` inside the BED.
