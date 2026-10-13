# PREREG — population-fixed PSVs for O2 on the human Y (2026-10-12 plan; written 2026-10-13, before any run below)

**Why.** The YAG O2 test (`docs/PREREG_yag_o2_2026-10-12.md`) found O2 calibrated on the Y ampliconic genes but its
split-half concordance on real A119b testis reads was 0.544 (31/57, FAIL). The reading offered there (post hoc): the
testis donor's Y is not the reference man's (CHM13 Y = HG002), so a read's two halves can point to different REFERENCE
copies with no error by O2. The user asked for the fallback when genome and reads come from different men: keep only
PSVs that are FIXED in the population, and label each call by whether it rests on them
(`docs/superpowers/specs/2026-10-12-population-fixed-psvs-design.md`, plan `docs/superpowers/plans/2026-10-12-population-fixed-psvs.md`).

**Built before this prereg (no O2 run used it yet).**
- Population: the 48 HGSVC3 / Hallast et al. 2023 chrY assemblies (`/mnt/c/Users/jfris/Desktop/data/hgsvc_sigY/`, md5-verified).
- `bench/popfix/` -> `/mnt/linuxdisk/tmp/popfix/psv_population.tsv` from the A119b catalog's 110 chrY families with <= 100
  copies (`utilities psv-discover`, default semiglobal engine): 586,019 sites. On chrY: FIXED 6,495 (alignment 404 /
  k-mer 6,091), POLYMORPHIC 1,311, UNCOVERED 1,288, NOT_ASSESSED 27,381 (columns with a non-Y copy); every non-Y site
  NOT_ASSESSED. Where both layers give a verdict they agree on 194 / 195 columns.
- Rust (`src/family.rs`, `src/copy_assign_bin.rs`): `copy_assign --psv-population FILE` writes
  `<out>.population_tier.tsv` (per contested molecule: the run's final call, the pairwise read-star certificate re-run on
  columns where BOTH copies' genome positions are FIXED, tier PORTABLE / REFERENCE_ONLY / UNASSIGNED); `--psv-fixed-only`
  makes the fixed-only certificate the decision; `--dump-star` gains a `positions` column.

**Observed while building (controller, 2026-10-12), which this prereg must be able to separate.** In the YAG real-read
star dump (`/mnt/linuxdisk/tmp/yag_o2/yag_human/run.star_reads.tsv`) BPY2 (MCL1210) molecules carry up to 121 read-star
columns in dense runs (e.g. 55 columns in ~230 bp, copy 0 against copies 1+2 at nearly every base), while the catalog has
ONE PSV column for that family (copies 99.975% identical). Those columns are the read forced across non-homologous
sequence (a different exon / shifted frame), not PSVs. That is a SECOND explanation of the 0.544, competing with the
donor-versus-reference one. Restricting to FIXED PSVs removes both kinds of column, so arm (c) alone cannot tell them
apart; arm (c') is added for that.

## Arms and bars

Families: the YAG families of each catalog (`prep_yag.py`; A119b: MCL54 TSPY, MCL480 RBMY, MCL479/MCL377/MCL815 HSFY,
MCL475 VCY, MCL478 CDY, MCL635 PRY, MCL193/MCL2142 DAZ, MCL1210 BPY2). A column is "fixed for the pair (a, b)" iff both
candidates have a genome position there and both `(chrom, pos)` are FIXED in the table (the Rust rule).

**(a) Table sanity — descriptive, no bar.** Status x method counts per YAG family (`bench/popfix/sanity.py`); k-mer vs
alignment agreement (above). Expectation stated in the spec: palindrome-arm PSVs often POLYMORPHIC, unique-flank PSVs
mostly FIXED.

**(b) Simulation with truth, `--psv-fixed-only`.** The 10-06 chrY simulation exactly as YAG Arm A (`simY2` command:
`sim.bam`, roster `asg.fam.copies.tsv/.fa`, `--skip-families asg.skip.families.txt`, `--skip-poa-diagnostic
--region-threads 4 --dump-star`) plus `--psv-population psv_population.tsv --psv-fixed-only`; scored by
`bench/yag_o2/arm_a.py` with `YAG_FAMILIES` (YAG families of the roster). Reads come from the reference copies, so this
checks that the strict path keeps O2's calibration, not portability.
- B1: 0 wrong among assigned (upper 95% bound reported; UNDERPOWERED if < 30 assigned).
- B2: 0 of the NOT-distinguishable molecules assigned (Arm A's definition, full read-star columns).
- Reported: assigned count with and without the flag (the without-flag numbers are YAG Arm A's).

**(c) Real A119b testis, YAG families — the pre-specified test.** The YAG Arm B command (`A119b.t2t.bam`, `yag_human`
copies/regions, `--skip-poa-diagnostic --tied-prefilter --dump-star --threads 4`) plus `--psv-population` (no strict).
`POPULATION=psv_population.tsv bench/yag_o2/split_half.py`: the same certificate as Arm B on columns fixed for the pair.
- C1: split-half concordance >= 0.98 among molecules whose two halves both assign, with >= 30 such molecules (else
  UNDERPOWERED).
- Reported: PORTABLE / REFERENCE_ONLY / UNASSIGNED counts from `population_tier.tsv`, per family; molecules with a
  PORTABLE single-copy call; the fraction of contested molecules with >= 1 assigning half (the both-halves denominator
  is conditioned on the certificate — never read C1 without it).

**(c') Discriminating arm (added before any run).** The same split-half on columns where both candidates' positions are
CATALOG PSV SITES — present in `psv_population.tsv` with ANY status (`POPULATION_ANY=1`).
- C'1: split-half concordance on catalog-PSV columns, same bar and power rule as C1.
- Reported per family: the fraction of read-star columns (contested YAG molecules) whose positions are catalog PSV sites
  for the pair, and the fraction FIXED.
- Reading, fixed now: C1 PASS and C'1 PASS -> the 0.544 came from columns that are not catalog PSVs (alignment across
  non-homologous sequence), and the donor explanation is not needed; C1 PASS and C'1 FAIL -> disagreement sits at
  catalog PSVs that are not population-FIXED (the donor explanation); C1 FAIL -> the disagreement is in O2's evidence,
  not the donor; any UNDERPOWERED -> no verdict for that comparison.

**(d) No regression.** Without the flags, `assignments`, `quant`, `families`, `psv_cols`, `psv_copies` byte-identical to
the default baseline `/mnt/linuxdisk/tmp/psv_sun/run_default/` (gorilla fibroblast dev subset); `cargo test --release
--no-fail-fast` passes except the 3 known `tests/pipeline_driver_dispatch.rs` failures.

**(e) Implementation check (bar = exact).** In run (c), molecules whose every read-star column is fixed for every pair
of candidates, and whose final status was not changed by a later gate (origin rejection, tie outside the catalog): where
either certificate says `assigned`, both say `assigned` with the same copy. Any mismatch is a bug: stop and fix before
reading (b)-(c'). (Exact ties among non-assigned molecules may pick a different best copy — ruling R16 — and are not
compared.)

Outputs `/mnt/linuxdisk/tmp/popfix/runs/`; results and verdicts appended below, nothing above edited after a run.

## Results

(pending)
