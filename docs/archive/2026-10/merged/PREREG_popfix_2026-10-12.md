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

## Amendment 1 (2026-10-13, before any run)

- **Gate P0-fixed.** `split_half.py` under `POPULATION` re-implements the fixed-only certificate in Python. On the full
  read it must reproduce `population_tier.tsv`'s `status_fixed` (and `copy_fixed` when `assigned`) for >= 0.99 of tier
  rows with `status_fixed != NA`. Fail -> (c) and (c') are not read until the mismatch is explained.
- **The mask is per pair**, exactly the Rust rule above. For candidates a and b, a column counts iff both have a
  genome position there and both `(chrom, pos)` are in the set: FIXED for (c), any status for (c'). Each candidate's chrom
  comes from the run's `family_join.tsv`. The PSV-score tie key of a candidate counts columns where its own position is
  in the set.
- **Arms (c), (c') and (e) read contested rows only** (`contested` = 1). Uncontested rows are the aligner's MAPQ-60
  placement and are not certificate calls. Under `--psv-fixed-only` they keep the aligner's call (`--psv-fixed-only`
  help text).
- **Strict-mode semantics fixed in code before any run** (rulings R20, R24, R25): a family with no genome positions has
  no fixed column (-> Tied); under strict no sole-candidate promotion; the read-through credit is withdrawn when the
  fixed-only best differs from the full-column best.

## Results (2026-10-13; binary ba525453, scorers 91fc1635; outputs `/mnt/linuxdisk/tmp/popfix/runs/`)

**Gate P0-fixed: 9,670 / 9,670 = 1.0, PASS.** The Python fixed-only certificate reproduces Rust's `status_fixed` and
`copy_fixed` on every tier row.

**(e) Implementation check: VACUOUS.** 9 molecules are eligible (every column FIXED for every candidate, no later gate),
plus 2,560 with no read-star column. Neither certificate assigns any of them, so the exact-agreement bar was not
exercised. 0 mismatches.

**(a) Table sanity (descriptive).** Rows are attributed to the strictest column's family.

| Family | Genes | FIXED / POLYMORPHIC / UNCOVERED / NOT_ASSESSED |
|---|---|---|
| MCL479 | HSFY | 783 / 0 / 18 / 0 |
| MCL635 | PRY | 472 / 12 / 0 / 0 |
| MCL478 | CDY | 144 / 8 / 4 / 0 |
| MCL54 | TSPY | 230 / 276 / 368 / 0 |
| MCL480 | RBMY | 12 / 24 / 0 / 0 |
| MCL2142 | DAZ | 2 / 6 / 0 / 0 |
| MCL193 | DAZ | 0 / 4 / 0 / 0 |
| MCL1210 | BPY2 | 0 / 3 / 0 / 0 |
| MCL377 | HSFY + chr22 | 0 / 0 / 0 / 432 |
| MCL475 | VCY + chrX | 0 / 0 / 0 / 227 |
| MCL815 | HSFY | no sites |

The palindrome-array genes (DAZ, BPY2, RBMY) are mostly POLYMORPHIC, as the spec expected.

Share of the contested molecules' read-star columns at catalog PSV sites (any status / FIXED):

| Family | Any status | FIXED |
|---|---|---|
| BPY2 | 0.09% | 0% |
| MCL193 | 2.8% | 0% |
| PRY | 11.5% | 11.0% |
| MCL479 | 12.2% | 9.0% |
| RBMY | 21.6% | 7.3% |
| MCL377 | 69.7% | 0% |
| TSPY | 84.7% | 23.0% |
| CDY | 93.3% | 80.4% |
| MCL2142 | 99.4% | 32.4% |

**(b) Simulation, `--psv-fixed-only`.**
- B1 UNDERPOWERED: 0 assigned.
- **B2 PASS: 0 of 1,458 not-distinguishable molecules assigned.**
- Identical to the no-flag YAG Arm A: 0 assigned among the 2,815 contested chrY-roster rows either way.

**(c) Real A119b testis, FIXED columns.**
- **C1 UNDERPOWERED.** 0 molecules have both halves assigning. 2 of 9,661 contested molecules (0.02%) have at least one
  assigning half.
- Full-read fixed-only certificate: 2 assigned, 224 ambiguous, 9,435 tied.
- Tiers: **PORTABLE 0**, REFERENCE_ONLY 4 (3 RBMY, 1 VCY), UNASSIGNED 12,027 (2,361 of them orphans with no family
  candidate, which have no tier).
- `assignments.tsv` is byte-identical to the YAG Arm B run.

**(c') Real reads, all catalog-PSV-site columns.**
- **C'1 UNDERPOWERED.** 0 molecules have both halves assigning. 124 / 9,661 (1.3%) have at least one assigning half:
  RBMY 110, MCL2142 11, PRY 2, CDY 1.
- Full-read certificate on these columns: 125 assigned (RBMY 110, MCL2142 11, PRY 3, CDY 1), 292 ambiguous, 9,244 tied.
- Reference, same run on ALL read-star columns: YAG Arm B reproduced exactly. 57 molecules have both halves assigning,
  31 of them agree (0.544), and 386 have at least one assigning half.

**(d) No regression: PASS.**
- Without the flags, `assignments`, `quant`, `families`, `psv_cols` and `psv_copies` are byte-identical to
  `run_default` (5 / 5).
- `params.tsv` gains two rows (`psv_population`, `psv_fixed_only`), the same precedent as `psv_engine`.
- `cargo test --release --no-fail-fast` passes everything except the 3 known `pipeline_driver_dispatch` failures (lib
  1,044 passed).

**Verdicts.**
- P0-fixed PASS.
- (b): B1 UNDERPOWERED, B2 PASS.
- (c): C1 UNDERPOWERED.
- (c'): C'1 UNDERPOWERED. So there is **no registered verdict** on donor-versus-artifact.
- (d) PASS.
- (e) VACUOUS.

**Correction to Amendment 1 (wording only, no bar affected).** The tier and the strict decision apply to every
molecule that goes through the read-star certificate. Under the default AS-tied gate that includes a few uncontested
molecules: 9 in run (c), 1 of them assigned. The sentence saying uncontested rows keep the aligner's call under
`--psv-fixed-only` is wrong. The final review found that even on the `--no-as-tied-only` placement path molecules go
through the certificate and get tier rows, and a strict Assigned call takes precedence over the placement; the
`copy_assign` help text states the exact behaviour. The arms above read contested rows only, as registered.

**What (b) does and does not show (final review).** Strict mode assigned 0 molecules in (b), as did the no-flag run,
so B2 PASS could not have failed. It shows strict mode adds no wrong calls there; it does not exercise a strict
Assigned call. That path is tested end to end only on a committed fixture (`tests/copy_assign_popfix.rs` on `o2_origin_resolve`):
the 1 baseline Assigned call is reproduced under an all-FIXED table and none is made with no FIXED site. It is not
tested on real data.

**Reading (post hoc, not tested here).**
1. **On the YAGs, population-FIXED PSVs are far too sparse per read to assign anything to a single copy.**
   - 0 molecules are PORTABLE.
   - The fallback "assign at least something across donors" gives nothing at copy level on these genes.
   - That is the case for sub-project A: assign to copy GROUPS under a flag.
2. **The YAG split-half failure was carried entirely by columns the catalog's PSV discovery does not report.**
   - All 57 molecules whose two halves both assigned lose that property once only catalog-PSV sites count, whatever
     their population status.
   - For BPY2, 99.9% of read-star columns are not catalog PSV sites.
   - Such columns have three possible sources:
     - the read aligned across non-homologous sequence of two copies;
     - differences the default engine's significant-block mask drops;
     - PSVs between two non-reference copies in stretches copy 0 lacks, which discovery never reports because it is
       anchored on copy 0 (final review).

     None of them is population variation at PSV sites.
   - This favours the alignment-artifact explanation over donor-versus-reference mosaicism for the 0.544. It is
     descriptive only, because the registered comparison is UNDERPOWERED.
3. **On catalog-PSV columns alone, the full-read certificate assigns 125 contested YAG molecules, against 3 on all
   columns** (mostly RBMY).
   - Non-PSV read-star columns appear to make O2 abstain on reads that catalog PSVs alone would assign.
   - This is a candidate O2 change, to be tested with truth (Arm A simulation) before any claim. It is not adopted here.
