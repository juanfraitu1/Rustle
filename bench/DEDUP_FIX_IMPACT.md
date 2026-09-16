# Dedup Fix Impact Measurement

Measures the effect of the Task 1 `copy_assign` primary-read dedup fix (commit `b7bcd63e`) against the
pre-fix binary (base commit `d1cbc468`) and the `RUSTLE_LEGACY_PLACEMENT_DEDUP=1` escape hatch, on real
data: an O2 `--families` run, a de novo run, and a human chr20 `--gtf` run.

## The bug

`src/bin/copy_assign.rs` ~2403-2410 (pre-fix): the primary-read dedup meant to drop a read returned by two
adjacent windows keyed on placement `(chrom, ref_start, ref_end, introns)` and ran **within each window**, so
distinct molecules that happen to share a placement collapse to one. On chr20 (a single window): 25,341
primaries -> 11,755 survive (53.6% of real molecules dropped; only 474 of the dropped set are same-UMI
duplicates). The code's own comment stated the intended key was per-molecule.

## The fix (commit `b7bcd63e`)

Primary reads are now deduplicated **across windows only** (the intended behavior: drop a read genuinely
re-emitted by an adjacent window, never a distinct same-placement molecule inside one window).
`RUSTLE_LEGACY_PLACEMENT_DEDUP=1` restores the old within-window placement dedup byte-for-byte, as an escape
hatch.

## Binaries

| Arm | Binary | Commit |
|---|---|---|
| base | `/mnt/linuxdisk/home/juanfraitu/dedupfix/copy_assign.base` | `d1cbc468` (pre-fix) |
| fixed | `/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign` | `b7bcd63e` (this worktree HEAD), default |
| legacy | same binary as fixed | `b7bcd63e`, run with `RUSTLE_LEGACY_PLACEMENT_DEDUP=1` |

Build:
```
CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target cargo build --release --bin copy_assign
```
Exit 0.

## Commands

`bench/dedup_fix_impact.sh` (restructured to take optional `SUBSTRATE [ARM]` args so each `(substrate, arm)`
job is one bash call under a 10-minute cap; no-arg form still runs everything). Nine jobs run serially,
foreground, one at a time:

```
bash bench/dedup_fix_impact.sh o2_families base
bash bench/dedup_fix_impact.sh o2_families fixed
bash bench/dedup_fix_impact.sh o2_families legacy
bash bench/dedup_fix_impact.sh denovo base      # 20m52s
bash bench/dedup_fix_impact.sh denovo fixed     # 17m43s
bash bench/dedup_fix_impact.sh denovo legacy    # 21m51s
bash bench/dedup_fix_impact.sh chr20_gtf base   # 11s
bash bench/dedup_fix_impact.sh chr20_gtf fixed  # 9s
bash bench/dedup_fix_impact.sh chr20_gtf legacy # 8s
```
All nine exit 0.

Underlying per-substrate commands (from the script):
- `o2_families`: `--bam npip3.bam --fasta GGO.fasta --families batch.copies.tsv --copies-fa cat.copies.fa --regions regions.txt --dump-psv`
- `denovo`: `--bam npip3.bam --fasta GGO.fasta --regions regions.txt`
- `chr20_gtf`: `--gtf --bam chr20.bam --fasta chr20.fa --region chr20:1-66210255`

## Byte comparison

`cmp -s` on every `run.*` output file except `run.stdout`/`run.stderr`.

| Substrate | base vs legacy | fixed vs legacy |
|---|---|---|
| o2_families | **same** on all 9 files (assignments, famcn_readonly, families, family_join, params, psv_cols, psv_copies, psv_reads, quant) | **same** on all 4 pre-registered files: assignments, quant, families, family_join |
| denovo | **same** on all 5 files (assignments, famcn_readonly, families, params, quant) | not required to match; see deltas below (no family_join / psv_* files in this arm — no `--families`/`--dump-psv`) |
| chr20_gtf | **same** on all 6 files (assignments, famcn_readonly, families, gtf, params, quant) | `run.gtf` **DIFFERS** (the measured effect); assignments/famcn_readonly/families/params/quant same |

No `run.params.tsv` difference occurred anywhere (byte-identical in every base-vs-legacy pair), so there is no
Task 1 defect: **base and legacy agree exactly on every substrate and every output file.**

## O2 prediction outcome (HARD STOP check)

Pre-registered prediction: `run.assignments.tsv`, `run.quant.tsv`, `run.families.tsv`, `run.family_join.tsv`
identical between `o2_families/fixed` and `o2_families/legacy`. **CONFIRMED** — all four files byte-identical
(`cmp -s` exit 0 on each). No stop condition triggered; de novo and chr20 arms proceeded.

## De novo deltas

Primary-read counts (from `run.stderr`, `[detect_and_assign] N primary -> ...`, one line per swept region):

| Arm | Region 1 primaries | Region 2 primaries |
|---|---|---|
| legacy | 90,913 | 143,666 |
| fixed | 246,800 | 438,595 |

Direction matches the bug: legacy over-drops distinct same-placement molecules within a window, so fixed
retains far more primary reads feeding family/transcript construction.

Families and copies:

| Arm | families rows | quant rows (copies) |
|---|---|---|
| legacy | 8 | 16 |
| fixed | 7 | 14 |

Assignment status counts (`run.assignments.tsv`, column `status`, AS-tied-multimapper rows only — O2's default
scope):

| Status | legacy | fixed |
|---|---|---|
| ambiguous | 465 | 457 |
| tied | 78 | 51 |
| assigned | 16 | 3 |
| **total** | 559 | 511 |

Families differing between arms, matched by copy coordinates (`run.quant.tsv`, `copy_chrom:copy_start-copy_end`;
family IDs are not stable across runs — [[project_o1_holdout_mcl166]] applies):

- 6 of legacy's 8 families have an exact coordinate match in fixed (legacy CAFAM0-4, CAFAM6 == fixed CAFAM1-6
  respectively, byte-identical copy spans and `n_reads_hard`/`anchored_reads`/`n_reads_soft`).
- **legacy CAFAM5** (copy0 `NC_073242.2:34634146-35184412`, copy1 `35469209-35541577`, n_reads=58,
  assigned_j=13) and **legacy CAFAM7** (copy0 `35115121-35117262`, copy1 `35472289-35474430`, n_reads=31, all
  PSV/junction stats 0) are both **absent** from fixed by exact coordinate match.
- In their place, **fixed CAFAM0** has copy0 `35112041-35184412`, copy1 `35469209-35541577` (n_reads=41,
  assigned_j=0) — copy1 exactly matches legacy CAFAM5's copy1; copy0 shares CAFAM5's end coordinate
  (35184412) and sits close to CAFAM7's copy0 start (35115121 vs 35112041, ~3 kb).
- Corroborating evidence from stderr's `SharedAcrossFamilies` warnings: legacy flags CAFAM5 vs CAFAM7 as
  sharing genomic sequence at both copy pairs (`recip=0.00` and `recip=0.03`); that warning is **absent** in
  fixed's output — the two overlapping legacy families collapse into fixed's single CAFAM0. (Both arms
  still separately flag the CAFAM1/CAFAM2-legacy == CAFAM2/CAFAM3-fixed pair as `SharedAcrossFamilies`,
  unaffected by the fix.)

Interpretation offered as observation, not a causal claim beyond what the coordinates show: the pre-fix
within-window dedup fragmented one locus into two overlapping low-support "families"; the fix's retained
reads let that locus resolve as one family.

## chr20 `--gtf` anchor outcome

```
cd /mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
gffcompare -r chr20_ref.gtf -o $O/chr20_gtf/$arm/gffc $O/chr20_gtf/$arm/run.gtf
```

| Metric | fixed (measured) | fixed (expected, spec fidelity anchor A1) | legacy (measured) | legacy (expected) |
|---|---|---|---|---|
| Query mRNAs | 1022 | 1022 | 976 | 976 |
| Matching transcripts | 350 | 350 | 347 | 347 |
| Transcript level (sn/pr) | 7.7 / 34.2 | 7.7 / 34.2 | 7.6 / 35.6 | 7.6 / 35.6 |
| Intron chain level (sn/pr) | 8.1 / 42.9 | 8.1 / 42.9 | 8.0 / 44.6 | 8.0 / 44.6 |

**All six anchor numbers match exactly** for both arms. Query loci: fixed 468, legacy 456 (812 vs 773
multi-exon transcripts) — not part of the pre-registered anchor set but consistent with fixed retaining more
transcript-supporting reads.

## Default decision pending user review before merge.
