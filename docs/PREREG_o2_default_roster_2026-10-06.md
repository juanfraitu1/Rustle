# PREREG: O2 with read-level truth on the default families roster, at NPIP, TBC1D3 and the Y-ampliconic genes (2026-10-06)

Written and committed before any roster, simulation or assignment product of this protocol exists. The only thing run before this registration is the report script
(`bench/o2_default_roster/report.py`) on the RECORDED legacy chr16 simulation (`o2sim/human`, 2026-09-25): it reproduces that run's tallies exactly (OWN 163 correct / 0 wrong /
1,084 abstain / 16 other; PRIMARY 112 / 205 / 3 conflict; ANY 60 / 432 / 156; union 0 assigned; n = 1,263 MAPQ-0 reads) and passes the cross-check of `figures/_o2.py`
against the ALL rows of `score.py reads`. That is a check of the instrument on recorded data, not a reading of any arm below.

## Question

Real reads carry no copy of origin, so O2's accuracy comes from simulation (`docs/PREREG_o2_read_truth_2026-09-23.md`). On the legacy chr16 catalog the certificate was exact within the
true family (157/157 correct, then 163/163) but the readings a consumer can use were not (PRIMARY accuracy .35, ANY .09 among assigned), and the table `assign` receives in the shipped
pipeline, the DEFAULT families roster (`mcl_families --emit-units`: one copy per member locus), has never been scored with read truth (`docs/PREREG_families_copy_table_2026-09-25.md`
holds engineering checks only). At the three families of this evaluation, on that roster, how many contested reads does O2 assign, and are its assignments right?

## Design (the registered protocol, on a different roster)

For each contig C in chr16 (target NPIP), chr17 (TBC1D3), chrY (the Y-ampliconic genes, "Y"), human A119b, HEAD binaries (`rustle_target_m2/release`):

1. **Roster R_C**: the default pipeline's families on C: `copy_assign --assemble-only --region C` with the driver's `assemble` flags (strict junctions, shipped polish, `--bridge-regroup f1v2`,
   secondaries >= 0.98 of the best AS from the stored genome-wide table), then the driver's `families` stage (`--min-cov-shorter 0.70`, most-reads representative). Its copy table
   (`PREFIX.fam.copies.tsv/.fa/.regions`: one copy per member locus of every multi-copy family) is what `assign` is handed. The products are compared byte for byte with the stored e163d955
   products where those exist (reported).
2. **Reads**: `bench/sim.py copies` (seed 20260925, the figures' seed): every copy of every family of >= 2 copies and >= 300 bp spliced sequence gets min(100, max(10, real read count))
   HiFi-model reads (substitution .001, indel .0003, up to 10% trimmed from a random end, +-0-30 bp end jitter), named `family|copy|i` (the truth), mapped against the WHOLE CHM13 genome with the
   shipped command (`minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, prebuilt splice index), so every paralogue competes as in production. The read model, depth rule and mapping are not changed after any number is seen.
3. **Catalog**: R_C minus the copies that received no primary alignment over their span (`copy_assign --families` refuses a copy with no reads; `figures/_o2.py derive_catalog`).
4. **Assignment**: the driver's `assign` stage on the simulated BAM (regions = the families' hulls, shipped defaults; arm O2) and `copy_assign --union-certificate` with the same inputs (arm U2, beside).
5. **Scoring**: `bench/score.py reads --per-read` (the registered scorer) joined with the BAM by `figures/_o2.py per_read` (its tallies must equal the scorer's ALL rows; a failure stops the report).
   Readings as registered: **OWN** (the true family's row: the certificate's intrinsic accuracy), **PRIMARY** (rows with `primary_local = 1`), **ANY** (any assigned row; two loci = conflict), and UNION (arm U2).
   Outcomes per read: correct / wrong / conflict / abstain / other (no row). Assigned = correct + wrong + conflict; accuracy = correct / assigned; coverage = assigned / MAPQ-0 reads.
   Population: the **MAPQ-0 reads** (the tied set O2 exists for; MAPQ > 0 is the aligner's). Identity band of a source copy = `figlib.identity_band` of the roster's own `max_family_identity`
   (identical, 99.5-100, 99-99.5, 98-99, < 98%).

## Target sets (defined now; a roster copy, same chromosome and strand, whose exons overlap >= 1 bp)

- **NPIP**: the exon union of any of the 25 CAT/Liftoff NPIP copies (`copy_recovery_tools_cat/ann`: `copies.hsa.tsv`, `truth.hsa.gtf`).
- **TBC1D3**: the exon union of any of the 16 CAT/Liftoff TBC1D3 copies (same files).
- **Y**: the body of a RefSeq chrY gene (`families_gw/species/human/genes_only.gff`) named DAZ*, RBMY*, TSPY*, BPY2*, HSFY*, VCY*, CDY*, PRY* or XKRY*.
Every table is given for the target's MAPQ-0 reads and for all MAPQ-0 reads of the roster (context). A source copy that overlaps no target truth is not in the target.

## Rules (per target, MAPQ-0 reads of the target's copies; DEV, first measurement; nothing is tuned)

- **B1 (certificate, OWN):** 0 wrong and 0 conflict among OWN-assigned reads, with >= 30 OWN-assigned reads: PASS. Any wrong or conflict: FAIL. Fewer than 30 assigned: UNDERPOWERED (no verdict; the count and,
  with 0 wrong, the rule-of-three upper bound 3 / n on the error rate are reported).
- **B2 (the readings a consumer can use, registered bar of `PREREG_o2_read_truth`):** for PRIMARY and for ANY separately: accuracy among assigned >= 0.95 overall, and >= 0.90 in each of the bands 98-99% and < 98%
  that holds >= 10 assigned reads, with >= 30 assigned reads overall: PASS. Otherwise FAIL, naming the part that fails. Fewer than 30 assigned: UNDERPOWERED.
- **B3 (union arm, beside, no bar):** the B1 test on the UNION reading.
- **No rescue:** a FAIL or UNDERPOWERED outcome is the finding. No read is added, no depth changed, no family dropped. A deeper design would be a new registration.

## Predicted before looking

- P1: OWN has 0 wrong wherever it assigns (157/157 and 163/163 on the legacy roster); NPIP and TBC1D3 reach >= 30 assigned.
- P2: the consumer readings (PRIMARY, ANY) stay below 0.95 at NPIP and TBC1D3 (legacy roster: .35 and .09); one family per NPIP gene array removes some foreign-family claims, not all of them. No prediction at Y.
- P3: Y is UNDERPOWERED on B1 (real chrY reads: 12 of 3,641 contested molecules assigned; the arms are 99.9% identical).
- P4: the union arm assigns almost nothing (0 of 1,263 on the legacy roster).

## Checks

GS1: the simulated BAM holds exactly one primary-or-unmapped record per simulated read (`figures/_o2.py`). GS2: the per-read tallies equal the ALL rows of `score.py reads` (OWN, ANY, UNION). GS3: the roster's byte comparison with the stored products is reported.
The scorer, the report and the runner are `bench/score.py`, `bench/o2_default_roster/report.py` and `bench/o2_default_roster/run.sh` (sha1 in every log).

## Not claimed

This is a simulation. Reads are drawn from the roster's own spliced exon sums (CHM13 sequence) with one error model, so there is no allelic variation, no readthrough, no copy the roster missed, and no difference between the sequenced individual and the
reference (real reads at NPIP copies often match no catalog copy, `docs/PSV_CEILING_2026-10-03.md`): the numbers are an upper bound on real-data accuracy, and the identity-band structure is what transfers. DEV only: chr16 and chr17 are the development blocks of NPIP
and TBC1D3, and no truth for the Y exists beyond annotation and this simulation. Nothing here measures O3 or the families themselves (O1).

## Amendment 1 (2026-10-06, written after the chr16 and chr17 results and before any chrY assignment product exists)

The driver's `assign` stage on the chrY simulation (458 copies, 64 families, one family of 173 copies, 10,224 reads) did not finish inside the 550 s cap of the runner (`timeout` exit 124, no tables written; the cap exists
because a foreground call is limited to 10 minutes). The chrY assignment is therefore run with the two levers the driver's merged stage already applies in its assign phase (`stage_assign_fam`) and documents as byte-identical:
`copy_assign --skip-poa-diagnostic --region-threads 4`, with the driver's own regions file and inputs. To support that equivalence, the same levers are run on chr16 and chr17 and `assignments.tsv`, `families.tsv` and `family_join.tsv` are
compared byte for byte with the default `assign` stage's tables (gate GS4). If they are identical, the three targets are read from the default `assign` tables (chr16, chr17) and the lever run (chrY) and the equivalence is reported; if
not, chr16 and chr17 are read from the default `assign` run, chrY from the lever run is flagged as a different arm, and the difference is reported. If the lever run on chrY also exceeds the cap, it is run detached and the arm is stated.
Rules B1 to B3, targets, read model and population are unchanged.

## Amendment 2 (2026-10-06, written after the chrY lever run was stopped and after the chr16 / chr17 skip tests, before the chrY skip arm was scored; the runner had printed "0 assigned rows of 2,815" for it)

The chrY assignment with the byte-identical levers (Amendment 1) was started detached with no time cap. A child `minimap2 -x splice -N 173` ("star" alignment of 614 sequences of the 173-copy family MCL0 against its 173 copies) used all five cores for 15 minutes without finishing and the run was stopped
(no table was written). None of the 40 chrY target copies is in MCL0. The chrY arm is therefore run with `copy_assign --skip-families` naming MCL0 (arms O2 `asg.skip`, U2 `asg.skipunion`), and every chrY table excludes MCL0's reads (`report.py --exclude-family MCL0`).
**This is NOT byte-identical to a full run in general** (gate GS5): skipping the largest family leaves the rows of all other families identical on chr17 (MCL0, 12 copies: 230 / 77 / 214 rows identical) but changes 19 of 3,219 `assignments.tsv` rows, 5 of 110 `families.tsv` rows and 3 of 342 `family_join.tsv` rows on chr16 (skipping MCL1, the NPIP
family, 30 copies; reads tied between the skipped family and others lose or gain a cross-family demotion). The chrY arm is stated as an approximation of that size; the full chrY run (MCL0 included) is not done and needs a long, unattended, five-core slot.
Rules B1 to B3, targets, read model and population are unchanged.

Note added 2026-10-06 after the user asked how a 173-copy family can exist: MCL0 is the Yq12 (GGAAT)n satellite array, not a gene family (161 of 173 members single-exon, median 9.2 kb, chrY 31-61 Mb, 33% GGAAT/ATTCC, zlib ratio 0.107 against 0.32 for the other families) and carries 51% of the simulated chrY reads; leaving it out is therefore principled, not only a time saving. Rules, targets and read model are unchanged.

