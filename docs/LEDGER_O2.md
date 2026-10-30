# Ledger: O2 (copy assignment of tied reads)

Merged 2026-10-24 from 53 dated docs, oldest first; each kept verbatim under its own heading (its headings one
level down). The originals are in `docs/archive/2026-10/merged/` (git history has them at their old paths).

## Contents

- O2_DEFAULT_ROSTER_2026-10-06.md
- PREREG_o2_default_roster_2026-10-06.md
- PREREG_asj_o2_2026-10-10.md
- PREREG_o2_population_2026-10-12.md
- PREREG_o2_pruning_and_eichler_2026-10-12.md
- PREREG_popfix_2026-10-12.md
- PREREG_psv_semiglobal_engine_2026-10-12.md
- PREREG_psv_significant_blocks_2026-10-12.md
- PREREG_psv_sun_pileup_2026-10-12.md
- PREREG_s1_o2_certificate_2026-10-12.md
- PREREG_yag_o2_2026-10-12.md
- PSV_SUN_PILEUP_2026-10-12.md
- PREREG_o2_anchor_2026-10-13.md
- PREREG_o2_exact_tie_2026-10-13.md
- PREREG_o2_tie_outside_gate_2026-10-13.md
- PREREG_o2_proof_2026-10-17.md
- PREREG_o2_hg002_tied_2026-10-18.md
- PREREG_o2_fair_heldout2_2026-10-19.md
- PREREG_o2_fair_heldout_2026-10-19.md
- PREREG_o2_per_family_2026-10-19.md
- PREREG_o2_regimen_2026-10-19.md
- PREREG_o2_single_genome_2026-10-19.md
- PREREG_samplecopy_haploid_2026-10-19.md
- PREREG_samplecopy_hg002_2026-10-19.md
- PREREG_samplecopy_par_2026-10-19.md
- PREREG_samplecopy_yag_testis_2026-10-19.md
- O2_VAULT_METHODS_2026-10-20.md
- O2_WINDOW_SWEEP_2026-10-20.md
- PREREG_o2_corrector_2026-10-20.md
- PREREG_o2_corrector_families_2026-10-20.md
- PREREG_o2_fair_a119b_fresh_2026-10-20.md
- PREREG_o2_known_family_smn_2026-10-20.md
- PREREG_o2_personal_copies_2026-10-20.md
- PREREG_o2_second_round_2026-10-20.md
- PREREG_o2_second_round_donor_2026-10-20.md
- PREREG_o2_set_null_band_2026-10-20.md
- PREREG_o2_smn_showcase_2026-10-20.md
- PREREG_o2_subs_2026-10-20.md
- PREREG_o2_window_set_2026-10-20.md
- PREREG_o2_corrector_idealsim_2026-10-21.md
- O2O3_EXPLORE_2026-10-22.md
- O2VC_DEV_2026-10-22.md
- O2_NAIVE_TIED_2026-10-22.md
- PREREG_o2_variant_calls_2026-10-22.md
- PREREG_o2_variant_calls_expected_failures_2026-10-22.md
- PREREG_o2_forced_heldout_2026-10-23.md
- PREREG_o2_p9_abstain_2026-10-23.md
- PREREG_o2_semi_p11_2026-10-23.md
- PREREG_o2_v2_heldout_2026-10-23.md
- PREREG_o2_guard_p13_2026-10-24.md
- PREREG_o2_guard_p14_2026-10-24.md
- PREREG_o2_semi_sites_p12_2026-10-24.md
- REPORT_o2_refabsent_2026-10-24.md


---

## O2_DEFAULT_ROSTER_2026-10-06.md

*(original: `docs/archive/2026-10/merged/O2_DEFAULT_ROSTER_2026-10-06.md`)*

## O2 with read-level truth on the default families roster: NPIP, TBC1D3, Y (2026-10-06)

**Status: chr16 (NPIP), chr17 (TBC1D3) and chrY (approximation: the 173-copy family left out, Amendment 2) scored. The full chrY run (that family included) is NOT done.**

Protocol: `docs/PREREG_o2_default_roster_2026-10-06.md` (committed c8753249, amended 9e2deefa for chrY before any chrY assignment product existed). Runner and report: `bench/o2_default_roster/`.
Products: `/mnt/linuxdisk/tmp/o2_default_roster_2026-10-06/{chr16,chr17,chrY}/` (`report.txt`, `report.json`, `score/`). HEAD build in `rustle_target_m2/release` (`copy_assign` 87824d91, `mcl_families` a6308244).
Human A119b, CHM13 v2.0, seed 20260925. **Simulation, DEV, first measurement; no rule was changed after a number was seen.**

### Roster provenance and checks

The default roster (HEAD `assemble` + `families` on the contig) is **byte-identical** to the 2026-09-30 e163d955 products wherever those exist: chr16 and chr17 (assembled GTF, `families.gtf`, clusters, loci, copy table), chrY (clusters, loci, copy table).
Rosters: chr16 372 copies in 111 families, chr17 226 copies in 78 families, chrY 458 copies in 64 families (one family of 173). Simulated reads: 8,225 / 5,843 / 10,224.
Mapping parts (read-disjoint, they only bound one call's wall time): chr16 4, chr17 3, chrY 8.
GS1 holds (one primary-or-unmapped record per simulated read: 8,225 and 5,843). GS2 holds (the per-read tallies equal the ALL rows of `score.py reads` for OWN, ANY and UNION).
**GS4 (Amendment 1):** `copy_assign --skip-poa-diagnostic --region-threads 4` gives `assignments.tsv`, `families.tsv` and `family_join.tsv` byte-identical to the default `assign` stage's on chr16 and on chr17.

### NPIP (chr16 roster; target = roster copies overlapping a CAT/Liftoff NPIP copy: 41 copies)

The aligner decides almost every read of these copies: of 1,099 simulated reads, 371 map at MAPQ 60, 695 at MAPQ 1-59 and **33 at MAPQ 0**. The 33 contested reads come from 5 copies: 3 reads from the NPIP family itself (MCL1: copies 11 and 15), 30 from three other families
that overlap NPIP exons (10 reads each, the depth floor). One of the 41 target copies (MCL1|25) enters through the territory of NPIPB13, the one CAT copy with no exons in the truth GTF; it has no MAPQ-0 read.

| reading (33 MAPQ-0 reads) | correct | wrong | conflict | abstain | no row | assigned | accuracy |
|---|---|---|---|---|---|---|---|
| OWN (the true family's row) | 1 | 0 | 0 | 30 | 2 | 1 | 1.00 |
| PRIMARY (rows with `primary_local = 1`) | 1 | 11 | 0 | 19 | 2 | 12 | 0.08 |
| ANY (any assigned row) | 0 | 28 | 1 | 2 | 2 | 29 | 0.00 |
| UNION (arm U2, beside) | 0 | 0 | 0 | 31 | 2 | 0 | n/a |

**Bars: B1 UNDERPOWERED (1 assigned); B2 PRIMARY UNDERPOWERED (12 assigned); B2 ANY UNDERPOWERED (29 assigned, one under the registered 30), with 0 of 29 correct; B3: nothing assigned.** Aligner-primary on the same 33 reads: 0.36.

All 931 MAPQ-0 reads of the chr16 roster (every family; context with power):

| reading | correct | wrong | conflict | abstain | assigned | accuracy | coverage |
|---|---|---|---|---|---|---|---|
| OWN | 27 | 0 | 0 | 902 | 27 | 1.000 | 0.029 |
| PRIMARY | 48 | 78 | 4 | 799 | 130 | 0.369 | 0.140 |
| ANY | 11 | 393 | 106 | 419 | 510 | 0.022 | 0.548 |
| UNION | 0 | 0 | 0 | 929 | 0 | n/a | 0.000 |

OWN: all 27 assigned reads are in the 99.5-100% identity band, none wrong; 188 of the 931 reads come from copies with an identical sibling and none of them is assigned. The recorded legacy-roster run (2026-09-25) shows the same pattern: OWN 163/0, PRIMARY 112/205, ANY 60/432/156, union 0.

### TBC1D3 (chr17 roster; target = roster copies overlapping a CAT/Liftoff TBC1D3 copy: 22 copies)

Of 887 simulated reads from the 22 target copies, 105 map at MAPQ 60, 782 at MAPQ 1-59 and **none at MAPQ 0**: there is no contested read to assign, so B1, B2 and B3 are UNDERPOWERED with 0 assigned.
All 134 MAPQ-0 reads of the chr17 roster (other families): OWN 6 correct / 0 wrong (coverage .045); PRIMARY 26 / 0; ANY 26 / 0; union 0 assigned. On this contig the readings a consumer sees are right wherever they assign (26 of 26).

### Y (chrY roster; target = roster copies overlapping the body of a RefSeq DAZ, RBMY, TSPY, BPY2, HSFY, VCY, CDY, PRY or XKRY gene: 40 copies in 12 families)

**What the 173-copy family is (user question, 2026-10-06): not a gene family but the Yq12 satellite array.** MCL0's 173 members are 161 single-exon loci (median 9.2 kb, up to 57 kb) tiling chrY 31.4-60.8 Mb; their sequence is 33% GGAAT/ATTCC pentamers (other families: 1%) and compresses 3x better (zlib ratio 0.107 against 0.32): the (GGAAT)n DYZ1 / HSat3 heterochromatin. The default `assemble` turns the reads that align over the array into single-exon pseudo-genes and `mcl_families` joins them (96-99.9% identical; a single exon meets the shared-exon test trivially). The simulation is faithful to that roster and is dominated by it: 5,194 of the 10,224 simulated reads (51%) come from MCL0. The other large Y families are real arrays (MCL1 46, MCL2 36, MCL3 35 copies, all at 5.9-10.2 Mb, the TSPY / FAM197Y region; none single-exon-satellite-like).
**Approximation (Amendment 2): the 173-copy family MCL0 is left out of the assignment run** (it holds none of the 40 target copies, and it is not a gene family). The full assignment stalled in one minimap2 child (`-x splice -N 173`, 614 sequences against 173 copies, five cores, 15 minutes without finishing) and was stopped; skipping a family is byte-identical for the other families on chr17
and changes 0.6% of the assignment rows on chr16, so the numbers below carry that uncertainty.

Of 10,224 simulated reads, 609 MAPQ-0 reads come from MCL0 (excluded) and 2,425 from the other 63 families; **452 MAPQ-0 reads come from the 40 target copies** (134 from copies with an identical sibling, 201 from the 99.5-100% band, 117 from the 98-99% band).
**O2 assigns none of them under any reading** (OWN, PRIMARY, ANY and the union arm all give 0 assigned; 438 abstain, 14 have no row), and none of the 2,425 MAPQ-0 reads of the 63 families either. B1, B2 and B3 are UNDERPOWERED (0 assigned). The aligner-primary baseline puts 41% of the 452 reads on their source copy.
Even the 117 reads of the 98-99% band (copies 1-2% diverged from their closest sibling, where chr16 coverage was 0.80) are not assigned: on this roster the certificate has no decisive column for them.

### Reading, with its limits

- Under the registered read model the contested reads of these families are rare (NPIP 33 of 1,099, TBC1D3 0 of 887): reads simulated from a copy's own sequence with 0.13% error align best to that copy, so the aligner decides them. The registered bars therefore cannot be tested at NPIP and TBC1D3 on this
  design, and "UNDERPOWERED" is the outcome. Real A119b reads at NPIP are contested far more often (24% of primaries at MAPQ 0, `docs/archive/2026-10/PSV_CEILING_2026-10-03.md`): that excess comes from sequence variation that the simulation does not have.
- Where the simulation does have contested reads (all 931 on chr16), the certificate is exact (27 of 27) but assigns 2.9% of them, and the consumer readings are poor (PRIMARY .37, ANY .02): the cross-family claims of `docs/archive/2026-09/PREREG_o2_read_truth_2026-09-23.md` persist on the default roster. At NPIP, 0 of 29 ANY-assigned reads are right.
- P1 (OWN 0 wrong, >= 30 assigned at NPIP and TBC1D3): 0 wrong holds, the >= 30 does not (1 and 0). P2 (consumer readings below .95 at NPIP and TBC1D3): holds at NPIP (0.08 and 0.00 on 12 and 29 reads), untestable at TBC1D3. P3 (Y UNDERPOWERED on B1): holds (0 assigned, approximation). P4 (union assigns almost nothing): holds (0 of 931, 0 of 134, 0 of 2,425).
- Upper bound on real-data accuracy: no allelic variation, no readthrough, no uncatalogued copies, one error model, CHM13 sequence as the source.


---

## PREREG_o2_default_roster_2026-10-06.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_default_roster_2026-10-06.md`)*

## PREREG: O2 with read-level truth on the default families roster, at NPIP, TBC1D3 and the Y-ampliconic genes (2026-10-06)

Written and committed before any roster, simulation or assignment product of this protocol exists. The only thing run before this registration is the report script
(`bench/o2_default_roster/report.py`) on the RECORDED legacy chr16 simulation (`o2sim/human`, 2026-09-25): it reproduces that run's tallies exactly (OWN 163 correct / 0 wrong /
1,084 abstain / 16 other; PRIMARY 112 / 205 / 3 conflict; ANY 60 / 432 / 156; union 0 assigned; n = 1,263 MAPQ-0 reads) and passes the cross-check of `figures/_o2.py`
against the ALL rows of `score.py reads`. That is a check of the instrument on recorded data, not a reading of any arm below.

### Question

Real reads carry no copy of origin, so O2's accuracy comes from simulation (`docs/archive/2026-09/PREREG_o2_read_truth_2026-09-23.md`). On the legacy chr16 catalog the certificate was exact within the
true family (157/157 correct, then 163/163) but the readings a consumer can use were not (PRIMARY accuracy .35, ANY .09 among assigned), and the table `assign` receives in the shipped
pipeline, the DEFAULT families roster (`mcl_families --emit-units`: one copy per member locus), has never been scored with read truth (`docs/archive/2026-09/PREREG_families_copy_table_2026-09-25.md`
holds engineering checks only). At the three families of this evaluation, on that roster, how many contested reads does O2 assign, and are its assignments right?

### Design (the registered protocol, on a different roster)

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

### Target sets (defined now; a roster copy, same chromosome and strand, whose exons overlap >= 1 bp)

- **NPIP**: the exon union of any of the 25 CAT/Liftoff NPIP copies (`copy_recovery_tools_cat/ann`: `copies.hsa.tsv`, `truth.hsa.gtf`).
- **TBC1D3**: the exon union of any of the 16 CAT/Liftoff TBC1D3 copies (same files).
- **Y**: the body of a RefSeq chrY gene (`families_gw/species/human/genes_only.gff`) named DAZ*, RBMY*, TSPY*, BPY2*, HSFY*, VCY*, CDY*, PRY* or XKRY*.
Every table is given for the target's MAPQ-0 reads and for all MAPQ-0 reads of the roster (context). A source copy that overlaps no target truth is not in the target.

### Rules (per target, MAPQ-0 reads of the target's copies; DEV, first measurement; nothing is tuned)

- **B1 (certificate, OWN):** 0 wrong and 0 conflict among OWN-assigned reads, with >= 30 OWN-assigned reads: PASS. Any wrong or conflict: FAIL. Fewer than 30 assigned: UNDERPOWERED (no verdict; the count and,
  with 0 wrong, the rule-of-three upper bound 3 / n on the error rate are reported).
- **B2 (the readings a consumer can use, registered bar of `PREREG_o2_read_truth`):** for PRIMARY and for ANY separately: accuracy among assigned >= 0.95 overall, and >= 0.90 in each of the bands 98-99% and < 98%
  that holds >= 10 assigned reads, with >= 30 assigned reads overall: PASS. Otherwise FAIL, naming the part that fails. Fewer than 30 assigned: UNDERPOWERED.
- **B3 (union arm, beside, no bar):** the B1 test on the UNION reading.
- **No rescue:** a FAIL or UNDERPOWERED outcome is the finding. No read is added, no depth changed, no family dropped. A deeper design would be a new registration.

### Predicted before looking

- P1: OWN has 0 wrong wherever it assigns (157/157 and 163/163 on the legacy roster); NPIP and TBC1D3 reach >= 30 assigned.
- P2: the consumer readings (PRIMARY, ANY) stay below 0.95 at NPIP and TBC1D3 (legacy roster: .35 and .09); one family per NPIP gene array removes some foreign-family claims, not all of them. No prediction at Y.
- P3: Y is UNDERPOWERED on B1 (real chrY reads: 12 of 3,641 contested molecules assigned; the arms are 99.9% identical).
- P4: the union arm assigns almost nothing (0 of 1,263 on the legacy roster).

### Checks

GS1: the simulated BAM holds exactly one primary-or-unmapped record per simulated read (`figures/_o2.py`). GS2: the per-read tallies equal the ALL rows of `score.py reads` (OWN, ANY, UNION). GS3: the roster's byte comparison with the stored products is reported.
The scorer, the report and the runner are `bench/score.py`, `bench/o2_default_roster/report.py` and `bench/o2_default_roster/run.sh` (sha1 in every log).

### Not claimed

This is a simulation. Reads are drawn from the roster's own spliced exon sums (CHM13 sequence) with one error model, so there is no allelic variation, no readthrough, no copy the roster missed, and no difference between the sequenced individual and the
reference (real reads at NPIP copies often match no catalog copy, `docs/archive/2026-10/PSV_CEILING_2026-10-03.md`): the numbers are an upper bound on real-data accuracy, and the identity-band structure is what transfers. DEV only: chr16 and chr17 are the development blocks of NPIP
and TBC1D3, and no truth for the Y exists beyond annotation and this simulation. Nothing here measures O3 or the families themselves (O1).

### Amendment 1 (2026-10-06, written after the chr16 and chr17 results and before any chrY assignment product exists)

The driver's `assign` stage on the chrY simulation (458 copies, 64 families, one family of 173 copies, 10,224 reads) did not finish inside the 550 s cap of the runner (`timeout` exit 124, no tables written; the cap exists
because a foreground call is limited to 10 minutes). The chrY assignment is therefore run with the two levers the driver's merged stage already applies in its assign phase (`stage_assign_fam`) and documents as byte-identical:
`copy_assign --skip-poa-diagnostic --region-threads 4`, with the driver's own regions file and inputs. To support that equivalence, the same levers are run on chr16 and chr17 and `assignments.tsv`, `families.tsv` and `family_join.tsv` are
compared byte for byte with the default `assign` stage's tables (gate GS4). If they are identical, the three targets are read from the default `assign` tables (chr16, chr17) and the lever run (chrY) and the equivalence is reported; if
not, chr16 and chr17 are read from the default `assign` run, chrY from the lever run is flagged as a different arm, and the difference is reported. If the lever run on chrY also exceeds the cap, it is run detached and the arm is stated.
Rules B1 to B3, targets, read model and population are unchanged.

### Amendment 2 (2026-10-06, written after the chrY lever run was stopped and after the chr16 / chr17 skip tests, before the chrY skip arm was scored; the runner had printed "0 assigned rows of 2,815" for it)

The chrY assignment with the byte-identical levers (Amendment 1) was started detached with no time cap. A child `minimap2 -x splice -N 173` ("star" alignment of 614 sequences of the 173-copy family MCL0 against its 173 copies) used all five cores for 15 minutes without finishing and the run was stopped
(no table was written). None of the 40 chrY target copies is in MCL0. The chrY arm is therefore run with `copy_assign --skip-families` naming MCL0 (arms O2 `asg.skip`, U2 `asg.skipunion`), and every chrY table excludes MCL0's reads (`report.py --exclude-family MCL0`).
**This is NOT byte-identical to a full run in general** (gate GS5): skipping the largest family leaves the rows of all other families identical on chr17 (MCL0, 12 copies: 230 / 77 / 214 rows identical) but changes 19 of 3,219 `assignments.tsv` rows, 5 of 110 `families.tsv` rows and 3 of 342 `family_join.tsv` rows on chr16 (skipping MCL1, the NPIP
family, 30 copies; reads tied between the skipped family and others lose or gain a cross-family demotion). The chrY arm is stated as an approximation of that size; the full chrY run (MCL0 included) is not done and needs a long, unattended, five-core slot.
Rules B1 to B3, targets, read model and population are unchanged.

Note added 2026-10-06 after the user asked how a 173-copy family can exist: MCL0 is the Yq12 (GGAAT)n satellite array, not a gene family (161 of 173 members single-exon, median 9.2 kb, chrY 31-61 Mb, 33% GGAAT/ATTCC, zlib ratio 0.107 against 0.32 for the other families) and carries 51% of the simulated chrY reads; leaving it out is therefore principled, not only a time saving. Rules, targets and read model are unchanged.


---

## PREREG_asj_o2_2026-10-10.md

*(original: `docs/archive/2026-10/merged/PREREG_asj_o2_2026-10-10.md`)*

## Allele-specific junctions as an O2 extension: haplotype assignment by the O2 certificate, allele vs copy by the haplotype assemblies, 2026-10-10

**Pre-registered before any read is phased or any junction tested.** User (2026-10-10): bring the dropped allele-specific-junction result (ASJ, June:
475 calls / 235 genes on `GGO.bam`, OR6737 testis, read-anchored, `project_allele_specific_junctions`) back as a sub-result of O2: re-run it on KB3781's own
fibroblast reads, phase reads with O2's rule instead of the June 0.8-agreement 2-means, and use mGorGor1's maternal / paternal assemblies (the same
animal) to tell an allele from a copy.

**Why O2.** O2 assigns a read to one of a set of sequences that differ at known columns (PSVs), or abstains. A haplotype is the same object one level down:
two sequences of one locus that differ at heterozygous SNVs. The ASJ test is then a downstream comparison of junction usage between the two assigned groups.
Junctions are NOT used to phase (they are what is tested).

### Inputs (fixed)

- Reads: `fibroblasts/GCA_029281585.2_flnc_mm.bam` (KB3781 fibroblast Iso-Seq on the primary, `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1`). Primary
  alignments only (no secondary / supplementary / unmapped); no MAPQ filter (paralogous reads are part of what is measured; per-call MAPQ-0 share reported).
- Genes: `GGO_genomic.gff` gene records on the 23 autosomes of `rna_allele/chrmap.tsv` (and chrX for the control), with >= 2 exons in some transcript and
  >= 20 primary reads overlapping the gene span.
- Haplotypes: each primary chromosome is byte-identical to one haplotype chromosome (`same_hap`); the other haplotype's chromosome is `B_name`
  (`gorilla_haps/{mat,pat}.fa`). Het sites of a gene = SNVs between the primary gene span (+ 1 kb flanks) and its alignment to `B_name`
  (`minimap2 -x asm5 -c --cs`, best alignment on `B_name` covering >= 90% of the span); SNVs within 5 bp of an indel in that alignment are dropped.
  chrX has no `B_name` (KB3781 is male): no het sites, no phasing.

### Arm A — assembly-phased, O2 rule (the new method)

Per read: at each het site it covers with a base (not deleted / not in an intron), the observation is the read base. **O2's certificate as in
`copy_pair_significance`** (`src/family.rs`): for the MLE haplotype `best` (log-likelihood, per column ln(1 - eps) if the read base equals that
haplotype's allele, else ln eps) vs the other, eps_j = 0.003 / 3 per distinguishing column (O2's `error_rate` default, no base qualities, editing filter
off), k = columns where the read matches `best`, p = P(X >= k), X ~ Poisson-binomial(eps). **Assign iff `best` is the strict MLE and p < alpha / (n - 1) =
1e-3** (n = 2); else abstain. One site never assigns (p = eps = alpha); two matching sites do.

### Arm R — read-anchored (June method, replicated on KB3781)

Per gene, candidate anchors from the reads (`--eqx`: X operations): a position with two alleles, each on >= 5 reads, minor fraction 0.35-0.65, depth >= 12;
the anchor = the candidate with the highest depth (ties: smallest position). Reads split by their base at the anchor.

### Junction test (both arms, June constants)

Per gene and group split: introns = the exact `N` operations of the reads, kept when >= 3 reads carry them. A read SPANS an intron when its alignment covers
both ends; PSI_g = reads of group g using it / reads of group g spanning it. Junctions with PSI = 1 in both groups are not tested. Fisher exact (two-sided)
on the 2 x 2 table; Benjamini-Hochberg over all tests of the arm (all chromosomes in the half scored); **ASJ = q < 0.05 and |dPSI| >= 0.3.**

### Allele vs copy (the haplotype assemblies as referee)

Each arm-R anchor is **ALLELIC** when it is an arm-A het site and the two read alleles are the primary and `B_name` bases; **NON-ALLELIC** when the `B_name`
alignment covers it (not within 5 bp of an indel) with the SAME base as the primary (homozygous in the animal: the second read allele comes from a paralog,
an edit or an error); **UNRESOLVED** otherwise (no qualifying `B_name` alignment, near an indel, or a third base).

### Tests

**Split.** Chromosomes by `chrmap.tsv` number: odd (1, 3, ..., 23) = DEV, scored and written first; even = HELD-OUT, computed after; at most ONE amendment
between them; the held-out half decides T1 and T2. Nothing is tuned: every constant is O2's default or the June ASJ value.

- **T1 (O2's bound holds on haplotypes).** Each gene's het sites split by site index parity; a read assigned by BOTH halves separately is concordant when
  the two agree. **Bar: discordance <= 0.005** over such reads (the certificate implies about 0.002 for two independent assignments at alpha 1e-3).
- **T2 (the assemblies separate allele from copy).** For each arm-R anchor with >= 10 reads assigned by arm A using the gene's OTHER het sites
  (anchor excluded), concordance = max over the two pairings of the share of assigned reads whose anchor allele matches the pairing. **Bar: median
  concordance >= 0.95 over ALLELIC anchors AND <= 0.80 over NON-ALLELIC anchors** (each with >= 10 anchors; fewer = UNDERPOWERED).
- **T3 (descriptive, all autosomes).** Arm-A ASJ calls (genes, calls, median |dPSI|); arm-R ASJ calls by referee class (the copy / edit-confounded share
  = NON-ALLELIC / resolved, beside June's 44 / 120); overlap of arm-R ALLELIC calls with arm-A calls (same gene and junction); per class the share of
  reads with MAPQ 0 and, among NON-ALLELIC anchors, transitions A/G vs others (editing vs paralog); abstention rate of arm A.
- **Control (descriptive).** chrX >= 5 Mb from either end (male, hemizygous, PAR excluded): arm-R candidate anchors per tested gene and arm-R ASJ calls,
  beside the autosomal rates; every such call is non-allelic by biology.

**Claim if T1 and T2 hold:** "O2's certificate assigns reads to haplotypes with its stated error bound (T1), and the animal's haplotype assemblies separate
allele-specific from copy-specific junction usage (T2); N junctions are allele-specific in KB3781, and X% of the read-anchored calls are copy- or
edit-confounded." If T1 fails, the bound does not transfer to haplotypes; if T2 fails, the referee is not established. Not claimed either way: a mechanism
(splice-site disruption; retracted for the June flagships), or anything about OR6737.

Code: `bench/asj_o2/` (test-first). Compute: per chromosome, foreground, `tools/rlock.sh heavy` for the alignment step.

### DEV result (2026-10-10; code ce1769f0 + two speed-ups shown byte-identical on chr21, scorer 23 tests; held-out NOT computed)

12 odd chromosomes, 6,807 expressed multi-exon genes (6,562 with an allelic alignment; median 27 het sites per gene).
**Arm A phasing:** 1,270,013 reads to the primary's haplotype, 1,230,894 to the other, 5,211,093 abstained (67.6%: reads covering < 2 het sites or tied).
- **T1 PASSES on dev: 872 discordant of 883,715 reads assigned by both halves = 0.00099** (bar 0.005; the certificate implies about 0.002).
- **T2 PASSES on dev:** median concordance **1.000 over 2,585 ALLELIC anchors** (bar >= 0.95) and **0.642 over 29 NON-ALLELIC anchors** (bar <= 0.80);
  UNRESOLVED 1.000 over 68.
- T3 (dev): **arm A 289 ASJ / 137 genes** (median |dPSI| 0.552; 69,633 tests). **Arm R 485 ASJ / 172 genes: ALLELIC 187, NON-ALLELIC 78, UNRESOLVED 220**;
  NON-ALLELIC = **29.4% of the resolved calls** (June, other animal and criteria: 44 / 120). 141 of the 187 ALLELIC arm-R calls are also arm-A calls.
  Anchors: ALLELIC 3,999, NON-ALLELIC 135, UNRESOLVED 822; transitions 0.75 of ALLELIC anchors vs **0.90 of NON-ALLELIC anchors** (an editing signature);
  median MAPQ-0 share 0 in every class (the non-allelic anchors are not low-MAPQ reads). Examples: ATXN3 two donors 7 bp apart switch with the allele
  (dPSI +0.89 / -0.97); AGTRAP two acceptors 21 bp apart (+0.68 / -0.68).
**No amendment.** The held-out half and the chrX control are run with the same code.

### HELD-OUT result, verdict and T3 (2026-10-10; no amendment; dev chromosomes run with the code before three speed-ups, each shown byte-identical on chr21)

**HELD-OUT (11 even chromosomes, 6,841 genes):**
- **T1 PASSES: 1,850 discordant of 1,285,832 = 0.00144** (bar 0.005).
- **T2 PASSES: median concordance 1.000 over 2,455 ALLELIC anchors and 0.577 over 18 NON-ALLELIC anchors** (bars >= 0.95 / <= 0.80); UNRESOLVED 1.000 over 88.
**Verdict: the claim holds.** O2's certificate assigns reads to haplotypes within its stated bound on both halves (discordance 0.0010 dev / 0.0014
held-out, against about 0.002 implied), and the haplotype assemblies separate allele-tracking anchors from non-allelic ones.

**T3, all 23 autosomes (13,648 genes):** arm A assigns 2,680,826 + 2,652,380 reads, abstains on 66.1%; **636 allele-specific junctions in 286 genes**
(median |dPSI| 0.567). Arm R 879 calls / 352 genes: ALLELIC 442, NON-ALLELIC 97, UNRESOLVED 340; **NON-ALLELIC = 18.0% of resolved calls**; 353 of the 442
ALLELIC arm-R calls are also arm-A calls. Anchors: ALLELIC 7,887, NON-ALLELIC 252, UNRESOLVED 1,723; **transitions 0.91 of NON-ALLELIC vs 0.75 of
ALLELIC anchors**, median MAPQ-0 share 0 in both: in KB3781 the read-anchored method's non-allelic anchors look like RNA editing more than paralogous
reads (June's confound on OR6737 was named paralogous, LOC*; not compared directly).
**Control (chrX >= 5 Mb from the ends, 424 genes, male):** 0.059 arm-R anchors per gene (autosomes 0.723), 92% transitions; **1 arm-R call in 424 genes
(0.0024 per gene; autosomes 0.064)**, UNRESOLVED (no other haplotype). The read-anchored method's floor on a hemizygous chromosome is about 4% of its
autosomal call rate.
**Not claimed:** a splicing mechanism; anything about OR6737; that UNRESOLVED arm-R calls (39%) are allelic (their anchors track the haplotype, median 1.000
over 156, but the referee did not classify them).


---

## PREREG_o2_population_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_population_2026-10-12.md`)*

## How many multi-placement reads can O2 assign, and how many does it move off the aligner's primary? 2026-10-12

**Pre-registered before the script runs.** User (2026-10-12): "how common are assignable reads to copies, do we have a lot of cases in
which it is impossible or do we have enough cases of a decent number of reads correctly assigned that StringTie and FLAIR cannot do
since they don't use secondary alignments".

**Population** (per species, benchmark BAMs `GGO_mm.bam` / `A119b.t2t.bam`, all chromosomes): molecules with another alignment within
the good window — `as_table`'s second-best AS >= 0.98 x best AS. **Sample**: 20,000 molecules per species, uniform, seed 1
(`random.Random(1).sample` over the population in table order). **Per molecule**: the primary record vs every non-supplementary
record with AS >= 0.98 x best (O2's certificate exactly as `bench/tool_difference/s1_o2.py`). Classes, first that applies:
**moved** (some rival certified over the primary: the aligner's — and so StringTie's / FLAIR's — placement is wrong and O2 can
correct it), **confirmed** (the primary certified over every rival), **identical** (no distinguishing column against any rival:
no read-level method can assign it), **undecided** (otherwise). Reported with 95% Wilson intervals and extrapolated to the
population; also split by the primary's AS rank (true tie primary = best = second; primary below the best; primary best with a close
rival), the mosaic share, and, for moved molecules, whether the winning placement is spliced where the primary is not (the
retrocopy direction). Descriptive; nothing adopted. Expected (written now): identical and confirmed dominate; moved is a small share.

**Instrument** (committed before running): `bench/tool_difference/o2_population.py` sha1 1a1180dcdf3b5358c2d254fa67a155da8c69dfe1 (2 tests).

### Result (2026-10-12)

Outputs `/mnt/linuxdisk/tmp/asm_gap/op/{gorilla,human}.o2pop.{txt,molecules.tsv}`.

| | gorilla | human |
|---|---|---|
| molecules / multi-placement population | 4,404,440 / 49,444 (1.1%) | 21,391,966 / 532,876 (2.5%) |
| **moved** (O2 certifies a rival over the aligner's primary) | 1.0% [0.9, 1.2] ~ 504 reads | 2.7% [2.5, 3.0] ~ 14,494 reads |
| **confirmed** (primary certified over every rival) | 76.0% ~ 37,565 | 62.9% ~ 335,046 |
| **identical** (no distinguishing column: unassignable by any read-level method) | 8.2% ~ 4,035 | 10.2% ~ 54,380 |
| **undecided** (too few columns) | 14.8% ~ 7,340 | 24.2% ~ 128,956 |
| true ties (primary = best = second): moved / confirmed / identical / undecided | 0.6 / 12.3 / 59.2 / 28.0% | 1.6 / 3.2 / 48.1 / 47.1% |
| primary below the best AS | 0.5 / 96.7 / 0.8 / 2.0% | 4.5 / 70.8 / 0.8 / 23.9% |
| primary best, close rival | 1.4 / 78.5 / 0.9 / 19.2% | 2.4 / 77.2 / 2.7 / 17.7% |
| mosaic (>= 2 columns for each of a pair) | 18.8% | 39.7% |
| moved where the winner has more introns (retrocopy reversal) | 8 / 204 | 24 / 544 |

**Reading.** Most multi-placement reads are assignable and the aligner's primary is already right (confirmed 76 / 63%); a further
8-10% are truly unassignable (identical copies) and 15-24% lack evidence. The reads StringTie and FLAIR misplace and O2 can move are
1.0% / 2.7% of the multi-placement population — about 500 gorilla and 14,500 human reads genome-wide (0.01% / 0.07% of all reads).
"Primary below the best AS" is confirmed 97% / 71% of the time: the aligner's primary, not the best-AS placement, is usually the
right one (the retrocopy pattern of register 1270). Mosaic reads are common (19 / 40%): pairs where the read carries each copy's
allele at >= 2 columns — the O3 regime (variants absent from both reference copies). Expected outcome (written before) held.
Open: where the moved reads land (how many copies get >= 2 of them, i.e. copies the tools cannot assemble) needs the full population.

### Amendment 1 (2026-10-12, user: "yes lets do that"; before any full-population number)

The same classification over the FULL population (`--n 0`; names processed in batches of 100,000 per BAM pass to bound memory;
identical per-molecule rule), recording for every moved molecule the winning placement (chromosome, aligned span). Aggregated by
`bench/tool_difference/o2_moved_loci.py`: each winning placement is mapped to the reference genes whose exons its aligned span
overlaps (none = unannotated); reported: genes receiving >= 1 / 2 / 5 / 10 moved reads, and for genes with >= 2, whether OP_good,
StringTie and FLAIR (benchmark GTFs) emit any transcript touching the gene's exons. A gene with >= 2 moved reads and no StringTie /
FLAIR transcript is a copy the tools miss and O2's reads could support. Descriptive; nothing adopted.

**Instruments for amendment 1** (committed before the full run; the batched script reproduces the 20,000 sample exactly with 3 batches): `o2_population.py` sha1 f8d7435e304215df6c0c510c87179c3140a1ebec, `o2_moved_loci.py` sha1 82bf3325f87839dca3452370d16abe704fa1c958 (4 tests).

### Amendment 1 result (2026-10-12): the full population

Outputs `/mnt/linuxdisk/tmp/asm_gap/op/{gorilla,human}.o2full.{txt,molecules.tsv,loci.txt}`. The full populations reproduce the
sample estimates.

| | gorilla (49,444) | human (532,876) |
|---|---|---|
| moved / confirmed / identical / undecided | 519 (1.0%) / 76.4% / 7.9% / 14.7% | 14,945 (2.8%) / 62.5% / 10.6% / 24.2% |
| moved reads on annotated genes / unannotated | 306 / 213 | 4,697 / 10,248 |
| genes receiving >= 1 / 2 / 5 / 10 moved reads | 188 / 66 / 22 / 16 | 1,115 / 588 / 262 / 125 |
| of the genes with >= 2: OP_good / StringTie / FLAIR emit a transcript there | 41 / 38 / 33 of 66 | 521 / 511 / 528 of 588 |
| genes with >= 2 moved reads and NO StringTie or FLAIR transcript | 25 (390 reads; ours touches 1) | 44 (227 reads; ours touches 10) |

(Gene counts count a read once per overlapping gene: the gorilla top entries are one 47-read placement over a tRNA / LOC cluster of 7
annotations; the human top entries are LOC / ncRNA / snRNA / miRNA genes.)
**Reading.** The reads O2 moves mostly land where the tools already emit transcripts (they misplace some reads but still assemble the
gene from others) or on unannotated sequence (41% gorilla, 69% human). Copies that the tools miss entirely and that >= 2 O2-moved
reads support are few — 25 / 44 genes, mostly tRNA / ncRNA / LOC annotations — and our current assembly touches 1 / 10 of them.
"Reads only O2 can place" does not translate into many "copies only we can assemble". Open: the 10,248 human (213 gorilla) moved reads
on UNANNOTATED sequence — clustered into loci with >= 2 reads they are candidate unannotated copies (O3's detection question).


---

## PREREG_o2_pruning_and_eichler_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_pruning_and_eichler_2026-10-12.md`)*

## O2 on the coin-toss genes, O2-pruned good pool, and Eichler's AS-margin rule, 2026-10-12

**Pre-registered before any S2 read is tested and before any pruned arm exists.** Follows register 1269 (secondary-only genes S1:
1 O2-supported gene per species; most are contradicted — paralog shadows). User (2026-10-12): "lets try on coin toss genes, and also
drop the transcripts in contradicted ones; also can we compare this with Eichler's rule that the copy belongs to the alignment that
is within 10 units of the AS?"

### 1. The certificate on S2 (coin-toss-starved genes)

`bench/tool_difference/s1_o2.py --stratum S2`, otherwise exactly as `docs/PREREG_s1_o2_certificate_2026-10-12.md` (same pairs: each
good secondary record overlapping the gene's exons vs the molecule's primary; same certificate; same gene classes). The S1 numbers
must reproduce exactly with the extended script (`--stratum S1`, the default).

### 2. Eichler's AS-margin rule beside O2 (per record, both strata)

As implemented 2026-09-21 (`docs/archive/2026-09/EICHLER_COMPARISON_2026-09-21.md`, register 935-936): a molecule is assigned to its
best-AS alignment iff no other alignment of it scores within 10 AS units (best - second >= 10, or no second alignment); otherwise it
is discarded. Per tested record (secondary in the stratum gene = copy B, primary = copy A): **B** if B's AS is the molecule's best and
the margin >= 10; **A** if the primary's AS is the best and the margin >= 10; **other** if the best is a third placement with margin
>= 10; **discard** if the margin < 10. Eichler gene classes, parallel to O2's: supported >= 2 records to B; contradicted 0 to B and
>= 2 to A or other; undecidable otherwise. Reported: per-record agreement where both rules call (O2 A/B vs Eichler A/B), and the
gene-class cross-table.

### 3. Pruned good-pool arms

For a rule R (O2 or Eichler), `OP_good_R` = `OP_good` minus every transcript that (a) touches an exon of a gene R calls contradicted
(S1 or S2) and (b) is GOOD-POOL-ONLY: its intron chain is absent from `OP_primary` on that chromosome (mono-exon: no `OP_primary`
transcript overlaps it). Transcripts built from the gene's own primary reads stay. Scored with `bench/tool_difference/score_op.py`
(same reference, strata and development-chromosome exclusion as register 1268), beside the existing arms.

### Claims and expectations (both species)

- **E1 (S2 certificate).** Written now: tied secondaries in S2 are certified to the primary copy more often than to the S2 copy.
- **E2 (pruning, annotation metrics).** Expected: genome-wide metrics move by < 0.5 point; stratum precision moves little, because a
  shadow can match the annotation as well as a real transcript (register 1269: reference-exact S1 transcripts sit mostly in
  contradicted genes). An annotation-only truth cannot score copy-of-origin; the number of shadow transcripts removed is the effect.
- **E3 (O2 vs Eichler).** Expected: where both call, they agree mostly (09-21: 97.1%); Eichler discards most of these records
  (margins are small among good secondaries), so it contradicts fewer genes and prunes fewer shadows.
Nothing is adopted: the outcome is descriptive, for the advisor's figure discussion. A test of copy-of-origin needs a truth with known
source copies (the family simulator, `bench/famsim`), named here as the follow-up.

**Instruments** (committed before running): `s1_o2.py` sha1 a79218c0bc6aa231f40cb7eaca7bca5da1ad1d12 (adds `--stratum` and the Eichler columns; 9 tests), `prune_good.py` sha1 4d17bbe843d8e99c4f406de683e9019abdfe37ca (2 tests).

### Result (2026-10-12)

Outputs `/mnt/linuxdisk/tmp/asm_gap/op/{gorilla,human}.v2.{S1,S2}.{txt,records.tsv,genes.tsv}`, `score_{gorilla,human}_pruned.txt`. The
extended `s1_o2.py` reproduces register 1269's S1 records exactly in both species.

**1. Certificate on S2 (coin-toss-starved genes).**

| | gorilla S2 | human S2 |
|---|---|---|
| records tested | 56,905 | 597,412 |
| AS == best ("tied"): to S2 copy / to primary copy / abstain | 75 / 41,569 / 5,527 | 386 / 82,882 / 240,068 (235,673 with no distinguishing column) |
| strictly worse: to S2 copy / to primary copy / abstain | 24 / 8,193 / 1,517 | 902 / 66,925 / 206,249 |
| genes supported / contradicted / undecidable | 9 / 200 / 77 | 89 / 424 / 374 |

**E1 holds** (tied S2 secondaries certified to the primary copy 41,569 vs 75; 82,882 vs 386).

**2. Eichler's AS-margin rule beside O2.** Records both rules call A or B: agreement S1 0.356 / 0.540, S2 0.141 / 0.300 (gorilla / human).
**E3's "agree mostly" FAILS**: Eichler sends 5,780 / 2,485 (S1) and 35,649 / 75,063 (S2) records to the stratum copy, O2 almost
none. The cause, checked on 300 sampled disagreeing gorilla reads per stratum (Eichler -> stratum copy, O2 -> primary copy): the
best-AS alignment has NO intron in 95% (S1) / 92% (S2), fewer introns than the primary in 99% / 100%, more mismatches in 100% / 99%
(median +23 in S1). These are spliced parent-gene reads that also align unspliced to a processed (intronless) copy; minimap2's AS
is higher there because no splice penalty is paid, so an AS-margin rule assigns them to the retrocopy, while the PSV bases put them
back on the parent. E3's other half holds: Eichler contradicts fewer genes (187 vs 600 gorilla; 563 vs 776 human) and prunes fewer
transcripts (214 vs 416; 804 vs 969).

**3. Pruned good-pool arms** (genome-wide minus the development chromosome; stratum = reference transcripts matched exactly / precision).

| arm | gorilla chain | gorilla S1 | gorilla S2 | human chain | human S1 | human S2 |
|---|---|---|---|---|---|---|
| OP_good | 29.4 / 30.0 | 34 / 19.9% | 84 / 24.5% | 23.8 / 15.0 | 36 / 16.5% | 210 / 15.9% |
| OP_good_o2 | 29.4 / 30.1 | 16 / 34.8% | 72 / 34.6% | 23.7 / 15.0 | 0 / 0.0% (16 tx) | 119 / 15.0% |
| OP_good_eichler | 29.4 / 30.1 | 20 / 29.4% | 77 / 30.4% | 23.7 / 15.0 | 5 / 10.6% | 119 / 13.9% |
| StringTie | 24.4 / 34.4 | 0 | 51 / 39.0% | 18.1 / 15.0 | 0 | 77 / 13.3% |
| FLAIR | 25.5 / 17.3 | 0 | 42 / 36.8% | 22.4 / 6.2 | 0 | 72 / 8.0% |

**E2**: genome-wide metrics move <= 0.3 point (holds). Stratum precision is NOT flat as predicted: in gorilla O2 pruning raises it
(S1 19.9 -> 34.8%, S2 24.5 -> 34.6%); in human it removes every S1 match (all 36 sat in contradicted genes) and lowers S2 precision
slightly (15.9 -> 15.0%) while dropping 91 matches. After O2 pruning, the S2 comparison with the tools: gorilla 72 matches vs 51 / 42
at 34.6% vs 39.0 / 36.8% precision; human 119 vs 77 / 72 at 15.0% vs 13.3 / 8.0%.
**Reading.** (a) O2's certificate and Eichler's AS-margin rule disagree mostly because of retrocopies: AS rewards unspliced
alignments to processed copies, the PSV bases do not; on this population O2's call is the biologically right one. (b) Neither rule
makes the secondary-only stratum a win: after pruning, gorilla keeps 16 S1 matches and human none. (c) In human S2 most records have
no distinguishing column (identical copies, 380,372 of 597,412): the per-gene "contradicted" rule (>= 2 certified elsewhere, 0 here)
then rests on a minority of reads, and pruning removes transcripts that match the annotation 20% of the time (221 / 1,107) — the
pruning rule is too coarse where copies are identical. Copy-of-origin needs a simulation truth (`bench/famsim`) to decide.


---

## PREREG_popfix_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PREREG_popfix_2026-10-12.md`)*

## PREREG — population-fixed PSVs for O2 on the human Y (2026-10-12 plan; written 2026-10-13, before any run below)

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

### Arms and bars

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

### Amendment 1 (2026-10-13, before any run)

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

### Results (2026-10-13; binary ba525453, scorers 91fc1635; outputs `/mnt/linuxdisk/tmp/popfix/runs/`)

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


---

## PREREG_psv_semiglobal_engine_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PREREG_psv_semiglobal_engine_2026-10-12.md`)*

## PREREG — a log-odds, semi-global, anchored-band engine for O2's copy-vs-copy PSV alignment (2026-10-12)

**Why.** `PREREG_psv_significant_blocks_2026-10-12.md` (B2 FAILED, mask opt-in) showed that `banded_msa_pair`'s cost
model (mismatch 1, gap open 32, extend 1; global) makes a ~30-base misalignment cheaper than one gap: inside stretches
an independent aligner supports, 9,703 of 16,183 of O2's apparent copy differences (60%) sit at the wrong position, and
14% of all PSV pairs come from non-homologous sequence forced into the global path. The user also asked that this part
be VERY FAST.

### Engine (fixed before any held-out data is touched; may be tuned on DEV only, each change logged as an amendment)

`semiglobal_msa_pair(a, b)` — same 2-row gapped-MSA output as today's engines, so `discover_psvs` is unchanged downstream:
- **Scores** (maximised): match +2, mismatch or non-ACGT -3; **two-piece affine gaps** cost min(5 + 2k, 24 + k) for a gap
  of k (minimap2-style: a long gap — an exon one copy lacks — costs ~1 / base, less than the -1.75 / base expectation of
  aligning unrelated sequence).
- **Free end gaps** in both sequences (overlap alignment): an unshared 5' / 3' extent is left unaligned, not forced.
- **Anchored band**: unique shared 11-mers, longest colinear chain; row i between anchors (i1, j1), (i2, j2) is
  computed for j in [i + min(d1, d2) - W, i + max(d1, d2) + W] (d = j - i), W = 100; before the first / after the last
  anchor the band follows that anchor's diagonal. No anchor -> no alignment (no shared 11-mer: no detectable homology).
- **Guard**: `significant_block_mask` (Karlin-Altschul, alpha 1e-3) on the result.
- Selected by `RUSTLE_PSV_ENGINE=semiglobal` (opt-in) until this prereg passes; default unchanged.

### Harness (no BAM)

`utilities psv-discover --families COPIES.tsv --copies-fa COPIES.fa --out PREFIX [--engine global|semiglobal]`
runs `discover_psvs_oriented` per family and writes `PREFIX.psv_cols.tsv` / `.psv_copies.tsv` (the dump's formats),
`PREFIX.psv_pos.tsv` (every copy's spliced offset and base per column) and `PREFIX.time.tsv` (wall per family).
Validated first: `--engine global` reproduces the 10-10 `copy_assign --dump-psv` columns on the dev subset (a column
the RNA-editing pre-pass removed in copy_assign may be extra in the harness; nothing else may differ).

### Yardstick (`bench/psv_sun/engine_eval.py`)

Per family, copy 0 vs each copy X in copy 0's frame, the colinear chain of `minimap2 -cx asm20` pieces (independent of
both engines). A **difference pair** = (X, column) where X's allele differs from copy 0's.
- In a chained piece: **V** if the engine's X offset equals minimap2's mapping of copy 0's offset, else **M**.
- Outside every piece: Smith-Waterman (`swscore`, BLASTN scores) of the surrounding unaligned stretches vs 20 shuffles:
  **H** (beats all; divergent homolog) or **N** (does not; unrelated).

### Substrates

- DEV (may be looked at freely): the 56 KB3781 families of `PREREG_psv_sun_pileup_2026-10-12.md` (`run/S.txt`).
- **HELD-OUT 1**: the 269 other KB3781 families with 2-20 copies (`rustle_figures/runs/gorilla_KB3781/gorilla_KB3781.fam.copies.tsv`).
- **HELD-OUT 2** (another species and sample): the 2,101 human A119b families with 2-20 copies
  (`rustle_figures/runs/human_A119b/human_A119b.fam.copies.tsv`).

### Bars (each held-out substrate separately; semiglobal vs global)

- **P1 precision**: V / (V + M) >= 0.95 for semiglobal.
- **P2 no loss**: V(semiglobal) >= 0.99 x V(global).
- **P3 junk**: N(semiglobal) <= 0.10 x N(global).
- **P4 speed**: total harness wall (same threads) semiglobal <= global; and no family slower than 2 s.
- **P5**: library tests pass; on DEV `copy_assign` default assignments stay byte-identical with the engine switched.
- Reported, no bar: H counts, total columns, per-family time distribution.
All of P1-P4 on BOTH held-out substrates -> the engine becomes the default. Otherwise it stays opt-in and the failing
bar is reported.

### Dev amendments (2026-10-12, DEV only; written before any held-out family is aligned or scored)

- **E1 engine (dev test failure):** free end gaps let the path start on the first row OR column only, so two unshared
  prefixes could not both be skipped and 60 unrelated bases were aligned (unit test). The engine is now LOCAL at its ends
  (Smith-Waterman: start / end anywhere) and the long-gap piece is a CONSTANT, cost min(5 + 2k, 30) (like an intron in a
  spliced alignment): a missing exon of any length costs 30, and unrelated stretches on both sides can be skipped by two
  gaps. Scores, anchors, band and the significant-block guard unchanged. 6 engine tests + 1,029 library tests pass.
  Anchor limit documented: a PERIODIC 20% difference pattern shares no exact 11-mer and gets no alignment.
- **Harness validated:** `psv-discover --engine global` reproduces the dev `copy_assign --dump-psv` tables exactly
  (sorted: 33,016 / 209 lines) in 6 s instead of a 5-min BAM run.
- **A1 yardstick (dev P1 = 0.736):** 1,567 of the 1,725 M pairs sat > 20 bp from minimap2's position, in copies with
  INTERNAL repeats: MCL8's copy 0 is a ~42-bp-period tandem array (323 off-diagonal self-alignments), MCL37 an 84-bp
  shift, MCL32's copy X a duplicated segment matching copy 0's offset equally well at two places (local identity 0.84
  vs 0.84). Positional homology is undefined there for ANY aligner. New class **R**: a pair inside a minimap2 piece whose
  copy-0 offset lies in an off-diagonal self-alignment of copy 0, or whose minimap2 X position lies in one of X
  (`minimap2 -cx asm20 -X -N 50 -p 0.1`; engine-independent). R is reported and excluded from V and M. Bars unchanged.
- Dev with A1: global V 3,192 / M 7,543 (0.297), R 5,448, H 27,138, N 15,058; semiglobal V 4,239 / M 29 (**0.993**),
  R 2,259, H 1,501, N 10; harness 10.16 s vs 1.90 s summed, slowest family 2.27 s vs 0.72 s.

### Held-out result 1 (2026-10-12) and Amendment E2 + fresh held-out 3 (written before any chimp family is aligned)

**HELD-OUT 1 (gorilla, 269 families) — P1-P4 PASS:** global V 16,445 / M 22,848 (precision 0.419), R 7,090, H 63,784,
N 86,288; semiglobal V **18,317** / M **31** (**0.998**), R 2,187, H 2,939, N **119**; harness 34.53 s vs **5.15 s**
summed, 0 semiglobal families > 2 s (`eng/ho1_eval/engine_eval.json`).

**HELD-OUT 2 (human, 2,101 families) — P4 FAILS for semiglobal:** 10 families > 2 s (MCL1887 34.1 s: a 265-kb copy vs an
11-kb one; MCL1904 24.5 s). The GLOBAL engine could not be run as comparator: it aborted allocating 135 GB after 158 s
(whole set) and did not finish 2,099 families (without MCL1887 / MCL1904) within 590 s on 2 threads — semiglobal did all
2,101 in 115.7 s wall. **P1 on human (semiglobal alone, 7 chunks): V 377,664 / M 1,128 = 0.997 PASS** (R 49,458, H 89,694,
N 4,505; `eng/ho2_eval/engine_eval_semi_total.json`). P2 / P3 are not computable on human (no comparator output). **Verdict as registered: not all bars pass on both held-out substrates -> the engine stays
OPT-IN.**

**Amendment E2 (band; motivated by held-out 2, so judged only on held-out 3):** between anchors the diagonal envelope
spans every column across a huge one-sided insertion (250k rows x 11k columns). Alignment paths are monotone, so each
row's band is intersected with the columns between the flanking anchors (+- W): before the first anchor columns
<= its end, after the last >= its start. New unit test: a 30-kb one-sided insertion — banded score equals the full DP
and < 10 M cells.

**HELD-OUT 3 (fresh; another species, never aligned or scored):** the 393 chimp PTR families with 2-20 copies
(`rustle_figures/runs/chimp_PTR/chimp_PTR.fam.copies.tsv`). Bars P1-P4 exactly as above (P4: total semiglobal <= global
and no semiglobal family > 2 s); yardstick with A1. If the global comparator cannot finish a family, that family is
excluded from P2/P3 and listed, and counts against global in P4. **All of P1-P4 pass on held-out 3 (with held-out 1
already passed) -> the engine becomes the default; otherwise it stays opt-in.**

### Held-out result 3 (chimp PTR, 393 families) and P5 — 2026-10-12

| bar | global (current) | semiglobal (+E2) | verdict |
|---|---|---|---|
| P1 V/(V+M) | 0.422 (V 38,970 / M 53,395) | **0.993** (V 42,955 / M 296) | PASS |
| P2 V(semi) >= 0.99 V(global) | 38,970 | 42,955 | PASS |
| P3 N(semi) <= 0.10 N(global) | 120,013 | 766 | PASS |
| P4 total time | 155.83 s summed (392 families) + **MCL36 aborted** (13.8 GB allocation, 86-kb copy) | **7.71 s** (all 393) | PASS |
| P4 no family > 2 s | MCL9 105.8 s | **MCL9 3.64 s** (15 copies, 424 kb) | **FAIL (1 family)** |

R 15,753 / 5,391; H 142,063 / 12,081 (global / semi). MCL36 excluded from P2/P3 (comparator could not run it).
**P5 PASS:** `copy_assign` on DEV with `RUSTLE_PSV_ENGINE=semiglobal`: assignments and quant byte-identical; families'
PSV stats, PSV columns (33,015 -> 4,298; the RNA-editing pre-pass removes some) and copy_conversions (5 -> 0) change.

**Verdict as registered: P4's per-family 2-s bar fails on held-out 2 (10 families, pre-E2) and held-out 3 (1 family,
post-E2) -> the engine stays OPT-IN.** P1-P3 pass on every substrate where they are computable; total time is 4.5-20x
lower than the current engine, which itself aborts on 86-kb (chimp) and longer (human) copies. Adopting it as the
default despite the per-family bar is a decision for the user, not this prereg.

### Adoption (user decision, 2026-10-12)

The user chose to make the semiglobal engine the DEFAULT despite the per-family 2-s bar (P4) failing, and to cap the
previous engine's memory. Implemented: `psv_engine_from` (default semiglobal + guard; `RUSTLE_PSV_ENGINE=global` or a
legacy variable selects the global DP), `GLOBAL_DP_MAX_CELLS = 2^29` with pairs above it routed to the semiglobal
engine, `params.tsv` rows `psv_engine` / `psv_blocks`, `utilities psv-discover --engine` default semiglobal. Checks:
library 1,032 pass (2 new tests); full suite 1,257 pass / 3 fail (pre-existing `pipeline_driver_dispatch`, unrelated);
capped global on chimp: all 393 families in 37 s wall, 429 MB (MCL36 / MCL9 / MCL233 pairs routed); default on human
with E2: 66.4 s summed / 39.7 s wall (was 212.9 / 115.7), 2 families > 2 s (MCL1887 5.1 s, MCL1904 2.4 s);
`copy_assign` default on DEV: assignments and quant byte-identical, `psv_engine semiglobal`, columns = the opt-in run.


---

## PREREG_psv_significant_blocks_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PREREG_psv_significant_blocks_2026-10-12.md`)*

## PREREG — PSV columns only inside significant local blocks of the copy-vs-copy alignment (2026-10-12)

**Problem** (`docs/PSV_SUN_PILEUP_2026-10-12.md`; register 683). `discover_psvs` aligns each copy's spliced sequence to
copy 0's GLOBALLY (banded Gotoh; a 150-bp gap costs 182, 150 unrelated bases as mismatches ~112), so stretches with no
counterpart are forced into columns. Diagnostic `bench/psv_sun/gap_diag.py` on the 56-family KB3781 subset (HEAD
a37399a5 dump, 126,499 decisive (copy, column) pairs): 54.5% lie in stretches a chained `minimap2 asm20` alignment
supports ("aligned"); 31.6% in stretches asm20 misses but Smith-Waterman beats all 20 shuffles (divergent homologs:
`both_homologous` 21.3%, `no_alignment_homologous` 10.3%); 14.0% in stretches SW does not (`both_unrelated` 7.8%,
`terminal_overhang` 4.3%, `no_alignment_unrelated` 1.9%). Labels per pair: `/mnt/linuxdisk/tmp/psv_sun/gap_diag/pair_labels.tsv`
(frozen before the fix is written). Chained local alignment alone (register 577) would lose the 31.6%.

**Fix (decided before implementation).** Keep the global alignment; a column becomes a PSV column only inside a block
of the alignment path chosen by
  maximize  sum over kept blocks (S_block - T),   T = ln(L / alpha) / lambda
- column scores BLASTN-style: match +2, mismatch (or non-ACGT) -3, gap run of k columns -(5 + 2k);
- L = alignment length (columns), alpha = 1e-3 (the certificate's alpha), lambda solves
  0.25 e^{2 lambda} + 0.75 e^{-3 lambda} = 1 (ungapped Karlin-Altschul, uniform bases; K taken as 1);
- exact O(L) DP over prefix sums; columns outside kept blocks contribute no alignment (the copy is a gap there).
A block must pay its own significance cost, so an unrelated exon between two homologous blocks is cut out when its
negative score exceeds T, while a few PSV mismatches inside a block (a small dip) stay. Applies to every
`discover_psvs` engine (banded DP, poasta, minimap2).

**Acceptance (56-family subset, binary with vs without the fix, identical arguments):**
- **B1** assignments, quant, posterior/abundance byte-identical (the default read-star certificate does not read these
  columns). Otherwise stop and explain.
- **B2** pairs labelled `aligned`: >= 99% retained (same family, same copy-0 genome position, same copy, same allele).
- **B3** pairs labelled `both_unrelated` + `terminal_overhang` + `no_alignment_unrelated`: >= 90% removed.
- **B4** reported, no bar: retention of `both_homologous` and `no_alignment_homologous`.
- **B5** library tests pass (the 3 known `tests/pipeline_driver_dispatch.rs` failures are pre-existing and unrelated).
Fail of B2 or B3 -> the fix is not adopted as is.

### Result (2026-10-12) — B2 FAILED: not adopted as default; kept OPT-IN (`RUSTLE_PSV_BLOCKS=1`)

Run: binary with the mask vs a37399a5 behaviour, same 56-family subset and arguments (`/mnt/linuxdisk/tmp/psv_sun/
{run_fix2,run_blk}/`, scorer `bench/psv_sun/block_eval.py` -> `run_blk/block_eval.json`). PSV columns 33,015 -> 4,377.
- **B1 PASS**: assignments, quant byte-identical (copy_conversions: 5 -> 0 calls; 4 are MCL58, where each of its 4 copies
  is "converted" from the others — circular; not tested further).
- **B2 FAIL**: `aligned` pairs retained **0.257** (bar 0.99).
- **B3 PASS**: non-homologous pairs removed **0.9987** (both_unrelated 9,807/9,810, terminal_overhang 5,380/5,399,
  no_alignment_unrelated 2,386/2,386).
- **B4**: both_homologous retained 0.068, no_alignment_homologous 0.014.
- Library tests 1,023 pass (3 new: unrelated middle exon / unrelated terminal overhang -> no columns; 80%-identity
  homolog keeps all 80 columns).

**Why B2 failed (post hoc, not a re-scoring of the bar).** (1) The pair yardstick counts every copy with an allele at a
decisive column, including copies equal to copy 0, whose column disappears when the only differing copy was junk.
Restricted to pairs that DIFFER from copy 0: `aligned` retained 0.488. (2) The `aligned` label is about the STRETCH
(copy 0's offset inside a minimap2 block), not about O2's pairing: by the independent position verifier, 9,703 of
16,183 such differing pairs were the global DP's own MISALIGNMENTS (O2's column != minimap2's position) — retained
0.212 — while the 5,697 position-verified ones were retained **0.956**. So the mask mostly removes misalignment and
junk, but it also loses 4.4% of correct columns and cannot RECOVER the misaligned ones.

**Root cause exposed:** `banded_msa_pair` costs (mismatch 1, gap open 32, extend 1) make a ~30-base misalignment
cheaper than one gap, so inside homologous regions 60% of O2's apparent copy differences are misaligned columns.
The next step is the DP's scoring (match reward, log-odds mismatch, a cheap long-gap piece, free end gaps), with the
mask as a final guard — pre-registered with a position-verified yardstick on families held out from this subset.
Default behaviour: byte-identical to a37399a5 on all 7 outputs (`run_optin/`).


---

## PREREG_psv_sun_pileup_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PREREG_psv_sun_pileup_2026-10-12.md`)*

## PREREG — SUNs vs our PSVs (M1) and a per-copy pileup check of PSVs (M2), gorilla KB3781 (2026-10-12)

**Why.** The advisor (a) calls PSVs "SUNs" and (b) argues pileup-style variant evidence would tell real variants from
errors better than our alignment-derived PSVs. M1 measures how our O2 PSVs relate to SUNs; M2 measures whether a
per-copy pileup AFTER assignment finds places where this animal does not carry the reference PSV allele, with
KB3781's own other haplotype as truth. Descriptive study; bars below say what each outcome means.

### Substrate and run (fixed before any output is read)

- Reads: KB3781 fibroblast IsoSeq `fibroblasts/GCA_029281585.2_flnc_mm.bam`; genome `winloci_scratch/GGO.fasta`
  (mGorGor1 v2 primary, RefSeq names); catalog = the 2026-10-10 de novo run
  `/mnt/linuxdisk/tmp/rustle_figures/runs/gorilla_KB3781/` (`fam.copies.tsv`, `fam.copies.fa`, 326 families).
- **Subset S (56 families):** rows of that run's `assign.families.tsv` with `psv_cols >= 1`, `n_reads >= 20`,
  `n_copies <= 20` (chosen for resources — the full assign is 9.4 min / 5.3 GB; no outcome of this study is known).
- **Binary:** snapshot of `rustle_target_m2/release/copy_assign` (md5 `8b0026b9b82de0878743814be1b01297`, the binary
  of the 10-10 run), same arguments as that run's driver call plus `--dump-psv --threads 2`, families filtered to S,
  regions = the 10-10 `regions.txt` lines overlapping an S copy. Run under `tools/rlock.sh light` with a 4-GB
  virtual-memory cap (no interference with other sessions: no cargo build, own output dir `/mnt/linuxdisk/tmp/psv_sun/`).
- **Sanity S1:** the subset run's (read, family) rows agree with the 10-10 run on status and copy for >= 95% of shared
  rows (differences allowed only for reads shared with families outside S). Fail -> stop and diagnose.
- **Gate V (coordinates):** `psv_col_pos` gives one copy's genome position per column; the other copies' positions
  are derived by aligning their `fam.copies.fa` spliced sequences to that copy's (`minimap2 -cx asm20 --cs`) and
  mapping offsets through each copy's exons and strand. A derived position is used only if the genome base there
  (complemented on `-`) equals the dumped allele. Gate: >= 95% of (copy, column) pairs with a non-gap allele verify;
  else stop and fix before reading M1/M2.

### M1 — how many of our decisive PSVs are SUNs

- **Private allele:** at a decisive column, an allele carried by exactly one copy of the family (gaps excluded).
- **SUN** (our operational form of Sudmant 2010): a private allele with >= 1 overlapping genomic 30-mer (taken from
  that copy's genome position, canonical) that occurs exactly once in GGO.fasta. Counted with a purpose-built k-mer
  scan (candidate k-mers only; no genome-wide k-mer database — disk).
- **Reported (per family and totals, with n):** decisive columns; columns with >= 1 private allele; SUN alleles /
  private alleles; **G = private alleles with NO genome-unique 30-mer / private alleles** (the allele+context exists
  elsewhere in the genome = an outside-the-catalog sequence our catalog-relative PSV definition cannot see);
  **R = O2-assigned reads that carry >= 1 SUN allele of their assigned copy / O2-assigned reads** (what a SUN-only
  method could assign; n ~ 529 rows on 10-10, so R is indicative only).
- **Meaning:** G >= 0.05 -> recommend annotating O2's PSVs with genome-uniqueness (a catalog-private allele is often
  not genome-private); G < 0.05 -> the catalog-relative definition loses little. R is descriptive (how much the
  non-private PSVs add over SUN-only evidence).

### M2 — per-copy pileup check of PSVs

- For copy X and decisive column j (verified position g): the per-copy pileup = (a) reads PLACED on X by the aligner —
  primary alignments (`-F 2308`) with MAPQ >= 20 covering g (`samtools mpileup -q 20 -Q 0 --ff 2308 -B`, per copy
  window, bases in genome orientation; ref-skips and deletions not counted) — plus (b) reads O2 ASSIGNED to X in the
  dump (only 529 such rows in S on 10-10, hence (a); the default AS-tied gate never certifies uncontested reads).
  Pools are reported separately; flags use their sum.
- **Flagged:** >= 10 covering reads, a single non-X allele in >= 3 reads and >= 0.20 of them (thresholds of the
  RNA allele caller, register row 1194). **Control:** >= 10 covering reads, X's allele in >= 0.95.
- **Truth:** X's locus (genomic span +-2 kb) aligned to the OTHER haplotype's chromosome (`rna_allele/chrmap.tsv`
  `B_name`; the primary chromosome is identical to the same-name haplotype) with `minimap2 -cx asm5 --cs`; counterpart
  = one alignment covering >= 95% of the query. Classes: **B_ALT** (other haplotype's base = the reads' alt allele),
  **B_OTHER** (differs from X's reference allele, not the alt), **B_SAME** (= X's reference allele), **NO_COUNTERPART**.
  One minimap2 call per chromosome, each as a light job.
- **Bar (pileup check carries real signal):** among flagged pairs with a counterpart, B_ALT >= 0.50 AND the flagged
  B_ALT+B_OTHER rate >= 5x the control rate. Otherwise flags are mostly not the animal's own sequence (mis-assignment,
  editing, error).
- **Negative control:** copies on chrX (male, no second X: `B_name` empty) — reported flag rate per covered pair vs
  autosomes; het-driven flags cannot exist there.
- Also reported: transitions share among flags (A<->G editing; O2 already filters within-copy A<->G columns), and the
  share of O2-assigned reads crossing a flagged column.

### Not claimed

One animal, one tissue, 56 families; SUN = our operationalization (30-mer, primary assembly only); M2 truth tests only
the other haplotype (the primary IS one haplotype). Code `bench/psv_sun/`; outputs `/mnt/linuxdisk/tmp/psv_sun/`.

### Amendment 1 (2026-10-12, after gate V failed; before any M1/M2 output existed)

- **S1 PASS** (26,721 / 26,821 = 99.6% of the 10-10 rows identical). **Gate V FAIL**: 118,169 / 216,593 (54.6%)
  (copy, column) pairs verify (`gate_v.json`), after the mapper was fixed to use every alignment piece on either strand.
- **Diagnosis:** (i) 10 of 56 families contain copies whose spliced transcripts are REVERSE-COMPLEMENT homologs of the
  anchor's (inconsistent strand calls on 1-2-exon de novo copies). Their dumped alleles match the homologous base 26%
  and its complement 29% (chance level): O2's column frame pairs NON-HOMOLOGOUS bases there. Median decisive-column
  density 723 / kb vs 9.8 / kb in the 46 strand-consistent families. (ii) In strand-consistent families 37% of pairs
  are unmapped: no local alignment covers the anchor's position (non-homologous segments forced into the global frame,
  or divergence beyond `asm20`). Reported as a separate finding; not fixed here (O2 code is another session's).
- **Amended population:** a decisive column is **homologous** iff every copy carrying a non-gap allele has a VERIFIED
  position. M1 and M2 are computed over homologous columns only, and the share of O2's decisive columns that are
  homologous is reported (overall and by family strand-consistency). Gate V is replaced by this restriction; no other
  definition or bar changes.


---

## PREREG_s1_o2_certificate_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PREREG_s1_o2_certificate_2026-10-12.md`)*

## Are the secondary-only copies expressed? O2's certificate on their reads, 2026-10-12

**Pre-registered before the script exists and before any read is tested.** User (2026-10-12), after register 1268 (our good pool
assembles transcripts on secondary-only genes S1, where StringTie and FLAIR have none): "for the loci we get from secondary-only
alignments, does SQANTI prove that they are real?" Answer given: no — SQANTI checks structure, not which copy the reads came from;
a paralog's reads can rebuild an identically structured annotated transcript on a silent copy. The test is O2's certificate.

### Population

The S1 genes of `docs/PREREG_inbetween_operating_point_2026-10-11.md` (strata `/mnt/linuxdisk/tmp/asm_gap/op/<species>.strata`,
benchmark samples gorilla OR6737 / human A119b, development chromosome excluded). For each S1 gene: every GOOD secondary record
(AS >= 0.98 x the molecule's genome-wide best AS) whose aligned blocks overlap one of the gene's exons by >= 1 bp — the records the
good pool used there. Each is paired with the same molecule's PRIMARY record (its placement elsewhere, copy A; the S1 locus is copy B).

### The certificate (O2's, unchanged: `copy_pair_significance` as ported in `bench/asj_o2/asj_o2.py`)

Columns = read positions aligned (`=` or `X`) in BOTH records whose reference bases differ between the two placements (read
orientation; bases from the FASTA): the PSVs this read spans. `k_A` / `k_B` = columns where the read matches A's / B's base (`=` in
that record). Strict maximum (`k_A != k_B`), p = P(X >= k_best), X ~ Binomial(columns, 0.003 / 3); assigned when p < 1e-3 (two
copies), else abstain. Indel columns are not used. Every (secondary record, primary record) pair is tested once.

### Per-gene class (fixed now)

- **supported**: >= 2 records certified to the S1 copy (B) — the read floor of the assembler;
- **contradicted**: 0 records certified to B and >= 2 certified to A;
- **undecidable**: otherwise (mostly: every record abstains — identical copies or too few PSVs).

### What is reported (both species; nothing is adopted or tuned)

1. Records tested; certified to B / to A / abstained; split by TIED (AS == the molecule's best) vs STRICTLY WORSE secondaries.
2. S1 genes by class, all S1 genes and those where OP_good emits a transcript.
3. OP_good's transcripts touching S1 genes, by their gene's class, and how many match a reference transcript exactly (tmap `=`).

**Expected, written now.** Strictly worse secondaries: certified to A far more often than to B (the aligner already prefers A by
>= 1 AS point). Tied secondaries: mostly abstain (identical copies). So most S1 genes will be undecidable or contradicted, and only
"supported" S1 genes may be called expressed copies; that, not register 1268's 34 / 36 matched transcripts, is the claim we would
make for the advisor.

**Instrument** (committed before running): `bench/tool_difference/s1_o2.py` sha1 6b6831e392a2349c7fb5f51b6eb9113b08484117, tests `test_s1_o2.py` (7 pass).

### Result (2026-10-12): the secondary-only copies are NOT supported by O2 — their reads belong to the other copy

Outputs `/mnt/linuxdisk/tmp/asm_gap/op/{gorilla,human}.s1o2.{txt,records.tsv,genes.tsv}`. Instrument self-check (gorilla read
m64076_221110_210557/127338300/ccs): every `=` column's reference base equals the read base on both placements (0 / 1,827 and
0 / 1,457), so the columns and bases are right.

| | gorilla | human |
|---|---|---|
| good secondary records tested (S1 genes) | 14,698 | 18,613 |
| AS == molecule best ("tied"): to S1 copy / to primary copy / abstain | 3 / 7,531 / 414 | 2 / 3,883 / 5,947 (5,859 with no distinguishing column) |
| AS < best (strictly worse): to S1 copy / to primary copy / abstain | 16 / 6,259 / 475 | 24 / 7,292 / 1,465 |
| read matches the PRIMARY copy at every distinguishing column | 9,780 (66.5%) | 7,107 (38.2%) |
| read matches the S1 copy at every distinguishing column | 13 | 16 |
| mosaic (>= 2 columns for each copy) | 1,805 (12.3%) | 3,655 (19.6%) |
| S1 genes: supported / contradicted / undecidable | 1 / 400 / 131 | 1 / 352 / 419 |
| ... with an OP_good transcript | 1 / 151 / 49 | 1 / 202 / 134 |
| OP_good transcripts there (reference-exact): supported / contradicted / undecidable | 4 (2) / 193 (38) / 45 (16) | 1 (0) / 267 (44) / 185 (3) |

(The transcript counts here use every S1 gene; register 1268's 171 / 224 used the S1 genes with a multi-exon reference transcript.)
**Reading.** In both species exactly ONE secondary-only gene has >= 2 reads certified to its own copy. Most S1 genes are
contradicted (gorilla 400 / 532) or undecidable because the copies are identical over the reads (human 419 / 772, almost all with no
distinguishing column). The transcripts our good pool builds there — including most of the reference-exact ones — are SHADOWS of the
paralog the reads came from, rebuilt on a copy with the same exon structure. Register 1268's C1 ("wins on secondary-only loci")
stands as scored but does not mean those copies are expressed. Corrections to the wording of the earlier preregistrations: (1) "tied"
(AS == the molecule's best) is not "tied with the primary" — the primary's AS is below the best in 3.8% (gorilla) / 3.9% (human) of
multi-record molecules; true ties (primary == best == second) are 0.5% / 1.6%; (2) S2 ("coin-toss-starved") used the same "tied" and
deserves the same certificate check before it is read as O2's regime. Useful side result: O2's certificate separates the shadows
(contradicted) from the undecidable — the step that should prune good-pool transcripts, as the user's 10-11 rule puts it (a primary or
secondary is not trusted as the origin until O2 says so).


---

## PREREG_yag_o2_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PREREG_yag_o2_2026-10-12.md`)*

## PREREG — does O2 work on the Y ampliconic genes (YAGs)? (2026-10-12)

**Why.** The user: "YAGs are definitely an important test" for O2. The Y is HAPLOID, so every difference between two Y
copies is a PSV, never a heterozygous allele. Known before this prereg (10-06 audit, `rustle_eval_2026-10-06/03_yag.md`):
real human Y — O2 assigns 12 of 3,641 contested molecules; 10-06 simulation (human chrY roster, 0.13% error) — O2
assigns 0 of 2,425 MAPQ-0 reads, including the 117 reads whose copies are only 98-99% identical (that last fact is
unexplained: such copies should usually be distinguishable).

### Arm A — simulation with truth (human chrY; reuses the 10-06 products, nothing re-simulated)

Reads `o2_default_roster_2026-10-06/chrY/sim.bam` (10,224 reads from the 458-copy chrY roster, names = `FAM|copy|i`,
genome-mapped), families `asg.fam.copies.tsv` / `.fa`, MCL0 (the Yq12 satellite) skipped — the exact 10-06 command,
with the CURRENT binary and `--dump-star`. Per molecule O2 decided on (any status), from `star_reads.tsv`:
- **distinguishable** = the truth copy t is a candidate and every other candidate c has >= 2 columns where t's and c's
  aligned bases differ and the read covers both (O2's own m >= 2 step; register: "m = 1 never assigns").
- **A-P1 precision:** 0 wrong among assigned (upper 95% bound reported; UNDERPOWERED if < 30 assigned).
- **A-P2 recall:** of distinguishable molecules, >= 0.80 assigned to t.
- **A-P3 calibration:** of NOT distinguishable molecules, <= 0.05 assigned.
- Reported: end-to-end placement (O2's copy where assigned, else the aligner's primary) vs truth; by identity band.

### Arm B — real testis reads, split-half concordance (human A119b, CHM13 chrY)

`copy_assign --dump-star` on `A119b.t2t.bam` restricted to the A119b catalog's chrY families (the default pipeline's
own copies; satellite/array families with > 100 copies skipped as in Arm A). The read-star certificate uses PSV columns
only here (`read_star_junctions false`), with per-column error e/3 = 0.001: for the best copy b and each rival c,
p = P(Binomial(n_bc, 0.001) >= k_bc) (n = columns where b and c differ and the read covers both, k = those matching b);
Assigned iff every rival has p < alpha / (#rivals) and LLR > 0, Tied if some rival has n = 0, else Ambiguous; b = the
maximin pairwise-LLR candidate. Python reimplementation `bench/yag_o2/split_half.py`.
- **B-P0 reproduction:** on the full columns it reproduces O2's own status and copy for >= 0.99 of molecules whose final
  status was not changed by a later gate (origin rejection, tie outside the catalog). Fail -> stop, nothing else read.
- Halves: columns at read positions < L/2 vs >= L/2 (L = read length from the BAM) — disjoint evidence from disjoint
  parts of the same molecule; on a haploid Y no allele can make them disagree.
- **B-P1 split-half concordance:** among molecules where BOTH halves are Assigned, the same copy in >= 0.98 (>= 30 such
  molecules, else UNDERPOWERED).
- **B-P2:** where O2 (full read) and a half both assign, the same copy in >= 0.98.
- Reported: per family; the fraction of contested molecules where at least one half assigns.

### Side check (no bar, report only)

SUN calls vs an external SUN set (human: the QuicK-mer2 k-mer set reproduced in register row 1171), if its files are
on disk.

Gorilla (OR6737 testis) needs a families stage first; added as an amendment with the same bars only if run.
Code `bench/yag_o2/`; outputs `/mnt/linuxdisk/tmp/yag_o2/`.

### Amendments made while running (2026-10-12)

- **Dump bug (found, fixed, TDD):** `--dump-star` kept ONE proof per read name, so a molecule evaluated in two families
  got family B's columns next to family A's verdict. First-pass numbers built on it are RETRACTED (53 "distinguishable"
  simulated reads, a "17 wrong" re-certification). Fix: every registration kept with its candidates' tids; the writer
  picks each row's own family proof, once per (molecule, family). Assignments byte-identical to the unfixed binary.
- **Scope (user, mid-run): "we only need this for the YAGs, so DAZ, PRY, etc."** Families = catalog families with a
  copy overlapping a RefSeq CHM13 chrY gene DAZ / RBMY / TSPY / VCY / XKRY / CDY / HSFY / BPY2 / PRY / PRORY
  (`bench/yag_o2/prep_yag.py`): A119b catalog 11 families / 86 copies (TSPY 46, RBMY 6, HSFY 6+3+2, VCY 6, DAZ 4+2,
  PRY 4, CDY 4, BPY2 3); simulation roster 11 families / 88 copies. (A first Y-only run is not reported as a result.)

### Results

**Arm A (simulation, YAG families):** of 1,460 contested molecules whose true copy is a candidate, 1,416 have a rival
IDENTICAL over every column the read covers and 42 differ at one column; 2 are distinguishable.
- A-P1 precision: 0 assigned -> UNDERPOWERED. A-P2 recall: 0 / 2 -> FAIL as registered (n = 2). **A-P3 calibration:
  0 / 1,458 undecidable assigned -> PASS.** Reproduction of O2 by the Python certificate: 2,404 / 2,404.

**Arm B (A119b testis reads, YAG families):** 7,844 molecules; O2 assigns 4. **B-P0 reproduction 3,445 / 3,445 = 1.0
PASS.** **B-P1 split-half concordance 31 / 57 = 0.544 FAIL** (HSFY 11/23, BPY2 16/22, PRY 2/9, other HSFY 2/3; DAZ and
TSPY: no half assigns). B-P2 3 / 3 -> UNDERPOWERED.

**Reading (post hoc; not tested here).** The haploid argument assumed the reads come from the reference man's Y. They
do not: CHM13's Y is HG002's, the testis RNA is A119b's. Real reads carry MORE distinguishing columns than reads
simulated from the reference copies (386 molecules with an assignable half vs 37), as expected if A119b's YAG copies
differ from HG002's (gene conversion between palindrome arms differs between men), so a read's halves can point to
different REFERENCE copies with no error by O2. O2's full-read certificate abstains on all of those. The split-half
test therefore cannot separate O2 error from donor-vs-reference mosaicism without a matched donor: a T2T Y and testis
Iso-Seq from the SAME man (YAGs are testis-specific; the T2T Ys come from cell lines). The defensible statement for
the YAGs is calibration, not accuracy.


---

## PSV_SUN_PILEUP_2026-10-12.md

*(original: `docs/archive/2026-10/merged/PSV_SUN_PILEUP_2026-10-12.md`)*

## SUNs vs our PSVs (M1) + per-copy pileup check (M2) — results, gorilla KB3781 (2026-10-12)

Prereg `docs/PREREG_psv_sun_pileup_2026-10-12.md` (3be23dc3) + Amendment 1 (b658a58e, before any M1/M2 output).
Code `bench/psv_sun/` (6 unit tests); outputs `/mnt/linuxdisk/tmp/psv_sun/` (`sanity.json`, `gate_v.json`, `m1.json`,
`m1_families.tsv`, `m2_pairs.tsv`, `m2.json`, `truth/`). Subset S = 56 families / 208 copies of the 10-10 de novo catalog.
copy_assign snapshot (md5 8b0026b9) `--dump-psv`, 290 s / 9.5 GB under the heavy lock; every other step a light job.

### Finding 0 (gate V, unplanned) — many of O2's PSV columns are not homologous positions

- S1 PASS: 99.6% of the 10-10 rows reproduced.
- Gate V FAIL: only 54.6% of 216,593 (copy, column) pairs verify by independent local alignment (minimap2 asm20,
  either strand, all pieces) + genome base check.
- **10 / 56 families contain copies whose spliced transcripts are reverse-complement homologs** (strand calls on
  1-2-exon de novo copies). O2's own dump shows it: at O2's columns these copy pairs agree at median 0.263 (20 pairs,
  chance level) vs 0.508 for same-strand pairs (95 pairs); decisive-column density 723 / kb vs 9.8 / kb. O2's global
  copy-vs-copy alignment pairs NON-HOMOLOGOUS bases there, so those "PSVs" are alignment artifacts (they inflate
  n_decisive and certificate confidence).
- In strand-consistent families 37% of pairs are not covered by any local alignment (forced global alignment of
  non-homologous segments, or divergence beyond asm20 — not separated here).
- Only **5,235 of 47,134 (11.1%)** decisive columns are homologous in the strict sense (every carrying copy
  verified); strictness is part of that number (one unverified copy excludes the column).

### M1 — PSVs vs SUNs (homologous columns)

| quantity | value |
|---|---|
| homologous decisive columns | 5,235 |
| ... with >= 1 copy-private allele | 5,069 |
| ... with >= 1 SUN (private allele with a genome-unique 30-mer) | 1,594 (30% of homologous PSV columns) |
| private alleles / SUNs | 6,122 / 2,000 (33%) |
| G = private alleles with no genome-unique 30-mer | pooled 0.673; **per-family median 0.125** (43 families with >= 5 private) |

- G is concentrated: MCL116 (1,131/1,162), MCL316 (860/860), MCL38 (654/678), MCL90 (508/609) = 76% of the
  non-unique ones; mostly NOT soft-masked (3,471 unmasked vs 651 masked) — the allele + context exists elsewhere in
  the genome outside the expressed catalog (unexpressed copies the de novo catalog cannot contain).
- Prereg reading: G >= 0.05 -> recommend annotating O2's PSVs with genome-uniqueness. Holds per family (median .125)
  and pooled.
- **R not measurable:** none of the 529 O2-assigned rows covers a homologous column (reported as undefined, not 0).

### M2 — per-copy pileup at PSV columns vs KB3781's other haplotype

- 18,168 verified (copy, column) pairs; 5,370 covered by >= 10 reads (MAPQ>=20 primaries; O2-assigned reads add
  none here). **66 flagged**, 5,270 control.
- Truth: flagged 61/62 with a counterpart = B_ALT (**0.984**); control B differs 10/4,734 (**0.0021**); fold 469.
  **Bar PASS.** NO_COUNTERPART 4 flagged / 456 control.
- 58 of the 61 confirmed flags: the alternative base IS ANOTHER COPY'S PSV ALLELE (3 novel) — on one haplotype copy X
  carries copy Y's allele (median alt fraction 0.47 = heterozygous). Reads from that haplotype give O2 evidence for
  the wrong copy at that column. 30 copies in 17 families.
- No-second-haplotype chromosomes (chrX 4 pairs, chrY 2, mito 74): 0 flags / 80 (chrX alone underpowered).
- Flag transitions 68% (germline-like Ts/Tv; the genome confirms them, so not editing).

### Answer to the advisor

- Pileups do not replace alignment for FINDING PSVs, but a per-copy pileup after placement is a real, cheap check:
  1.2% of covered PSV columns in this animal are heterozygous (one haplotype carries another copy's allele), and the
  check is right 98% of the time against the animal's own haplotype assembly.
- SUN != old name for PSV: SUNs are the genome-unique, copy-private subset (~30% of homologous PSV columns here).
- The bigger problem the check exposed is upstream: O2's PSV columns include non-homologous artifacts (reverse-strand
  copies, forced global alignment). Proposed register rows (not added — register shared with another session):
  (1) "O2 PSV columns are homologous positions" — 54.6% of pairs verify; 10/56 families reverse-complement artifacts.
  (2) "A per-copy pileup flags real deviations from the reference PSV allele" — 61/62 vs 0.2% control, PASS.

### Fix applied (2026-10-12, same day) and corrected numbers

**Fix** (`src/family.rs`): `discover_psvs_oriented` orients each copy against copy 0 before the copy-vs-copy alignment
(`is_reverse_homolog`: 15-mers shared forward vs reverse-complemented; flip iff more reverse AND sign test p < 1e-3,
the certificate's alpha). A flipped copy is aligned reverse-complemented, its column entry is (its own genomic position,
its base in copy 0's frame), `FamilyProfiles::strand` / `FamilyDetail::copy_frame_strand` / the dump's `copy_strand`
carry its FRAME strand so reads and alleles stay comparable. Opt-in intron columns follow the same frame. 4 new tests
(TDD, each watched failing first: 240 artifact columns for a 300-bp reverse pair; a 20-bp inverted stretch flipped by a
bare majority). Library suite 1,020 passed / 0 failed; full `cargo test --no-fail-fast` 1,244+ passed, 3 FAILED =
`tests/pipeline_driver_dispatch.rs` (pre-existing and unrelated: the driver's `all)` now calls `run_all`, which the
test's stubs do not define).

**Before/after on the same subset** (HEAD 3e0f3fcd built without vs with the fix, identical arguments):

| output | change |
|---|---|
| assignments / quant / copy_conversions | **byte-identical** (default certificate = read-star; it never read these columns) |
| strand-consistent families (46) | PSV columns **byte-identical** (0 changed) |
| reverse-pair families (10) | median 723 -> **10.0** columns/kb (e.g. MCL116 3,586 -> 5, MCL90 1,864 -> 19, MCL316 1,192 -> 0) |
| PSV columns total | 47,134 -> 33,015 |
| reverse-pair positions verified | 26% -> **77%** (9,233 / 12,030) |
| gate V overall | 54.6% -> 56.2% (the rest = register 683: different exon composition; `--psv-genomic` is its opt-in remedy) |

A first version flipped on a bare majority and changed 2 strand-consistent families (MCL55 1 vs 3 k-mers, MCL59 60 vs
61); the sign test removed both.

**Corrected M1/M2 (fixed dump; the first-pass numbers above were inflated by artifact columns that verified by chance):**

- M1: homologous decisive columns 1,799; SUN columns **1,409 = 78%** (was 30%); private alleles 2,116, SUN 1,815 = 86%;
  **G pooled 0.142** (was 0.673; MCL116/316/38/90 were artifacts), per-family median 0.10 (43 families). G >= 0.05 still
  -> genome-uniqueness annotation recommended. R still undefined.
- M2: 7,246 pairs, 4,711 covered; **258 flags**, 4,419 control. Flags B_ALT **222/224 = 0.991**; control differs
  10/4,111 = 0.0024; fold 415; **bar PASS**. Flags concentrate in the repaired families (MCL8 152, MCL47 36): their real
  PSV columns now exist and many are heterozygous in KB3781. Outputs `/mnt/linuxdisk/tmp/psv_sun/{run_head,run_fix2,eval_fix2}/`.


---

## PREREG_o2_anchor_2026-10-13.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_anchor_2026-10-13.md`)*

## PREREG — does O2 recover the right copy? Truth from the read's own unique anchor (2026-10-13)

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

### Arm 1 — anchor trim (gorilla = dev; human = held-out)

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

### Results — gorilla (dev), 2026-10-13

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

### Amendment 1 (2026-10-13, after the gorilla dev run, before the human held-out run)

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

### Results — human A119b (held-out)

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

### Amendment 2 (2026-10-13, before it runs): arm 2, junction truth from rival copies that cannot splice

Both arms run on the final binary: exact-tie fix 7ac557e7 plus the tie-outside gate fix d5d982a3, default flags plus
`--dump-star`.

**Idea.** A spliced read cannot come from a copy whose genome cannot produce its splice junction. The evidence sits in
the INTRON (the copy's splice-site dinucleotides), which the read does not contain. O2 decides on exonic PSV columns
only (`read_star_junctions` is off by default), so the two kinds of evidence are disjoint.

A parent gene against a processed retrocopy is NOT used as truth. The retrocopy's own transcript has the parent's
spliced structure, so a spliced read does not exclude it.

**Population.**
- gorilla_KB3781 (dev): the full assign with the final binary.
- Human A119b (held-out): the anchor run's 129 families, original molecules only, no `|trim`.
- Molecules: contested rows (`contested` = 1) with a read-star row, whose primary alignment has >= 1 intron (`N` >= 30 bp).

**Per molecule and candidate (`bench/o2_truth/junction_truth.py`).**
- The read, in its original orientation, is aligned to the candidate copy's genomic span +- 1 kb
  (`minimap2 -c -x splice:hq -uf`, secondaries kept); the best hit by AS is used.
- For each read junction q (read position of an `N` in the genome primary alignment):
  - the candidate **breaks** q if its alignment has an intron within 5 bp of q whose motif, in transcript orientation,
    is not GT-AG, GC-AG or AT-AC;
  - it **makes** q if that intron is canonical;
  - otherwise q is undetermined for it (no hit, unspliced there, or the intron sits elsewhere).
- A candidate is **EXCLUDED** iff it breaks >= 1 read junction, and **PRODUCES** iff it makes every read junction.
- **Informative:** >= 1 candidate EXCLUDED and >= 1 not.
- **Determined:** the non-excluded candidates form one locus (their spans overlap) and it PRODUCES; truth = that locus.

**Bars** (each on the dev and the held-out set separately):
- **J1 (falsification).** Among O2-assigned informative molecules, wrong = O2's copy is EXCLUDED. PASS iff
  wrong / assigned <= 0.01 with >= 30 assigned. UNDERPOWERED if < 30.
- **J2 (determined truth).** Among determined molecules that O2 assigns, PASS iff wrong / assigned <= 0.01 with >= 30.
- **Reported:** the aligner primary's rate of sitting on an EXCLUDED copy (informative molecules); recall on
  determined molecules; per family; every wrong call and up to 20 right ones, with columns and motifs.
- **Caveat:** about 1% of real introns are non-canonical. A wrong call against an EXCLUDED copy is examined
  individually, and the motif is reported.

### Amendment 3 (2026-10-13, before it runs): arm 3, split-half on the reference animal itself

**Why.** The YAG split-half (0.544) used reads from a man who is not the reference. KB3781's fibroblast reads come from
the animal whose genome IS the reference, so that donor confound is absent.

**Run.**
- The gorilla full assign star dump (as in Amendment 2), contested rows, all families. Read lengths come from the BAM.
- `bench/yag_o2/split_half.py`, unmasked (all read-star columns).
- `split_half.py` with `POPULATION_ANY` on a gorilla catalog-PSV site table: every `(chrom, genome_pos)` of
  `utilities psv-discover` (default engine) on the gorilla_KB3781 catalog.

**Bars.**
- **Gate P0.** The Python certificate reproduces Rust's status and copy (unmasked, rows no later gate changed) for
  >= 0.99. Fail: stop.
- **S1.** Split-half concordance, all columns: >= 0.98 among molecules whose two halves both assign, with >= 30 such
  molecules (else UNDERPOWERED).
- **S2.** The same on catalog-PSV columns.
- **Reported:** the fractions with >= 1 assigning half, per family, and every disagreeing molecule with its columns.

**Reading, fixed now.**
- **S1 PASS:** O2's evidence is self-consistent where the genome is the animal's own. The YAG failure was specific
  to the donor or to artefacts.
- **S1 FAIL and S2 PASS:** non-PSV read-star columns (alignment artefacts) drive the disagreement, as on the YAGs.
  Restricting O2 to catalog-PSV columns becomes a candidate change.
- **Both FAIL:** mosaic molecules (diploid alleles or gene conversion), or O2 error. Telling them apart needs haplotype
  resolution (mat/pat), which is not attempted here.
- **Caveat:** the primary assembly carries one haplotype per copy. Reads from the other haplotype can legitimately
  match different copies at different PSVs.

### Results — arm 2 (junction truth), 2026-10-13

Final binary. `bench/o2_truth/junction_truth.py` (8 tests, committed before scoring). Outputs `/mnt/linuxdisk/tmp/o2_arms/`.

| | gorilla (dev) | human (held-out) |
|---|---|---|
| spliced contested molecules scored | 34,161 | 14,437 |
| informative (>= 1 candidate EXCLUDED, >= 1 not) | 90 | 1,395 |
| **J1** O2-assigned informative / wrong (O2's copy EXCLUDED) | 4 / 2: UNDERPOWERED | 80 / 51: **FAIL as registered** |
| determined / **J2** O2-assigned (right) | 3 / 2 (1): UNDERPOWERED | 47 / 1 (1): UNDERPOWERED |
| aligner primary on an EXCLUDED copy (informative) | 26 / 90 | 458 / 1,395 |
| aligner right on determined | 1 / 3 | 9 / 47 |

**The examination the prereg required: the truth is invalid.**
- **The motifs are shifted canonical sites.** The exclusions behind the "wrong" calls are one-base shifts of
  canonical sites: human GT-TA 21, AT-AG 10, others <= 3; gorilla GT-AA 2.
- **They sit at the read ends.** 51 of the 57 human exclusions are at the read's first or last intron. That is
  minimap2 placing terminal junctions imperfectly, not a copy that cannot splice.
- **The truth fails its own sanity checks.** The aligner, which scores splice signals, sits on an "excluded" copy
  for 29-33% of informative molecules. On "determined" molecules it is right only 19% of the time.
- **So J1's FAIL is not evidence against O2.**
- **Arm 2 is uninformative by construction.** Copies that differ in splice-site canonicity get different alignment
  scores, so their reads rarely become the exact ties O2 handles: 0.26% of gorilla's spliced contested molecules are
  informative.

A sound version (internal introns only; no canonical pair within +-3 bp in the copy's genome) would need a fresh
substrate.

### Results — arm 3 (split-half on the reference animal), 2026-10-13

Gorilla full assign, final binary (`/mnt/linuxdisk/tmp/o2_arms/gorilla/`). Catalog-PSV table: 72,705 sites
(`psv-discover`, 326 families).

- **Gate P0 PASS:** 18,632 / 18,636 = 0.9998.
- **S1 FAIL:** all columns, 107 / 477 = 0.224 (6,566 molecules with >= 1 assigning half).
- **S2 FAIL:** catalog-PSV columns, 37 / 382 = 0.097 (6,061 / 59,423 = 10.2% with >= 1 assigning half).

**Reading as registered: both FAIL**, meaning mosaic molecules or O2 error. Telling them apart needs haplotype
resolution.

**Descriptive.**
- **Disagreement is concentrated:** MCL90 117 and MCL274 76 of the 370 disagreeing molecules. The two copies usually
  sit on the same chromosome within 1 Mb.
- **O2's full-read certificate abstains on 342 of the 370.**
- **On O2's own 135 contested calls:**
  - 51 have one decisive half, which agrees with O2 in 50;
  - 32 have both halves decisive: 4 agree and 28 disagree.
- **The 28 that disagree:** 27 are MCL26 with one pattern. Six 5' columns favour copy 1 (O2's call); the 3' half
  favours copy 0 by 2 columns to 1.
- **Hypothesis (not tested above):** an allele of copy 1 carrying gene-converted, copy-0-like sites. Checked next
  against the animal's own haplotype assemblies (mat / pat).

**Haplotype diagnostic of arm 3's 28 conflicting O2 calls (post hoc, descriptive).** The 28 reads and the MCL26 / MCL59
copies were aligned to mGorGor1 mat and pat (`splice:hq -uf`).
- **The primary assembly mixes haplotypes:** MCL26 copy 1 (the 3-exon parent) is mat chr17; copy 0 (single exon, a
  retrocopy) is pat chr7.
- **The parent's alleles differ:** its pat allele is 18 edits from its mat allele and 13 from the retrocopy.
- **All 27 MCL26 reads hit the parent's locus best,** never the retrocopy: pat 19 at NM <= 1, 4 at NM 2-5,
  4 at NM > 5; mat all NM >= 2.
- **So O2's 27 calls (copy 1) are right.** The split-half "conflict" is the parent's pat allele carrying
  retrocopy-like alleles at the 3' sites (allelic gene conversion).
- **MCL59 (1 read):** unresolved (NM 3 on mat; tied on pat).

### Amendment 4 (2026-10-13, before it runs): arm 4, truth from the animal's own haplotypes

**Idea.** A KB3781 read must match its true locus on mat or pat up to sequencing error. Each catalog copy (primary
coordinates) maps to its own locus on one haplotype at NM 0, and to its homolog on the other haplotype.

**Population.**
- **For precision:** every contested O2-assigned molecule of the gorilla final-binary run (135).
- **For recall and calibration:** a random 3,000 contested molecules (seed 7 on the md5 of the name).

**Per molecule (`bench/o2_truth/haplotype_truth.py`).**
- **Alignment:** `minimap2 -ax splice:hq -uf --eqx -N 10` to `mGorGor1.mat` and `.pat`.
- **Best hit:** the minimum NM across both haplotypes, accepted iff NM <= max(2, 0.003 x aligned bases), and no hit
  at a different locus is within NM + 2.
- **Copy loci:** each catalog copy's spliced sequence is aligned the same way, and its loci on each haplotype are kept
  (all hits within NM + 2 of its best there).
- **Truth:** the catalog copies whose loci overlap the best hit. With one copy (or one overlapping set of nested
  units), the truth is that copy. With none, the truth is "outside the catalog". With several unrelated copies, the
  truth is ambiguous and the molecule is not scored.

**Bars.**
- **H1 (precision):** among O2-assigned molecules with a one-copy truth, wrong <= 0.01, with >= 30.
- **H2 (reported):** the aligner primary's copy on the same molecules, right or wrong (O2 against the aligner).
- **H3 (reported):** on the 3,000-sample molecules with a one-copy truth: O2's recall, and O2's assigned-wrong count.
  Molecules whose truth is "outside the catalog" are counted (the O3 regime).

**Caveat.** This is single-substrate: no other animal has both a T2T diploid assembly and RNA from the same
individual here.

**Amendment 4a (before it runs; code review of the scorer).** "No hit at a different locus within NM + 2" would call
every read ambiguous whose two alleles both match: its mat and pat loci are different loci of the SAME copy. Uniqueness
is therefore judged at the copy level:
- **Near hits:** the read's hits within NM + 2 of its best are kept.
- **Truth:** the copies whose loci overlap any near hit.
- **Ambiguous:** only when those copies are not one locus in primary coordinates.

The best hit must still satisfy NM <= max(2, 0.003 x aligned bases).

### Results — arm 4 (haplotype truth), 2026-10-13

`bench/o2_truth/haplotype_truth.py` (5 tests, committed before running). Outputs `/mnt/linuxdisk/tmp/o2_arms/arm4_gorilla/`.
4,825 contested rows considered (135 O2-assigned molecules + a 3,000-molecule sample).

**Truth coverage.**
- **One-locus truth:** 29 assigned rows and 14 abstained rows.
- **Ambiguous truth:** 64 assigned and 2,985 abstained. Near-identical copies tie on the haplotypes too.
- **Outside the catalog:** 2 assigned and 137 abstained (the O3 regime).
- **No near-perfect hit:** 775. **No primary record in the windows:** 819.

**H1: 28 / 29 right, 1 wrong: UNDERPOWERED** (one short of 30). The aligner's primary is right on 28 / 29 (H2).
- **The right calls are dominated by MCL26** (all 20 stored examples). Precision rests largely on one family.
- **The wrong call** is O2's MCL3 copy 12, while the read matches MCL132 copy 1 (another family) at NM 4 on pat. This is
  the known per-family blind spot: the certificate never sees a copy of another family. The aligner was wrong too.
- **The 2 assigned rows with truth outside the catalog** are one molecule (`SRR27438212.7269840`), assigned in MCL265
  and in MCL267 by sole-candidate promotion (0 decisive columns). Its true locus is an uncatalogued mat chr16 locus
  at NM 0, so both calls are wrong. Counting them: 28 right / 31 verifiable calls.

**H3 (sample):** 19 molecules have a one-locus truth. O2 assigns 5 (all right; recall 0.26); the aligner is right on 15.

**Reading.**
- Where the animal's own haplotypes give a truth, O2's PSV calls are right: 28 / 29 on one-locus truth.
- The failures are the two documented blind spots, not PSV evidence gone wrong:
  - the true copy in another family (per-family certificate);
  - an uncatalogued origin promoted by the sole-candidate rule.
- The test is underpowered and concentrated in MCL26. Haplotype truth resolves few tied molecules, because exact ties
  are mostly between copies that are near-identical on both haplotypes.


---

## PREREG_o2_exact_tie_2026-10-13.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_exact_tie_2026-10-13.md`)*

## PREREG — O2 exact ties decided by rounding noise: the fix (2026-10-13)

**Found by** the anchor-trim truth test (`docs/PREREG_o2_anchor_2026-10-13.md`, human held-out). Three MCL1623
molecules match both candidates at exactly 5 of 10 PSV columns. They were `assigned`, with margin 0.000 and p 2.5e-13,
and all three were wrong.

**Root cause.** Every O2 certificate decides on log-likelihood differences built as floating-point sums in column order:
- the read-star pair LLR, which drives both the maximin best-copy pick and the "LLR > 0" assign rule;
- the population `fixed_only_certificate`;
- the legacy per-read `read_copy_evidence` logl, which drives the argmax and the `margin > 0` rule.

When two copies' evidence is exactly balanced, the true difference is 0. The computed difference can be off by one
rounding error: the order `+++--+--+-` of ±lr sums to 1.78e-15. So "> 0" and the best-copy pick are decided by rounding
noise. The design itself says exact twins must fall through to the tie key, in the read-star maximin comment: "twins
(0 both ways) fall through to the old key".

In the shipped runs, the assigned contested calls with margin 0.000 are: gorilla_KB3781 7 of 140, chimp_PTR 2 of 46.

**Fix.**
- **Read-star certificate (and the maximin pick) and `fixed_only_certificate`:** the pair LLR is computed from integer
  counts, LLR(a, b) = (k - m) x lr + (kj - mj) x lrj. Here k / m are the decisive columns matching a / b, and kj / mj
  the junction terms for and against a. It is exactly 0 iff the evidence is balanced.
- **Legacy `read_copy_evidence`:** the decision uses each copy's terms summed in a canonical order (sorted), so copies
  with the same multiset of terms get bit-equal totals. The printed numbers keep the old sums.

No statistical rule changes: alpha, e/3, lr, lrj and the tie keys are untouched.

**Bars.**
- **F1 (unit, exactness).** An exactly balanced read is not Assigned and has margin exactly 0.0, in every column order.
  This holds for the read-star pair LLR helper, `fixed_only_certificate` and the legacy `assign_read`. Each test must
  fail on the unfixed code (red first).
- **F2 (real data, gorilla_KB3781).** The full assign stage is run with the before binary (`copy_assign.before`, md5
  15ba26eb, HEAD b84ddf84) and the fixed binary, same command (the pipeline's), and `assignments.tsv` is compared
  row-for-row.
  - Every row whose status changes must have had an old margin with |margin| < 5e-4 (printed `0.000`), and must move
    from `assigned` to `ambiguous` or `tied`.
  - Expected: those 7 rows.
  - Rows whose `assigned_copy` changes without a status change are reported. These are tie-key picks for exact twins.
  - Any other difference is a FAIL: stop and explain.
- **F3 (the anchor test).** The human arm-1 batches are re-run with the fixed binary and re-scored. The 3 MCL1623 wrong
  calls must no longer be `assigned`. A1 / A2 are reported, descriptive only (A1 stays UNDERPOWERED).
- **F4.** `cargo test --release --no-fail-fast` passes, except the 3 known `pipeline_driver_dispatch` failures. Golden
  tests that pin bit-level sums are updated only if their decisions are unchanged, and each such update is listed.

Outputs `/mnt/linuxdisk/tmp/o2_tie/`.

### Results (2026-10-13; fix 7ac557e7; outputs `/mnt/linuxdisk/tmp/o2_tie/`)

**F1 PASS.** The three tests were red first: legacy `assign_read` assigned the 5:5 tie, and `fixed_only_certificate`
had margin 1.78e-15. Now green; lib 1,047 / 0.

**F2 (gorilla_KB3781 full assign, before vs after): FAIL as registered, explained.** 60,214 rows; 60,189 identical.
- **4 status changes,** all assigned -> ambiguous, with old margin `0.000`, 8-10 decisive columns, p 7e-11 to 2.5e-13,
  two candidates. These are the fix.
- **The registered expectation of 7 was wrong.** The other 3 zero-margin assigned rows are sole-candidate calls
  (n_candidates 1, 0 decisive columns), whose margin is 0 by design. The fix correctly leaves them.
- **3 rows changed only the printed margin,** `-0.000` -> `0.000`. These are exact ties whose old rounding residue
  was negative. This is the same mechanism, with no status or copy change, but the bar did not foresee it, hence the
  FAIL.
- **18 rows changed only the copy:** exact twins now fall through to the tie key, as the prereg allows.
- `families.tsv` and `quant.tsv` change only through those rows.
- Nothing outside the exact-tie mechanism.

**F3 PASS (human anchor test, held-out).** The 3 MCL1623 wrong calls are now `ambiguous`.
- Anchor-trimmed reads: A1 (amended) 1 assigned, 1 right, 0 wrong (UNDERPOWERED). A2 0 / 42.
- Whole human run (129 families, 44,491 rows): assigned 804 -> 675; contested assigned **691 -> 562 (-19%)**.
  - 127 tie-noise calls removed;
  - 185 exact-twin copy changes;
  - 82 margin-sign rows;
  - 5 other status changes:
    - 3 tied -> ambiguous: exact twins; the tie key picks another best copy and both statuses abstain;
    - 2 assigned -> tied (MCL23, real margins 6.9 / 13.8): **a separate, pre-existing race**, explained below.

**F4 PASS.** `cargo test --release --no-fail-fast`: 1,286 passed. The only failures are the 3 known
`pipeline_driver_dispatch` ones.

**Second bug found (not fixed here): the tie-outside gate depends on region order.**
- **The mechanism.** A molecule with an AS-tied placement outside every catalog unit must never be Assigned. The
  registry of such molecules is filled while each region is swept, but read when a family's own region is
  certified. If the region holding the outside placement is swept later, the gate misses and the call stays
  `assigned`. The write-time column `tie_outside_catalog` still says 1.
- **The proof.** The before binary on 1 region thread gives `tied` for both MCL23 rows; on 4 threads it gave
  `assigned`. So "`--region-threads` is byte-identical by construction" is false here.
- **In the shipped runs,** assigned rows with `tie_outside_catalog` = 1: gorilla_KB3781 1 / 624, chimp_PTR 1 / 54, and
  the human anchor run 1 / 675.
- **Fix candidate:** apply the demotion after all sweeps, at write time, where the column is computed.

**Reading.** O2's "LLR > 0" rule was decided by rounding noise on exactly balanced molecules. In human testis these
are common: mosaic molecules matching both copies at several columns each. They made up about a fifth of the
contested assigned calls in these families. Every human O2 number from before 7ac557e7 overstates assignments by that
much.


---

## PREREG_o2_tie_outside_gate_2026-10-13.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_tie_outside_gate_2026-10-13.md`)*

## PREREG — the tie-outside gate depends on region sweep order: the fix (2026-10-13)

**Found by** the exact-tie fix's held-out check (`docs/PREREG_o2_exact_tie_2026-10-13.md`, F3; register 1274). Two
human MCL23 molecules were `assigned` by the binary before the fix and `tied` after it, with an unchanged copy and real
margins (6.9 / 13.8). The same before binary on 1 region thread also gives `tied`.

**Root cause.** The §6gz gate says a molecule with an AS-tied placement outside every supplied family unit is never
`Assigned`.
- The registry of such molecules (`register_tie_outside`) is filled by each region's worker as it sweeps its records.
- The gate reads that registry inside the worker that certifies the molecule's family (family.rs, the read-star verdict).

So when the region holding the outside placement is swept after the family's region (serially in region order, or
concurrently under `--region-threads`), the gate misses. The call stays `assigned`, while the write-time column
`tie_outside_catalog` says 1. Which calls leak depends on order and thread timing. Assigned rows with
`tie_outside_catalog` = 1 in the shipped runs: gorilla_KB3781 1 / 624, chimp_PTR 1 / 54; human anchor run 1 / 675.

**Fix.**
- The same demotion as the worker's gate (`tie_row`: `Tied`, not resolvable, uniform posterior), plus the
  `families.tsv` recount (`recount_after_union`), is re-applied to every region's families once ALL regions have been
  computed: `denovo_pipeline::apply_tie_outside_gate`, called in `copy_assign` before the cross-family pass and the
  drain.
- Every emitted status therefore sees the complete registry: assignments, posterior, psv_reads, phased reads, the
  placement fallback.
- Quant is unaffected: abundance comes from the per-family EM over PSV observations, and hard counts ignore status.
- Skipped under `--union-certificate`, whose pass scores outside placements itself.

**Not fixed (documented residual).** With `--families --gtf`, GTF routing by O2 status runs inside the worker. It can
still see a leaked `Assigned`. The pipeline never uses this combination: its assemble stage is `--assemble-only`.

**Bars.**
- **G1 (unit).** `apply_tie_outside_gate` demotes exactly the `assigned` rows whose read is outside, sets the worker
  gate's fields, recounts the family, and leaves other rows and families untouched. Red first.
- **G2 (determinism, real data).** The human anchor batch b4 (MCL23's batch) with the fixed binary on 1 and on 4
  region threads gives byte-identical `assignments.tsv`, and both MCL23 rows are `tied`. Before the fix these
  differed.
- **G3 (no other change, real data).** gorilla_KB3781 full assign, fixed binary against the exact-tie binary, same
  args (`--region-threads 4`): `assignments.tsv` differs only in rows that go assigned -> tied and have
  `tie_outside_catalog` = 1. Their count is reported (the exact-tie run had 0 such rows at 4 threads, so 0 changes
  are expected there). Then the fixed binary on 1 region thread (the pipeline default) must give 0 assigned rows
  with `tie_outside_catalog` = 1. The shipped run had 1.
- **G4.** `cargo test --release --no-fail-fast` passes, except the 3 known `pipeline_driver_dispatch` failures.

Outputs `/mnt/linuxdisk/tmp/o2_tie/`.

### Amendment 1 (2026-10-13, after G2's first run, before G3)

**G2's first run.** `assignments.tsv` was identical on 1 and 4 threads, but `families.tsv` and `quant.tsv` were not,
in MCL23 only. The statement "quant is unaffected" above was **wrong**, for two reasons:
- On the read-star path the abundance is `star_soft_abundance`, the normalised sum of the molecules' posteriors.
  The worker computes it after its own gate, so an in-time demotion (uniform posterior) moves it.
- `uniq` / `uniq_agree` are counted in the worker from the then-current status, and the recount does not rebuild them.

**The fix now also:**
- records each row's unique-mapper contribution, `FamilyAssignment::uniq_contrib`, and subtracts it when a row is
  demoted late;
- recomputes the read-star soft abundance for families with a demoted row. This uses the same function
  (`star_soft_abundance_of`, which `star_soft_abundance` now calls), the same row order and the same width, so it is
  bit-identical to what the worker computes when its gate fires in time. It is skipped under the VG re-align leg,
  whose own EM replaces the abundance.

**Bars, extended.**
- **G2:** every output identical on 1 and 4 threads.
- **G3 adds:** the fixed binary on the full gorilla run gives byte-identical `assignments` / `families` / `quant` on
  1 and 4 threads.

### Results (2026-10-13; fix 014126be + d5d982a3; outputs `/mnt/linuxdisk/tmp/o2_tie/`)

- **G1 PASS.** `apply_tie_outside_gate` was red first. It now demotes only the `assigned` outside rows, recounts the
  family, takes the `uniq` contribution out and recomputes the abundance bit-identically. lib 1,048 / 0.
- **G2 PASS** (human anchor batch b4, the fixed binary).
  - 1 and 4 region threads give byte-identical `assignments`, `families`, `quant`, `star_reads`, `famcn_readonly` and
    `family_join`.
  - Both MCL23 rows are `tied` in both runs. On 4 threads the post-pass demoted exactly those 2 rows; on 1 thread the
    worker's own gate caught them.
  - Before the fix, the same binary gave `assigned` on 4 threads and `tied` on 1.
- **G3 PASS** (gorilla_KB3781, full assign).
  - The fixed binary on 4 threads against the exact-tie binary on 4 threads: 60,214 / 60,214 rows identical, and
    `quant` / `families` identical. No late demotion was needed in this run.
  - The fixed binary on 1 thread (the pipeline default): the post-pass demoted **1** row. That is the same
    deterministic leak as the shipped run (624 -> 623).
  - 1 and 4 threads: `assignments` / `families` / `quant` / `famcn_readonly` / `family_join` byte-identical.
    0 assigned rows with `tie_outside_catalog` = 1 on either.
- **G4 PASS.** 1,287 passed. The only failures are the 3 known `pipeline_driver_dispatch` ones. Warnings unchanged
  (58).

**Reading.**
- The §6gz gate is now applied once the outside-tie registry is complete, so no `Assigned` call can carry
  `tie_outside_catalog` = 1.
- `--region-threads` is again byte-identical across every assign output checked.
- The remaining in-worker use of the order-dependent status is GTF routing under `--families --gtf` (not a pipeline
  combination). It is documented, not fixed.


---

## PREREG_o2_proof_2026-10-17.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_proof_2026-10-17.md`)*

## Pre-registration: does O2 work? (held-out near-identical families, 2026-10-17)

**User (2026-10-17):** "In a previous session we found that GH, USP17L, TSPY, RBMY, BPY2, CDY, OR4F and LCE1; six of these
eight have near-identical copies; can we use these to prove that O2 works?" Then, on the proposed design: "yes please".

**Why those eight cannot be the proof:** they have been development set 9 since round 8, so the integrated O2's rules were
tuned on them. Five of the eight test only abstention (copies at divergence 0 to 0.0005). They are used here as development
evidence. The claim is decided once, on the held-out families below.

### The question

The integrated O2 (the O1 -> O2 table `graph/o2.integv.m2.tsv`, pipeline v10) assigns each read of a family to one of the
family's copies or abstains. **"O2 works" means three things, each with a bar:**
- **precise:** it rarely assigns a read to the wrong copy;
- **calibrated:** it does not guess on reads whose sequence cannot tell their copy apart;
- **powerful:** it assigns the reads whose sequence does tell.

### The oracle (`bench/o3_graph/o2proof.py`, frozen with this document)

Every read of a catalogued copy (truth kind `present` or `silent`; absent and uncatalogued copies are O3's) is classed by
what its own sequence can tell.
- **Costs:**
  - The read's alignments in `E0.bam` (primary and secondary; supplementary pieces left out) are placed on the copy whose
    locus each overlaps most. Off every copy, an alignment goes to its 10-kb genome bin.
  - A locus costs its best alignment's NM plus the read bases that alignment leaves unaligned.
- **Margin** = the best family rival's cost minus the true copy's. Family rivals are the other copies of the read's family in
  the reference genome, catalogued or not.
- **Classes:**
  - `off_family`: a locus outside the family costs no more than the true copy. Which loci form the family is O1's
    question, so these reads are reported only.
  - `undecidable`: margin <= 0. A family rival explains the read at least as well.
  - `weak`: margin 1. Reported only.
  - `decidable`: margin >= 2 (O2's own margin), or no family rival aligned.
  - `true_unaligned`: the aligner never placed the read on its true copy. Reported only.
- **Outcomes:** `correct` (the read's own copy), `wrong` (another copy, another family, or a locus that is no copy),
  `abstain`.
- **Baselines** (reported, no bar):
  - `mm2`: always the primary alignment's copy.
  - `mm2_q1`: the primary at MAPQ >= 1, else abstain.
  - `current`: copy_assign's O2.

**Development evidence** (spent: the 16 v10 development runs, sets 8-10 and rounds 10-13, both profiles;
`/mnt/linuxdisk/tmp/o2proof/*.json`):
- **Wrong / assigned:** 0 in 13 runs. The others are 1/13,150 (set 8 A119b), 1/13,456 (set 8 public) and 16/12,909
  (round 12 A119b: the DUX4 / KRTAP2 coincidences no sequence rule separates).
- **Wrong on undecidable reads:** at most 11 of 14,203 (round 12 A119b). `mm2` errs on 23-89% of them.
- **Correct on decidable reads:** at least 0.991 in every run.
- **Set 9 alone (the eight families above):** 0 wrong of 5,023 / 4,611 assigned; decidable 4,099/4,137 and 3,730/3,776.
  776 / 788 reads that tie with a locus outside the family are placed on their true copy.

### Families (held out; picked before any build)

The candidates were listed in a fixed order before any build (`families.py` GROUPS `o2proof`, c647219b). No family is in any
earlier list; VCY, HSFY and PRY were scored once before, for copy_assign's O2 only, by another simulator (2026-10-12).
- **Build rule:** round 12's, `families.py build families_o2p_all.json --cap 1000 --link locus --min-cover 0.8`, every
  group (111 s).
- **The set** is the first 8 families with >= 2 copies, in list order. HBG gave no family.

| family | copies | isoforms | genes | chromosome |
|---|---|---|---|---|
| HBA | 2 | 8 | HBA1, HBA2 | chr16 |
| HSPA1 | 2 | 2 | HSPA1A, HSPA1B | chr6 |
| TUBA3 | 3 | 5 | TUBA3C, TUBA3D, TUBA3E | chr2, chr13 |
| VCY | 2 | 2 | VCY, VCY1B | chrY |
| TPSB | 2 | 7 | TPSAB1, TPSB2 | chr16 |
| C4 | 2 | 8 | C4A, C4B | chr6 |
| HSFY | 2 | 8 | HSFY1, HSFY2 | chrY |
| SERPINB34 | 2 | 6 | SERPINB3, SERPINB4 | chr18 |

17 copies. Files: `/mnt/linuxdisk/tmp/o3_graph/families_o2p.json` (this set), `families_o2p_all.json`, `families_o2p.log`.

### Worlds and runs

- **Worlds:** six world types (INTACT, CROSS, OUTSIDE, REMOVED, BRANCH, SILENT) under `profiles/A119b_testis_v2.json` and
  `profiles/public_testis_v2.json`, with a diploid donor (`--sim-het-draw`).
- **Seeds:** 401, 402, 403, fixed in advance with no seed selection. That makes 6 cells (seed x profile) and 36 worlds.
- **Pipeline:** v10, round 14's flags plus `--haplotype-check-genome`.

```
F=/mnt/linuxdisk/tmp/o3_graph/families_o2p.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/o2proof
V="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome --sim-het-draw --heavy-current --step-timeout 1800 --het 0 --jobs 2 --threads 2"
for S in 401 402 403; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F o2proof $S --local $O/local $V --sim-profile bench/o3_graph/profiles/${P}_testis_v2.json
done; done
python3 bench/o3_graph/o2proof.py score $O/<cell>.json $O/<cell>/{INTACT,CROSS,OUTSIDE,REMOVED,BRANCH,SILENT}
```

### Bars (integrated O2), per profile pooled over its 3 seeds and 6 worlds, and in each of the 6 cells

- **B1, precise:** wrong / assigned <= 0.005.
- **B2, calibrated:** wrong on `undecidable` reads <= 0.01 of the `undecidable` reads.
- **B3, powerful:** correct on `decidable` reads >= 0.95 of the `decidable` reads.

**Power floors:**
- A pooled profile needs >= 500 `undecidable` and >= 500 `decidable` reads.
- A cell needs >= 100 of each.
- A bar below its floor is UNDERPOWERED for that profile or cell, not failed.

**Validity:**
- B2 is a test only if the set is hard: `mm2` must err on >= 25% of the pooled `undecidable` reads. Otherwise B2 is
  UNDERPOWERED (an easy set).
- B1 needs >= 500 assigned reads per pooled profile.

**Verdict:**
- **PASS:** B1, B2 and B3 hold in both pooled profiles and in every powered cell.
- **FAIL:** any powered bar fails.
- **UNDERPOWERED:** otherwise.

**Reported with no bar:** copy_assign's O2 on the same metrics, `mm2` and `mm2_q1`, the `weak`, `off_family` and
`true_unaligned` classes, and each world and family.

**Rules:**
- Nothing is changed after the first run.
- A bug found in the oracle or the pipeline is recorded as an amendment, with the reason, before any rerun.
- Development sets are never rescored to choose between versions of the oracle.

### What a PASS would and would not claim

**A PASS would claim** that on simulated reads of 8 held-out near-identical families, under two read profiles and three
seeds, the integrated O2 is precise, calibrated and powerful in the senses above.

**It would not claim:**
- **Accuracy on real reads.** For the Y families there is no matched donor: the complete Y assemblies come from cell lines,
  and these genes are expressed in testis only.
- **Anything about off-family ties.** In the simulation an unlisted paralog never expresses; in real data it can.
- **Anything about the families' scope.** That is O1's question.

### Results (2026-10-17, run once as pre-registered; `/mnt/linuxdisk/tmp/o3_graph/fast/o2proof/`)

The six cells ran in 11 min in total, with no failure. The oracle was unchanged from b3b5df59.

**Verdict: PASS.** B1, B2 and B3 hold in both pooled profiles and in all six cells, every one of them powered. The
validity condition holds too: `mm2` errs on 53% of the undecidable reads in both profiles.

| cell | assigned | wrong (B1) | undecidable: wrong (B2) | decidable: correct (B3) |
|---|---|---|---|---|
| A119b 401 | 3,938 | 0 | 0 / 1,236 | 3,640 / 3,649 (0.998) |
| A119b 402 | 2,317 | 0 | 0 / 645 | 2,052 / 2,083 (0.985) |
| A119b 403 | 5,237 | 0 | 0 / 910 | 5,151 / 5,170 (0.996) |
| public 401 | 3,851 | 0 | 0 / 1,432 | 3,531 / 3,540 (0.997) |
| public 402 | 2,027 | 0 | 0 / 849 | 1,724 / 1,762 (0.978) |
| public 403 | 5,261 | 0 | 0 / 1,018 | 5,198 / 5,233 (0.993) |

**Pooled per profile, every assigner (catalogued reads of the 6 worlds x 3 seeds):**

| assigner | profile | assigned | wrong | undecidable: wrong / correct | decidable correct | weak: correct / wrong |
|---|---|---|---|---|---|---|
| integrated O2 | A119b | 11,492 | 0 | 0 / 381 of 2,791 | 0.9946 of 10,902 | 268 / 0 of 2,032 |
| integrated O2 | public | 11,139 | 0 | 0 / 375 of 3,299 | 0.9922 of 10,535 | 311 / 0 of 2,013 |
| copy_assign O2 | A119b | 8,428 | 0 | 0 / 0 | 0.7729 | 2 / 0 |
| copy_assign O2 | public | 8,054 | 0 | 0 / 0 | 0.7642 | 3 / 0 |
| mm2 at MAPQ >= 1 | A119b | 12,955 | 21 | 21 / 0 | 1.0000 | 2,032 / 0 |
| mm2 at MAPQ >= 1 | public | 12,582 | 34 | 34 / 0 | 1.0000 | 2,013 / 0 |
| mm2, always | A119b | 15,729 | 1,483 | 1,479 / 1,312 | 1.0000 | 2,032 / 0 |
| mm2, always | public | 15,847 | 1,761 | 1,761 / 1,538 | 1.0000 | 2,013 / 0 |

**Reading:**
- **The integrated O2 made no wrong assignment in 22,631.** On the reads whose sequence cannot tell their copy apart,
  always-assigning errs on 53%.
- **It also placed 756 reads the edit-distance oracle calls undecidable, all correctly.** Its gap-compressed distance
  counts an indel run once and aligns clipped ends. That separates copies where NM plus clipping ties: homopolymer errors
  inflate NM equally on both copies.
- **copy_assign is as precise but much less powerful,** with 76-77% of decidable reads against 99.2-99.5%.
- **A plain MAPQ >= 1 filter is the strong baseline here.** Every decidable read has MAPQ >= 1, so it recovers all of
  them, and all weak ones. It pays for that with 55 wrong assignments, all on undecidable reads.
  - The integrated O2 trades recall on weak reads (13-15% against 100%) and 0.5-0.8% of decidable reads for zero errors.
  - The claim is precision and calibration at near-MAPQ power, not more power than MAPQ.
- **The `off_family` class is nearly empty here** (4 reads): these families' paralogs are all in the family.

**Not part of this test, but seen in the same runs (the round-8 end-to-end bars, O1 -> O2 -> O3):**
- **OUTSIDE fails B2 in all six cells:** integrated 138-171 errors, against 4-77 for the current method.
- **The cause is O1 missing the uncatalogued twin.**
  - O1 seeds the uncatalogued copy as a new locus in 6 of the 8 families, but not in C4 or HSFY, whose twins are
    near-identical.
  - Those copies' reads then have a single candidate in the catalog, and O2 places them on the catalogued twin. The
    current method abstains on them.
- **A119b REMOVED also fails B2 in two cells** (3-4 errors against 0-1).
- **These are O1's family-scope misses, which this test excludes by construction.** They are the open item: a read whose
  genome alignment ties with a locus O1 did not seed should be abstained on.


---

## PREREG_o2_hg002_tied_2026-10-18.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_hg002_tied_2026-10-18.md`)*

## Pre-registration: does O2 put AS-tied reads on the right copy? (HG002, the donor's own genome as truth, 2026-10-18)

**User (2026-10-18):**
- "can we prove that we O2 and O3 work? ... can we prove we can effectively assign reads to copies and point out possible
  cases in which there is a missing copy that would explain the PSVs?"
- "minimap2 already assigns them if they have mapq>0 the problem is the multimappers that have equal AS so we don't know to
  which copy do they truly belong and the distinction between primary and secondary is meaningless (coin toss)"
- "Ok we need to prove it in that way, any dataset it just needs to assign reads to copies based on PSVs for tied AS reads"
- "lets ensure that O2 can only work for tied AS reads, everything else is already resolved by minimap, then O3 is just an
  extension in which if we notice the PSVs nothing in the genome explains them"

**The split this tests:**
- **Resolved:** a read whose best-AS alignments all lie on one copy belongs to minimap2.
- **O2:** an AS-tied read. O2 may choose only among the copies it ties on, or abstain.
- **O3:** reads whose PSVs no genome copy explains. An O3 call overrides any placement.

### Why HG002

- **Real reads with their donor's genome.** PacBio Kinnex FLNC of HG002 (37,200,651 reads; Revio, 2023; the BAM header
  names no tissue) and HG002's diploid T2T assembly (Q100 v1.1, maternal and paternal), with the Q100 chain to CHM13 v2.0.
- **The 10-17 proof was simulation only** and said so (`docs/PREREG_o2_proof_2026-10-17.md`: "never claim O2 on real
  reads from this").
- **Held out.** No O2 rule was developed or tuned on HG002 reads. They were used once before, on 2026-09-29, to test a DNA
  presence gate. That gate failed and is not in the pipeline (register 1180).

### The population: reads minimap2 can only toss a coin for

1. **Candidates.** Every FLNC read with a minimap2 `map-hifi` hit to a transcript of the adopted multigene pool, at >= 0.90
   identity over >= 200 bp: 633,939 reads.
   - The pool is `seqfamilies.py build --rule multigene` (register 1328): 439 families, 1,685 copies, 3,583 transcripts
     (`r19seq/families19cal_full2_all.json`, set `seq19`).
2. **CHM13.** The candidates mapped to the whole of CHM13 v2.0 with
   `minimap2 -a -x splice:hq -uf --secondary=yes -N 100 -p 0.5`.
3. **Tied.** A read is tied when its best AS is shared by alignments on two or more copies of one family. A copy is matched
   by exon overlap.
   - Reads tied across two families are counted, not scored.
   - These reads' primary alignment is minimap2's coin toss among the tied copies. (Final review: 1,034 of the 23,424
     primaries are not on a tied copy; 947 lie off every copy and 87 on another copy.)
   - **Result:** 23,424 reads tied within one family, plus 543 across two.
     - **Size of the tie:** 23,266 tie on 2 copies, 141 on 3, 15 on 4, 1 on 5 and 1 on 7.
     - **Primary MAPQ:** 0 for 23,131; 1-29 for 82; 30 or more for 211.
     - **Families:** 72 have tied reads.
4. **Pseudoautosomal families are a separate stratum.** A family is pseudoautosomal when every copy lies in CHM13
   v2.0's PAR1 or PAR2 on chrX or chrY (`hg002truth.pseudoautosomal`): 16 of the 439 families, 12 of them with tied reads.
   - **Why apart:** in this male donor, the chrX and chrY copies are the maternal and paternal alleles of one recombining
     locus. Telling them apart is phasing, not copy assignment. Also, CHM13 v2.0's chrY is HG002's own chrY.
   - **Size:** they hold 9,410 of the 23,424 tied reads. SLC25A6 alone has 7,338.
   - **The verdict is decided on the other 60 families with tied reads (14,014 tied reads).** The pseudoautosomal
     stratum runs separately, scored the same way, and is reported with no bar.

### The truth: the donor's own haplotypes (`bench/o3_graph/hg002truth.py`, frozen with this document)

- **Mapping.** Each tied read is mapped to each HG002 haplotype with the same settings.
- **Lift.** Its best-scoring placements, over both haplotypes, are lifted to CHM13 through the chain. Each placement is put
  on the copy its lifted blocks overlap most.
- **Classes:**
  - `unique`: every best placement lands on one copy, and CHM13 tied the read on it. **Truth = that copy.** Two
    haplotypes landing on one copy are the copy's two alleles.
  - `unique_elsewhere`: as `unique`, but on a copy of the family that CHM13 did not tie the read on. **Truth = that copy.**
  - `own_tie`: the best placements land on two or more copies. The donor's copies tie too, so the read is undecidable.
  - `elsewhere`: the best placements lift onto CHM13, but onto no copy of the family. The read comes from another CHM13
    locus.
  - `no_counterpart`: no best placement lifts onto CHM13 at all. This is a donor copy that CHM13 lacks: O3's case.
  - `poor_match`: the donor's best placement has NM > max(2, 0.3% of its aligned bases), arm 4's floor
    (`bench/o2_truth/haplotype_truth.py`). Not scored.
  - `unmapped`: no alignment to either haplotype. Not scored.
- **The donor margin** of a `unique` read is reported but not used: the best AS minus the best AS of a placement off the
  truth copy.

**Counts (the truth step only; no assigner's output had been seen):**

| stratum | tied | `unique` | `unique_elsewhere` | `own_tie` | `elsewhere` | `no_counterpart` | `poor_match` |
|---|---|---|---|---|---|---|---|
| non-pseudoautosomal | 14,014 | 466 | 118 | 9,948 | 948 | 31 | 2,503 |
| pseudoautosomal | 9,410 | 624 | 0 | 7,151 | 0 | 0 | 1,635 |

- **No read was unmapped.**
- **`poor_match`:** NM quartiles 6, 9 and 13 (the floor is kept as arm 4 set it).
- **Donor margin of the `unique` reads:** quartiles 10, 15 and 20 AS points over both strata. The non-pseudoautosomal
  ones alone are 5, 10 and 25 (final review); the 624 pseudoautosomal ones are all 15-20. 44 reads have no
  placement off their copy at all.
- **So most AS ties are ties in the donor too.** 71% of the non-pseudoautosomal tied reads are `own_tie`, as register 1180
  found for PAR and EIF3C reads.

### The world and the run

- **Families:** every non-pseudoautosomal `seq19` family with at least one scored tied read, set `hg002`
  (`hg002truth.py world truth.tsv ... seq19 1 world --par exclude`): 54 families, 280 copies.
  > AC0064533 AC1275021 AGAP4 ANKRD36C BMS1P1 BOLA2B CBWD2 CCT6P3 CFC1B CSAG2 CTAGE15 DMRTC1 DUSP22 EIF3C F8A2 FAM156B
  > FCGR2A GOLGA2P10 GOLGA6L10 GOLGA8A GTF2H2C GTF2IRD2 H2AB2 H2AC18 H4C14 HSFX2 HSFY1 HSPA6 IKBKG LRRC37A2 MAGED4B
  > MSTRG3150 NBPF1 NCF1 NOMO1 NPIPB2 NUTM2B OR4G2P PDPR PMS2P11 RASA4B RBMY1A1 RBMY2BP RPL6 SERF1B SLC35E2B SMG1 SMN2
  > SPDYE1 SULT1A1 TBC1D3B ULK4P2 USP17L11 ZC3H11A
- **Reads:** every candidate whose CHM13 primary lies within 1 kb of a copy of those families, as sequenced
  (`realreads.py pool`): 83,146 reads.
  - **Truth reads in the world:** 582 of the 584.
  - **Scored tied reads in the world:** 10,839 of 11,511. A read whose coin-toss primary lies off every copy is not
    pooled. The scorer reports these reads as `not_in_world`.
- **The world's truth** is the tied reads with a donor copy (`unique` and `unique_elsewhere`). Every other read is in the
  world but not scored.
- **One INTACT world of real reads**, with no simulation. The seed is 1 and unused.
- **Pipeline:** v14 as in round 19 (register 1321, not changed since), with `--step-timeout 21600`. The timeout does not
  change outputs.
  ```
  P=/mnt/linuxdisk/tmp/o3_graph/hg002proof; O=/mnt/linuxdisk/tmp/o3_graph/fast/hg002
  V14="--o2-distance gapcomp --member-genome-check --member-attribution --union-exact --o2-genome veto --per-read-detectability --twins --error-floor --o1-hosts --branch-genome-novel --mosaic-shared --haplotype-check --haplotype-check-genome --branch-linear --aln-cache --veto-abstain --locus-pairwise --o2-significance --structure-check --struct-rescue --short-path-align --ref-copies --sim-het-draw --heavy-current --step-timeout 21600 --het 0 --jobs 2 --threads 2"
  python3 bench/o3_graph/fastset.py run $O/nonpar $P/world/families.json hg002 1 --worlds INTACT --real-reads $P/world/real --local $O/local_nonpar $V14
  python3 bench/o3_graph/o2proof.py score $O/nonpar_o2proof.json $O/nonpar/INTACT --tied $P/tied.tsv
  python3 bench/o3_graph/hg002truth.py score $P/truth.tsv $P/tied.tsv $O/nonpar/INTACT $O/nonpar_score.json
  python3 bench/o3_graph/hg002truth.py bars $O/nonpar_o2proof.json $O/nonpar_score.json $O/nonpar_bars.json
  ```
  - **The pseudoautosomal stratum** runs the same way after the premise test, with `--par only`, into `$O/par`.
  - **Known risk.** v14's genome veto trusts a primary at MAPQ >= 30 even when a secondary ties its AS (round 20,
    register 1327). 211 of the 23,424 tied reads have such a primary. The scoped systems keep such a call when it is on a tied copy.
- **Systems under test, each judged on its own:**
  - **Integrated O2:** the O1 -> O2 table `graph/o2.integv.m2.tsv`.
  - **Current O2:** `copy_assign` (`/mnt/linuxdisk/tmp/o3_fb_build/copy_assign.fb`, built 2026-10-15 from a tree with the
    exact-tie fix 7ac557e7 and the gate fix d5d982a3).
  - **The user's scoped O2,** `integrated_tied` and `current_tied`. These are the same calls, but a tied read's call stands
    only on a copy it ties on, and is otherwise an abstention; an untied read takes minimap2's copy. On this population
    every read is tied, so only the first rule applies. It is computed from the same run, from `tied.tsv`.
- **Baselines, with no bar:**
  - `mm2`: always the primary's copy (the coin toss).
  - `mm2_q1`: the primary's copy when MAPQ >= 1, else abstain.

### Scoring

- **O2 level** (`o2proof.py score ... --tied tied.tsv`; the oracle is unchanged since b3b5df59, and `--tied` only adds
  the scoped systems): the world's truth reads.
  - Each read is classed by what its own sequence can tell on the world's `E0.bam`: `decidable`, `weak`, `undecidable`,
    `off_family`, `true_unaligned`.
  - Each assigner's outcome is `correct` (the donor's copy), `wrong` (another copy, another family, or a locus that is no
    copy) or `abstain`. O3 is ignored at this level.
- **By donor class** (`hg002truth.py score`): every scored tied read in the world, for each assigner, with O3.
  - An absent-copy call takes precedence over a placement.
  - `unique` / `unique_elsewhere`: right only on the donor's copy. An O3 call on it is an error, since its copy is in
    CHM13.
  - `own_tie`: a placement on a tied copy is a guess; anywhere else, or an O3 call, it is an error.
  - `no_counterpart`: any placement is an error; an O3 call is right.
  - `elsewhere`: a placement is an error; an O3 call is reported apart.

### The premise: minimap2 resolves the untied reads

This test is run after the O2 run, with the same truth code and settings.
- **Reads:** the candidates whose best-AS CHM13 alignments all lie on one copy (`hg002truth.py untied`).
- **Truth:** mapped to both haplotypes (`truth_rest.sh`).
- **P1:** minimap2's copy is the donor's (`unique`) on >= 0.99 of the untied reads with a donor copy (`unique` +
  `unique_elsewhere`), with a floor of 100 reads.
- **Reported:**
  - `own_tie`: the donor's copies tie where CHM13's do not, so minimap2's choice cannot be checked.
  - `no_counterpart` and `elsewhere`: O3's domain.

### Bars

**O2 (each system), on the world's truth reads (`unique` + `unique_elsewhere`), O2 placements:**
- **T1, precise:** wrong / placed <= 0.02.
- **T2, better than the coin toss:** the upper Wilson 95% bound of the system's wrong / placed is below the lower bound of
  `mm2`'s wrong / all.
- **T3, calibrated:** placements on `own_tie` reads <= 0.05 of them, from the by-class score.
- **T4, uses the evidence:** correct on `decidable` reads >= 0.90 of them.

**Power floors:**
- T1 and T2 need >= 100 placed reads.
- T3 needs >= 100 `own_tie` reads.
- T4 needs >= 100 `decidable` reads.
- A bar below its floor is UNDERPOWERED, not failed.

**Validity:** `mm2` must be wrong on >= 25% of the truth reads. Otherwise the population is not a coin toss, and T2 is
void.

**Systems judged:** `integrated`, `current`, `integrated_tied`, `current_tied` (`hg002truth.py bars`, frozen with this
document: `FLOOR, T1_MAX, T3_MAX, T4_MIN, T5_MAX, VALID_MIN = 100, 0.02, 0.05, 0.90, 0.05, 0.25`).

**O2 verdict (per system):**
- **PASS:** T1 and T2 are powered and hold, and T3 and T4 hold or are underpowered.
- **FAIL:** any powered bar fails.
- **UNDERPOWERED:** otherwise.

**O3 (each system), from the by-class score:**
- **T5, safe:** O3 calls on `unique`, `unique_elsewhere` and `own_tie` reads <= 0.05 of them.
- **Reported:** O3 calls on `no_counterpart` reads, the "missing copy that would explain the PSVs" (recall, no bar), and
  on `elsewhere` reads.

**Reported with no bar:**
- Each assigner's correct / unique-truth reads (recall): the share of the coin-toss reads a system rescues.
- The `weak`, `undecidable`, `off_family` and `true_unaligned` classes.
- The donor margin, MAPQ and family breakdowns.

**Rules:**
- Nothing is changed after the run.
- A bug found in the truth, the scorers or the pipeline is recorded as an amendment, with the reason, before any rerun.

### What a PASS would and would not claim

**A PASS would claim** that on real HG002 reads whose CHM13 alignment score ties between copies, the system's placements
agree with the donor's own genome at >= 98%, better than minimap2's coin toss. It abstains on reads the donor's genome
cannot resolve, and places the reads its sequence decides.

**It would not claim:**
- **Families HG002's library barely expresses.** The testis-only Y families give 5 tied reads (HSFY1 3, RBMY1A1 1,
  RBMY2BP 1).
- **Power where CHM13 has no evidence.** Register 1182: a tied read often carries the donor's own alleles, which CHM13's
  copies lack equally. A CHM13-based assigner can only abstain on it.
- **Gorilla, or other donors.**

### Results: the non-pseudoautosomal world (2026-10-19, run once as pre-registered; `/mnt/linuxdisk/tmp/o3_graph/fast/hg002/`)

The run took 48 min (map_union 747 s, current 726 s, O2 650 s, O3 543 s), peak RSS 14.4 GB, with no failure. The code was
unchanged from 7e865cba.

**Verdict: UNDERPOWERED for all four systems.**
- **No system placed 100 of the 584 truth reads,** so T1 and T2 are below their floor. Their point estimates are in the
  table below.
- **Validity holds:** `mm2` is wrong on 341 of 582 placed (58.6%), or 341 of all 584 (58.4%, as the code computes it).
- **T3 holds for every system** (powered).
- **T4 is underpowered:** only 8 of the truth reads are `decidable`.
- **T5 holds for the integrated system and FAILS for the current one.**

| system | placed | wrong (T1) | T2: upper / coin toss's lower | own_tie placed (T3) | O3 calls on present-copy reads (T5) |
|---|---|---|---|---|---|
| integrated | 57 | 7 (0.123) | 0.232 / 0.544 | 96 / 9,607 (0.010) | 40 / 10,189 (0.004) |
| integrated_tied | 57 | 7 (0.123) | 0.232 / 0.544 | 96 / 9,607 (0.010) | 40 / 10,189 (0.004) |
| current | 41 | 5 (0.122) | 0.255 / 0.544 | 40 / 9,607 (0.004) | **6,999 / 10,189 (0.687)** |
| current_tied | 39 | 3 (0.077) | 0.203 / 0.544 | 37 / 9,607 (0.004) | **6,999 / 10,189 (0.687)** |
| `mm2`, always | 582 | 341 (0.586) | | | |
| `mm2` at MAPQ >= 1 | 22 | 12 (0.545) | | | |

**Reported breakdowns (no bar):**
- **Most AS ties are ties in the donor too, and both O2s abstain on them.** 9,607 of the world's scored tied reads are
  `own_tie`; the systems place 0.4-1.0% of them.
- **CHM13 holds almost no evidence for the rest.** The oracle finds 560 of the 584 truth reads `undecidable`, 14
  `off_family`, 8 `decidable` and 2 `true_unaligned`. These reads carry the donor's own alleles, which both CHM13 copies
  lack (register 1182).
- **MAPQ means nothing on a tie.** A tied read's primary at MAPQ >= 1 is wrong 12 times in 22, as the user said. That is
  the local remap's primary (`E0.bam`, the o2proof baseline). On `chm13.bam`, the alignment that defines the tie, it is 9
  in 17 (final review).
- **Placements by donor margin** (AS points between the truth copy and the best placement off it):

  | donor margin | reads | integrated: right / wrong | current: right / wrong | `mm2`: right / wrong |
  |---|---|---|---|---|
  | no placement off the truth copy | 44 | 42 / 0 | 0 / 0 | 22 / 22 |
  | <= 5 (about one site) | 169 | 1 / 5 | 3 / 0 | 83 / 86 |
  | 6-10 | 192 | 2 / 0 | 4 / 1 | 66 / 126 |
  | > 10 | 179 | 5 / 2 | 29 / 4 | 70 / 107 |

- **Each system's right placements come mostly from one family.**
  - The integrated system's 50 are 42 IKBKG (the donor has no second placement for them), 3 NBPF1, 3 GTF2H2C, 1 RPL6
    and 1 ZC3H11A.
  - Where the donor has a rival copy, the integrated system is 8 right and 7 wrong, no better than the coin toss.
  - The current system's 36 are 26 NBPF1, then 2 each of SPDYE1, GTF2H2C and GOLGA8A, and 1 each of EIF3C, RPL6,
    ZC3H11A and AC1275021.
  - Its PSV tie-break is right 36 times in 41, and 36 in 39 when scoped.
- **The integrated system's 7 wrong placements:**
  - Five are reads the donor separates by 4-5 AS points, about one site: EIF3C four times, SULT1A1 once.
  - Two are AC1275021 reads whose CHM13 truth copy aligns 228-492 bases worse than its rival. The donor's copy has a
    structure CHM13's lacks.
- **O3 on donor copies CHM13 lacks.** The world holds 30 `no_counterpart` reads: NPIPB2 12, GTF2H2C 9, SPDYE1 5, USP17L11
  3 and BMS1P1 1.
  - The integrated O3 flags 10 of them, all NPIPB2: a missing NPIPB copy. It places 7 and abstains on 13.
  - The current O3 flags none of them and places 16. Instead it flags 6,907 EIF3C reads, all but 22 of which tie in the donor too:
    allelic variation called a missing copy (register 1181).
- **The round-8 end-to-end bars on the truth reads** (fastset's own): integrated 50 correct and 12 errors; current 33 and
  37.

**What it means:**
- **On real AS-tied reads, CHM13's PSVs rarely carry the answer, and a CHM13-based O2 can mostly only abstain.** Both do:
  the calibration bar T3 holds on real data.
- **Where they place, they beat the coin toss (88% against 41%) but are far from the 98% bar,** on too few reads to
  decide. The current system's PSV tie-break is the better of the two on reads with a real rival (36/41).
- **The answer lies in the donor's own copies.** The reads that decide them are the tied reads' own alleles, so the lever
  is to learn the sample's copies from its uniquely placed reads and assign against those. That is a design question; it
  was not tested here.

### Results: the premise (2026-10-19, `hg002proof/truth_rest.sh`, as registered)

- **Reads:** 599,678 of the 633,939 candidates are resolved by minimap2, meaning their best-AS CHM13 alignments all lie on
  one copy.
- **Truth classes:** `unique` 215,514, `unique_elsewhere` 6,342, `own_tie` 3,112, `elsewhere` 1,353, `no_counterpart`
  259,889, `poor_match` 113,468.

**P1: FAILS.** minimap2's copy is the donor's on 215,514 of 221,856 untied reads with a donor copy (0.971, against the
0.99 bar; powered).

**Reported breakdowns (no bar):**

| families | `unique` | `unique_elsewhere` | P1 |
|---|---|---|---|
| not pseudoautosomal, no chrM copy | 207,998 | 4,481 | 0.979 |
| pseudoautosomal | 7,515 | 1,861 | 0.802 |
| with a chrM copy | 1 | 0 | |

- **`no_counterpart` here is a truth artifact.** The Q100 chain has no chrM chain. So the donor's chrM placements of
  258,334 mitochondrial-transcript reads (families MTND6, MTRNR2L8 and MTCO2, whose other copies are nuclear) do not lift.
  Outside them, 1,555 reads are `no_counterpart`. The tied reads were barely affected (31 `no_counterpart`).
- **Where minimap2 misses.** Its misses concentrate in copy-number-variable segmental duplications.
  - The largest are P2RY8 1,149 (pseudoautosomal), GPAT2 1,045, SMN2 606, AKAP17A 473 (pseudoautosomal), THOC3 416,
    FCGR2A 400, WASHC1 344, CCDC144A 210, TBC1D3B 199 and PLCXD1 178 (pseudoautosomal).
  - Donor margins of the misses: quartiles 5, 10 and 191 pooled, driven by the 1,861 pseudoautosomal misses, all at
    margin 5. The non-pseudoautosomal misses alone are **10, 20 and 466** (final review), so most of the CNV
    segmental-duplication misses are not one-site differences.
- **What it means.** "Everything else is already resolved by minimap2" holds for about 98% of the untied reads in
  ordinary families, not 99%.
  - The rest are reads whose own copy carries donor alleles that make another CHM13 copy fit better. The user's split
    gives them to minimap2 and gets them wrong.
  - Only a check of the sample's own variation, the O3 extension "PSVs nothing in the genome explains", can notice them.

### Results: the pseudoautosomal stratum (2026-10-19, `fast/hg002/run_par.sh`; reported, no bar)

- **The world:** 12 families, 24 copies, 22,302 reads and 624 truth reads (all `unique`). The run took 6 min.
- **Neither system attempts phasing.** The integrated system places 1 truth read (right) and the current one none. Both
  place 0 of the 7,151 `own_tie` reads.
- **O3** calls 3 (integrated) and 7 (current) of the 7,775 present-copy reads.
- **The coin toss is wrong on 287 of 624.**
- **CHM13 cannot separate them.** The oracle finds 620 of the 624 truth reads `undecidable` on CHM13: their X and Y alleles
  differ by donor variants that neither CHM13 copy carries.

### Corrections after the final review (2026-10-19)

A fresh reviewer re-derived all 23,424 truth classes from the BAMs and the chain, with 0 mismatches. It also rebuilt the
o2proof card and recomputed the score, bars and premise. No verdict changes. The text above is corrected where it pooled
strata or named the wrong baseline. The other findings are below; the code is unchanged, as the rules require.
- **P1's "minimap2's copy" is the copy holding all its best-AS alignments** (`hg002truth.resolved_copy`), as run.
  - On 10,137 untied reads (1.7%) the primary lies elsewhere: 6,731 off every copy and 3,406 on another copy.
  - Scored on the primary, P1 is 0.967 overall and 0.974 outside the pseudoautosomal families. It FAILS either way.
- **`own_tie` also holds 48 reads (39 RASA4B)** whose best placements land on a single copy plus a placement that lifts
  nowhere or off the copies (`len(found) > 1`). The prereg defined `own_tie` as two or more copies.
  - The integrated system abstains on all 48. The current system's O3 calls 40 of them.
  - Without them, its T5 is 6,959 / 10,141 = 0.686. It still FAILS.
- **A truth copy can belong to another family**, since the exon index spans all 439 families.
  - Tied reads: 0 cases.
  - Premise: 294 of the 6,342 misses, for example WASHC1 -> DDX11L9 (132 reads). P1 counting only within-family misses is
    0.973.
- **The two scorers use different coin tosses.**
  - o2proof's `mm2` is the primary of the pipeline's local remap. `hg002truth score`'s is the primary in `chm13.bam`.
  - On the 582 truth reads they are wrong 341 and 350 times.
- **The chain lift is one-to-one here.** The Q100 chain has no overlapping target blocks, so `chainlift`'s search start
  cannot miss a lift. The strand math held on 365 random lifts on both strands.
- **Register 1329's "73% own_tie" pools both strata.** The non-pseudoautosomal stratum is 71%.
- **Not judged:** whether the chain pairs true orthologs inside inverted or CNV segmental duplications. One maternal EIF3C
  placement lifts onto both CHM13 copies, and TBC1D3B.8 holds 78 of the 118 `unique_elsewhere` truths.


---

## PREREG_o2_fair_heldout2_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_fair_heldout2_2026-10-19.md`)*

## Pre-registration: O2 against minimap2's trusted calls on a deep new donor (2026-10-19)

**Why:** register 1341 held O2 out on the public testis donor and was UNDERPOWERED: 30 source reads placed, 0 wrong,
against minimap2's MAPQ >= 1 calls at 0.894. That library is shallow, so a PSV rarely reaches 3 validating reads.

This repeats the same test on a deep donor never used in O2's work.

### The donor

PacBio's public Kinnex full-length RNA release, `DATA-Revio-SCRI-Sample5-Cerebellum`
(https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/; released with the collaborator's permission).
- **Sample:** one prenatal cerebellum from Seattle Children's Research Institute (one donor). Revio, Kinnex full-length
  RNA kit.
- **Reads:** 44,962,471 FLNC reads (`2-FLNC/flnc.bam`).

### Design

- **Steps 1-2, as HG002's proof (`hg002proof/screen.sh`):**
  - the FLNC reads are screened against the pool's transcripts (`family_tx.fa`, `minimap2 -x map-hifi`); a hit at
    >= 0.90 identity over >= 200 bp makes a candidate;
  - the candidates are aligned to whole CHM13 v2.0 with A119b's parameters:
    `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`.

  That BAM is this donor's library. Its reads are the sources and the validators.
- **Step 3, exactly as register 1341** (`docs/PREREG_o2_fair_heldout_2026-10-19.md`):
  - **sources:** up to 500 per family, seed 1, MAPQ-60 primaries of >= 1 kb;
  - **fragments:** 3 per source, 300-1,200 bp, realigned with the same minimap2 parameters;
  - **truth:** the source's copy when its full read is clear by 30 AS, the fragment scored only when that copy is one of
    its tied copies;
  - **O2:** `o2psv.py svp --k 3 --validated-only`, a fragment's own source left out;
  - **scorer:** `o2perfamily.py fair` on the 281 qualifying non-pseudoautosomal families, counted by source read,
    against minimap2's primary at MAPQ >= 1.

### Bars (unchanged)

- **F1:** at least 100 source reads placed.
- **F2:** at most 2% of them wrong.
- **F3:** O2's precision at least 5 points above the baseline's.
- **F4:** more than 10 families with at least 5 placed and none wrong.

**Verdict:**
- **PASS:** F1 to F4 hold.
- **FAIL:** F1 holds and another bar fails.
- **UNDERPOWERED:** F1 fails.

**Rules:** code frozen at f2d19688 (regimen 7bf68347, scorer fad7e43b). Run script: `scri5/run.sh`. A bug is recorded as an
amendment before any rerun.


---

## PREREG_o2_fair_heldout_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_fair_heldout_2026-10-19.md`)*

## Pre-registration: O2 against minimap2's own trusted calls, held out on a new donor (2026-10-19)

**Why:** register 1340's audit voided the per-family claim on two counts.
- The data was not held out.
- O2's precision was set against minimap2's on all reads, and minimap2's primary at MAPQ >= 1 passed that rule too.

O2's measured value is precision on almost-tied reads. On HG002 and A119b (development), O2's placements were 1.1% and 0.6%
wrong, while minimap2's MAPQ >= 1 calls on the same reads were 27% and 14% wrong. This test states that claim against the
equally selective baseline, on a donor never used in O2's development.

### Design

- **The regimen:** `o2psv.py svp --k 3 --validated-only` (7bf68347), unchanged.
- **The donor:** the public testis library (ENA ERR13885926, GENCODE CapTrap; `_from_wsl/human_val/human_testis.t2t.bam`).
- **Sources:** per pool family (`r19seq`), up to 500 MAPQ-60 primaries of at least 1 kb on a copy (`o2frag.py sources
  --per-family 500 --seed 1`).
- **Fragments:** 3 per source, 300-1,200 bp. Realigned with the library's own parameters:
  `minimap2 -ax splice:hq --eqx -Y -N 50 -p 0.1 --secondary=yes`.
- **Truth:** the source's copy when the realigned full read is at least 30 AS clear of any other locus. A fragment is
  scored only when that copy is one of its tied copies.
- **Validation:** the library's own reads. A fragment's own source never validates it.
- **Families:** the 281 qualifying families (every member pair at identity [0.90, 1.00)), with the 17 pseudoautosomal
  families left out (`famident/qualify_nonpar.json`).
- **Counting:** by source read. A source counts as right for a method when every fragment the method places is on its
  copy (`o2perfamily.py fair`, e1acb42d).
- **The baseline:** minimap2's primary for the same fragments, taken only at MAPQ >= 1.

### Bars

- **F1, powered:** O2 places at least 100 source reads.
- **F2, precise:** at most 2% of them are wrong.
- **F3, adds precision:** O2's precision is at least 5 points above the MAPQ >= 1 baseline's on the same fragments.
- **F4, the user's handful:** more than 10 families where O2 places at least 5 source reads and none wrong.

**Verdict:**
- **PASS:** F1 to F4 hold.
- **FAIL:** F1 holds and another bar fails.
- **UNDERPOWERED:** F1 fails. The library is shallow (28,662 MAPQ-60 reads on the pool's copies), so this is possible.

**Development, for reference only (A119b fragments, already seen):** O2 placed 232 source reads, 1 wrong (99.6%); the
MAPQ >= 1 baseline placed 773, 103 wrong (86.7%); 18 families pass F4's rule.

**Simulation (one genome):** register 1339's test A stands. On 18 fresh worlds, version 1 placed 349 reads, 0 wrong.
Version 2 placed 1,104, 0 wrong, over rounds 23-25. In one genome minimap2's MAPQ >= 1 calls are already 99.86% right,
so the simulation claim is no harm, not a precision gain.

**Rules:** code frozen at e1acb42d. A bug is recorded as an amendment before any rerun.

### Results (`o2frag_public/run_real.log`, `fair.json`)

**UNDERPOWERED: F1 fails, with 30 source reads placed against 100 needed.**
- F2 holds: 0 wrong.
- F3 holds: precision 1.000, against the MAPQ >= 1 baseline's 0.894 (160 placed, 17 wrong) on the same fragments.
- F4 fails: 1 family (ROPN1, 6 of 6).

| public testis library | value |
|---|---|
| sources drawn / kept (clear by 30 AS) | 6,674 from 215 families / 6,070 |
| fragments with a truth | 18,210 |
| O2's decisions on almost-tied fragments | validated 49, sets 153, unexplained 178 |

The library is about 45 times shallower than A119b. Few PSV positions reach 3 validating reads, so O2 places little.

**Where minimap2's trusted calls fail, O2 abstains or is right:**

| family | minimap2 at MAPQ >= 1: right / wrong | O2: right / wrong |
|---|---|---|
| TPSB2 | 13 / 7 | abstained |
| DHRS4 | 5 / 4 | 2 / 0 |
| HSPA1A | 19 / 2 | 2 / 0 |
| THOC3 | 9 / 1 | 2 / 0 |

The direction matches development (A119b: 99.6% against 86.7%), but this donor cannot carry the claim.


---

## PREREG_o2_per_family_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_per_family_2026-10-19.md`)*

## Pre-registration: O2 works family by family in the near-identical regime (2026-10-19)

**User (2026-10-19):** "for now I only need to prove it works in a handful of families, as long as it works on more than 10
families of the regimen >90% similarity of all members but less than 100% (since these would be impossible to place)".

### The families

A pool family (`r19seq`, 439 families) qualifies when **every** pair of its members' main transcripts aligns
(`minimap2 -X -c -x asm20`, best alignment covering >= 50% of the shorter) at identity >= 0.90 and < 1.00
(`famident/family_identity.json`).
- **294 families qualify.**
- 69 have an identical pair and are left out: no read can be placed between identical copies.
- 47 have a pair that does not align.

### The regimen

Version 2: `o2psv.py svp --k 3 --validated-only`, code 7bf68347.

A read almost tied (within 10 AS) on two or more copies is placed only by a PSV the sample's own reads validate, and only
when its own bases contradict the call nowhere. Otherwise it keeps the set of copies the validated PSVs leave.

### The data

**Simulation (one genome):** every CHM13-only world of rounds 23-25.
- Profiles `*_v4_ref`; worlds INTACT, CROSS and SILENT.
- Round 23: 37 families. Round 24: 21. Round 25: 14, fresh.
- Truth: each read's simulated copy.

**Real reads:**
- **HG002**, natively almost-tied reads: truth is the donor's own assembly (`hg002proof`). poor_match is left unscored.
  HG002 was the regimen's development substrate.
- **A119b testis fragments:** truth is the full read's copy, scored among the fragment's tied copies (`o2frag_a119b`).
  These held out version 1 and served as version 2's development.

### A family works in a data set when

- O2 places at least 5 of the family's almost-tied reads;
- at most 5% of those placements are wrong;
- the family's sets hold the true copy for at least 95% of its almost-tied reads;
- O2's precision on those reads is at least minimap2's primary's.

### The claim

**More than 10 qualifying families work in simulation, and more than 10 work on real reads** (in either real data set).

**Reported:** every qualifying family with its result in each data set, including the ones that fail.

### Results (`famident/perfamily.tsv`, scorer bbf73d5e)

**By the rule: 20 qualifying families work in simulation, and 26 on real reads** (HG002 7, A119b fragments 21; CCT6P3 and
NOMO1 work in both). Only three qualifying families have any wrong real placement:
- SERF1B, 7 on HG002;
- AC0064533, 4 on HG002;
- TCAF2, 2 on A119b.

All three fail the rule.

### Audit (metric-auditor, same day): the verdict does not hold

The auditor recomputed every count independently, and they match. Three findings void the claim.

1. **Not held out.** Most scored tables existed before this pre-registration's commit:
   - HG002 version 2;
   - round 24 version 2, byte-identical;
   - the A119b k = 3 run, which already gives 21 families.

   Round 25 is the only fresh data: 7 qualifying families have reads there, and 3 work (CR3822851, DUXAP9, RABL2A).
2. **The rule does not single out O2.** It compared O2's precision on the reads O2 places with minimap2's on all reads. The
   baseline "minimap2's primary when MAPQ >= 1" passes the same rule in more families (simulation 30, real 33). In 46 of
   the 48 working cells every O2 call equals minimap2's primary.
3. **A floor of 5 placements cannot show <= 5% wrong.** With 0 wrong of 5, the 95% upper bound is 45%.

Smaller points:
- CD99, IL3RA and ASMTL are pseudoautosomal (phasing, not paralogs).
- A119b's 3 fragments per source read are not independent.
- The among-candidates rule (7bf68347) was added after the A119b runs.

**Pooled against the fair baseline (reads within 10 AS; minimap2's primary taken only at MAPQ >= 1):**

| data | reads (MAPQ-0) | O2 placed (wrong) | O2 on MAPQ-0 reads | O2 differs from minimap2 (right) | MAPQ >= 1 baseline placed (wrong) |
|---|---|---|---|---|---|
| simulation, rounds 23-25 | 12,585 (8,269) | 1,104 (0%) | 0 | 0 | 4,316 (0.14%) |
| HG002, donor truth | 27,199 (17,247) | 1,381 (1.09%) | 2 | 9 (8) | 8,958 (27.1%) |
| A119b fragments | 2,394 (1,108) | 319 (0.63%) | 4 | 5 (5) | 1,286 (13.8%) |

**Reading:**
- **O2 does not resolve exact ties (MAPQ 0).** Nothing in those reads separates the copies (0 / 2 / 4 placed), so they stay
  sets.
- **On near-ties O2 rarely overrules minimap2.** Its value is precision: on real reads, minimap2's MAPQ >= 1 near-tie calls
  are 14-27% wrong, and the subset O2 validates is about 1% wrong.
- **In one genome minimap2 is already right** (0.14% wrong), so O2 adds only certification and sets there.

The bars that would test this claim:
- held out on a donor never used here (the public testis library, ERR13885926);
- a minimum number of placements;
- O2's precision above the MAPQ >= 1 baseline's on the same reads;
- source reads counted, pseudoautosomal families apart.


---

## PREREG_o2_regimen_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_regimen_2026-10-19.md`)*

## Pre-registration: O2's regimen on almost-tied reads, in simulation and on real reads (round 24 + A119b fragments, 2026-10-19)

**User (2026-10-19, /goal):** "establish a working O2 regimen from both simulations and real reads, O2 must get almost tied
multi-mapping reads and help assign via PSVs from which copy did they actually come from".

### The regimen (`bench/o3_graph/o2psv.py svp`, code 26262631, defaults: k = 3, M = 10, fit 0.1%)

- **Domain (almost tied):** a read whose alignments within 10 AS of its best (Eichler) lie on two or more copies of one
  family and on no locus off it. Each alignment belongs to one copy: its strand, then exon overlap.
- **Evidence (PSVs):** per tied copy, the substitutions of its best tied alignment. Clips and indels are left out.
- **The read is placed** when either:
  - **two PSVs:** one copy has at least 2 fewer substitutions than every other; or
  - **one sample-validated PSV:** the read carries copy X's base at a PSV column, and the sample's own reads placed on the
    rival copy Y (AS margin >= M = 10 over every other copy: Eichler's well placed; the read itself left out) cover Y's
    position at least k = 3 times
    and none shows that base. The one copy no validated PSV excludes takes the read, provided the read's bases contradict
    it at no PSV column.
- **Otherwise** the read keeps the set of copies neither rule removes; disagreeing rules abstain. A read its best copy does
  not explain is O3's (`psv_unexplained`): more than 0.1% substitutions, HiFi's Q30 rate.

Why both rules:
- in one genome, almost-tied reads differ by at most about one PSV, which is enough there (round 23: 2,045 right, 0 wrong);
- on a real donor, a single PSV is the donor's own allele often enough to matter (HG002: one-PSV placements 93.4% right).

The sample-validated PSV is the bridge: it is used only where the sample's own copy shows it holds.

### Development (spent)

**HG002, donor truth:** all reads within 10 AS on >= 2 copies, poor_match unscored.

| HG002 setting | right | wrong | precision |
|---|---|---|---|
| **the regimen** | **3,169** | **35** | **98.9%** (37 own_tie placements) |
| two PSVs alone | 2,635 | 25 | 99.06% |
| validated alone | 534 | 10 | 98.1% |
| fit at the sample's own p90 / p99 instead of 0.1% | | | 98.3% / 97.9% (more placed) |
| validator margin 1 instead of 10 | | | the same |
| minimap2's primary on the same reads | 7,091 | 3,097 | 69.6% |

**Simulation:** three spent CHM13-only round-23 worlds.
- 124 placed, 0 wrong; every set holds the true copy (1,227 of 1,227); copies 2.01 -> 1.85-2.00.
- minimap2's primary is 55-64% right on the same reads.
- The fixed 0.1% fit costs power under the noisier public profile: 178 reads unexplained in one world.

### Held-out test A: simulation, one genome

- **Families:** round 24's draw: 21 fresh families, 90 copies (`families24.json`, set `r24`, sha1 8611cffce41b; the pool
  less round 23's).
- **Worlds:** INTACT, CROSS and SILENT, under `profiles/{A119b,public}_testis_v4_ref.json` (no donor alleles), seeds 1401,
  1402 and 1403. That makes 6 cells and 18 worlds.
- **Run:** `fastset.py run` with round 22's `$V17` builds the worlds. Then, per world, `o2psv.py copies-world` and
  `o2psv.py svp E0.bam E0.bam copies genome.ref.fa --k 3` (the defaults), scored against
  `reads.truth.E0.tsv`.
- **Bars (pooled):**
  - **A1:** wrong at most 0.5% of placements; floor 200 placements.
  - **A2:** sets hold the true copy for at least 99% of the domain's reads.
  - **A3:** O2's precision is above minimap2's primary's on the same reads.

### Held-out test B: real reads, A119b testis (a donor the regimen was not developed on)

- **Sources:** per pool family (`r19seq`, 439 families), 60 reads by a seeded hash (`o2frag.py sources`, seed 1). Each is
  a MAPQ-60 primary on a copy (strand and exon overlap), at least 1 kb long. Scripts:
  `fast/heldout24/run_sim.sh` and `o2frag_a119b/run_real.sh`.
- **Fragments:** 3 per source, 300-1,200 bp, seeded (`o2frag.py fragments`). Sources and fragments are realigned to whole
  CHM13 with A119b's own parameters: `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`.
- **Truth:** a fragment's truth is its source's copy, when the realigned full read lies on that copy at least 30 AS clear of
  any other locus (`o2frag.py truth`). On HG002, reads that clear agree with the donor's assembly 98.4% of the time, so
  this proxy carries 1-2% noise.
- **O2:** `o2psv.py svp FRAG.bam A119b.t2t.bam pool_copies target.fa --names truth.tsv --exclude sources.tsv` (a source
  never validates its own fragment), same k, M and quantile.
- **Bars (pooled):**
  - **B1:** wrong at most 2% of placements; floor 100.
  - **B2:** O2's precision is above minimap2's primary's on the same fragments.
  - **B3:** sets hold the true copy for at least 95%.

### Verdict

- **PASS:** A and B pass, each powered.
- **FAIL:** a powered bar fails.
- **UNDERPOWERED:** otherwise.

### Not claimed

- That identical copies can be split: they stay a set.
- Copies the reference lacks (O3).
- Truth beyond B's proxy. HG002's donor truth was development.

**Rules:** code frozen at 26262631. A bug is recorded as an amendment before any rerun.

### Results, test A (simulation, one genome; `fast/heldout24/run_sim.log`)

**Test A passes.**
- **A1:** 0 wrong of 349 placements (floor 200).
- **A2:** sets hold the true copy for 2,454 of 2,454 reads.
- **A3:** precision 1.000, against minimap2's primary's 0.789 on the same reads.

| cell | almost-tied reads | placed: right / wrong | sets holding the true copy | mean copies: tied -> set | minimap2 primary right |
|---|---|---|---|---|---|
| A119b s1401 | 569 | 102 / 0 | 569 / 569 | 4.80 -> 4.60 | 0.775 |
| A119b s1402 | 613 | 69 / 0 | 613 / 613 | 4.35 -> 4.18 | 0.856 |
| A119b s1403 | 438 | 68 / 0 | 438 / 438 | 5.52 -> 5.30 | 0.799 |
| public s1401 | 415 | 48 / 0 | 415 / 415 | 3.00 -> 2.88 | 0.677 |
| public s1402 | 217 | 25 / 0 | 217 / 217 | 3.10 -> 2.94 | 0.820 |
| public s1403 | 202 | 37 / 0 | 202 / 202 | 3.63 -> 3.45 | 0.802 |
| **pooled (18 worlds)** | **2,454** | **349 / 0** | **2,454 / 2,454** | **4.26 -> 4.09** | **0.789** |

**Reasons:** psv_validated 335, psv_two 14, psv_set 1,571, psv_unexplained 534.
- SILENT worlds place almost nothing, because a silent copy has no reads to validate a PSV against.
- AC0061583's 43 copies dominate the mean copies.
- The fixed 0.1% fit sets 534 reads apart as unexplained: the cost of the real-read calibration.

### Results, test B (real reads, A119b testis fragments; `o2frag_a119b/run_real.log`)

**Test B fails on B1, so the verdict is FAIL. Test A passes.**
- **B1 fails:** 7 wrong of 312 placements, 2.24% against the 2% bar (powered: floor 100).
- **B2 holds:** precision 0.978, against minimap2's primary's 0.656 on the same fragments.
- **B3 holds:** sets hold the true copy for 2,389 of 2,459 (97.2%).

| A119b fragments | value |
|---|---|
| sources: kept / drawn (clear by 30 AS after realignment) | 15,376 / 15,573, 397 families |
| fragments with a truth / almost-tied ones | 46,128 / 2,459 |
| placed by validated PSVs: right / wrong | 279 / 3 |
| placed by two PSVs: right / wrong | 26 / 4 |
| sets / unexplained (fit 0.1%) | 1,021 / 1,126 |
| mean copies, tied -> set | 2.34 -> 2.18 |

**The seven wrong placements:**
- In all seven, minimap2's primary for the fragment agrees with O2's call. The fragment's own PSVs point away from its full
  read's copy.
- Five are within one family:
  - TRIM64C, FAM153B and DHX40, by two PSVs;
  - TCAF2 (two fragments of one read), by validated PSVs.

  Two or more PSVs inside a short window disagree with the full read's placement: a local conversion tract in the donor's
  copy, or a heterozygous PSV that the 3 validating reads missed (all of them from the other haplotype has probability 1/8).
- Two fall in another family than their source read's (AC1423814 vs AC1378001, NOS2P2 vs MSTRG2518): the window is a
  segment the source shares with another family.

Unvalidated two-PSV placements were wrong 4 times in 30, validated ones 3 in 282.


---

## PREREG_o2_single_genome_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_single_genome_2026-10-19.md`)*

## Pre-registration: O2 in one genome -- the fewest-edit rule on CHM13-only worlds (round 23, 2026-10-19)

**User (2026-10-19):**
- "in O2 we only need gene copies, the multi-mapping reads they carry that cannot be resolved by minimap2 so they are given
  an arbitrary coin toss assignation of which one is the primary. In this case I just need to get a way for O2 to try to
  find what is the most likely real placement (or the real primary) or if there are five multimappers only two could be
  actually informative and the others an artifact, so that is why I would prefer that until we make O2 work we stay with
  only 1 genome, in this case CHM13".
- Then, after the proposal below: "yes please".

### The rule (code d63df788, `fastset.py ... --fewest-edits` = `o2graph.py run W intega 2 ... --fewest-edits`)

For a read tied within 10 AS (or by sequence) on two or more copies of its family:
- **O2's own call stands**, its structural checks included.
- **Where O2 abstained for want of its two-edit margin** (`thin`, `tied`, `locus_thin`, `locus_tied`), the copy with the
  fewest edits takes the read when it is the only one. Edits count every unaligned base (`veto_distance`).
- **Otherwise** the read keeps the set of its fewest-edit copies. The other placements explain it worse and are dropped.
- **Any other abstention** (a structural conflict, a call off the tied copies) keeps every tied copy.

Output: `o2.intege.m2.tsv` (call, reason, v17's call as `base_call`, `candidates`, `set`), beside an unchanged
`o2.intega.m2.tsv` (v17).

### One genome

The worlds simulate every read from its CHM13 copy with no donor alleles (`profiles/*_v4_ref.json`: v4 with
`het_per_base` = `hom_per_base` = 0). The reads keep the profile's errors, truncations, retained introns and adapters.

Only the single-genome world types are run:
- **INTACT:** every copy in the genome and the catalogue;
- **CROSS:** copies also express paralogs' isoform structures;
- **SILENT:** a third of the copies express nothing.

REMOVED, OUTSIDE and BRANCH involve a copy the reference or catalogue lacks, which is O3's question.

**Development (spent: round 22's families, A119b profile, seed 1201, the three worlds; `fast/ref_dev`):**
- Of 846 reads tied in their own family, the rule places **428 right, 0 wrong**, against v17's 109 / 0.
- Every set holds the true copy (846 / 846), and the mean set is 1.65 copies against 2.46 tied.
- minimap2's primary is wrong on 232 of the 846.

### Families (held out; drawn by the fixed rule)

**The pool:** round 22's unspent pool (98) less every family sharing a gene with round 22's draw, leaving 72 families with
213 copies (`r23seq/families23seq_all.json`). None of them was run by the sample-copy step on real reads.

**The draw:** the rule of rounds 17-22 with seed 23 (`seqfamilies.py draw`): **37 families, 90 copies**.
- **Copies per family:** {2: 32, 3: 3, 4: 1, 13: 1}.
- **Copy biotypes:** 43 protein-coding, 22 transcribed unprocessed, 15 transcribed processed, 5 unprocessed, 2 processed,
  2 lncRNA, 1 polymorphic.
- **File:** `/mnt/linuxdisk/tmp/o3_graph/families23.json` (set `r23`, sha1 0e46c2d3b8d8).

The families: AC1387761, AL5913791, AP0005524, AP0065871, CASP1, CKMT1B, CTRB2, CXorf51B, CYB5B, CYP4Z1, FAM153B, KRT16,
LRRC37B, MAGEA9B, MRPL36, MST1, MSTRG5154, NOS2P2, NPY4R, OR1D5, OR4C13, OR5M10, OR5T3, OTOAP1, PRDM7, PRSS40B, REXO1L1P (13
copies), RPL38, RPS27AP5, RPS4Y1, SETSIP, SMPD4, THOC3, TPTE2P2, TRIM77, ZNF33A, ZNF852 (`S19` prefix).

### Runs

- **Cells:** 2 profiles (A119b, public) x 3 seeds (1301, 1302, 1303), with no seed selection.
- **Worlds:** 3 per cell, so 18 in all.

```
F=/mnt/linuxdisk/tmp/o3_graph/families23.json; O=/mnt/linuxdisk/tmp/o3_graph/fast/heldout23
for S in 1301 1302 1303; do for P in A119b public; do
  python3 bench/o3_graph/fastset.py run $O/${P}_s$S $F r23 $S --local $O/local $V17 --fewest-edits \
    --sim-profile bench/o3_graph/profiles/${P}_testis_v4_ref.json --worlds INTACT,CROSS,SILENT
done; done
python3 bench/o3_graph/summary.py sets $O/sets.json $O/*_s13*/{INTACT,CROSS,SILENT} --off-family-ties 9 --seq-tie
```

`$V17` is round 22's flag set, unchanged.

### Bars (pooled over the 18 worlds; reads tied in their own family, off-family ties apart)

- **F1, precise:** the rule's wrong placements are at most 0.5% of its placements. Floor: 200 placements; below it,
  UNDERPOWERED.
- **F2, sets keep the truth:** the set holds the true copy for at least 99% of these reads. A placement counts as a set of
  one.
- **F3, more placed:** the rule places more of these reads than v17.
- **F4, no harm end to end:** fastset's verdict, now scored on the rule's table, is PASS in every cell (B1 and B2 in every
  world).

**Verdict:**
- **PASS:** F1 is powered, and F1 to F4 hold.
- **FAIL:** a powered bar fails.
- **UNDERPOWERED:** otherwise.

**Reported:**
- per cell and world: placements right / wrong, v17's, the sets' sizes against the tied copies, and minimap2's primary
  right / wrong;
- the reasons in `o2.intege`.

**Rules:** code frozen at d63df788. A bug is recorded as an amendment before any rerun.

### What a PASS would and would not claim

**A PASS would claim:** in one genome (reads from CHM13's own copies), O2 places most of the reads minimap2 leaves within
10 AS on their true copy at no more than 0.5% error. For the rest, it returns the smaller set of equally explained copies,
and that set holds the true one.

**It would not claim:**
- real samples, whose copies differ from CHM13's (O3: on HG002 a one-edit margin picked the wrong copy 46.5% of the time);
- copies the reference or catalogue lacks;
- that identical copies can be split. No rule can split them; the set says so.

### Amendment 1 (2026-10-19, after the run and before rescoring: a fresh code review of d63df788)

The review found no defect in the rule. It found two defects in the scorer's population and one reported number, plus two
minor gaps.

- **minimap2's primary was credited to the first copy whose span it overlaps** (`summary.primary_copies`). Two families
  have overlapping copies:
  - CYB5B's single-exon AC026464.2 lies inside CYB5B's last exon;
  - RPL38's antisense MGC16275 overlaps RPL38's end.

  **Fix:** the copy on the alignment's strand whose exons its blocks overlap most (`summary.attribute`).
- **A read minimap2 placed once counted as tied** when its one alignment spans both overlapping copies. That is
  `o2graph.best_as_loci` crediting an alignment to every span it overlaps, as v17 always has. Such a read is no
  multimapper.

  **Fix:** reads with a single tied alignment (`summary.single_tied_alignment`) are set apart, like the off-family ties.
  The rule's table is unchanged, so this changes the scored population only.
- **Not changed** (minor, and recorded):
  - `insignificant` (`--o2-significance`) is a margin-type abstention that FILL leaves out. It occurred 0 times in rounds
    22-23. It is kept out on purpose: that gate prices each edit type by its error rate, and a raw one-edit count would
    override it.
  - A set's `dropped` entry stands for every copy of the family O1 dropped. It is exposed on 52 of 11,609 tied rows and
    changes no `set holds truth`.

The bars and the run are as frozen. The same tables are rescored with the fixed scorer (`sets_amd1.json`). The original
scoring stays below for the record.

### Results under Amendment 1 (the verdict; `fast/heldout23/sets_amd1.json`, `sets_amd1.txt`)

**PASS: F1 to F4 hold.**
- **F1:** 0 wrong of 2,045 placements.
- **F2:** sets hold the true copy for 7,450 of 7,450 multimapped reads.
- **F3:** 2,045 placed, against v17's 430 (4.8 times).
- **F4:** fastset's verdict is PASS in all 6 cells; this is unaffected by the scorer.

| cell | multimapped tied reads | rule: right / wrong | v17: right / wrong | sets holding the true copy | mean copies: tied -> set | minimap2 primary wrong |
|---|---|---|---|---|---|---|
| A119b s1301 | 1,290 | 321 / 0 | 53 / 0 | 1,290 | 2.01 -> 1.75 | 483 (37.4%) |
| A119b s1302 | 1,231 | 398 / 0 | 114 / 0 | 1,231 | 2.05 -> 1.68 | 424 (34.4%) |
| A119b s1303 | 1,329 | 345 / 0 | 36 / 0 | 1,329 | 2.01 -> 1.74 | 525 (39.5%) |
| public s1301 | 1,091 | 237 / 0 | 23 / 0 | 1,091 | 2.01 -> 1.78 | 402 (36.8%) |
| public s1302 | 1,229 | 432 / 0 | 141 / 0 | 1,229 | 2.00 -> 1.65 | 415 (33.8%) |
| public s1303 | 1,280 | 312 / 0 | 63 / 0 | 1,280 | 2.00 -> 1.76 | 474 (37.0%) |
| **pooled** | **7,450** | **2,045 / 0** | **430 / 0** | **7,450** | **2.02 -> 1.73** | **2,723 (36.6%)** |

The 910 reads set apart were placed once. 529 of v17's 959 original placements were among them, so on true multimappers
v17 places 430.

### Results as first scored (superseded by Amendment 1; kept for the record; `sets.json`, `sets.txt`)

**PASS: F1 to F4 hold.**
- **F1:** 0 wrong of 2,574 placements (limit 0.5%, floor 200).
- **F2:** sets hold the true copy for 8,360 of 8,360 reads.
- **F3:** 2,574 placed, against v17's 959.
- **F4:** fastset's verdict is PASS in all 6 cells.

| cell | tied reads | rule: right / wrong | v17: right / wrong | sets holding the true copy | mean copies: tied -> set | minimap2 primary wrong |
|---|---|---|---|---|---|---|
| A119b s1301 | 1,398 | 421 / 0 | 153 / 0 | 1,398 | 2.01 -> 1.70 | 519 (37.1%) |
| A119b s1302 | 1,583 | 592 / 0 | 308 / 0 | 1,583 | 2.04 -> 1.63 | 598 (37.8%) |
| A119b s1303 | 1,477 | 404 / 0 | 95 / 0 | 1,477 | 2.01 -> 1.73 | 621 (42.0%) |
| public s1301 | 1,149 | 285 / 0 | 71 / 0 | 1,149 | 2.01 -> 1.75 | 402 (35.0%) |
| public s1302 | 1,429 | 526 / 0 | 235 / 0 | 1,429 | 2.00 -> 1.63 | 445 (31.1%) |
| public s1303 | 1,324 | 346 / 0 | 97 / 0 | 1,324 | 2.00 -> 1.74 | 480 (36.3%) |
| **pooled** | **8,360** | **2,574 / 0** | **959 / 0** | **8,360** | **2.01 -> 1.69** | **3,065 (36.7%)** |

**Reading:**
- On reads minimap2 leaves within 10 AS, its primary is wrong 36.7% of the time.
- The rule places 31% of them (2.7 times v17), all on their true copy.
- The other 69% keep a set of equally explained copies, and the set always holds the true copy. Most tied sets are pairs
  (mean 2.01), so a set of two means the two copies are identical over the read; no read placement can split them.

**Reasons on every tied row (all families, all reads):**
- placed: fewest_edits 1,847, O2's own assigned 962 and short_path 20 -- all 2,829 on the true copy, in the read's own
  family;
- kept as a set: fewest_edit_set 6,042;
- every tied copy kept: structure_only 58, unexplained 15.

**End to end:** the only integrated error is one read in two public_s1303 worlds (s0000683, FAM153B). minimap2 resolved it
by 10 AS or more onto the wrong copy (`as_resolved`), outside O2's domain, and v17 makes the same call.


---

## PREREG_samplecopy_haploid_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_samplecopy_haploid_2026-10-19.md`)*

## Pre-registration: the haploid rule of the O3 -> O2 sample-copy step, held out on HG002's chrX families (2026-10-19)

**User (2026-10-19):** "yes lets implement the rule", after register 1336. On A119b testis (male), two RBMY1 clusters were
labelled with one RBMY1F copy, and CDY1's two copies formed one cluster. On a haploid chromosome that cannot be one copy
per cluster, yet every agreed YAG call followed such a label.

### The rule (code 2adc91d7, `samplecopy.py run ... --male`)

A male sample carries its chrX and chrY copies once outside the pseudoautosomal regions (`hg002truth.PAR`), so each of
those copies gives one group of reads. A read tied only on such copies is placed only when each copy it ties on labels
**exactly one** cluster (`one_copy_one_group`):
- **two or more clusters** on one copy: at least one of them is not that copy;
- **no cluster**: the copy's reads may sit unseen inside another copy's cluster.

Reads tied on any diploid copy are untouched, so autosomal and pseudoautosomal families are unchanged. That was checked
byte-identical on PDPR, and on CDY2B without `--male`. The flags name, for O3, the haploid copies that label several
clusters or none.

**Development (A119b YAGs, spent):** the rule blocks all 208 placements (RBMY1 134, CDY 59, DAZ 15); split-half then calls
nothing. On HG002's two haploid families that have donor-copy reads (IKBKG 47, DMRTC1 5; spent, register 1333), the step
placed none, so the rule has nothing to remove there.

### The held-out families (by a fixed rule, never run by the step)

Every HG002 family whose copies are all haploid in a male, with at least 5 tied reads in `hg002proof/tied.tsv`, other
than the spent IKBKG and DMRTC1. Tied = equal best AS, the tables of registers 1329-1336.

| family | copies | tied | own_tie | elsewhere | poor_match | untied |
|---|---|---|---|---|---|---|
| F8A2 | F8A2, F8A3 (chrX) | 2,335 | 1,043 | 918 | 374 | 5 |
| FAM156B | FAM156B, FAM156A (chrX) | 634 | 532 | 0 | 102 | 36 |
| MAGED4B | MAGED4B, MAGED4 (chrX) | 50 | 39 | 0 | 11 | 2 |
| H2AB2 | H2AB2, H2AB3 (chrX) | 6 | 4 | 2 | 0 | 3 |

**No tied read here has a donor copy in its family** (no `unique`):
- own_tie reads tie in the donor too, so placing one is unsupported;
- elsewhere reads come from a donor copy off the family, so placing one is wrong;
- poor_match reads match the donor poorly and are reported but not scored.

The scored denominator is **2,538** (1,618 own_tie + 920 elsewhere).

### The run

```
H=/mnt/linuxdisk/tmp/o3_graph/hg002proof; O=/mnt/linuxdisk/tmp/o3_graph/samplecopy_haploid
for F in S19F8A2 S19FAM156B S19MAGED4B S19H2AB2; do
  python3 bench/o3_graph/samplecopy.py run /mnt/linuxdisk/tmp/o3_graph/r19seq/families19cal_full2_all.json seq19 $F $H/chm13.bam \
    $H/tied.tsv $H/untied.tsv /mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa $O/plain
  python3 bench/o3_graph/samplecopy.py run <same arguments, output $O/male> --male
done
python3 bench/o3_graph/samplecopy.py score $H/truth.tsv $O/plain <families>; ... score $H/truth.tsv $O/male <families>
```

### Bars

- **HS, calibrated on haploid families (with the rule):** own_tie placed + elsewhere placed <= 0.05 x 2,538, so at most
  126 reads.
- **HB, the rule's benefit:** the step without the rule places at least 50 own_tie + elsewhere reads (the floor), and
  the rule places fewer.

**Verdict:**
- **PASS:** HS and HB hold.
- **UNDERPOWERED:** HS holds, but the step without the rule places fewer than 50 such reads, so there was little to remove.
- **FAIL:** HS fails.

**Reported:** per family, placements by class with and without the rule, the clusters (labels, untied counts) and the
O3 flags.

**Rules:** code frozen at 2adc91d7. A bug is recorded as an amendment before any rerun.

### What a PASS would and would not claim

**A PASS would claim:** on real reads where the donor itself cannot place them, the haploid rule stops the sample-copy
step from placing reads that the step without it places.

**It would not claim:**
- that the rule keeps right placements: these families have none to keep. The no-harm side rests on construction, since
  diploid families are untouched;
- anything about the YAGs, which have no truth;
- anything about female samples.

### Results (2026-10-19, run as frozen; `samplecopy_haploid/run.log`, `score_plain.txt`, `score_male.txt`)

**UNDERPOWERED.**
- **HS holds:** 0 own_tie or elsewhere reads are placed with the rule, against a limit of 126.
- **HB is not reached:** the step without the rule also places **none**, below the floor of 50, so the rule had nothing
  to remove.

| family | reads in the model | informative columns | clusters (labelled) | tied reads placed, without / with the rule |
|---|---|---|---|---|
| F8A2 | 1,568 | 1 | 1 (0) | 0 / 0 of 1,555 |
| FAM156B | 675 | 0 | 1 (1) | 0 / 0 of 634 |
| MAGED4B | 52 | 0 | 1 (0) | 0 / 0 of 50 |
| H2AB2 | 1,416 | 0 | 1 (0) | 0 / 0 of 5 |

The donor's copies carry no variation that separates them in these reads (0-1 informative columns). Each family is one
cluster, and the step abstains whatever the rule. That is the right outcome for reads the donor itself ties (own_tie
1,283 scored), and for reads from a copy off the family (elsewhere 596). It also means these families cannot show what
the rule adds.

The rule therefore rests on two things:
- **construction:** it only turns placements into abstentions, and diploid families are untouched (byte-identical,
  PDPR and CDY2B);
- **the inconsistency it detects on the A119b YAGs:** all 208 placements blocked.

No real dataset in hand has truth for a haploid family where the step places reads.


---

## PREREG_samplecopy_hg002_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_samplecopy_hg002_2026-10-19.md`)*

## Pre-registration: the sample's own copies place O2's ambiguous reads (HG002, held-out families, 2026-10-19)

**User (2026-10-19):**
- "can we show clear cases of O2 working?"
- "Can we try to find families with 95-99% similarity members, those should have tied AS reads"
- "lets stop this and focus on building on those families that are more likely to work"
- On the design: "yes".

### Why

**Where ties live.** On HG002's reads, AS-tied reads are 0.0-0.1% of reads in families whose closest copies are 95-99%
identical, 3.5% at 99.5-99.9%, and 28% at 99.9% or more.

**Why O2 placed so few.** The donor-truth test (register 1329) found that most tied reads carry the donor's own alleles,
which both CHM13 copies lack. O2's PSVs, taken from CHM13, cannot place them.

**What `samplecopy.py` does instead.** It clusters each family's reads by the alleles they share, read-backed phasing over
the sample's own copies, after SDA / Chaisson 2017 / IsoCon in the thesis vault. A cluster takes the copy its untied reads
resolve to, by a confident majority. A tied read takes its cluster's copy when it carries two decisive sites; otherwise it
abstains.

### Development (spent): TBC1D3B and PDPR

- **PDPR:** 14 right, 0 wrong of 60 tied reads with a donor copy. v14 places none.
- **TBC1D3B:** abstains on all 97. Its reads hold no within-sample variation that separates the copies: only 4 informative
  columns, and the donor separates its tied reads at 2 sites that every one of its reads shares.
- **Rules fixed in development** (each was a development finding):
  - The error rate is the reads' background mismatch rate to the anchor.
  - The allele floor is half of one copy-haplotype's expected share.
  - A tied read needs two decisive sites.
  - Labels need a confident majority: the 95% Wilson lower bound of the share above 1/2.

### The held-out families (by a fixed rule, before any held-out run)

Every HG002 family with at least 5 AS-tied reads that have a donor copy (`hg002proof/truth.tsv`: `unique` or
`unique_elsewhere`), not pseudoautosomal, with no chrM copy, and other than the two development families.

**19 families, 399 tied reads with a donor copy, 7,228 donor-tied (`own_tie`):** AC0064533 14, AC1275021 6, ANKRD36C 29,
DMRTC1 5, DUSP22 7, EIF3C 24, GOLGA8A 9, GTF2H2C 5, GTF2IRD2 35, HSPA6 34, IKBKG 47, MSTRG3150 60, NBPF1 27, SERF1B 13,
SLC35E2B 18, SMG1 5, SMN2 22, SPDYE1 34, SULT1A1 5.

### The run

```
H=/mnt/linuxdisk/tmp/o3_graph/hg002proof; O=/mnt/linuxdisk/tmp/o3_graph/samplecopy_heldout
for F in <the 19 families>; do
  python3 bench/o3_graph/samplecopy.py run /mnt/linuxdisk/tmp/o3_graph/r19seq/families19cal_full2_all.json seq19 $F $H/chm13.bam \
    $H/tied.tsv $H/untied.tsv /mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa $O
done
python3 bench/o3_graph/samplecopy.py combine $H/truth.tsv $O /mnt/linuxdisk/tmp/o3_graph/fast/hg002/nonpar/INTACT <families> $O/combine.json
```

- **Truth:** the donor's own haplotypes (register 1329, `hg002truth.py`).
- **v14:** the HG002 world's `o2.integv.m2.tsv`.
- **The combined system:** v14's placement stands, and the sample copies fill v14's abstentions.

### Bars (pooled over the 19 families)

- **S1, precise:** the sample copies' wrong / placed <= 0.02 on tied reads with a donor copy. Floor: 50 placed; below it,
  UNDERPOWERED.
- **S2, calibrated:** the sample copies place <= 0.05 of the `own_tie` reads.
- **S3, adds without harm:** v14 + sample copies has more right placements than v14 and no more wrong ones.

**Verdict:**
- **PASS:** S1 is powered and holds, and S2 and S3 hold.
- **FAIL:** a powered bar fails.
- **UNDERPOWERED:** otherwise.

**Reported:** per family, recall (right / tied reads with a donor copy), and the clusters (columns, sizes, labels).

**Rules:** the code is frozen at the commit that adds this file. A bug is recorded as an amendment before any rerun.

### What a PASS would and would not claim

**A PASS would claim:** on real reads of 19 held-out families, clustering the sample's reads by their own alleles places
AS-tied reads that CHM13's PSVs cannot place, at >= 98% precision. It abstains on donor-tied reads, and it adds to v14
without errors.

**It would not claim:** testis-only families (the YAGs, absent from this library), power where the sample's reads carry no
separating variation (TBC1D3B), or other donors.

### Amendment 1 (2026-10-19, 13:05, during the run and before any score was computed)

**The bug.** The scorer counted only the tied reads present in `<family>.calls.tsv`. A tied read that never aligned to the
anchor copy has no row there and so dropped out of every count. MSTRG3150, for example, has 94 tied reads and 30 rows.

**The fix.** `samplecopy.py combine ... --tied hg002proof/tied.tsv` scores every tied read of the family that has a truth
row. A read with no call counts as an abstention, for the sample copies and in v14 + sample copies alike.

**What does not change.** No call. S1 (wrong / placed) is unaffected. S2's and S3's denominators now cover all the
family's tied reads. The pre-registered score is the `--tied` one.

### Results (2026-10-19, run once; `/mnt/linuxdisk/tmp/o3_graph/samplecopy_heldout/`, scored with Amendment 1's `--tied`)

**Verdict: FAIL.** S1 and S3 are powered and fail; S2 holds.
- **S1:** the sample copies are wrong on 4 of 86 placements (0.047) against the 0.02 bar.
- **S2:** they place 7 of the 7,228 `own_tie` reads (0.001).
- **S3:** v14 + sample copies makes 126 right placements against v14's 48, but 11 wrong against 7.

| pooled, 19 families (399 tied reads with a donor copy) | right | wrong | abstain | own_tie placed |
|---|---|---|---|---|
| sample copies | 82 | 4 | 313 | 7 |
| v14 | 48 | 7 | 344 | 52 |
| v14 + sample copies | 126 | 11 | 262 | 59 |

**Clear cases** (sample copies against v14, on the same tied reads):

| family | sample copies right / wrong | v14 right / wrong |
|---|---|---|
| HSPA6 | **34 / 0** | 0 / 0 |
| GTF2IRD2 | **31 / 0** | 0 / 0 |
| MSTRG3150 | 7 / 0 | 0 / 0 |
| GTF2H2C | 4 / 0 | 3 / 0 |
| SMG1 | 2 / 0 | 0 / 0 |
| AC1275021 | 2 / 0 | 0 / 2 |

The development family PDPR adds 14 / 0 against v14's 0 / 0.

**Other families:**
- **The 4 errors:** GOLGA8A 3 and AC0064533 1.
- **Abstained on every tied read:** IKBKG (no informative column; v14 is right on 42 by structure), SMN2, SPDYE1, SERF1B,
  SULT1A1, NBPF1, DUSP22, SLC35E2B, ANKRD36C, DMRTC1 and EIF3C.

**The errors, read post hoc.** Each sits in a cluster labelled by a majority that is not pure:
- **GOLGA8A:** 829 of 907 untied reads (91%) and 152 of 162 (94%).
- **AC0064533:** 35 of 37.
- **Their donor margins:** 5-15, i.e. one to three sites.

A tied read assigned by its cluster inherits that cluster's impurity. The majority test (Wilson lower bound above 1/2) admits
clusters that are 6-9% another copy.

**Correction (same day).** Purity does not separate the errors from the clear cases. The clear cases' clusters are no purer
by their untied reads:

| family | untied-read purity of its labelled clusters |
|---|---|
| HSPA6 | 0.61 and 1.00 |
| GTF2IRD2 | 0.94 and 0.74 |
| PDPR | 0.87 |
| GOLGA8A | 0.91 and 0.94 |

So a purity bar would remove the clear cases with the errors.

**What does separate them.** Every untied read teaches the labels, including near-ties that minimap2 places at random
(46-49% wrong at an AS margin under 10, register 1330). The candidate fix is to teach labels only with untied reads resolved
by 10 AS or more, Eichler's rule. It needs its own held-out test.


---

## PREREG_samplecopy_par_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_samplecopy_par_2026-10-19.md`)*

## Pre-registration: the O3 -> O2 sample-copy step, a second held-out test (HG002 pseudoautosomal families, 2026-10-19)

**User (2026-10-19):** "Second held-out test (Recommended)", on the plan to run the frozen original rule with the same bars
on the six unused pseudoautosomal families.

### Why

**The first held-out test failed** (register 1333). On 19 families, 82 placements were right and 4 wrong (0.047 against the
0.02 bar), with clear wins in six families.

**Development found no rule that separates the 4 errors from the clear cases:**
- cluster purity (the register 1333 correction);
- teaching labels only with reads resolved by >= 10 AS (register 1334).

So the original rule is tested again, unchanged, for an independent estimate of its precision.

**The framing** (user, the same day): an O3 step feeding O2. Reads that CHM13 cannot explain are clustered by their own
variants and handed back to a CHM13 copy, or left as candidate copies CHM13 lacks. Phased data is never an input; HG002's
haplotypes are the truth only.

### The method

`bench/o3_graph/samplecopy.py` as in the first held-out test (2f23fe47), run without `--margins`. Its rule:
- the anchor is the copy most untied reads resolve to;
- the allele floor is 0.5 / kmax;
- the error rate is the reads' background mismatch rate;
- tied reads need two decisive sites and a posterior of at least 0.99 on their own evidence;
- labels need a confident majority of all untied reads.

### The held-out families (by the same fixed rule, never run by this method)

The pseudoautosomal HG002 families with at least 5 AS-tied reads that have a donor copy.

**6 families, 622 tied reads with a donor copy, 6,681 donor-tied:** AKAP17A 10, CD99 8, CD99P1 24, DHRSX 544, PLCXD1 20 and
SLC25A6 16.

**Here the copies are chrX and chrY,** HG002's maternal and paternal alleles of one recombining locus. The clusters are
those two haplotypes. CHM13 v2.0's chrY is HG002's own Y.

### The run

```
H=/mnt/linuxdisk/tmp/o3_graph/hg002proof; O=/mnt/linuxdisk/tmp/o3_graph/samplecopy_par
for F in S19AKAP17A S19CD99 S19CD99P1 S19DHRSX S19PLCXD1 S19SLC25A6; do
  python3 bench/o3_graph/samplecopy.py run /mnt/linuxdisk/tmp/o3_graph/r19seq/families19cal_full2_all.json seq19 $F $H/chm13.bam \
    $H/tied.tsv $H/untied.tsv /mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa $O
done
python3 bench/o3_graph/samplecopy.py combine $H/truth.tsv $O /mnt/linuxdisk/tmp/o3_graph/fast/hg002/par/INTACT <families> $O/combine.json --tied $H/tied.tsv
```

### Bars (pooled; as the first held-out test)

- **S1, precise:** wrong / placed <= 0.02. Floor: 50 placed.
- **S2, calibrated:** `own_tie` placed <= 0.05.
- **S3, adds without harm:** v14 + sample copies has more right placements than v14 and no more wrong ones. v14 is the PAR
  world's `o2.integv.m2.tsv`.

**Verdict:**
- **PASS:** S1 is powered and holds, and S2 and S3 hold.
- **FAIL:** a powered bar fails.
- **UNDERPOWERED:** otherwise.

**Rules:** the code is frozen at the commit that adds this file. A bug is recorded as an amendment before any rerun.

### Results (2026-10-19, run once; `/mnt/linuxdisk/tmp/o3_graph/samplecopy_par/`)

**Verdict: PASS.** S1, S2 and S3 hold, and S1 is powered.

| pooled, 6 families (622 tied reads with a donor copy) | right | wrong | abstain | own_tie placed |
|---|---|---|---|---|
| sample copies | **508** | **0** | 114 | 0 / 6,681 |
| v14 | 1 | 0 | 621 | 0 |
| v14 + sample copies | 509 | 0 | 113 | 0 |

- **S1:** 0 of 508 placements are wrong.
- **S2:** 0 of the 6,681 donor-tied reads are placed.
- **S3:** v14 + sample copies makes 509 right placements against v14's 1, with no wrong one.

**One family carries it.**
- **DHRSX:** 852 reads, 41 informative columns and 2 clusters, both labelled. 508 of its 544 tied reads with a donor copy
  are placed right; v14 places 1.
- **The other five abstain throughout:**
  - AKAP17A, CD99 and SLC25A6 have a single cluster (SLC25A6 has 4 informative columns over 7,535 reads, CD99 has 35);
  - CD99P1's clusters are unlabelled;
  - PLCXD1's 4 labelled clusters place nothing (none of its tied reads passed the two-site and posterior rules).

**What it claims.** On a held-out stratum of real reads, the O3 -> O2 step resolves AS-tied reads without error when the
sample's own variants separate the copies. Here the copies are HG002's chrX and chrY alleles of DHRSX. Where nothing
separates them, it abstains.

**Limits.**
- The pass rests on one family.
- Separating chrX from chrY in the PAR is read-backed phasing of one recombining locus, not paralog assignment. It tests the
  method's machinery: clustering by the sample's own alleles and labelling by confidently placed reads.
- Together with register 1333 (82 / 4 on paralog families), the step's real-data precision is 590 right / 4 wrong over the
  two held-out tests.


---

## PREREG_samplecopy_yag_testis_2026-10-19.md

*(original: `docs/archive/2026-10/merged/PREREG_samplecopy_yag_testis_2026-10-19.md`)*

## Pre-registration: the O3 -> O2 sample-copy step on the Y-ampliconic families in testis, by split-half agreement (2026-10-19)

**User (2026-10-19):** "ok lets do the next step". That step is DAZ and the other YAGs in testis, checked by split-half,
since no donor assembly exists.

### Why split-half, and what HG002 says it means

There is no truth for a testis donor's Y. The check is reproducibility, using `samplecopy.py splithalf`:
- the family's reads are split into two disjoint halves by a seeded hash (seed 1);
- each half learns its own clusters and labels;
- every AS-tied read is called by both models.

**Calibration on HG002** (spent families, the donor's haplotypes as truth; `fast/../samplecopy_split_cal/`), tied reads with
a donor copy:
- **Both halves call the same copy:** 519 right, 0 wrong.
- **Only one half calls:** 38 right, 3 wrong.
- **The halves disagree:** never.

So agreement is the calibrated signal. The calibration comes from HG002's families, not from YAGs, and that limit is stated
below.

### Data

- **Reads:** A119b testis (`winloci_data/A119b.t2t.bam`, CHM13 v2.0, minimap2 splice:hq `-N 50 -p 0.1`). They are the
  4,953 reads whose primary lies within 1 kb of a copy of the 10 Y-ampliconic families of the adopted pool, with all their
  chrX / chrY alignments (`yag_testis/yag.bam`).
- **Limit:** secondaries off chrX / chrY are not seen. DAZL and the TSPY-like genes are autosomal but divergent.
- **Tables:** `hg002truth.py tied` / `untied`. 3,270 reads are AS-tied within one family and 1,021 are resolved by minimap2.

| family | copies | tied | untied |
|---|---|---|---|
| HSFY1 | 2 | 1,797 | 30 |
| RBMY1A1 | 6 | 590 | 293 |
| DAZ1 | 4 | 579 | 346 |
| TSPY2 | 48 | 136 | 123 |
| CDY2B | 4 | 105 | 84 |
| BPY2 | 3 | 47 | 5 |
| PRY | 4 | 10 | 11 |
| XKRY / VCX3A / RBMY2BP | 2 / 6 / 3 | 3 / 2 / 1 | 0 / 123 / 4 |

### The run (the frozen original rule; no change since be1edda6)

```
Y=/mnt/linuxdisk/tmp/o3_graph/yag_testis
for F in S19HSFY1 S19RBMY1A1 S19DAZ1 S19TSPY2 S19CDY2B S19BPY2 S19PRY; do
  python3 bench/o3_graph/samplecopy.py run /mnt/linuxdisk/tmp/o3_graph/r19seq/families19cal_full2_all.json seq19 $F $Y/yag.bam \
    $Y/tied.tsv $Y/untied.tsv /mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa $Y/out
  python3 bench/o3_graph/samplecopy.py splithalf <same arguments>
done
```

The families with at least 10 tied reads are run. XKRY, VCX3A and RBMY2BP have 1-3 each and are left out.

### Bar and reports

- **Y1, reproducible:** pooled over the families, the tied reads both half-models call agree on >= 0.95 of them. Floor: 50
  such reads; below it, UNDERPOWERED.
- **Reported:**
  - per family, the full model's placements, both-halves calls and agreements, and the agreed calls per copy (DAZ1-4, ...);
  - the clusters (labels, unlabelled clusters as candidate copies CHM13 lacks: O3's flag).
- **The claim (if Y1 holds):** testis tied reads of the YAGs are placed on CHM13 copies reproducibly from the sample's own
  variants.
- **Not claimed:** correctness. The HG002 calibration suggests agreed calls are rarely wrong, but YAG copies, testis
  expression and a donor whose Y is not CHM13's are untested by it.

### Results (2026-10-19, run as frozen; `yag_testis/run.log`, `out/*.calls.tsv`, `*.split.tsv`, `*.clusters.json`)

**Y1 holds: 187 of 187 tied reads called by both half-models agree (1.000), above the floor of 50.** The claim the bar allows
stands: the YAG tied reads that are placed are placed reproducibly. Correctness is not claimed, and the post-hoc reading
below gives a reason to doubt that these placements are provenance.

| family | reads / informative columns | clusters (labelled) | tied reads | full model placed | A / B / both / agree |
|---|---|---|---|---|---|
| HSFY1 | 2,191 / 6 | 1 (1) | 1,797 | 0 | 0 / 0 / 0 / 0 |
| RBMY1A1 | 942 / 101 | 7 (5) | 589 * | 134 | 163 / 218 / 134 / 134 |
| DAZ1 | 1,053 / 326 | 3 (1) | 579 | 15 | 0 / 19 / 0 / 0 |
| TSPY2 | 263 / 275 | 12 (6) | 136 | 0 | 0 / 0 / 0 / 0 |
| CDY2B | 207 / 38 | 2 (1) | 105 | 59 | 59 / 53 / 53 / 53 |
| BPY2 | 57 / 0 | 1 (0) | 47 | 0 | 0 / 0 / 0 / 0 |
| PRY | 61 / 214 | 2 (0) | 10 | 0 | 0 / 0 / 0 / 0 |

\* 590 in the table above; one read has only 73 bp inside the copy window (its other 1,442 bp lie in a downstream exon)
and does not reach the model: an abstention.

**The agreed calls, by copy:** RBMY1F@chrY:23,239k 132, RBMY1HP 2, CDY1@chrY:26,434k 53. DAZ1-4, HSFY1/2, TSPY, BPY2 and PRY: none.

**What each family's reads show** (untied reads = the copy minimap2's unique best AS picks):
- **HSFY1/2:** every read ties on the pair (1,797 of 1,827 reads), and the reads have 6 variable columns and form one group.
  The two copies are the same sequence over what is expressed. Nothing can be assigned, and nothing is.
- **DAZ1-4** (DAZ1 / DAZ4 and DAZ2 / DAZ3 are inverted pairs): 393 of the 579 tied reads tie on all four copies. The reads
  form three groups:
  - a DAZ1 group, labelled (untied DAZ1 79, DAZ4 6, DAZ2 3, DAZ3 2);
  - a DAZ2 + DAZ3 group, unlabelled (untied DAZ3 78, DAZ2 76): the two copies are not separable by this sample's reads;
  - a mixed group, unlabelled (DAZ1 32, DAZ3 27, DAZ4 18, DAZ2 11).
  The full model places 15 tied reads on DAZ1, and no read is placed by both halves. **DAZ is not resolved.**
- **RBMY1:** 7 groups. Five are labelled: RBMY1HP, RBMY1E, RBMY1A1, and **two groups (115 and 78 reads) both labelled
  RBMY1F@23,239k**. Two groups are unlabelled (73 and 44 reads). All 132 RBMY1F calls are reads tied between the two RBMY1F
  copies (23,004k on the minus strand, 23,239k on the plus strand). Of the family's untied reads, 62 resolve to 23,239k and
  8 to 23,004k.
- **CDY:** two groups.
  - A CDY1 group, labelled 26,434k by 76 of its 84 untied reads; the other 8 resolve to the other CDY1 copy (24,814k).
    All 53 agreed calls are CDY1-pair tied reads placed on 26,434k.
  - A CDY2A + CDY2B group of 47 reads, none untied (the pair is one sequence over what is expressed), unlabelled. Its 46
    tied reads abstain.
- **TSPY** (48 copies, a tandem array): 12 groups, 6 labelled, and no tied read meets the placement bar. **BPY2** (no
  variable column) and **PRY** (no labelled group) abstain.

**Post-hoc reading (not a bar; it does not change Y1).** The Y is haploid, so one copy gives one group of reads. In RBMY1,
two groups carry one label: the donor has two RBMY1F-like sequences that both look like 23,239k at the sites where CHM13's
two RBMY1F copies differ. Each read's call then follows its group's label, and nothing in the read tells the two CHM13
copies apart. Two explanations fit:
- an arm-to-arm conversion that left the donor's 23,004k with 23,239k's alleles, in which case about one group's reads are
  placed on the wrong copy;
- a donor duplication of RBMY1F.

CDY1 shows the converse: one group for two copies, labelled by a majority of 76 to 8. Its calls are the majority copy
for every CDY1 tied read, which is right only if the donor's other CDY1 is silent or absent. (Losing one CDY1 copy is
common: the gr/gr deletion removes one.)

Split-half agreement cannot see either problem, because both halves learn the same groups and the same labels. **The 187
agreed calls are reproducible, not provenance.** No YAG call is offered as a clear case of O2 working. The clear cases
remain those with donor truth:
- HSPA6, GTF2IRD2 and PDPR (register 1333);
- DHRSX (register 1335).

**O3 flags (groups no CHM13 copy labels):** RBMY1 has two unlabelled groups plus the duplicate RBMY1F label; DAZ has the
DAZ2 + DAZ3 group and the mixed group; CDY has the CDY2A + CDY2B group; PRY has two groups. These are candidates for donor
copies or alleles that CHM13 lacks. On these RNA reads alone, a group that is an extra copy cannot be told from one that
is a copy pair identical over what is expressed.


---

## O2_VAULT_METHODS_2026-10-20.md

*(original: `docs/archive/2026-10/merged/O2_VAULT_METHODS_2026-10-20.md`)*

## The ThesisVault's methods, assessed for O2 (2026-10-20)

**Goal (user, 2026-10-20):** "look into the @ThesisVault/ for pertinent methods that could be used for O2 or add them to the
negative results if not applicable".

**O2, as the user defines it** (`feedback_o2_thesis_core_2026-10-20`): for a long read whose primary and secondary flags are
a coin toss between near-identical copies, find its true copy from PSVs, BAM tags, statistics or other means, assign or
abstain (identical copies: the set is the answer), and use the result to improve assemblies that ignore secondaries.

**Source:** all 24 paper notes in `ThesisVault/Papers/` (their "Method in brief" sections and O2 paragraphs), checked
against this repo's negative-results register. The two methods with no register row (rescue / EM; facility location) were
tested exploratorily today with `bench/o3_graph/o2prior.py` (`personal/window/prior/eval.json`); nothing here is a
pre-registered verdict.

### Verdicts

| paper (vault note) | method, as it would apply to O2 | status for O2 | evidence |
|---|---|---|---|
| Vollger 2019 SDA; Chaisson 2017 | PSVs, read × PSV matrix, nearest copy by Hamming distance; correlation clustering of reads | **Pertinent, in use.** Reference PSV columns: tested, fail on fresh data. PSVs observed in the sample's own confidently placed reads: the corrector, the only vault-lineage method with held-out gains | r1348 (Pareto set), r1350 (substitution count) fail; corrector r1343 H1 130 fixes / 3 breaks, 40 families better / 1 worse; r1346 408 / 40; r1347 82 / 0 |
| Deschamps 2020 (review): set / clustering; uncertainty reporting | report the ambiguous group instead of a copy | Tested: the set holds the truth, but a PSV-free null does as well; the decided / coin-toss certificate equals an AS null | r1348, r1350 |
| Deschamps 2020: rescue, EM, MMR (uniquely-mapped-read ratio) | give a tied read to the copy with the most uniquely placed reads | **Not applicable (exploratory today)** | see "Statistics" below; row 1351 |
| Canzar 2016: facility location | reads may only go to "open" copies (parsimony) | **Not applicable (exploratory today)** | see below; row 1352 |
| Dishuck 2025 (NPIP): best hit beats next by >= 1 mismatch on 169 haplotypes; Guitart 2024 (TBC1D3): primary beats all by >= 10 AS, else drop | a margin rule on mismatches or AS | Tested in our forms: margin on the reference = substitution / NM count (fails); population haplotypes = population PSVs (nothing portable); dropping the ambiguous = abstention (cleans one copy, loses the other's reads) | r1350, r1273, r1349; row 1353 |
| Huang 2026 longcallR; Hosseini 2025 pHapCompass; Shaw 2025 devider | phasing / haplotype reconstruction from read co-occurrence | Not applicable: phasing is outside O2 (user, 10-20); their per-read scoring against known copies = the Pareto / likelihood set already tested; devider refuted as an O3 arm | r1348; devider arm REFUTED on all six clauses (branch `devider-arm`, 6994b228, `docs/DEVIDER_ARM_*_2026-10-0x.md`); row 1354 |
| Zheng 2025 Clair3-RNA; longcallR's editing class | separate variants, RNA editing and artefacts | Not applicable: the allele-blind rules fail on genomic donor alleles, not editing (18 / 26 gorilla breaks confirmed on his haplotype; S2 break bases shared by >= 30% of reads) | r1350 audit; row 1355 |
| Sahlin 2018 IsoCon; Tomaszkiewicz 2023 | reference-free clustering + consensus, reads to the closest candidate | Not applicable to O2 (no reference copies; cannot separate identical copies); used for O3 | r1195-r1229; row 1356 |
| Pal 2026 | per-read variant fingerprint against one representative copy | Not applicable as is (never resolved to a copy); its idea, read the sample's variants, is the corrector | row 1357 |
| Kovaka 2019 StringTie2; Tang 2020 FLAIR; Pardo-Palacios 2024 LRGASP | assemblers and their benchmark | Not applicable: baselines. FLAIR's MAPQ >= 1 realignment to its isoforms re-introduces the coin toss | r1344, r1349; row 1358 |
| Greshnova 2024 (own); Makova 2024; Rhie 2023; Yoo 2025; Bailey 2002; Bosch 2007; Cechova 2019; Soto 2025 | no read-to-copy assignment (catalogs, assemblies, SD detection, family definition) | Not applicable to O2 (O1 / O3 papers) | row 1359 |

### Statistics for the reads the bases cannot decide (exploratory, `o2prior.py`)

Reads O2's set leaves with more than one copy (window 30), truth determined. **Rescue**: the kept copy with the most
uniquely placed reads (primary, MAPQ >= 20, not tied; a fragment's own source left out). **Parsimony**: assign only when
exactly one kept copy has any uniquely placed read.

| substrate | reads | minimap2 right | rescue right / decided (minimap2 on those) | parsimony right / decided (minimap2 on those) |
|---|---|---|---|---|
| S1 one genome, round 27 (unique reads) | 1,838 | 997 | 619 / 989 (539) | 111 / 229 (104) |
| HG002 Kinnex on CHM13 | 6,526 | 4,962 | 2,865 / 6,508 (4,959) | 264 / 375 (144) |
| A119b fragments, seed 2 | 1,265 | 966 | 808 / 1,251 (963) | 52 / 54 (37) |
| gorilla held-out | 1,038 | 1,019 | 1,027 / 1,034 (1,015) | 129 / 130 (129) |
| *equal-AS ties only:* S1 | 1,669 | 831 | 539 / 861 (413) | 103 / 217 (92) |
| HG002 | 184 | 68 | 103 / 182 (67) | 45 / 76 (20) |
| A119b fragments | 396 | 184 | 254 / 385 (184) | 26 / 27 (11) |

- On reads the set keeps ambiguous, rescue is **harmful** on real data (HG002: 2,094 fewer right): minimap2's primary is
  right on 76% of them, using evidence the substitution set ignores.
- On true equal-AS ties, where minimap2 is a coin toss, rescue is right 57-66% against minimap2's 37-48%: some signal,
  far from provenance. Its precision does not rise cleanly with the support ratio (S1 0.54 / 0.71 / 0.73 / 0.57 at ratio
  < 2 / >= 2 / >= 4 / >= 8; HG002 no trend; A119b fragments up to 0.90 on 63 reads).
- Parsimony fails on copies that are expressed but never uniquely placed (S1: 103 / 217 right).
- These methods are quantification devices (expected counts), not per-read provenance.

### What is left for O2

- The only vault-lineage method with held-out gains on real reads is the **sample-aware corrector** (SDA / Chaisson's
  idea, PSVs read from the sample's own confidently placed reads, no phasing). Whether it is inside O2's scope is the
  user's decision (asked 2026-10-20).
- For reads identical in the sample too (DAZ), the set is the answer; no vault method recovers provenance there.
- Not in the vault and not tried: none of the papers models the read's *structure* against copy-specific gene models
  (ends, junction sets) as evidence for exact ties; that would need a copy-specific annotation difference to exist.


---

## O2_WINDOW_SWEEP_2026-10-20.md

*(original: `docs/archive/2026-10/merged/O2_WINDOW_SWEEP_2026-10-20.md`)*

## O2's window past Eichler's 10 AS: a development sweep on gorilla KB3781 (2026-10-20)

**User:** "we could also go past the 10 AS difference from Eichler, if that is too close to the fully tied impossible
cases"; "Sounds good, please proceed".

**Status.** Development families of `PREREG_o2_personal_copies_2026-10-20` (140 families), descriptive and exploratory. No
verdict; the held-out half was not used. Code: `bench/o3_graph/o2window.py` (+ `test_o2window.py`). Work dir:
`/mnt/linuxdisk/tmp/o3_graph/personal/window/`.

### What was swept

- **The window.** A read is tied in window W when its alignments within W of its best `ms` lie on 2 or more copies of one
  family, none off them, minimap2's primary among them. `ms` is minimap2's intron-free score; on these reads it tracks the
  edit count, one substitution = 5 points. W = 10, 20, 30, 50: **89,777** development reads are tied in some window.
- **Stitched alignments.** An alignment laying >= 50 bp on two copies' exons (SM392's chimeras) is no placement for O2.
- **O2's pick:** the tied copy with the unique fewest edits over the common window (Amendment 1's rule).
- **O2's set:** the tied copies that no other tied copy beats at every substitution column where the two differ
  (`pareto_set`). It drops a copy only when the read contradicts it and never supports it against that rival.
- **The answer key.** Every read was aligned to his maternal and paternal assemblies (6 batches, about 5 min each).
  - The registered edit key (NM + clipped bases, strict minimum) **over-penalizes clipped read ends**. In the windows up to
    20, 344 of minimap2's 754 "errors" are decided only by minimap2's alignment clipping more bases (5-19 bp in 241 of
    them). Example: SM375.0 matches the read with 0 mismatches and clips 4 end bases; SM375.1 aligns all of it with 2
    mismatches. The edit key gives the read to SM375.1.
  - The sweep therefore also uses an **intron-free key**: the copy with the highest `ms` in his genome, at least one
    substitution (5 points) above any other copy or locus (`truth_ms`).

### Results under the intron-free key (each band = the reads first tied in that window)

| band (`ms` points from the best) | scored | minimap2 wrong | O2 pick fixes / breaks | set holds the truth | set of one copy (right) | mean pool → set |
|---|---|---|---|---|---|---|
| < 10 | 519 | 225 | 3 / 60 | 517 | 218 (216) | 2.01 → 1.59 |
| 10-19 | 1,257 | 191 | 1 / 22 | 1,255 | 980 (980) | 2.03 → 1.24 |
| 20-29 | 7,700 | 97 | 0 / 11 | 7,698 | 7,155 (7,153) | 2.01 → 1.08 |
| 30-49 | 80,151 | 18 | 1 / 55 | 80,150 | 79,520 (79,519) | 2.00 → 1.01 |

**By origin, bands < 30:**

| band | reads whose copy matches the reference: scored / minimap2 wrong / set of one (right) | reads carrying the other parent's copy: scored / minimap2 wrong / set holds truth / set of one (right) |
|---|---|---|
| < 10 | 199 / 1 / 188 (187) | 320 / 224 / 319 / 30 (29) |
| 10-19 | 1,033 / 3 / 976 (976) | 224 / 188 / 222 / 4 (4) |
| 20-29 | 7,583 / 0 / 7,146 (7,146) | 117 / 97 / 115 / 9 (7) |

minimap2's errors in these three bands are almost all **SM251**: 212, 188 and 94. Its reads come from his maternal SM251.0,
which is about 41 edits from the reference's own copy.

### Reading

- **Widening adds reads minimap2 already gets right.** Past 20-30 points minimap2 is almost never wrong (18 of 80,151 in
  the 30-49 band). Nearly every read it gets wrong is inside 30.
- **The set does what the user asked.** It shrinks the pool without dropping the true copy:
  - it holds the truth on 99.6-100% of reads in every band;
  - it narrows to one copy for about 95% of reads whose copy matches the reference, and that copy is right;
  - for reads from his divergent maternal SM251.0 it keeps both copies (sets of one: 43, 40 of them right). An allele-blind
    rule cannot name that copy, and the set does not guess.
- **O2's single pick (fewest edits) breaks more than it fixes in every band.** The set is the usable output.
- **Caveat: the reference-matching reads are close to self-scoring.** For them the intron-free key on his own haplotype is
  the reference's own score ordering. minimap2's near-perfect score there, and the set holding the truth, are largely by
  construction. The independent evidence is the other-parent reads: the set holds the truth on 656 of 661 of them, but
  narrows to one copy on only 43.

### Next

Pre-register "intron-free window 30 + the set rule" and test it on fresh reads:
- the mirror reference (his other haplotype as the reference), or HG002 with its own assembly;
- simulations with a known origin.

Score the set by:
- truth held;
- set size against the pool;
- singletons right.


---

## PREREG_o2_corrector_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_corrector_2026-10-20.md`)*

## Pre-registration: O2 as a corrector of minimap2 on held-out reads (2026-10-20)

**User (2026-10-20):** "if minimap can already find the best alignment there is no need for this thesis. The entire point is
to use the tied and nearly tied AS alignments and find where they truly belong as an addendum to minimap2 since it gave
everything a coin toss. It is fine the ones that are more than >99% similar are impossible but there must be a way with
PSVs, de, nm, md tags to find the true placement of the multimappers".

So O2 is judged as a **corrector**, not a certifier. For every read whose alignments within 10 AS of its best lie on two
or more copies of one family (Eichler), O2 gives a copy. It is scored against minimap2's primary on the same reads:
- **fixes:** O2 right where minimap2 was wrong;
- **breaks:** O2 wrong where minimap2 was right;
- **accuracy:** of both, over the same reads.

Version 2 (registers 1341-1342) only certified a subset of minimap2's calls and never overruled one.

### The corrector (`bench/o3_graph/o2correct.py donor`, frozen at the commit that adds this file)

**Domain.** A read whose alignments within 10 AS of its best lie on two or more copies of one family and on no locus off it,
minimap2's primary among them.

**The sample's copies.** Per copy, the sample's confident reads on it: the primary is on the copy and every alignment on
another copy of its family is more than 9 AS lower. From them, each copy's **events**: substitutions, and indels shifted
left as far as the reference repeats them.

**Scoring.** All tied alignments of a read are compared over the read bases every one of them aligns. A tied copy scores:
- the read's events its own reads do not explain. An event is explained when at least 2 of them carry it, making up at least
  10% of their depth, at a depth of at least 3;
- plus the events at least 80% of its reads carry that the read lacks.

**Decision.** The lowest score takes the read. It overrules minimap2's primary only if it is still lower after the primary
is granted every event its reads do not cover (dominance). A shared lowest keeps minimap2's pick. A fragment's own source
read is left out of the sample that judges it.

**In one genome** a copy's reads explain only their shared errors, so the rule is minimap2's ordering by events.

### Development: every number seen before this file

- **HG002** (`hg002proof`, donor-assembly truth, PAR out). The rule and its constants were chosen here:
  - substitutions only (the MD/NM rule);
  - then the sample-aware rule, with HET 0.2 / 0.1 / 0.05 crossed with dominance on/off.

  | HG002, chosen rule | reads | corrector | minimap2 | fixes / breaks |
  |---|---|---|---|---|
  | near ties | 5,488 | 87.3% | 66.4% | 1,161 / 15 |
  | exact ties | 576 | 52.6% | 39.6% | 97 / 22 |

  15 families gain and 4 lose. The substitution-only rule was right in two families (SMN, THOC3) and a coin flip elsewhere.
- **Simulation, heldout23** (one genome). Breaks in CTRB2 led to comparing alignments over common read bases. After that,
  near ties were 2,016 reads with 0 fixes and 0 breaks.
- **Simulation, heldout21** (donor alleles, PAR out):

  | heldout21 | exact ties fixes / breaks | near ties fixes / breaks |
  |---|---|---|
  | worlds without BRANCH | 31 / 0 | 0 / 0 |
  | BRANCH worlds | 38 / 35 | 0 / 7 |

  A BRANCH world's reference lacks one of the donor's copies (O3's case). Its reads, confidently placed on the nearest
  reference copy, teach that copy alleles it does not carry.
- **A119b fragments, seed 1** (`o2frag_a119b`, qualifying non-PAR families). This set also exercised the library pass.
  1,518 fragments: corrector 80.0% vs minimap2 72.1%, 125 fixes / 5 breaks; 45 families gain, 2 lose.

### Held-out sets: none has been run through the corrector

| set | what | truth |
|---|---|---|
| **H1** | A119b fragments, **seed 2** (`o2frag_a119b2`: 11,280 sources, 33,840 fragments). Used only by version 2 in register 1342, never by this rule. Held out by reads, not by donor. | the source read's copy, clear by 30 AS |
| **H2** | public testis fragments (`o2frag_public`, library ERR13885926): **another donor**, shallow | the same |
| **S1** | simulation, one genome: heldout24 + heldout25 (36 worlds: CROSS, INTACT, SILENT) | `reads.truth.E0.tsv` |
| **S2** | simulation, donor alleles: heldout22 (36 worlds: BRANCH, CROSS, INTACT, OUTSIDE, REMOVED, SILENT) | the same |

**Common scoring rules:**
- A read is scored only when its truth is one of its tied copies (O2 cannot pick an absent copy).
- PAR families are out (`famident/par_families.json`).
- Fragments are counted one by one; their source reads are counted for power.
- Real sets are scored on the qualifying non-PAR families: every member pair at [0.90, 1.00) identity,
  `famident/qualify_nonpar.json`.
- Commands: `o2correct.py score OUT truth.tsv --among-candidates --families qualify_nonpar.json` and
  `o2correct.py score-worlds OUT.json WORLDS... --exclude par_families.json`.

### Bars

**H1 and H2, each, pooled over exact and near ties:**
- **P (power):** at least 50 overruled fragments, from at least 20 source reads.
- **G (gain):** corrector accuracy at least 2 points above minimap2's.
- **S (safety):** fixes at least 5 times breaks.
- **F (families):** at least 11 qualifying families where it fixes more than it breaks (user: "more than 10 families"). No
  more than a fifth as many families, rounded up, where it breaks more than it fixes.

**S1, without BRANCH worlds (none exist in S1):**
- near ties: corrector accuracy no more than 0.2 points below minimap2's, and breaks no more than fixes + 2;
- exact ties: corrector accuracy no more than 0.5 points below minimap2's.

**S2, the worlds without BRANCH:**
- corrector accuracy at least minimap2's in each band;
- breaks at most fixes / 5 + 2.

BRANCH worlds are reported apart. That is O3's case, and development already showed its harm.

### Verdict

- **PASS:** H1 passes P, G, S and F, and S1 and S2 pass.
- **FAIL:** H1 has P and fails G, S or F; or S1 or S2 fails.
- **UNDERPOWERED:** H1 fails P.
- **H2** gets its own PASS / FAIL / UNDERPOWERED on the same bars. A FAIL there limits the claim to A119b's donor.

**A PASS would claim:** among a read's tied and near-tied copies, the sample's own reads of each copy overrule minimap2's coin
toss far more often right than wrong. This holds on reads the rule never saw, in more than 10 near-identical families, and
does no harm in one genome.

**It would not claim:**
- accuracy against a donor assembly on a new donor (H1 and H2 have fragment truth);
- copies absent from the reference (O3);
- PAR;
- exact ties over identical sequence.

**Known bias:** a fragment's truth is minimap2's call on its full read. A systematic full-read error therefore counts against
the corrector, never for it.

**Rules:**
- Run script: `/mnt/linuxdisk/tmp/o3_graph/o2correct_heldout/run.sh`. The library pass resumes in foreground batches (`--cache`,
  `--budget`).
- A bug found during the run is fixed and recorded as an amendment before any rerun.
- No constant changes after this file.

### Results (`/mnt/linuxdisk/tmp/o3_graph/o2correct_heldout/run1.log`, `run2.log`; frozen copy `frozen/`)

**Overall: FAIL. S2 fails. H1 and S1 pass. H2 is UNDERPOWERED.**

| set | reads scored | corrector | minimap2 | fixes / breaks | bars |
|---|---|---|---|---|---|
| **H1** A119b seed 2, qualifying non-PAR | 1,241 fragments | **0.8171** | 0.7147 | **130 / 3** | P ✓ (133 overruled from 118 sources), G ✓ (+10.2 points), S ✓, F ✓ (40 families better, 1 worse; the limit was 8) |
| H1 exact ties | 440 | 0.6091 | 0.4523 | 69 / 0 | — |
| H1 near ties | 801 | 0.9313 | 0.8589 | 61 / 3 | — |
| H1, every non-PAR family | 1,727 | 0.7823 | 0.6775 | 186 / 5 | 55 families better, 2 worse |
| **H2** public testis, qualifying | 276 | 0.7826 | 0.7464 | 11 / 1 | **P fails** (12 overruled, from 12 sources): UNDERPOWERED |
| **S1** one genome (36 worlds), near ties | 2,309 | 0.9991 | 0.9991 | 0 / 0 | ✓ |
| S1 exact ties | 2,824 | 0.4770 | 0.4763 | 2 / 0 | ✓ |
| **S2** donor alleles, without BRANCH (30 worlds), near ties | 3,447 | 0.9985 | **1.0000** | **0 / 5** | **✗** (accuracy below minimap2's; 5 breaks > 14 / 5 + 2) |
| S2 exact ties | 3,079 | 0.4303 | 0.4258 | 14 / 0 | ✓ |
| S2 BRANCH (reported apart) | 786 exact, 812 near | — | — | 4 / 1 exact, 0 / 2 near | — |

**The S2 failure is the BRANCH mechanism, which development showed.** All 5 breaks are in one REMOVED world,
`A119b_s1201`, family ANKRD30BL.
- The world masks copy .3 out of the reference while the donor still expresses it: 111 reads.
- Those reads sit confidently on copy .4 and make up most of its sample.
- True .4 reads then look foreign on .4, and the corrector moves them to .0 (scores .0:4 vs .4:5).
- 10 more overrules move .3 reads themselves; they are not scored.

BRANCH breaks the same way in 3 of its 6 worlds. Copies the reference lacks are O3's case; O2's corrector is not safe on a
family until O3 has flagged them.

**Erratum in the development section.** "Worlds without BRANCH" in heldout21 was 24 worlds, not 30. All 6 REMOVED runs
crashed: their `genome.ref.fa` points to a deleted masked genome. `score-worlds` dropped them without saying so; it now lists
skipped worlds. In heldout22 all 36 worlds ran.

### Audit (metric-auditor, 2026-10-20): holds with caveats

Every number recomputes exactly from the frozen files with an independent script, and the overall FAIL stands. The audit also
confirmed:
- **Hold-out:** 0 of 15,573 seed-1 sources are reused. The prereg was committed at 12:50:27, before the first held-out output
  at 12:50:36.
- **Leave-one-out:** 1,329 of the 1,349 domain sources were confident reads, and each was left out of its own fragments'
  judgement.

Caveats:
1. **H1's truth and the corrector's teachers come from the same library's confident calls.** When the donor has a copy the
   reference lacks, both are wrong the same way, and that is scored as a fix. So this file's "never for it" is false in that
   case. H1 "right" means "agrees with minimap2's full-read call". The top fix families (NBPF6, TBC1D3B, FAM205A) vary in
   copy number between people. Only donor-assembly truth settles this; HG002, the one donor with it, was development.
2. **H1 replicates development rather than testing a new donor.** It shares seed 1's library, its 161,037 teacher reads and
   its families; only the source reads are disjoint. H2, the other donor, is underpowered.
3. **Fixes are spread out.** 130 fixes come from 115 sources. The top family has 10.8% of them, the top 10 have 59.2%, and
   17 families have at least 3.
4. **Filters do not decide a bar.** 14 non-PAR overrules have no scored truth. S still holds if all 14 were breaks
   (130 vs 17). G is still +10.0 points if the 33 rows whose truth is not among the tied copies count against it.
5. Excluding REMOVED from S2 after the fact (near 0/0, exact 11/0) would be choosing which worlds to score. It is not a
   result.

**What this shows.** Reading PSVs against the sample's own copies, learned from its confidently placed reads, corrects
minimap2 among tied copies:
- HG002 development, donor truth: near ties 66% → 87%;
- fresh A119b fragments: 72% → 82%, in 40 families;
- one genome: nothing to correct and nothing broken.

It is **not adopted**. The pre-registered donor simulation shows the corrector moves reads wrongly whenever the reference lacks
one of the donor's copies (BRANCH, REMOVED), a case O3 must flag first. A next version needs that guard and fresh held-out
data, ideally a second donor with an assembly.


---

## PREREG_o2_corrector_families_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_corrector_families_2026-10-20.md`)*

## Pre-registration: does the sample-aware corrector recover a share of minimap2's errors across many families? (2026-10-20)

**User (2026-10-20):** after the vault survey found the sample-aware corrector to be the only O2 method with held-out gains:
"ok can we try to prove that it works in other families, maybe not DAZ but as long as it can recover a small percent.."

### The rule (frozen; nothing changes)

`bench/o3_graph/o2correct.py donor`, version 2, unchanged since bfa15f12, default flags (no `--v1`, no `--guard`):
- **Domain:** a read whose alignments within 10 AS of its best (Eichler) lie on two or more copies of one family and on no
  locus off it, minimap2's primary among them.
- **The sample's copies:** per copy, its confident reads (primary on it, every alignment on another copy of its family more
  than 9 AS lower) give its events (substitutions, left-shifted indels): PSVs as this sample carries them. No phasing, no
  haplotypes, no assembly.
- **Decision:** the tied copy whose alignment has the fewest events its sample reads do not explain takes the read, and
  overrules minimap2 only when it dominates; otherwise minimap2's pick stays.

### Substrates, none ever seen by this rule

- **D (primary): a never-used donor, simulated from its own assembly.** HG03516 (HPRC year 1, Esan, Nigeria; hifiasm
  phased `year1_f1_assembly_v2_genbank`, maternal and paternal), never mentioned in this project. The recipe is register
  1346's, unchanged:
  - every CHM13 copy of the 281 qualifying non-PAR families (`famident/qualify_nonpar.json`), ± 20 kb, aligned to each
    haplotype (`minimap2 -x asm20 -c --secondary=yes -N 10`);
  - `donorsim.py orthologs`, then `donorsim.py build ... A119b_testis_v4_ref.json --per-copy 150 --seed 3516`: reads
    simulated from the donor's own transcripts, so its alleles, indels and missing copies are real;
  - aligned to CHM13 v2.0 (`minimap2 -a -x splice:hq -uf --secondary=yes -N 100 -p 0.5`).
  - Truth: the CHM13 copy whose donor ortholog produced the read.
- **A (replication, real reads): a third A119b fragment cut.** `o2frag.py sources` with `--seed 3`, excluding every
  source read of seeds 1 and 2; fragments `--per-read 3 --min 300 --max 1200 --seed 3`; truth = the source read's copy,
  kept when clear by 30 AS (register 1343's design and known bias). The corrector learns from the whole library
  (`--teach A119b.t2t.bam`), each fragment's source read left out.
- **N (no harm, one genome):** simulation round 27's six CROSS worlds (`fast/heldout27/*/CROSS`; INTACT and SILENT repeat
  their reads), never run through the corrector.

### Scoring

`o2correct.py score OUT.tsv TRUTH.tsv --among-candidates --families qualify_nonpar.json` (D also `--donor`): a read is
scored when its truth is among its tied copies. Per substrate: scored reads, minimap2 right, corrector right, fixes
(corrector right where minimap2 was wrong), breaks (the reverse), and per family whether it fixes more than it breaks.
N: `o2correct.py score-worlds` on the six CROSS worlds, PAR families excluded.

### Bars (D and A each)

- **P (power):** at least 300 scored reads from at least 11 families, minimap2 wrong on at least 50.
- **R (recovers a share):** fixes − breaks >= 2% of minimap2's wrong calls, and >= 20 reads.
- **S (safety):** fixes >= 5 × breaks.
- **F (families):** at least 11 families where it fixes more than it breaks, and no more than a fifth as many (rounded up)
  where it breaks more.

**N (no harm):** on the one-genome worlds, breaks <= fixes + 2.

**Verdict:** **D PASS** = P, R, S, F on D, and N. FAIL = P holds and a bar fails (or N fails). UNDERPOWERED = P fails.
A gets its own verdict on P, R, S, F; a claim across a real library needs A to PASS as well.

**What a PASS would claim:** for reads tied within 10 AS between near-identical copies, reading each copy's PSVs from the
sample's own unambiguous reads recovers a share of minimap2's wrong placements, many times more than it breaks, in more
than 10 families of a donor it never saw. **It would not claim:** reads identical in the donor too (they stay minimap2's
coin toss: no fixes, no breaks), copies the reference lacks (O3), the assembly step, or anything past 10 AS.

### Rules

- Each substrate is built and scored once; bugs → amendment + rerun.
- Heavy jobs (the donor download and alignments, the CHM13 splice alignments, the library pass) in the foreground, one at
  a time, each command under 10 minutes; the corrector's library pass resumes from its cache.
- The metric-auditor checks the verdicts on a frozen copy before they are reported. Human numbers only.
- Push only when asked.

### Amendments

#### Amendment 1 (2026-10-20, before any data): gorilla KB3781 only — replaces the substrates, bars D / A / N

**User:** "wait no, I think this is one of the over-reaching problems my advisor complains, can we stick with gorilla
fibroblasts and its T2T genome?"

- **Withdrawn, unrun:** D (HG03516), A (A119b seed 3), N (round 27). The interrupted command had started HG03516's
  maternal download; the 84 MB partial file was deleted unread. Nothing else was built.
- **The substrate (G2):** KB3781's fibroblast Iso-Seq on the GGO primary (`GCA_029281585.2_flnc_mm.bam`), in **the 383
  GWFAM families never scored** (`rna_allele/fibro_cat.copies.tsv` less the 283 of `personal/qual.json`, used by registers
  1347-1350). Every alignment of every read touching their copies is cut from the library, as `personal/sub.bam` was.
  DAZ-like reads (identical on two copies in his own genome too) are undetermined and not scored.
- **The rule:** the frozen corrector, as above, with that family table.
- **The answer key (evaluation only; the corrector never sees it):** each domain read aligned to his maternal and paternal
  assemblies; its copy is the one whose alignment has the best intron-free score (`ms`), at least one substitution (5
  points) ahead of any other copy or locus (`o2window.py truth-ms`), through each copy's ALT ortholog (`o2personal.py
  queries` / `orthologs`, unchanged). **refcons** = his reference-haplotype copy is among the best (his copy equals the
  reference); **altonly** = only the other haplotype's.
- **Scoring:** the frozen `o2correct.py score --among-candidates` on the determined reads (truth as `read, copy`); then the
  same on the refcons reads alone.
- **Bars (G2):** P (>= 300 scored reads from >= 11 families, minimap2 wrong on >= 50), R (fixes − breaks >= 2% of
  minimap2's wrong calls and >= 20), S (fixes >= 5 × breaks), F (>= 11 families fixing more than breaking, <= a fifth as
  many rounded up the reverse), **H (no harm where his copy equals the reference: on refcons reads, breaks <= fixes + 2)**.
  PASS = all five; FAIL = P holds and another fails; UNDERPOWERED = P fails.
- **Known caveats, stated now:** on refcons reads the key is the reference's own ordering (near self-scoring; H is a harm
  check, not a gain); in the 283 families tested before, the corrector's gorilla gains were concentrated (r1347: 82 fixes,
  76 in SM392).

#### Amendment 2 (2026-10-20, before any corrector run on these families): the thesis's mainstay families only

**User:** "I don't need the GWFAM, my advisor is very picky and he does not believe any of the method works, we need to stick
to the YAGs or mainstays that we have tried earlier, unless those families don't have nearly AS tied multi-mapping reads";
on the census below and the proposed test: "Yes".

**Census (substrate only; no rule run):** KB3781's whole fibroblast library, reads with a primary on a copy / reads whose
alignments on two or more of the family's copies lie within 10 AS:
- YAGs on his Y (RBMY 13, TSPY 6, HSFY 5, BPY2 2, DAZ 2, CDY 1 copies): **0 reads** (testis-only genes): untestable here.
- TBC1D3 (14 copies) 5 / 3; SMN (RefSeq, 3) 1,478 / 0; GOLGA6 (2) 0 / 0; FAM90A (1) 4 / 0.
- **NPIP** (the curated 25-copy set, `fibro_cat` family `NPIP`): 1,040 / **182** (44 exact).
- **LRPAP1** (its 22-copy family, `fibro_cat` `SM22`; the 10-04 case-study copies): 2,555 / **38**.
- SM21 (204 tied) is left out: genes riding in the NPIP duplication (PLA2G10, PARN, SRRM2-like), not a mainstay.

**Replaces Amendment 1's family set and bars.** Families: NPIP and LRPAP1 (47 copies; neither among the 283 scored in
registers 1347-1350). Everything else in Amendment 1 holds (his library, the frozen corrector, his two assemblies as the
answer key by the intron-free score, refcons / altonly).

**Bars (pooled over the two families; each family reported):**
- **P:** >= 50 scored reads and minimap2 wrong on >= 10; else UNDERPOWERED.
- **R:** fixes > breaks, and fixes − breaks >= 2% of minimap2's wrong calls.
- **S:** fixes >= 5 × breaks.
- **H:** on refcons reads, breaks <= fixes + 2.

PASS = P, R, S, H; FAIL = P holds and another fails.

### Results (code 67383d28; run once; `o3_graph/mainstay/`): **UNDERPOWERED**

KB3781's fibroblast reads touching NPIP's 25 and LRPAP1's 22 copies: 5,821 reads, 77,090 alignments. The corrector's
domain (tied within 10 AS): **185 reads** (NPIP 164, LRPAP1 21); it learns from 3,037 confident reads on 33 copies and
overrules minimap2 on 4.

**What his genome says about the 185** (each copy's counterpart on his other haplotype lifted for 20 of the 47 copies):

| family | tied | determined, his copy = the reference | determined, other parent's copy | identical on two copies in his genome | no traceable copy |
|---|---|---|---|---|---|
| NPIP | 164 | 75 (1 outside the tie) | 6 | 48 | 35 |
| LRPAP1 | 21 | 0 | 1 | 2 | 18 |

**Scored (truth among the tied copies): 81 reads; minimap2 right 79, corrector right 79; fixes 0, breaks 0.**
- P fails (minimap2 wrong on 2, bar 10): **UNDERPOWERED.** R and S cannot be met with 0 fixes; H holds (refcons 74: 0
  breaks).
- minimap2's 2 errors are NPIP reads (one his reference copy NPIP.14 placed on NPIP.16; one his other parent's NPIP.18
  placed on NPIP.17); the corrector kept minimap2's pick on both.
- The corrector's 4 overrules are all LRPAP1 reads whose best place in his genome has no traceable reference copy, so
  the key cannot score them. Descriptive look (not a registered result): two are reads whose best match is his maternal
  copy of SM22.1, minimap2's pick, by one substitution (the overrule would be a break); two leave a copy his genome
  scores 5 and 30 points worse for one of two equally best copies.

**Reading.** In the matched sample, for the thesis's mainstay families, minimap2 is already right on 79 of the 81 tied
reads his own genome can place. The rest are reads identical on two copies in his genome too (50: no method can place
them) or in duplications his assemblies cannot trace (53). There is almost nothing left for O2 to recover here. YAGs are
not expressed in fibroblasts and cannot be tested in this sample. Audit pending.

### Audit (metric-auditor, 2026-10-21): UNDERPOWERED holds with caveats; the reading above is withdrawn

Every number reproduces (0 status mismatches over 185 reads). UNDERPOWERED survives every orthology-respecting repair of
the lift once a key artefact is removed (minimap2 errors 4-9 against the bar of 10). Corrections:

- **Self-scored majority:** 74 of the 81 scored reads are refcons, keyed by the reference's own best alignment (73 / 74).
  The only informative reads are the **7 from his other parent's copy: minimap2 right on 6.**
- **Key artefact:** one of minimap2's two "errors" (SRR27178663.1441099) and the "outside the tie" read are truths built on
  stitched alignments (272-418 kb gaps; `truth_ms` never applied `stitched()`). **Genuine scored errors: 1**
  (SRR27438212.7310158, NPIP.17 vs NPIP.18).
- **The 53 off-copy reads are outcome-linked:** all have minimap2's copy >= 5 ms below their best placement (median 262).
  - **29 NPIP reads come from one maternal-only locus, mat chr18 (CM054600.2) 6.72-6.82 Mb**: 0-3 substitutions there,
    18-167 at their best paternal placement. A copy the reference lacks (O3's case).
  - 24 are traceable but the lift refused their copies: 11 NPIP copies tripped the 50-bp exon guard (MAPQ 31-60,
    99.77-99.97% identity), 11 SM22 dropped as ambiguous, 6 SM22 refused as interchromosomal. Repaired lifts give 105-113
    scored reads and 4-9 genuine minimap2 errors.
- **The corrector could not act on 66 of 185 reads** (minimap2's copy has no confident reads, so dominance never passes),
  including the one genuine error. NPIP.3's sample carries 1.35% events, consistent with a missing copy's reads teaching
  their nearest copy foreign alleles (unverified).
- **The 4 overrules (LRPAP1):** two are exact reference ties where the maternal locus is 99.83% / 99.94% identical to the
  two candidates (undecidable); in the other two his genome scores the corrector's copy 4-5 ms below minimap2's. The
  descriptive "two likely breaks" and "one of two equally best" were wrong.
- "50 identical" is 29 exact ties, 5 within 1-4 ms, 16 tied with an unlifted locus. "YAGs 0 reads" is 0 primary
  alignments (107 secondary, from 51 reads at 2.7-16.7% divergence). The census counts were not saved and reproduce
  only approximately (NPIP 175-180 tied; LRPAP1 22-61).
- Wording not to use: "minimap2 is right on 79 of 81"; "almost nothing left to recover"; "53 his assemblies cannot trace";
  "no method can place the 50"; "H holds" as safety; "minimap2 wrong on 2" without the artefact.

**What stands:** too few minimap2 errors among NPIP / LRPAP1 tied reads that his genome places on a reference copy to test
the corrector (1 genuine; 4-9 under repaired lifts). **What the data does show:** the largest group of minimap2's
mis-placements here is 29 NPIP reads from a maternal NPIP copy absent from the reference — O3, not a choice among
reference copies.

#### Amendment 3 (2026-10-21, before rescoring): the repaired answer key

**User:** "lets try 1" (fix the answer key and rescore the frozen corrector output; no new alignments).

**Already seen, stated plainly:** the r1360 audit's exploratory bracket. Under several lift repairs the scored set was
105-113 reads with 4-9 genuine minimap2 errors (UNDERPOWERED in each). This rescore applies one convention, chosen by
principle and fixed here; it is a confirmation, not a blind test. **Expected: UNDERPOWERED.**

**The repaired key** (`--repaired` on `o2personal.py orthologs` and `o2window.py truth-ms`; without it both reproduce
r1360's files byte for byte):
1. **Key alignments:** an alignment of a read to his assemblies is no placement when it is stitched across two copies'
   exons (`stitched()`, >= 50 bp on each) or bridges a gap of >= 100 kb (`no_placement`, the audit's recipe).
2. **Counterparts (identity convention):** a copy's counterpart on his other haplotype is the primary alignment (MAPQ
   >= 20) of its ± 20 kb query **on any contig that is not the reference's own** — the homologous chromosome or another;
   its exons are lifted without the 50-bp growth guard.
3. **Overlaps:** when two copies' counterparts overlap, the copy whose query aligned with the lower divergence (`de`)
   keeps it; the other gets none (instead of both dropped).

Unchanged: the frozen corrector output (`arm.tsv`), the scorer, the bars (P: minimap2 wrong on >= 10 of >= 50 scored; R,
S, H) and the verdict rules.

### Results under Amendment 3 (rescore, run once; `o3_graph/mainstay/repaired/`): **UNDERPOWERED**

- **Repaired counterparts:** 36 of 47 copies (11 lose an overlap to a copy with lower divergence). **Repaired key:** 51
  stitched or >= 100-kb-gap alignments dropped; statuses of the 185 tied reads: refcons 84, altonly 25, undetermined 42,
  off-copy 34 (NPIP 94 determined / 40 / 30; LRPAP1 15 / 2 / 4).
- **Scored 107:** minimap2 right 104, corrector right 106; **fixes 2, breaks 0**; families better 1 (LRPAP1) / worse 0.
  P fails (minimap2 wrong on 3, bar 10): **UNDERPOWERED**. R, S and H would hold.
- **Where his copy equals the reference (82 reads):** both 82 / 82 — scored against the reference's own ordering.
- **Reads from his other parent's copy (25, the informative ones):** minimap2 right 22, corrector 24; 2 fixes, 0 breaks.
  - The 2 fixes are LRPAP1 reads, SM22.1 -> SM22.18: exact AS ties on the reference, decided in his genome by one
    substitution. They are fixes under the identity convention fixed in Amendment 3; under synteny they would be breaks.
  - The remaining error, NPIP.17 for NPIP.18, is one the corrector cannot act on (no confident reads on either copy).

Audit pending.

#### Amendment 4 (2026-10-21, before the run it names): counterparts by position — final key for this experiment

**Why.** The audit of the Amendment 3 rescore found the identity convention chosen after the evidence was seen, and it
sets the direction: identity 2 fixes / 0 breaks (UNDERPOWERED); positional variants 0 / 2, one of them with 11 minimap2
errors (FAIL). The two LRPAP1 reads hinge on whether a maternal chr12 locus (25.85-25.91 Mb) is SM22.1's counterpart (same
chromosome) or SM22.18's (closer sequence, chr16). The project's own precedent judges counterparts by position (registers
1181, 1207).

**The one change:** each copy's counterpart is found with **± 500 kb** of reference context (`o2personal.py queries /
orthologs --flank 500000`), past the duplication, so unique sequence fixes its position; a counterpart on another
chromosome then means a real rearrangement, not a closer paralog. Everything else in Amendment 3 stays (longest primary
with MAPQ >= 20 on a non-reference contig, no growth guard, overlaps to the lower divergence, stitched / >= 100-kb key
alignments dropped), as do the frozen corrector output, the scorer and the bars.

**Seen before this amendment:** the audit's variants listed above. **The verdict this run gives is final for this
experiment; the key is not changed again.**

### Results under Amendment 4 (run once; `o3_graph/mainstay/positional/`): **UNDERPOWERED** — final for this experiment

- **Counterparts by position (± 500 kb):** 40 of 47 copies; 7 LRPAP1 copies lose an overlap to a copy with lower
  divergence (SM22.2, .5, .7, .12-.15). **Key:** 56 stitched or >= 100-kb-gap alignments dropped; statuses of the 185 tied
  reads: refcons 85, altonly 24, undetermined 44, off-copy 32 (NPIP 96 determined / 39 / 29; LRPAP1 13 / 5 / 3).
- **Scored 107** (NPIP 96, LRPAP1 11): minimap2 right 104, corrector right 106; **fixes 2, breaks 0**; families better 1
  (LRPAP1) / worse 0. P fails (minimap2 wrong on 3, bar 10): **UNDERPOWERED**. R, S and H would hold.
- **Where his copy equals the reference (83 reads):** both 83 / 83 (scored against the reference's own ordering).
- **Reads from his other parent's copy (24, the informative ones):** minimap2 right 21, corrector 23; 2 fixes, 0 breaks.
- **The 2 fixes** (LRPAP1 SRR27438212.2337814 and .5312667): exact AS ties on the reference between SM22.1 and SM22.18
  (ms 1541 on both); minimap2 SM22.1, the corrector SM22.18 (one event fewer that the sample's reads do not explain).
  - In his genome each read fits best at maternal chr12 25.87 Mb, by one substitution (5 ms points, the key's minimum
    margin); single-copy alignments, not stitched.
  - With ± 500 kb of context, SM22.1's maternal counterpart is at chr12 27.43 Mb (MAPQ 60), where the reads fit 82-85 ms
    points worse: **minimap2's pick is wrong under position.**
  - The maternal 25.87-Mb locus is claimed by two reference blocks: the chr16 block of SM22.18-.20 (946 kb aligned,
    divergence 0.0020) and the chr12 block of SM22.12-.15 (670 kb, 0.0057). The overlap rule (Amendment 3, kept) gives it
    to the lower divergence, SM22.18. Had it gone to SM22.14, the truth would lie outside the reads' tied copies:
    unscored, no fix and no break. **No position reading makes these breaks.**
  - The reference's (paternal) chr16 block sits on his maternal chr12: a rearrangement between his parents' haplotypes or
    a mis-join in one assembly. The key cannot tell which.
- **The remaining error** (NPIP SRR27438212.7310158): truth NPIP.18 by 68 ms points (his maternal copy); minimap2 and the
  corrector both NPIP.17. With no confident reads on either copy, the corrector counts the read's own differences
  (22 vs 23), which lean the wrong way by one.
- **Key artefact, no effect on the verdict:** 2 LRPAP1 reads (SRR27438212.2142599, .3445690) span two adjacent copies
  (SM22.1 -> SM22.11, SM22.20 -> SM22.18, and the maternal pair). The key drops those alignments as stitched and takes a
  worse one (176 and 100 ms points below their best) on SM22.14 as the truth. SM22.14 is outside their tied copies, so they are
  unscored; a read whose best alignments are all stitched should be undetermined.

Audit: below (holds with caveats).

### Audit of Amendment 4 (metric-auditor, 2026-10-21): UNDERPOWERED holds with caveats; the direction is not established

Every number reproduces: rebuilt counterparts equal `orthologs.tsv` (40 / 47), the rebuilt key equals `truth_ms.tsv` on
all 185 reads, the headline matches to the digit. UNDERPOWERED survives every ± 500-kb reading (minimap2 wrong on 1-3
against the bar of 10). Corrections:

- **The direction rests on divergence at one locus, not on position.**
  - Maternal chr12 25.82-25.89 Mb ("X") decides 9 of the 11 scored LRPAP1 reads: the 2 fixes, and 7 where both arms pick
    SM22.18.
  - ± 500 kb does not reach unique sequence there. The chr12 windows cover X too (SM22.14's: 682 kb at divergence
    0.0057; SM22.1's inverted upstream flank: 0.0078), so the lower-divergence rule picks the owner.
  - SM22.18's 0.0020 is closer than this region's own allelic windows (0.0054-0.0077 on chr12 and chr14).
  - X to the same-chromosome block: 98 scored, 1 minimap2 error, 0 fixes / 0 breaks. Overlap dropped: 97, 1, 0 / 0.
  - **Never 0 / 2 at ± 500 kb.** SM22.1's two flanks anchor at mat chr12 27.39-27.93 Mb (MAPQ 60), and the reads fit
    there 82-85 ms points worse. .5312667 matches X exactly (NM 0).
- **The 2 fixes are 2 reads from one locus, each by one substitution:** an anecdote, not a direction.
- **The key artefact reaches 2 scored reads.** The note above names only 2 unscored ones. These are SRR27178663.1441099
  (truth NPIP.16, 43 ms below its dropped best: the r1360 artefact read, now counted "right") and SRR27438212.1546517
  (NPIP.5, 27 below). Both count as right for both arms. `best_ms` is computed after the drop (`o2window.py` truth-ms),
  which hides them. With such reads undetermined: 105 scored, minimap2 102, corrector 104; verdict unchanged. None of the
  3 minimap2 errors rests on a leftover.
- **Self-scoring and exclusion:**
  - 83 of the 107 scored reads are refcons (both 83 / 83). Only the 24 altonly reads inform (minimap2 21, corrector 23).
  - 76 of the 185 tied reads are unscored. All 32 off-copy reads have minimap2's copy >= 5 ms below their best placement;
    the 29 maternal-only NPIP reads are 104-987 below.
- **Forking paths:** the key was amended twice after the data were seen. Amendment 4 was written knowing that X sets the
  direction. It was committed before its outputs; the only code change is the flank; the corrector output, BAMs, copies
  and chrmap are byte-identical.
- **Wording:** "ms 1541 on both" holds for .2337814 only (.5312667: 1578 on both). Drop "R, S and H would hold": P failed,
  and those bars would rest on 2 reads.

**What stands (final for this experiment):** UNDERPOWERED. Among the NPIP / LRPAP1 reads tied on the GGO reference that
his own genome places on a reference copy, minimap2 is wrong on 1-3: too few to test the corrector. The corrector broke
nothing. Its 2 fixes rest on one locus whose owner is set by divergence, not position. Settling that locus needs HiFi
reads spanning the maternal chr12 25.4-26.4 Mb junctions (not run). **The key is not changed again.**


---

## PREREG_o2_fair_a119b_fresh_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_fair_a119b_fresh_2026-10-20.md`)*

## Pre-registration: O2 against minimap2's trusted calls on fresh A119b reads (2026-10-20)

**User (2026-10-20):** "no lets focus on redoing it on CHM13 with A119b". That replaces the cerebellum donor of
`PREREG_o2_fair_heldout2_2026-10-19.md`: its whole-genome alignment was stopped by the session for low memory.

**What is held out:** the reads, not the donor. A119b was O2's development library:
- version 2's validated-only rule was chosen after register 1339's errors on A119b fragments;
- k = 3 was chosen on the same fragments;
- the fair scorer's development numbers came from them (232 placed, 1 wrong, against the baseline's 103 wrong of 773).

This draw takes **none of those 15,573 source reads**. Every fragment, and every validating read's role, is new.

### Design (register 1341's, unchanged except for the draw)

- **Library:** A119b testis, `winloci_data/A119b.t2t.bam` (CHM13 v2.0; minimap2 splice:hq -uf --eqx -Y -N 50 -p 0.1).
- **Sources:** up to 60 per pool family, MAPQ-60 primaries of >= 1 kb, **seed 2**, excluding `o2frag_a119b/sources.tsv`
  (`o2frag.py sources --exclude`, 4a57339f).
- **Fragments:** 3 per source, 300-1,200 bp, seed 2.
- **Realignment:** with A119b's parameters plus `-K 100M` (batch size: memory only). Written unsorted, then sorted.
- **Truth:** the source's copy when its full read is clear by 30 AS; scored only when that copy is one of the fragment's
  tied copies.
- **O2:** `o2psv.py svp --k 3 --validated-only` (7bf68347). The validators are A119b's own reads, a fragment's own source
  left out.
- **Scorer:** `o2perfamily.py fair` (fad7e43b) on the 281 qualifying non-pseudoautosomal families, by source read, against
  minimap2's primary at MAPQ >= 1.

### Bars (unchanged)

- **F1:** at least 100 source reads placed.
- **F2:** at most 2% wrong.
- **F3:** precision at least 5 points above the baseline's.
- **F4:** more than 10 families with at least 5 placed and none wrong.

**Verdict:**
- **PASS:** F1 to F4 hold.
- **FAIL:** F1 holds and another bar fails.
- **UNDERPOWERED:** F1 fails.

**A PASS would claim:** on A119b reads O2 has never seen, O2 places almost-tied reads more precisely than minimap2's own
trusted calls, in more than 10 near-identical families.

**It would not claim:** another donor (registers 1341 and the cerebellum attempt), or exact ties. Exact ties carry no
separating PSV; register 1340.

**Rules:** code frozen at 4a57339f. Run script: `o2frag_a119b2/run.sh`. A bug is recorded as an amendment before any rerun.

### Results (`o2frag_a119b2/run.log`, `fair.json`)

**PASS: F1 to F4 hold.**
- **F1:** 197 source reads placed.
- **F2:** 0 wrong. With 0 of 197, the 95% upper bound on the error rate is 1.9%.
- **F3:** precision 1.000, against the MAPQ >= 1 baseline's 0.845 (671 placed, 104 wrong) on the same fragments.
- **F4:** 18 families with at least 5 placed and none wrong.

| fresh A119b draw | value |
|---|---|
| sources drawn / kept (clear by 30 AS) | 11,421 / 11,280, from 212 families |
| fragments with a truth | 33,840 |
| O2's decisions on almost-tied fragments | validated 277, sets 656, unexplained 859 |

**The 18 families (all members at 90-100% identity):**

| family | member identity | O2 right / wrong | minimap2 at MAPQ >= 1: right / wrong |
|---|---|---|---|
| AC1143242 | 0.991 | 6 / 0 | 8 / 2 |
| AC1378001 | 0.998 | 8 / 0 | 12 / 0 |
| BMS1P1 | 0.985-0.998 | 5 / 0 | 9 / 5 |
| CCT6P3 | 0.988 | 6 / 0 | 9 / 3 |
| CLEC18A | 0.971-0.9996 | 9 / 0 | 13 / 0 |
| GPAT2 | 0.992 | 10 / 0 | 18 / 0 |
| GPR89A | 0.971 | 11 / 0 | 16 / 0 |
| GTF2IRD2 | 0.988-0.993 | 8 / 0 | 23 / 5 |
| HSPA1A | 0.996 | 10 / 0 | 36 / 2 |
| HYDIN | 0.985 | 8 / 0 | 14 / 0 |
| MSTRG596 | 0.986 | 6 / 0 | 8 / 0 |
| NOMO1 | 0.996-0.998 | 5 / 0 | 13 / 1 |
| NOS2P2 | 0.997 | 5 / 0 | 9 / 8 |
| PDE4DIP | 0.911-0.998 | 10 / 0 | 14 / 1 |
| SPATA31A1 | 0.929-0.996 | 10 / 0 | 22 / 0 |
| SPATA31D4 | 0.906-0.977 | 6 / 0 | 11 / 1 |
| UGT2B17 | 0.966 | 6 / 0 | 7 / 0 |
| ZNF658 | 0.988 | 8 / 0 | 19 / 1 |

No qualifying family has an O2 error.

### Audit (metric-auditor, same day): holds with caveats

**Recomputed independently:**
- Truth and O2's rows match exactly: 197 placed, 0 wrong; F4's 18 families.
- The baseline is 670 placed, 103 wrong (0.846). That is one source off, because the row is filed under the primary's
  family. It does not flip.
- By fragment: O2 220 / 220 right; minimap2 at MAPQ >= 1 689 / 801 (0.860).
- The seed-1 and seed-2 sources share no read. The code was unchanged and the pre-registration preceded every output.
  No validation leakage.

**Caveats:**
1. **O2 certifies; it does not overrule.** Its call equals minimap2's primary on all 220 placed fragments, and every one of
   the 112 fragments minimap2 got wrong is one O2 abstained on. So O2 certifies about 29% of minimap2's trusted calls, at
   no error. Nothing simpler does that:
   - no MAPQ cut matches it (MAPQ >= 30: 91 placed, 0.879);
   - neither does a substitution count without the sample's validation (360 placed, 0.908).
2. **Held out by reads only.** The donor, families, PSV sites and validating reads are A119b's, where k = 3 and
   validated-only were chosen. 14 of the 18 F4 families also passed in development. This is a replicate. The other-donor
   test (register 1341) is still UNDERPOWERED.
3. **F4 does not single out O2.** minimap2's MAPQ >= 1 calls meet the same per-family rule in 19 families, and a floor of
   5 cannot bound a family's error rate.
4. **Proxy truth.** The truth is minimap2's call on the full read. A wrong full read would be invisible here; on HG002
   this truth agrees with the donor 98.4-99.8% at this margin.


---

## PREREG_o2_known_family_smn_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_known_family_smn_2026-10-20.md`)*

## Pre-registration: a known family, SMN1/SMN2: does O2's pick of the real primary give better loci than minimap2's, StringTie's and FLAIR's? (2026-10-20)

**User (2026-10-20):** "before moving to O3, can we ensure it works on a known gene family? One that has members that are
similar to each other by >90% but less than 100%, can we outperform minimap2 in picking which is the real primary or can we
reduce the pool from several multimappers to just a handful of more likely candidate reads that we can use to make new loci
and thus outperform stringtie and flair?"

### The family and why it was chosen

**SMN1 / SMN2** (pool family S19SMN2: `.0` = SMN2, chr5:70,808,816-70,837,687 −; `.1` = SMN1, chr5:71,381,817-71,409,805 +).
The two copies are 99.8% identical.

The biology gives a truth that does not depend on any aligner. SMN2 skips exon 7 in most transcripts (Δ7, its c.840C>T), and
SMN1 includes it. A read minimap2 gives to the wrong copy therefore mixes the two copies' exon-7 patterns.

**Not held out by family.** On HG002, which was development data for the corrector, its version-1 and version-2 calls on
SMN's tied reads were already seen: 53% → 99% against the donor assembly.

**What is new:**
- what that does to **loci**: StringTie and FLAIR run on each read's chosen alignment;
- the **exon-7 readout**;
- **A119b's native SMN reads**, never used: the A119b fragment sets cut other reads.

### Methods: which alignment is each read's primary

**Reads:** every read whose minimap2 primary lies on SMN1 or SMN2 (`summary.attribute`).

**Arms:**
- **mm2:** minimap2's primary.
- **O2:** for a read in the corrector's domain (tied within 10 AS on both copies, nothing off the family), the corrector's
  copy; otherwise minimap2's primary. The corrector is `o2correct.py donor`, version 2, commit bfa15f12.
- **truth:** HG002 only. The donor-assembly copy (`hg002proof` truth classes unique / unique_elsewhere), else minimap2's
  primary.

**Building each arm's BAM.** One record per read: its best-AS alignment on the chosen copy, given the read's sequence and
qualities and made primary. A read whose chosen copy has no alignment keeps minimap2's primary (counted).

**Assemblers**, on each arm's BAM (the lab's recipes, `bench/entangled/run.sh`):
- **StringTie** 3.0.3: `-L`, no `-G`;
- **FLAIR** 3.0.0: `bam2bed` + `collapse --trust_ends --generate_map` (the genome `chm13v2.0.fa`).

### Readout: exon-7 inclusion per copy

Exons are taken from the copy table (0-based, half-open):

| copy | exon 6 | exon 7 | exon 8 |
|---|---|---|---|
| SMN2 | 70,815,947-70,816,056 | 70,810,763-70,810,817 | 70,808,816-70,810,319 |
| SMN1 | 71,403,494-71,403,603 | 71,408,728-71,408,782 | 71,409,227-71,409,805 |

**Classifying a read or transcript on a copy.** It **spans** when its blocks overlap exon 6 and exon 8. Then it is:
- **FL** when a block overlaps exon 7 by at least 20 bp;
- **Δ7** otherwise.

**Per copy,** inclusion = FL / (FL + Δ7):
- reads: counted one each;
- StringTie: transcripts weighted by `cov`;
- FLAIR: isoforms weighted by their mapped reads.

### Bars

- **L1 (HG002, reads):** for SMN1 and for SMN2, O2's inclusion is closer to the truth arm's than mm2's is.
- **L2 (HG002, loci):** each assembler is compared with itself on the truth arm. On each copy, the assembler on O2's BAM must
  be no farther from that than the assembler on mm2's BAM, and strictly closer on at least 2 of the 4 (tool, copy) pairs.
- **L3 (A119b, held out, biology):**
  - SMN1's inclusion minus SMN2's is at least as large under O2 as under mm2, for the reads and for both assemblers;
  - under O2, SMN2's inclusion is below SMN1's.
  - The minimap2 coin toss is expected to shrink that gap.
  - **UNDERPOWERED** if fewer than 20 of the corrector's domain reads span exons 6 to 8.

**Verdict:**
- **PASS:** L1, L2 and L3 hold.
- **FAIL:** any of them fails.
- **UNDERPOWERED:** L3 lacks power; L1 and L2 are then reported alone.

**A PASS would claim:** on a known 99.8%-identical family, O2's pick of the real primary turns minimap2's coin-tossed SMN reads
into each copy's own isoforms, and StringTie and FLAIR built on O2's picks get closer to the donor's truth than on
minimap2's.

**It would not claim:** other families, or new copies.

**Rules:**
- The tool is `bench/o3_graph/o2loci.py`, written and tested before any arm is scored.
- Run script: `/mnt/linuxdisk/tmp/o3_graph/o2loci_smn/run.sh`.
- A bug found during the run is fixed and recorded as an amendment before any rerun.

### Amendment 1 (2026-10-20, after the first scoring, before any bar was read)

Two bugs, both found in the first scoring (`o2loci_smn/run_before_amendment1.log`). That scoring classified only 53 of 3,651
HG002 reads and gave no FLAIR numbers at all.

1. **Exon 6 was the wrong exon.** The table above gives the 109-bp exon (SMN2 70,815,947-70,816,056; SMN1
   71,403,494-71,403,603). That is **exon 6B**, the Alu-derived alternative exon in intron 6. SMN exon 6 is the 111-bp exon:
   - SMN2: 70,816,585-70,816,696;
   - SMN1: 71,402,854-71,402,965.

   Exons 7 and 8 are unchanged.
2. **FLAIR's read map** keys each isoform as `isoform_gene`. The scorer now strips the gene suffix to match the GTF's
   `transcript_id`; test `test_o2loci.py`.

The arms' BAMs and both assemblers' outputs do not depend on either fix and are not rerun. Only `score` is rerun. The bars
are unchanged.

### Results (`o2loci_smn/score_amended.log`, `hg002/score.json`, `a119b/score.json`)

**Overall: FAIL.**
- **L1 PASSES:** O2 gives each copy its own isoforms at the read level.
- **L2 FAILS:** StringTie and FLAIR do not keep the assignment they are given.
- **L3 is UNDERPOWERED:** A119b barely expresses SMN2.

| HG002 exon-7 inclusion (FL / (FL + Δ7)) | minimap2 | O2 | truth |
|---|---|---|---|
| **reads, SMN2** | 0.653 (485 / 258) | **0.625** (890 / 533) | 0.617 (836 / 519) |
| **reads, SMN1** | 0.887 (2,279 / 291) | **0.992** (1,874 / 16) | 0.985 (1,928 / 30) |
| StringTie, SMN2 (`cov`) | 0.746 | 0.705 | 0.707 |
| StringTie, SMN1 (`cov`) | 0.929 | 1.000 | 0.424 |
| FLAIR, SMN2 (mapped reads) | 0.763 | 0.765 | 0.700 |
| FLAIR, SMN1 (mapped reads) | 0.897 | 0.880 | 0.934 |

- **L1.** minimap2 hands about 700 of SMN2's reads to SMN1, many of them Δ7. SMN1 then seems to skip exon 7 in 11% of its
  reads, against 1.5% in truth. O2 errs by 0.008 (SMN2) and 0.007 (SMN1), against minimap2's 0.036 and 0.098.

  Per-copy read counts:

  | copy | O2 | truth | minimap2 |
  |---|---|---|---|
  | SMN2 | 1,423 | 1,355 | 743 |
  | SMN1 | 1,890 | 1,958 | 2,570 |
- **L2.** Each assembler is compared with its own output on the truth arm:
  - StringTie on SMN2: O2 is closer (0.002 vs 0.039).
  - StringTie on SMN1: O2 is farther (0.576 vs 0.505).
  - FLAIR on SMN2: O2 is farther (0.065 vs 0.063).
  - FLAIR on SMN1: O2 is farther (0.054 vs 0.037).

  Neither assembler carries a read's copy through:
  - **FLAIR** realigns every read to its own isoforms, which reintroduces the coin toss between SMN1's and SMN2's
    near-identical isoforms. Its per-copy totals barely move between arms.
  - **StringTie** gives SMN1, with about 2,000 reads, transcripts of `cov` 1-19. On the truth BAM it even reports 42%
    inclusion. `longcov` gives the same.
- **L3.** Of A119b's 297 domain reads, 15 span exons 6 to 8 (bar: 20). SMN2 has one spanning read in the whole library, and
  the corrector overruled no read.

**Not pre-registered, descriptive.** Against the truth arm's read-level inclusion (0.617 / 0.985), the absolute errors on
SMN2 / SMN1 are:

| arm | SMN2 | SMN1 |
|---|---|---|
| O2's reads | 0.008 | 0.007 |
| StringTie on minimap2's BAM | 0.129 | 0.056 |
| FLAIR on minimap2's BAM | 0.146 | 0.088 |

To beat them as loci, the copy assignment must survive assembly: an assembler that counts the reads it is given per copy.
Feeding O2's primaries to StringTie or FLAIR is not enough.


### Amendment 2 (2026-10-20, after the audit): StringTie was scored on `longcov`, not `cov`

**The bug.** The audit found that the scorer's `cov "..."` pattern matched inside StringTie's `longcov "..."`, which StringTie
writes first. Every StringTie number above is therefore `longcov`. So are the statements "`cov` 1-19" and "0.424 even on
truth", and the claim that `longcov` gives the same. All of these are wrong.

**The fix.** Attributes are now read by name, and `test_o2loci.py` has a `longcov` line. Only `score` is rerun
(`o2loci_smn/score_amended2.log`).

| HG002 exon-7 inclusion, StringTie (`cov`) | minimap2 | O2 | truth |
|---|---|---|---|
| SMN2 | 0.654 | **0.628** | 0.618 |
| SMN1 | 0.894 | **1.000** | 0.989 |

**L2, corrected:**
- **StringTie on O2's BAM is closer to its truth arm on both copies** (SMN2 0.010 vs 0.036; SMN1 0.011 vs 0.095). StringTie
  keeps the copy it is given.
- **FLAIR is farther on both** (0.065 vs 0.063; 0.054 vs 0.037), so **L2 still FAILS, on FLAIR alone**. FLAIR re-maps reads
  to its own isoforms: it puts 3,285 of 3,443 reads on the same copy in the O2 and minimap2 arms, and returns 367 of O2's 682
  moves to SMN1. By design it cannot register an input assignment.
- **The verdict stays FAIL.**

**Also disclosed by the audit:**
- **L1 is development data** (O2 1,289 / 1,299 right on SMN's tied reads with donor truth). In A119b the corrector overruled
  no read, so this test has no held-out evidence for SMN.
- **The first scoring (wrong exon 6) failed L1 on SMN1** (0.110 vs 0.096). Amendment 1 came 54 s after that run.
- **The truth arm falls back to minimap2's primary for 751 of 3,651 reads** (633 poor_match, 114 own_tie, 4 without a row).
  This biases it toward minimap2. On the 2,900 donor-truth reads alone, L1 holds tighter: O2 0.599 / 0.994, truth
  0.598 / 0.996, minimap2 0.622 / 0.876.
- **At c.840**, the truth puts all 756 T reads on SMN2 and all 1,544 C reads on SMN1.
- **The post-hoc line above mixed units.** Like for like, StringTie (`cov`) errs on SMN2 / SMN1 by:
  - 0.037 / 0.091 on minimap2's BAM;
  - 0.011 / 0.015 on O2's BAM;
  - 0.001 / 0.004 on the truth BAM.

  FLAIR errs 0.083 / 0.051 even on the truth BAM.
- **Minor:**
  - The per-copy counts above are reads spanning exons 6-8. All SMN2 reads number 1,582 / 1,518 / 902 (O2 / truth /
    minimap2).
  - A119b's SMN domain is 275 reads (22 more are SERF1B's).
  - SMN2 skips exon 7 in 38% of its truth reads, not "most".


---

## PREREG_o2_personal_copies_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_personal_copies_2026-10-20.md`)*

## Pre-registration: personal copies against the corrector, on gorilla KB3781's own reads with exact truth (2026-10-20)

**User (2026-10-20):**
- "lets go back to O2, and assume we have a family with members with different degrees of similarity (90 to slighly less
  than 100%); what is still missing to be able to disambibuate (ambiguous AS) to know the real provenance and use that for
  assembly"
- On the proposal to test per-copy sample sequences against the corrector's pileup: "yes please, also maybe we can use
  gorilla for this since we have the fibroblasts and T2T genome?"
- On the design below: "yes please proceed".

### The question

The corrector (`o2correct.py donor`, registers 1343-1346) learns each copy's sample alleles as a **pileup**: per copy, how
many of its confident reads carry each substitution or indel. A tied read's events are then checked one by one against
that pileup.

A **personal copy** is the copy's own sequence in the sample: the reference copy with the events its confident reads
carry written in, split into at most two haplotypes by the reads that carry them. A tied read is realigned to the personal
copies of its tied copies.

The personal copy keeps what the pileup loses: which events sit together on one haplotype, and the read's alignment to the
sample's sequence rather than to the reference.

The experiment asks:
1. Do personal copies place tied reads better than the pileup?
2. How close does either get to the animal's real copies?

### Substrate: one animal, his reads and his own genome

- **The reads.** KB3781 ("Jim", the animal of the mGorGor1 T2T assembly): fibroblast Iso-Seq, SRR27178662, 13,516,443
  molecules and 34,762,361 records. They were aligned to the RefSeq primary GCF_029281585.2 (`GGO.fasta`) by minimap2 2.31
  `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`
  (`/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam`).
- **The reference is a mosaic of his two haplotypes** (`rna_allele/chrmap.tsv`, sha1 e9182791, each chromosome checked
  identical in sequence):
  - his father's: chromosomes 1, 2 (NC_086017.1), 3, 4, 7, 8, 9, 12, 13, 14, 15, 16, 18, 19, 21 and Y;
  - his mother's: 5, 6, 10, 11, 17, 20, 22, 23 (NC_086018.1) and X.
  - This corrects the design note in chat, which said the autosomes were all his father's.
- **The ALT haplotype** is the other parent's copy of each autosome, plus every unplaced contig of either assembly. His
  assemblies are GCA_028885475.2 (paternal) and GCA_028885495.2 (maternal), in `gorilla_haps/{pat,mat}.fa`.
- **Two built-in groups of reads.** A read from the haplotype the reference was taken from has a sample copy identical to
  the reference: the **one-genome control**. A read from the ALT haplotype carries the other parent's alleles: the **test**.

### The families

- **Catalog.** The KB3781 fibroblast catalog GWFAM: `rna_allele/fibro_cat.copies.tsv`, sha1 7521851f, 666 families, 3,981
  copies.
- **Rule.** A family qualifies when every pair of its members' main transcripts aligns at identity in [0.90, 1.00), covering
  at least 50% of the shorter. The test is `minimap2 -X -c -x asm20` with `o2perfamily.qualifies`, the human list's rule.
  - It is run one family at a time (`o2personal.py ava`), since only same-family pairs count.
  - A family of more than 20 members is first tested on 20 of them; a pair failing there fails the family.
- **Result.** 283 families qualify, with 716 copies (`personal/qual.json`, sha1 f4bab888; `qual.copies.tsv`, sha1 78a3afa8).
  NPIP does not qualify.
- **The hold-out is by family.** The parity of sha1(family name), fixed before any result was seen (`o2personal.half`):
  - development: 140 families;
  - held-out: 143 families.

  | minimum pair identity | development | held-out |
  |---|---|---|
  | [0.90, 0.95) | 15 | 11 |
  | [0.95, 0.98) | 50 | 64 |
  | [0.98, 0.99) | 37 | 33 |
  | [0.99, 1.00) | 38 | 35 |

### The domain (the corrector's, unchanged)

A read is in the domain when:
- its alignments within 10 AS of its best lie on two or more copies of one qualifying family;
- none of them lies on a locus off those copies;
- minimap2's primary is among them.

The reads are every alignment of the 450,984 reads with any alignment on a qualifying copy (`personal/sub.bam`, 4,055,491
records, `names.txt` sha1 26dd6d3c).

### The truth: the read's best placement in the animal's own genome

**Alignment.** Each domain read is aligned to his maternal and his paternal assembly (`mGorGor1.{mat,pat}.splice.mmi`):
minimap2 2.30, `-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`.

**Each alignment maps to a reference copy:**
- On a haplotype contig the reference is identical to, it maps to the copy it overlaps there, in the reference's
  coordinates (`summary.attribute`, the corrector's own attribution).
- On an ALT contig, it maps to the copy whose ALT ortholog it overlaps. A copy's ALT ortholog:
  - its reference span ± 20 kb, aligned with `minimap2 -x asm20 -c` to the other haplotype's whole assembly;
  - the longest primary alignment with MAPQ >= 20 through which its exons lift (`donorsim.pick_ortholog`);
  - on the homologous chromosome or an unplaced contig;
  - two copies whose orthologs overlap are both dropped (`donorsim.drop_ambiguous`).

**The truth** (`o2personal.truth_of`) is the one copy that every alignment within 5 AS of the read's best lies on. 5 AS is
one substitution under splice:hq (match 1 + mismatch 4).

| status | meaning |
|---|---|
| determined, **refcons** | The reference's own copy is among them: the sample's copy equals the reference over the read. |
| determined, **altonly** | Only the ALT haplotype's copy: the read carries the other parent's alleles. |
| undetermined | Two copies, or a copy and an off-copy locus, lie within 5 AS. No method can place the read. |
| off-copy | The best placement is on no reference copy: a copy the reference lacks, an unannotated locus, or an ortholog that did not lift. |

**Scored reads (S):** determined, with the truth among the read's tied copies.

### The arms

Each arm gives one copy per read. An arm that abstains keeps minimap2's primary.

1. **minimap2**: its primary.
2. **Corrector**: `o2correct.py donor`, version 2 without `--guard`, unchanged since bfa15f12.
3. **Personal copies, phased** (`o2personal.py`). Per copy:
   - **Reads:** the corrector's own confident reads (primary on the copy, every alignment on another copy of its family more
     than 9 AS lower).
   - **Events:** those the corrector calls explained (>= 2 reads, >= 10% of a depth >= 3). Events at >= 80% of the depth are
     **carried**.
   - **Haplotypes:** the other events split into at most two haplotypes by a 2-means over the reads (`o2personal.phase`).
   - **Sequence:** each haplotype is the reference span ± 5 kb with its events written in.
   - **Decision:** each tied read is realigned to the personal copies (the same minimap2 options, `-N 200`). The tied copy
     with the unique best AS takes it. A shared best, or a tied copy the read does not align to, keeps minimap2's pick.
4. **Oracle**: the same decision over the copies' true sequences: the reference span ± 5 kb and its ALT ortholog ± 5 kb.
5. **Exploratory, no bar**: personal copies unphased, one sequence per copy (the events more than half the depth carries).

### Scoring (`o2personal.py score`)

**Per arm, on S:**
- right and accuracy;
- against minimap2: **fixes** (minimap2 wrong, arm right) and **breaks** (minimap2 right, arm wrong);
- reads decided, and decided right.

**Strata:**
- refcons and altonly;
- exact ties (AS margin 0 between the two best tied copies) and near ties;
- the family's minimum pair identity, in the four bins above;
- each family.

**Also reported:** undetermined, off-copy, and truth outside the tie set.

### Bars, held-out half

- **G0 (gate):** oracle accuracy on S >= 0.99. Failing it means a pipeline error: INVALID.
- **P (power):** |S| >= 300, from >= 11 families, with minimap2 wrong on >= 50 of them.

**V1. Personal copies beat the corrector:**
- **G1:** personal accuracy >= the corrector's + 1.0 point.
- **S1:** personal fixes >= 5 × personal breaks.
- **F1:** >= 11 families where personal is right more often than the corrector, and no more than a fifth as many (rounded
  up) where it is right less often.
- **H1:** on refcons reads, personal accuracy >= minimap2's − 0.2 points.

**V2. The corrector beats minimap2** (secondary; its first real-read test with exact truth):
- **G2:** corrector accuracy >= minimap2's + 2.0 points.
- **S2:** fixes >= 5 × breaks.
- **F2:** >= 11 families where it fixes more than it breaks, and no more than a fifth as many (rounded up) where it breaks
  more than it fixes.
- **H2:** on refcons reads, corrector accuracy >= minimap2's − 0.2 points.

### Verdict

- **INVALID:** G0 fails.
- **UNDERPOWERED** (both V1 and V2): G0 holds and P fails.
- **PASS:** all four of that verdict's bars hold. **FAIL:** any of them does not.

**A V1 PASS would claim:** on one gorilla's own fibroblast reads, with exact truth from his own assemblies, the sample's
phased copy sequences place tied reads better than the corrector's pileup. This holds in more than 10 near-identical
families and leaves reads whose copy equals the reference unharmed.

**It would not claim:**
- other tissues or animals;
- copies absent from the reference (off-copy: O3);
- reads tied in the sample itself (undetermined);
- a better assembly (no second round here).

### Development, then the freeze

The development half (140 families) may be run and looked at freely. On it, I may:
- fix bugs;
- change the personal copies' construction, the realignment's minimap2 options and FLANK_P;
- change the truth's mechanics (contig names, ortholog lift).

Fixed for the whole experiment:
- the truth margin (5 AS);
- the domain;
- the corrector (bfa15f12);
- the bars.

Every change goes into Amendment 1, with its reason, at a commit. The held-out half is scored once, with the code at that
commit.

### Rules

- Gorilla numbers are never pooled with human ones.
- Heavy jobs run in the foreground, one at a time: the truth alignments (about 14 GB of index each) and the ALT ortholog
  alignments.
- The metric-auditor checks the held-out numbers before any verdict is reported.
- Not run, and kept as a later replication: the mirror, with reference = the other haplotype.

### Amendments

#### Amendment 1 (2026-10-20): allele-blind O2, and an edit-based truth. Written before any arm was scored against it

**Why.**

1. **The user stopped the phased arm:** "I think O2 should not care about alleles, phasing, paternal and maternal
   haplotypes. It should only be about finding the true place of multimappers based on PSVs, BAM tags and other info".
   - Dropped: the personal copies (arms 3 and 5) and the oracle (arm 4, his assembled copies as targets), with gate G0.
   - His two assemblies stay as the **answer key only**. No arm sees them.
2. **AS was the wrong yardstick for the truth.** AS charges every intron, so a copy with fewer introns (a retrocopy) can
   outscore the spliced gene the read matches perfectly. This was seen on the development half only:
   - Of its 3,484 domain reads, 2,592 have minimap2's primary below the highest-AS tied copy. In 2,539 of them the primary
     has fewer edits (NM); in 2,540 it has more introns.
   - Example: SRR27178662.117798 matches SM438.1 (minimap2's primary, AS 1,196) with NM 0, and SM438.0 (AS 1,204, fewer
     introns) with NM 24.
   - The AS truth (5-AS margin) therefore called minimap2 wrong on 740 development refcons reads like these.
   - minimap2's primary has the fewest NM of the tied copies on 3,405 of the 3,484.

**What had been seen before this amendment:**
- **Development:**
  - minimap2's primary against the AS truth, by origin and margin;
  - the AS / NM / intron comparison above;
  - the corrector's output (`arm2.tsv`), scored against no truth.
- **Held-out:** only the AS truth's status counts: determined altonly 404, determined refcons 215, off-copy 7,
  undetermined 270. No arm's calls were compared with any truth.
- **Personal copies:** built (`personal.*.fa`), never realigned or scored. `o2personal._depths` was sped up (bisect),
  without changing its output.

**The truth** (replaces the section above). The alignments, and their mapping to reference copies, are unchanged
(`truth.mat.bam`, `truth.pat.bam`, `orthologs.tsv`).
- **An alignment's edits:** its NM plus its clipped read bases (soft and hard), i.e. the read bases it gets wrong or does
  not explain. Introns are free.
- **The truth** (`o2personal.truth_of`, now on edits): the one copy whose alignments hold the fewest edits, strictly fewer
  than every alignment on another copy or off every copy. Edits are whole numbers, so that is at least one edit fewer.
- **Statuses:**
  - **determined, refcons:** a reference-haplotype alignment holds the fewest edits;
  - **determined, altonly:** only an ALT one does;
  - **undetermined:** two copies, or a copy and an off-copy locus, share the fewest edits;
  - **off-copy:** the fewest edits lie off every reference copy.
- **S** is unchanged: held-out, determined, the truth among the read's tied copies.

**The arms** (replace the section above). Each gives one copy per read; a shared best keeps minimap2's pick.
1. **minimap2:** its primary.
2. **O2, fewest edits (the main arm):**
   - each tied copy is represented by its highest-AS alignment (the corrector's choice);
   - its edits are its substituted, inserted and deleted bases (`o2correct.alignment_events`, indels shifted left), counted
     over the read bases every tied alignment aligns (`o2correct.common_window`);
   - the unique fewest takes the read.

   Errors and alleles every tied copy shares cancel, so the difference is the PSVs.
3. **Fewest substitutions:** `o2correct.py run` (X ops, or NM less its indel bases, over the whole alignment).
4. **Highest AS (negative control):** the unique highest AS among the tied copies.
5. **Corrector:** `o2correct.py donor`, bfa15f12. A comparison only, with no bar.

**Strata** (one added): **structure**, i.e. whether the tied copies' alignments (each copy's highest-AS one) differ in
intron count ("intron-differ") or not ("same structure").

**Bars** (replace the section above). On the held-out half, O2 (arm 2) against minimap2:
- **P:** |S| >= 300, from >= 11 families, with minimap2 wrong on >= 50.
- **G:** accuracy >= minimap2's + 2.0 points.
- **S:** fixes >= 5 × breaks.
- **F:** >= 11 families where it fixes more than it breaks, and no more than a fifth as many (rounded up) where it breaks
  more than it fixes.
- **H:** on refcons reads, accuracy >= minimap2's − 0.2 points.

**Verdict:** PASS = all five hold; FAIL = P holds and another bar fails; UNDERPOWERED = P fails. This one verdict replaces
V1 and V2. Arms 3-5 are reported with the same numbers, without verdicts.

**Development and freeze** are unchanged: changes go into Amendment 2, at the freeze commit. Fixed from here on:
- the edit rule of the truth;
- the domain;
- the arms' definitions;
- the bars.

#### Amendment 2 (2026-10-20): the development half, one secondary arm, and the freeze. Written before the held-out half was scored

**The development half, scored under Amendment 1** (`personal/score_dev.json`, 140 families):
- Domain statuses: determined refcons 1,805; determined altonly 1,661; truth outside the tie 7; off-copy 1; undetermined 10.
- Scored reads S: 3,466 from 31 families.
- The truth is robust: the runner-up copy has >= 5 more edits for 3,270 of them, 2-4 for 106, exactly 1 for 90.

| arm | accuracy | fixes / breaks against minimap2 |
|---|---|---|
| minimap2 | 0.9833 (58 wrong) | – |
| O2, fewest edits | 0.9645 | 4 / 69 |
| fewest substitutions | 0.9786 | 0 / 16 |
| highest AS (control) | 0.2726 | 4 / 2,467 |
| corrector | 0.9815 | 0 / 6 |

**What it showed:**
- **AS is the wrong yardstick between copies.** 3,220 of the 3,466 scored reads are ties between copies whose alignments
  differ in intron count. minimap2 is right on 99.97% of them; the highest-AS control is right on 24.7%.
- **minimap2's 58 errors are almost all beyond any allele-blind rule.** For 48 of them the true copy has more edits than
  minimap2's pick against the reference, and for 6 it has as many. Those reads carry his own alleles. Only 4 have the true
  copy with fewer edits; O2 fixes those 4.
- **O2's 69 breaks are one family's reads, not a bug.** 65 are SM251 (chr1, 3 copies at [0.93, 0.98]):
  - The reads match his maternal SM251.0 with 1-2 edits, inside SM251.0's ALT ortholog (MAPQ 60).
  - That maternal copy differs from the reference's own SM251.0 by about 41 edits per 2.5-kb read, so against the reference
    the reads are closer to SM251.1 (34 edits against 43).
  - It is a copy the reference effectively lacks (O3's case), which no allele-blind rule can place. minimap2's pick of
    SM251.0 is right.

**Changes:**
1. **Fewest substitutions ran on a smaller BAM.** `o2correct.py run` on `sub.bam` exceeded the 10-minute limit (it
   attributes all 4.06M records). It was run unchanged on `domain.bam`, every alignment of the 4,380 domain reads cut from
   `sub.bam`. All 4,380 reads got the same tie sets and primaries as `arm2.tsv`.
2. **Truth output.** The printout shows the development half only. A descriptive `runner_up` column was added (the fewest
   edits on any other copy or locus).
3. **A secondary arm, designed after the development result: O2, consistent.**
   - It takes O2's pick only where the read's bases contradict it nowhere in the common window. A contradiction is a read
     base that O2's copy mismatches (an X in its `--eqx` alignment) and that a rival tied copy's alignment matches (an =).
   - Otherwise it keeps minimap2's pick.
   - This is the O2 regimen's own principle, "its own bases contradict the call nowhere" (register 1342). In SM251 the
     maternal copy is a mosaic of SM251.0 and SM251.1 bases, so it vetoes those moves.
   - It gets the same five bars and its own verdict, secondary to O2's.
   - Development: 2 overrules, 0 fixes / 1 break, accuracy 0.9830. It vetoes O2's 4 fixes too.

Unchanged: O2's main arm, the truth's edit rule, the domain, the bars, and the primary verdict (O2, fewest edits, against
minimap2).

**Freeze.** The held-out half is scored once with the code at the commit that adds this amendment. Before scoring, every
output is rebuilt with that code and compared with the development run: `truth_edits.tsv`, `arms.tsv` and `run.tsv`.

### Results: the held-out half, scored once at 00edc953 (`personal/score_heldout.json`)

**Domain:** 896 held-out reads tied within 10 AS. 887 are determined, all with the truth among the tied copies (refcons 479,
altonly 408); 3 are undetermined, 6 off-copy.

**S:** 887 reads from 35 families; minimap2 is wrong on 91, so P holds.

| arm | accuracy | fixes / breaks against minimap2 | families fixed > broken / broken > fixed |
|---|---|---|---|
| minimap2 | 0.8974 | – | – |
| **O2, fewest edits (primary)** | **0.8985** | **2 / 1** | **2 / 1** |
| O2, consistent (secondary) | 0.8985 | 1 / 0 | 1 / 0 |
| fewest substitutions | 0.8884 | 0 / 8 | 0 / 1 |
| highest AS (control) | 0.6077 | 68 / 325 | 1 / 22 |
| corrector (descriptive, no bar) | 0.9899 | 82 / 0 | 2 / 0 |

**O2's bars:**
- P holds;
- **G fails:** +0.11 points, against +2.0;
- **S fails:** 2 fixes against 5 × 1 break;
- **F fails:** 2 families, against 11;
- H holds: −0.21 points.

**Verdict: FAIL by the pre-registered bars.** O2 changed 3 of minimap2's 887 calls: SM556 and SM456 fixed, SM463 broken.
All three are noisy reads, with 1.8-4.6% edits against a median of 0.29%. Substantively that is **no measurable effect**.
O2, consistent: FAIL (G, F).

**Two caveats on the bars** (audit):
- **No arm could have passed F.** minimap2's 91 errors lie in 5 families (SM392 82, SM502 6, SM556, SM456, SM277 1 each), so
  even an oracle would be better in 5 families against the 11 F needs. On the development half an oracle would gain +1.67
  points, short of G. Amendment 1 dropped the oracle arm and its gate, which would have shown this before the freeze. It is a
  design flaw of the bars, as in register 950.
- **P holds only through SM392.** Without it, S is 805 reads with 9 minimap2 errors: UNDERPOWERED.

**Where minimap2 errs:**
- **SM392 (82 errors):** two copies 4 kb apart on chr5 (0.989 identity), on a chromosome where the reference is his mother's.
  - Every scored read comes from his paternal SM392.1, whose ortholog is sound: order and spacing are kept (4,179 vs
    4,172 bp), and 0.913 of the reference SM392.1's 31-mers are in it. It differs from the reference's SM392.1 by a median 21
    edits per read.
  - On the reference, every one of the 82 primaries is a **two-copy chimera**: about 2.0 kb on SM392.0, a 12,220-bp "intron",
    then 475 bp on SM392.1. `summary.attribute` credits it to SM392.0.
  - minimap2's own AS favours SM392.1 on 67 of the 82, and the highest-AS control fixes those 67.
  - O2 does not fix them. Its edit count treats the 12-kb gap as a free intron, and the chimera stands for SM392.0.
  - So these reads are not beyond an allele-blind rule: AS, which charges for the gap, places most of them.
- **SM502 (6 errors):** fixed by the corrector only.
- **The other 3 errors:** one read in each of 3 families.
- **By origin:** 89 of the 91 errors are altonly reads (minimap2 0.7819 on them), 2 are refcons (0.9958).

**Structure of the ties.** 841 of the 887 scored reads tie copies whose alignments differ in intron count. Ranking by AS is
wrong on 39.2% of the scored reads (development: 72.7%), yet right on SM392's chimeras. Exact ties (AS margin 0): 29 reads,
right for every arm.

**The refcons stratum is close to self-scoring** (audit). On all 479 held-out refcons reads (and 1,805/1,805 in
development), the truth equals the copy with the fewest edits among the library's own reference alignments, which is what O2
minimises. minimap2's 0.9958 there, and H, hold largely by construction. The two assemblies add information on about 96 of
the 887 scored reads: SM392's 82 and 14 ties.

**How many reads his genome cannot place.** 3 of the 896 domain reads tie exactly under the strict-minimum rule (development:
10 of 3,484):
- only one (SM491) is a tie between two catalogued copies;
- the SM357 read ties with SM357.0's own maternal ortholog, which the ortholog lift missed (the maternal assembly carries
  0.95 of SM357.0's 31-mers at CM054600.2 ~5.07-5.11 Mb); 5 of the 6 off-copy reads are the same gap;
- the SM339 read ties with an uncatalogued paternal locus.

Another 24 scored reads are decided by a single edit and 12 by two. In register 1338, one-edit margins were 46.5% wrong on
HG002. Requiring a margin of at least 2 or at least 5 leaves the verdict unchanged (O2 2 / 1; corrector 82 / 0 and 76 / 0).

**The corrector is descriptive only.** It has no bar after Amendment 1, and it learns the sample's alleles, which is outside
the user's allele-blind definition of O2. 76 of its 82 fixes are SM392. Its 0 breaks depend on clipped bases counting as
edits: with an NM-only truth it is 78 / 3.

**The reads.** The domain reads come from four runs of the same animal: SRR27178662 (147), SRR27178663 (306), SRR27438212
(3,804) and SRR27438213 (123), all SAMN04003007 = KB3781. The substrate section names only the first.

**Reading.** By the pre-registered bars, O2's fewest-edits rule FAILS on held-out reads, with no measurable effect: 3 calls
changed in 887. Three things follow:
- For reads whose copy matches the reference, minimap2's pick is already the fewest-edit copy almost always, so an
  allele-blind rule has little to change there.
- minimap2's held-out errors are mostly one family, where its primary is a chimera spanning two copies. AS (which charges for
  the spurious intron) places most of them; the edit count, which makes every intron free, does not.
- The bars could not be passed by any arm (F), and the power bar holds only through that family.

### Audit (metric-auditor, 2026-10-20): holds with caveats

- **Every number reproduces to the digit**, from an independent rebuild of the truth (4,380 / 4,380 statuses, truths and
  origins agree) and an independent scorer. FAIL is the pre-registered verdict.
- **The hold-out order is clean:** prereg 725aafb2 → Amendment 1 99595983 → development score 19:00 → freeze 00edc953 at
  19:01:11 → held-out score 19:01:57. The frozen code matches 00edc953, and the frozen TSVs are byte-identical to the
  development-run files.
- **The secondary arm came from development only.** One point the auditor could not settle: whether any run printed held-out
  statuses under the edit truth before the freeze. The session shows the edit truth was first run with the
  development-only printout, so they were first seen at the held-out scoring (counts only in any case).
- **Corrections taken into the Results above:**
  - F was out of reach for any arm, and P holds only through SM392;
  - SM392's minimap2 errors are two-copy chimeric primaries that AS places, so not "beyond any allele-blind rule" and not
    O3's case;
  - "No allele-blind arm fixes any altonly read" was false: highest AS is +67 / −93 on altonly reads;
  - refcons is close to self-scoring;
  - "3 truly ambiguous" became 3 exact ties under the strict-minimum rule, of which 1 is a two-copy tie, plus 24 one-edit
    decisions;
  - the corrector's 0 breaks depend on counting clips;
  - four sequencing runs, not one.
- **Wording not to use:**
  - "O2 is refuted on real reads";
  - "P holds, so the test was powered";
  - "minimap2's errors are beyond any allele-blind rule";
  - "only the corrector fixed them";
  - "exact truth" without qualification;
  - any pass wording for the corrector.
- **Unclear, low risk:** the splice indexes were built before `gorilla_haps` was downloaded. Contigs map by unique length (24
  / 24 paternal, 223 / 225 maternal, no alignments on the 2 ambiguous ones), but sequence identity is not proven.
- The auditor's scripts are in `/mnt/linuxdisk/tmp/audit_2026-10-20_o2personal/`; the frozen copy is
  `personal/audit_frozen_00edc953/`.


---

## PREREG_o2_second_round_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_second_round_2026-10-20.md`)*

## Pre-registration: the second round -- loci rebuilt where O2's corrected multimappers have a say (2026-10-20)

**User (2026-10-20):** "it needs to make loci where the corrected multimappers now have a say, that is why I suggest it as
a post-assembly second round". Then, to the proposal of a held-out test on other known families: "yes please".

### The second round

This is the closed loop of `docs/archive/2026-09/PREREG_tied_read_loop_2026-09-25.md`, with O2's corrector in place of the
union certificate. The certificate assigned 0 reads, so the loop never changed anything.

- **Pass 1** is the shipped `assemble` stage. Loci are seeded from primaries plus secondaries with AS >= 0.98 x the read's
  best. A read tied between two copies therefore echoes at both.
- **The home table** (`bench/loop_home.py o2 --every-read`, b15b5948) gives every read whose primary lies on a family copy one
  home:
  - O2's copy, when the corrector (`bench/o3_graph/o2correct.py donor`, version 2) decided it;
  - minimap2's primary copy, when the read is not tied within 10 AS (Eichler: minimap2's call stands);
  - no home for a tie the corrector could not break: it keeps pass 1's echo.
- **Pass 2** is the same assembly with `RUSTLE_READ_HOME_TABLE`. A homed read enters the read pool only through its records at
  home.

### Development: every number seen before this file

All of it is HG002, chr5, SMN1/SMN2, scored as exon-7 inclusion of the assembled loci against the donor truth (0.617 / 0.985),
StringTie read on `cov` (Amendment 2 of register 1344):

| arm | SMN2 | SMN1 |
|---|---|---|
| pass 2 with O2 homes | 0.585 | 0.969 |
| primaries only | 0.602 | 0.801 |
| pass 1 | 0.743 | 0.744 |
| StringTie | 0.654 | 0.894 |
| FLAIR | 0.763 | 0.897 |

SMN is excluded below. The corrector's read-level calls were seen in development for every HG002 family. What is held out
here is the **loci** of the families below: no second-round output for any of them has been looked at.

### Families (fixed by a power rule, not by results)

Every qualifying non-PAR family (`famident/qualify_nonpar.json`: every member pair at [0.90, 1.00) identity) with at least 20
of the corrector's domain reads carrying a donor truth (unique / unique_elsewhere), SMN excluded. There are 15:
AC0064533, AC1275021, DUSP22, FCGR2A, GOLGA8A, GTF2H2C, GTF2IRD2, HSPA6, NCF1, NOMO1, PDPR, SERF1B, TBC1D3B, THOC3 and
ZC3H11A.

File: `/mnt/linuxdisk/tmp/o3_graph/famident/second_round_families.json`, sha1 6a6624394d9c.

### Arms (HG002 Kinnex, `hg002proof/chm13.bam`; CHM13 v2.0)

| arm | what |
|---|---|
| **LOOP-O2** | pass 2 with the home table above |
| **P** | primaries only (`--seed-pool primary`): minimap2's one home for every read |
| **GOOD** | pass 1, the shipped assembly |
| **ST** | StringTie 3.0.3 `-L` on the primaries (`-F 2308`) over the 15 families' copy spans |
| **FL** | FLAIR 3.0.0 `bam2bed` + `collapse --trust_ends --generate_map` on the same primaries |
| ORACLE | pass 2 with the donor truth's homes (`--every-read --truth`, never O2's calls); reported, not in a bar |

- **Assembler.** The Rustle arms use `tools/rustle_pipeline.sh assemble` genome-wide with one best-AS table, binary
  `rustle_target/release`. The binary is unchanged in `src/` since it was built.
- **Corrector.** Its output is regenerated at this commit and must equal `o2psv_dev/hg002.v2.tsv` (audit of register 1344).

### Metric

**Truth shares.** Per family, the family's reads whose donor-truth class is unique / unique_elsewhere, counted by their true
copy (`truth.tsv` + `truth_untied.tsv`).

**Arm shares.** Per family, the arm's assembled transcripts on each copy (`summary.attribute`: strand, then exon overlap),
weighted:
- `reads` for the Rustle arms;
- `cov` for StringTie, read by name;
- mapped reads for FLAIR.

**Per-family score.** The total variation distance (TVD) between the arm's shares and the truth's. An arm with no transcript
on the family scores 1. Scorer: `bench/o3_graph/o2second.py score`.

### Bars

- **S1:** the pooled TVD (the mean over families, weighted by truth reads) of LOOP-O2 is below that of P, GOOD, ST and FL,
  each.
- **S2:** LOOP-O2's TVD is lower than P's in more families than it is higher, and the same against GOOD.
- **S3 (power):** at least 5 families have truth reads.

**Verdict:**
- **PASS:** S1 and S2 hold.
- **FAIL:** one of them fails.
- **UNDERPOWERED:** S3 fails.

The unweighted family mean is reported beside the pooled one (THOC3 carries the most reads).

**A PASS would claim:** on 15 near-identical families the corrector never had its loci judged on, rebuilding loci after O2
gives each copy its own reads. The loci are closer to the donor's truth than the shipped assembly, minimap2's primaries,
StringTie or FLAIR.

**It would not claim:**
- another sample: HG002 is the only one with per-read truth, and a simulation follow-up is separate;
- transcript structures, since the metric is the read share per copy;
- copies the reference lacks (O3).

**Rules:**
- Run script: `/mnt/linuxdisk/tmp/o3_graph/o2second/run.sh`.
- A bug is recorded as an amendment before any rerun.

### Results (`o2second/run.log`, `o2second/hg002/score.json`)

**FAIL: S1 and S2 fail. S3 holds (15 families, 29,965 truth reads).** The corrector's output at this commit is identical to the
development file (`cmp`). Pass 2 kept 597,640 records at home and dropped 57,294 (O2 moved 2,220 reads off minimap2's copy). Erratum (audit): this line first quoted ORACLE's totals, 591,438 / 56,395.

| arm | pooled TVD (truth-read weighted) | unweighted mean TVD | LOOP-O2 better / worse, families |
|---|---|---|---|
| **LOOP-O2** (pass 2, O2 homes) | 0.112 | 0.266 | — |
| P (primaries only) | 0.127 | 0.283 | 7 / 4 |
| GOOD (pass 1, shipped) | 0.168 | 0.273 | 7 / 7 |
| ST (StringTie) | **0.091** | 0.223 | 6 / 8 |
| FL (FLAIR) | 0.097 | 0.251 | 7 / 7 |
| ORACLE (pass 2, true homes) | 0.108 | 0.215 | 4 / 7 |

**S1.** LOOP-O2 beats P and GOOD but not ST or FL.

**S2.** LOOP-O2 beats P (7 / 4) but ties GOOD (7 / 7).

**Even perfect placement (ORACLE, 0.108) does not reach StringTie (0.091), so the gap is in the assembled loci, not in O2's
assignment.** LOOP-O2 is within 0.004 of the ORACLE pooled, and equal to it in 5 families. Read from the per-family shares:

- **Our assembler drops a locus StringTie keeps.** FCGR2A.1 (chr1:160,843,885-160,862,793) has a StringTie transcript of `cov`
  395 and no transcript in any Rustle arm. Every arm scores about 0.97 there anyway, because the truth puts the family's reads
  on FCGR2A.2.
- **The donor's copies are not CHM13's.** In FCGR2A and TBC1D3B the truth copy is not where any aligner puts the reads. ORACLE
  scores 0.97 and 0.51; every other arm, 0.89-0.98. These are copy-number-variable families: O3's case.
- **The echo is real, and the second round removes it.** DUSP22: truth 95 / 5; pass 1 echoes at 50 / 50 (TVD 0.448); pass 2
  0.018. ZC3H11A: pass 1 0.169, pass 2 0.000. NOMO1: 0.223 → 0.142. PDPR: 0.405 → 0.238.
- **Where pass 1 happens to match the truth ratio** (THOC3 0.004, SERF1B 0.027, AC1275021 0.037), the second round moves it
  away: THOC3 0.026, SERF1B 0.349, AC1275021 0.167.
- **The Rustle arms' read weights are a fraction of the truth reads.** THOC3 has 1,256 reads in transcripts against 5,962
  truth reads, so its shares rest on fewer reads than StringTie's `cov`.

**Not claimed:** that the second round's loci share a family's reads among its copies better than StringTie or FLAIR. It
does better than the shipped pass 1 and than minimap2's primaries alone, pooled.


### Audit (metric-auditor, 2026-10-20): holds with caveats; the explanation above is corrected

Every per-family TVD, both means and the better/worse counts recompute exactly. The tie counts are stable under 16 hash seeds,
and the power rule regenerates the same 15 families. **The FAIL stands, and no variant flips S2: against GOOD it is 7 / 7 or
7 / 8.** The paragraph "the gap is in the assembled loci, not in O2's assignment" does not hold as written:

1. **The weights, not the loci, decide S1.** The share of a family's reads each arm credits to its transcripts is:
   - Rustle's `reads` attribute: 0.22-0.23 (GOOD 0.32);
   - StringTie's `cov`: 0.67;
   - FLAIR: 0.96.

   The Rustle shortfall is uneven across copies (NCF1.1: 417 credited against 3,620 truth reads; StringTie 3,146). With one
   common weight, each arm's own read placement counted at the copies where it built a transcript, the pooled TVD is:

   | arm | TVD |
   |---|---|
   | LOOP-O2 | 0.055 |
   | ORACLE | 0.052 |
   | ST | 0.066 |
   | P | 0.068 |
   | FL | 0.068 (its own read map: 0.097) |
   | GOOD | 0.126 |

   S1 would pass; S2 still fails. This is post hoc, not a result. **The S1 gap is how our assembler credits reads to its
   transcripts, not which loci it builds.**
2. **ORACLE is not perfect placement, and it fails both bars itself**, so the test had no feasibility gate. Its homes cover
   only reads with a donor-truth class: 8-78% of each family's reads stay at minimap2's copy or echo. 114 reads have no
   record at their true copy.
3. **Most of the weight is reads O2 never touches.** 99.0% of truth rows come from `truth_untied`. Only 3,978 truth reads
   (13.3%) are in the corrector's domain: 573 fixed, 7 broken. THOC3 gives 65% of the P → LOOP-O2 gain. Three of S2's
   "better than P" families had no truth read moved.
4. **LOOP-O2 vs ORACLE.** The pooled gap is 0.0047 (not 0.004), and the two are equal in 4 families (not 5). At read level,
   LOOP-O2 is 0.031 against ORACLE's 0.000. Assembly hides the difference: no FCGR2A.2 transcript is built even with true
   homes.
5. **"Donor CNV, O3's case" is not supported.** FCGR2A and TBC1D3B have 0 `no_counterpart` reads. Their truth is chain-lifted
   `unique_elsewhere` (395 / 414 and 289 / 297 reads) and may itself be mislabelled.
6. **The echo removal is single-homing, not O2.** Primaries only does the same: DUSP22 0.017, ZC3H11A 0.000, NOMO1 0.142.
7. **Confirmed:** our assembler builds no transcript over FCGR2A.1 or .2. The overlap is at most 189 exonic bp, against 460
   primaries there and StringTie `cov` 394.9. Attribution variants (same strand only, exon overlap required, own family only)
   move the pooled scores by at most 0.007.

**What a valid test needs:**
- a read → transcript map from our assembler, scored like FLAIR's;
- ORACLE passing the bars as a pre-registered gate;
- a metric over the reads the corrector can move, the domain, rather than the 87% minimap2 already placed.


---

## PREREG_o2_second_round_donor_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_second_round_donor_2026-10-20.md`)*

## Pre-registration: the second round on a held-out donor (HG00733, simulated from its own assembly) (2026-10-20)

**User (2026-10-20):**
- After the audit of register 1345: "lets try to implement all of that, I think that is paramount to prove that O2 works".
- For the held-out data: "Donor simulation first (Recommended)".

### What is new against register 1345

1. **One read-crediting rule for every arm** (`bench/o3_graph/o2second.py score-reads`, 7c72b036). Each read contributes the
   records its arm's read pool admits:
   - primaries only;
   - `good`: primaries + secondaries with AS >= 0.98 x the read's best, the shipped pool;
   - `home`: a homed read's records at home, otherwise `good`.

   A record is credited to the copy of the arm's transcript covering at least half of its aligned bases. No tool's own read
   counting is used.
2. **Only the reads O2 can move are scored:** the corrector's domain reads (tied within 10 AS on two or more copies of one
   family), each scored when its true copy is among its tied copies. Each read is:
   - **right:** credited to its true copy alone;
   - **wrong:** credited to another copy only;
   - **echo:** credited to its true copy and another;
   - **uncredited.**
3. **A feasibility gate.** ORACLE, every read homed at its true copy, must beat the others before O2 is judged. Here the truth
   is exact for every read, so ORACLE is real perfect placement.

### The held-out donor

**HG00733** (HPRC year 1, phased hifiasm assemblies, maternal and paternal). It has never been used in this project.

Its reads were simulated from its own sequence (`bench/o3_graph/donorsim.py`):
- **Orthologs.** Each CHM13 copy of the 281 qualifying non-PAR families (`famident/qualify_nonpar.json`; every member pair at
  [0.90, 1.00) identity), with ±20 kb of flank, was aligned to each haplotype (minimap2 -x asm20 -c). The ortholog is the
  longest primary with MAPQ >= 20 through which every isoform lifts. Copies whose orthologs overlap on a haplotype were
  dropped (122 pairs).
  - Result: 1,252 copy-haplotype orthologs; 663 of the 724 copies have one or two. Median identity 0.998 over the whole
    flanked region.
- **Reads.** Each expressed copy's CHM13 isoforms were lifted onto its orthologs. 150 reads per copy, times realworld's
  expression multipliers, were split over haplotypes and isoforms (main 2.0, others 0.5). Each read is the donor
  transcript's own sequence, carrying the `A119b_testis_v4_ref` profile's:
  - retained introns;
  - truncation;
  - end jitter;
  - HiFi errors.

  The donor's alleles come from its genome, not from the profile.
  - Result: 83,700 reads (41,437 maternal, 42,263 paternal); 61 copies have no ortholog and 8 are silent.
- **Truth.** A read's truth is the CHM13 copy whose ortholog produced it.
- **Alignment.** The reads were aligned to CHM13 v2.0 as HG002's were: `-ax splice:hq -uf --secondary=yes -N 100 -p 0.5`,
  giving 271,686 records.

**Seen before this file:**
- the ortholog and build counts above;
- the alignment's flagstat.

No corrector call, assembly or score on these reads has been looked at.

The rules were set on HG002 (development): the corrector (version 2, bfa15f12), the home table (b15b5948) and the scorer
(7c72b036, development 0.805 vs 0.681 on 15 HG002 families).

### Arms (Rustle arms genome-wide on the simulated BAM)

| arm | what | read pool for the score |
|---|---|---|
| **LOOP-O2** | pass 2, every read homed (`loop_home.py o2 --every-read`): the corrector's copy when it decided, else minimap2's when not tied within 10 AS; ties O2 left echo | home |
| **P** | primaries only (`--seed-pool primary`) | primary |
| **GOOD** | pass 1, the shipped assembly | good |
| **ST** | StringTie 3.0.3 `-L` on the primaries | primary |
| **FL** | FLAIR 3.0.0 `bam2bed` + `collapse --trust_ends --generate_map` on the primaries | primary |
| ORACLE | pass 2, every read homed at its true copy (`--every-read --truth`) | home |

The corrector runs on the simulated BAM with its own confident reads as the sample (`o2correct.py donor`).

### Bars

- **G (gate):** ORACLE's right rate is at least 0.05 above the best of P, GOOD, ST and FL.
- **B1:** LOOP-O2's right rate is above P's, GOOD's, ST's and FL's.
- **B2:** LOOP-O2's wrong + echo rate is no higher than P's.
- **B3 (power):** at least 100 scored reads.
- **B4 (the user's bar, "more than 10 families"):** LOOP-O2's right rate is above P's in more than 10 families with at least
  10 scored reads.

**Verdict:**
- **UNDERPOWERED:** B3 fails.
- **INFEASIBLE:** G fails. No verdict on O2.
- **PASS:** B1, B2 and B4 hold.
- **FAIL:** otherwise.

**A PASS would claim:** on a donor never used in development, with real donor-vs-CHM13 copy differences and exact truth,
rebuilding loci after O2 gives the reads minimap2 leaves tied to their true copy's locus. That happens more often than with
minimap2's primaries, the shipped assembly, StringTie or FLAIR, and in more than 10 near-identical families.

**It would not claim:**
- real reads: expression and error profile are simulated, though calibrated on A119b;
- copies the donor has and CHM13 lacks: no reads are simulated from donor loci without a CHM13 ortholog.

**Rules:**
- Run script: `/mnt/linuxdisk/tmp/o3_graph/donorsim/HG00733/run.sh`.
- A bug is fixed and recorded as an amendment before any rerun.

### Results (`donorsim/HG00733/run.log`, `sim/second/score_reads.json`)

**FAIL. The gate holds, so the test was feasible (ORACLE 0.970). B1 and B2 fail; B3 and B4 hold.**

| HG00733, 11,596 scored reads (281 families) | right | wrong | echo | uncredited |
|---|---|---|---|---|
| **LOOP-O2** | 0.477 | 0.079 | **0.428** | 0.016 |
| P (primaries only) | 0.619 | 0.350 | 0 | 0.031 |
| GOOD (pass 1) | 0.077 | 0.023 | 0.889 | 0.012 |
| ST (StringTie) | 0.619 | 0.346 | 0 | 0.036 |
| FL (FLAIR) | 0.630 | 0.349 | 0 | 0.021 |
| ORACLE | 0.970 | 0.002 | 0 | 0.028 |

**Bars:**
- **B1 fails:** 0.477 < P's 0.619.
- **B2 fails:** wrong + echo is 0.507, against P's 0.350.
- **B4 holds:** LOOP-O2 credits more reads right than P in 20 families with at least 10 reads, but fewer in 71.

**Why it fails: the pre-registered home rule.** The corrector decided 6,440 of the 11,970 domain reads. The other 5,530
(4,937 `tie`, 593 `undominated`) have no home, so they keep pass 1's echo at every tied copy, and an echo is not "right".
Primaries only gives each of those reads minimap2's single pick, which is right about 62% of the time.

**On the reads it moved, the corrector is right.** Of the 455 reads O2 took off minimap2's copy:

| arm | right | wrong |
|---|---|---|
| LOOP-O2 | 406 | 47 |
| P, ST, FL | about 41 | 381-413 |

**Descriptive: the corrector itself on this held-out donor** (`o2correct.py score`, printed by the pre-registered run script;
not a bar):

| HG00733 | reads | corrector right | minimap2 right | fixes / breaks |
|---|---|---|---|---|
| all | 11,596 | 0.6695 | 0.6377 | 408 / 40 |
| exact ties | 4,975 | 0.488 | 0.449 | 225 / 29 |
| near ties | 6,621 | 0.806 | 0.780 | 183 / 11 |

52 families have more fixes than breaks; 7 have more breaks than fixes.

**Not claimed:** that loci rebuilt after O2 credit tied reads to their true copy more often than StringTie, FLAIR or
minimap2's primaries. Under this file's rule they do not, because echoed ties are not right.

**Post hoc, not a result:** with minimap2's primary as the home of the ties O2 leaves, LOOP-O2 would be P plus O2's moves:
about 0.619 + (406 - 41) / 11,596 = 0.650 right. That variant needs its own pre-registration on another donor.

### Audit (metric-auditor, 2026-10-20): holds with caveats

**What reproduces.** All six arms recompute to the digit with independent code: the home table rebuilt (77,226 rows),
independent crediting, the corrector line, and a hash-seed rerun. The FAIL stands. The descriptive corrector line is
genuinely held out: version 2 (13:44) predates the HG00733 download (14:47) and takes no truth input.

**Errata:**
1. **The undecided reads.** "Primaries only gives each of those reads minimap2's single pick, which is right about 62% of
   the time" is wrong. On the 5,381 scored reads the corrector left undecided, minimap2 is right 42.1% (P credits 40.7%);
   62% is P's rate over all reads.

   The mechanism holds. LOOP-O2 loses 2,022 reads on the undecided ones and gains 374 on the decided ones, exactly the
   1,648-read B1 gap. B1 counts P's coin toss on exact ties (0.433) as right, while an echo scores 0.

   LOOP-O2 also trails on near ties (0.748, against P's 0.760 and FL's 0.770). Its wrong rate is 0.079, against P's 0.350.
2. **The 455-moved-reads table is selected by O2's own decision.** Its unselected effect is the +365 reads (+3.1 points) of
   the corrector line, not extra evidence.
3. **The post-hoc 0.650 is arithmetically right** (0.651). Nothing was re-assembled, and its ceiling is P plus the corrector's
   read-level gain: the single-homing of register 1345's audit.
4. **Counts.**
   - 172 families have scored reads, not 281; the scorer printed the size of the family list.
   - "122 pairs" means 122 copy-haplotype entries.
   - The top 10 families hold 51% of the scored reads; CXorf51B alone holds 12.9%.

**The truth checks out**, with these exceptions:
- Strands: at least 99.9% of reads align on their gene's strand, and 98-99% of junctions are annotated.
- Copy identity: copy IDs match `pool_copies.tsv` for all 724 copies.
- Paralog-closer orthologs: in 29 cases a paralog's gene fits the called ortholog with fewer exon edits than the positional
  copy (587 scored reads). Excluding them leaves the verdict unchanged (LOOP-O2 0.499, P 0.637). These may be gene
  conversions or HG00733 assembly errors.
- SMN is simulated from SMN2 only: SMN1 is collapsed on the maternal haplotype and has MAPQ 14 on the paternal.
- TBC1D3 lacks 2 of its 11 copies.
- **Simulator bug:** `lift_exons` never checked exon length. 4 isoforms grew by 500 bp or more (PRAMEF11.5 paternal:
  1,021 → 43,526 bp), giving 284 reads, 26 of them scored. Fixed after the audit: an exon may change by at most 50 bp
  (`MAX_EXON_GROWTH`), with a test. This run is not regenerated; the bug does not change the verdict.

**The families the corrector worsens** (user, 2026-10-20: "are they real families or artifacts?"). All 7 are real annotated
families:

| family | breaks | fixes |
|---|---|---|
| AC114324.2 pair, DTX2 / DTX2P1, H4C11 / H4C12, KRT17P1 / P2 / P6 | 1 each | 0 |
| TRIM43B | 2 | — |
| PRAMEF7 | 7 | — |
| PRAMEF10 | 22 | 3 |

The 31 PRAMEF and TRIM43B breaks are all on the paternal haplotype, where HG00733's copies do not match CHM13's:
- **PRAMEF8 has no confident HG00733 ortholog.** 125 of the reads minimap2 puts on CHM13's PRAMEF8 come from the donor's
  PRAMEF7. They become PRAMEF8's sample, and the corrector moves 7 of PRAMEF7's tied reads there. This is the
  missing-copy case (register 1343): O3's.
- **The paternal PRAMEF cluster** is split over two contigs and inverted against the maternal one; its orthologs have MAPQ
  37-38.


---

## PREREG_o2_set_null_band_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_set_null_band_2026-10-20.md`)*

## Pre-registration: the Pareto set against a null arm, on the near-tie band alone, fresh one-genome simulations (2026-10-20)

**User (2026-10-20):** after register 1348's audit (a PSV-free null passes the same bars): "ok, lets try both, in order" —
first this test, then the FCGR2A lift check.

### Why

Register 1348 showed that the bars of `PREREG_o2_window_set` certify the 30-point window, not the set: a rule using no read
bases (minimap2's primary when no rival is within 10 points, every copy otherwise) passes them too. The set can show value
only where that null gives nothing: the **near-tie band**, reads whose closest rival is within 10 intron-free points. There
the null keeps every copy, so the question is whether the PSV columns narrow the pool, and at what cost in dropped truths.

### The rule (unchanged from `PREREG_o2_window_set`, code b4594ea2; `bench/o3_graph/o2window.py`)

- Window: alignments within 30 points of the best `ms` on two or more copies of one family, none off every copy, the
  primary among them; stitched alignments excluded.
- **Band under test:** reads first tied within 10 points (`window_copies` at W = 10, i.e. `ms > best − 10`).
- The set: `pareto_set` over the substitution columns of the common window.
- Single pick (fewest edits): descriptive.

### Arms, all scored on the band's reads

1. **The set.**
2. **Null, keep-all:** every copy in the 30-point window. It holds the truth whenever the truth is in the window and
   never narrows.
3. **Null, primary:** minimap2's primary alone (a one-copy set on every read).

### Substrate: round 26, fresh one-genome simulations

- `families26.json` (set `r26`, sha1 1826986addb5, 90 copies): 30 families drawn with `random.seed(2026)` from the 130 qualifying non-PAR families of the
  pool (`r19seq/families19cal_full2_all.json`) that no earlier round (families2-25) used. 2-9 copies each.
- Worlds INTACT, CROSS and SILENT under `profiles/{A119b,public}_testis_v4_ref.json` (no donor alleles), seeds 1601-1603:
  6 cells, 18 worlds, built by `fastset.py run` with round 25's `$V17` (`fast/heldout25/run_sim.sh`).
- Truth: each read's simulated copy. Every copy is in the reference.
- Scored read: determined, truth among the window's copies, and in the band. Pooled over the 18 worlds, reads keyed by
  world.

### Bars (the band only)

- **P:** at least 300 scored reads from at least 11 families.
- **H:** the set holds the truth on at least 99% of scored reads.
- **N (beats keep-all):** the set is one copy for at least 25% of scored reads. (Gorilla development: 42%; S1's band,
  seen after the fact: 38.6%. The keep-all null scores 0.)
- **Q (beats the primary):** of the set's singletons, at least 98% are right, and more are right than minimap2's primary
  on those same reads.

**Verdict:** PASS = P, H, N and Q hold; FAIL = P holds and another bar fails; UNDERPOWERED = P fails.

**A PASS would claim:** in one genome, among reads whose rival copy is within 10 intron-free points, the PSV columns
narrow at least a quarter of them to a single copy that is right at least 98% of the time and more often than minimap2's
coin toss, while keeping the true copy on at least 99% of all such reads. It would not claim anything about donors whose
copies differ from the reference (register 1348's S2), copies absent from the reference, or exact ties.

### Reported alongside, no bars

- the three arms on the whole 30-point window and on the 10-20 and 20-30 bands;
- the single pick's fixes and breaks;
- mean pool and set sizes;
- per family: scored, held, singletons, singletons right.

### Rules

- Round 26 is built and scored once with the code at this file's commit; bugs go into an amendment with the rerun.
- Heavy jobs foreground, one at a time, under 10 minutes each (one `fastset.py run` per command).
- The metric-auditor checks the verdict on a frozen copy before it is reported.
- Human numbers only; nothing pooled with gorilla.

### Amendments

#### Amendment 1 (2026-10-20): substrate G, gorilla KB3781's held-out families — written before any gorilla band was scored

User: "I thought that for this particular experiment we would use gorilla and not use the paternal and maternal
haplotypes till O3, can we check if it holds there too?"
- **Reads and reference:** his fibroblast library on GCF_029281585.2 (`personal/sub.bam`), the 143 held-out families of
  `PREREG_o2_personal_copies` (`qual.json`, half `heldout`). No arm uses his haplotypes.
- **Answer key only:** the intron-free key (`o2window.py truth-ms`: the copy with the best `ms` in his two assemblies, by
  at least 5 points, through `orthologs.tsv`), from alignments of the band's reads to `mGorGor1.{mat,pat}.splice.mmi`.
- **Same band, arms and bars** (P, H, N, Q). Verdict reported separately from round 26; gorilla and human never pooled.
- **Caveats fixed in advance:**
  - These families are not wholly fresh: their 10-AS core (896 reads) was scored once for the single pick in register
    1347, and that run's SM392 chimeras motivated the stitched rule. No set or window result was ever seen on them.
  - A read whose copy equals the reference cannot be decided by the key when the reference itself has the near-tie
    (both copies score alike on his identical haplotype); such reads will be undetermined and unscored. The scored band
    will therefore be weighted to reads carrying the other parent's copy, as register 1348's S2 was. Reported by origin.

### Results (code 7a7be665; round 26 built and scored once; `personal/window/r26/band.json`)

**The near-tie band (first tied within 10 points): PASS** — 8,955 scored reads from 22 families.

| arm | holds the truth | one copy | singletons right | minimap2's primary right on those reads |
|---|---|---|---|---|
| **the set** | 8,953 (99.98%) | 5,045 (56.3%) | 5,043 (99.96%) | 5,041 |
| null, keep-all | 8,955 (100%) | 0 | – | – |
| null, primary | 6,972 (77.9%) | 8,955 | 6,972 (77.9%) | – |

P, H, N and Q hold. **Q holds by two reads** (5,043 against 5,041): the set almost never corrects minimap2.

**What the set does, read by read:**
- On the 5,045 reads it narrows to one copy, minimap2's primary was already right on 5,041 (99.9%).
- On the 3,910 reads it keeps two or more copies, minimap2's primary is right on 1,931 (49.4%): a coin toss. Those are
  reads identical to their copies over the aligned bases (pool 2 → set 2: 1,341 reads; 3 → 3: 401; 4 → 4: 455).
- So in one genome the set **separates minimap2's safe picks from its coin tosses** rather than fixing the coin tosses:
  it says which tied reads are decidable from their bases (and decides them, 99.96% right) and which are not.
- Single pick: carries no bar; dropped truths 2 (both S19PRAMEF7).

**Other bands (no bars):** 10-20: set holds 6,098 / 6,098, one copy on 6,077, all right; minimap2's primary right on all
6,098. 20-30: 5,024 / 5,024, singletons 5,014 all right; primary right on all 5,024. Whole window: 20,077 scored, held
20,075, singletons 16,136 (right 16,134), primary right 18,094.

### Results, substrate G (gorilla KB3781, 143 held-out families; code 92f54122; `personal/window/g/band.json`)

Answer key: the intron-free key from his two assemblies (statuses over the 2,755 window-30 reads: refcons 1,514, altonly
1,212, undetermined 26, off-copy 3). O2 never saw the assemblies.

**The near-tie band (first tied within 10 points): UNDERPOWERED** — 175 scored reads from 12 families (P needs 300 / 11).
Descriptive: set holds 173, one copy for 84 (right 82), minimap2's primary right on those 77 (5 fixes, all SM392, 0 breaks);
on the 91 reads kept ambiguous minimap2 is right on 80.

**Other bands (no bars; the whole window passes P, H, N, Q descriptively):**

| band | scored | set holds | one copy (right) | minimap2 right on those | set fixes / breaks | kept > 1: minimap2 right |
|---|---|---|---|---|---|---|
| 10-20 | 1,336 | 1,335 | 802 (801) | 598 | 203 / 0 | 527 of 534 |
| 20-30 | 1,200 | 1,199 | 787 (786) | 770 | 16 / 0 | 412 of 413 |
| window 30, all | 2,711 | 2,707 | 1,673 (1,669) | 1,445 | 224 / 0 | 1,019 of 1,038 |

- **Every one of the 224 fixes is SM392**, the family whose minimap2 primaries are two-copy chimeras (register 1347's
  audit): the stitched primary is credited to SM392.0; O2 takes no stitched alignment as a placement and the PSV columns of
  the plain alignments decide SM392.1, right 224 / 224 by his paternal assembly.
- Outside SM392 the set decides only reads minimap2 already had right (1,445 of 1,449), and keeps ambiguous reads on which
  minimap2 is right 98% (unlike the one-genome simulations, where those were a coin toss): on gorilla the key cannot break
  a near-tie the reference itself has, so the band's refcons reads are mostly undetermined and unscored.
- Human and gorilla numbers are not pooled. Audit pending.

#### Audit of substrate G (metric-auditor, 2026-10-20): holds with caveats — the mechanism stated above is wrong

- Every number reproduces. UNDERPOWERED is the right verdict; the band also fails H (98.86%) and Q (97.6%), so with power it
  would have been FAIL.
- **No PSV column decided any of the 224 SM392 reads.** Each has exactly one non-stitched alignment in the window
  (`candidates` = SM392.1 alone), so the set is a singleton by elimination of the chimeric primary; `pareto_set` never ran.
  SM392.0 entered the window only through the chimera that was then excluded. This is an alignment-defect repair (the plain
  SM392.1 secondary already has the higher AS on 217 / 225), not a copy decision, and the stitched rule was designed on
  this family in this half: a named diagnostic, not a held-out finding.
- **Without SM392: fixes 0 / breaks 0** (1,445 = 1,445), the same picture as register 1348's null.
- **The narrowing is self-scored.** Singletons are 1,442 / 1,503 refcons reads (keyed by the reference's own alignment)
  and 7 / 983 altonly reads outside SM392; on those altonly reads minimap2 is right 98.0%.
- The SM392 truth itself is sound (paternal SM392.1 exon 3, margin >= 19; both copies present on the paternal haplotype in
  order).
- Wording not to use: "the PSV columns decide SM392.1"; "O2 fixes 224 reads"; "held-out replication" for SM392; "the set
  separates safe picks from coin tosses" on gorilla.


---

## PREREG_o2_smn_showcase_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_smn_showcase_2026-10-20.md`)*

## Pre-registration: the SMN1 / SMN2 showcase — PSV-decided reads, assign-or-abstain, and an assembler that honours them (2026-10-20)

**User (2026-10-20):** "O2 ... is to be able to offer a way to find the true provenance of multi-mapping long reads where the
flag as primary or secondary is just a coin toss and meaningless ... among the [families] there must be a family in which
we can tell based on PSVs, tags ... and use that to improve the assemblies that just ignore secondary alignments." Asked
which tried families' members the PSV columns tell apart, and offered SMN1/SMN2 as the showcase: "yes lets do that".

### The family and the reads

**SMN1 / SMN2** on CHM13 (pool family `S19SMN2`: `.0` = SMN2 chr5:70,808,816-70,837,687 −; `.1` = SMN1 chr5:71,381,817-
71,409,805 +), 99.8% identical. **HG002 Kinnex** reads (`hg002proof/chm13.bam`), the one real donor whose own assembly gives
every read's copy (`hg002truth`: classes unique / unique_elsewhere). SMN2 skips exon 7 in most transcripts; SMN1 keeps it.

**Seen before, stated plainly.** On these reads: register 1344 (the corrector's picks, StringTie/FLAIR, exon-7 readout) and
register 1348 (the set rule's read-level accuracy on SMN: 2,885 tied reads with a copy in the window, 2,258 decided, 2,249
right; 627 ambiguous, on which minimap2 is right 28). This is a demonstration with pre-registered bars on a matched
donor, not a held-out test. What has not been seen: what any assembler builds from the set rule's homes.

### O2's rule for this showcase (code at this file's commit)

For every read whose minimap2 primary lies on SMN1 or SMN2:
- **tied** when a rival copy's alignment is within 30 points of its best intron-free score (`ms`; `o2window.py`,
  window 30, stitched alignments excluded);
- a tied read is **decided** when, at every column where its two alignments differ by a substitution, the read agrees
  with one copy and not the other (`pareto_set` leaves one copy); it goes to that copy;
- a tied read that is not decided is **abstained**: identical over its bases on both copies, or mixed; it is given to
  neither copy;
- a read that is not tied keeps minimap2's primary.

No thresholds beyond the window; no use of the assembly.

### Arms (one record per read, `o2loci.py arms`)

| arm | each read's primary |
|---|---|
| **mm2** | minimap2's primary (what StringTie, FLAIR and a primaries-only assembler see) |
| **O2** | decided → its copy; abstained → **left out of the family**; untied → minimap2's |
| **O2-coin** | as O2, but abstained → minimap2's primary (reported; shows what the coin toss costs) |
| **truth** | HG002's assembly copy (unique / unique_elsewhere), else minimap2's |

**Assemblers on each arm's BAM:** Rustle (`tools/rustle_pipeline.sh assemble --seed-pool primary`, the thesis's own:
it builds loci from the records it is given), StringTie 3.0.3 `-L`, FLAIR 3.0.0 `collapse --trust_ends`.

### Readout (register 1344's, unchanged): exon-7 inclusion per copy

A read or transcript on a copy **spans** when its blocks overlap exon 6 and exon 8; it is FL when a block overlaps exon 7
by >= 20 bp, else Δ7. Inclusion = FL / (FL + Δ7): reads counted one each, StringTie transcripts by `cov`, FLAIR isoforms by
mapped reads, Rustle transcripts by `reads`. Exon coordinates: `o2loci_smn/smn_spec.json`.

**The reference values** are the truth arm's reads: SMN2 0.617, SMN1 0.985 (register 1344).

### Bars

- **R1 (reads):** on each copy, the O2 arm's read inclusion is closer to the truth arm's than mm2's is.
- **R2 (loci, O2 + Rustle against the tools that ignore secondaries):** on each copy, |Rustle on O2 − truth reads| is
  smaller than each of |StringTie on mm2 − truth reads|, |FLAIR on mm2 − truth reads| and |Rustle on mm2 − truth reads|.
- **R3 (the assembler honours the homes):** on each copy, |Rustle on O2 − Rustle on truth| <= 0.03.
- **R4 (abstaining is not free of cost but is honest):** per copy, the O2 arm's spanning-read count is at least 60% of
  the truth arm's (the abstained reads are the price; report their number).

**Verdict:** PASS = R1-R4 hold; FAIL otherwise. No power bar: SMN has thousands of reads.

**A PASS would claim:** on a 99.8%-identical family in a real donor, reading each tied read's PSV columns gives most reads
their copy and leaves the undecidable ones out, and an assembler that keeps those homes recovers each copy's own exon-7
pattern where StringTie, FLAIR and a primaries-only assembly mix them. It would not claim other families or donors.

### Rules

- Code written and tested before any arm is assembled; bugs found while running go into an amendment with the rerun.
- Foreground, one assembler at a time. Human numbers only.
- The metric-auditor checks the result on a frozen copy before it is reported.

### Amendments

(none yet)

### Results (code eca1c331; `o2smn/score.log`, `bars.json`; run once): **FAIL (R2 on SMN2, R4 on SMN2)**

3,651 reads with a primary on SMN1/SMN2; 3,623 tied within 30 points; the set decides 2,861 and abstains on 762. O2 moves
8 reads against minimap2; the truth moves 616.

| exon-7 inclusion | mm2 | O2 | O2-coin | truth |
|---|---|---|---|---|
| reads SMN2 | 0.653 | 0.635 | 0.646 | **0.617** |
| reads SMN1 | 0.887 | 0.988 | 0.890 | **0.985** |
| Rustle SMN2 | 0.602 | 0.590 | 0.602 | 0.585 |
| Rustle SMN1 | 0.801 | 0.969 | 0.801 | 0.941 |
| StringTie SMN2 | 0.654 | 0.634 | 0.649 | 0.618 |
| StringTie SMN1 | 0.894 | 1.000 | 0.896 | 0.989 |
| FLAIR SMN2 | 0.763 | 0.710 | 0.763 | 0.700 |
| FLAIR SMN1 | 0.897 | 0.929 | 0.897 | 0.934 |

- **R1 holds** on both copies (O2 0.018 / 0.003 from the truth reads; minimap2 0.036 / 0.098).
- **R2 holds on SMN1** (Rustle on O2 0.016 from truth; StringTie, FLAIR and Rustle on primaries 0.091 / 0.088 / 0.184)
  and **against StringTie and FLAIR on SMN2** (0.027 vs 0.037 / 0.146), but **fails against Rustle on primaries on SMN2**
  (0.015): minimap2's SMN2 reads are mostly right, the coin toss damages SMN1.
- **R3 holds** (0.005 / 0.028): Rustle keeps the homes it is given; FLAIR does not (O2-coin = mm2 for it).
- **R4 fails on SMN2:** abstaining leaves 723 of its 1,355 spanning reads (53%); SMN1 keeps 96%.
- O2-coin (abstained reads at the primary) is minimap2's picture again: the abstained reads are where minimap2 errs.

**Reading.** The PSV columns give SMN1 its own isoforms (exon 7 kept in 98.8% of reads, as in truth, against 88.7% under
minimap2), and an assembler that keeps the homes reproduces that; StringTie and FLAIR on primaries do not. The price is
abstention: nearly half of SMN2's reads are identical to SMN1 over their span and are left out. Audit pending.

### Audit (metric-auditor, 2026-10-20): holds with caveats

- Every number reproduces. The FAIL stands.
- **The SMN1 gain is removal, not assignment:** of the 616 reads the truth moves, O2 fixes 8 and abstains on 599; the
  other 2,241 decisions confirm a primary minimap2 already had right.
- **349 of the abstained FL reads carry SMN2's c.840 T on an SMN1 primary** while matching CHM13's SMN1 at other columns:
  HG002's SMN2 differs from CHM13's (the donor-vs-reference case), and the every-column rule abstains. "Identical over
  their span" was wrong for them.
- **R2 on SMN2 is infeasible:** Rustle on the truth arm is itself 0.033 from the truth reads; Rustle credits `reads` to
  16% of the family's reads, and R3 on SMN1 turns on one read (8 vs 7). The FLAIR "O2" column re-admitted the abstained
  reads from reads.fq and is not FLAIR on the O2 arm. R4's inputs pre-dated the prereg (r1348's 627 ambiguous).
- Wording not to use: "the PSV columns give most reads their copy"; "O2 recovers SMN1's isoforms"; "FLAIR on O2".

### Tags (descriptive, `o2tags.py`; `o2smn/tags.log`, `o2smn/tags/*.log`)

Per tied read, which statistic of its best alignment on each copy picks the donor's copy:

| statistic | SMN: decides / right / wrong (2,885 tied) | on the 627 abstained |
|---|---|---|
| AS, ms | 2,863 / 2,257 / 606 | 606 / 17 / 589 |
| NM | 2,509 / 2,252 / 257 | 262 / 14 / 248 |
| de | 2,511 / 2,256 / 255 | 263 / 16 / 247 |
| **substitutions only (NM − CIGAR indel bases = MD mismatches)** | **2,604 / 2,603 / 1** | **354 / 354 / 0** |
| indel bases | 2,683 / 1,476 / 1,207 | 603 / 11 / 592 |

Substitutions alone also on THOC3 5,482 / 5; DUSP22 1,927 / 1; GTF2H2C 1,408 / 3; NOMO1 348 / 0; GOLGA8A 257 / 10
(minimap2 256 / 282). Where the donor's copies differ from CHM13 every statistic fails: TBC1D3B 0 / 75 (AS 0 / 147), PDPR
179 / 98, GTF2IRD2 157 / 37.


---

## PREREG_o2_subs_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_subs_2026-10-20.md`)*

## Pre-registration: O2 by substitutions only, on fresh substrates (2026-10-20)

**User (2026-10-20):** "how about the analysis of CIGAR and bam tags such as md, nm, de, ms etc" — and, on the proposal to
test the rule it revealed on fresh data: "yes please".

### The rule

For a read tied within 30 intron-free points on two or more copies of one family (`o2window.py` window 30, stitched
alignments excluded, each copy by its highest-`ms` alignment): count each alignment's **substitutions** over the read bases
every tied alignment aligns (`alignment_events` "S" events in the common window: what the MD tag counts, less indels). The
copy with the **unique fewest** takes the read; a shared minimum keeps minimap2's primary (abstain). Indels, clips, AS,
ms, NM and de are not used.

**Why (seen on HG002, descriptive, `PREREG_o2_smn_showcase` "Tags"):** substitutions pick the true copy 2,603 / 1 on SMN,
5,482 / 5 on THOC3, and decide 354 / 0 of the reads the Pareto set abstained on; AS, ms, NM and de are misled by indels.
Register 1343's preamble said the same rule was a coin flip outside SMN/THOC3 on the 10-AS domain, with the edit answer
key that over-penalized clips; this test settles it.

### Substrates, never used for this rule

- **S1, one genome:** simulation round 27, `families27.json` (sha1 ebedbf33389e, 30 of the 100 qualifying non-PAR pool
  families no earlier round drew, seed 2027, 83 copies), the same build as round 26 (`heldout26/run_cell.sh` with
  families27, seeds 1701-1703, profiles A119b / public, INTACT / CROSS / SILENT). Truth = the simulated copy. Note: 3 drawn
  families (DUSP22, GTF2H2C, HSPA6) appeared in the HG002 tag table above; the rule was not tuned on them.
- **S2, real reads, another donor:** A119b testis fragments, seed 2 (`o2frag_a119b2/frag.bam`, 33,840 fragments of 11,280
  reads, `--eqx`); truth = the source read's copy (`truth.tsv`), as in register 1343. Known bias: a systematic full-read
  error counts against the rule. Qualifying families only (`famident/qualify_nonpar.json`).
- **S3, replication, seen families:** gorilla KB3781's 143 held-out families (`personal/window/g`), intron-free key from
  his assemblies (answer key only). Reported, no verdict.

### Scoring (`o2window.py pick-bars`, column `subs_pick`)

Scored read: determined truth among the window's copies. Per substrate: reads, families, the rule's decisions, right,
minimap2's primary right, fixes (rule right, primary wrong), breaks (the reverse), families fixing more than breaking and
the reverse; by band (first tied within 10 / 20 / 30).

### Bars (S1 and S2 each; as register 1343's)

- **P:** >= 300 scored reads from >= 11 families, minimap2 wrong on >= 50.
- **G:** accuracy >= minimap2's + 2.0 points.
- **S:** fixes >= 5 × breaks.
- **F:** >= 11 families fixing more than breaking, no more than a fifth as many (rounded up) the reverse.

**Verdict per substrate:** PASS = all four; FAIL = P holds and another fails; UNDERPOWERED = P fails. S1 and S2 are not
pooled; a claim needs both to PASS.

### Rules

Code committed with this file; each substrate prepared and scored once; bugs → amendment + rerun; foreground, one heavy
job at a time; audit on a frozen copy before reporting; human and gorilla never pooled.

### Amendments

(none yet)

### Results (code a45c91bc; each substrate scored once; frozen copy `personal/window/audit_frozen_subs_a45c91bc/`)

| substrate | scored / families | minimap2 right | substitution rule right | fixes / breaks | families better / worse | verdict |
|---|---|---|---|---|---|---|
| **S1** one genome, round 27 (18 worlds) | 23,145 / 27 | 21,125 | 21,125 | 0 / 0 | 0 / 0 | **FAIL** (G, F) |
| **S2** A119b fragments, seed 2 | 3,646 / 107 | 3,259 | 3,251 | 4 / 12 | 2 / 11 | **FAIL** (G, S, F) |
| S3 gorilla held-out (replication) | 2,711 / 29 | 2,464 | 2,439 | 1 / 26 | 1 / 6 | (no verdict; bars fail) |

**S1, read by read.** Of 23,145 scored reads the substitution count decides 18,713, and on **every one** it picks minimap2's
own copy (right on 18,710; both wrong on 3). The other 4,432 are substitution ties — 3,086 with equal counts on every copy —
and there minimap2 is right on 2,415 and wrong on 2,017: a coin toss. Every minimap2 error in S1 is a read whose bases
carry no substitution that tells its copies apart. Bands: first tied within 10, 11,501 scored, minimap2 right 9,483;
10-20, 7,394 / 7,392; 20-30, 4,250 / 4,250.

**S2.** Fixes 4, breaks 12, in 107 families (2 better, 11 worse). The truth is each fragment's source read's own call
(register 1343's known bias: a systematic full-read error counts against the rule).

**S3.** 1 fix, 26 breaks: 23 on reads carrying his other parent's alleles (SM233 10, SM502 8, SM556 5), 3 on reads whose
copy equals the reference. SM392's stitched reads abstain (one plain candidate, nothing to compare).

**Reading (pending the audit).** The HG002 SMN tag result (seen data) does not generalize. In one genome the substitution
count adds no correction — minimap2's primary already is the fewest-substitution copy whenever one exists — but it does
say exactly which primaries are coin tosses. In real donors the count is misled by the donor's own alleles, as every
allele-blind statistic is.

### Audit (metric-auditor, 2026-10-20): holds with caveats — the FAILs stand, the reading is narrowed

All numbers reproduce to the digit (own scorer; substitution counts re-derived from `E0.bam` X ops, 792 / 792). S1 FAIL
holds with caveats, S2 FAIL holds, S3 (no verdict) holds with caveats. Corrections to the Results above:

- **S1's FAIL was built in.** The profiles carry no donor alleles (`het_per_base = hom_per_base = 0`), so the true copy has
  both the best AS and the fewest substitutions wherever a substitution PSV is in the window; minimap2's unique errors are
  exact ties (837 of 839 equal AS; 834 with identical reference under the read). A read-sequence rule could fix at most
  ~7 reads: G and F were unreachable. "minimap2's primary is the fewest-substitution copy" describes this simulator.
- **S1 rows are read-world rows:** INTACT and SILENT reads are a subset of CROSS. 23,145 rows = 9,658 distinct reads
  (minimap2 right 8,815; decided 7,795, all equal to minimap2, 7,793 right; ties 1,863, minimap2 right 1,022).
- **Errors:** 2,020 (not 2,017): 2,017 on substitution ties, 3 rows (2 reads) decided with both picks wrong. "3,086
  with equal counts" is 3,020 equal plus 66 single-candidate IFNA21 reads; the coin-toss rate is 53.8% (2,349 / 4,366).
- **The decided / coin-toss split adds nothing over minimap2's own tags in one genome:** accepting the primary when its
  AS strictly beats every other window copy accepts 19,100 rows at 99.99%, the rest 50.1% right. This repeats register
  1342's certifier framing.
- **S2 is not fresh for this rule:** register 1342's audit measured a substitution count on it and register 1343 held it
  out. Its truth is AS-set but does not penalize the rule by construction: the full source read's substitution count
  favours the truth on 10 of the 12 breaks and 4 of 4 fixes; at most 2 breaks (ZNF705E, DHX40) may be truth artefacts.
  Families: 106, not 107.
- **S3:** 18 of the 26 breaks are confirmed donor alleles (his haplotype carries the read's base at every column where
  the read sides with the other copy), from 3 family-level haplotype differences. Families: 27, not 29. The code abstains
  when one non-stitched candidate remains (225 reads, SM392 224); had that copy taken the read, S3 would be 225 fixes / 26
  breaks. Without SM392, S3 is underpowered (22 minimap2 errors).
- **Register:** this re-tests register 1347's refuted fewest-substitutions arm (0 / 8 held-out, 0 / 16 development on the
  same gorilla families) and overlaps registers 707 and 1343; the prereg cited only 1343.
- **HG002 (seen data):** SMN **and THOC3** both show real gains (at least 334 and 236 reads where substitutions are right
  and minimap2 wrong), consistent with a donor indel misleading AS; DUSP22, GTF2H2C and NOMO1 show precision only. In
  S2 / S3 the donor differences at the deciding columns are substitutions, which mislead the count.
- Wording not to use: "every minimap2 error is a read with no distinguishing substitution"; "O2's honest output is a
  certificate"; "minimap2's primary is the fewest-substitution copy" without "in a no-allele simulation"; "SMN is the
  exception"; "S2 never used for this rule".

**Exploratory, after the verdicts (S1 development use, unique reads = CROSS worlds):** on the 1,834 two-or-more-candidate
coin-toss reads, giving each to the tied copy with more base-decided reads in the same world is right 1,095 / 1,824
(minimap2 988 on the same reads); at a >= 4:1 ratio 263 / 306 (86%) against minimap2's 161 (53%). It depends on the
simulator's expression model; real reads have no truth for exact ties.


---

## PREREG_o2_window_set_2026-10-20.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_window_set_2026-10-20.md`)*

## Pre-registration: O2 as a pool reducer in a 30-point intron-free window, on fresh substrates (2026-10-20)

**User (2026-10-20):**
- "the problem with nearly tied multimappers is that we cannot trust the primary, the primary and all the secondaries are a
  coin toss so that is why I am trying to either define a primary via PSVs and BAM tags or at least reduce the pool of
  possible confounded placements";
- "we could also go past the 10 AS difference from Eichler, if that is too close to the fully tied impossible cases";
- on the proposed test of "intron-free window 30 + the set rule" on fresh reads and simulations: "yes, please proceed".

### The rule under test (frozen at the commit that adds this file; `bench/o3_graph/o2window.py`)

1. **The window.** A read is in the domain when its alignments within 30 points of its best `ms` lie on two or more copies
   of one qualifying family, none off every copy, and minimap2's primary among them.
   - `ms` is minimap2's DP score without the intron charges. One substitution costs 5 points; "within 30" is `ms > best − 30`
     (as Eichler's 10 AS is a slack of 9).
   - Eichler's 10-AS domain is inside this one.
2. **Placements.** Each copy in the window is represented by its highest-`ms` alignment there. An alignment laying 50 or more
   bases on the exons of two copies (`stitched`) is no placement.
3. **The set** (`pareto_set`): the copies in the window that no other copy beats at every substitution column where the two
   differ, over the read bases every placement aligns. Copy r beats copy c when the read matches r and mismatches c at some
   column, and never the reverse. The set drops a copy only when the read contradicts it and never supports it against that
   rival.
4. **The single pick** (the tied copy with the fewest edits over the common window) is reported but carries no bar: on the
   gorilla development families it broke more than it fixed in every band.

### Development evidence (gorilla KB3781, development families; `docs/O2_WINDOW_SWEEP_2026-10-20.md`)

In the 30-point window, under the intron-free answer key: 9,482 scored reads; the set holds the truth on 9,476 (99.94%); it
is one copy for 8,340 (88.0%), right on 8,336 of them (99.95%); minimap2's primary is wrong on 529. The rule's constants
(30 points, 50 bp for a stitch, the Pareto comparison) were chosen there and are fixed here.

### Substrates: two, never seen by this rule

**S1, one genome (the user's regime for O2): simulations heldout24 + heldout25.** 36 worlds: samples A119b and public
× seeds 1401-1403 / 1501-1503 × INTACT, CROSS, SILENT; CHM13-only profiles (`*_v4_ref`), 21 + 14 families. Every copy is in
the reference. Truth = each read's simulated copy (`reads.truth.E0.tsv`). Copies: each world's `o2psv.copies.tsv`; genome:
its `genome.ref.fa`; reads: its `E0.bam`. These worlds were used before only for the corrector's S1 (register 1343) and the
regimen rounds, never for a window or a set rule.

**S2, real reads of another donor: HG002 Kinnex on CHM13** (`/mnt/linuxdisk/tmp/o3_graph/hg002proof/chm13.bam`, 633,939
candidate reads, `minimap2 -x splice:hq -uf --secondary=yes -N 100 -p 0.5`, no `--eqx`; match columns are read from the M
ops against the reference). Families: the 281 qualifying non-PAR families (`famident/qualify_nonpar.json`) of the pool
(`o2psv_dev/pool_copies.tsv`). Truth = the HG002 assembly's placement lifted to CHM13 (`hg002truth`: `truth.tsv` for the
10-AS-tied reads, `truth_untied.tsv` for the rest): class `unique` and `unique_elsewhere` → determined; `own_tie` →
undetermined; `no_counterpart`, `elsewhere` → off-copy; `poor_match` → left out.
- Caveat known in advance: HG002 is not the reference's donor, so every read may carry its own alleles (the O3 confound).
  That is the case where a set that does not guess should earn its keep.
- HG002 was development data for the corrector (register 1343) and the substrate of register 1329; neither used a window
  or a set.

Gorilla and human numbers are never pooled. The gorilla held-out half is not used: its SM392 chimeras motivated the
stitched rule, so it is not fresh for this test.

### Scoring (`o2window.py bars`, one window = 30)

A read is **scored** when its truth is determined and among the window's copies. Reported: scored reads, families, the
truth's statuses, reads whose truth lies outside the window, minimap2's primary right, the set holding the truth, singleton
sets and how many are right, minimap2's primary right on the singleton reads, the single pick's fixes and breaks, and the
same per band (reads first tied within 10, 20, 30 points).

### Bars, per substrate

- **P (power):** at least 300 scored reads from at least 11 families (S1: pooled over its 36 worlds, families distinct).
- **H (holds):** the set holds the truth on at least 99% of scored reads.
- **R (reduces):** the set is one copy for at least 50% of scored reads.
- **Q (quality):** of the singleton sets, at least 98% are right, and at least as many as minimap2's primary is right on
  those same reads.

**Verdict per substrate:** PASS = P, H, R and Q hold; FAIL = P holds and another bar fails; UNDERPOWERED = P fails. S1 is
the primary verdict (one genome); S2 is reported with its own verdict and never pooled with S1.

**A PASS on S1 would claim:** in one genome, O2 reduces the pool of a read tied within 30 intron-free points to the copies
its bases do not rule out; that pool keeps the true copy at least 99% of the time, is a single copy for at least half the
reads, and that single copy is right at least 98% of the time and at least as often as minimap2's primary. **A PASS on S2**
would extend that to a real donor whose copies differ from the reference's.

**It would not claim:** a single correct primary for every tied read (the set may keep several copies); copies absent from
the reference (off-copy: O3); exact ties (undetermined); anything about Eichler's 10-AS reads alone (they are a band here).

### Rules

- S1 and S2 are each prepared and scored once, with the code at this file's commit. Bugs found while running are logged in
  an amendment with the rerun; no constant changes.
- Heavy jobs in the foreground, one at a time, under 10 minutes per command.
- The metric-auditor checks both verdicts on a frozen copy before they are reported.
- Push only when asked.

### Amendments

#### Amendment 1 (2026-10-20): two pooling bugs, fixed before any bar was read

Both surfaced on the first S1 scoring and touch no rule, constant or bar:
1. **Simulated read names repeat across worlds** (`s0000001` in every world), so the pooled truth table held 3,779 reads
   instead of 60,345 and pooled calls collided. Truth and calls tables of pooled worlds now key each read by its world
   (`keyed`, `world_tag`, `calls_tag`); a single calls table (S2) is unkeyed.
2. **The `bars` parser** read `--window 30` as two truth files when it followed `--truth`; it is now parsed first.

The first S1 bars output (3,779-read truth) was discarded unread beyond the traceback; S2's bars had not run. The
`prepare` step and its outputs are untouched. Both substrates were then scored once with the fixed code.

### Results (code b4594ea2; `personal/window/{s1,s2}/bars.json`; each substrate scored once)

#### S1, one genome (36 simulated worlds, 29 families): **PASS**

| | scored | minimap2 right | set holds the truth | set of one copy | singletons right | minimap2 right on those reads |
|---|---|---|---|---|---|---|
| all, window 30 | 9,549 | 8,165 (85.5%) | 9,549 (100%) | 6,846 (71.7%) | 6,846 (100%) | 6,832 |
| first tied within 10 | 4,347 | 2,979 (68.5%) | 4,347 (100%) | 1,679 (38.6%) | 1,679 | 1,677 |
| first tied within 20 | 3,255 | 3,239 (99.5%) | 3,255 (100%) | 3,242 (99.6%) | 3,242 | 3,230 |
| first tied within 30 | 1,947 | 1,947 (100%) | 1,947 (100%) | 1,925 (98.9%) | 1,925 | 1,925 |

- P holds (9,549 reads, 29 families); H 100%; R 71.7%; Q 100% and 6,846 >= 6,832.
- The single pick: 12 fixes, 0 breaks.
- No read had its truth outside the window; no primary was stitched.
- In the band first tied within 10, the set is one copy for only 38.6% of reads: those are the near-identical copies, where
  the read often matches both copies over its whole length. The set keeps both rather than guessing, and minimap2 is wrong
  on 1,368 of those 4,347 reads.

#### S2, HG002 Kinnex on CHM13 (61 families): **FAIL (H, Q)**

| | scored | minimap2 right | set holds the truth | set of one copy | singletons right | minimap2 right on those reads |
|---|---|---|---|---|---|---|
| all, window 30 | 17,881 | 15,798 (88.4%) | 17,353 (97.05%) | 11,355 (63.5%) | 10,852 (95.6%) | 10,836 |
| first tied within 10 | 5,272 | 3,341 (63.4%) | 4,771 (90.5%) | 1,379 (26.2%) | 900 (65.3%) | 896 |
| first tied within 20 | 2,865 | 2,770 (96.7%) | 2,843 (99.2%) | 1,971 (68.8%) | 1,952 (99.0%) | 1,952 |
| first tied within 30 | 9,744 | 9,687 (99.4%) | 9,739 (99.95%) | 8,005 (82.2%) | 8,000 (99.94%) | 7,988 |

- P holds (17,881 reads, 61 families); **H fails** at 97.05%; R holds at 63.5%; **Q fails** at 95.6% (the set's singletons
  are still right at least as often as minimap2 on the same reads, 10,852 >= 10,836).
- The single pick: 18 fixes, 59 breaks.
- Domain reads left unscored: undetermined 1,851 (HG002's own copies tie), poor_match 6,023, off-copy 376, missing 36;
  truth outside the window 93.
- **Both failures sit in the band first tied within 10.** There the set drops the truth on 501 reads (9.5%) and its
  singletons are right on 65%: a read carrying HG002's own allele of its copy, where that allele equals the other copy's
  base, contradicts its true copy at that column and sides with the rival. That is the donor-vs-reference confound (O3's
  case) the pre-registration named in advance.
- Bands 10-20 and 20-30 pass every bar on HG002 too.

### Audit (metric-auditor, 2026-10-20): S1 holds with caveats, S2 holds; the reading does not

- **Every number reproduces to the digit** (independent scorer, `/mnt/linuxdisk/tmp/audit_2026-10-20_o2window/`); the
  frozen code is byte-identical to b4594ea2; Amendment 1 changes only the pooling keys and the parser.
- **The bars do not test the set.** A PSV-free null — minimap2's primary when no rival is within 10 points, every window
  copy otherwise — passes both substrates (recomputed by me as well): S1 held 99.83%, singletons 54.5%, right 99.69%;
  S2 held 99.15%, singletons 70.5%, right 98.79%. The set beats that null on S1 (100% / 71.7% / 100%) and is **worse than
  it on S2 on every bar**. So S1's PASS certifies the 30-point window, not the Pareto set; a null arm was missing.
- **R is carried by the window.** In the Eichler band (first tied within 10) the set is one copy for 38.6% (S1) and 26.2%
  (S2) of reads; the 10-30 bands supply R, and there minimap2 is already right 99.5-100%.
- **The set's gain over minimap2 on S1 is 14 reads** (2 in the <10 band, 12 in 10-20; ZNF658 and GPR89A), 0 the reverse;
  on the other 1,370 reads minimap2 gets wrong the set keeps more than one copy.
- **S2's H failure is mostly one family, FCGR2A:** 383 of 528 dropped truths and 379 of 479 wrong singletons, all one
  pattern (truth .2, minimap2 and set .1; 381 are `unique_elsewhere` with donor_nm 0-1 and donor margin 274-321). Whether
  HG002's locus carries .1 alleles (the O3 confound) or the chain lift lands on the wrong paralog of the FCGR2A/2C/2B
  region is not established. Without FCGR2A, S2's "all" would pass; its <10 band would still fail (held 97.5%, singletons
  right 90.0%). SMN2 accounts for 9 drops.
- **"Holds the truth" is conditioned on the window:** 93 S2 reads whose truth lies outside the 30-point copies are
  excluded; 87 of them are singletons (the primary), all wrong. Counting them, S2 H = 96.55%.
- **Unscored singletons:** the set is one copy on 774 of the 1,851 reads HG002's own genome calls undecidable and on 4,288
  of the 6,023 poor-match reads. Q is on determinable reads only.
- S2's 10-20 band passes by 6 reads (22 failures against 28 allowed); minimap2 alone passes the 20-30 band.
- Stitched primaries: S1 0, S2 59 (11 among the dropped truths). S1 reads are not error-free (16% NM 0, mean de 0.0023),
  so H = 100% is favoured, not forced.
- S1 and S2 are fresh for this rule (constants from gorilla); but r1329 had shown that 73% of HG002's ties are donor ties
  before the bars were chosen.
- **Wording not to use:** "never drops the truth"; "the set does not guess"; "O2 reduces the pool" as a claim distinct
  from "abstain when a rival is within 10 points"; "the 10-30 bands pass on HG002 too" without the null; "the failures
  are the donor's own alleles" without naming FCGR2A and its unverified lift.

**What a next test needs:** a pre-registered null arm (the rule above), scored on the Eichler band alone, where the set
must beat the null on held truth and singletons; and the FCGR2A lift settled by aligning CHM13's FCGR2A.1 / .2 to HG002's
haplotypes.

#### Amendment 2 (2026-10-20, post hoc, after the audit): the FCGR2A truth was a lift artifact — an erratum, not a verdict change

User: "ok, lets try both, in order" (second item: settle FCGR2A's lift).
- CHM13's FCGR2A.1 (chr1:160,843,885-160,862,793) and FCGR2A.2 (160,925,676-160,944,582), each ± 20 kb, were aligned to
  HG002's haplotypes (`minimap2 -x asm20 -c`, `s2/fcgr2a/orth.paf`). **Both land on the same 59-kb locus** of each
  haplotype (chr1_MATERNAL 156,468,027-156,527,020; chr1_PATERNAL 164,682,294-164,741,283): HG002 carries one copy where
  CHM13 carries two. That locus matches copy .1 far better (357 / 362 mismatches over 59 kb) than copy .2 (959 / 964).
- The 379 disputed reads (truth .2, set .1, class `unique_elsewhere`) all lie inside that locus (maternal 156,487,403-
  156,505,290), 337 of them with 0-1 mismatches. Their truth of ".2" came from the chain lift of a collapsed copy; by
  sequence they are HG002's .1-like copy, and the set's call was the better one.
- **Recount with every FCGR2A ".2" truth corrected to ".1" (401 reads), post hoc:** whole window H 99.18%, R 63.5%, Q
  98.91% (11,231 singletons right against the primary's 11,215) — the pre-registered bars would hold. **The near-tie band
  still fails:** held 97.67%, singletons right 92.68% (bars 99% / 98%). The remaining drops there are TBC1D3B, OR2T4 and
  SMN2 reads: the donor-vs-reference confound stands, smaller than reported.
- The registered S2 verdict (FAIL) is not changed by this amendment; it is recorded as an erratum of the answer key.


---

## PREREG_o2_corrector_idealsim_2026-10-21.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_corrector_idealsim_2026-10-21.md`)*

## Pre-registration: the frozen corrector on ideal-expression reads simulated from both of KB3781's haplotypes (2026-10-21)

**User (2026-10-21):** after register 1361 (his real fibroblast reads: UNDERPOWERED, minimap2 wrong on 3 of 107 tied reads),
asked whether reads simulated at optimal expression from every known gene would work better; on the proposal to simulate
from his non-reference haplotype: "ok lets do that".

**Question.** When every known copy is expressed from both of his haplotypes, does the frozen sample-aware corrector fix a
share of minimap2's wrong placements among tied reads, many times more than it breaks? The answer is a ceiling under ideal
expression, from a simulation, not a real-read result.

### The rule (frozen; nothing changes)

`bench/o3_graph/o2correct.py donor`, version 2, unchanged since bfa15f12, default flags, run as for register 1360
(`o2correct.py donor BAM COPIES.tsv GGO.fasta OUT.tsv --cache DIR --budget 500`). Domain, sample and decision are those of
`PREREG_o2_corrector_families_2026-10-20.md`. The corrector sees only the reads aligned to the GGO reference and the copies
table, never his haplotype assemblies.

### Families and copies

- **NPIP** (25 copies, catalog family `NPIP`), **LRPAP1** (22, `SM22`), **TBC1D3** (14, `TBC1D3`): 61 copies, each with its
  catalog exon chain (`rna_allele/fibro_cat.copies.tsv`, one model per copy).
- **SMN is dropped.** The census's "SMN (RefSeq, 3)" (Amendment 2 of PREREG_o2_corrector_families) matched the words "motor
  neuron": MNX1, SMNDC1 and SMN2. The gorilla annotation has one SMN gene (SMN2), in no catalog family, so SMN is not a
  multi-copy family here. This corrects that census line.
- **The YAGs stay out.** His Y has no second haplotype, so in his own genome every YAG tie is a true tie.

### Simulation (reads only)

- **Copy instances.**
  - **Reference instance:** each copy on the reference haplotype. Per chromosome the GGO reference is identical to one of his
    haplotype contigs (`rna_allele/chrmap.tsv`).
  - **Alternate instance:** the copy's counterpart on his other haplotype, by Amendment 4's positional key of
    PREREG_o2_corrector_families: ± 500 kb, `o2personal.py queries / orthologs --repaired --flank 500000`; the longest primary
    with MAPQ >= 20 on a non-reference contig, no exon-growth guard, overlaps to the lower divergence. It is recomputed for the
    61 copies; the 47 NPIP / LRPAP1 copies must reproduce `mainstay/positional/orthologs.tsv`.
- **Contested counterparts.** A counterpart that overlapped another copy's counterpart before the lower-divergence rule resolved
  it, winner and losers alike (for example maternal chr12 25.8-25.9 Mb, register 1361). Its reads are simulated, because his
  genome has the copy, but they are **never scored**. A copy that lost the overlap has a reference instance only.
- **Ideal expression.** Every instance gets the same **100 reads** over its one model: no expression spread, no silent copies
  (profile `A119b_testis_v4_ref` with `expr_log_sd` 0 and without the paralog-expression keys).
- **Reads.** `donorsim.py build`, the simulator of register 1346: the instance's own sequence (the reference, or his
  haplotype assembly through the lifted exons), with the profile's retained introns, truncation, end jitter and HiFi errors.
  - The read model is A119b testis Iso-Seq, the audited profile of register 1346 and the one-genome rounds. It is not fitted
    to his fibroblast library: `calibrate.py` is bound to CHM13 and its annotation.
  - His alleles, indels and conversions are real; only expression and sequencing are simulated.
- **Not simulated.** Copies his haplotypes carry but the reference lacks, such as the maternal-only NPIP locus of register
  1360. They are not known genes. Leaving them out favours the corrector: no foreign alleles teach its sample.
- **Two replicates:** seeds 3781 and 3782.
- **Alignment:** the real library's command, `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`, against the
  GGO reference (prebuilt `winloci_data/GGO.splice.mmi`), with minimap2 2.30 (the library used 2.31).

### Truth and scoring

- **Truth:** the reference copy whose instance produced the read. It is exact; no answer key is built from alignments.
- **Scoring:** `o2correct.py score ARM.tsv TRUTH.tsv --donor --among-candidates`. A read is scored when its truth is among its
  tied copies; reads from contested instances are left out.
- **Strata:**
  - **ref:** reads from reference instances, one genome. This is the no-harm control: minimap2's errors there should be exact
    sequence ties.
  - **alt:** reads from his other haplotype, the informative ones.
  - **all:** both together.

### Bars (each replicate, pooled over the three families)

- **P (power):** >= 50 scored reads and minimap2 wrong on >= 10; else UNDERPOWERED.
- **R (recovers a share):** fixes > breaks, and fixes − breaks >= 2% of minimap2's wrong calls.
- **S (safety):** fixes >= 5 × breaks.
- **H (no harm):** on ref reads, breaks <= fixes + 2.
- **F (families):** fixes > breaks in at least 2 of the 3 families, and breaks > fixes in none.

**Verdict:**
- **PASS:** P, R, S, H and F hold in both replicates.
- **FAIL:** a replicate where P holds and another bar fails.
- **UNDERPOWERED:** anything else.

### Expectations, written before any read is built

- **P:** expected to hold. The real library had 24 informative reads; ideal expression gives 100 reads to every alternate
  instance.
- **R, S, F:** no prediction. The corrector's earlier donor simulation (HG00733, register 1346) gave 408 fixes / 40 breaks;
  his haplotypes differ from the reference far less than a stranger's differ from CHM13.

### What a PASS would claim

With every known copy expressed, the sample's own unambiguous reads carry enough of his copy-specific variants for the frozen
corrector to recover a share of minimap2's wrong placements among tied reads, in more than one family, many times more than
it breaks. That is a ceiling for NPIP, LRPAP1 and TBC1D3 in KB3781.

**It would not claim:**
- real-read performance (register 1361: UNDERPOWERED);
- reads from copies the reference lacks (O3);
- contested loci;
- the YAGs;
- other read models or expression levels.

### Rules

- Each replicate is built and scored once; bugs → amendment + rerun.
- Heavy commands in the foreground, one at a time, each under 10 minutes.
- The metric-auditor checks the verdict on a frozen copy before it is reported.
- Gorilla only; never pooled with human numbers. Push only when asked.
- Code added for this run, behind new options with defaults unchanged:
  - `o2personal.py orthologs` writes the contested counterparts beside its output;
  - `donorsim.py build --per-instance` gives each haplotype instance its own read count;
  - `bench/o3_graph/o2idealsim.py` builds the copies, the simulation's instances, the truth strata and the bars.

### Results (code 9acb0233; run once per replicate; `/mnt/linuxdisk/tmp/o3_graph/idealsim/`): **FAIL**

- **Counterparts:** 53 of 61 copies (TBC1D3 13 of 14). The 47 NPIP / LRPAP1 copies reproduce Amendment 4's
  `orthologs.tsv` exactly. 14 LRPAP1 counterparts are contested: 7 simulated but unscored, 7 without an alternate instance.
- **Instances:** 61 reference, 46 alternate scored, 7 alternate contested; 100 reads each = 11,400 reads per replicate, all
  aligned (11,400 primaries, 0 unmapped).
- **Corrector:** 1,544 / 1,595 reads tied within 10 AS; its sample has 7,958 / 7,742 confident reads on 57 / 56 copies.

| replicate | scored | minimap2 wrong | fixes | breaks | families better / worse | P | R | S | H | F | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 3781 | 1,319 | 471 | 8 | 7 | 1 / 1 | ✓ | ✗ | ✗ | ✓ | ✗ | FAIL |
| 3782 | 1,365 | 503 | 10 | 4 | 3 / 0 | ✓ | ✗ | ✗ | ✓ | ✓ | FAIL |

The corrector changes ~1% of the tied calls and recovers a net 1 and 6 of minimap2's 471 and 503 errors (0.2%, 1.2%).

**Where minimap2's errors are (descriptive; 3781 / 3782):**
- **Reference instances, one genome:** 147 / 142 errors, all but one exact AS ties (identical copies over the read; the
  coin toss the set answers). Near ties: minimap2 right on 386 / 386 and 392 / 393.
- **His other haplotype, exact ties:** 188 / 192.
- **His other haplotype, near ties:** 136 / 169. Here the corrector's own evidence points to minimap2's pick: the true copy
  has more unexplained events than minimap2's copy on 116 / 153 reads, equal on 11 / 8, fewer on 9 / 8.
  - They concentrate in four instances (120 of 136, 129 of 169), each sending its reads to one neighbour:
    - TBC1D3.1 → .0: 65 / 73;
    - TBC1D3.7 → .6: 20 / 28;
    - TBC1D3.0 → .1: 16 / 16;
    - NPIP.1 → .2 (and .13): 19 / 12.
  - TBC1D3.0 and .1 are an inverted pair 105 kb apart. The reference segment spanning both (31.88-32.03 Mb) aligns reversed
    on his maternal chr4 (32.28-32.42 Mb), between collinear flanks: an inversion between his haplotypes with its ends in
    the two copies.
  - TBC1D3.6 and .7 lie inside a region whose whole ± 500 kb aligns reversed on his maternal chr4.
  - In these pairs his copy at a copy's position is closer in sequence to the neighbouring reference copy. Position and
    sequence disagree on provenance, and no read-level evidence follows position there.

Audit: below (FAIL holds with caveats; the TBC1D3 reading above is withdrawn).

### Audit (metric-auditor, 2026-10-21): FAIL holds with caveats; the TBC1D3 diagnosis above is withdrawn

Every number reproduces: scored, minimap2 wrong, fixes, breaks, families and bars, plus the diagnosis counts (116 / 11 / 9
and 153 / 8 / 8). Checked and found sound:
- **Order:** design 09:04, code 09:07, first output 09:07; the design text is unchanged.
- **Simulation:** 114 instances × 100 reads; the ideal profile changes only the expression keys.
- **Blindness:** no truth leaks through read names (a seeded permutation); the corrector reads only the BAM, the copies table
  and GGO.
- **Read sources:** 13 spot-checked reads match their own instance best.
- **Robustness:** no reasonable reading turns FAIL into PASS. Adding the contested reads gives 11 / 8 and 10 / 5. The
  31 / 36 reads whose truth lies outside the tied copies are minimap2 errors either way.

Corrections:
- **The power is positional only.** On just 3 / 471 and 9 / 503 of minimap2's errors does the true copy have strictly fewer
  edits (NM + clipped bases) than minimap2's pick. Under a fewest-edit key the test would be UNDERPOWERED. Of the
  corrector's scored overrules (19 / 16), one per replicate moves a read to the copy with fewer edits. **FAIL means the
  reads do not carry the positional provenance minimap2 misses; no read-based rule can recover it.**
- **TBC1D3: "closer to the neighbouring copy" is withdrawn.** His maternal copies are mosaics of the pair's paralog
  differences, equidistant from the two reference copies over the copy models.
  - Edits to its own copy vs the neighbour: maternal .1 4 / 4, .0 5 / 5, .7 2 / 2, .6 21 / 21.
  - His maternal .0 and .1 differ at one base (exon 14); reads that miss the 5′ exons look like .0.
  - These two copies give 86 / 471 and 97 / 503 of minimap2's errors.
  - The .0 / .1 inversion is supported: a 144.6-kb reversed piece between collinear pieces of 526 and 411 kb. Its ends are
    not pinned to the copies.
  - The .6 / .7 whole-region reversal keeps the positional correspondence and explains nothing.
- **"No read-level evidence follows position" overreaches.** On the alternate-haplotype near-tie errors, NM is higher on the
  true copy for 96 / 129 reads, equal for 39 / 36, lower for 1 / 4.
- **F fails in 3781 on LRPAP1.** All 3 SM22 breaks are reads of his maternal SM22.1 moved to SM22.8. That copy is closer to
  SM22.8 (25-mer containment 0.865 vs 0.792), and 72 / 74 of its reads sit confidently on SM22.8, where they enter the
  corrector's sample (register 1343's mechanism).
- **Small corrections:**
  - The corrector overrules 37 / 1,544 and 18 / 1,595 tied calls (2.4% / 1.1%), not "~1%".
  - In 3782 the fourth instance is TBC1D3.5 → .3 (13 reads), not NPIP.1 (12).
- **Outside O2's domain:** 526 / 518 untied reads from his other haplotype sit confidently on the wrong copy by the
  positional key. Lifted NPIP exons change by up to −5.7 kb, because `--repaired` has no growth guard.
- **Not checked:** whether the near-identity of his maternal TBC1D3.0 / .1 is an assembly collapse. Checking it needs his
  HiFi reads at the exon-14 base, a heavy run.

**What stands:** FAIL. With every known NPIP, LRPAP1 and TBC1D3 copy expressed from both of his haplotypes, minimap2 is
wrong by position on ~500 tied reads per replicate. On all but 3-9 of them, though, the read is at least as close to
minimap2's pick as to its own copy, and the frozen corrector recovers a net 1 and 6. In his own genome, tied reads carry
almost no sequence evidence for the copy minimap2 misses: ideal expression adds power, not information.


---

## O2O3_EXPLORE_2026-10-22.md

*(original: `docs/archive/2026-10/merged/O2O3_EXPLORE_2026-10-22.md`)*

## EXPLORATORY (post hoc, not pre-registered)

### Can a truth-free O3 gate make O2-VC abstain on bad copies? (2026-10-22)

**Status.** This is a post hoc analysis of the frozen outputs (code `07441e7c`, prereg `eccbc4e0`). It changes no
pre-registered verdict: E stays FAIL, A for NBPF1 stays FAIL, and A for GTF2IRD2 and PDPR stays PASS.

**Inputs.** `/mnt/linuxdisk/tmp/o3_graph/o2vc/audit_frozen_eccbc4e0/{hg002,hg00733}`, plus the truth tables
(`hg002proof/truth*.tsv`, `donorsim/HG00733/sim/truth.tsv`). Nothing was re-called and no code was changed.

**Scripts** (scratch, not in the repo): `/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/e542c6ba-28cf-42e7-becf-fd99d54f830f/scratchpad/ex/`.
- They reproduce the scorer exactly: HG002 pooled VC 271 / 4,579 (1,137 fixes, 24 breaks); HG00733 VC 176 / 2,801 (279
  fixes, 30 breaks).
- A re-implementation of the assign rule over `columns.tsv` reproduces `calls.tsv`'s `vc` on all 6,333 and 12,015 scored reads.

**Weighting.** HG00733 is primary: it is simulated, so its truth is exact. HG002 is reported with all scored classes and
without the chain-dependent `elsewhere` / `no_counterpart`, which the controller's audit found unreliable (NPIPB2#0's 16
reads are a chain gap over the donor's 57-bp insertion in NPIPB2.2; 29 of 133 `no_counterpart` sit in gaps aligning to a
CHM13 copy; GTF2H2C's 62 are undecided; NBPF1's 3 were probably right). Without them HG002's pooled VC is **85 / 4,393
(1.93%)**, not 271 / 4,579 (5.92%).

### 0. Bottom line

1. **On HG00733, no truth-free gate tested brings E1 under 2%.**
   - The three pre-specified copy gates remove 63-84 of the 176 wrong placements and leave E1 at 4.1-4.7%.
   - The cheapest is G3 (HomRef at AF ≥ 0.10): 64 wrong removed for 31 right lost (1 fix).
2. **The largest failure is invisible to every VCF-level signal.** PRAMEF7 makes 69 of the 176. PRAMEF7.1 has no donor
   counterpart (all 216 simulated reads come from .0), yet it has 28 teachers and only four PASS 1/1 records: no het,
   LowQual, HomRef, Dense or flag.
3. **On HG002, every gate removes most wrong placements (G1 94%, G2 93%, G3 98%; without the chain classes 81%, 81%, 94%),
   but each one lands its cost on a family VC gets right.**
   - G1 and G3 abstain on all of SMN2 (1,221 right, 612 fixes).
   - G2 abstains on all of PDPR (36 placements, 21 fixes) and DUSP22 (115 right).
   - GTF2IRD2 is untouched by every gate.
4. **The costs do not travel between donors.** G2 loses 323 right placements on HG002 and 759 on HG00733; G3 loses 1,526
   and 31. Two read-level patterns (§3b) flip sign between the donors.
5. **The best-separating signal is a `HomRef` record on the copy:** 2/8 bad vs 0/79 good copies on HG00733, 9/9 vs 6/21
   on HG002. It is specific on the simulated donor and sensitive on the real one, never both.
6. **An oracle bound (truth-side).** Marking every copy with no simulated donor read, plus G3, would give HG00733
   42 / 2,599 = 1.62%. This is the room that a real absent-copy detector could take (§4).

### 1. VC's wrong placements

**Key.**
- `placed←truth` gives copy suffixes.
- `kinds` are the columns voting for the placed copy (Q ≥ 10): `call` means ≥ 1 `separating_call` column; `ref` means only
  `separating_ref` columns.
- `erased` means the read also has an `erased_ref` column.

#### 1a. HG00733: 176 wrong, all of class `unique` (exact truth)

| family | wrong / placed | fixes | breaks | placed←truth (n) | vc_cols | kinds, erased |
|---|---|---|---|---|---|---|
| PRAMEF7 | 69 / 104 | 24 | 6 | .1←.0 (69) | 1: 23, ≥2: 46 | call 66 (58 with erased); ref+erased 3 |
| PRAMEF10 | 64 / 86 | 6 | 17 | .1←.4 (57), .3←.2 (7) | 1: 7, 2: 57 | ref+erased 64 |
| THOC3 | 9 / 17 | 4 | 4 | .1←.0 (9) | 1 | call+erased 9 |
| TBC1D3B | 7 / 18 | 1 | 0 | .4←.0 (5), .2←.1 (2) | 2-3 | call+ref (5 erased) |
| RHOXF2B | 6 / 188 | 6 | 0 | .1←.0 (6) | 1 | ref, no erased |
| NBPF6 | 6 / 21 | 1 | 0 | .1←.0 (6) | 1-2 | ref (2 erased) |
| TRIM43B | 5 / 9 | 1 | 2 | .1←.0 (5) | 1-3 | mixed |
| 7 others (VCX3A 3, GOLGA6L6 2, ...) | 10 / 197 | 14 | 1 | scattered | mostly 1 | mostly ref |

**Pooled HG00733 wrong placements.** 87 have a `separating_call` vote, 89 rest on reference PSVs only, 148 have an
`erased_ref` column, and 58 rest on a single column.

**Mechanisms, matched to the expected-failures signatures:**
- **PRAMEF7 is F3/F5: one donor copy teaches two reference copies, split by haplotype.** PRAMEF7.1 never appears in the
  truth, yet its 28 teachers give it hom-alt calls (e.g. chr1:12,725,212 T>C, 1-based; A>G in the read's frame; 1/1, DP 28)
  where .0 keeps {ref} from its 63 teachers. Every .0 read carrying that haplotype then votes `separating_call` for .1.
- **PRAMEF10 is F1/F2.** .1 carries six `LowQual` unphased 0/1 at AF 0.37-0.68 and a `HomRef` at AF 0.24; the truth copy .4
  has 240 simulated reads against .1's 74, consistent with .4's alleles in .1's teachers. The hets erase the reference PSVs
  that would have voted .4; the surviving reference PSVs vote .1.

#### 1b. HG002: 271 wrong; the chain-dependent classes are marked

| family | wrong / placed | `elsewhere` + `no_counterpart` | fixes | breaks | placed←truth | kinds, erased |
|---|---|---|---|---|---|---|
| GTF2H2C | 145 / 166 | 145 (83 + 62) | 2 | 0 | .1←(none) | ref+erased, ≥ 2 cols (143) |
| SPDYE1 | 59 / 77 | 24 (0 + 24) | 0 | 6 | .18←(none) 24; .20←.18 33; .20←.2/.19 2 | mostly call+erased |
| ANKRD36C | 21 / 64 | 0 | 9 | 14 | .1←.2 | ref, 1 col (20 of 21), 15 erased |
| AC0064533 | 11 / 21 | 0 | 0 | 0 | .0←.1 | ref, no erased |
| SLC35E2B | 9 / 47 | 0 | 5 | 4 | .1←.0 | ref+erased, 1 col (7 of 9) |
| NBPF1 | 9 / 14 | 9 (6 + 3) | 0 | 0 | .7←(none) 6; .3←(none) 3 | ref+erased / call |
| NPIPB2 | 5 / 6 | 4 (4 + 0) | 0 | 0 | .9←(none) 4; .18←.19 1 | ref+erased |
| CBWD2 | 3 / 46 | 3 | 0 | 0 | .2/.4←(none) | ref+erased |
| 6 others (LRRC37A2, H2AC18, FCGR2A 2 each; MSTRG3150, BMS1P1, AC1275021 1 each) | 9 | 1 | — | 0 | — | — |

**Pooled HG002 wrong placements.**
- 212 of the 271 rest on reference PSVs only, and 239 have an `erased_ref` column.
- Without the chain classes: 85 wrong, 60 of them with an `erased_ref` column.

### 2. Truth-free signals per copy: bad vs good

**Signals** (frozen VCFs, `teachers.tsv`, `flags.tsv`, `columns.tsv`): teachers; FILTER counts; `het_mid` / `het_low` =
0/1 at AF 0.2-0.8 / < 0.2; `het_paralog_mid` = a 0/1 at AF 0.2-0.8 whose ALT is another tied copy's reference base at the
homologous column; `alt_paralog` = the copy lacks its own reference base and holds a rival's; `swap` = complementary swaps
(F3's follow-up); `teacher_flag_reads`; `tied_flag_reads` = tied-read flags whose alleles touch the copy.

**Groups.**
- Bad copies: VC placed ≥ 3 wrong reads on them. Good copies: ≥ 5 right placements and 0 wrong.
- The family-level groups (bad = ≥ 3 VC wrong; good = ≥ 10 placed, 0 wrong) give the same picture. HG00733 has 8 bad and
  37 good families; HG002 has 8 bad and 9 good (4 and 12 without the chain classes).
- The fraction of tied reads with an `erased_ref` column is > 0 in 7/8 bad vs 25/37 good families on HG00733, and in
  7/8 vs 6/9 on HG002.

**Copies with signal > 0:**

| signal | HG00733 bad (8) | HG00733 good (79) | HG002 bad (9) | HG002 good (21) | HG002 no-chain bad (4) | good (23) |
|---|---|---|---|---|---|---|
| LowQual | 1 | 1 | 6 | 4 | 2 | 6 |
| **HomRef** | **2** | **0** | **9** | **6** | 4 | 8 |
| RnaEdit | 1 | 6 | 8 | 11 | 4 | 13 |
| Dense | 0 | 2 | 2 | 3 | 1 | 4 |
| het_mid | 5 | 39 | 9 | 16 | 4 | 18 |
| het_low | 0 | 1 | 7 | 9 | 3 | 11 |
| het_paralog_mid | 2 | 21 | 7 | 6 | 3 | 8 |
| alt_paralog | 5 | 30 | 4 | 9 | 1 | 11 |
| swap | 0 | 0 | 0 | 0 | 0 | 0 |
| teacher_flag_reads | 0 | 10 | 7 | 12 | 4 | 13 |
| tied_flag_reads | 0 | 2 | 2 | 7 | 1 | 8 |

**Reading.**
- **Dense, both flags and RnaEdit do not separate.** They mark good copies as often as bad ones, which fits expected-failures
  §0.5: Dense marks artifacts.
- **`het_paralog_mid` separates on HG002 only.** That is the "two donor copies on one reference copy" signature: 7/9 vs
  6/21 on HG002, but 2/8 vs 21/79 on HG00733.
- **The F3 swap test fires on no copy.** PRAMEF7 is not a swap: .0 stays at its reference base, and .1 goes hom-alt.
- **On HG00733 the bad copies are thinly taught:** PRAMEF7.1 has 28 teachers, THOC3.1 30, PRAMEF10.3 16 and NBPF6.1 11.
  156 of the 176 wrong placements land on the tied copy with the fewest teachers, against 551 of 2,625 right placements.
- **On HG002 the same pattern runs the other way.** 3,913 of 4,308 right placements land on the least-taught tied copy.

### 3. Gates, written down BEFORE their effect was computed (file first written 20:37 PDT, before `gates.py` ran)

Each gate marks copies from truth-free outputs of the frozen run. VC then abstains on every tied read whose tied set holds a
marked copy; nothing else changes. The gates were chosen after reading the step-2 signal tables and before any read-level
effect was computed.

- **G1 UNPHASED.** A copy is marked when its VCF has ≥ 1 `LowQual` record with GT 0/1: a het longcallR could not phase into
  the copy's two haplotypes. No new threshold.
- **G2 COLLAPSE.** A copy is marked when it has ≥ 1 0/1 record (any FILTER) with FORMAT AF in [0.2, 0.8] whose ALT is another
  tied copy's reference base at the homologous column (`columns.tsv`): two donor copies on one reference copy.
  - Threshold: the AF band 0.2-0.8, taken from the brief.
  - It was set by eye and not tuned.
- **G3 HOMREF.** A copy is marked when it has ≥ 1 `HomRef` record with FORMAT AF ≥ 0.10. The caller saw the allele at its own
  candidate floor and still wrote 0/0 (expected-failures F15's follow-up).
  - Threshold: 0.10, the caller's own constant.

#### 3a. Effect

"Right lost" counts right placements that become abstentions; "fixes lost" is the subset that minimap2 had wrong. Every
gate also removes breaks (HG00733: 30 → 16 / 9 / 16; HG002: 24 → 0).

| arm | copies marked | wrong removed | right lost | fixes lost | VC wrong / placed after | rate |
|---|---|---|---|---|---|---|
| **HG00733, ungated** | — | — | — | — | 176 / 2,801 | 6.28% |
| HG00733 G1 | 12 | 63 | 89 | 5 | 113 / 2,649 | 4.27% |
| HG00733 G2 | 65 | 84 | 759 | 146 | 92 / 1,958 | 4.70% |
| HG00733 G3 | 7 | 64 | 31 | 1 | 112 / 2,706 | **4.14%** |
| **HG002 pooled, ungated** | — | — | — | — | 271 / 4,579 | 5.92% |
| HG002 G1 | 40 | 254 | 1,403 | 632 | 17 / 2,922 | 0.58% |
| HG002 G2 | 35 | 251 | 323 | 61 | 20 / 4,005 | 0.50% |
| HG002 G3 | 43 | 266 | 1,526 | 651 | 5 / 2,787 | 0.18% |
| **HG002 no-chain, ungated** | — | — | — | — | 85 / 4,393 | 1.93% |
| HG002 no-chain G1 | 40 | 69 | 1,403 | 632 | 16 / 2,921 | 0.55% |
| HG002 no-chain G2 | 35 | 69 | 323 | 61 | 16 / 4,001 | 0.40% |
| HG002 no-chain G3 | 43 | 80 | 1,526 | 651 | 5 / 2,787 | 0.18% |

**HG00733, by family:**
- G1 removes PRAMEF10 57, VCX3A 3, GOLGA6L6 2 and FOXD4L1 1. It loses VCX3A 32, GOLGA6L6 28, PRAMEF10 12 and PRAMEF18 12.
- G2 removes PRAMEF10 57, THOC3 9, NBPF6 6, TRIM43B 5, ... It loses OR11H12 203, ZC3H11A 80, ANXA8L1 43, ...
- G3 removes PRAMEF10 57, NBPF6 6 and FOXD4L1 1. It loses NBPF6 15, PRAMEF10 12, OR2L8 3 and FOXD4L1 1.
- **No gate touches PRAMEF7's 69.**

**HG002 showcase families** (ungated: GTF2IRD2 0 / 44 with 24 fixes; PDPR 0 / 36 with 21 fixes):
- G1 and G3 leave both unchanged.
- G2 leaves GTF2IRD2 unchanged and takes **PDPR to 0 / 0, losing all 21 fixes.**
- NBPF1 (9 / 14 ungated) goes to 0 / 2 under G1, 3 / 6 under G2 and 0 / 0 under G3.

**Why G1 and G3 cost so much on HG002.** The cost is mostly SMN2 (a dev family). S19SMN2.0 carries the prereg's AF-0.14
`LowQual` 0/1 (chr5:70823918) and a `HomRef` at AF 0.25 (chr5:70811068). Its leakage is real but harmless: VC is 0 wrong on
its 1,221 placements.

#### 3b. Variants NOT pre-specified (hindsight; recorded so that nobody re-runs them as new ideas)

| variant | HG00733 after (wrong removed / right lost / fixes lost) | HG002 pooled after (same) | verdict |
|---|---|---|---|
| G3c: HomRef AF ≥ 0.10 makes the copy unknown at that column only (F15(f) verbatim) | 176 / 2,801 (0 / 0 / 0) | 282 / 4,602 (−11 / −12 / −4) | useless: no HG00733 read changes; on HG002 it only ADDS placements, since unknown-ing removes counter-votes |
| H4: mark copies whose records are all 1/1 (made to catch PRAMEF7.1) | 91 / 1,693 (85 / 1,023 / 168) | 268 / 1,966 (3 / 2,610 / 439) | far too costly |
| R1 (read-level): abstain on any read with an `erased_ref` column | 28 / 1,966 = 1.42% (148 / 687 / 171) | 32 / 252 (239 / 4,088 / 1,119) | flips sign between donors |
| R2 (read-level): abstain when VC picks the least-taught tied copy | 20 / 2,094 = 0.96% (156 / 551 / 91) | 67 / 462 (204 / 3,913 / 1,058) | flips sign; a 1/k-like expression prior |
| ORACLE (truth-side): mark the 50 copies with no simulated donor read (includes copies donorsim dropped as ambiguous) | 106 / 2,694 (70 / 37 / 24) | — | upper bound for an absent-copy detector |
| ORACLE + G3 | 42 / 2,599 = **1.62%** (134 / 68 / 25) | — | upper bound |

**What the read-level rules show.**
- R1 and R2 meet 2% on HG00733 but destroy HG002. HG00733's dominant failure is leakage onto the less-taught copy, which
  donorsim's even haplotype split and lognormal expression favour (F13). On real HG002, the right placements sit on exactly
  the reads those rules discard. Neither is an O3 signal; neither should be pre-registered.
- The oracle rows use the truth, so they are not a gate. They only bound what a correct absent-copy mark could add.

### 4. Recommendation

**1. Do not pre-register any VCF-level copy gate as an abstention rule.** On exact truth the best (G3) leaves E1 at 4.14%;
each gate's cost moves between donors; and the largest failure class leaves no VCF trace. RNA alone cannot tell one
heterozygous donor copy split over two reference copies from two homozygous copies (expected-failures Row 7, IsoCon), so no
teacher-allele test closes it either.

**2. Pre-register next: G3 plus a DNA-side absent-copy mark (G5), as one O3 gate.**
- G5 marks a reference copy whose copy-specific PSV k-mers are absent from the sample's own HiFi WGS: the planned O3 DNA
  step (`docs/PENDING_2026-09-23.md` §3), used here only to make O2 abstain.
- G3 covers the collapse/leakage class (PRAMEF10), G5 the absent-copy class (PRAMEF7). Freeze G3 exactly as in §3.
- The bound here for G3 plus a perfect G5 is 1.62% at 25 fixes lost.
- Pre-state G3's known real-read cost: on HG002 it loses 651 of 1,137 fixes (57%), 612 of them on SMN2.

**3. Fresh substrate.** A never-used HPRC year-1 donor (e.g. HG01978 or HG02257) with phased assemblies and public HiFi
WGS; only HG00733 exists under `donorsim/`. Simulate its RNA with `bench/o3_graph/donorsim.py`, adding beta-binomial
allelic ratios (F13's follow-up), because even splits flatter the leakage rules. Feed G5 the donor's real HiFi WGS. HG002 is
not fresh: it was used, and the caller was trained on it.

**4. What falsifies it on the fresh donor** (proposed bars, set by eye here). Any one of:
- gated E1 > 2%;
- > 10% of VC's fixes lost;
- no more wrong placements removed than right ones;
- G5 marks fewer than half of the copies that have no ortholog on either donor haplotype in donorsim's `orthologs` table.

The last bar uses that table rather than the oracle set: the oracle set (copies with no simulated read) also holds copies
donorsim dropped as ambiguous.

**5. Real-read check.** On HG002's verified-label reads (no chain classes), the same gate must keep GTF2IRD2's and PDPR's
fixes and lose < 10% of the pooled fixes. G3 alone fails that today.

### 5. RNA-only gates from F3 (appended 2026-10-22; the user rules out DNA reads, so §4.2's G5 is withdrawn)

**User (via the controller):** "assign reads to copies and find in those patterns the ones that could be explained by copies
not in the genome". O3 must be RNA-only, so the DNA-side mark proposed in §4.2 is dropped. The names G4 and G5 below are new
and replace it.

#### 5.1 Gates, written down BEFORE their effect was computed (appended 21:00 PDT, before `g45.py` ran at 21:00:34)

Both gates read only `columns.tsv`: per tied read, every column that separates two copies under VC or N3, with each tied
copy's reference base and called genotype in the read's frame.
- A column is counted once per homologous position pair, (A's chrom:pos, B's chrom:pos), whichever reads carry it.
- A copy whose genotype is `.` (unknown) at a column never counts there.
- The one constant is **≥ 2 columns**, set by eye.

- **G4 SWAP PAIR.**
  - **Mark.** A pair of tied copies (A, B) is a swap pair when ≥ 2 columns meet all of these:
    - VC separates A and B there (both known, disjoint genotypes);
    - A's genotype excludes A's reference base;
    - B's genotype excludes B's reference base.
  - **Merge.** Swap pairs are merged into groups, the connected components of the swap-pair graph. Each group is one copy in
    this sample.
  - **Effect.** A read VC places on a group member abstains when another member of that group is in its tied set: VC cannot
    choose between them. A read placed on a copy outside the group already beats every group member, so it is unchanged.
- **G5 TAUGHT-BY-PARALOG.**
  - **Mark.** Copy A is possibly absent when, for some rival B, ≥ 2 columns meet all of these:
    - the reference separates A and B (ref_A ≠ ref_B, both copies known);
    - A's genotype is exactly {ref_B}, so ref_A is absent from it.
  - **Effect.** A read VC places on a possibly-absent copy abstains, whatever else is tied.
- **Combinations.** G3+G4, G3+G5 and G3+G4+G5 abstain on the union of the reads each gate abstains on. G3 is the gate of
  §3, unchanged.

#### 5.2 Effect (computed after 5.1)

The "wrong removed / right lost (fixes lost) → after" columns follow §3a. G3 is repeated from §3a for reference.

| gate | marks (HG00733 / HG002) | HG00733 → E1 | HG002 pooled → rate | HG002 no-chain → rate | GTF2IRD2 / PDPR | PRAMEF7 / PRAMEF10 (HG00733) |
|---|---|---|---|---|---|---|
| ungated | — | 176 / 2,801 = 6.28% | 271 / 4,579 = 5.92% | 85 / 4,393 = 1.93% | 0/44 (24 fixes) / 0/36 (21) | 69/104 / 64/86 |
| G3 | 7 / 43 copies | 64 / 31 (1) → 4.14% | 266 / 1,526 (651) → 0.18% | 80 / 1,526 (651) → 0.18% | unchanged / unchanged | 69/104 / 7/17 |
| **G4** | **0 / 0 pairs** | 0 / 0 (0) → 6.28% | 0 / 0 (0) → 5.92% | 0 / 0 (0) → 1.93% | unchanged / unchanged | not caught / not caught |
| **G5** | 21 / 24 copies | 12 / 224 (65) → **6.39%** | 211 / 300 (89) → 1.47% | 38 / 300 (89) → 1.16% | **0/0, 24 fixes lost / 0/0, 21 fixes lost** | not caught / 57/63 |
| G3+G4 | = G3 | 4.14% | 0.18% | 0.18% | unchanged / unchanged | 69/104 / 7/17 |
| G3+G5 (= G3+G4+G5) | — | 76 / 249 (66) → 4.04% | 266 / 1,715 (730) → 0.19% | 80 / 1,715 (730) → 0.19% | lost / lost | 69/104 / 0/0 |

**G4 is inert.** On either donor, no tied pair has a single column where both copies' genotypes exclude their own
reference base, let alone 2. F3's pre-registered "complementary swap" signature does not occur in these data.

**G5 marks the source copy, not the sink.**
- PRAMEF10.4 equals .1's reference allele at 4 columns, so G5 marks .4. But .4 is the truth copy, with 240 simulated reads,
  of the 57 reads VC wrongly places on .1. So G5 abstains on .4's own right placements and leaves the 57 in place.
- It also marks PRAMEF10.3 (via .2, 2 columns), which removes the 7 .3←.2 errors.
- On HG002 it marks the showcase copies GTF2IRD2.1 and PDPR.0 and abstains on all of both families. Its good pooled rate
  there comes from GTF2H2C's 145 and SPDYE1's 59, both chain-dependent or unverified.
- On HG00733 it makes E1 worse than ungated.

**PRAMEF7.1 is caught by neither gate.** The (.0, .1) columns read:
- 4 sample-private hom-alt calls on .1 at reference-identical positions (VC separates; .0 keeps its reference base);
- 3 reference PSVs intact on both copies;
- 1 donor allele on .1 at a reference PSV;
- 1 column each way where one copy holds the other's reference base (G5 needs 2).

.0 keeps its own reference base at every column but one, so there is no swap. Read this way, the donor's second PRAMEF7.0
haplotype matches CHM13 .1 at the reference PSVs and adds four private SNVs. Each reference copy is taught by the haplotype
that resembles it, so in RNA the phantom looks exactly like a present, homozygous .1 carrying 4 private SNVs (IsoCon's
alleles-vs-copies problem; expected-failures Row 7).

#### 5.3 One more variant, NOT pre-specified (designed after reading PRAMEF7's columns above; dev-set numbers)

**H5 (single-haplotype pair).** Tied copies A and B are one sample copy when both conditions hold:
- ≥ 2 columns where ref_A = ref_B, VC separates them, and exactly one of them has left its reference base;
- neither copy's VCF has any 0/1 record.

VC then abstains on a read placed on A or B when the other is tied. This is the RNA-only reading of "a reference copy this
sample lacks": two homozygous haplotypes of one donor copy spread over two reference copies.

| arm | marks | HG00733 → E1 | HG002 pooled | HG002 no-chain | GTF2IRD2 / PDPR | PRAMEF7.1 |
|---|---|---|---|---|---|---|
| H5 | 13 / 2 pairs | 77 / 66 (31) → 3.72% | 3 / 108 (35) → 6.00% | 3 / 108 (35) → 1.91% | unchanged | **caught** |
| G3+H5 | — | 141 / 97 (32) → **1.37%** | 269 / 1,634 (686) → 0.07% | 83 / 1,634 (686) → 0.07% | unchanged | caught |

**Where H5 fires.**
- HG00733: PRAMEF7 (.0, .1); 7 TBC1D3B pairs; 2 PRAMEF11 pairs; 2 OR13C5 pairs; AL5913791.
- HG002: FCGR2A (.0, .1) and MSTRG3150 (.5, .6), costing 108 right placements for 3 wrong removed.

**Its cost.** Present homozygous copies with private SNVs look the same, so H5 trades right placements for wrong ones
wherever such copies are expressed. G3's SMN2 cost on HG002 (1,221 right, 612 fixes) is unchanged.

#### 5.4 Revised recommendation (RNA-only; replaces §4.2-4.4)

1. **Neither F3 gate is worth pre-registering.** G4 never fires, because the swap signature is absent from both donors. G5
   abstains on the source copy and destroys both HG002 showcases.
2. **The RNA-only candidate is G3+H5**, pre-registered as stated in §3 and §5.3 with "≥ 2" as its one constant. Its 1.37% on
   HG00733 was found on the data that inspired H5, so it is a dev number. Pre-state both known costs:
   - G3 on SMN2-like families (HG002: 1,526 right lost);
   - H5 on homozygous copies with private SNVs (HG00733: 66 right lost; HG002: 108).
3. **Substrate.** A never-used HPRC year-1 donor simulated with `donorsim.py`, with beta-binomial allelic ratios, so that allelic
   imbalance tests H5's "no het" condition. Truth: donorsim's per-copy ortholog table, used only for scoring.
4. **Falsification.** Any one of:
   - gated E1 > 2%;
   - > 10% of VC's fixes lost;
   - wrong removed ≤ right lost;
   - H5 marking fewer than half of the copy pairs where one reference copy has no ortholog in the donor while its tied rival
     has one. This is the PRAMEF7 class, and these thresholds were set by eye.

   On HG002's verified-label reads, G3+H5 already fails the fixes bar (686 of 1,137 lost). That cost is pre-stated, not
   discovered.


---

## O2VC_DEV_2026-10-22.md

*(original: `docs/archive/2026-10/merged/O2VC_DEV_2026-10-22.md`)*

## O2 variant calls: the development run on HG002 (S19SMN2, S19MSTRG3150), 2026-10-22

Task 6 of `docs/superpowers/plans/2026-10-22-o2-variant-calls.md`, under the controller rulings
(`.superpowers/sdd/2026-10-22-o2-variant-calls/task-6-rulings.md`). **These are the numbers seen before the pre-registration
(Task 8); the prereg quotes them as already seen.** Only the two development families were run and scored: no panel family,
no other donor, no world.

Code: branch `o3-fast`, commits 643bf8ae (the `--weak-het` policy, the exit code, the A2/C1 evaluability ruling,
`copies-from-sets`) and fdb41f25 (the weak-het default frozen at P1 `keep`).

### 1. Commands and timings

```bash
T=/mnt/linuxdisk/tmp/o3_graph; G=/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa
## P1 (keep): the one longcallR run
RLOCK_TIMEOUT=5400 bash tools/rlock.sh heavy python3 bench/o3_graph/o2vc.py run $T/hg002proof/chm13.bam \
    $T/o2psv_dev/pool_copies.tsv $G $T/o2vc/dev_P1 --families S19SMN2,S19MSTRG3150 --exclude-chrM
## P2 (unknown) and P3 (ref) reuse P1's VCFs (copied into dev_P2/vcf, dev_P3/vcf; byte-identical, 8 files)
... run ... $T/o2vc/dev_P2 --families S19SMN2,S19MSTRG3150 --exclude-chrM --weak-het unknown
... run ... $T/o2vc/dev_P3 --families S19SMN2,S19MSTRG3150 --exclude-chrM --weak-het ref
## each policy scored twice: no panel, then the panel A view of SMN2
python3 bench/o3_graph/o2vc.py score $T/o2vc/dev_PX $T/hg002proof/truth.tsv $T/hg002proof/truth_untied.tsv \
    --mode donor --n4 $T/o2psv_dev/hg002.v2.tsv [--panel A --focus S19SMN2]
```

| step | wall | peak RSS | exit |
|---|---|---|---|
| run P1 (scan of the 1.73 M-record BAM ~61 s, 8 longcallR calls) | 80.6 s | 0.84 GB | 0 |
| run P2 (8 VCFs reused, 0 calls attempted) | 69.9 s | 0.36 GB | 0 |
| run P3 (8 VCFs reused, 0 calls attempted) | 69.3 s | 0.36 GB | 0 |
| score, each of 6 | ~2 s | - | 0 |

Logs: `$T/o2vc/logs/` (`dev_P*.stderr`, `dev_P*.time`, `score_P*_{nopanel,panelA}.stdout`).

### 2. The run (identical for the three policies except the weak-het treatment)

- `run.json`: longcallR 2.0.0, preset hifi-masseq, min_mapq 0; families S19MSTRG3150 and S19SMN2; no chrM family; 0 failed
  copies; 0 Q0 teacher bases; 0 tied reads without QUAL.
- **depth_agreement: 41 of 41 records (100%)** have |DP - teacher_depth| <= 1, in all three runs (>= 95% required: no
  teacher_depth change needed).
- weak hets (a het whose GT alt has FORMAT AF < 0.15): **1**, `chr5:70823918 A>G LowQual 0/1 GQ 92 DP 79 AF 0.14` on S19SMN2.0
  (SMN2, minus strand), in exon 2a (`[70823848, 70823920)`), so T>C in the read frame. Against it S19SMN2.1 (SMN1) has C at
  71395632: a reference PSV (SMN2.0 T / SMN1 C), on **1,589 of the 1,655** SMN2 tied reads (860 carry C, 729 T).
- Tied reads: 1,849 = S19SMN2 1,655 + S19MSTRG3150 194, all kept (0 skipped). Teachers: 4,779 on 8 of 9 copies
  (S19MSTRG3150.2 none; .1 and .4 have 2 each, so they stay unknown).
- Hand checks (all three policies): every VC placement has `vc_cols >= 1` (0 exceptions); no flagged read is placed (4 flagged,
  0 placed); no tied read's primary lies outside its tied copies (`outside` 0).
- **The exon-7 c.840 site is in `S19SMN2.columns.tsv`**: S19SMN2.0 70810811 (T) vs S19SMN2.1 71408733 (C), the 6th base of the
  54-bp exon (`[70810763, 70810817)` minus / `[71408728, 71408782)` plus): `separating_ref` on 852 reads (830 T, 22 C). The
  calls agree with the reference there (SMN2.0 T, DP 72; SMN1 C, DP 1,860).

### 3. Truth classes of the tied reads (the same under every policy)

| family | tied | scored | unique | unique_elsewhere | own_tie | poor_match | no_truth_row |
|---|---|---|---|---|---|---|---|
| S19SMN2 | 1,655 | 1,299 | 693 | 606 | 113 | 243 | **0** |
| S19MSTRG3150 | 194 | 141 | 123 | 18 | 15 | 37 | **1** |
| pooled | 1,849 | 1,440 | 816 | 624 | 128 | 280 | **1** |

No elsewhere and no no_counterpart read. Decidable (truth copy among the tied copies): 1,439 of 1,440 scored.

### 4. Per policy (pooled over both dev families, all arms on the same 1,440 scored reads)

P1 `keep` = the record's GT alleles (the spec); P2 `unknown` = the copy unknown at a weak het; P3 `ref` = the weak alt dropped.

| policy | VC placed | right | wrong | wrong/placed | fixes | breaks | both wrong | abstained | K1 (own_tie placed) | K2 (flagged) | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|
| P1 keep | 1,330 | 1,329 | 1 | 1/1330 = 0.0008 | 647 | 0 | 1 | 110 | 0/128 = 0 | 3/1568 = 0.0019 | 2 (4 reads) |
| P2 unknown | 1,330 | 1,329 | 1 | 1/1330 = 0.0008 | 647 | 0 | 1 | 110 | 0/128 = 0 | 3/1568 = 0.0019 | 2 (4 reads) |
| P3 ref | 785 | 782 | 3 | 3/785 = 0.0038 | 43 | 0 | 3 | 655 | **98/128 = 0.766** | 3/1568 = 0.0019 | 2 (4 reads) |

VC placements on reads that are not scored: P1 181 and P2 181 (all poor_match); P3 271 (98 own_tie, 173 poor_match).

Erased-dependent VC placements (VC places, N3 abstains or places elsewhere, the read has a voting `erased_ref` column):

| policy | placements | on scored reads | right | wrong |
|---|---|---|---|---|
| P1 keep | 797 (SMN2 688, MSTRG3150 109) | 706 | 706 | 0 |
| P2 unknown | 109 (MSTRG3150 only) | 92 | 92 | 0 |
| P3 ref | 109 (MSTRG3150 only) | 92 | 92 | 0 |

Every arm on the same reads (placed / right / wrong / fixes / breaks; N1 always places, so it has no fixes or breaks of its own):

| arm | P1 keep | P2 unknown | P3 ref |
|---|---|---|---|
| N1 minimap2 | 1440 / 784 / 656 | same | same |
| N2 MAPQ >= 1 | 1358 / 735 / 623 / 0 / 0 | same | same |
| N3 reference only | 692 / 677 / 15 / 8 / 0 | **1237 / 1224 / 13 / 612 / 0** | 692 / 677 / 15 / 8 / 0 |
| N4 corrector (holds 1,849 of 1,849) | 1440 / 1421 / 19 / 638 / 1 | same | same |
| VC | 1330 / 1329 / 1 / 647 / 0 | 1330 / 1329 / 1 / 647 / 0 | 785 / 782 / 3 / 43 / 0 |

Fix evidence (VC): P1 30 sample call, 609 reference PSV corrected by the sample, 8 reference PSV; P2 30 / 5 / 612; P3 30 / 5 / 8.
Read by read: P1 and P2 give the same `vc` and `flag` on all 1,849 tied reads (their `n3` differs on 922); P1 and P3 differ on
921 reads' `vc` (their `n3` is identical).

#### Per family under P1 (the chosen policy)

| family | arm | placed | right | wrong | abstained | fixes | breaks |
|---|---|---|---|---|---|---|---|
| S19SMN2 | N1 | 1299 | 683 | 616 | 0 | - | - |
| | N2 | 1277 | 671 | 606 | 22 | 0 | 0 |
| | N3 | 677 | 675 | 2 | 622 | 8 | 0 |
| | N4 | 1299 | 1289 | 10 | 0 | 607 | 1 |
| | VC | 1221 | 1221 | 0 | 78 | 612 | 0 |
| S19MSTRG3150 | N1 | 141 | 101 | 40 | 0 | - | - |
| | N2 | 81 | 64 | 17 | 60 | 0 | 0 |
| | N3 | 15 | 2 | 13 | 126 | 0 | 0 |
| | N4 | 141 | 132 | 9 | 0 | 31 | 0 |
| | VC | 109 | 108 | 1 | 32 | 35 | 0 |

K1 per family: SMN2 0/113, MSTRG3150 0/15. K2: SMN2 1/1412, MSTRG3150 2/156. Flags: S19SMN2#0 (2 reads: 1 unique,
1 poor_match), S19MSTRG3150#0 (2 unique reads). C1: no no_counterpart read (not evaluable).

The one wrong VC placement (all policies): MSTRG3150 read `m84036_230422_223801_s1/41616922/ccs/69_1474`, class
unique_elsewhere (truth S19MSTRG3150.4, outside its tied copies .5/.6; .4 has 2 teachers, so it is unknown), placed on .5 by one
sample-call column; minimap2 also wrong (.5).

#### The panel A bars on SMN2 under each policy (a dev view, not a panel run)

| policy | A1 VC wrong/placed | A2 fixes >= 10 | A3 VC fixes > N3 fixes | K1 | K2 | verdict |
|---|---|---|---|---|---|---|
| P1 keep | 0/1221 PASS | 612 PASS | 612 > 8 PASS | 0/128 PASS | 3/1568 PASS | PASS |
| P2 unknown | 0/1221 PASS | 612 PASS | **612 > 612 FAIL** | 0/128 PASS | 3/1568 PASS | FAIL |
| P3 ref | 2/676 PASS | **8 FAIL** (616 minimap2-wrong: evaluable) | **8 > 8 FAIL** | **98/128 FAIL** | 3/1568 PASS | FAIL |

### 5. The choice (frozen rule, ruling B8)

1. Exclude a policy with dev K1 > 2%: **P3 excluded** (98/128 = 76.6%). P1 and P2: 0/128.
2. Lowest pooled VC wrong/placed: P1 1/1330 = P2 1/1330, a tie.
3. Ties broken by more fixes: P1 647 = P2 647, still a tie (their VC placements are the same read by read).

The rule does not separate P1 from P2. **The default stays P1 `keep`**, the spec's genotype rule ("the record's GT alleles,
whatever its FILTER"), since the rule gives no reason to change it (commit fdb41f25, pinned by
`test_the_weak_het_default_is_the_policy_the_frozen_dev_rule_chose`).

What P1 vs P2 changes is N3, and with it A3. N3 knows only what VC knows. Under P1, VC knows SMN2.0 at the exon-2a column
({C, T}), so N3 votes there with the reference PSV and abstains on the reads it contradicts. Under P2, both are blind there,
and N3 equals VC on SMN2 (612 fixes each).

**Finding: on SMN2 the whole of VC's gain over the reference alone is the erasure of one reference PSV by one LowQual record
(AF 0.14, DP 79).** 604 of SMN2's 612 VC fixes rest on it ("reference PSV corrected by the sample"). Its 688 erased-dependent
placements (614 scored) are all right. Restoring the PSV (P3) places 98 of 113 own_tie reads and drops the fixes to 8.

Interpretation (not tested here): a copy's teachers are chosen by AS margin, which favours reads that carry the copy's own
reference base at a PSV. A teacher allele fraction at a reference PSV is therefore biased low. SMN2.0's teachers show 14% C,
but the tied reads show 860 C and 729 T.

### 6. HG002 BAM provenance

- `hg002proof/chm13.bam`: **no @RG**. @PG: `minimap2 2.30-r1287 -a -x splice:hq -uf --secondary=yes -N 100 -p 0.5 -t 4
  .../idx/target.splice.mmi candidates.fq`, then samtools 1.22.1 sort. The reads were realigned from FASTQ, so the read group
  was lost. Every read name, checked over the SMN2 span, comes from movie **m84036_230422_223801_s1**.
- The source FLNC, `/mnt/linuxdisk/home/juanfraitu/hg002_kinnex/flnc.bam` (63.2 GB, 37.2 M FLNC, size-verified download from
  `https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/` on 2026-10-13):
  - @RG: `PU:m84036_230422_223801_s1 PM:REVIO PL:PACBIO LB:HG02_8mer_1.05X SM:BioSample_1
    DS:READTYPE=SEGMENT;SOURCE=CCS;...;BASECALLERVERSION=5.0`, barcodes Iso-Seq_v2_Barcoded_cDNA_Primers.
  - @PG chain: ccs 7.0.0, primrose 1.4.0, lima 2.8.99, **skera split 1.1.99** (Kinnex/MAS-seq segmentation), lima 2.9.0
    `--isoseq`, refine 4.1.0 `--require-polya`, under SMRT Link 13.0 at PacBio's `/pbi/dept/...` (2023 "honeybee" project).
- **It is PacBio's public HG002 Kinnex full-length RNA release:** the same movie as every `chm13.bam` read, Revio, skera
  segmentation. The header carries no "HG002a:MAS:cDNA" label, so I cannot confirm that it is the very dataset longcallR's MAS
  model was trained on ("the Kinnex full-length RNA dataset (HG002a:MAS:cDNA) from PacBio"). It is very likely the same
  sample, and possibly the same reads. Any HG002 result may then be optimistic for the caller, a train/test overlap.

### 7. Bugs found

None in the dev run. The code items (weak-het policy, exit code, A2/C1 evaluability, copies-from-sets) were built TDD-first.
Three pre-ruling test expectations were updated to the A2 ruling: A2 at 9 and at 3 minimap2-wrong reads, and the
"VC places nothing" case.

### 8. Also done: the DAZ copy table

`/mnt/linuxdisk/tmp/o3_graph/o2vc/inputs/daz_copies.tsv` was written by `o2vc.py copies-from-sets
.../r19seq/families19cal_full2_all.json seq19 S19DAZ1`. It has **4 copies** (DAZ1-4, chrY), byte-identical to the S19DAZ1 rows
of `o2psv_dev/pool_copies.tsv`. Nothing was run on it.


---

## O2_NAIVE_TIED_2026-10-22.md

*(original: `docs/archive/2026-10/merged/O2_NAIVE_TIED_2026-10-22.md`)*

EXPLORATORY (post hoc, not pre-registered)
## Naive tied-read policies (NAIVE-BEST, NAIVE-ALL) against minimap2's primary: read level and copy level, 2026-10-22

Question (user): as a naive O2, impute that the right copy is "the one within 10 AS", and test whether that improves on the primary placement and on loci determination. Read-only analysis of the frozen audit `/mnt/linuxdisk/tmp/o3_graph/o2vc/audit_frozen_eccbc4e0/{hg002,hg00733}` (calls.tsv, `<family>.reads.tsv`, `<family>.case.tsv`, score.C.json, score.E.json); no decision rule was fixed beforehand, so every number is descriptive. Scripts (not committed): `/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/e542c6ba-28cf-42e7-becf-fd99d54f830f/scratchpad/naive_tied.py` (analysis) and `.../build_md.py` (this layout).

### Definitions and checks
- Tied read: alignments within 10 AS of the best on >= 2 copies of one family. On every read the AS spread across its tied copies is <= 10 (checked: 0 violations).
- **NB-abst**: place on the tied copy with the strictly highest AS, abstain on an exact tie at the max (scored reads: 622 of 6333 in HG002, 5132 of 12015 in HG00733). **NB-prim**: same, but at an exact tie keep minimap2's primary (N1, even when it is not one of the max-tied copies: 156 / 25 tied reads, HG002 / HG00733).
- **NAIVE-ALL**: the read is credited to every tied copy (`candidates`). Not a single placement, so it has no wrong/placed; it is scored by set containment and, at copy level, by credits.
- Scoring as in `cmd_score`: scored = tied reads of class unique, unique_elsewhere, elsewhere, no_counterpart; right iff the placement is the truth copy (elsewhere / no_counterpart have none, so any placement is wrong); fixes = right where N1 is wrong, breaks = wrong where N1 is right; PAR families (par_families.json) excluded (15 HG002 families with tied reads, none in HG00733). Of the 20100 non-PAR tied HG002 reads 6333 are scored (not scored, no definite truth: own_tie 9642, poor_match 3800, no truth row 325); all 12015 HG00733 tied reads are.
- Checks: my N1, N2, N4, VC placed/right/wrong/abstained/fixes/breaks equal score.C.json and score.E.json exactly (asserted); NB-abst / NB-prim / N1 were re-derived a second way from the `case.tsv` truth columns and agree to the read.

### 1. Read level

**HG002**, left block all four scored classes (N = 6333), right block unique + unique_elsewhere only, i.e. without elsewhere / no_counterpart (N = 6083); abstained = N - placed.

| policy | placed | right | wrong | wrong/placed | fixes | breaks | placed | right | wrong | wrong/placed | fixes | breaks |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N1 (minimap2 primary) | 6333 | 3875 | 2458 | 0.388 | - | - | 6083 | 3875 | 2208 | 0.363 | - | - |
| N2 (MAPQ >= 1) | 5725 | 3647 | 2078 | 0.363 | 0 | 0 | 5507 | 3647 | 1860 | 0.338 | 0 | 0 |
| VC | 4579 | 4308 | 271 | 0.059 | 1137 | 24 | 4393 | 4308 | 85 | 0.019 | 1137 | 24 |
| N4 (corrector) | 6333 | 5096 | 1237 | 0.195 | 1252 | 31 | 6083 | 5096 | 987 | 0.162 | 1252 | 31 |
| **NB-abst** | 5711 | 3645 | 2066 | 0.362 | 36 | 38 | 5505 | 3645 | 1860 | 0.338 | 36 | 38 |
| **NB-prim** | 6333 | 3873 | 2460 | 0.388 | 36 | 38 | 6083 | 3873 | 2210 | 0.363 | 36 | 38 |

**HG00733** (N = 12015, every read class unique):

| policy | placed | right | wrong | wrong/placed | fixes vs N1 | breaks vs N1 |
|---|---|---|---|---|---|---|
| N1 (minimap2 primary) | 12015 | 7420 | 4595 | 0.382 | - | - |
| N2 (MAPQ >= 1) | 6953 | 5223 | 1730 | 0.249 | 0 | 0 |
| VC | 2801 | 2625 | 176 | 0.063 | 279 | 30 |
| N4 (corrector) | 12015 | 7788 | 4227 | 0.352 | 408 | 40 |
| **NB-abst** | 6883 | 5410 | 1473 | 0.214 | 342 | 112 |
| **NB-prim** | 12015 | 7650 | 4365 | 0.363 | 342 | 112 |

wrong/placed is not comparable across arms that abstain (N2, VC, NB-abst place fewer reads): compare NB-abst with N2, which abstains on the same exact ties, and the others by `right` out of N.

**NAIVE-ALL** (the set of tied copies as the answer):

| dataset | scored reads | reads with a truth copy | set contains truth | fraction of scored | fraction of truth-copy reads | mean set size | credits on the truth copy / all credits | N1 right / scored |
|---|---|---|---|---|---|---|---|---|
| HG002, four classes | 6333 | 6083 | 5892 | 0.930 | 0.969 | 2.06 | 0.470 | 0.612 |
| HG002, unique + unique_elsewhere | 6083 | 6083 | 5892 | 0.969 | 0.969 | 2.06 | 0.470 | 0.637 |
| HG00733 | 12015 | 12015 | 11607 | 0.966 | 0.966 | 2.37 | 0.407 | 0.618 |

Truth outside the tied copies (never containable): 191 reads in HG002, 408 in HG00733. Set sizes in HG002: {2: 5973, 3: 337, 4: 15, 5: 5, 6: 1, 7: 2}; HG00733: {2: 9566, 3: 1270, 4: 850, 5: 110, 6: 49, 7: 71, 8: 42, 9: 56, 10: 1}.

### 2. Copy level (loci determination)
Credits per family over the scored tied reads. TRUTH = truth copy of unique / unique_elsewhere reads; N1 = primary; VC = placed reads only; VC+N1 = VC, else N1; NB-prim; NAIVE-ALL = +1 to every tied copy. Expressed = credited >= 2 reads; sens / prec pooled over (family, copy) against the copies with >= 2 TRUTH reads. A primary on no copy counts as a credit (false expression, TVD mass) but not as a called copy. TVD = 0.5 * sum |policy fraction - TRUTH fraction| per family, over families with a TRUTH read; weighted by the family's TRUTH reads; VC (and NB-abst) are undefined where they place nothing, so the family count is shown.

**HG002, tied reads only** (TRUTH: 6083 reads on 117 copies of 61 families, 87 expressed copies); left block four classes (6333 reads, 250 credited but with no truth copy), right block unique + unique_elsewhere (6083 reads):

| policy | called | sens | prec | zero-truth credits | TVD med | TVD wtd | called | sens | prec | zero-truth credits | TVD med | TVD wtd |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N1 | 105 | 0.874 | 0.724 | 743/6333 = 0.117 | 0.187 | 0.347 | 97 | 0.874 | 0.784 | 700/6083 = 0.115 | 0.187 | 0.347 |
| VC | 48 | 0.471 | 0.854 | 17/4579 = 0.004 | 0.083 | 0.126 | 45 | 0.471 | 0.911 | 7/4393 = 0.002 | 0.049 | 0.129 |
| VC+N1 | 105 | 0.885 | 0.733 | 728/6333 = 0.115 | 0.100 | 0.166 | 97 | 0.885 | 0.794 | 683/6083 = 0.112 | 0.060 | 0.165 |
| NB-abst | 99 | 0.839 | 0.737 | 530/5711 = 0.093 | 0.153 | 0.338 | 91 | 0.839 | 0.802 | 514/5505 = 0.093 | 0.132 | 0.338 |
| NB-prim | 110 | 0.897 | 0.709 | 731/6333 = 0.115 | 0.180 | 0.345 | 102 | 0.897 | 0.765 | 700/6083 = 0.115 | 0.164 | 0.344 |
| NAIVE-ALL | 147 | 0.954 | 0.565 | 1658/13062 = 0.127 | 0.500 | 0.469 | 135 | 0.954 | 0.615 | 1598/12543 = 0.127 | 0.500 | 0.468 |

TVD on the 35 four-class families where every policy has credits (median / read-wtd mean): N1 0.182/0.327, VC 0.083/0.126, VC+N1 0.079/0.130, NB-abst 0.125/0.315, NB-prim 0.164/0.324, NAIVE-ALL 0.499/0.464. VC places a scored read in 35 of the 61 families with a truth read, NB-abst in 59.

**HG00733, tied reads only** (12015 reads; TRUTH: 12015 reads on 352 copies of 174 families, 294 expressed copies):

| policy | copies called (>= 2 reads) | sens | prec | credits on zero-truth copies | TVD median (families) | TVD read-wtd mean |
|---|---|---|---|---|---|---|
| N1 | 330 | 0.939 | 0.836 | 762/12015 = 0.063 | 0.136 (174) | 0.248 |
| VC | 158 | 0.531 | 0.987 | 71/2801 = 0.025 | 0.064 (116) | 0.304 |
| VC+N1 | 331 | 0.956 | 0.849 | 729/12015 = 0.061 | 0.093 (174) | 0.236 |
| NB-abst | 272 | 0.833 | 0.901 | 398/6883 = 0.058 | 0.086 (168) | 0.226 |
| NB-prim | 329 | 0.932 | 0.833 | 772/12015 = 0.064 | 0.143 (174) | 0.227 |
| NAIVE-ALL | 399 | 0.997 | 0.734 | 2741/28533 = 0.096 | 0.373 (174) | 0.358 |

TVD on the 114 HG00733 families where every policy has credits (median / read-wtd mean): N1 0.138/0.254, VC 0.069/0.304, VC+N1 0.090/0.234, NB-abst 0.085/0.173, NB-prim 0.143/0.213, NAIVE-ALL 0.327/0.327.

**Expressed copies that N1 misses versus copies NAIVE-ALL adds** ("beyond N1" = called by NAIVE-ALL, not by N1; "expressed" = >= 2 TRUTH reads):

| analysis | TRUTH-expressed copies | N1 misses | NAIVE-ALL finds | VC+N1 finds | NAIVE-ALL calls beyond N1 | of which expressed | with 1 truth read | with 0 truth reads |
|---|---|---|---|---|---|---|---|---|
| HG002 tied, four classes | 87 | 11 | 7 | 1 | 44 | 7 | 2 | 35 |
| HG002 tied, unique + unique_elsewhere | 87 | 11 | 7 | 1 | 39 | 7 | 4 | 28 |
| HG00733 tied | 294 | 18 | 17 | 5 | 69 | 17 | 17 | 35 |
| HG002 tied + untied, unique + unique_elsewhere | 407 | 22 | 4 | 0 | 15 | 4 | 1 | 10 |
| HG00733 tied + untied | 656 | 2 | 1 | 0 | 14 | 1 | 0 | 13 |

**Tied + untied reads** (untied = 205524 HG002 reads of unique / unique_elsewhere in truth_untied.tsv, placed on minimap2's primary, which is right on 98.9% of them; HG00733: 71685 reads, placed on the truth copy by assumption). TRUTH HG002: 211607 reads on 486 copies of 222 families, 407 expressed; HG00733: 83700 reads on 657 copies of 274 families, 656 expressed.

| policy | HG002 sens | prec | zero-truth credits | TVD wtd | HG00733 sens | prec | zero-truth credits | TVD wtd |
|---|---|---|---|---|---|---|---|---|
| N1 | 0.946 | 0.941 | 1927/211607 = 0.009 | 0.020 | 0.997 | 0.965 | 717/83700 = 0.009 | 0.036 |
| VC | 0.919 | 0.964 | 1267/209917 = 0.006 | 0.015 | 0.980 | 0.998 | 70/74486 = 0.001 | 0.039 |
| VC+N1 | 0.946 | 0.944 | 1915/211607 = 0.009 | 0.014 | 0.997 | 0.965 | 699/83700 = 0.008 | 0.034 |
| NB-abst | 0.941 | 0.950 | 1752/211029 = 0.008 | 0.020 | 0.994 | 0.979 | 381/78568 = 0.005 | 0.034 |
| NB-prim | 0.951 | 0.942 | 1921/211607 = 0.009 | 0.019 | 0.997 | 0.965 | 725/83700 = 0.009 | 0.033 |
| NAIVE-ALL | 0.956 | 0.920 | 2514/218067 = 0.012 | 0.021 | 0.998 | 0.947 | 2522/100218 = 0.025 | 0.072 |

With the 260772 untied no_counterpart / elsewhere reads kept (HG002 four classes), every policy puts ~55% of its credits on zero-truth copies (labels partly liftover artifacts) and prec differs by only 0.034 between N1 and NAIVE-ALL, so that variant is dominated by those reads.

### Caveats
- HG002 `elsewhere` / `no_counterpart` labels partly come from liftover-chain artifacts; hence both variants above. HG00733 is a simulation from the donor haplotypes, aligned to CHM13, not a real-read result.
- The >= 2-reads rule is lenient for NAIVE-ALL: two reads tied on the same pair make both copies "expressed" by construction. Tied-only copy sets are small (HG002 61 families with a truth copy, HG00733 174).
- Untied reads dilute every difference (the policies agree on them); HG00733 untied reads are credited to the truth copy because the frozen outputs hold no minimap2 placement for them, an upper bound for N1.

### Reading
- NAIVE-BEST gives no gain over the primary in HG002 (NB-prim wrong/placed 0.388 vs N1 0.388, 36 fixes / 38 breaks; NB-abst 0.362 vs N2 0.363) and a modest one in the simulation (0.363 vs 0.382, 342 / 112; NB-abst 0.214 vs N2 0.249), far from VC (0.059 / 0.063): the highest tied AS carries little beyond the primary's own pick.
- NAIVE-ALL does recover expressed copies the primaries miss (HG002 7 of 11, HG00733 17 of 18; VC+N1 only 1 and 5), its set holds the truth for 96.9% / 96.6% of reads, at 2.06 / 2.37 copies per read.
- It also invents expression: 39 / 69 extra copies called beyond N1, only 7 / 17 expressed (18% / 25%) and 28 / 35 with zero truth reads (~4 / ~2 per recovered copy); prec 0.615 vs N1 0.784 and 0.734 vs 0.836, TVD median 0.500 / 0.373 vs 0.187 / 0.136: a recall-only lever (HG002 / HG00733).


---

## PREREG_o2_variant_calls_2026-10-22.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_variant_calls_2026-10-22.md`)*

## Pre-registration: O2 from the sample's own variant calls (2026-10-22)

**User (2026-10-22):** "it just needs to be something that passes from beginning to end, like if I posit that I have a gene
family I need to show it (lets say DAZ or so) then show how if there are PSVs, bam tags or other differences among the possible
placements and show that we got the right placement instead of minimap2 or that we found a possibly missing copy based on the
variants" — one family picked by hand (GTF2IRD2), with backups.

**Design:** `docs/superpowers/specs/2026-10-22-o2-variant-calls-design.md`. The spec's "How reads are scored" and "Not
evaluable" amendments are in commit c63bbd6a. Plan: `docs/superpowers/plans/2026-10-22-o2-variant-calls.md`.

**Frozen code:** `07441e7c` on branch `o3-fast`, holding `bench/o3_graph/o2vc.py`, `o2vc_data.py`, `o2vc_rule.py`. The
weak-het policy is `P1 `keep` (the record's GT alleles whatever FILTER; `--weak-het` default, commit fdb41f25)`, chosen on the development families by the rule in "What was seen before".

**Expected failures, pre-registered separately:** `docs/PREREG_o2_variant_calls_expected_failures_2026-10-22.md` (0b31688c). It
holds 16 failure modes with their signatures, and blind panel predictions.

### The rule (verbatim from the spec's method, summarised)

- **Domain.** A tied read is one whose alignments within 10 AS of its best lie on ≥ 2 copies of one family and on no locus off it.
  minimap2's primary may lie outside the tied set.
- **Teachers.**
  - A copy's teachers are reads primary on it whose every alignment on another copy of the family is ≥ 10 AS lower. A tied read
    never teaches.
  - longcallR 2.0.0 (`-p hifi-masseq --min-mapq 0 --no-bam-output`, defaults otherwise) calls SNVs on each copy's teachers.
- **Genotypes.**
  - **Known:** the copy is known at a position when depth ≥ 6, counted as longcallR's DP counts it. The genotype is the GT
    alleles whatever FILTER, with the weak-het policy `P1 `keep` (the record's GT alleles whatever FILTER; `--weak-het` default, commit fdb41f25)` applied.
  - **Unknown:** `RnaEdit`, `Dense` and depth > 50,000 are unknown.
  - **No record:** {reference base}.
- **Assign or abstain.**
  - **Separates:** a column separates two copies when both are known and their allele sets are disjoint.
  - **Votes:** the read's base, at Q ≥ 10, votes for the copy whose set holds it.
  - **Beats:** c beats d when ≥ 1 separating column votes for c and none votes for d.
  - **Decision:** assign c iff c beats every other tied copy; else abstain with the set of unbeaten copies.
  - **Frame:** all comparisons are made in the read's frame, so a copy aligned opposite to the read has its alleles complemented.
- **Possible missing copy.**
  - **Foreign:** an allele is foreign at a column where every tied copy is known and none holds it.
  - **Flag:** a flag is a connected component of recurrent foreign alleles (≥ 2 reads each) that are linked (≥ 2 reads carry
    both, at different columns).
  - **Editing:** A>G changes in the transcript frame never count toward a flag.
  - **Effect:** a flagged read is not placed.

### Arms (all on the same reads)

| arm | definition |
|---|---|
| N1 | minimap2's primary |
| N2 | N1 at MAPQ ≥ 1 |
| N3 | the same rule, every known genotype = {reference base} |
| N4 | the frozen corrector (`o2correct.py donor` v2), where its output holds the read |
| VC | this design |

### Panels and exact commands (run once each, in this order, after this file is committed)

```bash
cd /mnt/linuxdisk/home/juanfraitu/rustle_soto_subset
T=/mnt/linuxdisk/tmp/o3_graph; G=/mnt/linuxdisk/home/juanfraitu/npip_ladder/idx/target.fa; P=$T/famident/par_families.json
## A + C-human: HG002 real Kinnex reads, every pool family (chrM families out); truth = HG002's own assembly (answer key only)
python3 bench/o3_graph/o2vc.py run $T/hg002proof/chm13.bam $T/o2psv_dev/pool_copies.tsv $G $T/o2vc/hg002 --exclude-chrM
for F in S19GTF2IRD2 S19NBPF1 S19PDPR; do
  python3 bench/o3_graph/o2vc.py score $T/o2vc/hg002 $T/hg002proof/truth.tsv $T/hg002proof/truth_untied.tsv --mode donor \
    --n4 $T/o2psv_dev/hg002.v2.tsv --panel A --focus $F --exclude $P
done
python3 bench/o3_graph/o2vc.py score $T/o2vc/hg002 $T/hg002proof/truth.tsv $T/hg002proof/truth_untied.tsv --mode donor \
  --n4 $T/o2psv_dev/hg002.v2.tsv --panel C --exclude $P
## B: DAZ, A119b testis, no truth (descriptive)
python3 bench/o3_graph/o2vc.py run $T/yag_testis/yag.bam $T/o2vc/inputs/daz_copies.tsv $G $T/o2vc/daz --families S19DAZ1
python3 bench/o3_graph/o2vc.py score $T/o2vc/daz --mode none
## E: HG00733 simulated from its own haplotypes, held out (never run or scored during development)
S=$T/donorsim/HG00733/sim
python3 bench/o3_graph/o2vc.py run $S/chm13.bam $T/o2psv_dev/pool_copies.tsv $G $T/o2vc/hg00733
python3 bench/o3_graph/o2vc.py score $T/o2vc/hg00733 $S/truth.tsv --mode donor --n4 $S/second/o2correct.tsv --panel E --exclude $P
```

**Disclosed precompute.**
- **What was run early:** the longcallR calls (`OUTDIR/vcf/*.vcf`) for `hg002` and `hg00733` were computed BEFORE this file was
  committed. They come from a frozen checkout at a2098a96, in batches, while the policy was being chosen. This was to save time:
  the user has a meeting.
- **Why it is safe:** calls do not depend on the weak-het policy, and `run` reuses saved VCFs. The development-time changes after
  a2098a96 did not touch calling (`git diff a2098a96..07441e7c` touches no calling code. It changes only `parse_vcf(with_af)`, the weak-het interpretation of saved VCFs, the attempt counter and exit code, the scorer and the CLI.).
- **What was not looked at:** the truth-free placements those batches wrote were neither opened nor scored. The commands above
  regenerate every output except the VCFs with the frozen code.

**Deferred (fast path, user 2026-10-22):**
- panel C-gorilla (KB3781 NPIP unit C);
- panel D (one-genome no-harm worlds).

Each gets its own amendment before it runs.

### Scoring (spec, as amended in c63bbd6a)

- **Scored reads.** These are tied reads with a definite truth: classes unique, unique_elsewhere, elsewhere and no_counterpart.
  - A placement is right iff the placed copy is the truth copy. Every other placement on a scored read is wrong. That includes a
    truth copy outside the tied set, `elsewhere`, and `no_counterpart`.
  - Fixes, breaks and wrong/placed are counted over all scored reads, for every arm, on the same reads. A break is a read N1 has
    right that the arm places wrong; this includes outside-primary reads VC moves.
  - The decidable subset (truth among the tied copies) is reported for comparison with earlier registers.
- **Calibration counts.**
  - K1 = VC placements among `own_tie` tied reads (HG002 only).
  - K2 = flagged tied reads whose truth copy is present (unique, unique_elsewhere, own_tie; in E every read).
  - C1 = flags with ≥ 2 `no_counterpart` reads.
- **Exclusions.**
  - PAR families (`famident/par_families.json`) are out of every pooled number and reported apart.
  - chrM families are out of the HG002 run.
- **Not evaluable.** A bar is not evaluable when there is nothing to measure, or when no outcome could meet it for lack of reads:
  - a zero denominator;
  - A2 with < 10 minimap2-wrong scored reads on the focus family;
  - C1 with < 2 `no_counterpart` tied reads;
  - E2 with no VC placement.
- **Verdict precedence.** FAIL if any evaluable bar fails; else UNDERPOWERED if any bar is not evaluable; else PASS.
- **The "fix possible" bound** for A2, A3 and E2 counts every scored read that N1 places wrong, as implemented. That includes
  `elsewhere` / `no_counterpart` reads and reads whose truth copy lies outside the tied set, none of which any arm can fix. It is
  therefore an upper bound on the possible fixes (Task 6 review, minor 1).
  - At scoring, the decidable N1-wrong count (`decidable.arms.N1.wrong`) is reported next to it.
  - If a bar is evaluable only by the upper bound, the result says so.

### Bars (verbatim from the spec)

- **A (focus GTF2IRD2; the backups NBPF1 and PDPR get the same bars and are reported):**
  - A1 VC wrong ≤ 2% of VC placements;
  - A2 fixes ≥ 10;
  - A3 VC fixes > N3 fixes;
  - K1 ≤ 2% and K2 ≤ 2%, pooled over HG002's pool families, PAR out.
- **B:** descriptive.
- **C-human:**
  - C1: at least one flag holds ≥ 2 `no_counterpart` tied reads;
  - C2 (reported): every flag with its reads' truth classes. The audit checks each counted flag's donor locus for a CHM13 copy.
- **E:**
  - E1 wrong ≤ 2% of placements;
  - E2 fixes > breaks;
  - K2 ≤ 2%.

**Decision.**
- The showcase passes iff panel A for GTF2IRD2 passes.
- A failed GTF2IRD2 is a FAIL. The backups are reported, never substituted.
- C-human and E each get their own verdict, never pooled. Human and gorilla are never pooled.

### What was seen before this file

- **The showcase choice.** The corrector's HG002 per-family numbers on tied reads with donor truth chose it: GTF2IRD2 19 fixes /
  0 breaks, NBPF1 11 / 0, PDPR 14 / 0. So did samplecopy's held-out GTF2IRD2 31 / 0 (register 1333). The showcase is hand-picked,
  and hand-picked sets overstate (metric trap 09-18).
- **Development, HG002 SMN2 and MSTRG3150 only** (`docs/O2VC_DEV_2026-10-22.md`): pooled over both development families (1,440 scored tied reads, every arm on the same reads): **VC 1,330 placed, 1 wrong, 647 fixes, 0 breaks**; N1 (minimap2) 656 wrong of 1,440; N2 623 wrong of 1,358; N3 (reference only) 692 placed, 15 wrong, 8 fixes; N4 (corrector) 19 wrong, 638 fixes, 1 break; K1 0/128; K2 3/1,568; depth agreement 41/41. Policy choice: P3 (`ref`) excluded (K1 98/128); P1 (`keep`) and P2 (`unknown`) tie on every term (identical VC on all 1,849 reads), so the spec's P1 stays. ⚠On SMN2, VC's whole gain over N3 rests on ONE LowQual 0/1 record (chr5:70823918, AF 0.14, DP 79, exon 2a) that erases one reference PSV; the audit re-calls SMN2 and GTF2IRD2 five times each.
- **The SMN2 smoke** (Task 4): one LowQual 0/1 at AF 0.14 erased a reference PSV for 1,589 of 1,655 tied reads.
- **P1 vs P2.** The frozen rule does not separate them. VC is identical on all 1,849 dev tied reads; only N3 differs, on 922
  SMN2 reads. P1 is kept as the spec's own rule, the status quo; it is not derived from the rule.
- **Policy choice rule (frozen before scoring the dev runs):**
  1. Exclude any policy whose dev K1 is above 2%.
  2. Take the lowest dev wrong/placed.
  3. Break ties by more fixes.

### Caveats carried into every claim

- **The caller is not held out on HG002.** longcallR's PacBio MAS-seq model was trained on PacBio's HG002 Kinnex cDNA set.
  Provenance of our HG002 reads: the `chm13.bam` header and run IDs (m84036_230422…) look like PacBio's public HG002 Kinnex release, so the reads may be longcallR's own MAS-seq training data. The headers cannot prove it. Any HG002 result is therefore caller-optimistic..
- **SNVs only.** Indel PSVs are not used.
- **Absorbed extra copies are not flagged.** An extra donor copy absorbed by its nearest reference copy is never flagged, and is
  visible only through `no_counterpart` placements.
- **E is simulated.** It has even allele expression, no editing and no donor-extra copies, so it is optimistic.
- **Non-determinism.** longcallR is not deterministic at weak sites (4 of 41 identical runs). The saved VCFs are the record. The
  audit re-calls GTF2IRD2 five times and reports how many placements change.
- **Novelty.** This is DuploMap's reliable-PSV principle (Prodanov & Bansal 2020, DNA) carried to long-read RNA, with an
  assign-or-abstain certificate and an RNA-side missing-copy flag.

### After the runs

1. Append "Results (code `07441e7c`; run once)":
   - one table per panel, with each bar's value and verdict;
   - GTF2IRD2's worked reads from its `summary.md`.
2. Audit with the `metric-auditor` agent on a frozen copy of the outputs.
3. Add register row 1374 and update memory.

### Results (code 07441e7c; run once, 19:54–20:10 PDT; audit 20:31)

Outputs are frozen in `/mnt/linuxdisk/tmp/o3_graph/o2vc/audit_frozen_eccbc4e0/`. The audit (metric-auditor) recomputed every
headline from `calls.tsv` and the truth tables, without the scorer; its scripts are in `/mnt/linuxdisk/tmp/audit_2026-10-22_o2vc/`.

| panel | verdict (scorer) | audit | numbers |
|---|---|---|---|
| A GTF2IRD2 (showcase) | PASS | HOLDS WITH CAVEATS | VC 44 placed / 0 wrong / 24 fixes / 0 breaks of 61 scored; N1 37 wrong; N3 9 placed, all wrong, 1 break; N4 7 wrong, 30 fixes. All 44 placements are on copy .1, so the 95% upper bound on the error rate is 6.6%. The decidable N1-wrong count is 37, so A2 is evaluable without the upper bound. |
| A PDPR (backup) | PASS | HOLDS WITH CAVEATS | VC 36 / 0 / 21 fixes. All placements are on PDPR.0; VC abstained on all 82 reads whose truth is PDPR.3. |
| A NBPF1 (backup) | FAIL | HOLDS | VC 14 placed, 9 wrong (6 elsewhere + 3 no_counterpart). The 3 are probably right, since a chain gap lines them up with NBPF1.3. That still leaves 6/14 wrong. |
| K1 / K2 (HG002 pooled, PAR out) | PASS | HOLDS | K1 49/9,642; K2 12/15,725 |
| C-human | PASS | **DOES NOT HOLD** | The 16 S19NPIPB2#0 reads are the donor's NPIPB2.2 with a 57-bp insertion, in a gap of the HG002→CHM13 chain. The insertion pushes .2 out of the 10-AS tie. Relabelled, C1 has 0 counted flags and **FAILS** (it is evaluable: 107 no_counterpart reads remain). |
| B DAZ | descriptive | HOLDS | 622 tied reads; VC and N3 place 0; N2 places 15; 5 reads flagged in 2 flags (no truth). |
| E HG00733 (held out) | **FAIL** (E1) | HOLDS | VC 176 wrong / 2,801 placed (6.3%); E2 279 fixes > 30 breaks; K2 43/12,015. PRAMEF7 + PRAMEF10 = 133 of the 176 wrong and 23 of the 30 breaks (pre-registered F3). |

**Robustness.**
- Five fresh longcallR re-calls of S19GTF2IRD2 and S19SMN2 gave byte-identical VCFs and 0 placement changes.
- The weak SMN2 record (chr5:70823918) was kept every time.
- Depth agreement: 2,571 of 2,571 records.

**Pooled HG002** (reported only; no bar):
- VC 271 wrong / 4,579 placed, 1,137 fixes, 24 breaks; N1 2,458 / 6,333 wrong; N4 1,237 wrong.
- Without the development families: 270 / 3,249 = 8.3%.
- With the 801 poor_match reads scored against the donor's best placement: 339 / 5,366 = 6.3%.
- Of the 18 families with ≥ 10 placements, 10 are at ≤ 2% wrong. The worst are:
  - GTF2H2C 145/166 (all elsewhere/no_counterpart labels, which depend on the chain and are unverified);
  - SPDYE1 59/77;
  - ANKRD36C 21/64 (14 breaks);
  - AC0064533 11/21;
  - SLC35E2B 9/47.

**Reading.**
- **Showcase.** It passes, but on the development donor, with a family chosen from earlier HG002 results.
- **Held-out donor.** It fails E1, mainly on PRAMEF, the pre-registered donor-lacks-a-copy case.
- **Missing-copy flag.** The pass is a liftover artifact. The rule did correctly refuse to place those reads on the wrong copies.
- **Against the other methods.** The method is far more precise than minimap2, the corrector and N2 on the same reads (HG002 5.9% vs
  38.8% / 19.5% / 36.3% wrong; HG00733 6.3% vs 38.2% / 35.2% / 24.9%). It does not reach the 2% bar on the held-out donor or across
  families.
- **Where the errors are.** They concentrate in copy-number-variable families: PRAMEF, GTF2H2C, SPDYE1, NBPF1 and ANKRD36C. That
  is the case for an O3 step before O2.
- **Blind predictions.** They were too pessimistic for A and too optimistic for C and E.


---

## PREREG_o2_variant_calls_expected_failures_2026-10-22.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_variant_calls_expected_failures_2026-10-22.md`)*

## Pre-registered expected failures: O2 from the sample's own variant calls (2026-10-22)

**Status.** Committed BEFORE the development run (Task 6) and before any panel run (Task 9), at the user's request (2026-10-22):
"explore [the ThesisVault] and see if we can corroborate what we are doing plus pre-register possible expected failures".

**How it was made.** A blind literature survey by a subagent:
- read the spec, the SDD ledger and memory, all 28 vault paper notes and the cited full texts;
- read no panel output and no truth file.

The controller adopted it unchanged below the line, apart from this header.

**Already addressed when this was committed** (so these predictions are about the code as it will run):
- **§0.7, F2 and F16(3) ("placed `no_counterpart` / `elsewhere` reads are invisible to A1/E1").** Fixed by spec amendment
  c63bbd6a and the Task 5 fix round:
  - every tied read with a definite truth is scored;
  - those placements count as wrong;
  - every arm is scored on the same reads.
- **F16(4): T6's policy choice ignored K1.** T6 will also report K1 for each genotype policy (P1/P2/P3). The choice rule stays
  dev precision, then fixes. A policy whose dev K1 exceeds 2% is excluded.
- **§2.18 reported-only counts.** Adopted as a diagnostics step in Task 6/9, reported only, no new bars.

**Caveat to carry into every claim (§1 row 2, F15).** longcallR's PacBio MAS-seq model was trained on PacBio's HG002 Kinnex
cDNA set. That set is very likely the HG002 reads of panels A and C, so HG002 is not held out for the caller (to verify).

**Novelty framing (§3.3).** "DuploMap's reliable-PSV principle (Prodanov & Bansal 2020, DNA) carried to long-read RNA, with
an assign-or-abstain certificate and an RNA-side missing-copy flag." It is not a new principle. Search for an RNA analogue
before claiming "first".

---

## O2 from the sample's own variant calls: vault survey, pre-registered expected failures, novelty delta

**Date:** 2026-10-22.
**Destination:** to be committed as `docs/PREREG_o2_variant_calls_expected_failures_2026-10-22.md` BEFORE T6 and T9 run.
**Blind:** I read no panel output, no truth file and no directory under `o3_graph/o2vc/`. The SMN2 dev-smoke facts come only from
the ledger and the Task-4 report. These predictions must not be tuned toward.

**What I read:**
- the spec, including "Context: related work and known gaps";
- the SDD ledger, the Task-4 and Task-5 rulings, and the Task-4 report;
- memory `project_o2_variant_calls_2026-10-22` and `reference_paraphase_parascopy`, plus the O2 notes they rest on (r1329, r1330,
  r1333, r1334, r1336, r1343-r1348, the corrector note, `feedback_metric_traps`);
- the docstring of `bench/o3_graph/donorsim.py`;
- every vault paper note (28) and the Topics notes;
- the full text behind every claim cited here.

**Citation key** (short name = vault paper note):
- DM = Prodanov_2020_DuploMap; PSC = Prodanov_2022_Parascopy;
- PP23 = Chen_2023_Paraphase_SMN; PP25 = Chen_2025_Paraphase;
- LR = Huang_2026_LongRNA_SNP_Phasing (longcallR); C3R = Zheng_2025_Clair3-RNA;
- SDA = Vollger_2019_SDA; CH = Chaisson_2017_Multicopy_Duplications_Polyploid_Phasing; IC = Sahlin_2018_IsoCon;
- DV = Shaw_2025_devider; PH = Hosseini_2025_pHapCompass;
- CZ = Canzar_2016_Multimapping_Conflicts; DS = Deschamps_2020_MultiMapped_Review;
- PAL = Pal_2026_Y_Multicopy_Survival; DSH = Dishuck_2025_NPIP; GU = Guitart_2024_TBC1D3;
- GR = Greshnova_2024_YAG_T2T_Isoforms; RH = Rhie_2023_T2T_Y.

`§` names the full-text section. `rNNNN` is a negative-results-register row, quoted from memory. "(inf.)" marks my own inference,
not a paper's claim.

### 0. What the literature adds to the spec (read first)

1. **The core loop is published, for DNA.** DM genotypes reference PSVs per copy from the sample's confidently placed long reads. It
   decides only with PSVs the sample confirms at both copies, and gives MAPQ 0 when no location beats all the others. O2-VC is that
   loop moved to RNA:
   - restricted to AS-tied reads;
   - with a hard certificate and a linked-allele flag;
   - without DM's indel PSVs and without DM's iteration.
2. **Teacher error is about 2%, measured twice.**
   - DM §Results: after PSV re-scoring, 1.9% of HG002 HiFi reads with MAPQ ≥ 30 fall below 10.
   - r1330: 2.1% of non-PAR untied HG002 reads lie off the donor's copy, concentrated in CNV SDs.

   The leaked allele's allele fraction is εr/(1+εr), where ε is the teacher error and r how much more the rival copy is expressed
   (inf., §2.0). With ε = 0.02 it crosses longcallR's 0.10 floor once the rival copy has about 5.6 times the teachers. SMN2's
   AF-0.14 het is the expected product (r ≈ 8).
3. **longcallR is used off-label.** LR filters MAPQ < 10. It trains and benchmarks only inside GIAB high-confidence regions, at
   coverage ≥ 20. Its PacBio MAS-seq model was trained on PacBio's HG002 Kinnex full-length RNA set, very likely panels A and C's own
   reads (verify that `hg002_kinnex` is that PacBio release). C3R's stratified benchmark ranks SDs among the two worst strata, and there longcallR averages about 66% F1.
4. **Caller misses become reference votes.** A missed call or a HomRef record at depth ≥ 6 becomes {reference base}, which is an N3
   vote. In SDs, caller false negatives therefore re-import exactly the reference-PSV errors VC exists to remove. Substitutions alone
   were a coin flip outside SMN/THOC3: 17 fixes, 18 breaks (memory, corrector).
5. **"Dense = collapsed copy" is unsupported for the copies O2 sees.** Reads tie only between copies ≥ 99.5% identical (memory, 10-19).
   A collapsed copy of that kind adds ≤ 0.5 candidates per 100 bp, against the Dense trigger of ≥ 5 (inf.). Dense will mark alignment
   artifacts.
6. **Teachers are selected by their alleles.** A copy's reads that carry its rival's allele are the reads that tie, so they rarely
   teach. DV names the mechanism: reference bias. The genotypes are pushed toward separating alleles (false homozygotes), and a copy's
   second haplotype can look "foreign" (false flags).
7. **Three scoring facts the bars hide.**
   - With fewer than 50 VC placements, A1 ≤ 2% means zero wrong.
   - Placed `no_counterpart` reads are invisible to A1 and E1 (metric trap 09-22).
   - On `outside = 1` reads N1 is wrong by construction, which inflates fixes. T6's criterion (dev precision) ignores K1.

### 1. Corroboration

| # | Element | Supports | Contradicts or qualifies | Does the literature change anything? |
|---|---|---|---|---|
| 1 | Teachers by an AS margin (≥ 10 over every family rival) | GU §Methods "Iso-Seq and transcript analyses", the same rule: reads "with primary minimap2 alignment scores of 10 or greater for a given paralog cluster relative to all other cluster mappings were retained, while other mappings were marked as ambiguous and ignored". DSH §Results: the best mapping must beat the next "by at least one additional mismatch". DM §Methods step 2: PSV genotypes come from reads "aligned with high confidence". | GU and DSH measure the margin against panels of assembled human haplotypes (GU: "all validated human haplotypes"; DSH: 169), not against one reference. A margin against CHM13 copies does not certify the donor's copy. DM §Results: 1.9% of MAPQ-≥ 30 HiFi reads demoted. PP25 §Results: "even reads with high MAPQ can be misaligned due to reference genome artifacts, common CNVs, and high rates of gene conversion". r1330: 2.1% off the donor's copy. r1334: GTF2IRD2's and HSPA6's sample clusters were labelled by near-tie reads, because margin-≥ 10 reads were too few to label them. | Keep the rule, since it is the field's. Treat the teachers as about 2% contaminated, most where expression is skewed and copy number varies (§2 F1-F3, F5). GTF2IRD2's thin margin-≥ 10 support is a known risk for A2. |
| 2 | Per-copy genotypes from an off-the-shelf RNA caller | DM §Methods "Identifying reliable PSVs" makes the same move in DNA: per-copy PSV genotypes from confident reads, and only PSVs that are reference-consistent at both copies are used. LR §Results: PacBio precision 98.5-99.0%, recall 91.0-95.3%, with an A-to-I class. | Off-label on four counts. (i) LR §Methods "candidate SNP sites": MAPQ ≥ 10, secondary and supplementary reads dropped; the design runs `--min-mapq 0`. (ii) LR §Methods "Training dataset": candidates outside GIAB high-confidence regions are removed, and the PacBio MAS model was trained on "the Kinnex full-length RNA dataset (HG002a:MAS:cDNA) from PacBio", so HG002 is not held out for the caller. (iii) LR §Methods "Benchmarking": accuracy is measured only where coverage ≥ 20 inside GIAB high-confidence regions. (iv) C3R §Results "Performance analysis by genomic context": "Low complexity" and "Segmental duplications" have the lowest F1, and longcallR averages about 66% over PacBio's difficult strata (CDS 91.74%, homopolymers 37.07%). | The headline accuracy does not carry over to SD copies at depth 6. Report caller QC on the teachers (§2 F15), and treat HG002 as caller-training data. |
| 3 | Disjoint-allele vote with assign-or-abstain | DM §Methods "Read assignment using PSVs" has the same beats-every-rival structure: the best location b must satisfy P_bi ≥ P_ib for every i, else "we keep the original alignment of the read and assign it a mapping quality of 0". DV §Methods: "If two or more equally good paths exist for a read, we do not assign the read to any path". PH §1: "haplotype mapping uncertainty in the presence of identical haplotype segments". DS §3: the field's defaults are fractional; rescue is "equivalent to a single iteration of the expectation maximization algorithm". A shared allele separating nothing matches DM's constant term for unreliable PSVs and PSC's f ≥ 0.95 rule. | "One separating column is enough" is the most permissive evidence rule in the vault. DM caps each PSV's weight (ŝ = 10⁻³, "to reduce impact of a single PSV on read mapping"). CH (Abstract) assumes "a sizable number of paralog-specific variants". SDA's attraction edges need reads that span ≥ 2 PSVs. samplecopy (r1333) needed ≥ 2 decisive sites. The no-counter-vote veto makes one erased column decisive: on the SMN2 smoke, one AF-0.14 het produced 688 placements. | The shape of the rule is well supported, but its precision rests on single-column genotypes. Report the outcomes of `vc_cols == 1` placements separately (§2 F1, F10). |
| 4 | SNVs only, indels ignored | LR §Results: Clair3-RNA's indel precision is 82.7% (HG004 MAS cDNA) and 60.5% (ONT cDNA), "much lower than the precision of SNP calling". C3R §Discussion: indel accuracy "remains challenging". DM drops low-complexity PSVs, the homopolymer class. | DM §Methods "identifying PSVs": a PSV is any differing substring between unique anchors, so indels are kept unless low-complexity. PP23 §Results: haplotypes are "phased variants (SNPs or indels)". PP25 §Results: a common "PMS2CL-specific" duplication sits on 72.8% of PMS2CL and 31.1% of PMS2 haplotypes, so indel PSVs are polymorphic too. Memory (corrector): AS charges about 8 per indel against 5 per mismatch; HG002's SMN and THOC3 copies differ from CHM13 by a 1-bp insertion; substitutions alone are right on SMN+THOC3 (588/605) but a coin flip elsewhere (17 fixes, 18 breaks). | Yes. Indels create many real ties, so the SNV-only rule both abstains (where copies differ only by indels) and overrides indel-driven AS preferences (§2 F8). DM's high-complexity-indel filter is a ready middle path. |
| 5 | RNA-editing, Dense and low-depth sites treated as unknown | LR §Introduction: editing events "may be mistaken as SNPs". C3R §Results: tagging REDIportal sites raises SNP F1. IC §Discussion: "RNA editing may generate transcripts that have substitutions but originate from the same copy, leading to overestimating the true copy number", which supports excluding A>G from flags. Depth 6 is LR's own candidate floor. Using "unknown" rather than a reference stand-in is more conservative than DM's prior p(0/0) = 0.95. | Dense does not detect collapsed copies at O2's identities. SDA §Results resolves collapses below 0.5% divergence, and O2's tied copies are ≥ 99.5% identical: ≤ 0.5 candidates per 100 bp against Dense's ≥ 5 (inf.). Depth 6 is below every vault benchmark (PP25 §Results, 10× per haplotype; LR, ≥ 20; C3R, zygosity errors fall as depth rises). RnaEdit labels flip between identical runs (ledger, 4/41), and on the smoke 21 of S19SMN2.0's 26 records were RnaEdit. | Keep "unknown". Drop the claim that Dense marks collapsed copies: it will mark artifacts (§2 F2). Expect run-to-run flips (T10). |
| 6 | Weak hets (SMN2, AF 0.14) | Reasons to doubt them: C3R §Methods "Variant zygosity flipping" uses VAF_HET 0.5 with tolerance ε 0.2, so AF 0.14 lies outside the het band (0.3-0.7) and within ε of VAF_Artifact 0.08. LR's candidate floor is alt AF ≥ 0.1 with ≥ 2 alt reads, so 0.14 is barely a candidate. Row 1's ≈ 2% error predicts exactly such AFs (§2.0). | Reasons to believe them: C3R §Results says zygosity flipping is "common in RNA data", through allelic imbalance and isoforms. PP23's "c" haplotypes and DM §Results, where 16.5% of HG001's two-copy PSVs are (0/0, 0/1), make real hets at PSVs common. No vault method separates leakage from allele-specific expression at a single site. | P1 risks placements created by erasure (A1, K1). P3 risks false homozygotes from allele-specific expression (§2 F4). T6's choice should also score K1 (§2 F16). |
| 7 | Foreign-allele missing-copy flag | DM §Methods "DuploMap algorithm": reads discordant at ≥ 5 PSVs get low MAPQ, as possibly from "missing duplicated sequences in the reference genome". In SDA, linking PSVs needs reads that carry both. PP25 Fig. 1: "an extra haplotype … indicates a copy number gain". PSC §Results: OTOA's aggregate copy number is 6 where the reference has 4. SDA: 51.5% of resolved collapses "diverged". | IC §Discussion: "no functionality … that would distinguish allele-specific variants from distinct gene copies". PP23 phased haplotypes into alleles only with pedigrees. SDA tells PSVs from alleles by single-copy depth, which RNA lacks. DV §Results: minority "haplotypes" were mixes "well supported by 'chimeric' reads". DM §Discussion: "does not currently account for missing sequences or copy number changes", so a copy absorbed into a reference copy is never foreign. | The flag means "unexplained haplotype". In a diploid that includes the untaught allele of a present copy (§2 F5), and it misses absorbed copies (§2 F2). C1 is weak evidence in either direction. |
| 8 | No iteration or EM | DS §3: EM and rescue are fractional. Thesis: r579 (EM gained nothing), r1351-1359 (rescue/EM harmful on real reads), r1345/r1346 (echo). | DM §Methods iterates placement and PSV reliability "until the read assignments do not change"; it is DM's only teacher-cleaning step. PSC estimates reliability by EM across samples. | No change. But O2-VC has no way to correct its ≈ 2% wrong teachers. The middle path is one non-fractional teacher check: DM's discordance demotion, applied to teachers. |
| 9 | No facility location | CZ §1: reads are "assigned (uniquely) … such that the total weight of the assignment is maximized", with no abstention. Its 0.19 bound is on the total weight, not on any read's correctness (CZ note §Caveats). Among AS-tied reads the weights tie, so the optimum is indifferent (inf.). Thesis: facility location was 47-59% precise. | Nothing in the vault. | No. |
| 10 | Reference-only null arm (N3) | PSC §Results: when population-unreliable PSVs are present, GATK and Freebayes precision fall from > 0.96 to 0.56/0.59 and recall from 0.52/0.55 to 0.25/0.29. That is the all-reference-PSV against reliable-PSV contrast. DM §Introduction: unreliable PSVs reduce "sensitivity and specificity". r1348: an abstaining metric needs a null arm. | The literature predicts that the sample calls mainly cut N3's breaks (precision). A3 instead tests fixes, the less-supported direction. N3 shares VC's known/unknown gating, so it isolates what the genotypes contain, not the gating. | Report N3's and VC's breaks next to A3 (reported only). |
| 11 | A hand-picked showcase plus a held-out donor | PP25 (259 people, 36 trios, 47 HPRC assemblies), PSC (2,504 + 698 genomes, MLPA, trios) and DM (genome-wide simulation, GIAB) all validate far beyond one locus, so holding out a donor is the right instinct. | Metric trap 09-18: hand-picked sets overstate F by about 0.27. GTF2IRD2 was picked after the corrector (19/0) and samplecopy (31/0) had succeeded on it. Panel E is donorsim output: reads "split evenly over its haplotypes' orthologs", no editing, and no donor copy without a CHM13 ortholog. Its PRAMEF/TRIM43B failure is already known (r1346). DM's real-data gains are checked "only indirectly" (DM note §Caveats). | Panel A shows the method can work, not how accurate it is. Panel E is optimistic for allele-specific expression, editing and missing copies (§2 F13). |

### 2. Pre-registered expected failures

#### 2.0 Conventions and derived quantities

**Fields:**
- `calls.tsv`: `vc vc_set vc_cols n3 n3_set flag outside margin mapq`.
- `<family>.columns.tsv`: `kind` ∈ `separating_call`, `separating_ref`, `erased_ref`, `foreign`, `editing`. Per copy it gives copy,
  chrom, pos, ref, genotype and depth, with `.` meaning unknown.
- `<family>.teachers.tsv`: FILTER counts per copy and `teacher_flag_reads`.
- `<family>.flags.tsv`; `vcf/<copy>.vcf` (FILTER, and FORMAT GT, DP, AF, PS).
- `<family>.case.tsv`: `truth_class`, `truth`, `vc_outcome`.
- `<family>.summary.md`: the `vc_cols == 1` count, the evidence class of each fix, and the placements that depend on erased columns.
- `score.<panel>[.<focus>].json`: arms, K1, K2, flags and bars; `run.json`.

**Leakage arithmetic (inf.).**
- Suppose a fraction ε of rival B's confident reads lands on copy A, and B has r times A's teachers at a site. The leaked allele's AF
  on A is then εr/(1+εr).
- It passes LR's floor of 0.10 at εr ≥ 0.111, and the preset's 0.15 (ledger) at εr ≥ 0.176. Two alt reads need ε·n_B ≥ 2.
- With ε = 0.02 (DM 1.9%, r1330 2.1%): r ≥ 5.6 for 0.10, r ≥ 8.8 for 0.15, and n_B ≥ 100.
- SMN2's AF 0.14 at DP 79 (≈ 11 alt reads) fits ε = 0.02 with n_B ≈ 550 and r ≈ 8.

**Small-n A1.** With fewer than 50 placements, wrong/placed ≤ 0.02 means zero wrong. GTF2IRD2 had about 35 decidable tied reads in
the corrector's domain, and minimap2 was wrong on 23 of them.

**Bars:**
- A1: wrong placements ≤ 2% of VC's placements. A2: ≥ 10 fixes. A3: VC fixes more reads than N3.
- K1: ≤ 2% of donor-tied (`own_tie`) reads placed. K2: flagged reads with a present copy ≤ 2%.
- C1: at least one flag holds ≥ 2 `no_counterpart` reads.
- E1: wrong ≤ 2% of placements. E2: fixes > breaks.

#### F1. Teacher contamination by a referenced rival (cross-copy leakage)
| | |
|---|---|
| (a) mechanism | About 2% of rival B's confident reads land on A with margin ≥ 10. On a copy expressed 5-9 times less than its rival, they push B's allele past the caller's floor at reference PSVs. The resulting record (usually a LowQual 0/1 whose ALT is B's reference base) makes A = {a, b}, and the column becomes `erased_ref`. Every tied read whose only counter-vote sat there is now placed: on SMN2, one site gave 688 placements. Such a placement is a break when the erased column was the read's correct vote and its remaining votes come from a biased genotype (F5). |
| (b) literature | DM §Results: 1.9% of MAPQ-≥ 30 reads demoted. PP25 §Results: misalignment despite high MAPQ. r1330: 2.1% off the donor's copy, in CNV SDs (GPAT2, SMN2, THOC3, FCGR2A, WASHC1, TBC1D3B). LR §Methods: AF ≥ 0.1 and ≥ 2 alt reads. Arithmetic in §2.0. |
| (c) signature | In `vcf/<copy>.vcf`, a 0/1 record (mostly LowQual) with FORMAT AF 0.10-0.25 on the copy with fewer teachers (`teachers.tsv`), whose ALT equals the paired copy's reference base. In `columns.tsv`, `erased_ref` at that position across most of the family's tied reads. In `summary.md`, "VC placements that depend on erased columns" > 0. In `case.tsv`, breaks over-represented among erased-dependent placements. |
| (d) prediction | Pushes A1 toward FAIL (one break fails GTF2IRD2, §2.0), raises K1, and raises E1 (donorsim's lognormal expression often gives r ≥ 5). Falsifiable: (i) in HG002, at least one non-dev pool family carries such a weak het on its lower-expressed copy; (ii) wherever ≥ 5 breaks exist, the break rate among erased-dependent placements is at least 3 times the rate among other VC placements. |
| (e) likelihood | High: it already happened on SMN2, and ε ≈ 2% is measured twice. |
| (f) mitigation | Partly, through T6's frozen policy: P3 restores the counter-vote, but only below AF 0.15. Follow-up: let a het erase a column only if the copy's own teachers carry the weak allele without the rival's other alleles (a teacher-side linkage test, computed already for `teacher_flag_reads`). Alternatively, demote teachers that disagree with their copy at ≥ 2 separating columns, once, without iterating (DM's discordance step). |

#### F2. Teacher contamination by an unreferenced copy (absorption)
| | |
|---|---|
| (a) mechanism | A donor copy C* with no CHM13 counterpart teaches its nearest reference copy A whenever its reads beat every rival there by ≥ 10. C*'s private alleles then enter A's genotype (0/1 at an AF equal to C*'s share, phased as A's "second haplotype" when A is otherwise homozygous). C*'s tied reads therefore carry no foreign allele, are never flagged, and are placed on A or abstain. The scorer cannot see this: a placed `no_counterpart` read is not decidable, so it is in neither the numerator nor the denominator of A1 or E1 (metric trap 09-22). Recent duplicates of a reference copy, the commonest extra copies, are absorbed by construction (inf.). |
| (b) literature | DM §Discussion: "does not currently account for missing sequences or copy number changes". PP25 §Discussion: Paraphase "is not influenced by the CN difference between an individual and the reference", the property O2-VC lacks. PP25 §Results: three CN-variable groups look like false duplications "due to missing genes in the CHM13 T2T assembly". PSC §Results: 15 loci above the reference copy number in more than half the samples (OTOA, 6 vs 4). SDA: 51.5% of resolved collapses diverged. r1343 and r1346: BRANCH/REMOVED worlds and PRAMEF. |
| (c) signature | `case.tsv` rows with `truth_class = no_counterpart` and a non-empty `vc`. In the VCF, ≥ 2 PASS 0/1 records sharing one PS on a copy, at AF 0.15-0.35, with ALTs that are no family copy's reference base (a collapsed haplotype). `teacher_flag_reads` > 0 in `teachers.tsv` only when C*'s share is below the caller's floor. Dense stays at 0-2 records per copy and away from PSVs. |
| (d) prediction | Pushes C1 toward FAIL; A1 and E1 are blind to it by construction. Falsifiable: in HG002, VC places a non-zero number of `no_counterpart` tied reads in at least one pool family. Panel E cannot show F2, because donorsim simulates only orthologs of CHM13 copies. |
| (e) likelihood | High that it occurs; medium for its effect on C1. |
| (f) mitigation | None: it is the spec's "Known limit". Follow-up: a Paraphase-style haplotype count from longcallR's phase sets, where a copy carrying a phased minor haplotype of private alleles becomes unknown and is flagged on the teacher side. Report `no_counterpart` placements. |

#### F3. The donor lacks a reference copy, or one donor copy sits under two reference copies (the PRAMEF case)
| | |
|---|---|
| (a) mechanism | The donor's copy D teaches two reference copies A and A′, each read going to the one it resembles more. If D is homogeneous, A and A′ learn the same alleles, the columns become `erased_ref`, and VC abstains, as the spec expects. If D's two haplotypes resemble different reference copies, hap1 teaches A and hap2 teaches A′. One donor copy then yields two disjoint genotypes, a "phantom separation", and D's tied reads are placed by haplotype: half of them land on the copy that is not D's lifted counterpart. |
| (b) literature | r1346: 31 of the corrector's 40 breaks are PRAMEF/TRIM43B; HG00733 lacks PRAMEF8, and the donor's PRAMEF7 reads fill PRAMEF8. r1336 (Y ampliconic genes): two RBMY1 groups share one label, alleles taken for copies. PP23, Table 1 note: the only copy-number error was "two of the three haplotypes being identical in sequence". PP25 §Results: in 22 GRCh38 groups, > 95% of people have copy number 2. PAL: HG002 carries 45 TSPY copies against a median of 32. |
| (c) signature | `separating_call` columns in `columns.tsv` where both copies' genotypes differ from their own reference base, each holding the other's base or a donor allele (complementary swaps). In `case.tsv`, breaks in both directions inside one family. |
| (d) prediction | Pushes E1 toward FAIL, since E contains PRAMEF/TRIM43B; risk to A1 in HG002 in r1330's CNV families. Falsifiable: in E, VC abstains on ≥ 90% of PRAMEF/TRIM43B tied reads and makes ≤ 5 breaks there (corrector v1: 31). |
| (e) likelihood | Medium. |
| (f) mitigation | Erasure covers the homogeneous case; nothing covers teaching split by haplotype. Follow-up: when two copies' genotypes are swapped relative to the reference at ≥ 2 columns, treat them as one copy and abstain. |

#### F4. Allele-specific or monoallelic expression
| | |
|---|---|
| (a) mechanism | Take a donor het at copy A, a/b, where b is rival B's base. If b is expressed below the caller's floor, A is called {a}; the column separates A = {a} from B = {b}, and A's minority b-reads vote for B. A read with no other vote for A is placed on B (a break), and own_tie reads from that haplotype get placed (K1). Fully monoallelic a is harmless, since no b-reads exist. Fully monoallelic b makes A = {b} = B, which erases the column (also harmless). P3 (weak alt → {ref}) creates this case on purpose. |
| (b) literature | C3R §Results "Zygosity flipping analysis": "also common in RNA data"; ignoring zygosity raised F1 by 2.43 points (ONT) and 2.88 (PacBio), and removed 4,334 false positives and 4,331 false negatives in BCM Iso-Seq HG004. C3R §Methods: a het "appears as a homozygous reference or a homozygous variant due to allelic imbalance". LR §Introduction: "skewed allele fraction caused by allele-specific expression". |
| (c) signature | A break whose deciding column in `columns.tsv` (often `vc_cols == 1`) has the true copy homozygous, while that copy's teacher BAM carries the read's base in 2-15% of reads. This needs an audit pileup; HomRef records expose FORMAT AF. None in E, because donorsim splits each copy's reads evenly over its haplotypes. |
| (d) prediction | Raises A1 and K1 in HG002; E is unaffected, so part of any precision gap between HG002 and E is allele-specific expression. Falsifiable: if VC breaks at least once in HG002, at least one break is decided at a column where the truth copy's teachers carry the read's base at AF < 0.15. |
| (e) likelihood | Medium: it needs a het at a PSV, strong imbalance and a sole vote. |
| (f) mitigation | None: the GT is taken whatever the AF. Follow-up: no single-site rule separates allele-specific expression from leakage (Row 6). Use linkage among the teachers instead. If the weak allele co-occurs with the rival's alleles at other columns on the same teacher reads, it is leakage; if it is independent of them, it is a real het. |

#### F5. Teachers selected by their alleles (reference bias of the teacher set; intrinsic to the design)
| | |
|---|---|
| (a) mechanism | Whether a read teaches depends on its alleles. Copy A's reads that carry B's allele at a PSV lose about two mismatches' worth of margin and tend to tie, so A's teachers under-sample exactly the haplotype (or conversion tract) of A that resembles B. Three effects follow. (i) False homozygotes at donor-het PSVs: F4 without any allele-specific expression. (ii) The private SNVs of A's B-like haplotype are missing from A's genotype, so that haplotype's tied reads carry linked "foreign" alleles and are flagged although A is present (K2). (iii) In a conversion tract, A's teacher depth falls below 6, so A is unknown there and VC abstains. |
| (b) literature | DV §Methods: "rtoa and ator differ owing to reference bias: If a read comes from a haplotype with a true alternate allele, it may systematically align with the reference allele incorrectly". PP25 §Methods realigns every read of a group to one archetype, so no read is chosen by where it was placed. GU and DSH measure margins against panels of assembled haplotypes, so copies are compared against real copies, not one reference. DM re-places reads iteratively. r1329: 73% of HG002 tied reads also tie in the donor. r1334. |
| (c) signature | In `teachers.tsv`, a copy's teacher count far below its share of tied reads. In `flags.tsv`, flags whose reads have `truth_class` unique or own_tie on one copy (the K2 numerator) and whose alleles all sit on one haplotype. Separating columns where about half of a copy's tied reads carry the rival's base: a haplotype, not sequencing error. |
| (d) prediction | Raises K2 in HG002 and E, and lowers A2 through effect (iii). Falsifiable: if any flag holds ≥ 3 reads of a present copy, its alleles will be heterozygous in the donor's copy (checked against the donor assembly on the truth side), not errors. |
| (e) likelihood | Medium-high: the mechanism is structural; its size is unknown. |
| (f) mitigation | None. The disjoint rule confines effect (i) to reads without another vote. Follow-up: genotype each copy from all of its family's reads realigned to one archetype (Paraphase), and use the teachers only to name the haplotypes. |

#### F6. RNA editing
| | |
|---|---|
| (a) mechanism | A-to-I editing (A>G on the transcript strand) at an A/G reference PSV makes the A-copy's edited reads look like the G-copy's. It is harmless where longcallR labels the site RnaEdit (unknown) or calls it 0/1 (erased). It is harmful where the editing level in the teachers stays below the floor, so the A-copy = {A}, and an edited tied read has that column as its only vote. Weak T>C and A>G records flip FILTER between identical runs, so a column's known/unknown status, and the placements that depend on it, can flip too. |
| (b) literature | C3R §Results: 8,049 (Iso-Seq) and 15,659 (MAS-Seq) high-quality editing sites; about half of the PacBio sites are intronic; REDIportal tagging raises F1. LR §Introduction. IC §Discussion: editing inflates copy counts. Ledger: 21 of 26 records RnaEdit on S19SMN2.0, and flips in 4 of 41 runs. |
| (c) signature | RnaEdit counts in `teachers.tsv`; `editing` rows in `columns.tsv`; breaks decided by an A↔G column (in the transcript frame) whose A-copy genotype is {A}; differences at RnaEdit/LowQual sites between T10's five re-calls. |
| (d) prediction | At most 1 HG002 break attributable to editing. E has none, since no editing is simulated, so E is optimistic. T10: if S19GTF2IRD2 has any weak A>G or T>C record, at least one VC placement or flag changes across the five re-calls. |
| (e) likelihood | Low for errors; high for power lost at A/G PSVs; medium for run-to-run instability. |
| (f) mitigation | Treating RnaEdit as unknown and excluding editing from flags covers the main cases; editing below the caller's floor is not covered. Follow-up: an A↔G column (transcript frame) may not be a read's only vote. |

#### F7. Gene conversion and mosaic copies
| | |
|---|---|
| (a) mechanism | A donor copy converted toward its paralog. A homozygous tract gives it the paralog's base, the column becomes `erased_ref`, and VC abstains or decides elsewhere while N3 votes wrong. An allelic (het) tract is also erased. Mosaic reads vote for both copies and abstain. VC errs here only through F5, when the converted haplotype is the untaught one. |
| (b) literature | PP25 §Results: more than 75% of African PMS2/PMS2CL haplotypes are partly or fully converted, and 36 trios show 4 de novo conversions, 2 of them non-allelic. PP23 §Results: "c" haplotypes. PSC §Discussion: its constant-copy-number assumption breaks under "gene conversion events and hybrid alleles"; VC decides per column and makes no such assumption. DM §Discussion: "gene conversion is well known to result in overlap between PSVs and variants". PAL: 99.94% and 99.87% homogenization. |
| (c) signature | Runs of adjacent `erased_ref` columns on one copy, homozygous for the paralog's base. Fixes whose evidence class in `summary.md` is "reference PSV corrected by the sample". `n3` wrong on the same reads where `vc` is right or empty. |
| (d) prediction | Pushes A3 toward PASS. In E, where HG00733's conversions are real, N3 makes at least twice as many breaks as VC. |
| (e) likelihood | High that conversion occurs; low that it causes VC errors. |
| (f) mitigation | Mitigated by design (`erased_ref`); what remains goes through F5. |

#### F8. Real indel PSVs lost (SNV-only)
| | |
|---|---|
| (a) mechanism | Many real-read ties are made by an indel: AS charges about 8 per 1-bp indel and 5 per mismatch. VC sees only SNVs. Where the copies differ only by indels inside the common window, it abstains. Elsewhere it overrides an AS preference that came from an indel, on SNV evidence. That is right when the indel is a donor-vs-reference difference (SMN, THOC3) and wrong when it is a real PSV the donor carries. N4, which uses left-aligned indels, can decide reads VC cannot. |
| (b) literature | LR §Results: SNVs only; Clair3-RNA's indel precision is 82.7% and 60.5%. C3R §Results "Indel performance": PacBio indel F1 averages 87.46%. DM §Methods: PSVs of any length, with low-complexity ones filtered out. PP23: "SNPs or indels". PP25: the PMS2CL duplication on 72.8% and 31.1% of haplotypes. Memory (corrector): 588/605 right on SMN+THOC3, 17 fixes and 18 breaks elsewhere. |
| (c) signature | `reads.tsv` rows whose tied copies differ in `indels` but not in `substitutions`, giving an empty `vc` and `vc_cols` 0. N4 fixes, in `score.json`, on reads VC abstains on. Among VC placements that disagree with `minimap2`, the placed copy's alignment carries more indel bases, which follows from how the tie arises. VC breaks concentrated in families whose CHM13 copies differ by a non-homopolymer indel inside the common window. |
| (d) prediction | Lowers A2. On GTF2IRD2, VC fixes fewer reads than N4's 19; falsifiable: VC makes ≤ 15 fixes there. The sign of E2 is unaffected. A1 and E1 rise only in families with indel PSVs. |
| (e) likelihood | High for lost power; medium for errors. |
| (f) mitigation | Not mitigated: it is known gap 1. Follow-up: DM-style high-complexity indel PSVs as columns under the same disjoint rule, with their own prereg. |

#### F9. Low expression of a copy
| | |
|---|---|
| (a) mechanism | A tied copy whose teacher depth is below 6 at a column is unknown there, and an unknown copy can never be beaten. Every read whose tie set contains it abstains. The foreign-allele test needs every tied copy known, so the same copy also blocks flags. Many-copy families (NBPF1, NPIP, SPDYE1) are blocked by their least-expressed tied copy, which is also where leakage AF peaks (F1). |
| (b) literature | LR §Methods: a depth floor of 6, and a benchmark at ≥ 20. PP25 §Results: 10× per haplotype. C3R §Results: DP ≥ 8 instead of ≥ 4 adds 4.19 F1 points (ONT). IC §Results: DAZ reaches > 90% recall only at about 410×. DSH §Results: 55.6% of reads individually assignable, even with 169 haplotypes. r1334. |
| (c) signature | An empty `vc` in `calls.tsv` whose `vc_set` contains a copy with genotype `.` at every column of the read in `columns.tsv`. That copy's count in `teachers.tsv`; `failed_copies` in `run.json`. |
| (d) prediction | Pushes A2 toward FAIL: on GTF2IRD2 I guess 3-9 fixes, and on NBPF1 A2 needs 10 of about 11 minimap2-wrong reads. Pushes C1 toward FAIL: SPDYE1 is blocked, and NPIPB2 needs both .5 and .7 known. Some families in E place no read at all. |
| (e) likelihood | High. |
| (f) mitigation | Correctly conservative; there is no fix inside the rule. Follow-up: before scoring, report per family the tied reads whose tie set is fully known at ≥ 1 column, so an A2 miss is classed as a power failure, not an error. |

#### F10. 5′-truncated or short reads
| | |
|---|---|
| (a) mechanism | A column must be aligned in every tied alignment, and separating sites cluster in a few exons (for SMN, exons 7-8 and downstream). Truncated or short reads, or reads whose tied alignments clip or splice differently, cover no separating column and abstain. A read that covers a single column becomes a `vc_cols == 1` placement, resting on one genotype and one base. |
| (b) literature | GR §Discussion: multi-mapping reads "covered transcripts only partially and thus could not capture the small difference between gene copies". IC §Discussion: "impossible to separate transcripts originating from copies with identical exonic sequences". PP23 §Results, for exons 1-6: "read alignments are often ambiguous in this region, even for long reads". C3R: F1 drops from 90.31% to 83.98% within 2-16 bp of splice boundaries (Iso-Seq). |
| (c) signature | Abstentions dominated by `vc_cols = 0`; VC's placement rate rises with `ms` (`reads.tsv`); the `vc_cols == 1` share in `summary.md`. |
| (d) prediction | Lowers A2; precision is unaffected apart from single-column errors. Falsifiable: in pooled HG002, VC's placement rate among tied reads in the top `ms` tertile is at least twice the rate in the bottom tertile. |
| (e) likelihood | High for lost power; low for errors. |
| (f) mitigation | Correct abstention. The exon-only footprint is RNA's built-in ceiling compared with DNA methods. |

#### F11. Near-identical copies (DAZ)
| | |
|---|---|
| (a) mechanism | DAZ copies are near-identical, and DAZ's 2.4-kb repeat units lie inside one gene, so tied alignments can pair a read base with non-homologous repeat positions (spurious columns). Most DAZ reads tie on all copies, so teachers are scarce, copies are unknown, and VC abstains. The Y is haploid: a 0/1 or 1/2 record on a non-PAR Y copy cannot be a germline het, so every het on DAZ marks a copy mixture, a collapse or a repeat misalignment. A119b's Y is not CHM13's (CHM13 v2.0's chrY is HG002's T2T-Y), so A119b's AZFc copies can differ in number and structure. |
| (b) literature | PAL §Results: palindrome arms 99.94% identical, arrays 99.87%; DAZ has 2-19 repeats. RH: the T2T-Y is "from the HG002 genome", and the DAZ "entire repeat array, consisting of 2.4 kb composite units … falls within one gene annotation". GR §Discussion: "the DAZ gene includes almost identical repeats that could cause alignment issues". PP25: Y groups indistinguishable by sequence alone. IC: DAZ needs about 410× for > 90% recall. r1336: 393 of 579 DAZ tied reads tie on all four copies, and DAZ2 and DAZ3 differ at one expressed site. |
| (c) signature | Panel B's `score` (`--mode none`): placed, abstained, flagged, and `vc_cols == 1`. Deciding positions in `columns.tsv`; 0/1 and 1/2 records in the DAZ copies' VCFs; `flags.tsv`. |
| (d) prediction | B is descriptive. ≥ 95% of DAZ tied reads abstain. Every VC placement has `vc_cols == 1`, and the deciding columns fall at ≤ 2 sites (expected: r1336's DAZ2/DAZ3 site). N3 places at least as many reads as VC. At least one DAZ copy has a het record. A DAZ flag whose alleles lie inside the repeat units is an artifact. |
| (e) likelihood | High for the abstention; medium for a het record. |
| (f) mitigation | Mitigated by abstaining. Haploidy is a free control for contamination: report het records on the non-PAR X and Y copies of male samples. |

#### F12. Reference errors or minority alleles in CHM13 at PSVs
| | |
|---|---|
| (a) mechanism | A CHM13 base error at a PSV, or CHM13 carrying a minority allele: the copy's teachers carry the population base, the copy is called 1/1, and the column is erased (VC right, N3 wrong). A CHM13 false duplication: reads tie exactly, no teacher exists, the copy is unknown, and VC abstains. A copy missing from CHM13 is F2. |
| (b) literature | DM §Methods: "some PSVs may be assembly artifacts". PP25 §Results: 22 GRCh38 false duplications, most already flagged by CHM13, and three groups where CHM13 lacks a gene. PP23 §Results: GRCh38 carries the minority SMN1 haplotype S1-2. Yoo_2025 §Key Points: read-depth conclusions about copy number overturned by complete assemblies. |
| (c) signature | 1/1 records whose ALT equals the paralog's reference base at the same CHM13 site in HG002, HG00733 and A119b wherever they are covered: a non-reference homozygote shared across independent donors. |
| (d) prediction | Helps A3 and E1 and harms no bar; such sites will exist in at least one pool family. |
| (e) likelihood | Low that it causes harm. |
| (f) mitigation | Mitigated (`erased_ref`). |

#### F13. Simulation artifacts in HG00733 (panel E)
| | |
|---|---|
| (a) mechanism | donorsim splits each copy's reads evenly over its haplotypes (no allele-specific expression). It simulates no editing and no library chimeras. It drops copies whose orthologs overlap on a haplotype, and exons that change by > 50 bp. It never simulates a donor copy that lacks a CHM13 ortholog. So E cannot show F2, F4 or F6, and shows only part of F3. Its truth is exact, so a placed read from donor-identical copies is right or wrong by chance. E is also not pristine: r1346 ran the corrector on it, and the PRAMEF/TRIM43B failure shaped the spec's "Known limit". |
| (b) literature | The `donorsim.py` docstring, read as the method's definition. PP25 §Results: the simulation floor (10 kb, 10×, 0.05%) shows how simulations set favourable conditions. DM note §Caveats: real-data gains checked "only indirectly". CH note §Caveats: simulation only, with the sizable-PSV assumption true by construction. DV note §Caveats: accuracy simulated. |
| (c) signature | E's 0/1 records sit at AF ≈ 0.4-0.6, and RnaEdit records appear only at simulated errors. E has fewer erased-dependent placements per placement than HG002. Every flag in E counts toward K2. |
| (d) prediction | E's wrong fraction is below HG002's. Falsifiable: if both have ≥ 100 VC placements, HG002's wrong fraction exceeds E's. E says nothing about C1. |
| (e) likelihood | High. |
| (f) mitigation | Not mitigated. Follow-up: re-simulate with beta-binomial allelic ratios, editing at REDIportal sites, and donor-extra copies. |

#### F14. Truth-label artifacts (liftover)
| | |
|---|---|
| (a) mechanism | HG002's truth lifts the donor's copies to CHM13 through a chain. In SDs and CNVs a chain can fail: a gap gives `no_counterpart` although a CHM13 copy exists, so C1 counts a spurious hit and K2 is understated. A chain can also pair a donor copy with the paralogous CHM13 copy: the truth is swapped, so a right VC placement scores as a break and a wrong primary as right. chrM has no chain and is excluded. CHM13 v2.0's non-PAR Y is HG002's own, so any `no_counterpart` there is a lift artifact. |
| (b) literature | DSH note §Method in brief: paralogs are named by phylogeny, "not … synteny alone", because interlocus gene conversion breaks positional identity. GU: minigraph could not group TBC1D3 copies. PP25 §Results: most disagreements with HPRC assemblies were traced to assembly errors. RH: the T2T-Y comes from HG002. r1330: 258k chrM `no_counterpart` reads. r1348: the FCGR2A lift was never verified. |
| (c) signature | Breaks whose read aligns with fewer edits to the donor locus of VC's copy than to the truth copy's (an audit re-alignment, on the truth side). C1 flags whose donor locus aligns to a CHM13 copy (C2). Breaks in both directions inside one family. |
| (d) prediction | Spurious breaks in A1, spurious passes in C1, and an understated K2. Falsifiable: at least one of any HG002 panel-A breaks turns out to be a lift artifact on audit. |
| (e) likelihood | Medium. |
| (f) mitigation | Partly mitigated by C2 and the chrM exclusion. Follow-up: before reporting, re-align every break, and the reads of every counted flag, to both donor haplotypes. |

#### F15. Caller false negatives and off-label behaviour
| | |
|---|---|
| (a) mechanism | Where a copy's teachers carry a non-reference allele but longcallR writes HomRef, or the allele misses the 0.10-AF / 2-read floor, the copy is "known = {reference base}". That is an N3 vote, so caller misses in SDs re-import the reference-PSV errors VC exists to remove. Other quirks: non-determinism (4/41); unphased hets written as LowQual, whose GT is still used; approximate phasing; `--min-mapq 0`, which the authors never benchmarked; and HG002 Kinnex in the training set. |
| (b) literature | LR §Results: recall 91.0-95.3%, and only in high-confidence regions; it "may miss true SNPs due to phasing difficulties". LR §Methods: the phasing "may miss the global minimum"; trained on HG002a:MAS:cDNA, in high-confidence regions only. C3R §Results: SDs among the lowest-F1 strata, longcallR averaging about 66%. Memory (corrector): substitutions alone go 17/18 outside SMN/THOC3. |
| (c) signature | HomRef records with FORMAT AF ≥ 0.10 on a copy; `separating_ref` columns at those sites; breaks decided there; differences between T10's re-calls. |
| (d) prediction | Raises A1, E1 and K1. Falsifiable: at least one HG002 pool family has a HomRef record with AF ≥ 0.3 at a reference PSV. |
| (e) likelihood | Medium. |
| (f) mitigation | Not mitigated. Follow-up: treat a HomRef record whose AF reaches the caller's own candidate floor as unknown; this uses the caller's constant, not a new one. |

#### F16. Scoring artifacts
| | |
|---|---|
| (a) mechanism | (1) With fewer than 50 placements, A1 means zero wrong. (2) On `outside = 1` reads N1 is wrong by construction, so every right VC placement counts as a fix; NPIPB2's 10 reads show such reads exist. (3) Placed `no_counterpart`, `poor_match` and `elsewhere` reads fall outside A1 and E1 (a universe intersection). (4) T6 picks the genotype policy by dev precision on decidable reads, breaking ties by fixes. own_tie placements (K1) are not in its criterion, and P1 maximizes the placements created by erased counter-votes (688 on SMN2). (5) If SMN2 and MSTRG3150 count among "HG002's pool families", data used for T6 enters K1, K2 and C. (6) K1's denominator, about 9.6k non-PAR own_tie reads (r1329), means about 190 placements fail it, and one SMN2-like family can supply them. |
| (b) literature | Metric traps 09-18 (hand-picked sets) and 09-22 (a universe intersection inflates the worst arm's precision); r1348 (null arm). |
| (c) signature | `arms.VC.placed` in `score.A.S19GTF2IRD2.json`; fixes split by `outside`; the truth classes of placed reads in `case.tsv`; K1 per policy in the T6 doc. |
| (d) prediction | A1 is a knife-edge, `outside = 1` inflates A2, and one family can decide K1. |
| (e) likelihood | Certain: these are arithmetic. |
| (f) mitigation | The reported-only counts in §2.18. |

#### 2.17 Panel verdicts: best guesses, to be checked and never tuned toward

P is my probability that every bar of the panel passes.

- **A, GTF2IRD2: FAIL**, most likely on A2. I guess 3-9 fixes: the teachers are too thin at depth ≥ 6 (r1334), and VC sees only SNVs where
  N4 made 19 fixes. A1 holds at 0 wrong, and A3 passes if VC places anything. P ≈ 0.15 with K1 and K2 included; 0.3 for A1-A3 alone.
- **A backup, PDPR: the best of the three**, P ≈ 0.35. In r1334, margin-taught labels gave 37/0, and minimap2 is wrong on 34 reads. A1
  is a zero-error knife-edge.
- **A backup, NBPF1: FAIL, underpowered.** A2 needs 10 of about 11 minimap2-wrong reads, and its many copies leave unknown copies that
  block placements.
- **K1 (HG002 pooled): PASS**, P ≈ 0.7. It fails only through erased-counter-vote placements in an SMN2-like family under P1 or P2.
- **K2 (HG002 pooled): PASS**, P ≈ 0.85. The SMN2 smoke flagged 2 of 1,655 tied reads, and the denominator is large.
- **B, DAZ: as expected.**
  - ≥ 95% of tied reads abstain.
  - At most a few `vc_cols == 1` placements, at the one DAZ2/DAZ3 site.
  - At least one het record on a DAZ copy.
- **C-human, C1: FAIL**, P(pass) ≈ 0.4. Absorption, the every-copy-known rule and the linkage rule each remove candidates. A pass would
  come from GTF2H2C or NPIPB2, and must clear C2.
- **E, HG00733: lean PASS**, P ≈ 0.5. E2 ≈ 0.9 and K2 ≈ 0.9; E1 ≈ 0.6, failing (if at all) in families with skewed expression and in
  PRAMEF/TRIM43B.
- **D and C-gorilla: not run** (the fast path), so this round has no one-genome no-harm check.

#### 2.18 Reported-only counts that make these predictions checkable (no new bars)
1. VC placements of `no_counterpart`, `poor_match` and `elsewhere` reads (F2, F16).
2. Per family, before scoring: the tied reads whose tie set is fully known at ≥ 1 column (F9).
3. VC's and N3's breaks next to A3 (Row 10).
4. Fixes split by `outside` (F16).
5. Per copy: 0/1 records with AF < 0.25, and whether each ALT is the paired copy's reference base (F1); HomRef records with AF ≥ 0.10
   (F15).
6. Het records on haploid copies in male samples (F11).
7. K1 on the dev families under each T6 policy (F16).
8. Right and wrong among erased-dependent placements, wherever a truth exists (F1).

### 3. Missing prior art and the novelty delta

#### 3.1 Side by side

| | DuploMap (DM) | Paraphase (PP23, PP25) | Parascopy (PSC) | O2-VC |
|---|---|---|---|---|
| Substrate | Long DNA reads (CLR, HiFi, ONT) in SDs > 5 kb at ≥ 97% | HiFi DNA, 160 GRCh38 groups > 10 kb at > 99% | Short-read WGS, 167 loci, cohorts | Long RNA reads (Kinnex FLNC), AS-tied reads of named CHM13 families |
| What it decides | Location and MAPQ of every read whose primary overlaps an SD | Haplotypes (= copies), copy number, phased variants, fusions | Aggregate and per-paralog copy number per sample; each PSV's f | Each AS-tied read's copy, an abstention set, or a flag |
| Evidence | Reference PSVs of any length (high-complexity), pair-HMM windows, unique-k-mer LCS | Every variant site of all reads realigned to one archetype | Read depth plus PSV allele counts | SNV columns common to all tied alignments |
| Whose variants | The sample's: per-copy PSV genotypes from reads with MAPQ ≥ t; reliable = posterior (0/0, 0/0) | The sample's own haplotypes | The population's: f across ≥ 50 samples (f ≥ 0.95) | The sample's: per-copy GT from margin-≥ 10 teachers via longcallR |
| Sample-private alleles | Not used (reference PSVs only) | They define the haplotypes | Not used | Can separate (hom-alt) or erase (het) |
| Decision rule | Likelihood; the best location beats each rival pairwise; per-PSV cap | Graph assembly, then read → haplotype (tie handling not described) | HMM, then EM, with quality scores | Combinatorial: assign c iff c beats every rival (≥ 1 vote, no counter-vote) |
| Abstention | MAPQ 0, original alignment kept | Groups of low diversity left inseparable; no-calls | Low quality or no per-paralog estimate | An explicit set of unbeaten copies |
| Certificate | None (a probability) | Haplotype plus a diagnostic base | None | Deciding columns listed per placement |
| Missing copy | Reads discordant at ≥ 5 PSVs demoted; copy-number changes "not accounted for" | An extra haplotype = a gain; immune to the reference's copy number | Aggregate copy number above the reference (OTOA) | ≥ 2 linked recurrent foreign alleles in tied reads |
| Iteration | Yes, until assignments stop changing | Not described | EM | None |
| Indels | Yes (high-complexity) | Yes | Not described | No |
| Validation | Simulations (precision 0.995), GIAB HG002 calls | 259 people, 36 trios, 47 HPRC assemblies | 3,202 1kGP genomes, MLPA, trios | One hand-picked HG002 family; one simulated donor (HG00733) |

#### 3.2 What overlaps, what O2-VC adds and what it lacks

**DuploMap (the closest analogue).**
- Overlap:
  - it post-processes minimap2;
  - it genotypes per-copy PSVs from the sample's confidently placed reads (DM's teachers);
  - only PSVs the sample confirms at both copies decide;
  - a location must beat every rival pairwise, else MAPQ 0;
  - reads discordant at ≥ 5 PSVs are demoted as possibly coming from a duplicated sequence the reference lacks.
- O2-VC adds:
  - RNA: spliced, exon-only, expression-skewed and edited reads;
  - a scope of AS-tied reads only;
  - a stock RNA caller in place of a built-in PSV genotyper;
  - sample-private alleles that can separate, and shared alleles that erase (DM only rates reference PSVs as reliable or not);
  - a hard certificate that lists the deciding columns;
  - a linked foreign-allele flag (≥ 2 recurrent alleles, ≥ 2 reads) in place of a per-read discordance count;
  - a reference-only null (N3) and a selective baseline (N2) on the same reads.
- O2-VC lacks:
  - indel and multi-base PSVs;
  - DM's iteration, which is its teacher-cleaning loop;
  - a per-PSV influence cap;
  - the unique-k-mer LCS pre-filter for reads that reach past the SD;
  - a probabilistic MAPQ;
  - a genome-wide evaluation;
  - any copy-number handling (DM lacks it too).

**Paraphase.**
- Overlap:
  - the evidence is the sample's own copies, not reference positions;
  - copies are treated as haplotypes;
  - it abstains where copies are as similar as alleles;
  - gene conversion is visible per haplotype.
- O2-VC adds:
  - RNA;
  - per-read placement onto a named reference copy. The names come from where minimap2 put the teachers. PP25 §Results: "It is often
    not easy to assign haplotypes of a paralog group to individual genes without prior knowledge", so Paraphase names haplotypes by
    unique flanks, curated PSVs or c.840;
  - a per-read certificate;
  - the tied-read scope.
- O2-VC lacks:
  - phasing: multi-site haplotypes from all reads realigned to one archetype, immune both to placement errors and to copy-number
    differences;
  - copy number from the haplotype count plus depth (in RNA, depth is expression);
  - detection of identical haplotypes, which needs depth (PP23's only copy-number error);
  - fusion calls;
  - population and trio validation.

**Parascopy.**
- Overlap:
  - a reference PSV is a hypothesis;
  - only reliable PSVs count;
  - no estimate is given where PSVs are unreliable (C4A/B: 7 of 50 reliable).
- O2-VC adds: one sample at a time, per-copy genotypes in place of population frequencies, read assignment, and RNA.
- O2-VC lacks:
  - population-level reliability, so it cannot tell a donor het from leakage (Row 6);
  - copy number from depth, which finds missing copies such as OTOA directly;
  - replicate and Mendelian validation.

#### 3.3 Novelty, stated honestly

**Not new.**
- Genotyping PSVs per copy from the sample's confident long reads, and deciding only with the PSVs the sample confirms (DM: DNA,
  2020).
- Abstaining when copies cannot be told apart (DM, PP25, DV).

**New, as far as the vault shows.**
- Doing it for long RNA reads.
- With an off-the-shelf RNA caller.
- With a threshold-free combinatorial certificate restricted to AS-tied reads.
- Letting sample-private alleles separate copies.
- With a linked foreign-allele flag.
- With an evaluation design that holds a reference-only null, the selective baseline and own_tie calibration.

The claim to make is "DuploMap's reliable-PSV principle carried to long-read RNA, with a certificate and an RNA-side O3 flag". It is
not a new principle. The vault holds no RNA analogue of DuploMap, but a targeted search must come before any "first" claim.

#### 3.4 Methods from the vault to adopt, compare against or cite

**Adopt:**
- DM's high-complexity indel PSVs (gap 1).
- DM's discordance demotion, applied to teachers once (F1, F2).
- A Paraphase-style haplotype count from longcallR's phase sets, to mark collapsed copies (gap 2).
- PSC-style "no estimate" reporting: each family's power ceiling, before scoring (F9).

**Compare against (truth side or reported only):**
- Paraphase run on HG002 HiFi DNA, if GTF2IRD2, PDPR or NBPF1 are among its 160 GRCh38 groups, to cross-check copy counts. It must
  never become an input.
- Clair3-RNA as a second caller, if it installs without Docker; C3R reports it ahead of longcallR in SDs.
- DSH's best-plus-one-mismatch rule needs the donor assembly as input (an anti-goal), so it can only be cited.

**Cite:**
- GU, the source of the 10-AS rule.
- PAL §Limitations: transcriptome and assembly from different individuals prevent "a direct one-to-one transcript-to-copy mapping";
  O2-VC's same-sample premise answers this.
- GR §Discussion (own work): multi-mapping reads left unresolved.
- Tomaszkiewicz_2023_YAG_Isoforms: asks for such tools.
- IC: alleles versus copies, and editing.
- DV: abstaining on ties, and reference bias.
- CZ and DS as foils.
- LR and C3R for the caller's limits.
- Tang_2020_FLAIR: FLAIR's reassignment keeps MAPQ ≥ 1, the same selective baseline as N2.

#### 3.5 What the vault lacks
- **Allele-specific mapping bias.** WASP-style correction, and RNA alignment to personal genomes or pangenomes. This is the theory
  behind F5, and the alternative the anti-goals rule out.
- **How common allele-specific and monoallelic expression is** in lymphoblastoid lines and testis, needed to size F4.
- **RNA editing** in SDs and in testis (REDIportal itself, ADAR).
- **Kinnex/MAS-seq library artifacts**: concatemer splitting, chimeras, RT template switching. These are sources of foreign alleles.
- **The minimap2 paper.** The definitions of AS, MAPQ and secondary hits (`-N`, `-p`) that "tied" rests on.
- **How accurate liftover and chains are in SDs**, for the truth side.
- **Gene-conversion rates** (Dumont 2015; Vollger 2023). The vault cites them only through DM and PP25.
- **DM's and Paraphase's supplements**: the discordance filter, the haplotype graph, and how ties are handled.


---

## PREREG_o2_forced_heldout_2026-10-23.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_forced_heldout_2026-10-23.md`)*

## Pre-registration: O2 as forced haplotype-to-copy matches, simulated gorilla held-out (P7), 2026-10-23

**Status:** committed before any P7 family is run with the O2 code. Goal milestone 3 (`ThesisVault/Topics/Goal_O2_O3.md`); the
user said go on 2026-10-23.

**Method.** `o3port.py` as of 791c0bcc (default options). An AS-tied read (in >= 2 copies' piles within 10 AS) is handled as follows:
- it goes to the copy its haplotype is matched to in every maximum allele-slot matching (a forced edge);
- or it goes to the certified new copy it belongs to;
- otherwise it is ABSTAIN.

**Panel.** The fresh gorilla families whose present copies produce tied reads in the tie scan
(`bench/o3_graph/o2tiescan.py`, 10 reads per allele; ties only, no O2 result looked at): GWFAM168, GWFAM224, GWFAM328, GWFAM430.
GWFAM415, 407, 28 and 401 are development (P1) families and were used for development.

**Simulation.** `bench/o3_graph/o2sim.py`, 100 reads per allele, seed 1. Two cases per family:
- control: nothing erased;
- missing: the copy with the most tied reads in the control is masked and left out of the table.

**Truth and baseline.** Truth is each read's origin. The baseline is minimap2's primary.

### Decision rules (pooled over the 8 runs)
- **O2-1:** wrong <= 2% of placed tied reads (a read of a copy absent from the table is right only on a certified new copy).
- **O2-2:** fixes > breaks against minimap2.
- **Verdict.** O2 passes iff O2-1 and O2-2 pass. With nothing placed, O2-2 fails (0 > 0 is false).
- Run once; nothing is changed after looking.

**Development** (P1 families, same simulation): 0 wrong of 24 placed; 1 fix, 0 breaks; 1,196 of 1,220 tied reads abstained
(identical copies with free slots: nothing forced).

### Result (run once, 2026-10-23)
| run | tied | placed | wrong | abstain | fixes | breaks | minimap2 wrong |
|---|---|---|---|---|---|---|---|
| GWFAM168 control | 200 | 0 | 0 | 200 | 0 | 0 | 200 |
| GWFAM168 missing | 0 | | | | | | |
| GWFAM224 control | 200 | 0 | 0 | 200 | 0 | 0 | 200 |
| GWFAM224 missing | 0 | | | | | | |
| GWFAM430 control | 2 | 2 | 0 | 0 | 2 | 0 | 2 |
| GWFAM430 missing | 12 | 9 | 0 | 3 | 9 | 0 | 12 |
| GWFAM328 control (27 copies) | 42 | 42 | **24** | 0 | 15 | 0 | 39 |
| GWFAM328 missing | not run: exceeded the heavy-job time limit twice | | | | | | |

**Pooled over 7 runs:** 53 placed, 24 wrong (45%); 26 fixes, 0 breaks.
- **O2-1 FAILS** (bar <= 2%). The missing 8th run could only change this by placing more than 1,100 tied reads with no new error.
- **O2-2 passes.**
- **Verdict: O2 as forced matches DOES NOT PASS.**

**What it shows:**
- In two- and three-copy families the rule is safe but mostly abstains: identical copies with free slots force nothing.
- In the 27-copy family, forced edges were wrong 24 of 42 times. With many near-identical copies, haplotypes merged across piles
  (the read-overlap merge) can be forced onto the wrong copy.


---

## PREREG_o2_p9_abstain_2026-10-23.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_p9_abstain_2026-10-23.md`)*

## Pre-registration: O2 v2 on fresh gorilla fibroblast families with a synthetic near-identical duplicate (P9), 2026-10-23

**Status:** committed before any fresh family is run with `o2dup.py`. Goal milestone 3.

**Criterion 3, reframed by the user (2026-10-23, option i).** O2 places AS-tied reads with <= 2% wrong and abstains on ties that
RNA cannot decide; it never places reads on a copy the sample lacks. "More fixes than breaks" is no longer required. In a copy set
with free allele slots (homozygous copies), no tie is decidable from RNA (Theorem A), so abstaining there is the right answer.

**Method.** `o3port.py` as of ee6fb354, default options (O2 v2).

**Substrate** (`bench/o3_graph/o2dup.py`, `O2DUP_CONV=1`, 100 reads per allele, seed 1; development: GWFAM248, GWFAM415):
- copy X (most allele differences; synthetic B allele with 6 substitutions when < 2) duplicated as X' with d substitutions on its
  own contig, homozygous;
- the sample's X B allele carries X''s base at X''s first substitution (a conversion) when d > 0;
- case control (the sample has X'), case missing (the sample lacks X'; X' stays in the reference and the table).

**Panel.** Every fresh >= 3-copy fibroblast family with <= 13 copies (outside the 53 development families, refabsent's 34,
GWFAM248 and P7's four), sorted by name: GWFAM135, GWFAM264, GWFAM332, GWFAM34, GWFAM36, GWFAM419, GWFAM433, GWFAM441.
GWFAM2 (27 copies) and GWFAM1 (36) are left out for runtime (P7's 27-copy run exceeded the time limit). A family whose run fails
(no isoform, time limit) is reported and skipped. d in {0, 2}.

**Truth.** Each read's origin.

### Decision rules (pooled over every run that completes)
- **Q1:** wrong <= 2% of placed tied reads (with nothing placed, Q1 holds: 0 wrong).
- **Q2:** 0 tied reads placed on X' in the missing case (a copy the sample lacks).
- **Verdict.** Pass iff Q1 and Q2. Run once; nothing changed after looking.

**Reported, not decided:** tied reads, placed, abstained, fixes and breaks against minimap2, minimap2's wrong count, runtime.

### Result (run once, 2026-10-23; work dir `/mnt/linuxdisk/tmp/o2p9`): PASS, but only by abstaining
All 8 families x d {0, 2} x {control, missing} = 32 runs completed (GWFAM36 needed the synthetic B allele).

| | tied | placed | wrong | on X' (missing) | abstain | minimap2 wrong |
|---|---|---|---|---|---|---|
| d 0 control | 3,201 | 0 | 0 | | 3,201 | 1,612 |
| d 0 missing | 1,601 | 0 | 0 | 0 | 1,601 | 788 |
| d 2 control (conversion) | 3,201 | 0 | 0 | | 3,201 | 410 |
| d 2 missing (conversion) | 1,601 | 0 | 0 | 0 | 1,601 | 408 |
| **pooled** | **9,604** | **0** | **0** | **0** | **9,604** | **3,218 (33.5%)** |

- **Q1 PASS** (0 wrong of 0 placed). **Q2 PASS** (0 reads on the lacked copy X'; minimap2 put about half of X's tied reads there).
- **Verdict: PASS under the reframed criterion 3, and the pass is pure abstention.** O2 v2 placed nothing: X' is homozygous in
  every run, so a free allele slot remains and no haplotype-to-copy match is forced (Theorem A). The result shows O2 is safe where
  minimap2 is wrong a third of the time, not that it resolves any tie. Fixes 0, breaks 0.


---

## PREREG_o2_semi_p11_2026-10-23.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_semi_p11_2026-10-23.md`)*

## Pre-registration: O2 `--semi` (two-site lead) on fresh gorilla fibroblast families (P11), 2026-10-23

**Status:** committed before any P11 family is simulated or run. Goal milestone 3, under the user's rule "the totally tied ones are
hopeless" (keep semi-tied reads and MAPQ-0 secondaries; drop only exact ties).

**Method.** `o3port.py --semi` as of 319be331 (register 1389). A tied read (>= 2 copies within 10 AS):
- certified new copy or forced haplotype edge: as O2 v2;
- exact tie (same AS and edit distance on its two best copies): HOPELESS (abstain);
- semi-tied: its best AS / fewest-edit copy, only when the runner-up needs >= 2 more edits (two distinguishing sites); else abstain.

**Development** (GWFAM248, GWFAM415; d 0 / 1 / 2): 584 placed, 2 wrong (0.34%); minimap2 MAPQ >= 1 on the same reads 40.6% wrong.

**Panel.** The first 10 families (by name) of P10's panel: GWFAM117, 118, 119, 120, 132, 133, 139, 154, 171, 172. Two-copy gorilla
fibroblast families, used for O3 (P10) but never for O2.

**Substrate.** `bench/o3_graph/o2dup.py` with `O2DUP_CONV=1`, `O2DUP_PORT_ARGS=--semi`, 100 reads per allele, seed 1, d in {1, 2, 3}:
the synthetic duplicate X' of the family's most heterozygous copy, a conversion of X's B allele at X''s first substitution; cases
control and missing (the sample lacks X'). 60 runs. A run that fails (no isoform, time limit) is reported and skipped.

**Truth.** Each read's origin. **Baselines:** minimap2's primary; minimap2's primary where MAPQ >= 1.

### Decision rules (pooled over every completed run)
- **S1:** O2 wrong <= 2% of its placed tied reads.
- **S2:** 0 tied reads placed on X' in the missing case.
- **S3:** O2's wrong rate < minimap2 MAPQ >= 1's wrong rate on the same tied reads.
- **Power:** >= 100 tied reads placed by O2, else underpowered (no pass).
- **Verdict.** Pass iff S1, S2, S3 and power. Run once; nothing changed after looking.

**Reported, not decided:** fixes and breaks against minimap2's primary; hopeless and abstain counts; per d.

### Result (run once, 2026-10-23; work dir `/mnt/linuxdisk/tmp/o2p11`): FAILS S2
All 60 runs completed.

| | tied | placed | wrong | hopeless | minimap2 wrong | MAPQ >= 1 placed / wrong |
|---|---|---|---|---|---|---|
| d 1 | 6,000 | 0 | 0 | 120 | 2,015 | 5,882 / 1,955 |
| d 2 | 5,997 | 3,782 | 4 | 1,891 | 994 | 4,106 / 22 |
| d 3 | 2,599 | 804 | 0 | 41 | 28 | 2,558 / 6 |
| **pooled** | **14,596** | **4,586** | **4 (0.09%)** | 2,052 | 3,037 | 12,546 / 1,983 (15.8%) |

- **S1 PASS** (0.09% <= 2%). **S3 PASS** (0.09% < 15.8%). **Power PASS** (4,586 placed).
- **S2 FAILS: 2 reads placed on X' in the missing case** (GWFAM118 d2, GWFAM172 d2; both X|B reads). GWFAM118's read has edit
  distance 41 on X' and 43 on X (AS 3,771 vs 3,761): the "two-site" lead was one converted site plus one read error / alignment
  difference, because the rule counts edits, not distinguishing sites.
- **Verdict: DOES NOT PASS.** Fixes 0, breaks 0 against minimap2's primary.
- Next (development, not this held-out): count distinguishing sites (reference positions where the copies differ) instead of the
  edit-distance difference.


---

## PREREG_o2_v2_heldout_2026-10-23.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_v2_heldout_2026-10-23.md`)*

## Pre-registration: O2 v2 (forced matches, weighed on every tied copy) on fresh erased-copy simulations (P8), 2026-10-23

**Status:** committed before the families are picked or run. Goal milestone 3.

**Method.** `o3port.py` as of 4b6ba85f, default options. Its O2 output: a tied read goes to its haplotype's copy only when that
haplotype was weighed for every copy the read ties to and the edge is forced, or to the certified new copy it belongs to; otherwise
ABSTAIN.

**Panel.** `o3sim_one.py pick` (its rule: the simplest >= 3-copy family whose erased copy is >= 97% identical over >= 80% to >= 2
survivors and has >= half of its reads tied on >= 2 survivors).
- Excluded: GWFAM248 and every family used in P1-P7 or the tie scan's development families (O3SIM_EXCLUDE).
- Picked repeatedly, excluding each pick, until 3 families are found.
- Per family, the one-copy sweep's 40 configurations (expression, divergence, same vs other individual, controls), at the code's
  defaults.

**Truth.** Each read's origin. **Baseline:** minimap2's primary.

### Decision rules (pooled over every run)
- **O2-1:** wrong <= 2% of placed tied reads.
- **O2-2:** fixes > breaks against minimap2.
- **Verdict.** Pass iff both hold. Run once.

### Result (2026-10-23): cannot be run
No fresh gorilla family meets the pick rule once GWFAM248, P1-P7's families and the 53 development families are excluded. In the
remaining >= 3-copy families, no erased copy is both >= 97% identical over >= 80% to two survivors and hidden as tied multimappers.
So the fibroblast data hold no fresh held-out for O2 v2. Nothing was simulated or scored.

### After the user chose "another held-out substrate" (2026-10-23, development only; nothing fresh run)
- **Fresh pool.** 10 fresh >= 3-copy fibroblast families remain (GWFAM34, 36, 264, 332, 135, 419, 441, 433, 2, 1). The tie scan
  gives 0 AS-tied reads in every one with their real copies, so real fresh ties do not exist.
- **Synthetic duplicate substrate** (`bench/o3_graph/o2dup.py`): a fresh family plus X' = copy X with d substitutions on its own
  contig; cases control / missing; optional conversion (`O2DUP_CONV=1`: X's B allele carries X''s base at one substitution).
- **Development (GWFAM415, GWFAM248), d = 0-10, with and without conversion: O2 v2 places 0 tied reads.**
  - d = 0: the tie is symmetric (X and X' identical in the reference); RNA cannot decide; abstaining is right.
  - d >= 2 without conversion: minimap2 is already right (1 / 396 wrong); nothing to fix.
  - with conversion (minimap2 52-100 wrong): X' is homozygous, so its second slot is free and X's B haplotype fits it (Theorem A);
    no edge is forced, so v2 abstains, correctly.
- **So on this substrate O2-2 (fixes > breaks) cannot pass by construction** unless the substrate is tuned (heterozygous X', fewer
  free slots) until v2 places reads. That would be fitting the test to the method, so no P9 was pre-registered.


---

## PREREG_o2_guard_p13_2026-10-24.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_guard_p13_2026-10-24.md`)*

## Pre-registration: O3 -> O2 with `--semi --guard` on fresh gorilla fibroblast families (P13), 2026-10-24

**Status:** committed before any P13 family is simulated or run.

**Method.** `o3port.py --semi --guard` as of 500819a1 (registers 1391, 1394): P12's O2, plus O2 abstains on every tied read
touching a copy where O3 flags a putative copy (or a copy that putative haplotype could take a slot on).

**Development** (GWFAM248, GWFAM415, case `extra`, d 2 / 3): flagged run GWFAM248 d3 91% -> 0 wrong; the silent run GWFAM248 d2
stays 40% wrong (nothing flagged; Theorem A).

**Panel.** P10 panel families 21-30 by name (never used for O2): GWFAM259, 271, 277, 289, 29, 293, 295, 297, 302, 308.

**Substrate.** `o2dup.py`, `O2DUP_CONV=1`, `O2DUP_PORT_ARGS="--semi --guard"`, `O2DUP_CASES=control,missing,extra`, 100 reads per
allele, seed 1, d in {1, 2, 3}: 90 runs. `extra` = control plus a third near copy X* (d substitutions, seed 3) the reference lacks;
X*'s reads are right only on a certified new copy. A failed run is reported and skipped.

**Baseline:** minimap2's primary where MAPQ >= 1, on the same tied reads.

### Decision rules
- **S1:** control + missing pooled: wrong <= 2% of placed. **S2:** 0 tied reads on X' in missing runs.
- **S3:** control + missing pooled: wrong rate < minimap2 MAPQ >= 1's. **Power:** >= 100 placed in control + missing.
- **E1:** `extra` runs where O3 flags a putative copy on X or X' (the guarded runs): O2 wrong <= 2% of placed (0 placed passes).
- **Verdict:** pass iff S1, S2, S3, power and E1. Run once; nothing changed after looking.

**Reported, not decided (the limit):** the number of `extra` runs with no flag (silent), and O2's and minimap2's wrong rates there.

### Result (run once, 2026-10-24; work dir `/mnt/linuxdisk/tmp/o2p13`): FAILS S2
All 90 runs completed.

| arm | runs | tied | placed | wrong | MAPQ >= 1 placed / wrong |
|---|---|---|---|---|---|
| control | 30 | 9,857 | 3,373 | 2 | 8,682 / 891 |
| missing | 30 | 5,301 | 1,252 | 2 | 4,207 / 886 |
| extra, flagged (guarded) | 7 | 1,184 | 0 | 0 | 1,180 / 182 |
| extra, silent (limit) | 23 | 13,188 | 5,654 | 2,254 (39.9%) | 11,947 / 5,154 (43.1%) |

- **S1 PASS** (4 / 4,625 = 0.09%). **S3 PASS** (0.09% < 19.7%). **Power PASS.** **E1 PASS** (guarded runs: 0 placed).
- **S2 FAILS: 2 X|B reads on X' in missing runs** (GWFAM293 d2 r079, GWFAM297 d2 r098). Cause: an indel beside a distinguishing
  site shifts the read's two alignments by 1 bp, so adjacent positions compare non-homologous reference bases and count as sites.
- **The limit (reported):** O3 flags the extra copy in 7 / 30 `extra` runs; in the 23 silent runs O2 is about as wrong as minimap2
  (39.9% vs 43.1%).
- **Verdict: DOES NOT PASS.**


---

## PREREG_o2_guard_p14_2026-10-24.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_guard_p14_2026-10-24.md`)*

## Pre-registration: O3 -> O2 with `--semi --guard` and the site homology check (P14), 2026-10-24

**Status:** committed before any P14 family is simulated or run.

**Method.** `o3port.py --semi --guard` with SITE_FLANK = 5 (register 1396): P13's method, plus a distinguishing site counts only
when both of the read's alignments are gap-free over 5 bases on each side (P13 failed S2 on two reads where an indel shifted one
alignment).

**Panel.** P10 panel families 31-40 by name (never used for O2): GWFAM31, GWFAM312, GWFAM349, GWFAM35, GWFAM363, GWFAM378, GWFAM380, GWFAM387, GWFAM389, GWFAM395.

**Substrate, baseline, decision rules (S1, S2, S3, power, E1) and the reported limit:** identical to P13
(`docs/PREREG_o2_guard_p13_2026-10-24.md`, sections Substrate to Reported; 90 runs). Run
once; nothing changed after looking.

### Result (run once, 2026-10-24; work dir `/mnt/linuxdisk/tmp/o2p14`): PASS
All 90 runs completed.

| arm | runs | tied | placed | wrong | MAPQ >= 1 placed / wrong |
|---|---|---|---|---|---|
| control | 30 | 9,553 | 2,496 | 4 | 8,521 / 885 |
| missing | 30 | 5,116 | 935 | 0 | 4,162 / 880 |
| extra, flagged (guarded) | 8 | 1,476 | 0 | 0 | 1,457 / 314 |
| extra, silent (limit) | 22 | 12,556 | 3,988 | 1,585 (39.7%) | 11,458 / 4,965 (43.3%) |

- **S1 PASS** (4 / 3,431 = 0.12%). **S2 PASS** (0 on X'). **S3 PASS** (0.12% < 13.9%). **Power PASS.** **E1 PASS** (0 placed).
- **Verdict: PASSES.**
- **The limit (reported):** O3 flags the extra copy in 8 / 30 `extra` runs; in the 22 silent runs O2 is about as wrong as minimap2
  (39.7% vs 43.3%): an extra copy that fits a free allele slot is invisible to RNA (Theorem A) and misleads both.


---

## PREREG_o2_semi_sites_p12_2026-10-24.md

*(original: `docs/archive/2026-10/merged/PREREG_o2_semi_sites_p12_2026-10-24.md`)*

## Pre-registration: O2 `--semi` with supported distinguishing sites on fresh gorilla fibroblast families (P12), 2026-10-24

**Status:** committed before any P12 family is simulated or run. Goal milestone 3, user rule "only exact ties are hopeless".

**Method.** `o3port.py --semi` as of 3b11426e (register 1391). A tied read nothing forces: HOPELESS on an exact tie (same AS and
edit distance on its two best copies); otherwise placed on the leader only when it leads by >= 2 distinguishing sites (reference
bases of the two copies differ), each counted only if the read's base is carried by >= 6 reads of the leader's pile at that column.

**Development** (GWFAM248, GWFAM415; d 0/1/2, conversion): 582 placed, 1 wrong (0.17%), 0 on the lacked copy.

**Panel.** P10 panel families 11-20 by name (never used for O2): GWFAM176, 178, 189, 193, 194, 226, 233, 255, 257, 258.

**Substrate, truth, baselines, decision rules:** identical to P11 (`docs/PREREG_o2_semi_p11_2026-10-23.md`): o2dup.py,
`O2DUP_CONV=1`, `O2DUP_PORT_ARGS=--semi`, 100 reads per allele, seed 1, d in {1, 2, 3}, cases control / missing, 60 runs.
- **S1:** wrong <= 2% of placed tied reads. **S2:** 0 tied reads on X' in missing runs. **S3:** wrong rate < minimap2 MAPQ >= 1's.
- **Power:** >= 100 placed. **Verdict:** pass iff all four. Run once; nothing changed after looking.

### Result (run once, 2026-10-24; work dir `/mnt/linuxdisk/tmp/o2p12`): PASS
All 60 runs completed.

| | tied | placed | wrong | hopeless | MAPQ >= 1 placed / wrong |
|---|---|---|---|---|---|
| d 1 | 6,000 | 0 | 0 | 57 | 5,943 / 1,987 |
| d 2 | 6,000 | 3,959 | 1 | 1,686 | 4,314 / 10 |
| d 3 | 2,240 | 603 | 0 | 12 | 2,228 / 4 |
| **pooled** | **14,240** | **4,562** | **1 (0.02%)** | 1,755 | 12,485 / 2,001 (16.0%) |

- **S1 PASS** (0.02%). **S2 PASS** (0 reads on X' in missing runs). **S3 PASS** (0.02% < 16.0%). **Power PASS** (4,562).
- **Verdict: PASSES.** At d 1 (one distinguishing site) every read abstains, by design.


---

## REPORT_o2_refabsent_2026-10-24.md

*(original: `docs/archive/2026-10/merged/REPORT_o2_refabsent_2026-10-24.md`)*

## Report-only: O3 -> O2 on KB3781's real reference-absent copies (refabsent panel), 2026-10-24

**Status:** written before the port is run on these families. The refabsent panel (34 families with haplotype-only loci, truth
from KB3781's own haplotypes, registers 1205-1207) has been used before (O3 stage, 2026-10-01), so this is a REPORT, not a test:
no verdict, nothing tuned on it.

**Method.** `o3port.py --semi --guard` as of 11779a31 (P14's method), on `refabsent/R0.bam` (fibroblast reads on `_pri`), copies =
`refabsent/copies_genes.tsv`, read set = `refabsent/fam/<family>.fa`. Script: `bench/o3_graph/o2refabsent.py`.

**Truth.** Each read's best record on the haplotype the family's copies lift to: on a lifted copy -> that copy; on a B-only locus
-> absent; a tie between two loci -> unknown (dropped). Within-delta B-only loci (114 of 127) are near-identical to a `_pri` locus,
so their truth is close to allelic and is reported apart from beyond-delta loci.

**Reported.** (1) Per family: expressed B-only loci vs the port's certified / putative flags. (2) O2 on tied reads from reference
copies: placed, wrong, vs minimap2 MAPQ >= 1. (3) O2 on tied reads from B-only loci: placed anywhere but a certified new copy =
wrong; vs minimap2 MAPQ >= 1 (always wrong there).

### Result (run 2026-10-24; work dir `/mnt/linuxdisk/tmp/o2refabsent`)
Truth: 16,236 reads from reference copies, 3,408 from B-only loci, 1,108 tied between loci (dropped). A first truth build mapped
haplotype contigs by their `chrN` prefix, so `chrN_*_random_utig*` overwrote `chrN` and 87 correct GWFAM64 placements scored as
wrong; fixed before the numbers below (the port's output was not rerun or changed).

**O3 flags vs expressed B-only loci** (8 families have one; reads in parentheses):

| family | expressed B-only loci | certified | putative |
|---|---|---|---|
| GWFAM175 | B0 (281), B1 (6) | 0 | 0 |
| GWFAM205 | B0 (77) | 1 | 0 |
| GWFAM227 | B0 (13) | 0 | 0 |
| GWFAM26 | B2 (13), B3 (9) | 0 | 0 |
| GWFAM382 | B0 (4) | 0 | 0 |
| GWFAM390 | B1 (20) | 1 | 0 |
| GWFAM4 | B0 (5) | 0 | 0 |
| GWFAM70 | B1 (16), B4 (24) | 0 | 0 |
| GWFAM256 | none | 0 | 1 |
| GWFAM360 | none | 1 | 0 |

- O3 flags 2 of the 8 families with an expressed haplotype-only locus (GWFAM205, GWFAM390); 2 of the 26 families without one are
  flagged (GWFAM256 putative, GWFAM360 certified). Below MIN_SUPPORT (6 reads) a locus cannot be flagged: GWFAM382, GWFAM4, GWFAM175_B1.

**O2 on tied reads:**
- From reference copies: 432 tied, **157 placed, 0 wrong**; minimap2 MAPQ >= 1 418 placed, 29 wrong (6.9%).
- From B-only loci: 18 tied, 2 placed, **2 wrong** (GWFAM175, onto copy 1: the silent free-slot case of P14, now in real data);
  minimap2 MAPQ >= 1 16 placed, 16 wrong.
- 308 tied reads have no truth (tied between haplotype loci).

### Addendum: O3 with the port's IsoCon side (the stage, `candidate_copies`, run on the 34 families; work dir `.../o2refabsent/h`)
The run above had no IsoCon side (no stage products). With it (`o3port_p4.sh` stage / concat / nets / port c cc): 10 certified
candidates + 1 putative. Each contig was aligned to KB3781's two haplotypes and to `_pri`; a contig the haplotypes hold at ~1.0
but `_pri` only far lower is a copy the reference lacks (the truth table lists only loci found from the families' copies, so it
misses some):

| contig | best haplotype identity | best `_pri` identity | reading |
|---|---|---|---|
| cand_GWFAM4_1 | 0.9978 | 0.8673 | **reference-absent copy** |
| cand_GWFAM103_1 | 0.9993 | 0.8980 | **reference-absent copy** |
| cand_GWFAM175_0 | 1.0000 | 0.9378 | **reference-absent copy** (pat chr5 40.10-40.13 Mb, in the tandem array next to B1) |
| cand_GWFAM70_0 | 0.9996 | 0.9750 | **reference-absent copy** |
| putative_GWFAM256_0 | 0.9996 | 0.9908 | allele-level gap (0.9%): undecided |
| cand_GWFAM205_0 | 0.9996 | 0.9958 | allele-level (0.4%) |
| cand_GWFAM4_0 | 0.9905 | 0.9865 | allele-level / noisy consensus |
| cand_GWFAM360_0 | 1.0000 | 0.9971 | allele-level (0.3%) |
| cand_GWFAM205_1 | 1.0000 (70% covered) | 0.9988 | not absent |
| cand_GWFAM169_0 | 0.9974 | 0.9974 | not absent (false flag) |
| cand_GWFAM103_0 | 0.7148 | 0.7148 | junk contig (false flag) |

- **4 clearly reference-absent copies found, in 4 families** (GWFAM4, 103, 175, 70), each held by a haplotype at >= 0.998 and by
  `_pri` at <= 0.975. Only GWFAM175 was in the truth table's beyond-delta expressed list; the other three were unlisted.
- 4 flags are at an allele-level gap (0.3-0.9%) and cannot be called absent or allelic from this check; 3 are false flags.
- The IsoCon side is what finds them: without it the port flagged GWFAM205 / GWFAM390 / GWFAM360 / GWFAM256 and none of the four.

---

## The AS gap is a splicing artefact, not a copy signal (2026-10-24, register 1403)

**Question (advisor's reading of the Eichler rule):** gather every alignment of a read, and assign the read to a secondary whose AS is >=
10 above the primary's. The user's position: the rule was a window for flagging ambiguity, not for assignment, and the primary is the best
alignment anyway. Tested on the real gorilla fibroblast BAM `fibroblasts/GCA_029281585.2_flnc_mm.bam` (minimap2 2.31 `-ax splice:hq -uf
--eqx -Y -N 50 -p 0.1`; 13,516,443 primaries, 4,190,262 with a secondary). Script `bench/as_tie_splicing.py`; data `/mnt/linuxdisk/tmp/as_flip`.

**A secondary with a HIGHER AS than the primary is common** (so "the primary has the highest AS" is false), and almost always a worse alignment:

| best secondary AS - primary AS | reads | secondary has fewer introns (unspliced vs spliced) | secondary has more edits | same introns, fewer NM |
|---|---|---|---|---|
| 0 (exact tie) | 27,063 | 7.3% (7.1%) | 7.5% | 23 |
| +1 to +9 | 15,405 | 97.9% (95.9%) | 97.9% | 231 |
| >= +10 (the advisor's rule fires) | 399,865 | 99.9% (96.4%) | 99.4% | 193 |

Mechanism: splice-mode AS charges every intron (a long gap, ~32), so an unspliced alignment to a processed copy, with a dozen mismatches,
outscores the read's true spliced alignment; minimap2 still makes the spliced one primary (so it does not pick the primary by AS). Example
SRR27438212.8997843: primary 5 introns, NM 0, MAPQ 60, AS 592; secondary unspliced on another chromosome, NM 19, MAPQ 0, AS 626. The rule
would move ~400k reads (9.5% of reads with a secondary), ~99% of them onto an alignment with more mismatches.

**Splice-aware windows** (reads with SOME secondary in the window):

| window | reads | a candidate with > primary NM + 5 | primary MAPQ >= 20 |
|---|---|---|---|
| A raw AS within 10 (or higher) | 469,468 | 409,592 (87%) | 410,262 (87%) |
| B splice-neutral AS (AS + 32 per intron) within 10 | 44,778 | 525 (1.2%) | 3,962 (8.8%) |
| C same intron count, then AS within 10 | 40,978 | 162 (0.4%) | 2,510 (6.1%) |

A same-intron-count rule with fewest edits (NM <= primary + 1) instead of AS is NOT usable: 112,183 reads, 66% with a confident primary,
because a clipped, shorter secondary has fewer edits. **Conclusion:** compare copies only between alignments of the same exon structure
(C; B is its single-score form); a spliced-vs-unspliced difference is a structure call (processed copy or not), not a copy tie. Only 142
reads have both the advisor's trigger and a same-structure tie. Not yet run through the O2 regimen's held-out truth (counts only, no labels).

## PREREG: splice-aware tie window on simulated gorilla reads with known origin (P15), 2026-10-24

**Status:** committed before any read is simulated. Rules frozen in `bench/as_tie_splicing.py` (register 1403): windows A (raw AS within 10
or higher), B (AS + 32 per intron, within 10), C (same intron count, then AS within 10); the advisor's rule R (move the read to a secondary
with AS >= primary AS + 10, the highest such).

**Substrate (held out by construction: no truth has been seen for any rule).** Gorilla only (never pooled with human). Genes: the GGO_genomic.gff
genes (any type with mRNA / lnc_RNA / transcript children) whose span holds >= 5 primary starts of real fibroblast reads in window A
(`/mnt/linuxdisk/tmp/as_flip`); a seeded draw (seed 20261024) of 300. Per gene up to 5 of its transcripts (seeded), 20 reads each: the
spliced transcript sequence from GGO.fasta (strand-aware), full length, with 0.05% substitutions, 0.025% insertions, 0.025% deletions per base
(seeded). Aligned genome-wide with the real BAM's command: `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` against
`GGO.splice.mmi`. Truth: an alignment is RIGHT when it lies on the source transcript's chromosome and overlaps its span.

**Hypotheses and bars** (bars set by eye before the run):
- **H1 (the advisor's rule harms):** among reads where R fires, the primary is right in >= 90% and R's chosen secondary is right in <= 10%.
- **H2 (dropping structure differences is safe):** among reads in A but not in C, the primary is right in >= 98%.
- **H3 (C keeps the real ties):** of the RECOVERABLE reads A catches (primary wrong, a right secondary in A's window), C's window holds a right
  secondary for >= 90%.
- Reported, no bar: the same for B; read counts per window; primary accuracy inside each window.

**Fails if** any bar is missed; a failed bar is reported as is (no re-run with other constants).

### P15 result (2026-10-24): ✅ PASS on all three bars (register 1404)
300 genes (of 532 eligible), 842 transcripts, 16,840 reads; minimap2 2.30, 480 s, 20 GB. Primary right for 16,550 / 16,840 (98.3%);
15,831 reads have a secondary. Scripts `bench/as_tie_p15_sim.py`, `bench/as_tie_p15_score.py`; data `/mnt/linuxdisk/tmp/p15`.

| | reads | primary right | recoverable (primary wrong, a right secondary in the window) |
|---|---|---|---|
| A raw AS window | 2,240 | 1,967 (87.8%) | 270 |
| B splice-neutral AS | 1,419 | 1,146 (80.8%) | 270 |
| C same intron count + AS within 10 | 1,384 | 1,114 (80.5%) | 270 |

- **H1 PASS:** R (move to the secondary with AS >= primary + 10) fires on 812 reads; the primary is right for 809 (99.6%), R's secondary for
  16 (2.0%); 802 of R's secondaries have fewer introns than the primary. Applied, R would make ~98% of the reads it moves wrong.
- **H2 PASS:** the 856 reads in A but not C have a right primary in 853 (99.65%): the structure-difference cases are not ambiguity.
- **H3 PASS:** C holds a right secondary for 270 / 270 of A's recoverable reads (B too), with 38% fewer reads in the window.
- Caveat: reads come from reference transcripts with 0.1% error, no allelic variation and no expressed processed copies (by construction the
  parent is the truth); whether a processed copy itself is expressed is a separate call (spliced vs unspliced), outside this test.

### `--tie-structure` wired into the tie-read step (2026-10-24, register 1405)
`copy_phasing` and `o3port.py` take `--tie-structure` (opt-in; default unchanged): a record is in a copy's pile, and counts as
near-best in the landing loci, when it has the primary's intron count and an AS within 10 of the primary's (`pile::Prim::near`,
`o3port.near`). By default: an AS within 10 of the read's best, as before.
- Flag off: both ports byte-identical to the frozen goldens (81/81, Python and Rust). Flag on: Python = Rust on 13 new golden cases
  (`*_tie`), at 1 and 4 threads.
- Effect: the family-local panels (o2p14, o2p12) are unchanged, since their references hold no processed copies. Refabsent `h`: 8 fewer
  tied rows, 2 fewer assignments. P10 control: flags 3 -> 4. The new flag (cand_GWFAM257_0, 18 reads) is a haplotype whose reads the
  old rule dropped from their own copy's pile, because a secondary with one more intron aligned ~250 more read bases and won on AS.
  The contig matches KB3781's maternal haplotype at 3 / 4,989 differences and the reference / paternal at 40 (a 0.74% gap, the
  allele-level band of register 1401). So the rule recovers a real divergent maternal allele, which the copy finder then certifies as a
  copy: a false flag in the control arm, from the allele-vs-copy blind spot, not from the window.
- Why it is not the default yet: P15 validated the window on reads with no allelic variation; on real reads it feeds the copy finder
  alleles the old rule hid. Flip the default once the copy finder separates a divergent allele from an extra copy.

## PREREG: O2 placement with and without `--tie-structure` on fresh simulated gorilla families (P16), 2026-10-24

**Status:** committed before any read is simulated or any family is run. Frozen code: main 54984246 (`copy_phasing`).

**Panel.** Families of `GGO_gwcat` with 2-8 copies, at least one copy multi-exon, and a family id named nowhere in docs/, bench/ or the
earlier panel copy tables (297 ids excluded); a seeded draw (seed 20261025) of 40. Gorilla only.

**Reads.** Every copy equally expressed: 30 reads per copy, its full spliced sequence from `GGO_gwcat.copies.fa` in transcript orientation,
with 0.05% substitutions, 0.025% insertions, 0.025% deletions per base (seeded). Aligned genome-wide: `minimap2 -ax splice:hq -uf --eqx -Y
-N 50 -p 0.1 --secondary=yes` against `GGO.splice.mmi`. Truth = the source copy.

**Runs.** `copy_phasing run --semi --guard` (P14's configuration) with the default window (arm D) and with `--tie-structure` (arm T), on the
same inputs: copies table = the catalog spans, net = every read with a primary or secondary record overlapping a copy of the family (all
used), no IsoCon side. Baseline M = minimap2's primary. A read's final placement = O2's assignment when it assigns, else the primary. Right
= the assigned copy is the source / the primary overlaps the source copy's span.

**Hypotheses** (bars by eye before the run):
- **H1:** T's O2 assignments are right in >= 98%.
- **H2:** T makes no more wrong assignments than D.
- **H3:** final per-read accuracy over all panel reads: T >= D and T >= M.
- Reported, no bar: tied rows, assigned / abstain / hopeless per arm; T - D differences per family.
**Fails if** any bar is missed; reported as is.

### P16 result (2026-10-24): ❌ FAIL by the pre-registered rule (H1 undefined: O2 assigns nothing); uninformative (register 1406)
40 fresh families (203 eligible), 90 copies, 2,700 reads; data `/mnt/linuxdisk/tmp/p16`, scripts `bench/o2_p16_sim.py`, `bench/o2_p16_score.py`.
- minimap2's primary is right for 2,671 / 2,700; all 29 wrong primaries are among the 360 tied rows (6 families x 60).
- Arms D and T: the same 360 tied rows, all HOPELESS (exact ties), 0 assigned; final placements identical (2,671 for D, T and M).
  H1 (T precision >= .98) is 0 / 0: FAIL by rule. H2, H3 pass trivially.
- The tied rows: in 5 of the 6 families the two catalog "copies" are neighbouring opposite-strand genes with overlapping spans
  (one record in both piles); GWFAM462 looks like a near-identical inverted pair. Exact ties, so abstaining is right.
- Why the window cannot matter here: copy_phasing's O2 compares only the family's own copies; a processed copy elsewhere never
  enters a pile. On ideal reads of reference copies there are no alleles, and minimap2 resolves every non-exact tie.
- The only product difference: loci.tsv for GWFAM421. 30 reads land unspliced on NC_073245.2:55889087-55890310 (1.2 kb, a processed
  copy); D marks them near-best (an expressed copy's own reads), T marks them weaker, which is right.
- Lesson: an O2 test of the window needs tied reads whose rival is ANOTHER COPY OF THE FAMILY with a different structure (as GWFAM257,
  register 1405), with allelic variation; equal expression from reference copies does not produce them.

## REPORT (no verdict): `--tie-structure` vs the default window on real KB3781 reads with haplotype truth (P5, P10 control arms), 2026-10-24
**Status:** written before the runs. No fresh family remains under P4's eligibility (all 158 used), so this is report-only on spent
O3 panels; the window itself was validated held-out in P15 (register 1404). Nothing is tuned on it.
**Method.** `copy_phasing run --semi --guard --isocon k c/cand` (main 69ddbe3f) on `R0.bam` (real fibroblast reads on `_pri`), the
`c.copies.tsv` copies and `kc/net_reads.tsv`, with the default window (D) and with `--tie-structure` (T). Truth = `labels.tsv` (each
read's copy from KB3781's own haplotypes). Reported per panel: tied rows, assigned / abstain / hopeless, right / wrong among assigned
reads with a label, minimap2's primary on the same rows, reads that change pile, and the T - D differences per family.

### Result (report-only; register 1407): the window changes nothing in O2 on real reads
Work dir `/mnt/linuxdisk/tmp/p17`; scorer `bench/o2_tie_window_real_score.py`.

| panel | arm | tied rows | assigned | right / wrong (labelled) | abstain | hopeless |
|---|---|---|---|---|---|---|
| P5 control | D | 149 | 1 | 1 / 0 | 0 | 148 |
| P5 control | T | 149 | 1 | 1 / 0 | 0 | 148 |
| P10 control | D | 107 | 18 | 17 / 1 | 69 | 20 |
| P10 control | T | 107 | 18 | 17 / 1 | 69 | 20 |

- The same tied rows and the same assignments in both arms, on both panels (0 reads change). P10: minimap2's primary is right for 99 of
  the 105 labelled tied reads; O2's one wrong placement breaks a right primary, in both arms.
- P5: 144 of the 149 tied rows are GWFAM279 readthrough reads whose primary runs from copy :0 into copy :1, so one record is in
  both piles: HOPELESS is right (a readthrough, not a copy tie).
- So in copy_phasing the window reaches O2's assignments only through a same-family rival with a different structure. On these two
  panels that never decides a tied read; its real-data effects are in the copy finder (GWFAM257, register 1405) and the landing
  loci (GWFAM421, register 1406).

### REPORT (register 1409): what minimap2 ranks by is `ms`, not AS
Source (minimap2 v2.30 `hit.c` `mm_hit_sort`, `align.c` `mm_update_extra`): after base-level alignment, hits are sorted by `dp_max`
(the `ms` tag: the best local DP score along the CIGAR, with no intron cost), with ties broken by a hash of the read name. For reads
>= 500 bp (`rank_min_len`) whose top two hits are within `rank_frac` 0.9, `mm_update_dp_max` first recomputes `dp_max` with a
mismatch penalty set from the best hit's divergence. AS (`dp_score`) includes the intron and splice-site costs and is never used to rank.
On the real fibroblast BAM (minimap2 2.31 splice:hq, 4,190,262 reads with >= 1 secondary; `bench/mm2_primary_rank.py` on
`/mnt/linuxdisk/tmp/as_flip/grp*.tsv`), the primary has the highest ms for 4,162,202 (99.3%); 26,639 (0.6%) tie on ms (the hash decides:
the coin toss); 1,421 (0.03%) are below the best secondary's ms (the divergence recalibration). AS disagrees with the
primary for 415,270 (9.9%), which is why the AS + 10 rule fired so often. Of the ms ties, 26,384 have the same NM and 25,176 the same intron
count: exact ties, which no sequence score can break.

## P18 pre-registration (2026-10-25): does an error-model score beat minimap2's `ms` on near-tie reads?
**Question (user, after register 1409):** minimap2 ranks by `ms`, which weighs every mismatch alike and every 1-bp gap alike. HiFi errors
are mostly 1-bp homopolymer indels. Does a likelihood under the real read error profile place near-tie reads better?
**Scorer (fixed before the test):** log-likelihood of each alignment under per-event rates measured on 24,000 real fibroblast MAPQ-60
primaries (`bench/o2_p18_errprof.py`; substitution 7.9e-4, 1-bp homopolymer deletion 7.6e-4 / insertion 5.6e-4, non-homopolymer
indels 2-5e-5; real allele differences included). Matched column log(1 - total rate); soft-clipped base = one substitution.
Candidates = the primary + every secondary with ms >= primary ms - 10 ("band"). LL exact ties are left to the primary.
**Substrate:** simulated reads from BOTH of Jim's alleles of each GGO_gwcat copy. The reference chromosome is one haplotype's chromosome
(exact length match). The other allele is the copy's transcript hit on the other haplotype's chromosome between its 20-kb flanks'
asm5 anchors (`bench/o2_p18_alleles.py`: 1,193 built; 53 tandem copies with several hits in the window and 27 others excluded; 142 on X/Y,
which are haploid, so 20 reference-allele reads). 10 reads per allele; errors from badread's pacbio2021 error and qscore models (independent of the
scorer), identity beta(99.8, max 100, sd 0.15). Aligned as the real BAM (`splice:hq -uf --eqx -Y -N 50 -p 0.1`). Right = overlaps
the true copy's span.
**Split:** families by sha1("p18:"+family) mod 3 (0 = dev, else test); P16's 40 families excluded. Dev (153 families, 10,000 reads)
was used only to build the pipeline. Dev result (already seen): 453 band reads, 308 LL exact ties; LL agrees with the primary on
142 / 145 others (3 switches: 1 right, 2 wrong); at posterior >= 0.99, 136 assigned, 16 wrong (15 of them reads of the other
allele) vs the primary's 14 wrong on the same reads.
**Test hypotheses (held-out families, run once):**
- H1: switches (LL argmax != primary) are right more often than wrong: right > wrong and one-sided binomial p < 0.05.
- H2: the LL argmax's wrong rate on non-tied band reads < the primary's on the same reads.
- H3: posterior >= 0.99 assignments are <= 1% wrong and no worse than the primary on the same reads.
- PASS = H1 and H2 and H3. FAIL means `ms` stays the ranking; report the allele share of the wrong calls.

### P18 result (2026-10-25): ❌ FAIL (H1, H3); `ms` stays the ranking (register 1410)
Test: 294 held-out families, 14,960 reads (8,170 from the reference allele, 6,790 from the other allele); minimap2's primary is right
for 14,742 (98.5%). Band reads 471: 249 LL exact ties (the copies are identical over the read), 222 decidable.
- H1 FAIL: 2 switches, both right (p = 0.25). The error model agrees with `ms` on 220 / 222.
- H2 PASS (barely): LL wrong 40 / 222 vs primary 42 / 222.
- H3 FAIL: posterior >= 0.99 assigns 220, 40 wrong (18%), the same 40 the primary gets wrong.
- All 40 confident wrong calls are reads of Jim's OTHER allele, whose sequence matches another copy better than its own (its allele
  carries that copy's base at the distinguishing sites). No read-error model can fix that; it is the donor-allele problem O3 -> O2 is for.
- Where the 218 wrong primaries are: 146 exact ties in the band (coin toss), 42 near-ties (40 other-allele reads), 30 outside the band.
- Conclusion: sequencing errors are not what makes minimap2 wrong on near-ties; HiFi errors fall on both alignments alike. What
  makes it wrong is (a) copies identical over the read and (b) the sample's own alleles. Both need evidence from outside the read
  (the sample's copies and alleles: O3 and O2's sample-aware lane), not a better per-read score.

### REPORT (register 1411, dev only): no other tag or CIGAR feature separates the true copy on near-ties
P18 dev reads (153 families), band candidates (ms within 10), true copy vs each wrong candidate (`bench/o2_p18_features.py`,
`bench/o2_p18_features_cond.py`). Features: ms, AS, NM, de, cm, s1 (tags); and from the CIGAR + reference: introns, annotated junctions
(GGO GFF), non-canonical junctions, soft clip, aligned length, A count in the 20 bp past the 3' end.
- Primary wrong, LL-tied (147 reads, 135 of them reference-allele reads): every feature is EQUAL between the primary and the true copy.
  The copies are identical over the read, structure included.
- Primary wrong, decidable (30 reads, 29 other-allele): ms, LL, NM, de, cm and s1 favour the primary (the wrong copy) in 25-29 of the 30;
  introns, junction motifs and clips are equal. Annotated junctions split 9 / 10. Downstream-A favours the truth 10 / 0, but all 10
  are one family at 3 vs 2 A's (noise).
- So within a single read's alignments there is nothing left: ties are identical everywhere, and the near-ties are the sample's
  alleles. Not tested here, because the simulation lacks it: real read ends (TSS / polyA site per copy) and evidence pooled over the
  sample's reads (unique reads of each copy, its alleles).

### REPORT (register 1412, dev only): abundance from the sample's other reads as a tie-breaker
Option 2 of register 1411: evidence pooled over the sample's reads. Units = catalog copies (else 10-kb locus bins); each read's
candidates = its band (ms within 10) with the P18 error-model likelihood; abundances from all reads. Scripts `bench/o2_p20_*.py`;
work dir `/mnt/linuxdisk/tmp/p18`. P18's family split; dev families only; test not touched.
- **Even expression, full-length reads (P18 dev):** no information. 146 / 151 tied reads have two tied copies and neither has a
  single unique read; the copies are identical over every read.
- **Uneven expression (30% silent, else lognormal around 10 reads) + Iso-Seq-like 5' truncation (half the reads keep 20-100% from the
  3' end), dev (`p20`, 5,727 reads, 142 families):** 180 LL-tied reads; minimap2's primary is right for 80 (44%, coin toss).
  - EM point estimate, assign at >= 0.99: 40 assigned, 9 wrong (22%, all in GWFAM276 with 0 unique reads on either copy). EM pushes a
    copy's abundance to ~0 on weak soft evidence, so its "0.99" is not a real posterior.
  - Collapsed Gibbs with a flat Dirichlet prior (the abundance uncertainty enters the posterior): 0 tied reads reach 0.99, so as
    a certifier it adds nothing. As a best guess, its pick is right for 150 / 180 tied reads (83%) vs the primary's 80 (44%), at
    stated posteriors mostly 0.5-0.8: right but underconfident.
- **Real fibroblast reads (no truth, 20-kb loci, `bench/o2_p20_emreal.py`):** 26,353 exactly tied reads (same ms and NM on >= 2 loci).
  Several tied loci have unique reads: 90.2%. Exactly one: 8.1%. None: 1.7%. The top locus holds >= 99% of the tied loci's unique
  reads for 2,420.
- So: pooled abundance can put the right copy first for tied reads (83% vs 44% in dev), but it does not certify single reads at 99%.
  Its natural output is a per-copy expression split with uncertainty, or a probability per read, not an O2 assignment under
  assign-or-abstain. Sensitive to the simulated expression and truncation model; not yet held out.

## PREREG P21 (2026-10-25): does copy_phasing (O3 -> O2, read clusters on alignment columns) fix the other-allele errors of P18?
**Status:** committed before copy_phasing is run on these reads. Frozen code: main HEAD (`copy_phasing` binary as built from 69ddbe3f;
no code change since).
**Substrate:** the P18 TEST reads and BAM as they are (294 held-out families, 14,960 reads from both of Jim's alleles; `test.bam`
in `/mnt/linuxdisk/tmp/p18`). P18 measured minimap2's primary here: 218 wrong, of which 40 are confident near-tie errors on
other-allele reads and 146 are exact ties. P18's error-model score was the only thing run on them; copy_phasing was not.
**Run:** `copy_phasing run --semi --guard` (P14's configuration, default window, no IsoCon side) with the copies table = the
catalog spans of the 294 families, genome = GGO.fasta, net = every read with a record overlapping a copy of the family (all used).
**Scoring:** a read's final placement = O2's assignment when it assigns, else the primary. Right = the assigned copy is the source
copy / the primary overlaps the source copy's span (P18's rule).
**Hypotheses** (bars set before the run):
- H1: O2's assignments are right in >= 98%.
- H2: O2 fixes more primaries than it breaks (fixed > broken), and fixes >= 1 of the 40 confident other-allele errors.
- H3: final per-read accuracy >= minimap2's primary (14,742 / 14,960).
- Reported, no bar: tied rows, assigned / abstain / hopeless; the same split for other-allele reads; haplotypes, certified and
  putative copies found per family (every extra haplotype here is an allele of a reference copy, so certified copies are false).
**Fails if** any bar is missed.

### P21 result (2026-10-25): ❌ FAIL (H2); O2 stays safe but fixes none of the allele errors (register 1413)
Work dir `/mnt/linuxdisk/tmp/p21`; scorer `bench/o2_p21_score.py`; copy_phasing 28 s.
- 1,190 tied rows: 9 assigned (all right, all where the primary was already right), 180 ABSTAIN, 1,001 HOPELESS. H1 PASS (9/9),
  H3 PASS (final 14,742 = primary), H2 FAIL (fixed 0, broken 0).
- The 40 confident other-allele errors: 30 are tied rows and all ABSTAIN (GWFAM4 20; GWFAM405 10, where the other allele's cluster
  is flagged as a putative copy and `--guard` abstains); 10 (GWFAM65) are not tied among the family's copies.
- O3 side, reported: 846 haplotypes in 294 families; 2 putative copies (GWFAM318, GWFAM405); **2 certified copies, both false**:
  GWFAM197 and GWFAM384, where every haplotype is one of Jim's two alleles. In both, the divergent copy-0 allele (11 and 15
  differences) yields more clusters than the family has allele slots (5 haplotypes for 2-3 copies), so one is forced to count as
  extra. Badread's errors depend on sequence context, so they recur across reads; repeated errors that look like alleles are the
  suspect, not yet checked.
- Reading: clustering the reads by their alignment columns makes the right choice on the allele errors, which is to abstain, but
  it does not turn them into correct placements. To place them it would have to join the other allele's cluster to its own
  copy's reference-allele cluster, and nothing in the reads links the two. On a diploid sample, the certified tier can be fooled
  by a divergent allele plus errors that repeat.

### REPORT (register 1414): the two false certified copies of P21 are not recurrent errors
Dumps: `/mnt/linuxdisk/tmp/p21/X.dump` (default window), `XT.dump` (`--tie-structure`), `PY.dump` (Python port, identical).
- **GWFAM197 = the AS intron artefact (r1403).** Copy 1's 20 reads (both alleles) are spliced on copy 1 (AS 698, ms 866) and have an
  unspliced secondary on copy 0, a processed copy, with AS 793. In the raw AS window they join copy 0's pile, and copy 1's own pile
  is empty, so their haplotype (11 private columns vs copy 0) has no slot: certified. With `--tie-structure` copy 1 gets its pile
  back and nothing is certified. Whole test panel under `--tie-structure` (arm T): O2 unchanged (9 assigned, 9 right); certified 1.
- **GWFAM384 = reads outside their own copy's window.** Copy 0's other allele (15 edits vs copy 0, 65 vs copy 1 over the whole
  transcript) aligns better locally to copy 1's locus (primary on copy 1, AS 3185; on copy 0 AS 3160, NM 27 vs 22). Its 10 reads are
  more than 10 AS behind on copy 0, so they never enter copy 0's pile. Copy 0 has a free slot (one haplotype seen), but the reads
  don't have copy 0 as an option, so copy 1's pile shows a third haplotype: certified.
- Rule this points to (not implemented): a certificate must hold over every copy the aligner placed the reads on, not only those in
  the 10-AS window. Certify only if no copy where the reads have any record has a free slot. GWFAM384's reads have a record on copy
  0, which has a free slot, so they would not be certified. The rule needs a held-out check that certified recall (P5) survives.
- The original suspect, badread's context errors repeating like alleles, is not the cause in either family.

## PREREG P22 (2026-10-25): `--certify-all-records` held out: precision and recall of certified copies
**Status:** committed before the masked genome is built or any run on it. Code: `--certify-all-records` in `o3port.py` and
`copy_phasing` (opt-in; Rust = Python on 28 flag-on golden cases at 1 and 4 threads, the 94 earlier cases unchanged).
**Why held out:** the rule was built from P21's GWFAM197 / GWFAM384 (P18 test families). The P18 DEV families (153, 10,000 reads
from both of Jim's alleles) have never been used for certification.
**Panel:** dev families with >= 2 copies whose copies all have built alleles or sit on X/Y. A seeded draw (seed 20261025) puts half
of them in arm EXTRA: one copy each (seeded) is replaced by N in a copy of GGO.fasta. Its reads become reads of a copy the
reference lacks. The other half are CONTROL (nothing masked). The dev reads are realigned to the masked genome with the real-BAM
command. The copies table drops the masked copies; net = every read with a record on a remaining copy of the family.
**Runs:** `copy_phasing run --semi --guard` with and without `--certify-all-records`, each with and without `--tie-structure` (4 arms).
**Truth:** a certified copy is TRUE if most of its reads come from the masked copy of its family (either allele), else FALSE.
**Hypotheses** (bars before the run), comparing the flag arm with the same window without the flag:
- H1 (precision): false certified copies with the flag <= without, and the flag arm has <= 2 false certified copies over all
  families (default window).
- H2 (recall): masked copies certified with the flag >= 90% of those certified without it (default window).
- Reported, no bar: the same under `--tie-structure`; putative flags; O2 assigned / right on tied rows.
**Fails if** H1 or H2 is missed.

### P22 result (2026-10-25): FAIL H2 (register 1415)
Run: `/mnt/linuxdisk/tmp/p22` (panel.tsv, mask.bed, dev.bam on the masked genome, arms W / WC / T / TC, scorer `score.py`).
141 eligible dev families: 70 EXTRA (one copy masked, 3.37 Mb of N), 71 CONTROL; 7,860 reads realigned, 7,178 net rows.

| arm | certified | true | false | masked copies found (of 70) | putative (EXTRA / CONTROL) |
|---|---|---|---|---|---|
| default window (W) | 14 | 13 | 1 | 11 | 11 / 2 |
| + `--certify-all-records` (WC) | 10 | 9 | 1 | 8 | 11 / 2 |
| `--tie-structure` (T) | 17 | 16 | 1 | 13 | 10 / 2 |
| `--tie-structure --certify-all-records` (TC) | 12 | 11 | 1 | 9 | 10 / 2 |

- H1 PASS: 1 false certified copy in every arm, the same one (GWFAM256 EXTRA: copy 1 masked, its two alleles and copy 0's other
  allele compete for copy 0's one free slot; the matching certifies a copy-0 allele). 0 false in the CONTROL families.
- H2 FAIL: the flag keeps 8 / 11 masked copies (73% < 90%). Losses (diagnosed from `--dump`, after scoring):
  GWFAM345 and GWFAM401 - a copy whose pile has no variable column yields no haplotype, so its slots look free to the extended
  matching and the masked copy's haplotype takes one; GWFAM204 - probably a reshuffle: the extended matching leaves a different
  haplotype unmatched, one O2's matching placed, so nobody is certified (not checked line by line).
- O2: 0 assignments on tied rows in every arm (400 HOPELESS, 97-99 ABSTAIN).
- Recall of the certified tier is low on this panel: 11-13 / 70 masked copies (16-19%).
- Decision: `--certify-all-records` stays opt-in, not adopted. `--tie-structure` found 2 more masked copies, no extra false one.

## PREREG O2-TIE-MS (register 1472, 2026-10-30): O2's tie window defined on `ms`, the score minimap2 ranks by

**User:** "does it work for tied reads, not only AS but also ms since that is the true tag that determines primaries from
secondaries?" then "lets run it". Written before any P12 run directory is re-run under the new flag (the P12 substrate is being
copied from the other machine to `~/work/jfim/rustle_data/reads/o2p12/` as this is committed; `--tie-ms` is implemented and
unit-tested, commit pending).

**Why.** Register 1409: minimap2 ranks hits by `ms` (best local DP score, no intron cost) and breaks exact `ms` ties by a hash of
the read name; AS (with intron and splice costs) is never used to rank. O2's pile window (`o3port.near`) is AS within 10 of the
read's best AS, so "AS-tied" is O2's own definition, not minimap2's coin toss.

**Measurement on real BAMs (no bar; `bench/o3_graph/tie_ms_vs_as.py`, `~/runs/o2_ties/`), KB3781 fibroblast on GGO, reads with
>= 2 non-supplementary records (4,190,262):**
- Genome-wide: 26,669 exact `ms` ties (26,405 same NM, 25,210 same intron count); 25,172 (94.4%) fall inside O2's AS window,
  1,497 outside, nearly all (1,236 of 1,497) with an intron-count difference >= 5 and an AS gap >= 100: a spliced and an unspliced
  placement with the same `ms` (the processed-copy artefact of register 1403). O2's AS window holds 82,070 reads, 56,683 of them
  not `ms` ties: reads minimap2 did rank (primary = max `ms` in 56,467).
- Restricted to reads whose records all fall on one family's copies (kbcols, 165,931 reads): 16,949 `ms` ties, **16,940 (99.95%)
  inside the AS window**, 9 outside (intron difference 3-5). The AS window holds 27,401 reads, 10,461 of them (38%) not `ms` ties.
- Reading: within families O2 already sees every coin-toss read; the AS window is wider than the coin toss, not narrower. The
  question for the bar is whether narrowing it to `ms` changes placements or error.

**Arms (P12 substrate: 10 fresh gorilla families x d 1,2,3 x control/missing = 60 runs, `o2dup.py` products as copied; the port
re-run per directory with `o3port.py run --semi [--tie-ms W]`, scored by `o2eval_sim.py`, nothing re-simulated):**
- A (reproduction): default window. Must reproduce P12 pooled: tied 14,240, placed 4,562, wrong 1.
- M0 (registered): `--tie-ms 0`, the window = records with the read's best `ms` (minimap2's exact ties; the hopeless rule then
  compares `ms` and NM).
- M10 (reported): `--tie-ms 10`.

**Bars (M0 vs A, pooled over the 60 runs; run once, nothing tuned):**
- R: A reproduces P12 (placed 4,562, wrong 1). Otherwise the comparison is void and reported as such.
- B1: M0 wrong <= 2% of M0 placed, and M0 placed >= 100 (P12's own S1 and power bars).
- B2: M0 placed >= 0.95 x A placed (the `ms` window loses at most 5% of the placements).
- B3: M0 wrong rate <= A wrong rate + 0.5 percentage points.
- Reported: tied-set sizes per arm; reads tied in one window and not the other; hopeless counts; M10 the same.
- PASS = R, B1, B2, B3. Then `--tie-ms 0` is proposed as the default window (a separate decision, since P14's `--guard` and the
  refabsent panel were not re-run here). FAIL = any bar fails; `--tie-ms` stays opt-in and the page keeps "AS-tied".

**Not claimed:** anything about reads outside families (the 1,497 genome-wide `ms` ties outside the AS window are spliced vs
unspliced placements and are the copy finder's problem, register 1403); anything about the Rust port (`copy_phasing` has no
`--tie-ms` yet).
