# O3 on real reads: copies removed from the genome, then added back (2026-10-15)

Pre-registration, written before any input of this test was built or run. Code: `bench/o3_realdata/` at `dd37275e`
(15 tests). Branch `o3-realdata`. Method under test: `copy_assign --psv-clusters` v2 (merged into `main` as `a493a66a`,
`docs/PREREG_o3_psv_clusters_v2_2026-10-15.md`), unchanged.

## Question

The advisor's O3 (2026-10-15): cluster the reads O2 places on a family by their PSVs; a cluster that no copy in the genome
explains better is a copy absent from the genome. In simulation it recovered 425 / 426 absent copies and flagged 0 / 152 intact
families. On real gorilla reads:

1. When a copy is removed from the GENOME and its reads are re-mapped, do they form an `unexplained` cluster on a surviving
   paralog?
2. Does the cluster go away when the copy is back in the genome?
3. Do single-copy genes, with the real heterozygous allele background, stay free of false `copy` calls?

This is step 6 of `docs/HANDOFF_laptop_to_main_2026-10-15.md`.

## Known before running (recorded now)

- **Simulation** (laptop, register 1292-1293): recovery 425 / 426 at ≥ 0.6% divergence, 0 / 152 control families flagged, RNA
  editing kept out. With exact AS ties only, the method sees twins only (80 / 310), so `--no-as-tied-only` is required.
- **The development substrate** (A13, `/mnt/linuxdisk/tmp/rna_allele/`) was built and used in 2026-10 for earlier O3 methods
  (IsoCon chain, `o3_candidates`; Amendments 7-15 of `docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md`):
  - 53 gorilla families of an older catalog (GWFAM), each with one copy hard-masked to N in the genome;
  - 59,013 real reads (up to 500 per copy) re-mapped to the masked genome (`linktest/R.bam`) and to the unmasked one
    (`control/R0.bam`);
  - the copy of origin of each read (`linktest/labels.tsv`).

  PSV clustering has never been run on it. Seen before: those methods' results (A13 detection 25 / 53 families; A14's control
  35 / 53 families with false flags), the families and their read counts.
- **Not seen**: any `--psv-clusters` output on real reads; which catalog families the held-out rule picks.
- **Allele cut-off**: delta = 0.00958, the 99th percentile of exonic divergence between KB3781's two haplotypes over its
  single-copy genes. So about 1% of single-copy genes carry alleles more divergent than delta. If expressed and clustered, such
  an allele is labelled `copy` by construction: that is the method's known limit, and the control measures its rate.

## Method (fixed)

**Removal and re-mapping.**
- The masked span is a copy's locus hull, hard-masked to N (`bedtools maskfasta`) in `GGO.fasta`.
- The reads are re-mapped with the baseline's exact command:
  `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes` (A13's `R.bam` header), sorted and indexed.
- The masked genome gets its own index (`minimap2 -x splice:hq -d`). The added-back arm maps the same reads to the unmasked
  genome with `winloci_data/GGO.splice.mmi`.

**The O3 run** (every arm):
`copy_assign --bam BAM --fasta GENOME --regions regions.txt --families COPIES.tsv --copies-fa COPIES.fa --no-as-tied-only
--psv-clusters --region-threads 1 --skip-poa-diagnostic --out PREFIX`.
The catalog is the surviving copies in the removed arm, and every copy in the added-back arm.

**Held-out families** (`prepare.py heldout`, seed 1):
- from the shipped KB3781 catalog (`runs/gorilla_KB3781/gorilla_KB3781.fam.copies.tsv`), families with ≥ 3 copies;
- every copy has ≥ 20 primary reads with an aligned base on its exons (a read on the exons of two copies counts for neither);
- every locus hull is ≤ 200 kb;
- no copy overlaps a copy of A13's 201;
- the copy last by (chrom, start) is masked, and a family whose masked hull overlaps any other catalog copy, its own siblings
  included, is dropped;
- if more than 60 families remain, a seeded sample of 60 is taken.

Reads: up to 500 per copy (seeded), labelled D (the masked copy) or S (a survivor), written in their original orientation
for `-uf`.

**Single-copy control** (`prepare.py control`, seed 1):
- candidates are KB3781's `S_single` genes with ≥ 500 exonic bp lifted ≥ 90% to the other haplotype; heterozygosity =
  (mismatches + indels) / lifted bp;
- the R3 set is a random sample of 300 expressed genes (≥ 20 primary reads on the exons), taken in seeded random order;
- also reported: up to 100 expressed genes with heterozygosity > delta (`high_het`), where `copy` calls are expected by
  construction;
- each gene is a one-copy family (its exons spliced from the genome), run on the baseline BAM and the unmasked genome.

**Scoring** (`score.py`; the simulation's rules, `bench/o3_unexplained/score.py`):
- A read is *placed* on a family when it has an assignment row of that family.
- **Recovered**: the `unexplained` cluster holding the most of the masked copy's placed reads holds ≥ 0.8 of them, with purity
  ≥ 0.9. Its label is `copy` when est_div > delta, else `allele_or_copy`.
- **Cleared** (added-back arm): no `unexplained` cluster of the family holds ≥ 0.5 of the masked copy's placed reads.
- **False copy** (control): an `unexplained` cluster with est_div > delta in a one-copy family.
- A masked copy is eligible when ≥ 6 of its reads are placed on its family in the removed arm.

## Bars (held-out families and the control; species: gorilla only)

| id | bar |
|---|---|
| R1 | recovered in ≥ 0.70 of the eligible masked copies |
| R2 | cleared in ≥ 0.90 of the recovered families |
| R3 | false `copy` in ≤ 0.02 of the 300 random single-copy genes |

**Development (A13's 53 families):** the same runs and scoring, reported without a bar. Its purpose is to exercise the
pipeline on a substrate that was built for this kind of test. If a defect shows up there, the fix is a recorded amendment
made before the held-out runs.

## Readings (fixed now)

1. **R1, R2 and R3 pass**: claimed, "O3 (PSV clusters) detects copies removed from the genome on real gorilla reads, the
   detection goes away when the copy is back, and single-copy genes stay below 2% false `copy` calls". The claim is quoted with
   the eligible and recovered counts.
2. **R1 passes, R2 or R3 fails**: detection reported; specificity not established, and the failing side is named.
3. **R1 fails**: not claimed; what happens to the masked copy's reads is reported (absorbed by a paralog, placed on another
   family, not placed), with the clusters they form.

Reported, not barred:
- the labels of the recovered clusters and the divergence of the masked copy from its nearest survivor;
- the masked copy's reads by fate;
- unexplained clusters per family in both arms;
- the control by heterozygosity bin, and the `high_het` stratum;
- the development substrate's numbers.

## Amendment 1 (2026-10-15, before any run; user decision)

Building the held-out inputs (`prepare.py heldout`, the rule above) picked only 14 families. Of the catalog's 90 families with
>= 3 copies, 62 fail the rule that every copy has >= 20 reads, 7 overlap A13 and 3 have a masked hull over another copy.
That leaves R1 coarse: 10 of 14 would pass. Seen at this point: the counts of families each rule removes, and nothing else.

**Change (the user chose it from three options):** the held-out families are those with **>= 2 copies** under the same
rules (`prepare.py heldout --min-copies 2`). 74 families pass, and a seeded sample of 60 is taken (seed 1). In a 2-copy
family, the masked copy's reads have one surviving paralog.

R1-R3 and the readings are unchanged. R1 and R2 are also reported separately for 2-copy and >= 3-copy families
(`score.py verdict --picked`). The 14-family build stays on disk (`/mnt/linuxdisk/tmp/o3_realdata/heldout/`) and is not run.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.
