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

## Results (2026-10-15)

Run: code `39e333c1` (Amendment 1), `copy_assign` from `main` `a493a66a`. Outputs are in `/mnt/linuxdisk/tmp/o3_realdata/`
(`heldout2/`, `dev/`, `control/`).

**Cost.** The masked genome's index took 195 s at 21.2 GB. The 24,920 held-out reads were mapped in 4 parts per arm (about
220 s per part; one call of all of them hit the 10-minute limit). Each `copy_assign --psv-clusters` run took 166-478 s at
≤ 2.4 GB.

**Held out (Amendment 1's 60 families: 48 with 2 copies, 12 with >= 3).**
- Reads: 9,893 masked-copy (D) and 15,027 survivor (S). Masked arm: 24,904 primary and 16 unmapped; added-back arm: 24,920
  primary.

| | held out | development (A13, 53 families) |
|---|---|---|
| eligible masked copies (>= 6 reads placed on the family) | 12 (all >= 3-copy) | 25 |
| **R1: recovered** | **4 / 12 = 0.333, FAIL** (bar 0.70) | 15 / 25 = 0.600 |
| **R2: cleared after the copy is added back** | **3 / 4 = 0.750, FAIL** (bar 0.90) | 9 / 15 = 0.600 |
| labels of the recovered clusters | `copy` 3, `allele_or_copy` 1 | `copy` 10, `allele_or_copy` 5 |
| unexplained clusters, removed / added-back arm | 58 / 138 | 177 / 210 |

**One-copy families get no placed reads.** `copy_assign` writes no assignment row for a family with one copy, so
`--psv-clusters` has no reads to cluster there. Two consequences follow, both found while scoring:
- **The 48 two-copy held-out families are untested.** After masking they have one copy, and none of their masked-copy reads
  is placed. They are not eligible, so R1 and R2 stand on the 12 families with >= 3 copies.
- **R3 is not measurable.** The single-copy control (300 random genes and 26 `high_het`, each a one-copy family) placed no
  read, so its 0 / 300 false `copy` calls are vacuous and are not counted as a pass.

Amendment 1 was chosen on the understanding that a one-copy family is clustered like any other. That understanding came from a
check of mine that confirmed the run finished, not that it placed reads, so it was wrong.

### Verdict (2026-10-15)

**Reading 3: not claimed.**
- R1 fails: 4 of 12 eligible masked copies are recovered held out, against a bar of 0.70; development gives 15 of 25.
- R2 fails as well.
- R3 cannot be measured.

### Post hoc (descriptive, not pre-registered; `bench/o3_realdata/posthoc.py`)

- **About half the misses are fragmentation.** The masked copy's placed reads often sit in several PURE unexplained clusters
  (purity >= 0.9), none of which holds 80% of them alone. Taken together, those pure clusters hold >= 80% of the reads in
  8 of 12 held-out copies (4 by the registered single-cluster rule) and 18 of 25 development copies (15).
- The remaining misses:
  - **unlinked:** the reads form a cluster the linkage test sets aside (MCL20, GWFAM4);
  - **impure:** the reads share an unexplained cluster with survivor reads (MCL34, MCL69, GWFAM123, GWFAM402, GWFAM407,
    GWFAM415);
  - **explained or unresolved:** MCL27;
  - **in no cluster:** GWFAM47;
  - **just under the bar:** GWFAM331, with combined recall 0.795.
- **Specificity on real reads is poor.** In the added-back arm every copy is present, yet unexplained clusters with
  est_div > delta (`copy`) remain in 17 of 60 held-out families and 16 of 53 development families.

  | | held out | development |
  |---|---|---|
  | unexplained clusters labelled `copy` / `allele_or_copy` | 57 / 81 | 27 / 183 |
  | `copy` clusters made mostly of survivor reads | 41 of 57 | 21 of 27 |
  | `copy` clusters made mostly of the masked copy's reads | 16 | 6 |

  The simulation flagged 0 of 152 intact families. So on real reads, a `copy` cluster in a family whose copies are all in
  the genome is common.
- **What a fix would need** (for a new prereg, not done here):
  - place a one-copy family's reads on its copy when `--psv-clusters` is on, which gives 2-copy losses and the single-copy
    control a test;
  - decide whether fragments of one source should count together (they did in 4 held-out families);
  - find what the `copy` clusters of intact families are before relying on the label: collapsed copies in the assembly,
    isoforms or exons the catalog copy does not carry, or reads from other loci.

### Follow-up (post hoc, 2026-10-15): what are the `copy` clusters of intact families?

The added-back arms have every copy in the genome, yet they hold 57 `copy`-labelled unexplained clusters (held out) and 27
(development). Each was traced (`bench/o3_realdata/copyclusters.py`, `copyclusters_summary.py`; outputs
`/mnt/linuxdisk/tmp/o3_realdata/{heldout2,dev}/classified.tsv`, `copyclusters_summary.json`) with four kinds of evidence:
- each read's copy of origin (its baseline primary);
- the prototype's clustering re-run on each host. It matched every Rust cluster exactly (Jaccard 1.0, 84 / 84), and splits
  each cluster's unexplained sites into those no catalog copy covers and those covered;
- the reads aligned to KB3781's own maternal and paternal assemblies, with the family's copies lifted onto both;
- for clusters matching neither haplotype, the reads aligned to human (CHM13), chimpanzee (mPanTro3) and orangutan
  (mPonPyg2).

Rules, applied in order:
1. **other family's copy**: most reads match a haplotype (>= 0.995 identity, >= 0.9 coverage) at another family's copy.
2. **not in the genome**: most reads match neither haplotype at >= 0.995. These are then given the species with the
   highest median identity if it is >= 0.99. Below that, they are gorilla-like when gorilla is closest by >= 0.02, else
   unresolved.
3. **missing from the reference**: the reads match a haplotype locus that is no copy of the family.
4. Own-family reads:
   - **catalog lacks the exons**: >= half the unexplained sites lie where no catalog copy has sequence;
   - **minority segment**: the sites are covered by a median < 0.2 of the cluster's reads;
   - **own-copy variant**: otherwise.

| what the cluster is | held out: clusters / reads | development: clusters / reads |
|---|---|---|
| reads of another gorilla family's copy | **35 / 2,292** | 0 / 0 |
| other apes' transcripts (chimpanzee; orangutan) | 4 / 103; 1 / 5 | 2 / 13; 2 / 8 |
| gorilla-like sequence in neither haplotype assembly | 3 / 1,985 | 0 / 0 |
| a copy missing from the reference assembly | 0 / 0 | **3 / 469** (GWFAM175) |
| the catalog copy lacks the exons | 4 / 320 | 2 / 16 |
| a few reads' divergent segment flips the cluster consensus | 0 / 0 | 5 / 930 |
| allele-like variant of the copy's own reads | 5 / 129 | 4 / 39 |
| unresolved (no species >= 0.99) | 5 / 183 | 9 / 38 |
| total | 57 / 5,017 | 27 / 1,513 |

**Reading.**
- **Reads from other loci** cause most of these clusters. On the shipped catalog, 35 of 57 are reads of another
  family's copy. They sit in 7 families (MCL69 10 clusters, MCL126 9, MCL43 6, MCL288 5, MCL296 3, MCL80 1, MCL286 1).
  Their reads match that other locus closely (each cluster's median read identity to its best haplotype hit: median
  0.999, range 0.994-1.0). `--no-as-tied-only` places them on this family, whose catalog lacks their source.
- **Other apes' transcripts** (chimpanzee, orangutan) remain in the library, as found before (A44).
- **Collapsed or missing copies are rare.**
  - Only GWFAM175's reads match a haplotype locus that is not a copy of the family (cluster medians 0.9988-0.9996): a real copy missing from
    the reference, and the case O3 is for.
  - MCL20's and MCL33's large clusters are gorilla sequence 2-3% from any locus in either haplotype. They may be copies in
    assembly gaps, but the assemblies cannot say.
- **Catalog and consensus artifacts** make the rest:
  - exons the catalog copy lacks are unexplained by construction;
  - a few reads with a locally divergent segment can turn a cluster's consensus to `alt` where only they cover the locus
    (GWFAM47: 6 of 176 reads, a ~70 bp window; the median read has 0 mismatches).

**What would cut the false `copy` calls** (each a method change, needing its own prereg):
- keep only reads whose primary or best alignment lies in the family, or cluster related families together;
- screen out other species (the A44 read-share screen);
- count only sites that a catalog copy covers;
- require the unexplained alleles to be carried by most of the cluster's reads.
