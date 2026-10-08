# Figure 6 — Rustle's default de novo family definition across the paralogue identity spectrum

> **Status (2026-09-25).** The genome-wide version is pre-registered (`docs/archive/2026-09/PREREG_genome_wide_families_2026-09-25.md`,
> claims F6.1–F6.5 as amended by Amendment 1, written before any genome-wide family number) and not built yet: it
> needs the `families` stage of both human samples with its copy table (a `mcl_families` that writes
> `<id>.fam.copies.tsv`) and each sample's genome-wide spectrum. Until its tables exist (`fig6_gw_recall`,
> `fig6_gw_groups`, `fig6_gw_precision`), `make.py plot fig6` draws the development tables (human A119b chr16)
> summarised at the end of this caption. After the build, `python3 figures/fig_family_spectrum.py summary` prints every
> number marked ‹…›; replace them, state each claim as the pre-registration's table says, and delete this box.

**The families shown (stated once).** Rustle has one default de novo family definition, at the RNA level:
IsoSeq reads → loci assembled from primary alignments plus secondary alignments scoring at least 98% of the read's
best alignment score anywhere in the genome → one representative per locus, the transcript with the most reads (its
exons: coordinates from the reads, bases from the genome) → loci linked when their genomic sequences align (minimap2
asm20; ≥ 300 bp at ≥ 70% identity, covering ≥ 30% of the longer locus's exonic length, with one alignment joining
exons to exons over ≥ 60% of the smaller locus's exonic length) → Markov clustering (inflation 2.8), families of ≥ 2
loci. The same families are the copies that copy assignment uses (Figs 4–5). No protein is used anywhere in this
figure.

**Claim (pre-registered; wording decided by the outcome).** Per human sample (A119b, testis), Ensembl Compara pairs
across the whole genome, Rustle's default families:
- F6.1: at ≥ 90% protein identity the families recover ‹k of m› of the pairs they can recover (both genes are family
  members);
- F6.2: at 60–90% a direct alignment between the genes and the families are ‹not separable / separable in band …›;
- F6.3: the largest paralogue group holds ‹k of n› of the 60–90% pairs;
- F6.4: at ≥ 90%, pairs whose genes lie on two chromosomes are recovered ‹as well as / less than› pairs on one
  chromosome (‹k/n› vs ‹k/n›);
- F6.5: within-family precision is ‹k/n› on multi-exon copies and ‹k/n› over all copies.

## Definitions (each used once below)

- **Rustle's default families:** as stated above (the driver's `families` stage, `mcl_families --from-gtf
  --min-exonic-bp 1 --min-shared-exon-frac 0.60`); a family's **copies** are its member loci, each represented by its
  most-read transcript (the copy table `<id>.fam.copies.tsv`).
- **Ensembl Compara paralogue pairs** (release 116, human, exported from BioMart for every chromosome). A pair's
  **protein identity** is the larger of Compara's two percent identities. Pairs across chromosomes are kept. Genes are
  matched to Rustle's loci by gene symbol through the RefSeq CHM13 annotation.
- **Pairs scored:** Compara pairs whose two genes both have a locus in the sample's genome-wide assembly (one
  representative transcript per locus, the one with the most reads, ≥ 200 bp). This is conditioned on the
  assembler's loci, not on read counts, and it is the same for the direct alignment and for the families.
- **Paralogue group:** a connected component of the scored pairs (built here; not an Ensembl gene tree).
- **Direct alignment** (grey): some locus of one gene aligns to some locus of the other with minimap2 (asm20, ≥ 80%
  nucleotide identity; or asm20 `-k11 -w5`, ≥ 60%), over ≥ 50% of the shorter representative (query span ÷ shorter
  length). Nucleotide identity here is matching bases ÷ alignment block length. Nucleotide only.
- **Same family** (blue + navy): both genes have copies in one default family. Blue: their loci also align directly
  (the grey bar's test, same pair); navy: in one family without a direct alignment (joined through other loci of the
  family, or by the family rule's genomic alignment where the representatives' direct alignment falls below its floor).
  **Most the families can recover** (dashed outline): both genes are members of some family.
- **Precision over judgeable pairs:** a predicted pair is judged only when both genes have Compara data; pairs with a
  gene Compara does not list are not counted, so **precision is an upper bound**.
- **Genome minus development contigs:** human chr16 (the family rules were developed there) removed.

## Panels

**a** Human A119b (T2T-CHM13 v2.0; tissue not recorded) and human testis (T2T-CHM13 v2.0; aligned without `-uf`), whole
genome. Sensitivity of Compara pairs by the pair's protein-identity band: left bar = direct alignment, right bar = same
default family (blue with a direct alignment, navy without one), dashed outline = the most the families can recover.
Numbers: pairs recovered by the direct alignment; k/m = pairs in one family / pairs the families can recover.

**b** Per human sample: the three bands from 60 to 90% pooled, the largest paralogue group vs the other groups; and the
≥ 90% band split into pairs whose genes lie on two chromosomes vs on one. Square = direct alignment, circle = same
family; the numbers at the right are the recovered pairs (alignment · family).

**c** Per human sample: precision against Compara over judgeable pairs, for the direct alignment (asm20 `-k11 -w5`) and
for the families' within-family gene pairs (all copies; multi-exon copies only). Open circle = without the family that
holds the most judged pairs; k/n at the right with that value in brackets. No interval is drawn: pairs of one family
are not independent.

**Numbers shown** (‹…› = from `fig_family_spectrum.py summary` after the build): per band and sample, pairs, direct
alignment, families (with / without a direct alignment), the most the families can recover; the 60–90% largest group;
≥ 90% across vs within chromosomes; precision rows.

## Supplementary figures (same module)

**Supplementary Fig. 6s-seeding — seeding loci with secondary alignments, at the family level.** The default
(secondary alignments within 2% of the read's best score also seed loci) against primary alignments only, each through
the same families stage. Genome-wide version (pre-registered claims F6.6–F6.8, Amendment 1): per sample, genome minus
development contigs, (a) Compara pairs at ≥ 90% protein identity whose two genes both have ≥ 2 reads whose primary
alignment has an aligned block on the gene's annotated exons (human samples), (b) the same families' Compara precision,
(c) Liftoff copy pairs (a record's annotated placement and one of its extra copies, ≥ 95% identical, both loci with
≥ 2 primary reads; every sample), recovered when one family covers both loci (Liftoff's ≥ 50% of exon bases); open
circle = primaries only, filled = default, black tick = the default without its largest family. The annotation's
protein-homology families appear only as a secondary panel, and only when built (`fig6s_protein_homology 1`). The
panel moved out of the main figure because it compares two read pools, not two family definitions; Fig. 3 shows the
seeding gain genome-wide against the annotation, and on the development contig the whole family-level gain was one
tandem array. Its six `families_primary` runs are optional. Development version (drawn until the genome-wide table
exists): gorilla OR6737 chr20 (NC_073244.2, the contig used to choose the seeding rule), against that contig's
protein-homology families (a SECONDARY reference, labelled so), by annotated-mRNA identity band; genes with ≥ 2
primary reads on their exons.

**Supplementary Fig. 6s-protein — the translated protein search as a comparator.** Human A119b chr16 (development):
the direct alignment's recovered Compara pairs by band with and without a translated protein search (mmseqs
`--search-type 2`, ≥ 30% protein identity, ≥ 50% coverage), and each tier's precision against Compara. The
translated search is not part of Rustle's family rule and is never run genome-wide (hours and ≥ 13 GB per sample).

## Methods

**Build.** `python3 figures/make.py data fig6` (genome scope, the default; `--set fig6_scope=dev` rebuilds the
development tables, with `--set families_bin=<dir>` pointing the chr16 families run at a `mcl_families` that writes
the copy table until `bin` has one). Each call does at most `fig6_budget_s` (default 540 s) of heavy work and exits 75
while work remains; repeat under `flock /mnt/linuxdisk/tmp/rustle_heavy.lock` until it exits 0. It never runs a
pipeline stage (assemblies and families come from `make.py runs`); a missing sample is listed in the table notes and
the tables are marked provisional.
1. **Spectrum** (human samples; HEAVY, resumable):
   `python3 bench/score.py spectrum --gtf <sample>.gtf --ref exons.gtf --fasta chm13v2.0.fa --chrom ALL --compara
   human_e116.tsv --out ${work}/families_gw/samples/<sample>/spectrum/genome --threads 4 --minimap2 tools/mm2_shard.sh
   --skip-t3`. Both nucleotide all-vs-alls run through the shard wrapper (3 Mb of query per shard).
2. **Families scoring:** `score.py pairs --members <sample>.fam.copies.tsv --genes exons.gtf --chrom ALL --truth
   compara:human_e116.tsv --universe <spectrum>.truth_pairs.tsv`; the paralogue-group, chromosome and
   direct-alignment attribution reuses score.py's loaders and is checked band by band against the scorer's totals
   and against spectrum.tsv.
3. **Supplement fig6s-seeding:** Compara universe (`_o1.compara_universe`: `_o1.primary_counts_gw`, one BAM pass per
   contig over every gene and pseudogene record) and `score.py pairs --universe` on both configurations' copy tables;
   Liftoff pairs from the Fig. 8 self-lift and read-support tables (`_liftoff.copy_pairs`, `pair_families`), 2,000
   bootstrap resamples of source records.
4. **Supplement fig6s-protein:** the chr16 development spectrum run with T3 (`score.py spectrum --mmseqs mmseqs`).

## Caveats

- **Genome-wide pairs include the development chromosome.** The tables also give the pairs with no chr16 gene
  (substrate S1 of `fig6_gw_recall`).
- **Pairs are not independent.** Pairs within one family or tandem array share their evidence; Wilson intervals in the
  tables treat them as independent and are descriptive only. Read the group and family counts beside every pair count.
- **Same loci, different alignments.** The direct alignment and the families use the same assembly's loci; the direct
  alignment compares spliced representatives (≥ 200 bp) and labels a locus with the RefSeq gene of largest exonic
  overlap, while the families align genomic spans and `score.py pairs` labels a copy with the gene of largest span
  overlap.
- **Estimators.** The direct alignment uses score.py spectrum's default identity (matching bases ÷ block length) and
  coverage (query span ÷ shorter length, which can exceed 1), not the family rule's.
- **Compara below about 50%.** Compara separates old paralogues by gene trees; below about 50% protein identity,
  sensitivity measures agreement with a gene-tree reference that an RNA-level definition does not claim to reproduce.
- **The seeding default was chosen against a protein referee** (`docs/archive/2026-09/PREREG_locus_read_pool_2026-09-22.md`, register
  1100, closed); the supplement rescored it on external references, and a new verdict on them needs a new
  pre-registration.
- **The legacy copy catalog is not shown.** It was the copy-assignment roster before 2026-09-25 (development tables of
  2026-09-24/25: chr16, 1,418 copies in 290 families, 85% single-exon); the families stage replaced it.
- **Human testis was aligned without `-uf`**, so its loci carry that caveat.

---

## Development tables (drawn until the genome-wide tables exist)

Human A119b chr16 (the development chromosome; Compara same-chromosome pairs; the default families built by the
families stage on the chr16 restriction of the genome-wide assembly). Tables `fig6_chr16_recall`, `fig6_chr16_groups`,
`fig6_chr16_precision`; supplements `fig6_gorilla_recall`, `fig6_gorilla_clusters` (seeding, development) and
`fig6s_protein_tiers`.

‹development numbers: filled from `python3 figures/fig_family_spectrum.py summary`›
