# Figure 7 — Family recovery by Rustle's two modes, de novo and guided (a Rustle-internal comparison)

> **Status (2026-09-25).** The genome-wide version is pre-registered (`docs/PREREG_genome_wide_families_2026-09-25.md`,
> claims F7.1–F7.4 as amended by Amendment 1, written before any genome-wide family number) and not built yet: it
> needs the families stage of every sample (with its copy table), the guided run of every species, and, for the
> Liftoff rows, the Fig. 8 self-lift and read-support tables. Until its tables exist (`fig7_gw_summary`,
> `fig7_gw_contigs`, `fig7_gw_per_family`), `make.py plot fig7` draws the per-chromosome development tables
> summarised at the end of this caption. After the build, `python3 figures/fig_family_recovery.py summary` prints
> every number marked ‹…›; replace them, state each claim as the pre-registration's table says, and delete this box.

**This is not a comparison with other tools.** It compares the two modes of Rustle with each other. Both modes group
loci with the same family rule; they differ in where the loci come from. The de novo mode is Rustle's one default de
novo family definition, the families copy assignment uses.

**Claim (pre-registered; wording decided by the outcome).**
- F7.1: guided recovers the Ensembl Compara families (duplications within primates) at least as well as de novo, in
  both human samples (genome minus development contigs): ‹bipartite F, de novo vs guided, per sample›.
- F7.2: the sign of the guided − de novo gap is the same without every contig used for a family decision: ‹…›.
- F7.3 (descriptive): the de novo families place ‹k of n› read-supported Liftoff copy pairs in one family, per sample
  (every species; the guided mode cannot hold an unannotated copy).
- F7.4: no single contig carries the gap: ‹…›.

## Definitions (each used once below)

- **De novo mode = Rustle's default de novo families.** IsoSeq reads → loci assembled from primary alignments plus
  secondary alignments that score at least 98% of the read's best alignment score anywhere in the genome → one
  representative per locus (the transcript with the most reads; its exons, with the genome's bases) → the family rule
  below, on the loci's genomic spans (`minimap2 -x asm20 -c -X -N 50 -p 0.1 --secondary=yes`, genome-wide, so a family
  may span chromosomes). One run per sample.
- **Guided mode.** Loci are the annotated gene and pseudogene bodies of the species' RefSeq annotation (loci whose
  exons overlap are folded into one). They are aligned all-vs-all **with the same minimap2 flags as de novo**, so the
  two modes differ only in their loci. Guided reads no RNA: one run per species, drawn as the lower dot on every
  sample of that species.
- **Family rule (both modes).** Two loci are linked when an alignment of ≥ 300 bp at ≥ 70% identity covers ≥ 30% of
  the longer locus's exonic length, touches exonic bases on both sides, and one record covers ≥ 60% of the smaller
  locus's exonic length exon to exon; Markov clustering (MCL, inflation 2.8); a predicted family has ≥ 2 loci. No
  protein is used.
- **Reference families** (external; how independent each is, stated with it).
  - *Ensembl Compara families* (human; release 116; the headline): protein-coding genes connected by Compara paralogue
    pairs whose duplication node lies within the primates (Homo sapiens … Primates in Compara's species tree). Within
    one Compara gene tree this is a clade below a duplication, not a chain of similarity; no identity threshold.
  - *Soto et al. 2025* (human; Table S1C, segmental-duplication gene families; each gene's first family). **Not
    independent:** the 60% shared-exon threshold of the family rule was chosen against these families.
  - *NPIP reference set* (human chr16 only; an inset): Soto's NPIP families extended with RefSeq NPIP genes and genes
    aligning over ≥ 95% of their length to a Soto NPIP gene.
  - *Liftoff copy pairs* (every species; genome-wide version): the Fig. 8 self-lift pairs each annotated gene or
    pseudogene with the unannotated copies (≥ 95% identical, exons ≥ 200 bp) Liftoff finds for it in the same genome;
    a pair is scored when both loci have ≥ 2 reads of the sample whose primary alignment has an aligned block on the
    locus's exons, and recovered when one de novo family's loci cover ≥ 50% of both loci's exon bases (Liftoff's own
    criterion). Pairwise sensitivity only (the relation certifies pairs; it is not a partition). **The guided mode is
    not scored on it:** its loci are the annotation, and an extra copy is unannotated by construction.
  - The annotation's protein-homology families are a SECONDARY reference, in Supplementary Fig. 7s only.
- **Scoring** (`family_score`). A locus is labelled with the one annotated gene it overlaps most. Predicted families
  are intersected with the genes of the reference families; genes without a reference family are not scored, so
  **precision is an upper bound**, the same bound for both modes.
  - *One-to-one (bipartite) matching:* each reference family is matched to at most one predicted family, maximising
    the shared genes. Sensitivity = matched genes ÷ reference genes; precision = matched genes ÷ genes of the matched
    predicted families; **F (pooled)** = their harmonic mean.
  - *Pairwise:* sensitivity = same-family gene pairs placed in one predicted family ÷ all same-family pairs;
    precision = those pairs ÷ all within-family predicted pairs.
  - *Exact / partial / missed:* a reference family whose matched predicted family equals it (F = 1) / overlaps it /
    shares no gene with it.
- **Substrates.** One genome-wide run per mode, restricted afterwards (a reference family keeps its members on the
  kept contigs, a predicted family its loci; a family with fewer than 2 members left is not scored):
  - *genome minus development contigs* (the headline; human chr16, gorilla chr20 = NC_073244.2);
  - *genome minus every contig used for a family decision* (the ring in panel b; human chr16, chr5, chr7, chr21,
    chr2, chr8, chr10; gorilla chr20);
  - *per contig* (panel d): the genome-wide run restricted to one contig. This is a breakdown, not a per-chromosome
    run: a family built genome-wide may have members elsewhere.
- **Contig exposure** (human; the gorilla development contig is chr20, which chose the seeding rule):
  development = chr16 (the early family rules); threshold selection = chr5, chr7, chr21 (the 60% exon threshold was
  chosen there against Soto families); reused verdict set = chr2, chr8, chr10 (about 30 guided-mode tests since
  2026-09-20); scored once, no decision = human chr6 and gorilla chr10 (NC_073234.2; the per-chromosome figure);
  every other contig was never used for a family decision. Human testis, gorilla KB3781, chimpanzee and orangutan
  were never used for any family decision.

## Panels

Samples: human A119b (T2T-CHM13 v2.0; tissue not recorded), human testis (T2T-CHM13 v2.0; public library, aligned
without `-uf`), gorilla OR6737 testis and gorilla KB3781 fibroblast cell line (mGorGor1), chimpanzee (mPanTro3) and
orangutan (mPonPyg2). Species are never pooled.

**a–c** One row per sample and reference, grouped by species; upper dark dot = de novo (that sample), lower light dot =
guided (that species' annotation; the same value on both samples of a species). Genome minus development contigs. The
NPIP reference set is an inset on chr16, the development chromosome. Liftoff rows carry the de novo pairwise
sensitivity only (panel c gives k of n pairs).
- **a** Pairwise sensitivity and precision.
- **b** One-to-one sensitivity, precision and F (pooled). In the F column, the ring is the same F without every contig
  used for a family decision.
- **c** Reference families, exact · partial · missed, per mode.

**d** Bipartite F per contig against the Compara families (human samples and the human guided run). Marker shape:
triangle = development, diamond = threshold selection or reused verdict set, open circle = scored once, filled circle
= never used for a family decision. The black tick is the genome-minus-development value.

**e–g** Per-family F, de novo (x) against guided (y), genome minus development contigs; one point per reference
family that at least one mode recovers (jitter ≤ 0.015): **e** human A119b, Compara families; **f** human testis,
Compara families; **g** human A119b, Soto 2025 (not independent). Flagship families are named.

**Numbers shown** (‹…› = from `fig_family_recovery.py summary` after the build):

| sample | reference | families · genes | sensitivity (dn / g) | precision (dn / g) | F (dn / g) | F gap (g − dn) | exact·partial·missed (dn / g) |
|---|---|---|---|---|---|---|---|
| human A119b | Compara, primates | ‹…› | ‹…› | ‹…› | ‹…› | ‹…› | ‹…› |
| human A119b | Soto 2025 | ‹…› | ‹…› | ‹…› | ‹…› | ‹…› | ‹…› |
| human testis | (same two) | ‹…› | ‹…› | ‹…› | ‹…› | ‹…› | ‹…› |
| every sample | Liftoff copy pairs | ‹pairs› | ‹de novo only› | – | – | – | – |
| human, chr16 inset | NPIP reference set | ‹…› | ‹…› | ‹…› | ‹…› | ‹…› | ‹…› |

**Supplementary Fig. 7s-protein-homology.** The same panels against the annotation's protein-homology families (one
protein per gene, the longest CDS; pseudogenes and immunoglobulin / T-cell-receptor gene segments excluded; all-vs-all
BLASTP, E ≤ 10⁻⁵; genes linked when non-overlapping aligned segments cover ≥ 30% of the longer protein; MCL inflation
2.8), a SECONDARY reference: most of their same-family pairs have no nucleotide alignment, they mix fold-sharing genes
with paralogues, and they exclude pseudogenes. Genome-wide only with `fig7_protein_homology 1`; the development
version draws the per-chromosome protein-homology rows (human and gorilla).

## Methods

**Build.** `python3 figures/make.py data fig7` (genome scope, the default; `--set fig7_scope=dev` rebuilds the
development tables, `--set fig7_dev_cached=1` rescores the cached per-contig runs without re-running them). Each call
does at most `fig7_budget_s` (default 540 s) of heavy work and exits 75 while work remains; repeat under `flock
/mnt/linuxdisk/tmp/rustle_heavy.lock` until it exits 0. It never runs a pipeline stage: the de novo families are the
run-cache stage `families` (`make.py runs --sample S --stage families`, with `RUSTLE_MINIMAP2=tools/mm2_shard.sh`); a
missing sample is listed in the table notes and the tables are marked provisional.
1. **Annotation caches** (one pass per species GFF, about 25–40 s): gene records, their exon unions, a
   gene-records-only GFF for the scorer (`${work}/families_gw/species/<species>/`).
2. **Guided** (per species, HEAVY, resumable): every gene / pseudogene record → `samtools faidx -r` →
   `tools/mm2_shard.sh paf OUT -x asm20 -c -X -N 50 -p 0.1 --secondary=yes -t 4 BODIES BODIES` (shards of 10 Mb of
   query, one shared index; `cmp`-identical to one minimap2 run, checked on gorilla chrY) → `mcl_families --paf OUT
   --gff <full GFF> --min-exonic-bp 1 --min-shared-exon-frac 0.60`.
3. **References:** Compara families from `bench/truth.py compara` (`compara_gw`) and the CHM13 RefSeq protein-coding
   genes; Soto S1C and the NPIP set rewritten with each gene's contig; Liftoff pairs from the Fig. 8 tables
   (`_liftoff.copy_pairs`, `pair_families`, on the sample's `<id>.fam.copies.tsv`).
4. **Scoring:** `family_score --clusters C --gff genes_only.gff --soto REF --chrom ALL --per-family P --pairwise`, on
   each substrate and each contig; pooled values are recomputed from the per-family rows and checked against the
   printed line. Units are cached (`${work}/fig7/gw/`).
5. **Bridging measurement** (pre-registered; `fig7_gw_bridge`, supplementary, decides nothing): on human chr16 and
   gorilla chr20, per-chromosome guided families with the recorded flags and with the de novo flags, scored against
   the contig's own protein-homology families and (human) Soto 2025.
6. **Compara level sweep** (`fig7_gw_compara_levels`, table only): Compara families at Hominidae, Catarrhini,
   Primates (the headline), Eutheria, Vertebrata and all levels.

## Caveats

- **Precision is an upper bound** in both modes: genes without a reference family are not scored.
- **Sensitivity mixes read coverage, alignability and the definition.** Guided loci are the annotation's own genes,
  the same genes that label every reference, so guided has a representational advantage by construction. The gap is
  node coverage plus definition; the figure does not separate the two.
- **Liftoff copy pairs score de novo only.** They measure what the annotation lacks; they are not a mode comparison.
- **Soto 2025 is not independent** of the 60% exon threshold, and it is a cover (6.4% of its genes carry more than
  one family; the first is kept).
- **Development and reused contigs** are included in "whole genome" and removed in the headline; the stricter
  substrate removes every contig used for a family decision. The rules were developed on human families (NPIP,
  TBC1D3, AMY, the Soto families) that have homologues in the apes. Human chr6 of A119b was scored in the
  per-chromosome development tables (below), so the Compara claims must also hold on human testis.
- **Human testis was aligned without `-uf`** (strand not forced); its de novo loci carry that caveat.
- **Pairwise numbers are dominated by the largest families** (pairs grow quadratically); the bipartite columns weight
  families by genes.
- **Display only.** Per-family labels are the longest common prefix of a family's named members; points are
  jittered; flagship labels are placed so that no label covers a point or another label.

---

## Development tables (drawn until the genome-wide tables exist)

Per-chromosome runs of both modes (2026-09-25; the genome-wide A119b assembly restricted to each chromosome for de
novo, the chromosome's annotation for guided, with the recorded guided flags `-x asm20 -c --eqx -P`; families within
one chromosome only): human A119b chr16 (development), chr2, chr8, chr10 (reused verdict set), chr6 (scored once, no
decision). Rescored 2026-09-25 from the cached runs (`fig7_dev_cached 1`, not re-run) against the Compara families at
Primates restricted to each chromosome, Soto 2025 and the NPIP set (tables `fig7_summary`, `fig7_per_family`). The
main development figure has no gorilla rows: the gorilla contigs' only development reference is the secondary one
(Supplementary Fig. 7s), and Liftoff copy pairs exist only in the genome-wide version. These are development numbers,
not the pre-registered claims.

- Compara families (primates), one-to-one F, de novo vs guided: chr6 (never used before) 0.444 vs 0.833 (+0.389; 15
  families, 35 genes; exact 3 vs 9); chr2 0.588 vs 0.785 (+0.197); chr8 0.364 vs 0.667 (+0.303); chr10 0.542 vs 0.642
  (+0.100); chr16 (development) 0.615 vs 0.875 (+0.260). Guided was higher on all five chromosomes.
- Soto 2025 (not independent): +0.178 (chr2), +0.274 (chr8), +0.193 (chr10), +0.212 (chr16); equal on chr6 (0.800,
  3 two-gene families). NPIP set (chr16): 0.576 vs 0.721 (+0.145).
- Per family, pooled over the five chromosomes (panel d): against Compara, guided higher in 36 families, de novo
  higher in 4, equal in 23 (16 missed by both); against Soto 2025, 33 / 7 / 23 (12 missed by both).
- Bipartite precision, an upper bound, was 0.82–1.00 in every Compara and Soto row (NPIP set: 0.68 de novo, 0.82 guided).
- Supplementary Fig. 7s (protein-homology families, secondary; unchanged from the earlier development tables):
  guided − de novo +0.081 on chr6 (0.115 vs 0.196), +0.047 to +0.260 on chr2, chr8, chr10, +0.087 on chr16; gorilla
  chr10 (NC_073234.2) +0.076, chr20 (NC_073244.2, development) +0.122; both modes missed most protein-homology
  families (321 of 469 human, 178 of 212 gorilla families missed by both).
