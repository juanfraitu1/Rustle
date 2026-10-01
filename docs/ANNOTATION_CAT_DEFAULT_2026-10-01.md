# Default human annotation: CAT/Liftoff on CHM13 v2.0 (2026-10-01)

**Decision (user, 2026-10-01):** human work uses the T2T-CHM13 **v2.0 CAT/Liftoff** annotation from now on, the annotation family
Soto et al. 2025 used (CAT/Liftoff v4 on v1.0 for genes; the v2.0 CAT/Liftoff transcriptome for RNA) and the one the T2T ape papers
combine with NCBI. Results already registered against RefSeq stay as registered; the key benchmarks are re-run on CAT and both are
reported. Gorilla keeps NCBI RefSeq (`GGO_genomic.gff`) until a CAT annotation of our gorilla reference is checked.

## Files (`/mnt/linuxdisk/home/juanfraitu/winloci_data/gencode_chm13/`)

| file | what |
|---|---|
| `chm13v2.0_gencode.gff3` | the source: 64,213 genes (60,239 CAT + 3,974 Liftoff), 234,903 transcripts, 1,444,511 exons, 4 GB (long attributes) |
| `chm13v2.0_CAT_Liftoff.slim.gff3.gz` (+ `.tbi`) | **use this**: gene / transcript / exon lines, sorted, tabix-indexed, 17 MB. Gene lines `ID=Name=<gene_id>` (unique; CAT names repeat across paralogs) + `gene_name=`; exon lines `Parent=<transcript>;gene=<gene_id>` (what our GFF readers attach exons by) |
| `chm13v2.0_CAT_Liftoff.genes.tsv` | one row per gene: id, name, biotype, CAT/Liftoff, span, strand, transcripts, exon-union blocks |
| `chm13v2.0_CAT_Liftoff.vs_v4.tsv` | per CAT v4 gene (v1.0): present in v2.0, same exon structure |
| `chm13v2.0_CAT_Liftoff.refseq_map.tsv` | per RefSeq gene: the CAT gene sharing the most exonic bases on the same strand (ties: Jaccard), shared bp, quality |

Built by `bench/annotation/cat_setup.py` (73 s). Never join the two annotations by name: CAT's NPIPB3 is RefSeq's NPIPB5.

## Checks

- **v2.0 vs v4 (Soto's annotation):** of 62,671 v4 genes, 62,617 have the same exon structure in v2.0, 1 differs, 53 are absent;
  1,595 genes are new in v2.0. **All 2,334 of Soto's genes are present with identical exon structure**, so the Soto replication
  carries over unchanged.
- **RefSeq -> CAT (58,563 RefSeq genes):** strong 28,654, partial 8,799, weak 1,810, none 19,300. The unmatched are mostly pseudogenes
  (11,119) and lncRNAs (6,835); 214 protein-coding. Protein-coding: strong 15,989, partial 3,650, weak 236, none 214.
- **NPIP / TBC1D3 RefSeq copies (41):** strong 18, partial 15, weak 1, none 7. Truth sets built from RefSeq copies need translating
  copy by copy before a CAT re-run.

## Coordinates

CAT v4 is CHM13 v1.0; this annotation, RefSeq, our reads and indexes are v2.0. v1.0 -> v2.0 differs on every chromosome (by sequence,
2026-09-30: chr1 -169, chr9 +8, chr16 -5, chr17 -291, ...; acrocentrics by up to 737 kb from the rDNA rebuild). Use v2.0 throughout.

## Results registered against RefSeq (re-run on CAT, report both)

The O1 NPIP / TBC1D3 truth (26 human NPIP copies), the copy-recovery comparison (`copy_recovery_tools`, truth from RefSeq), the
held-out family sets (PREREG 2026-09-20), the SQANTI3 polish table, and the layer-order / lattice tables (`light/work/refseq/`).
