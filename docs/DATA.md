# Datasets — what they are, where they came from, how to rebuild them

**No sequencing data is in this repository.** The aligned BAMs are 4–96 GB each. This file records the
provenance of every substrate a result in `docs/o1_ledger.md` depends on, precisely enough to rebuild it.
Paths under `/mnt/linuxdisk/...` are local to the machine the work was done on; treat them as labels.

⚠ Two human Iso-Seq libraries are used and **they are not interchangeable** (register row 867). A result
quoted on one cannot be compared with a tool arm run on the other.

## Reference genomes and annotations

| id | file | source |
|---|---|---|
| CHM13v2.0 | `chm13v2.0.fa` (3.1 GB) | T2T-CHM13 v2.0 |
| CHM13 RefSeq | `chm13v2.0_RefSeq_full.gff.gz` (118 MB) | RefSeq annotation of CHM13v2.0 |
| Gorilla | `GGO.fasta` | **GCF_029281585.2** (mGorGor1) |
| Gorilla RefSeq | `GGO_genomic.gff` (693 MB) | RefSeq annotation of GCF_029281585.2; 41,193 gene+pseudogene |

⚠ **Never use `Reference/HSA_genomic.gff`** — it drops 29.1% of loci. Use `chm13v2.0_RefSeq_full.gff.gz`.

Per-chromosome reference GTFs are built with the `gff_to_gtf` binary (a stand-in for `gffread -T`,
which is not installed on the work machine; validated to reproduce gffread's `chr20_ref.gtf` transcript
count exactly, 4,574 = 4,574):

```sh
target/release/gff_to_gtf chm13v2.0_RefSeq_full.gff.gz chr20 chr20_ref.gtf
target/release/gff_to_gtf GGO_genomic.gff NC_073244.2 ggo_ref.gtf
```

## Iso-Seq libraries

### 1. `human_testis.t2t.bam` — the six-chromosome panel (§6p8-§6q6)

Public human testis Iso-Seq. **Read names carry the ENA run accession: `ERR13885926`.** 188,864 records
on chr20.

```sh
minimap2 -ax splice:hq --eqx -Y -N 50 -p 0.1 --secondary=yes -t 8 \
         chm13v2.0.splice.mmi human_testis.fastq.gz | samtools sort -@4 -o human_testis.t2t.bam
samtools index human_testis.t2t.bam
```
(minimap2 2.30-r1287, samtools 1.22.1. ⚠ no `-uf` — unlike the two libraries below.)

### 2. `A119b.t2t.bam` — the deep human library the lab's tool arms use (§6q7)

96 GB, 68,026,217 alignment records, 9,841 unmapped; **1,104,846 records on chr20, ~6× library 1**.
PacBio movie `m64404e_240606_134033`; lab-internal, so ⚠**an accession still has to be filled in here
before anyone outside the group can rebuild it.** Pipeline: `ccs` → `skera split` → `lima` → `isoseq
refine --require-polya` (giving `A119b.flnc.bam`) → minimap2:

```sh
minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 28 chm13v2.0.fa - \
  | samtools sort -o A119b.t2t.bam ; samtools index A119b.t2t.bam
```
(minimap2 2.31-r1302. `A119b.flnc.bam` + its `.pbi` are what `isoseq collapse` needs for FL counts.)

### 3. `GGO_mm.bam` — gorilla, the thesis substrate (§6q7)

11.7 GB. FLNC reads `GGO_OR6737.IS.8a36218d3f23-filtered.fastq` against GCF_029281585.2; 473,231 records
on `NC_073244.2`. ⚠**Accession to be filled in.**

```sh
minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes GGO.fasta \
         GGO_OR6737.IS.8a36218d3f23-filtered.fastq -o GGO_mm.sam
samtools sort -o GGO_mm.bam GGO_mm.sam ; samtools index GGO_mm.bam
```
(minimap2 2.31-r1302.)

⚠ `GGO_ds.bam` is a **downsample** of a gorilla library and is NOT the same substrate as `GGO_mm.bam`.

## Comparison-tool outputs (not produced by this repo)

`~/Desktop/isoseq_upload/` (isoseq collapse, A119b + GGO_OR6737) and `~/Desktop/benchmark_collapse/`
(StringTie + FLAIR, both samples) hold cluster-produced, genome-wide transcript sets, with the exact
sbatch recipes beside them. §6q7 reuses these **as-is**. Their StringTie arms are 3.0.1.

## Layer-order / nested-lattice substrate (NPIP/TBC1D3, human CHM13, 2026-09-16)

`bench/layer_order/npip_tbc1d3.py` (library `lattice_common.py`) reproduces `bench/LAYER_ORDER_NPIP_TBC1D3.md` and
`bench/NESTED_LATTICE_NPIP_TBC1D3.md`. The results tree is `ROOT = /mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3`
(`--root` / `LO_ROOT`). Every stage **overwrites** its outputs under ROOT, so reruns go to a copy. The off-repo
`light/scripts/` and `heavy/scripts/` built the input tables (see the README in each); this repo does not rebuild them.
The code reads the following.

| where | files | used by |
|---|---|---|
| `ROOT/light/work/refseq/` | `genes.tsv`, `exons.tsv`, `cds.tsv`, `gene_dbxref.tsv` (RefSeq CHM13 gene records) | every stage |
| `ROOT/light/work/P/` | `proteins.index.tsv`, `blastp.tsv`, `searched.txt` (§6ko protein layer) | corrected-tables, lattice-edges |
| `ROOT/light/work/D/`, `work/S1/` | `c15_17_22.e1.*`, `c16_19_20.e1.*` (`graph.tsv`, `loci.tsv`, `clusters.tsv`); `*.e1s.graph.tsv` | corrected-tables, lattice-edges |
| `ROOT/light/work/cat/` | `cat_genes_exons.tsv` (CAT v2.0 exon unions; the Soto gene-id mapping) | Soto mapping |
| `ROOT/light/` | `members.tsv`, `P.groups.tsv`, `P.members_status.tsv`, `P.edges.tsv`, `D.groups.tsv`, `D.edges.tsv`, `C.groups.tsv`, `C.supported_clades.tsv`, `truth_soto.tsv`, `truth_soto_families.tsv`; plus the `*.corrected.tsv` that corrected-tables writes | corrected-tables, layer-order, lattice-edges |
| `ROOT/heavy/` | `EXPR.counts.tsv`, `S2.genes.tsv`, `S2.edges.tsv` | expr-recount, corrected-tables, lattice-edges |
| `ROOT/integrate_slim/` | `P_N2_clusters.tsv`: a **frozen input**, written by the archived `lo_p_variants.py` (tag `notebook-2026-09-20`) with the off-repo `light/scripts/layer_protein_bounded.py` | layer-order |
| `ROOT/lattice/pre_correction_1712/` | `nodes.tsv` (optional; only for the "V identical to the 17:03 build" log line) | lattice-edges |
| `/mnt/linuxdisk/home/juanfraitu/o1_falsemerge/` | `human2/genes.asm20.paf`, `human2/genes.regions`, `human2/guided.*`; `lit/aj_dev/refseq_e1.*`; `lit/aj_ho/refseq/{all.paf, nodes.tsv, nodes.tsv.names.tsv, e0.*, e1.*}` (the two E1 catalogs: gene-body PAFs, keys, E0/E1 loci and clusters) | corrected-tables, layer-order, lattice-edges, lattice-check-c2 |
| `winloci_data/` | `hgnc/hgnc_complete_set.txt`; `soto_replication/soto_gene_to_families.tsv`; `Reference/chm13v2.0_RefSeq_full.gff.gz` | truths, Soto mapping, expression |
| `_from_wsl/human_val/human_testis.t2t.bam` | the testis Iso-Seq library above (library 1) | expr-recount, lattice-expr |
| repo | `docs/lit_subclusters_npip_tbc1d3_truth.tsv`, `bench/soto/soto_famCN_S1C.tsv` | corrected-tables, Soto mapping |

Tools: `samtools` (on PATH); numpy and scipy (the scorers); the `mcl_port` Rust bin, called through `bench/lib.py`
(layer-order REFINE; `RUSTLE_MCL_PORT_BIN`); `bench/guided_pipeline.py` (`gene_body_chains`, lattice-check-c2);
`bench/truth.py` (`excluded`, `edges_from`, `pair_hsps`).

## Soto 2025 replication substrate (ledger §6ie–§6ip; `REPRODUCE.md` §5a)

`bench/soto/soto_replication.py` reads these. The in-repo inputs are `bench/soto/soto_famCN_S1C.tsv`,
`soto_parCN_S1E.tsv`, `acro_extra_anchors.tsv` and the frozen edge table `shared_exons_2334_finalhuman.tsv`. They are
enough for `genesets`, `cluster`, `dennislab` and `score`.

| where (`/mnt/linuxdisk/home/juanfraitu/winloci_data/`) | files | used by |
|---|---|---|
| `soto_replication/` | `final_human_clean.bed`: CHM13 v2.0 SEDEF, 34 columns, 88,756 rows. The user-supplied `final_human.bed` with its trailing header row stripped (§6ip) | `edges` |
| `soto_replication/` | `cat_v4.bed`: CAT v4 transcripts (CHM13 v1.0, 37 columns, gene id in column 19), the BED dump of the `cat_v4.bb` beside it (chr21 spot check: the same 2,912 records) | `edges` |
| `soto_replication/` | `soto_{1793,2334}_geneset.tsv` (hand-made before wave 7; `genesets` now derives the same gene/biotype sets from S1C), `shared_exons_1793_final.tsv` (§6if input to `dennislab`), and the frozen outputs `replicated_families_2334_{median,mean}_finalhuman.tsv`, `replicated_families_dennislab_{mean,median}.tsv` | checks |
| `soto_wssd/` | 271 per-sample SGDP `*_wssd.bb` WSSD copy-number tracks (CHM13 v1.0), fetched from the UCSC hub (`BASE_URL` in the module) by the `fetch.sh` beside them | `famcn --wssd-dir` |

Tools: scikit-learn, numpy and scipy for `score`, and pyBigWig (the miniforge python) or `bigBedToBed` for `famcn`.

## Derived working substrates

Built under `/mnt/linuxdisk/home/juanfraitu/bakeoff/`, one directory per chromosome, each holding
`chrN.bam`, `chrN.fa`, `chrN_ref.gtf` and the tool arms:

| directory | library | note |
|---|---|---|
| `human_chr{20,11,7,14,5,9}` | `human_testis.t2t.bam` | the §6p8-§6q6 six-chromosome panel |
| `a119b_chr20` | `A119b.t2t.bam` | §6q7 human |
| `ggo_NC_073244.2` | `GGO_mm.bam` | §6q7 gorilla |

```sh
samtools view -@3 -b <library>.bam <contig> -o chrN.bam && samtools index chrN.bam
samtools faidx <genome>.fa <contig> > chrN.fa && samtools faidx chrN.fa
```
