# Datasets — what they are, where they came from, how to rebuild them

**No sequencing data is in this repository.** The aligned BAMs are 4–96 GB each. This file records the
provenance of every substrate a result in `docs/archive/2026-11/o1_ledger.md` (the O1 ledger to 2026-10) depends on, precisely enough to rebuild it.
Paths under `/mnt/linuxdisk/...` are local to the laptop the early work was done on; treat them as labels. Since 2026-10 the work runs on the Fritz cluster: see the next section.

⚠ Two human Iso-Seq libraries are used and **they are not interchangeable** (register row 867). A result
quoted on one cannot be compared with a tool arm run on the other.

## Where things live on Fritz (2026-11-01)

| what | path |
|---|---|
| genomes | `~/data/chm13/chm13v2.0.fa`; `~/data/hg002/hg002.pat.fa` (HG002 T2T Q100 paternal, names 1..22,X,Y + JAKCWS); `~/work/jfim/rustle_data/genomes/GGO.fa` |
| Iso-Seq BAMs | `~/work/jfim/rustle_data/reads/`: `A119b.t2t.bam` (human testis, CHM13), `GCA_029281585.2_flnc_mm.bam` (KB3781 gorilla fibroblast), `GGO_mm.bam` (OR6737 gorilla testis) |
| HG002 Kinnex | `~/work/jfim/hg002_run/hg002.pat.bam` (aligned to the PATERNAL haplotype; 131 GB, 9 h to make, do not redo). `HG002_kinnex.t2t.bam` is CHM13-aligned and does not match the hg base run |
| HG002 truth | `~/data/hg002/genes.gff3`, `~/data/hg002/xref.tsv` (`DENOVO_HG002_GFF`, `DENOVO_HG002_XREF`; Ensembl 2022_07) |
| family truth | `~/data/denovo_eval/` (next-to-last section), CAT TSV `~/work/jfim/rustle_data/annotation/chm13v2.0_CAT_Liftoff.genes.tsv` |
| base runs | `~/runs/human_A119b/a1.*`, `~/runs/gorilla_KB3781/kb.*`, `~/runs/gorilla_OR6737/or.*`, `~/runs/hg002_kinnex/hg.*` (first three symlink into `~/work/jfim/runs/`); driver `tools/rustle_pipeline.sh` |
| PSV VCF report (register 1484) | `~/o2ties_real/psv/or.psv60.*` |
| storage | big outputs under `$WORK` (`~/work`; group FILE quota full since 2026-10-31, so few big files), intermediates in node-local `$TMPDIR`, kept data on `$HPCVAULT`; `$HOME` 200 GB hard |

## Reference genomes and annotations

| id | file | source |
|---|---|---|
| CHM13v2.0 | `chm13v2.0.fa` (3.1 GB) | T2T-CHM13 v2.0 |
| CHM13 RefSeq | `chm13v2.0_RefSeq_full.gff.gz` (118 MB) | RefSeq annotation of CHM13v2.0 |
| Gorilla | `GGO.fasta` | **GCF_029281585.2** (mGorGor1) |
| Gorilla RefSeq | `GGO_genomic.gff` (693 MB) | RefSeq annotation of GCF_029281585.2; 41,193 gene+pseudogene |

⚠ **Never use `Reference/HSA_genomic.gff`** — it drops 29.1% of loci. Use `chm13v2.0_RefSeq_full.gff.gz`.

Per-chromosome reference GTFs are built with the `utilities gff-to-gtf` subcommand (a stand-in for `gffread -T`,
which is not installed on the work machine; validated to reproduce gffread's `chr20_ref.gtf` transcript
count exactly, 4,574 = 4,574):

```sh
target/release/utilities gff-to-gtf chm13v2.0_RefSeq_full.gff.gz chr20 chr20_ref.gtf
target/release/utilities gff-to-gtf GGO_genomic.gff NC_073244.2 ggo_ref.gtf
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

## BAM tags and what the binaries read (2026-11-01)
All three working BAMs (`A119b.t2t.bam`, `GGO_mm.bam`, `hg002.pat.bam`) were made with the command above: `--eqx` (CIGAR `=`/`X`, so every
mismatch position and base is in the alignment; the MD and `--cs` tags would repeat it and are NOT added), `-Y` (soft clips, full query
in secondaries), `-N 50 -p 0.1 --secondary=yes` (the secondaries O2 needs). Tags present on every record: `NM`, `ms`, `AS`, `de`, `tp`,
`s1`, `cm`, `nn`, `rl`, `ts`; `s2` where a secondary chain existed. What reads what:
- `copy_assign` (assign stage): AS, CIGAR (and the 10-AS margin, `as_margin`).
- `copy_phasing` (O2 bench binary): AS, NM, CIGAR, and since register 1483 `ms` (`--tie-ms W`).
- `bench/stgaps/stratify.py`: `de` (minimap2 `-c` PAF has no `dv`).
- Unused, more informative than MD for tie decisions: `de` (gap-compressed divergence of the aligned part; equal-`ms` records can differ),
  `s1` / `s2` (chaining score of the record and of the best rival chain, minimap2's own pre-alignment confidence, what sets MAPQ; independent
  of the splice-aware AS that register 1477 showed misranks retrocopies). Adding them is a scorer column first (pre-register, run once,
  report), not an aligner change. No re-alignment is needed for anything the pipeline does today.

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
`soto_parCN_S1E.tsv`, `acro_extra_anchors.tsv`, the frozen edge tables `shared_exons_2334_finalhuman.tsv` (SEDEF
walk) and `shared_exons_5154_exon_mapback.tsv` (Soto's exon map-back, 2026-09-29), and the frozen DEV / HELD-OUT
split `soto_split_2026-09-29.tsv`. They are enough for `genesets`, `cluster` (both rules), `dennislab`, `score` and
the S1C rung of `ladder`.

| where (`/mnt/linuxdisk/home/juanfraitu/winloci_data/`) | files | used by |
|---|---|---|
| `soto_replication/` | `final_human_clean.bed`: CHM13 v2.0 SEDEF, 34 columns, 88,756 rows. The user-supplied `final_human.bed` with its trailing header row stripped (§6ip) | `edges` |
| `soto_replication/` | `final_v1.bed`: CHM13 v1.0 native SEDEF WITH CIGARs (col 33), 34 columns, 151,759 rows, identity in col 21. User-supplied 2026-09-28; strip the trailing `#chr1 ...` header row to `final_v1_clean.bed` before use. `edges --native-v1` skips the v2.0->v1.0 liftover for this file (docs/archive/2026-09/SOTO_REPLICATION_STATUS_2026-09-28.md) | `edges --native-v1` |
| `soto_replication/` | `cat_v4.bed`: CAT v4 transcripts (CHM13 v1.0, 37 columns, gene id in column 19), the BED dump of the `cat_v4.bb` beside it (chr21 spot check: the same 2,912 records) | `edges`, `curate`, `famcn --interval` |
| `soto_replication/` | `sd98_v1.bed`: Soto's SD98 regions, the merged UCSC CHM13 v1.0 `sedefSegDups` track at `$24 >= 0.98` (817 autosomal regions, 97,797,568 bp; no chrX/Y) | `curate --sd98-bed`, `edges --exon-mapback`, `famcn --interval` |
| `soto_replication/` | `t2t-chm13-v1.0.fa.gz` (+ `.fai`, `.gzi`): the CHM13 v1.0 genome (no chrY), samtools-indexed; the exon map-back's query sequences and mapping target. A prebuilt default-parameter `minimap2 -d` index (7.4 GB, 86 s / 12 GB to build) is at `/mnt/linuxdisk/tmp/rustle_figures_dev/soto_reconcile/exon/idx/t2t_v1_default.mmi` (`--mm2-index`; scratch, regenerable) | `edges --exon-mapback` |
| `soto_replication/` | `soto_{1793,2334}_geneset.tsv` (hand-made before wave 7; `genesets` now derives the same gene/biotype sets from S1C), `shared_exons_1793_final.tsv` (§6if input to `dennislab`), and the frozen outputs `replicated_families_2334_{median,mean}_finalhuman.tsv`, `replicated_families_dennislab_{mean,median}.tsv` | checks |
| `soto_replication/` | `famcn_ours_all.tsv` (2026-08-01: our WSSD famCN of the 5,154 autosomal SD98 genes over their merged CAT v4 exons, first 10 SGDP tracks, CRLF line ends; = `famcn --interval exons --samples 10`) and `famcn_ours_allwssd.tsv` (2026-09-29, sha1 11daa3ce: columns `famCN`/`famCN_mad`/`n_samples` = 268 tracks over merged exons, `famCN_269` = the outlier kept, `famCN_sotoiv`/`n_sotoiv_rows` = Soto's interval, gene body ∩ SD98; = the paste of `famcn --interval exons --samples all`, `... --outlier ''` and `--interval sd98 --samples all`, byte-identical) | `cluster --pair-mad --famcn`, `ladder --famcn-ours` |
| `soto_wssd/` | 271 per-sample `*_wssd.bb` WSSD copy-number tracks (CHM13 v1.0) = the UCSC hub listing: 269 SGDP (`LP…`/`SS…`, LP6005442-DNA_A08 among them — the sample Soto's code drops) + `chm13_wssd.bb` + `hg38_wssd.bb`; fetched by the `fetch.sh` beside them (`BASE_URL` in the module) | `famcn --wssd-dir` |

Tools: scikit-learn (optional; a stdlib ARI is used when absent), numpy and scipy for `score` / `ladder`; pyBigWig (the
miniforge python) for `famcn --interval` (`bigBedToBed` serves only the `--intervals` mode); samtools and minimap2
2.30 for `edges --exon-mapback`.

### Assembly parCN (`bench/soto/parcn_assembly.py`; `PREREG_soto_parcn_assembly_2026-09-29.md`, register 1171-1174)

| where | files | used by |
|---|---|---|
| `winloci_data/` | `chm13v2.0.fa` (soft-masked, `.fai`): the k-mer reference (chrY excluded at count time = v1.1 content); `Reference/chm13v2.0_RefSeq_full.gff.gz`: the 300 control genes | `regions`, `kmers`, `count --exclude chrY`, `edit-depth` |
| `winloci_data/` | `GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna` (chimpanzee), `GCF_028885625.2_NHGRI_mPonPyg2-v2.0_pri_genomic.fna` (orangutan): haploid primaries, counts x 2 | `count` |
| `~/gorilla_haps/` | `mat.fa`, `pat.fa` (mGorGor1 v2.0 haplotypes): counted separately, summed for the diploid famCN | `count` |
| `/mnt/linuxdisk/tmp/rustle_figures_dev/soto_parcn_asm/` (scratch; the 09-29 run) | `genomes/hg002v1.1.fasta.gz` (HG002 v1.1, MAT + PAT + chrM, 47 records); `work/regions.tsv`, `Q.u64` (17,538,361 canonical 30-mers), `pos.npz`, `cand{0,1,2}.u64`; `counts/{chm13noY,hg002,ggo_mat,ggo_pat,ptr,ppy}.i32` (exact counts per Q entry, the C counter `kc30`) and `counts/ed{0,1,2}.u32` (edit depth per candidate, `kn30`); `results_*.tsv`, `summary.json` | `analyze` (reads the counts as they are; `count` / `edit-depth` regenerate them, slowly) |

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

## DENOVOFAM truth files, re-downloaded on Fritz (2026-10-29)

`bench/denovo_eval/famtruth.py` reads them from the laptop paths unless `DENOVO_CAT_GFF`, `DENOVO_CAT_TSV`, `DENOVO_PARALOGS`
and `DENOVO_GGO_GFF` are set (plus `DENOVO_HGNC` for the HGNC gene groups). On Fritz they live in `~/data/denovo_eval/`:

| file | source | note |
|---|---|---|
| `chm13v2.0_gencode.genes.gff3` | `s3://human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3` (4.0 GB), gene lines only | the loader reads only `gene` lines (`gene_id`, `source_gene`, `collapsed_gene_ids`) |
| `chm13v2.0_CAT_Liftoff.genes.tsv` | the cluster bundle (`rustle_data/annotation/`) | unchanged |
| `paralogs.tsv` | Ensembl BioMart, release 116, dataset `hsapiens_gene_ensembl`, attributes `ensembl_gene_id`, `hsapiens_paralog_ensembl_gene`, `hsapiens_paralog_orthology_type`, `hsapiens_paralog_subtype`, one query per chromosome (`~/data/denovo_eval/fetch_paralogs.sh`) | the laptop export (3.3M rows) was an earlier release; row counts may differ |
| `GGO_genomic.gff` | the bundle's `GGO_genomic.gff.gz`, decompressed | |
| `hgnc_complete_set.txt` | `https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt` (2026-10-29) | `pairjudge.py`, `cdsseed.py` and `exon_rule_judge.py` read it from `DENOVO_HGNC` (laptop default `/mnt/linuxdisk/tmp/hgnc/`) |

## Docs archive

Closed studies: `docs/archive/2026-09/`, `docs/archive/2026-10/` (`merged/` removed 2026-11-08, see `MERGED_REMOVED.md`). `docs/archive/2026-11/` holds the verbatim pre-compaction copies (2026-11-01) of the ledgers, the register and the definition docs (`LEDGER_*.md`, `o1_ledger.md`, `NEGATIVE_RESULTS_REGISTER.md`, `seeded_family_definition.md`, `copy_assignment_definition.md`, and the living docs); every ledger citation of a dataset above resolves there in full.

## Assembly products (2026-11-02)
`PREFIX.gtf` is the assembly the later stages read (families, assign); `PREFIX.families.gtf` is it without bridge transcripts (the
families input). `PREFIX.final.gtf` is the assembly to look at and compare with other tools: the same run plus the three end layers
(median 5' ends, 3' UTR extension, single-exon models from end clusters; registers 1487-1489), produced by a parallel second run so the
families input is unchanged (registers 1490 / 1491; `docs/archive/2026-11/PREREG_assembly_ends_2026-11-01.md`). `RUSTLE_ASSEMBLY_LAYERS=off` skips it.
A fourth layer drops models whose intron chain is a strict sub-chain of another model (`--polish-subchain drop`, register 1498:
contained share 10.1% -> 5.4% vs RefSeq, complete matches -0.44%); `RUSTLE_FINAL_SUBCHAIN=off` keeps them.

## Scoring any species (2026-11-02)
`famtruth.py`, `pairscore.py`, `seqscore.py` load genes through `famtruth.load_genes`: `--species human` = CHM13 CAT (`DENOVO_CAT_GFF/TSV`);
presets `gorilla` (`DENOVO_GGO_GFF`) and `hg002` (`DENOVO_HG002_GFF` + `DENOVO_HG002_XREF`); any other name needs
`--annotation GTF|GFF3[.gz]` (format sniffed: GTF, Ensembl GFF3, RefSeq-style GFF3). `--ids symbol` (default) maps genes to human
Ensembl ids by symbol, so the human Compara table applies; `--ids own` uses the genes' own Ensembl ids with that species' Compara
paralogue export (`--paralogs`, any 4-column BioMart export) and `--level` as a comma list of taxa (`'Mus musculus,Murinae'`).
No script carries a machine path default: unset inputs exit naming the variable. Test: `bench/denovo_eval/test_famtruth_load.py`.

## Multi-sample runs (2026-11-02)
`tools/rustle_pipeline.sh all --bam a.bam,b.bam [--samples A,B] --fasta G --out PREFIX ...` (a comma in `--bam` hands the run to
`tools/rustle_samples.sh`; one BAM is unchanged). BAMs: same genome (@SQ), indexed, no shared read names. Products: PREFIX.pooled.bam,
the usual assembly / families on it, `PREFIX.samples/<S>.assign.*` per sample on the shared catalog (copy_assign
RUSTLE_SHARED_CATALOG=1: a copy silent in one sample is allowed), `PREFIX.counts.tsv` (transcript x sample read counts, a read
compatible with k transcripts adds 1/k; `utilities sample-counts`), `PREFIX.copy_counts.tsv` (copy x sample reads, abundance).
`pool` builds only the pooled BAM; `copy-table --samples A,B --out PREFIX` rebuilds the copy table. v1 refuses --candidates,
--legacy-catalog, --as-table, --to, --from and runs no flag stage. r1503 (A119b + HG002, chr17+chr19): pooled finds 29% more
exact RefSeq chains, chain precision 20.2 vs 22.0-25.9 alone: more sensitive, not more precise.

## O2 contested tables and haplotype truth (2026-11-03)
Workspace `/anvme/workspace/v138da13-rustle/` (expires 2027-01-31; `ws_extend`): `o2contested/{A119b.t2t,HG002_kinnex.t2t,GGO_mm,KB3781}.molecules.tsv` (`as_table` with columns 7-17: best/second/primary ms, MAPQ, splice shape, best-AS rival, positions; one `#as_table` line, no column row), `o2contested/KB3781.records.tsv` (every primary and secondary record: name, secondary, chrom, pos, end, AS, ms, spliced; 34.8M records, 2.2 GB) and its by-name sort, `o2contested/haps/KB3781.{mat,pat}.mmi` (Jim's haplotypes for `haplotype_truth.py`), runs `o2mstie/` (r1509/1510/1514 assignments, `truth_all10/hits.tsv`, `*.strict_dump.tsv`), `o2contested/r1511/`, `o2pool/` (r1513 arms, `GGO.splice.mmi`). KB3781 fibroblast reads are Jim's (the mGorGor1 individual), so his haplotype assemblies are the read-level truth. Sbatch files and logs: `~/o2contested/`.

## SEDEF SD calls on Fritz (2026-11-08; register 1580-1581, `--sd-damp-bed`)

All under `/anvme/workspace/v138da13-rustle/` (workspace expires 2027-01-31; `ws_extend`). SEDEF BEDs have 34 columns:
pair coordinates in 1-6, identity (`fracMatch`) in 21, CIGAR in 33, and a trailing `#chr1 start1 ...` header line, which
`SdPairs::from_bed_str` skips.

**Families default (user decision 2026-11-09):** SD damping is on whenever a genome's SEDEF BED is given:
`tools/rustle_pipeline.sh` with `RUSTLE_SEDEF_BED=<file below>` (factor `RUSTLE_SD_DAMP`, default 0.25; 1 = off), or with
the file linked as the sidecar `<FASTA>.sedef.bed`, used when the variable is unset (`RUSTLE_SEDEF_BED=` turns it off;
CHM13's is linked at `~/data/chm13/chm13v2.0.fa.sedef.bed`), or
`mcl_families --sd-damp-bed <file>` directly (`--sd-damp` defaults to 0.25 with a BED). Human CHM13 runs use
`final_human.bed`, gorilla runs `gorilla/final.bed`, HG002 runs `hg002pat/final.bed`.

| genome | file | provenance |
|---|---|---|
| CHM13 v2.0 | `sedef_dl/incoming/sedef_send/final_human.bed` (88,756 pairs, md5 2cd7d079...) | user-supplied (croc 2026-11-07), SEDEF v1.1 `sedef.sh` defaults on soft-masked chm13v2.0.fa (Penn State Roar); `sedef_send/sedef_provenance.txt` |
| CHM13 v2.0 (UCSC) | `sedef_dl/chm13v2.sedef.bed` (83,516 rows; identity col 24, partner cols 10-12) | UCSC hub bigBedToBed; NOT the same call set as above (the native file is the one used) |
| gorilla mGorGor1 v2.0 (GCF_029281585.2, KB3781) | `sedef_dl/gorilla/final.bed` (253,029 pairs, md5 991ed343...) | user-supplied, SEDEF v1.1 on GGO.fasta; contigs match `~/work/jfim/rustle_data/genomes/GGO.fa` |
| HG002 paternal (`~/data/hg002/hg002.pat.fa`) | `sedef_build/hg002pat/final.bed` | built here 2026-11-08: SEDEF 1.1 (github vpc-ccg/sedef tag 1.1) rebuilt per node in `$TMPDIR`, `sedef_build/sedef_notime.sh` (= `sedef.sh` with `/usr/bin/time`, `--bar` and the awk command generation removed, and `grep -r`/`find` instead of globs over 264,714 seed logs); seeding 56 min on 72 cores; 41 align buckets OOM-killed at 72-way and rerun 8-way (`job_retry_killed.sh`, `job_resume.sh`). HG002 is soft-masked less than CHM13 (50 Mb slice 37% vs 59% lowercase). This build's rows carry one extra field and no header line: identity (`fracMatch`) is column **22**, CIGAR 34 (pair coordinates still 1-6). 714,007 pairs; 4 chr9-satellite candidates > 1-5 Mb dropped (row 1580 amendments 2-3) |

## Repeat and protein annotation for tie groups (2026-11-10; register 1607, 1611)
Under `/anvme/workspace/v138da13-rustle/annot/`:

| file | source | used by |
|---|---|---|
| `dl_rmsk/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed` | T2T CHM13 RepeatMasker | `tie_groups.py build ... RMSK_BED:0.5` (human) |
| `dl_rmsk_ggo/GGO.repeatMasker.bed` | UCSC GenArk `GCF_029281585.2.repeatMasker.out.gz`, converted (NCBI ships no `rm.out` for this assembly) | same, gorilla |
| `dl_gencode/gencode.v47.pc_translations.fa.gz` | GENCODE v47 protein translations | `bench/famcmp/protein_truth.py` -> `r1611/ens_map.tsv` (`DENOVO_GGO_ENS_MAP`) |

Run directories: `r1603/` (A119b arms, `recs.bam` 26.6 GB, `recs.slim.tsv`), `r1610/`, `r1611/` (gorilla OR6737). SEDEF files: see the section above.
