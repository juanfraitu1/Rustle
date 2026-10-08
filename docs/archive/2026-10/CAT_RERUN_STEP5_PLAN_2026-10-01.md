# CAT re-run, steps 5 and 6: plan (written before any CAT result was computed)

2026-10-01. Follows `docs/archive/2026-10/CAT_RERUN_PROTOCOL_2026-10-01.md` (rulings R1-R6 and Amendment 1, binding; R6: nothing is re-tuned
on CAT) and `docs/archive/2026-10/ANNOTATION_CAT_DEFAULT_2026-10-01.md`. Step 5 rebuilds the NPIP/TBC1D3 layer-order / nested-lattice study
(`bench/LAYER_ORDER_NPIP_TBC1D3.md`, `bench/NESTED_LATTICE_NPIP_TBC1D3.md`) with CAT/Liftoff v2.0 in place of RefSeq at every
input. Step 6 recomputes the family certificates (`bench/FAMILY_CERTIFICATES_NPIP_TBC1D3.md`, D1) on CAT evidence. The report
comparing both is `docs/LAYER_ORDER_CAT_2026-10-01.md`.

## Where things go

- `ROOT` = `/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3` (frozen, never written). `CAT` =
  `/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3_cat`, a copy of `ROOT/{light,heavy,integrate_slim,lattice}`; every
  input below is overwritten inside `CAT`. The file names stay those the code reads (e.g. `CAT/light/work/refseq/genes.tsv`
  holds CAT genes; the directory name `refseq` is kept because the code hardcodes it).
- Off-repo scripts get `_cat` copies next to the copied originals, inside `CAT/light/scripts/` and `CAT/heavy/scripts/`; the
  originals in `ROOT` are untouched. CAT logic shared by them lives in one new repo file, `bench/annotation/cat_layer_order.py`.
- Step 6 works in `/mnt/linuxdisk/home/juanfraitu/family_cert_cat/` (new; `_cat` copies of the `family_cert/` scripts).
  `family_cert/` is never written.
- Repo code (`bench/layer_order/lattice_common.py`, `npip_tbc1d3.py`): CAT behaviour only behind `LO_ANNOT=cat`
  (default `refseq` = unchanged). Proof of no RefSeq change: `npip_tbc1d3.py all --with-check-c2` run on a fresh copy of
  ROOT with the code BEFORE the edit and on another fresh copy AFTER it; every output file compared byte for byte (timing
  tokens excepted, as in the wave-7 check).

## Amendment 2 of the protocol (applied everywhere below)

Protocol Amendment 2 (commit 695be9eb, received while this plan was being written, before any CAT result): for a RefSeq
read-through record `A-B` that a truth table names as the B copy, R1 is applied to the exons outside A's RefSeq record. On
these families this changes one record: member `PKD1P6-NPIPP1` (chr16:15,105,353-15,141,806, −) maps to CAT NPIPP1
`CHM13_G0020725` (15,105,363-15,124,458, −), not to R1's pick `LOFF_G0001012` (Liftoff PKD1P6). It is applied in the member
image (row 7), in every truth re-keyed here (rows 18, 27: NPIPA = Dishuck A4), in the S2 seeds if the record were a seed (it is
not: seeds were selected by name) and in the name map used for hardcoded names. The step-1 TSVs stay as R1 produced them;
the override is applied downstream, as the amendment says.

## Node naming on CAT (needed by every row)

CAT gene names repeat across paralogs (1,198 names), and the code keys several tables by name and uses
`gene_id = "gene-" + name`. Each CAT gene gets a unique **label**: its `gene_name` when that name is unique among the 64,213
CAT genes, else `gene_name~<CAT gene_id>` (e.g. `NPIPB15~LOFF_G0001213`). `gene_id` = `gene-<label>`. The CAT gene id, gene
name, source and Ensembl `source_gene` are kept in extra columns. Labels are display names only: every RefSeq-to-CAT
correspondence goes through R1 (`chm13v2.0_CAT_Liftoff.refseq_map.tsv`; exon-less RefSeq records by span Jaccard,
Amendment 1), never by name. The ~70 RefSeq names hardcoded in `npip_tbc1d3.py` (family anchors NPIPB2/TBC1D3, chaining
"genes of interest", report examples, the clause-5 TBC1D3 groups) are translated in CAT mode to the label of their R1 image
(`name_map.tsv`); a name without an image is reported as absent, never substituted by a same-named CAT gene.

## Inputs, one row each

| # | input (RefSeq, as read by the code) | CAT replacement | command / how | ruling |
|---|---|---|---|---|
| 1 | `light/work/refseq/genes.tsv` (every gene/pseudogene record) | every CAT gene of `chm13v2.0_CAT_Liftoff.genes.tsv` (64,213), columns kept (`gene_id, name, biotype, description, chrom, start0, end, strand`) + `cat_gene_id, cat_gene_name, cat_source, source_gene`. `description` = `CAT <id> <gene_name> (<source>)`, plus `; readthrough (R3)` for R3 read-throughs (the code's `"readthrough" in description` test) | `python3 CAT/light/scripts/refseq_tables_cat.py` (calls `cat_layer_order.py tables`) | R3 |
| 2 | `light/work/refseq/exons.tsv` (exon union per gene; no exon -> span) | exon union of the gene's transcripts' exons (the slim GFF3; every CAT gene has exons) | same | — |
| 3 | `light/work/refseq/cds.tsv` (one protein per gene: the transcript with the longest total CDS, ties: first transcript id in sorted order, `annotation_nodes.longest`) | the same rule on the CDS lines of the full `chm13v2.0_gencode.gff3` (streamed, never loaded whole), grouped by transcript `Parent` -> gene | same | same rule as `refseq_tables.py` |
| 4 | `light/work/refseq/gene_dbxref.tsv` (HGNC id per gene) | CAT `source_gene` (Ensembl id, version stripped) -> HGNC `ensembl_gene_id` -> `HGNC:<id>`; genes without a match listed. In CAT mode `hgnc_lookup` does NOT fall back to the symbol (R4 has no symbol route) | same | R4 |
| 5 | `light/work/cat/cat_genes_exons.tsv` (CAT exon unions for the RefSeq -> Soto map) | kept (it is already CAT v2.0); unused by the CAT Soto join | — | R5 |
| 6 | Soto mapping (`soto_map_gene`: RefSeq exon overlap -> Soto gene) | identity: a CAT node's Soto gene is its own CAT id when that id is in Soto's 2,334-gene universe (`soto_gene_to_families.tsv`); quality `strong`; ambiguity flag from Soto | `LO_ANNOT=cat` branch in `lattice_common.soto_map_gene` | R5 |
| 7 | `light/members.tsv` (44) and `members.corrected.tsv` (46, description rule) | R1 image of the 46 corrected members (family carried over; Amendment 2 for PKD1P6-NPIPP1). Drops (no CAT gene, e.g. RefSeq NPIPB3) and merges (two RefSeq members -> one CAT gene) listed. CAT genes whose `gene_name` starts with NPIP/TBC1D3 and has no `-` are listed beside, not substituted. In CAT mode `corrected-tables` takes the member set from `members.tsv` instead of the RefSeq description rule | `cat_layer_order.py members`; `LO_ANNOT=cat` branch | R1, R2 |
| 8 | Read-through records overlapping a member on the same strand (expr-recount `--ignore`, 6 RefSeq names) | R3 read-throughs (name `A-B` with A and B CAT gene names, or the R1 image of a RefSeq read-through) overlapping a CAT member on the same strand | `cat_layer_order.py members` -> `CAT/light/work/refseq/readthrough_over_members.txt` | R3 |
| 9 | `light/work/P/proteins.faa`, `proteins.index.tsv`, `protdb.*` (20,088 proteins) | CAT proteome by the same rule (`layer_protein.py prep`: r2 biotype filter `truth.excluded(bt, 2)`, translation of the longest CDS, >= 10 aa) | `python3 CAT/light/scripts/layer_protein_cat.py prep` | — |
| 10 | `light/work/P/blastp.tsv`, `searched.txt`, `S.round*.txt`, `todo.round*.txt` (closure rounds 1-8 complete, round 9 stopped after 600 of its queue; 4,430 searched) | the same closure on the CAT proteome: `plan N` / `advance N` for N = 1..8 to completion, then round 9's queue searched for its first 600 proteins (queue order = pid order, as in RefSeq) and stopped. Same blastp command (`-evalue 1e-5 -max_target_seqs 100000`, outfmt 9 columns). **The HSPs come from the step-6 all-vs-all** (row 24): `chunk N K` takes a query's rows from it instead of re-running blastp. Valid because a query's HSPs do not depend on its batch, thread count or on `-dbsize` set to the true residue count (RefSeq check, family_cert A1.2); re-checked on CAT on >= 25 queries run with the closure's own flags. If that check fails, the closure runs real blastp per round instead | `layer_protein_cat.py plan/chunk/advance` | same rounds, E-value, stopping point |
| 11 | `light/P.groups.tsv`, `P.edges.tsv`, `P.stability.tsv`, `P.members_status.tsv` | `layer_protein_bounded_cat.py hops 9`, `probe 5 4`, `final 2 3 4` (same k set; the stopping rule "first k in {2,3} whose member clusters equal those at k+1", unchanged) | as listed | same k rule |
| 12 | `light/work/D/c15_17_22.e1.*`, `c16_19_20.e1.*` (E1 catalogs: graph, loci, clusters) | CAT gene-body catalogs on the same trios. **Reuse:** `o1_falsemerge/lit/annot_gencode/` (chr15/17/22) and `lit/aj_ho/gencode/` (chr16/19/20) are the CAT/GENCODE gene-span catalogs of prereg AI/AJ, built with `bench/node_graph_mcl.py` (`minimap2 -x asm20 -c -X -N 50 -p 0.1`, index `-x asm20`). Their node tables equal, as multisets of (chrom, start, end, strand, exon union), every CAT v2.0 gene on those chromosomes (7,388 and 7,300; checked before any result), and their `all.paf` = the concatenation of all chunk PAFs (checked). Before reuse: spans re-extracted from the genome and compared, and one chunk re-mapped with the recorded command and diffed. `mcl_families` is re-run on them with the recorded E1 flags `--min-exonic-bp 1 --min-shared-exon-frac 0.0 --dump-graph`, plus `--min-cov-shorter 0` (the containment escape became default 0.70 on 2026-09-29; 0 = the behaviour every RefSeq catalog was built with). Check first: the current binary with these flags reproduces the RefSeq E1 `clusters.tsv` md5s aef97aa0 / 8e753e3c | `mcl_families --paf <reused all.paf> --gff <its nodes.gff> ...` | same flags |
| 13 | `light/work/S1/*.e1s.*` (validation dumps only) | same, `--min-shared-exon-frac 0.30` | as 12 | — |
| 14 | E0 catalogs (`human2/guided`, `lit/aj_ho/refseq/e0`; truth_guided E0 column) | same PAFs, `--min-exonic-bp 0 --min-shared-exon-frac 0.0 --min-cov-shorter 0` (checked to reproduce the RefSeq E0 md5s first) | as 12 | — |
| 15 | catalog record keys -> genes (`genes.regions` / `nodes.tsv` + `.names.tsv`) | both CAT catalogs read as `nodes` kind; `nodes.tsv.names.tsv` rewritten with CAT labels (idx -> label of the CAT gene with those coordinates and exon union) | `cat_layer_order.py catalogs` | — |
| 16 | `light/D.groups.tsv`, `D.edges.tsv` | `layer_dna_cat.py D D c15_17_22=...:nodes:... c16_19_20=...:nodes:...` on the CAT catalogs | as listed | — |
| 17 | `light/C.groups.tsv`, `C.supported_clades.tsv` (§6js IQ-TREE trees of the RefSeq literature records) | **RefSeq-only, relabelled.** The trees are a separate pipeline (`o1_falsemerge/lit/guided_t`, reference-projected alignments + IQ-TREE) not rebuilt here. Leaves are relabelled to the R1 image of each RefSeq leaf; a leaf without a CAT image (NPIPB3) is dropped from every split. C_tree / literature-C rows are therefore RefSeq trees on CAT labels, flagged as such in the report | `layer_c_cat.py` | R1 |
| 18 | `docs/lit_subclusters_npip_tbc1d3_truth.tsv` (31 literature records, by RefSeq name) | the same 31 rows re-keyed through the step-1 tables (`docs/lit_subclusters_npip_dishuck_check.CAT.tsv`, `docs/tbc1d3_members.CAT.tsv`): name = CAT label of the R1 image, level1/level2 carried over; NPIPB3 dropped. `LIT` points to `CAT/light/work/refseq/lit_truth.cat.tsv` in CAT mode | `cat_layer_order.py lit` | step 1 |
| 19 | `light/truth_soto.tsv`, `truth_soto_families.tsv`, `truth_hgnc.tsv`, `truth_guided.tsv`, `truth_literature_subfamilies.tsv`, `universe_light.tsv` | `truths_universe_cat.py`: Soto by Gene ID (R5; `best_refseq_gene_id` column = the CAT node of each Soto gene id), HGNC by R4, literature from row 18, guided from the CAT E0/E1 (rows 12, 14) | as listed | R4, R5 |
| 20 | `heavy/S2.*` (SD98 closure; seeds = 39 RefSeq genes named NPIP*/TBC1D3* with exons) | same closure (`s2_round_cat.py plan/collect`, rounds until a round adds nothing, cap 5 + check round as recorded) with CAT genes/exons (distinct exon intervals; every CAT gene has exons) and seeds = R1 images of the 39 RefSeq seeds (R2). Same SD98 input and minimap2 command (`-c --end-bonus 5 --eqx -N 50 -p 0.5 -t 4`, `chm13v2.0.fa.mmi`). A side already mapped in the RefSeq run reuses its PAF records (a side's records do not depend on the other queries of its chunk); new sides are mapped. Then `s2_finalize_cat.py` | `heavy/scripts/*_cat` | R2 |
| 21 | `heavy/EXPR.counts.tsv` | `expr_counts_cat.py` + `expr_write_cat.py`: same count rule; "unique" = exactly one CAT gene genome-wide | as listed | — |
| 22 | expr-recount `EXTRA` ids (RefSeq U genes EXPR lacks) and the expression GFF (`chm13v2.0_RefSeq_full.gff.gz` in `gff_exon_index`) | EXTRA = CAT U genes absent from EXPR (computed from members ∪ P ∪ D ∪ C before the run); GFF = a CAT GFF3 derived from the slim file with gene `ID=gene-<label>` | `cat_layer_order.py exprgff`; env | — |
| 23 | `integrate_slim/P_N2_clusters.tsv` (+ `P_stability_plain.tsv`, `P_aa50.tsv`, `P_variant_labels.tsv`), frozen | rebuilt with the archived `lo_p_variants.py` (`git show notebook-2026-09-20:archive/bench/layer_order/lo_p_variants.py`), paths pointed at `CAT`, on `layer_protein_bounded_cat.py` | `python3 CAT/light/scripts/lo_p_variants_cat.py` | — |
| 24 | (step 6) `family_cert/protein/hsps.tsv`, `pairs.tsv` (all 20,088 x database blastp, `-dbsize`) | CAT proteome (row 9) all-vs-all with the D2/R3 flags: `-evalue 1e-5 -max_target_seqs 100000 -dbsize <CAT total residues>`, outfmt 15 columns; **4 threads** instead of 5 (threads do not change HSPs, R3); batches of 400, each run under `rlock.sh heavy`, < 9 min per call; TTN-like very long queries searched standalone (A1.3) | `family_cert_cat/protein/run_batches_cat.sh`, `build_tables_cat.py` | D2 |
| 25 | (step 6) `family_cert/dna/nodes.tsv` (RefSeq records with >= 1 exon, representative transcript, body) | every CAT gene (all have exons): exon union of transcripts clipped to span; representative transcript = among `MANE_Select`-tagged transcripts if any (CAT has no "RefSeq Select" tag; 15,532 MANE_Select transcripts), else all, most exonic bp, ties smallest id. No exon-less records exist, so the `exonless_span` variant is empty on CAT | `build_nodes_cat.py` | D3 |
| 26 | (step 6) DNA witnesses (`dna/batches`, two-hop expansion, `witnesses.tsv`) | same queries (tx: `-x splice -uf -c -N 50 -p 0.1`; body: `-x asm20 -c -N 50 -p 0.1`; prebuilt `npip_ladder/idx` indexes), same two-hop expansion from the CAT member nodes, same witness rules; query reuse by sequence md5 from `npip_ladder/union` and from `family_cert/dna/batches` (same command) | `dna_cert_cat.py reuse/queries/touch/witnesses` | D3 |
| 27 | (step 6) certificate sets | NPIP / TBC1D3 = CAT members (row 7); NPIPA/NPIPB and Iso-Seq groups from the step-1 Dishuck table (`dishuck_group` field, CAT id) | `certify_cat.py` | D4, step 1 |
| 28 | `lattice/pre_correction_1712/nodes.tsv` (only for a log line) | removed from the CAT copy (comparing CAT V with the RefSeq 17:03 V has no meaning) | `rm` | — |
| 29 | `light/work/xchrom/*`, `heavy/diag_p0`, `heavy/scripts/check_bedtools.sh`, `diag_p0.py` (diagnostics, read by no stage) | not rebuilt; reported RefSeq-only | — | — |
| 30 | `human_testis.t2t.bam`, `hgnc_complete_set.txt`, `soto_gene_to_families.tsv`, `bench/soto/soto_famCN_S1C.tsv`, `final_human_clean.bed` | unchanged (annotation-free) | — | — |

## Unchanged by construction

Every threshold, identity cut, coverage cut, MCL inflation, k set, rounds count and stopping rule of the RefSeq run (R6).
The lattice tests themselves (`lattice_common.tests`) are not touched.

## What is reported RefSeq-only (not rebuilt)

- The clause-5 / literature C layer trees (row 17): leaves relabelled, trees not re-inferred.
- Diagnostics of row 29.
- Anything in this table that fails its pre-use check (rows 10, 12, 14, 20): the report says so and falls back as stated.

## Order of execution

1. Byte-identity baseline: copy ROOT -> `/mnt/linuxdisk/tmp/step56_tmp/refA`, run the CURRENT CLI `all --with-check-c2` there.
2. CAT tables (rows 1-8, 15, 18, 22), heavy tables (row 20 input), copy ROOT -> `CAT`.
3. Step-6 protein all-vs-all (row 24) -> P closure (rows 10-11) -> `P_N2_clusters` (row 23).
4. Catalog checks and `mcl_families` (rows 12-16), C relabel (17), truths (19), S2 + EXPR (20-21).
5. Code edit behind `LO_ANNOT=cat`; RefSeq byte-identity proof on `/mnt/linuxdisk/tmp/step56_tmp/refB`.
6. `LO_ANNOT=cat npip_tbc1d3.py --root CAT all --with-check-c2`.
7. Step 6 DNA (rows 25-27) and certificates.
8. `docs/LAYER_ORDER_CAT_2026-10-01.md`.
