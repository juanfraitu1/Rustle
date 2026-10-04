# Reproduce

What to install, what to run, and what number should come back. Substrate provenance is in
[`docs/DATA.md`](docs/DATA.md); what each file in the repo is for is in
[`docs/ACTIVE_WORKING_SET.md`](docs/ACTIVE_WORKING_SET.md).

⚠ **This is not an assembler project.** The assembly-only mode below is the substrate the thesis
objectives stand on, not the contribution — see `README.md` and `docs/METHOD_PSEUDOCODE.md`.

⚠⚠ **`RUSTLE_JUNCTION_MAJORITY` default flipped to ON 2026-09-21** (register row 960;
`docs/IDEAL_CHROMOSOME_SIM_2026-09-21.md` §10–11; `bench/CHR16_JUNCTION_MAJORITY_ARM.md`). A non-canonical
splice junction no longer discards an otherwise-canonical, well-supported transcript outright. Every
number quoted below that does NOT explicitly set `RUSTLE_JUNCTION_MAJORITY=1` in its own command (i.e.
most of §4's "Expected numbers") was measured under the OLD strict default and has not been re-run under
the new one — set `RUSTLE_JUNCTION_MAJORITY=0` to reproduce those numbers exactly.

## 1. Build

Rust 1.93.1, edition 2021, no submodules.

```sh
git clone <this repo> && cd Rustle
cargo build --release            # binaries land in target/release/
cargo test  --release --lib --bins
```
Expected: **883 passed / 0 failed** (lib) and **28 / 0** (bins); 19 ignored.

On a WSL2 machine with a small VHDX, build onto the big disk instead:
```sh
export CARGO_TARGET_DIR=/mnt/<bigdisk>/rustle_target TMPDIR=/mnt/<bigdisk>/tmp
```

## 2. Third-party tools

None are vendored. Only what you actually want to compare against:

```sh
conda install -c bioconda stringtie gffcompare samtools minimap2
conda create -n flair -c bioconda flair          # optional, FLAIR 3.0.0
conda install -c bioconda isoseq                 # optional
```
⚠ `gffread` is not required — the `gff_to_gtf` binary builds the per-chromosome reference GTFs.

## 3. Data

Nothing ships with the repo. Follow `docs/DATA.md` to obtain a genome, its RefSeq annotation and one
aligned Iso-Seq BAM, then slice one chromosome:

```sh
samtools view -@3 -b <library>.bam chr20 -o chr20.bam && samtools index chr20.bam
samtools faidx <genome>.fa chr20 > chr20.fa && samtools faidx chr20.fa
target/release/gff_to_gtf chm13v2.0_RefSeq_full.gff.gz chr20 chr20_ref.gtf
```

⚠ Two human libraries appear in the ledger and **their numbers are not comparable** (register row 867):
the six-chromosome panel uses `human_testis.t2t.bam`, the lab's tool comparison uses the ~6× deeper
`A119b.t2t.bam`. Our chr20 output is 658 transcripts on one and 5,844 on the other.

## 4. The assembly-only mode, and the polish

`--assemble-only` runs reads → skeletons → gate → loci → GTF and skips family detection and assignment
entirely. `--assembly-polish` adds the §6p8–§6q6 filters, which use **only** the emitted `reads "N"`
attribute — no reference, no annotation, so they are legal de novo.

```sh
# raw assembly (default: --assembly-polish none, byte-identical to the pre-polish emit)
target/release/copy_assign --assemble-only \
  --bam chr20.bam --fasta chr20.fa --region chr20:1-66210255 --out raw

# the shipped polish (2026-09-23 defaults: strict canonical junctions for the transcript product and the
# retained-intron filter at ratio 10 are ON by default under --assemble-only; both flags are shown explicitly)
target/release/copy_assign --assemble-only --assembly-junctions strict \
  --assembly-polish full --polish-isoform-fraction 0.02 \
  --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10 \
  --bam chr20.bam --fasta chr20.fa --region chr20:1-66210255 --out polished

# ⭐ 2026-09-23 (§6zb): genome-wide in ONE process — the assemble-only path streams the BAM (no reads held,
# O(distinct chains) memory per contig) and polishes per contig, so no batch script is needed:
#   whole human A119b genome ≈ 2 GB peak; gorilla ≈ 1 GB. `--materialize-reads` restores the old path.
target/release/copy_assign --assemble-only --genome-wide \
  --assembly-polish full --polish-isoform-fraction 0.02 \
  --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10 \
  --bam A119b.t2t.bam --fasta chm13v2.0.fa --out genome

# the 2026-09-22 output, byte-for-byte (majority-tolerated junctions, no retained-intron filter)
target/release/copy_assign --assemble-only --assembly-junctions majority --polish-retained-ratio 0 \
  --assembly-polish full --polish-isoform-fraction 0.02 \
  --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 \
  --bam chr20.bam --fasta chr20.fa --region chr20:1-66210255 --out polished_0922

# high-recall levers (use on deep libraries; see the caveat below). Without the polish flags this is
# the raw high-recall arm; add them back for the polished one -- the two give very different counts.
RUSTLE_JUNCTION_MAJORITY=1 target/release/copy_assign --assemble-only --read-isoform-k 3 \
  --bam chr20.bam --fasta chr20.fa --region chr20:1-66210255 --out recall_raw

RUSTLE_JUNCTION_MAJORITY=1 target/release/copy_assign --assemble-only --read-isoform-k 3 \
  --assembly-polish full --polish-isoform-fraction 0.02 \
  --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 \
  --bam chr20.bam --fasta chr20.fa --region chr20:1-66210255 --out recall_polished
```

### Abundance and locus outputs (§6r5)

```sh
# add count-based cov + TPM to the GTF (default off; unset the GTF is byte-identical)
target/release/copy_assign --assemble-only --gtf-tpm ... --out run

# loci as BED, plus a one-to-one match against the annotation showing size agreement
target/release/locus_bed run.gtf --out run --ref chr20_ref.gtf
#   run.loci.bed / run.ref_loci.bed / run.locus_match.tsv
```
`TPM_i = reads_i / Σ reads * 1e6` — **not** length-normalised, because a long read is one molecule.
Measured against StringTie's TPM on chr20: count-based ρ 0.879, length-normalised 0.714. On chr20 our
median locus size ratio against the annotation is **0.994**, with 55.6% inside ±10% (StringTie 53.7%,
FLAIR 47.2%). Details: `bench/TPM_AND_LOCUS_BED.md`. Both tools are Rust binaries (§6r9), so nothing in this file needs Python.

Score any of them:
```sh
gffcompare -r chr20_ref.gtf -o cmp polished.gtf && cat cmp.stats
```

### Expected numbers

**`human_testis.t2t.bam`, chr20** (188,864 records) — `polished` should give **658 mRNAs, 336 matching
intron chains, intron-chain 7.8/51.6, transcript 7.4/51.2**.

**`A119b.t2t.bam`, chr20** (1,104,846 records) — `polished` (2026-09-23 defaults) gives **5,522 mRNAs,
1,059 chains, 24.7/21.0, 23.2/19.2**; `polished_0922` gives 5,844 mRNAs, 1,064 chains, 24.8/19.9, 23.3/18.2
(`docs/PREREG_assembly_precision_levers_2026-09-23.md`: held-out gorilla precision 33.1 → 35.6 for −0.46% chains). `recall_raw` gives **20,699 mRNAs and 1,259 chains at 29.4/8.2** — past isoseq's
1,253 — while `recall_polished` gives 7,437 mRNAs and 1,101 chains at 25.7/16.4. Against the lab's arms
on the same BAM: StringTie 861 chains at 20.1/16.8, FLAIR 1,026 at 23.9/7.6, isoseq collapse 1,253 at
29.2/**3.0** from 64,384 transcripts.

**`GGO_mm.bam`, `NC_073244.2`** (473,231 records) — `polished` gives **4,063 mRNAs, 1,575 chains,
28.4/38.9, 26.6/38.8**, beating StringTie (1,374 at 24.7/37.0) and FLAIR (1,393 at 25.1/23.6) on every
metric. `recall_polished` gives **5,298 mRNAs, 1,688 chains at 30.4/31.9, 28.5/31.9 — which beats isoseq
collapse outright on every metric** (20,643 mRNAs, 1,655 chains at 29.8/9.0, 28.1/8.1).

Full tables: `bench/ASSEMBLY_POLISH.md`, `bench/LAB_DATASET_BAKEOFF.md`, `bench/SQANTI3_POLISH.md`.

⚠ **The polish is depth-sensitive.** `--polish-isoform-fraction` costs ~2 matching chains on the shallow
library and **65 on A119b chr20**; the mono floor and shadow rule are free on both. It was tuned on the
shallow library, and making it depth-aware is the top open item.

## 5. Family definition (the actual objectives)

```sh
# de novo family stage in ONE command (2026-09-23): loci from the assembled GTF, all-vs-all, MCL
target/release/mcl_families --from-gtf run.gtf --fasta GENOME.fa --min-exonic-bp 1 --min-shared-exon-frac 0.60 --out run.fam
#   (writes run.fam.loci.gff3 / .loci.fa / .loci.paf, then run.fam.clusters.tsv; the older two-step form:)
target/release/mcl_families --paf all_vs_all.paf --gff loci.gff3 --min-exonic-bp 1 --min-shared-exon-frac 0.60 ...
```
Score the resulting clusters against a family truth (sensitivity / precision / one-to-one bipartite F /
collapse — the standing reporting rule) with the native scorer, byte-identical to the retired
`bench/mode_family_score.py` including scipy's assignment tie-breaking:

```sh
target/release/family_score --clusters run.clusters.tsv --gff chr20_ref.gff --soto bench/soto/soto_famCN_S1C.tsv --chrom chr20 [--family NPIP]
```
⚠ `--min-shared-exon-frac` is **inert without `--min-exonic-bp 1`**. The RNA-level definition is
components of L3 at `w_98 >= 0.985`, with L4 (0.995) applied selectively; see
`docs/seeded_family_definition.md` §0★★ and ledger §6p0/§6p1/§6p5.

⚠⚠ **`--min-cov-shorter` default flipped 0 → 0.70 on 2026-09-29** (the user's decision; shipped opt-in in f2144faf,
register 1006/1014, `docs/PREREG_cov_shorter_adoption_2026-09-22.md`): a pair whose coverage of the LONGER locus fails
also passes at coverage ≥ 0.70 of the SHORTER locus's exonic length, with that coverage as its edge weight. Every
`mcl_families` number quoted in this file that does not name the flag was measured at 0 — **`--min-cov-shorter 0`
(driver `RUSTLE_MIN_COV_SHORTER=0`) reproduces it byte for byte** (proven on human_testis against a 3007c3d4 build).
Known regressions of the new default: NPIP in GUIDED mode, Soto F .833 → .800; semi-guided SD-region nodes, precision
.973 → .833 (register 1007/1009) — never use it with region nodes.

### 5a. Soto 2025 family replication (concordance, not independent: register T15 / 858; ledger §6ie–§6ip)

The chain uses Soto's own famCN (S1C), CAT v4 genes and gene universe, so it measures **concordance** with Soto, and
the famCN leg is circular (register 858). In repo: `bench/soto/soto_famCN_S1C.tsv`, `soto_parCN_S1E.tsv`,
`acro_extra_anchors.tsv`, and `shared_exons_2334_finalhuman.tsv` (the frozen output of `edges`). Not in repo:
`final_human_clean.bed` (CHM13 v2.0 SEDEF, 88,756 rows, header stripped) and `cat_v4.bed` (CAT v4, CHM13 v1.0), under
`winloci_data/soto_replication/`.

```sh
S=bench/soto/soto_replication.py
python3 $S genesets --out-eligible g1793.tsv --out-full g2334.tsv
python3 $S edges --sedef final_human_clean.bed --geneset g2334.tsv --cat-bed cat_v4.bed \
    --extra-anchors bench/soto/acro_extra_anchors.tsv --out-shared shared.tsv        # 4,192 edges
#   (or skip it: bench/soto/shared_exons_2334_finalhuman.tsv is this step's frozen output)
python3 $S cluster --shared shared.tsv --geneset g1793.tsv --full-geneset g2334.tsv \
    --famcn bench/soto/soto_famCN_S1C.tsv --mad-statistic median --out rep_median.tsv
python3 $S score --predicted rep_median.tsv          # --truth defaults to bench/soto/soto_famCN_S1C.tsv
```
Expected (median / mean MAD): ARI **0.6959 / 0.6862**, exact 241/491 / 264/491, pair P/R/F1 0.841/0.595/0.697 /
0.906/0.554/0.687; bipartite MICRO 0.784/0.709 / 0.812/0.721, MACRO 0.731/0.718 / 0.770/0.739, undetected 99/491 /
88/491. `rep_{median,mean}.tsv` are byte-identical to the frozen `replicated_families_2334_{median,mean}_finalhuman.tsv`
(re-verified 2026-09-29; `genesets` + `cluster` + `score` for both statistics take about 4 s). `famcn` (WSSD famCN
at arbitrary intervals) needs pyBigWig (the miniforge python) or `bigBedToBed`; `score` uses scikit-learn's ARI when
importable (a stdlib one otherwise), numpy and scipy.

**Reconciled recipe (2026-09-29; `docs/PREREG_soto_reconciliation_2026-09-29.md`, register 1162-1166;
`docs/SOTO_REPLICATION_STATUS_2026-09-28.md` §1).** Two choices of Soto's *released code* — map SD98 exons back
(not regions) and gate each shared-exon pair by famCN MAD < 1 then grow families through coding genes (not the
component split) — close the gap to their Table S1C. Both are opt-in flags; the chain above is unchanged.

```sh
S=bench/soto/soto_replication.py; W=/mnt/linuxdisk/home/juanfraitu/winloci_data/soto_replication
python3 $S genesets --out-eligible g1793.tsv --out-full g2334.tsv
bash tools/rlock.sh heavy python3 $S edges --exon-mapback --cat-bed $W/cat_v4.bed --sd98-bed $W/sd98_v1.bed \
    --genome $W/t2t-chm13-v1.0.fa.gz --threads 5 --out-shared exon_edges.tsv     # 12,231 edges; ~8 GB, 2-8 min
#   (--mm2-index IDX.mmi reuses a `minimap2 -d` index; or skip the step: bench/soto/shared_exons_5154_exon_mapback.tsv
#    is its frozen output, byte-identical on the 2026-09-29 re-run)
python3 $S cluster --pair-mad --shared bench/soto/shared_exons_5154_exon_mapback.tsv --geneset g1793.tsv \
    --full-geneset g2334.tsv --famcn bench/soto/soto_famCN_S1C.tsv --out pair_s1c.tsv --out-cover pair_s1c.cover.tsv
python3 $S score --predicted pair_s1c.tsv
python3 $S score --predicted pair_s1c.tsv --split bench/soto/soto_split_2026-09-29.tsv --half heldout --only pairs
```
Expected: ARI **0.9698**, exact **479/491**, pair P/R/F1 1.000/0.942/0.970, MICRO 1.000/0.980, MACRO 0.999/0.992,
undetected 0/491, 504 predicted families (158 genes in ≥ 2, `pair_s1c.cover.tsv`); held-out **0.9681**, 263/266
(DEV 0.9708, 216/225). The `gene_id`/`family_id` projection of `pair_s1c.tsv` over the 2,334 genes is byte-identical
to the frozen reconcile partition. This is concordance with Soto's tables (famCN, universe and curated list are
theirs; register 858 / 1085).

**The ladder** — what the copy-number gate buys on the same edges (register 1169 / 1170; quote exact families and
the ARI without FAM90A beside every famCN rung, the ALL-491 ARI moves ±0.035 on that one family):
```sh
bash tools/rlock.sh heavy /home/juanfra/miniforge3/bin/python3 $S famcn --interval exons --samples all \
    --cat-bed $W/cat_v4.bed --sd98-bed $W/sd98_v1.bed --wssd-dir /mnt/linuxdisk/home/juanfraitu/winloci_data/soto_wssd \
    --matrix famcn_matrix.npz --jobs 4 --out famcn_exons_268.tsv        # 269 tracks once, ~130 s, 0.6 GB; then instant:
/home/juanfra/miniforge3/bin/python3 $S famcn --interval sd98  --samples all --matrix famcn_matrix.npz --out famcn_sd98_268.tsv
/home/juanfra/miniforge3/bin/python3 $S famcn --interval exons --samples all --outlier '' --matrix famcn_matrix.npz --out famcn_exons_269.tsv
/home/juanfra/miniforge3/bin/python3 $S famcn --interval exons --samples 10  --matrix famcn_matrix.npz --out famcn_exons_10.tsv
python3 $S ladder --shared bench/soto/shared_exons_5154_exon_mapback.tsv --geneset g1793.tsv --full-geneset g2334.tsv \
    --famcn-ours $W/famcn_ours_allwssd.tsv --famcn-ours10 $W/famcn_ours_all.tsv \
    --split bench/soto/soto_split_2026-09-29.tsv --drop-family ID_356
```
Expected: `famcn_exons_268.tsv` = columns 1-4 of `$W/famcn_ours_allwssd.tsv`, `famcn_exons_269.tsv` column 2 = its
`famCN_269`, `famcn_sd98_268.tsv` columns 2 and 5 = its `famCN_sotoiv` / `n_sotoiv_rows` (the paste of the three is
byte-identical, sha1 11daa3ce), `famcn_exons_10.tsv` = `famcn_ours_all.tsv` modulo its CRLF line ends. Ladder (ARI all
/ DEV / HELD-OUT, exact, ARI without ID_356): sequence only 0.7307 / .6418 / .8693, 345, 0.7057; our famCN 10
samples exons 0.9198 / .9096 / .9317, 373, 0.9131; 268 samples exons 0.8855 / .9039 / .8610, 375, 0.9089; **268
samples, Soto's interval 0.9277 / .9227 / .9343, 411, 0.9251**; S1C 0.9698 / .9708 / .9681, 479, 0.9650. (Rungs 2-3
are 0.9197 / 0.8853 in the prereg, whose scorer ordered families sharing a smallest member by Python-set order;
the module's order is deterministic and equals the frozen reconcile output; exact counts are identical.)

**Assembly parCN** (`bench/soto/parcn_assembly.py`; `docs/PREREG_soto_parcn_assembly_2026-09-29.md`, register
1171-1174; QuicK-mer2 itself needs ~52 GB, register 1167). From the frozen k-mer count tables (`docs/DATA.md`):
```sh
A=/mnt/linuxdisk/tmp/rustle_figures_dev/soto_parcn_asm
P=bench/soto/parcn_assembly.py
python3 $P analyze --regions $A/work/regions.tsv --q $A/work/Q.u64 --pos $A/work/pos.npz \
    --counts CHM13=$A/counts/chm13noY.i32 HG002=$A/counts/hg002.i32 GGO=$A/counts/ggo_mat.i32+$A/counts/ggo_pat.i32 \
             GGOmat=$A/counts/ggo_mat.i32 GGOpat=$A/counts/ggo_pat.i32 PTR=$A/counts/ptr.i32 PPY=$A/counts/ppy.i32 \
    --edit-depth $A/counts/ed0.u32,$A/counts/ed1.u32,$A/counts/ed2.u32 --out-prefix parcn_      # ~45 s, 1.6 GB
python3 bench/soto/test_parcn_assembly.py                                                        # 8 tests, ~6 s
```
Expected (`parcn_summary.json`, all 113 values equal to the frozen run): controls HG002 parCN = 2 for 297/299;
resolved 1,163/1,831; **Fixed 321/322 = 0.997** within 0.5 of S1E (mean statistic 0.994), Nearly-Fixed 408/629 =
0.649, Polymorphic 83/212 = 0.392, Spearman 0.414; H1 105/109 = 0.963, H2 105/132 = 0.795, H3 131/631 = 0.208,
H4 105/118 = 0.890; F-cal 0.617 / 0.886. Regenerating the counts (`regions` → `kmers` → `count --meryl --exclude chrY`
per genome → `edit-depth`) is heavy (a whole-genome meryl DB per assembly) and was not re-run; `count` and
`edit-depth` are verified against a brute force by the unit tests.

## 6. Before proposing anything

`docs/NEGATIVE_RESULTS_REGISTER.md` (870 rows) records what has already been refuted, with the reason.
`docs/o1_ledger.md` is the running measurement record. Several claims in older sections are explicitly
retracted in later ones — the ledger is append-only, so **the latest section wins**.

## Known reproducibility caveats

- **Tool builds differ across the six-chromosome panel**: chr20's StringTie arm is 3.0.1, chr11/7/14/5/9
  are 3.0.3 (register row 870).
- **A second, broken FLAIR install can shadow a working one**; `flair collapse` then dies with
  `ModuleNotFoundError: No module named 'flair'`. Fix with per-script shims (§6q6).
- **`isoseq collapse` needs PacBio-style read names AND a relaxed `--min-aln-coverage`** together; with
  SRA-style names it silently skips every read (register row 865, itself a retraction).
- Two datasets (`A119b`, `GGO_OR6737`) still need public accessions filled into `docs/DATA.md`.

## Missing copies from RNA alone: flag, characterise, screen, hand to DNA (thesis O3; §6ze, 2026-09-23)

`docs/PREREG_o3_rna_only_2026-09-23.md`. One BAM (primaries with `de:f` and `--eqx`), a locus set, the primary
genome and its minimap2 splice index; optional annotation GFF (IG/TR screen), `--confirm` genomes (a
haplotype assembly: DNA confirmation) and `--foreign` genomes (another species: contamination screen).

```bash
# scan in contig batches (a laptop-sized foreground job each), then align once per genome
missing_copy_flag --bam READS.bam --fasta GENOME.fa --loci GENES.gff --gff GENES.gff --index x \
            --contigs chr1,chr2,... --out run_b1 --scan-only
missing_copy_flag --bam READS.bam --fasta GENOME.fa --loci GENES.gff --index GENOME.splice.mmi \
            --confirm pat=PAT.splice.mmi --confirm mat=MAT.splice.mmi --foreign human=CHM13.splice.mmi \
            --out run --from-scan run_b1,run_b2
```

Output `run.missing_copy.tsv` (one row per expressed locus; `class` ∈ divergent / structural / both — the structural
detector (exon-order rearrangements carried as ≥ 50 bp insertions, addendum 2) is on by default; verdict ∈
contamination / foreign_species / hypermutation / rna_editing / scattered / unannotated_paralogue /
reference_absent_candidate, plus
`expected_dna_depth_ratio` and per-genome confirmation) and `run.consensus.fa` (the hidden copy's spliced
consensus — the probe for a DNA k-mer/depth check). Positive control:
`python3 bench/sim.py missing-copy genomic REF.gtf REF.fa chr20 40 0.02 simB SEED` (40/40 at 2% divergence, all
confirmed — measured with the out-of-repo `/mnt/linuxdisk/tmp/gw22/o3/simB.py`, whose mutations were seeded by Python's
per-process `hash()`; `sim.py` seeds with `stable_seed()` since wave 7, so re-measure before quoting).

---

> **Moved documents (wave 4, 2026-09-23).** Files this document cites that were pruned from the working tree — `docs/o1_investigations.md`, `docs/OBJECTIVES_AND_VERIFICATION.md`, `docs/o3_missing_copy_evidence.md`, `docs/NUMBERS.md`, `docs/ONE_METHOD.md`, `docs/METHOD_PSEUDOCODE.md`, `docs/OPEN_ITEMS_2026-09-09.md`, `docs/superpowers/` — are at git tag `notebook-2026-09-23` (`git checkout notebook-2026-09-23 -- <path>`) and in `~/Desktop/Rustle_attic/2026-09-23/` (see its `MANIFEST.tsv`). The citations above are provenance and were left as written.

> **Renamed bench scripts (wave 7, 2026-09-24).** The 21 top-level `bench/*.py` scripts became `bench/lib.py`
> (shared helpers) and three subcommand scripts: `bench/score.py` (scorers), `bench/sim.py` (simulators) and
> `bench/truth.py` (truth builders); `bench/guided_pipeline.py` and `bench/mcl_port.py` stay. Each new module's
> docstring maps old command -> new command, and the old files are at git tag `notebook-2026-09-24`
> (`git show notebook-2026-09-24:bench/<old>.py`). Every replaced scorer and the two hash-free simulators were checked
> byte-identical on recorded inputs; the O2/O3 simulators (`sim.py copies`, `sim.py missing-copy`) now use stable
> seeds, so their reads differ from every earlier run.
> The 8 `bench/soto/*.py` Soto-replication scripts became `bench/soto/soto_replication.py` (subcommands `genesets`,
> `edges`, `cluster`, `dennislab`, `famcn`, `score`; §5a), checked byte-identical on the recorded inputs, and
> `bench/soto/rustlib.py` (0 importers) left the tree; both are at the same tag.
> The 12 `bench/layer_order/*.py` files became `lattice_common.py` (library) and `npip_tbc1d3.py` (one subcommand per
> stage; recipes in `bench/LAYER_ORDER_NPIP_TBC1D3.md` and `bench/NESTED_LATTICE_NPIP_TBC1D3.md` §11).
> `bench/README.md` has the full old-name → new-command table, including the names these scripts had before 2026-09-24.

## Copy-assignment accuracy with read-level truth (thesis O2; §6zf, 2026-09-23)

```bash
python3 bench/sim.py copies CAT.copies.tsv CAT.copies.fa GENOME.splice.mmi run 20260923   # reads from every copy, mapped genome-wide
copy_assign --bam run.bam --fasta GENOME.fa --regions whole_chromosomes.txt --families CAT.copies.tsv --copies-fa CAT.copies.fa --out run_o2
python3 bench/score.py reads --catalog CAT.copies.tsv run run_o2      # OWN / PRIMARY / ANY readings, per divergence bin
```
Expected (human chr16 catalog `chr16_arm/on`, copies < 300 bp dropped): OWN 157/157 correct, 0 wrong, 1,088 abstain of 1,259 MAPQ-0 reads.
⚠ That run's per-copy read seeds came from Python's per-process `hash()` (wave-7 defect B2), so it cannot be regenerated
read for read; `score.py reads` reproduces the numbers exactly on the recorded run (`/mnt/linuxdisk/tmp/gw22/o2sim/h16`,
`h16_o2`). `sim.py copies` now seeds with `stable_seed()`, so a fresh simulation draws different reads.
**Re-measured 2026-09-25 with the stable seed 20260925** (`python3 figures/make.py data fig4`, which runs exactly the
three commands above, mapping in 8 read-disjoint parts — `sim.py copies --parts 8`, identical records): human 28,453
reads, 1,263 MAPQ-0; OWN 163 correct, 0 wrong, 0 conflict, 1,084 abstain, 16 not scored; ANY 60 correct / 432 wrong /
156 conflict; `--union-certificate` 0 assigned (1,255 of the 1,263 tied reads have an NM-identical genomic twin; row 1103 measured all 990 scored molecules of the hash-seeded run). Gorilla
(`hom_c234`, 11,448 reads): 30 MAPQ-0, 0 assigned. Tables: `figures/data/fig4_assignability_upset.tsv`,
`fig5_assign_accuracy_bands.tsv`.

## The whole pipeline in one driver (2026-09-23)

```bash
tools/rustle_pipeline.sh all --bam READS.bam --fasta GENOME.fa --out run --index GENOME.splice.mmi --gff ANNOT.gff \
    [--confirm pat=PAT.splice.mmi --confirm mat=MAT.splice.mmi] [--foreign human=CHM13.splice.mmi] [--threads 4]
# stages, each also runnable alone: assemble -> families (mcl_families --from-gtf; its copy table run.fam.copies.* is what
#   assign reads since 2026-10-02) -> assign (copy_assign --families) -> flag (missing_copy_flag scan + align). OPT-IN:
#   candidates (o3_candidates + augmentation + patch realignment, between families and assign; `--candidates` runs it in
#   `all` and makes assign and flag use it; ruling R14: its first acceptance failed, docs/O3_CANDIDATES_ACCEPTANCE_2026-10-02.md (the re-run passed: A13, docs/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md), and
#   Amendment 14's no-deletion control failed, docs/O3_CANDIDATES_CONTROL_A14_2026-10-03.md).
#   LEGACY: catalog (gw_family_catalog), which `--legacy-catalog` builds in `all` and assign then reads (refused together
#   with --candidates). Products all carry the --out prefix.
# Intermediates are cached in run.cache/ (default; --no-cache off): the catalog's collapsed representatives and every
# all-vs-all PAF, keyed by the binary, the BAM/FASTA and every upstream setting. A re-run that changes only the edge
# rule replays them (human chr16 catalog 358 s cold -> 0.9 s warm; gorilla families 225 s -> 0.6 s; byte-identical).
tools/rustle_pipeline.sh catalog --bam READS.bam --fasta GENOME.fa --out run --inspect   # + edge tables, collapse stats
tools/rustle_pipeline.sh cache-ls --out run                                               # what run.cache holds
# bridge-aware regrouping (bench/ASSEMBLY_POLISH.md 2026-09-29 addendum 3): f1v2 is THE DEFAULT since 2026-09-29 (the
# user's decision), so the line above already runs it. On the BAMs and best-AS tables of the 2026-09-25 runs, `assemble`
# writes the held-out products of docs/PREREG_f1_bridge_locus_2026-09-28.md (f1: gorilla OR6737, KB3781) and
# docs/PREREG_f1v2_readshare_2026-09-29.md (f1v2: human A119b, testis) byte for byte: run.gtf, the families input
# run.families.gtf and the side tables. Their FUSED counts are readthrough_eval `a.fused` on the latter.
RUSTLE_BRIDGE_REGROUP=f1 tools/rustle_pipeline.sh all --bam READS.bam --fasta GENOME.fa --out run ...   # the F1 arm
# THE PRE-FLIP PIPELINE, byte for byte (every driver number in this file dated before 2026-09-29 was measured this way; its
# assign ran on the legacy catalog, hence --legacy-catalog since 2026-10-02; assign on a catalog with cross-chromosome
# families differs since then, see `copy_assign --help`, --families):
RUSTLE_BRIDGE_REGROUP=off RUSTLE_MIN_COV_SHORTER=0 tools/rustle_pipeline.sh all --legacy-catalog --bam READS.bam --fasta GENOME.fa --out run ...
```

## Tandem-copy simulations: what the aligner and the pipeline do with near-identical adjacent copies (§6zg, 2026-09-24)

```bash
python3 bench/sim.py tandem --fasta chr20.fa --gtf chr20_ref.gtf --out t --layout tandem --copies 2 \
    --sweep 0.9,0.95,0.98,0.99,0.995,1.0 --pipeline --bin target/release      # also --layout interleaved, --copies 3
```
Per read: same/other copy, cross-copy chain, MAPQ, AS tie; per condition: assembler transcripts (chimeric), catalog
copies, assignment. Expected (`docs/PREREG_tandem_copy_sim_2026-09-24.md`): 0 cross-copy chains below identity 1.0.

## Identity spectrum against Ensembl Compara (§6zh, 2026-09-24)

```bash
# truth: BioMart paralogue table for one chromosome (useast mirror), then the three edge tiers on the expressed loci
python3 bench/score.py spectrum --gtf chr16_assembled.gtf --ref chr16_ref.gtf --fasta chm13v2.0.fa --chrom chr16 \
    --compara compara_chr16.tsv --out chr16 --mmseqs mmseqs      # recall by Compara identity band, precision by ours
# the same truth at the FAMILY level: score a gw_family_catalog copies.tsv (transitive families, not direct edges);
# --universe fixes the recall denominator to the tier run's expressed pairs (§6zi, rows 1097-1098)
python3 bench/score.py pairs --members chr16.copies.tsv --genes chr16_ref.gtf --chrom chr16 \
    --truth compara:compara_chr16.tsv --universe chr16.truth_pairs.tsv     # numbers unchanged (byte-identical, wave 7)
# seeding loci with GOOD secondaries (rows 1060/1100): one pass over the BAM, then the assembler reads the table
as_table --bam reads.bam --out reads.molecules.tsv --threads 4          # 66 s / 0.8 GB on an 11.7 GB gorilla BAM
RUSTLE_GTF_SECONDARY=1 RUSTLE_GTF_SECONDARY_AS_RATIO=0.98 RUSTLE_GTF_SECONDARY_AS_TABLE=reads.molecules.tsv \
    copy_assign --assemble-only ...                 # the pipeline driver does this by default (--no-seed-secondaries turns it off)
# copy assignment with ONE certificate per tied read over every placement it touches (catalog copies across
# families + outside loci built from the genome); opt-in — on the chr16 truth sim it removes every foreign claim
# (844 -> 0) and every wrong row, and ties the reads whose tied partner is an identical genomic twin (row 1103)
copy_assign --families cat.copies.tsv --copies-fa cat.copies.fa --union-certificate --bam reads.bam --fasta ref.fa --out o2
# deep human libraries whose node set is dominated by unspliced Alu stubs (chr16: 85% single-exon copies, all-pair
# precision 0.27): RUSTLE_ER_COVERAGE_LONGER_FLOOR=0.30 doubles precision (0.53) for 8/103 pairs, but costs 24% of
# referee pairs on gorilla — opt-in, never the default (§6zj, row 1099)
```

## Publication figures (2026-09-25)

`figures/` builds every figure of the paper from tidy tables with provenance headers (`figures/README.md`): gffcompare
intron chains (fig 1), SQANTI3 (fig 2), transcripts built from tied secondary alignments (fig 3), the copy-assignability
UpSet (fig 4), copy-assignment accuracy and the hard-locus benchmark (fig 5), family recall across the identity
spectrum (fig 6) and de novo vs guided family recovery (fig 7).

```bash
python3 figures/make.py data figN      # regenerate one figure's tables (foreground, cached; costs in figures/README.md)
python3 figures/make.py plot all       # render figures/out/*.{pdf,png,svg}
python3 figures/make.py check          # every table present with provenance; every figure renders
```
Every number in `figures/captions/*.md` is read from `figures/data/*.tsv`.
