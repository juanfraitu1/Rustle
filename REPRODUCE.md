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
consensus — the probe for a DNA k-mer/depth check). Positive control: `/mnt/linuxdisk/tmp/gw22/o3/simB.py`
(40/40 at 2% divergence, all confirmed).

---

> **Moved documents (wave 4, 2026-09-23).** Files this document cites that were pruned from the working tree — `docs/o1_investigations.md`, `docs/OBJECTIVES_AND_VERIFICATION.md`, `docs/o3_missing_copy_evidence.md`, `docs/NUMBERS.md`, `docs/ONE_METHOD.md`, `docs/METHOD_PSEUDOCODE.md`, `docs/OPEN_ITEMS_2026-09-09.md`, `docs/superpowers/` — are at git tag `notebook-2026-09-23` (`git checkout notebook-2026-09-23 -- <path>`) and in `~/Desktop/Rustle_attic/2026-09-23/` (see its `MANIFEST.tsv`). The citations above are provenance and were left as written.

## Copy-assignment accuracy with read-level truth (thesis O2; §6zf, 2026-09-23)

```bash
python3 bench/copy_assign_read_truth.py sim CAT.copies.tsv CAT.copies.fa GENOME.splice.mmi run 20260923   # reads from every copy, mapped genome-wide
copy_assign --bam run.bam --fasta GENOME.fa --regions whole_chromosomes.txt --families CAT.copies.tsv --copies-fa CAT.copies.fa --out run_o2
CATALOG_TSV=CAT.copies.tsv python3 bench/copy_assign_read_truth.py score run run_o2      # OWN / PRIMARY / ANY readings, per divergence bin
```
Expected (human chr16 catalog `chr16_arm/on`, copies < 300 bp dropped): OWN 157/157 correct, 0 wrong, 1,088 abstain of 1,259 MAPQ-0 reads.

## The whole pipeline in one driver (2026-09-23)

```bash
tools/rustle_pipeline.sh all --bam READS.bam --fasta GENOME.fa --out run --index GENOME.splice.mmi --gff ANNOT.gff \
    [--confirm pat=PAT.splice.mmi --confirm mat=MAT.splice.mmi] [--foreign human=CHM13.splice.mmi] [--threads 4]
# stages, each also runnable alone: assemble -> families (mcl_families --from-gtf) -> catalog (gw_family_catalog)
#   -> assign (copy_assign --families) -> flag (missing_copy_flag scan + align). Products all carry the --out prefix.
```

## Tandem-copy simulations: what the aligner and the pipeline do with near-identical adjacent copies (§6zg, 2026-09-24)

```bash
python3 bench/tandem_copy_sim.py --fasta chr20.fa --gtf chr20_ref.gtf --out t --layout tandem --copies 2 \
    --sweep 0.9,0.95,0.98,0.99,0.995,1.0 --pipeline --bin target/release      # also --layout interleaved, --copies 3
```
Per read: same/other copy, cross-copy chain, MAPQ, AS tie; per condition: assembler transcripts (chimeric), catalog
copies, assignment. Expected (`docs/PREREG_tandem_copy_sim_2026-09-24.md`): 0 cross-copy chains below identity 1.0.
