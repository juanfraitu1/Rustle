# Rustle

A long-read transcript assembler written in Rust, with a variation graph mode for gene family assembly.

## Overview

Rustle assembles transcripts from long-read RNA-seq alignments (PacBio, ONT) using a splice-graph and max-flow decomposition pipeline. It includes a variation graph (VG) mode that links gene family copies via multi-mapping reads and jointly resolves read assignments across paralogs.

## Features

- **Transcript assembly**: splice graph construction, Edmonds-Karp max-flow, seeded path extraction
- **Long-read optimized**: poly-A/T detection, junction correction, hard boundary inference
- **Variation graph mode** (`--vg`): gene family discovery, EM/MFLP/flow-based multi-mapping resolution
- **SNP-based copy assignment** (`--vg-snp`): use sequence variants to distinguish gene copies
- **Novel copy discovery** (`--vg-discover-novel`): find missing paralogs from unmapped reads via k-mer matching
- **Phased assembly scaffold** (`--vg-phase`): haplotype-aware assembly using HP tags
- **Splice consensus validation**: genome-based GT-AG/GC-AG motif checking
- **Guided and expression-only modes**: `-G` for annotation-guided, `-e` for quantification

## Quick Start

```bash
# Build
cargo build --release

# De novo assembly (long reads)
./target/release/rustle -L -o output.gtf input.bam

# With variation graph mode (gene families)
./target/release/rustle -L --vg --vg-report families.tsv -o output.gtf input.bam

# Guided assembly
./target/release/rustle -L -G reference.gtf -o output.gtf input.bam

# With genome for splice consensus validation
./target/release/rustle -L --genome-fasta genome.fa -o output.gtf input.bam
```

## VG Mode: Gene Family Assembly

The `--vg` flag enables variation graph mode for multi-copy gene families:

```bash
# EM solver (default) — junction-based compatibility scoring
./target/release/rustle -L --vg --vg-solver em -o output.gtf input.bam

# MFLP solver — linear programming via good_lp
./target/release/rustle -L --vg --vg-solver mflp -o output.gtf input.bam

# With SNP-based copy assignment
./target/release/rustle -L --vg --vg-snp --genome-fasta genome.fa -o output.gtf input.bam

# Novel copy discovery from unmapped reads
./target/release/rustle -L --vg --vg-discover-novel --genome-fasta genome.fa -o output.gtf input.bam
```

Family groups are discovered automatically from multi-mapping read patterns (supplementary alignments) and exonic sequence similarity (k-mer Jaccard). The `--vg-report` flag outputs a TSV with per-family details.

## Key Options

| Flag | Description | Default |
|------|-------------|---------|
| `-L` | Long-read mode | off |
| `-G <GTF>` | Guided assembly with reference annotation | — |
| `-e` | Expression-only (quantify known transcripts) | off |
| `-o <GTF>` | Output GTF path | required |
| `-p <N>` | Threads | auto |
| `-f <F>` | Minimum isoform fraction | 0.01 |
| `-c <F>` | Minimum coverage per bp | 1.0 |
| `--genome-fasta <FA>` | Genome for splice consensus | — |
| `--vg` | Enable variation graph mode | off |
| `--vg-solver {em,mflp,flow}` | Multi-mapping solver | em |
| `--vg-snp` | SNP-based copy assignment | off |
| `--vg-phase` | Phased assembly (HP tags) | off |
| `--vg-discover-novel` | Find novel gene copies | off |
| `--vg-report <TSV>` | Family group report | — |
| `--vg-min-shared <N>` | Min shared reads to link bundles | 3 |

## Output

Standard GTF format with additional attributes:

```
gene_id "RSTL.1"; transcript_id "RSTL.1.1"; cov "12.5"; FPKM "0.5"; TPM "150.0";
source "flow"; longcov "15.0";
```

In VG mode, transcripts from gene families include:
```
family_id "FAM_0"; copy_id "1"; family_size "3";
```

## Installation

Requires Rust toolchain (1.70+):

```bash
git clone https://github.com/juanfraitu1/Rustle.git
cd Rustle
cargo build --release
# Binary: ./target/release/rustle
```

## License

MIT License — see [LICENSE](LICENSE).

Transcript assembly pipeline inspired by StringTie (Pertea et al., Johns Hopkins University). VG mode, multi-mapping resolution, SNP assignment, and novel copy discovery are original contributions.
