#!/usr/bin/env bash
# Subset a genome FASTA to only the contigs that have mapped reads in a BAM, so a
# `rustle --vg --genome-fasta` run loads a tiny reference instead of the whole
# genome. On a region-scoped BAM this is a ~100-300x speedup (e.g. DAZ: 1.2 s vs
# ~6 min) with an identical result — the coordinates are unchanged.
#
# Usage:  subset_genome_for_vg.sh <genome.fa> <reads.bam> <out_subset.fa>
# Then:   rustle --vg --genome-fasta <out_subset.fa> -L <reads.bam> -o out.gtf
#
# Requires: samtools on PATH. Uses the BAM index (idxstats, instant) + the FASTA
# .fai (faidx, seeks directly to the needed contigs); neither reads the whole
# genome.
set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "usage: $0 <genome.fa> <reads.bam> <out_subset.fa>" >&2
    exit 2
fi
GEN=$1; BAM=$2; OUT=$3

command -v samtools >/dev/null || { echo "error: samtools not on PATH" >&2; exit 1; }
[ -f "$GEN" ] || { echo "error: genome FASTA not found: $GEN" >&2; exit 1; }
[ -f "$BAM" ] || { echo "error: BAM not found: $BAM" >&2; exit 1; }

# Ensure indexes exist (cheap; created once).
[ -f "$GEN.fai" ] || samtools faidx "$GEN"
if [ ! -f "$BAM.bai" ] && [ ! -f "${BAM%.bam}.bai" ] && [ ! -f "$BAM.csi" ]; then
    samtools index "$BAM"
fi

# Contigs with >=1 mapped read (idxstats col 1 = contig, col 3 = mapped reads),
# intersected with contigs actually present in the genome FASTA.
genome_contigs=$(cut -f1 "$GEN.fai" | sort -u)
mapped_contigs=$(samtools idxstats "$BAM" | awk '$3 > 0 {print $1}' | sort -u)
contigs=$(comm -12 <(echo "$genome_contigs") <(echo "$mapped_contigs") | tr '\n' ' ')

if [ -z "${contigs// /}" ]; then
    echo "error: no mapped contigs in $BAM are present in $GEN" >&2
    exit 1
fi

# Extract just those contigs (faidx seeks via the .fai; no full-genome scan).
# shellcheck disable=SC2086
samtools faidx "$GEN" $contigs > "$OUT"
samtools faidx "$OUT"

n=$(echo $contigs | wc -w)
echo "subset $n contig(s) [$contigs] -> $OUT ($(du -h "$OUT" | cut -f1))" >&2
