#!/usr/bin/env bash
# Prove that a gene-family copy PRESENT in the reads but ABSENT from the reference genome
# is detectable. A hidden copy's reads mismap to the closest sibling, carrying their private
# SNPs as a COHERENT second haplotype the reference (1 copy at that locus) cannot explain.
# Signal = count of positions where 20-60% of reads carry a non-reference allele (an alt
# haplotype, not sequencing error). HIDDEN: many; ALL-VISIBLE: ~0.
set -euo pipefail
GEN=bench/tandem_attribution/gen_synthetic.py
WORK=${1:-/tmp/hidden_copy_probe}; rm -rf "$WORK"; mkdir -p "$WORK"
python3 "$GEN" --identity 0.97 --copies 3 --reads-per-copy 30 --hidden-copies "2" --spacing 200000 --seed 8 --out "$WORK/hidden"  >/dev/null
python3 "$GEN" --identity 0.97 --copies 3 --reads-per-copy 30                       --spacing 200000 --seed 8 --out "$WORK/visible" >/dev/null
for c in hidden visible; do
  minimap2 -ax splice:hq -uf -N20 -p0.5 "$WORK/$c/genome.fa" "$WORK/$c/reads.fa" 2>/dev/null | samtools sort -o "$WORK/$c/r.bam" - 2>/dev/null
  samtools index "$WORK/$c/r.bam"
  n=$(samtools mpileup -f "$WORK/$c/genome.fa" -r synTANDEM:1000-12600 "$WORK/$c/r.bam" 2>/dev/null | \
    awk '{d=$4; if(d<10)next; b=$5; gsub(/[.,]/,"",b); gsub(/[\^$].|[+-][0-9]+/,"",b); f=length(b)/d; if(f>=0.2&&f<=0.6)c++} END{print c+0}')
  dp=$(samtools view -c -F 0x900 "$WORK/$c/r.bam" synTANDEM:1000-12600 2>/dev/null)
  echo "$c: copy0-locus depth=$dp  alt-haplotype-positions=$n"
done
