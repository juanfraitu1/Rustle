#!/bin/bash
# Substrate prep for a human ordinary-chromosome assembler bakeoff (parameterized sibling of
# bench/prep_chr20_ref.sh -- DO NOT overwrite that one): extracts CHROM from the genome-wide human testis
# IsoSeq BAM, the T2T-CHM13 genome FASTA, and the RefSeq annotation, then converts the annotation GFF3 -> GTF
# for gffcompare/SQANTI3. Run once, before bakeoff_chrom_{ours,stringtie,flair}.sh and chrom_score.sh.
# usage: prep_chrom_ref.sh CHROM   (e.g. prep_chrom_ref.sh chr17)
set -euo pipefail
CHROM=${1:?usage: prep_chrom_ref.sh CHROM}
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_${CHROM}
SRC_BAM=/mnt/linuxdisk/home/juanfraitu/_from_wsl/human_val/human_testis.t2t.bam
GENOME=/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0.fa
GFF=/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz
mkdir -p "$W"; cd "$W"

# 1. CHROM BAM (coordinate-sorted subset, re-indexed).
if [ ! -s "${CHROM}.bam" ]; then
  samtools view -b "$SRC_BAM" "$CHROM" > "${CHROM}.bam"
  samtools index "${CHROM}.bam"
fi
samtools quickcheck "${CHROM}.bam" && echo "${CHROM}.bam OK: $(samtools view -c "${CHROM}.bam") records"

# 2. CHROM genome FASTA (+.fai).
if [ ! -s "${CHROM}.fa" ]; then
  samtools faidx "$GENOME" "$CHROM" > "${CHROM}.fa"
  samtools faidx "${CHROM}.fa"
fi
cat "${CHROM}.fa.fai"

# 3. CHROM reference annotation, GFF3 subset.
if [ ! -s "${CHROM}_ref.gff3" ]; then
  zcat "$GFF" | awk -F'\t' -v c="$CHROM" '$1==c' > "${CHROM}_ref.gff3"
fi
wc -l "${CHROM}_ref.gff3"

# 4. RE-SORT before GFF3->GTF conversion. The upstream RefSeq GFF3 lists a feature's leftmost EXON
#    before its own GENE/TRANSCRIPT record whenever they share the same start coordinate (confirmed:
#    this ordering is already present in the un-filtered genome-wide chm13v2.0_RefSeq_full.gff.gz, not
#    introduced by our chromosome filter). gffread's single-pass GFF3 parser needs the parent seen before the
#    child to propagate the true gene ID into the GTF's `gene_id`; without this resort, EVERY converted
#    transcript's gene_id silently collapses to its own transcript_id (verified on the chr20 instance of this
#    exact file: a minimal 2-line extract of one gene+transcript converts correctly in isolation, but the same
#    lines convert WRONG when left in the full 99,049-row file in original order). Stable sort key: (chrom,
#    start, feature-rank[gene/pseudogene=0, mRNA/transcript/other top-level=1, exon/CDS/UTR=2], original
#    line order) -- this is a topological fix for the common "same-start" case, not a full Parent-graph
#    topological sort.
if [ ! -s "${CHROM}_ref.sorted.gff3" ]; then
  python3 - "${CHROM}_ref.gff3" "${CHROM}_ref.sorted.gff3" << 'PYEOF'
import sys
RANK0 = {"gene", "pseudogene"}
RANK2 = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR", "start_codon", "stop_codon"}
def rank(feat):
    if feat in RANK0: return 0
    if feat in RANK2: return 2
    return 1
lines = []
with open(sys.argv[1]) as f:
    for i, line in enumerate(f):
        if line.startswith('#') or not line.strip():
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 5:
            continue
        lines.append((cols[0], int(cols[3]), rank(cols[2]), i, line))
lines.sort(key=lambda t: (t[0], t[1], t[2], t[3]))
with open(sys.argv[2], 'w') as out:
    for _, _, _, _, line in lines:
        out.write(line)
PYEOF
fi

# 5. GFF3 -> GTF (gffread, from the sqanti3 conda env -- SQANTI3's --refGTF requires GTF; gffcompare
#    accepts GFF3 directly but we use the same GTF for both tools for consistency).
source /home/juanfra/miniforge3/etc/profile.d/conda.sh; conda activate sqanti3
if [ ! -s "${CHROM}_ref.gtf" ]; then
  gffread "${CHROM}_ref.sorted.gff3" -T -o "${CHROM}_ref.gtf" 2> gffread.log
fi
echo "${CHROM}_ref.gtf: $(wc -l < "${CHROM}_ref.gtf") rows, $(awk -F'\t' '$3=="transcript"' "${CHROM}_ref.gtf" | wc -l) transcripts"
# sanity: gene_id must not still be silently equal to transcript_id for the general case (a small residual
# is expected and correct -- genes with no explicit mRNA/transcript child, where the exon's Parent points
# directly at the gene; gffread then uses the gene's own ID as both).
awk -F'\t' '$3=="transcript"' "${CHROM}_ref.gtf" | perl -ne '$tid=$1 if /transcript_id "([^"]+)"/; $gid=$1 if /gene_id "([^"]+)"/; print if $tid eq $gid' | wc -l
