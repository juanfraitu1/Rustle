#!/bin/bash
# Substrate prep for the human chr20 ordinary-chromosome assembler bakeoff
# (bench/CHR20_ASSEMBLER_COMPARISON.md): extracts chr20 from the genome-wide human testis IsoSeq BAM,
# the T2T-CHM13 genome FASTA, and the RefSeq annotation, then converts the annotation GFF3 -> GTF for
# gffcompare/SQANTI3. Run once, before bakeoff_chr20_{ours,stringtie,flair}.sh and chr20_score.sh.
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
SRC_BAM=/mnt/linuxdisk/home/juanfraitu/_from_wsl/human_val/human_testis.t2t.bam
GENOME=/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0.fa
GFF=/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz
mkdir -p "$W"; cd "$W"

# 1. chr20 BAM (coordinate-sorted subset, re-indexed).
if [ ! -s chr20.bam ]; then
  samtools view -b "$SRC_BAM" chr20 > chr20.bam
  samtools index chr20.bam
fi
samtools quickcheck chr20.bam && echo "chr20.bam OK: $(samtools view -c chr20.bam) records"

# 2. chr20 genome FASTA (+.fai).
if [ ! -s chr20.fa ]; then
  samtools faidx "$GENOME" chr20 > chr20.fa
  samtools faidx chr20.fa
fi
cat chr20.fa.fai

# 3. chr20 reference annotation, GFF3 subset.
if [ ! -s chr20_ref.gff3 ]; then
  zcat "$GFF" | awk -F'\t' '$1=="chr20"' > chr20_ref.gff3
fi
wc -l chr20_ref.gff3

# 4. RE-SORT before GFF3->GTF conversion. The upstream RefSeq GFF3 lists a feature's leftmost EXON
#    before its own GENE/TRANSCRIPT record whenever they share the same start coordinate (confirmed:
#    this ordering is already present in the un-filtered genome-wide chm13v2.0_RefSeq_full.gff.gz, not
#    introduced by our chr20 filter). gffread's single-pass GFF3 parser needs the parent seen before the
#    child to propagate the true gene ID into the GTF's `gene_id`; without this resort, EVERY converted
#    transcript's gene_id silently collapses to its own transcript_id (verified on this exact file: a
#    minimal 2-line extract of one gene+transcript converts correctly in isolation, but the same lines
#    convert WRONG when left in the full 99,049-row file in original order). Stable sort key: (chrom,
#    start, feature-rank[gene/pseudogene=0, mRNA/transcript/other top-level=1, exon/CDS/UTR=2], original
#    line order) -- this is a topological fix for the common "same-start" case, not a full Parent-graph
#    topological sort.
if [ ! -s chr20_ref.sorted.gff3 ]; then
  python3 - chr20_ref.gff3 chr20_ref.sorted.gff3 << 'PYEOF'
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
if [ ! -s chr20_ref.gtf ]; then
  gffread chr20_ref.sorted.gff3 -T -o chr20_ref.gtf 2> gffread.log
fi
echo "chr20_ref.gtf: $(wc -l < chr20_ref.gtf) rows, $(awk -F'\t' '$3=="transcript"' chr20_ref.gtf | wc -l) transcripts"
# sanity: gene_id must not still be silently equal to transcript_id for the general case (a small residual
# is expected and correct -- genes with no explicit mRNA/transcript child, where the exon's Parent points
# directly at the gene; gffread then uses the gene's own ID as both).
awk -F'\t' '$3=="transcript"' chr20_ref.gtf | perl -ne '$tid=$1 if /transcript_id "([^"]+)"/; $gid=$1 if /gene_id "([^"]+)"/; print if $tid eq $gid' | wc -l
