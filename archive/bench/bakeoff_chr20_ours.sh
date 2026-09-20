#!/bin/bash
# "Ours" arm of the human chr20 ordinary-chromosome assembler bakeoff
# (advisor wants a standard gffcompare+SQANTI3 comparison on a NOT-multicopy-targeted human chromosome,
# see bench/CHR20_ASSEMBLER_COMPARISON.md). Sibling of bakeoff_flair.sh/bakeoff_stringtie.sh but targets
# human T2T-CHM13 chr20 IsoSeq instead of the gorilla NPIP substrate those target -- DO NOT overwrite those.
#
# Pure de novo detection: NO --families (no catalog supplied), --gtf only. The GTF's own `multicopy`
# attribute (true/false) comes from copy_assign's own de novo family/copy detection at --min-copies
# (default 2), not from any external catalog -- this is what bench/chr20_score.sh stratifies on.
set -euo pipefail
W=/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20
BAM="$W/chr20.bam"
FA="$W/chr20.fa"
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
mkdir -p "$W/ours"
cd "$W/ours"

[ -s "$BAM" ] || { echo "missing $BAM" >&2; exit 1; }

"$BIN" --gtf --bam "$BAM" --fasta "$FA" --region chr20:1-66210255 --out ours \
  > ours.stdout.log 2> ours.stderr.log
echo "exit=$?"
tail -20 ours.stderr.log
ls -la "$W/ours"
grep -vc '^#' ours.gtf || true
