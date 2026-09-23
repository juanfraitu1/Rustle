#!/usr/bin/env bash
# ⚠ 2026-09-23 (§6zb): for `--assemble-only` this script is SUPERSEDED by `copy_assign --genome-wide` — one
# process, streaming pass-1, per-contig polish, ~2 GB peak for the whole human genome and byte-identical
# to this sweep's concatenated output apart from the O2-only `matched_reads` attribute. Kept for the
# non-streaming configurations (`--read-isoform-k`, `RUSTLE_FOOTPRINT_NODES`, `--materialize-reads`).
#
# Genome-wide `--assemble-only` sweep, one process per contig, with bounded concurrency.
#
# WHY PROCESSES AND NOT `--region-threads` (§6r8, register 880/881):
#   * `--region-threads 4` measured 25m44s against a 26m11s serial baseline -- a 1.02x speedup -- because
#     the parallel `compute` map is a small share of wall time and the serial drain dominates.
#   * Worse, it is NOT byte-identical despite its docstring saying so: one molecule flipped `tied` ->
#     `assigned` at an exact tie (margin 0.000), and `families.tsv` moved with it.
#   * One process per contig gave 1.76x on 4 chromosomes (426s -> 242s) with every GTF byte-identical.
#
# MEMORY WAS THE BINDING CONSTRAINT AND IS NOW MOSTLY GONE. Peak RSS used to be 8.4-11.8 GB for ONE
# chromosome -- read bases and qualities held for every alignment -- so 4 concurrent chromosomes
# OOM-killed a 25 GB box (exit 137). `--assemble-only` now drops both at parse time (assembly never reads
# them; only O2 does), taking peak RSS to 3.5-4.5 GB. 4 concurrent chromosomes then fit and give 3.74x
# (426s -> 114s, byte-identical). Budget ~4.5 GB per slot.
#
# usage: genome_wide_sweep.sh --bam IN.bam --fasta GENOME.fa --outdir DIR [--jobs 4] [--contigs "chr1 chr2"]
set -euo pipefail

BAM=""; FA=""; OUT=""; JOBS=4; CONTIGS=""
while [ $# -gt 0 ]; do
    case "$1" in
        --bam) BAM=$2; shift 2 ;;
        --fasta) FA=$2; shift 2 ;;
        --outdir) OUT=$2; shift 2 ;;
        --jobs) JOBS=$2; shift 2 ;;
        --contigs) CONTIGS=$2; shift 2 ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done
[ -n "$BAM" ] && [ -n "$FA" ] && [ -n "$OUT" ] || { sed -n '2,20p' "$0"; exit 2; }

BIN=${COPY_ASSIGN:-target/release/copy_assign}
[ -x "$BIN" ] || BIN=$(command -v copy_assign) || { echo "copy_assign not found; set \$COPY_ASSIGN" >&2; exit 127; }
# 2026-09-23 (§6za): strict junctions + retained-intron ratio 10 are the binary defaults under --assemble-only; listed here explicitly.
POLISH=${POLISH:-"--assembly-junctions strict --assembly-polish full --polish-isoform-fraction 0.02 --polish-mono-shadow --polish-mono-quantile 0.82 --polish-ism-ratio 0.7 --polish-retained-ratio 10"}

mkdir -p "$OUT"
# contigs with at least one alignment, longest first so the stragglers start early
if [ -z "$CONTIGS" ]; then
    CONTIGS=$(samtools idxstats "$BAM" | awk '$1!="*" && $3>0 {print $2"\t"$1}' | sort -k1,1nr | cut -f2 | tr '\n' ' ')
fi
declare -A LEN
while read -r name len _; do [ -n "$name" ] && LEN[$name]=$len; done < <(samtools idxstats "$BAM" | awk '{print $1, $2, $3}')

echo "sweep: $(echo "$CONTIGS" | wc -w) contigs, ${JOBS} at a time, into $OUT/"
start=$(date +%s)
running=0
for c in $CONTIGS; do
    l=${LEN[$c]:-0}; [ "$l" -gt 0 ] || continue
    ( "$BIN" --assemble-only $POLISH --bam "$BAM" --fasta "$FA" \
        --region "$c:1-$l" --out "$OUT/$c" > "$OUT/$c.log" 2>&1 \
      && echo "  done $c ($(awk -F'\t' '$3=="transcript"' "$OUT/$c.gtf" 2>/dev/null | wc -l) transcripts)" \
      || echo "  FAILED $c (see $OUT/$c.log)" ) &
    running=$((running+1))
    if [ "$running" -ge "$JOBS" ]; then wait -n 2>/dev/null || wait; running=$((running-1)); fi
done
wait
echo "sweep wall: $(( $(date +%s) - start ))s"
echo "merge with: cat $OUT/*.gtf > $OUT/genome.gtf"
