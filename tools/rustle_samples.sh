#!/bin/bash
# rustle_samples.sh — multi-sample runs: pool, then split (docs/superpowers/specs/2026-11-02-multi-sample-design.md).
# Called by tools/rustle_pipeline.sh when --bam holds a comma; the same arguments plus --samples:
#   tools/rustle_pipeline.sh all  --bam a.bam,b.bam [--samples A,B] --fasta G --out PREFIX [options]
#   tools/rustle_pipeline.sh pool --bam a.bam,b.bam [--samples A,B] --out PREFIX            (the pooled BAM only)
# copy-table --samples A,B --out PREFIX rebuilds PREFIX.copy_counts.tsv from PREFIX.samples/.
# all = pool -> assemble + families on PREFIX.pooled.bam -> assign once per sample on the shared catalog
# (PREFIX.samples/<S>.assign.*) -> PREFIX.counts.tsv (utilities sample-counts on PREFIX.final.gtf, else PREFIX.gtf) ->
# PREFIX.copy_counts.tsv. BAMs must share the genome (@SQ lines), be indexed, hold mapped reads and share no read name.
# v1 limits: no --candidates, --legacy-catalog, --as-table, --to / --from; the flag stage is not run.
set -euo pipefail
STAGE=${1:?stage}; shift
say() { echo "[rustle_samples] $*" >&2; [ -n "${OUT:-}" ] && echo "$(date '+%F %T') $*" >> "$OUT.samples.log"; return 0; }
die() { echo "[rustle_samples] $*" >&2; exit 2; }
BAMS=""; NAMES=""; OUT=""; BIN="$(dirname "$0")/../target/release"; THREADS=4; REST=()
while [ $# -gt 0 ]; do
  case "$1" in
    --bam) BAMS=$2; shift 2;; --samples) NAMES=$2; shift 2;; --out) OUT=$2; shift 2;;
    --bin) BIN=$2; REST+=("$1" "$2"); shift 2;; --threads) THREADS=$2; REST+=("$1" "$2"); shift 2;;
    --candidates|--legacy-catalog|--to|--from|--as-table) die "$1 is not supported in multi-sample runs (v1; the as-table is built from the pooled BAM)";;
    *) REST+=("$1"); shift;;
  esac
done
case "$STAGE" in all|pool|copy-table) ;; *) die "multi-sample runs use the stages all or pool (got '$STAGE')";; esac
[ -n "$OUT" ] || die "need --out"
# PREFIX.copy_counts.tsv from PREFIX.samples/<S>.assign.quant.tsv: family_id copy_index copy_tid, then <S>_reads <S>_abundance
# (n_reads_hard, abundance); a sample without a quant file (or a copy it lacks) gets 0, a family on that sample's assign skip
# list (PREFIX.samples/<S>.assign.skipped.txt: --max-window-records on ITS reads, --low-mem, RUSTLE_ANNOTATED_ONLY,
# RUSTLE_EXCLUDE_BED) gets NA: it was not assigned, which is not zero reads.
copy_table() {
  local names=$1 SD=$OUT.samples tmp=$OUT.copy_counts.tsv.tmp n i=0
  : > "$tmp.ids"
  IFS=',' read -r -a NM <<< "$names"
  for n in "${NM[@]}"; do
    i=$((i + 1)); f=$SD/$n.assign.quant.tsv
    if [ -s "$f" ]; then
      awk -F'\t' -v OFS='\t' -v k=$i 'FNR == 1 { for (j = 1; j <= NF; j++) c[$j] = j; next }
        { print $c["family_id"], $c["copy_index"], $c["copy_tid"], k, $c["n_reads_hard"], $c["abundance"] }' "$f" >> "$tmp.ids"
    else
      say "copy table: no $f, sample $n counted 0"
    fi
    if [ -s "$SD/$n.assign.skipped.txt" ]; then awk -v OFS='\t' -v k=$i 'NF { print "#skip", $1, k }' "$SD/$n.assign.skipped.txt" >> "$tmp.ids"; fi
  done
  awk -F'\t' -v OFS='\t' -v m=$i -v names="$names" '
    $1 == "#skip" { sk[$2, $3] = 1; next }
    { id = $1 OFS $2 OFS $3; if (!(id in seen)) { seen[id] = 1; order[++n] = id; fam[id] = $1 }; r[id, $4] = $5; a[id, $4] = $6 }
    END { split(names, nm, ","); h = "family_id" OFS "copy_index" OFS "copy_tid"
          for (j = 1; j <= m; j++) h = h OFS nm[j] "_reads" OFS nm[j] "_abundance"; print h
          for (x = 1; x <= n; x++) { l = order[x]
            for (j = 1; j <= m; j++) { z = ((fam[order[x]], j) in sk) ? "NA" : 0
              l = l OFS ((order[x], j) in r ? r[order[x], j] : z) OFS ((order[x], j) in a ? a[order[x], j] : z) }
            print l } }' "$tmp.ids" > "$tmp"
  rm -f "$tmp.ids"; mv -f "$tmp" "$OUT.copy_counts.tsv"
}
if [ "$STAGE" = copy-table ]; then [ -n "$NAMES" ] || die "copy-table needs --samples"; copy_table "$NAMES"; exit 0; fi
command -v samtools > /dev/null || die "samtools not on PATH"
IFS=',' read -r -a B <<< "$BAMS"
if [ -n "$NAMES" ]; then IFS=',' read -r -a S <<< "$NAMES"
else S=(); for b in "${B[@]}"; do n=$(basename "$b"); S+=("${n%.bam}"); done; fi
[ "${#S[@]}" = "${#B[@]}" ] || die "--samples names ${#S[@]} samples for ${#B[@]} BAMs"
declare -A SEEN=()
for n in "${S[@]}"; do
  [[ "$n" =~ ^[A-Za-z0-9._-]+$ ]] || die "sample name '$n' must match [A-Za-z0-9._-]+"
  [ -z "${SEEN[$n]:-}" ] || die "sample name '$n' is used twice (BAM basenames collide: pass --samples)"
  SEEN[$n]=1
done
mkdir -p "$(dirname "$OUT")"
sq0=""
for i in "${!B[@]}"; do
  b=${B[$i]}
  [ -r "$b" ] || die "cannot read BAM '$b'"
  # the binaries open BAM.bai and nothing else (a .csi or x.bai passed here, then failed in every per-sample stage)
  [ -f "$b.bai" ] || die "BAM '$b' has no $b.bai: samtools index it"
  sq=$(samtools view -H "$b" | grep '^@SQ' | md5sum | cut -d' ' -f1)
  [ -z "$sq0" ] || [ "$sq" = "$sq0" ] || die "BAM '$b' (${S[$i]}) has other @SQ lines than '${B[0]}': samples must share the genome"
  sq0=$sq
  [ "$(samtools idxstats "$b" | awk '{m += $3} END {print m + 0}')" -gt 0 ] || die "BAM '$b' (${S[$i]}) holds no mapped reads"
done
key=$(for b in "${B[@]}"; do printf '%s\t%s\t%s\n' "$(readlink -f "$b")" "$(stat -c %s "$b")" "$(stat -c %Y "$b")"; done)
POOLED=$OUT.pooled.bam
if [ -s "$POOLED" ] && [ -f "$POOLED.bai" ] && [ "$(cat "$OUT.pool.key" 2>/dev/null)" = "$key" ]; then
  say "pool: $POOLED reused (same BAMs)"
else
  say "pool: read-name check over ${#B[@]} samples"
  # `|| exit 3` per BAM: the loop's status is its LAST iteration's, and set -e is off inside $(...)
  dup=$(for i in "${!B[@]}"; do samtools view -F 0x900 "${B[$i]}" | cut -f1 | LC_ALL=C sort -u -T "${TMPDIR:-/tmp}" | awk -v s="${S[$i]}" '{print $1 "\t" s}' || exit 3; done |
        LC_ALL=C sort -k1,1 -T "${TMPDIR:-/tmp}" | awk -F'\t' '!hit && $1 == p && $2 != q {hit = $1 " in " q " and " $2} {p = $1; q = $2} END {print hit}') \
    || die "read-name check failed (samtools view / sort; see above)"
  [ -z "$dup" ] || die "read name $dup: samples must not share read names"
  rm -f "$OUT.pool.key"
  # `a && b` is exempt from set -e: without the die a failed merge went on to index and key the OLD pool
  { samtools merge -f -@ "$THREADS" -o "$POOLED.tmp" "${B[@]}" && mv -f "$POOLED.tmp" "$POOLED"; } || { rm -f "$POOLED.tmp"; die "samtools merge failed: $POOLED not updated"; }
  samtools index -@ "$THREADS" "$POOLED"
  printf '%s\n' "$key" > "$OUT.pool.key"
  say "pool: $POOLED from ${S[*]}"
fi
[ "$STAGE" = pool ] && exit 0

P=$(dirname "$0")/rustle_pipeline.sh
say "all: assemble + families on the pooled reads"
# the stages one by one, not `all --to families`: run_all times each stage with /usr/bin/time, absent on some compute nodes
for st in assemble families; do bash "$P" "$st" --bam "$POOLED" --out "$OUT" "${REST[@]}"; done
SD=$OUT.samples; mkdir -p "$SD"
for i in "${!B[@]}"; do
  n=${S[$i]}
  rm -f "$SD/$n".gtf "$SD/$n".families.gtf "$SD/$n".fam.* "$SD/$n".splice.mmi
  for f in "$OUT".gtf "$OUT".families.gtf "$OUT".fam.* "$OUT".splice.mmi; do
    [ -e "$f" ] && ln -sfn "$(readlink -f "$f")" "$SD/$n${f#"$OUT"}"
  done
  say "assign: $n"
  # the catalog comes from the pooled reads: a copy silent in this sample is expected (copy_assign RUSTLE_SHARED_CATALOG)
  RUSTLE_SHARED_CATALOG=1 bash "$P" assign --bam "${B[$i]}" --out "$SD/$n" "${REST[@]}"
done
GTF=$OUT.final.gtf; [ -s "$GTF" ] && [ ! "$GTF" -ot "$OUT.gtf" ] || GTF=$OUT.gtf
say "counts: per-sample reads on $GTF"
"$BIN/utilities" sample-counts --gtf "$GTF" --bam "$(IFS=,; echo "${B[*]}")" --samples "$(IFS=,; echo "${S[*]}")" \
  --out "$OUT.counts.tsv" 2>> "$OUT.samples.log"
copy_table "$(IFS=,; echo "${S[*]}")"
say "done: $OUT.counts.tsv, $OUT.copy_counts.tsv, per-sample assign in $SD/ (flag stage not run in multi-sample mode)"
