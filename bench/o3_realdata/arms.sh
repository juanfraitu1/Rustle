#!/bin/bash
# arms.sh — the masked and added-back BAMs of a real-data O3 build (prereg docs/archive/2026-10/merged/PREREG_o3_realdata_2026-10-15.md, Method),
# one step per call so each heavy step fits the 10-minute window under tools/rlock.sh:
#
#   arms.sh mask  DIR              DIR/masked.fa (+ .fai): GGO.fasta with DIR/mask.bed hard-masked to N; checks every span is N
#   arms.sh index DIR              DIR/masked.splice.mmi (minimap2 -x splice:hq -d)
#   arms.sh split DIR K            DIR/reads.part{0..K-1}.fa, round-robin over DIR/reads.fa
#   arms.sh map   DIR ARM PART     ARM R (masked index) or R0 (GGO.splice.mmi): DIR/ARM.pPART.bam, the baseline's exact command
#   arms.sh merge DIR ARM K        DIR/ARM.bam (+ .bai) from the K parts
set -euo pipefail
GENOME=/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta
UNMASKED=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi
step=${1:?mask|index|split|map|merge}; D=${2:?DIR}
case $step in
  mask)
    bedtools maskfasta -fi "$GENOME" -bed "$D/mask.bed" -fo "$D/masked.fa"
    samtools faidx "$D/masked.fa"
    while IFS=$'\t' read -r c s e f; do
      n=$(samtools faidx "$D/masked.fa" "$c:$((s + 1))-$e" | tail -n +2 | tr -d 'Nn\n' | wc -c)
      [ "$n" -eq 0 ] || { echo "span $f $c:$s-$e not fully masked ($n bases)" >&2; exit 1; }
    done < "$D/mask.bed"
    echo "masked $(wc -l < "$D/mask.bed") spans, all N" ;;
  index) minimap2 -x splice:hq -t 4 -d "$D/masked.splice.mmi" "$D/masked.fa" ;;
  split)
    k=${3:?K}
    awk -v k="$k" -v d="$D" '/^>/{n++} {print > (d "/reads.part" ((n - 1) % k) ".fa")}' "$D/reads.fa"
    grep -c '^>' "$D"/reads.part*.fa ;;
  map)
    arm=${3:?R|R0}; p=${4:?PART}
    case $arm in R) idx=$D/masked.splice.mmi ;; R0) idx=$UNMASKED ;; *) echo "ARM must be R or R0" >&2; exit 2 ;; esac
    minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4 "$idx" "$D/reads.part$p.fa" 2> "$D/$arm.p$p.bam.mm2.log" \
      | samtools sort -@2 -o "$D/$arm.p$p.bam" - ;;
  merge)
    arm=${3:?R|R0}; k=${4:?K}
    parts=(); for ((i = 0; i < k; i++)); do parts+=("$D/$arm.p$i.bam"); done
    samtools merge -f -@2 "$D/$arm.bam" "${parts[@]}"
    samtools index "$D/$arm.bam"
    echo "$arm: $(samtools view -c -F 2308 "$D/$arm.bam") primary, $(samtools view -c -f 4 "$D/$arm.bam") unmapped" ;;
  *) echo "arms.sh: unknown step $step" >&2; exit 2 ;;
esac
