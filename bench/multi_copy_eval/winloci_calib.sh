#!/usr/bin/env bash
# Calibration: on a fixed region, run StringTie / baseline / several VG configs
# and report how many RefSeq transcripts each matches (class '='). Used to pin
# the correct VG "win" config before fanning out across candidates.
# Usage: winloci_calib.sh CHROM START END LABEL
set -uo pipefail
CHROM="$1"; START="$2"; END="$3"; LABEL="${4:-calib}"
ROOT=/mnt/c/Users/jfris/Desktop/Rustle
SCRATCH="${SCRATCH:-/home/juanfra/winloci_scratch}"
GGO="${GGO:-$SCRATCH/GGO.bam}"
GFF="${GFF:-$SCRATCH/GGO_genomic.gff}"
FASTA="${FASTA:-$SCRATCH/GGO.fasta}"
RUSTLE="$ROOT/target/release/rustle"; ST="$ROOT/tools/stringtie/stringtie"
W="/tmp/calib/$LABEL"; rm -rf "$W"; mkdir -p "$W"
samtools view -b "$GGO" "$CHROM:$START-$END" 2>/dev/null | samtools sort -o "$W/l.bam" - 2>/dev/null
samtools index "$W/l.bam"
echo "region $CHROM:$START-$END  primary=$(samtools view -c -F0x900 $W/l.bam) sec=$(samtools view -c -f0x100 $W/l.bam)"

run() { # name  gtf  [env...]
  local name="$1"; local gtf="$2"; shift 2
  env "$@" "$RUSTLE" --vg --vg-snp --genome-fasta "$FASTA" -L "$W/l.bam" -o "$gtf" >/dev/null 2>&1 || true
}
eqn() { local f; f=$(ls $W/$1.*.tmap 2>/dev/null|head -1); [ -n "$f" ] && awk -F'\t' 'NR>1&&$3=="="{print $2}' "$f"|sort -u; }

"$ST" -L "$W/l.bam" -o "$W/st.gtf" >/dev/null 2>&1 || true
"$RUSTLE" -L "$W/l.bam" -o "$W/base.gtf" >/dev/null 2>&1 || true
run vgtan   "$W/vgtan.gtf"   RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1
run vgtandg "$W/vgtandg.gtf" RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1 RUSTLE_VG_DECISIVE_GATE=1
run vgplain "$W/vgplain.gtf"
env RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1 "$RUSTLE" --vg --vg-snp --read-chain --genome-fasta "$FASTA" -L "$W/l.bam" -o "$W/vgtanrc.gtf" >/dev/null 2>&1 || true

for t in st base vgtan vgtandg vgplain vgtanrc; do
  [ -s "$W/$t.gtf" ] && gffcompare -r "$GFF" "$W/$t.gtf" -o "$W/$t" >/dev/null 2>&1 || true
  eqn "$t" > "$W/$t.eq" 2>/dev/null || true
done
txn() { grep -cP '\ttranscript\t' "$W/$1.gtf" 2>/dev/null || echo 0; }
printf "%-10s %6s %6s\n" tool tx ref_eq
for t in st base vgtan vgtandg vgplain vgtanrc; do
  printf "%-10s %6s %6s\n" "$t" "$(txn $t)" "$(wc -l < $W/$t.eq 2>/dev/null||echo 0)"
done
echo "-- vgtan GAIN vs base:"; comm -23 "$W/vgtan.eq" "$W/base.eq" | sed 's/^/   +/'
echo "-- vgtan LOSS vs base:"; comm -13 "$W/vgtan.eq" "$W/base.eq" | sed 's/^/   -/'
echo "-- vgtanrc GAIN vs base:"; comm -23 "$W/vgtanrc.eq" "$W/base.eq" | sed 's/^/   +/'
