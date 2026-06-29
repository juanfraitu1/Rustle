#!/usr/bin/env bash
# Config ladder on an already-extracted locus.bam: localize which VG stage drops
# transcripts. Prints, per config: #tx emitted and #RefSeq matched (class '='),
# plus the matched-ref set so callers can diff.
# Usage: winloci_ladder.sh GENE   (uses /tmp/winloci/GENE/locus.bam)
set -uo pipefail
G="${1:?gene}"
SC=/home/juanfra/winloci_scratch
ROOT=/mnt/c/Users/jfris/Desktop/Rustle
R="$ROOT/target/release/rustle"; ST="$ROOT/tools/stringtie/stringtie"
BAM="/tmp/winloci/$G/locus.bam"; W="/tmp/ladder/$G"; rm -rf "$W"; mkdir -p "$W"
[ -s "$BAM" ] || { echo "no bam $BAM"; exit 1; }

run() { local name="$1"; shift; local last="${!#}"; local gtf="$W/$name.gtf"
  env "${@:1:$#-1}" "$last" >/dev/null 2>&1 || true; }

# config -> command (last arg is the rustle/ST invocation as a string via bash -c)
declare -A CFG
CFG[0_base]="$ST -L $BAM -o $W/0_base.gtf"
CFG[1_rustlebase]="$R -L $BAM -o $W/1_rustlebase.gtf"
CFG[2_vg]="$R --vg --vg-snp --genome-fasta $SC/GGO.fasta -L $BAM -o $W/2_vg.gtf"
CFG[3_vgtandem]="env RUSTLE_VG_TANDEM=1 $R --vg --vg-snp --genome-fasta $SC/GGO.fasta -L $BAM -o $W/3_vgtandem.gtf"
CFG[4_vgtandem_gate]="env RUSTLE_VG_TANDEM=1 RUSTLE_VG_TANDEM_PRIMARY_JUNCTIONS=1 $R --vg --vg-snp --genome-fasta $SC/GGO.fasta -L $BAM -o $W/4_vgtandem_gate.gtf"

eqn(){ local f; f=$(ls $W/$1.*.tmap 2>/dev/null|head -1); [ -n "$f" ]&&awk -F'\t' 'NR>1&&$3=="="{print $2}' "$f"|sort -u; }
txn(){ grep -cP '\ttranscript\t' "$W/$1.gtf" 2>/dev/null||echo 0; }

for k in 0_base 1_rustlebase 2_vg 3_vgtandem 4_vgtandem_gate; do
  bash -c "${CFG[$k]}" >/dev/null 2>&1 || true
  [ -s "$W/$k.gtf" ] && gffcompare -r "$SC/GGO_tx.gff" "$W/$k.gtf" -o "$W/$k" >/dev/null 2>&1 || true
  eqn "$k" > "$W/$k.eq"
done
echo "config              tx  refeq   (matched set vs 1_rustlebase: -=lost +=gained)"
base="$W/1_rustlebase.eq"
for k in 0_base 1_rustlebase 2_vg 3_vgtandem 4_vgtandem_gate; do
  printf "%-18s %4s %5s   " "$k" "$(txn $k)" "$(wc -l < $W/$k.eq)"
  if [ "$k" != "1_rustlebase" ] && [ "$k" != "0_base" ]; then
    lost=$(comm -23 "$base" "$W/$k.eq"|paste -sd,); gain=$(comm -13 "$base" "$W/$k.eq"|paste -sd,)
    printf "lost[%s] gained[%s]" "$lost" "$gain"
  fi
  echo
done
