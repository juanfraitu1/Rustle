#!/bin/bash
# Anti-overfit validation for mis-chain read salvage. GATE (byte-identity + known families) runs FIRST;
# target recovery is measured LAST. Controller-run (heavy binary, WSL crash rule): outputs to winloci_scratch.
set -u
BIN=/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog
FA=/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa
SAM=/home/juanfra/miniforge3/bin/samtools
SCR=/home/juanfra/winloci_scratch; CACHE=$SCR/soto_cache; PC=$CACHE/perchrom
OUT=$CACHE/salvage_val; mkdir -p "$OUT"
ncopies(){ echo $(($(wc -l < "$1" 2>/dev/null || echo 1)-1)); }
overlaps(){ awk -F'\t' -v c="$2" -v s="$3" -v e="$4" 'NR>1 && $4==c && !($5>e||$6<s){print $5"-"$6" ("$9"r)";f=1} END{exit !f}' "$1"; }

echo "######## PART A — OFF byte-identity (chr9) ########"
timeout 600 "$BIN" --bam "$PC/chr9.bam" --fasta "$FA" --out "$OUT/chr9_off1" >/dev/null 2>&1
timeout 600 "$BIN" --bam "$PC/chr9.bam" --fasta "$FA" --out "$OUT/chr9_off2" >/dev/null 2>&1
RUSTLE_MISCHAIN_SALVAGE=1 timeout 600 "$BIN" --bam "$PC/chr9.bam" --fasta "$FA" --out "$OUT/chr9_on" >/dev/null 2>&1
echo "OFF run1 md5: $(md5sum < "$OUT/chr9_off1.copies.tsv")"
echo "OFF run2 md5: $(md5sum < "$OUT/chr9_off2.copies.tsv")   (must equal run1 = deterministic OFF)"
echo "ON      md5: $(md5sum < "$OUT/chr9_on.copies.tsv")   (should DIFFER = salvage active)"
echo "OFF copies: $(ncopies "$OUT/chr9_off1.copies.tsv")   ON copies: $(ncopies "$OUT/chr9_on.copies.tsv")"

echo; echo "######## PART B — known families: ON must EQUAL OFF (no over-splitting) ########"
for f in GSTM MAGEA DAZ RBMY TSPY PCDHB; do
  bam=$SCR/reinv_${f}_orig.bam; [ -f "$bam" ] || { echo "$f: NO BAM"; continue; }
  timeout 400 "$BIN" --bam "$bam" --fasta "$FA" --out "$OUT/${f}_off" >/dev/null 2>&1
  RUSTLE_MISCHAIN_SALVAGE=1 timeout 400 "$BIN" --bam "$bam" --fasta "$FA" --out "$OUT/${f}_on" >/dev/null 2>&1
  off=$(ncopies "$OUT/${f}_off.copies.tsv"); on=$(ncopies "$OUT/${f}_on.copies.tsv")
  echo "  $f: OFF=$off ON=$on  $([ "$off" = "$on" ] && echo 'SAME (ok)' || echo '*** DIFF — INVESTIGATE ***')"
done

echo; echo "######## PART C — target recovery on confirmed cases (ON vs OFF) ########"
# CDH12 (chr5, both loci same chrom)
$SAM view -b "$PC/chr5.bam" chr5:70800000-71700000 -o "$OUT/cdh12.bam" 2>/dev/null && $SAM index "$OUT/cdh12.bam"
# AC134878.2(chr22)+TEKT4P2(chr21) cross-chrom
$SAM view -b "$PC/chr22.bam" chr22:5780000-5880000 -o "$OUT/_a.bam" 2>/dev/null
$SAM view -b "$PC/chr21.bam" chr21:5670000-5770000 -o "$OUT/_b.bam" 2>/dev/null
$SAM merge -f "$OUT/ac13.bam" "$OUT/_a.bam" "$OUT/_b.bam" 2>/dev/null && $SAM index "$OUT/ac13.bam"
# AC126603.1(chr15)+AC142384.1(chr16) cross-chrom
$SAM view -b "$PC/chr15.bam" chr15:18100000-18200000 -o "$OUT/_c.bam" 2>/dev/null
$SAM view -b "$PC/chr16.bam" chr16:34860000-34950000 -o "$OUT/_d.bam" 2>/dev/null
$SAM merge -f "$OUT/ac12.bam" "$OUT/_c.bam" "$OUT/_d.bam" 2>/dev/null && $SAM index "$OUT/ac12.bam"
# NCF1/B/C (chr7)
$SAM view -b "$PC/chr7.bam" chr7:74000000-76500000 -o "$OUT/ncf1.bam" 2>/dev/null && $SAM index "$OUT/ncf1.bam"

check(){ # name bam xflag "chrom:s-e:label ..."
  local name=$1 bam=$2 xf=$3; shift 3
  timeout 600 "$BIN" --bam "$bam" --fasta "$FA" $xf --out "$OUT/${name}_off" >/dev/null 2>&1
  RUSTLE_MISCHAIN_SALVAGE=1 timeout 600 "$BIN" --bam "$bam" --fasta "$FA" $xf --out "$OUT/${name}_on" >/dev/null 2>&1
  echo "-- $name: OFF=$(ncopies "$OUT/${name}_off.copies.tsv")cp  ON=$(ncopies "$OUT/${name}_on.copies.tsv")cp"
  for t in "$@"; do IFS=: read c s e lab <<< "$t"
    o=$(overlaps "$OUT/${name}_off.copies.tsv" "$c" "$s" "$e" && echo Y || echo N)
    n=$(overlaps "$OUT/${name}_on.copies.tsv"  "$c" "$s" "$e" && echo Y || echo N)
    echo "     $lab: OFF=$o ON=$n  $([ "$o" = N ] && [ "$n" = Y ] && echo '<< RECOVERED' || true)"
  done
}
check cdh12 "$OUT/cdh12.bam" "" "chr5:70901860:70902156:CDH12P1" "chr5:71589693:71589989:CDH12P3"
check ac13  "$OUT/ac13.bam"  "--cross-chrom" "chr22:5801017:5862418:AC134878.2" "chr21:5692610:5754003:TEKT4P2"
check ac12  "$OUT/ac12.bam"  "--cross-chrom" "chr15:18136928:18186041:AC126603.1" "chr16:34884759:34929924:AC142384.1"
check ncf1  "$OUT/ncf1.bam"  "" "chr7:75976253:75991692:NCF1" "chr7:74420835:74436153:NCF1B" "chr7:76360590:76375995:NCF1C"
echo; echo "######## DONE ########"
