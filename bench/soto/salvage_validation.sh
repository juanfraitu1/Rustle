#!/bin/bash
# Anti-overfit validation for mis-chain read salvage. GATE (byte-identity + known families) runs FIRST;
# target recovery is measured LAST. Controller-run (heavy binary, WSL crash rule): outputs to winloci_scratch.
# §6am: `set -u` alone let a 127 from a missing binary fall straight through to the verdicts below.
set -euo pipefail
# The mandated build dir is CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target; the repo-local
# ./target/release DOES NOT EXIST, so the old path here invoked NOTHING. Overridable for A/B-ing binaries.
BIN=${BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target/release/gw_family_catalog}
FA=/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa
SAM=/home/juanfra/miniforge3/bin/samtools
SCR=/home/juanfra/winloci_scratch; CACHE=$SCR/soto_cache; PC=$CACHE/perchrom
OUT=$CACHE/salvage_val
die(){ echo "ABORT: $*"; exit 2; }
vacuous(){ echo "ABORT (VACUOUS EVIDENCE): $*"; exit 4; }

# §6am PREFLIGHT — every input checked BEFORE any run, mkdir or truncation: a missing binary used to
# return 127 into >/dev/null, after which leftover outputs were md5'd and a byte-identity PASS printed.
[ -x "$BIN" ] || die "gw_family_catalog not executable at $BIN (build with CARGO_TARGET_DIR=/mnt/linuxdisk/home/juanfraitu/rustle_target --release, or pass BIN=)"
[ -x "$SAM" ] || die "samtools not executable at $SAM"
[ -s "$FA" ]  || die "reference FASTA missing/empty at $FA"
[ -d "$PC" ]  || die "per-chrom BAM cache missing at $PC"
for c in chr5 chr7 chr9 chr15 chr16 chr21 chr22; do
  [ -s "$PC/$c.bam" ] || die "$PC/$c.bam missing/empty — Parts A/C cannot be evaluated"
done

mkdir -p "$OUT"
# §6am: drop this script's own prior outputs, so no future failed run can be certified against leftovers
# (the Jul-24 chr9_off1/off2/on trio sitting here is exactly what the broken run md5-compared).
rm -f "$OUT"/*.copies.tsv "$OUT"/*.copies.fa "$OUT"/*.families.tsv "$OUT"/*.log "$OUT"/*.bam "$OUT"/*.bam.bai
ncopies(){ echo $(($(wc -l < "$1" 2>/dev/null || echo 1)-1)); }
overlaps(){ awk -F'\t' -v c="$2" -v s="$3" -v e="$4" 'NR>1 && $4==c && !($5>e||$6<s){print $5"-"$6" ("$9"r)";f=1} END{exit !f}' "$1"; }
# §6am: a comparison must never run against a file the current run did not produce.
need_out(){ [ -s "$1" ] || die "expected output $1 was not produced — any verdict below would score leftovers"; }

echo "######## PART A — OFF byte-identity (chr9) ########"
timeout 600 "$BIN" --bam "$PC/chr9.bam" --fasta "$FA" --out "$OUT/chr9_off1" >/dev/null 2>&1 || die "chr9 OFF run 1 failed (exit $?)"
timeout 600 "$BIN" --bam "$PC/chr9.bam" --fasta "$FA" --out "$OUT/chr9_off2" >/dev/null 2>&1 || die "chr9 OFF run 2 failed (exit $?)"
RUSTLE_MISCHAIN_SALVAGE=1 timeout 600 "$BIN" --bam "$PC/chr9.bam" --fasta "$FA" --out "$OUT/chr9_on" >/dev/null 2>&1 || die "chr9 ON run failed (exit $?)"
need_out "$OUT/chr9_off1.copies.tsv"; need_out "$OUT/chr9_off2.copies.tsv"; need_out "$OUT/chr9_on.copies.tsv"
# §6am: with 0 copies all three files are bare headers, so OFF1==OFF2 and ON==OFF are true by construction
# and the "should DIFFER" check cannot fail — an empty evidence set scoring as a byte-identity PASS.
[ "$(ncopies "$OUT/chr9_off1.copies.tsv")" -gt 0 ] || vacuous "chr9 OFF produced 0 copies; the md5 triple compares bare headers"
echo "OFF run1 md5: $(md5sum < "$OUT/chr9_off1.copies.tsv")"
echo "OFF run2 md5: $(md5sum < "$OUT/chr9_off2.copies.tsv")   (must equal run1 = deterministic OFF)"
echo "ON      md5: $(md5sum < "$OUT/chr9_on.copies.tsv")   (should DIFFER = salvage active)"
echo "OFF copies: $(ncopies "$OUT/chr9_off1.copies.tsv")   ON copies: $(ncopies "$OUT/chr9_on.copies.tsv")"

echo; echo "######## PART B — known families: ON must EQUAL OFF (no over-splitting) ########"
scored=0
for f in GSTM MAGEA DAZ RBMY TSPY PCDHB; do
  bam=$SCR/reinv_${f}_orig.bam; [ -f "$bam" ] || { echo "$f: NO BAM"; continue; }
  timeout 400 "$BIN" --bam "$bam" --fasta "$FA" --out "$OUT/${f}_off" >/dev/null 2>&1 || die "$f OFF run failed (exit $?)"
  RUSTLE_MISCHAIN_SALVAGE=1 timeout 400 "$BIN" --bam "$bam" --fasta "$FA" --out "$OUT/${f}_on" >/dev/null 2>&1 || die "$f ON run failed (exit $?)"
  need_out "$OUT/${f}_off.copies.tsv"; need_out "$OUT/${f}_on.copies.tsv"
  off=$(ncopies "$OUT/${f}_off.copies.tsv"); on=$(ncopies "$OUT/${f}_on.copies.tsv")
  # §6am: 0 copies in BOTH arms is not "SAME (ok)" — an over-split is undetectable there, so label it
  # and keep it out of the evidence count instead of printing a pass.
  if [ "$off" = 0 ] && [ "$on" = 0 ]; then
    echo "  $f: OFF=0 ON=0  *** NO COPIES IN EITHER ARM — VACUOUS, NOT A PASS ***"
  else
    scored=$((scored+1))
    echo "  $f: OFF=$off ON=$on  $([ "$off" = "$on" ] && echo 'SAME (ok)' || echo '*** DIFF — INVESTIGATE ***')"
  fi
done
# §6am: with no family carrying copies, Part B's "no over-splitting" claim rests on nothing at all.
[ "$scored" -gt 0 ] || vacuous "Part B scored 0 families with copies (all BAMs absent, or every arm empty)"
echo "  Part B evidence: $scored/6 families scored with >0 copies"

echo; echo "######## PART C — target recovery on confirmed cases (ON vs OFF) ########"
# CDH12 (chr5, both loci same chrom)
{ $SAM view -b "$PC/chr5.bam" chr5:70800000-71700000 -o "$OUT/cdh12.bam" 2>/dev/null && $SAM index "$OUT/cdh12.bam"; } || die "could not build/index $OUT/cdh12.bam from $PC/chr5.bam"
# AC134878.2(chr22)+TEKT4P2(chr21) cross-chrom
$SAM view -b "$PC/chr22.bam" chr22:5780000-5880000 -o "$OUT/_a.bam" 2>/dev/null || die "could not subset chr22:5780000-5880000"
$SAM view -b "$PC/chr21.bam" chr21:5670000-5770000 -o "$OUT/_b.bam" 2>/dev/null || die "could not subset chr21:5670000-5770000"
{ $SAM merge -f "$OUT/ac13.bam" "$OUT/_a.bam" "$OUT/_b.bam" 2>/dev/null && $SAM index "$OUT/ac13.bam"; } || die "could not merge/index $OUT/ac13.bam"
# AC126603.1(chr15)+AC142384.1(chr16) cross-chrom
$SAM view -b "$PC/chr15.bam" chr15:18100000-18200000 -o "$OUT/_c.bam" 2>/dev/null || die "could not subset chr15:18100000-18200000"
$SAM view -b "$PC/chr16.bam" chr16:34860000-34950000 -o "$OUT/_d.bam" 2>/dev/null || die "could not subset chr16:34860000-34950000"
{ $SAM merge -f "$OUT/ac12.bam" "$OUT/_c.bam" "$OUT/_d.bam" 2>/dev/null && $SAM index "$OUT/ac12.bam"; } || die "could not merge/index $OUT/ac12.bam"
# NCF1/B/C (chr7)
{ $SAM view -b "$PC/chr7.bam" chr7:74000000-76500000 -o "$OUT/ncf1.bam" 2>/dev/null && $SAM index "$OUT/ncf1.bam"; } || die "could not build/index $OUT/ncf1.bam from $PC/chr7.bam"

check(){ # name bam xflag "chrom:s-e:label ..."
  local name=$1 bam=$2 xf=$3; shift 3
  # §6am: an absent/empty region subset would make every locus score OFF=N ON=N — a silent non-result.
  need_out "$bam"
  timeout 600 "$BIN" --bam "$bam" --fasta "$FA" $xf --out "$OUT/${name}_off" >/dev/null 2>&1 || die "$name OFF run failed (exit $?)"
  RUSTLE_MISCHAIN_SALVAGE=1 timeout 600 "$BIN" --bam "$bam" --fasta "$FA" $xf --out "$OUT/${name}_on" >/dev/null 2>&1 || die "$name ON run failed (exit $?)"
  need_out "$OUT/${name}_off.copies.tsv"; need_out "$OUT/${name}_on.copies.tsv"
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
