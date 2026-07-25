#!/bin/bash
# Test whether dna_family_fallback recovers each TARGET family, ONE AT A TIME, writing results INCREMENTALLY
# (survives terminal closure — a killed run only loses the in-flight target). Each target = a small region run
# with --homology-primary + RUSTLE_DNA_FAMILY_FALLBACK=1 (~1-2 min). A member is recovered if a main copy OR a
# dna_family locus (its own coords OR a projection locus) overlaps the member span.
set -u
CACHE=/home/juanfra/winloci_scratch/soto_cache; PC=$CACHE/perchrom
BIN=/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog
FA=/mnt/wsl/PHYSICALDRIVE0p2/home/juanfraitu/winloci_data/chm13v2.0.fa
SAM=/home/juanfra/miniforge3/bin/samtools
OUT=$CACHE/dna_targets; mkdir -p "$OUT"
RES=$OUT/results.tsv
[ -f "$RES" ] || echo -e "family\tmember\tregion\trecovered\tvia\tdetail" > "$RES"

# member recovered? check copies.tsv + dna_family.tsv (own coords + projection loci) for overlap of chrom:s-e
recovered() { # $1=out_prefix $2=chrom $3=s $4=e  -> prints "Y<TAB>via<TAB>detail" or "N\t-\t-"
  local p=$1 c=$2 s=$3 e=$4
  # main copies
  local m; m=$(awk -F'\t' -v c="$c" -v s="$s" -v e="$e" 'NR>1 && $4==c && !($5>e||$6<s){print $5"-"$6;exit}' "$p.copies.tsv" 2>/dev/null)
  [ -n "$m" ] && { echo -e "Y\tRNA-copy\t$m"; return; }
  # dna_family own coords OR projection loci
  if [ -f "$p.dna_family.tsv" ]; then
    awk -F'\t' -v c="$c" -v s="$s" -v e="$e" 'NR>1{
      if($2==c && !($3>e||$4<s)){print "Y\tdna-own\t"$2":"$3"-"$4; found=1; exit}
      n=split($8,ps,";"); for(i=1;i<=n;i++){ split(ps[i],a,"@"); split(a[1],b,":"); split(b[2],r,"-");
        if(b[1]==c && !(r[1]>e||r[2]<s)){print "Y\tdna-proj\t"ps[i]; found=1; exit} }
    } END{if(!found)print "N\t-\t-"}' "$p.dna_family.tsv"
  else echo -e "N\t-\t-"; fi
}

run_target() { # $1=fam $2=member $3=member_chrom:s-e $4=region_bed_spec (space-sep chrom:s-e list to extract)
  local fam=$1 member=$2 mspec=$3; shift 3
  grep -q "^$fam	$member	" "$RES" && { echo "[$fam/$member] already done"; return; }
  local pre="$OUT/${fam}_${member}"
  # build region BAM from per-chrom cache mini-BAMs
  local parts=() ; for reg in "$@"; do local c=${reg%%:*}; $SAM view -b "$PC/${c}.bam" "$reg" -o "$pre.$c.bam" 2>/dev/null; parts+=("$pre.$c.bam"); done
  $SAM merge -f "$pre.bam" "${parts[@]}" 2>/dev/null && $SAM index "$pre.bam"
  RUSTLE_DNA_FAMILY_FALLBACK=1 timeout 600 "$BIN" --bam "$pre.bam" --fasta "$FA" --homology-primary --out "$pre" >"$pre.log" 2>&1
  local mc=${mspec%%:*}; local mr=${mspec#*:}; local ms=${mr%-*}; local me=${mr#*-}
  local rec; rec=$(recovered "$pre" "$mc" "$ms" "$me")
  echo -e "$fam\t$member\t$mspec\t$rec" >> "$RES"
  echo "[$fam/$member] -> $rec"
  rm -f "$pre".*.bam "$pre.bam" "$pre.bam.bai" 2>/dev/null   # keep .copies/.dna_family/.log, drop bulky BAMs
}

# targets: fam  member  member_chrom:s-e  region(s)-to-extract
run_target ID_313 CDH12P1    chr5:70901860-70902156  chr5:70800000-71700000
run_target ID_313 CDH12P3    chr5:71589693-71589989  chr5:70800000-71700000
run_target ID_402 NCF1       chr7:75976253-75991692  chr7:74300000-76500000
run_target ID_402 NCF1B      chr7:74420835-74436153  chr7:74300000-76500000
run_target ID_148 AC126603.1 chr15:18136928-18186041 chr15:18080000-18200000 chr16:34860000-34950000
run_target ID_148 AC142384.1 chr16:34884759-34929924 chr15:18080000-18200000 chr16:34860000-34950000
run_target ID_302 BOLA2B     chr16:29735379-29736755 chr16:29700000-29780000
run_target ID_481 UBE2Q2P8   chr15:82169947-82170050 chr15:82150000-82190000
run_target ID_222 AC243829.6 chr17:37415116-37415223 chr17:37380000-37440000
run_target ID_407 TRIM74     chr7:74187708-74197691  chr7:74100000-76800000
echo "=== DONE. results: ==="; column -t "$RES"
