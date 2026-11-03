#!/usr/bin/env bash
# Single-exon models whose span has 1-kb bins with no PRIMARY read depth (r1518). `samtools depth` keeps its default flag
# filter (UNMAP,SECONDARY,QCFAIL,DUP), so "zero depth" = no primary record over the bin.
# usage: mono_gaps.sh GTF BAM OUT.tsv [MIN_SPAN=20000]   -> OUT.tsv: tid reads n_bins zero_bins max_zero_run ; summary on stdout
set -eu
G=$1; BAM=$2; OUT=$3; MIN=${4:-20000}
python3 - "$G" "$MIN" <<'PY' > "$OUT.regions"
import re,sys
nex={}; tr={}
for l in open(sys.argv[1]):
    f=l.rstrip('\n').split('\t')
    if len(f)<9: continue
    t=re.search(r'transcript_id "([^"]+)"',f[8]).group(1)
    if f[2]=='exon': nex[t]=nex.get(t,0)+1
    elif f[2]=='transcript': tr[t]=(f[0],int(f[3]),int(f[4]),int((re.search(r'reads "(\d+)"',f[8]) or [0,0])[1]))
for t,(c,s,e,r) in tr.items():
    if nex.get(t)==1 and e-s+1>=int(sys.argv[2]): print(f"{c}\t{s-1}\t{e}\t{t}\t{r}")
PY
: > "$OUT"
while IFS=$'\t' read c s e t r; do
  samtools depth -a -r "$c:$((s+1))-$e" "$BAM" | awk -v t="$t" -v r="$r" -v s="$s" -v L=$((e-s)) '{b=int(($2-1-s)/1000); d[b]+=$3; n[b]++} END{z=0;mx=0;run=0;nb=0; for(b=0;b<=int((L-1)/1000);b++){nb++; if(n[b]==0 || d[b]/n[b]<1){z++; run++; if(run>mx)mx=run} else run=0}; printf "%s\t%d\t%d\t%d\t%d\n", t, r, nb, z, mx}' >> "$OUT"
done < "$OUT.regions"
awk -F'\t' '{n++; if($4>0) g++; if($5>=3) g3++; if($5>=10) g10++} END{printf "models %d, with a zero-primary 1-kb bin %d, zero-run >= 3 kb %d, >= 10 kb %d\n", n, g, g3, g10}' "$OUT"
