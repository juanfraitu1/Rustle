#!/usr/bin/env bash
# A/B exhibit: a copy starved of PRIMARY reads (2 < gate 3) — recovered only when the multimapping (AS-tied
# secondary) reads are fed in via --recover-copies.
set -euo pipefail
OUT=/home/juanfra/winloci_scratch
PY=/home/juanfra/miniforge3/bin/python
MM2=/home/juanfra/miniforge3/bin/minimap2
SAM=/home/juanfra/miniforge3/bin/samtools
REL="$(cd "$(dirname "$0")/../target/release" && pwd)"

"$PY" "$(dirname "$0")/sim_starved.py"
"$MM2" -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes "$OUT/simsc.fasta" "$OUT/simsc_reads.fastq" 2>/dev/null \
  | "$SAM" sort -o "$OUT/simsc.bam"; "$SAM" index "$OUT/simsc.bam"

echo; echo "### primary vs secondary read counts per planted copy (minimap2) ###"
"$PY" - <<'PY'
import pysam
from collections import defaultdict
loci=[(f.split("\t")[1],f.split("\t")[2],int(f.split("\t")[3]),int(f.split("\t")[4]),f.split("\t")[6])
      for f in open("/home/juanfra/winloci_scratch/simsc_truth.tsv").read().splitlines()[1:]]
b=pysam.AlignmentFile("/home/juanfra/winloci_scratch/simsc.bam","rb"); P=defaultdict(int); S=defaultdict(int)
for a in b.fetch():
    if a.is_unmapped: continue
    m=(a.reference_start+a.reference_end)//2
    for ci,c,s,e,role in loci:
        if a.reference_name==c and s<=m<=e:
            (S if (a.is_secondary or a.is_supplementary) else P)[ci]+=1; break
for ci,c,s,e,role in loci: print(f"  copy{ci} ({role:16}): primary={P[ci]:3}  secondary={S[ci]:3}")
PY

for mode in base recover; do
  echo; echo "### copy_assign  ${mode}  ###"
  flag=""; [ "$mode" = recover ] && flag="--recover-copies --as-ratio 0.98"
  "$REL/copy_assign" --bam "$OUT/simsc.bam" --fasta "$OUT/simsc.fasta" \
    --region sc:0-70000 --min-copies 2 $flag --out "$OUT/simsc_$mode" 2>/dev/null || true
  echo "-- families.tsv (n_copies / rescued_copies) --"
  if [ -f "$OUT/simsc_$mode.families.tsv" ]; then
    awk -F'\t' 'NR==1||$3>=1{print "   "$1, "n_copies="$3, "rescued="$4}' "$OUT/simsc_$mode.families.tsv" | head
  else echo "   (no family detected)"; fi
done
echo; echo "DONE. base = primary-only; recover = +AS-tied secondaries (the multimapping-read rescue)."
