#!/bin/bash
# align_driver.sh: for each autosome, minimap2 -x asm5 -c --cs of the _pri chromosome (10 Mb chunks, query) against its B partner
# (target), via tools/mm2_shard.sh (resumable). Stops with exit 75 at the call's deadline; rerun to resume. Writes out/<chrom>.paf.
set -uo pipefail
D=/mnt/linuxdisk/tmp/rna_allele; H=/mnt/linuxdisk/home/juanfraitu/gorilla_haps
PRI=/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta
SHARD=/mnt/c/Users/jfris/Desktop/Rustle/tools/mm2_shard.sh
export TMPDIR=$D/tmp; mkdir -p $TMPDIR $D/q $D/t $D/out
DEADLINE=$(( $(date +%s) + ${BUDGET:-540} ))
while IFS=$'\t' read pri chrom same src chk bh bname; do
  [ -z "$bname" ] && continue
  [ -s $D/out/chr$chrom.paf.done ] && continue
  if [ ! -s $D/q/chr$chrom.fa ]; then
    python3 - "$PRI" "$pri" $D/q/chr$chrom.fa <<'PY'
import subprocess, sys, os
CH = int(os.environ.get("CHUNK", "10000000"))
fa, name, out = sys.argv[1:4]
seq = "".join(l.strip() for l in subprocess.run(["samtools","faidx",fa,name],capture_output=True,text=True,check=True).stdout.splitlines()[1:])
with open(out + ".tmp", "w") as f:
    for i in range(0, len(seq), CH):
        f.write(f">{name}:{i}\n"); s = seq[i:i+CH]
        for j in range(0, len(s), 80): f.write(s[j:j+80] + "\n")
import os; os.replace(out + ".tmp", out)
PY
  fi
  [ -s $D/t/chr$chrom.fa ] || samtools faidx $H/$bh.fa $bname > $D/t/chr$chrom.fa
  now=$(date +%s); [ $now -ge $DEADLINE ] && { echo "deadline before chr$chrom"; exit 75; }
  MM2_SHARD_DEADLINE=$DEADLINE bash $SHARD paf $D/out/chr$chrom.paf -x asm5 -c --cs -t 4 $D/t/chr$chrom.fa $D/q/chr$chrom.fa 2> $D/out/chr$chrom.shard.log
  rc=$?
  if [ $rc -eq 0 ]; then echo "$(wc -l < $D/out/chr$chrom.paf) records" > $D/out/chr$chrom.paf.done; echo "chr$chrom done $(cat $D/out/chr$chrom.paf.done)";
  else echo "chr$chrom rc=$rc $(grep -c 'done:' $D/out/chr$chrom.shard.log) shards done this call"; exit 75; fi
done < <(tail -n +2 $D/chrmap.tsv)
echo ALL_DONE
