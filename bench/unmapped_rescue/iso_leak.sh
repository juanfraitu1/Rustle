#!/bin/bash
# resumable: IsoCon per locus until the deadline; skips loci with final_candidates.fa or a .failed marker
H=/mnt/linuxdisk/tmp/o3_rescue/mattruth/isocon_leak; ISO=/home/juanfra/miniforge3/envs/isocon/bin/IsoCon
DEADLINE=$(( $(date +%s) + ${BUDGET:-540} ))
mkdir -p $H/iso
for fa in $(ls -S -r $H/reads/*.fa); do
  f=$(basename $fa .fa); out=$H/iso/$f
  [ -s $out/final_candidates.fa ] && continue
  [ -e $out.failed ] && continue
  now=$(date +%s); [ $now -ge $(( DEADLINE - 30 )) ] && { echo "deadline"; exit 75; }
  rm -rf $out
  timeout $(( DEADLINE - now )) $ISO pipeline -fl_reads $fa -outfolder $out --nr_cores 4 > $out.log 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then [ $rc -eq 124 ] && { rm -rf $out; echo "timeout on $f"; exit 75; }; echo "$f rc=$rc"; touch $out.failed; fi
done
echo ALL_DONE
