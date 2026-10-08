#!/bin/bash
# resumable: IsoCon per family until the deadline; skips families with final_candidates.fa. H = the work dir (env, required).
# Copy of TRUTH/refabsent/iso_batch.sh with H from the environment.   H=W/isoc BUDGET=540 tools/rlock.sh heavy bash iso_batch.sh
H=${H:?work dir}; ISO=/home/juanfra/miniforge3/envs/isocon/bin/IsoCon
DEADLINE=$(( $(date +%s) + ${BUDGET:-540} ))
for fa in $(ls -S -r $H/fam/*.fa); do
  f=$(basename $fa .fa); out=$H/iso/$f
  [ -s $out/final_candidates.fa ] && continue
  [ -e $out.empty ] && continue
  [ $(grep -c ">" $fa) -lt 2 ] && { mkdir -p $H/iso; touch $out.empty; continue; }
  now=$(date +%s); [ $now -ge $(( DEADLINE - 30 )) ] && { echo "deadline"; exit 75; }
  rm -rf $out; mkdir -p $H/iso
  timeout $(( DEADLINE - now )) $ISO pipeline -fl_reads $fa -outfolder $out --nr_cores 4 > $out.log 2>&1
  rc=$?
  if [ $rc -ne 0 ]; then echo "$f rc=$rc"; [ $rc -eq 124 ] && { rm -rf $out; exit 75; }; touch $out.empty; fi
done
echo ALL_DONE
