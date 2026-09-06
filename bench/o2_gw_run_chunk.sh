#!/bin/bash
# Foreground chunk driver: runs not-yet-done families for at most BUDGET seconds of launches (each family capped
# at PER_FAM seconds; a capped family records exit=124 and can be rerun singly with a longer cap), then exits.
B=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign
BAM=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_ds.bam; FA=/home/juanfra/winloci_scratch/GGO.fasta
BUDGET=${BUDGET:-300}; PER_FAM=${PER_FAM:-240}; t0=$(date +%s); n=0
for d in $(ls -d fam_* | sort); do
  [ -s $d/A.done ] && continue
  [ $(( $(date +%s) - t0 )) -ge $BUDGET ] && break
  ( ulimit -v 12000000; /usr/bin/time -f "TIME elapsed=%e rss_kb=%M" timeout $PER_FAM $B --bam $BAM --fasta $FA --families $d/copies.tsv --copies-fa $d/copies.fa --regions $d/regions --out $d/A > $d/A.log 2> $d/A.err )
  echo "exit=$? copies=$(($(wc -l < $d/copies.tsv)-1))" > $d/A.done; n=$((n+1))
done
echo "chunk: $n families in $(( $(date +%s) - t0 )) s; done $(ls */A.done | wc -l)/$(ls -d fam_* | wc -l); non-zero exits $(grep -L 'exit=0' */A.done | wc -l)"
