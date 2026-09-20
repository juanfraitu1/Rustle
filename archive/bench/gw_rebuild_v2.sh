#!/bin/bash
# Genome-wide gorilla catalog under the L1/L2 defaults (§6fh): dropped members emitted (`member_status`), the
# read-supported locus extent (`locus_start`/`locus_end`) — same inputs as bench/gw_rebuild.sh (gw_units_v1).
# Output: mcl_ann/gw_units_v2.*  Then: bench/o2_sweep_split.py gw_units_v2 sweep_gw_v2 and sweep_gw_v2/run_sweep.sh.
set -u
M=/mnt/linuxdisk/home/juanfraitu/mcl_ann; W=/mnt/linuxdisk/home/juanfraitu/winloci_data
MF=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/mcl_families
cd $M && ( ulimit -v 20000000; /usr/bin/time -f "TIME %e s rss %M kB" $MF --paf allgenes_gw.asm20.paf --gff $W/GGO_genomic.gff --min-exonic-bp 1 \
  --sedef $W/GGO_sedef_final.bed --rmsk $W/rmsk/GCF_029281585.2.repeatMasker.out \
  --bam $W/GGO_ds.bam --fasta /home/juanfra/winloci_scratch/GGO.fasta --out $M/gw_units_v2 ) > $M/adj/gw_units_v2.log 2>&1
echo "gw exit=$?"; tail -3 $M/adj/gw_units_v2.log
