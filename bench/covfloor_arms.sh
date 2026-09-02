#!/bin/bash
# §6bp — the ASYMMETRIC two-sided coverage clause, END TO END through the real binary.
# Env is arm_f2's exactly (jm_arms.sh's OFF arm), so any difference is charged to the code change alone.
# K1 OFF must be BYTE-IDENTICAL to arm_f2 (else the clause is not additive)
# K2 the params certificate must distinguish the arms (defect M2)
# K3 ON at 0.30 must hold NPIP recall at 14/31
set -o pipefail
cd /mnt/linuxdisk/home/juanfraitu/npip_cat
BIN=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/gw_family_catalog
FA=/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta
run () {
  ARM=$1; FLOOR=$2
  unset RUSTLE_FLAGFREE_SITES RUSTLE_TIER2_ADMIT RUSTLE_COLLAPSE_EXONIC RUSTLE_COLLAPSE_UNSTRANDED \
        RUSTLE_JUNCTION_MAJORITY RUSTLE_READ_STRAND RUSTLE_ER_WEIGHTED_PARTITION RUSTLE_FOOTPRINT_NODES \
        RUSTLE_ER_COVERAGE_LONGER RUSTLE_ER_SUM_COVERAGE RUSTLE_ER_CORE_COVERAGE \
        RUSTLE_ER_COVERAGE_LONGER_FLOOR
  export RUSTLE_GATE_MIN_READS=2
  [ -n "$FLOOR" ] && export RUSTLE_ER_COVERAGE_LONGER_FLOOR=$FLOOR
  mkdir -p arm_$ARM/dump
  export RUSTLE_ER_EDGE_DUMP=/mnt/linuxdisk/home/juanfraitu/npip_cat/arm_$ARM/dump/e
  date +"%T start $ARM floor=${FLOOR:-unset}"
  "$BIN" --bam npip3.bam --fasta "$FA" --out arm_$ARM/cat --homology-primary --threads 4 \
    > arm_$ARM/run.log 2> arm_$ARM/run.err
  echo "EXIT=$? arm=$ARM"
  grep -E "^min_coverage_longer" arm_$ARM/dump/e.params.tsv
  wc -l < arm_$ARM/dump/e.edges.tsv
}
run cfoff ""
run cf030 0.30
