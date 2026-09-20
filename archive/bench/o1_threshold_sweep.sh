#!/bin/bash
# G3 grid (PREREG adj/g3): one-at-a-time from the defaults, gorilla 3-contig (clusters only) and the Soto slice.
set -u
MF=/mnt/linuxdisk/home/juanfraitu/rustle_target/release/mcl_families
M=/mnt/linuxdisk/home/juanfraitu/mcl_ann; S=/mnt/linuxdisk/home/juanfraitu/soto_mcl; OUT=$M/adj/g3; mkdir -p $OUT
run () { tag=$1; shift
  (cd $M && $MF --paf allgenes.asm20.paf --gff /mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff --min-exonic-bp 1 --no-emit-units "$@" --out $OUT/ggo_$tag > $OUT/ggo_$tag.log 2>&1); echo "ggo $tag exit=$?"
  (cd $S && $MF --paf allgenes.asm20.paf --gff hsa.gff --min-exonic-bp 1 --no-emit-units "$@" --out $OUT/hsa_$tag > $OUT/hsa_$tag.log 2>&1); echo "hsa $tag exit=$?"
}
run default
for v in 0.60 0.65 0.75 0.80 0.85 0.90; do run id$v --min-identity $v; done
for v in 0.10 0.20 0.40 0.50 0.60; do run cov$v --min-cov-longer $v; done
for v in 100 200 500 1000 2000; do run bp$v --min-bp $v; done
echo G3_DONE
