#!/bin/bash
# identity_check.sh — the consolidation safety harness (plan 2026-10-05, Phases 0-4 of the .rs merge work).
#   tools/identity_check.sh golden [DIR]   run everything, SAVE the product set as the golden copy
#   tools/identity_check.sh check [DIR]    run everything, cmp every product against the golden copy
#   DIR defaults to /mnt/linuxdisk/home/juanfraitu/identity_golden
# What it proves after any merge/refactor: (1) the full release test suite; (2) an end-to-end mcl_families
# --from-gtf --emit-units --emit-relations run on a synthetic fixture (exercises the MCL cluster: fam_from_gtf,
# annotation_families, family_container, family_relations, run_cache PAF caching OFF); (3) with FULL=1 also a
# real-data slice: copy_assign --assemble-only + mcl_families on a 5-Mb gorilla BAM slice (the Phase-3
# hash-swap order-leak detector — minimap2 runs for real, ~6 min, heavy lock).
set -uo pipefail
MODE=${1:?golden|check}
GOLD=${2:-/mnt/linuxdisk/home/juanfraitu/identity_golden}
REPO=$(cd "$(dirname "$0")/.." && pwd)
cd "$REPO"
export PATH="$HOME/miniforge3/bin:$PATH"
export CARGO_TARGET_DIR=${CARGO_TARGET_DIR:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2}
BIN=$CARGO_TARGET_DIR/release
mkdir -p "$GOLD"

echo "[identity] 1/3 release test suite"
bash tools/rlock.sh heavy cargo test --release > /tmp/identity_test.log 2>&1 || { echo "[identity] TEST SUITE FAILED"; tail -5 /tmp/identity_test.log; exit 1; }
passed=$(grep -E "test result" /tmp/identity_test.log | awk -F'[.;] ' '{s+=$2} END{print s}')
echo "[identity] suite: $passed tests passed"

echo "[identity] 2/3 synthetic MCL end-to-end"
W=/tmp/identity_syn; rm -rf $W; mkdir -p $W
python3 - <<'EOF'
import random
rnd = random.Random(7)
def seq(n): return ''.join(rnd.choice("ACGT") for _ in range(n))
# two loci on c1 at ~97% identity (one edge), a third diverged locus on c2
a = seq(3000)
def mut(s, r):
    l = list(s)
    for i in range(len(l)):
        if rnd.random() < r: l[i] = rnd.choice("ACGT")
    return ''.join(l)
b = mut(a, 0.03)
c = seq(2500)
with open("/tmp/identity_syn/genome.fa", "w") as f:
    f.write(f">c1\n{a}{b}\n>c2\n{c}\n")
def gtf(gene, chrom, off, ln, txs):
    out = []
    ex = [(off+1, off+300), (off+501, off+ln)]
    for i, (tid, reads) in enumerate(txs):
        (s, e) = (ex[0][0], ex[1][1])
        out.append(f'{chrom}\trustle\ttranscript\t{s}\t{e}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{tid}"; reads "{reads}";')
        for (x, y) in ex:
            out.append(f'{chrom}\trustle\texon\t{x}\t{y}\t.\t+\t.\tgene_id "{gene}"; transcript_id "{tid}";')
    return "\n".join(out)
rows = [gtf("G1", "c1", 0, 1200, [("T1", 30), ("T2", 12)]),
        gtf("G2", "c1", 3000, 1100, [("T3", 25)]),
        gtf("G3", "c2", 0, 900, [("T4", 20)])]
open("/tmp/identity_syn/in.gtf", "w").write("\n".join(rows) + "\n")
EOF
( cd $W && env -u RUSTLE_CACHE_DIR "$BIN/mcl_families" --from-gtf in.gtf --fasta genome.fa --threads 2 \
    --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --emit-relations --out out.fam > mcl.log 2>&1 ) || { echo "[identity] MCL RUN FAILED"; cat $W/mcl.log; exit 1; }
if [ "$MODE" = golden ]; then
  cp $W/out.fam.copies.tsv $W/out.fam.copies.fa $W/out.fam.clusters.tsv $W/out.fam.params.tsv $W/out.fam.relations.tsv $W/out.fam.members_by_locus.tsv $W/out.fam.loci.gff3 $GOLD/
  echo "[identity] golden saved to $GOLD"
else
  for f in out.fam.copies.tsv out.fam.copies.fa out.fam.clusters.tsv out.fam.params.tsv out.fam.relations.tsv out.fam.members_by_locus.tsv out.fam.loci.gff3; do
    cmp -s $W/$f "$GOLD/$f" && echo "[identity] IDENTICAL $f" || { echo "[identity] DIFFERENT $f"; exit 1; }
  done
fi

if [ "${FULL:-0}" = 1 ]; then
  echo "[identity] 3/3 real-data slice (assemble + families on a 5-Mb gorilla slice)"
  BAM=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_mm.bam
  FASTA=/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta
  R=$W/real; mkdir -p $R
  if [ "$MODE" = golden ]; then
    bash tools/rlock.sh heavy bash -c "samtools view -b $BAM NC_073244.2:1-5000000 > $R/slice.bam && samtools index $R/slice.bam"
    bash tools/rlock.sh heavy bash -c "cd $R && '$BIN/copy_assign' --assemble-only --genome-wide --assembly-junctions strict --bam slice.bam --fasta '$FASTA' --out slice --threads 4 > assemble.log 2>&1"
    bash tools/rlock.sh heavy bash -c "cd $R && env -u RUSTLE_CACHE_DIR '$BIN/mcl_families' --from-gtf slice.families.gtf --fasta '$FASTA' --threads 4 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --out slice.fam > families.log 2>&1"
    cp $R/slice.gtf $R/slice.fam.copies.tsv $R/slice.fam.clusters.tsv $GOLD/real_
    echo "[identity] real-data golden saved"
  else
    bash tools/rlock.sh heavy bash -c "samtools view -b $BAM NC_073244.2:1-5000000 > $R/slice.bam && samtools index $R/slice.bam"
    bash tools/rlock.sh heavy bash -c "cd $R && '$BIN/copy_assign' --assemble-only --genome-wide --assembly-junctions strict --bam slice.bam --fasta '$FASTA' --out slice --threads 4 > assemble.log 2>&1"
    bash tools/rlock.sh heavy bash -c "cd $R && env -u RUSTLE_CACHE_DIR '$BIN/mcl_families' --from-gtf slice.families.gtf --fasta '$FASTA' --threads 4 --min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units --out slice.fam > families.log 2>&1"
    for f in slice.gtf slice.fam.copies.tsv slice.fam.clusters.tsv; do
      cmp -s $R/$f "$GOLD/real_$f" && echo "[identity] IDENTICAL real $f" || { echo "[identity] DIFFERENT real $f"; exit 1; }
    done
  fi
fi
echo "[identity] ALL GREEN ($MODE)"
