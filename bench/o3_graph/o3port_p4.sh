#!/bin/bash
# o3port_p4.sh — P4, the two-tier O3 output on a fresh gorilla fibroblast held-out (docs/PREREG_o3_two_tier_p4_2026-10-23.md).
# The registered harness's commands (bench/rna_allele/accept_o3_candidates.sh, ACC=a13 for the deletion arm, ACC=a15h for the
# control arm) with P4's paths; that harness stays as registered. Every heavy step is ONE foreground call under rlock heavy.
#   o3port_p4.sh <step> [arg]
# Substrate: panel, mask, reads, split, index (masked splice index), ralign <p> / rmerge (R.bam: scored reads on the masked genome),
#   r0align <p> / r0merge (R0.bam: on the unmasked `_pri`), copies, batches
# Arms (d = deletion, c = control; k = the stage, p = the port): stage <d|c> <g>, concat <d|c>, nets <d|c>, port <d|c>
# Deletion scoring (X = kd | pd): contigs X, mindex X, malign X <p>, mmerge X, label X, score X, clean X
# Control scoring (X = kc | pc): contigs X, hap X, lift X, clabel X, classify X, cindex X, calign X <p>, cmerge X, cscore X, clean X
# verdict: H4-1 / H4-2 / H4-3 (o3port_p4.py verdict)
set -euo pipefail
REPO=${REPO:-$(cd "$(dirname "$0")/../.." && pwd)}
W=${W:-/mnt/linuxdisk/tmp/rna_allele/p4}
BIN=${BIN:-/mnt/linuxdisk/tmp/rna_allele/p4/bin/candidate_copies}   # frozen copy, sha1 2e3cd592 (both preregs)
GGO_FA=/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta
GGO_MMI=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi
MAT_MMI=/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.mat.splice.mmi
PAT_MMI=/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.pat.splice.mmi
HEAVY="bash $REPO/tools/rlock.sh heavy"
LIGHT="bash $REPO/tools/rlock.sh light"
ACC_PY="python3 $REPO/bench/rna_allele/accept_o3_candidates.py"
CT="python3 $REPO/bench/rna_allele/control_test.py"
MT="python3 $REPO/bench/rna_allele/merge_test.py"
LT="python3 $REPO/bench/rna_allele/link_test.py"
P4="python3 $REPO/bench/o3_graph/o3port_p4.py"
MM_R="-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4"      # Amendment 7 / 9's alignment command

arm_in() {   # the stage's / port's inputs of arm d (masked) or c (unmasked)
  case $1 in
    d) echo "--bam $W/R.bam --fasta $W/masked.fa --copies $W/d.copies.tsv --copies-fa $W/d.copies.fa --index $W/masked.splice.mmi" ;;
    c) echo "--bam $W/R0.bam --fasta $GGO_FA --copies $W/c.copies.tsv --copies-fa $W/c.copies.fa --index $GGO_MMI" ;;
    *) echo "arm must be d or c" >&2; exit 2 ;;
  esac
}
mkdir -p $W/logs
step=${1:?step}; shift
case $step in
  panel)   $LIGHT $P4 panel --w $W ;;
  subpanel) $LIGHT $P4 subpanel --w $W --n ${1:-20} ;;   # P5: the 20 of P4's draw and their reads
  mask)    $HEAVY $LT mask --w $W ;;
  reads)   $HEAVY $LT reads --w $W ;;
  split)   $LIGHT $P4 split --w $W ;;
  index)   $HEAVY /usr/bin/time -v minimap2 -x splice -t 4 -d $W/masked.splice.mmi $W/masked.fa > $W/logs/index.log 2>&1
           grep -E "Elapsed \(wall|Maximum resident" $W/logs/index.log ;;
  ralign|r0align)
    p=${1:?part}
    if [ $step = ralign ]; then IDX=$W/masked.splice.mmi; OUT=$W/R.$p.bam; else IDX=$GGO_MMI; OUT=$W/R0.$p.bam; fi
    $HEAVY /usr/bin/time -v bash -c "set -o pipefail; minimap2 $MM_R $IDX $W/scored.part$p.fa 2> $W/logs/$step$p.mm2.log \
      | samtools sort -@ 1 -m 500M -o $OUT -" > $W/logs/$step$p.log 2>&1
    grep -E "Elapsed \(wall|Maximum resident" $W/logs/$step$p.log ;;
  rmerge)  $HEAVY bash -c "samtools merge -f -@ 2 $W/R.bam $W/R.0.bam $W/R.1.bam $W/R.2.bam && samtools index $W/R.bam" ;;
  r0merge) $HEAVY bash -c "samtools merge -f -@ 2 $W/R0.bam $W/R0.0.bam $W/R0.1.bam $W/R0.2.bam && samtools index $W/R0.bam" ;;
  copies)  # the deletion arm's copies are the survivors on the masked genome (A12's way); the control's are every copy (A14's --all)
    $LIGHT python3 $REPO/bench/rna_allele/panel_to_copies.py --panel $W/panel.json --bam $W/R.bam --fasta $W/masked.fa --out $W/d
    $LIGHT python3 $REPO/bench/rna_allele/panel_to_copies.py --all --panel $W/panel.json --bam $W/R0.bam --fasta $GGO_FA --out $W/c ;;
  batches) $LIGHT $P4 batches --w $W ;;
  stage)   # the harness's stage step: stderr timestamped, the binary's sha1 first (accept_o3_candidates.py nets reads these logs)
    arm=$1; g=${2:?batch}; X=$W/k$arm; mkdir -p $X/logs
    fams=$(sed -n "$((g + 1))p" $W/batches.txt); [ -n "$fams" ] || { echo "no batch $g" >&2; exit 2; }
    printf '%s\tbinary sha1 %s\n' "$(date +%s)" "$(sha1sum $BIN | cut -c1-40)" > $X/logs/stage_g$g.log
    env -u RUSTLE_CACHE_DIR $HEAVY /usr/bin/time -v $BIN $(arm_in $arm) --out $X/cand_g$g --threads 4 --families "$fams" 2>&1 >/dev/null \
      | gawk '{ print systime() "\t" $0; fflush() }' >> $X/logs/stage_g$g.log
    grep -E "Elapsed \(wall|Maximum resident|pass B:|done:" $X/logs/stage_g$g.log ;;
  concat)  $LIGHT $ACC_PY concat --w $W/k$1 --prefix cand ;;
  nets)    arm=$1; X=$W/k$arm
           B=""; [ $arm = c ] && B="--bam $W/R0.bam"
           [ $arm = d ] && B="--bam $W/R.bam"
           $HEAVY $ACC_PY nets --w $X --linktest $W --prefix cand --copies $W/$arm.copies.tsv --copies-fa $W/$arm.copies.fa $B > $X/attrib.out
           cat $X/attrib.out ;;
  port)    # the port on the stage's net; its IsoCon side = the stage's products. arg 2 = cc: correlation clustering (dir p<arm>cc)
    arm=$1; CL=${2:-}; X=$W/p$arm$CL; mkdir -p $X/logs
    CLA=""; [ "$CL" = cc ] && CLA="--clusterer cc"
    if [ $arm = d ]; then BAMG="--bam $W/R.bam --genome $W/masked.fa"; else BAMG="--bam $W/R0.bam --genome $GGO_FA"; fi
    for f in net_reads.tsv nets.tsv attrib.json attrib.out; do cp $W/k$arm/$f $X/$f; done
    $HEAVY /usr/bin/time -v python3 $REPO/bench/o3_graph/o3port.py run $BAMG --copies $W/$arm.copies.tsv --net $X/net_reads.tsv \
      --out $X/cand --isocon $W/k$arm/cand --dump $X/cand.dump.jsonl $CLA > $X/logs/port.log 2>&1
    grep -E "Elapsed \(wall|Maximum resident" $X/logs/port.log; echo "flags $(($(wc -l < $X/cand.candidates.tsv) - 1))" ;;
  contigs) $LIGHT $ACC_PY contigs --w $W/$1 --prefix cand ;;
  mindex)  X=$W/$1
    $HEAVY /usr/bin/time -v bash -c "cat $W/masked.fa $X/iso.contigs.fa > $X/M.fa && minimap2 -x splice -t 4 -d $X/M.splice.mmi $X/M.fa" > $X/logs/mindex.log 2>&1
    grep -E "Elapsed \(wall|Maximum resident" $X/logs/mindex.log ;;
  malign)  X=$W/$1; p=${2:?part}
    $HEAVY /usr/bin/time -v bash -c "set -o pipefail; minimap2 $MM_R -K 100M $X/M.splice.mmi $W/scored.part$p.fa 2> $X/logs/malign$p.mm2.log \
      | samtools sort -@ 1 -m 500M -o $X/M.$p.bam -" > $X/logs/malign$p.log 2>&1
    grep -E "Elapsed \(wall|Maximum resident" $X/logs/malign$p.log ;;
  mmerge)  X=$W/$1; $HEAVY bash -c "samtools merge -f -@ 2 $X/RIL.bam $X/M.0.bam $X/M.1.bam $X/M.2.bam && samtools index $X/RIL.bam" ;;
  label|clabel)  # contigs vs the unmasked genome -> contigs.tsv sources; clabel: every copy is in the reference (--all-copies)
    X=$W/$1; ALLO=""; [ $step = clabel ] && ALLO=--all-copies
    $HEAVY bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $GGO_MMI $X/iso.contigs.fa > $X/iso.base.paf 2> $X/logs/label.mm2.log"
    $LIGHT $ACC_PY label --w $X --linktest $W $ALLO ;;
  score)   X=$W/$1
    $LIGHT $MT components --w $X > $X/logs/components.out
    $HEAVY $MT score --w $X > $X/score.out 2>&1; grep -E "^\[M\]|^M1|^OVERALL" $X/score.out || tail -5 $X/score.out ;;
  hap)     X=$W/$1
    $HEAVY bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $MAT_MMI $X/iso.contigs.fa > $X/iso.mat.paf 2> $X/logs/hap.mat.log && \
      minimap2 -c -x splice:hq -uf -N 20 -t 4 $PAT_MMI $X/iso.contigs.fa > $X/iso.pat.paf 2> $X/logs/hap.pat.log"
    ln -sfn iso.mat.paf $X/contigs_L.mat.paf; ln -sfn iso.pat.paf $X/contigs_L.pat.paf; wc -l $X/iso.mat.paf $X/iso.pat.paf ;;
  lift)    $LIGHT $CT lift --w $W/$1 --l $W ;;
  classify) X=$W/$1; $LIGHT $CT classify --w $X --l $W > $X/classify.out; tail -6 $X/classify.out ;;
  cindex)  X=$W/$1
    $HEAVY /usr/bin/time -v bash -c "cat $GGO_FA $X/iso.contigs.fa > $X/C.fa && minimap2 -x splice -t 4 -d $X/C.splice.mmi $X/C.fa" > $X/logs/cindex.log 2>&1
    grep -E "Elapsed \(wall|Maximum resident" $X/logs/cindex.log ;;
  calign)  X=$W/$1; p=${2:?part}
    $HEAVY /usr/bin/time -v bash -c "set -o pipefail; minimap2 $MM_R $X/C.splice.mmi $W/scored.part$p.fa 2> $X/logs/calign$p.mm2.log \
      | samtools sort -@ 1 -m 500M -o $X/C.$p.bam -" > $X/logs/calign$p.log 2>&1
    grep -E "Elapsed \(wall|Maximum resident" $X/logs/calign$p.log ;;
  cmerge)  X=$W/$1; $HEAVY bash -c "samtools merge -f -@ 2 $X/C.bam $X/C.0.bam $X/C.1.bam $X/C.2.bam && samtools index $X/C.bam" ;;
  cscore)  X=$W/$1; for f in R0.bam R0.bam.bai; do ln -sfn $W/$f $X/$f; done
    $HEAVY $CT score --w $X --l $W > $X/score.out; tail -8 $X/score.out ;;
  clean)   rm -f $W/$1/M.splice.mmi $W/$1/M.fa $W/$1/C.splice.mmi $W/$1/C.fa ;;
  verdict) $LIGHT $P4 verdict --w $W ;;
  verdict5) $LIGHT $P4 verdict5 --w $W ;;
  *) echo "unknown step $step" >&2; exit 2 ;;
esac
