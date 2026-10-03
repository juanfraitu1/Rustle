#!/bin/bash
# accept_o3_candidates.sh — Amendment 12 (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md): `o3_candidates` on Amendment 7's 53-family
# held-out, scored by `merge_test.py score` (arm M = masked genome + one union per flagged candidate, each candidate its own component).
#
#   accept_o3_candidates.sh <step> [run] [arg]    run = reg (delta 0.00958, the registered run, work dir $A) | half (0.00479, $A/half) |
#                                                 double (0.01916, $A/double)
#
# Steps, in order (every heavy step is ONE foreground call under tools/rlock.sh heavy, each < 10 min; light ones under rlock light):
#   copies               panel.json -> A12.copies.tsv / A12.copies.fa / A12.regions (panel_to_copies.py)              light
#   nets                 the stage's read nets recomputed (pass A / pass B / cap) -> nets.tsv, net_reads.tsv          light
#   plan [n] [k] [file]  families not in the first k lines of the batch file -> n new batches on an estimated cost     light
#   split                each scored part in two halves -> $A/parts/scored.part{0,1,2}{a,b}.fa (shorter malign calls)   light
#   stage <run> <g>      o3_candidates on batch g (line g+1 of $W/batches.txt, else $A/batches.txt), time -v -> cand_g<g>.* heavy
#   concat <run>         the batches' products -> cand.* (one header; families in the stage's order)                 light
#   contigs <run>        flagged unions renamed iso_<family>_<k> -> iso.contigs.fa, iso_names.tsv                    light
#   mindex <run>         M.fa = masked.fa + iso.contigs.fa; minimap2 -x splice -d M.splice.mmi (~5 min, ~21 GB)      heavy
#   malign <run> <p>     scored.part<p>.fa (p = 0..2, or a half 0a..2b) realigned to M (R / RIL arms' flags + -K 100M)  heavy
#   mmerge <run>         M.{0,1,2}.bam or M.{0,1,2}{a,b}.bam -> RIL.bam (+ .bai)                                      heavy
#   label <run>          iso.contigs.fa vs the UNMASKED genome (minimap2 -c -x splice:hq -uf -N 20) -> iso.base.paf;
#                        contigs.tsv (source D / S:<copy> / elsewhere / none), empty merge/paf/*.paf, links            heavy + light
#   score <run>          merge_test.py components (singletons: no pair alignments) + merge_test.py score -> score.out  light + heavy
#   keep <run>           A12-2: candidates' reads vs unions and cluster consensus sequences -> keep.out, keep.json     heavy
#   report <run>         stage counts, labels, detection, floors, causes -> report.out, report.json                  light
#   decompose <run>      post hoc (no registered rule): arm M read by read vs IsoCon's Amendment 8 arm M, D reads split by
#                        where the R arm put them -> decompose.out, decompose.json                                     light
#   clean <run>          rm M.splice.mmi and M.fa (the BAMs and tables stay)
#
# As run on 2026-10-02: copies; nets; batches.txt = batch 0 written by hand, then `plan 4 1` (calibrated on it); stage reg 0..4; concat,
# contigs, label, mindex reg; malign reg 0..2; mmerge reg; score reg; keep reg; report reg; decompose reg; clean reg. Reruns (the heavy
# lock was shared with another session, so shorter calls): `plan 9 0 $A/half/batches.txt`, copied to $A/double/; split; stage half|double
# 0..8; concat, contigs, label, mindex; malign half 0a..2b / double 0..2; mmerge; score; report; clean.
#
# The stage never sees RUSTLE_CACHE_DIR (every run is computed, so A12-3's wall times are real). Batching does not change any result:
# nets, clusters and candidates are per family, the unmapped-read attribution indexes every family's copies in every batch, and the
# genome hits are per consensus; only nets.fa's de-duplication across families (R9, not used here) depends on the batch.
set -euo pipefail
REPO=${REPO:-$(cd "$(dirname "$0")/../.." && pwd)}   # this checkout unless the caller sets REPO
L=/mnt/linuxdisk/tmp/rna_allele/linktest
A=/mnt/linuxdisk/tmp/rna_allele/a12
BIN=${BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release/o3_candidates}   # the binary of the 2026-10-02 run unless set
GGO_MMI=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi
PY="python3 $REPO/bench/rna_allele/accept_o3_candidates.py"
HEAVY="bash $REPO/tools/rlock.sh heavy"
LIGHT="bash $REPO/tools/rlock.sh light"

dir_of() { case $1 in reg) echo $A ;; half) echo $A/half ;; double) echo $A/double ;; *) echo "run must be reg|half|double" >&2; exit 2 ;; esac; }
delta_args() { case $1 in reg) echo "" ;; half) echo "--delta 0.00479" ;; double) echo "--delta 0.01916" ;; esac; }

step=${1:?step}; shift
case $step in copies|nets|plan|split) run=reg ;; *) run=${1:?run}; shift ;; esac
W=$(dir_of "$run"); mkdir -p "$W/logs"
case $step in
  copies) $LIGHT python3 $REPO/bench/rna_allele/panel_to_copies.py --panel $L/panel.json --bam $L/R.bam --fasta $L/masked.fa --out $A/A12 ;;
  nets)   $LIGHT $PY nets --w $A --linktest $L ;;
  plan)   $LIGHT $PY plan --w $A --batches "${1:-7}" --keep "${2:-0}" --batch-file "${3:-$A/batches.txt}" ;;
  split)  # each scored part in two halves (alternate records): shorter heavy calls when the lock is contended; per-read results unchanged
    mkdir -p $A/parts
    for p in 0 1 2; do
      $LIGHT gawk -v a=$A/parts/scored.part${p}a.fa -v b=$A/parts/scored.part${p}b.fa \
        '/^>/ { n++ } { if (n % 2) print > a; else print > b }' $L/scored.part$p.fa
    done
    grep -c ">" $A/parts/*.fa ;;
  stage)
    BF=$W/batches.txt; [ -s $BF ] || BF=$A/batches.txt          # a run may carry its own batch plan (same families, other grouping)
    g=${1:?batch}; fams=$(sed -n "$((g + 1))p" $BF); [ -n "$fams" ] || { echo "no batch $g in $BF" >&2; exit 2; }
    # stderr (the stage's log and time -v) timestamped line by line (epoch s) for the per-family timing
    env -u RUSTLE_CACHE_DIR $HEAVY /usr/bin/time -v $BIN --bam $L/R.bam --fasta $L/masked.fa --copies $A/A12.copies.tsv --copies-fa $A/A12.copies.fa \
      --index $L/masked.splice.mmi --out $W/cand_g$g --threads 4 $(delta_args "$run") --families "$fams" 2>&1 >/dev/null \
      | gawk '{ print systime() "\t" $0; fflush() }' > $W/logs/stage_g$g.log
    grep -E "Elapsed \(wall|Maximum resident|done:" $W/logs/stage_g$g.log ;;
  concat)  $LIGHT $PY concat --w $W --prefix cand ;;
  contigs) $LIGHT $PY contigs --w $W --prefix cand ;;
  mindex)
    $HEAVY /usr/bin/time -v bash -c "cat $L/masked.fa $W/iso.contigs.fa > $W/M.fa && minimap2 -x splice -t 4 -d $W/M.splice.mmi $W/M.fa" > $W/logs/mindex.log 2>&1
    grep -E "Real time|Elapsed \(wall|Maximum resident" $W/logs/mindex.log ;;
  malign)
    p=${1:?part}; Q=$L/scored.part$p.fa; [ -s $A/parts/scored.part$p.fa ] && Q=$A/parts/scored.part$p.fa     # 0..2, or a half 0a..2b
    $HEAVY /usr/bin/time -v bash -c "set -o pipefail; minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -K 100M -t 4 $W/M.splice.mmi \
      $Q 2> $W/logs/malign$p.mm2.log | samtools sort -@ 1 -m 500M -o $W/M.$p.bam -" > $W/logs/malign$p.log 2>&1
    grep -E "Real time|Elapsed \(wall|Maximum resident" $W/logs/malign$p.mm2.log $W/logs/malign$p.log ;;
  mmerge) # the three parts, or the six halves; never both
    if [ -s $W/M.0.bam ]; then parts="$W/M.0.bam $W/M.1.bam $W/M.2.bam"; else parts=$(ls $W/M.[012][ab].bam | tr '\n' ' '); fi
    [ $(echo $parts | wc -w) -eq 3 ] || [ $(echo $parts | wc -w) -eq 6 ] || { echo "expected 3 parts or 6 halves: $parts" >&2; exit 2; }
    $HEAVY bash -c "samtools merge -f -@ 2 $W/RIL.bam $parts && samtools index $W/RIL.bam" ;;
  label)
    $HEAVY /usr/bin/time -v bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $GGO_MMI $W/iso.contigs.fa > $W/iso.base.paf 2> $W/logs/label.mm2.log" > $W/logs/label.log 2>&1
    $LIGHT $PY label --w $W --linktest $L ;;
  score)
    $LIGHT python3 $REPO/bench/rna_allele/merge_test.py components --w $W > $W/logs/components.out
    $HEAVY python3 $REPO/bench/rna_allele/merge_test.py score --w $W > $W/score.out 2>&1
    cat $W/score.out ;;
  keep)   $HEAVY $PY keep --w $W --linktest $L --prefix cand > $W/keep.out 2> $W/logs/keep.log; cat $W/keep.out ;;
  report) $LIGHT $PY report --w $W --linktest $L --prefix cand --nets $A > $W/report.out; cat $W/report.out ;;
  decompose) $LIGHT $PY decompose --w $W --linktest $L > $W/decompose.out; cat $W/decompose.out ;;
  clean)  rm -f $W/M.splice.mmi $W/M.fa ;;
  *) echo "unknown step $step" >&2; exit 2 ;;
esac
