#!/bin/bash
# accept_o3_candidates.sh — the `o3_candidates` acceptance on Amendment 7's 53-family held-out (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md),
# scored by `merge_test.py score` (arm M = masked genome + one union per flagged candidate, each candidate its own component).
#   ACC=a13 (default): Amendment 13 (+ 13b-13e) — the A13 stage (alignment attribution, structural template), work dir a13/; the copies table,
#                      FASTA, regions and batch plan of A12 are reused (`link`); A13-1 is decided on the re-registered comparator C (`comparator`).
#   ACC=a12:           Amendment 12, as run on 2026-10-02 into a12/ (its k-mer net recomputation `nets` and binary a3564999 are at commit fde90c0a;
#                      the current stage has no k-mer rule, so only its scoring steps re-run here).
#
#   [ACC=a12|a13] accept_o3_candidates.sh <step> [run] [arg]    run = reg (delta 0.00958, the registered run, work dir $A) | half (0.00479,
#                                                               $A/half) | double (0.01916, $A/double)
#
# Steps, in order (every heavy step is ONE foreground call under tools/rlock.sh heavy, each < 10 min; light ones under rlock light):
#   copies               (a12) panel.json -> A12.copies.tsv / A12.copies.fa / A12.regions (panel_to_copies.py)         light
#   link                 (a13) A12's copies table / FASTA / regions / batches.txt and R.bam (+ .bai) linked into $A           light
#   plan [n] [k] [file]  (a13) families not in the first k lines of the batch file -> n new batches on an estimated cost;
#                        refuses to write through a link (a13/batches.txt links to A12's registered plan)               light
#   split                each scored part in two halves -> $A/parts/scored.part{0,1,2}{a,b}.fa (shorter malign calls)   light
#   stage <run> <g>      (a13) o3_candidates on batch g (line g+1 of $W/batches.txt, else $A/batches.txt), time -v ->
#                        cand_g<g>.* (ACC=a12 refused: A12's registered products stay as they are)                    heavy
#   concat <run>         the batches' products -> cand.* (one header; families in the stage's order)                 light
#   nets <run>           (a13) the nets as the stage built them: pass A from R.bam and pass B's attribution re-run per batch
#                        (minimap2 map-ont), checked against the stage's own products (every batch log's pass-B counts,
#                        families.tsv n_net / n_used, reads.tsv, nets.fa) -> nets.tsv, net_reads.tsv, attrib.json, attrib.out   heavy
#   contigs <run>        flagged unions renamed iso_<family>_<k> -> iso.contigs.fa, iso_names.tsv                    light
#   mindex <run>         M.fa = masked.fa + iso.contigs.fa; minimap2 -x splice -d M.splice.mmi (~5 min, ~21 GB)      heavy
#   malign <run> <p>     scored.part<p>.fa (p = 0..2, or a half 0a..2b) realigned to M (R / RIL arms' flags + -K 100M)  heavy
#   mmerge <run>         M.{0,1,2}.bam or M.{0,1,2}{a,b}.bam -> RIL.bam (+ .bai)                                      heavy
#   label <run>          iso.contigs.fa vs the UNMASKED genome (minimap2 -c -x splice:hq -uf -N 20) -> iso.base.paf;
#                        contigs.tsv (source D / S:<copy> / elsewhere / none), empty merge/paf/*.paf, links            heavy + light
#   score <run>          merge_test.py components (singletons: no pair alignments) + merge_test.py score -> score.out  light + heavy
#   keep <run>           A12-2 / A13-2: candidates' reads vs unions and cluster consensus sequences -> keep.out       heavy
#   report <run>         stage counts, labels, detection, floors, causes -> report.out, report.json                  light
#   comparator <run>     (a13) Amendment 13b's C: IsoCon's Amendment 8 arm M right D reads over the truth-free attainable D
#                        reads (a record on a surviving copy in R.bam, or attributed into their own family's net by this run)
#                        and the A13-1 verdict -> comparator.out, comparator.json                                        light
#   decompose <run>      post hoc (no registered rule): arm M read by read vs IsoCon's Amendment 8 arm M, D reads split by
#                        where the R arm put them -> decompose.out, decompose.json                                     light
#   clean <run>          rm M.splice.mmi and M.fa (the BAMs and tables stay)
#
# As run on 2026-10-02 (ACC=a12): copies; nets; batches.txt = batch 0 written by hand, then `plan 4 1` (calibrated on it); stage reg 0..4;
# concat, contigs, label, mindex reg; malign reg 0..2; mmerge reg; score reg; keep reg; report reg; decompose reg; clean reg. Reruns (the
# heavy lock was shared with another session, so shorter calls): `plan 9 0 $A/half/batches.txt`, copied to $A/double/; split; stage
# half|double 0..8; concat, contigs, label, mindex; malign half 0a..2b / double 0..2; mmerge; score; report; clean.
# As run on 2026-10-03 (ACC=a13, the o3_candidates of 0f5824a7): link; stage reg 0..4 (A12's five batches); concat, nets, contigs, label,
# mindex reg; malign reg 0..2; mmerge; score; comparator; keep; report; decompose; clean. Reruns: stage half|double 0..4 (the same five
# batches), concat, nets, contigs, label; the flagged contig sets differed from the registered run's (34 of 102 and 43 of 68 unions
# byte-identical), so each got its own mindex, malign 0..2, mmerge, score, report, comparator, clean.
#
# The stage never sees RUSTLE_CACHE_DIR (every run is computed, so the wall times are real). Batching changes the A13 result only through
# ruling R18: the poorly placed reads (Amendment 13b) are those in no net of THIS run's families, and the attribution targets are this run's
# nets + every family's copies, so a read netted in one batch may be attributed in another (`nets` counts the reads in nets of two batches);
# clusters, candidates and genome hits are per family / per consensus; nets.fa's de-duplication across families (R9) is per batch (unused).
set -euo pipefail
REPO=${REPO:-$(cd "$(dirname "$0")/../.." && pwd)}   # this checkout unless the caller sets REPO
ACC=${ACC:-a13}
case $ACC in a12|a13) ;; *) echo "ACC must be a12 or a13" >&2; exit 2 ;; esac
L=/mnt/linuxdisk/tmp/rna_allele/linktest
A12=/mnt/linuxdisk/tmp/rna_allele/a12                 # A12's work dir: the copies table, FASTA, regions and batch plan A13 reuses
A=/mnt/linuxdisk/tmp/rna_allele/$ACC
BIN=${BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release/o3_candidates}   # rebuilt from the commit under test before each run
GGO_MMI=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi
PY="python3 $REPO/bench/rna_allele/accept_o3_candidates.py"
HEAVY="bash $REPO/tools/rlock.sh heavy"
LIGHT="bash $REPO/tools/rlock.sh light"

dir_of() { case $1 in reg) echo $A ;; half) echo $A/half ;; double) echo $A/double ;; *) echo "run must be reg|half|double" >&2; exit 2 ;; esac; }
delta_args() { case $1 in reg) echo "" ;; half) echo "--delta 0.00479" ;; double) echo "--delta 0.01916" ;; esac; }
only() { [ "$ACC" = "$1" ] || { echo "step $step is ACC=$1 only" >&2; exit 2; }; }

step=${1:?step}; shift
case $step in copies|link|plan|split) run=reg ;; *) run=${1:?run}; shift ;; esac
W=$(dir_of "$run"); mkdir -p "$W/logs"
case $step in
  copies) only a12; $LIGHT python3 $REPO/bench/rna_allele/panel_to_copies.py --panel $L/panel.json --bam $L/R.bam --fasta $L/masked.fa --out $A/A12 ;;
  link)   # Amendment 13: substrate and scoring unchanged — the same copies table, FASTA and regions; the same five batches
    only a13
    for f in A12.copies.tsv A12.copies.fa A12.regions batches.txt; do ln -sfn $A12/$f $A/$f; done
    for f in R.bam R.bam.bai; do ln -sfn $L/$f $A/$f; done
    ls -l $A | grep -- "->" ;;
  plan)   # A12's plan is a registered product: never rewritten (under a13 the default batch file is a link to it, so name another)
    only a13; BF=${3:-$A/batches.txt}
    [ ! -L "$BF" ] || { echo "plan: $BF is a link (to A12's registered batch plan): name another batch file" >&2; exit 2; }
    $LIGHT $PY plan --w $A --batches "${1:-7}" --keep "${2:-0}" --batch-file "$BF" ;;
  split)  # each scored part in two halves (alternate records): shorter heavy calls when the lock is contended; per-read results unchanged
    mkdir -p $A/parts
    for p in 0 1 2; do
      $LIGHT gawk -v a=$A/parts/scored.part${p}a.fa -v b=$A/parts/scored.part${p}b.fa \
        '/^>/ { n++ } { if (n % 2) print > a; else print > b }' $L/scored.part$p.fa
    done
    grep -c ">" $A/parts/*.fa ;;
  stage)
    only a13                                                     # A12's registered stage products (a12/cand_g*) are never overwritten
    BF=$W/batches.txt; [ -s $BF ] || BF=$A/batches.txt          # a run may carry its own batch plan (same families, other grouping)
    g=${1:?batch}; fams=$(sed -n "$((g + 1))p" $BF); [ -n "$fams" ] || { echo "no batch $g in $BF" >&2; exit 2; }
    # stderr (the stage's log and time -v) timestamped line by line (epoch s) for the per-family timing
    env -u RUSTLE_CACHE_DIR $HEAVY /usr/bin/time -v $BIN --bam $L/R.bam --fasta $L/masked.fa --copies $A/A12.copies.tsv --copies-fa $A/A12.copies.fa \
      --index $L/masked.splice.mmi --out $W/cand_g$g --threads 4 $(delta_args "$run") --families "$fams" 2>&1 >/dev/null \
      | gawk '{ print systime() "\t" $0; fflush() }' > $W/logs/stage_g$g.log
    grep -E "Elapsed \(wall|Maximum resident|pass B:|done:" $W/logs/stage_g$g.log ;;
  concat)  $LIGHT $PY concat --w $W --prefix cand ;;
  nets)    only a13; $HEAVY $PY nets --w $W --linktest $L --prefix cand --copies $A/A12.copies.tsv --copies-fa $A/A12.copies.fa > $W/attrib.out
           cat $W/attrib.out ;;
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
  report) # A12: the registered run's k-mer nets for every run; A13: each run's own nets (pass A and n_net do not depend on delta; `nets` checks)
    if [ "$ACC" = a12 ]; then NETS=$A; else NETS=$W; fi
    $LIGHT $PY report --w $W --linktest $L --prefix cand --nets $NETS > $W/report.out; cat $W/report.out ;;
  comparator) only a13; $LIGHT $PY comparator --w $W --linktest $L --prefix cand > $W/comparator.out; cat $W/comparator.out ;;
  decompose) $LIGHT $PY decompose --w $W --linktest $L > $W/decompose.out; cat $W/decompose.out ;;
  clean)  rm -f $W/M.splice.mmi $W/M.fa ;;
  *) echo "unknown step $step" >&2; exit 2 ;;
esac
