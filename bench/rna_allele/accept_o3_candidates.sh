#!/bin/bash
# accept_o3_candidates.sh — the `o3_candidates` acceptance on Amendment 7's 53-family held-out (docs/PREREG_rna_allele_haplotype_count_2026-10-01.md),
# scored by `merge_test.py score` (arm M = masked genome + one union per flagged candidate, each candidate its own component).
#   ACC=a13 (default): Amendment 13 (+ 13b-13e) — the A13 stage (alignment attribution, structural template), work dir a13/; the copies table,
#                      FASTA, regions and batch plan of A12 are reused (`link`); A13-1 is decided on the re-registered comparator C (`comparator`).
#   ACC=a12:           Amendment 12, as run on 2026-10-02 into a12/ (its k-mer net recomputation `nets` and binary a3564999 are at commit fde90c0a;
#                      the current stage has no k-mer rule, so only its scoring steps re-run here).
#   ACC=a14:           Amendment 14 — the no-deletion control of the stage (ruling R22): Amendment 9's substrate (the UNMASKED `_pri` GGO.fasta +
#                      GGO.splice.mmi, control/R0.bam, a copies table of all 201 copies), the same five batches, work dir a14/; every flagged union
#                      classified against KB3781's mat / pat by Amendment 9's rule (`hap`, `overlap`, `classify`: C1'), arm C = `_pri` + the flagged
#                      unions (`cindex`, `calign`, `cmerge`, `cscore`: C2'), `creport`; and ruling R23's whole-BAM cost (`wcopies`, `wplan`,
#                      `wstage`, `wreport`) in a14/wholebam/. Only run `reg`.
#
#   [ACC=a12|a13|a14] accept_o3_candidates.sh <step> [run] [arg]    run = reg (delta 0.00958, the registered run, work dir $A) | half (0.00479,
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
#   clean <run>          rm M.splice.mmi and M.fa (the BAMs and tables stay); a14: C.splice.mmi and C.fa
# ACC=a14 only (the steps above that also run under a14: copies, link, stage, concat, nets, contigs, label, report is replaced by creport):
#   copies               panel.json -> A14.copies.tsv / .fa / .regions: EVERY copy (--all), sequences from the unmasked `_pri`, n_reads by R0.bam light
#   hap reg              iso.contigs.fa vs mat and pat (minimap2 -c -x splice:hq -uf -N 20, Amendment 9's command) -> iso.mat.paf, iso.pat.paf  heavy
#   overlap reg          iso.contigs.fa vs A9's new-copy contigs (control/contigs_L.fa) and vs the A13 run's survivor-derived unions
#                        (minimap2 -c -x asm20 -N 50 -p 0.1, Amendment 9's overlap command) -> overlap_a9.paf, overlap_a13.paf          heavy (tiny)
#   classify reg         C1' (Amendment 9's a / b / c rule; self-check on A9's own candidates) -> classify.out, candidates_classified.tsv  light
#   cindex reg           C.fa = GGO.fasta + iso.contigs.fa; minimap2 -x splice -d C.splice.mmi (~3 min, ~21 GB)                 heavy
#   calign reg <p>       scored.part<p>.fa (p = 0..2) realigned to C with R0's flags (Amendment 9's arm C command) -> C.<p>.bam      heavy
#   cmerge reg           C.{0,1,2}.bam -> C.bam (+ .bai)                                                                       heavy
#   cscore reg           C2' (control_test.py score's classes, each flagged union its own locus; self-check: R0 = A9's) -> score.out   heavy
#   creport reg          wall time, stage totals, counters, pass B, the joined reads' fate, per-family table -> report.out        light
#   wcopies              ruling R23: refabsent/panel.json (378 families, 915 copies) -> wholebam/W.copies.tsv / .fa (--all, `_pri`,
#                        n_reads by the full fibroblast BAM)                                                                  heavy
#   wplan [size]         wholebam/batches.txt: the families in table order, `size` (50) per batch                              light
#   wstage <g>           the stage on the FULL fibroblast BAM, batch g of wholebam/batches.txt, time -v -> wholebam/cand_g<g>.*;
#                        stopped by SIGINT at WSTAGE_LIMIT (570 s; time -v still reports) when it cannot finish in a 10-min call  heavy
#   wreport              per batch: Elapsed, peak RSS, attribution set and targets, flagged -> wholebam/wreport.out             light
#
# As run on 2026-10-02 (ACC=a12): copies; nets; batches.txt = batch 0 written by hand, then `plan 4 1` (calibrated on it); stage reg 0..4;
# concat, contigs, label, mindex reg; malign reg 0..2; mmerge reg; score reg; keep reg; report reg; decompose reg; clean reg. Reruns (the
# heavy lock was shared with another session, so shorter calls): `plan 9 0 $A/half/batches.txt`, copied to $A/double/; split; stage
# half|double 0..8; concat, contigs, label, mindex; malign half 0a..2b / double 0..2; mmerge; score; report; clean.
# As run on 2026-10-03 (ACC=a14, the o3_candidates rebuilt from 75827d5d = 0f5824a7's logic): copies; link; stage reg 0..4; concat, nets,
# contigs, label, hap, overlap, classify; cindex; calign reg 0..2; cmerge; cscore; creport; clean; wcopies; wplan 50; wstage 0 twice (stopped at
# the time limit both times: the first by the lock's SIGTERM, which /usr/bin/time does not survive, the second by the SIGINT limit; the
# measurement stopped there, ruling R23); wreport.
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
case $ACC in a12|a13|a14) ;; *) echo "ACC must be a12, a13 or a14" >&2; exit 2 ;; esac
L=/mnt/linuxdisk/tmp/rna_allele/linktest
CTRL=/mnt/linuxdisk/tmp/rna_allele/control            # Amendment 9's control: R0.bam (the scored reads on the unmasked `_pri`), copies_lift.tsv, contigs_L.fa
A13=/mnt/linuxdisk/tmp/rna_allele/a13                 # the A13 run (its survivor-derived flags, for Amendment 14's overlap)
GGO_FA=/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta          # the unmasked `_pri`
MAT_MMI=/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.mat.splice.mmi
PAT_MMI=/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.pat.splice.mmi
FIBRO=/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam       # ruling R23: the full gorilla fibroblast BAM (23 GB)
RA=/mnt/linuxdisk/tmp/rna_allele
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
case $step in copies|link|plan|split|wcopies|wplan|wstage|wreport) run=reg ;; *) run=${1:?run}; shift ;; esac
W=$(dir_of "$run"); mkdir -p "$W/logs"
[ "$ACC" != a14 ] || [ "$run" = reg ] || { echo "ACC=a14 runs only reg" >&2; exit 2; }
# the masked run's arm-M steps would mix masked.fa / R.bam with the control's unions: refused under a14 (its arm C: cindex .. creport)
case $step in split|plan|mindex|malign|mmerge|score|keep|report|comparator|decompose) [ "$ACC" != a14 ] || { echo "step $step is the masked run's (ACC=a12/a13); ACC=a14 uses cindex / calign / cmerge / cscore / creport" >&2; exit 2; } ;; esac
WB=$A/wholebam
case $step in
  copies)
    if [ "$ACC" = a14 ]; then   # Amendment 14: every copy of the 53 families (mask + keep), from the unmasked `_pri`; n_reads by R0.bam
      $LIGHT python3 $REPO/bench/rna_allele/panel_to_copies.py --all --panel $L/panel.json --bam $CTRL/R0.bam --fasta $GGO_FA --out $A/A14
    else
      only a12; $LIGHT python3 $REPO/bench/rna_allele/panel_to_copies.py --panel $L/panel.json --bam $L/R.bam --fasta $L/masked.fa --out $A/A12
    fi ;;
  link)   # Amendment 13: substrate and scoring unchanged — the same copies table, FASTA and regions; the same five batches
          # Amendment 14: the same five batches (A12's plan); the BAM is Amendment 9's R0.bam (the stage reads it from $CTRL)
    if [ "$ACC" = a14 ]; then
      ln -sfn $A12/batches.txt $A/batches.txt
    else
      only a13
      for f in A12.copies.tsv A12.copies.fa A12.regions batches.txt; do ln -sfn $A12/$f $A/$f; done
      for f in R.bam R.bam.bai; do ln -sfn $L/$f $A/$f; done
    fi
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
    [ "$ACC" != a12 ] || only a13                                # A12's registered stage products (a12/cand_g*) are never overwritten
    BF=$W/batches.txt; [ -s $BF ] || BF=$A/batches.txt          # a run may carry its own batch plan (same families, other grouping)
    g=${1:?batch}; fams=$(sed -n "$((g + 1))p" $BF); [ -n "$fams" ] || { echo "no batch $g in $BF" >&2; exit 2; }
    if [ "$ACC" = a14 ]; then    # Amendment 14: nothing masked — the scored reads on the unmasked `_pri` (Amendment 9's R0), every copy
      IN="--bam $CTRL/R0.bam --fasta $GGO_FA --copies $A/A14.copies.tsv --copies-fa $A/A14.copies.fa --index $GGO_MMI"
    else
      IN="--bam $L/R.bam --fasta $L/masked.fa --copies $A/A12.copies.tsv --copies-fa $A/A12.copies.fa --index $L/masked.splice.mmi"
    fi
    # stderr (the stage's log and time -v) timestamped line by line (epoch s) for the per-family timing; the first line names the
    # binary's sha1 (final review 2026-10-03, finding 6: a `cargo test` can rebuild the shared target between runs)
    printf '%s\tbinary sha1 %s\n' "$(date +%s)" "$(sha1sum $BIN | cut -c1-40)" > $W/logs/stage_g$g.log
    env -u RUSTLE_CACHE_DIR $HEAVY /usr/bin/time -v $BIN $IN --out $W/cand_g$g --threads 4 $(delta_args "$run") --families "$fams" 2>&1 >/dev/null \
      | gawk '{ print systime() "\t" $0; fflush() }' >> $W/logs/stage_g$g.log
    grep -E "Elapsed \(wall|Maximum resident|pass B:|done:" $W/logs/stage_g$g.log ;;
  concat)  $LIGHT $PY concat --w $W --prefix cand ;;
  nets)    [ "$ACC" != a12 ] || only a13
           if [ "$ACC" = a14 ]; then NA="--copies $A/A14.copies.tsv --copies-fa $A/A14.copies.fa --bam $CTRL/R0.bam"
           else NA="--copies $A/A12.copies.tsv --copies-fa $A/A12.copies.fa"; fi
           $HEAVY $PY nets --w $W --linktest $L --prefix cand $NA > $W/attrib.out
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
  label)  # a14: the unmasked genome IS the reference, so the copy Amendment 7 masked is one more S:<copy> (Amendment 9's copy order)
    $HEAVY /usr/bin/time -v bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $GGO_MMI $W/iso.contigs.fa > $W/iso.base.paf 2> $W/logs/label.mm2.log" > $W/logs/label.log 2>&1
    $LIGHT $PY label --w $W --linktest $L $([ "$ACC" = a14 ] && echo --all-copies) ;;
  hap)    only a14
    $HEAVY /usr/bin/time -v bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $MAT_MMI $W/iso.contigs.fa > $W/iso.mat.paf 2> $W/logs/hap.mat.mm2.log && \
      minimap2 -c -x splice:hq -uf -N 20 -t 4 $PAT_MMI $W/iso.contigs.fa > $W/iso.pat.paf 2> $W/logs/hap.pat.mm2.log" > $W/logs/hap.log 2>&1
    grep -E "Real time|Elapsed \(wall|Maximum resident" $W/logs/hap.*.log ;;
  overlap) only a14   # Amendment 9's overlap command (asm20) against A9's 76 new-copy contigs and the A13 run's survivor-derived unions
    gawk -F'\t' 'NR == FNR { if (FNR > 1 && $8 ~ /^S:/) s[$1] = 1; next } /^>/ { keep = (substr($1, 2) in s) } keep' $A13/contigs.tsv $A13/iso.contigs.fa > $W/a13_S.fa
    $HEAVY bash -c "minimap2 -c -x asm20 -N 50 -p 0.1 -t 2 $CTRL/contigs_L.fa $W/iso.contigs.fa > $W/overlap_a9.paf 2> $W/logs/overlap_a9.log && \
      minimap2 -c -x asm20 -N 50 -p 0.1 -t 2 $W/a13_S.fa $W/iso.contigs.fa > $W/overlap_a13.paf 2> $W/logs/overlap_a13.log"
    echo "A13 survivor-derived unions: $(grep -c '>' $W/a13_S.fa); PAF lines: A9 $(wc -l < $W/overlap_a9.paf), A13 $(wc -l < $W/overlap_a13.paf)" ;;
  classify) only a14; $LIGHT $PY classify --w $W --linktest $L --control $CTRL --a13 $A13 > $W/classify.out; cat $W/classify.out ;;
  cindex) only a14
    $HEAVY /usr/bin/time -v bash -c "cat $GGO_FA $W/iso.contigs.fa > $W/C.fa && minimap2 -x splice -t 4 -d $W/C.splice.mmi $W/C.fa" > $W/logs/cindex.log 2>&1
    grep -E "Real time|Elapsed \(wall|Maximum resident" $W/logs/cindex.log ;;
  calign) only a14     # Amendment 9's arm C command (= R0's flags)
    p=${1:?part}
    $HEAVY /usr/bin/time -v bash -c "set -o pipefail; minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4 $W/C.splice.mmi \
      $L/scored.part$p.fa 2> $W/logs/calign$p.mm2.log | samtools sort -@ 1 -m 500M -o $W/C.$p.bam -" > $W/logs/calign$p.log 2>&1
    grep -E "Real time|Elapsed \(wall|Maximum resident" $W/logs/calign$p.mm2.log $W/logs/calign$p.log ;;
  cmerge) only a14; $HEAVY bash -c "samtools merge -f -@ 2 $W/C.bam $W/C.0.bam $W/C.1.bam $W/C.2.bam && samtools index $W/C.bam" ;;
  cscore) only a14; $HEAVY $PY cscore --w $W --linktest $L --control $CTRL > $W/score.out; cat $W/score.out ;;
  creport) only a14; $LIGHT $PY creport --w $W --linktest $L --prefix cand > $W/report.out; cat $W/report.out ;;
  wcopies) only a14; mkdir -p $WB
    $HEAVY python3 $REPO/bench/rna_allele/panel_to_copies.py --all --panel $RA/refabsent/panel.json --bam $FIBRO --fasta $GGO_FA --out $WB/W ;;
  wplan)  only a14; $LIGHT $PY wplan --w $WB --copies $WB/W.copies.tsv --batch-file $WB/batches.txt --batch-size "${1:-50}"; cat -A $WB/batches.txt | cut -c1-120 ;;
  wstage) only a14   # one sequential sweep of the whole 23-GB BAM per call (pass B). A call that cannot finish is stopped by the inner `timeout -s INT`
                     # (WSTAGE_LIMIT, default 570 s): /usr/bin/time ignores SIGINT, so it still reports Elapsed and the peak RSS of the stage (and of the
                     # minimap2 runs it had reaped); the lock's RLOCK_TIMEOUT (590 s) is the backstop. Exit 124 = stopped.
    g=${1:?batch}; fams=$(sed -n "$((g + 1))p" $WB/batches.txt); [ -n "$fams" ] || { echo "no batch $g in $WB/batches.txt" >&2; exit 2; }
    mkdir -p $WB/logs
    printf '%s\tbinary sha1 %s\n' "$(date +%s)" "$(sha1sum $BIN | cut -c1-40)" > $WB/logs/wstage_g$g.log
    set +e
    env -u RUSTLE_CACHE_DIR RLOCK_TIMEOUT=${RLOCK_TIMEOUT:-590} $HEAVY timeout -s INT -k 10 ${WSTAGE_LIMIT:-570} /usr/bin/time -v $BIN --bam $FIBRO \
      --fasta $GGO_FA --copies $WB/W.copies.tsv \
      --copies-fa $WB/W.copies.fa --index $GGO_MMI --out $WB/cand_g$g --threads 4 --families "$fams" 2>&1 >/dev/null \
      | gawk '{ print systime() "\t" $0; fflush() }' >> $WB/logs/wstage_g$g.log
    rc=${PIPESTATUS[0]}; set -e
    echo "exit $rc"; grep -E "Elapsed \(wall|Maximum resident|exit status|signal|BAM:|pass B:|done:" $WB/logs/wstage_g$g.log ;;
  wreport) only a14; $LIGHT $PY wreport --w $WB --copies $WB/W.copies.tsv --copies-fa $WB/W.copies.fa > $WB/wreport.out; cat $WB/wreport.out ;;
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
  clean)  rm -f $W/M.splice.mmi $W/M.fa $W/C.splice.mmi $W/C.fa ;;
  *) echo "unknown step $step" >&2; exit 2 ;;
esac
