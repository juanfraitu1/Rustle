#!/bin/bash
# accept_o3_candidates.sh — the `o3_candidates` acceptance on Amendment 7's 53-family held-out (docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md),
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
#   ACC=a15h:          Amendment 15's held-out (docs/archive/2026-10/PREREG_rna_allele_haplotype_count_2026-10-01.md, Amendment 15, registered before
#                      any re-run): Amendment 10's read set with nothing deleted — refabsent/R0.bam (32,219 scored reads of 34 families
#                      on the unmasked `_pri`), the 30 disjoint families (refabsent/bonly.tsv's 34 minus GWFAM4 / GWFAM169 / GWFAM175 /
#                      GWFAM402, dev-overlapping), the a14/wholebam 915-copy table restricted to them (H.copies.*), refabsent/panel.json,
#                      the unmasked `_pri` + GGO.splice.mmi; 3 batches of 10 (rlock heavy's timeout is 600 s and ~11 families of this
#                      read set take ~5 min, so A12/A14's five batches are too fine-grained — H carries its own plan), work dir a15h/;
#                      classification as Amendment 9 by control_test.py (`lift`, `classify`; H1: families with >= 1 class b/c/pri flag
#                      <= 4 of 30), arm C = `_pri` + the flagged unions with the read set realigned (`cindex`, `calign`, `cmerge`,
#                      `cscore`; H3: false moves <= 5% of all reads of the 30 families). Reported, not decided: H2 (the expressed
#                      beyond-delta locus GWFAM175_B0 — GWFAM175 is dev-overlapping, so it is NOT among the 30). Only run `reg`.
#
#   [ACC=a12|a13|a14|a15h] accept_o3_candidates.sh <step> [run] [arg] run = reg (delta 0.00958, the registered run, work dir $A) | half (0.00479,
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
# ACC=a15h only (the shared steps stage / concat / nets / contigs run as under a14 on the H inputs; `copies` is refused — the H
#                copies derive from a14/wholebam's table — and so are a14's overlap / wcopies / wplan / wstage / wreport):
#   hfamilies            refabsent/bonly.tsv's 34 families minus the 4 dev-overlapping -> H.families.txt (the registered 30); labels.tsv
#                        restricted to the 30 -> H.labels.tsv (H3's denominator: all reads of the 30 families)                  light
#   hcopies              a14/wholebam/W.copies.{tsv,fa} (915 copies, 378 families) restricted to the 30 -> H.copies.{tsv,fa};
#                        verified: rows > 0, every family present                                                                        light
#   link                 H's batch plan: 3 batches of 10 in H.families.txt order -> batches.txt (not a link: H's own plan); R0.bam
#                        (+ .bai) linked (control_test.py score's arm R0 reads {w}/R0.bam); htruth/ = panel.json + the restricted
#                        labels (control_test.py's --l for lift / classify / score)                                            light
#   hap reg              (as a14) iso.mat.paf / iso.pat.paf; contigs_L.{mat,pat}.paf linked (control_test.classify's file names) heavy
#   lift                 control_test.py lift, reused exactly (its --l = htruth: the refabsent panel; chrmap / the asm5 PAFs stay its
#                        hard-coded TRUTH dir): every copy interval -> its B-haplotype interval -> copies_lift.tsv (classify's b / a) light
#   classify reg         control_test.py classify (Amendment 9's a / b / c / pri rule; its prints keep A9's hard-coded /53 denominator
#                        — noted, not modified: the JSON is the record) -> classify.json; H1 = delta.families_false <= 4 of 30  light
#   cindex reg           (as a14) C.fa = GGO.fasta + iso.contigs.fa; minimap2 -x splice -d C.splice.mmi                        heavy
#   calign reg <p>       (as a14, but the read set = every record of refabsent/R0.bam as FASTA, samtools fasta -0, in 3 parts, one per
#                        heavy call) realigned to C with R0's flags (Amendment 9's arm C command) -> C.<p>.bam               heavy
#   cmerge reg           (as a14) C.{0,1,2}.bam -> C.bam (+ .bai)                                                               heavy
#   cscore reg           control_test.py score (Amendment 9's classes, each flagged union its own locus; arms R0 = refabsent/R0.bam,
#                        C = this run's; no A9 self-check exists for this read set, so control_test.py runs as-is) -> score.json;
#                        H3 = false moves <= 5% of the 30 families' reads                                                    heavy
#   creport reg          (as a14, --linktest htruth) -> report.out                                                              light
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
# As PLANNED, NOT yet run (ACC=a15h, Amendment 15's held-out; registered before any re-run): hfamilies; hcopies; link; stage reg 0..2;
# concat, nets, contigs, label, hap, lift, classify; cindex; calign reg 0..2; cmerge; cscore; creport; clean.
#
# The stage never sees RUSTLE_CACHE_DIR (every run is computed, so the wall times are real). Batching changes the A13 result only through
# ruling R18: the poorly placed reads (Amendment 13b) are those in no net of THIS run's families, and the attribution targets are this run's
# nets + every family's copies, so a read netted in one batch may be attributed in another (`nets` counts the reads in nets of two batches);
# clusters, candidates and genome hits are per family / per consensus; nets.fa's de-duplication across families (R9) is per batch (unused).
set -euo pipefail
REPO=${REPO:-$(cd "$(dirname "$0")/../.." && pwd)}   # this checkout unless the caller sets REPO
ACC=${ACC:-a13}
case $ACC in a12|a13|a14|a15h) ;; *) echo "ACC must be a12, a13, a14 or a15h" >&2; exit 2 ;; esac
L=/mnt/linuxdisk/tmp/rna_allele/linktest
CTRL=/mnt/linuxdisk/tmp/rna_allele/control            # Amendment 9's control: R0.bam (the scored reads on the unmasked `_pri`), copies_lift.tsv, contigs_L.fa
A13=${A13:-/mnt/linuxdisk/tmp/rna_allele/a13}   # the A13 run (its survivor-derived flags, for Amendment 14's overlap); overridable so a re-run never overwrites it
GGO_FA=/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta          # the unmasked `_pri`
MAT_MMI=/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.mat.splice.mmi
PAT_MMI=/mnt/linuxdisk/home/juanfraitu/winloci_data/mGorGor1.pat.splice.mmi
FIBRO=/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam       # ruling R23: the full gorilla fibroblast BAM (23 GB)
RA=/mnt/linuxdisk/tmp/rna_allele
A12=/mnt/linuxdisk/tmp/rna_allele/a12                 # A12's work dir: the copies table, FASTA, regions and batch plan A13 reuses
A=${A:-/mnt/linuxdisk/tmp/rna_allele/$ACC}   # the run's work dir; overridable (Amendment 15 re-runs use a13a15/a14a15, leaving the registered a13/a14 products in place)
BIN=${BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release/o3_candidates}   # rebuilt from the commit under test before each run
GGO_MMI=/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO.splice.mmi
PY="python3 $REPO/bench/rna_allele/accept_o3_candidates.py"
CT="python3 $REPO/bench/rna_allele/control_test.py"
A14W=${A14W:-$RA/a14/wholebam}   # ruling R23's copies table (915 copies, 378 families); Amendment 15's held-out restricts it to the 30
HEAVY="bash $REPO/tools/rlock.sh heavy"
LIGHT="bash $REPO/tools/rlock.sh light"

dir_of() { case $1 in reg) echo $A ;; half) echo $A/half ;; double) echo $A/double ;; *) echo "run must be reg|half|double" >&2; exit 2 ;; esac; }
delta_args() { case $1 in reg) echo "" ;; half) echo "--delta 0.00479" ;; double) echo "--delta 0.01916" ;; esac; }
only() { [ "$ACC" = "$1" ] || { echo "step $step is ACC=$1 only" >&2; exit 2; }; }

step=${1:?step}; shift
case $step in copies|link|plan|split|wcopies|wplan|wstage|wreport|hfamilies|hcopies) run=reg ;; *) run=${1:?run}; shift ;; esac
W=$(dir_of "$run"); mkdir -p "$W/logs"
case $ACC in a14|a15h) [ "$run" = reg ] || { echo "ACC=$ACC runs only reg" >&2; exit 2; } ;; esac
# the masked run's arm-M steps would mix masked.fa / R.bam with the control's unions: refused under a14 / a15h (their arm C:
# cindex .. cscore; the masked run's `report` is refused, a14 / a15h report with creport)
case $step in split|plan|mindex|malign|mmerge|score|keep|report|comparator|decompose)
  case $ACC in a14|a15h) echo "step $step is the masked run's (ACC=a12/a13); ACC=a14/a15h use cindex / calign / cmerge / cscore / creport" >&2; exit 2 ;; esac ;; esac
WB=$A/wholebam
case $step in
  copies)
    case $ACC in
      a15h) echo "step copies is ACC=a12/a14 only; ACC=a15h derives its copies from $A14W with hfamilies / hcopies" >&2; exit 2 ;;
      a14)  # Amendment 14: every copy of the 53 families (mask + keep), from the unmasked `_pri`; n_reads by R0.bam
        $LIGHT python3 $REPO/bench/rna_allele/panel_to_copies.py --all --panel $L/panel.json --bam $CTRL/R0.bam --fasta $GGO_FA --out $A/A14 ;;
      *)    only a12; $LIGHT python3 $REPO/bench/rna_allele/panel_to_copies.py --panel $L/panel.json --bam $L/R.bam --fasta $L/masked.fa --out $A/A12 ;;
    esac ;;
  hfamilies)  # Amendment 15's held-out: the 34 families of refabsent/bonly.tsv minus the 4 dev-overlapping = the registered 30
    only a15h
    $LIGHT python3 - "$RA/refabsent/bonly.tsv" "$RA/refabsent/labels.tsv" "$A/H.families.txt" "$A/H.labels.tsv" <<'EOF'
import sys
bonly, labels, fam_out, lab_out = sys.argv[1:5]
DEV = {"GWFAM4", "GWFAM169", "GWFAM175", "GWFAM402"}   # among the 53 development families (Amendment 15)
fams = []
for ln in open(bonly):
    f = ln.split("\t")
    if f[0] == "family" or not f[0].strip():
        continue
    if f[0] not in fams:
        fams.append(f[0])
h = [f for f in fams if f not in DEV]
assert len(fams) == 34 and len(h) == 30, (len(fams), len(h))
open(fam_out, "w").write("".join(f + "\n" for f in h))
lines = open(labels).read().splitlines()      # H3's denominator: all reads of the 30 families
sel = [ln for ln in lines[1:] if ln and ln.split("\t", 2)[1] in set(h)]
assert sel, "no labels of the 30 families"
open(lab_out, "w").write(lines[0] + "\n" + "\n".join(sel) + "\n")
print(f"{len(fams)} families in bonly.tsv; minus the 4 dev-overlapping: {len(h)} -> {fam_out}; their {len(sel)} scored reads -> {lab_out}")
# H2 (reported, not decided): the expressed beyond-delta locus GWFAM175_B0 (Amendment 10's D1) — GWFAM175 is dev-overlapping,
# so it is NOT among the 30 and this harness does not decide H2.
EOF
    ;;
  hcopies)  # ruling R23's 915-copy / 378-family table restricted to the 30 families -> H.copies.tsv / H.copies.fa (the stage's --copies)
    only a15h
    $LIGHT python3 - "$A/H.families.txt" "$A14W/W.copies.tsv" "$A14W/W.copies.fa" "$A/H.copies.tsv" "$A/H.copies.fa" <<'EOF'
import sys
fam_file, wtsv, wfa, htsv, hfa = sys.argv[1:6]
h = {ln.strip() for ln in open(fam_file) if ln.strip()}
lines = open(wtsv).read().splitlines()
hdr, rows = lines[0], [ln for ln in lines[1:] if ln and ln.split("\t", 1)[0] in h]
have = {ln.split("\t", 1)[0] for ln in rows}
assert rows and have == h, f"rows {len(rows)}, families {len(have)} != {len(h)}"
open(htsv, "w").write("\n".join([hdr] + rows) + "\n")
out, keep = [], False
for ln in open(wfa):
    if ln.startswith(">"):
        keep = ln[1:].split("|", 1)[0] in h
    if keep:
        out.append(ln)
open(hfa, "w").writelines(out)
nrec = sum(1 for ln in out if ln.startswith(">"))
assert nrec == len(rows), (nrec, len(rows))
print(f"{len(rows)} rows / {nrec} FASTA records of the {len(h)} families -> {htsv}, {hfa}")
EOF
    ;;
  link)   # Amendment 13: substrate and scoring unchanged — the same copies table, FASTA and regions; the same five batches
          # Amendment 14: the same five batches (A12's plan); the BAM is Amendment 9's R0.bam (the stage reads it from $CTRL)
          # Amendment 15's held-out: H's OWN plan (3 batches of 10; the heavy timeout is 600 s and ~11 families of this read set take
          # ~5 min, so A12/A14's five batches are too fine-grained), written here, not linked
    if [ "$ACC" = a14 ]; then
      ln -sfn $A12/batches.txt $A/batches.txt
    elif [ "$ACC" = a15h ]; then
      [ -s $A/H.families.txt ] || { echo "link: run hfamilies first" >&2; exit 2; }
      [ -s $A/H.copies.tsv ] || { echo "link: run hcopies first" >&2; exit 2; }
      [ ! -L $A/batches.txt ] || rm -f $A/batches.txt
      gawk 'NR % 10 == 1 { if (s != "") print s; s = $1; next } { s = s "," $1 } END { print s }' $A/H.families.txt > $A/batches.txt
      for f in R0.bam R0.bam.bai; do ln -sfn $RA/refabsent/$f $A/$f; done   # control_test.py score's arm R0 reads {w}/R0.bam
      mkdir -p $A/htruth   # control_test.py's --l for lift / classify / score: the panel + the labels restricted to the 30
      ln -sfn $RA/refabsent/panel.json $A/htruth/panel.json               # (H3's denominator: all reads of the 30 families)
      ln -sfn $A/H.labels.tsv $A/htruth/labels.tsv
      wc -l $A/batches.txt; cat $A/batches.txt
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
    elif [ "$ACC" = a15h ]; then # Amendment 15's held-out: Amendment 10's scored reads on the unmasked `_pri`, the H copies
      IN="--bam $RA/refabsent/R0.bam --fasta $GGO_FA --copies $A/H.copies.tsv --copies-fa $A/H.copies.fa --index $GGO_MMI"
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
           LT=$L; [ "$ACC" != a15h ] || LT=$A/htruth   # a15h: labels of the refabsent read set (restricted to the 30)
           if [ "$ACC" = a14 ]; then NA="--copies $A/A14.copies.tsv --copies-fa $A/A14.copies.fa --bam $CTRL/R0.bam"
           elif [ "$ACC" = a15h ]; then NA="--copies $A/H.copies.tsv --copies-fa $A/H.copies.fa --bam $RA/refabsent/R0.bam"
           else NA="--copies $A/A12.copies.tsv --copies-fa $A/A12.copies.fa"; fi
           $HEAVY $PY nets --w $W --linktest $LT --prefix cand $NA > $W/attrib.out
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
  label)  # a14 / a15h: the unmasked genome IS the reference, so every copy is one more S:<copy> (Amendment 9's copy order); a15h's
          # panel and labels are the refabsent read set's (restricted to the 30 families)
    LT=$L; ALLO=""
    if [ "$ACC" = a15h ]; then LT=$A/htruth; ALLO=--all-copies; elif [ "$ACC" = a14 ]; then ALLO=--all-copies; fi
    $HEAVY /usr/bin/time -v bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $GGO_MMI $W/iso.contigs.fa > $W/iso.base.paf 2> $W/logs/label.mm2.log" > $W/logs/label.log 2>&1
    $LIGHT $PY label --w $W --linktest $LT $ALLO ;;
  hap)
    [ "$ACC" = a14 ] || [ "$ACC" = a15h ] || { echo "step hap is ACC=a14/a15h only" >&2; exit 2; }
    $HEAVY /usr/bin/time -v bash -c "minimap2 -c -x splice:hq -uf -N 20 -t 4 $MAT_MMI $W/iso.contigs.fa > $W/iso.mat.paf 2> $W/logs/hap.mat.mm2.log && \
      minimap2 -c -x splice:hq -uf -N 20 -t 4 $PAT_MMI $W/iso.contigs.fa > $W/iso.pat.paf 2> $W/logs/hap.pat.mm2.log" > $W/logs/hap.log 2>&1
    if [ "$ACC" = a15h ]; then   # a15h's classifier is control_test.classify, which reads A9's names (accept's classify reads iso.*.paf)
      ln -sfn iso.mat.paf $W/contigs_L.mat.paf; ln -sfn iso.pat.paf $W/contigs_L.pat.paf
    fi
    grep -E "Real time|Elapsed \(wall|Maximum resident" $W/logs/hap.*.log ;;
  overlap) only a14   # Amendment 9's overlap command (asm20) against A9's 76 new-copy contigs and the A13 run's survivor-derived unions
    gawk -F'\t' 'NR == FNR { if (FNR > 1 && $8 ~ /^S:/) s[$1] = 1; next } /^>/ { keep = (substr($1, 2) in s) } keep' $A13/contigs.tsv $A13/iso.contigs.fa > $W/a13_S.fa
    $HEAVY bash -c "minimap2 -c -x asm20 -N 50 -p 0.1 -t 2 $CTRL/contigs_L.fa $W/iso.contigs.fa > $W/overlap_a9.paf 2> $W/logs/overlap_a9.log && \
      minimap2 -c -x asm20 -N 50 -p 0.1 -t 2 $W/a13_S.fa $W/iso.contigs.fa > $W/overlap_a13.paf 2> $W/logs/overlap_a13.log"
    echo "A13 survivor-derived unions: $(grep -c '>' $W/a13_S.fa); PAF lines: A9 $(wc -l < $W/overlap_a9.paf), A13 $(wc -l < $W/overlap_a13.paf)" ;;
  lift)   # a15h: every copy interval of the panel lifted to its B-haplotype interval through the frozen truth's asm5 alignments —
          # control_test.py's lift reused exactly (its --l = htruth: refabsent/panel.json; chrmap / the PAFs are its hard-coded TRUTH
          # dir); classify's b_allele / a_haplotype_only rule reads the result -> copies_genes.tsv, copies_lift.tsv
    only a15h; $LIGHT $CT lift --w $A --l $A/htruth ;;
  classify)
    [ "$ACC" = a14 ] || [ "$ACC" = a15h ] || { echo "step classify is ACC=a14/a15h only" >&2; exit 2; }
    if [ "$ACC" = a15h ]; then
      # H1 (Amendment 15): families with >= 1 class b/c/pri flag <= 4 of the 30 = classify.json's delta.families_false.
      # control_test.py's prints keep A9's hard-coded /53 denominator — noted, not modified (the JSON is the record).
      $LIGHT $CT classify --w $W --l $A/htruth > $W/classify.out
    else
      $LIGHT $PY classify --w $W --linktest $L --control $CTRL --a13 $A13 > $W/classify.out
    fi; cat $W/classify.out ;;
  cindex)
    [ "$ACC" = a14 ] || [ "$ACC" = a15h ] || { echo "step cindex is ACC=a14/a15h only" >&2; exit 2; }
    $HEAVY /usr/bin/time -v bash -c "cat $GGO_FA $W/iso.contigs.fa > $W/C.fa && minimap2 -x splice -t 4 -d $W/C.splice.mmi $W/C.fa" > $W/logs/cindex.log 2>&1
    grep -E "Real time|Elapsed \(wall|Maximum resident" $W/logs/cindex.log ;;
  calign) # Amendment 9's arm C command (= R0's flags)
    [ "$ACC" = a14 ] || [ "$ACC" = a15h ] || { echo "step calign is ACC=a14/a15h only" >&2; exit 2; }
    p=${1:?part}
    if [ "$ACC" = a15h ]; then
      # the read set = EVERY record of Amendment 10's R0.bam (mapped and unmapped) as FASTA, in 3 parts (one per heavy call);
      # samtools fasta -0 takes the unpaired reads (all of them here). (refabsent/scored.part*.fa is the same read set in the
      # same 3 parts, written by Amendment 10 — the parts are rebuilt from the BAM here so this step stands alone.)
      mkdir -p $W/parts
      if [ ! -s $W/parts/H.part2.fa ]; then
        $LIGHT samtools fasta -0 $W/parts/H.all.fa $RA/refabsent/R0.bam
        $LIGHT gawk -v d=$W/parts '/^>/ { n++ } { print > (d "/H.part" (n % 3) ".fa") }' $W/parts/H.all.fa
        grep -c ">" $W/parts/H.part*.fa
      fi
      Q=$W/parts/H.part$p.fa
    else
      Q=$L/scored.part$p.fa
    fi
    $HEAVY /usr/bin/time -v bash -c "set -o pipefail; minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4 $W/C.splice.mmi \
      $Q 2> $W/logs/calign$p.mm2.log | samtools sort -@ 1 -m 500M -o $W/C.$p.bam -" > $W/logs/calign$p.log 2>&1
    grep -E "Real time|Elapsed \(wall|Maximum resident" $W/logs/calign$p.mm2.log $W/logs/calign$p.log ;;
  cmerge)
    [ "$ACC" = a14 ] || [ "$ACC" = a15h ] || { echo "step cmerge is ACC=a14/a15h only" >&2; exit 2; }
    $HEAVY bash -c "samtools merge -f -@ 2 $W/C.bam $W/C.0.bam $W/C.1.bam $W/C.2.bam && samtools index $W/C.bam" ;;
  cscore)
    [ "$ACC" = a14 ] || [ "$ACC" = a15h ] || { echo "step cscore is ACC=a14/a15h only" >&2; exit 2; }
    if [ "$ACC" = a15h ]; then
      # H3: control_test.py score's C2 over the 30 families' reads (--l htruth: labels restricted to the 30) — arms R0 (the linked
      # refabsent/R0.bam) and C (this run's); each flagged union its own locus (merge/paf/*.paf are label's empty files, as a14);
      # no A9 self-check exists for this read set, so control_test.py runs as-is -> score.json
      $HEAVY $CT score --w $W --l $A/htruth > $W/score.out
    else
      $HEAVY $PY cscore --w $W --linktest $L --control $CTRL > $W/score.out
    fi; cat $W/score.out ;;
  creport)
    [ "$ACC" = a14 ] || [ "$ACC" = a15h ] || { echo "step creport is ACC=a14/a15h only" >&2; exit 2; }
    LT=$L; [ "$ACC" != a15h ] || LT=$A/htruth
    # a15h's classify / cscore are control_test.py's (A9's layouts): candidates_classified.tsv / calls.tsv do not exist; creport
    # guards on both (accept_o3_candidates.py) and reports the stage / nets / pass-B numbers without them
    $LIGHT $PY creport --w $W --linktest $LT --prefix cand > $W/report.out; cat $W/report.out ;;
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
