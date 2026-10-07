#!/bin/bash
# bench/seed_pool/run.sh — the runner of docs/PREREG_seed_pool_real_reads_2026-10-07.md: which alignments seed the assembly (primary / good / all)
# and which transcript represents a locus (most reads / primary-first), on real reads, one substrate at a time. Every arm is the pipeline driver
# (tools/rustle_pipeline.sh assemble, then families) with one --seed-pool setting, so the comparison is the switch the advisor can flip.
#
#   run.sh arm SUB ARM        assemble + families of one POOL arm (P, G100, ..., A; an unknown name or an *_R1 name is refused). Resumable: a step is skipped only
#                             when its DONE marker (PREFIX.assemble.done / PREFIX.families.done: the arm, its flags and the substrate) matches
#   run.sh r1 SUB ARM         Rule 1 on top of the pool arm ARM (needs the P arm): locus_units.py --rule1, then the driver's families on it -> ARM_R1
#   run.sh all SUB [ARMS]     every primary arm, then the descriptive ones, in order; stops with exit 75 when this call's time is used or the
#                             all-vs-all of an A-type arm needs another bounded call: run the same command again
#   run.sh g4 SUB             G4: the P arm clustered again through the shard wrapper (under P_SH), compared with the single-process run
#   run.sh gates SUB          G0 (human_chr16 only), G1, G1b, G2, G3 from the products
#   run.sh score SUB          composition (own nodes, node classes), copy_support (twice, PYTHONHASHSEED 0/1), family_score, for every finished arm
#   run.sh table SUB          the comparison tables, the per-copy matrix, the non-dominated set
#   run.sh adopt SUB          write the DONE markers for arms made before the markers existed (the human arms and the first gorilla arms), from their names
#   run.sh clean SUB          delete the PAF and loci FASTA of finished arms (the shard cache of an arm that is still running is never touched)
# SUB = human_chr16 (NPIP) | human_chr17 (TBC1D3) | gorilla_npip (NC_073241.2 + NC_073242.2) | gorilla_tbc1d3 (NC_073228.2 + NC_073224.2); the gorilla
# substrates and their instruments are fixed by Amendment 1 of the prereg (no cap signal: found = Amendment A; no family truth: M4 not scored).
# ARM = P | G100 | G995 | G98 | G95 | G90 | A   (G<rho>: --seed-pool good --seed-as-ratio 1 / .995 / .98 / .95 / .90; A: --seed-pool all)
# Environment: RS_BIN (release dir), RS_WORK (products), RS_CALL_S (seconds of one `all` call, default 330). Nothing is read from a RUSTLE_* variable
# of the calling shell. The work disk is nearly full: a step refuses to start with less than 8 GB free.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
BIN=${RS_BIN:-/mnt/linuxdisk/home/juanfraitu/rustle_target_m2/release}
W=${RS_WORK:-/mnt/linuxdisk/tmp/seed_pool_2026-10-07}
HB=/mnt/linuxdisk/home/juanfraitu
ANN=/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat/ann
FROZEN=/mnt/linuxdisk/tmp/readpool_npip
DEF=/mnt/linuxdisk/tmp/rescore_2026-10-06/DEF
export TMPDIR=$W/tmp; mkdir -p "$TMPDIR"
T_CALL0=$(date +%s)    # a call must end inside the 10-minute tool limit: the shard wrapper's deadline is anchored here
export RLOCK_TIMEOUT=${RLOCK_TIMEOUT:-570}
leak=$(env | grep '^RUSTLE_' || true)
[ -z "$leak" ] || { echo "[seed_pool] refusing: RUSTLE_* set in the calling environment: $leak" >&2; exit 2; }
heavy() { bash "$REPO/tools/rlock.sh" heavy "$@"; }
light() { bash "$REPO/tools/rlock.sh" light "$@"; }
say() { echo "[seed_pool] $(date +%H:%M:%S) $*" >&2; }

cmd=${1:?arm|r1|all|g4|gates|score|table|adopt|clean}; shift
SUB=${1:?substrate}; shift
case "$SUB" in
human_chr16) CONTIG=chr16; FAMILY=NPIP;   BAM=$HB/winloci_data/A119b.t2t.bam; FA=$HB/winloci_data/chm13v2.0.fa; TRUTHS=compara,u2,soto
             COPIES=$ANN/copies.hsa.tsv; TRUTH=$ANN/truth.hsa.gtf; EXONS=(--exons-json "$FROZEN/npip_read_pool.json")
             TARGETS="compara:CF153;u2:ID_154,ID_149,ID_151;soto:ID_154,ID_149"      # the truth families that hold the NPIP copies (rows of the family scores)
             TABLE_SRC=/mnt/linuxdisk/tmp/rustle_figures/runs/human_A119b/human_A119b.molecules.tsv;;
human_chr17) CONTIG=chr17; FAMILY=TBC1D3; BAM=$HB/winloci_data/A119b.t2t.bam; FA=$HB/winloci_data/chm13v2.0.fa; TRUTHS=compara,soto
             COPIES=$ANN/copies.hsa.tsv; TRUTH=$ANN/truth.hsa.gtf; EXONS=()
             TARGETS="compara:CF185;soto:ID_468,ID_469"                              # the truth families that hold the TBC1D3 copies
             TABLE_SRC=/mnt/linuxdisk/tmp/rustle_figures/runs/human_A119b/human_A119b.molecules.tsv;;
gorilla_npip|gorilla_tbc1d3)
             GANN=/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools/ann
             BAM=$HB/winloci_data/GGO_mm.bam; FA=$HB/_from_wsl/winloci_scratch/GGO.fasta; TRUTHS=""; FOUND=ann; EXONS=(); TRUTH=$GANN/truth.ggo.gtf
             TABLE_SRC=/mnt/linuxdisk/tmp/rustle_figures/runs/gorilla_OR6737/gorilla_OR6737.molecules.tsv; MOLNAME=mol_gorilla.tsv
             case "$SUB" in gorilla_npip) CONTIG=NC_073241.2,NC_073242.2; FAMILY=NPIP;; *) CONTIG=NC_073228.2,NC_073224.2; FAMILY=TBC1D3;; esac
             # only the copies on the assembled contigs are the family's truth here (Amendment 1: NC_073244.2 is not assembled)
             mkdir -p "$W/$SUB"
             python3 - "$GANN/copies.ggo.tsv" "$FAMILY" "$CONTIG" "$W/$SUB/copies.tsv" <<'PY'
import csv, sys
src, fam, contigs, out = sys.argv[1:5]
keep = set(contigs.split(","))
rows = [r for r in csv.DictReader(open(src), delimiter="\t") if r["family"] == fam and r["chrom"] in keep]
with open(out, "w") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(rows)
PY
             COPIES=$W/$SUB/copies.tsv;;
*) echo "unknown substrate $SUB" >&2; exit 2;;
esac
S=$W/$SUB; mkdir -p "$S" "$W/mol"
FOUND=${FOUND:-tc}; MOLNAME=${MOLNAME:-mol.tsv}; TARGETS=${TARGETS:-}
[ -e "$W/mol/$MOLNAME" ] || ln -s "$TABLE_SRC" "$W/mol/$MOLNAME"      # the .asbin sidecar lands here, not in the stored run

pool_args() {   # ARM -> the driver's pool flags
  case $1 in
    P)    echo "--seed-pool primary";;
    A)    echo "--seed-pool all";;
    G100) echo "--seed-pool good --seed-as-ratio 1";;
    G995) echo "--seed-pool good --seed-as-ratio 0.995";;
    G98)  echo "--seed-pool good --seed-as-ratio 0.98";;
    G95)  echo "--seed-pool good --seed-as-ratio 0.95";;
    G90)  echo "--seed-pool good --seed-as-ratio 0.90";;
    *) echo "unknown arm $1" >&2; return 1;;
  esac
}
prefix() { echo "$S/$1/run"; }
# DONE markers: a step is finished only if its marker names this arm, these flags and this substrate (a non-empty file is not proof: a crashed or foreign run leaves one)
marker_text() { echo "arm=$1 flags=$2 substrate=$SUB contig=$CONTIG"; }
done_ok() { [ -s "$1" ] && [ "$(head -1 "$1")" = "$2" ]; }
r1_text() {   # ARM: Rule 1 over the pool arm ARM with the primaries of P; the sha1s tie the marker to the two inputs
  echo "arm=${1}_R1 flags=rule1(base=$1,primary=P) substrate=$SUB contig=$CONTIG base=$(sha1sum "$(prefix "$1").families.gtf" | cut -c1-12) primary=$(sha1sum "$(prefix P).families.gtf" | cut -c1-12)"
}
need_disk() { local free; free=$(df -BG --output=avail "$W" | tail -1 | tr -dc 0-9); [ "$free" -ge 8 ] || { echo "[seed_pool] only ${free} GB free on the work disk" >&2; exit 2; }; }
stamp() { for b in copy_assign mcl_families family_score as_table; do echo "$b	$(sha1sum "$BIN/$b" | cut -d' ' -f1)"; done; for f in tools/rustle_pipeline.sh bench/entangled/locus_units.py bench/copy_support.py bench/seed_pool/composition.py bench/seed_pool/run.sh bench/seed_pool/table.py bench/seed_pool/gates.py; do echo "$f	$(sha1sum "$REPO/$f" | cut -d' ' -f1)"; done; }
# A-type arms always go through the shard wrapper; so does every gorilla arm that widens the pool (G98 on gorilla NPIP took 6 minutes in one process, P 41 s)
sharded() { case $1 in A|A_R1) return 0;; P|P_R1) return 1;; *) case "$SUB" in gorilla_*) return 0;; *) return 1;; esac;; esac; }

# the driver's families stage on PREFIX (assembled there); the all-vs-all of A-type arms goes through the shard wrapper in a bounded call
BUDGET_RE="budget: shard|hit the deadline|budget exhausted before the index|is gone: stopping"
do_families() {   # ARM PREFIX
  local arm=$1 P=$2 fam_env=() rc=0 wl="" before=0 after=0 use_shard=0
  if sharded "$arm" || [ -n "${FORCE_SHARD:-}" ]; then
    use_shard=1
    # mcl_families executes RUSTLE_MINIMAP2 directly and tools/mm2_shard.sh is not an executable file in the repo (mode 644): a shim
    mkdir -p "$W/bin"; printf '#!/bin/bash\nexec bash %q "$@"\n' "$REPO/tools/mm2_shard.sh" > "$W/bin/mm2_shard.sh"; chmod +x "$W/bin/mm2_shard.sh"
    fam_env=(RUSTLE_MINIMAP2="$W/bin/mm2_shard.sh" MM2_SHARD_DIR="$W/mm2_shard" MM2_SHARD_BUDGET_S=530 MM2_SHARD_DEADLINE=$(( T_CALL0 + 555 )))
    [ -z "${FORCE_SHARD:-}" ] || fam_env+=(MM2_SHARD_MIN_BYTES=1000)
  fi
  # THIS arm's wrapper.log: the newest one whose start line names this arm's loci FASTA (the cache is shared by every arm and substrate)
  wrapper_log() { local f; for f in $(ls -t "$W"/mm2_shard/*/*/q*/t*/wrapper.log 2>/dev/null); do if grep -qF "target=$P.fam.loci.fa " "$f"; then echo "$f"; return 0; fi; done; return 0; }
  progress() { local f; f=$(wrapper_log); if [ -n "$f" ]; then grep -cE ' done: |index: built' "$f" || true; else echo 0; fi; }
  [ "$use_shard" = 0 ] || before=$(progress)
  rm -f "$P".fam.* "$P.families.log"
  /usr/bin/time -v env "${fam_env[@]}" bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" families --bam "$BAM" --fasta "$FA" --out "$P" \
    --bin "$BIN" --threads 4 --no-cache > "$P.families.driver.log" 2> "$P.families.driver.stderr" || rc=$?
  if [ "$rc" != 0 ]; then
    # a budget stop keeps its progress and is worth another call; a call that made none (a shard that fits no call: exit 76, or any other failure) is a failure
    if [ "$use_shard" = 1 ]; then
      wl=$(wrapper_log); after=$(progress)
      if [ -n "$wl" ] && tail -2 "$wl" | grep -qE "$BUDGET_RE" && [ "$after" -gt "$before" ]; then
        say "$arm: the all-vs-all needs another bounded call ($(tail -1 "$wl" | sed 's/^[0-9-]* [0-9:]* //')): run the same command again"
        exit 75
      fi
    fi
    say "$arm: families failed (exit $rc), see $P.families.driver.stderr and $P.families.log"; exit "$rc"
  fi
  # the aligner's intermediates are re-creatable and the work disk is nearly full; the record count and the wrapper's log are kept
  [ ! -s "$P.fam.loci.paf" ] || wc -l < "$P.fam.loci.paf" > "$P.paf_records"
  rm -f "$P.fam.loci.paf" "$P.fam.loci.fa"
  if [ "$use_shard" = 1 ]; then
    wl=$(wrapper_log)
    if [ -n "$wl" ]; then cp "$wl" "$P.wrapper.log"; rm -rf "${wl%/*/*/*/*}"; fi          # only this arm's key directory
  fi
  say "$arm: families $(grep Elapsed "$P.families.driver.stderr" | awk '{print $NF}'), $(grep -h 'families:' "$P.families.driver.stderr" | sed 's/.*families: //' | head -1)"
}

step_arm() {   # ARM
  local arm=$1 pa P mt
  case "$arm" in *_R1) echo "[seed_pool] arm: $arm is a Rule-1 arm: use 'r1 SUB ${arm%_R1}'" >&2; exit 2;; esac
  pa=$(pool_args "$arm") || exit 2
  P=$(prefix "$arm"); mt=$(marker_text "$arm" "$pa"); mkdir -p "$(dirname "$P")"; need_disk
  if [ ! -s "$P.gtf" ] || [ ! -s "$P.families.gtf" ] || ! done_ok "$P.assemble.done" "$mt"; then
    rm -f "$P.assemble.done" "$P.families.done"
    { echo "date	$(date -Is)"; echo "arm	$arm	$pa"; echo "substrate	$SUB"; echo "head	$(git -C "$REPO" rev-parse --short HEAD)"; stamp; } > "$P.run.log"
    # shellcheck disable=SC2086
    /usr/bin/time -v bash "$REPO/tools/rlock.sh" heavy bash "$REPO/tools/rustle_pipeline.sh" assemble --bam "$BAM" --fasta "$FA" --out "$P" --bin "$BIN" --threads 4 \
      --no-cache --contig "$CONTIG" --as-table "$W/mol/$MOLNAME" $pa > "$P.assemble.driver.log" 2> "$P.assemble.driver.stderr"
    echo "$mt" > "$P.assemble.done"
    say "$arm: assemble $(grep Elapsed "$P.assemble.driver.stderr" | awk '{print $NF}'), $(grep -h 'assemble:' "$P.assemble.driver.stderr" | grep transcripts | tail -1 | sed 's/.*assemble: //')"
  else
    say "$arm: assemble skipped (finished)"
  fi
  if [ ! -s "$P.fam.clusters.tsv" ] || [ "$P.fam.clusters.tsv" -ot "$P.families.gtf" ] || ! done_ok "$P.families.done" "$mt"; then
    do_families "$arm" "$P"; echo "$mt" > "$P.families.done"
  else
    say "$arm: families skipped (finished)"
  fi
}

step_r1() {   # ARM (the pool arm Rule 1 is applied to)
  local arm=$1 B Pp R rt
  case "$arm" in *_R1) echo "[seed_pool] r1: give the POOL arm (P, G98, ...), not $arm" >&2; exit 2;; esac
  pool_args "$arm" > /dev/null || exit 2
  B=$(prefix "$arm"); Pp=$(prefix P); R=$(prefix "${arm}_R1"); mkdir -p "$(dirname "$R")"; need_disk
  [ -s "$Pp.families.gtf" ] && [ -s "$B.families.gtf" ] || { echo "[seed_pool] r1 $arm: run the P arm and the $arm arm first" >&2; exit 2; }
  rt=$(r1_text "$arm")
  if [ ! -s "$R.families.gtf" ] || [ "$R.families.gtf" -ot "$B.families.gtf" ] || ! done_ok "$R.rule1.done" "$rt"; then
    rm -f "$R.rule1.done" "$R.families.done"
    cp "$B.gtf" "$R.gtf"        # the driver's guard wants the families input to be no older than the assembly
    sleep 1
    python3 "$REPO/bench/entangled/locus_units.py" --base "$B.families.gtf" --out "$R.families.gtf" --primary-gtf "$Pp.families.gtf" --rule1 > "$R.rule1.log" 2>&1
    { echo "date	$(date -Is)"; echo "arm	${arm}_R1	Rule 1 over $arm, primaries of P"; echo "substrate	$SUB"; stamp; } > "$R.run.log"
    echo "$rt" > "$R.rule1.done"
    say "${arm}_R1: Rule 1 applied ($(tail -1 "$R.rule1.log" | cut -c1-100))"
  fi
  if [ ! -s "$R.fam.clusters.tsv" ] || [ "$R.fam.clusters.tsv" -ot "$R.families.gtf" ] || ! done_ok "$R.families.done" "$rt"; then
    do_families "${arm}_R1" "$R"; echo "$rt" > "$R.families.done"
  else
    say "${arm}_R1: families skipped (finished)"
  fi
}

case "$cmd" in
arm) step_arm "${1:?ARM}";;
r1)  step_r1 "${1:?ARM}";;
all)
  t0=$(date +%s); limit=${RS_CALL_S:-330}
  case "$SUB" in
    # Amendment 1: the sweep ends G100 / G90 and their Rule-1 arms are not run on gorilla. Amendment 2 (user, 2026-10-07): the all-secondaries
    # arms are deferred; run them later with `run.sh all gorilla_npip "A_R1"` and `run.sh all gorilla_tbc1d3 "A A_R1"` (gorilla NPIP's A had finished)
    gorilla_*) list=${1:-"P P_R1 G98 G98_R1 G995 G95"};;
    *)         list=${1:-"P P_R1 G98 G98_R1 A A_R1 G100 G995 G95 G90 G100_R1 G995_R1 G95_R1 G90_R1"};;
  esac
  for a in $list; do
    el=$(( $(date +%s) - t0 ))
    # an A-type arm (its assembly and a sharded all-vs-all of up to 9 minutes) starts only at the start of a call
    if sharded "$a" && [ "$el" -gt 15 ]; then say "$a needs a call of its own; run the same command again"; exit 75; fi
    if [ "$el" -gt "$limit" ]; then say "this call's time is used; run the same command again"; exit 75; fi
    case "$a" in *_R1) step_r1 "${a%_R1}";; *) step_arm "$a";; esac
  done
  say "all arms of $SUB finished";;
g4)
  # G4: the all-vs-all through the shard wrapper equals a single minimap2 run (P is clustered again, sharded, under another prefix)
  Pp=$(prefix P); Q=$S/P_SH/run; mkdir -p "$(dirname "$Q")"; [ -s "$Pp.fam.clusters.tsv" ] || { echo "run the P arm first" >&2; exit 2; }
  cp "$Pp.gtf" "$Q.gtf"; sleep 1; cp "$Pp.families.gtf" "$Q.families.gtf"
  FORCE_SHARD=1 do_families P_SH "$Q"
  if cmp -s "$Pp.fam.clusters.tsv" "$Q.fam.clusters.tsv"; then echo "G4 PASS: sharded and unsharded clusters are identical ($SUB)"; else echo "G4 FAIL: the clusters differ ($SUB)"; exit 1; fi;;
gates)  python3 "$HERE/gates.py" --work "$S" --sub "$SUB" --def "$DEF" --bin "$BIN" --copies "$COPIES" --truth "$TRUTH" --family "$FAMILY" "${EXONS[@]}" ;;
score)
  arms=(); for d in "$S"/*/; do a=$(basename "$d"); case "$a" in *_SH) continue;; esac; [ -s "$d/run.fam.clusters.tsv" ] && arms+=("$a"); done
  [ ${#arms[@]} -gt 0 ] || { echo "no finished arm in $S" >&2; exit 2; }
  armspec=(); lociarg=()
  for a in "${arms[@]}"; do P=$(prefix "$a"); armspec+=(--arm "$a=$P.fam.loci.gff3,$P.fam.clusters.tsv"); lociarg+=(--loci "$a=$P.fam.loci.gff3,$P.families.gtf"); done
  python3 "$HERE/composition.py" --copies "$COPIES" --truth "$TRUTH" --family "$FAMILY" "${EXONS[@]}" "${armspec[@]}" --out "$S/nodes.json" --report "$S/comp.json" \
    --emit-node-loci "$S/nodeloci"
  for seed in 0 1; do
    PYTHONHASHSEED=$seed light python3 "$REPO/bench/copy_support.py" --copies "$COPIES" --truth "$TRUTH" --bam "$BAM" --family "$FAMILY" --out "$S/support_s$seed" \
      "${lociarg[@]}" --nodes "$S/nodes.json" > "$S/support_s$seed.stdout"
  done
  cp "$S/support_s0.json" "$S/support.json"; cp "$S/support_s0.copies.tsv" "$S/support.copies.tsv"
  # M2n: the same instrument over the NODE loci only (the instrument's own M2 looks at every locus that overlaps the copy)
  nodearg=(); for a in "${arms[@]}"; do nodearg+=(--loci "$a=$S/nodeloci/$a.nodes.gff3"); done
  PYTHONHASHSEED=0 light python3 "$REPO/bench/copy_support.py" --copies "$COPIES" --truth "$TRUTH" --bam "$BAM" --family "$FAMILY" --out "$S/support_nodeonly" \
    "${nodearg[@]}" --nodes "$S/nodes.json" > "$S/support_nodeonly.stdout"
  for a in "${arms[@]}"; do
    [ -n "$TRUTHS" ] || continue          # no family truth on this substrate: M4 is not scored
    light python3 "$REPO/bench/default_rescore/score_families.py" --fs "$BIN/family_score" --clusters "$(prefix "$a").fam.clusters.tsv" --contig "${CONTIG%%,*}" --label "$a" \
      --out "$S/fs_$a.json" --work "$S/fs" --truths "$TRUTHS" > "$S/fs_$a.stdout"
  done
  say "scored ${#arms[@]} arms: ${arms[*]}";;
table) python3 "$HERE/table.py" --work "$S" --family "$FAMILY" --truths "$TRUTHS" --found "$FOUND" --targets "$TARGETS" ;;
adopt)
  # markers for arms that were made before the markers existed: the name decides what they are claimed to be; nothing is checked beyond the files being there
  for a in P G100 G995 G98 G95 G90 A; do
    P=$(prefix "$a"); [ -s "$P.fam.clusters.tsv" ] || continue
    mt=$(marker_text "$a" "$(pool_args "$a")"); echo "$mt" > "$P.assemble.done"; echo "$mt" > "$P.families.done"; say "$a: adopted"
  done
  for a in P G100 G995 G98 G95 G90 A; do
    R=$(prefix "${a}_R1"); [ -s "$R.fam.clusters.tsv" ] || continue
    rt=$(r1_text "$a"); echo "$rt" > "$R.rule1.done"; echo "$rt" > "$R.families.done"; say "${a}_R1: adopted"
  done;;
clean)
  # only finished arms; the shard cache of an arm that is still running is never touched (do_families removes a finished arm's own key directory)
  for d in "$S"/*/; do [ -s "$d/run.fam.clusters.tsv" ] && rm -f "$d"run.fam.loci.paf "$d"run.fam.loci.fa; done; say "cleaned";;
*) echo "usage: run.sh arm|r1|all|g4|gates|score|table|adopt|clean SUB ..." >&2; exit 2;;
esac
