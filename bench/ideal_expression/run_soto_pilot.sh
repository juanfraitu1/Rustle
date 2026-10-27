#!/bin/bash
# bench/ideal_expression/run_soto_pilot.sh -- the runner of docs/PREREG_ideal_arm_pilot_2026-10-07.md and its amendment: the ideal arm of the expression ladder on the Soto genes, at two depths.
# One stage per call, every call foreground and at most 590 s (rlock.sh kills a command at RLOCK_TIMEOUT = 560 s and gives up on a held lock after RLOCK_WAIT = 30 s; use a Bash timeout of 600000).
#   run_soto_pilot.sh [--dry-run] STAGE [DEPTH]      DEPTH = reads per transcript, one of RS_DEPTHS (default "10 30"); --dry-run prints each command with its lock class and runs nothing
#
#   truth        sim_genes.py at every depth, seed 20261008 (RS_SEED), PYTHONHASHSEED=0, ONE read pool: depth 10 is the first 10 reads of every transcript of depth 30 (light, or heavy if RS_PROCS > 1)
#   truthqc      G1 manifest + canonical introns and G2 read ends at every depth, and G-depth (the shallow arm is the deep pool's first reads, read for read) (light)
#   map          the DEEPEST pool, once: minimap2 (the shipped command, -N 50) against the whole-genome splice index, read-disjoint parts, parts = ceil(read Mb / 35), ONE PART PER CALL (heavy);
#                re-run until it prints "mapping complete"
#   cut          the shallower arms' BAMs = the pool's BAM cut by read name (samtools view -N), checked: one primary-or-unmapped record per read (light)
#   bamqc D      G4 from the BAM alone (light); required BEFORE the pipeline
#   assemble D   the driver's `assemble`, HEAD release binaries, defaults, a NEW --out per arm (heavy); prints the ASSEMBLY POLISH counts of the log
#   families D   the driver's `families` (heavy) through the resumable shard wrapper (an executable shim and a shard directory per arm); exit 75 = the budget is spent, call again with the
#                same command; exit 1 with "no progress" = a call finished no more steps (index, shards) than the one before: one step does not fit in the budget, do not call again
#   anno D       the controls: the annotation as nodes and GTF (ctrl0) and shifted by 50 kb (ctrl50k) (light)
#   score D      score_soto_nodes.py, five variants (node = primary, transcript, gene_id locus, compete soto, raw GTF); the calibration depth (RS_CUT_DEPTH, default the shallowest) FREEZES the
#                single-exon cut in $W/single_cut.json, the others read it (light)
#   score-ctrl D the two controls through the same scorer (light)
#   score-real   the real-library products (old A119b GTF, A119b f1v2 GTF, current-default testis run) through the same scorer, on the ideal arm's truth (light)
#   report D     the primary score again with every other run beside it (COMPARATORS) (light)
#   gates-sim    G3, simulation half: tables re-simulated under PYTHONHASHSEED=1 and a 120-gene FASTQ under 0 and 1 (light)
#   gates D      G1-G5 recomputed from the files, G6 (positive control), G7 (negative control), G-depth, the scorer re-run under two hash seeds; prints VALID or INVALID (light)
#   cleanup      deletes the pool's mapping part files by exact name, once its BAM is complete (light)
#   status       what exists, mapping progress (help and status only read: they create nothing)
# Environment: RS_WORK (products), RS_DEPTHS, RS_SEED, RS_PAD (25000), RS_PROCS (1), RS_PARTS, RS_FORCE=1 (let `truth` overwrite), RS_SHARD=0 (no shard wrapper), RS_SHARD_BUDGET (480), RS_CUT_DEPTH,
# RS_BIN (HEAD release dir), RS_SRC_OK (commits known not to change the binaries, for G5), RS_BASELINE_PER_GENE (the 2026-09-22 per-gene table, names), RS_GENOME/RS_GFF/RS_INDEX/RS_GENES/RS_S1C/
# RS_SPLIT (inputs), RS_REAL_A119B_OLD, RS_REAL_A119B_F1V2, RS_REAL_TESTIS (prefix), RS_A119B_READS (comparators). Nothing is read from a RUSTLE_* variable of the calling shell: the runner
# refuses to start if one is set. Heavy stages run under `tools/rlock.sh heavy`, light ones under `tools/rlock.sh light`; no background job, no pkill, files deleted only by exact name.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
HB=${RS_HOME:-/mnt/linuxdisk/home/juanfraitu}
BIN=${RS_BIN:-$HB/rustle_target_m2/release}
W=${RS_WORK:-/mnt/linuxdisk/tmp/ideal_arm_pilot_2026-10-07}
FA=${RS_GENOME:-$HB/winloci_data/chm13v2.0.fa}
GFF=${RS_GFF:-$HB/winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz}
INDEX=${RS_INDEX:-$HB/npip_ladder/idx/target.splice.mmi}
GENES=${RS_GENES:-$REPO/bench/soto/soto_famCN_S1C.tsv}
S1C=${RS_S1C:-$GENES}
SPLIT=${RS_SPLIT:-$REPO/bench/soto/soto_split_2026-09-29.tsv}
BASE_FID=${RS_BASELINE_PER_GENE-/mnt/linuxdisk/tmp/gw22/soto_fid.per_gene.tsv}
REAL_OLD=${RS_REAL_A119B_OLD:-/mnt/linuxdisk/tmp/gw22/ours_genome.gtf}
REAL_F1V2=${RS_REAL_A119B_F1V2:-/mnt/linuxdisk/tmp/rustle_figures_dev/f1_port/runs/human_A119b/new.f1v2.families.gtf}
REAL_TESTIS=${RS_REAL_TESTIS:-/mnt/linuxdisk/tmp/rustle_figures_dev/defaults_flip/runs/new/human_testis}
A119B_READS=${RS_A119B_READS:-/mnt/linuxdisk/tmp/expr_ladder/prior/a119b_reads_per_soto_gene.tsv}
SEED=${RS_SEED:-20261008}; PAD=${RS_PAD:-25000}; PROCS=${RS_PROCS:-1}
DEPTHS=${RS_DEPTHS:-"10 30"}
POOL=$(printf '%s\n' $DEPTHS | sort -n | tail -1); FIRST=$(printf '%s\n' $DEPTHS | sort -n | head -1)
CUT_DEPTH=${RS_CUT_DEPTH:-$FIRST}
DRIVER="$REPO/tools/rustle_pipeline.sh"
export RLOCK_TIMEOUT=${RLOCK_TIMEOUT:-560} RLOCK_WAIT=${RLOCK_WAIT:-30}
leak=$(env | grep '^RUSTLE_' || true)
[ -z "$leak" ] || { echo "[soto_pilot] refusing: RUSTLE_* set in the calling environment: $leak" >&2; exit 2; }
DRY=0
if [ "${1:-}" = "--dry-run" ]; then DRY=1; shift; fi
cmd=${1:-help}; DEPTH=${2:-}
export TMPDIR=$W/tmp
PD=$W/d$POOL; PP=$PD/reads          # the deepest pool: truth, FASTQ and the one BAM
CUTF=$W/single_cut.json
sha() { if [ -f "$1" ]; then sha1sum "$1" | cut -d' ' -f1; else echo missing; fi; }
say() { echo "[soto_pilot] $*" >&2; }
heavy() { if [ "$DRY" = 1 ]; then echo "DRY heavy: $*"; else bash "$REPO/tools/rlock.sh" heavy "$@"; fi; }
light() { if [ "$DRY" = 1 ]; then echo "DRY light: $*"; else bash "$REPO/tools/rlock.sh" light "$@"; fi; }
logto() { if [ "$DRY" = 1 ]; then cat; else tee "$1"; fi; }
logq() { if [ "$DRY" = 1 ]; then cat; else tee "$1" > /dev/null; fi; }   # log to a file, print nothing
need() { # file message
  [ "$DRY" = 1 ] && return 0
  [ -s "$1" ] || { echo "[soto_pilot] $2" >&2; exit 2; }
}
gate_ok() { # json key-or-empty message: the gate recorded in the file must have passed
  [ "$DRY" = 1 ] && return 0
  python3 - "$1" "$2" <<'PY' || { echo "[soto_pilot] $3" >&2; exit 2; }
import json, sys
d = json.load(open(sys.argv[1]))
ok = d[sys.argv[2]]["ok"] if sys.argv[2] else d["ok"]
sys.exit(0 if ok else 1)
PY
}
stamp() {
  for b in copy_assign mcl_families as_table family_score; do echo "$b	$(sha "$BIN/$b")"; done
  for f in sim_genes.py sim_windows.py score_soto_nodes.py map_reads.py run_soto_pilot.sh; do echo "$f	$(sha "$HERE/$f")"; done
  for f in rustle_pipeline.sh rlock.sh mm2_shard.sh; do echo "$f	$(sha "$REPO/tools/$f")"; done
  echo "head	$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo none)"
  echo "minimap2	$(minimap2 --version 2>/dev/null || echo none)"
  echo "samtools	$(samtools --version 2>/dev/null | head -1 || echo none)"
  echo "python	$(python3 --version 2>&1)"
}
begin_log() { # file mode(w|a) extra-lines...
  [ "$DRY" = 1 ] && return 0
  local f=$1 mode=$2; shift 2
  { echo "date	$(date -Is)"; for l in "$@"; do echo "$l"; done; stamp; } | { if [ "$mode" = a ]; then cat >> "$f"; else cat > "$f"; fi; }
}
need_depth() {
  [ -n "$DEPTH" ] || { echo "[soto_pilot] stage $cmd needs a depth, one of: $DEPTHS" >&2; exit 2; }
  case " $DEPTHS " in *" $DEPTH "*) ;; *) echo "[soto_pilot] depth $DEPTH is not one of: $DEPTHS" >&2; exit 2;; esac
}
arm() { DD=$W/d$1; P=$DD/reads; A=$DD/asm; C0=$DD/ctrl0; C5=$DD/ctrl50k; }
sim_args=(--genes "$GENES" --gff "$GFF" --genome "$FA" --s1c "$S1C" --split "$SPLIT" --pad "$PAD" --seed "$SEED")
case "$cmd" in help|-h|--help|status) ;; *) [ "$DRY" = 1 ] || mkdir -p "$W" "$TMPDIR" ;; esac     # help and status only read

# the ids of the baseline's 1,078 RefSeq-name population, joined to the S1C ids by Gene Name (1,345 ids): a labelled sub-population, never the denominator
baseline_ids() {
  [ -n "$BASE_FID" ] && [ -s "$BASE_FID" ] || return 1
  [ -s "$W/baseline_ids.txt" ] || python3 - "$BASE_FID" "$S1C" "$W/baseline_ids.txt" <<'PY'
import csv, sys
names = {r["gene"] for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
ids = sorted({r["Gene ID"] for r in csv.DictReader(open(sys.argv[2], newline=""), delimiter="\t") if r["Gene Name"] in names})
open(sys.argv[3], "w").write("Gene ID\n" + "\n".join(ids) + "\n")
PY
}
score_gate() { # the cut must be frozen before any arm but the calibration arm
  if [ "$DRY" != 1 ] && [ ! -e "$CUTF" ] && [ "$DEPTH" != "$CUT_DEPTH" ]; then
    echo "[soto_pilot] the single-exon cut is not frozen yet: freeze it on the calibration arm first (run: score $CUT_DEPTH)" >&2; exit 2
  fi
}
score_needs() {
  need "$A.fam.loci.gff3" "families first (the nodes are the families stage's fam.loci.gff3)"
  need "$A.families.gtf" "families first"
  need "$A.gtf" "assemble first"
  need "$P.bam" "map first"
}
score_common() { # sets SC (array): truth, bam, assemble log, cut, subset
  SC=(--truth "$P" --bam "$P.bam" --assemble-log "$A.assemble.log")
  if [ -e "$CUTF" ]; then SC+=(--single-cut-file "$CUTF"); elif [ "$DEPTH" = "$CUT_DEPTH" ]; then SC+=(--freeze-single-cut "$CUTF"); fi
  if [ "$DRY" != 1 ] && baseline_ids; then SC+=(--subset "$W/baseline_ids.txt" --subset-label "baseline name-join population (names -> ids)"); fi
}
score_one() { # name show(show|quiet) args...
  local name=$1 show=$2; shift 2
  if [ "$show" = show ]; then light python3 "$HERE/score_soto_nodes.py" --out "$DD/$name" "$@" | logto "$DD/$name.log"
  else light python3 "$HERE/score_soto_nodes.py" --out "$DD/$name" "$@" | logq "$DD/$name.log"; fi
}

case "$cmd" in
help|-h|--help)
  awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0" ;;
truth)
  for d in $DEPTHS; do
    if [ "$DRY" != 1 ] && [ -e "$W/d$d/reads.fq" ] && [ "${RS_FORCE:-0}" != 1 ]; then
      echo "[soto_pilot] refusing to overwrite $W/d$d/reads.fq: the replicate exists (everything downstream is keyed to it); RS_FORCE=1 overrides" >&2; exit 2
    fi
  done
  begin_log "$W/truth.run.log" w "seed	$SEED" "pad	$PAD" "depths	$DEPTHS" "procs	$PROCS"
  # one process is a light job (the full set at depth 10 was measured at 95 MB and under a minute); more processes take the heavy lock
  if [ "$PROCS" -gt 1 ]; then lock=heavy; else lock=light; fi
  for d in $DEPTHS; do
    [ "$DRY" = 1 ] || mkdir -p "$W/d$d"
    $lock env PYTHONHASHSEED=0 python3 "$HERE/sim_genes.py" --out "$W/d$d/reads" "${sim_args[@]}" --reads-per-tx "$d" --procs "$PROCS" | logto "$W/d$d/sim.log"
  done ;;
truthqc)
  need "$PP.fq" "truth first"
  for d in $DEPTHS; do
    light python3 "$HERE/sim_genes.py" gates --out "$W/d$d/reads" --genome "$FA" --genes "$GENES" | logto "$W/d$d/g12.log"
  done
  for d in $DEPTHS; do
    [ "$d" = "$POOL" ] && continue
    light python3 "$HERE/sim_genes.py" depth-subset --json "$W/depth_subset_$d.json" "$W/d$d/reads" "$PP" | logto "$W/depth_subset_$d.log"
  done
  if [ "$DRY" != 1 ] && [ "$POOL" != "$FIRST" ]; then
    python3 - "$W" "$POOL" $DEPTHS <<'PY'
import json, sys
w, pool = sys.argv[1], sys.argv[2]
parts = [json.load(open(f"{w}/depth_subset_{d}.json")) for d in sys.argv[3:] if d != pool]
json.dump(dict(ok=all(p["ok"] for p in parts), lines=[ln for p in parts for ln in p["lines"]]), open(f"{w}/depth_subset.json", "w"), indent=1)
PY
  fi ;;
map)
  need "$PP.fq" "truth first"
  for d in $DEPTHS; do
    need "$W/d$d/reads.gates.json" "truthqc first (G1 and G2 are checked before the reads are mapped)"
    gate_ok "$W/d$d/reads.gates.json" G1 "G1 did not pass at depth $d: fix the simulation before mapping (see $W/d$d/g12.log)"
    gate_ok "$W/d$d/reads.gates.json" G2 "G2 did not pass at depth $d: fix the simulation before mapping (see $W/d$d/g12.log)"
  done
  if [ "$POOL" != "$FIRST" ]; then gate_ok "$W/depth_subset.json" "" "G-depth did not pass: the shallow arms are not the deep pool's first reads (see $W/depth_subset_*.log)"; fi
  [ "$DRY" = 1 ] || [ -s "$INDEX" ] || { echo "[soto_pilot] no splice index at $INDEX" >&2; exit 2; }
  if [ -n "${RS_PARTS:-}" ]; then PARTS=$RS_PARTS
  elif [ -s "$PD/map.parts" ]; then PARTS=$(cat "$PD/map.parts")
  elif [ -s "$PP.sim.json" ]; then PARTS=$(python3 -c "import json,math,sys; print(max(1, math.ceil(json.load(open(sys.argv[1]))['read_bases'] / 35e6)))" "$PP.sim.json")
  else PARTS=1; fi
  if [ "$DRY" != 1 ] && [ ! -s "$PD/map.parts" ]; then echo "$PARTS" > "$PD/map.parts"; fi
  begin_log "$PP.map.run.log" a "parts	$PARTS"
  heavy python3 "$HERE/map_reads.py" "$PP" "$INDEX" --parts "$PARTS" --threads 4 --max-parts-per-call 1 2>&1 | tail -2 ;;
cut)
  need "$PP.bam" "map first (the pool's BAM is the only mapping)"
  for d in $DEPTHS; do
    [ "$d" = "$POOL" ] && continue
    arm "$d"
    n=0; [ "$DRY" = 1 ] || n=$(( $(wc -l < "$P.fq") / 4 ))
    if [ "$DRY" != 1 ] && [ -s "$P.bam" ] && [ -s "$P.bam.bai" ] && [ "$(light samtools view -c -F 2304 "$P.bam")" = "$n" ]; then say "depth $d: the cut BAM exists and is complete"; continue; fi
    light python3 "$HERE/sim_genes.py" readnames --out "$P" --to "$P.names.txt"
    light samtools view -b -N "$P.names.txt" -o "$P.cut.tmp.bam" "$PP.bam"
    if [ "$DRY" != 1 ]; then mv "$P.cut.tmp.bam" "$P.bam"; fi
    light samtools index "$P.bam"
    if [ "$DRY" != 1 ]; then
      got=$(light samtools view -c -F 2304 "$P.bam")
      [ "$got" = "$n" ] || { echo "[soto_pilot] the cut BAM holds $got primary-or-unmapped records for the $n reads of depth $d" >&2; exit 1; }
      say "depth $d: $got reads cut from the pool by name"
    fi
  done ;;
bamqc)
  need_depth; arm "$DEPTH"
  need "$PP.bam" "map first"
  need "$P.bam" "cut first"
  light python3 "$HERE/sim_genes.py" bamqc --out "$P" --bam "$P.bam" | logto "$DD/g4.log" ;;
assemble)
  need_depth; arm "$DEPTH"
  need "$P.bam" "map first"
  need "$P.g4.json" "bamqc first (G4 is checked before the pipeline)"
  gate_ok "$P.g4.json" "" "G4 did not pass: the BAM is not a valid mapping of the reads (see $DD/g4.log)"
  begin_log "$A.run.log" w "bam	$P.bam" "depth	$DEPTH"
  if [ "$DRY" = 1 ]; then echo "DRY heavy: bash $DRIVER assemble --bam $P.bam --fasta $FA --out $A --bin $BIN --threads 4"
  else
    # the driver reuses any best-AS table whose header names the same BAM path, whatever the BAM now holds: a stale table of this prefix is removed (exact names), so the table is built from this BAM
    rm -f "$A.molecules.tsv" "$A.molecules.tsv.asbin"
    /usr/bin/time -v bash "$REPO/tools/rlock.sh" heavy bash "$DRIVER" assemble --bam "$P.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 > "$A.assemble.driver.log" 2> "$A.assemble.driver.stderr"
    say "assemble: $(grep -h 'assemble:' "$A.assemble.driver.stderr" | tail -1 | sed 's/.*assemble: //')"
    grep -h 'ASSEMBLY POLISH' "$A.assemble.log" > "$DD/polish.txt" || true
    cat "$DD/polish.txt" >&2
  fi ;;
families)
  need_depth; arm "$DEPTH"
  need "$A.gtf" "assemble first"
  need "$A.families.gtf" "assemble first (the families input is asm.families.gtf)"
  if [ "$DRY" != 1 ]; then
    for s in clusters.tsv loci.tsv loci.gff3 loci.fa loci.paf copies.tsv copies.fa copies.regions copies.merged.tsv params.tsv; do rm -f "$A.fam.$s"; done
    rm -f "$A.families.log"
  fi
  begin_log "$DD/families.run.log" a "shard	${RS_SHARD:-1}" "shard_budget_s	${RS_SHARD_BUDGET:-480}"
  shard=()
  if [ "${RS_SHARD:-1}" = 1 ]; then
    # the PAF cache key holds the minimap2 command path: one shim path (and one shard directory) per arm, the same in every call
    if [ "$DRY" != 1 ]; then
      mkdir -p "$DD/mm2"
      printf '#!/bin/bash\nexec bash "%s/tools/mm2_shard.sh" "$@"\n' "$REPO" > "$DD/mm2_shard_exec.sh"; chmod +x "$DD/mm2_shard_exec.sh"
    fi
    shard=(env "RUSTLE_MINIMAP2=$DD/mm2_shard_exec.sh" "MM2_SHARD_DIR=$DD/mm2" "MM2_SHARD_BUDGET_S=${RS_SHARD_BUDGET:-480}" "MM2_SHARD_DEADLINE=$(( $(date +%s) + ${RS_SHARD_DEADLINE:-500} ))")
  fi
  if [ "$DRY" = 1 ]; then echo "DRY heavy: ${shard[*]:-} bash $DRIVER families --bam $P.bam --fasta $FA --out $A --bin $BIN --threads 4"
  else
    t0=$(date +%s.%N); rc=0
    /usr/bin/time -v bash "$REPO/tools/rlock.sh" heavy ${shard[@]+"${shard[@]}"} bash "$DRIVER" families --bam "$P.bam" --fasta "$FA" --out "$A" --bin "$BIN" --threads 4 > "$A.families.driver.log" 2> "$A.families.driver.stderr" || rc=$?
    if [ "$rc" != 0 ] || [ ! -s "$A.fam.clusters.tsv" ]; then
      # a spent shard budget: the driver exits 1 with no clusters file, and the newest shard wrapper log, written in THIS call, ends in a budget / deadline line. Call again -- unless this call
      # finished no more steps (the index, a shard) than the previous one: then one step does not fit in the budget, the wrapper's exit 76 looks like 75 from here, and calling again would spin.
      prc=0
      python3 - "$DD/mm2" "$t0" "$DD/families.progress" <<'PY' || prc=$?
import os, re, sys
root, t0, prog = sys.argv[1], float(sys.argv[2]), sys.argv[3]
logs = [os.path.join(d, f) for d, _, fs in os.walk(root) for f in fs if f == "wrapper.log"]
logs = [p for p in logs if os.path.getmtime(p) >= t0]
if not logs:
    sys.exit(1)
tail = [ln for ln in open(max(logs, key=os.path.getmtime)).read().splitlines() if ln.strip()][-3:]
if not any(re.search(r"budget|deadline", ln) for ln in tail):
    sys.exit(1)
done = sum(1 for d, _, fs in os.walk(root) for f in fs if f == "idx.mmi" or f.endswith(".ok"))
try:
    before = int(open(prog).read().strip())
except (OSError, ValueError):
    before = None
with open(prog, "w") as fh:
    fh.write(f"{done}\n")
sys.exit(3 if before == done else 0)
PY
      case $prc in
        0) say "families: the shard budget is spent (exit $rc); the progress is kept: call again with the same command"; exit 75 ;;
        3) say "families: no progress since the previous call: one step (the index or a single shard) does not fit in the budget of ${RS_SHARD_BUDGET:-480} s; calling again would repeat it. Use smaller shards (export MM2_SHARD_BP, bases per shard) or a larger RS_SHARD_BUDGET"; exit 1 ;;
        *) say "families failed (exit $rc; 124 = the $RLOCK_TIMEOUT s limit): see $A.families.driver.stderr"; exit 1 ;;
      esac
    fi
    rm -f "$DD/families.progress"
    say "families: $(grep -h 'families:' "$A.families.driver.stderr" | head -2 | sed 's/.*families: //' | tr '\n' ';')"
  fi ;;
anno)
  need_depth; arm "$DEPTH"
  need "$P.transcripts.tsv" "truth first"
  light python3 "$HERE/sim_genes.py" control --out "$P" --to "$C0"
  light python3 "$HERE/sim_genes.py" control --out "$P" --to "$C5" --shift 50000 ;;
score)
  need_depth; arm "$DEPTH"; score_gate; score_needs; score_common
  score_one score show "${SC[@]}" --unit node --loci-gff3 "$A.fam.loci.gff3"
  score_common      # the primary run froze the cut if this is the calibration arm: the other variants read it
  score_one score_tx quiet "${SC[@]}" --unit transcript --gtf "$A.families.gtf"
  score_one score_locus quiet "${SC[@]}" --unit locus --gtf "$A.families.gtf"
  score_one score_soto quiet "${SC[@]}" --unit node --loci-gff3 "$A.fam.loci.gff3" --compete soto
  score_one score_raw quiet "${SC[@]}" --unit transcript --gtf "$A.gtf" ;;
score-ctrl)
  need_depth; arm "$DEPTH"
  need "$C0.loci.gff3" "anno first"
  need "$CUTF" "the single-exon cut is not frozen yet: run score $CUT_DEPTH first"
  score_one score_ctrl0 quiet --truth "$P" --single-cut-file "$CUTF" --unit node --loci-gff3 "$C0.loci.gff3"
  score_one score_ctrl0_locus quiet --truth "$P" --single-cut-file "$CUTF" --unit locus --gtf "$C0.gtf"
  score_one score_ctrl50k quiet --truth "$P" --single-cut-file "$CUTF" --unit node --loci-gff3 "$C5.loci.gff3" ;;
score-real)
  arm "$CUT_DEPTH"
  need "$P.transcripts.tsv" "truth first"
  need "$CUTF" "the single-exon cut is not frozen yet: run score $CUT_DEPTH first"
  [ "$DRY" = 1 ] || mkdir -p "$W/real"
  real_one() { # label args...
    local name=$1; shift
    light python3 "$HERE/score_soto_nodes.py" --truth "$P" --out "$W/real/$name" --real-library --single-cut-file "$CUTF" "$@" | logq "$W/real/$name.log"
  }
  rd=()
  if [ "$DRY" = 1 ] || [ -s "$A119B_READS" ]; then rd=(--reads-per-gene "$A119B_READS" --reads-column prim_on_exons); fi
  if [ "$DRY" = 1 ] || [ -s "$REAL_OLD" ]; then
    real_one old_tx --unit transcript --gtf "$REAL_OLD" ${rd[@]+"${rd[@]}"}; real_one old_locus --unit locus --gtf "$REAL_OLD" ${rd[@]+"${rd[@]}"}
  else say "no old A119b GTF at $REAL_OLD: skipped"; fi
  if [ "$DRY" = 1 ] || [ -s "$REAL_F1V2" ]; then
    real_one f1v2_tx --unit transcript --gtf "$REAL_F1V2" ${rd[@]+"${rd[@]}"}; real_one f1v2_locus --unit locus --gtf "$REAL_F1V2" ${rd[@]+"${rd[@]}"}
  else say "no A119b f1v2 GTF at $REAL_F1V2: skipped"; fi
  if [ "$DRY" = 1 ] || [ -s "$REAL_TESTIS.families.gtf" ]; then
    real_one testis_tx --unit transcript --gtf "$REAL_TESTIS.families.gtf" --assemble-log "$REAL_TESTIS.assemble.log"
    real_one testis_locus --unit locus --gtf "$REAL_TESTIS.families.gtf" --assemble-log "$REAL_TESTIS.assemble.log"
    real_one testis_node --unit node --loci-gff3 "$REAL_TESTIS.fam.loci.gff3" --assemble-log "$REAL_TESTIS.assemble.log"
  else say "no testis run at $REAL_TESTIS: skipped"; fi ;;
report)
  need_depth; arm "$DEPTH"; score_gate; score_needs
  need "$CUTF" "the single-exon cut is not frozen yet: run score $CUT_DEPTH first (report never freezes it)"
  score_common
  cmps=()
  for pair in "ctrl0 node=$DD/score_ctrl0" "ctrl50k node=$DD/score_ctrl50k" "ctrl0 locus=$DD/score_ctrl0_locus" "old A119b tx=$W/real/old_tx" "old A119b locus=$W/real/old_locus" \
              "f1v2 A119b tx=$W/real/f1v2_tx" "f1v2 A119b locus=$W/real/f1v2_locus" "testis node=$W/real/testis_node" "testis tx=$W/real/testis_tx" "testis locus=$W/real/testis_locus"; do
    if [ "$DRY" = 1 ] || [ -s "${pair#*=}.summary.json" ]; then cmps+=(--compare "${pair%%=*}=${pair#*=}.summary.json"); fi
  done
  for d in $DEPTHS; do
    if [ "$d" != "$DEPTH" ] && { [ "$DRY" = 1 ] || [ -s "$W/d$d/score.summary.json" ]; }; then cmps+=(--compare "depth $d node=$W/d$d/score.summary.json"); fi
  done
  score_one report show "${SC[@]}" --unit node --loci-gff3 "$A.fam.loci.gff3" ${cmps[@]+"${cmps[@]}"} ;;
gates-sim)
  need "$PP.transcripts.tsv" "truth first"
  [ "$DRY" = 1 ] || mkdir -p "$W/g3"
  light env PYTHONHASHSEED=1 python3 "$HERE/sim_genes.py" --out "$W/g3/tables" "${sim_args[@]}" --reads-per-tx "$POOL" --no-reads | logq "$W/g3/tables.log"
  for s in 0 1; do light env PYTHONHASHSEED=$s python3 "$HERE/sim_genes.py" --out "$W/g3/sub$s" "${sim_args[@]}" --reads-per-tx "$POOL" --first 120 | logq "$W/g3/sub$s.log"; done
  rc=0
  light python3 "$HERE/sim_genes.py" compare --json "$W/g3/cmp_tables.json" --only .transcripts.tsv,.targets.tsv,.windows.tsv,.skipped.tsv,.identical_seq.tsv "$PP" "$W/g3/tables" || rc=1
  light python3 "$HERE/sim_genes.py" compare --json "$W/g3/cmp_subset.json" "$W/g3/sub0" "$W/g3/sub1" || rc=1
  if [ "$DRY" != 1 ]; then
    python3 - "$W" <<'PY'
import json, sys
w = sys.argv[1]
parts = [json.load(open(f"{w}/g3/cmp_tables.json")), json.load(open(f"{w}/g3/cmp_subset.json"))]
json.dump(dict(ok=all(p["ok"] for p in parts), lines=[ln for p in parts for ln in p["lines"]]), open(f"{w}/g3_sim.json", "w"), indent=1)
PY
  fi
  exit $rc ;;
gates)
  need_depth; arm "$DEPTH"
  need "$P.fq" "truth first"
  need "$A.fam.loci.gff3" "families first (the scorer half of G3 re-runs the scorer)"
  need "$W/g3_sim.json" "gates-sim first (the simulation half of G3)"
  need "$DD/score_ctrl0.summary.json" "score-ctrl first (G6 and G7 read the controls)"
  # the verdicts are read from the files the steps below write; those files are removed first, so a step that crashes leaves its gate FAIL, never a stale PASS
  if [ "$DRY" != 1 ]; then
    rm -f "$DD/gates.recheck.json"
    for s in 0 1; do for x in per_gene.tsv by_family.tsv summary.json report.txt; do rm -f "$DD/g3_s$s.$x"; done; done
  fi
  light python3 "$HERE/sim_genes.py" gates --out "$P" --genome "$FA" --genes "$GENES" --json "$DD/gates.recheck.json" | logq "$DD/g12.recheck.log" || true
  score_gate; score_common
  for s in 0 1; do
    light env PYTHONHASHSEED=$s python3 "$HERE/score_soto_nodes.py" --out "$DD/g3_s$s" "${SC[@]}" --unit node --loci-gff3 "$A.fam.loci.gff3" | logq "$DD/g3_s$s.log" || true
  done
  if [ "$DRY" = 1 ]; then exit 0; fi
  NOW_STAMP=$(mktemp "$TMPDIR/stamp.XXXXXX"); stamp > "$NOW_STAMP"
  rcp=0
  python3 - "$W" "$DD" "$P" "$A" "$BIN" "$REPO" "$NOW_STAMP" "${RS_SRC_OK:-}" "$DEPTH" "$POOL" <<'PY' || rcp=$?
import json, os, subprocess, sys
W, D, P, A, BIN, REPO, NOW, SRC_OK, DEPTH, POOL = sys.argv[1:11]
res = {}
try:
    g12 = json.load(open(D + "/gates.recheck.json"))
    res["G1"] = (g12["G1"]["ok"], g12["G1"]["lines"])
    res["G2"] = (g12["G2"]["ok"], g12["G2"]["lines"])
except (OSError, ValueError, KeyError):
    res["G1"] = (False, ["G1 FAIL: the recheck of the simulation wrote no result (see g12.recheck.log)"])
    res["G2"] = (False, ["G2 FAIL: the recheck of the simulation wrote no result (see g12.recheck.log)"])
sim = json.load(open(W + "/g3_sim.json"))
try:
    same = all(open(f"{D}/g3_s0.{x}", "rb").read() == open(f"{D}/g3_s1.{x}", "rb").read() for x in ("per_gene.tsv", "by_family.tsv", "summary.json", "report.txt"))
except OSError:
    same = False
ok3 = bool(sim["ok"]) and same
res["G3"] = (ok3, ["G3 PASS: the simulation (tables in full, a 120-gene FASTQ) and the scorer's tables are byte-identical under PYTHONHASHSEED 0 and 1" if ok3 else
                   "G3 FAIL: " + ("; ".join(sim["lines"]) if not sim["ok"] else "") + ("" if same else " the scorer's tables differ between hash seeds 0 and 1")])
g4 = json.load(open(P + ".g4.json"))
res["G4"] = (g4["ok"], g4["lines"][:1])
# G5: the binaries are those of the assemble run, the scripts those of the truth stage, and no commit touched src/ after the oldest binary was built (except the whitelisted ones)
def stamps(path):
    return dict(ln.rstrip("\n").split("\t", 1) for ln in open(path) if "\t" in ln)
rec, first, now = stamps(A + ".run.log"), stamps(W + "/truth.run.log"), stamps(NOW)
bins = ["copy_assign", "mcl_families", "as_table", "family_score"]
scripts = ["sim_genes.py", "sim_windows.py", "score_soto_nodes.py", "map_reads.py", "run_soto_pilot.sh", "rustle_pipeline.sh", "rlock.sh", "mm2_shard.sh"]
fails = [f"G5 FAIL binaries: {b} changed since the assemble run ({rec.get(b, 'unrecorded')[:8]} -> {now[b][:8]})" for b in bins if rec.get(b) != now[b]]
fails += [f"G5 FAIL scripts: {f} changed since the truth stage ({first.get(f, 'unrecorded')[:8]} -> {now[f][:8]})" for f in scripts if first.get(f) != now[f]]
mt = [os.path.getmtime(f"{BIN}/{b}") for b in bins if os.path.exists(f"{BIN}/{b}")]
late, white = [], []
try:
    log = subprocess.run(["git", "-C", REPO, "log", "--format=%h %ct", "--", "src", "Cargo.toml", "Cargo.lock"], capture_output=True, text=True, check=True).stdout.split("\n")
    for ln in log:
        if ln.strip() and mt and int(ln.split()[1]) > min(mt):
            (white if any(ln.split()[0].startswith(w) or w.startswith(ln.split()[0]) for w in SRC_OK.split() if w) else late).append(ln.split()[0])
except (subprocess.CalledProcessError, FileNotFoundError):
    fails.append("G5 FAIL binaries: the repository history could not be read, so the source commits newer than the binaries could not be listed")
if late:
    fails.append(f"G5 FAIL binaries: commit(s) {' '.join(late)} touched src/ or Cargo.* after the oldest binary was built; RS_SRC_OK=<commit> whitelists a commit known not to change the binaries")
res["G5"] = (not fails, fails or [f"G5 PASS binaries: {len(bins)} binaries identical to the assemble-run stamps and {len(scripts)} scripts to the truth-stage stamps; head {now.get('head')} (assemble run {rec.get('head')})"
                                  + (f"; whitelisted source commits: {' '.join(white)}" if white else "; no source commit newer than the oldest binary")])
# G6, G7: the controls. Perfect nodes cannot show fused excess, a wrong extent or a gene without a node; a perfect annotation in the wrong place must recover fewer genes than in the right one
try:
    c0, c5 = json.load(open(D + "/score_ctrl0.summary.json")), json.load(open(D + "/score_ctrl50k.summary.json"))
    bad = {k: c0["classes"][k] for k in ("fused_excess", "extent_wrong", "no_locus") if c0["classes"][k]}
    res["G6"] = (not bad, ["G6 PASS positive control: the annotation scored as nodes has no fused excess, no wrong extent and no gene without a node"] if not bad else
                 [f"G6 FAIL positive control: the annotation scored as nodes loses genes it cannot lose: {bad}"])
    ok7 = c5["classes"]["node_ok"] < c0["classes"]["node_ok"]
    res["G7"] = (ok7, ["G7 PASS negative control: the annotation shifted by 50 kb recovers fewer genes than the annotation in place"] if ok7 else
                 ["G7 FAIL negative control: the annotation shifted by 50 kb does not recover fewer genes than the annotation in place"])
except (OSError, ValueError, KeyError):
    res["G6"] = (False, ["G6 FAIL: the control summaries are missing or unreadable"])
    res["G7"] = (False, ["G7 FAIL: the control summaries are missing or unreadable"])
order = ["G1", "G2", "G3", "G4", "G5", "G6", "G7"]
if DEPTH != POOL:
    try:
        ds = json.load(open(f"{W}/depth_subset_{DEPTH}.json"))
        res["G-depth"] = (ds["ok"], ds["lines"][:1] if ds["ok"] else ds["lines"])
    except (OSError, ValueError, KeyError):
        res["G-depth"] = (False, ["G-depth FAIL: the depth check wrote no result"])
    order.append("G-depth")
valid = all(ok for ok, _ in res.values())
for g in order:
    ok, lines = res[g]
    print("\n".join(lines if not ok else lines[:1]))
print("VALID" if valid else "INVALID")
json.dump(dict(valid=valid, gates={g: dict(ok=ok, lines=lines) for g, (ok, lines) in res.items()}), open(D + "/gates.json", "w"), indent=1)
sys.exit(0 if valid else 1)
PY
  rm -f "$NOW_STAMP"
  exit $rcp ;;
cleanup)
  need "$PP.bam" "map first"
  if [ "$DRY" != 1 ]; then
    [ -s "$PD/map.parts" ] || { say "no mapping parts recorded: nothing to clean"; exit 0; }
    n=$(( $(wc -l < "$PP.fq") / 4 ))
    got=$(light samtools view -c -F 2304 "$PP.bam")
    [ "$got" = "$n" ] || { echo "[soto_pilot] the pool's BAM holds $got primary-or-unmapped records for $n reads: not complete, the part files are kept" >&2; exit 1; }
    PARTS=$(cat "$PD/map.parts")
  else PARTS=${RS_PARTS:-1}; fi
  for i in $(seq 0 $((PARTS - 1))); do
    for ext in .fq .bam .bam.tmp .mm2.log; do
      f="$PP.part$i$ext"
      if [ "$DRY" = 1 ]; then echo "DRY rm -f $f"; else rm -f "$f"; fi
    done
  done
  say "the part files of the pool's mapping are deleted (exact names); the BAM and the FASTQ stay" ;;
status)
  echo "work dir: $W   depths: $DEPTHS   pool: $POOL   calibration depth: $CUT_DEPTH"
  for f in truth.run.log depth_subset.json single_cut.json g3_sim.json baseline_ids.txt; do
    if [ -e "$W/$f" ]; then echo "  $f	$(stat -c '%s bytes, %y' "$W/$f" | cut -c1-40)"; else echo "  $f	-"; fi
  done
  for d in $DEPTHS; do
    echo "  d$d:"
    for f in reads.fq reads.transcripts.tsv reads.targets.tsv reads.gates.json reads.bam reads.g4.json asm.gtf asm.families.gtf asm.fam.loci.gff3 asm.fam.clusters.tsv ctrl0.loci.gff3 score.summary.json score_ctrl0.summary.json gates.json; do
      if [ -e "$W/d$d/$f" ]; then echo "    $f	$(stat -c '%s bytes, %y' "$W/d$d/$f" | cut -c1-40)"; else echo "    $f	-"; fi
    done
  done
  n=0; np=0; if [ -s "$PD/map.parts" ]; then np=$(cat "$PD/map.parts"); for i in $(seq 0 $((np - 1))); do if [ -s "$PP.part$i.bam" ]; then n=$((n + 1)); fi; done; fi
  echo "  parts mapped: $n / $np (pool d$POOL)" ;;
*) echo "usage: run_soto_pilot.sh [--dry-run] truth|truthqc|map|cut|bamqc D|assemble D|families D|anno D|score D|score-ctrl D|score-real|report D|gates-sim|gates D|cleanup|status|help" >&2; exit 2 ;;
esac
