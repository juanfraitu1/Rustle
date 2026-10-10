#!/usr/bin/env bash
# tools/mm2_shard.sh -- sharded, resumable, cached minimap2 all-vs-all (drop-in RUSTLE_MINIMAP2 wrapper).
#
# WHY. mcl_families --from-gtf and gw_family_catalog run ONE all-vs-all, `minimap2 <flags> X.fa X.fa`, which
# genome-wide takes hours in a single process. This wrapper builds the index of X.fa once (`-d`, same flags),
# maps contiguous query shards against that one index (each shard only if its PAF is missing), and prints
# the shards' PAFs concatenated IN INPUT ORDER. minimap2 writes records in query order and every per-query
# decision (-X/--dual=no/-D by name, -N/-p/--secondary, mid_occ from the index, tie-break hash from qname and
# qlen) is independent of the other queries, so the concatenation equals a single run byte for byte
# (checked with cmp on chr16, see docs below). Two conditions are enforced because they WOULD break that:
#   * one index part: minimap2 splits an index every -I bases (default 8G) and a multi-part index writes
#     part-major output, so the wrapper refuses (exit 2) when the target has >1 part;
#   * the index is SHARED: never --split-prefix, never a per-shard index.
#
# USAGE (wrapper mode; argv exactly as minimap2 would get it):
#   RUSTLE_MINIMAP2=tools/mm2_shard.sh gw_family_catalog ...          # or mcl_families --from-gtf ...
#   tools/mm2_shard.sh <minimap2 flags...> TARGET.fa QUERY.fa > out.paf
# Only an all-vs-all (TARGET and QUERY the same file, plain FASTA, >= MM2_SHARD_MIN_BYTES) is sharded;
# everything else (--version, -d, -o, -a/SAM, --split-prefix, target != query, .mmi/.gz/FASTQ inputs,
# small inputs, MM2_SHARD_DISABLE=1) is exec'd to the real minimap2 unchanged.
#
# CLI mode (for Python / shell drivers; always shards, target may differ from query):
#   tools/mm2_shard.sh paf OUT.paf <minimap2 flags...> TARGET.fa QUERY.fa   # writes OUT.paf atomically
#   tools/mm2_shard.sh guided BODIES.fa OUT.paf [THREADS]                    # = paf OUT.paf -x asm20 -c --eqx -P -t THREADS BODIES.fa BODIES.fa
#   tools/mm2_shard.sh status <minimap2 flags...> TARGET.fa QUERY.fa        # progress; exit 0 complete / 75 not yet / 3 no layout
#   tools/mm2_shard.sh help
# Python: loop `subprocess.run(["flock", LOCK, "tools/mm2_shard.sh", "guided", bodies, paf, "4"],
#         env={**os.environ, "MM2_SHARD_BUDGET_S": "540"})` until returncode 0; 75 means "call again".
#
# ENVIRONMENT
#   MM2_SHARD_MINIMAP2   the real minimap2 (default: `minimap2` on PATH; must not resolve to this script)
#   MM2_SHARD_DIR        cache root (default /mnt/linuxdisk/tmp/mm2_shard_cache)
#   MM2_SHARD_READS      records per shard            } first one set wins; default MM2_SHARD_BP=10000000.
#   MM2_SHARD_BP         bases per shard (a shard     } Changing the shard size starts a new layout (old
#                        closes once it reaches BP)   } shards are kept, not reused).
#   MM2_SHARD_N          number of shards (equal-bp contiguous split)
#   MM2_SHARD_BUDGET_S   wall budget in seconds from this call's start; MM2_SHARD_DEADLINE = absolute epoch s.
#                        The earlier of the two wins. A step that cannot finish in time is not started (when
#                        earlier shards give a rate) or is killed at the deadline; its .tmp is removed.
#   MM2_SHARD_MIN_BYTES  wrapper mode: all-vs-alls smaller than this are exec'd unchanged (default 4000000)
#   MM2_SHARD_DISABLE=1  wrapper mode: always exec the real minimap2
#   MM2_SHARD_ON_BUDGET  auto|exit|kill-parent (wrapper mode only; default auto). gw_family_catalog,
#                        copy_assign and missing_copy_flag treat a non-zero minimap2 exit as an EMPTY edge set
#                        and carry on (denovo_pipeline.rs "a non-zero minimap2 exit is still SILENTLY an empty
#                        edge set"), which would write a wrong catalog. `auto` therefore sends SIGTERM to the
#                        parent when the parent's executable is one of those three; mcl_families already
#                        fails on a non-zero exit and is left alone.
#
# EXIT CODES: 0 done (PAF on stdout / OUT.paf written) · 75 budget exhausted, re-run the same command to resume
# (progress is kept) · 76 one step (index or a single shard) did not fit in a whole budget: use smaller shards
# or a larger budget · 2 usage / refused input · 70 minimap2 or verification failure (see the key dir's logs).
#
# LAYOUT: $MM2_SHARD_DIR/<md5(target)>/a<sha1(flags without -t, minimap2 version)>/
#            args.txt  idx.mmi  idx.log
#            q<md5(query)>/t<threads>.<spec>/  layout.tsv  shard_0000.{fa,names,paf,ok,log}  wrapper.log
# A shard's .paf appears only by rename after minimap2 exited 0, its log ends with "[M::main] Real time",
# every record has >= 12 columns and its query names are the shard's own names in input order.
set -uo pipefail

SELF=$(readlink -f "${BASH_SOURCE[0]}")
T0=$(date +%s)
MM2=${MM2_SHARD_MINIMAP2:-minimap2}
ROOT=${MM2_SHARD_DIR:-/mnt/linuxdisk/tmp/mm2_shard_cache}
EX_BUDGET=75; EX_TOOBIG=76; EX_USAGE=2; EX_FAIL=70

log() { echo "[mm2_shard] $*" >&2; [ -n "${WLOG:-}" ] && echo "$(date '+%F %T') $*" >> "$WLOG"; return 0; }

resolve_mm2() {
  local p
  p=$(command -v "$MM2" 2>/dev/null) || { echo "[mm2_shard] real minimap2 '$MM2' not found" >&2; exit $EX_USAGE; }
  p=$(readlink -f "$p")
  [ "$p" = "$SELF" ] && { echo "[mm2_shard] MM2_SHARD_MINIMAP2 resolves to this wrapper (recursion)" >&2; exit $EX_USAGE; }
  MM2=$p
}

passthrough() { resolve_mm2; exec "$MM2" "$@"; }

# ---------------------------------------------------------------- deadline
DEADLINE=""
if [ -n "${MM2_SHARD_BUDGET_S:-}" ]; then DEADLINE=$(( T0 + MM2_SHARD_BUDGET_S )); fi
if [ -n "${MM2_SHARD_DEADLINE:-}" ]; then
  if [ -z "$DEADLINE" ] || [ "$MM2_SHARD_DEADLINE" -lt "$DEADLINE" ]; then DEADLINE=$MM2_SHARD_DEADLINE; fi
fi
remaining() { if [ -z "$DEADLINE" ]; then echo 1000000000; else echo $(( DEADLINE - $(date +%s) )); fi; }

# ---------------------------------------------------------------- input classification
is_plain_fasta() {  # regular, non-empty, not gz, not an .mmi, first byte '>'
  local f=$1 c
  [ -f "$f" ] && [ -s "$f" ] || return 1
  case "$f" in *.gz|*.mmi) return 1;; esac
  c=$(head -c 1 "$f")
  [ "$c" = ">" ]
}

# split "$@" (minimap2 argv) into FLAGS (all but the last two), TARGET, QUERY, THREADS, KEYFLAGS (FLAGS minus -t)
parse_argv() {
  local n=$#
  [ "$n" -ge 2 ] || return 1
  local a=("$@")
  TARGET=${a[$((n-2))]}; QUERY=${a[$((n-1))]}
  FLAGS=("${a[@]:0:$((n-2))}")
  THREADS=3; KEYFLAGS=(); IBASES=8000000000
  local i=0 m=${#FLAGS[@]} f
  while [ $i -lt $m ]; do
    f=${FLAGS[$i]}
    case "$f" in
      -t) THREADS=${FLAGS[$((i+1))]:-3}; i=$((i+2)); continue;;
      -t[0-9]*) THREADS=${f#-t}; i=$((i+1)); continue;;
      -I) IBASES=$(num_bases "${FLAGS[$((i+1))]:-8G}");;
      -I[0-9]*) IBASES=$(num_bases "${f#-I}");;
    esac
    KEYFLAGS+=("$f"); i=$((i+1))
  done
  return 0
}

num_bases() {  # minimap2 mm_parse_num: k/K 1e3, m/M 1e6, g/G 1e9
  awk -v s="$1" 'BEGIN{ n=s+0; u=substr(s,length(s)); m=1;
    if(u=="k"||u=="K")m=1e3; else if(u=="m"||u=="M")m=1e6; else if(u=="g"||u=="G")m=1e9; printf "%.0f", n*m }'
}

refused_flag() {  # flags whose output cannot be sharded/concatenated as-is
  local f
  for f in "$@"; do
    case "$f" in
      --version|-V|--help|-h|-d|-o|-a|--split-prefix|--split-prefix=*) return 0;;
    esac
  done
  return 1
}

md5of() { md5sum "$1" | cut -c1-32; }

# ---------------------------------------------------------------- the engine
# engine MODE(wrap|cli|status) ; uses TARGET QUERY FLAGS KEYFLAGS THREADS IBASES ; on success cats shards to fd 1
engine() {
  local mode=$1
  resolve_mm2
  local ver md5t md5q argsha
  ver=$("$MM2" --version 2>/dev/null)
  md5t=$(md5of "$TARGET")
  if [ "$(readlink -f "$TARGET")" = "$(readlink -f "$QUERY")" ]; then md5q=$md5t; else md5q=$(md5of "$QUERY"); fi
  argsha=$( { printf '%s\n' "${KEYFLAGS[@]}"; printf 'minimap2 %s\n' "$ver"; } | sha1sum | cut -c1-16)
  local spec
  if [ -n "${MM2_SHARD_READS:-}" ]; then spec="R${MM2_SHARD_READS}"
  elif [ -n "${MM2_SHARD_BP:-}" ]; then spec="B${MM2_SHARD_BP}"
  elif [ -n "${MM2_SHARD_N:-}" ]; then spec="N${MM2_SHARD_N}"
  else spec="B10000000"; fi
  KEYDIR=$ROOT/$md5t/a$argsha
  local qdir=$KEYDIR/q$md5q/t$THREADS.$spec
  if [ "$mode" = status ]; then
    status_report "$qdir"; return $?
  fi
  mkdir -p "$qdir" || { echo "[mm2_shard] cannot create $qdir" >&2; return $EX_FAIL; }
  WLOG=$qdir/wrapper.log
  # one wrapper per key at a time (never the shared heavy lock: the caller holds it)
  exec 9>"$KEYDIR/.lock"
  flock 9
  [ -f "$KEYDIR/args.txt" ] || { printf '%s\n' "${KEYFLAGS[@]}" > "$KEYDIR/args.txt.tmp" && printf 'minimap2 %s\n' "$ver" >> "$KEYDIR/args.txt.tmp" && mv "$KEYDIR/args.txt.tmp" "$KEYDIR/args.txt"; }
  log "start mode=$mode target=$TARGET ($md5t) query=$QUERY ($md5q) flags='${FLAGS[*]}' dir=$qdir deadline_in=$(remaining)s"

  # ---- layout (split the query once)
  if [ ! -f "$qdir/layout.tsv" ]; then
    local tmp=$qdir/.layout.tmp.$$
    rm -rf "$tmp"; mkdir -p "$tmp"
    local R=${MM2_SHARD_READS:-0} BP=${MM2_SHARD_BP:-0} N=${MM2_SHARD_N:-0}
    if [ -z "${MM2_SHARD_READS:-}${MM2_SHARD_BP:-}${MM2_SHARD_N:-}" ]; then BP=10000000; fi
    [ -n "${MM2_SHARD_READS:-}" ] && { BP=0; N=0; }
    [ -z "${MM2_SHARD_READS:-}" ] && [ -n "${MM2_SHARD_BP:-}" ] && N=0
    # pass 1: record lengths (bases, newlines excluded) -> $tmp/recs (name \t bp)
    LC_ALL=C awk -v OUT="$tmp/recs" '
      /^>/ { if (n) print name "\t" bp > OUT; n++; name=substr($1,2); bp=0; next }
      { bp += length($0) }
      END { if (n) print name "\t" bp > OUT }' "$QUERY" || { log "layout pass 1 failed"; return $EX_FAIL; }
    local total nrec ndup
    read -r nrec total ndup < <(LC_ALL=C awk -F'\t' '{n++; t+=$2; if (seen[$1]++) d++} END{printf "%d %.0f %d\n", n, t, d}' "$tmp/recs")
    # shard assignment per record -> $tmp/assign (shard index, one per record)
    LC_ALL=C awk -F'\t' -v R="$R" -v BP="$BP" -v N="$N" -v T="$total" '
      BEGIN { s=0; cnt=0; acc=0; cum=0; kcur=0 }
      {
        if (N > 0) { k=int(cum * N / (T > 0 ? T : 1)); if (k >= N) k=N-1; if (k > kcur && cnt > 0) { s++; kcur=k; cnt=0; acc=0 } }
        else if (R > 0) { if (cnt >= R) { s++; cnt=0; acc=0 } }
        else { if (acc >= BP && cnt > 0) { s++; cnt=0; acc=0 } }
        print s; cnt++; acc+=$2; cum+=$2
      }' "$tmp/recs" > "$tmp/assign"
    # pass 2: write the shards verbatim (header and sequence lines byte-identical to the query)
    LC_ALL=C awk -v A="$tmp/assign" -v D="$tmp" '
      /^>/ { if ((getline s < A) <= 0) { print "assign underflow" > "/dev/stderr"; exit 3 }
             if (s != cur) { if (fn) { close(fn); close(nn) } cur=s; fn=sprintf("%s/shard_%04d.fa", D, s); nn=sprintf("%s/shard_%04d.names", D, s) }
             print substr($1,2) > nn }
      { print > fn }
      BEGIN { cur=-1 }' "$QUERY" || { log "layout pass 2 failed"; rm -rf "$tmp"; return $EX_FAIL; }
    # layout.tsv: shard, records, bases, first, last
    LC_ALL=C paste "$tmp/assign" "$tmp/recs" | awk -F'\t' '
      { if (!($1 in n)) { first[$1]=$2; order[m++]=$1 } n[$1]++; b[$1]+=$3; last[$1]=$2 }
      END { print "#shard\trecords\tbases\tfirst\tlast"; for (i=0;i<m;i++){s=order[i]; printf "%d\t%d\t%.0f\t%s\t%s\n", s, n[s], b[s], first[s], last[s]} }' > "$tmp/layout.tsv.part"
    # sanity: the shards reassemble to the query exactly
    local nshard; nshard=$(( $(wc -l < "$tmp/layout.tsv.part") - 1 ))
    if ! cat "$tmp"/shard_*.fa | cmp -s - "$QUERY"; then log "layout: shards do not reassemble to the query"; rm -rf "$tmp"; return $EX_FAIL; fi
    { echo "# query=$QUERY md5=$md5q records=$nrec bases=$total duplicate_names=$ndup spec=$spec"; cat "$tmp/layout.tsv.part"; } > "$tmp/layout.tsv.new"
    mv "$tmp"/shard_* "$qdir"/ && mv "$tmp/layout.tsv.new" "$qdir/layout.tsv"
    rm -rf "$tmp"
    log "layout: $nrec records, $total bases -> $nshard shards ($spec)"
  fi
  local nshard; nshard=$(grep -vc '^#' "$qdir/layout.tsv")
  local dupnames; dupnames=$(head -1 "$qdir/layout.tsv" | sed -n 's/.*duplicate_names=\([0-9]*\).*/\1/p')

  # ---- target bases vs -I (one index part), from the layout when target == query
  local tbases
  if [ "$md5t" = "$md5q" ]; then tbases=$(head -1 "$qdir/layout.tsv" | sed -n 's/.* bases=\([0-9]*\).*/\1/p')
  else tbases=$(LC_ALL=C awk '!/^>/{b+=length($0)} END{printf "%.0f", b}' "$TARGET"); fi
  if [ "$tbases" -ge "$IBASES" ]; then
    log "REFUSED: target has $tbases bases >= -I $IBASES: a multi-part index writes part-major output, shards would not concatenate to a single run"
    return $EX_USAGE
  fi

  # ---- the shared index
  if [ ! -f "$KEYDIR/idx.mmi" ]; then
    local rem; rem=$(remaining)
    if [ "$rem" -le 0 ]; then log "budget exhausted before the index build"; return $EX_BUDGET; fi
    log "index: building $KEYDIR/idx.mmi (-t $THREADS, <= ${rem}s)"
    local t1; t1=$(date +%s)
    run_bounded "$rem" "$KEYDIR/idx.log" /dev/null "$MM2" "${KEYFLAGS[@]}" -t "$THREADS" -d "$KEYDIR/idx.mmi.tmp.$$" "$TARGET"
    local rc=$?
    if [ $rc -eq 124 ] || [ $rc -eq 137 ]; then
      rm -f "$KEYDIR/idx.mmi.tmp.$$"; log "index build hit the deadline"
      whole_budget "$rem" && return $EX_TOOBIG
      return $EX_BUDGET
    fi
    if [ $rc -ne 0 ]; then rm -f "$KEYDIR/idx.mmi.tmp.$$"; log "index build failed (rc $rc, $KEYDIR/idx.log)"; return $EX_FAIL; fi
    local parts; parts=$(grep -c 'loaded/built the index for' "$KEYDIR/idx.log")
    if [ "$parts" != 1 ]; then rm -f "$KEYDIR/idx.mmi.tmp.$$"; log "REFUSED: index has $parts parts (need exactly 1)"; return $EX_USAGE; fi
    mv "$KEYDIR/idx.mmi.tmp.$$" "$KEYDIR/idx.mmi"
    log "index: built in $(( $(date +%s) - t1 ))s"
  fi

  # ---- shards in order
  local i s paf base
  for (( i=0; i<nshard; i++ )); do
    base=$(printf '%s/shard_%04d' "$qdir" "$i")
    paf=$base.paf
    if [ -f "$paf" ] && [ -f "$base.ok" ]; then continue; fi
    local bp; bp=$(awk -F'\t' -v i="$i" '!/^#/ && $1==i {print $3}' "$qdir/layout.tsv")
    if ! parent_alive; then log "parent (pid $PPID) is gone: stopping before shard $i"; return $EX_BUDGET; fi
    local rem; rem=$(remaining)
    local est; est=$(estimate "$qdir" "$bp")
    if [ "$rem" -le 0 ] || { [ -n "$est" ] && [ "$est" -gt "$rem" ]; }; then
      log "budget: shard $i/$nshard not started (remaining ${rem}s, estimate ${est:-?}s); $(done_count "$qdir" "$nshard")/$nshard shards done"
      [ "$rem" -gt 0 ] && whole_budget "$rem" && return $EX_TOOBIG
      return $EX_BUDGET
    fi
    local t1; t1=$(date +%s)
    run_bounded "$rem" "$base.log" "$paf.tmp" "$MM2" "${KEYFLAGS[@]}" -t "$THREADS" "$KEYDIR/idx.mmi" "$base.fa"
    local rc=$?
    local dt=$(( $(date +%s) - t1 ))
    if [ $rc -eq 124 ] || [ $rc -eq 137 ]; then
      rm -f "$paf.tmp"
      log "shard $i/$nshard hit the deadline after ${dt}s; $(done_count "$qdir" "$nshard")/$nshard done"
      whole_budget "$rem" && return $EX_TOOBIG
      return $EX_BUDGET
    fi
    if [ $rc -ne 0 ]; then rm -f "$paf.tmp"; log "shard $i: minimap2 failed (rc $rc, $base.log)"; return $EX_FAIL; fi
    if ! grep -q '\[M::main\] Real time' "$base.log"; then rm -f "$paf.tmp"; log "shard $i: log has no completion line"; return $EX_FAIL; fi
    local chk
    chk=$(LC_ALL=C awk -F'\t' -v NAMES="$base.names" -v ORD="$([ "${dupnames:-0}" = 0 ] && echo 1 || echo 0)" '
      BEGIN { while ((getline n < NAMES) > 0) if (!(n in idx)) idx[n]=++k; last=0 }
      NF < 12 { print "short record at line " NR; exit 1 }
      { j=idx[$1]; if (!j) { print "foreign query " $1; exit 1 } if (ORD && j < last) { print "query out of order at line " NR; exit 1 } last=j }
      END { }' "$paf.tmp")
    if [ -n "$chk" ]; then rm -f "$paf.tmp"; log "shard $i: verification failed: $chk"; return $EX_FAIL; fi
    local bytes lines; bytes=$(stat -c %s "$paf.tmp"); lines=$(wc -l < "$paf.tmp")
    printf 'bytes\t%s\nlines\t%s\nseconds\t%s\nbases\t%s\n' "$bytes" "$lines" "$dt" "$bp" > "$base.ok"
    mv "$paf.tmp" "$paf"
    log "shard $i/$nshard done: ${dt}s, $lines records"
  done

  # ---- all present: verify sizes, then emit in input order
  for (( i=0; i<nshard; i++ )); do
    base=$(printf '%s/shard_%04d' "$qdir" "$i")
    local want; want=$(awk -F'\t' '$1=="bytes"{print $2}' "$base.ok")
    [ "$(stat -c %s "$base.paf")" = "$want" ] || { log "shard $i: PAF size differs from its .ok record"; rm -f "$base.ok"; return $EX_FAIL; }
  done
  log "complete: $nshard shards; emitting $(( $(date +%s) - T0 ))s after start"
  local files=(); for (( i=0; i<nshard; i++ )); do files+=("$(printf '%s/shard_%04d.paf' "$qdir" "$i")"); done
  cat "${files[@]}"
}

# the caller still exists (our parent pid unchanged); an orphaned wrapper stops instead of mapping for nobody
parent_alive() { local pp; pp=$(awk '{print $4}' /proc/$$/stat 2>/dev/null); [ "${pp:-$PPID}" = "$PPID" ]; }

# did a step that started with $1 seconds left have (nearly) a whole budget? Then a deadline hit means it cannot
# fit in ANY call of this budget (exit 76, not 75), so a resume loop stops instead of spinning.
whole_budget() { [ -n "${MM2_SHARD_BUDGET_S:-}" ] && [ "$1" -ge $(( MM2_SHARD_BUDGET_S * 3 / 4 )) ]; }

done_count() { local n=0 i; for (( i=0; i<$2; i++ )); do [ -f "$(printf '%s/shard_%04d.ok' "$1" "$i")" ] && [ -f "$(printf '%s/shard_%04d.paf' "$1" "$i")" ] && n=$((n+1)); done; echo $n; }

estimate() {  # seconds for a shard of $2 bases, from the finished shards' seconds per base (+20% margin)
  cat "$1"/shard_*.ok 2>/dev/null | awk -F'\t' -v b="$2" '$1=="seconds"{s+=$2} $1=="bases"{x+=$2}
    END{ if (b == 0) print 0; else if (x > 0 && s > 0) printf "%d\n", s / x * b * 1.2 + 5 }'
}

# run_bounded SECONDS LOG OUT cmd... : foreground-equivalent child, killable, bounded by `timeout`
CHILD=""
run_bounded() {
  local secs=$1 lg=$2 out=$3; shift 3
  timeout --kill-after=15 "$secs" "$@" > "$out" 2> "$lg" &
  CHILD=$!
  wait "$CHILD"; local rc=$?
  CHILD=""
  return $rc
}
on_signal() { [ -n "$CHILD" ] && kill "$CHILD" 2>/dev/null; rm -f "${KEYDIR:-/nonexistent}"/idx.mmi.tmp.$$ 2>/dev/null; exit 143; }
trap on_signal TERM INT HUP

status_report() {
  local qdir=$1
  echo "key_dir	$KEYDIR"
  echo "shard_dir	$qdir"
  echo "index	$([ -f "$KEYDIR/idx.mmi" ] && echo built || echo missing)"
  if [ ! -f "$qdir/layout.tsv" ]; then echo "layout	missing"; return 3; fi
  local n d; n=$(grep -vc '^#' "$qdir/layout.tsv"); d=$(done_count "$qdir" "$n")
  echo "shards	$n"
  echo "done	$d"
  local rest; rest=$(awk -F'\t' -v q="$qdir" '!/^#/ { f=sprintf("%s/shard_%04d.ok", q, $1); if ((getline x < f) <= 0) b+=$3; close(f) } END{printf "%.0f", b}' "$qdir/layout.tsv")
  echo "bases_left	$rest"
  echo "est_seconds_left	$(estimate "$qdir" "$rest")"
  [ "$d" = "$n" ] && return 0 || return $EX_BUDGET
}

kill_parent_maybe() {
  local policy=${MM2_SHARD_ON_BUDGET:-auto} exe
  [ "$policy" = exit ] && return 0
  exe=$(basename "$(readlink -f "/proc/$PPID/exe" 2>/dev/null)" 2>/dev/null)
  if [ "$policy" = kill-parent ] || { [ "$policy" = auto ] && case "$exe" in gw_family_catalog|copy_assign|missing_copy_flag) true;; *) false;; esac; }; then
    log "stopping the parent ($exe, pid $PPID) with SIGTERM: it would treat this non-zero exit as an empty edge set. Re-run the same command to resume."
    kill -TERM "$PPID" 2>/dev/null
  fi
}

# ---------------------------------------------------------------- dispatch
# Everything runs inside main(), called from the LAST line: bash has parsed the whole file before any work
# starts, so editing this script while a long run uses it cannot change that run.
main() {
case "${1:-}" in
  help|--mm2-shard-help)
    sed -n '2,/^set -uo/p' "$SELF" | sed '$d' | sed 's/^# \{0,1\}//'; exit 0;;
  paf)
    shift; OUT=${1:?paf OUT.paf <flags> TARGET QUERY}; shift
    parse_argv "$@" || { echo "usage: mm2_shard.sh paf OUT.paf <minimap2 flags> TARGET QUERY" >&2; exit $EX_USAGE; }
    refused_flag "${FLAGS[@]}" && { echo "[mm2_shard] flags not shardable: ${FLAGS[*]}" >&2; exit $EX_USAGE; }
    is_plain_fasta "$TARGET" && is_plain_fasta "$QUERY" || { echo "[mm2_shard] TARGET and QUERY must be plain FASTA files" >&2; exit $EX_USAGE; }
    mkdir -p "$(dirname "$OUT")"
    engine cli > "$OUT.tmp.$$"; rc=$?
    if [ $rc -eq 0 ]; then mv "$OUT.tmp.$$" "$OUT"; else rm -f "$OUT.tmp.$$"; fi
    exit $rc;;
  guided)
    shift; B=${1:?guided BODIES.fa OUT.paf [THREADS]}; OUT=${2:?guided BODIES.fa OUT.paf [THREADS]}; TH=${3:-4}
    exec "$SELF" paf "$OUT" -x asm20 -c --eqx -P -t "$TH" "$B" "$B";;
  status)
    shift
    parse_argv "$@" || { echo "usage: mm2_shard.sh status <minimap2 flags> TARGET QUERY" >&2; exit $EX_USAGE; }
    engine status; exit $?;;
esac

# wrapper mode
[ "${MM2_SHARD_DISABLE:-0}" = 1 ] && passthrough "$@"
parse_argv "$@" || passthrough "$@"
refused_flag "${FLAGS[@]}" && passthrough "$@"
is_plain_fasta "$TARGET" && is_plain_fasta "$QUERY" || passthrough "$@"
[ "$(readlink -f "$TARGET")" = "$(readlink -f "$QUERY")" ] || passthrough "$@"
[ "$(stat -L -c %s "$TARGET")" -ge "${MM2_SHARD_MIN_BYTES:-4000000}" ] || passthrough "$@"
engine wrap; local rc=$?
if [ $rc -ne 0 ]; then kill_parent_maybe; fi
exit $rc
}
main "$@"; exit $?
