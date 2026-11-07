#!/bin/bash
# rlock.sh — run one command under the machine's job locks (written for a 5-core, 26 GB workstation; RLOCK_DIR = lock dir).
#
#   tools/rlock.sh heavy <cmd...>   one at a time: assemblies, minimap2, cargo build/test, anything > 2 GB or > 3 min
#   tools/rlock.sh light <cmd...>   up to LIGHT_SLOTS (default 2) at a time, ALSO while a heavy job runs: gffcompare,
#                                   python scoring, small assemblies (< 2 GB RSS, < 3 min) — never a build or minimap2
#
# Both classes wait up to RLOCK_WAIT s (default 900) and run the command under `timeout RLOCK_TIMEOUT` (default 600),
# in the foreground. Light slots are separate lock files, so light jobs never block a heavy one and vice versa; the
# only shared resource is the 5 cores, which is why LIGHT_SLOTS stays at 2. Exit code = the command's (124 = timeout,
# 73 = the lock wait ran out and the command never ran; not 75, which mm2_shard.sh and the bench drivers use for "run again"). shortcut: when every light slot is busy the wait is on slot 1
# only, so a free slot 2 is not taken; fine at 2 slots, use a semaphore if LIGHT_SLOTS grows.
set -u
class=${1:?heavy|light}; shift
D=${RLOCK_DIR:-${TMPDIR:-/tmp}}; W=${RLOCK_WAIT:-900}; T=${RLOCK_TIMEOUT:-600}; N=${LIGHT_SLOTS:-2}
case $class in
  heavy) flock -E 73 -w "$W" "$D/rustle_heavy.lock" timeout "$T" "$@"; rc=$?; lock="$D/rustle_heavy.lock" ;;
  light)
    # hold a slot atomically: take the first free slot now, else wait on slot 1
    for i in $(seq 1 "$N"); do
      exec 9>"$D/rustle_light$i.lock"
      if flock -n 9; then exec timeout "$T" "$@"; fi
      exec 9>&-
    done
    flock -E 73 -w "$W" "$D/rustle_light1.lock" timeout "$T" "$@"; rc=$?; lock="$D/rustle_light1.lock" ;;
  *) echo "rlock.sh: class must be heavy or light" >&2; exit 2 ;;
esac
# a lock wait that ran out must not read as the command's own failure
[ "$rc" = 73 ] && echo "rlock.sh: waited ${W} s for $lock and gave up; the command did not run (RLOCK_WAIT)" >&2
exit "$rc"
