#!/bin/bash
# rlock.sh — run one command under the machine's job locks (WSL2, 5 cores, 26 GB; see the crash rule in memory).
#
#   tools/rlock.sh heavy <cmd...>   one at a time: assemblies, minimap2, cargo build/test, anything > 2 GB or > 3 min
#   tools/rlock.sh light <cmd...>   up to LIGHT_SLOTS (default 2) at a time, ALSO while a heavy job runs: gffcompare,
#                                   python scoring, small assemblies (< 2 GB RSS, < 3 min) — never a build or minimap2
#
# Both classes wait up to RLOCK_WAIT s (default 900) and run the command under `timeout RLOCK_TIMEOUT` (default 600),
# in the foreground. Light slots are separate lock files, so light jobs never block a heavy one and vice versa; the
# only shared resource is the 5 cores, which is why LIGHT_SLOTS stays at 2. Exit code = the command's (124 = timeout).
set -u
class=${1:?heavy|light}; shift
D=/mnt/linuxdisk/tmp; W=${RLOCK_WAIT:-900}; T=${RLOCK_TIMEOUT:-600}; N=${LIGHT_SLOTS:-2}
case $class in
  heavy) exec flock -w "$W" $D/rustle_heavy.lock timeout "$T" "$@" ;;
  light)
    # hold a slot atomically: take the first free slot now, else wait on slot 1
    for i in $(seq 1 "$N"); do
      exec 9>$D/rustle_light$i.lock
      if flock -n 9; then exec timeout "$T" "$@"; fi
      exec 9>&-
    done
    exec flock -w "$W" $D/rustle_light1.lock timeout "$T" "$@" ;;
  *) echo "rlock.sh: class must be heavy or light" >&2; exit 2 ;;
esac
