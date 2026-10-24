#!/usr/bin/env bash
# tools/minimap2_gate/build.sh -- build `minimap2-gate`: minimap2 2.30 (r1287) + gate.patch, the exon gate that
# `mcl_families --from-gtf` uses for its all-vs-all (docs/archive/2026-10/merged/O1_FAST_2026-10-21.md).
#
#   tools/minimap2_gate/build.sh [BIN_DIR]        default BIN_DIR = ~/.local/bin; needs git, make, a C compiler, zlib headers
#
# The patch adds `--gate FILE` (`name<TAB>start<TAB>end` per exon, 0-based half-open, in the sequence's own coordinates) and
# `--gate-margin INT`: the chains between a query and a target are aligned at base level only when one of them touches a
# listed interval on both sides (+-margin), and `--dp-memo`: a chain whose window bases and seeds equal an earlier chain's
# reuses that alignment, and a DP call whose arguments and sequences equal an earlier call's reuses that result
# (overlapping loci present the same windows many times). On CPUs with AVX2 the dual-affine DP runs an AVX2 build of
# minimap2's SSE4.1 kernel with identical results (MM_KSW_NO_AVX2=1 forces the SSE kernel); the build checks the two
# kernels against each other on 20,000 random cases. With --gate, each batch's longest queries are mapped first and a
# query whose chains span >= 10 Mb has its chains aligned on all threads (results identical, written in input order).
# Without the Rustle options the binary writes the stock output byte
# for byte; its version is 2.30-r1287-gate6. mcl_families finds it as `minimap2-gate` on PATH (or RUSTLE_MINIMAP2), and
# tools/mm2_shard.sh runs it for any command line carrying --gate.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
BIN=${1:-$HOME/.local/bin}
COMMIT=79c9cc186b95f50bd899f69b48eba995ced810c6   # tag v2.30
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
git clone -q https://github.com/lh3/minimap2 "$WORK/minimap2"
git -C "$WORK/minimap2" checkout -q "$COMMIT"
git -C "$WORK/minimap2" apply "$HERE/gate.patch"
make -C "$WORK/minimap2" -j"${JOBS:-4}" minimap2 > "$WORK/build.log" 2>&1 || { tail -20 "$WORK/build.log" >&2; exit 1; }
if grep -qw avx2 /proc/cpuinfo 2>/dev/null; then   # the AVX2 DP kernel must return exactly what the SSE4.1 kernel returns
	make -C "$WORK/minimap2" ksw2_avx2_test >> "$WORK/build.log" 2>&1 || { tail -20 "$WORK/build.log" >&2; exit 1; }
	"$WORK/minimap2/ksw2_avx2_test" 20000 11 || { echo "AVX2 DP kernel differs from the SSE4.1 kernel; not installed" >&2; exit 1; }
fi
mkdir -p "$BIN"
install -m 755 "$WORK/minimap2/minimap2" "$BIN/minimap2-gate"
echo "installed $BIN/minimap2-gate ($("$BIN/minimap2-gate" --gate /dev/null --version))"
