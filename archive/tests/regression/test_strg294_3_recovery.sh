#!/bin/bash
# Regression test: STRG.294.3 recovery (16-exon alt-TTS isoform)
#
# StringTie emits STRG.294.3 at NC_073243.2:41949957-42027248 (16 exons,
# cov=24.79) — an alt-TTS sibling of STRG.294.1. Rustle currently
# produces a 17-exon variant (extra exon via splice to 42032760) with
# cov=2.43 that gets killed by the isofrac filter.
#
# This test asserts that Rustle emits a tx whose intron chain matches
# STRG.294.3 AND whose last exon ends at 42027248 (alt-TTS) with 16 exons.
#
# Usage: ./test_strg294_3_recovery.sh [extra-rustle-flags...]
# Exit 0 on pass, 1 on fail. Prints diagnostic summary.

set -u
# Find the repository root relative to the script location
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
BAM=$REPO/GGO_19.bam
REF=$REPO/GGO_19_stringtie.gtf
RUSTLE=$REPO/target/release/rustle
OUT=/tmp/strg294_test.gtf

BUNDLE="NC_073243.2:41947910-42033006"

"$RUSTLE" -L --only-debug-bundle --debug-bundle "$BUNDLE" "$@" "$BAM" -o "$OUT" 2>/dev/null

if [ ! -s "$OUT" ]; then
  echo "FAIL: Rustle produced no output"
  exit 1
fi

python3 <<PYEOF
import re, sys
path = "$OUT"
target_exons = [
    (41949957, 41951275), (41952343, 41952513), (41953196, 41953264),
    (41953360, 41953536), (41953762, 41953904), (41954606, 41954731),
    (41954957, 41955137), (41955241, 41955339), (41956085, 41956187),
    (41967491, 41967570), (41967769, 41967970), (41968767, 41968898),
    (41974487, 41974561), (41975018, 41975110), (42021694, 42022413),
    (42025055, 42027248),
]
target_introns = [(target_exons[i][1], target_exons[i+1][0])
                  for i in range(len(target_exons) - 1)]

tx_exons = {}
for line in open(path):
    if line.startswith('#'):
        continue
    p = line.rstrip().split('\t')
    if len(p) < 9 or p[2] != 'exon':
        continue
    m = re.search(r'transcript_id "([^"]+)"', p[8])
    if not m:
        continue
    tid = m.group(1)
    tx_exons.setdefault(tid, []).append((int(p[3]), int(p[4])))

def intron_chain(exons):
    e = sorted(exons)
    return [(e[i][1], e[i+1][0]) for i in range(len(e) - 1)]

TOL = 5
best = None
for tid, exs in tx_exons.items():
    ic = intron_chain(exs)
    if len(ic) != len(target_introns):
        continue
    ok = all(abs(a[0]-b[0]) <= TOL and abs(a[1]-b[1]) <= TOL
             for a, b in zip(ic, target_introns))
    if ok:
        best = (tid, len(exs), sorted(exs)[-1][1])
        break

if best:
    tid, n_exons, last_end = best
    print(f"PASS: tx {tid} matches STRG.294.3 intron chain "
          f"({n_exons} exons, last_end={last_end})")
    sys.exit(0)

closest_tid = None
closest_matches = 0
for tid, exs in tx_exons.items():
    ic = intron_chain(exs)
    matches = 0
    for a, b in zip(ic, target_introns):
        if abs(a[0]-b[0]) <= TOL and abs(a[1]-b[1]) <= TOL:
            matches += 1
    if matches > closest_matches:
        closest_matches = matches
        closest_tid = tid
emitted = len(tx_exons)
print(f"FAIL: no tx matches STRG.294.3's 16-exon intron chain")
print(f"  emitted tx in locus: {emitted}")
if closest_tid:
    ce = sorted(tx_exons[closest_tid])
    print(f"  closest: {closest_tid} ({len(ce)} exons, last_end={ce[-1][1]}, "
          f"{closest_matches}/{len(target_introns)} introns match)")
sys.exit(1)
PYEOF
rc=$?
exit $rc
