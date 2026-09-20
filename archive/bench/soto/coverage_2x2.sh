#!/bin/bash
# The 2x2: coverage DENOMINATOR x SUBSTRATE. Four runs, strictly SERIAL (WSL2 crash rule).
#
# The question: should the denominator depend on whether a length difference between two sequences is
# INFORMATIVE or NOISE?
#
#   EXON-SUM   length is whichever isoform assembled, so two genuine copies can differ purely by isoform
#              choice (measured: both-spliced coverage 0.259 -> 0.914 when the substrate switched to
#              genomic span). max() charges that mismatch against the pair as if it were divergence.
#              PREDICTION: max() is WRONG here. This cell is the expected-failure control that makes the
#              comparison interpretable -- without it, a bad max() result is unattributable.
#
#   GEN-SPAN   both sequences assert "this locus runs from here to here". A 3x disagreement means either
#              they are not copies or one boundary is wrong, and either way the edge rule should notice.
#              max() acts as a stub filter: a rep at 0.03x its true extent can no longer buy an edge by
#              aligning fully into a complete sibling.
#              PREDICTION: this is the interesting cell.
#
# ⚠ COUNTER-EVIDENCE, stated up front so the result is not read selectively: the duplicated unit is
#   size-invariant while annotated spans are not (NPIP spans 10.6-49.4 kb around a ~16 kb cassette, block
#   length vs max(span) r=+0.196). Only 134/171 NPIP pairs can reach max-coverage 0.50 AT ALL. So a recall
#   cost is expected in BOTH gen-span cells, and the question is whether precision pays for it.
#
# ⚠ A is also the re-established BASELINE: the committed 76.2% was measured on --cross-chrom, which is now
#   deprecated, so it is not the reference for any of these.
set -uo pipefail
CACHE=${SOTO_CACHE:-/home/juanfra/winloci_scratch/soto_cache}
BIN=/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog
FA=${SOTO_FASTA:-/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa}
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=/home/juanfra/winloci_scratch/cov2x2
mkdir -p "$OUT"
[ -f "$CACHE/soto_regions.bam" ] || { echo "cache missing"; exit 1; }

run () {   # $1=name  $2=flags  $3=env
  local name=$1 flags=$2 extra=$3
  echo "=================================================================="
  echo "ARM $name | flags: $flags | env: ${extra:-<none>}"
  date '+  start %H:%M:%S'
  local t0=$SECONDS
  env $extra "$BIN" --bam "$CACHE/soto_regions.bam" --fasta "$FA" $flags \
      --threads 8 --out "$OUT/$name" > "$OUT/$name.log" 2>&1
  local rc=$?
  echo "  exit=$rc wall=$((SECONDS-t0))s"
  if [ $rc -ne 0 ]; then echo "  FAILED:"; tail -20 "$OUT/$name.log"; return 1; fi
  echo "  copies=$(($(wc -l < "$OUT/$name.copies.tsv")-1)) families=$(($(wc -l < "$OUT/$name.families.tsv")-1))"
  python3 "$HERE/soto_cache_score.py" "$OUT/$name.copies.tsv" 2>&1 | tail -14
  echo
}

run A_min_exonsum  "--refine"                            ""                            || exit 1
run B_max_exonsum  "--refine"                            "RUSTLE_ER_COVERAGE_LONGER=1" || exit 1
run C_min_genspan  "--refine --homology-genomic-span"    ""                            || exit 1
run D_max_genspan  "--refine --homology-genomic-span"    "RUSTLE_ER_COVERAGE_LONGER=1" || exit 1

echo "=================================================================="
echo "2x2 SUMMARY"
echo "=================================================================="
python3 - "$OUT" <<'PY'
import sys, os
from collections import defaultdict
O=sys.argv[1]
print(f"{'arm':<16}{'copies':>8}{'families':>10}{'med_fam':>9}{'largest':>9}{'med_span':>12}")
for a in ("A_min_exonsum","B_max_exonsum","C_min_genspan","D_max_genspan"):
    p=f"{O}/{a}.copies.tsv"
    if not os.path.exists(p): print(f"{a:<16}  (missing)"); continue
    fam=defaultdict(list); n=0
    for i,ln in enumerate(open(p)):
        if i==0: continue
        f=ln.rstrip("\n").split("\t")
        if len(f)<6: continue
        fam[f[0]].append(int(f[5])-int(f[4])); n+=1
    sz=sorted(len(v) for v in fam.values()); sp=sorted(x for v in fam.values() for x in v)
    print(f"{a:<16}{n:>8}{len(fam):>10}{sz[len(sz)//2]:>9}{max(sz):>9}{sp[len(sp)//2]:>12,}")
print("\nPREDICTION on record: B (max x exon-sum) is the expected failure; D (max x gen-span) is the test.")
PY
echo "DONE"
