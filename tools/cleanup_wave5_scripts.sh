#!/bin/bash
# Wave 5 (2026-09-23): reduce bench/ to the scripts that produce the CURRENT O1 / O2 / O3 results, the Soto
# replication, and the shared infrastructure; every other analysis script (refuted arms, one-off probes, ported
# parity oracles, fixture generators whose imports were archived in wave 2) moves to
# ~/Desktop/Rustle_attic/2026-09-23b/<same path>, tracked ones after tagging HEAD as notebook-2026-09-23b.
#   git checkout notebook-2026-09-23b -- <path>     recovers any of them.
# The KEEP list below is the reviewable inventory; its import closure is computed and kept automatically.
# usage: tools/cleanup_wave5_scripts.sh [--apply]    (dry run by default)
set -euo pipefail
REPO=/mnt/c/Users/jfris/Desktop/Rustle
ATTIC=/mnt/c/Users/jfris/Desktop/Rustle_attic/2026-09-23b
APPLY=${1:-}
cd "$REPO"

KEEP=$(cat <<'EOF'
# --- shared infrastructure
bench/soto/rustlib.py
bench/sim_reads.py
# --- O1: the RNA-level definition machinery (nested edge-test lattice, NPIP/TBC1D3) and the guided mode + truths
bench/layer_order/lattice_common.py
bench/layer_order/lattice_edges.py
bench/layer_order/lattice_expr.py
bench/layer_order/lattice_levels.py
bench/layer_order/lattice_truth.py
bench/layer_order/lattice_filtration.py
bench/layer_order/lattice_report_tables.py
bench/layer_order/lattice_check_c2.py
bench/layer_order/lo_analysis.py
bench/layer_order/lo_corrected_tables.py
bench/layer_order/lo_expr_recount.py
bench/layer_order/soto_map.py
bench/guided_pipeline.py
bench/adjudicated_truth.py
bench/annotation_nodes.py
bench/protein_families.py
bench/locus_reads.py
bench/heldout_family_score.py
bench/soto_vs_us_referee.py
bench/rna_truth_from_protein.py
bench/mcl_port.py
bench/ideal_chromosome_sim.py
# --- O2: the read-level truth, the excision robustness test, the hard-locus tool bakeoff, the Eichler comparator
bench/o2_read_truth_sim.py
bench/o2_read_truth_score.py
bench/copy_assign_excision.py
bench/hard_locus_bakeoff.py
bench/tool_bakeoff.py
bench/eichler_compare.py
# --- O3: the simulations behind the RNA-only chain (one script, three modes)
bench/missing_copy_sim.py
# --- Soto 2025 replication chain (ARI 0.6959) and its scorer
bench/soto/soto_replicate_from_sedef.py
bench/soto/soto_replicate_clustering.py
bench/soto/soto_cluster_from_shared.py
bench/soto/soto_cluster_dennislab_algorithm.py
bench/soto/soto_attach_noncoding_members.py
bench/soto/famcn_from_wssd.py
bench/soto/soto_score_against_truth.py
bench/soto/soto_bipartite_match_score.py
EOF
)
KEEPLIST=$(echo "$KEEP" | grep -v '^#' | grep -v '^$')

# import closure over bench/ python: a kept script's `import x` / `from x import` of a sibling keeps x
python3 - "$KEEPLIST" > /tmp/wave5_closure.txt <<'PY'
import sys, re, os, subprocess
keep=set(sys.argv[1].split())
py=[p for p in subprocess.run(['git','ls-files','bench'],capture_output=True,text=True).stdout.split() if p.endswith('.py')]
stems={}
for p in py: stems.setdefault(os.path.basename(p)[:-3],[]).append(p)
def imports(p):
    s=open(p,errors='replace').read(); out=set()
    for m in re.finditer(r'^\s*(?:from\s+([A-Za-z_]\w*)\s+import|import\s+([A-Za-z_]\w*))',s,re.M):
        for q in stems.get(m.group(1) or m.group(2),[]):
            if q!=p: out.add(q)
    return out
added=True
while added:
    added=False
    for p in list(keep):
        if not os.path.exists(p): continue
        for q in imports(p):
            if q not in keep: keep.add(q); added=True; print(f'closure: {p} imports {q}')
print('KEEPSET='+' '.join(sorted(keep)))
PY
grep '^closure:' /tmp/wave5_closure.txt || true
KEEPSET=$(grep '^KEEPSET=' /tmp/wave5_closure.txt | cut -d= -f2-)

# data files referenced by a kept script stay; everything else under bench/ that is not .md moves
MOVE=()
for f in $(git ls-files bench | grep -v '\.md$'); do
  case " $KEEPSET " in *" $f "*) continue;; esac
  b=$(basename "$f")
  if [[ "$f" != *.py && "$f" != *.sh ]]; then
    if grep -qF -- "$b" $KEEPSET 2>/dev/null; then echo "keep (data used by a kept script): $f"; continue; fi
  fi
  MOVE+=("$f")
done
# untracked leftovers under bench (caches) go too
for f in $(git status --short --ignored bench | grep '^!!' | cut -c4-); do MOVE+=("$f"); done

echo "kept scripts: $(echo $KEEPSET | wc -w)   moving: ${#MOVE[@]}"
# references from the load-bearing places to anything that moves
for f in "${MOVE[@]}"; do
  b=$(basename "$f")
  hits=$(grep -l -F -- "$b" REPRODUCE.md tests/*.rs tools/*.sh tools/*.py 2>/dev/null | tr '\n' ' ' || true)
  [ -n "$hits" ] && echo "⚠ referenced by live places: $f <- $hits"
done

if [ "$APPLY" = "--apply" ]; then
  git tag -f notebook-2026-09-23b HEAD >/dev/null
  mkdir -p "$ATTIC"; MAN="$ATTIC/MANIFEST.tsv"
  [ -f "$MAN" ] || printf "path\ttracked\trecover_with\n" > "$MAN"
  for f in "${MOVE[@]}"; do
    [ -e "$f" ] || continue
    if [ -n "$(git ls-files "$f")" ]; then tracked=yes; rec="git checkout notebook-2026-09-23b -- $f"; else tracked=no; rec="attic only"; fi
    mkdir -p "$ATTIC/$(dirname "$f")"; mv "$f" "$ATTIC/$f"
    [ $tracked = yes ] && git rm -r -q --cached "$f"
    printf "%s\t%s\t%s\n" "$f" "$tracked" "$rec" >> "$MAN"
  done
  find bench -type d -empty -delete
  echo "APPLIED: $(wc -l < "$MAN") manifest rows"
else
  printf '%s\n' "${MOVE[@]}" | sed 's|^|  move: |' | head -200
  echo "[DRY RUN]"
fi
