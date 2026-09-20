#!/usr/bin/env bash
# Wave-1 cleanup: the highest-confidence classes from tools/audit_cleanup_candidates.py
#   TEMP, LEGACY-ASSEMBLER, SUPERSEDED-VERSION, SUPERSEDED-PORTED, SUPERSEDED-CITED
#
# DRY RUN BY DEFAULT -- prints what it would do and touches nothing. Pass --apply to act.
#
# Three dispositions, so nothing irreplaceable is ever destroyed:
#   tracked            -> git mv into archive/   (recoverable from git either way)
#   untracked SOURCE   -> mv into archive/untracked/   (*.py/*.sh/*.md/*.tsv/... exist ONLY here)
#   caches and logs    -> rm -rf   (regenerable: __pycache__, .pytest_cache, *.log, *.err, *.out)
#
# Regenerate the lists first: python3 tools/audit_cleanup_candidates.py
# Rules: docs/CLEANUP_CANDIDATES.md. What is LIVE and must not be touched: docs/ACTIVE_WORKING_SET.md.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

APPLY=0
[ "${1:-}" = "--apply" ] && APPLY=1

TSV=docs/cleanup_candidates.tsv
[ -s "$TSV" ] || { echo "missing $TSV -- run: python3 tools/audit_cleanup_candidates.py" >&2; exit 2; }

CLASSES='TEMP|LEGACY-ASSEMBLER|SUPERSEDED-VERSION|SUPERSEDED-PORTED|SUPERSEDED-CITED'

# regenerable = caches, or log-ish files that nothing cites
is_regenerable () {
    case "$1" in
        *__pycache__*|*.pytest_cache*|*.pyc) return 0 ;;
        *.log|*.err|*.out|*.stderr.log|*.stdout) return 0 ;;
        *) return 1 ;;
    esac
}

mapfile -t ROWS < <(awk -F'\t' -v c="$CLASSES" 'NR>2 && $5 ~ "^("c")$" {print $2"\t"$1}' "$TSV" | sort -k2)

n_git=0; n_arch=0; n_rm=0
for row in "${ROWS[@]}"; do
    st=${row%%$'\t'*}; p=${row#*$'\t'}
    [ -e "$p" ] || continue
    if [ "$st" = tracked ]; then
        dest="archive/$p"
        if [ $APPLY -eq 1 ]; then mkdir -p "$(dirname "$dest")"; git mv -k -- "$p" "$dest" || { echo "  SKIP git mv: $p" >&2; continue; }
        else echo "  [git mv]  $p -> $dest"; fi
        n_git=$((n_git+1))
    elif is_regenerable "$p"; then
        if [ $APPLY -eq 1 ]; then rm -rf -- "$p"; else echo "  [rm -rf]  $p"; fi
        n_rm=$((n_rm+1))
    else
        dest="archive/untracked/$p"
        if [ $APPLY -eq 1 ]; then mkdir -p "$(dirname "$dest")"; mv -n -- "$p" "$dest"
        else echo "  [mv]      $p -> $dest   (untracked: this is the only copy)"; fi
        n_arch=$((n_arch+1))
    fi
done

echo
echo "wave 1: $n_git tracked -> archive/ , $n_arch untracked source -> archive/untracked/ , $n_rm caches/logs -> deleted"
[ $APPLY -eq 1 ] || { echo "(dry run -- pass --apply to act)"; exit 0; }
echo "next: cargo test --release --lib --bins ; python3 tools/audit_cleanup_candidates.py ; git add -A && git commit"
