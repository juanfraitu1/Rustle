#!/bin/bash
# Wave 4 (2026-09-23): move superseded documents, earlier logs/data tables and caches out of the repo into
# ~/Desktop/Rustle_attic/2026-09-23/<same relative path>. Tracked files are `git rm --cached` (deletion STAGED,
# not committed) after tagging HEAD as notebook-2026-09-23, so every tracked file stays recoverable with
#   git checkout notebook-2026-09-23 -- <path>
# Untracked files (caches, run outputs) exist only in the attic afterwards.
# usage: tools/cleanup_wave4_attic.sh [--apply]   (dry run by default)
set -euo pipefail
REPO=/mnt/c/Users/jfris/Desktop/Rustle
ATTIC=/mnt/c/Users/jfris/Desktop/Rustle_attic/2026-09-23
APPLY=${1:-}
cd "$REPO"

# class<TAB>path
LIST=$(cat <<'EOF'
superseded-doc	README.md
superseded-doc	AGENT_HANDOFF.md
superseded-doc	docs/OBJECTIVES_AND_VERIFICATION.md
superseded-doc	docs/o1_investigations.md
superseded-doc	docs/NUMBERS.md
superseded-doc	docs/ONE_METHOD.md
superseded-doc	docs/METHOD_PSEUDOCODE.md
superseded-doc	docs/RETIREMENT_AND_MIGRATION.md
superseded-doc	docs/o1_catalog_provenance.md
superseded-doc	docs/OPEN_ITEMS_2026-09-09.md
superseded-doc	docs/o3_missing_copy_evidence.md
superseded-doc	docs/REPRODUCE.md
superseded-doc	docs/METHOD_EMAIL.txt
superseded-doc	docs/CLEANUP_CANDIDATES.md
superseded-doc	docs/cleanup_candidates.tsv
superseded-doc	docs/archive
superseded-doc	docs/artifacts
superseded-doc	docs/experiments
superseded-doc	docs/superpowers
superseded-doc	bench/ASJ.md
superseded-doc	bench/OBJECTIVES_FLOW.md
superseded-doc	bench/PANELS_AND_NOTES.md
superseded-doc	bench/REFERENCE_ABSENT_AND_UNMAPPED.md
earlier-data	docs/audit46_2026-09-04.tsv
earlier-data	docs/rna_bp1_p9_cores_2026-09-04.tsv
earlier-data	docs/gw_units_v1_2026-09-05.tsv
earlier-data	docs/gw_units_v1_params_2026-09-05.tsv
earlier-data	docs/npip_truth_audit_2026-09-05.tsv
earlier-data	docs/o2_expected_sites_2026-09-05.txt
earlier-data	docs/o2_unexplained_v8_families_2026-09-05.tsv
earlier-data	docs/rna_admit_v8_2026-09-05.tsv
earlier-data	docs/rna_units_v8_blocks_2026-09-05.tsv
earlier-data	docs/rna_units_v9_2026-09-05.tsv
earlier-data	docs/rna_units_v9_merged_2026-09-05.tsv
earlier-data	docs/roster_admission_summary_2026-09-05.tsv
earlier-data	docs/soto_fragmentation_diagnosis_2026-09-05.txt
earlier-data	docs/soto_rna_admitted_2026-09-05.tsv
earlier-data	docs/sweep_v10_families_2026-09-05.tsv
earlier-data	docs/sweep_v10_fast_families_2026-09-05.tsv
earlier-data	docs/sweep_v2_units_2026-09-05.tsv
earlier-data	docs/sweep_v3_nofilter_families_2026-09-05.tsv
earlier-data	docs/sweep_v3_unbounded_families_2026-09-05.tsv
earlier-data	docs/sweep_v3_unbounded_nofilter_families_2026-09-05.tsv
earlier-data	docs/sweep_v7_star_families_2026-09-05.tsv
earlier-data	docs/sweep_v9_final_families_2026-09-05.tsv
earlier-data	docs/sweep_v9_gpair_families_2026-09-05.tsv
earlier-data	docs/sweep_v9_pair_final_2026-09-05.tsv
earlier-data	docs/excision_sweep_v2_2026-09-09.txt
earlier-data	docs/stricter_options_catalog.tsv
earlier-data	test_data/vg_hmm
superseded-script	tools/cleanup_wave1.sh
superseded-script	tools/cleanup_wave2_bench.py
superseded-script	tools/cleanup_wave3_outputs.py
superseded-script	tools/family_vg_report.py
cache	bench/__pycache__
cache	tools/demo
cache	tests/fixtures/same_chrom_supplement/out_conflict_enum.copies.fa
cache	tests/fixtures/same_chrom_supplement/out_conflict_enum.copies.tsv
cache	tests/fixtures/same_chrom_supplement/out_conflict_enum.families.tsv
cache	tests/fixtures/same_chrom_supplement/out_conflict_enum.pairs.tsv
cache	tests/fixtures/same_chrom_supplement/out_famcn.copies.fa
cache	tests/fixtures/same_chrom_supplement/out_famcn.copies.tsv
cache	tests/fixtures/same_chrom_supplement/out_famcn.famcn.tsv
cache	tests/fixtures/same_chrom_supplement/out_famcn.families.tsv
cache	tests/fixtures/same_chrom_supplement/out_famcn.pairs.tsv
cache	tests/fixtures/same_chrom_supplement/out_hom.copies.fa
cache	tests/fixtures/same_chrom_supplement/out_hom.copies.tsv
cache	tests/fixtures/same_chrom_supplement/out_hom.families.tsv
cache	tests/fixtures/same_chrom_supplement/out_hom.pairs.tsv
cache	tests/fixtures/same_chrom_supplement/out_hom_protein_qc.copies.fa
cache	tests/fixtures/same_chrom_supplement/out_hom_protein_qc.copies.tsv
cache	tests/fixtures/same_chrom_supplement/out_hom_protein_qc.families.tsv
cache	tests/fixtures/same_chrom_supplement/out_hom_protein_qc.pairs.tsv
EOF
)

if [ "$APPLY" = "--apply" ]; then
  git tag -f notebook-2026-09-23 HEAD >/dev/null
  mkdir -p "$ATTIC"
  MAN="$ATTIC/MANIFEST.tsv"
  [ -f "$MAN" ] || printf "path\ttracked\tclass\trecover_with\n" > "$MAN"
fi
n=0; nt=0; bytes=0
while IFS=$'\t' read -r cls p; do
  [ -z "$p" ] && continue
  [ -e "$p" ] || { echo "skip (absent): $p"; continue; }
  if git ls-files --error-unmatch "$p" >/dev/null 2>&1 || [ -n "$(git ls-files "$p" 2>/dev/null)" ]; then tracked=yes; else tracked=no; fi
  sz=$(du -sb "$p" | cut -f1)
  n=$((n+1)); bytes=$((bytes+sz)); [ $tracked = yes ] && nt=$((nt+1))
  printf "%-18s tracked=%-3s %9d  %s\n" "$cls" "$tracked" "$sz" "$p"
  if [ "$APPLY" = "--apply" ]; then
    mkdir -p "$ATTIC/$(dirname "$p")"
    mv "$p" "$ATTIC/$p"
    if [ $tracked = yes ]; then git rm -r -q --cached "$p"; rec="git checkout notebook-2026-09-23 -- $p"; else rec="attic only (never tracked)"; fi
    printf "%s\t%s\t%s\t%s\n" "$p" "$tracked" "$cls" "$rec" >> "$MAN"
  fi
done <<< "$LIST"
echo "entries: $n (tracked $nt), bytes: $bytes  [$([ "$APPLY" = "--apply" ] && echo APPLIED || echo DRY RUN)]"
