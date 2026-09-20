#!/bin/bash
# Extract family counts and per-drop-reason from sweep logs.
cd /scratch/jxi21/Assembler/Rustle/bench/jaccard_sweep

printf "%-8s %-8s %-8s %-8s %-8s %-8s %-15s\n" "thr" "raw" "kept" "low_pj" "low_sh" "megafam" "high_cv"
printf "%-8s %-8s %-8s %-8s %-8s %-8s %-15s\n" "----" "----" "----" "------" "------" "------" "-------"
for thr in 0.05 0.10 0.15 0.20 0.25 0.30 0.40 0.50; do
  log="sweep_${thr}.log"
  [ -f "$log" ] || continue
  raw=$(grep -oP "Discovered \K\d+" "$log" | head -1)
  kept=$(grep "family-quality filter" "$log" | grep -oP "→ \K\d+" | head -1)
  drop_pj=$(grep "family-quality filter" "$log" | grep -oP '"low_primitive_jaccard": \K\d+' || echo 0)
  drop_sh=$(grep "family-quality filter" "$log" | grep -oP '"low_shared": \K\d+' || echo 0)
  drop_mf=$(grep "family-quality filter" "$log" | grep -oP '"megafamily": \K\d+' || echo 0)
  drop_cv=$(grep "family-quality filter" "$log" | grep -oP '"high_exon_cv": \K\d+' || echo 0)
  printf "%-8s %-8s %-8s %-8s %-8s %-8s %-15s\n" "$thr" "$raw" "$kept" "$drop_pj" "$drop_sh" "$drop_mf" "$drop_cv"
done
