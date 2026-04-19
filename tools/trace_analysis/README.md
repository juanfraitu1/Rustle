# StringTie trace analysis tools

Diff-based divergence analysis between Rustle and StringTie using
`stringtie_debug/trace_GGO_19_deep.log` (910MB full instrumented trace).

## Workflow

```bash
# 1. Rank divergent loci by |st_tx - rustle_tx|
python3 tools/trace_analysis/rank_divergent_loci.py

# 2. Slice the full trace per target bundle
python3 tools/trace_analysis/slice_trace_by_locus.py
# writes /tmp/locus_traces/*.log

# 3. Generate per-locus summary (StringTie vs Rustle output counts + event histograms)
python3 tools/trace_analysis/locus_divergence_summary.py

# 4. Event counts per bundle (COLOR_BREAK, EJUNC_DELETE, etc.)
python3 tools/trace_analysis/per_locus_event_stats.py

# 5. Extract every FILTER / PRED_FATE / FINAL_KEEP line per locus
python3 tools/trace_analysis/extract_key_decisions.py
# writes /tmp/locus_traces/KEY_DECISIONS.md (all decisions for every locus)

# 6. Compact tabular summary
python3 tools/trace_analysis/generate_compact_report.py
# writes /tmp/locus_traces/COMPACT_SUMMARY.md
```

## Output files

- `COMPACT_SUMMARY.md` — kill-reason histogram + FINAL keep/drop per locus
- `KEY_DECISIONS.md` — full FILTER lines with KILLED_BY pred IDs and reasons

## Reference

The full trace file `stringtie_debug/trace_GGO_19_deep.log` contains:
- ~992 bundle entries (PASS_ID=1..992)
- ~20 function categories (build_graphs, good_junc, create_graph, print_predcluster, etc.)
- ~91,552 unique event patterns
- Full decision lineage per bundle: read processing → junction decisions →
  graph construction → path extraction → filter fate

This is the canonical source for "what StringTie actually decided" when
investigating divergence.
