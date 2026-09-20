# Figures for advisor presentation

Three diagrams explaining the work. Each has a PNG (presentation) and PDF
(print-quality) version.

---

## `fig1_golga6l7_recovery.{png,pdf}` — The motivating result

Shows the GOLGA6L7 tandem triple on chr19:
- **Row 1**: reference annotation — 3 expressed paralogs
- **Row 2**: StringTie `-L` output — only L7_1 recovered (1/3)
- **Row 3**: Rustle `-L --vg` output — all three recovered (3/3)

Use this to open: "same BAM, same reads, Rustle recovers 2 paralogs StringTie
misses entirely."

## `fig2_stringtie_flow.{png,pdf}` — How StringTie works (for context)

Three-step walkthrough of StringTie's long-read assembly:
1. Reads align to exons
2. Build a splice graph (DAG) with read-support capacities on edges;
   black = consecutive junctions, red = alternative (skip) junctions
3. Max-flow decomposition: extract the heaviest path, subtract that flow,
   repeat — each extracted path is one predicted isoform

Use this when the advisor asks "what does StringTie actually do?"

## `fig3_rustle_vg_extension.{png,pdf}` — What Rustle adds

Four-panel diagram showing how we extend the StringTie model:
1. Three tandem paralog copies (same splice topology, shifted genomic coords)
2. StringTie's view: copies 2 and 3 starve because secondary alignments are
   dropped and 5-6 primary reads don't reach the TSS
3. Rustle VG view: cross-copy multi-mapping reads with 1/NH fractional
   weight + chain-homology detection of paralogous splice topology
4. Final result: SNP/alignment-score-based copy assignment + flow
   redistribution + family-graph projection → all 3 assembled as `=` matches

Use this when the advisor asks "what's new compared to StringTie?"

## Reproducing

Scripts that generate each figure live in `tools/`:
- `tools/diagram_golga6l7_recovery.py`
- `tools/diagram_stringtie_flow.py`
- `tools/diagram_rustle_vg.py`

Run each directly with `python3`; requires matplotlib.
