# Glass-box demo for Rustle multi-copy gene-family work — advisor meeting prep

**Date:** 2026-04-26 (meeting on 2026-04-27)
**Author:** advisor-meeting prep with Claude
**Status:** draft for sign-off

## Problem

Advisor objections to address tomorrow:
1. "No working assembler"
2. "No understanding of complexity"
3. "Results are luck or hardcoded"

Advisor refuses to read source code on principle. Pre-rendered narrative figures (`figures/fig1-3`) have already been seen and discounted. We need a demo that earns credibility from inputs and outputs alone.

## Goal

Produce a defensible **glass-box live demo** the advisor can watch end-to-end on data he can verify is unmodified, with every algorithmic stage emitting a TSV/PNG he can inspect, culminating in the multi-copy gene family result.

Out of scope: implementing Phase 3 variation-graph alignment in Rust. We sketch it conceptually only.

## Tier structure (Combo D, recommended)

Three escalating tiers, each defeats one objection:

| Tier | Dataset | Purpose | Defeats |
|---|---|---|---|
| 1 | `synthetic_family/` (95 reads, deterministic seed=42, ground-truth GTF) | Hand-traceable walkthrough | "Hardcoded" — generator + truth are independent files |
| 2 | `GGO_19.bam` chr19 GOLGA6L triple at NC_073243.2:104.79–104.88 Mb (LOC115931294, LOC134757625, LOC101137218 — all + strand, all 8 introns, "golgin subfamily A member 6-like protein 7") | Live run on real data with all stages dumped | "No assembler" + "No understanding" |
| 3 | Cross-chromosome GOLGA family table (precomputed from `MULTI_COPY_FAMILY_PROOF.md`) | Closing slide showing single-locus isn't a fluke | "Single locus is luck" |

Tier 2 is the centerpiece. Tier 1 builds credibility for Tier 2. Tier 3 is one slide.

## Five visuals to produce

All generated mechanically from data dumps — no hand-laid layouts.

### V1 — Splice graph from `RUSTLE_PARITY_GRAPH_TSV` + new edges dump

For one bundle (synthetic Tier 1, then chr19 Tier 2):
- Nodes labeled with exon coordinates and coverage
- Edges labeled with type (collinear / junction) and capacity
- Layout: graphviz `dot` (left-to-right, mechanical)
- Output: `graph_dump.tsv` (nodes) + `graph_edges.tsv` (edges) → PNG

### V2 — Three side-by-side isomorphic graphs of the chr19 triple

- Each of LOC115931294, LOC134757625, LOC101137218 rendered with V1 pipeline
- NetworkX VF2 isomorphism algorithm computes node bijection
- Bijection drawn as colored arcs between the three layouts
- "Similar" defined precisely as: graph isomorphism preserving (a) node count, (b) edge multiset by type (collinear/junction), (c) source/sink topology, with genomic coordinates normalized by offset within locus

### V3 — Per-iteration flow decomposition

- Requires new `RUSTLE_FLOW_ITER_TSV` dump (Edmonds-Karp augmenting path + bottleneck + residual-after each iteration)
- Render: small grid of graph snapshots, one per iteration, residual capacities color-coded
- Annotated with: extracted path = transcript

### V4 — Multi-mapper read pinning

- Pick one specific NH≥2 read crossing 2-3 paralogs in chr19 triple
- Diagram: read aligned at each copy, with weight under (a) StringTie 1/NH uniform = 0.33 each, (b) Rustle EM after convergence (one copy gets ~1.0, others ~0.0)
- Source: existing VG trace stderr (parse it; do NOT add new dump for this)

### V5 — Variation graph for Phase 3 (conceptual)

- Synthetic minimal example: 2 copies → variation graph with bubble nodes at SNP positions
- One read with copy-discriminating SNPs aligned to the graph, settling on one path
- Pure illustration, hand-built (matplotlib or graphviz), no Rust code
- Caption clearly labels: "Future direction, not yet implemented"

## Code/dump additions required

Two genuinely new Rust dumps; everything else uses existing instrumentation.

| New dump | Env var | Where | Cost |
|---|---|---|---|
| Graph edges + capacities | `RUSTLE_PARITY_GRAPH_EDGES_TSV` | path_extract.rs near existing graph nodes dump (~5435) | Low — iterate edge map alongside existing node loop |
| Flow per-iteration state | `RUSTLE_FLOW_ITER_TSV` | max_flow.rs, inside Edmonds-Karp augmenting-path loop | Medium — need to capture path + residual without breaking the algorithm |

Skip for now (time pressure):
- VG EM per-iteration dump (V4 derives from existing stderr trace)
- Final transcripts dump (the GTF is already final)

## Existing infrastructure we reuse

- `RUSTLE_PARITY_READ_DUMP_TSV` (reads in)
- `RUSTLE_PARITY_PARTITION_TSV` + `RUSTLE_PARITY_BUNDLENODE_TSV` (bundle)
- `RUSTLE_PARITY_GRAPH_TSV` (graph nodes)
- `RUSTLE_READ_NODE_TSV` (reads → graph nodes)
- `RUSTLE_PARITY_TF_TSV` (transfrags)
- `RUSTLE_DEBUG_STAGE_TSV` (per-bundle summaries with 30+ counters)
- `tools/trace_analysis/slice_trace_by_locus.py` (filter dumps to a locus)

## New tools to build (all Python, in `Rustle/tools/demo/`)

| Tool | Input | Output |
|---|---|---|
| `render_splice_graph.py` | nodes TSV + edges TSV + locus | PNG (graphviz) |
| `compare_graphs_isomorphism.py` | 3× nodes+edges TSVs | bijection JSON + 3-panel PNG with arcs |
| `render_flow_iterations.py` | flow_iter TSV + nodes/edges TSV | grid PNG |
| `render_multimapper_pinning.py` | read ID + VG trace stderr | side-by-side PNG |
| `render_phase3_vg_sketch.py` | synthetic params (hand-coded) | PNG |
| `run_demo.sh` | (none) | orchestrates all of the above end-to-end on both Tier 1 and Tier 2 |

## Demo runbook (meeting flow)

1. **Open Terminal** — empty workspace, all data files visible.
2. **Tier 1: Synthetic walkthrough** (~5 min)
   - `cat truth.gtf` (independent ground truth, advisor reads first)
   - `cat generate_isoseq.py | head -30` (deterministic generator, seed=42)
   - `samtools view reads_sorted.bam | head -5` (5 reads visible)
   - Run rustle with all dumps enabled, `cat` each TSV in order: reads → bundle → graph → flow → transcripts
   - Compare last to truth.gtf
3. **Tier 2: Real data** (~10 min)
   - Same workflow on the chr19 GOLGA6L triple region from `GGO_19.bam`
   - Show V1, V2 (the headline isomorphism), V3, V4
   - Compare StringTie output side-by-side: 1/3 vs 3/3 recovery
4. **Tier 3: Closing slide** (~3 min)
   - Cross-chromosome GOLGA family table, total recovery deltas
5. **Phase 3 sketch** (~3 min) — V5 conceptual diagram + 2-minute pitch.

## What we explicitly don't claim

- We do NOT claim a working production variation-graph aligner — that's Phase 3 future.
- We do NOT claim isomorphism for ALL paralogs in a family — the GOLGA broader family has heterogeneous intron counts, so isomorphism is sub-clade-restricted. We make this scoping explicit when discussing GOLGA.
- We do NOT claim novel-paralog discovery is finished — the scaffold exists in `vg.rs:discover_novel_copies` but is incomplete.

## Risks and mitigations

| Risk | Mitigation |
|---|---|
| Rust dump additions don't compile in time | Skip them: V1/V2/V3 still possible from existing graph node TSV (V3 uses node-coverage-residual approximation; less precise but defensible) |
| Live demo crashes mid-run | Pre-record terminal session as backup video; have all dumps pre-saved on disk |
| Advisor disputes a TSV value | Show the source file/line that emitted it AND the BAM record(s) that produced it — verifiability via `samtools view` |
| Chr19 triple isomorphism fails to hold (graphs not actually isomorphic) | This would be a real finding; we report it honestly. Expected to hold based on annotation (all 8 introns same orientation), but assembled graphs may differ slightly. We pre-run before meeting and adjust narrative if needed |
| Time runs out today | Priority order: V1 → V2 → run scripts on Tier 1 → V3 → Tier 2 dumps → V4 → V5. Drop V5 first, then V4 if needed. |

## Acceptance for this design

User sign-off needed on:
1. Tier structure (combo D)
2. Five visuals as scoped
3. Two dump additions only (graph edges, flow iters); skip VG EM dump
4. New tooling lives in `tools/demo/`
5. Spec for "similar" = graph isomorphism on (node count, edge multiset by type, source/sink topology) with normalized coords

## Next step after sign-off

Invoke `superpowers:writing-plans` to produce ordered implementation steps with verification checkpoints, time budgets, and a critical-path graph. Then `superpowers:executing-plans` (or `subagent-driven-development` for parallel pieces) to actually build.
