# Network Flow in Rustle — A Walkthrough

This document complements `docs/ALGORITHMS.md` §3–4 with a hands-on,
code-level walkthrough of how Rustle turns aligned long reads into
transcripts via network flow. Read this when you want to:

- Trace a specific path decision in the source
- Extend the flow logic safely
- Diagnose why Rustle emits / misses a particular transcript

Goal: after this page you should be able to open `src/rustle/max_flow.rs`
and `src/rustle/path_extract.rs` and know what each function is doing
and where each decision actually lives.

---

## 1. The object at the center: the splice graph

A **splice graph** is the per-bundle data structure built after read
alignment and junction filtering. Definitions:

- **Node** — a contiguous exon segment between two adjacent "events"
  (a bundlenode boundary, a junction donor/acceptor, or a
  coverage-cliff split point). See `GraphNode` in
  `src/rustle/graph.rs`.
- **Edge** — either
  - a **splice edge** (node `u` ends before node `v` starts with a
    gap `v.start > u.end + 1` → intron between), or
  - a **contiguous edge** (`v.start == u.end + 1`, no intron).
- **Source** (`source_id`, always node 0) — a virtual node with an
  outgoing edge to every transcript start (hardstart or coverage
  cliff entry). Every assembled transcript begins with source.
- **Sink** (`sink_id`, always the last node) — virtual, with
  incoming edges from every transcript end. Every transcript ends
  at sink.

So each real transcript is a path `source → n₁ → n₂ → … → sink`.

The splice graph is built in `create_graph_inner` (graph_build.rs:468).
Longtrim splits (coverage-cliff-based node splits with source/sink
edges) happen DURING iteration via `longtrim_inline`
(graph_build.rs:1128), so transfrags mapped later reference the final
node layout.

---

## 2. Transfrags — the "capacity" carriers

A **transfrag** is a subpath of the splice graph carrying some
abundance. There are three origins:

1. **Read-derived**: each long read is a multi-exon chain → becomes
   one transfrag with `abundance = read_weight`. See
   `map_reads_to_graph_bundlenodes` in `map_reads.rs:645`.
2. **Synthetic source/sink**: `add_coverage_source_sink_edges`
   (graph_build.rs:68) emits `source → n` and `n → sink` transfrags
   at coverage-cliff boundaries with abundance proportional to the
   cliff size. These are the "entry/exit tickets" for the flow.
3. **Longtrim-synthetic**: when `longtrim_inline` splits a node at a
   read-boundary, it emits a `source → new_start_node` and a
   `cur → new_start_node` continuity transfrag so flow still threads
   through the split.

The `GraphTransfrag` struct (graph.rs:100+) stores:
- `node_ids: Vec<usize>` — the path through the graph
- `pattern: GBitVec` — node-IDs + edge-IDs as a bitvec, used for
  fast "is this path compatible with that path" checks
- `abundance: f64` — the carried capacity
- `longstart/longend: u64` — furthest read coverage span, used to
  restrict path extension

Edge capacity is NOT a simple scalar — it's derived by summing
transfrag abundances that traverse each `(u, v)` edge. See
`build_capacity_network` in `max_flow.rs`.

---

## 3. Max flow — what and why

### The formal object

Given the splice graph `G = (V, E)` with capacities `c : E → ℝ≥0`,
max-flow finds `f : E → ℝ≥0` such that:

1. **Capacity**: `0 ≤ f(u, v) ≤ c(u, v)` for every edge.
2. **Conservation**: for every non-source, non-sink node `v`,
   `Σ_{(u,v)} f(u, v) = Σ_{(v,w)} f(v, w)`.
3. **Maximum**: the total flow out of source (equivalently, into
   sink) is as large as possible.

In Rustle (and StringTie) we use **Edmonds-Karp**, the BFS variant of
Ford-Fulkerson: repeatedly find a shortest augmenting path in the
residual graph and push its bottleneck.

### Code pointer

`bfs_augmenting_path` (max_flow.rs:88) is the BFS. Given current
`capacity`, `flow`, adjacency `link`, it writes a `pred` array
pointing each reached node at its BFS parent, and returns true if
sink is reachable.

The driver loop (max_flow.rs:1187-1220) repeatedly:

```
while bfs_augmenting_path(...) {
    reconstruct path from pred
    bottleneck = min residual along path
    for each edge on path: flow += bottleneck; reverse_flow -= bottleneck
    total_flow += bottleneck
}
```

### Why this is the right object

The edge capacity on splice edge `(u, v)` equals the number of reads
that crossed that junction (times per-base normalization). A
transcript that uses junction `(u, v)` consumes part of that
capacity. Max flow asks "how much transcript mass can I push through
this graph?" — and the answer is always less than or equal to "how
much read mass I observed". The gap (observed − assembled) is the
unexplained/noise coverage.

### What Rustle does differently from textbook Ford-Fulkerson

1. **Node capacity** — not just edges. A node's capacity = average
   per-base coverage × length. When `weighted_node_cap = true`
   (max_flow.rs:108), flow through a node is also bounded by its own
   capacity, not just edges. Enforced via the self-loop
   `capacity[u][u]`.

2. **Not one max-flow, but one max-flow per seed.** Classical
   max-flow returns *a* flow decomposition; Rustle needs the
   decomposition to MAP to real isoforms. It iterates: pick highest-
   abundance seed, run max-flow through that seed, emit as a
   transcript, subtract, repeat. See `long_max_flow_seeded`
   (max_flow.rs:2269).

---

## 4. Seed-driven decomposition

This is where Rustle's flow usage diverges most sharply from
textbook flow problems.

### The top-level loop (path_extract.rs, `parse_trflong` family)

```
sort transfrags by abundance desc
for each seed s in order:
    if s's abundance < threshold: skip
    path = [s.node_ids]                # start from seed's chain
    extend_left(path, seed)            # back_to_source_fast_long
    extend_right(path, seed)           # fwd_to_sink_fast_long
    if path is complete (source → sink):
        emit as transcript
        subtract seed's abundance along path from transfrag pool
    else:
        stash as checktrf candidate    # may still be rescued later
```

### Extending left: `back_to_source_fast_long`

Located in `path_extract.rs:4145`. Given the current path's
leftmost node `i` and the seed's pathpat (bitvec of allowed nodes/
edges), choose a parent `p` of `i` that maximizes transfrag
compatibility. Repeat until `source` is reached or no parent has
support.

Key heuristics:
- **Fast path**: if `i-1` is already on pathpat, pick it (maintains
  path continuity without a full scan).
- **Coverage drop exclusion**: if `i-1` is coordinate-adjacent
  (touching) to `i` with a sharp coverage drop, reject — likely an
  intergenic bridging artifact.
- **Witness check**: `has_lr_witness_two_splices` requires some
  transfrag to span two consecutive splice edges on the path — a
  proxy for "at least one real read supports this stretch".

### Extending right: `fwd_to_sink_fast_long`

Symmetric (path_extract.rs:3518). Given current rightmost node `i`,
pick a child `c` by scoring each child's transfrag support
(`childcov`) and picking the highest. Near-ties use a `nodecov`
tiebreak.

### Why "seed-driven"?

A long read that spans 15 exons *directly witnesses* a 15-exon
isoform. Using it as a seed sidesteps the "which flow decomposition"
ambiguity — we anchor on the chain that actual reads exhibited, then
only the source-side and sink-side tails need path decisions.

### The BFS shortest-augmenting-path inside the seed loop

`long_max_flow_seeded_with_used` (max_flow.rs:2282) runs Edmonds-Karp
*within the seed's pathpat*. Because the seed restricts which nodes
can carry flow, the max-flow returned is the maximum mass the seed's
specific chain can carry — that's the transcript's abundance.

---

## 5. A concrete walkthrough: a 3-isoform locus

Simplify: suppose a locus has exons A, B, C, D, E and three real
isoforms:
- **Isoform 1** (dominant): A → B → C → D → E (cov = 100)
- **Isoform 2** (alt-splice): A → B → D → E (skips C, cov = 20)
- **Isoform 3** (alt-TSS): B → C → D → E (short, cov = 5)

### Graph after construction

Nodes: `src`, `A`, `B`, `C`, `D`, `E`, `sink`.

Edges with capacities:
```
src → A: 100       (isoform 1 start)
src → B: 5         (isoform 3 hardstart TSS)
A → B: 100+20=120  (iso 1 + iso 2)
B → C: 100+5=105   (iso 1 + iso 3)
B → D: 20          (iso 2 skip-splice)
C → D: 100+5=105   (iso 1 + iso 3)
D → E: 125         (all three)
E → sink: 125
```

Read-derived transfrags (simplified — imagine one transfrag per
observed isoform, abundance proportional to read count):

- t1: A–B–C–D–E, abund 100 (long reads covering full iso 1)
- t2: A–B–D–E, abund 20 (long reads covering iso 2)
- t3: B–C–D–E, abund 5 (long reads covering iso 3)

### Seed selection

Sort by abundance desc → t1, t2, t3.

### Iteration 1: seed = t1

`pathpat = {src, A, B, C, D, E, sink}` (all nodes t1 touches plus
source/sink).

- `back_to_source_fast_long` from A: source is a parent → path
  extends `src → A → B → C → D → E`.
- `fwd_to_sink_fast_long` from E: sink is a child → path extends to
  sink.
- Emit isoform 1 with abundance = min edge capacity along path = 100.
- Subtract 100 from every edge on the path.

After subtraction:
```
src → A: 0; A → B: 20; B → C: 5; C → D: 5; D → E: 25; E → sink: 25
(src → B: 5 and B → D: 20 untouched)
```

### Iteration 2: seed = t2

`pathpat = {src, A, B, D, E, sink}`.

- Back from A: src → A capacity is now 0. Could also use reverse
  flow from the first iteration (residual) — in practice the
  decomposition is monotone and t2's own transfrag capacity on
  `A → B → D → E` was preserved (B → D is disjoint from t1's path).
- Path: `src → A → B → D → E → sink`.
- Emit isoform 2 with abundance 20.

### Iteration 3: seed = t3

- Back from B: `src → B` capacity = 5. Path ends in source.
- Forward from E: sink. Path ends.
- Emit isoform 3 with abundance 5.

All three isoforms recovered with correct abundances.

### What could go wrong

- If isoform 2's seed (t2) was short — say only 2 exons — the back
  extension would need to guess A via transfrag compatibility; if
  the dominant t1 has already depleted `A → B`'s residual, the path
  might not find a way back to source and get stashed in checktrf.
- If an extra "noise" edge `B → E` existed with capacity 3, `fwd_to_sink`
  from B might pick the direct `B → E` route over `B → C → D → E`
  for a seed that doesn't explicitly include C, producing a 2-exon
  rather than 4-exon transcript.

The latter is exactly the STRG.294.3 class — see
`memory/project_strg294_flow_path_choice.md` for a real example and
why a simple "prefer contiguous" tiebreak doesn't fix it.

---

## 6. Where to look in the source

| Concern | File / function |
|---|---|
| Graph construction | `graph_build.rs::create_graph_inner` |
| Inline longtrim splits | `graph_build.rs::longtrim_inline` |
| Read → transfrag mapping | `map_reads.rs::map_reads_to_graph_bundlenodes` |
| Transfrag patterns & edge bits | `graph.rs::GraphTransfrag` |
| BFS augmenting path | `max_flow.rs::bfs_augmenting_path` |
| Seeded max-flow | `max_flow.rs::long_max_flow_seeded_with_used` |
| Left extension | `path_extract.rs::back_to_source_fast_long` |
| Right extension | `path_extract.rs::fwd_to_sink_fast_long` |
| Witness check (two-splice) | `killed_junctions.rs::has_lr_witness_two_splices` |
| Seed ordering / main loop | `path_extract.rs::parse_trflong` |
| checktrf rescue/redistribute | `path_extract.rs` (search for `checktrf`) |

Tracing knobs (set at runtime via env vars):

- `RUSTLE_TRACE_LOCUS=start-end` — scope all traces to a locus
- `RUSTLE_TRACE_FWD`, `RUSTLE_TRACE_BACK` — per-child extension
  decisions in fwd/back
- `RUSTLE_TRACE_LONGTRIM=1` — longtrim split events
- `RUSTLE_FILTER_DEBUG=1` — dump per-tx filter decisions
- `RUSTLE_GRAPH_DOT_DIR=/tmp/dots` — write per-bundle GraphViz
  DOT files for the splice graph

---

## 7. Two properties that matter for debugging

### Flow decomposition is not unique

Two different path sets can realize the same max-flow value. Rustle's
greedy seed-ordered decomposition picks *a* decomposition, which is
usually biologically sensible because dominant isoforms have more
reads and get processed first. But "usually" is not "always" — see
STRG.294.3.

Concretely: if at some forking node both a splice edge and a
contiguous edge have sufficient residual capacity to carry the
current seed's abundance, the choice is heuristic (
highest-transfrag-cov-support in `fwd_to_sink_fast_long`). StringTie
sometimes takes the opposite choice, which is how the same graph +
same reads + same filters can produce a different final GTF.

### Pattern bits must stay consistent with node IDs

`GraphTransfrag::pattern` is a `GBitVec` indexed by node ID + edge ID
(edge IDs live in `graph.gpos`). If you change the graph (add nodes,
split nodes, prune nodes) after transfrags are built, you MUST
rebuild the patterns via `rebuild_transfrag_patterns`. Otherwise
witness checks, onpath checks, and max-flow compatibility will read
stale bits and silently mis-decide.

This is why `apply_longtrim_direct` is currently run before
`map_reads_to_graph_bundlenodes` — so transfrags are built on the
final graph. The opt-in post-hoc split path
(`RUSTLE_LONGTRIM_APPLY_ON=1`) has a documented F1 regression
precisely because of stale patterns.

---

## 8. Mental model summary

1. The splice graph is a factored representation of reads: nodes =
   exon segments, edges = observed junctions.
2. Reads become transfrags carrying abundance. Transfrag abundance
   sums give edge capacities.
3. Max-flow finds the total isoform mass the graph can support.
   Flow decomposition turns that total into individual
   transcript-path assignments.
4. Rustle's decomposition is seed-driven: pick highest-abundance
   transfrag, extend it source-ward and sink-ward, emit, subtract,
   repeat.
5. The quality of the output hinges on (a) how well the graph
   reflects real biology (→ graph construction), (b) how well
   seeds are ordered and extended (→ path_extract), and (c) how
   many spurious paths survive the final filter stack.

If an observed transcript is missing from Rustle's output, the
question to ask is: which of those three levels failed?
