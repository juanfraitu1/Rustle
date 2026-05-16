# Option C — faithful port of StringTie read→transfrag (one transfrag per long-read pattern)

The deferred maximal-blast-radius rewrite. This doc grounds it in a precise
cross-tool side-by-side (rlink.cpp vs map_reads.rs, both fully read) and
gives a staged, env-gated, A/B-able implementation plan.

## The reframed root cause (corrects READ_TO_TRANSFRAG_DIVERGENCE.md)

Earlier we said the divergence was "upstream in read→node mapping
(`collect_read_nodes_exact`)". The full side-by-side **refutes that**:

- **rustle `collect_read_nodes_exact` (map_reads.rs:250)** marks a node a
  member on *any* `> 0` half-open overlap — **no bp slop, no
  junction-correction in the membership test** (map_reads.rs:320,328).
- **StringTie `get_read_pattern` (rlink.cpp:4456)** marks a node a member
  on *any* `overlapLen > 0` — **no bp slop** (rlink.cpp:4486,4522/4526).

→ The two tools build the **same `unique_nodes` chain** for a read. Node
membership is *already faithful*. The divergence is **entirely
downstream of membership**:

| | StringTie | rustle |
|---|---|---|
| read → #transfrags | **exactly 1 per valid strand** (rlink.cpp:5170, single `update_abundance`) | **1..N** via 5 fragmentation sites |
| junction with no graph edge | node bits set, **edge bit unset, NOT split, NOT dropped** (rlink.cpp:4507-4512) | **SPLIT_BADJUNC** splits the read (map_reads.rs:1152-1262) |
| coverage valley | (no such concept) | **split_chimeric_transfrags** splits (map_reads.rs:672, called :947) |
| single-node read | dropped → NULL (rlink.cpp:4812) | single-node-fragment skip (map_reads.rs:920-922) |
| leading/trailing node into intron/src/sink | **trim only** (shrink, never split), guarded by hardstart/hardend (rlink.cpp:4693-4809) | `trim_longread_path_for_update_abundance` (map_reads.rs:1647) — analogous |

**Precision in StringTie is NOT read fragmentation.** A spurious/killed
junction simply is not a graph *edge*; the transfrag still carries both
node bits, but downstream path/flow cannot traverse a non-existent edge.
rustle currently enforces the same precision by *fragmenting the read* at
that junction — a different mechanism with a different (worse) side
effect: it also destroys the *valid* long-range spanning chain that
happens to cross one killed junction, which is exactly the STRG.15.1-class
ST_ONLY / NEVER_CONSTRUCTED transfrag.

Why Option A failed (9→12, F1-neg): it suppressed **only**
SPLIT_BADJUNC, and its `whole_chain_graph_valid` guard required
`a.children.contains(b)` — *false at precisely the killed junctions* — so
it rarely fired; chimeric split + single-node skip + trim still mangled
the chain. Option C removes *all* fragmentation for long reads and moves
precision to edge-absence, exactly as ST does.

## Target semantics (what Stage 0 must reproduce)

For each long read, on its valid strand:

1. `unique_nodes` = `collect_read_nodes_exact` chain (UNCHANGED — already
   faithful; pure overlap).
2. Emit **one** transfrag whose node set = `unique_nodes`, pattern =
   node bits for all of them + edge bits **only where a real graph edge
   exists** (`ensure_edges_for_read_path` already refuses killed-junction
   edges, map_reads.rs:222-228 — verify). NO split at killed junctions,
   NO chimeric split, NO single-node-fragment skip.
3. Apply the ST `update_abundance` trim (leading/trailing nodes that
   over-run into intron/source/sink), guarded so a `hardstart`/`hardend`
   terminal node is never trimmed (rustle
   `trim_longread_path_for_update_abundance` already does this — keep on).
4. A genuine single-node read → drop (matches ST NULL,
   rlink.cpp:4812).
5. Dedup by exact pattern, accumulate abundance (rustle
   `add_or_update_transfrag` already does this; ST `findtrf_in_treepat`).

## Staged plan (each stage A/B-able, env-gated, default byte-identical)

**Stage 0 — scaffold `RUSTLE_ST_READ_PATTERN=1` (default off).** In
`map_reads_to_graph_bundlenodes` (map_reads.rs:774): when the env is set
and the read is long, **skip `split_read_segments`** and instead form a
single segment = full `unique_nodes` (carry `orphan=false`); **skip the
single-node-fragment skip** (:920-922); after the read loop, **skip
`split_chimeric_transfrags`** (:947). Keep `ensure_edges_for_read_path`,
keep `add_or_update_transfrag` + its trim. ~30-60 lines, localized,
default path untouched → byte-identical 1746/1948.

**Stage 1 — measure (task #103).** `transfrag_construction_diff.py --pre`
on STRG.15.1 (expect the 9 ST_ONLY / NEVER_CONSTRUCTED → ~0); full chr19
GGO matching/query vs **1746/1948 F1=92.21**. Two outcomes:
- F1 ≥ baseline → precision held by edge-absence; promote toward default.
- F1 < baseline → a precision leak exists: path/flow is traversing a
  node-adjacency that has **no graph edge** (the thing read-fragmentation
  was papering over). Go to Stage 2.

**Stage 2 — close the precision leak at its real site.** Confirm
rustle's path-extraction / `long_max_flow` only ever steps along real
graph edges (not "both nodes are in this transfrag"). If it does step on
non-edges, that — not read fragmentation — is the precision bug; fix it
there so killed junctions can't form spurious chains while valid
long-range chains survive. (This is the ST invariant: transfrag node
membership ⊋ traversable edges.)

**Stage 3 — reconcile trim constants to ST exactly.** Map
`trim_longread_path_for_update_abundance` thresholds to ST
`longintronanchor=25` (rlink.h:33), `CHI_WIN=100`/`CHI_THR=50`
(rlink.h:17-18), `DROP=0.5` (rlink.h:13), hardstart/hardend skip
(rlink.cpp:4698-4699). Only after Stage 1/2 are green.

## Risk / gate

Highest blast radius of the arc: changes the transfrag set feeding
seeding → `long_max_flow` → every predcluster filter → all of chr19.
Mitigation: strict env gate (default byte-identical, verified each
stage); STRG.15.1 `--pre` as the direct construction readout; 1746/1948
F1 as the regression gate; staged so a red Stage 1 localizes the leak
(Stage 2) rather than forcing a full revert. StringTie changes: none
(reference only, local build untouched).

## Stage 0 IMPLEMENTED + measured (2026-05-16) — decisive redirect

`RUSTLE_ST_READ_PATTERN=1` landed (default off, both
`map_reads_to_graph_bundlenodes` and `map_reads_to_graph`,
map_reads.rs): long reads bypass `split_read_segments` +
single-node-fragment skip + `split_chimeric_transfrags`; emit one
transfrag = full `unique_nodes` chain; single-node read dropped (ST
NULL).

- **Default (env unset): byte-identical 1746/1948 F1=92.21%** (verified
  via gffcompare vs stringtie_GGO_19.gtf: matching 1746, query 1948).
- **Gate ON, full chr19: F1-NEGATIVE.** matching **1746 → 1663**, query
  1890 → F1 **89.19%**.
- **Gate ON, STRG.15.1 `--pre`: ST_ONLY 9 → 15 (WORSE), all still
  NEVER_CONSTRUCTED.** Removing *all* transfrag-level fragmentation did
  **not** construct the missing chains.

**→ Fragmentation (`split_read_segments` / `split_chimeric_transfrags`)
is NOT the read→transfrag divergence.** Decisive evidence from the
`[BUNDLEMAP]` exon→node dump (gate on): rustle *does* assemble the long
in-graph chain for these reads (e.g. a read with far-upstream exons
16479596… whose in-graph chain is `16601687-16601787, 16601926-16602024,
16602259-16602330, …, 16611955-16612061`). Its **interior introns
exactly match** ST_ONLY `16601927-16612061 nI18`
(16602025-16602259, 16602331-16602408, …). The only difference: rustle
keeps an **extra 5′ leading node `16601687-16601787`** (extra intron
16601788-16601925) that StringTie's `update_abundance` **trims off**
(rlink.cpp:4693-4809, leading node abutting intron / coverage-drop).
Different first node → different canonical intron-chain key →
`transfrag_construction_diff.py` buckets it ST_ONLY / NEVER_CONSTRUCTED
even though it is essentially the same read's chain.

**Real Option C locus = leading/trailing-node TRIM reconciliation**
(rustle `trim_longread_path_for_update_abundance` map_reads.rs:1647 vs ST
`update_abundance` rlink.cpp:4693-4809), and secondarily the candidate
node-set scoping (`ordered_nodes_from_readgroup`). NOT Stage 2
(edge-absence precision). Stage 3 is promoted to **Stage 1 of the real
work**: make rustle's boundary trim byte-faithful to ST
(`longintronanchor=25`, `CHI_WIN=100`, `DROP=0.5`,
hardstart/hardend-skip), then re-key STRG.15.1 `--pre`.

Stage-0 fragmentation-bypass is kept env-gated **default-off** (it is
F1-negative and now known not to be the lever) — useful only as the
isolation harness that produced this redirect.

## Status

Stage 0 FALSIFIED as the lever; redirected to boundary-trim
reconciliation (above). Scaffold `RUSTLE_ST_READ_PATTERN` retained
env-gated default-off, byte-identical 1746/1948. Baseline to beat:
**1746/1948 F1=92.21%**.
