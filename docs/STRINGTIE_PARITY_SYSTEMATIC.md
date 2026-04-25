# Systematic StringTie / rlink parity: traces, coverage, and chain selection

This document ties together **trace outputs**, **`rlink.cpp` reasoning**, and **Rustle** so we can separate (a) missing paths StringTie would emit from (b) extra paths StringTie would not emit, with **coverage consistency** as the first-class invariant.

For operational handoff (commands, env flags, current batch status), see:
`docs/STRINGTIE_LONGREAD_PARITY_CONTINUATION_GUIDE.md`.

## 1. What the trace files are for

| Artifact | Role |
|----------|------|
| `--trace-reference REF.gtf` + `--trace-output` | Per-reference transcript: **NotExtracted** (never built) vs **filtered** (built then removed), with **blocker stage** (`src/rustle/trace_reference.rs`, `find_print_predcluster_killing_stage`). |
| `trace_blockers_*.tsv` | Aggregated blocker counts from trace/bootstrap runs (e.g. `isofrac`, `not_extracted:junctions_present`). |
| `GGO_19_*_stage_explain.lineage.*_trace.tsv` | Bundle-level lineage / stage dumps for deep parity debugging. |
| `GGO_19_trace_replay_gap_audit*.md` | Coverage of **trace event families** vs replay harness (which C++ stages are reconstructible). |

**Workflow:** Run rustle on the same BAM as StringTie reference GTF → trace → classify misses by blocker → map blocker to pipeline stage → compare that stage’s **coverage inputs** to `rlink.cpp`.

## 2. rlink.cpp: coverage vocabulary (single source of truth)

### 2.1 Bundle-level read coverage (`bpcov`)

- **`cov_edge_add` / `add_read_to_cov`**: Strand-aware cumulative coverage; neutral reads also update `bpcov[1]` (all-strand track). Pair overlap is counted once where noted in comments.
- **`get_cov(s, start, end, bpcov)`**: Integrates the cumulative array over `[start,end]` on strand track `s` (see `rlink.cpp` ~1830–1851).

This is the **ground truth for intron retention / incomplete / strand-neutral** checks inside **`print_predcluster`** (e.g. introncov vs exoncov with `intronfrac`, long-read vs short-read branches ~18445–18590).

### 2.2 Transfrag abundance vs flow-derived `pred->cov`

- **Long path extraction (`parse_trflong`)**: `longcov = transfrag[t]->abundance` before flow (~10094). Final prediction **`cov`** comes from **flux / node depletion** along the chosen path (store_transcript / exon accumulation logic; see Rust port comments in `path_extract.rs` ~6340+).
- **`print_predcluster`**:
  - Priority / max-interval structure uses **`abs(tlen) * cov`** for ordering (~18424, 18539).
  - Exon-interval bookkeeping uses **per-bp `pred->cov`** for long reads (`excov = pred->cov`, not multiplied by length ~18463–18466), vs short reads where `excov *= abs(tlen)` before adding to intervals.
  - **Longunder** (dominant-isoform fraction) uses **`pred[n]->cov`** against **`usedcov` / `multicov`** (~18975–19011), with **`isofraclong`** scaled by `CHI_WIN` / `CHI_THR` using the same **`pred[n]->cov`** (~18991–18994), not `longcov`.
  - **Short-read branch** can **recompute** coverages from `bpcov` in the `else` at ~19017+ (“recompute coverages -> mostly for short reads; long reads should be more reliable…”).

**Implication:** StringTie’s “correct” chain choice is defined by **same graph + same flux bookkeeping + same predcluster ordering and thresholds**. Any Rust drift in **nodeflux**, **noderate**, or **which number is stored as `coverage` vs `longcov`** shows up as **wrong survivors** (extra or missing isoforms), not only as gffcompare boundary noise.

### 2.3 Readthr / rawreads

- **`pred[n]->cov < readthr`** drops predictions in **rawreads** printing (~20090+). Rustle’s **`readthr_gate`** uses flow **`coverage`** (`pred->cov`) by default; set **`RUSTLE_READTHR_LONGCOV_FALLBACK`** to also retain long-read transcripts when **`longcov >= threshold`** but **`coverage`** is low (legacy heuristic).

## 3. Rustle: where coverage is defined (audit checkpoints)

| Concept | Rustle location | Should match |
|--------|-----------------|--------------|
| Bundle `bpcov` | `bam.rs` / bundle build | `add_read_to_cov` + `get_cov` semantics |
| `Transcript.coverage` | Flow / parse_trflong port (`path_extract.rs`) | C++ `pred->cov` after flow (per-bp or flux-derived as in port comments) |
| `Transcript.longcov` | Seed abundance; preserved through flow | C++ `longcov` / `transfrag->abundance` |
| `print_predcluster` | `transcript_filter.rs` | `rlink.cpp` `print_predcluster` ordering, longunder, pairwise, runoff, readthr |
| Junction-path / assist transcripts | `pipeline.rs` + filter “protected” sources | Same **predcluster** inputs; bottleneck cov must not use a different unit than flow transcripts |

**Optional Rust extensions (parity drift risk):** `RUSTLE_READTHR_LONGCOV_FALLBACK` restores longcov-based readthr survival; unanchored multi-exon threshold bump; rescue / guide protection. Prefer default **rlink-aligned** `pred->cov` / `coverage` for gates.

### 3.1 Where extra transcripts come from (`source` tags + histogram)

- **Max-flow / `parse_trflong` de novo** transcripts are tagged **`flow`** in GTF (`path_extract.rs`). **`short_flow`** is the short-read analogue.
- **`checktrf_rescue`**: checktrf rescue paths; **`junction_path`**, **`ref_chain`**, **`ref_chain_*`**: see `pipeline.rs` (`emit_junction_paths`, reference assist).
- After a run, **`RUSTLE_SOURCE_HISTOGRAM=1`** prints a stderr table: count and % per `source`.
- Offline: **`scripts/gtf_transcript_source_histogram.py out.gtf`**.

Use the histogram to see whether overproduction is dominated by **`flow`** (graph + seeds: tune caps, readthr, isofrac) vs **`junction_path`** (lower `max_junction_paths`, bottleneck cov) vs **reference assists** (`ref_chain*`, traces).

**Optional emit suppressors (parity experiments):** `RUSTLE_DISABLE_TERMINAL_ALT_ACCEPTOR`, `RUSTLE_DISABLE_MICRO_EXON_RESCUE`, `RUSTLE_CHECKTRF_ABUNDANCE_FLOOR_FRAC` (see `bench/full_bam_parity_run/SOURCE_HISTOGRAM.txt`).

### 3.2 Bundle coverage TSV (`RUSTLE_BUNDLE_COV_DUMP`)

With **`RUSTLE_TRACE_LOCUS=start-end`** and **`RUSTLE_BUNDLE_COV_DUMP=1`**, stderr emits **`RUSTLE_BUNDLE_COV_TSV`** rows at predcluster stages **`entry`**, **`after_pairwise`**, **`after_isofrac`**, **`after_readthr`** (transcripts overlapping the locus). Use for diffing **`cov` / `longcov` / `bpcov_cov` / `abs(tlen)*cov`** against a C++ trace on the same region.

## 4. Systematic error taxonomy (ties trace → cause)

### 4.1 “Miss StringTie path” (false negative)

| Trace / gffcompare signal | Likely mechanism | Check first |
|---------------------------|------------------|-------------|
| `not_extracted:junctions_present` | Graph/seed didn’t represent chain, or early filter dropped transfrag | Seeds, `process_transfrags`, junction graph vs `rlink` |
| `not_extracted` + no junctions | Bundle bounds, strand, or trim removed evidence | Bundle merge/runoff, `-t` / coverage trim |
| Filter: **junction / pairwise / dedup** | Right transcript beaten by wrong **cov ordering** | Compare `cov`/`longcov` vs C++ for same bundle |
| gffcompare **`j`** | Splice boundary shift | `-E`, ref assist, longtrim; not always cov—but cov errors can change which path wins |

### 4.2 “Emit StringTie would not” (false positive)

| Signal | Likely mechanism | Check first |
|--------|------------------|-------------|
| Many extra loci / isoforms | Low `-f`, extra path enumeration, weak flow paths kept | `transcript_isofrac`, `emit_junction_paths`, caps |
| **`c`** class | Containment / subset vs dominant | Pairwise + longunder + `dedup_intron_chain_subsets` |
| High query count, decent Sn | Predcluster keeps low-**cov** paths | **Coverage units** in predcluster vs rlink |

### 4.3 Preset tuning (`apply_stringtie_parity_preset`)

- **Default (balanced):** `-f` / pairwise / lowisofrac floor **0.02**, `max_junction_paths` cap **600**, no forced `filter_contained`, no automatic per-splice graph pruning.
- **Strict:** `RUSTLE_STRINGTIE_PRESET=strict` — floor **0.025**, junction cap **400**, `filter_contained`, per-splice **≥ 0.02** (precision sweep; Sn drops sharply on full BAM).

## 5. Work plan (ordered)

1. **Coverage invariants (highest leverage)**  
   - For a fixed bundle, use **`RUSTLE_BUNDLE_COV_DUMP`** + **`RUSTLE_TRACE_LOCUS`** (see §3.2) or `[BUNDLE_TRACE]` lines; diff **`coverage` / `longcov`** against C++ **`pred->cov` / `longcov`** on the same slice.

2. **Stage-by-stage predcluster parity**  
   - Already partially ported (`transcript_filter.rs` cites rlink line ranges). Re-verify **longunder** uses the same cov field as C++ **for long reads** (18991–18999 use `pred[n]->cov`, not `longcov`).

3. **Trace-driven fixes**  
   - Sort `trace_blockers_*.tsv` by frequency; for top blockers, add a **minimal replay** or unit test from saved bundle JSON if available.

4. **Anti-goals**  
   - Do not tune gffcompare thresholds until **cov + predcluster** match on shared traces; otherwise Pr/Sn moves are opaque.

## 6. Pointers into this repo

- **rlink:** `print_predcluster` ~18402+; longunder ~18975+; short recomputation ~19017+; `get_cov` ~1830+.  
- **Rustle:** `transcript_filter::print_predcluster_with_summary`; `path_extract` long-read cov ~6340+; `trace_reference` for blocker classification.  
- **Scripts:** `scripts/trace_replay_gap_audit.py`, `scripts/compare_trace_bootstrap_to_reference.py`, `scripts/run_full_bam_stringtie_trace.sh`.

This file is the living checklist for **“equivalent coverage throughout the pipeline”** and **correct chain choice** relative to StringTie’s C++ reference implementation.
