//! Parallel "StringTie-exact" mode scaffold.
//!
//! When `RUSTLE_STRINGTIE_EXACT=1` is set, Rustle disables every known
//! Rustle-specific relaxation at once and enables every identified
//! StringTie-parity gate. The goal is to see where Rustle's accumulated
//! divergence from StringTie's algorithm actually lives, instead of
//! chasing divergences one at a time.
//!
//! This is a BUILDING MODE: as new divergences are discovered via trace
//! comparison (e.g. `RUSTLE_PATHPAT_OR_TRACE=1 RUSTLE_TRACE_LOCUS=...`),
//! they should be gated on `stringtie_exact()` here and wired into the
//! relevant site. Over time the mode accumulates StringTie-parity until
//! it converges.
//!
//! # Usage
//!
//! ```bash
//! RUSTLE_STRINGTIE_EXACT=1 ./target/release/rustle -L input.bam -o out.gtf
//! ```
//!
//! ## Multi-stage alignment check (StringTie vs Rustle)
//!
//! After a paired run, compare instrumentation TSVs with:
//!
//! ```bash
//! PARITY_STAGE_TSV=st_parity.tsv ./stringtie -L in.bam -o st.gtf
//! ./target/release/rustle in.bam -o ru.gtf --debug-stage-tsv ru.stage.tsv
//! RUSTLE_PARITY_SHADOW=1 RUSTLE_PARITY_SHADOW_TSV=ru.shadow.tsv ...   # optional
//!
//! python3 scripts/parity_stage_alignment_check.py \
//!   --stringtie st_parity.tsv --rustle ru.stage.tsv \
//!   --parity-shadow ru.shadow.tsv --strict
//! ```
//!
//! Tunables: `--read-max-delta`, `--junction-tolerance`, shadow `max_allowed` rates inside the script.
//!
//! ## Partition **decision** dump (geometry signature)
//!
//! Same bundle key on both sides; one canonical line per bundle:
//!
//! - StringTie: `PARITY_PARTITION_TSV=st_partition.tsv` (segment coords normalized to Rustle space).
//! - Rustle: `RUSTLE_PARITY_PARTITION_TSV=ru_partition.tsv`
//!
//! ```bash
//! python3 scripts/compare_partition_geometry.py st_partition.tsv ru_partition.tsv
//! python3 scripts/compare_partition_geometry.py st.tsv ru.tsv --mismatch-tsv mismatches.tsv
//! ```
//!
//! Each segment’s **right** coordinate in the signature uses the same `-1` as StringTie’s
//! `parity_partition_emit_stringtie` (`cur->end - 1`); see `parity_partition_dump::encode_partition_signature`.
//! Chains are built **one per `SubBundleResult` linked list** (StringTie: one per `CBundle` walk),
//! not from the post-merge color component split, so `RUSTLE_PARITY_PARTITION_TSV` matches
//! `PARITY_PARTITION_TSV` bundle geometry on the same BAM (remaining diffs are mostly ±few bp
//! bundlenode boundaries, not missing chains).
//!
//! An **empty** Rustle signature means no bundlenode list was produced for that bundle (real
//! divergence vs StringTie’s partition); fixing that is follow-on assembly work, not the dump.
//!
//! ## Junction-set dump (multi-stage, `good_junc` alignment)
//!
//! StringTie’s de novo path applies **internal** steps that are hard to mirror line-for-line
//! (synthetic transfrags, replacing weak junctions, exon/intron extension or coverage trims).
//! Rustle still exposes **several bundle-local junction sets** so you can compare against
//! StringTie’s output **or** an earlier internal table once StringTie emits the same schema:
//!
//! - `RUSTLE_PARITY_JUNCTION_TSV=ru_junctions.tsv` — see `parity_junction_dump` for stages:
//!   `read_union`, `pre_canonical`, `graph_input`, `cgroup_feed`.
//!
//! ```bash
//! RUSTLE_PARITY_JUNCTION_TSV=ru_junctions.tsv ./target/release/rustle -L in.bam -o ru.gtf
//! python3 scripts/compare_junction_parity.py st_junctions.tsv ru_junctions.tsv \
//!   --ref-stage graph_input --cand-stage cgroup_feed
//! python3 scripts/compare_junction_parity.py a.tsv b.tsv --slop 2
//! python3 scripts/parity_fidelity_report.py --partition-st st_part.tsv --partition-ru ru_part.tsv \\
//!   --junction-st st_junc.tsv --junction-ru ru_junc.tsv --json parity_summary.json
//! ```
//!
//! When StringTie’s export corresponds to a checkpoint **before** final graph input, pick the
//! Rustle `stage` row that best matches (e.g. `pre_canonical` vs ST pre-canonicalize), and use
//! `--slop` for a few bp of splice-site leeway while bisecting algorithmic differences.
//!
//! ### StringTie instrumentation (`rlink.cpp`, same row schema for junctions)
//!
//! - `PARITY_JUNCTION_TSV=st_junctions.tsv` — `junction_set_v1` rows with `source=stringtie`.
//!   Bundle coordinates match `PARITY_PARTITION_TSV` / Rustle: **normalized** `start`/`end`
//!   (`BundleData::start-1`, `end-1`). Stages (in pipeline order):
//!   | `stage` | When |
//!   |---------|------|
//!   | `read_union_at_infer_entry` | First line of `infer_transcripts` (after `bundle_input` stage emit) |
//!   | `read_union_pre_count_good` | After LR exon trim + read sort, before `count_good_junctions` |
//!   | `read_union_post_count_good` | After `count_good_junctions` returns |
//!   | `aggregate_stranded_post_count_good` | Global junction table, `strand != 0`, index 0 skipped |
//!   | `aggregate_graph_ready_after_read_pass` | After per-read `good_junc` / strand pass in `build_graphs`, `mm >= 0`, `nreads_good >= junctionthr` or `guide_match` |
//!   | `aggregate_graph_ready_post_bundle_partition` | After `post_bundle_partition` / partition emit |
//! - `PARITY_BUNDLE_ORDER_TSV=st_order.tsv` — **global monotonic** `seq` across all bundles; one
//!   row per `parity_stage_emit` call (does **not** require `PARITY_STAGE_TSV`). Use to align
//!   processing order with Rustle logs when debugging multi-bundle races or batching.
//!
//! ```bash
//! PARITY_JUNCTION_TSV=st_junctions.tsv PARITY_BUNDLE_ORDER_TSV=st_order.tsv \
//!   PARITY_STAGE_TSV=st_stage.tsv PARITY_PARTITION_TSV=st_part.tsv \
//!   ./stringtie -L in.bam -o st.gtf
//! ```
//!
//! # Currently gated divergences
//!
//! | Site | Default behavior | StringTie-exact behavior |
//! |------|------------------|--------------------------|
//! | `junctions::apply_junction_filters_and_canonicalize` | Coalesce, nearby-weak, per-donor isofrac, canonicalize after `good_junc` | **Minimal pass**: drop `strand=0` kills + `filter_weak_junctions` only (matches StringTie shape: no extra merge/remap before the graph) |
//! | `RunConfig::apply_compat_preset` LR junction merge | `junction_canonical_tolerance` forced to ≥10bp | **Skipped** when `RUSTLE_STRINGTIE_EXACT=1`; then `apply_stringtie_exact_overrides` forces **0** |
//! | `pipeline` alt-junction demotion (`demote_alt_junctions_to_canonical`) | Optional merge of weak alt splice sites into canonical | **Skipped** (Rustle-only; no StringTie equivalent at this stage) |
//! | `path_extract::fwd_to_sink_fast_long` `past_seed` gate | Skip non-sink children when `i_coord > maxpath_coord` (Rustle-specific) | Gate disabled; StringTie has no equivalent |
//! | `path_extract::onpath_long` reach relaxation (line ~3142) | Allow tf past maxp when `node_can_reach` returns true | Strict: `trnode.last() > maxp` → false |
//! | `transfrag_process::process_transfrags` `novel_splice_rescue` | Promote contained long-read tfs with novel junctions to keeptrf | Disabled; StringTie strictly sets weak=1 |
//! | `path_extract::extract_transcripts` strict_redist filter | Already off-by-default (2026-04-21 flip) | Same — loose best_trf_match like StringTie |
//!
//! # Adding a new divergence gate
//!
//! 1. Identify a Rustle-specific deviation from StringTie's rlink.cpp behavior.
//! 2. Add an opt-in env flag for the specific site (e.g. `RUSTLE_FOO_OFF=1`).
//! 3. In the site, check `stringtie_exact() || specific_flag_set()`.
//! 4. Add an entry to the table above.
//! 5. Re-run with `RUSTLE_STRINGTIE_EXACT=1` to see the compound effect.
//!
//! ## Toward **100%** `partition_geometry` parity (same BAM, deterministic)
//!
//! After bundle-key alignment and segment end normalization, remaining mismatches are almost
//! always **different cgroup / bundlenode decompositions** (different `|` / `+` groupings) or
//! **±few bp exon boundaries** on the same intron graph. Closing the last ~20% is not a single
//! knob; it is a staged port of StringTie's **pre-graph** pipeline:
//!
//! 1. **Junction table identity** — `RUSTLE_PARITY_JUNCTION_TSV` multi-stage rows vs StringTie
//!    (compare `cgroup_feed` to ST post-`good_junc`, or an earlier ST stage to `pre_canonical` /
//!    `graph_input`; see `parity_junction_dump` module docs).
//! 2. **Read exon arrays** — byte-identical exon `(start,end)` lists per read after correction
//!    (junction snap, `junction_correction`, any trim that runs before cgroup).
//! 3. **`bundle_builder` vs `add_read_to_group`** — align `merge_fwd_groups`, `cgroup_keep_exon_*`,
//!    strand-colored union-find, and `merge_close_groups` with `rlink.cpp` / `bundle.cpp`.
//! 4. **3-strand projection** — confirm `SubBundleResult` strand filtering matches how StringTie
//!    fills `bundle[3]` / `bnode[3]` before `post_bundle_partition`.
//!
//! Use `compare_partition_geometry.py --mismatch-tsv` to drive bisection; add stage TSVs (like
//! `PARITY_STAGE_TSV`) for each checkpoint above so mismatches narrow to a single stage.
//!
//! ## cgroup / graph topology trace (Rustle-first)
//!
//! `RUSTLE_PARITY_TRACE_TSV=ru_trace.tsv` emits multi-kind rows (`cgroup_summary_v1`,
//! `graph_topology_v1` with stages `post_graph_build` / `pre_read_map`, `read_map_summary_v1`,
//! optional `read_exons_v1`) for cgroup counts, graph fingerprints before/after synthesis, mapped
//! read-transfrag counts, and per-read exon chains when `RUSTLE_PARITY_READ_EXONS=1`.
//! StringTie does not yet emit the same schema; add hooks in `rlink.cpp` for side-by-side graph parity.
//!
//! ```bash
//! RUSTLE_PARITY_TRACE_TSV=ru_trace.tsv RUSTLE_PARITY_PARTITION_TSV=ru_part.tsv \\
//!   RUSTLE_PARITY_JUNCTION_TSV=ru_junc.tsv ./target/release/rustle -L in.bam -o ru.gtf
//! PARITY_PARTITION_TSV=st_part.tsv PARITY_JUNCTION_TSV=st_junc.tsv stringtie -L in.bam -o st.gtf
//! python3 scripts/disambiguate_parity_bundles.py \\
//!   --partition-st st_part.tsv --partition-ru ru_part.tsv \\
//!   --junction-st st_junc.tsv --junction-ru ru_junc.tsv --rustle-trace ru_trace.tsv \\
//!   -o bundle_disambig.tsv
//! ```
//!
//! Event streams for near full-function tracing (high volume):
//! - StringTie: `PARITY_TRACE_TSV=st_trace.tsv PARITY_TRACE_LEVEL=1|2|3`
//! - Rustle: `RUSTLE_PARITY_TRACE_TSV=ru_trace.tsv RUSTLE_PARITY_TRACE_LEVEL=1|2|3`
//! and pass `--stringtie-trace st_trace.tsv` / `--rustle-trace ru_trace.tsv` to the
//! disambiguation script.

/// Returns true when the caller has opted into the full StringTie-exact mode.
/// Individual sites can still have their own finer-grained opt-out flags —
/// this is the meta-flag that enables everything in sync.
#[inline]
pub fn stringtie_exact() -> bool {
    std::env::var_os("RUSTLE_STRINGTIE_EXACT").is_some()
}

/// Returns true if StringTie-exact mode is active OR the site's specific
/// opt-out flag is set. Use at each divergence site to check whether the
/// Rustle-specific relaxation should be disabled.
#[inline]
pub fn parity_requested(specific_env_var: &str) -> bool {
    stringtie_exact() || std::env::var_os(specific_env_var).is_some()
}

/// One-time startup banner. Call from the top of the main pipeline entry
/// when stringtie_exact is active, so runs are clearly labeled.
pub fn maybe_emit_banner() {
    if stringtie_exact() {
        eprintln!(
            "[STRINGTIE_EXACT] RUSTLE_STRINGTIE_EXACT=1 — disabling Rustle-specific relaxations. \
             Gated sites: junction post-good_junc merge filters, alt-junction demotion, \
             fwd_to_sink past_seed gate, onpath_long reach relaxation, novel_splice_rescue; \
             RunConfig: junction_canonical_tolerance forced to 0 after compat preset. \
             See src/rustle/stringtie_parity.rs for details."
        );
    }
}
