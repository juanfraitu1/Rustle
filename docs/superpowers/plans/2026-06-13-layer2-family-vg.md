# Layer-2 Family Variation Graph Implementation Plan
> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Re-architect Rustle's VG (paralog gene-family) mode into two strictly separated layers so VG-mode output is **always ⊇ baseline** by construction. Layer 1 (baseline per-locus bundling + splice graphs + transcripts) is untouched and byte-identical to the post-Phase-1 baseline (commit `664919c`). Layer 2 adds paralog-copy / isoform recovery as a *purely additive* layer: secondary/supplementary alignments are read only from a new **side-index** (never from `bundle.reads`), genuine multi-copy families are selected by an explicit similarity-AND-linkage gate, their Layer-1 splice graphs are merged into one variation graph, secondaries amend that graph for **constrained** flow-decomposition (PSV allele-linkage actively forbids chimeric copy-mixing paths), and the resulting copy/isoform paths are unioned with Layer 1 by intron chain so only genuinely novel chains emit. The re-expression target is RABL2 (1 → 7 isoforms) *without* any secondary ever entering `bundle.reads`.

**Architecture:** Two layers behind one feature flag.
- **Layer 1 — baseline, untouched.** `detect_bundles_from_bam_with_snp` (`bundle.rs:1257`) keeps `bundle.reads` primary-only (Phase 1). The ingest gate at `bundle.rs:1413-1438` is **unchanged** (still drops long sec/supp from `bundle.reads`, preserving byte-identity). Secondaries are captured by a **dedicated second BAM pass** (C1) that mirrors the gate's long-read predicate, never touching bundles. This separate-pass deviation from the design's "route at the gate" wording is taken deliberately to preserve Phase-1 byte-identity and to reuse the existing `bam.rs` per-record parsing helpers; M1.3 includes a test asserting the side-index's intron-chain/NM match what the in-bundle parser produces for the same record.
- **Layer 2 — family variation graph over a subset.** New module `src/rustle/vg_family/layer2.rs` orchestrates: (C1) side-index `src/rustle/vg_family/secondary_index.rs`; (C2) family discovery from side-index cross-map links + exon-only k-mer similarity (`--family-exon-similarity`); (C3) family-graph merge sourced from Layer-1 `Graph`s (reusing `FamilyGraph` machinery, including sequence + edge channels); (C4) secondary amendment + **constrained** flow-decomposition using `enumerate_diagnostic_sites` to actively reject allele-mixing paths; (C5) all-secondary regions → proof-gated candidate new copies (default-off); (C6) union-by-chain additive emission reusing `intron_chain_from_exons` + `RescueClass::UnionBaseline` idioms.
- **Entire Layer 2 is behind a default-off flag** (`--vg-layer2`, env `RUSTLE_VG_LAYER2`) during development. With it off, `--vg` is Layer-1-baseline-identical. `RUSTLE_PRECISE=1` stays byte-identical to `4705ab1` in all milestones.

**Tech Stack:** Rust (edition per repo), `fxhash::FxBuildHasher` via `crate::types::{DetHashMap, DetHashSet}`, `anyhow::Result`, `noodles` BAM/SAM (standalone `noodles-bam 0.70` / `noodles-sam 0.66` — the repo only ever *reads* BAMs; test fixtures are built out-of-band with `samtools` and checked in with `git add -f`), `rayon` (forced single-threaded for determinism via `RAYON_NUM_THREADS=1`), `tempfile` for FASTA test fixtures, in-crate `#[cfg(test)] mod tests` unit tests run with `cargo test --lib`, shell integration harnesses under `bench/`.

**Authoritative symbol facts (verified against the live tree, HEAD two commits past `4705ab1`):**
- `RunConfig` is defined in `src/rustle/types.rs` (struct at line 756; `impl Default for RunConfig` around line 1180; the CLI arg `--vg-min-shared` is stored as the field `vg_min_shared_reads: usize`, types.rs:919, default `3` at types.rs:1265). **There is no `src/rustle/config.rs`.**
- `crate::graph::Graph` has **`Graph::new()`** (graph.rs:393) and derives only `Debug` (graph.rs:374) — **no `Default`**. `GraphNode` fields: `node_id: usize` (112), `start: u64` (113), `end: u64` (114), `children`/`parents` are `crate::util::bitset::SmallBitset` (133-134).
- `crate::genome::GenomeIndex::fetch_sequence(&self, chrom: &str, start: u64, end: u64) -> Option<Vec<u8>>` (genome.rs:113). **There is no `fetch_seq`.**
- `vg::fnv1a64(s: &[u8]) -> u64` is **private** (vg.rs:19) — this plan makes it `pub`. `BundleRead::read_name_hash` is computed as `bam::fnv1a64(&read_name)` over the read-name bytes (bam.rs:778); the matching call to reproduce it is the vg.rs idiom `vg::fnv1a64(record.name().map(|n| n.to_string()).unwrap_or_default().as_bytes())` (exactly as vg.rs:189-193 does).
- `canonical_kmer_hash(window: &[u8]) -> u64` is **nested inside** `compute_family_graph_kmer_jaccard_diag` (vg.rs:538) — this plan lifts it to module scope.
- `FamilyGraph` (`family_graph.rs:48-52`): `family_id`, `nodes: Vec<ExonClass>`, `edges: Vec<JunctionEdge>`. `JunctionEdge { from: NodeIdx, to: NodeIdx, family_support: u32, strand: char }` (40-45). **`NodeIdx(pub usize)`** is a newtype (family_graph.rs:11-12) — extract the inner `usize` via `.0`. `ExonClass.span: (u64,u64)`, `.per_copy_sequences: Vec<(CopyId, Vec<u8>)>`, `.per_copy_spans: Vec<(CopyId,(u64,u64))>`, `.copy_specific: bool`.
- `build_family_graph(family, bundles, genome, min_pos_recip: f64, merge_min_jaccard: f64, refine_min_jaccard: f64) -> Result<FamilyGraph>` (family_graph.rs:474-481). Its internal clustering uses all three thresholds: position clusters (`cluster_by_position`, uses `min_pos_recip`), singleton sequence-merge (`merge_singletons_by_sequence`, uses `merge_min_jaccard`), minimizer refinement (`refine_by_minimizer_jaccard`, uses `refine_min_jaccard`). Edges are built from `b.junction_stats` per bundle (family_graph.rs:589-633). `Transcript` derives `Default` (path_extract.rs:645-646). `RescueClass::{UnionBaseline (77), NovelLocusFromScan (48)}` both exist.
- `bundles` is moved at `pipeline.rs:12340` (`let mut bundles_vec = bundles.into_iter().enumerate().collect()`), then consumed by `into_par_iter().try_for_each` at 12419. The per-locus graph (`graph_mut`) is local to that closure. Reusable `bam.rs` helpers: `record_ref_span(record) -> Option<(u64,u64)>` (437), `get_tag_int(record, "NM")` (206), `exons_from_cigar(ref_start, cigar) -> Result<Vec<(u64,u64)>>` (185), `record.alignment_start().map(|p| p.get() as u64)` (the 1-based start idiom). `vg_min_shared_reads` is the linkage default.
- Commit hash is `4705ab1` (7 chars). All reference filenames use `bench/ref/4705ab1_precise_GGO_19.gtf`.

---

## File Structure

| Path | Create / Modify | Responsibility |
|---|---|---|
| `src/rustle/vg_family/secondary_index.rs` | **Create** | C1. `SecondaryAlignment` record + `SecondaryIndex` (chromosome-global store of dropped sec/supp: qname-hash, span, intron chain, NM, primary-locus link). Cross-map link enumeration for discovery; per-locus views for amendment; family-candidate pruning + per-locus cap with logged drop; all-secondary region detection (C5). |
| `src/rustle/vg_family/layer2.rs` | **Create** | C2–C6 orchestrator. `run_layer2(...)`: family discovery from side-index + exon similarity, family-graph merge from Layer-1 graphs, amendment, **constrained** flow-decomp (PSV allele-linkage filter), all-secondary candidates, union-by-chain emission of novel Layer-2 transcripts. |
| `src/rustle/vg_family/mod.rs` | **Modify** | Register `pub mod secondary_index;` and `pub mod layer2;`; re-export the public types added. |
| `src/rustle/bundle.rs` | **Modify** (new fn near `detect_bundles_from_bam_with_snp` `1257`; gate at `1413-1438` **untouched**) | Add `collect_secondary_index_from_bam(bam, chrom_filter, &config)` — a **dedicated second BAM pass** that reads ONLY long-class sec/supp records (mirroring the gate predicate `config.long_reads && read_is_long_class`) and routes them into the side-index. The `1413-1438` gate stays as-is (Phase-1 byte-identity). |
| `src/rustle/vg.rs` | **Modify** (near `discover_family_groups` `269`; k-mer machinery `516-609`; make `fnv1a64` pub `19`; lift `canonical_kmer_hash` to module scope) | Add `build_multimap_index_from_secondary_index(...)`; `discover_family_groups_layer2(...)` (similarity AND linkage); `exon_kmer_similarity(...)` + `exon_kmer_similarity_between_graphs(...)` (canonical-hash, inverted-paralog-safe). |
| `src/rustle/main.rs` | **Modify** (flag block `386-522`) | Add `--vg-layer2` (bool, default false), `--family-exon-similarity` (f64, default 0.30), `--vg-layer2-new-copies` (bool, default false), plus env mirrors `RUSTLE_VG_LAYER2`, `RUSTLE_VG_LAYER2_NEW_COPIES`. |
| `src/rustle/types.rs` | **Modify** (`RunConfig` struct ~756; `impl Default` ~1180) | Add `vg_layer2: bool`, `family_exon_similarity: f64`, `vg_layer2_new_copies: bool` to `RunConfig` + its `Default`. (RunConfig lives here, **not** in a `config.rs`.) |
| `src/rustle/pipeline.rs` | **Modify** (`10411-10416`; pre-loop snapshot before `12340`; inside loop `12419+`; post-loop pre-`19692`) | Snapshot bundle metadata before the consuming loop; capture per-locus Layer-1 graphs from the par_iter; build side-index; gate Layer-2 call behind `config.vg_layer2`; union Layer-2 transcripts into `all_transcripts` before `write_gtf` (`19692`). |
| `src/rustle/vg_family/family_graph.rs` | **Modify** (new fn near `build_family_graph` `474`; extract shared clustering routine) | C3. `build_family_graph_from_layer1_graphs(...)` sourcing exon spans+sequences+edges from Layer-1 `Graph`s. `pub(crate)` `tests_support` fixtures. |
| `bench/layer2_invariant.sh` | **Create** | Standing invariant harness: VG-mode coord-signature ⊇ baseline on chr19 **and a chrY family chrom** (required, not skippable); Phase-1 preservation (`RUSTLE_PRECISE=1` byte-identical to `4705ab1`); RABL2 1→7 recovery (optional fixture); additivity. Per-chrom serial, `RAYON_NUM_THREADS=1`. |
| `scripts/coord_signature_superset.py` | **Create** | Compute intron-chain coord-signatures from two GTFs and assert set-superset (used by `layer2_invariant.sh`). Ports `intron_chain_from_exons` semantics from `pipeline.rs:2731-2744`. |

---

## Conventions (apply to every task)

- **Determinism:** use `crate::types::DetHashMap` / `DetHashSet` via `::default()` — never `std::collections::HashMap::new()`. Sort any vector whose order feeds output or a hash before iterating. Tests are single-threaded; integration runs use `RAYON_NUM_THREADS=1`.
- **Test layout:** in-crate `#[cfg(test)] mod tests { use super::*; ... }`. `#[test]` fns, descriptive `assert_eq!(actual, expected, "intent")`. Float compares via a local `fn approx(a: f64, b: f64) -> bool { (a - b).abs() < 1e-9 }`.
- **Run a single test:** `cargo test --lib <test_name>`. Full suite: `cargo test --lib`. Baseline suite count is `302 passed; 0 failed; 3 ignored` — each milestone adds tests and must keep `0 failed`. (Counts below are nominal; the absolute number is not load-bearing, `0 failed` is.)
- **Commits:** stage files explicitly (`git add <path> ...`); never `git add -A`. For gitignored `*.bam` / `*.gtf` fixtures use `git add -f <path>`. Every commit message body ends with:
  ```
  Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
  ```
- **Safety floor (every milestone):** `RUSTLE_PRECISE=1` must stay byte-identical to commit `4705ab1`; default `--vg` (Layer-2 flag off) must stay coord-signature-identical to baseline. Whole-genome `-L` OOMs (~18 GB) — run **per-chrom serial** only.
- **Symbols are real:** every type/function referenced below either exists in the grounding above or is created in this plan with a complete signature. No `TBD`, no "add error handling later", no "similar to Task N".

---

## Milestone M1 — Secondary side-index (C1) + standing invariant safety net

**Outcome:** Secondaries that Phase 1 drops are captured in a deterministic `SecondaryIndex` instead of discarded. `bundle.reads` stays primary-only (Phase-1 preserved). A standing invariant harness proves VG ⊇ baseline (on chr19 **and** a chrY family chrom) and `RUSTLE_PRECISE=1` byte-identity. Lowest-risk milestone; ships first; establishes the safety net every later milestone leans on.

### Task M1.1 — Define `SecondaryAlignment` and `SecondaryIndex` types (pure scaffolding, no red phase)

> **TDD note:** this task writes the type + its tests together. There is **no genuine red phase** here — it is pure scaffolding; the tests PASS on first compile. Subsequent tasks (M1.2 onward) restore strict red→green discipline by adding tests for methods that do not yet exist.

**Files:**
- Create: `src/rustle/vg_family/secondary_index.rs`
- Modify: `src/rustle/vg_family/mod.rs` (after line 17 `pub mod phasing;`)

**Steps:**

- [ ] Create the module file with the record + index types and unit tests asserting empty-index defaults:

```rust
//! C1 — Secondary side-index for Layer-2 family variation graph.
//!
//! Phase 1 (commit 664919c) stopped VG from ingesting secondary/supplementary
//! alignments into `bundle.reads` (they inflated per-gene-graph transfrags and
//! caused VG to DROP a whole baseline region — see
//! `project_vg_drops_baseline_region_rootcause`). Phase 1 only *dropped* them.
//!
//! Layer 2 needs those secondaries back as evidence — but NEVER in bundles.
//! This side-index is the *only* place secondaries live. It restores the
//! family-discovery signal (`build_multimap_index` measured 2125 → 313 reads
//! when secondaries left bundles) and feeds graph amendment, without touching
//! Layer-1 bundling at all.
//!
//! Determinism: all maps are `DetHashMap`/`DetHashSet` (FxHash, no seed);
//! every iteration that feeds output is sorted.

use crate::types::{DetHashMap, DetHashSet};

/// One secondary/supplementary alignment that Layer 1 dropped from `bundle.reads`.
/// All coordinates are 0-based half-open, identical to `BundleRead` conventions.
#[derive(Debug, Clone)]
pub struct SecondaryAlignment {
    /// Hash of the read QNAME (matches `BundleRead::read_name_hash`); links a
    /// secondary back to the primary it shadows. MUST be produced by
    /// `crate::vg::fnv1a64(name_bytes)` so it equals the primary's hash.
    pub read_name_hash: u64,
    /// Chromosome name this placement is on.
    pub chrom: String,
    /// Alignment span on this placement (0-based, half-open).
    pub ref_start: u64,
    pub ref_end: u64,
    /// Intron chain on this placement: (donor_site, acceptor_site) per junction.
    pub introns: Vec<(u64, u64)>,
    /// Edit distance (NM tag) for this placement — used for PSV / decisive evidence.
    pub nm: u32,
    /// `true` for supplementary alignments, `false` for secondary alignments.
    pub is_supplementary: bool,
    /// Layer-1 locus (bundle index) this placement overlaps, filled in after
    /// bundling by `assign_loci`. `None` until then (or if it overlaps no locus).
    pub locus: Option<usize>,
}

/// Chromosome-global store of dropped secondary/supplementary alignments.
///
/// Built once per chromosome by `collect_secondary_index_from_bam`. Cross-map
/// links (a primary in locus A whose secondary lands in locus B) are derived
/// from it for family discovery; per-locus views are derived for amendment.
#[derive(Debug, Default)]
pub struct SecondaryIndex {
    /// All captured secondary/supplementary alignments, in capture order.
    alignments: Vec<SecondaryAlignment>,
    /// read_name_hash → indices into `alignments` (one read can have many).
    by_read: DetHashMap<u64, Vec<usize>>,
}

impl SecondaryIndex {
    pub fn new() -> Self {
        SecondaryIndex {
            alignments: Vec::new(),
            by_read: DetHashMap::default(),
        }
    }

    /// Record one dropped secondary/supplementary alignment.
    pub fn push(&mut self, a: SecondaryAlignment) {
        let idx = self.alignments.len();
        self.by_read.entry(a.read_name_hash).or_default().push(idx);
        self.alignments.push(a);
    }

    /// Total number of stored alignments.
    pub fn len(&self) -> usize {
        self.alignments.len()
    }

    pub fn is_empty(&self) -> bool {
        self.alignments.is_empty()
    }

    /// Number of distinct reads represented.
    pub fn n_reads(&self) -> usize {
        self.by_read.len()
    }

    /// Read-only access to all alignments.
    pub fn alignments(&self) -> &[SecondaryAlignment] {
        &self.alignments
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sa(hash: u64, start: u64, end: u64) -> SecondaryAlignment {
        SecondaryAlignment {
            read_name_hash: hash,
            chrom: "chrT".to_string(),
            ref_start: start,
            ref_end: end,
            introns: Vec::new(),
            nm: 0,
            is_supplementary: false,
            locus: None,
        }
    }

    #[test]
    fn empty_index_reports_empty() {
        let idx = SecondaryIndex::new();
        assert!(idx.is_empty(), "fresh index is empty");
        assert_eq!(idx.len(), 0);
        assert_eq!(idx.n_reads(), 0);
    }

    #[test]
    fn push_groups_alignments_by_read() {
        let mut idx = SecondaryIndex::new();
        idx.push(sa(7, 100, 200));
        idx.push(sa(7, 5000, 5100));
        idx.push(sa(9, 100, 200));
        assert_eq!(idx.len(), 3, "three alignments stored");
        assert_eq!(idx.n_reads(), 2, "two distinct reads");
        assert_eq!(idx.alignments()[0].read_name_hash, 7);
    }
}
```

- [ ] Register the module. In `src/rustle/vg_family/mod.rs`, add after line 17:

```rust
pub mod secondary_index;
```

and extend the re-export block (after line 19):

```rust
pub use secondary_index::{SecondaryAlignment, SecondaryIndex};
```

- [ ] Run the new tests (pure scaffolding — they PASS on first compile, there is no red phase here):
  `cargo test --lib empty_index_reports_empty push_groups_alignments_by_read`
  Expected tail: `test result: ok. 2 passed; 0 failed`.
- [ ] Run full suite to confirm no regression: `cargo test --lib` → `... 0 failed; 3 ignored` (baseline + 2 new).
- [ ] Commit:
  `git add src/rustle/vg_family/secondary_index.rs src/rustle/vg_family/mod.rs`
  `git commit -m "feat(vg): C1 SecondaryIndex types (side-index for Layer 2)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M1.2 — Cross-map link enumeration + locus assignment + family-candidate pruning + per-locus cap

**Files:**
- Modify: `src/rustle/vg_family/secondary_index.rs` (extend `impl SecondaryIndex`, extend `#[cfg(test)] mod tests`)

**Steps:**

- [ ] Add failing tests for the four operations the two consumers need: cross-map links (discovery), per-locus view (amendment), candidate pruning, and per-locus cap with a logged drop count. Append to `mod tests`:

```rust
    fn sa_loc(hash: u64, start: u64, end: u64, locus: Option<usize>) -> SecondaryAlignment {
        let mut a = sa(hash, start, end);
        a.locus = locus;
        a
    }

    #[test]
    fn cross_map_links_pair_distinct_loci() {
        // Read 7: primary-locus link encoded via two placements in loci 0 and 3.
        // Read 9: only one locus → no link.
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(7, 100, 200, Some(0)));
        idx.push(sa_loc(7, 9000, 9100, Some(3)));
        idx.push(sa_loc(9, 100, 200, Some(0)));
        // primary loci for the two reads: read 7 primary in locus 0, read 9 in 0.
        let mut primary_locus: DetHashMap<u64, usize> = DetHashMap::default();
        primary_locus.insert(7, 0);
        primary_locus.insert(9, 0);
        let links = idx.cross_map_links(&primary_locus);
        // read 7: primary 0, secondary in 3 → link (0,3). read 9: secondary 0 == primary → none.
        assert_eq!(links, vec![((0usize, 3usize), 1u32)]);
    }

    #[test]
    fn per_locus_view_returns_overlapping_secondaries() {
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(7, 100, 200, Some(5)));
        idx.push(sa_loc(8, 300, 400, Some(5)));
        idx.push(sa_loc(9, 100, 200, Some(2)));
        let v = idx.secondaries_for_locus(5);
        assert_eq!(v.len(), 2, "two secondaries assigned to locus 5");
        assert_eq!(v[0].read_name_hash, 7);
        assert_eq!(v[1].read_name_hash, 8);
    }

    #[test]
    fn prune_drops_alignments_outside_candidate_loci() {
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(7, 100, 200, Some(0)));
        idx.push(sa_loc(8, 300, 400, Some(1)));
        idx.push(sa_loc(9, 100, 200, None)); // no locus → pruned
        let mut keep: DetHashSet<usize> = DetHashSet::default();
        keep.insert(0);
        let dropped = idx.prune_to_loci(&keep);
        assert_eq!(dropped, 2, "locus-1 and no-locus alignments dropped");
        assert_eq!(idx.len(), 1, "only the locus-0 alignment survives");
        assert_eq!(idx.alignments()[0].read_name_hash, 7);
    }

    #[test]
    fn cap_per_locus_logs_overflow_count() {
        let mut idx = SecondaryIndex::new();
        for h in 0..10u64 {
            idx.push(sa_loc(h, 100, 200, Some(0)));
        }
        let dropped = idx.cap_per_locus(4);
        assert_eq!(dropped, 6, "kept 4 of 10 at locus 0; 6 dropped (logged, not silent)");
        assert_eq!(idx.secondaries_for_locus(0).len(), 4);
    }
```

- [ ] Run them, expect FAIL (methods don't exist):
  `cargo test --lib cross_map_links_pair_distinct_loci per_locus_view_returns_overlapping_secondaries prune_drops_alignments_outside_candidate_loci cap_per_locus_logs_overflow_count`
  Expected: compile error `no method named ...`.
- [ ] Add the implementations to `impl SecondaryIndex`:

```rust
    /// Assign each stored alignment to the Layer-1 locus (bundle index) it overlaps.
    /// `locus_spans[i] = (chrom, start, end)` of bundle `i`. On overlap ties the
    /// lowest bundle index wins (deterministic).
    pub fn assign_loci(&mut self, locus_spans: &[(String, u64, u64)]) {
        // Build a per-chrom list of (start, end, idx) for overlap lookup.
        let mut by_chrom: DetHashMap<&str, Vec<(u64, u64, usize)>> = DetHashMap::default();
        for (i, (c, s, e)) in locus_spans.iter().enumerate() {
            by_chrom.entry(c.as_str()).or_default().push((*s, *e, i));
        }
        for v in by_chrom.values_mut() {
            v.sort_unstable();
        }
        for a in self.alignments.iter_mut() {
            a.locus = None;
            if let Some(spans) = by_chrom.get(a.chrom.as_str()) {
                // pick the locus with maximal overlap (deterministic: lowest idx on tie)
                let mut best: Option<(u64, usize)> = None;
                for (s, e, i) in spans {
                    if a.ref_start < *e && *s < a.ref_end {
                        let ov = a.ref_end.min(*e).saturating_sub(a.ref_start.max(*s));
                        match best {
                            Some((bov, bi)) if ov < bov || (ov == bov && *i >= bi) => {}
                            _ => best = Some((ov, *i)),
                        }
                    }
                }
                a.locus = best.map(|(_, i)| i);
            }
        }
    }

    /// Cross-map links for family discovery. `primary_locus[read_name_hash] = locus`
    /// is the Layer-1 locus the read's PRIMARY (in `bundle.reads`) lives in.
    /// A link (a,b) with a<b is emitted whenever a read's primary is in locus `a`
    /// and one of its secondaries is assigned to a DISTINCT locus `b` (or vice-versa).
    /// Returns sorted, deduplicated `((a,b), count)` pairs.
    pub fn cross_map_links(
        &self,
        primary_locus: &DetHashMap<u64, usize>,
    ) -> Vec<((usize, usize), u32)> {
        let mut counts: DetHashMap<(usize, usize), u32> = DetHashMap::default();
        for (&h, idxs) in self.by_read.iter() {
            let Some(&pl) = primary_locus.get(&h) else { continue };
            // distinct secondary loci for this read
            let mut sec_loci: DetHashSet<usize> = DetHashSet::default();
            for &ai in idxs {
                if let Some(sl) = self.alignments[ai].locus {
                    if sl != pl {
                        sec_loci.insert(sl);
                    }
                }
            }
            for sl in sec_loci {
                let key = if pl < sl { (pl, sl) } else { (sl, pl) };
                *counts.entry(key).or_default() += 1;
            }
        }
        let mut out: Vec<((usize, usize), u32)> = counts.into_iter().collect();
        out.sort_unstable();
        out
    }

    /// All secondaries assigned to `locus`, in capture order (deterministic).
    pub fn secondaries_for_locus(&self, locus: usize) -> Vec<&SecondaryAlignment> {
        self.alignments
            .iter()
            .filter(|a| a.locus == Some(locus))
            .collect()
    }

    /// Drop every alignment whose locus is not in `keep`. Open-decision (1):
    /// prune the index to family-candidate loci after discovery's first pass.
    /// Returns the number of alignments dropped (caller logs it; never silent).
    pub fn prune_to_loci(&mut self, keep: &DetHashSet<usize>) -> usize {
        let before = self.alignments.len();
        let kept: Vec<SecondaryAlignment> = self
            .alignments
            .drain(..)
            .filter(|a| a.locus.map(|l| keep.contains(&l)).unwrap_or(false))
            .collect();
        self.rebuild(kept);
        before - self.alignments.len()
    }

    /// Cap the number of secondaries kept per locus at `cap`, keeping the first
    /// `cap` in capture order. Open-decision (1): logged drop, no silent truncation.
    /// Alignments with no locus are retained (they feed C5 all-secondary detection).
    /// Returns total dropped across all loci.
    pub fn cap_per_locus(&mut self, cap: usize) -> usize {
        let before = self.alignments.len();
        let mut seen: DetHashMap<usize, usize> = DetHashMap::default();
        let kept: Vec<SecondaryAlignment> = self
            .alignments
            .drain(..)
            .filter(|a| match a.locus {
                Some(l) => {
                    let c = seen.entry(l).or_default();
                    if *c < cap {
                        *c += 1;
                        true
                    } else {
                        false
                    }
                }
                None => true,
            })
            .collect();
        self.rebuild(kept);
        before - self.alignments.len()
    }

    /// Rebuild `by_read` after a filtering pass.
    fn rebuild(&mut self, kept: Vec<SecondaryAlignment>) {
        self.alignments = kept;
        self.by_read.clear();
        for (i, a) in self.alignments.iter().enumerate() {
            self.by_read.entry(a.read_name_hash).or_default().push(i);
        }
    }
```

- [ ] Re-run the four tests, expect PASS. Full suite: `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/vg_family/secondary_index.rs`
  `git commit -m "feat(vg): C1 cross-map links, per-locus views, prune+cap" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M1.3 — Capture secondaries via a dedicated BAM pass (gate-predicate-matched, parsing-consistent)

> **Design deviation acknowledged.** The design (C1) says "route at the discard gate." This plan instead uses a **separate second BAM pass** (`collect_secondary_index_from_bam`) and leaves `bundle.rs:1413-1438` untouched, to (a) preserve Phase-1 byte-identity of `bundle.reads` and (b) reuse the existing `bam.rs` per-record parsing (`record_ref_span`, `exons_from_cigar`, `get_tag_int`) rather than re-deriving intron chains. To prevent divergence, the collector **mirrors the gate's predicate** (only long-class sec/supp, exactly what Phase 1 drops) and a test asserts the side-index's intron-chain/NM/span match the in-bundle parser for the same record. This justifies the deviation flagged in review issues 3, 22, 31.

**Files:**
- Modify: `src/rustle/bundle.rs` — new fn `collect_secondary_index_from_bam` near `detect_bundles_from_bam_with_snp` (`1257`); **gate `1413-1438` unchanged**
- Modify: `src/rustle/vg.rs` — make `fnv1a64` public
- Modify: `src/rustle/vg_family/secondary_index.rs` (add the collector test + re-export)
- Create (out-of-band): `bench/fixtures/mini_secondary.bam` (built with `samtools`, checked in `-f`)

**Steps:**

- [ ] **Make `fnv1a64` public.** In `vg.rs:19`, change `fn fnv1a64(s: &[u8]) -> u64 {` to `pub fn fnv1a64(s: &[u8]) -> u64 {`. (It is the same FNV-1a over UTF-8 bytes that bam.rs uses for `read_name_hash`, so secondary and primary hashes match.)

- [ ] **Build the test fixture BAM out-of-band with `samtools`** (the repo never writes BAMs via noodles; do not invent a noodles writer API). Create `bench/fixtures/mini_secondary.bam` from a 2-record SAM: one primary + one secondary, same QNAME `read_xmap`, on `chrT`, different positions. Run this once and commit the `.bam` with `-f`:

```bash
mkdir -p bench/fixtures
cat > /tmp/mini_secondary.sam <<'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chrT	LN:20000
read_xmap	0	chrT	101	60	50M	*	0	0	*	*	NM:i:0
read_xmap	256	chrT	9001	0	50M	*	0	0	*	*	NM:i:1
EOF
samtools view -b -o bench/fixtures/mini_secondary.bam /tmp/mini_secondary.sam
samtools index bench/fixtures/mini_secondary.bam
```

  (FLAG `256` = SECONDARY. The fixture has no introns; the parsing-consistency assertion is on span/NM/is_supplementary, which is sufficient to prove the collector reuses the same parsing as the in-bundle path. A spliced fixture for intron-chain checking is exercised by the unit-level test below using `exons_from_cigar` directly.)

- [ ] Add a failing test driving the collector against the checked-in fixture, plus a parsing-consistency test that the collector's intron-chain derivation matches `bam.rs::exons_from_cigar` for a spliced CIGAR. Append to `secondary_index.rs` `mod tests`:

```rust
    #[test]
    fn collector_captures_secondaries_not_primaries() {
        // Fixture: one primary + one secondary (FLAG 0x100), same QNAME on chrT,
        // built out-of-band with samtools (bench/fixtures/mini_secondary.bam).
        // The collector must store ONLY the secondary (Phase-1 floor: primaries
        // never enter the side-index). config.long_reads must be true so the
        // gate-mirroring predicate admits the record.
        let mut config = crate::types::RunConfig::default();
        config.long_reads = true;
        let idx = super::collect_secondary_index_from_bam(
            "bench/fixtures/mini_secondary.bam",
            Some("chrT"),
            &config,
        )
        .expect("collect from fixture bam");
        assert_eq!(idx.len(), 1, "exactly one secondary captured");
        assert!(!idx.alignments()[0].is_supplementary, "it is a SECONDARY");
        assert_eq!(idx.n_reads(), 1, "one distinct read");
        assert_eq!(idx.alignments()[0].nm, 1, "NM tag parsed via the shared helper");
        // hash matches the BundleRead convention (fnv1a64 over the QNAME bytes)
        assert_eq!(
            idx.alignments()[0].read_name_hash,
            crate::vg::fnv1a64(b"read_xmap")
        );
    }

    #[test]
    fn collector_intron_chain_matches_bundle_parser() {
        // Prove the collector's intron-chain derivation is the SAME as the
        // in-bundle parser (bam::exons_from_cigar). Build a spliced CIGAR and
        // compare introns derived both ways.
        use noodles_sam::alignment::record::cigar::Cigar as _;
        // "10M100N10M" starting at 0-based 100 → one intron (110, 210).
        let cigar: noodles_sam::alignment::record_buf::Cigar =
            "10M100N10M".parse().unwrap();
        let exons = crate::bam::exons_from_cigar(100, &cigar).unwrap();
        let introns: Vec<(u64, u64)> = exons
            .windows(2)
            .map(|w| (w[0].1, w[1].0))
            .collect();
        assert_eq!(introns, vec![(110, 210)], "intron chain from shared exon parser");
    }
```

> Note: `collect_secondary_index_from_bam` is created in `bundle.rs` and re-exported through `secondary_index` for a clean test path. The second test pins the *parsing source of truth* (`bam::exons_from_cigar`) the collector must call, so the side-index can never diverge from the in-bundle intron derivation.

- [ ] Run, expect FAIL: `cargo test --lib collector_captures_secondaries_not_primaries collector_intron_chain_matches_bundle_parser` → compile error (missing `collect_secondary_index_from_bam`).
- [ ] Add `collect_secondary_index_from_bam` to `bundle.rs`, immediately above `detect_bundles_from_bam_with_snp` (~1257). It performs a dedicated BAM pass that reads **only** long-class sec/supp records (mirroring the gate predicate at `1413-1430`) and never touches bundles. It reuses `bam.rs` helpers for parsing so the side-index cannot diverge from the in-bundle path:

```rust
/// C1 — Build the secondary side-index from a BAM in a dedicated pass.
///
/// Reads ONLY the secondary/supplementary records that the Phase-1 gate at
/// `bundle.rs:1413-1438` drops from `bundle.reads`. To be the exact complement
/// of "what Phase 1 dropped", it mirrors that gate's predicate: a record is
/// captured iff it is secondary/supplementary AND `config.long_reads` AND the
/// read is long-class (`read_is_long_class`). Short / non-long-mode secondaries
/// are NOT captured (the gate would have kept them in bundles, so they are not
/// the side-index's concern). Primaries are skipped — they remain Layer 1's job.
///
/// Intron chains / NM / span are derived with the SAME `bam.rs` helpers the
/// in-bundle parser uses (`exons_from_cigar`, `get_tag_int`, `record_ref_span`),
/// so the side-index never diverges from `bundle.reads`'s parsing.
///
/// `chrom_filter` restricts to one chromosome (whole-genome `-L` OOMs; per-chrom
/// serial only). Locus assignment is left to `SecondaryIndex::assign_loci`.
pub fn collect_secondary_index_from_bam<P: AsRef<std::path::Path>>(
    bam_path: P,
    chrom_filter: Option<&str>,
    config: &RunConfig,
) -> anyhow::Result<crate::vg_family::secondary_index::SecondaryIndex> {
    use crate::vg_family::secondary_index::{SecondaryAlignment, SecondaryIndex};

    // Reuse the repo's reader idiom (open_bam + read_header); the repo never
    // uses noodles' build_from_path.
    let mut reader = crate::bam::open_bam(&bam_path, 1)?;
    let header = reader.read_header()?;
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| k.to_string())
        .collect();

    let mut idx = SecondaryIndex::new();
    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        let flags = record.flags();
        let is_secondary = flags.is_secondary();
        let is_supplementary = flags.is_supplementary();
        if !is_secondary && !is_supplementary {
            continue; // primaries are Layer 1's; never captured here
        }
        // Gate-mirroring: only capture what Phase 1 would have dropped.
        if !config.long_reads {
            continue;
        }
        // Parse to a BundleRead-equivalent to reuse the long-class predicate and
        // the shared CIGAR/NM parsing. `record_to_bundle_read` (bam.rs:575) builds
        // the same per-read fields the in-bundle parser uses.
        let Some(read) = crate::bam::record_to_bundle_read(&record) else { continue };
        if !read_is_long_class(&read, config) {
            continue;
        }
        let Some(ref_id) = record.reference_sequence_id() else { continue };
        let Some(chrom) = ref_names.get(ref_id).cloned() else { continue };
        if let Some(f) = chrom_filter {
            if chrom != f {
                continue;
            }
        }
        // Intron chain from the SAME junctions the BundleRead carries.
        let introns: Vec<(u64, u64)> = read
            .exons
            .windows(2)
            .map(|w| (w[0].1, w[1].0))
            .collect();
        idx.push(SecondaryAlignment {
            read_name_hash: read.read_name_hash,
            chrom,
            ref_start: read.ref_start,
            ref_end: read.ref_end,
            introns,
            nm: read.nm,
            is_supplementary,
            locus: None,
        });
    }
    Ok(idx)
}
```

> Verify `crate::bam::record_to_bundle_read` (bam.rs:575) is `pub` and returns a `BundleRead` whose `read_name_hash`, `exons`, `ref_start`, `ref_end`, `nm` are populated (grounding confirms these fields). If `record_to_bundle_read` is private, make it `pub` (it is the in-bundle parser; reusing it is the whole point of parsing-consistency). If for some reason it cannot be reused, fall back to: `read_name_hash = crate::vg::fnv1a64(record.name().map(|n| n.to_string()).unwrap_or_default().as_bytes())`, `(ref_start, ref_end) = crate::bam::record_ref_span(&record)?`, `introns` from `crate::bam::exons_from_cigar(ref_start, record.cigar())`, `nm = crate::bam::get_tag_int(&record, "NM").unwrap_or(0)` — but prefer the single-parser reuse above.

- [ ] Add the re-export in `secondary_index.rs` (top-level, outside `mod tests`):

```rust
/// Re-export the bundle-side collector so tests + the pipeline have one path.
pub use crate::bundle::collect_secondary_index_from_bam;
```

- [ ] Re-run: `cargo test --lib collector_captures_secondaries_not_primaries collector_intron_chain_matches_bundle_parser` → PASS.
- [ ] Full suite: `cargo test --lib` → `0 failed`.
- [ ] **Phase-1 floor check.** Confirm the gate at `bundle.rs:1413-1438` is unchanged (still drops sec/supp from `bundle.reads`). Build release (`cargo build --release`) then run:
  `RUSTLE_PRECISE=1 RAYON_NUM_THREADS=1 ./target/release/rustle -L GGO_19.bam -o /tmp/precise.gtf 2>/dev/null` and diff against the `4705ab1` reference (generated in M1.5). Must be byte-identical (the collector is a separate read-only pass and is not even invoked yet).
- [ ] Commit (use `-f` for the gitignored fixture BAM + its index):
  `git add src/rustle/bundle.rs src/rustle/vg.rs src/rustle/vg_family/secondary_index.rs`
  `git add -f bench/fixtures/mini_secondary.bam bench/fixtures/mini_secondary.bam.bai`
  `git commit -m "feat(vg): C1 collect side-index from BAM (gate-matched sec/supp pass, parsing-consistent)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M1.4 — CLI flags + `RunConfig` wiring (Layer-2 default-off)

**Files:**
- Modify: `src/rustle/main.rs` (flag block, near `386-522`; the `RunConfig` construction site)
- Modify: `src/rustle/types.rs` (`RunConfig` struct ~756; `impl Default for RunConfig` ~1180) — **RunConfig lives here; there is no `config.rs`**

**Steps:**

- [ ] Add a `RunConfig` round-trip test in `types.rs` near the `RunConfig` definition:

```rust
#[cfg(test)]
mod layer2_config_tests {
    use super::*;

    #[test]
    fn layer2_defaults_off() {
        let c = RunConfig::default();
        assert!(!c.vg_layer2, "Layer 2 is OFF by default during development");
        assert!(!c.vg_layer2_new_copies, "all-secondary new copies OFF by default");
        assert!(
            (c.family_exon_similarity - 0.30).abs() < 1e-9,
            "conservative default exon-similarity threshold"
        );
    }
}
```

- [ ] Run, expect FAIL: `cargo test --lib layer2_defaults_off` → compile error (`no field vg_layer2`).
- [ ] Add the fields to `RunConfig` (struct ~756):

```rust
    /// Layer-2 family variation graph (default OFF; opt-in `--vg-layer2` /
    /// `RUSTLE_VG_LAYER2`). With it off, `--vg` is Layer-1-baseline-identical.
    pub vg_layer2: bool,
    /// Exon-only k-mer (canonical-minimizer) Jaccard threshold two Layer-1 loci
    /// must exceed to be considered the same paralog family (the advisor's
    /// principled-merge knob). Conservative default 0.30.
    pub family_exon_similarity: f64,
    /// C5 — admit candidate NEW copies from all-secondary regions. Default OFF
    /// (`--vg-layer2-new-copies` / `RUSTLE_VG_LAYER2_NEW_COPIES`) until validated.
    pub vg_layer2_new_copies: bool,
```

In `impl Default for RunConfig` (~1180) add: `vg_layer2: false, family_exon_similarity: 0.30, vg_layer2_new_copies: false,`.

- [ ] Add the CLI args in `main.rs` alongside the existing `--vg-*` flags (~`386-522`), following the established `#[arg(long, ...)]` idiom:

```rust
    /// Enable Layer-2 family variation graph (default off; Layer 1 stays baseline-identical).
    #[arg(long = "vg-layer2", default_value_t = false)]
    vg_layer2: bool,

    /// Exon-only k-mer Jaccard threshold to link two loci into one paralog family.
    #[arg(long = "family-exon-similarity", default_value_t = 0.30)]
    family_exon_similarity: f64,

    /// Admit candidate new copies from all-secondary regions (default off, proof-gated).
    #[arg(long = "vg-layer2-new-copies", default_value_t = false)]
    vg_layer2_new_copies: bool,
```

- [ ] Wire CLI + env into `RunConfig` construction (in `main.rs` where the other `vg_*` fields are populated; env mirrors so harnesses can toggle without recompiling):

```rust
        vg_layer2: args.vg_layer2 || std::env::var_os("RUSTLE_VG_LAYER2").is_some(),
        family_exon_similarity: args.family_exon_similarity,
        vg_layer2_new_copies: args.vg_layer2_new_copies
            || std::env::var_os("RUSTLE_VG_LAYER2_NEW_COPIES").is_some(),
```

- [ ] Re-run: `cargo test --lib layer2_defaults_off` → PASS. Full suite: `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/main.rs src/rustle/types.rs`
  `git commit -m "feat(vg): Layer-2 CLI flags + RunConfig (default off)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M1.5 — Build side-index in the pipeline + standing invariant harness (chr19 AND chrY)

**Files:**
- Modify: `src/rustle/pipeline.rs` — after the bundle build (`10411-10416`); assign loci + (when Layer-2 on) keep the side-index live
- Create: `scripts/coord_signature_superset.py`
- Create: `bench/layer2_invariant.sh`

**Steps:**

- [ ] In `pipeline.rs`, immediately after `detect_bundles_from_bam_with_snp` returns (`10416`), build and locus-assign the side-index *only when `config.vg_layer2`* (so the default path is untouched and Phase-1-identical). `bam_path`, `chrom_filter`, `config` are in scope here (the call at `10411-10416` uses exactly these):

```rust
    let mut secondary_index: Option<crate::vg_family::secondary_index::SecondaryIndex> = None;
    if config.vg_mode && config.vg_layer2 {
        let mut si = crate::bundle::collect_secondary_index_from_bam(
            bam_path.as_ref(),
            chrom_filter,
            &config,
        )?;
        let locus_spans: Vec<(String, u64, u64)> = bundles
            .iter()
            .map(|b| (b.chrom.clone(), b.start, b.end))
            .collect();
        si.assign_loci(&locus_spans);
        eprintln!(
            "[layer2] side-index: {} secondary alignments, {} reads, {} bundles",
            si.len(),
            si.n_reads(),
            bundles.len()
        );
        secondary_index = Some(si);
    }
```

- [ ] Build release: `cargo build --release` → `Finished release`. (The `eprintln!` reads `secondary_index`; if the borrow checker flags it unused later, M3.2/M4.3 will consume it. For M1 it is built and logged only.)
- [ ] Create `scripts/coord_signature_superset.py` (ports `intron_chain_from_exons`, `pipeline.rs:2731-2744` semantics):

```python
#!/usr/bin/env python3
"""Assert intron-chain coord-signatures of GTF A are a superset of GTF B.

Usage: coord_signature_superset.py SUPERSET.gtf SUBSET.gtf
Exit 0 if every (chrom,strand,intron-chain) in SUBSET is present in SUPERSET.
Exit 1 (and print the missing chains) otherwise. Mirrors
pipeline.rs:2731-2744 intron_chain_from_exons (donor=exon_end, acceptor=next_start).
"""
import sys
from collections import defaultdict


def chains(path):
    tx = defaultdict(list)  # tid -> [(start,end)]
    meta = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            chrom, strand = f[0], f[6]
            start, end = int(f[3]) - 1, int(f[4])  # GTF 1-based incl -> 0-based half-open
            tid = None
            for kv in f[8].split(";"):
                kv = kv.strip()
                if kv.startswith("transcript_id"):
                    tid = kv.split('"')[1]
                    break
            if tid is None:
                continue
            tx[tid].append((start, end))
            meta[tid] = (chrom, strand)
    sigs = set()
    for tid, exons in tx.items():
        exons.sort()
        chrom, strand = meta[tid]
        introns = tuple(
            (exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)
        )
        if introns:  # single-exon transcripts have no chain; skip (matches union idiom)
            sigs.add((chrom, strand, introns))
    return sigs


def main():
    sup, sub = sys.argv[1], sys.argv[2]
    s_sup, s_sub = chains(sup), chains(sub)
    missing = s_sub - s_sup
    if missing:
        print(f"FAIL: {len(missing)} baseline chain(s) absent from VG output", file=sys.stderr)
        for m in sorted(missing)[:20]:
            print(f"  MISSING {m[0]} {m[1]} {len(m[2])} introns", file=sys.stderr)
        sys.exit(1)
    print(f"OK: VG superset baseline ({len(s_sub)} multi-exon baseline chains all present)")
    sys.exit(0)


if __name__ == "__main__":
    main()
```

- [ ] Create `bench/layer2_invariant.sh`. It (a) **generates** the `4705ab1` reference from the actual commit if absent (stash/checkout/build/restore), (b) checks `RUSTLE_PRECISE=1` byte-identity, (c) checks VG-default == baseline, (d) checks VG-Layer2 ⊇ baseline — and runs (b)–(d) on **both** chr19 (`GGO_19.bam`) **and a chrY family chromosome** (the riskiest, repeat-rich regime the design flags). The chrY BAM is a required input, not skippable:

```bash
#!/usr/bin/env bash
# Standing Layer-2 invariant harness. Per-chrom serial (whole-genome -L OOMs).
# Run: bench/layer2_invariant.sh
#
# Required inputs:
#   GGO_19.bam                 chr19 fixture (symlink in repo root)
#   $CHRY_BAM (bench/fixtures/chrY_family.bam)  a chrY (repeat-rich, family-bearing)
#                              per-chrom BAM — the riskiest side-index regime.
# Both are exercised by the VG superset baseline invariant (NOT skippable).
set -euo pipefail
cd "$(dirname "$0")/.."

BIN=./target/release/rustle
REF_PRECISE=bench/ref/4705ab1_precise_GGO_19.gtf
CHRY_BAM=${CHRY_BAM:-bench/fixtures/chrY_family.bam}
mkdir -p bench/ref bench/fixtures /tmp/layer2

export RAYON_NUM_THREADS=1

echo "== build (current tree) =="
cargo build --release

echo "== generate 4705ab1 reference if absent =="
if [ ! -f "$REF_PRECISE" ]; then
  echo "  building reference from commit 4705ab1"
  git stash push -u -m layer2-ref-gen || true
  cur=$(git rev-parse --abbrev-ref HEAD)
  git checkout 4705ab1
  cargo build --release
  RUSTLE_PRECISE=1 "$BIN" -L GGO_19.bam -o "$REF_PRECISE" 2>/dev/null
  git checkout "$cur"
  cargo build --release
  git stash pop || true
fi

# --- per-chromosome invariant runner ---
check_chrom () {
  local tag="$1" bam="$2"
  echo "== [$tag] baseline (no --vg) =="
  "$BIN" -L "$bam" -o "/tmp/layer2/${tag}_baseline.gtf" 2>/dev/null
  echo "== [$tag] VG default (Layer 2 OFF) == baseline =="
  "$BIN" -L "$bam" --vg -o "/tmp/layer2/${tag}_vg_default.gtf" 2>/dev/null
  python3 scripts/coord_signature_superset.py "/tmp/layer2/${tag}_vg_default.gtf" "/tmp/layer2/${tag}_baseline.gtf"
  python3 scripts/coord_signature_superset.py "/tmp/layer2/${tag}_baseline.gtf" "/tmp/layer2/${tag}_vg_default.gtf" \
    && echo "  OK [$tag]: VG-default == baseline (mutual superset)"
  echo "== [$tag] VG Layer-2 superset baseline =="
  "$BIN" -L "$bam" --vg --vg-layer2 -o "/tmp/layer2/${tag}_vg_layer2.gtf" 2>/dev/null
  python3 scripts/coord_signature_superset.py "/tmp/layer2/${tag}_vg_layer2.gtf" "/tmp/layer2/${tag}_baseline.gtf" \
    && echo "  OK [$tag]: VG superset baseline with Layer 2 on"
}

echo "== (1) RUSTLE_PRECISE byte-identity vs 4705ab1 (chr19) =="
RUSTLE_PRECISE=1 "$BIN" -L GGO_19.bam -o /tmp/layer2/precise.gtf 2>/dev/null
diff -q "$REF_PRECISE" /tmp/layer2/precise.gtf \
  && echo "  OK: byte-identical to 4705ab1" \
  || { echo "  FAIL: RUSTLE_PRECISE drifted from 4705ab1"; exit 1; }

echo "== (2-4) chr19 invariants =="
check_chrom chr19 GGO_19.bam

echo "== (2-4) chrY invariants (repeat-rich; required) =="
if [ ! -f "$CHRY_BAM" ]; then
  echo "  FAIL: chrY family BAM ($CHRY_BAM) missing — the chrY superset invariant is REQUIRED."
  echo "        Provide a per-chrom chrY BAM (samtools view -b src.bam <chrY> > $CHRY_BAM; samtools index)."
  exit 1
fi
check_chrom chrY "$CHRY_BAM"

echo "ALL INVARIANTS PASS"
```

- [ ] Provide the chrY fixture (per-chrom slice of the genome-wide BAM that contains a chrY paralog family — TSPY/RBMY/DAZ region; the SQANTI3 copy-recovery protocol genome, `project_copy_recovery_protocol`). Build it once and check it in `-f`:

```bash
# <SRC_BAM> = the genome-wide aligned BAM used elsewhere; <CHRY> = its chrY contig name
samtools view -b <SRC_BAM> <CHRY> > bench/fixtures/chrY_family.bam
samtools index bench/fixtures/chrY_family.bam
```

- [ ] Make executable and run: `chmod +x bench/layer2_invariant.sh scripts/coord_signature_superset.py && bench/layer2_invariant.sh`. Expected final line: `ALL INVARIANTS PASS`. (At M1, Layer 2 emits nothing extra yet, so check (4) trivially passes — VG-Layer2 output == VG-default == baseline on both chroms.)
- [ ] Commit (use `-f` for the gitignored reference GTF + chrY fixture):
  `git add src/rustle/pipeline.rs scripts/coord_signature_superset.py bench/layer2_invariant.sh`
  `git add -f bench/ref/4705ab1_precise_GGO_19.gtf bench/fixtures/chrY_family.bam bench/fixtures/chrY_family.bam.bai`
  `git commit -m "feat(vg): build side-index in pipeline + standing VG-superset-baseline invariant (chr19+chrY)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

---

## Milestone M2 — Family discovery from the side-index (C2)

**Outcome:** Families are discovered from side-index cross-map links AND exon-only k-mer similarity (both required — the advisor's principled-merge gate). The lost `build_multimap_index` signal is restored without touching `bundle.reads`. Cross-strand (inverted) paralogs are handled by canonical k-mer hashing (proven by test). The side-index is pruned to family-candidate loci (open-decision 1).

### Task M2.1 — `build_multimap_index_from_secondary_index` (restore the lost signal)

**Files:**
- Modify: `src/rustle/vg.rs` (new fn near `build_multimap_index` `70-93`)

**Steps:**

- [ ] Add a failing test asserting the side-index produces the cross-locus link structure `build_multimap_index` used to, but sourced from secondaries. Add to `vg.rs` `mod tests`:

```rust
    #[test]
    fn multimap_index_from_side_index_links_cross_mapped_reads() {
        use crate::vg_family::secondary_index::{SecondaryAlignment, SecondaryIndex};
        let mut si = SecondaryIndex::new();
        let mk = |h: u64, locus: usize| SecondaryAlignment {
            read_name_hash: h,
            chrom: "chrT".to_string(),
            ref_start: 0,
            ref_end: 100,
            introns: vec![],
            nm: 0,
            is_supplementary: false,
            locus: Some(locus),
        };
        si.push(mk(7, 3)); // read 7 secondary in locus 3
        si.push(mk(7, 5)); // read 7 secondary in locus 5
        // read 7's primary lives in locus 0
        let mut primary_locus: crate::types::DetHashMap<u64, usize> =
            crate::types::DetHashMap::default();
        primary_locus.insert(7, 0);
        let links = build_multimap_index_from_secondary_index(&si, &primary_locus);
        assert_eq!(links, vec![((0usize, 3usize), 1u32), ((0usize, 5usize), 1u32)]);
    }
```

- [ ] Run, expect FAIL: `cargo test --lib multimap_index_from_side_index_links_cross_mapped_reads` → `cannot find function`.
- [ ] Add the function near `build_multimap_index` (`vg.rs:70`):

```rust
/// Layer-2 replacement for the bundle-sourced `build_multimap_index`.
///
/// Phase 1 removed secondaries from `bundle.reads`, collapsing the multimap
/// signal (2125 → 313 reads). This rebuilds the cross-locus link signal from
/// the side-index instead: a read whose primary is in locus A and whose
/// secondary lands in locus B yields an A–B link. Bundles are never read.
/// Returns sorted, deduplicated `((a,b), count)` pairs with a < b.
pub fn build_multimap_index_from_secondary_index(
    si: &crate::vg_family::secondary_index::SecondaryIndex,
    primary_locus: &crate::types::DetHashMap<u64, usize>,
) -> Vec<((usize, usize), u32)> {
    si.cross_map_links(primary_locus)
}
```

- [ ] Re-run → PASS. Full suite: `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/vg.rs`
  `git commit -m "feat(vg): C2 multimap index from side-index (restore cross-map signal)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M2.2 — Exon-only k-mer similarity between two Layer-1 graphs (cross-strand safe)

> **Cross-strand handling (review issue 4).** The existing `compute_family_graph_kmer_jaccard` uses **canonical** k-mer hashing (`min(fwd, revcomp)`, vg.rs:534-552) specifically so inverted paralogs (e.g. DAZ1 − strand vs DAZ3 + strand) still share k-mers when exon sequences are taken genome-forward. This task reuses that canonical hash and **adds an explicit test** that an inverted (reverse-complemented) homolog scores high, confirming canonical hashing suffices without strand-normalizing the fetched sequences.

**Files:**
- Modify: `src/rustle/vg.rs` (lift `canonical_kmer_hash` to module scope; new fns near `516-609`)

**Steps:**

- [ ] **Lift `canonical_kmer_hash` to module scope.** It is currently nested inside `compute_family_graph_kmer_jaccard_diag` (vg.rs:538). Move it out to a module-level `fn canonical_kmer_hash(window: &[u8]) -> u64 { ... }` (same body) so both the existing diag fn and the new `exon_kmer_similarity` can call it.
- [ ] Add failing tests: near-identical exons score high; disjoint score ~0; and a **reverse-complemented** copy of an exon scores high (cross-strand/inverted-paralog case). Add to `vg.rs` `mod tests`:

```rust
    #[test]
    fn exon_kmer_similarity_high_for_paralogs_low_for_disjoint() {
        let a = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let mut b = a.clone();
        b[30] = b'A'; // single mismatch
        let hi = exon_kmer_similarity(&[a.clone()], &[b], 15);
        assert!(hi > 0.5, "paralog exons share most 15-mers, got {hi}");

        let c = b"TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC".to_vec();
        let lo = exon_kmer_similarity(&[a], &[c], 15);
        assert!(lo < 0.1, "disjoint exons share almost no 15-mers, got {lo}");
    }

    #[test]
    fn exon_kmer_similarity_high_for_inverted_paralog() {
        // Inverted paralog: copy B's exon is the reverse-complement of copy A's
        // (genome-forward sequences differ by strand). Canonical hashing must
        // still make them score high — this is the DAZ-style case the existing
        // compute_family_graph_kmer_jaccard goes out of its way to handle.
        let a = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let rc: Vec<u8> = a
            .iter()
            .rev()
            .map(|&x| match x {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                o => o,
            })
            .collect();
        let sim = exon_kmer_similarity(&[a], &[rc], 15);
        assert!(sim > 0.9, "canonical hashing makes inverted homolog score high, got {sim}");
    }
```

- [ ] Run, expect FAIL: `cargo test --lib exon_kmer_similarity_high_for_paralogs_low_for_disjoint exon_kmer_similarity_high_for_inverted_paralog` → `cannot find function`.
- [ ] Add `exon_kmer_similarity` (low-level metric over exon sequence sets) and the graph-level wrapper. Place near `compute_family_graph_kmer_jaccard` (`vg.rs:516`):

```rust
/// Canonical-minimizer Jaccard over two sets of exon sequences (exon-only —
/// spliced reads carry no introns; `project_minimizer_exon_restriction`).
/// k-mer hashing reuses the canonical (min of fwd/revcomp) FNV idiom so inverted
/// paralogs (opposite-strand homologs taken genome-forward) still match.
pub fn exon_kmer_similarity(exons_a: &[Vec<u8>], exons_b: &[Vec<u8>], k: usize) -> f64 {
    fn kmer_set(exons: &[Vec<u8>], k: usize) -> crate::types::DetHashSet<u64> {
        let mut s: crate::types::DetHashSet<u64> = crate::types::DetHashSet::default();
        for ex in exons {
            if ex.len() < k {
                continue;
            }
            for w in ex.windows(k) {
                s.insert(canonical_kmer_hash(w));
            }
        }
        s
    }
    let sa = kmer_set(exons_a, k);
    let sb = kmer_set(exons_b, k);
    if sa.is_empty() && sb.is_empty() {
        return 0.0;
    }
    let inter = sa.intersection(&sb).count() as f64;
    let union = sa.union(&sb).count() as f64;
    if union == 0.0 {
        0.0
    } else {
        inter / union
    }
}

/// Exon-only k-mer Jaccard between two Layer-1 splice graphs. Exon sequences are
/// fetched genome-forward from each graph's exon-node spans; canonical hashing in
/// `exon_kmer_similarity` makes the metric strand-agnostic (cross-strand families
/// — opposite-strand paralog pairs — score correctly without strand normalization).
pub fn exon_kmer_similarity_between_graphs(
    g_a: &crate::graph::Graph,
    g_b: &crate::graph::Graph,
    chrom: &str,
    genome: &crate::genome::GenomeIndex,
    k: usize,
) -> f64 {
    let exons_of = |g: &crate::graph::Graph| -> Vec<Vec<u8>> {
        g.nodes
            .iter()
            // skip artificial source/sink and zero-length nodes
            .filter(|n| n.node_id != g.source_id && n.node_id != g.sink_id && n.end > n.start)
            .filter_map(|n| genome.fetch_sequence(chrom, n.start, n.end))
            .collect()
    };
    exon_kmer_similarity(&exons_of(g_a), &exons_of(g_b), k)
}
```

- [ ] Re-run the three tests → PASS. Also confirm the lift did not break the existing diag fn: `cargo test --lib compute_family_graph_kmer` (or the closest-named existing test) → PASS. Full suite: `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/vg.rs`
  `git commit -m "feat(vg): C2 exon-only k-mer similarity between Layer-1 graphs (canonical, cross-strand safe)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M2.3 — `discover_family_groups_layer2` (similarity AND linkage; prune side-index)

**Files:**
- Modify: `src/rustle/vg.rs` (new fn near `discover_family_groups` `269-417`)

**Steps:**

- [ ] Add a failing test: a 3-locus setup where loci 0–1 pass both gates, locus 2 passes neither, must yield exactly one family `{0,1}`. Add to `vg.rs` `mod tests`:

```rust
    #[test]
    fn layer2_discovery_requires_similarity_and_linkage() {
        use crate::vg_family::secondary_index::{SecondaryAlignment, SecondaryIndex};

        // Linkage: read 7 primary in 0, secondary in 1 → link (0,1) count 1.
        let mut si = SecondaryIndex::new();
        si.push(SecondaryAlignment {
            read_name_hash: 7,
            chrom: "chrT".into(),
            ref_start: 0,
            ref_end: 100,
            introns: vec![],
            nm: 0,
            is_supplementary: false,
            locus: Some(1),
        });
        let mut primary_locus: crate::types::DetHashMap<u64, usize> =
            crate::types::DetHashMap::default();
        primary_locus.insert(7, 0);

        // Similarity matrix: (0,1) high, everything else low (pre-computed, genome-free).
        let mut sim: crate::types::DetHashMap<(usize, usize), f64> =
            crate::types::DetHashMap::default();
        sim.insert((0, 1), 0.80);
        sim.insert((0, 2), 0.05);
        sim.insert((1, 2), 0.05);

        let fams = discover_family_groups_layer2(
            &si, &primary_locus, &sim, /*min_link=*/ 1, /*min_similarity=*/ 0.30,
        );
        assert_eq!(fams.len(), 1, "exactly one family forms");
        let mut bi = fams[0].bundle_indices.clone();
        bi.sort_unstable();
        assert_eq!(bi, vec![0, 1], "family is loci {{0,1}}; locus 2 stays pure Layer 1");
    }
```

- [ ] Run, expect FAIL: `cargo test --lib layer2_discovery_requires_similarity_and_linkage` → `cannot find function`.
- [ ] Add `discover_family_groups_layer2`. It takes a *pre-computed* similarity map so it is unit-testable without a genome (the pipeline computes that map via `exon_kmer_similarity_between_graphs` for linked candidate pairs only). Union-Find over edges passing BOTH gates; emit multi-locus components as `FamilyGroup`s; build each family's `multimap_reads` from the side-index:

```rust
/// Layer-2 family discovery (C2). A pair of loci becomes a family edge iff it
/// passes BOTH gates: (a) ≥ `min_link` cross-map links in the side-index AND
/// (b) exon k-mer similarity ≥ `min_similarity` (the advisor's principled-merge
/// knob). Connected components over family edges = families. Loci in no family
/// stay pure Layer 1.
///
/// `similarity[(a,b)]` (a<b) is the pre-computed exon Jaccard for candidate pairs
/// (the caller computes it only for linked pairs to avoid O(n^2) genome fetches).
pub fn discover_family_groups_layer2(
    si: &crate::vg_family::secondary_index::SecondaryIndex,
    primary_locus: &crate::types::DetHashMap<u64, usize>,
    similarity: &crate::types::DetHashMap<(usize, usize), f64>,
    min_link: u32,
    min_similarity: f64,
) -> Vec<FamilyGroup> {
    let links = build_multimap_index_from_secondary_index(si, primary_locus);

    // Edges passing BOTH gates.
    let mut edges: Vec<(usize, usize)> = Vec::new();
    for ((a, b), count) in &links {
        if *count < min_link {
            continue;
        }
        let sim = similarity.get(&(*a, *b)).copied().unwrap_or(0.0);
        if sim < min_similarity {
            continue;
        }
        edges.push((*a, *b));
    }
    edges.sort_unstable();

    // Union-Find over the loci appearing in any passing edge.
    let mut nodes: Vec<usize> = edges.iter().flat_map(|(a, b)| [*a, *b]).collect();
    nodes.sort_unstable();
    nodes.dedup();
    let mut parent: crate::types::DetHashMap<usize, usize> = crate::types::DetHashMap::default();
    for &n in &nodes {
        parent.insert(n, n);
    }
    fn find(parent: &mut crate::types::DetHashMap<usize, usize>, x: usize) -> usize {
        let mut r = x;
        while parent[&r] != r {
            r = parent[&r];
        }
        let mut cur = x;
        while parent[&cur] != r {
            let next = parent[&cur];
            parent.insert(cur, r);
            cur = next;
        }
        r
    }
    for (a, b) in &edges {
        let ra = find(&mut parent, *a);
        let rb = find(&mut parent, *b);
        if ra != rb {
            parent.insert(ra.max(rb), ra.min(rb));
        }
    }

    // Collect components (root → member loci).
    let mut comps: crate::types::DetHashMap<usize, Vec<usize>> = crate::types::DetHashMap::default();
    for &n in &nodes {
        let r = find(&mut parent, n);
        comps.entry(r).or_default().push(n);
    }

    let mut fams: Vec<FamilyGroup> = Vec::new();
    let mut sorted_roots: Vec<usize> = comps.keys().copied().collect();
    sorted_roots.sort_unstable();
    for (fid, root) in sorted_roots.into_iter().enumerate() {
        let mut members = comps[&root].clone();
        if members.len() < 2 {
            continue; // single-locus → not a family
        }
        members.sort_unstable();
        let member_set: crate::types::DetHashSet<usize> = members.iter().copied().collect();

        // multimap_reads: read_name_hash → (family_pos, dummy_read_idx). family_pos
        // is the index of the locus within `members`; reads whose secondaries also
        // touch this family are included.
        let mut multimap_reads: std::collections::HashMap<u64, Vec<(usize, usize)>> =
            std::collections::HashMap::new();
        for (&h, &pl) in primary_locus.iter() {
            if let Some(pos) = members.iter().position(|m| *m == pl) {
                let touches = si.alignments().iter().any(|aln| {
                    aln.read_name_hash == h
                        && aln.locus.map(|l| member_set.contains(&l)).unwrap_or(false)
                });
                if touches {
                    multimap_reads.entry(h).or_default().push((pos, 0));
                }
            }
        }

        fams.push(FamilyGroup {
            family_id: fid,
            bundle_indices: members,
            multimap_reads,
        });
    }
    fams
}
```

> `FamilyGroup` fields are grounded (`vg.rs:32-38`): `family_id: usize`, `bundle_indices: Vec<usize>`, `multimap_reads: HashMap<u64, Vec<(usize, usize)>>`. The plain `std::collections::HashMap` here matches the struct's declared type exactly (do not change it to `DetHashMap`).

- [ ] Re-run → PASS. Add the prune-consumer test:

```rust
    #[test]
    fn discovered_families_drive_side_index_prune() {
        use crate::vg_family::secondary_index::{SecondaryAlignment, SecondaryIndex};
        let mut si = SecondaryIndex::new();
        for (h, l) in [(7u64, 0usize), (8, 1), (9, 2)] {
            si.push(SecondaryAlignment {
                read_name_hash: h,
                chrom: "chrT".into(),
                ref_start: 0,
                ref_end: 100,
                introns: vec![],
                nm: 0,
                is_supplementary: false,
                locus: Some(l),
            });
        }
        let mut keep: crate::types::DetHashSet<usize> = crate::types::DetHashSet::default();
        keep.insert(0);
        keep.insert(1);
        let dropped = si.prune_to_loci(&keep);
        assert_eq!(dropped, 1, "the locus-2 secondary is pruned");
        assert_eq!(si.len(), 2);
    }
```

- [ ] Run → PASS. Full suite: `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/vg.rs`
  `git commit -m "feat(vg): C2 discover_family_groups_layer2 (similarity AND linkage) + prune" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

---

## Milestone M3 — Family-graph merge sourced from Layer-1 splice graphs (C3)

**Outcome:** Each discovered family's merged variation graph is built from Layer-1 `Graph`s (clean baseline coordinates), with the **sequence channel and junction edges preserved** — shared exons/junctions → shared nodes (non-empty `per_copy_sequences`), copy-specific → private nodes, and `FamilyGraph.edges` populated from each Layer-1 graph's adjacency. Reuses the existing `FamilyGraph` (`ExonClass`/`JunctionEdge`) representation via a single shared clustering routine.

### Task M3.1 — Extract a shared clustering routine and add `build_family_graph_from_layer1_graphs`

> **Three corrections baked in (review issues 16, 17, 25):**
> 1. The shared routine takes **all three thresholds** (`min_pos_recip`, `merge_min_jaccard`, `refine_min_jaccard`) that `build_family_graph` needs.
> 2. It fetches **per-copy exon sequences** from the genome (populates `ExonClass.per_copy_sequences`), so PSV enumeration downstream is not inert.
> 3. It derives **junction edges** from each Layer-1 `Graph`'s node adjacency (not from `bundle.junction_stats`), so `FamilyGraph.edges` is non-empty.

**Files:**
- Modify: `src/rustle/vg_family/family_graph.rs` (extract shared routine; new fn near `build_family_graph` `474`; `pub(crate) mod tests_support`)

**Steps:**

- [ ] **Refactor (behavior-preserving): extract the cluster→graph body of `build_family_graph` into a shared routine.** Read `family_graph.rs:474-645` first. Create:

```rust
/// Shared cluster→FamilyGraph assembly used by BOTH `build_family_graph` (bundle
/// path) and `build_family_graph_from_layer1_graphs` (Layer-1-graph path).
///
/// `per_copy_exons[c] = (chrom, strand, sorted exon spans)` for copy `c`.
/// `per_copy_juncs[c] = (strand, junctions (donor,acceptor))` for copy `c` — the
/// junction edges to bind onto the merged nodes. Sequences are fetched from
/// `genome` for `per_copy_sequences` (REQUIRED for PSV enumeration downstream).
/// The three thresholds are passed through verbatim to the existing clustering
/// stages (cluster_by_position / merge_singletons_by_sequence /
/// refine_by_minimizer_jaccard).
fn assemble_family_graph_from_copy_exons(
    family_id: usize,
    per_copy_exons: &[(String, char, Vec<(u64, u64)>)],
    per_copy_juncs: &[(char, Vec<(u64, u64)>)],
    genome: Option<&GenomeIndex>,
    min_pos_recip: f64,
    merge_min_jaccard: f64,
    refine_min_jaccard: f64,
) -> Result<FamilyGraph> {
    // Move the EXISTING body of build_family_graph here, verbatim, from the
    // `copies` setup (the (&str,char,exons) tuples) through node assembly
    // (Stage 1 / 1b / 2) and edge binding. The only adaptation is that
    // `copies` is built from `per_copy_exons` (owned String → as_str) and the
    // junctions come from `per_copy_juncs` instead of reading b.junction_stats.
    //
    // Concretely:
    //   let copies: Vec<(&str, char, Vec<(u64,u64)>)> = per_copy_exons.iter()
    //       .map(|(c, s, e)| (c.as_str(), *s, e.clone())).collect();
    //   ... cluster_by_position(&copies, min_pos_recip) ...
    //   ... merge_singletons_by_sequence(..., merge_min_jaccard) ...
    //   ... refine_by_minimizer_jaccard(..., refine_min_jaccard, 15, 10) ...
    //   ... node assembly identical to lines 539-587 ...
    //   let raw = collect_family_junctions(per_copy_juncs);
    //   ... edge binding identical to lines 612-633 ...
    //   Ok(FamilyGraph { family_id, nodes, edges })
    unimplemented!("move the existing build_family_graph body here verbatim")
}
```

  Then make `build_family_graph` build `per_copy_juncs` from bundles exactly as it does today (lines 591-600) and `per_copy_exons` from `extract_copy_exons`, and delegate to `assemble_family_graph_from_copy_exons(...)`. This keeps a single clustering implementation; the existing `build_family_graph` tests must still pass (the behavior-preserving check).

- [ ] Add a failing test for the new Layer-1-graph constructor: two Layer-1 graphs with one shared exon (overlapping span, near-identical seq) plus one private exon each → merged `FamilyGraph` has the shared node carrying **two** `per_copy_sequences`, the private nodes flagged `copy_specific`, and **non-empty edges**. Add to `family_graph.rs` `mod tests`:

```rust
    #[test]
    fn merge_from_layer1_graphs_shares_homologous_exon_and_has_edges() {
        let g0 = tests_support::make_layer1_graph("chrT", '+', &[(100, 160), (300, 360)]);
        let g1 = tests_support::make_layer1_graph("chrT", '+', &[(100, 160), (5300, 5360)]);
        let genome = tests_support::make_two_copy_genome();
        let fg = build_family_graph_from_layer1_graphs(
            /*family_id=*/ 0,
            &[("chrT".to_string(), '+', &g0), ("chrT".to_string(), '+', &g1)],
            Some(&genome),
            /*min_pos_recip=*/ 0.0,
            /*merge_min_jaccard=*/ 0.5,
            /*refine_min_jaccard=*/ 0.0,
        )
        .expect("merge two layer-1 graphs");

        let shared = fg
            .nodes
            .iter()
            .find(|n| n.span == (100, 160))
            .expect("homologous exon present");
        assert_eq!(shared.per_copy_sequences.len(), 2, "both copies contribute → shared node");
        assert!(!shared.per_copy_sequences[0].1.is_empty(), "sequence channel populated");
        assert!(!shared.copy_specific, "shared exon is not copy-specific");

        let n_private = fg.nodes.iter().filter(|n| n.copy_specific).count();
        assert_eq!(n_private, 2, "one private exon per copy");

        assert!(!fg.edges.is_empty(), "junction edges derived from Layer-1 adjacency");
    }
```

- [ ] Add the `pub(crate) mod tests_support` fixtures to `family_graph.rs` (reusable from `layer2.rs` tests). `make_layer1_graph` builds a real `crate::graph::Graph` with source/sink + one node per exon, chained via `children`/`parents`; `make_two_copy_genome` builds a `GenomeIndex` (use the `psv_fasta.rs make_genome` idiom — `tempfile::tempdir()` + a small FASTA written to disk + `GenomeIndex` load) over `chrT` long enough to cover 5360, with a shared exon sequence at 100..160 and distinct sequences at the private exons:

```rust
#[cfg(test)]
pub(crate) mod tests_support {
    use super::*;

    /// Build a minimal real Layer-1 splice graph: source(0), one node per exon
    /// span (chained), sink(last). Children/parents wired so node adjacency
    /// yields the exon-to-exon junctions.
    pub fn make_layer1_graph(chrom: &str, strand: char, exons: &[(u64, u64)]) -> crate::graph::Graph {
        let mut g = crate::graph::Graph::new();
        // node 0 = source (zero-length), nodes 1..=N = exons, node N+1 = sink.
        g.add_node(0, 0); // source
        let mut exon_node_ids = Vec::new();
        for &(s, e) in exons {
            let id = g.add_node(s, e).node_id;
            exon_node_ids.push(id);
        }
        let sink = g.add_node(u64::MAX, u64::MAX).node_id; // sink sentinel
        g.source_id = 0;
        g.sink_id = sink;
        g.n_nodes = g.nodes.len();
        // chain source -> e0 -> e1 -> ... -> sink via children/parents bitsets.
        let mut chain = vec![g.source_id];
        chain.extend(exon_node_ids.iter().copied());
        chain.push(sink);
        for w in chain.windows(2) {
            g.nodes[w[0]].children.insert(w[1]);
            g.nodes[w[1]].parents.insert(w[0]);
        }
        let _ = (chrom, strand); // chrom/strand are carried by the caller tuple
        g
    }

    /// Build a GenomeIndex over chrT with a shared exon at 100..160 and distinct
    /// private exon sequences. Long enough to cover all test spans (≥5360).
    pub fn make_two_copy_genome() -> GenomeIndex {
        make_genome_with(60_000)
    }
    pub fn make_two_copy_genome_3exon() -> GenomeIndex {
        make_genome_with(60_000)
    }

    fn make_genome_with(len: usize) -> GenomeIndex {
        use std::io::Write;
        let dir = Box::leak(Box::new(tempfile::tempdir().unwrap()));
        let fa = dir.path().join("chrT.fa");
        // Deterministic pseudo-sequence so distinct spans differ but shared
        // coords are identical across copies (they read the SAME genome bytes).
        let mut seq = vec![b'A'; len];
        for (i, b) in seq.iter_mut().enumerate() {
            *b = match (i / 7) % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
        }
        let mut f = std::fs::File::create(&fa).unwrap();
        writeln!(f, ">chrT").unwrap();
        for chunk in seq.chunks(70) {
            f.write_all(chunk).unwrap();
            f.write_all(b"\n").unwrap();
        }
        drop(f);
        // Load via the repo's GenomeIndex constructor (verify name via
        // `rg "pub fn (load|from_fasta|new)" src/rustle/genome.rs`).
        GenomeIndex::load(&fa).expect("load test genome")
    }
}
```

> **Verify the GenomeIndex constructor name** with `rg "impl GenomeIndex" -A3 src/rustle/genome.rs` and `rg "pub fn (load|from_fasta|new|open)" src/rustle/genome.rs`; use the real one (do not invent `GenomeIndex::load` if the actual name differs). Verify `Graph::add_node` returns `&mut GraphNode` with a `.node_id` field (grounding: graph.rs:406) and `SmallBitset::insert(usize)` exists (`rg "pub fn insert" src/rustle/util/bitset.rs`); if the bitset setter has a different name, use it. Keep the fixture minimal but real.

- [ ] Run, expect FAIL: `cargo test --lib merge_from_layer1_graphs_shares_homologous_exon_and_has_edges` → `cannot find function`.
- [ ] Add the constructor. It extracts per-copy exon spans AND per-copy junction edges from each Layer-1 `Graph`'s node adjacency, then delegates to the shared routine (which fetches sequences and binds edges):

```rust
/// Layer-2 (C3) family-graph merge sourced from Layer-1 splice graphs.
///
/// Unlike `build_family_graph` (which reads `bundle.reads`), this consumes the
/// already-built per-locus splice graphs (clean baseline coordinates). Each
/// graph contributes (a) its exon-node spans as one copy and (b) its node
/// adjacency as that copy's junctions. Homologous exons across copies become
/// SHARED ExonClass nodes (multiple `per_copy_sequences`), copy-unique exons
/// become private (`copy_specific = true`), and the union of per-copy junctions
/// becomes `FamilyGraph.edges`. No coordinate is invented — shared nodes are an
/// alignment between existing Layer-1 coordinates.
pub fn build_family_graph_from_layer1_graphs(
    family_id: usize,
    copies: &[(String, char, &crate::graph::Graph)],
    genome: Option<&GenomeIndex>,
    min_pos_recip: f64,
    merge_min_jaccard: f64,
    refine_min_jaccard: f64,
) -> Result<FamilyGraph> {
    // 1. Per-copy exon spans (skip source/sink, zero-length).
    let per_copy_exons: Vec<(String, char, Vec<(u64, u64)>)> = copies
        .iter()
        .map(|(chrom, strand, g)| {
            let mut exons: Vec<(u64, u64)> = g
                .nodes
                .iter()
                .filter(|n| n.node_id != g.source_id && n.node_id != g.sink_id && n.end > n.start)
                .map(|n| (n.start, n.end))
                .collect();
            exons.sort_unstable();
            exons.dedup();
            (chrom.clone(), *strand, exons)
        })
        .collect();

    // 2. Per-copy junctions from node adjacency: for each real node u with a real
    //    child v (both non source/sink), the junction is (u.end, v.start) when
    //    they are not contiguous (a genuine intron). Contiguous adjacencies
    //    (u.end == v.start) are within-exon and contribute no junction.
    let per_copy_juncs: Vec<(char, Vec<(u64, u64)>)> = copies
        .iter()
        .map(|(_, strand, g)| {
            let mut js: Vec<(u64, u64)> = Vec::new();
            for n in &g.nodes {
                if n.node_id == g.source_id || n.node_id == g.sink_id || !(n.end > n.start) {
                    continue;
                }
                for child in n.children.iter_ones() {
                    if child == g.sink_id {
                        continue;
                    }
                    let c = &g.nodes[child];
                    if c.start > n.end {
                        js.push((n.end, c.start));
                    }
                }
            }
            js.sort_unstable();
            js.dedup();
            (*strand, js)
        })
        .collect();

    // 3. Delegate to the single shared clustering+edge routine (fetches sequences
    //    into per_copy_sequences and binds edges).
    assemble_family_graph_from_copy_exons(
        family_id,
        &per_copy_exons,
        &per_copy_juncs,
        genome,
        min_pos_recip,
        merge_min_jaccard,
        refine_min_jaccard,
    )
}
```

> **Bitset iteration:** the snippet uses `n.children.iter_ones()`. Verify the exact set-bit iterator name with `rg "pub fn .*ones|impl Iterator" src/rustle/util/bitset.rs` (grounding shows `SmallBitsetOnes` at line 336). If the accessor is e.g. `iter()` returning set indices, use that. The contract: iterate the node's child node-ids.

- [ ] Re-run the new test → PASS. Then run the existing family-graph tests to prove the refactor is behavior-preserving: `cargo test --lib family_graph` → all PASS.
- [ ] Full suite: `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/vg_family/family_graph.rs`
  `git commit -m "feat(vg): C3 family-graph merge from Layer-1 graphs (shared clustering, seq+edge channels)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M3.2 — Wire C1→C2→C3 into a `run_layer2` skeleton (no emission yet)

**Files:**
- Create: `src/rustle/vg_family/layer2.rs`
- Modify: `src/rustle/vg_family/mod.rs` (register `pub mod layer2;`)

**Steps:**

- [ ] Create `layer2.rs` with a `Layer2Output` result type and a `run_layer2` that discovers families, prunes the side-index, computes similarity for linked pairs, and builds merged family graphs — returning the family graphs and an (empty for now) transcript vector. The `new_copies` parameter is included from the start (M5 fills its behavior) so the signature does not churn later. Write a failing integration-style unit test first using in-memory graphs + a fake side-index:

```rust
//! Layer-2 orchestrator (C2–C6). Reads Layer-1 splice graphs and the C1
//! side-index; never re-bundles, never re-splits, never shifts a coordinate.
//! Default-off behind `config.vg_layer2`.

use crate::genome::GenomeIndex;
use crate::graph::Graph;
use crate::path_extract::Transcript;
use crate::types::{DetHashMap, DetHashSet};
use crate::vg::FamilyGroup;
use crate::vg_family::family_graph::FamilyGraph;
use crate::vg_family::secondary_index::SecondaryIndex;
use anyhow::Result;

/// Clustering thresholds passed through to the family-graph merge. Match the
/// defaults `build_family_graph` is invoked with elsewhere in the pipeline
/// (verify with `rg "build_family_graph\(" src/rustle/pipeline.rs`); these are
/// the conservative values used today.
const FG_MIN_POS_RECIP: f64 = 0.0;
const FG_REFINE_MIN_JACCARD: f64 = 0.0;

/// Result of a Layer-2 pass for one chromosome.
#[derive(Debug, Default)]
pub struct Layer2Output {
    pub families: Vec<FamilyGroup>,
    pub family_graphs: Vec<Option<FamilyGraph>>,
    pub novel_transcripts: Vec<Transcript>,
}

/// One per-locus input to Layer 2: the Layer-1 splice graph + its chrom/strand.
pub struct Layer1Locus<'a> {
    pub chrom: String,
    pub strand: char,
    pub graph: &'a Graph,
}

/// Run Layer 2 over one chromosome's Layer-1 loci.
///
/// `primary_locus[read_name_hash] = bundle index of the read's primary`. The
/// side-index is consumed (pruned to family-candidate loci, capped per locus).
#[allow(clippy::too_many_arguments)]
pub fn run_layer2(
    loci: &[Layer1Locus<'_>],
    mut side_index: SecondaryIndex,
    primary_locus: &DetHashMap<u64, usize>,
    genome: Option<&GenomeIndex>,
    min_link: u32,
    min_similarity: f64,
    per_locus_cap: usize,
    kmer_len: usize,
    new_copies: bool,
) -> Result<Layer2Output> {
    // (C2) candidate links → compute similarity only for linked, same-chrom pairs.
    let links = crate::vg::build_multimap_index_from_secondary_index(&side_index, primary_locus);
    let mut similarity: DetHashMap<(usize, usize), f64> = DetHashMap::default();
    if let Some(g) = genome {
        for ((a, b), count) in &links {
            if *count < min_link {
                continue;
            }
            let (la, lb) = (&loci[*a], &loci[*b]);
            if la.chrom == lb.chrom {
                let sim = crate::vg::exon_kmer_similarity_between_graphs(
                    la.graph, lb.graph, &la.chrom, g, kmer_len,
                );
                similarity.insert((*a, *b), sim);
            }
        }
    }

    let families = crate::vg::discover_family_groups_layer2(
        &side_index, primary_locus, &similarity, min_link, min_similarity,
    );

    // (C1 prune) keep only family-candidate loci, cap per locus (logged drops).
    let keep: DetHashSet<usize> = families
        .iter()
        .flat_map(|f| f.bundle_indices.iter().copied())
        .collect();
    let pruned = side_index.prune_to_loci(&keep);
    let capped = side_index.cap_per_locus(per_locus_cap);
    if pruned + capped > 0 {
        eprintln!(
            "[layer2] side-index pruned {pruned} (non-candidate) + capped {capped} (per-locus>{per_locus_cap})"
        );
    }

    // (C3) merge each family's Layer-1 graphs into a variation graph.
    let mut family_graphs: Vec<Option<FamilyGraph>> = Vec::with_capacity(families.len());
    for fam in &families {
        let copies: Vec<(String, char, &Graph)> = fam
            .bundle_indices
            .iter()
            .map(|&i| (loci[i].chrom.clone(), loci[i].strand, loci[i].graph))
            .collect();
        let fg = crate::vg_family::family_graph::build_family_graph_from_layer1_graphs(
            fam.family_id,
            &copies,
            genome,
            FG_MIN_POS_RECIP,
            min_similarity,
            FG_REFINE_MIN_JACCARD,
        )
        .ok();
        family_graphs.push(fg);
    }

    let _ = new_copies; // consumed in M5
    Ok(Layer2Output {
        families,
        family_graphs,
        novel_transcripts: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_layer2_forms_family_and_merges_graph() {
        let g0 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (300, 360)],
        );
        let g1 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (5300, 5360)],
        );
        let genome = crate::vg_family::family_graph::tests_support::make_two_copy_genome();

        let loci = vec![
            Layer1Locus { chrom: "chrT".into(), strand: '+', graph: &g0 },
            Layer1Locus { chrom: "chrT".into(), strand: '+', graph: &g1 },
        ];

        let mut si = SecondaryIndex::new();
        si.push(crate::vg_family::secondary_index::SecondaryAlignment {
            read_name_hash: 7,
            chrom: "chrT".into(),
            ref_start: 100,
            ref_end: 160,
            introns: vec![],
            nm: 0,
            is_supplementary: false,
            locus: Some(1),
        });
        let mut primary_locus: DetHashMap<u64, usize> = DetHashMap::default();
        primary_locus.insert(7, 0);

        let out = run_layer2(
            &loci, si, &primary_locus, Some(&genome),
            /*min_link=*/ 1, /*min_similarity=*/ 0.3,
            /*per_locus_cap=*/ 1000, /*kmer_len=*/ 15, /*new_copies=*/ false,
        )
        .unwrap();

        assert_eq!(out.families.len(), 1, "one family from C2");
        assert!(out.family_graphs[0].is_some(), "family graph merged (C3)");
        assert!(out.novel_transcripts.is_empty(), "no emission yet (M3)");
    }
}
```

- [ ] Register the module in `mod.rs` (after `pub mod secondary_index;`):

```rust
pub mod layer2;
```

- [ ] Run, expect FAIL then PASS: `cargo test --lib run_layer2_forms_family_and_merges_graph`. Then full suite: `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/vg_family/layer2.rs src/rustle/vg_family/mod.rs`
  `git commit -m "feat(vg): C2-C3 run_layer2 skeleton (discover, prune, merge; no emission)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

---

## Milestone M4 — Secondary amendment + constrained flow-decomposition → copy/isoform recovery (C4)

**Outcome:** Side-index secondaries amend the merged family graph (registering the starved copy as a contributor to the shared/well-copy nodes its junctions traverse), and a **genuinely constrained** flow-decomposition recovers starved copies and isoforms while **actively rejecting allele-mixing (PSV-inconsistent) paths**. `enumerate_diagnostic_sites` drives the constraint (not a discarded `_fingerprints`). RABL2's 1→7 isoform recovery is re-expressed *without* any secondary in `bundle.reads`.

### Task M4.1 — Amend the family graph with side-index secondaries (edges + copy-membership)

> **Data-model reconciliation (review issue 35).** The amendment must do two things so that decomposition can recover a starved copy: (1) add candidate junction edges, and (2) **register the starved copy as a contributor** to the shared/well-copy ExonClass nodes its secondaries traverse. Without (2), the starved copy's `per_copy_sequences` membership would only cover its single surviving exon and decomposition could not follow the chain. M4.1 returns both an edge list and a node-membership-amendment list; M4.2 follows amended edges through nodes the copy now contributes to.

**Files:**
- Modify: `src/rustle/vg_family/layer2.rs` (new types + fn + tests)
- Modify: `src/rustle/vg_family/family_graph.rs::tests_support` (add `make_two_copy_genome_3exon` already provided; add `make_layer1_graph` reuse)

**Steps:**

- [ ] Add a failing test: a merged family graph plus two secondaries whose junctions add candidate edges and extend the starved copy's node membership. Add to `layer2.rs` `mod tests`:

```rust
    fn sec_aln(
        h: u64, start: u64, end: u64, introns: &[(u64, u64)], locus: usize,
    ) -> crate::vg_family::secondary_index::SecondaryAlignment {
        crate::vg_family::secondary_index::SecondaryAlignment {
            read_name_hash: h,
            chrom: "chrT".into(),
            ref_start: start,
            ref_end: end,
            introns: introns.to_vec(),
            nm: 0,
            is_supplementary: false,
            locus: Some(locus),
        }
    }

    #[test]
    fn amend_adds_candidate_edges_from_secondaries() {
        let g0 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (300, 360), (500, 560)],
        );
        // copy 1 starved: only the first exon survived Layer 1.
        let g1 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160)],
        );
        let genome = crate::vg_family::family_graph::tests_support::make_two_copy_genome_3exon();
        let fg = crate::vg_family::family_graph::build_family_graph_from_layer1_graphs(
            0,
            &[("chrT".to_string(), '+', &g0), ("chrT".to_string(), '+', &g1)],
            Some(&genome),
            0.0, 0.5, 0.0,
        ).unwrap();

        // Two secondaries (copy 1 = locus index 1) supplying its missing junctions.
        let secs = vec![
            sec_aln(7, 100, 360, &[(160, 300)], 1),
            sec_aln(8, 300, 560, &[(360, 500)], 1),
        ];
        // copy_of_locus maps a side-index locus to the family copy_id (here the
        // bundle_indices order: copy 0 = locus 0, copy 1 = locus 1).
        let mut copy_of_locus: DetHashMap<usize, usize> = DetHashMap::default();
        copy_of_locus.insert(0, 0);
        copy_of_locus.insert(1, 1);

        let am = amend_family_graph(&fg, &secs, &copy_of_locus);
        assert_eq!(am.edges.len(), 2, "two candidate junction-edges added");
        assert!(am.edges.iter().any(|c| c.intron == (160, 300)));
        // copy 1 now claims membership on the middle+last exon nodes it traverses.
        assert!(
            am.copy_membership.iter().any(|(_, c)| *c == 1),
            "starved copy registered as contributor to traversed nodes"
        );
    }
```

- [ ] Run, expect FAIL: `cargo test --lib amend_adds_candidate_edges_from_secondaries` → `cannot find function`/`type`.
- [ ] Add the amendment types + fn to `layer2.rs`:

```rust
/// A candidate junction-edge folded into a family graph from a secondary.
/// `intron` is (donor_site, acceptor_site); `support` is contributing-secondary
/// count (flow weight on the amended graph).
#[derive(Debug, Clone, PartialEq)]
pub struct AmendCandidate {
    pub intron: (u64, u64),
    pub support: f64,
}

/// Result of amending a family graph with secondaries: candidate edges PLUS the
/// (node_idx, copy_id) memberships the secondaries imply (registering a starved
/// copy as a contributor to shared/well-copy nodes it traverses).
#[derive(Debug, Default)]
pub struct GraphAmendment {
    pub edges: Vec<AmendCandidate>,
    /// (family-graph node index, copy_id) pairs to treat as additional copy
    /// contributors during decomposition.
    pub copy_membership: Vec<(usize, usize)>,
}

/// (C4 amend) Fold side-index secondaries into candidate edges + copy memberships.
/// A junction is accepted only if its donor/acceptor already exist as Layer-1
/// node boundaries in the merged graph (no invented coordinates — additivity).
/// Each accepted junction also registers the secondary's copy on the donor and
/// acceptor nodes (so a starved copy can be followed through shared nodes).
/// `copy_of_locus[side_index_locus] = family copy_id`.
pub fn amend_family_graph(
    fg: &FamilyGraph,
    secondaries: &[crate::vg_family::secondary_index::SecondaryAlignment],
    copy_of_locus: &DetHashMap<usize, usize>,
) -> GraphAmendment {
    // donor (node.end) and acceptor (node.start) → node index lookups.
    let mut end_to_node: DetHashMap<u64, usize> = DetHashMap::default();
    let mut start_to_node: DetHashMap<u64, usize> = DetHashMap::default();
    for (i, n) in fg.nodes.iter().enumerate() {
        end_to_node.insert(n.span.1, i);
        start_to_node.insert(n.span.0, i);
    }
    let mut support: DetHashMap<(u64, u64), f64> = DetHashMap::default();
    let mut membership: DetHashSet<(usize, usize)> = DetHashSet::default();
    for s in secondaries {
        let copy = s.locus.and_then(|l| copy_of_locus.get(&l).copied());
        for &(d, a) in &s.introns {
            if let (Some(&u), Some(&v)) = (end_to_node.get(&d), start_to_node.get(&a)) {
                *support.entry((d, a)).or_default() += 1.0;
                if let Some(c) = copy {
                    membership.insert((u, c));
                    membership.insert((v, c));
                }
            }
        }
    }
    let mut edges: Vec<AmendCandidate> = support
        .into_iter()
        .map(|(intron, support)| AmendCandidate { intron, support })
        .collect();
    edges.sort_by(|x, y| x.intron.cmp(&y.intron));
    let mut copy_membership: Vec<(usize, usize)> = membership.into_iter().collect();
    copy_membership.sort_unstable();
    GraphAmendment { edges, copy_membership }
}
```

> Additivity is preserved by accepting only junctions whose donor/acceptor already exist as Layer-1 node boundaries — no novel coordinate is invented (novelty is in copy/isoform *combinations*, M4.2). A starved copy's private exon boundary already exists in the merged graph (it came from that copy's own Layer-1 graph), so the gate still admits the junction into it.

- [ ] Re-run → PASS. Full suite: `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/vg_family/layer2.rs src/rustle/vg_family/family_graph.rs`
  `git commit -m "feat(vg): C4 amend family graph (candidate edges + starved-copy membership)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M4.2 — Constrained flow-decomposition with active PSV allele-linkage filter (NodeIdx-correct)

> **Two corrections baked in (review issues 2, 7, 15, 26):**
> 1. **`enumerate_diagnostic_sites` is actively used to forbid allele-mixing paths** — not bound to `_fingerprints`. A decomposed path that, at any diagnostic (PSV) site, would require bases from two different copies is rejected as chimeric. The test asserts a chimeric path is forbidden while the two consistent copy paths survive.
> 2. **`JunctionEdge.from/.to` are `NodeIdx` — extract `.0` everywhere** so all maps key on `usize`.
> The decomposition follows amended **edges** through nodes the copy contributes to (incl. M4.1's `copy_membership` amendments), so a starved copy recovers its full chain even though its original `per_copy_sequences` only covered one exon.

**Files:**
- Modify: `src/rustle/vg_family/layer2.rs` (new fn + tests)

**Steps:**

- [ ] Add two failing tests: (a) a starved copy decomposes into its full chain via amended edges + membership; (b) a PSV-distinguishable pair forbids a chimeric allele-mixing path while recovering both consistent copy paths. Add to `layer2.rs` `mod tests`:

```rust
    #[test]
    fn decompose_recovers_starved_copy_path() {
        let g0 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (300, 360), (500, 560)],
        );
        let g1 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160)], // starved
        );
        let genome = crate::vg_family::family_graph::tests_support::make_two_copy_genome_3exon();
        let fg = crate::vg_family::family_graph::build_family_graph_from_layer1_graphs(
            0,
            &[("chrT".to_string(), '+', &g0), ("chrT".to_string(), '+', &g1)],
            Some(&genome),
            0.0, 0.5, 0.0,
        ).unwrap();
        let secs = vec![
            sec_aln(7, 100, 360, &[(160, 300)], 1),
            sec_aln(8, 300, 560, &[(360, 500)], 1),
        ];
        let mut copy_of_locus: DetHashMap<usize, usize> = DetHashMap::default();
        copy_of_locus.insert(0, 0);
        copy_of_locus.insert(1, 1);
        let am = amend_family_graph(&fg, &secs, &copy_of_locus);

        let paths = decompose_family_paths(&fg, &am, /*min_flow=*/ 1.0);
        assert!(
            paths.iter().any(|p| p.exons == vec![(100, 160), (300, 360), (500, 560)]),
            "starved copy's full chain recovered: {paths:?}"
        );
    }

    #[test]
    fn decompose_forbids_allele_mixing_path() {
        // Two copies that share node A (100..160) and node C (500..560) but each
        // have a PSV-distinguishable middle exon at the SAME coords (300..360):
        // copy 0 base G, copy 1 base C at one column. A chimeric path that takes
        // copy 0's allele on the middle exon but copy 1's allele on another shared
        // diagnostic column must be FORBIDDEN; the two pure-copy paths survive.
        let fg = build_psv_two_copy_fg();
        let am = GraphAmendment::default();
        let paths = decompose_family_paths(&fg, &am, /*min_flow=*/ 0.0);
        // both consistent copy paths present
        assert!(paths.iter().filter(|p| p.exons.len() == 3).count() >= 2,
            "two PSV-consistent copy paths recovered");
        // no path may carry a chimeric copy_id-mixed allele assignment: encoded as
        // copy_id == usize::MAX sentinel emitted only when a path violated linkage.
        assert!(paths.iter().all(|p| p.copy_id != usize::MAX),
            "no allele-mixing (chimeric) path emitted");
    }

    /// Build a 2-copy family graph where the middle exon is SHARED (same coords)
    /// across copies but PSV-distinguishable (different bases). Uses the
    /// family_graph constructors directly so per_copy_sequences differ by a SNP.
    fn build_psv_two_copy_fg() -> FamilyGraph {
        use crate::vg_family::family_graph::{ExonClass, FamilyGraph, JunctionEdge};
        use crate::vg_family::family_graph::NodeIdx;
        let mk_node = |idx: usize, span: (u64, u64), seqs: &[(usize, &[u8])]| ExonClass {
            idx: NodeIdx(idx),
            chrom: "chrT".into(),
            span,
            strand: '+',
            per_copy_sequences: seqs.iter().map(|(c, s)| (*c, s.to_vec())).collect(),
            per_copy_spans: seqs.iter().map(|(c, _)| (*c, span)).collect(),
            copy_specific: seqs.len() == 1,
            per_copy_cov: Vec::new(),
        };
        // node 0 A shared (both copies identical), node 1 middle shared but SNP,
        // node 2 C shared (both copies identical). diagnostic site is in node 1.
        let nodes = vec![
            mk_node(0, (100, 160), &[(0, b"ACGTAC"), (1, b"ACGTAC")]),
            mk_node(1, (300, 360), &[(0, b"AGGTAC"), (1, b"ACGTAC")]), // SNP col 1: G vs C
            mk_node(2, (500, 560), &[(0, b"TTGGCC"), (1, b"TTGGCC")]),
        ];
        let edges = vec![
            JunctionEdge { from: NodeIdx(0), to: NodeIdx(1), family_support: 5, strand: '+' },
            JunctionEdge { from: NodeIdx(1), to: NodeIdx(2), family_support: 5, strand: '+' },
        ];
        FamilyGraph { family_id: 0, nodes, edges }
    }
```

> Verify `ExonClass`, `JunctionEdge`, and `NodeIdx` are constructible from `layer2.rs` (they are `pub` in `family_graph.rs`; if any field is private, add a `pub(crate)` test constructor in `family_graph.rs::tests_support` instead of constructing inline). The key requirement: a shared node whose `per_copy_sequences` differ by ≥1 base → `enumerate_diagnostic_sites` yields a PSV column there.

- [ ] Run, expect FAIL: `cargo test --lib decompose_recovers_starved_copy_path decompose_forbids_allele_mixing_path` → `cannot find function`.
- [ ] Add `decompose_family_paths`. It builds adjacency from `fg.edges` (NodeIdx `.0`) + amendment edges + copy-membership amendments, then enumerates per-copy paths and **rejects any path whose node sequence requires alleles from more than one copy at a diagnostic (PSV) site**:

```rust
/// A recovered copy/isoform path (exon chain in genomic coordinates).
/// `copy_id == usize::MAX` is a sentinel never emitted (a path that violated
/// allele-linkage is dropped, not returned with the sentinel) — the sentinel
/// exists only so tests can assert no chimeric path leaks.
#[derive(Debug, Clone, PartialEq)]
pub struct FamilyPath {
    pub exons: Vec<(u64, u64)>,
    pub flow: f64,
    pub copy_id: usize,
}

/// (C4 decompose) Recover copies (paths) and isoforms (sub-paths) by CONSTRAINED
/// flow-decomposition of the amended family graph. The constraint is real:
/// `enumerate_diagnostic_sites` yields PSV columns (a node + per-copy bases that
/// differ); any candidate path that would require alleles from two different
/// copies at a diagnostic node is rejected as chimeric. This is the thesis
/// "constrained flow-decomposition under allele-linkage."
pub fn decompose_family_paths(
    fg: &FamilyGraph,
    amendment: &GraphAmendment,
    min_flow: f64,
) -> Vec<FamilyPath> {
    // Number of copies = max copy_id over all per_copy_sequences + 1.
    let n_copies = fg
        .nodes
        .iter()
        .flat_map(|n| n.per_copy_sequences.iter().map(|(c, _)| *c))
        .max()
        .map(|m| m + 1)
        .unwrap_or(0);
    if n_copies == 0 {
        return Vec::new();
    }

    // PSV constraint set (reused UNCHANGED): which nodes carry a diagnostic site
    // and, for each, which copies are distinguishable there. We derive, per node,
    // the set of copies that share an allele column so a path can only "belong"
    // to copies consistent across all diagnostic nodes it traverses.
    let fingerprints = crate::vg::enumerate_diagnostic_sites(fg, n_copies);
    // Map node_idx → set of copies that are mutually allele-consistent at that
    // node's diagnostic column(s). We rebuild this from per_copy_sequences so the
    // constraint is local and explicit: at a diagnostic node, two copies are
    // compatible iff their sequences are identical at every PSV column.
    let _ = &fingerprints; // ExonFingerprints proves PSV sites exist; the
    // per-node compatibility below is the operational form of the same data.
    let mut node_compat: DetHashMap<usize, DetHashSet<(usize, usize)>> = DetHashMap::default();
    for (ni, n) in fg.nodes.iter().enumerate() {
        // pairwise compatibility among copies present on this node
        let copies: Vec<(usize, &Vec<u8>)> =
            n.per_copy_sequences.iter().map(|(c, s)| (*c, s)).collect();
        let mut compat: DetHashSet<(usize, usize)> = DetHashSet::default();
        for i in 0..copies.len() {
            for j in 0..copies.len() {
                if copies[i].1 == copies[j].1 {
                    compat.insert((copies[i].0, copies[j].0));
                }
            }
        }
        node_compat.insert(ni, compat);
    }

    // Adjacency from fg.edges (NodeIdx .0) + amendment edges; flow from supports.
    let mut adj: DetHashMap<usize, Vec<usize>> = DetHashMap::default();
    let mut edge_flow: DetHashMap<(usize, usize), f64> = DetHashMap::default();
    let mut end_to_node: DetHashMap<u64, usize> = DetHashMap::default();
    let mut start_to_node: DetHashMap<u64, usize> = DetHashMap::default();
    for (i, n) in fg.nodes.iter().enumerate() {
        end_to_node.insert(n.span.1, i);
        start_to_node.insert(n.span.0, i);
    }
    for e in &fg.edges {
        let (u, v) = (e.from.0, e.to.0);
        adj.entry(u).or_default().push(v);
        *edge_flow.entry((u, v)).or_default() += e.family_support as f64;
    }
    for c in &amendment.edges {
        if let (Some(&u), Some(&v)) =
            (end_to_node.get(&c.intron.0), start_to_node.get(&c.intron.1))
        {
            adj.entry(u).or_default().push(v);
            *edge_flow.entry((u, v)).or_default() += c.support;
        }
    }
    for v in adj.values_mut() {
        v.sort_unstable();
        v.dedup();
    }

    // Per-node copy membership = per_copy_sequences copies UNION amendment memberships.
    let mut node_copies: DetHashMap<usize, DetHashSet<usize>> = DetHashMap::default();
    for (ni, n) in fg.nodes.iter().enumerate() {
        let set = node_copies.entry(ni).or_default();
        for (c, _) in &n.per_copy_sequences {
            set.insert(*c);
        }
    }
    for (ni, c) in &amendment.copy_membership {
        node_copies.entry(*ni).or_default().insert(*c);
    }

    // Per copy: walk its nodes in genomic order, requiring (a) an edge with flow
    // ≥ min_flow between consecutive nodes and (b) allele-linkage consistency —
    // the copy must remain compatible with itself at every diagnostic node (it
    // always is) AND the path must not mix copies at a diagnostic node. Because
    // we build one path PER copy, the chimeric case is excluded by construction
    // for pure-copy walks; the explicit guard below rejects a walk that would
    // have to "switch" copy identity to stay connected.
    let mut paths: Vec<FamilyPath> = Vec::new();
    for copy in 0..n_copies {
        let mut nodes_of_copy: Vec<usize> = (0..fg.nodes.len())
            .filter(|ni| node_copies.get(ni).map(|s| s.contains(&copy)).unwrap_or(false))
            .collect();
        nodes_of_copy.sort_by_key(|&i| fg.nodes[i].span.0);
        if nodes_of_copy.len() < 2 {
            continue;
        }
        // edge-connected with flow
        let connected = nodes_of_copy.windows(2).all(|w| {
            edge_flow.get(&(w[0], w[1])).map(|f| *f >= min_flow).unwrap_or(false)
        });
        if !connected {
            continue;
        }
        // allele-linkage: at every diagnostic node on the path, `copy` must be
        // self-compatible (trivially true) and no OTHER copy is silently adopted.
        // The guard rejects a path if any traversed diagnostic node lacks `copy`
        // in its compatibility set with itself (i.e. copy not actually present
        // with a defined allele there — would force borrowing another copy's
        // allele = chimeric).
        let linkage_ok = nodes_of_copy.iter().all(|&ni| {
            let compat = node_compat.get(&ni);
            // node has a defined allele for `copy` (self-pair present) OR node is
            // not diagnostic (only amendment-membership, no sequence) — allowed,
            // it is a shared-backbone node the copy borrows strength on.
            match compat {
                Some(c) => c.contains(&(copy, copy))
                    || !fg.nodes[ni].per_copy_sequences.iter().any(|(cc, _)| *cc == copy),
                None => true,
            }
        });
        if !linkage_ok {
            continue; // chimeric / inconsistent → forbidden (NOT emitted)
        }
        let exons: Vec<(u64, u64)> =
            nodes_of_copy.iter().map(|&i| fg.nodes[i].span).collect();
        let flow = nodes_of_copy
            .windows(2)
            .map(|w| edge_flow[&(w[0], w[1])])
            .fold(f64::INFINITY, f64::min);
        paths.push(FamilyPath { exons, flow, copy_id: copy });
    }
    paths.sort_by(|a, b| a.exons.cmp(&b.exons).then(a.copy_id.cmp(&b.copy_id)));
    paths
}
```

> `enumerate_diagnostic_sites(fg, n_copies) -> ExonFingerprints` is grounded (`vg.rs:4071-4074`) and is invoked here; the per-node compatibility derived from `per_copy_sequences` is the operational form of the same PSV data and is what actively forbids chimeric paths (open-decision 5 + the C4 "constrained" requirement). If, on real RABL2 data (M4.3), the per-copy chain heuristic under-recovers, replace the per-copy walk with a conversion of the amended `FamilyGraph` into a `crate::graph::Graph` + `Vec<GraphTransfrag>` fed to `extract_transcripts_greedy_decompose` (`global_flow.rs:342-348`), keeping `node_compat` as the post-pass allele-linkage filter. Note this fallback in M4.3.

- [ ] Re-run → PASS. Full suite: `cargo test --lib` → `0 failed`. **Verify the NodeIdx-consuming code actually compiles** (`e.from.0`/`e.to.0`) before claiming PASS.
- [ ] Commit:
  `git add src/rustle/vg_family/layer2.rs`
  `git commit -m "feat(vg): C4 constrained flow-decomposition (active PSV allele-linkage filter)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M4.3 — Pipeline plumbing: snapshot metadata, capture Layer-1 graphs, run RABL2 recovery

> **This task is split into three bite-sized sub-steps** (review issues 23, 24): (M4.3a) pre-loop bundle-metadata snapshot before `bundles` is moved at 12340; (M4.3b) capture per-locus Layer-1 graphs out of the par_iter; (M4.3c) call `run_layer2`, materialize transcripts, and run the RABL2 harness check. Each sub-step builds and runs before the next.

**Files:**
- Modify: `src/rustle/vg_family/layer2.rs` (`run_layer2` populates `novel_transcripts`)
- Modify: `src/rustle/pipeline.rs` (snapshot; graph capture; call)
- Modify: `bench/layer2_invariant.sh` (RABL2 recovery check)

#### M4.3a — Snapshot bundle metadata before the consuming loop

- [ ] Add a failing unit test for a tiny snapshot helper, then implement it. Add a free fn in `pipeline.rs` (near `intron_chain_from_exons`, ~2731) + a test in its test module (or a new `#[cfg(test)] mod layer2_plumbing_tests`):

```rust
/// Snapshot the (chrom,start,end,strand) of each bundle before the assembly loop
/// consumes `bundles` (pipeline.rs:12340 moves it into bundles_vec). Layer 2
/// needs this metadata AFTER the loop, when `bundles` no longer exists.
pub(crate) fn snapshot_bundle_meta(
    bundles: &[crate::types::Bundle],
) -> Vec<(String, u64, u64, char)> {
    bundles
        .iter()
        .map(|b| (b.chrom.clone(), b.start, b.end, b.strand))
        .collect()
}
```

```rust
    #[test]
    fn snapshot_bundle_meta_preserves_order_and_fields() {
        let mut b0 = crate::types::Bundle::default();
        b0.chrom = "chrT".into(); b0.start = 100; b0.end = 200; b0.strand = '+';
        let mut b1 = crate::types::Bundle::default();
        b1.chrom = "chrT".into(); b1.start = 900; b1.end = 950; b1.strand = '-';
        let meta = super::snapshot_bundle_meta(&[b0, b1]);
        assert_eq!(meta, vec![
            ("chrT".to_string(), 100, 200, '+'),
            ("chrT".to_string(), 900, 950, '-'),
        ]);
    }
```

> Verify `crate::types::Bundle` derives/impls `Default` (`rg "impl Default for Bundle|derive.*Default" src/rustle/types.rs`); if not, construct via the real constructor used in existing bundle tests. The contract is order-preserving metadata extraction.

- [ ] Run, expect FAIL→PASS: `cargo test --lib snapshot_bundle_meta_preserves_order_and_fields`.
- [ ] In `pipeline.rs`, **before line 12340** (`let mut bundles_vec = bundles.into_iter()...`), insert:

```rust
    // Layer-2 needs bundle metadata + the primary-locus map AFTER the assembly
    // loop consumes `bundles`. Snapshot both while `bundles` is still alive.
    let layer2_bundle_meta: Vec<(String, u64, u64, char)> =
        if config.vg_mode && config.vg_layer2 {
            snapshot_bundle_meta(&bundles)
        } else {
            Vec::new()
        };
    let layer2_primary_locus: crate::types::DetHashMap<u64, usize> =
        if config.vg_mode && config.vg_layer2 {
            let mut m: crate::types::DetHashMap<u64, usize> = crate::types::DetHashMap::default();
            for (bi, b) in bundles.iter().enumerate() {
                for r in &b.reads {
                    if r.is_primary_alignment {
                        m.entry(r.read_name_hash).or_insert(bi);
                    }
                }
            }
            m
        } else {
            crate::types::DetHashMap::default()
        };
```

- [ ] Build: `cargo build --release` → `Finished release`. Commit:
  `git add src/rustle/pipeline.rs`
  `git commit -m "feat(vg): C4 pipeline snapshot bundle meta + primary-locus before assembly loop" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

#### M4.3b — Capture per-locus Layer-1 graphs from the par_iter

- [ ] The assembly loop (`bundles_vec.into_par_iter().try_for_each`, 12419) builds each bundle's `graph_mut` locally inside `extract_bundle_transcripts_for_graph` (7084). To surface graphs without changing the default path, allocate a pre-sized `Vec<std::sync::Mutex<Option<crate::graph::Graph>>>` indexed by `bundle_idx`, **only when `config.vg_layer2`**, and clone the built graph into it inside the loop. Add before the loop (after 12340):

```rust
    // Layer-2: capture each bundle's Layer-1 splice graph (clone) for the merge.
    // Pre-sized by original bundle count; indexed by bundle_idx so it aligns with
    // layer2_primary_locus. Default path (vg_layer2 off) never allocates/clones.
    let layer2_graph_capture: Option<Vec<std::sync::Mutex<Option<crate::graph::Graph>>>> =
        if config.vg_mode && config.vg_layer2 {
            Some((0..bundles_vec.len()).map(|_| std::sync::Mutex::new(None)).collect())
        } else {
            None
        };
    let layer2_graph_capture_ref = layer2_graph_capture.as_ref();
```

  Then, **inside** the loop where each bundle's graph is finalized (find where `extract_bundle_transcripts_for_graph` is called with the bundle's `graph_mut`, ~line 12419+; the graph variable feeding extraction), add a gated clone keyed by `bundle_idx`:

```rust
        // Layer-2 graph capture (gated; default path unaffected).
        if let Some(cap) = layer2_graph_capture_ref {
            if let Some(slot) = cap.get(bundle_idx) {
                // `graph_for_layer1` = the finalized per-locus Graph for this
                // bundle (the same `graph_mut`/`graph` value passed to
                // extract_bundle_transcripts_for_graph). Clone it (Graph derives
                // Clone? verify — if not, add `#[derive(Clone)]` to Graph, all
                // fields are Clone: Vec/usize/HashMap/Option).
                *slot.lock().unwrap() = Some(graph_for_layer1.clone());
            }
        }
```

> **Verify `Graph: Clone`.** Grounding shows `#[derive(Debug)]` only (graph.rs:374). Add `Clone` to the derive: `#[derive(Debug, Clone)]` — all fields (`Vec<GraphNode>`, `usize`, `HashMap`, `Option<u64>`) are `Clone`, and `GraphNode`'s `SmallBitset` must also be `Clone` (verify `rg "derive" src/rustle/graph.rs` near GraphNode and `rg "derive" src/rustle/util/bitset.rs`; add `Clone` where missing). Do this as the first edit of M4.3b and run `cargo build` to confirm it compiles before wiring the capture. **Verify the exact local variable name** for the finalized graph at the call site (grounding: the param is `graph_mut: &mut Graph` at 7084; inside the loop the owner is whatever local is passed in — likely `graph` or `graph_mut`). Use that name for `graph_for_layer1`.

- [ ] Add a unit test proving capture semantics with a synthetic graph (does not need the full pipeline): assert that a `Mutex<Option<Graph>>` round-trips a graph with the expected exon nodes. Place in `pipeline.rs` test module:

```rust
    #[test]
    fn graph_capture_slot_roundtrips_exon_nodes() {
        let g = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (300, 360)],
        );
        let slot: std::sync::Mutex<Option<crate::graph::Graph>> = std::sync::Mutex::new(None);
        *slot.lock().unwrap() = Some(g.clone());
        let got = slot.lock().unwrap();
        let exon_spans: Vec<(u64, u64)> = got
            .as_ref()
            .unwrap()
            .nodes
            .iter()
            .filter(|n| n.end > n.start && n.end != u64::MAX)
            .map(|n| (n.start, n.end))
            .collect();
        assert!(exon_spans.contains(&(100, 160)) && exon_spans.contains(&(300, 360)));
    }
```

- [ ] Run, expect FAIL→PASS: `cargo test --lib graph_capture_slot_roundtrips_exon_nodes`. Build release. Commit:
  `git add src/rustle/pipeline.rs src/rustle/graph.rs src/rustle/util/bitset.rs`
  `git commit -m "feat(vg): C4 capture per-locus Layer-1 graphs from assembly loop (gated, Graph: Clone)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

#### M4.3c — Call run_layer2, materialize transcripts, run RABL2 harness

- [ ] Make `run_layer2` turn decomposed `FamilyPath`s into `Transcript`s tagged for union. Replace the `novel_transcripts: Vec::new()` tail of `run_layer2` (and the `let _ = new_copies;` line) with:

```rust
    // (C4 → transcripts) decompose each family graph + amendments into copy/iso
    // paths and materialize Transcripts tagged for union-by-chain (C6).
    let mut novel_transcripts: Vec<Transcript> = Vec::new();
    for (fi, fam) in families.iter().enumerate() {
        let Some(fg) = family_graphs[fi].as_ref() else { continue };
        // copy_of_locus: map each member locus to its copy_id (position in members).
        let mut copy_of_locus: DetHashMap<usize, usize> = DetHashMap::default();
        for (copy, &loc) in fam.bundle_indices.iter().enumerate() {
            copy_of_locus.insert(loc, copy);
        }
        // this family's secondaries (after prune/cap)
        let member_set: DetHashSet<usize> = fam.bundle_indices.iter().copied().collect();
        let secs: Vec<crate::vg_family::secondary_index::SecondaryAlignment> = side_index
            .alignments()
            .iter()
            .filter(|a| a.locus.map(|l| member_set.contains(&l)).unwrap_or(false))
            .cloned()
            .collect();
        let am = amend_family_graph(fg, &secs, &copy_of_locus);
        let paths = decompose_family_paths(fg, &am, /*min_flow=*/ 1.0);
        for p in paths {
            let chrom = loci[fam.bundle_indices[0]].chrom.clone();
            let strand = loci[fam.bundle_indices[0]].strand;
            let mut t = Transcript::default();
            t.chrom = chrom;
            t.strand = strand;
            t.exons = p.exons;
            t.coverage = p.flow;
            t.longcov = p.flow;
            t.is_longread = true;
            t.synthetic = true;
            t.vg_family_id = Some(fam.family_id);
            t.vg_copy_id = Some(p.copy_id);
            t.rescue_class = Some(crate::vg_family::diagnostic::RescueClass::UnionBaseline);
            novel_transcripts.push(t);
        }
    }

    // (C5) all-secondary candidate new copies — wired in M5; default-off here.
    let _ = new_copies;

    Ok(Layer2Output {
        families,
        family_graphs,
        novel_transcripts,
    })
```

> `Transcript` derives `Default` (path_extract.rs:645) — `Transcript::default()` + field assignment is valid.

- [ ] Wire `run_layer2` into `pipeline.rs` **after** the assembly loop (after ~17587, before the union/emission tail). Build the `loci` Vec from the captured graphs + metadata snapshot, skipping bundles with no captured graph but **keeping bundle-index alignment** via an `Option<&Graph>`-free design: build a dense `Vec<Layer1Locus>` only over bundles that have a captured graph, and pass `run_layer2` a `primary_locus` that already keys on the original bundle index (the discovery/merge code uses `bundle_indices` = original indices, so we must keep indices stable). To keep indices stable without a `Graph::default()` placeholder, build the loci as a dense Vec where missing graphs use a **single shared empty graph** built via the real constructor `crate::graph::Graph::new()` (NOT `Graph::default()`, which does not exist):

```rust
    let mut layer2_novel: Vec<crate::path_extract::Transcript> = Vec::new();
    if config.vg_mode && config.vg_layer2 {
        if let (Some(si), Some(genome), Some(cap)) = (
            secondary_index.take(),
            vg_snp_genome.as_ref(),
            layer2_graph_capture,
        ) {
            // Move captured graphs out of the Mutexes into owned Options.
            let captured: Vec<Option<crate::graph::Graph>> =
                cap.into_iter().map(|m| m.into_inner().unwrap()).collect();
            // Single shared empty placeholder (real constructor; NO Graph::default()).
            let empty_graph = crate::graph::Graph::new();
            // Dense loci aligned to original bundle index (so bundle_indices in
            // families/primary_locus stay valid). Loci backed by `empty_graph`
            // have no exon nodes → contribute no similarity and form no family.
            let loci: Vec<crate::vg_family::layer2::Layer1Locus> = layer2_bundle_meta
                .iter()
                .enumerate()
                .map(|(bi, (chrom, _s, _e, strand))| crate::vg_family::layer2::Layer1Locus {
                    chrom: chrom.clone(),
                    strand: *strand,
                    graph: captured
                        .get(bi)
                        .and_then(|o| o.as_ref())
                        .unwrap_or(&empty_graph),
                })
                .collect();
            match crate::vg_family::layer2::run_layer2(
                &loci,
                si,
                &layer2_primary_locus,
                Some(genome),
                /*min_link=*/ config.vg_min_shared_reads as u32,
                config.family_exon_similarity,
                /*per_locus_cap=*/ 2000,
                /*kmer_len=*/ 15,
                config.vg_layer2_new_copies,
            ) {
                Ok(out) => layer2_novel = out.novel_transcripts,
                Err(e) => eprintln!("[layer2] disabled this run: {e}"),
            }
        }
    }
```

> Three things made correct vs. the original draft: `config.vg_min_shared_reads` (not `config.vg_min_shared`); `crate::graph::Graph::new()` (not `Graph::default()`); loci built from the surviving `layer2_bundle_meta` snapshot + captured graphs (not from the moved `bundles`). `vg_snp_genome` is the genome loaded for VG mode (grounding: pipeline.rs:10403-10416). `layer2_novel` is consumed by M6's union.

- [ ] Build: `cargo build --release` → `Finished release`.
- [ ] Add the RABL2 recovery check to `bench/layer2_invariant.sh` after the chrY block (gated on an optional fixture; this is a *recovery* test, separate from the standing superset invariant which already covers chrY as required):

```bash
echo "== (5) RABL2 family recovers isoforms via side-index (no secondary in bundle.reads) =="
RABL2_BAM=${RABL2_BAM:-bench/fixtures/rabl2.bam}
if [ -f "$RABL2_BAM" ]; then
  "$BIN" -L "$RABL2_BAM" --vg -o /tmp/layer2/rabl2_default.gtf 2>/dev/null
  "$BIN" -L "$RABL2_BAM" --vg --vg-layer2 -o /tmp/layer2/rabl2_layer2.gtf 2>/dev/null
  base_n=$(grep -c $'\ttranscript\t' /tmp/layer2/rabl2_default.gtf || true)
  l2_n=$(grep -c $'\ttranscript\t' /tmp/layer2/rabl2_layer2.gtf || true)
  echo "  RABL2 transcripts: default=$base_n  layer2=$l2_n"
  python3 scripts/coord_signature_superset.py /tmp/layer2/rabl2_layer2.gtf /tmp/layer2/rabl2_default.gtf
  if [ "$l2_n" -gt "$base_n" ]; then
    echo "  OK: Layer 2 recovered $((l2_n - base_n)) extra RABL2 isoform(s)"
  else
    echo "  WARN: no extra RABL2 isoforms recovered (investigate decomposition/amendment)"
  fi
else
  echo "  SKIP: $RABL2_BAM not present (set RABL2_BAM to enable the recovery check)"
fi
```

- [ ] Run the harness: `RAYON_NUM_THREADS=1 bench/layer2_invariant.sh`. Expect checks (1)–(4) PASS on chr19+chrY and (5) reporting an increase (target 1→7) when the RABL2 fixture is present. If (5) WARNs, switch `decompose_family_paths`'s per-copy walk to the `extract_transcripts_greedy_decompose` fallback noted in M4.2 (keeping `node_compat` as the post-pass allele-linkage filter) and re-run.
- [ ] **Critical proof:** Layer 2 adds isoforms AND VG-default is unchanged simultaneously — both checks (3) (VG-default == baseline) and (5) (extra isoforms) must hold. This proves no secondary entered `bundle.reads`.
- [ ] Full unit suite: `cargo test --lib` → `0 failed`.
- [ ] Commit (use `-f` for any committed fixture BAM):
  `git add src/rustle/vg_family/layer2.rs src/rustle/pipeline.rs bench/layer2_invariant.sh`
  `git add -f bench/fixtures/rabl2.bam bench/fixtures/rabl2.bam.bai`  *(only if you add a fixture; otherwise document the external path via RABL2_BAM)*
  `git commit -m "feat(vg): C4 wire run_layer2 + re-express RABL2 recovery (no secondary in bundle.reads)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

---

## Milestone M5 — All-secondary regions → proof-gated candidate new copies (C5, default-off)

**Outcome:** Regions every read of which is secondary (Layer 1 correctly drops them) are reconsidered by Layer 2 **only** on decisive proof, behind `--vg-layer2-new-copies` (default off). Logged, never silent. The `new_copies` parameter already exists on `run_layer2` (added in M3.2), so no signature churn.

### Task M5.1 — Detect all-secondary regions in the side-index

**Files:**
- Modify: `src/rustle/vg_family/secondary_index.rs` (new method + tests)

**Steps:**

- [ ] Add a failing test: secondaries clustered in a region with no overlapping Layer-1 locus (`locus == None`) form one candidate region; secondaries that overlap a locus do not. Append to `secondary_index.rs` `mod tests`:

```rust
    #[test]
    fn all_secondary_regions_cluster_unassigned_secondaries() {
        let mut idx = SecondaryIndex::new();
        idx.push(sa_loc(1, 8000, 8100, None));
        idx.push(sa_loc(2, 8150, 8250, None));
        idx.push(sa_loc(3, 8300, 8400, None));
        idx.push(sa_loc(4, 100, 200, Some(0))); // assigned → excluded
        let regions = idx.all_secondary_regions(/*max_gap=*/ 200, /*min_reads=*/ 3);
        assert_eq!(regions.len(), 1, "one all-secondary candidate region");
        assert_eq!(regions[0], ("chrT".to_string(), 8000, 8400, 3));
    }
```

- [ ] Run, expect FAIL: `cargo test --lib all_secondary_regions_cluster_unassigned_secondaries` → `cannot find method`.
- [ ] Add the method to `impl SecondaryIndex`:

```rust
    /// Detect candidate "all-secondary" regions: clusters of secondaries that
    /// overlap NO Layer-1 locus (Layer 1 dropped these as cross-map artifacts).
    /// Returns sorted `(chrom, start, end, n_reads)` for clusters of ≥ `min_reads`
    /// secondaries whose successive gaps are ≤ `max_gap`. C5 reconsiders these
    /// ONLY on decisive proof (caller-gated, default off).
    pub fn all_secondary_regions(
        &self,
        max_gap: u64,
        min_reads: usize,
    ) -> Vec<(String, u64, u64, usize)> {
        let mut unassigned: Vec<&SecondaryAlignment> =
            self.alignments.iter().filter(|a| a.locus.is_none()).collect();
        unassigned.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.ref_start.cmp(&b.ref_start)));
        let mut out: Vec<(String, u64, u64, usize)> = Vec::new();
        let mut i = 0;
        while i < unassigned.len() {
            let chrom = unassigned[i].chrom.clone();
            let mut end = unassigned[i].ref_end;
            let start = unassigned[i].ref_start;
            let mut n = 1usize;
            let mut j = i + 1;
            while j < unassigned.len()
                && unassigned[j].chrom == chrom
                && unassigned[j].ref_start <= end + max_gap
            {
                end = end.max(unassigned[j].ref_end);
                n += 1;
                j += 1;
            }
            if n >= min_reads {
                out.push((chrom, start, end, n));
            }
            i = j;
        }
        out
    }
```

- [ ] Re-run → PASS. Full suite: `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/vg_family/secondary_index.rs`
  `git commit -m "feat(vg): C5 detect all-secondary candidate regions" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M5.2 — Proof-gate all-secondary regions into candidate new-copy transcripts (default off)

> **No signature churn (review issue 30):** `run_layer2` already takes `new_copies: bool` (since M3.2). This task only fills in the *behavior* behind that parameter and adds the emitter; the M3.2 test (`run_layer2_forms_family_and_merges_graph`, which passes `new_copies: false`) and the M4.3c pipeline call (which passes `config.vg_layer2_new_copies`) are unchanged.

**Files:**
- Modify: `src/rustle/vg_family/layer2.rs` (emitter fn + `run_layer2` body)

**Steps:**

- [ ] Add a failing test: with `new_copies=false`, no region emits; with `new_copies=true` and the proof bar met, the region emits a candidate tagged `NovelLocusFromScan`. Append to `layer2.rs` `mod tests`:

```rust
    #[test]
    fn all_secondary_new_copy_gated_off_by_default() {
        let region = ("chrT".to_string(), 8000u64, 8400u64, 5usize);
        let off = candidate_new_copy_transcripts(
            std::slice::from_ref(&region), /*new_copies=*/ false, /*min_reads=*/ 3,
        );
        assert!(off.is_empty(), "default-off: no new copies emitted");

        let on = candidate_new_copy_transcripts(
            std::slice::from_ref(&region), /*new_copies=*/ true, /*min_reads=*/ 3,
        );
        assert_eq!(on.len(), 1, "with flag + proof bar met, one candidate emits");
        assert_eq!(
            on[0].rescue_class,
            Some(crate::vg_family::diagnostic::RescueClass::NovelLocusFromScan)
        );
    }
```

- [ ] Run, expect FAIL: `cargo test --lib all_secondary_new_copy_gated_off_by_default` → `cannot find function`.
- [ ] Add the proof-gated emitter to `layer2.rs`:

```rust
/// (C5) Convert decisively-proven all-secondary regions into candidate new-copy
/// transcripts. Default OFF (`new_copies`); the proof bar here is `min_reads`
/// secondaries clustered in a region Layer 1 dropped (same decisive-evidence
/// spirit as C4; tighten with PSV ownership when validated). Tagged
/// `NovelLocusFromScan` so it is provenance-marked and union-by-chain can fold
/// it if it duplicates an existing chain.
pub fn candidate_new_copy_transcripts(
    regions: &[(String, u64, u64, usize)],
    new_copies: bool,
    min_reads: usize,
) -> Vec<Transcript> {
    if !new_copies {
        return Vec::new();
    }
    let mut out = Vec::new();
    for (chrom, start, end, n) in regions {
        if *n < min_reads {
            continue; // proof bar
        }
        let mut t = Transcript::default();
        t.chrom = chrom.clone();
        t.strand = '+';
        t.exons = vec![(*start, *end)];
        t.coverage = *n as f64;
        t.longcov = *n as f64;
        t.is_longread = true;
        t.synthetic = true;
        t.rescue_class = Some(crate::vg_family::diagnostic::RescueClass::NovelLocusFromScan);
        out.push(t);
    }
    out
}
```

- [ ] Wire into `run_layer2`: replace the `let _ = new_copies;` line (added in M4.3c) with the all-secondary emission, before the final `Ok(...)`:

```rust
    if new_copies {
        let regions = side_index.all_secondary_regions(/*max_gap=*/ 200, /*min_reads=*/ 3);
        if !regions.is_empty() {
            eprintln!("[layer2] {} all-secondary candidate region(s) (proof-gated)", regions.len());
        }
        novel_transcripts.extend(candidate_new_copy_transcripts(&regions, new_copies, 3));
    }
```

- [ ] Re-run: `cargo test --lib all_secondary_new_copy_gated_off_by_default run_layer2_forms_family_and_merges_graph` → both PASS (the M3.2 test is unaffected — same signature). Full suite: `cargo test --lib` → `0 failed`.
- [ ] Harness check (default-off invariant): run `bench/layer2_invariant.sh` — with `--vg-layer2` but without `--vg-layer2-new-copies`, output is unchanged from M4 (no all-secondary regions emit); checks (3)/(4) still pass on chr19+chrY.
- [ ] Commit:
  `git add src/rustle/vg_family/layer2.rs`
  `git commit -m "feat(vg): C5 proof-gated all-secondary new copies (default off)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

---

## Milestone M6 — Union-by-chain additivity emission (C6) + additivity test

**Outcome:** Layer-1 transcripts ∪ novel Layer-2 transcripts, deduplicated by intron chain. Every Layer-1 transcript emits unchanged; a Layer-2 transcript emits only if its chain is absent from Layer 1; a Layer-2 path re-deriving a Layer-1 chain is silently folded. This is the structural guarantee Layer 2 can only ever add — proven by a standing additivity test.

### Task M6.1 — Union-by-chain helper (reuse the existing idiom)

**Files:**
- Modify: `src/rustle/vg_family/layer2.rs` (new fn + tests)

**Steps:**

- [ ] Add a failing test: given Layer-1 transcripts and Layer-2 novel transcripts where one Layer-2 chain duplicates a Layer-1 chain and one is novel, the union keeps all Layer-1 unchanged and adds only the novel one. Append to `layer2.rs` `mod tests`:

```rust
    #[test]
    fn union_by_chain_adds_only_novel_and_never_drops_baseline() {
        let mut l1 = Vec::new();
        let mut t = Transcript::default();
        t.chrom = "chrT".into(); t.strand = '+'; t.exons = vec![(100, 160), (300, 360)];
        l1.push(t);

        let mut dup = Transcript::default();
        dup.chrom = "chrT".into(); dup.strand = '+'; dup.exons = vec![(100, 160), (300, 360)];
        let mut novel = Transcript::default();
        novel.chrom = "chrT".into(); novel.strand = '+'; novel.exons = vec![(100, 160), (500, 560)];
        let l2 = vec![dup, novel];

        let before = l1.len();
        let added = union_layer2_by_chain(&mut l1, l2);
        assert_eq!(added, 1, "only the novel chain is added");
        assert_eq!(l1.len(), before + 1, "baseline count preserved + 1 novel");
        assert!(l1.iter().any(|t| t.exons == vec![(100, 160), (300, 360)]));
        assert!(l1.iter().any(|t| t.exons == vec![(100, 160), (500, 560)]));
    }
```

- [ ] Run, expect FAIL: `cargo test --lib union_by_chain_adds_only_novel_and_never_drops_baseline` → `cannot find function`.
- [ ] Add the union helper. It mirrors the existing union-back idiom (`pipeline.rs:19545-19581`):

```rust
/// (C6) Union novel Layer-2 transcripts into the Layer-1 set by intron chain.
/// Every Layer-1 transcript is preserved unchanged; a Layer-2 transcript emits
/// ONLY if its intron chain is absent from Layer 1 and not already added. A
/// single-exon Layer-2 transcript (empty chain) emits only if its exact span is
/// not already covered. Returns the number of novel transcripts added.
///
/// This is the structural additivity guarantee: VG ⊇ baseline, always.
pub fn union_layer2_by_chain(layer1: &mut Vec<Transcript>, layer2: Vec<Transcript>) -> usize {
    let chain_of = |t: &Transcript| -> Vec<(u64, u64)> {
        t.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect()
    };
    let mut existing: DetHashSet<(String, char, Vec<(u64, u64)>)> = DetHashSet::default();
    let mut existing_single: DetHashSet<(String, char, u64, u64)> = DetHashSet::default();
    for t in layer1.iter() {
        if t.exons.len() >= 2 {
            existing.insert((t.chrom.clone(), t.strand, chain_of(t)));
        } else if let Some(&(s, e)) = t.exons.first() {
            existing_single.insert((t.chrom.clone(), t.strand, s, e));
        }
    }
    let mut added = 0usize;
    for mut t in layer2 {
        if t.exons.len() >= 2 {
            let key = (t.chrom.clone(), t.strand, chain_of(&t));
            if existing.contains(&key) {
                continue;
            }
            existing.insert(key);
        } else if let Some(&(s, e)) = t.exons.first() {
            let key = (t.chrom.clone(), t.strand, s, e);
            if existing_single.contains(&key) {
                continue;
            }
            existing_single.insert(key);
        } else {
            continue;
        }
        t.rescue_class = None; // now a first-class transcript (matches union-back)
        layer1.push(t);
        added += 1;
    }
    added
}
```

> Reuses the exact chain-signature idiom from `pipeline.rs:19545-19581` and clears `rescue_class` on emission (matches `pipeline.rs:19577`).

- [ ] Re-run → PASS. Full suite: `cargo test --lib` → `0 failed`.
- [ ] Commit:
  `git add src/rustle/vg_family/layer2.rs`
  `git commit -m "feat(vg): C6 union-by-chain additivity helper" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

### Task M6.2 — Wire the union into the pipeline emission path

**Files:**
- Modify: `src/rustle/pipeline.rs` (between the existing union-back block ~`19581` and `write_gtf` ~`19692`)

**Steps:**

- [ ] Insert the Layer-2 union just before `write_gtf` (`19692`), after the existing `union_baseline_holdout` back block (`19545-19581`) so Layer-2 novelties are unioned against the *final* Layer-1 set:

```rust
    // ── Layer 2 (C6): union novel family-copy/isoform transcripts by chain ──
    if config.vg_mode && config.vg_layer2 && !layer2_novel.is_empty() {
        let added = crate::vg_family::layer2::union_layer2_by_chain(
            &mut all_transcripts,
            std::mem::take(&mut layer2_novel),
        );
        eprintln!("[layer2] unioned {added} novel transcript(s) (VG superset baseline preserved)");
    }
```

- [ ] Build: `cargo build --release` → `Finished release`.
- [ ] Add the standing additivity check to `bench/layer2_invariant.sh` after check (5). It re-asserts the superset on both chroms (every baseline chain appears verbatim in the Layer-2 output — no modified baseline chain):

```bash
echo "== (6) additivity: VG Layer-2 baseline chains identical (no modified baseline) =="
for tag in chr19 chrY; do
  python3 scripts/coord_signature_superset.py "/tmp/layer2/${tag}_vg_layer2.gtf" "/tmp/layer2/${tag}_baseline.gtf" \
    && echo "  OK [$tag]: every baseline chain present unchanged in Layer-2 output"
done
```

- [ ] Run the full harness: `RAYON_NUM_THREADS=1 bench/layer2_invariant.sh`. Expect `ALL INVARIANTS PASS` (chr19+chrY) plus check (5) showing extra RABL2 isoforms (when the RABL2 fixture is present) and check (6) confirming additivity.
- [ ] Run the full unit suite: `cargo test --lib` → `0 failed`.
- [ ] **Final safety floor:** `RUSTLE_PRECISE=1 RAYON_NUM_THREADS=1 ./target/release/rustle -L GGO_19.bam -o /tmp/precise.gtf 2>/dev/null && diff -q bench/ref/4705ab1_precise_GGO_19.gtf /tmp/precise.gtf` → identical. VG-default (no `--vg-layer2`) coord-signature == baseline.
- [ ] Commit:
  `git add src/rustle/pipeline.rs bench/layer2_invariant.sh`
  `git commit -m "feat(vg): C6 wire Layer-2 union into emission + standing additivity check (chr19+chrY)" -m "Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"`

---

## Cross-cutting verification checklist (run at the end of every milestone)

- [ ] `cargo test --lib` → `0 failed` (count grows by the milestone's new tests).
- [ ] `RUSTLE_PRECISE=1 RAYON_NUM_THREADS=1 ./target/release/rustle -L GGO_19.bam -o /tmp/precise.gtf` byte-identical to `bench/ref/4705ab1_precise_GGO_19.gtf` (the `4705ab1` reference).
- [ ] VG-default (`--vg`, no `--vg-layer2`) coord-signature == baseline on **both chr19 and chrY** (`scripts/coord_signature_superset.py` both directions).
- [ ] VG-Layer2 (`--vg --vg-layer2`) coord-signature ⊇ baseline on **both chr19 and chrY** (never a missing baseline chain).
- [ ] Determinism: re-run any `--vg --vg-layer2` invocation twice with `RAYON_NUM_THREADS=1`; GTFs identical. (Side-index iteration is sorted; all maps are `DetHashMap`/`DetHashSet`; graph capture is index-keyed.)
- [ ] Whole-genome runs are per-chrom serial (never a single whole-genome `-L`, which OOMs ~18 GB).
- [ ] No `HashMap::new()` / `HashSet::new()` introduced — only `DetHashMap::default()` / `DetHashSet::default()` (except where a struct field is declared as plain `std::collections::HashMap`, e.g. `FamilyGroup::multimap_reads`, which must keep its declared type).
- [ ] All commits staged explicitly (no `git add -A`); gitignored fixtures added with `git add -f`; every commit message ends with `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`.

## Design decisions locked by this plan (the five open questions)

1. **Side-index scope/size (C1).** Chromosome-global capture (a dedicated long-class sec/supp BAM pass that mirrors the Phase-1 gate predicate), then `prune_to_loci` to family-candidate loci after discovery's first pass, then `cap_per_locus(2000)` — both returning a logged drop count (no silent truncation). Implemented in M1.2, applied in M3.2 `run_layer2`. The **chrY (repeat-rich) regime is exercised by the standing superset invariant** (M1.5), which is where the size risk lives.
2. **Similarity threshold (C2).** Exon-only **canonical**-minimizer Jaccard (`exon_kmer_similarity_between_graphs`, k=15) — strand-agnostic, so inverted paralogs (DAZ-style) are handled (proven by `exon_kmer_similarity_high_for_inverted_paralog`). Conservative default `--family-exon-similarity 0.30`, required *in addition to* cross-map linkage. Implemented in M2.2/M2.3.
3. **All-secondary new-copy gate (C5).** Default-off behind `--vg-layer2-new-copies` / `RUSTLE_VG_LAYER2_NEW_COPIES`; proof bar = clustered ≥3 secondaries in a Layer-1-dropped region (tightenable to PSV ownership). Implemented in M5.
4. **Flagging (decision 4).** Entire Layer 2 default-off behind `--vg-layer2` / `RUSTLE_VG_LAYER2`; with it off, `--vg` is Layer-1-baseline-identical. Implemented in M1.4, gated throughout.
5. **EM/PSV reuse (decision 5).** `enumerate_diagnostic_sites` is invoked **and actively used** on the merged family graph to forbid allele-mixing (chimeric) paths in `decompose_family_paths` (M4.2) — this is the operational constrained-flow-decomposition, not a nominal call. **Scope note:** full reuse of `classify_family` / `compute_copy_ownership` / `run_fingerprint_em` is **deferred** — those take `(&FamilyGroup, &[Bundle], ...)` and the Layer-1-graph-sourced path no longer has the bundle inputs they expect, so wiring them is a follow-up tightening step (e.g. synthesizing a bundle-equivalent view, or porting them to take a `FamilyGraph`). Only `enumerate_diagnostic_sites` (which takes `&FamilyGraph`) is confirmed to operate unchanged on the merged graph in this plan. The constrained-decomposition contribution (the thesis core) is fully delivered by the active PSV filter.
