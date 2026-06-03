# Intra-Bundle Tandem VG (Local Variation Graph) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recover per-copy transcripts from a tandem paralog array that collapses into a single bundle (e.g. RBMY1), by decomposing the mega-bundle into a synthetic per-copy `FamilyGroup` and feeding the existing (already tandem-aware) family-graph + fingerprint-EM + completion machinery.

**Architecture:** The new code is narrow — detect a tandem mega-bundle (≥2 genomically-separated, mutually sequence-similar read clusters) and split it into synthetic per-copy sub-bundles forming a `FamilyGroup`. Everything downstream is an UNCHANGED call to existing functions: `build_family_graph` (family_graph.rs:411), `annotate_per_copy_exon_coverage` (vg.rs:1406), `run_fingerprint_em` (vg.rs), and the default-on completion/borrow. The family graph over the decomposed copies IS the local variation graph (shared ExonClasses = tandem cycle; copy-specific exons = bubbles). Abundance/abstention come from the flow-capacity EM. All `--vg`-gated behind `RUSTLE_VG_TANDEM` (default OFF); default de-novo stays byte-identical.

**Tech Stack:** Rust (rustle), `FamilyGraph`/fingerprint-EM (`src/rustle/vg.rs`, `src/rustle/vg_hmm/`), `GenomeIndex` (`src/rustle/genome.rs`), cargo test.

**Spec:** `docs/superpowers/specs/2026-06-02-intrabundle-tandem-vg-design.md` (draft 2026-06-02).

**Depends on:** the flow-capacity apportionment EM on this branch (`vg/flow-capacity-apportionment`) — it supplies per-copy `read.weight`, `em_anchored`, and `capacity_confidence` (the abstention channel).

**Task dependency order:** T1 `tandem_copy` field → T2 GTF emit (foundation) → T3 `cluster_reads_by_position` → T4 `cluster_kmer_jaccard` (pure helpers) → T5 `TandemConfig` env → T6 `detect_tandem_bundle` (uses T3/T4/T5) → T7 `decompose_tandem_to_family` → T8 pipeline wiring → T9 validation. Implement in order.

---

### Task 1: Add `tandem_copy` field to `Transcript`

**Files:**
- Modify: `src/rustle/path_extract.rs:620` (the `Transcript` struct; add field after `family_verdict`)
- Modify: 37 `Transcript` literal sites (compiler-enforced): `src/rustle/path_extract.rs`, `src/rustle/pipeline.rs`, `src/rustle/transcript_filter.rs`, `src/rustle/merge_mode.rs`, `src/rustle/single_exon_pileup.rs`, `src/rustle/vg.rs`

Notes (verified): `Transcript` has no `Default` impl and no `..spread` literals — the compiler enumerates every missing-field site. Every literal already contains `capacity_confidence:`. New field defaults to `None` everywhere; only Task 8 sets `Some(true)`.

- [ ] **Step 1: Add the field to the struct (after the `family_verdict` field)**

```rust
    /// O5/tandem VG (spec 2026-06-02): `Some(true)` when this transcript is a
    /// copy resolved from an intra-bundle tandem decomposition (the array
    /// collapsed into one bundle; VG split it into synthetic per-copy
    /// sub-bundles). Provenance flag — pairs with `capacity_confidence` for the
    /// borrowed/abstained portion. `None` for ordinary transcripts. VG-only.
    pub tandem_copy: Option<bool>,
```

- [ ] **Step 2: Build to enumerate the literals**

Run: `cargo build 2>&1 | grep -E "missing field tandem_copy" | head -50`
Expected: ~37 `error[E0063]: missing field tandem_copy` lines.

- [ ] **Step 3: At each reported site, insert `tandem_copy: None,`** (next to `family_verdict: …,`). All default `None`. Example:
```rust
                    capacity_confidence: cap_conf_opt, abundance_min: abund_min_opt, family_verdict: None, tandem_copy: None, intron_low: Vec::new(), synthetic: false, rescue_class: None,
```

- [ ] **Step 4: Build clean**

Run: `cargo build 2>&1 | grep -E "^error" | head`
Expected: no output.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "feat(vg): add Transcript.tandem_copy field"
```

---

### Task 2: Emit `tandem_copy` in GTF

**Files:**
- Modify: `src/rustle/gtf.rs:202` (after the `capacity_confidence` emission block)
- Test: `src/rustle/gtf.rs` test module

- [ ] **Step 1: Write the failing test** (mirror the existing `capacity_confidence` test style at gtf.rs:307)

```rust
    #[test]
    fn tandem_copy_attr_emitted_when_some_true() {
        let mut tx = make_test_transcript();
        tx.vg_copy_id = Some(0);
        tx.tandem_copy = Some(true);
        let line = format_transcript_attrs(&tx, 0.05);
        assert!(line.contains("tandem_copy \"true\";"), "got: {line}");
    }

    #[test]
    fn tandem_copy_attr_absent_otherwise() {
        let mut tx = make_test_transcript();
        tx.tandem_copy = None;
        assert!(!format_transcript_attrs(&tx, 0.05).contains("tandem_copy"));
    }
```

(If the cc tests use inline `tx` construction + a different formatter call, match that exact pattern instead.)

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p rustle --lib gtf::tests::tandem_copy -- --nocapture`
Expected: FAIL — attribute not emitted.

- [ ] **Step 3: Add the emission block (after the `capacity_confidence` block, gtf.rs:~202)**

```rust
        // O5/tandem provenance flag (--vg only; spec 2026-06-02). VG-only so
        // non-VG GTF is byte-identical.
        if tx.tandem_copy == Some(true) {
            tx_attrs.push_str(" tandem_copy \"true\";");
        }
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p rustle --lib gtf::tests::tandem_copy`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/gtf.rs
git commit -m "feat(vg): emit tandem_copy GTF attribute"
```

---

### Task 3: Pure helper `cluster_reads_by_position`

Splits a bundle's read footprints into genomically-separated clusters (the copy loci). Pure 1-D gap clustering — trivially unit-testable.

**Files:**
- Create: `src/rustle/vg_hmm/tandem.rs`
- Modify: `src/rustle/vg_hmm/mod.rs` (add `pub mod tandem;`)
- Test: `src/rustle/vg_hmm/tandem.rs` (inline)

- [ ] **Step 1: Create the module with the function + failing tests**

```rust
//! Intra-bundle tandem decomposition (spec 2026-06-02). Detects a tandem
//! paralog array that collapsed into one bundle and splits it into per-copy
//! clusters so the existing family-graph/EM machinery can resolve the copies.

/// Cluster `[ (ref_start, ref_end) ]` read footprints into genomically-ordered
/// groups separated by gaps `> min_gap`. Returns, per cluster, the index list
/// into `intervals` plus the cluster's (start,end) span. Genomic-ascending.
pub fn cluster_reads_by_position(
    intervals: &[(u64, u64)],
    min_gap: u64,
) -> Vec<(Vec<usize>, (u64, u64))> {
    if intervals.is_empty() {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..intervals.len()).collect();
    order.sort_by_key(|&i| intervals[i].0);

    let mut clusters: Vec<(Vec<usize>, (u64, u64))> = Vec::new();
    for &i in &order {
        let (s, e) = intervals[i];
        match clusters.last_mut() {
            Some((idxs, span)) if s <= span.1.saturating_add(min_gap) => {
                idxs.push(i);
                span.1 = span.1.max(e);
                span.0 = span.0.min(s);
            }
            _ => clusters.push((vec![i], (s, e))),
        }
    }
    clusters
}
```

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_two_separated_copies() {
        // two copies ~23kb apart, gap threshold 5kb -> 2 clusters.
        let ivs = [(1000, 14000), (1500, 14500), (37000, 50000), (37200, 50200)];
        let c = cluster_reads_by_position(&ivs, 5000);
        assert_eq!(c.len(), 2);
        assert_eq!(c[0].0.len(), 2);
        assert_eq!(c[1].0.len(), 2);
        assert_eq!(c[0].1, (1000, 14500));
        assert_eq!(c[1].1, (37000, 50200));
    }

    #[test]
    fn single_cluster_when_within_gap() {
        let ivs = [(1000, 14000), (16000, 30000)]; // 2kb gap < 5kb
        assert_eq!(cluster_reads_by_position(&ivs, 5000).len(), 1);
    }

    #[test]
    fn empty_input() {
        assert!(cluster_reads_by_position(&[], 5000).is_empty());
    }
}
```

- [ ] **Step 2: Wire the module**

Add to `src/rustle/vg_hmm/mod.rs`: `pub mod tandem;`

- [ ] **Step 3: Run**

Run: `cargo test -p rustle --lib vg_hmm::tandem::tests::`
Expected: PASS (3 tests).

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg_hmm/tandem.rs src/rustle/vg_hmm/mod.rs
git commit -m "feat(vg): cluster_reads_by_position (tandem decomposition)"
```

---

### Task 4: Pure helper `cluster_kmer_jaccard`

Exact k-mer Jaccard between two sequences — the within-bundle cross-copy similarity test. Pure, unit-testable.

**Files:**
- Modify: `src/rustle/vg_hmm/tandem.rs`
- Test: same file

- [ ] **Step 1: Add the function + failing tests**

```rust
use std::collections::HashSet;

/// Exact k-mer Jaccard between two DNA sequences (ACGT only; other bases
/// skipped). Returns 0.0 if either has no valid k-mers. Clusters are short
/// (one copy's exonic sequence), so an exact set is fine — no min-hash needed.
pub fn cluster_kmer_jaccard(a: &[u8], b: &[u8], k: usize) -> f64 {
    fn kset(s: &[u8], k: usize) -> HashSet<&[u8]> {
        let mut set = HashSet::new();
        if s.len() >= k {
            for w in s.windows(k) {
                if w.iter().all(|&c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
                    set.insert(w);
                }
            }
        }
        set
    }
    let sa = kset(a, k);
    let sb = kset(b, k);
    if sa.is_empty() || sb.is_empty() {
        return 0.0;
    }
    let inter = sa.intersection(&sb).count() as f64;
    let union = sa.union(&sb).count() as f64;
    inter / union
}
```

```rust
    #[test]
    fn identical_sequences_jaccard_one() {
        let s = b"ACGTACGTACGTACGTACGT";
        assert!((cluster_kmer_jaccard(s, s, 15) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn near_identical_high_jaccard() {
        let a = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut b = a.to_vec();
        b[16] = b'T'; // one substitution
        assert!(cluster_kmer_jaccard(a, &b, 15) > 0.2);
    }

    #[test]
    fn disjoint_sequences_jaccard_zero() {
        let a = b"AAAAAAAAAAAAAAAAAAAA";
        let b = b"CCCCCCCCCCCCCCCCCCCC";
        assert_eq!(cluster_kmer_jaccard(a, b, 15), 0.0);
    }
```

- [ ] **Step 2: Run**

Run: `cargo test -p rustle --lib vg_hmm::tandem::tests::`
Expected: PASS (6 tests total).

- [ ] **Step 3: Commit**

```bash
git add src/rustle/vg_hmm/tandem.rs
git commit -m "feat(vg): cluster_kmer_jaccard (tandem similarity)"
```

---

### Task 5: `TandemConfig` env reader

**Files:**
- Modify: `src/rustle/vg_hmm/tandem.rs`
- Test: same file

- [ ] **Step 1: Add config + failing test**

```rust
/// O5/tandem runtime config (env). `enabled` defaults OFF (opt-in prototype).
#[derive(Debug, Clone, Copy)]
pub struct TandemConfig {
    pub enabled: bool,
    pub min_gap: u64,       // genomic separation to split copies
    pub min_jaccard: f64,   // cross-copy sequence-similarity floor
    pub kmer: usize,
}

impl TandemConfig {
    pub fn from_env() -> Self {
        let enabled = std::env::var("RUSTLE_VG_TANDEM")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false);
        let min_gap = std::env::var("RUSTLE_VG_TANDEM_MIN_SEP")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(5_000);
        let min_jaccard = std::env::var("RUSTLE_VG_TANDEM_MIN_JACCARD")
            .ok().and_then(|v| v.parse().ok())
            .filter(|f: &f64| *f > 0.0 && *f <= 1.0).unwrap_or(0.20);
        TandemConfig { enabled, min_gap, min_jaccard, kmer: 15 }
    }
}
```

```rust
    #[test]
    fn tandemconfig_defaults_off() {
        let c = TandemConfig::from_env();
        assert!(!c.enabled);
        assert_eq!(c.min_gap, 5_000);
        assert!((c.min_jaccard - 0.20).abs() < 1e-9);
    }
```

- [ ] **Step 2: Run** (do NOT set the env vars in the test)

Run: `cargo test -p rustle --lib vg_hmm::tandem::tests::tandemconfig`
Expected: PASS.

- [ ] **Step 3: Commit**

```bash
git add src/rustle/vg_hmm/tandem.rs
git commit -m "feat(vg): TandemConfig env reader (RUSTLE_VG_TANDEM default OFF)"
```

---

### Task 6: `detect_tandem_bundle`

Orchestrates T3/T4 against a bundle + genome. Returns the per-copy spans when the bundle is a tandem array, else `None`.

**Files:**
- Modify: `src/rustle/vg_hmm/tandem.rs`
- Reference (read): `src/rustle/types.rs` (`Bundle`, `BundleRead.exons`/ref start-end), `src/rustle/genome.rs` (`GenomeIndex::fetch_sequence`)
- Test: same file (None-path units; full positive path is the RBMY integration test in Task 9)

- [ ] **Step 1: Add the type + function**

```rust
use crate::types::Bundle;
use crate::genome::GenomeIndex;

/// A tandem array detected inside one bundle: the per-copy genomic spans
/// (genomic-ascending) and, for each copy, the read indices (into `bundle.reads`)
/// that fall in that copy's span.
pub struct TandemDecomposition {
    pub copy_spans: Vec<(u64, u64)>,
    pub copy_read_idxs: Vec<Vec<usize>>,
}

/// Detect whether `bundle` is a collapsed tandem array. Requires ≥2 position
/// clusters that are mutually sequence-similar (≥ `cfg.min_jaccard`). Returns
/// `None` when disabled, single-cluster, or clusters are dissimilar (a normal
/// multi-exon gene — must NOT be split).
pub fn detect_tandem_bundle(
    bundle: &Bundle,
    genome: &GenomeIndex,
    cfg: TandemConfig,
) -> Option<TandemDecomposition> {
    if !cfg.enabled {
        return None;
    }
    // Read footprints = (ref_start, ref_end) per read.
    let ivs: Vec<(u64, u64)> = bundle.reads.iter()
        .map(|r| (r.ref_start, r.ref_end))
        .collect();
    let clusters = cluster_reads_by_position(&ivs, cfg.min_gap);
    if clusters.len() < 2 {
        return None;
    }
    // Fetch each cluster's reference sequence (the copy's genomic span).
    let seqs: Vec<Vec<u8>> = clusters.iter()
        .map(|(_, (s, e))| genome.fetch_sequence(&bundle.chrom, *s, *e).unwrap_or_default())
        .collect();
    // Require mutual similarity: every cluster must be ≥min_jaccard to at least
    // one other (so a lone unrelated cluster doesn't gate, but a single gene with
    // one weird exon block isn't promoted either).
    let n = clusters.len();
    let mut similar_to_any = vec![false; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if cluster_kmer_jaccard(&seqs[i], &seqs[j], cfg.kmer) >= cfg.min_jaccard {
                similar_to_any[i] = true;
                similar_to_any[j] = true;
            }
        }
    }
    // Keep only similar clusters; need ≥2 to call a tandem array.
    let kept: Vec<usize> = (0..n).filter(|&i| similar_to_any[i]).collect();
    if kept.len() < 2 {
        return None;
    }
    Some(TandemDecomposition {
        copy_spans: kept.iter().map(|&i| clusters[i].1).collect(),
        copy_read_idxs: kept.iter().map(|&i| clusters[i].0.clone()).collect(),
    })
}
```

(Confirm the exact field names `ref_start`/`ref_end` on `BundleRead` when reading types.rs in this step; adjust if they differ, e.g. `start`/`end`.)

- [ ] **Step 2: Add None-path tests** (positive path needs a genome → covered by Task 9)

```rust
    #[test]
    fn detect_disabled_returns_none() {
        // enabled=false short-circuits before touching genome; pass a dummy via
        // the existing test GenomeIndex constructor if available, else assert the
        // early return by constructing TandemConfig{enabled:false,..} and an empty
        // bundle through the existing vg.rs test bundle helper.
        let cfg = TandemConfig { enabled: false, min_gap: 5000, min_jaccard: 0.2, kmer: 15 };
        // build a minimal bundle via the existing test helper (see vg.rs tests)
        let b = crate::vg::tests::make_test_bundle_empty();
        let g = crate::genome::tests::empty_genome();
        assert!(detect_tandem_bundle(&b, &g, cfg).is_none());
    }
```

If `make_test_bundle_empty`/`empty_genome` helpers don't exist, add tiny `#[cfg(test)]` constructors in those modules (an empty `Bundle` literal and a `GenomeIndex` over a 1-contig in-memory sequence) and reference them. Keep them minimal.

- [ ] **Step 3: Build + run**

Run: `cargo test -p rustle --lib vg_hmm::tandem`
Expected: PASS (build clean; None-path test passes).

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "feat(vg): detect_tandem_bundle (collapsed tandem-array detector)"
```

---

### Task 7: `decompose_tandem_to_family`

Builds synthetic per-copy sub-bundles + a `FamilyGroup` from a detection.

**Files:**
- Modify: `src/rustle/vg_hmm/tandem.rs`
- Reference (read): `src/rustle/vg.rs` (`FamilyGroup`), `src/rustle/bundle.rs` (`recompute_junction_stats`)
- Test: same file

- [ ] **Step 1: Add the function**

```rust
use crate::vg::FamilyGroup;

/// Turn a detection into synthetic per-copy sub-bundles (appended to `bundles`)
/// + a `FamilyGroup` linking them. Each sub-bundle = the parent cloned, with
/// `reads` restricted to that copy's reads, span set to the copy span, and
/// `synthetic=true`. `multimap_reads` links reads whose name appears in >1 copy
/// (the secondaries that caused the collapse — now cross-copy evidence).
pub fn decompose_tandem_to_family(
    parent: &Bundle,
    decomp: &TandemDecomposition,
    bundles: &mut Vec<Bundle>,
    family_id: usize,
    config: &crate::types::RunConfig,
) -> (FamilyGroup, Vec<usize>) {
    let mut bundle_indices = Vec::new();
    // Track read_name_hash -> Vec<(copy_pos, read_idx_in_subbundle)> for multimap.
    let mut name_to_copies: std::collections::HashMap<u64, Vec<(usize, usize)>> =
        std::collections::HashMap::new();

    for (copy_pos, (span, read_idxs)) in
        decomp.copy_spans.iter().zip(decomp.copy_read_idxs.iter()).enumerate()
    {
        let mut sub = parent.clone();
        sub.start = span.0;
        sub.end = span.1;
        sub.synthetic = true;
        sub.vg_family_id = Some(family_id);
        sub.reads = read_idxs.iter().map(|&i| parent.reads[i].clone()).collect();
        for (ri, r) in sub.reads.iter().enumerate() {
            name_to_copies.entry(r.read_name_hash).or_default().push((copy_pos, ri));
        }
        crate::bundle::recompute_junction_stats(&mut sub, config);
        bundle_indices.push(bundles.len());
        bundles.push(sub);
    }

    // A read shared across ≥2 copies is a multimapper for the EM.
    let multimap_reads: std::collections::HashMap<u64, Vec<(usize, usize)>> =
        name_to_copies.into_iter().filter(|(_, v)| v.len() > 1).collect();

    (FamilyGroup { family_id, bundle_indices: bundle_indices.clone(), multimap_reads },
     bundle_indices)
}
```

- [ ] **Step 2: Add a unit test** (build the parent via the existing vg.rs test helpers)

```rust
    #[test]
    fn decompose_builds_one_subbundle_per_copy() {
        let parent = crate::vg::tests::make_test_bundle_two_copies(); // helper: reads in 2 position clusters
        let decomp = TandemDecomposition {
            copy_spans: vec![(1000, 14000), (37000, 50000)],
            copy_read_idxs: vec![vec![0, 1], vec![2, 3]],
        };
        let mut bundles: Vec<crate::types::Bundle> = Vec::new();
        let cfg = crate::types::RunConfig::default();
        let (fam, idxs) = decompose_tandem_to_family(&parent, &decomp, &mut bundles, 0, &cfg);
        assert_eq!(idxs.len(), 2);
        assert_eq!(bundles.len(), 2);
        assert_eq!(bundles[0].reads.len(), 2);
        assert!(bundles[0].synthetic);
        assert_eq!(fam.bundle_indices, idxs);
    }
```

If `make_test_bundle_two_copies` / `RunConfig::default()` aren't available, add a minimal `#[cfg(test)]` helper (4 reads via the existing `make_read`, 2 in each span) and build the `RunConfig` the way other vg.rs tests do.

- [ ] **Step 3: Run**

Run: `cargo test -p rustle --lib vg_hmm::tandem`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add -A
git commit -m "feat(vg): decompose_tandem_to_family (synthetic per-copy bundles)"
```

---

### Task 8: Wire the tandem pre-pass into the VG pipeline

Run detection+decomposition before family discovery; replace each tandem mega-bundle with its synthetic sub-bundles so the existing family path resolves the copies.

**Files:**
- Modify: `src/rustle/pipeline.rs` — just before the family-discovery block (`discover_family_groups`, ~10247) where `bundles` and the VG genome are in scope.

- [ ] **Step 1: Read the surrounding code to fix names/insertion point**

Run: `sed -n '10150,10360p' src/rustle/pipeline.rs`
Identify: the `bundles` vec, the loaded VG genome handle (`vg_genome_for_discovery`), `config`, and the point just AFTER `vg_families` is finalized (the `let vg_families = vg_families.0;` line, ~10338) where we append the tandem families.

- [ ] **Step 2: Add the pre-pass** (before `discover_family_groups`). Collect the synthetic `FamilyGroup`s — they are used DIRECTLY (not re-discovered), so the quality filter can't drop them.

```rust
        // ── O5/tandem: decompose collapsed tandem arrays before discovery ──
        // (spec 2026-06-02). Opt-in via RUSTLE_VG_TANDEM. A tandem paralog array
        // collapses into ONE bundle (secondaries bridge the copies), so the
        // inter-bundle discovery below never sees a family. Split such bundles
        // into synthetic per-copy sub-bundles and build the family DIRECTLY.
        let tandem_cfg = crate::vg_hmm::tandem::TandemConfig::from_env();
        let mut tandem_families: Vec<crate::vg::FamilyGroup> = Vec::new();
        let mut tandem_sub_bundles: std::collections::HashSet<usize> = std::collections::HashSet::new();
        if tandem_cfg.enabled {
            if let Some(genome) = vg_genome_for_discovery.as_ref() {
                let n_orig = bundles.len();
                let mut replaced: Vec<usize> = Vec::new();
                for bi in 0..n_orig {
                    if let Some(decomp) = crate::vg_hmm::tandem::detect_tandem_bundle(
                        &bundles[bi], genome, tandem_cfg,
                    ) {
                        let parent = bundles[bi].clone();
                        // family_id assigned later (after real discovery) to avoid
                        // colliding with discovered family ids; use a placeholder now.
                        let (fam, sub_idxs) = crate::vg_hmm::tandem::decompose_tandem_to_family(
                            &parent, &decomp, &mut bundles, /*family_id placeholder*/ 0, &config,
                        );
                        for &i in &sub_idxs { tandem_sub_bundles.insert(i); }
                        tandem_families.push(fam);
                        replaced.push(bi);
                        eprintln!("[VG-TANDEM] bundle {} ({}:{}-{}) → {} copy sub-bundles",
                                  bi, parent.chrom, parent.start, parent.end, sub_idxs.len());
                    }
                }
                // Neutralize replaced parents so they don't also assemble.
                for bi in replaced {
                    bundles[bi].reads.clear();
                    bundles[bi].junction_stats = Default::default();
                }
            }
        }
```

- [ ] **Step 3: Append the tandem families AFTER discovery+filter** (just after `let vg_families = vg_families.0;`, ~10338), re-numbering family ids to follow the discovered ones:

```rust
        // Inject directly-built tandem families (bypass the inter-bundle quality
        // filter — they are intra-bundle by construction). Re-number ids to
        // continue past the discovered families.
        let mut vg_families = vg_families;
        let base_fid = vg_families.iter().map(|f| f.family_id + 1).max().unwrap_or(0);
        for (k, mut fam) in tandem_families.into_iter().enumerate() {
            fam.family_id = base_fid + k;
            vg_families.push(fam);
        }
```

Notes (resolve in Step 1):
- `vg_families` is rebound from `vg_families.0` at ~10338; make that binding `mut` (or shadow as above) so the tandem families can be pushed. `vg_family_bundle_set` (built just after) then naturally includes the sub-bundles.
- Clearing the parent's reads is the simplest "replace." If empty bundles are NOT skipped by the assembly loop, also exclude `replaced` indices from `vg_family_bundle_set`/assembly — verify by reading the loop.
- With the tandem families in `vg_families`, the EXISTING `build_family_graph` / `run_fingerprint_em` / completion path (pipeline.rs:~10485+/10700+) processes them unchanged. `vg_genome_for_discovery` must be loaded — guaranteed since the tandem path requires the genome (and the earlier ungate makes it load whenever `--genome-fasta` is given).

- [ ] **Step 4: Tag tandem-copy transcripts** — after VG family assembly, where `vg_copy_id`/`capacity_confidence` are set, also set `tandem_copy = Some(true)` for transcripts whose source bundle index is in `tandem_sub_bundles` (built in Step 2; reuse the existing per-bundle→transcript association used for `vg_copy_id`).

- [ ] **Step 5: Build**

Run: `cargo build --release 2>&1 | grep -E "^error" | head`
Expected: no errors.

- [ ] **Step 6: Smoke-run (default OFF unaffected)**

Run: `target/release/rustle -L /tmp/tspy.bam --vg --genome-fasta /tmp/daz_subset.fa -o /tmp/tandem_off.gtf 2>&1 | grep -c VG-TANDEM`
Expected: `0` (no tandem pass when flag off).

- [ ] **Step 7: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "feat(vg): wire intra-bundle tandem decomposition pre-pass (opt-in)"
```

---

### Task 9: Validation

**Files:** none (validation). Record results in the design doc under `## Validation results`.

- [ ] **Step 1: Regression — de-novo byte-identical (hard gate)**

```bash
target/release/rustle -L bench/<GGO_19 bam> -G bench/GGO_19.gtf -o /tmp/tandem_denovo.gtf
diff <(grep -v '^#' /tmp/tandem_baseline.gtf) <(grep -v '^#' /tmp/tandem_denovo.gtf) | head
```
Expected: empty diff; 95.6/90.5 unchanged.

- [ ] **Step 2: RBMY primary success (flag ON)**

```bash
RUSTLE_VG_TANDEM=1 RUSTLE_VG_ANCHOR_PRIOR=1 \
  target/release/rustle -L /tmp/tspy.bam --vg --genome-fasta /tmp/daz_subset.fa -o /tmp/tandem_rbmy.gtf 2>/tmp/tandem_rbmy.log
grep -E "VG-TANDEM|Discovered|family group" /tmp/tandem_rbmy.log
# c6 region 19,717,578-19,730,926
awk '$3=="transcript" && $4>=19717000 && $4<=19731500' /tmp/tandem_rbmy.gtf
# copies now carry vg_family_id / tandem_copy
awk '$3=="transcript" && $4>=19600000 && $4<=19731500' /tmp/tandem_rbmy.gtf | grep -oE 'vg_family_id "[^"]*"|tandem_copy "[^"]*"' | sort | uniq -c
```
Expected: `[VG-TANDEM]` fires; a family is discovered; **c6 emits a transcript** tagged `tandem_copy "true"` with low `capacity_confidence`; c1–c5 gain `vg_family_id`/`vg_copy_id`.

- [ ] **Step 3: DAZ unchanged (not a tandem array)**

```bash
RUSTLE_VG_TANDEM=1 target/release/rustle -L /tmp/daz.bam --vg --genome-fasta /tmp/daz_subset.fa -o /tmp/tandem_daz.gtf 2>&1 | grep -c VG-TANDEM
```
Expected: `0` (DAZ is dispersed/inverted → `detect_tandem_bundle` returns None); DAZ GTF identical to the non-tandem `--vg` run.

- [ ] **Step 4: False-split guard**

Construct a single large multi-exon gene with an internal repeat (or pick one from GGO_19) and run with the flag; assert it is **not** split into copies.
Expected: no `[VG-TANDEM]` for it (separation/similarity thresholds reject it).

- [ ] **Step 5: No-fabrication**

Assert a copy cluster with reads emits at its own abundance and abstaining reads get `low_confidence`; a region with a single silent copy (0 reads) produces no transcript.

- [ ] **Step 6: Genome-wide spillover guard**

Re-run `bench/paralog_secondary_scan` with `RUSTLE_VG_TANDEM=1`; assert no spurious tandem families on non-array loci and precision flat-or-up.

- [ ] **Step 7: Record + commit**

Append `## Validation results (2026-06-XX)` to the design doc with numbers, then:
```bash
git add docs/superpowers/specs/2026-06-02-intrabundle-tandem-vg-design.md
git commit -m "docs(vg): record intra-bundle tandem validation results"
```

---

## Promotion criteria (flag → default-on-in-VG)

Flip `RUSTLE_VG_TANDEM` default-on only when all hold:
1. Step 1 byte-identical (no de-novo impact) — hard requirement.
2. RBMY: family discovered, c6 recovered with honesty labels, c1–c5 attributed (Step 2).
3. DAZ + false-split + no-fabrication + spillover guards clean (Steps 3–6).
4. On a real multi-copy benchmark, tandem decomposition adds recovered real copies without adding FPs.

Until then it stays opt-in, documented in `docs/vg_genome_scoping.md`.
