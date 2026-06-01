# Multi-Copy Family Classification (`classify_family`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the operational `classify_family` predicate from the spec — compute the M/H/X/I criteria on a discovered VG family and emit a `FamilyVerdict` (verdict + graded identifiability + locus arrangement) as GTF attributes on `--vg` transcripts, using the raw edit-distance (ΔNM) copy anchor.

**Architecture:** A pure, unit-tested metric core in `vg.rs` (per-read raw-ΔNM anchor + identifiability-partition union-find), then `classify_family(&FamilyGroup, &[Bundle], &Params) -> FamilyVerdict` that reuses the existing `multimap_reads`/`bundle.reads[ri]` placement infrastructure (the `compute_copy_independent_support` scaffold at `vg.rs:1527`). The pipeline computes one verdict per VG family at EM time (where `compute_copy_independent_support` already runs, bundles intact) and threads it to GTF emission. Everything is VG-gated; the default de-novo path produces no `vg_copy_id` transcripts and is byte-identical.

**Tech Stack:** Rust. Reference logic: `bench/paralog_secondary_scan/{nm_anchor.py,run_scan.py,validate_known_families.py}`.

**Conventions:** build `cargo build --release 2>&1 | tail -3` (ends "Finished"); whole-suite tests `cargo test --release --lib` (external test files have a pre-existing `junction_pair_stats` compile break — use `--lib`); commits end with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`; **never** `git add tools/stringtie`.

**Key correction baked in:** the discriminant is **raw ΔNM** (count of distinguishing mismatches), NOT the per-base NM rate used by the existing `compute_copy_independent_support`. Per-base rate mis-scored DAZ3 (NM 14 vs 26 = 12 distinguishing positions = decisive, yet <1% rate gap). This is a NEW classification path; it does **not** modify or remove the existing copy-support guard (leave `compute_copy_independent_support` intact).

---

## File Structure

- `src/rustle/vg.rs` — add: `ReadAnchor` enum + `anchor_read` (pure); `identifiability_partition` (pure, union-find); `FamilyClass`/`Identifiability`/`LocusRel` enums + `FamilyVerdict` struct; `classify_family`. In-module `#[cfg(test)]` tests.
- `src/rustle/pipeline.rs` — call `classify_family` per VG family where `compute_copy_independent_support` is already called (VgSolver::On arm, EM time); store `HashMap<(usize,usize), FamilyVerdict-attrs>` keyed by `(vg_family_id, vg_copy_id)`; apply at GTF-write time.
- `src/rustle/path_extract.rs` (the `Transcript` struct) + `src/rustle/gtf.rs` — carry + emit `family_verdict`, `family_identifiability`, `family_n_copies`, `family_n_expressed` on VG transcripts (mirror the existing `copy_independent_support` carrier/emitter).
- Tests: in-module unit tests (Tasks 1–3); synthetic-fixture + default-de-novo gates (Task 5).

---

### Task 1: Pure raw-ΔNM read anchor (TDD)

**Files:** `src/rustle/vg.rs` (add near `copy_support_fraction`, ~line 1486); in-module `#[cfg(test)]`.

- [ ] **Step 1: Write the failing tests.**
```rust
#[test]
fn anchor_owns_when_dnm_large() {
    // this copy NM=6, sibling NM=504 over comparable extent -> owns
    assert_eq!(anchor_read(6, 4600, &[(504, 4145)], 2, 0.7), ReadAnchor::Owns);
}
#[test]
fn anchor_sibling_when_dnm_negative() {
    // this copy NM=504, sibling NM=6 -> belongs to sibling
    assert_eq!(anchor_read(504, 4145, &[(6, 4600)], 2, 0.7), ReadAnchor::Sibling);
}
#[test]
fn anchor_tie_when_no_distinguishing_position() {
    // NM 2 vs 2 -> tie (no distinguishing position)
    assert_eq!(anchor_read(2, 2900, &[(2, 2900)], 2, 0.7), ReadAnchor::Tie);
    // NM 14 vs 26 -> 12 distinguishing mismatches -> NOT a tie, owns
    assert_eq!(anchor_read(14, 4018, &[(26, 3934)], 2, 0.7), ReadAnchor::Owns);
}
#[test]
fn anchor_owns_when_no_comparable_other() {
    // sibling alignment too short (extent guard) -> not a competitor -> owns
    assert_eq!(anchor_read(50, 4000, &[(2, 200)], 2, 0.7), ReadAnchor::Owns);
    // no other placement at all -> owns (uniquely this copy)
    assert_eq!(anchor_read(50, 4000, &[], 2, 0.7), ReadAnchor::Owns);
}
```
- [ ] **Step 2: Run to verify it fails.** Run: `cargo test --release --lib anchor_ 2>&1 | tail -5` → FAIL (`anchor_read`/`ReadAnchor` undefined).
- [ ] **Step 3: Implement.**
```rust
/// Where a multi-mapping read truly belongs, by raw edit-distance margin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadAnchor { Owns, Sibling, Tie }

/// Classify one read's copy assignment from raw edit distance (the count of
/// distinguishing mismatches — the fingerprint-EM log-likelihood ratio).
/// `nm_c`/`alen_c`: NM and aligned length of the read's placement at THIS copy.
/// `others`: (nm, alen) of the read's placements at OTHER copies.
/// Only placements with `alen >= extent_frac * alen_c` compete (a short partial
/// hit elsewhere is not a real competitor). `dnm = nm_other_best - nm_c`:
/// `Owns` if `dnm >= t`, `Sibling` if `dnm <= -t`, else `Tie`. With no
/// comparable competitor the read uniquely Owns this copy.
pub fn anchor_read(nm_c: u32, alen_c: u64, others: &[(u32, u64)], t: i64, extent_frac: f64) -> ReadAnchor {
    let best_other = others.iter()
        .filter(|&&(_, al)| (al as f64) >= extent_frac * (alen_c as f64))
        .map(|&(nm, _)| nm)
        .min();
    match best_other {
        None => ReadAnchor::Owns,
        Some(bo) => {
            let dnm = bo as i64 - nm_c as i64;
            if dnm >= t { ReadAnchor::Owns }
            else if dnm <= -t { ReadAnchor::Sibling }
            else { ReadAnchor::Tie }
        }
    }
}
```
- [ ] **Step 4: Run to verify it passes.** Run: `cargo test --release --lib anchor_ 2>&1 | tail -5` → all PASS.
- [ ] **Step 5: Commit.**
```bash
git add src/rustle/vg.rs
git commit -m "vg: anchor_read — raw-dNM per-read copy anchor (pure, tested)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 2: Pure identifiability partition (TDD)

**Files:** `src/rustle/vg.rs`; in-module tests.

- [ ] **Step 1: Write the failing tests.**
```rust
#[test]
fn partition_full_when_all_pairs_distinguishable() {
    // 3 copies, every shared pair has a distinguishing read -> 3 classes (full)
    let nonid_pairs: Vec<(usize, usize)> = vec![];
    assert_eq!(identifiability_partition(3, &nonid_pairs), 3);
}
#[test]
fn partition_none_when_all_pairs_nonidentifiable() {
    // 3 identical copies, every pair non-identifiable -> 1 class (none)
    let nonid_pairs = vec![(0, 1), (1, 2), (0, 2)];
    assert_eq!(identifiability_partition(3, &nonid_pairs), 1);
}
#[test]
fn partition_partial_when_some_merge() {
    // copies 0,1 identical; 2 distinct -> 2 classes (partial)
    let nonid_pairs = vec![(0, 1)];
    assert_eq!(identifiability_partition(3, &nonid_pairs), 2);
}
```
- [ ] **Step 2: Run to verify it fails.** Run: `cargo test --release --lib partition_ 2>&1 | tail -5` → FAIL.
- [ ] **Step 3: Implement** (union-find merging non-identifiable pairs; returns the number of identifiability classes).
```rust
/// Number of identifiability classes among `n_copies`: copies are merged when
/// they are NON-identifiable (`nonid_pairs`). classes == n_copies -> full;
/// classes == 1 -> none; otherwise partial.
pub fn identifiability_partition(n_copies: usize, nonid_pairs: &[(usize, usize)]) -> usize {
    let mut parent: Vec<usize> = (0..n_copies).collect();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while parent[r] != r { r = parent[r]; }
        let mut c = x;
        while parent[c] != c { let n = parent[c]; parent[c] = r; c = n; }
        r
    }
    for &(a, b) in nonid_pairs {
        if a < n_copies && b < n_copies {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            if ra != rb { parent[ra] = rb; }
        }
    }
    (0..n_copies).filter(|&i| find(&mut parent, i) == i).count()
}
```
- [ ] **Step 4: Run to verify it passes.** Run: `cargo test --release --lib partition_ 2>&1 | tail -5` → all PASS.
- [ ] **Step 5: Commit.**
```bash
git add src/rustle/vg.rs
git commit -m "vg: identifiability_partition — equivalence classes via non-identifiable pairs

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 3: `FamilyVerdict` types + `classify_family`

**Files:** `src/rustle/vg.rs`; in-module test.

- [ ] **Step 1: Add the result types** (place above `classify_family`).
```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FamilyClass {
    Family,                      // >=2 resolvable expressed copies
    FamilyNonIdentifiable,       // expressed + connected, copies indistinguishable
    GenePlusUnexpressedParalog,  // exactly 1 expressed copy
    Spillover,                   // 0 expressed copies (reads belong to siblings)
    NotConnected,                // copies do not share reads
    NotExpressedHere,            // no reads at the family
    SingleExonOutOfScope,        // all copies single-exon
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Identifiability { Full, Partial, None }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocusRel { Tandem, Distal, Trans, Overlapping, Single }

#[derive(Debug, Clone)]
pub struct FamilyVerdict {
    pub class: FamilyClass,
    pub n_copies: usize,            // M
    pub n_expressed: usize,         // X
    pub connectivity: f64,          // H = fraction of family reads placing at >=2 copies
    pub identifiability: Identifiability,
    pub n_id_classes: usize,
    pub locus_rel: LocusRel,
}

impl FamilyClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            FamilyClass::Family => "family",
            FamilyClass::FamilyNonIdentifiable => "family_nonidentifiable",
            FamilyClass::GenePlusUnexpressedParalog => "gene_plus_unexpressed_paralog",
            FamilyClass::Spillover => "spillover",
            FamilyClass::NotConnected => "not_connected",
            FamilyClass::NotExpressedHere => "not_expressed_here",
            FamilyClass::SingleExonOutOfScope => "single_exon_out_of_scope",
        }
    }
}
impl Identifiability {
    pub fn as_str(&self) -> &'static str {
        match self { Identifiability::Full => "full", Identifiability::Partial => "partial", Identifiability::None => "none" }
    }
}
impl LocusRel {
    pub fn as_str(&self) -> &'static str {
        match self {
            LocusRel::Tandem => "tandem", LocusRel::Distal => "distal",
            LocusRel::Trans => "trans", LocusRel::Overlapping => "overlapping", LocusRel::Single => "single",
        }
    }
}

/// Classification parameters (env-overridable in the pipeline; see Task 4).
pub struct FamilyParams {
    pub dnm_t: i64,         // distinguishing-mismatch threshold (default 2)
    pub extent_frac: f64,   // comparable-extent guard (default 0.7)
    pub k_expr: usize,      // reads to call a copy expressed (default 3)
    pub h_min: f64,         // connectivity floor (default 0.10)
    pub tie_none: f64,      // tie-fraction >= this => identifiability None (default 0.60)
}
impl Default for FamilyParams {
    fn default() -> Self { FamilyParams { dnm_t: 2, extent_frac: 0.7, k_expr: 3, h_min: 0.10, tie_none: 0.60 } }
}
```
- [ ] **Step 2: Write the failing test** (synthetic family from `Bundle`s — mirror how `compute_copy_independent_support` tests build fixtures; if Bundle construction is impractical in a unit test, assert via the Task 5 end-to-end gate instead and note it here). Minimal fixture test:
```rust
#[test]
fn classify_family_smoke() {
    // Build a 2-copy family where copy 0 has 5 reads owning it and copy 1 has 5.
    // (Construct two Bundles with reads carrying nm/exons; link via multimap_reads.)
    let (family, bundles) = crate::vg::test_fixtures::two_copy_resolvable();
    let v = classify_family(&family, &bundles, &FamilyParams::default());
    assert_eq!(v.class, FamilyClass::Family);
    assert_eq!(v.n_copies, 2);
    assert!(v.n_expressed >= 2);
    assert_eq!(v.identifiability, Identifiability::Full);
}
```
  Add `two_copy_resolvable()` to a `#[cfg(test)] pub mod test_fixtures` building two `Bundle`s (chrom "chr1", distinct start/end, each with ≥5 `BundleRead`s with `exons` set so `read_aligned_len > 0`, distinct `read_name_hash`, `nm` low at own copy) and a `FamilyGroup{ bundle_indices: vec![0,1], multimap_reads }` where shared reads have placements at both copies with the owning copy's nm lower by ≥2.
- [ ] **Step 3: Run to verify it fails.** Run: `cargo test --release --lib classify_family_smoke 2>&1 | tail -8` → FAIL.
- [ ] **Step 4: Implement `classify_family`** (reuse the placement-gathering structure of `compute_copy_independent_support`, but raw-ΔNM via `anchor_read`).
```rust
/// Classify a discovered VG family per the M/H/X/I definition (spec §3-§4).
/// Operates on the ORIGINAL FamilyGroup (pre-strand-split) + intact bundles,
/// exactly like `compute_copy_independent_support` (call at EM time).
pub fn classify_family(family: &FamilyGroup, bundles: &[Bundle], p: &FamilyParams) -> FamilyVerdict {
    let n_copies = family.bundle_indices.len();

    // (nm, alen) of a read's placement, or None if unresolved / no aligned length.
    let placement = |fam_pos: usize, ri: usize| -> Option<(u32, u64)> {
        let bi = *family.bundle_indices.get(fam_pos)?;
        let bundle = bundles.get(bi)?;
        let read = bundle.reads.get(ri)?;
        let alen = read_aligned_len(read);
        if alen == 0 { return None; }
        Some((read.nm, alen))
    };

    // scope: a copy is spliced if any read in its bundle has >=1 junction.
    let copy_spliced = |fam_pos: usize| -> bool {
        family.bundle_indices.get(fam_pos)
            .and_then(|&bi| bundles.get(bi))
            .map(|b| b.reads.iter().any(|r| r.exons.len() >= 2 || !r.junctions.is_empty()))
            .unwrap_or(false)
    };
    let any_spliced = (0..n_copies).any(copy_spliced);

    // X: per-copy Owns count, over family multimap reads + unique reads.
    let mut owns: Vec<usize> = vec![0; n_copies];
    // I: pairwise non-identifiable detection.
    let mut pair_has_shared = std::collections::HashSet::<(usize, usize)>::new();
    let mut pair_has_distinguishing = std::collections::HashSet::<(usize, usize)>::new();
    let mut n_reads_total = 0usize;     // distinct family reads seen (for H denom)
    let mut n_reads_shared = 0usize;    // reads placing at >=2 copies (H numer)

    for placements in family.multimap_reads.values() {
        n_reads_total += 1;
        // copies this read places at, with best (nm, alen) per copy
        let mut per_copy: std::collections::HashMap<usize, (u32, u64)> = std::collections::HashMap::new();
        for &(fp, ri) in placements {
            if let Some((nm, al)) = placement(fp, ri) {
                per_copy.entry(fp).and_modify(|e| { if nm < e.0 { *e = (nm, al); } }).or_insert((nm, al));
            }
        }
        if per_copy.len() >= 2 { n_reads_shared += 1; }
        // anchor each copy this read touches vs the others
        let copies: Vec<usize> = per_copy.keys().copied().collect();
        for &c in &copies {
            let (nm_c, al_c) = per_copy[&c];
            let others: Vec<(u32, u64)> = copies.iter().filter(|&&o| o != c).map(|&o| per_copy[&o]).collect();
            if let ReadAnchor::Owns = anchor_read(nm_c, al_c, &others, p.dnm_t, p.extent_frac) {
                owns[c] += 1;
            }
        }
        // pairwise distinguishing: for each pair the read touches, does dNM>=t either way?
        for i in 0..copies.len() {
            for j in (i + 1)..copies.len() {
                let (a, b) = (copies[i].min(copies[j]), copies[i].max(copies[j]));
                pair_has_shared.insert((a, b));
                let dnm = per_copy[&copies[i]].0 as i64 - per_copy[&copies[j]].0 as i64;
                if dnm.abs() >= p.dnm_t { pair_has_distinguishing.insert((a, b)); }
            }
        }
    }
    // unique reads (in a copy's bundle, not a multimap key) always Owns that copy.
    for c in 0..n_copies {
        if let Some(&bi) = family.bundle_indices.get(c) {
            if let Some(bundle) = bundles.get(bi) {
                let n_unique = bundle.reads.iter()
                    .filter(|r| !family.multimap_reads.contains_key(&r.read_name_hash)).count();
                owns[c] += n_unique;
                n_reads_total += n_unique;
            }
        }
    }

    let n_expressed = owns.iter().filter(|&&o| o >= p.k_expr).count();
    let connectivity = if n_reads_total > 0 { n_reads_shared as f64 / n_reads_total as f64 } else { 0.0 };

    // identifiability classes: non-identifiable pair = shared but no distinguishing read.
    let nonid: Vec<(usize, usize)> = pair_has_shared.iter()
        .filter(|pr| !pair_has_distinguishing.contains(pr)).copied().collect();
    let n_id_classes = identifiability_partition(n_copies, &nonid);
    let total_assessed: usize = owns.iter().sum::<usize>().max(1);
    // tie fraction across shared reads (reads that placed at >=2 copies but Owns nobody by dNM)
    let n_tie_reads = n_reads_shared.saturating_sub(
        family.multimap_reads.values().filter(|pl| {
            // a shared read "resolves" if it Owns some copy
            let mut pc: std::collections::HashMap<usize,(u32,u64)> = std::collections::HashMap::new();
            for &(fp,ri) in pl.iter() { if let Some(v)=placement(fp,ri){ pc.entry(fp).and_modify(|e|{if v.0<e.0{*e=v;}}).or_insert(v);} }
            if pc.len() < 2 { return false; }
            let cs: Vec<usize> = pc.keys().copied().collect();
            cs.iter().any(|&c| {
                let (nm,al)=pc[&c]; let o: Vec<(u32,u64)>=cs.iter().filter(|&&x|x!=c).map(|&x|pc[&x]).collect();
                anchor_read(nm,al,&o,p.dnm_t,p.extent_frac)==ReadAnchor::Owns
            })
        }).count()
    );
    let tie_frac = if n_reads_shared > 0 { n_tie_reads as f64 / n_reads_shared as f64 } else { 0.0 };
    let identifiability = if n_id_classes == n_copies && tie_frac < 0.15 { Identifiability::Full }
        else if n_id_classes <= 1 || tie_frac >= p.tie_none { Identifiability::None }
        else { Identifiability::Partial };

    // locus arrangement from the copies' bundle loci.
    let locus_rel = locus_relationship(family, bundles);

    // verdict
    let class = if n_copies < 2 { FamilyClass::NotConnected }
        else if !any_spliced { FamilyClass::SingleExonOutOfScope }
        else if total_assessed < p.k_expr || owns.iter().all(|&o| o == 0) { FamilyClass::NotExpressedHere }
        else if connectivity < p.h_min { FamilyClass::NotConnected }
        else if n_expressed >= 2 { FamilyClass::Family }
        else if matches!(identifiability, Identifiability::None) { FamilyClass::FamilyNonIdentifiable }
        else if n_expressed == 1 { FamilyClass::GenePlusUnexpressedParalog }
        else { FamilyClass::Spillover };

    FamilyVerdict { class, n_copies, n_expressed, connectivity, identifiability, n_id_classes, locus_rel }
}

/// Genomic arrangement of a family's copies from their bundle loci.
fn locus_relationship(family: &FamilyGroup, bundles: &[Bundle]) -> LocusRel {
    let loci: Vec<(&str, u64, u64)> = family.bundle_indices.iter()
        .filter_map(|&bi| bundles.get(bi)).map(|b| (b.chrom.as_str(), b.start, b.end)).collect();
    if loci.len() < 2 { return LocusRel::Single; }
    let same_chrom = loci.iter().all(|l| l.0 == loci[0].0);
    if !same_chrom { return LocusRel::Trans; }
    // overlapping if any two spans overlap; else tandem if all within 1Mb, else distal
    let mut overlap = false; let mut maxgap = 0u64;
    for i in 0..loci.len() { for j in (i+1)..loci.len() {
        let (a, b) = (loci[i], loci[j]);
        if a.1 <= b.2 && b.1 <= a.2 { overlap = true; }
        let gap = a.1.max(b.1).saturating_sub(a.2.min(b.2));
        maxgap = maxgap.max(gap);
    }}
    if overlap { LocusRel::Overlapping } else if maxgap < 1_000_000 { LocusRel::Tandem } else { LocusRel::Distal }
}
```
- [ ] **Step 5: Run to verify it passes.** Run: `cargo test --release --lib classify_family_smoke 2>&1 | tail -8` → PASS. Then `cargo build --release 2>&1 | tail -3` → Finished.
- [ ] **Step 6: Commit.**
```bash
git add src/rustle/vg.rs
git commit -m "vg: classify_family — M/H/X/I FamilyVerdict (raw-dNM anchor, identifiability partition, locus_rel)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 4: Pipeline wiring + GTF emission (VG-gated)

**Files:** `src/rustle/pipeline.rs`, `src/rustle/path_extract.rs` (Transcript struct), `src/rustle/gtf.rs`.

- [ ] **Step 1: Carry fields on `Transcript`.** In `path_extract.rs`, beside the existing `copy_independent_support: Option<f64>`, add:
```rust
    pub family_verdict: Option<String>,
    pub family_identifiability: Option<String>,
    pub family_n_copies: Option<usize>,
    pub family_n_expressed: Option<usize>,
```
  Initialize all four to `None` in every `Transcript { .. }` literal (compiler will list them; mirror how `copy_independent_support: None` was added).
- [ ] **Step 2: Compute verdicts at EM time.** In `pipeline.rs`, in the `VgSolver::On` arm where `compute_copy_independent_support` is invoked over `vg_families` with intact bundles, add alongside it:
```rust
    // Multi-copy family classification (spec 2026-06-01). VG-gated; annotation only.
    let fam_params = crate::rustle::vg::FamilyParams {
        dnm_t: std::env::var("RUSTLE_VG_DNM_T").ok().and_then(|v| v.parse().ok()).unwrap_or(2),
        extent_frac: std::env::var("RUSTLE_VG_EXTENT_FRAC").ok().and_then(|v| v.parse().ok()).unwrap_or(0.7),
        k_expr: std::env::var("RUSTLE_VG_K_EXPR").ok().and_then(|v| v.parse().ok()).unwrap_or(3),
        h_min: 0.10, tie_none: 0.60,
    };
    // family_idx here is the same index used to assign vg_family_id to transcripts.
    let verdict = crate::rustle::vg::classify_family(family, &bundles, &fam_params);
    for (copy_pos, _bi) in family.bundle_indices.iter().enumerate() {
        vg_family_verdict.insert((family_idx, copy_pos), verdict.clone());
    }
```
  Declare `let mut vg_family_verdict: std::collections::HashMap<(usize,usize), crate::rustle::vg::FamilyVerdict> = std::collections::HashMap::new();` next to the existing `vg_copy_support` declaration (~`pipeline.rs:10562`).
- [ ] **Step 3: Apply to transcripts before GTF write.** Where `copy_independent_support` is set on surviving VG transcripts (the late guard ~`pipeline.rs:18063`), set the four fields from `vg_family_verdict.get(&(t.vg_family_id?, t.vg_copy_id?))`:
```rust
    if let (Some(fid), Some(cid)) = (t.vg_family_id, t.vg_copy_id) {
        if let Some(v) = vg_family_verdict.get(&(fid, cid)) {
            t.family_verdict = Some(v.class.as_str().to_string());
            t.family_identifiability = Some(v.identifiability.as_str().to_string());
            t.family_n_copies = Some(v.n_copies);
            t.family_n_expressed = Some(v.n_expressed);
        }
    }
```
  **Do not suppress** any transcript on the verdict — this task is annotation only (the DAZ3-not-a-phantom correction: real copies are kept; acting on the verdict is out of scope here).
- [ ] **Step 4: Emit attributes.** In `gtf.rs`, where `copy_independent_support "<v>"` is emitted for VG transcripts, append (only when present):
```rust
    if let Some(ref s) = t.family_verdict { attrs.push_str(&format!(" family_verdict \"{}\";", s)); }
    if let Some(ref s) = t.family_identifiability { attrs.push_str(&format!(" family_identifiability \"{}\";", s)); }
    if let Some(n) = t.family_n_copies { attrs.push_str(&format!(" family_n_copies \"{}\";", n)); }
    if let Some(n) = t.family_n_expressed { attrs.push_str(&format!(" family_n_expressed \"{}\";", n)); }
```
  (Match the exact attribute-builder idiom used at the `copy_independent_support` emission site — `attrs` may be a different variable; mirror it.)
- [ ] **Step 5: Build + suite.** Run: `cargo build --release 2>&1 | tail -3` → Finished; `cargo test --release --lib 2>&1 | tail -5` → green.
- [ ] **Step 6: Commit.**
```bash
git add src/rustle/pipeline.rs src/rustle/path_extract.rs src/rustle/gtf.rs
git commit -m "vg: emit family_verdict/identifiability/n_copies/n_expressed on --vg transcripts (annotation only)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

### Task 5: Validation gate (controller-evaluated)

- [ ] **DAZ default `--vg`** (`RUSTLE_VG_SUPPORT_TRACE` not needed): run `./target/release/rustle --vg --genome-fasta ../GGO.fasta -L /tmp/daz.bam -o /tmp/daz_cls.gtf`. Confirm DAZ transcripts carry `family_verdict "family"` (or `family_nonidentifiable`), `family_identifiability "partial"` (DAZ is 11% ties), `family_n_copies >= 2`. Record values.
- [ ] **Synthetic fixture**: `rm -f /tmp/synth_assign.attr.tsv; python3 bench/multi_copy_eval/run_oracle.py --fast --check` → still ALL PASS (verdict emission must not perturb assignment/isoform scores).
- [ ] **Default de-novo isolation**: a non-`--vg` GGO_19 run produces output with **no** `family_verdict` attribute and is otherwise unchanged vs the pre-change build (the four fields are `None` without `vg_copy_id`). Confirm with a diff of the GTF (ignoring nothing — it must be byte-identical except no new attrs appear).
- [ ] **Known-family cross-check (analysis, not gated in CI)**: spot-check that the Rust verdict for DAZ matches the Python `validate_known_families.py` row (FAMILY ✓ / partial). Note any divergence (the Rust path uses runtime `multimap_reads`; the Python uses genome-wide placements — small count differences are expected, the *class* must agree).
- [ ] If DAZ comes back `spillover` or `single_exon_out_of_scope`, debug `classify_family` inputs (are `multimap_reads` populated at the call site? are `bundles` intact? do DAZ reads have `exons`/`junctions`?) before proceeding — mirror the `RUSTLE_VG_SUPPORT_TRACE` approach to add a one-line `eprintln!` of `(n_copies, n_expressed, connectivity, tie_frac, class)` per family under a `RUSTLE_VG_CLASSIFY_TRACE` env gate.

---

## Self-Review

- **Spec coverage:** §2 splice-graph object → modeled via `FamilyGroup` copies + `multimap_reads` (the runtime proxy; documented). §3 M/H/X + scope → `classify_family` (n_copies, connectivity, n_expressed/owns, any_spliced). §3 anchor (raw ΔNM, T=2, extent guard) → Task 1 `anchor_read`. §4 identifiability partition (full/partial/none) → Task 2 + Task 3. §5 locus_rel → `locus_relationship`. §9 verdict enum + fields → `FamilyVerdict`. Emission/VG-gating → Task 4. Validation gate (§11) → Task 5. **Certificate/aggregate for non-identifiable classes (§4)** is represented by the `family_nonidentifiable` verdict + `family_identifiability "none"` attribute (the abstention signal); emitting a full per-class abundance polytope is deliberately deferred (YAGNI — annotation first). Note this gap explicitly.
- **Placeholder scan:** all code blocks complete; no TBD. The one judgment call left to the implementer (exact `attrs` variable name at the GTF emission site, the `family_idx` variable in the VG arm) is integration, with the mirror site named.
- **Type consistency:** `anchor_read(u32,u64,&[(u32,u64)],i64,f64)->ReadAnchor` used identically in Tasks 1 and 3; `identifiability_partition(usize,&[(usize,usize)])->usize` consistent in Tasks 2/3; `FamilyVerdict`/`FamilyClass`/`Identifiability`/`LocusRel` field+variant names consistent across Tasks 3/4; GTF attr names match the `Transcript` field names.
