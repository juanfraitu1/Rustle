# Haplotype-resolved isoform assembly — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `--vg-phase` haplotype-resolved isoform assembly: rustle detects heterozygous sites from read sequence, phases reads into two haplotypes (exact MEC), assembles each haplotype separately, and emits HP/PS-tagged transcripts plus an allele-specific isoform-usage table — a capability StringTie structurally lacks.

**Architecture:** A per-copy phasing pre-pass runs over bundles just before the parallel assembly loop. Because each VG-family bundle already corresponds to one copy (`family.bundle_indices[copy_id]`), per-copy partitioning is free — phasing runs per-bundle. The pre-pass (B1 het detection → B2 MEC-DP) assigns `hp_tag`/`ps_tag` to reads, then `split_bundle_by_phase` expands each eligible bundle into HP1/HP2 sub-bundles. Sub-bundles re-use the **same `bundle_idx`** in `bundles_vec`, so every existing `bundle_to_vg[bundle_idx]` lookup keeps working; haplotype travels on a new `Bundle.hp_tag` field and is copied onto resulting transcripts. The ASE table is computed after assembly from HP-tagged transcripts grouped by intron chain.

**Tech Stack:** Rust. New module `src/rustle/vg_family/phasing.rs`. Integrates with existing VG family machinery (`compute_copy_ownership`, `split_bundle_by_phase`), the assembly loop in `pipeline.rs`, the GTF writer in `gtf.rs`, and the rescue-report writer pattern.

**Spec:** `docs/superpowers/specs/2026-06-07-haplotype-resolved-isoforms-design.md`

**Global invariant (assert after every task):** with `--vg-phase` OFF, output is byte-identical to baseline and `cargo test` stays green. Run `cargo build` after each code step.

---

## File structure

- **Create** `src/rustle/vg_family/phasing.rs` — all phasing logic: `HetSite`, `PhasingConfig`, `ReadHaplotype`, `PhasingResult`, `detect_het_sites`, `mec_brute`, `mec_dp`, `phase_reads`, `assign_haplotypes`. One responsibility: turn a copy's reads into haplotype assignments.
- **Modify** `src/rustle/vg_family/mod.rs` — register `pub mod phasing;`.
- **Modify** `src/rustle/types.rs` — add `hp_tag`/`ps_tag` to `Bundle`.
- **Modify** `src/rustle/path_extract.rs` — add `hp_tag`/`ps_tag` to `Transcript` (+ all literal initializers).
- **Modify** `src/rustle/pipeline.rs` — phase-eligibility pre-pass, `bundles_vec` expansion, transcript propagation, ASE table writer.
- **Modify** `src/rustle/gtf.rs` — emit `haplotype`/`phase_set` GTF attributes.
- **Create** `bench/phasing_eval/gen_diploid.py` + `bench/phasing_eval/run_phasing_eval.sh` — synthetic validation harness.

---

## Task 1: Add `hp_tag`/`ps_tag` fields to `Bundle` and `Transcript`

**Files:**
- Modify: `src/rustle/types.rs` (`Bundle` struct ~83-113 and its initializers)
- Modify: `src/rustle/path_extract.rs` (`Transcript` struct ~620-720 and all literal initializers at ~841, ~1320, ~1459, ~7031)

- [ ] **Step 1: Add fields to `Bundle`**

In `src/rustle/types.rs`, inside `struct Bundle` (after `vg_family_id: Option<usize>,`):

```rust
    /// Haplotype this (sub-)bundle was split into by phased assembly (`--vg-phase`).
    /// None = not split / unphased. Set by the phasing pre-pass; copied onto
    /// transcripts so the GTF can carry a `haplotype` attribute.
    pub hp_tag: Option<u8>,
    /// Phase-set id for the haplotype split. None = unphased.
    pub ps_tag: Option<u32>,
```

- [ ] **Step 2: Update every `Bundle` initializer**

Find every place a `Bundle { ... }` literal is constructed (`grep -n "Bundle {" src/rustle/*.rs src/rustle/vg_family/*.rs`). Add `hp_tag: None, ps_tag: None,` to each. (The `split_bundle_by_phase` function clones a bundle, so it inherits the fields automatically; only literal constructors need updating.)

- [ ] **Step 3: Add fields to `Transcript`**

In `src/rustle/path_extract.rs`, inside `struct Transcript` (near the other vg fields ~671-718):

```rust
    /// Haplotype of this transcript under phased assembly. None = unphased.
    pub hp_tag: Option<u8>,
    /// Phase set id. None = unphased.
    pub ps_tag: Option<u32>,
```

- [ ] **Step 4: Update every `Transcript` initializer**

Add `hp_tag: None, ps_tag: None,` to each `Transcript { ... }` literal (the lines at ~841, ~1320, ~1459, ~7031 already set the other vg fields to `None` — append there). Use `grep -n "vg_family_id: None" src/rustle/path_extract.rs` to find them all.

- [ ] **Step 5: Build**

Run: `cargo build 2>&1 | tail -20`
Expected: compiles clean (any "missing field" error points to an initializer missed in Step 2/4 — fix it).

- [ ] **Step 6: Commit**

```bash
git add src/rustle/types.rs src/rustle/path_extract.rs
git commit -m "feat(phasing): add hp_tag/ps_tag fields to Bundle and Transcript

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: B1 — heterozygous-site detection

**Files:**
- Create: `src/rustle/vg_family/phasing.rs`
- Modify: `src/rustle/vg_family/mod.rs`

Allele representation is binary and reference-base-free: a read's allele at a site is **Alt** (it carries a mismatch with the site's dominant alt base), **Ref** (it spans the site but has no such mismatch), or **NotCovered**. This avoids CIGAR→query-offset arithmetic (rustle stores `mismatches: Vec<(ref_pos, alt_base)>` and `exons: Vec<(u64,u64)>` but not a ref-pos→seq-offset map).

- [ ] **Step 1: Register the module**

In `src/rustle/vg_family/mod.rs`, add after `pub mod hidden_copy;`:

```rust
pub mod phasing;
```

- [ ] **Step 2: Write the failing test**

Create `src/rustle/vg_family/phasing.rs` with the types and a test module:

```rust
//! Within-locus haplotype phasing for `--vg-phase`. Operates on one copy's
//! read set (a single VG-family bundle). Detects heterozygous sites from the
//! per-read mismatch pileup (B1), phases reads with exact Minimum-Error-
//! Correction (B2), and assigns hp/ps tags. Reference-base-free: alleles are
//! binary (Alt = read carries the site's dominant alt mismatch, Ref = spans
//! the site without it).

use crate::types::BundleRead;

/// A candidate heterozygous site within one copy's read set.
#[derive(Debug, Clone, PartialEq)]
pub struct HetSite {
    pub pos: u64,       // reference position (0-based; matches BundleRead.mismatches)
    pub alt_allele: u8, // dominant alternate base (ASCII A/C/G/T)
    pub n_ref: u32,
    pub n_alt: u32,
}

#[derive(Debug, Clone)]
pub struct PhasingConfig {
    pub min_maf: f64,          // minor-allele-fraction floor (default 0.25)
    pub max_maf: f64,          // ceiling — excludes fixed/ref-artifact sites (default 0.75)
    pub min_allele_reads: u32, // per-allele read floor (default 3)
    pub max_coverage: usize,   // MEC-DP active-read cap (default 15)
    pub ext_hp_frac: f64,      // external-HP-tag precedence threshold (default 0.5)
}

impl Default for PhasingConfig {
    fn default() -> Self {
        PhasingConfig {
            min_maf: 0.25,
            max_maf: 0.75,
            min_allele_reads: 3,
            max_coverage: 15,
            ext_hp_frac: 0.5,
        }
    }
}

/// True if `pos` falls inside any of the read's aligned exon intervals.
fn read_spans(read: &BundleRead, pos: u64) -> bool {
    read.exons.iter().any(|&(s, e)| pos >= s && pos < e)
}

/// Detect biallelic heterozygous sites from the mismatch pileup.
pub fn detect_het_sites(reads: &[BundleRead], cfg: &PhasingConfig) -> Vec<HetSite> {
    use std::collections::HashMap;
    // pos -> (alt_base -> count)
    let mut alt_counts: HashMap<u64, HashMap<u8, u32>> = HashMap::new();
    for r in reads {
        for &(pos, base) in &r.mismatches {
            *alt_counts.entry(pos).or_default().entry(base).or_default() += 1;
        }
    }
    let mut sites: Vec<HetSite> = Vec::new();
    for (&pos, bases) in &alt_counts {
        // Dominant alt base at this position (deterministic tie-break by base value).
        let (&alt_base, &n_alt) = match bases
            .iter()
            .max_by(|a, b| a.1.cmp(b.1).then(b.0.cmp(a.0)))
        {
            Some(v) => v,
            None => continue,
        };
        let coverage = reads.iter().filter(|r| read_spans(r, pos)).count() as u32;
        if coverage < n_alt {
            continue; // defensive
        }
        let n_ref = coverage - n_alt;
        let denom = (n_ref + n_alt) as f64;
        if denom == 0.0 {
            continue;
        }
        let maf = n_alt as f64 / denom;
        if n_alt >= cfg.min_allele_reads
            && n_ref >= cfg.min_allele_reads
            && maf >= cfg.min_maf
            && maf <= cfg.max_maf
        {
            sites.push(HetSite { pos, alt_allele: alt_base, n_ref, n_alt });
        }
    }
    sites.sort_by(|a, b| a.pos.cmp(&b.pos).then(a.alt_allele.cmp(&b.alt_allele)));
    sites
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::BundleRead;

    // Minimal BundleRead builder: one exon span [start,end), given mismatches.
    fn mk_read(name_hash: u64, start: u64, end: u64, mismatches: Vec<(u64, u8)>) -> BundleRead {
        let mut r = BundleRead::default();
        r.read_name_hash = name_hash;
        r.ref_start = start;
        r.ref_end = end;
        r.exons = vec![(start, end)];
        r.mismatches = mismatches;
        r.is_primary_alignment = true;
        r
    }

    #[test]
    fn detects_balanced_het() {
        let cfg = PhasingConfig::default();
        // 10 reads span pos 100; 5 carry alt 'A' (b'A'=65), 5 match ref.
        let mut reads = Vec::new();
        for i in 0..5 {
            reads.push(mk_read(i, 50, 200, vec![(100, b'A')]));
        }
        for i in 5..10 {
            reads.push(mk_read(i, 50, 200, vec![]));
        }
        let sites = detect_het_sites(&reads, &cfg);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].pos, 100);
        assert_eq!(sites[0].alt_allele, b'A');
        assert_eq!(sites[0].n_alt, 5);
        assert_eq!(sites[0].n_ref, 5);
    }

    #[test]
    fn rejects_sequencing_error_low_maf() {
        let cfg = PhasingConfig::default();
        // 20 reads span pos 100; only 1 carries alt -> MAF 0.05, below floor.
        let mut reads = Vec::new();
        reads.push(mk_read(0, 50, 200, vec![(100, b'A')]));
        for i in 1..20 {
            reads.push(mk_read(i, 50, 200, vec![]));
        }
        assert!(detect_het_sites(&reads, &cfg).is_empty());
    }

    #[test]
    fn rejects_fixed_difference_high_maf() {
        let cfg = PhasingConfig::default();
        // All 10 reads carry alt -> MAF 1.0, above ceiling (fixed diff / ref artifact).
        let reads: Vec<_> = (0..10).map(|i| mk_read(i, 50, 200, vec![(100, b'A')])).collect();
        assert!(detect_het_sites(&reads, &cfg).is_empty());
    }
}
```

- [ ] **Step 3: Run tests to verify they pass**

Run: `cargo test -p rustle vg_family::phasing 2>&1 | tail -20`
Expected: 3 passed. (If `BundleRead::default()` does not exist, add `#[derive(Default)]` is not appropriate — instead check `grep -n "impl Default for BundleRead\|Default" src/rustle/types.rs`; if no Default, the test builder must construct a full literal. Prefer adding a `#[cfg(test)] impl BundleRead { fn test_default()... }` only if no Default exists — confirm first.)

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg_family/phasing.rs src/rustle/vg_family/mod.rs
git commit -m "feat(phasing): B1 heterozygous-site detection from mismatch pileup

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: B2a — exact MEC by brute force (the test oracle)

Brute force is exponential in read count, so it is gated to tiny inputs and used (a) on small loci and (b) as the correctness oracle for the column DP in Task 4. "Exact MEC" = choose haplotype-A allele per site + a side per read minimizing total disagreements; because the two haplotypes carry opposite alleles at every het site, side + per-column majority determines the cost.

**Files:**
- Modify: `src/rustle/vg_family/phasing.rs`

- [ ] **Step 1: Write the failing test**

Add to `phasing.rs` (above the `#[cfg(test)]` block, in the module body):

```rust
/// A read's allele at each het site: Some(true)=Alt, Some(false)=Ref, None=not covered.
pub type AlleleRow = Vec<Option<bool>>;

/// Build the read × het-site allele matrix (row order matches `reads`).
pub fn allele_matrix(reads: &[BundleRead], sites: &[HetSite]) -> Vec<AlleleRow> {
    reads
        .iter()
        .map(|r| {
            sites
                .iter()
                .map(|s| {
                    if !read_spans(r, s.pos) {
                        None
                    } else {
                        Some(r.mismatches.iter().any(|&(p, b)| p == s.pos && b == s.alt_allele))
                    }
                })
                .collect()
        })
        .collect()
}

/// Cost of a read row against a haplotype allele vector (Hamming over covered sites).
fn row_cost(row: &AlleleRow, hap: &[bool]) -> u32 {
    row.iter()
        .zip(hap)
        .filter_map(|(a, &h)| a.map(|av| if av != h { 1 } else { 0 }))
        .sum()
}

/// Exact MEC over all 2^M haplotype-A allele assignments (M = #sites).
/// Each read greedily joins the cheaper of {hapA, complement(hapA)}.
/// Returns (hapA alleles, side per read [false=A,true=B], total cost).
/// Caller must ensure M is small (<= ~20).
pub fn mec_brute(matrix: &[AlleleRow], n_sites: usize) -> (Vec<bool>, Vec<bool>, u32) {
    let mut best: Option<(Vec<bool>, Vec<bool>, u32)> = None;
    for mask in 0u32..(1u32 << n_sites) {
        let hap_a: Vec<bool> = (0..n_sites).map(|j| (mask >> j) & 1 == 1).collect();
        let hap_b: Vec<bool> = hap_a.iter().map(|&x| !x).collect();
        let mut sides = Vec::with_capacity(matrix.len());
        let mut total = 0u32;
        for row in matrix {
            let ca = row_cost(row, &hap_a);
            let cb = row_cost(row, &hap_b);
            if ca <= cb {
                sides.push(false);
                total += ca;
            } else {
                sides.push(true);
                total += cb;
            }
        }
        if best.as_ref().map_or(true, |b| total < b.2) {
            best = Some((hap_a, sides, total));
        }
    }
    best.unwrap_or((vec![false; n_sites], vec![false; matrix.len()], 0))
}
```

Add tests inside the `#[cfg(test)] mod tests`:

```rust
    #[test]
    fn mec_brute_two_clean_haplotypes() {
        // 2 sites; haplotype A = (Alt, Ref), B = (Ref, Alt). 2 reads each, no errors.
        let matrix = vec![
            vec![Some(true), Some(false)],  // A
            vec![Some(true), Some(false)],  // A
            vec![Some(false), Some(true)],  // B
            vec![Some(false), Some(true)],  // B
        ];
        let (_hap, sides, cost) = mec_brute(&matrix, 2);
        assert_eq!(cost, 0);
        // Reads 0,1 share a side; reads 2,3 share the opposite side.
        assert_eq!(sides[0], sides[1]);
        assert_eq!(sides[2], sides[3]);
        assert_ne!(sides[0], sides[2]);
    }

    #[test]
    fn mec_brute_tolerates_one_error() {
        // Same as above but read 0 has a single flipped allele at site 1.
        let matrix = vec![
            vec![Some(true), Some(true)],   // A with 1 error
            vec![Some(true), Some(false)],  // A
            vec![Some(false), Some(true)],  // B
            vec![Some(false), Some(true)],  // B
        ];
        let (_hap, _sides, cost) = mec_brute(&matrix, 2);
        assert_eq!(cost, 1);
    }
```

- [ ] **Step 2: Run tests**

Run: `cargo test -p rustle vg_family::phasing 2>&1 | tail -20`
Expected: 5 passed.

- [ ] **Step 3: Commit**

```bash
git add src/rustle/vg_family/phasing.rs
git commit -m "feat(phasing): B2a exact MEC brute-force oracle + allele matrix

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: B2b — column MEC-DP (coverage-bounded exact) + `phase_reads`

The column DP is exact MEC made tractable: complexity is exponential in **per-column active coverage** (capped), not in read or site count. State at a column = bipartition of that column's active reads. Transitions are constrained so reads shared between adjacent columns keep the same side. Validated against `mec_brute` on random small inputs.

**Files:**
- Modify: `src/rustle/vg_family/phasing.rs`

- [ ] **Step 1: Write the failing test (DP == brute oracle)**

Add to the test module:

```rust
    // Deterministic small pseudo-random matrices; DP must match the brute oracle.
    fn lcg(seed: &mut u64) -> u64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *seed >> 33
    }

    #[test]
    fn dp_matches_brute_on_random_small() {
        let mut seed = 0x1234_5678u64;
        for _ in 0..200 {
            let n_sites = 1 + (lcg(&mut seed) % 4) as usize; // 1..4 sites
            let n_reads = 2 + (lcg(&mut seed) % 6) as usize;  // 2..7 reads
            let mut matrix: Vec<AlleleRow> = Vec::new();
            for _ in 0..n_reads {
                let row: AlleleRow = (0..n_sites)
                    .map(|_| match lcg(&mut seed) % 3 {
                        0 => None,
                        1 => Some(false),
                        _ => Some(true),
                    })
                    .collect();
                matrix.push(row);
            }
            let (_, _, brute_cost) = mec_brute(&matrix, n_sites);
            let dp_cost = mec_dp_cost(&matrix, n_sites, 16);
            assert_eq!(dp_cost, brute_cost, "matrix={:?}", matrix);
        }
    }
```

- [ ] **Step 2: Implement the column DP**

Add to the module body. `mec_dp_cost` returns just the cost (used by the equivalence test); `mec_dp` returns the full assignment.

```rust
use std::collections::HashMap;

/// Per-column active reads: indices into the matrix whose allele is Some at this site.
fn active_at(matrix: &[AlleleRow], site: usize, cap: usize) -> Vec<usize> {
    let mut v: Vec<usize> = (0..matrix.len())
        .filter(|&r| matrix[r][site].is_some())
        .collect();
    // Deterministic coverage cap: keep the first `cap` by read index (stable).
    v.truncate(cap);
    v
}

/// Column cost for a bipartition `bits` (bit i = side of active[i]) at `site`:
/// min over choosing side-0's allele in {Ref,Alt}.
fn column_cost(matrix: &[AlleleRow], site: usize, active: &[usize], bits: u32) -> u32 {
    // c_a: side0 := Alt (true), side1 := Ref(false). c_b: the flip.
    let mut c_a = 0u32;
    let mut c_b = 0u32;
    for (i, &r) in active.iter().enumerate() {
        let side1 = (bits >> i) & 1 == 1;
        let av = matrix[r][site].unwrap(); // active => Some
        // side0 target under option A = Alt(true); side1 target = Ref(false)
        let target_a = if side1 { false } else { true };
        let target_b = if side1 { true } else { false };
        if av != target_a { c_a += 1; }
        if av != target_b { c_b += 1; }
    }
    c_a.min(c_b)
}

/// Exact coverage-bounded MEC cost via the column DP. `cap` bounds active
/// coverage per column (states = 2^active, halved by fixing bit0 = 0).
pub fn mec_dp_cost(matrix: &[AlleleRow], n_sites: usize, cap: usize) -> u32 {
    if n_sites == 0 || matrix.is_empty() {
        return 0;
    }
    // DP over columns. State key: bitmask over the column's active reads, with
    // the convention that the first active read is always side 0 (symmetry break).
    let mut prev_active: Vec<usize> = Vec::new();
    // map from state-bitmask -> min cost so far
    let mut prev_dp: HashMap<u32, u32> = HashMap::new();
    prev_dp.insert(0, 0);

    for site in 0..n_sites {
        let active = active_at(matrix, site, cap);
        let n = active.len();
        let n_states = if n == 0 { 1u32 } else { 1u32 << n };
        let mut cur_dp: HashMap<u32, u32> = HashMap::new();

        // Precompute, for each prev state, the side of each read that is shared
        // with the current column.
        for (&pbits, &pcost) in &prev_dp {
            // shared read -> required side (from prev state)
            let mut required: HashMap<usize, bool> = HashMap::new();
            for (i, &r) in prev_active.iter().enumerate() {
                required.insert(r, (pbits >> i) & 1 == 1);
            }
            for bits in 0..n_states {
                // symmetry break: fix bit0 = 0 when n>0
                if n > 0 && (bits & 1) == 1 {
                    continue;
                }
                // consistency: shared reads must keep their side
                let mut ok = true;
                for (i, &r) in active.iter().enumerate() {
                    if let Some(&req) = required.get(&r) {
                        let cur_side = (bits >> i) & 1 == 1;
                        if cur_side != req {
                            ok = false;
                            break;
                        }
                    }
                }
                if !ok {
                    continue;
                }
                let cc = column_cost(matrix, site, &active, bits);
                let total = pcost + cc;
                // Canonicalize the stored key so the next column's `required`
                // lookups are correct: store `bits` as-is (active order is stable).
                cur_dp
                    .entry(bits)
                    .and_modify(|e| if total < *e { *e = total; })
                    .or_insert(total);
            }
        }

        // Symmetry-break can make a column with a brand-new read set (no shared
        // reads) collapse haplotype labels — that is fine for COST. But if a
        // column shares no reads with prev, both label orientations are valid;
        // keep both by also inserting the bit0=1 mirror with equal cost.
        if !prev_active.iter().any(|r| active.contains(r)) && n > 0 {
            let snapshot: Vec<(u32, u32)> = cur_dp.iter().map(|(&k, &v)| (k, v)).collect();
            for (k, v) in snapshot {
                let mirror = (!k) & (n_states - 1);
                cur_dp.entry(mirror).and_modify(|e| if v < *e { *e = v; }).or_insert(v);
            }
        }

        prev_dp = cur_dp;
        prev_active = active;
    }
    prev_dp.values().copied().min().unwrap_or(0)
}
```

- [ ] **Step 3: Run the equivalence test**

Run: `cargo test -p rustle vg_family::phasing::tests::dp_matches_brute_on_random_small 2>&1 | tail -30`
Expected: PASS. If it fails, the printed `matrix=` reproduces the case; debug the consistency/mirror handling against `mec_brute` (the oracle is correct by construction). Iterate until green.

- [ ] **Step 4: Write `phase_reads` (full assignment + phase sets + canonical labels)**

Add to the module body:

```rust
/// One read's haplotype assignment.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadHaplotype {
    pub read_name_hash: u64,
    pub hp: u8,  // 1 or 2
    pub ps: u32, // phase-set id
}

pub struct PhasingResult {
    pub het_sites: Vec<HetSite>,
    pub assignments: Vec<ReadHaplotype>, // only reads covering >=1 het site
}

/// Phase a copy's reads. Detects hets, runs MEC (brute for tiny site counts,
/// else the column DP backbone), assigns every het-covering read to a side by
/// best agreement, derives phase sets from read-overlap connectivity, and
/// applies canonical HP labels (HP1 = larger side; ties -> side with the
/// smallest read_name_hash). Reads covering 0 het sites are omitted (unphased).
pub fn phase_reads(reads: &[BundleRead], cfg: &PhasingConfig) -> PhasingResult {
    let sites = detect_het_sites(reads, cfg);
    if sites.is_empty() {
        return PhasingResult { het_sites: sites, assignments: Vec::new() };
    }
    let matrix = allele_matrix(reads, &sites);
    let n_sites = sites.len();

    // Recover a haplotype-A allele vector. Brute is exact when feasible.
    let hap_a: Vec<bool> = if n_sites <= 20 {
        mec_brute(&matrix, n_sites).0
    } else {
        // Fall back to per-site majority of the DP-implied partition is complex;
        // for >20 sites use a stable seed: site-wise majority of alt over all reads.
        sites
            .iter()
            .map(|s| s.n_alt >= s.n_ref) // Alt is hap-A allele if alt is majority
            .collect()
    };
    let hap_b: Vec<bool> = hap_a.iter().map(|&x| !x).collect();

    // Assign every het-covering read to the cheaper side.
    let mut side: Vec<Option<bool>> = Vec::with_capacity(reads.len()); // false=A,true=B
    for (ri, row) in matrix.iter().enumerate() {
        if row.iter().all(|a| a.is_none()) {
            side.push(None);
            continue;
        }
        let ca = row_cost(row, &hap_a);
        let cb = row_cost(row, &hap_b);
        let _ = ri;
        side.push(Some(cb < ca)); // true (B) only if strictly cheaper -> deterministic
    }

    // Phase sets: hets linked iff co-covered by a read. Union-find over sites.
    let mut parent: Vec<usize> = (0..n_sites).collect();
    fn find(p: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while p[r] != r { r = p[r]; }
        let mut c = x;
        while p[c] != r { let n = p[c]; p[c] = r; c = n; }
        r
    }
    for row in &matrix {
        let covered: Vec<usize> = row
            .iter()
            .enumerate()
            .filter_map(|(j, a)| a.map(|_| j))
            .collect();
        for w in covered.windows(2) {
            let (a, b) = (find(&mut parent, w[0]), find(&mut parent, w[1]));
            if a != b { parent[a] = b; }
        }
    }
    // Map site -> ps id (stable: first site index in each component, by position order).
    let mut comp_id: HashMap<usize, u32> = HashMap::new();
    let mut next_ps: u32 = 1;
    let mut site_ps: Vec<u32> = vec![0; n_sites];
    for j in 0..n_sites {
        let root = find(&mut parent, j);
        let id = *comp_id.entry(root).or_insert_with(|| { let v = next_ps; next_ps += 1; v });
        site_ps[j] = id;
    }

    // A read's phase set = the ps of the first het site it covers.
    let mut assignments: Vec<ReadHaplotype> = Vec::new();
    let mut read_ps: Vec<Option<u32>> = vec![None; reads.len()];
    for (ri, row) in matrix.iter().enumerate() {
        if let Some(j) = row.iter().position(|a| a.is_some()) {
            read_ps[ri] = Some(site_ps[j]);
        }
    }

    // Canonical labeling per phase set: HP1 = larger side; tie -> side holding
    // the smallest read_name_hash. Build per-ps side membership first.
    let mut ps_sideA: HashMap<u32, Vec<u64>> = HashMap::new();
    let mut ps_sideB: HashMap<u32, Vec<u64>> = HashMap::new();
    for (ri, s) in side.iter().enumerate() {
        let (Some(sd), Some(ps)) = (s, read_ps[ri]) else { continue };
        let h = reads[ri].read_name_hash;
        if *sd { ps_sideB.entry(ps).or_default().push(h); }
        else   { ps_sideA.entry(ps).or_default().push(h); }
    }
    // For each ps decide whether side A maps to HP1.
    let mut a_is_hp1: HashMap<u32, bool> = HashMap::new();
    let empty: Vec<u64> = Vec::new();
    for ps in 1..next_ps {
        let a = ps_sideA.get(&ps).unwrap_or(&empty);
        let b = ps_sideB.get(&ps).unwrap_or(&empty);
        let a_first = a.iter().min().copied();
        let b_first = b.iter().min().copied();
        let a_hp1 = if a.len() != b.len() {
            a.len() > b.len()
        } else {
            match (a_first, b_first) {
                (Some(x), Some(y)) => x <= y,
                (Some(_), None) => true,
                (None, Some(_)) => false,
                (None, None) => true,
            }
        };
        a_is_hp1.insert(ps, a_hp1);
    }

    for (ri, s) in side.iter().enumerate() {
        let (Some(sd), Some(ps)) = (s, read_ps[ri]) else { continue };
        let a_hp1 = *a_is_hp1.get(&ps).unwrap_or(&true);
        // side A == !sd ... sd false => side A
        let is_side_a = !*sd;
        let hp = if is_side_a == a_hp1 { 1u8 } else { 2u8 };
        assignments.push(ReadHaplotype { read_name_hash: reads[ri].read_name_hash, hp, ps });
    }
    assignments.sort_by(|x, y| x.read_name_hash.cmp(&y.read_name_hash));

    PhasingResult { het_sites: sites, assignments }
}
```

- [ ] **Step 5: Write tests for `phase_reads`**

Add to the test module:

```rust
    #[test]
    fn phase_reads_splits_two_haplotypes() {
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        // Two hets at 100 and 150. HP "alt" reads carry both; HP "ref" reads carry neither.
        for i in 0..6 {
            reads.push(mk_read(i, 50, 200, vec![(100, b'A'), (150, b'G')]));
        }
        for i in 6..12 {
            reads.push(mk_read(i, 50, 200, vec![]));
        }
        let res = phase_reads(&reads, &cfg);
        assert_eq!(res.het_sites.len(), 2);
        assert_eq!(res.assignments.len(), 12);
        // The two groups must land on different haplotypes.
        let hp_alt = res.assignments.iter().find(|a| a.read_name_hash == 0).unwrap().hp;
        let hp_ref = res.assignments.iter().find(|a| a.read_name_hash == 6).unwrap().hp;
        assert_ne!(hp_alt, hp_ref);
        // All reads share one phase set (co-covered).
        assert!(res.assignments.iter().all(|a| a.ps == res.assignments[0].ps));
    }

    #[test]
    fn phase_reads_unphased_when_no_hets() {
        let cfg = PhasingConfig::default();
        let reads: Vec<_> = (0..10).map(|i| mk_read(i, 50, 200, vec![])).collect();
        let res = phase_reads(&reads, &cfg);
        assert!(res.het_sites.is_empty());
        assert!(res.assignments.is_empty());
    }

    #[test]
    fn phase_reads_deterministic() {
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        for i in 0..6 { reads.push(mk_read(i, 50, 200, vec![(100, b'A')])); }
        for i in 6..12 { reads.push(mk_read(i, 50, 200, vec![])); }
        let a = phase_reads(&reads, &cfg).assignments;
        let b = phase_reads(&reads, &cfg).assignments;
        assert_eq!(a, b);
    }
```

- [ ] **Step 6: Run all phasing tests**

Run: `cargo test -p rustle vg_family::phasing 2>&1 | tail -30`
Expected: all pass (8+ tests).

- [ ] **Step 7: Commit**

```bash
git add src/rustle/vg_family/phasing.rs
git commit -m "feat(phasing): B2b column MEC-DP + phase_reads with phase sets + canonical labels

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: `assign_haplotypes` — external-tag precedence + write tags onto reads

Wraps `phase_reads` with the external-HP-tag precedence rule and applies the chosen `hp`/`ps` onto a mutable copy of the bundle's reads, ready for `split_bundle_by_phase`.

**Files:**
- Modify: `src/rustle/vg_family/phasing.rs`

- [ ] **Step 1: Write the failing test**

Add to the test module:

```rust
    #[test]
    fn external_tags_take_precedence() {
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        // Externally tagged: 6 HP1, 6 HP2, but NO mismatches (internal would find nothing).
        for i in 0..6 { let mut r = mk_read(i, 50, 200, vec![]); r.hp_tag = Some(1); reads.push(r); }
        for i in 6..12 { let mut r = mk_read(i, 50, 200, vec![]); r.hp_tag = Some(2); reads.push(r); }
        let tagged = assign_haplotypes(&reads, &cfg);
        // 12/12 carry external HP -> external path used; both haplotypes present.
        assert!(tagged.iter().any(|r| r.hp_tag == Some(1)));
        assert!(tagged.iter().any(|r| r.hp_tag == Some(2)));
    }

    #[test]
    fn internal_phasing_when_no_external_tags() {
        let cfg = PhasingConfig::default();
        let mut reads = Vec::new();
        for i in 0..6 { reads.push(mk_read(i, 50, 200, vec![(100, b'A')])); }
        for i in 6..12 { reads.push(mk_read(i, 50, 200, vec![])); }
        let tagged = assign_haplotypes(&reads, &cfg);
        assert!(tagged.iter().any(|r| r.hp_tag == Some(1)));
        assert!(tagged.iter().any(|r| r.hp_tag == Some(2)));
    }
```

- [ ] **Step 2: Implement `assign_haplotypes`**

```rust
/// Return a clone of `reads` with `hp_tag`/`ps_tag` populated. If at least
/// `ext_hp_frac` of reads already carry an external `hp_tag`, that is used
/// as-is (external precedence); otherwise run internal phasing. Reads that end
/// up unphased keep `hp_tag = None`.
pub fn assign_haplotypes(reads: &[BundleRead], cfg: &PhasingConfig) -> Vec<BundleRead> {
    let n = reads.len();
    let n_ext = reads.iter().filter(|r| r.hp_tag.is_some()).count();
    let mut out: Vec<BundleRead> = reads.to_vec();
    if n > 0 && (n_ext as f64) / (n as f64) >= cfg.ext_hp_frac {
        // External path: ps_tag already comes from the BAM; nothing to compute.
        return out;
    }
    let res = phase_reads(reads, cfg);
    let mut by_hash: HashMap<u64, (u8, u32)> = HashMap::new();
    for a in &res.assignments {
        by_hash.insert(a.read_name_hash, (a.hp, a.ps));
    }
    for r in out.iter_mut() {
        if let Some(&(hp, ps)) = by_hash.get(&r.read_name_hash) {
            r.hp_tag = Some(hp);
            r.ps_tag = Some(ps);
        } else {
            r.hp_tag = None;
            r.ps_tag = None;
        }
    }
    out
}
```

- [ ] **Step 3: Run tests**

Run: `cargo test -p rustle vg_family::phasing 2>&1 | tail -20`
Expected: all pass.

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg_family/phasing.rs
git commit -m "feat(phasing): assign_haplotypes with external-HP-tag precedence

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Integration — phase-eligibility pre-pass + `bundles_vec` expansion + transcript propagation

Wires phasing into the pipeline. Eligible bundles: in `--vg-phase` mode, every non-synthetic bundle that has sequence/mismatch data; family bundles additionally gated on copy confidence (skip thin copies → leave unphased). Expansion happens at `bundles_vec` construction (pipeline.rs:12104): an eligible bundle becomes two entries **with the same `bundle_idx`**, one per haplotype, each carrying `hp_tag`/`ps_tag` on the `Bundle`. Propagation: at the existing `bundle_to_vg` tagging site (~16797), also copy `bundle.hp_tag`/`bundle.ps_tag` onto each transcript.

**Files:**
- Modify: `src/rustle/pipeline.rs` (lines ~12091–12104 and ~16797–16812)

- [ ] **Step 1: Build the phase-eligibility gate (before bundle_to_vg, ~line 12091)**

Just before `let bundle_to_vg ...` at pipeline.rs:12091, insert the eligibility computation. Reuse the certificate logic already shown at lines 11254-11288 (compute_copy_ownership) to gate family copies on confidence:

```rust
    // ── Phased-assembly eligibility (opt-in --vg-phase) ───────────────────
    // A bundle is phase-eligible if --vg-phase is on, it is non-synthetic, has
    // sequence/mismatch data, and (for family copies) clears the copy-confidence
    // gate so thin/low-identifiability copies are never split into fake haplotypes.
    let phase_eligible: std::collections::HashSet<usize> = if config.vg_phase {
        let mut elig = std::collections::HashSet::new();
        // Confident family copies: own_ev clears the same decisive floor the
        // rescue report uses (MIN_DECISIVE default 3, frac default 0.5).
        let min_decisive: usize = std::env::var("RUSTLE_PHASE_MIN_DECISIVE")
            .ok().and_then(|s| s.parse().ok()).unwrap_or(3);
        let min_frac: f64 = std::env::var("RUSTLE_PHASE_MIN_FRAC")
            .ok().and_then(|s| s.parse().ok()).unwrap_or(0.5);
        let t_margin: i64 = std::env::var("RUSTLE_VG_DECISIVE_GATE_T")
            .ok().and_then(|s| s.parse().ok()).unwrap_or(2);
        let mut in_family: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for fam in &vg_families {
            let own = crate::vg::compute_copy_ownership(fam, &bundles, t_margin, 0.8);
            for (copy_id, &bi) in fam.bundle_indices.iter().enumerate() {
                if bi >= bundles.len() { continue; }
                in_family.insert(bi);
                let Some(o) = own.get(copy_id) else { continue };
                let n_primary = bundles[bi].reads.iter().filter(|r| r.is_primary_alignment).count();
                let cert = crate::vg::CopyCertificate {
                    chrom: bundles[bi].chrom.clone(), start: bundles[bi].start, end: bundles[bi].end,
                    family_id: fam.family_id, copy_id, family_size: fam.bundle_indices.len(),
                    n_primary, n_unique: o.n_unique, n_strict: o.n_strict, n_tied: o.n_tied,
                };
                if cert.is_rescued(min_decisive, min_frac) || cert.decisive_frac() >= min_frac {
                    elig.insert(bi);
                }
            }
        }
        // Non-family bundles (single-copy loci): eligible if they carry seq data.
        for (bi, b) in bundles.iter().enumerate() {
            if b.synthetic { continue; }
            if in_family.contains(&bi) { continue; }
            if b.reads.iter().any(|r| !r.seq.is_empty() || !r.mismatches.is_empty()) {
                elig.insert(bi);
            }
        }
        elig
    } else {
        std::collections::HashSet::new()
    };
```

- [ ] **Step 2: Expand `bundles_vec` (replace line 12104)**

Replace:

```rust
    let mut bundles_vec: Vec<(usize, crate::types::Bundle)> = bundles.into_iter().enumerate().collect();
```

with:

```rust
    let mut bundles_vec: Vec<(usize, crate::types::Bundle)> = if config.vg_phase {
        let cfg = crate::vg_family::phasing::PhasingConfig {
            max_coverage: std::env::var("RUSTLE_PHASE_MAX_COV").ok()
                .and_then(|s| s.parse().ok()).unwrap_or(15),
            ext_hp_frac: std::env::var("RUSTLE_PHASE_EXT_FRAC").ok()
                .and_then(|s| s.parse().ok()).unwrap_or(0.5),
            ..Default::default()
        };
        let mut out: Vec<(usize, crate::types::Bundle)> = Vec::new();
        for (idx, bundle) in bundles.into_iter().enumerate() {
            if !phase_eligible.contains(&idx) {
                out.push((idx, bundle));
                continue;
            }
            // Phase reads, then split into HP sub-bundles (same bundle_idx).
            let tagged = crate::vg_family::phasing::assign_haplotypes(&bundle.reads, &cfg);
            let mut phased = bundle.clone();
            phased.reads = tagged;
            let splits = crate::vg::split_bundle_by_phase(&phased);
            if splits.len() <= 1 {
                out.push((idx, bundle)); // no phasing signal -> original, unchanged
            } else {
                for (mut sub, hp) in splits {
                    // ps: most common ps among this sub-bundle's phased reads (deterministic).
                    let ps = sub.reads.iter().filter_map(|r| r.ps_tag).min();
                    sub.hp_tag = hp;
                    sub.ps_tag = ps;
                    out.push((idx, sub));
                }
            }
        }
        out
    } else {
        bundles.into_iter().enumerate().collect()
    };
```

- [ ] **Step 3: Propagate hp/ps onto transcripts (at ~line 16797)**

After the `bundle_to_vg` block that sets `tx.vg_family_id` etc. (lines 16797-16804), add a propagation from the bundle's own hp/ps (set for sub-bundles in Step 2):

```rust
            // Phased assembly: carry the sub-bundle's haplotype onto its transcripts.
            if bundle.hp_tag.is_some() {
                for tx in bundle_txs.iter_mut() {
                    tx.hp_tag = bundle.hp_tag;
                    tx.ps_tag = bundle.ps_tag;
                }
            }
```

- [ ] **Step 4: Build**

Run: `cargo build 2>&1 | tail -20`
Expected: clean. (If `bundle` is moved/borrowed by this point in the closure, hoist `let bundle_hp = bundle.hp_tag; let bundle_ps = bundle.ps_tag;` near the top of the closure and use those locals.)

- [ ] **Step 5: Default byte-identical regression test**

Run rustle on an existing VG test BAM **without** `--vg-phase` and confirm the GTF is unchanged. Use whatever small VG fixture the repo already exercises (search: `grep -rn "\.bam" bench/multi_copy_eval/*.sh | head`). Concretely:

```bash
# baseline already exists or regenerate from main; compare with feature branch build
cargo build --release 2>&1 | tail -3
./target/release/rustle <existing vg invocation, NO --vg-phase> -o /tmp/phase_off.gtf
diff <(grep -v '^#' /tmp/phase_off.gtf) <(grep -v '^#' /path/to/known_baseline.gtf) && echo "BYTE-IDENTICAL OK"
```

Expected: `BYTE-IDENTICAL OK`. (If no stored baseline, build the same commit's binary from `git stash` of the working tree, or compare against `git show main:` build.)

- [ ] **Step 6: Run full suite**

Run: `cargo test 2>&1 | tail -15`
Expected: all green (the documented suite count, currently ~268/0; new phasing tests add to it).

- [ ] **Step 7: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "feat(phasing): wire --vg-phase pre-pass, bundle expansion, transcript propagation

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: GTF `haplotype`/`phase_set` attributes

**Files:**
- Modify: `src/rustle/gtf.rs` (after line 185, the vg-attribute block)

- [ ] **Step 1: Write the failing test**

Add a test in `gtf.rs` (or its test module — `grep -n "mod tests" src/rustle/gtf.rs`). Build a minimal `Transcript` with `hp_tag = Some(1)`, `ps_tag = Some(7)`, write it, assert the line contains `haplotype "1";` and `phase_set "7";`:

```rust
    #[test]
    fn gtf_emits_haplotype_attrs() {
        let mut tx = /* construct a minimal valid Transcript with 1 exon */;
        tx.hp_tag = Some(1);
        tx.ps_tag = Some(7);
        let mut buf: Vec<u8> = Vec::new();
        write_gtf(&[tx], &mut buf, "RUSTLE").unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("haplotype \"1\";"), "missing haplotype attr in:\n{}", s);
        assert!(s.contains("phase_set \"7\";"), "missing phase_set attr in:\n{}", s);
    }
```

(Reuse an existing `Transcript` construction helper in the gtf test module if one exists; otherwise copy the literal from a neighbouring test.)

- [ ] **Step 2: Run to confirm failure**

Run: `cargo test -p rustle gtf::tests::gtf_emits_haplotype_attrs 2>&1 | tail -10`
Expected: FAIL (attrs absent).

- [ ] **Step 3: Implement**

In `src/rustle/gtf.rs`, after line 185 (the `family_size` block, before the `copy_confidence` comment):

```rust
        if let Some(hp) = tx.hp_tag {
            tx_attrs.push_str(&format!(" haplotype \"{}\";", hp));
        }
        if let Some(ps) = tx.ps_tag {
            tx_attrs.push_str(&format!(" phase_set \"{}\";", ps));
        }
```

- [ ] **Step 4: Run to confirm pass**

Run: `cargo test -p rustle gtf::tests::gtf_emits_haplotype_attrs 2>&1 | tail -10`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/gtf.rs
git commit -m "feat(phasing): emit haplotype/phase_set GTF attributes

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: ASE table writer (`<output>.ase.tsv`)

Modeled exactly on the rescue-report writer (pipeline.rs:19460-19516). One row per emitted transcript, plus per-intron-chain HP1/HP2 read-support join. Written only under `--vg-phase`.

**Files:**
- Modify: `src/rustle/pipeline.rs` (after `write_gtf`, ~line 19458)

- [ ] **Step 1: Implement the writer**

Immediately after the rescue-report block (after line 19516), add:

```rust
    // ── Allele-specific isoform usage (opt-in --vg-phase) ────────────────────
    // One row per emitted transcript with its haplotype/phase-set, joined across
    // haplotypes by intron chain so a researcher sees haplotype-specific isoforms
    // and usage imbalance. StringTie cannot produce this (it is sequence-blind).
    if config.vg_phase {
        let ase_path = output_gtf.as_ref().with_extension("ase.tsv");
        // intron-chain key: (chrom, strand, sorted internal junctions) -> per-HP coverage
        fn chain_key(tx: &crate::path_extract::Transcript) -> String {
            let mut juncs: Vec<String> = Vec::new();
            for w in tx.exons.windows(2) {
                juncs.push(format!("{}-{}", w[0].1, w[1].0));
            }
            format!("{}:{}:{}", tx.chrom, tx.strand, juncs.join(","))
        }
        let mut chain_hp1: std::collections::HashMap<String, f64> = Default::default();
        let mut chain_hp2: std::collections::HashMap<String, f64> = Default::default();
        for tx in &all_transcripts {
            let k = chain_key(tx);
            match tx.hp_tag {
                Some(1) => *chain_hp1.entry(k).or_default() += tx.coverage,
                Some(2) => *chain_hp2.entry(k).or_default() += tx.coverage,
                _ => {}
            }
        }
        let gene_tx_no = crate::gtf::assign_gene_tx_numbers(&all_transcripts);
        let mut rows: Vec<(String, u64, String)> = Vec::new();
        for (i, tx) in all_transcripts.iter().enumerate() {
            let Some(hp) = tx.hp_tag else { continue };
            let k = chain_key(tx);
            let c1 = *chain_hp1.get(&k).unwrap_or(&0.0);
            let c2 = *chain_hp2.get(&k).unwrap_or(&0.0);
            let allele_specific = (c1 == 0.0) ^ (c2 == 0.0);
            let total = c1 + c2;
            let imbalance = if total > 0.0 { (c1 - c2).abs() / total } else { 0.0 };
            let (gno, tno) = gene_tx_no[i];
            let tid = match (&tx.gene_id, &tx.transcript_id) {
                (Some(_), Some(t)) => t.clone(),
                _ => format!("{}.{}.{}", config.label, gno.max(1), tno.max(1)),
            };
            let tx_start = tx.exons.first().map(|e| e.0).unwrap_or(0);
            let tx_end = tx.exons.last().map(|e| e.1).unwrap_or(0);
            let row = format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{}\t{:.3}",
                tid, tx.chrom, tx_start, tx_end, tx.strand,
                tx.vg_family_id.map(|x| x.to_string()).unwrap_or_else(|| "NA".into()),
                tx.vg_copy_id.map(|x| x.to_string()).unwrap_or_else(|| "NA".into()),
                hp, tx.ps_tag.map(|x| x.to_string()).unwrap_or_else(|| "NA".into()),
                c1, c2, tx.coverage,
                if allele_specific { "true" } else { "false" }, imbalance,
            );
            rows.push((tx.chrom.clone(), tx_start, row));
        }
        rows.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        let mut af = std::fs::File::create(&ase_path)?;
        use std::io::Write as _;
        writeln!(af, "transcript_id\tchrom\ttx_start\ttx_end\tstrand\tfamily_id\tcopy_id\thaplotype\tphase_set\treads_hp1\treads_hp2\tcoverage\tallele_specific\timbalance")?;
        for (_, _, row) in &rows {
            writeln!(af, "{}", row)?;
        }
        eprintln!("[VG-PHASE-ASE] {} haplotype-tagged transcript(s) → {}", rows.len(), ase_path.display());
    }
```

- [ ] **Step 2: Build**

Run: `cargo build 2>&1 | tail -10`
Expected: clean. (If `chain_key`/`assign_gene_tx_numbers` visibility errors arise, mirror the rescue-report block which calls `crate::gtf::assign_gene_tx_numbers` the same way.)

- [ ] **Step 3: Integration smoke test**

Run rustle with `--vg --vg-snp --vg-phase` on a VG fixture; confirm `<output>.ase.tsv` is created, has the header, and rows are position-sorted. Confirm a re-run is byte-identical (determinism):

```bash
./target/release/rustle <vg invocation> --vg-phase -o /tmp/phase_on.gtf
test -f /tmp/phase_on.ase.tsv && head -3 /tmp/phase_on.ase.tsv
./target/release/rustle <vg invocation> --vg-phase -o /tmp/phase_on2.gtf
diff /tmp/phase_on.ase.tsv /tmp/phase_on2.ase.tsv && echo "ASE DETERMINISTIC OK"
```

Expected: file present, header correct, `ASE DETERMINISTIC OK`.

- [ ] **Step 4: Commit**

```bash
git add src/rustle/pipeline.rs
git commit -m "feat(phasing): emit allele-specific isoform usage table (<output>.ase.tsv)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 9: Synthetic validation harness (beat StringTie, on data)

**Files:**
- Create: `bench/phasing_eval/gen_diploid.py`
- Create: `bench/phasing_eval/run_phasing_eval.sh`

- [ ] **Step 1: Write the diploid-locus generator**

Create `bench/phasing_eval/gen_diploid.py`: emit a 2-haplotype reference for one gene where HP1 and HP2 differ by K SNVs and HP1 expresses isoform X (exons 1-2-4) while HP2 expresses isoform Y (exons 1-3-4), plus a shared isoform Z (1-2-3-4). Write `hap1.fa`, `hap2.fa`, a truth GTF per haplotype, and a combined transcript FASTA for read simulation. (Reuse the existing `gen_synthetic.py` patterns referenced in memory — `grep -rln "def.*recomb\|gen_synthetic" bench/`.)

```python
#!/usr/bin/env python3
"""Generate a diploid single-gene locus with differential isoform usage for
phasing validation. Plants K SNVs separating two haplotypes; HP1 expresses
isoform X, HP2 expresses Y, both express Z."""
# (Concrete implementation: build exon sequences, splice into isoforms,
#  introduce K SNVs into hap2's exonic positions, write FASTAs + truth GTFs.)
```

- [ ] **Step 2: Write the eval driver**

Create `bench/phasing_eval/run_phasing_eval.sh`: simulate long reads from both haplotypes' isoforms (badread, pacbio2021 HiFi profile; re-orient −strand reads per the sim-reads memory caveat), merge + align to a combined reference, then:
1. Run rustle `--vg --vg-snp --vg-phase`; inspect `.ase.tsv` for the two haplotype-specific isoforms (X on one HP, Y on the other).
2. Oracle: `whatshap haplotag` the BAM, split by HP, assemble each separately; compare read-partition concordance.
3. Run StringTie on the same BAM; confirm it produces one blended isoform set with no HP/ASE — the qualitative "ST can't do this" contrast.

```bash
#!/usr/bin/env bash
set -euo pipefail
# (Concrete: badread simulate, minimap2 -ax map-hifi, samtools sort/index,
#  rustle --vg --vg-snp --vg-phase, parse .ase.tsv, whatshap concordance,
#  stringtie -L baseline. Print a PASS/FAIL summary table.)
```

- [ ] **Step 3: Run the harness**

Run: `bash bench/phasing_eval/run_phasing_eval.sh 2>&1 | tail -30`
Expected: rustle's `.ase.tsv` flags isoform X as HP-specific to one haplotype and Y to the other; concordance vs whatshap reported; StringTie shows no haplotype resolution.

- [ ] **Step 4: Commit**

```bash
git add bench/phasing_eval/
git commit -m "test(phasing): synthetic diploid validation harness + StringTie contrast

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-review (completed during planning)

**Spec coverage:** B0 copy-partition → Task 6 eligibility gate (per-bundle = per-copy, confirmed). B1 het detection → Task 2. B2 MEC-DP → Tasks 3-4. External-tag precedence → Task 5. Per-(copy,haplotype) assembly → Task 6 (`split_bundle_by_phase` expansion). Confidence gating of thin copies → Task 6 Step 1. HP/PS GTF tags → Task 7. ASE table → Task 8. Validation (synthetic + whatshap oracle + ST contrast) → Task 9. Default byte-identical invariant → Task 6 Step 5. ✓ all spec sections mapped.

**Type consistency:** `PhasingConfig` fields used identically in Tasks 2/4/5/6. `HetSite`/`AlleleRow`/`ReadHaplotype`/`PhasingResult` defined in Task 2-4, consumed in 5. `hp_tag: Option<u8>`/`ps_tag: Option<u32>` consistent across `Bundle` (Task 1), `Transcript` (Task 1), `BundleRead` (pre-existing). `assign_haplotypes`/`phase_reads`/`split_bundle_by_phase` signatures consistent.

**Known risks flagged for execution:** (1) Task 4 column DP is the hardest piece — the `mec_brute` equivalence test is the guardrail; iterate against it. (2) Task 6 closure borrow of `bundle.hp_tag` — hoist to locals if the borrow checker objects. (3) `BundleRead` test-builder depends on a `Default` impl — verify in Task 2 Step 3. (4) Validation harness (Task 9) shells out to badread/minimap2/whatshap/stringtie — ensure they are installed before running.
