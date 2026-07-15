# Generalized projection (`--project-all-families`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Add an opt-in `--project-all-families` path that projects EVERY resolved copy's consensus (not one per family) at id≥0.98 + ≥3-read support, recovering "detected-family, no-projection" missed members, emitted to a distinct `<out>.allproj.tsv` (DNA-localized parCN leg), byte-identical when off.

**Architecture:** New pure module `src/rustle/vg_family/project_all.rs` (consensus/known extraction, dedup, row format) + wiring in `src/bin/gw_family_catalog.rs` reusing `genome_projection::project_families_batch` and `denovo_assemble::reads_in_region`.

**Tech Stack:** Rust; minimap2 (via `project_families_batch`); clap.

## Global Constraints

- **Default OFF = byte-identical:** no `allproj.tsv`, no changed columns in `families.tsv`/`copies.tsv`/`famcn.tsv`. Flag `--project-all-families` (CLI) + env `RUSTLE_PROJECT_ALL_FAMILIES=1`; requires `--homology-primary` (like `--enumerate-copies`).
- **Projection params:** `project_families_batch(consensuses, fasta, &known, 0.98, 0.90, minimap2, threads)`. Consensus key = the family_id (`GWFAM{fi}`) REPEATED per copy, so the batch unions all copies' hits under one family key and applies per-family `known` self-exclusion. `known[fid]` = all of that family's copy spans.
- **Admission gate (exact):** a projected locus is admitted iff (dedup survivor, reciprocal-overlap ≥ 0.50) AND `n_support_reads ≥ 3` (primary reads overlapping it, via `reads_in_region`).
- **Output:** `<out>.allproj.tsv` header `family_id\tchrom\tstart\tend\tidentity\tn_support_reads\toverlaps_existing_copy`. `overlaps_existing_copy` = locus overlaps a copy already in that family. DNA-localized leg — NOT added to copies.tsv/famcn.tsv.
- **Interfaces (verified):** `project_families_batch(consensuses: &[(String,Vec<u8>)], genome_fasta: &str, known: &HashMap<String,Vec<(String,u64,u64)>>, min_identity: f64, min_cov: f64, minimap2: &str, threads: usize) -> anyhow::Result<HashMap<String,Vec<CopyLocus>>>`; `CopyLocus { chrom:String, start:u64, end:u64, identity:f64, cov:f64 }`; `reads_in_region(bam_path:&str, chrom:&str, lo:u64, hi:u64, threads:usize) -> anyhow::Result<(Vec<PrimaryRead>,Vec<BamRead>)>`; `PrimaryRead { chrom, ref_start, ref_end, introns }`. In `gw_family_catalog.rs`, `fams: Vec<Vec<DenovoTranscript>>`, `fid = format!("GWFAM{fi}")`, copy fields `c.seq/c.chrom/c.start/c.end`.
- TDD, focused tests only (`cargo test --lib vg_family::project_all`). The crate is large — never a full `cargo test`.

---

### Task 1: Module scaffold + pure extraction helpers

**Files:**
- Create: `src/rustle/vg_family/project_all.rs`
- Modify: `src/rustle/vg_family/mod.rs` (`pub mod project_all;`)

**Interfaces:**
- Produces: `pub struct CopyIn { pub seq: Vec<u8>, pub chrom: String, pub start: u64, pub end: u64 }`; `pub fn all_copy_consensuses(fams: &[(String, Vec<CopyIn>)]) -> Vec<(String, Vec<u8>)>`; `pub fn known_from_fams(fams: &[(String, Vec<CopyIn>)]) -> std::collections::HashMap<String, Vec<(String, u64, u64)>>`.

- [ ] **Step 1: Write the failing test:**

```rust
#[test]
fn extraction_repeats_fid_and_collects_spans() {
    let fams = vec![
        ("GWFAM0".to_string(), vec![
            CopyIn { seq: b"ACGT".to_vec(), chrom: "chr1".into(), start: 100, end: 200 },
            CopyIn { seq: b"ACGA".to_vec(), chrom: "chr1".into(), start: 500, end: 600 },
        ]),
        ("GWFAM1".to_string(), vec![
            CopyIn { seq: b"TTTT".to_vec(), chrom: "chr2".into(), start: 10, end: 20 },
        ]),
    ];
    let cons = all_copy_consensuses(&fams);
    // one entry per copy, family_id repeated per copy (so project_families_batch unions per family)
    assert_eq!(cons, vec![
        ("GWFAM0".to_string(), b"ACGT".to_vec()),
        ("GWFAM0".to_string(), b"ACGA".to_vec()),
        ("GWFAM1".to_string(), b"TTTT".to_vec()),
    ]);
    let known = known_from_fams(&fams);
    assert_eq!(known["GWFAM0"], vec![("chr1".to_string(),100,200), ("chr1".to_string(),500,600)]);
    assert_eq!(known["GWFAM1"], vec![("chr2".to_string(),10,20)]);
}
```

- [ ] **Step 2: Run, verify RED** — `cargo test --lib vg_family::project_all::tests::extraction_repeats_fid_and_collects_spans` → FAIL (module/types not found).

- [ ] **Step 3: Implement + register:**

```rust
//! Generalized projection (`--project-all-families`): project EVERY resolved copy's consensus (not one per
//! family) to localize members a single family-consensus projection misses. OPTIONAL, additive, DNA-localized
//! parCN leg — never changes the RNA-split catalog or the family definition. Pure extraction + dedup + row
//! formatting here; the minimap2 projection + read-support gate are wired in gw_family_catalog.

use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct CopyIn { pub seq: Vec<u8>, pub chrom: String, pub start: u64, pub end: u64 }

/// One `(family_id, consensus)` entry PER COPY, with the family_id repeated across its copies so
/// `project_families_batch` unions all copies' hits under that one family key.
pub fn all_copy_consensuses(fams: &[(String, Vec<CopyIn>)]) -> Vec<(String, Vec<u8>)> {
    fams.iter().flat_map(|(fid, copies)| copies.iter().map(move |c| (fid.clone(), c.seq.clone()))).collect()
}

/// Per-family copy spans, for the projection's `known` self-exclusion (a copy projecting back onto its own
/// catalogued locus is not a new localization).
pub fn known_from_fams(fams: &[(String, Vec<CopyIn>)]) -> HashMap<String, Vec<(String, u64, u64)>> {
    fams.iter().map(|(fid, copies)| (fid.clone(), copies.iter().map(|c| (c.chrom.clone(), c.start, c.end)).collect())).collect()
}

#[cfg(test)]
mod tests { use super::*; /* tests */ }
```
Register in `mod.rs`: `pub mod project_all; // OPTIONAL --project-all-families recall leg (generalized projection); never alters the RNA-split catalog.`

- [ ] **Step 4: Run, verify GREEN** — PASS.
- [ ] **Step 5: Commit** — `git add src/rustle/vg_family/project_all.rs src/rustle/vg_family/mod.rs && git commit -m "feat(project-all-families): module scaffold + consensus/known extraction"`

---

### Task 2: Dedup + overlap-check + row formatting

**Files:**
- Modify: `src/rustle/vg_family/project_all.rs`

**Interfaces:**
- Consumes: `genome_projection::CopyLocus`.
- Produces: `pub fn dedup_overlapping(loci: Vec<CopyLocus>) -> Vec<CopyLocus>` (reciprocal-overlap ≥ 0.50, keep highest identity); `pub fn overlaps_any(chrom: &str, start: u64, end: u64, spans: &[(String,u64,u64)]) -> bool`; `pub fn format_allproj_row(family_id: &str, l: &CopyLocus, n_support: usize, overlaps_existing: bool) -> String`.

- [ ] **Step 1: Write the failing test:**

```rust
#[test]
fn dedup_overlap_and_row_format() {
    use crate::vg_family::genome_projection::CopyLocus;
    let mk = |s: u64, e: u64, id: f64| CopyLocus { chrom: "chr7".into(), start: s, end: e, identity: id, cov: 0.95 };
    // two overlapping hits (id .994 vs .982) + one disjoint -> 2 survivors, higher id kept
    let out = dedup_overlapping(vec![mk(75976253,75991692,0.994), mk(75976300,75991600,0.982), mk(76360590,76375995,0.990)]);
    assert_eq!(out.len(), 2);
    assert!(out.iter().any(|l| (l.identity - 0.994).abs() < 1e-9 && l.start == 75976253));
    assert!(overlaps_any("chr7", 75976300, 75991600, &[("chr7".into(), 75976253, 75991692)]));
    assert!(!overlaps_any("chr7", 76360590, 76375995, &[("chr7".into(), 75976253, 75991692)]));
    let row = format_allproj_row("GWFAM7", &mk(75976253,75991692,0.994), 41, false);
    assert_eq!(row, "GWFAM7\tchr7\t75976253\t75991692\t0.994\t41\tfalse");
}
```

- [ ] **Step 2: Run, verify RED** — FAIL.

- [ ] **Step 3: Implement:**

```rust
use crate::vg_family::genome_projection::CopyLocus;

fn recip_overlap(a: &CopyLocus, b: &CopyLocus) -> f64 {
    if a.chrom != b.chrom { return 0.0; }
    let (lo, hi) = (a.start.max(b.start), a.end.min(b.end));
    if hi <= lo { return 0.0; }
    let ov = (hi - lo) as f64;
    (ov / (a.end - a.start).max(1) as f64).min(ov / (b.end - b.start).max(1) as f64)
}

/// Collapse reciprocal-overlap ≥ 0.50 loci (from different sibling consensuses hitting one genomic locus)
/// into one, keeping the highest-identity survivor.
pub fn dedup_overlapping(mut loci: Vec<CopyLocus>) -> Vec<CopyLocus> {
    loci.sort_by(|a, b| b.identity.partial_cmp(&a.identity).unwrap_or(std::cmp::Ordering::Equal));
    let mut kept: Vec<CopyLocus> = Vec::new();
    for l in loci {
        if !kept.iter().any(|k| recip_overlap(k, &l) >= 0.50) { kept.push(l); }
    }
    kept
}

pub fn overlaps_any(chrom: &str, start: u64, end: u64, spans: &[(String, u64, u64)]) -> bool {
    spans.iter().any(|(c, s, e)| c == chrom && !(*s > end || *e < start))
}

pub fn format_allproj_row(family_id: &str, l: &CopyLocus, n_support: usize, overlaps_existing: bool) -> String {
    format!("{family_id}\t{}\t{}\t{}\t{:.3}\t{}\t{}", l.chrom, l.start, l.end, l.identity, n_support, overlaps_existing)
}
```

- [ ] **Step 4: Run, verify GREEN** — PASS.
- [ ] **Step 5: Commit** — `git add src/rustle/vg_family/project_all.rs && git commit -m "feat(project-all-families): locus dedup + overlap check + row format"`

---

### Task 3: Wire `--project-all-families` into gw_family_catalog

**Files:**
- Modify: `src/bin/gw_family_catalog.rs`

**Interfaces:**
- Consumes: Task 1/2 (`project_all::*`), `genome_projection::project_families_batch`, `denovo_assemble::reads_in_region`.

- [ ] **Step 1: Add the flag** near the other recall flags in `Args`:

```rust
/// Generalized projection: project EVERY resolved copy consensus (not one per family) at id>=0.98 +
/// >=3-read support to localize members a single-consensus projection misses. Writes <out>.allproj.tsv
/// (DNA-localized parCN, NOT added to copies.tsv). Requires --homology-primary. Default off (byte-identical).
#[arg(long, default_value_t = false)]
project_all_families: bool,
```
And after args parse: `let project_all = (args.project_all_families || std::env::var("RUSTLE_PROJECT_ALL_FAMILIES").ok().as_deref() == Some("1")) && args.homology_primary;`

- [ ] **Step 2: Add the projection block** AFTER the existing `--enumerate-copies` block (so it's an additive sibling; do NOT modify that block). Guard the ENTIRE block on `if project_all {`:

```rust
if project_all {
    use rustle::vg_family::project_all::{CopyIn, all_copy_consensuses, known_from_fams, dedup_overlapping, overlaps_any, format_allproj_row};
    // Build (fid, copies) with the SAME fid the catalog uses.
    let fam_copies: Vec<(String, Vec<CopyIn>)> = fams.iter().enumerate().map(|(fi, copies)| {
        (format!("GWFAM{fi}"), copies.iter().map(|c| CopyIn { seq: c.seq.clone(), chrom: c.chrom.clone(), start: c.start, end: c.end }).collect())
    }).collect();
    let consensuses = all_copy_consensuses(&fam_copies);
    let known = known_from_fams(&fam_copies);
    let proj = rustle::vg_family::genome_projection::project_families_batch(
        &consensuses, &args.fasta, &known, 0.98, 0.90, &refine_params.minimap2, args.threads,
    ).unwrap_or_default();
    let spans_by_fam: std::collections::HashMap<&String, &Vec<(String,u64,u64)>> = known.iter().collect();
    let mut rows: Vec<String> = Vec::new();
    for (fid, locs) in &proj {
        for l in dedup_overlapping(locs.clone()) {
            // read-support gate: >=3 primary reads over the locus
            let n_support = rustle::vg_family::denovo_assemble::reads_in_region(&args.bam, &l.chrom, l.start, l.end, args.threads)
                .map(|(p, _)| p.len()).unwrap_or(0);
            if n_support < 3 { continue; }
            let overlaps = spans_by_fam.get(fid).map(|s| overlaps_any(&l.chrom, l.start, l.end, s)).unwrap_or(false);
            rows.push(format_allproj_row(fid, &l, n_support, overlaps));
        }
    }
    if !rows.is_empty() {
        let mut af = std::fs::File::create(format!("{}.allproj.tsv", args.out))?;
        writeln!(af, "family_id\tchrom\tstart\tend\tidentity\tn_support_reads\toverlaps_existing_copy")?;
        for r in &rows { writeln!(af, "{r}")?; }
        eprintln!("[gw-catalog] project-all-families: {} loci (id>=0.98, >=3 reads) -> {}.allproj.tsv", rows.len(), args.out);
    }
}
```
Confirm the exact names in scope (`fams`, `args.fasta`, `args.bam`, `args.threads`, `refine_params.minimap2`, `refine_params` may be named differently — check the enumerate block, which already uses `refine_params.minimap2`). Adapt if a name differs. Do NOT change the `--enumerate-copies` block or any existing output.

- [ ] **Step 3: Verify** — `cargo build --bin gw_family_catalog` → clean. `cargo test --lib vg_family::project_all` → PASS (unchanged). A quick manual OFF run must not create `allproj.tsv`.

- [ ] **Step 4: Commit** — `git add src/bin/gw_family_catalog.rs && git commit -m "feat(project-all-families): wire the flag + projection + read-support gate into gw_family_catalog"`

---

### Task 4: Validation (byte-identical OFF + Soto A/B) + doc

**Files:**
- Modify: `bench/soto/SOTO_A119B_RECOVERY.md` (update Panel 4 with the measured gain)

- [ ] **Step 1: Byte-identical OFF** — build release; run `gw_family_catalog --homology-primary` (no `--project-all-families`) on a small catalog input; md5 `families.tsv`/`copies.tsv`/`famcn.tsv` vs a pre-feature run → identical, and confirm NO `allproj.tsv` created.
- [ ] **Step 2: Soto A/B** — on `/mnt/linuxdisk/.../soto_reads.bam` (CHM13 fasta), run with `--cross-chrom --homology-primary --enumerate-copies` OFF vs `--project-all-families` ON (HEAVY/genome-wide → harness-tracked background, outputs to /mnt/linuxdisk). Overlap the new `allproj.tsv` loci (where `overlaps_existing_copy=false`) against the 147 missed `80_fams.chr.bed` members (same-chrom interval overlap) → count DISTINCT newly-covered members. Confirm precision: every admitted locus overlaps a real Soto member (0 spurious). Target +10–20 members toward ~72–74 %.
- [ ] **Step 3: Update Panel 4** in `SOTO_A119B_RECOVERY.md` with the measured `--project-all-families` gain (new members, new %, precision), labeled as the DNA-localized leg. Commit.

---

## Self-Review
- **Spec coverage:** flag+env+byte-identical-off (T3/Global) ✓; project all copies via repeated-fid batch keys (T1) ✓; known self-exclusion (T1/T3) ✓; dedup + id≥0.98 + ≥3-read gate (T2/T3) ✓; distinct allproj.tsv, copies.tsv/famcn.tsv untouched (T3) ✓; Soto A/B + byte-identical validation (T4) ✓. Non-goals (no catalog change, no definition change, no re-seeding) respected — the block is guarded and additive.
- **Placeholder scan:** complete code per step; the one adaptation point (`refine_params` name) is flagged with the check.
- **Type consistency:** `CopyIn`, `all_copy_consensuses`/`known_from_fams`, `dedup_overlapping`/`overlaps_any`/`format_allproj_row`, `CopyLocus`, `project_families_batch`, `reads_in_region` used with matching signatures across T1–T3.
