# Expressed-collapsed families (`--collapse-expressed`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** A PSV-free admission path that re-admits dropped exon-identical (K=0, 0-PSV) but heavily-expressed families as an `expressed_collapsed` copy-number class, gated by ≥2 read-supported projection loci (the EEF1A1 guard), behind opt-in `--collapse-expressed`, byte-identical off.

**Architecture:** Extend `src/rustle/vg_family/collapse_enumerate.rs` (pure gate + driver + row format), thread a 3rd return vec through `detect_homology_catalog_genome_wide` (`denovo_pipeline.rs`), and wire the flag + writer in `src/bin/gw_family_catalog.rs`. Reuses `project_family_copies` + `reads_in_region`.

**Tech Stack:** Rust; minimap2 (via `project_family_copies`); clap.

## Global Constraints

- **Default OFF byte-identical:** flag off (and env `RUSTLE_COLLAPSE_EXPRESSED != "1"`) ⇒ no `expressed_collapsed.tsv`, and `families.tsv`/`copies.tsv`/`famcn.tsv`/`collapsed.tsv` unchanged. Requires `--homology-primary`.
- **Gate EXACTLY:** admit a dropped `<2`-distinct-loci candidate iff `project_family_copies(consensus, fasta, &[(chrom,lo,hi)], 0.98, 0.90, ...)` yields ≥2 loci EACH with `n_support ≥ MIN_LOCUS_READS` (=3, primary reads via `reads_in_region(...).map(|(p,_)| p.len())`). NO `hidden_copy`/PSV requirement. `depth_cn` is NOT used.
- **Independent of `--collapse-enumerate`:** both flags default off; a dropped candidate may be admitted by either/both. The `collapse-enumerate` PSV path and its `collapsed.tsv` are UNCHANGED (byte-identical).
- **Copy-NUMBER only:** `expressed_collapsed.tsv` only; never `copies.tsv`/`famcn.tsv`. famCN seed-inclusive (`famcn_from_projection(n) = n+1`, reuse the existing fn).
- **Interfaces (verified):** `project_family_copies(consensus:&[u8], genome_fasta:&str, known:&[(String,u64,u64)], min_identity:f64, min_cov:f64, minimap2:&str, threads:usize) -> anyhow::Result<Vec<CopyLocus>>`; `CopyLocus { chrom, start, end, identity, cov }`; `reads_in_region(bam_path, chrom, lo, hi, threads) -> anyhow::Result<(Vec<PrimaryRead>, Vec<BamRead>)>`; `famcn_from_projection(usize)->usize` (in collapse_enumerate.rs). `detect_homology_catalog_genome_wide` currently returns `Result<(Vec<Vec<DenovoTranscript>>, Vec<CollapsedFamily>)>`; drop-point branch is `else if cfg.collapse_enumerate && loci.len() < 2`. `DenovoConfig` has `collapse_enumerate: bool` + `from_env()`.
- TDD, focused tests only (`cargo test --lib vg_family::collapse_enumerate`).

---

### Task 1: Pure gate + struct + row format

**Files:** Modify `src/rustle/vg_family/collapse_enumerate.rs`

**Interfaces:** Produces `pub const MIN_LOCUS_READS: usize = 3`; `pub fn admit_expressed_collapse(read_supported_loci: usize) -> bool`; `pub struct ExpressedCollapsedFamily { pub chrom: String, pub start: u64, pub end: u64, pub famcn: usize, pub min_locus_reads: usize, pub projection: Vec<CopyLocus> }`; `pub fn format_expressed_collapsed_row(family_id: &str, f: &ExpressedCollapsedFamily) -> String`.

- [ ] **Step 1: Write the failing tests** (in `collapse_enumerate.rs` `mod tests`):

```rust
#[test]
fn admit_expressed_needs_two_read_supported_loci() {
    assert!(!admit_expressed_collapse(0));
    assert!(!admit_expressed_collapse(1));
    assert!(admit_expressed_collapse(2));
    assert!(admit_expressed_collapse(5));
}

#[test]
fn expressed_collapsed_row_format() {
    let f = ExpressedCollapsedFamily {
        chrom: "chr2".into(), start: 97950885, end: 98048181, famcn: 3, min_locus_reads: 32,
        projection: vec![
            CopyLocus { chrom: "chr2".into(), start: 97950885, end: 98048181, identity: 0.998, cov: 0.95 },
            CopyLocus { chrom: "chr2".into(), start: 99100000, end: 99198000, identity: 0.994, cov: 0.93 },
        ],
    };
    assert_eq!(format_expressed_collapsed_row("GWFAMe0", &f),
        "GWFAMe0\tchr2\t97950885\t98048181\t3\t32\tK0_COLLAPSED_EXPRESSED\tchr2:97950885-98048181@0.998;chr2:99100000-99198000@0.994");
}
```

- [ ] **Step 2: Run, verify RED** — `cargo test --lib vg_family::collapse_enumerate::tests::admit_expressed_needs_two_read_supported_loci` → FAIL (not found).

- [ ] **Step 3: Implement:**

```rust
/// Minimum PRIMARY reads at a projected locus for it to count as an EXPRESSED copy (the EEF1A1 guard:
/// silent pseudogenes fall below this).
pub const MIN_LOCUS_READS: usize = 3;

/// Expressed-collapsed admission (PSV-free): a dropped candidate is a real multi-copy family iff it
/// projects to >= 2 genomic loci that are EACH read-supported. No hidden-copy witness required — these
/// families are exon-identical (0 PSVs) but their copies are all transcribed.
pub fn admit_expressed_collapse(read_supported_loci: usize) -> bool {
    read_supported_loci >= 2
}

#[derive(Debug, Clone)]
pub struct ExpressedCollapsedFamily {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub famcn: usize,           // seed-inclusive: read-supported projected loci + the seed
    pub min_locus_reads: usize, // weakest admitted locus's support (transparency)
    pub projection: Vec<CopyLocus>,
}

/// One `<out>.expressed_collapsed.tsv` row: family_id, chrom, start, end, famCN, min_locus_reads,
/// status (`K0_COLLAPSED_EXPRESSED`), projection_loci (`chrom:start-end@identity` joined by `;`).
pub fn format_expressed_collapsed_row(family_id: &str, f: &ExpressedCollapsedFamily) -> String {
    let proj = f.projection.iter()
        .map(|c| format!("{}:{}-{}@{:.3}", c.chrom, c.start, c.end, c.identity))
        .collect::<Vec<_>>().join(";");
    format!("{family_id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        f.chrom, f.start, f.end, f.famcn, f.min_locus_reads, "K0_COLLAPSED_EXPRESSED", proj)
}
```

- [ ] **Step 4: Run, verify GREEN** — PASS.
- [ ] **Step 5: Commit** — `git add src/rustle/vg_family/collapse_enumerate.rs && git commit -m "feat(collapse-expressed): PSV-free admission gate + ExpressedCollapsedFamily + row format"`

---

### Task 2: `readmit_locus_expressed` driver

**Files:** Modify `src/rustle/vg_family/collapse_enumerate.rs`

**Interfaces:** Consumes `project_family_copies`, `reads_in_region`, `famcn_from_projection`, Task 1. Produces `pub fn readmit_locus_expressed(bam_path: &str, chrom: &str, lo: u64, hi: u64, consensus: &[u8], fasta_path: &str, minimap2: &str, threads: usize) -> Option<ExpressedCollapsedFamily>`.

- [ ] **Step 1: Write the failing test** (unit on the read-support-filter → admit decision, using hand-built `CopyLocus` + a stubbed count is not possible since it does I/O; instead test the *pure* decision path by asserting the gate + a construction helper). Add:

```rust
#[test]
fn expressed_driver_builds_family_from_supported_loci() {
    // exercises the post-projection assembly logic: given read-supported loci, famCN is seed-inclusive and
    // min_locus_reads is the weakest. (The minimap2/BAM I/O path is covered by the Soto/EEF1A1 live runs.)
    let loci = vec![
        CopyLocus { chrom: "c".into(), start: 0, end: 9, identity: 0.99, cov: 0.95 },
        CopyLocus { chrom: "c".into(), start: 50, end: 59, identity: 0.995, cov: 0.95 },
    ];
    let supports = vec![32usize, 4usize];
    let fam = build_expressed_family("c", 0, 9, loci.clone(), &supports);
    assert!(fam.is_some());
    let fam = fam.unwrap();
    assert_eq!(fam.famcn, 3);            // 2 supported loci + seed
    assert_eq!(fam.min_locus_reads, 4);  // weakest
    // fewer than 2 supported -> None
    assert!(build_expressed_family("c", 0, 9, loci[..1].to_vec(), &supports[..1]).is_none());
}
```

- [ ] **Step 2: Run, verify RED** — FAIL (fns not found).

- [ ] **Step 3: Implement** a small pure builder (testable) + the I/O driver that calls it:

```rust
use crate::vg_family::denovo_assemble::reads_in_region;

/// Pure: assemble an `ExpressedCollapsedFamily` from already-read-supported projection loci (paired with
/// their support counts). Returns `None` unless `>= 2` loci are supported. Factored out so the admit/famCN
/// logic is unit-testable without minimap2/BAM I/O.
fn build_expressed_family(chrom: &str, lo: u64, hi: u64, loci: Vec<CopyLocus>, supports: &[usize]) -> Option<ExpressedCollapsedFamily> {
    if !admit_expressed_collapse(loci.len()) { return None; }
    let min_locus_reads = supports.iter().copied().min().unwrap_or(0);
    Some(ExpressedCollapsedFamily {
        chrom: chrom.to_string(), start: lo, end: hi,
        famcn: famcn_from_projection(loci.len()), min_locus_reads, projection: loci,
    })
}

/// I/O driver: project the dropped candidate's consensus, keep projection loci with >= MIN_LOCUS_READS
/// primary reads (the EEF1A1 guard), and admit if >= 2 remain. No PSV/hidden-copy requirement.
pub fn readmit_locus_expressed(
    bam_path: &str, chrom: &str, lo: u64, hi: u64, consensus: &[u8], fasta_path: &str, minimap2: &str, threads: usize,
) -> Option<ExpressedCollapsedFamily> {
    let known = vec![(chrom.to_string(), lo, hi)];
    let loci = project_family_copies(consensus, fasta_path, &known, 0.98, 0.90, minimap2, threads).ok()?;
    let mut supported = Vec::new();
    let mut supports = Vec::new();
    for l in loci {
        let n = reads_in_region(bam_path, &l.chrom, l.start, l.end, threads).map(|(p, _)| p.len()).unwrap_or(0);
        if n >= MIN_LOCUS_READS { supports.push(n); supported.push(l); }
    }
    build_expressed_family(chrom, lo, hi, supported, &supports)
}
```
(`reads_in_region` is already imported at the top of the module for `readmit_locus`; reuse that import.)

- [ ] **Step 4: Run, verify GREEN** — `cargo test --lib vg_family::collapse_enumerate` → PASS.
- [ ] **Step 5: Commit** — `git add src/rustle/vg_family/collapse_enumerate.rs && git commit -m "feat(collapse-expressed): readmit_locus_expressed (project + per-locus read-support gate)"`

---

### Task 3: Config flag + 3-tuple wiring + writer

**Files:** Modify `src/rustle/vg_family/denovo_pipeline.rs`, `src/bin/gw_family_catalog.rs`

- [ ] **Step 1: Add the config field.** In `DenovoConfig` add `pub collapse_expressed: bool` (default false next to `collapse_enumerate`); in `from_env()` add `collapse_expressed: std::env::var("RUSTLE_COLLAPSE_EXPRESSED").ok().as_deref() == Some("1")` (matching the `collapse_enumerate` line's functional-update style).

- [ ] **Step 2: 3-tuple return + drop-point.** In `detect_homology_catalog_genome_wide`: change the return type to `Result<(Vec<Vec<DenovoTranscript>>, Vec<CollapsedFamily>, Vec<crate::vg_family::collapse_enumerate::ExpressedCollapsedFamily>)>`; add `let mut expressed: Vec<...ExpressedCollapsedFamily> = Vec::new();`. Restructure the drop branch so BOTH paths can fire on a `loci.len() < 2` collapse:

```rust
} else if (cfg.collapse_enumerate || cfg.collapse_expressed) && loci.len() < 2 {
    let chrom = copies[0].chrom.clone();
    let lo = copies.iter().map(|c| c.start).min().unwrap_or(0);
    let hi = copies.iter().map(|c| c.end).max().unwrap_or(0);
    let consensus = copies.iter().max_by_key(|c| c.seq.len()).map(|c| c.seq.clone()).unwrap_or_default();
    if cfg.collapse_enumerate {
        if let Some(cf) = crate::vg_family::collapse_enumerate::readmit_locus(
            bam_path, &chrom, lo, hi, &consensus, &genome, fasta_path, &refine.minimap2, threads) {
            collapsed.push(cf);
        }
    }
    if cfg.collapse_expressed {
        if let Some(ef) = crate::vg_family::collapse_enumerate::readmit_locus_expressed(
            bam_path, &chrom, lo, hi, &consensus, fasta_path, &refine.minimap2, threads) {
            expressed.push(ef);
        }
    }
}
```
Update the final `Ok((out, collapsed))` to `Ok((out, collapsed, expressed))`. (Check the exact in-scope names — `bam_path`/`fasta_path`/`genome`/`refine.minimap2`/`threads` — from the existing `readmit_locus` call; adapt if any differ.)

- [ ] **Step 3: Update the caller in `gw_family_catalog.rs`.** The `let (raw, collapsed): (...) = if args.homology_primary { detect_homology_catalog_genome_wide(...)? } else { (<expr>, Vec::new()) };` becomes a 3-tuple `let (raw, collapsed, expressed): (Vec<Vec<DenovoTranscript>>, Vec<...CollapsedFamily>, Vec<...ExpressedCollapsedFamily>) = ...` — add `Vec::new()` for the expressed slot in EVERY non-homology branch (there are a few: homology_primary / cross_chrom / plain — update all). Add the flag to `Args`:

```rust
/// Re-admit EXON-IDENTICAL (0-PSV) but heavily-EXPRESSED families that collapse to <2 RNA loci as a
/// K0_COLLAPSED_EXPRESSED copy-number class (projection >=2 loci, each read-supported; needs DNA parCN for
/// per-read resolution). Writes <out>.expressed_collapsed.tsv. Default off; byte-identical when off.
#[arg(long, default_value_t = false)]
collapse_expressed: bool,
```
After parse: `cfg.collapse_expressed = args.collapse_expressed || cfg.collapse_expressed;`. Then, AFTER the existing `collapsed.tsv` writer block, add (guarded):

```rust
if cfg.collapse_expressed && !expressed.is_empty() {
    let mut ef = std::fs::File::create(format!("{}.expressed_collapsed.tsv", args.out))?;
    writeln!(ef, "family_id\tchrom\tstart\tend\tfamCN\tmin_locus_reads\tstatus\tprojection_loci")?;
    for (i, fam) in expressed.iter().enumerate() {
        writeln!(ef, "{}", rustle::vg_family::collapse_enumerate::format_expressed_collapsed_row(&format!("GWFAMe{i}"), fam))?;
    }
    eprintln!("[gw-catalog] collapse-expressed: {} K0_COLLAPSED_EXPRESSED families -> {}.expressed_collapsed.tsv", expressed.len(), args.out);
}
```

- [ ] **Step 4: Verify** — `cargo build --bin gw_family_catalog --bin copy_assign` clean (copy_assign also calls detect_and_assign, not detect_homology — confirm it doesn't destructure detect_homology; if any in-module test destructures the 2-tuple, update to 3). `cargo test --lib vg_family::collapse_enumerate` → PASS. `cargo test --lib vg_family::denovo_pipeline` → PASS (the homology_catalog fixture test destructures the tuple — update it to 3).

- [ ] **Step 5: Commit** — `git add -A src/rustle/vg_family/denovo_pipeline.rs src/bin/gw_family_catalog.rs && git commit -m "feat(collapse-expressed): --collapse-expressed flag + 3-tuple wiring + expressed_collapsed.tsv writer"`

---

### Task 4: Validation (byte-identical OFF + Soto A/B + EEF1A1 control) + doc

**Files:** Modify `bench/soto/SOTO_A119B_RECOVERY.md`

- [ ] **Step 1: Byte-identical OFF** — build release; `copy_assign`/`gw_family_catalog` on GSTM/PCDHB/MAGEA/DAZ with the flag off → md5-identical `families.tsv`/`copies.tsv`/`famcn.tsv`/`collapsed.tsv`; no `expressed_collapsed.tsv`.
- [ ] **Step 2: Soto A/B** — on the Soto BAM (CHM13), `--cross-chrom --homology-primary --enumerate-copies` OFF vs `--collapse-expressed` ON (background). Overlap `expressed_collapsed.tsv` projection loci with the ~29 dead-family members (ANKRD36B/LIMS1/TCAF2/…, `bench/soto/80_fams.chr.bed`); report members recovered as expressed-collapsed copy-number + precision (each family's loci overlap a real Soto region).
- [ ] **Step 3: EEF1A1 control** — `gw_family_catalog --homology-primary --collapse-expressed` on the EEF1A1-scoped BAM (`NC_073229.2` region of GGO_mm.bam) against full GGO.fasta → assert NO `expressed_collapsed.tsv` row (silent pseudogenes fail per-locus read-support). If admitted, raise `MIN_LOCUS_READS` and re-measure.
- [ ] **Step 4: Update `SOTO_A119B_RECOVERY.md`** with the measured dead-family recovery + EEF1A1 control result, labeled as the copy-number leg. Commit.

---

## Self-Review
- **Spec coverage:** PSV-free gate = ≥2 read-supported loci (T1) ✓; per-locus read-support driver, EEF1A1 guard (T2) ✓; flag/env/byte-identical-off, 3-tuple wiring, own file (T3) ✓; Soto A/B + EEF1A1 control + byte-identical validation (T4) ✓; depth_cn not a gate ✓; collapse-enumerate collapsed.tsv untouched ✓.
- **Placeholder scan:** complete code per step; the in-scope-name adaptation points (drop-point locals, caller branches, in-module tuple tests) are flagged.
- **Type consistency:** `admit_expressed_collapse`/`ExpressedCollapsedFamily`/`format_expressed_collapsed_row`/`build_expressed_family`/`readmit_locus_expressed`, `CopyLocus`, `famcn_from_projection`, `project_family_copies`, `reads_in_region` used with matching signatures; the 3-tuple return threaded to every caller branch.
