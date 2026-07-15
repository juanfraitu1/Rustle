# K=0-collapsed family enumeration (`--collapse-enumerate`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a default-off `--collapse-enumerate` flag that re-admits near-identical gene families which collapse to <2 RNA-distinct loci — reporting their genome-projected copy *number* + a `K0_COLLAPSED` flag — instead of silently dropping them.

**Architecture:** A pure three-signal admission decision (`hidden_copy` local witness + balanced alt-fraction + ≥2 genome-projection loci) in a new `vg_family/collapse_enumerate.rs` module, wired into `detect_and_assign` at the existing `<min_copies` drop point behind a config flag, emitting a new `<out>.collapsed.tsv`. Reuses `hidden_copy::detect_hidden_copy` and `genome_projection::project_family_copies` — no new detection algorithm.

**Tech Stack:** Rust, existing `vg_family` modules (`hidden_copy`, `genome_projection`, `denovo_pipeline`), minimap2 (already a dependency of projection), noodles BAM.

## Global Constraints

- **Default OFF, byte-identical when OFF.** With the flag off, `families.tsv`/`copies.tsv`/`assignments.tsv` and all existing output are md5-identical to the current release. No new file is written when off.
- **Flag surfaces:** env `RUSTLE_COLLAPSE_ENUMERATE=1` AND CLI `--collapse-enumerate` on both `gw_family_catalog` and `copy_assign`.
- **Admission = ALL THREE:** `ev.flagged` AND `ev.alt_read_fraction >= 0.30` AND `projection_loci >= 2`. Never re-admit on fewer.
- **Re-admitted families are copy-NUMBER only** — written to `<out>.collapsed.tsv` with `status = K0_COLLAPSED`; NEVER added to `copies.tsv`/`assignments.tsv` (no fabricated per-read copies).
- **Existing signatures (do not change):**
  - `hidden_copy::detect_hidden_copy(reads: &[ReadObs], p: &HiddenCopyParams) -> HiddenCopyEvidence`
  - `hidden_copy::ReadObs { start: u64, end: u64, alts: Vec<u64> }`
  - `hidden_copy::HiddenCopyEvidence { n_primary_reads, n_alt_positions, n_alt_reads, alt_read_fraction: f64, flagged: bool }`
  - `genome_projection::project_family_copies(consensus: &[u8], genome_fasta: &str, known: &[(String,u64,u64)], min_identity: f64, min_cov: f64, minimap2: &str, threads: usize) -> anyhow::Result<Vec<CopyLocus>>`
  - `genome_projection::CopyLocus { chrom: String, start: u64, end: u64, identity: f64, cov: f64 }`

---

### Task 1: Config flag + CLI + env wiring

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (`DenovoConfig` struct + its `Default`/`from_env`)
- Modify: `src/bin/gw_family_catalog.rs` (clap args + set `cfg.collapse_enumerate`)
- Modify: `src/bin/copy_assign.rs` (clap args + set `cfg.collapse_enumerate`)

**Interfaces:**
- Produces: `DenovoConfig.collapse_enumerate: bool` (default `false`), read by Task 3.

- [ ] **Step 1: Write the failing test** (in `denovo_pipeline.rs` `#[cfg(test)] mod tests`)

```rust
#[test]
fn collapse_enumerate_defaults_off_and_reads_env() {
    let d = DenovoConfig::default();
    assert!(!d.collapse_enumerate, "must default OFF for byte-identical behavior");
    std::env::set_var("RUSTLE_COLLAPSE_ENUMERATE", "1");
    assert!(DenovoConfig::from_env().collapse_enumerate);
    std::env::remove_var("RUSTLE_COLLAPSE_ENUMERATE");
    assert!(!DenovoConfig::from_env().collapse_enumerate);
}
```

- [ ] **Step 2: Run test, verify it fails**

Run: `cargo test --lib collapse_enumerate_defaults_off_and_reads_env`
Expected: FAIL — `no field collapse_enumerate on DenovoConfig`.

- [ ] **Step 3: Implement** — add the field (next to the existing `collapse_gate: bool` field), set `false` in `Default`, and in `DenovoConfig::from_env()` add:

```rust
// in the struct:
pub collapse_enumerate: bool,
// in Default::default():
collapse_enumerate: false,
// in from_env() (follow the existing env-flag pattern in this fn):
collapse_enumerate: std::env::var("RUSTLE_COLLAPSE_ENUMERATE").ok().as_deref() == Some("1"),
```

- [ ] **Step 4: Run test, verify it passes**

Run: `cargo test --lib collapse_enumerate_defaults_off_and_reads_env`
Expected: PASS.

- [ ] **Step 5: Add the CLI flags.** In `src/bin/gw_family_catalog.rs` and `src/bin/copy_assign.rs`, add a clap arg (follow the existing `--collapse-gate`/`--enumerate-copies` boolean pattern):

```rust
/// Re-admit near-identical families that collapse to <2 RNA loci as K0_COLLAPSED, reporting
/// genome-projected copy NUMBER (needs DNA parCN for per-read resolution). Writes <out>.collapsed.tsv.
/// Default off; when off, all existing output is byte-identical.
#[arg(long, default_value_t = false)]
collapse_enumerate: bool,
```

Then where `cfg` is built, add `cfg.collapse_enumerate = args.collapse_enumerate;` (OR-ed with the env: `cfg.collapse_enumerate = args.collapse_enumerate || cfg.collapse_enumerate;` after `from_env`).

- [ ] **Step 6: Build both bins, commit**

Run: `cargo build --bin gw_family_catalog --bin copy_assign`
Expected: builds clean.
```bash
git add src/rustle/vg_family/denovo_pipeline.rs src/bin/gw_family_catalog.rs src/bin/copy_assign.rs
git commit -m "feat(collapse-enumerate): add --collapse-enumerate flag (default off, env RUSTLE_COLLAPSE_ENUMERATE)"
```

---

### Task 2: Pure admission decision (`collapse_enumerate` module)

**Files:**
- Create: `src/rustle/vg_family/collapse_enumerate.rs`
- Modify: `src/rustle/vg_family/mod.rs` (add `pub mod collapse_enumerate;`)

**Interfaces:**
- Consumes: `hidden_copy::HiddenCopyEvidence`, `genome_projection::CopyLocus`.
- Produces: `pub const MIN_ALT_FRAC: f64 = 0.30;`, `pub fn admit_collapse(ev: &HiddenCopyEvidence, n_projection_loci: usize) -> bool`, `pub struct CollapsedFamily { pub chrom: String, pub start: u64, pub end: u64, pub famcn: usize, pub n_alt_reads: usize, pub alt_read_fraction: f64, pub projection: Vec<CopyLocus> }` — consumed by Task 3 and Task 4.

- [ ] **Step 1: Write the failing tests** (in the new file)

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_family::hidden_copy::HiddenCopyEvidence;
    fn ev(flagged: bool, frac: f64) -> HiddenCopyEvidence {
        HiddenCopyEvidence { n_primary_reads: 300, n_alt_positions: 40, n_alt_reads: (300.0*frac) as usize, alt_read_fraction: frac, flagged }
    }
    #[test]
    fn admits_only_when_all_three_signals_hold() {
        assert!(admit_collapse(&ev(true, 0.50), 2), "flagged + balanced + >=2 loci -> admit");
        assert!(!admit_collapse(&ev(false, 0.50), 2), "not flagged -> reject");
        assert!(!admit_collapse(&ev(true, 0.10), 2), "minor 2nd haplotype (het/edit-like) -> reject");
        assert!(!admit_collapse(&ev(true, 0.50), 1), "single projection locus -> reject");
    }
}
```

- [ ] **Step 2: Run, verify fail**

Run: `cargo test --lib vg_family::collapse_enumerate`
Expected: FAIL — module/function not found.

- [ ] **Step 3: Implement the module**

```rust
//! K=0-collapsed family re-admission (behind `--collapse-enumerate`). A near-identical family that
//! collapses to <2 RNA-distinct loci is re-admitted as copy NUMBER iff it shows a LOCAL collapse:
//! a `hidden_copy` second-haplotype witness that is BALANCED (co-equal depth) AND projects to >=2 genomic loci.
use crate::vg_family::hidden_copy::HiddenCopyEvidence;
use crate::vg_family::genome_projection::CopyLocus;

/// Min balanced-alt fraction: a co-equal collapsed 2nd copy (~0.5), not a minor het/edit (~<=0.1).
pub const MIN_ALT_FRAC: f64 = 0.30;

/// The three-signal gate. ALL must hold (see spec / Global Constraints).
pub fn admit_collapse(ev: &HiddenCopyEvidence, n_projection_loci: usize) -> bool {
    ev.flagged && ev.alt_read_fraction >= MIN_ALT_FRAC && n_projection_loci >= 2
}

#[derive(Debug, Clone)]
pub struct CollapsedFamily {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub famcn: usize,          // genome-projected copy number
    pub n_alt_reads: usize,    // hidden 2nd-haplotype depth
    pub alt_read_fraction: f64,
    pub projection: Vec<CopyLocus>,
}
```

- [ ] **Step 4: Register module + run tests**

Add `pub mod collapse_enumerate;` to `src/rustle/vg_family/mod.rs`.
Run: `cargo test --lib vg_family::collapse_enumerate`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**
```bash
git add src/rustle/vg_family/collapse_enumerate.rs src/rustle/vg_family/mod.rs
git commit -m "feat(collapse-enumerate): pure three-signal admission decision + CollapsedFamily"
```

---

### Task 3: Detect collapsed loci and run the gate (wiring into `detect_and_assign`)

**Files:**
- Modify: `src/rustle/vg_family/collapse_enumerate.rs` (add the locus→decision driver)
- Modify: `src/rustle/vg_family/denovo_pipeline.rs` (`detect_and_assign`: collect dropped `<min_copies` single-locus candidates when `cfg.collapse_enumerate`, run the driver, attach results to the returned struct)

**Interfaces:**
- Consumes: `admit_collapse`, `CollapsedFamily` (Task 2); `hidden_copy::{ReadObs, HiddenCopyParams, detect_hidden_copy}`; `genome_projection::project_family_copies`; `bam::extract_mismatches_vs_fasta` (per-read mismatch positions vs the reference — the same helper the `hidden_copy` O4 scan path uses; check its exact signature in `src/rustle/bam.rs` and adapt).
- Produces: `pub fn readmit_locus(bam_reads_at_locus: &[/* BamRead */], consensus: &[u8], chrom: &str, lo: u64, hi: u64, genome: &GenomeIndex, fasta_path: &str, minimap2: &str, threads: usize) -> Option<CollapsedFamily>`; and a `collapsed: Vec<CollapsedFamily>` field on `detect_and_assign`'s return type (`RegionResult`/the struct returned — inspect the actual return type at `denovo_pipeline.rs` ~1252 and add the field, defaulting empty).

- [ ] **Step 1: Write the failing test** (in `collapse_enumerate.rs`) — the driver decision path with a synthetic balanced-collapse locus vs a minor-het locus, using `detect_hidden_copy` directly on hand-built `ReadObs` (no BAM/minimap2 needed for this unit):

```rust
#[test]
fn readmit_decision_from_readobs_balanced_vs_het() {
    use crate::vg_family::hidden_copy::{ReadObs, HiddenCopyParams, detect_hidden_copy};
    // 20 candidate columns; ~half the reads carry every alt (a co-equal collapsed 2nd copy)
    let cols: Vec<u64> = (0..20).map(|i| 1000 + i*10).collect();
    let mk = |carry: bool| ReadObs { start: 1000, end: 1200, alts: if carry { cols.clone() } else { vec![] } };
    let mut collapse: Vec<ReadObs> = (0..150).map(|_| mk(true)).collect();
    collapse.extend((0..150).map(|_| mk(false)));           // 0.50 balanced
    let ev = detect_hidden_copy(&collapse, &HiddenCopyParams::default());
    assert!(ev.flagged && ev.alt_read_fraction >= MIN_ALT_FRAC);
    assert!(admit_collapse(&ev, 2), "balanced collapse + 2 loci admits");
    // minor het: only 8% carry the alts
    let mut het: Vec<ReadObs> = (0..24).map(|_| mk(true)).collect();
    het.extend((0..276).map(|_| mk(false)));                // 0.08
    let ev2 = detect_hidden_copy(&het, &HiddenCopyParams::default());
    assert!(!admit_collapse(&ev2, 2), "minor het does not admit");
}
```

- [ ] **Step 2: Run, verify fail (or wrong)**

Run: `cargo test --lib readmit_decision_from_readobs_balanced_vs_het`
Expected: FAIL to compile until the test module `use`s resolve / PASS logically once wired. (If `HiddenCopyParams::default()` thresholds differ, adjust the synthetic counts so the balanced case flags and the het case does not — confirm against `HiddenCopyParams::default()` in `hidden_copy.rs`.)

- [ ] **Step 3: Implement `readmit_locus`** — build `ReadObs` from the locus reads (per-read mismatch positions via the `bam` helper), `detect_hidden_copy`, short-circuit on `!flagged`, else `project_family_copies(consensus, fasta_path, &[(chrom,lo,hi)], 0.98, 0.90, minimap2, threads)`, then `admit_collapse`. On admit, return `Some(CollapsedFamily{ chrom, start: lo, end: hi, famcn: loci.len(), n_alt_reads: ev.n_alt_reads, alt_read_fraction: ev.alt_read_fraction, projection: loci })`; else `None`. Pass the locus's own span as `known` so its own genomic locus is excluded from the projected count.

- [ ] **Step 4: Wire into `detect_and_assign`** — where `colocated_families`/refine discards `<min_copies` candidates, and only `if cfg.collapse_enumerate`, collect each discarded single-locus candidate's reads (`reads_in_region` for its `chrom:lo-hi`) + its consensus (the candidate rep's `seq`), call `readmit_locus`, push any `Some` onto the new `collapsed` field. Guard the ENTIRE block behind `cfg.collapse_enumerate` so the off-path is untouched.

- [ ] **Step 5: Run the unit test + a byte-identical smoke**

Run: `cargo test --lib vg_family::collapse_enumerate` → PASS.
Run: `cargo build --bin copy_assign` → clean.

- [ ] **Step 6: Commit**
```bash
git add src/rustle/vg_family/collapse_enumerate.rs src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(collapse-enumerate): detect collapsed loci + run three-signal gate in detect_and_assign"
```

---

### Task 4: Emit `<out>.collapsed.tsv` from both binaries

**Files:**
- Modify: `src/bin/gw_family_catalog.rs` (write collapsed rows after `families.tsv`/`copies.tsv`)
- Modify: `src/bin/copy_assign.rs` (same)

**Interfaces:**
- Consumes: the `collapsed: Vec<CollapsedFamily>` from Task 3's `detect_and_assign` return.

- [ ] **Step 1: Write the failing test** — a small integration test asserting the writer produces the exact header + a row. In `tests/regression/collapse_enumerate.rs` (add a `[[test]]` target to `Cargo.toml`):

```rust
use rustle::vg_family::collapse_enumerate::CollapsedFamily;
use rustle::vg_family::genome_projection::CopyLocus;
#[test]
fn collapsed_tsv_row_format() {
    let f = CollapsedFamily { chrom: "chr2".into(), start: 108994973, end: 109147842, famcn: 2, n_alt_reads: 600, alt_read_fraction: 0.49,
        projection: vec![CopyLocus{chrom:"chr2".into(),start:108994973,end:109147842,identity:0.99,cov:0.95},
                         CopyLocus{chrom:"chr2".into(),start:110869109,end:110895544,identity:0.993,cov:0.92}] };
    let row = rustle::vg_family::collapse_enumerate::format_collapsed_row("GWFAMc0", &f);
    assert_eq!(row, "GWFAMc0\tchr2\t108994973\t109147842\t2\t600\t0.490\tK0_COLLAPSED\tchr2:108994973-109147842@0.990;chr2:110869109-110895544@0.993");
}
```

- [ ] **Step 2: Run, verify fail**

Run: `cargo test --test collapse_enumerate collapsed_tsv_row_format`
Expected: FAIL — `format_collapsed_row` not found.

- [ ] **Step 3: Implement `format_collapsed_row`** in `collapse_enumerate.rs`:

```rust
pub fn format_collapsed_row(family_id: &str, f: &CollapsedFamily) -> String {
    let proj = f.projection.iter()
        .map(|c| format!("{}:{}-{}@{:.3}", c.chrom, c.start, c.end, c.identity))
        .collect::<Vec<_>>().join(";");
    format!("{family_id}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
        f.chrom, f.start, f.end, f.famcn, f.n_alt_reads, f.alt_read_fraction, "K0_COLLAPSED", proj)
}
```

- [ ] **Step 4: Wire the writer** — in each bin, only when `args.collapse_enumerate` AND `!collapsed.is_empty()`, create `format!("{}.collapsed.tsv", args.out)`, write header `family_id\tchrom\tstart\tend\tfamCN\tn_alt_reads\talt_frac\tstatus\tprojection_loci`, then one `format_collapsed_row` per family (ids `GWFAMc0`, `GWFAMc1`, …). Do NOT write the file when the flag is off or the vec is empty.

- [ ] **Step 5: Run test + build**

Run: `cargo test --test collapse_enumerate` → PASS. `cargo build --bins` → clean.

- [ ] **Step 6: Commit**
```bash
git add src/bin/gw_family_catalog.rs src/bin/copy_assign.rs src/rustle/vg_family/collapse_enumerate.rs Cargo.toml tests/regression/collapse_enumerate.rs
git commit -m "feat(collapse-enumerate): emit <out>.collapsed.tsv (famCN + K0_COLLAPSED + projection loci)"
```

---

### Task 5: Validation — byte-identical OFF + Soto ON/OFF measurement

**Files:** none (measurement + a doc)
- Create: `bench/soto/COLLAPSE_ENUMERATE_MEASURE.md`

- [ ] **Step 1: Byte-identical OFF regression.** Build release; run `copy_assign` on GSTM/PCDHB/MAGEA/DAZ (regions in `/home/juanfra/winloci_scratch/sigval_regions.txt`, flags `--homology-primary --skip-poa-diagnostic --min-copies 2`) WITHOUT the flag; md5 `families.tsv`/`famcn_readonly.tsv`/`assignments.tsv` against the pre-feature release outputs (`/home/juanfra/winloci_scratch/rf_*_on.*`).
Expected: **md5-identical** for all four families. (This is the isolation guarantee.)

- [ ] **Step 2: Soto ON vs OFF.** On the Debian disk (`/mnt/linuxdisk/home/juanfraitu/winloci_data`, remount if needed): run `gw_family_catalog --cross-chrom --homology-primary --enumerate-copies` on `soto_reads.bam` (a) without and (b) with `--collapse-enumerate`. Intersect the new `soto_gw*.collapsed.tsv` projection loci with the 15 category-B Soto family regions (`bench/soto/80_fams.chr.bed`, families ID_127/22/226/240/251/338/386/402/458).
Record: how many category-B families are re-admitted (recall gain), and that EVERY re-admitted family's projection loci overlap a real Soto family region (precision = no false re-admissions).

- [ ] **Step 3: EEF1A1 control.** Run `copy_assign --collapse-enumerate` on the EEF1A1 region `NC_073229.2:97600000-97620000` against `GGO_mm.bam`/`GGO.fasta`. Assert **no** `collapsed.tsv` row is produced (EEF1A1 is single-copy; if a row IS produced, the three-signal gate is insufficient → add an intronless/structure filter on the projection loci as a 4th signal, per the spec's honest caveat, and re-measure).

- [ ] **Step 4: Write `bench/soto/COLLAPSE_ENUMERATE_MEASURE.md`** with the three results (byte-identical off ✓, category-B recall gain, precision = 0 false re-admits, EEF1A1 control), and commit.
```bash
git add bench/soto/COLLAPSE_ENUMERATE_MEASURE.md
git commit -m "bench(collapse-enumerate): byte-identical off + Soto on/off measurement + EEF1A1 control"
```

---

## Self-Review

- **Spec coverage:** flag+env (T1) ✓; three-signal gate incl balanced-alt-frac (T2) ✓; hidden_copy + projection wiring at the drop point (T3) ✓; collapsed.tsv + families marker (T4) ✓; byte-identical-off + Soto measurement + EEF1A1 control (T5) ✓. Category-A and per-read K=0 resolution are explicit non-goals — no task, correct.
- **Placeholder scan:** signatures are concrete; the two spots the implementer must confirm against the real code (`bam::extract_mismatches_vs_fasta` exact signature in T3; `HiddenCopyParams::default()` thresholds for the synthetic test counts in T3) are named with the file to check, not left vague.
- **Type consistency:** `CollapsedFamily`, `admit_collapse`, `MIN_ALT_FRAC`, `format_collapsed_row` used consistently across T2/T3/T4; `CopyLocus`/`HiddenCopyEvidence` field names match the read signatures.
