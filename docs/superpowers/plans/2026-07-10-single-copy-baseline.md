# Single-copy baseline & formal transcript — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or
> superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Emit single-copy loci as the copy-number baseline and internalize λ_global into the Rust binary, so
`depth_cn = E_fam/λ_global` no longer depends on an external Python-computed scalar.

**Architecture:** A pure `single_copy` module (filter reps not in any family → lightweight records; median). A
genome-wide baseline function that reuses the existing conflict-catalog traversal via an extracted DRY core. A
flag on `gw_family_catalog` to emit the table + λ. A `--lambda-file` on `copy_assign` to consume it.

**Tech Stack:** Rust (`src/rustle/vg_family/`, `src/bin/`), existing genome-wide traversal + `depth_cn`.

**Spec:** `docs/superpowers/specs/2026-07-10-single-copy-baseline-design.md`

## Global Constraints

- **Not an assembler.** Emit a TABLE (`.single_copy.tsv`), never a GTF. No isoform-level quantification. Do
  NOT touch the `--gtf` path or its k=1 vs k=2 divergence.
- **λ = median n_reads over single-copy loci.** A single-copy locus = a rep whose `tid` appears in **no** emitted
  family. Loci already passed the gate (`n_reads ≥ GATE_MIN_READS`).
- **Memory:** the genome-wide reps are already all held during `colocated_families`; single-copy records are
  **lightweight (no `seq`)**, so the accumulator is small even at ~100k loci.
- **DRY:** the genome-wide single-copy path and the conflict-catalog path share ONE traversal core.
- **λ precedence in `copy_assign`:** explicit `--lambda-global <f64>` > `--lambda-file <path>` > absent (`NaN`).
- Existing genome-wide catalog behaviour must stay byte-identical (the refactor in Task 2 is behavior-preserving).

## File Structure

- `src/rustle/vg_family/single_copy.rs` — CREATE. Pure: `SingleCopyLocus`, `single_copy_loci`, `lambda_global`.
- `src/rustle/vg_family/mod.rs` — MODIFY. Register the module.
- `src/rustle/vg_family/denovo_pipeline.rs` — MODIFY. Extract `gw_reps_and_catalog` core;
  `detect_single_copy_baseline_genome_wide`; transcript-definition doc.
- `src/bin/gw_family_catalog.rs` — MODIFY. `--single-copy-baseline` flag + emit.
- `src/bin/copy_assign.rs` — MODIFY. `--lambda-file` + wire `depth_cn`.

---

### Task 1: Pure `single_copy` module + the formal transcript definition

**Files:**
- Create: `src/rustle/vg_family/single_copy.rs`
- Modify: `src/rustle/vg_family/mod.rs`; `src/rustle/vg_family/denovo_assemble.rs` (doc only)
- Test: in-crate `mod tests` in `single_copy.rs`

**Interfaces:**
- Consumes: `DenovoTranscript { tid: String, chrom: String, start: u64, end: u64, n_reads: u32, strand: char,
  introns: Vec<(u64,u64)>, seq: Vec<u8> }` (`family_detect.rs`); `ColocatedFamily { family_id, chrom, start,
  end, copies: Vec<DenovoTranscript> }` (`denovo_pipeline.rs`).
- Produces:
  - `pub struct SingleCopyLocus { pub chrom: String, pub start: u64, pub end: u64, pub strand: char, pub n_reads: u32, pub n_exons: usize }`
  - `pub fn single_copy_loci(reps: &[DenovoTranscript], families: &[ColocatedFamily]) -> Vec<SingleCopyLocus>`
  - `pub fn lambda_global(loci: &[SingleCopyLocus]) -> Option<f64>`

- [ ] **Step 1: Write the failing tests**

Create `src/rustle/vg_family/single_copy.rs` with only `use` lines + this test module:

```rust
//! Single-copy baseline: the χ(H)=1 loci that calibrate copy number.
//!
//! A TRANSCRIPT is the shipped `DenovoTranscript` object, defined by the `assemble_gate` predicate
//! (`denovo_assemble.rs`): an exact-intron-chain read cluster whose LOCUS (junction-incidence component,
//! `locus_support`) carries >= GATE_MIN_READS reads, whose junctions are all canonical and consistent-strand,
//! and whose spliced length is in [MIN_SPLICED, MAX_SPLICED]. A locus is a copy family of size >= 1:
//! single-copy is the degenerate χ(H)=1 case (0 PSVs); its transcripts are the λ_global baseline that
//! calibrates `depth_cn = E_fam / λ_global`. This module is the copy-number baseline, NOT an isoform catalog.

use super::family_detect::DenovoTranscript;
use super::denovo_pipeline::ColocatedFamily;

#[cfg(test)]
mod tests {
    use super::*;

    fn tx(tid: &str, chrom: &str, start: u64, end: u64, n_reads: u32, strand: char, introns: Vec<(u64,u64)>) -> DenovoTranscript {
        DenovoTranscript { tid: tid.into(), chrom: chrom.into(), start, end, n_reads, strand, introns, seq: vec![] }
    }
    fn fam(id: &str, copies: Vec<DenovoTranscript>) -> ColocatedFamily {
        let chrom = copies[0].chrom.clone();
        ColocatedFamily { family_id: id.into(), chrom, start: 0, end: 0, copies }
    }

    #[test]
    fn single_copy_loci_are_the_reps_no_family_claims() {
        let a = tx("a", "c1", 100, 200, 10, '+', vec![(120,150)]);
        let b = tx("b", "c1", 300, 400, 20, '-', vec![]);
        let c = tx("c", "c1", 500, 600, 30, '+', vec![(520,540),(560,580)]);
        // a family claims a and b; c is single-copy.
        let families = vec![fam("F0", vec![a.clone(), b.clone()])];
        let sc = single_copy_loci(&[a, b, c], &families);
        assert_eq!(sc.len(), 1);
        assert_eq!(sc[0].chrom, "c1");
        assert_eq!(sc[0].start, 500);
        assert_eq!(sc[0].end, 600);
        assert_eq!(sc[0].strand, '+');
        assert_eq!(sc[0].n_reads, 30);
        assert_eq!(sc[0].n_exons, 3, "n_exons = introns.len() + 1");
    }

    #[test]
    fn single_copy_loci_membership_is_by_tid() {
        let a = tx("a", "c1", 100, 200, 10, '+', vec![]);
        // family copy shares the tid but is a distinct clone; a must still be recognised as claimed.
        let claimed = tx("a", "c1", 100, 200, 10, '+', vec![]);
        let other = tx("z", "c1", 100, 200, 10, '+', vec![]);
        let families = vec![fam("F0", vec![claimed, other])];
        assert!(single_copy_loci(&[a], &families).is_empty(), "a's tid is claimed -> not single-copy");
    }

    #[test]
    fn single_copy_locus_carries_no_seq() {
        // SingleCopyLocus has no seq field at all -- this is a compile-level guarantee, asserted by shape.
        let a = tx("a", "c1", 1, 2, 5, '+', vec![]);
        let sc = single_copy_loci(&[a], &[]);
        let _ = SingleCopyLocus { chrom: sc[0].chrom.clone(), start: sc[0].start, end: sc[0].end,
                                  strand: sc[0].strand, n_reads: sc[0].n_reads, n_exons: sc[0].n_exons };
    }

    #[test]
    fn lambda_global_is_the_median_n_reads() {
        let mk = |n: u32| SingleCopyLocus { chrom: "c".into(), start: 0, end: 1, strand: '+', n_reads: n, n_exons: 1 };
        assert_eq!(lambda_global(&[mk(10), mk(20), mk(30), mk(40)]), Some(25.0), "even -> mean of middle two");
        assert_eq!(lambda_global(&[mk(30), mk(10), mk(20)]), Some(20.0), "odd, unsorted -> middle after sort");
        assert_eq!(lambda_global(&[]), None, "empty -> None");
        assert_eq!(lambda_global(&[mk(7)]), Some(7.0));
    }
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --lib single_copy 2>&1 | tail -5`
Expected: FAIL — `cannot find type SingleCopyLocus` / `function single_copy_loci`.

- [ ] **Step 3: Implement**

Add above the test module in `single_copy.rs`:

```rust
/// One single-copy locus (χ(H)=1). Carries NO `seq` — this is a lightweight baseline record, so a genome-wide
/// accumulator of ~100k of these is cheap.
#[derive(Clone, Debug, PartialEq)]
pub struct SingleCopyLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub n_reads: u32,
    pub n_exons: usize,
}

/// Reps that NO family claims, as lightweight `SingleCopyLocus` records. Membership is by `tid`.
pub fn single_copy_loci(reps: &[DenovoTranscript], families: &[ColocatedFamily]) -> Vec<SingleCopyLocus> {
    let claimed: std::collections::HashSet<&str> =
        families.iter().flat_map(|f| f.copies.iter().map(|c| c.tid.as_str())).collect();
    reps.iter()
        .filter(|r| !claimed.contains(r.tid.as_str()))
        .map(|r| SingleCopyLocus {
            chrom: r.chrom.clone(),
            start: r.start,
            end: r.end,
            strand: r.strand,
            n_reads: r.n_reads,
            n_exons: r.introns.len() + 1,
        })
        .collect()
}

/// Median read count over the single-copy loci — the genome-wide single-copy expression floor `λ_global`.
/// `None` on an empty set.
pub fn lambda_global(loci: &[SingleCopyLocus]) -> Option<f64> {
    if loci.is_empty() {
        return None;
    }
    let mut v: Vec<u32> = loci.iter().map(|l| l.n_reads).collect();
    v.sort_unstable();
    let n = v.len();
    Some(if n % 2 == 1 {
        v[n / 2] as f64
    } else {
        (v[n / 2 - 1] as f64 + v[n / 2] as f64) / 2.0
    })
}
```

Register the module — in `src/rustle/vg_family/mod.rs`, next to the other `pub mod` lines:

```rust
pub mod single_copy; // O1 baseline: single-copy (χ(H)=1) loci + λ_global for the copy-number normalizer.
```

Add the formal-transcript doc — in `src/rustle/vg_family/denovo_assemble.rs`, prepend to the `assemble_gate` doc
comment (find `pub fn assemble_gate`):

```rust
/// THE FORMAL TRANSCRIPT DEFINITION. A transcript is an exact-intron-chain cluster of primary reads whose LOCUS
/// (the junction-incidence component, see `locus_support`) carries >= `min_reads` (GATE_MIN_READS) reads, whose
/// junctions are all canonical + consistent-strand, and whose spliced length is in `[min_spliced, max_spliced]`.
/// A locus is a copy family of size >= 1; single-copy is the χ(H)=1 boundary case (see `single_copy`).
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test --lib single_copy 2>&1 | tail -5`
Expected: `test result: ok. 4 passed`

- [ ] **Step 5: Commit**

```bash
git add src/rustle/vg_family/single_copy.rs src/rustle/vg_family/mod.rs src/rustle/vg_family/denovo_assemble.rs
git commit -m "feat(single_copy): SingleCopyLocus + single_copy_loci + lambda_global; formal transcript doc"
```

---

### Task 2: Genome-wide single-copy baseline (DRY core extraction)

**Files:**
- Modify: `src/rustle/vg_family/denovo_pipeline.rs`
- Test: in-crate (composition + refactor safety via the existing suite)

**Interfaces:**
- Consumes: Task 1's `single_copy::{SingleCopyLocus, single_copy_loci}`.
- Produces:
  - `fn gw_reps_and_catalog(bam_path, fasta_path, threads, win, min_copies, cfg) -> Result<(Vec<DenovoTranscript>, Vec<ColocatedFamily>)>`
  - `pub fn detect_single_copy_baseline_genome_wide(bam_path: &str, fasta_path: &str, threads: usize, win: u64, min_copies: usize, cfg: &DenovoConfig) -> Result<Vec<crate::vg_family::single_copy::SingleCopyLocus>>`

- [ ] **Step 1: Extract the core (behavior-preserving)**

In `denovo_pipeline.rs`, rename the CURRENT body of `detect_conflict_catalog_genome_wide` into a private helper
`gw_reps_and_catalog` that returns `(reps, catalog)` instead of `Ok(catalog)`. Concretely: change the function
signature's return type to `Result<(Vec<DenovoTranscript>, Vec<ColocatedFamily>)>`, rename it to
`gw_reps_and_catalog`, and change the final `Ok(catalog)` to `Ok((reps, catalog))`. Everything between is
UNCHANGED (the reps build, per-chrom edges, `colocated_families`). Note: `reps` is currently `let reps` (line
~24) — it is still in scope at the end, so returning it needs no new plumbing.

Then re-implement the public API as thin wrappers:

```rust
pub fn detect_conflict_catalog_genome_wide(
    bam_path: &str, fasta_path: &str, threads: usize, win: u64, min_copies: usize, cfg: &DenovoConfig,
) -> Result<Vec<ColocatedFamily>> {
    let (_reps, catalog) = gw_reps_and_catalog(bam_path, fasta_path, threads, win, min_copies, cfg)?;
    Ok(catalog)
}

/// Genome-wide single-copy baseline: the reps that no family claims, as lightweight records.
/// Same traversal as the conflict catalog (DRY via `gw_reps_and_catalog`); a locus is single-copy iff it is not
/// a copy of any emitted family. Feeds `λ_global` and the `.single_copy.tsv` baseline table.
pub fn detect_single_copy_baseline_genome_wide(
    bam_path: &str, fasta_path: &str, threads: usize, win: u64, min_copies: usize, cfg: &DenovoConfig,
) -> Result<Vec<crate::vg_family::single_copy::SingleCopyLocus>> {
    let (reps, catalog) = gw_reps_and_catalog(bam_path, fasta_path, threads, win, min_copies, cfg)?;
    Ok(crate::vg_family::single_copy::single_copy_loci(&reps, &catalog))
}
```

- [ ] **Step 2: Add a composition test**

Add to `denovo_pipeline.rs` tests (this pins the composition without needing a BAM — `single_copy_loci` is the
logic; the traversal is exercised by the existing gw tests + the Task 3 smoke run):

```rust
#[test]
fn single_copy_baseline_excludes_family_copies() {
    // rep set {a,b,c}; a 2-copy family {a,b}; baseline = {c}. (Direct on single_copy_loci -- the genome-wide
    // wrapper is the same call over the traversal's reps/catalog.)
    use crate::vg_family::single_copy::single_copy_loci;
    let mk = |tid: &str, s: u64| DenovoTranscript {
        tid: tid.into(), chrom: "c1".into(), start: s, end: s + 100, n_reads: 12, strand: '+',
        introns: vec![], seq: vec![],
    };
    let (a, b, c) = (mk("a", 0), mk("b", 1000), mk("c", 2000));
    let families = vec![ColocatedFamily { family_id: "F0".into(), chrom: "c1".into(), start: 0, end: 1100,
                                          copies: vec![a.clone(), b.clone()] }];
    let sc = single_copy_loci(&[a, b, c], &families);
    assert_eq!(sc.len(), 1);
    assert_eq!(sc[0].start, 2000);
}
```

- [ ] **Step 3: Run — refactor safety + new test**

Run: `cargo test --lib 2>&1 | tail -3`
Expected: all pass (the existing `detect_conflict_catalog_genome_wide` / gw tests are unchanged in behaviour;
the new composition test passes). Confirm no test that named the old private-body structure broke.

- [ ] **Step 4: Commit**

```bash
git add src/rustle/vg_family/denovo_pipeline.rs
git commit -m "feat(denovo_pipeline): gw_reps_and_catalog core; detect_single_copy_baseline_genome_wide"
```

---

### Task 3: `gw_family_catalog --single-copy-baseline` flag + emit

**Files:**
- Modify: `src/bin/gw_family_catalog.rs`
- Test: CLI/build + a real-data smoke run (foreground, per the crash rule)

**Interfaces:**
- Consumes: Task 2's `detect_single_copy_baseline_genome_wide`; Task 1's `single_copy::lambda_global`.

- [ ] **Step 1: Add the flag**

In `gw_family_catalog.rs` `struct Args`, add (next to the other flags):

```rust
    /// Emit the single-copy baseline instead of the family catalog: `<out>.single_copy.tsv` (one row per
    /// single-copy χ(H)=1 locus) + `<out>.lambda_global.tsv` (the genome-wide median n_reads = λ_global, the
    /// copy-number normalizer). A TABLE, not a GTF -- this is copy-number calibration, not an isoform catalog.
    #[arg(long, default_value_t = false)]
    single_copy_baseline: bool,
```

- [ ] **Step 2: Emit in `main`**

In `main`, before the existing catalog branch (near where `detect_conflict_catalog_genome_wide` is called),
add:

```rust
    if args.single_copy_baseline {
        use rustle::vg_family::single_copy::lambda_global;
        let loci = detect_single_copy_baseline_genome_wide(
            &args.bam, &args.fasta, args.threads, args.win, args.min_copies, &cfg,
        )?;
        let mut sc = std::fs::File::create(format!("{}.single_copy.tsv", args.out))?;
        writeln!(sc, "chrom\tstart\tend\tstrand\tn_reads\tn_exons\tchi_h\tn_psv")?;
        for l in &loci {
            writeln!(sc, "{}\t{}\t{}\t{}\t{}\t{}\t1\t0", l.chrom, l.start, l.end, l.strand, l.n_reads, l.n_exons)?;
        }
        let lam = lambda_global(&loci);
        let mut lf = std::fs::File::create(format!("{}.lambda_global.tsv", args.out))?;
        writeln!(lf, "lambda_global\tn_single_copy_loci")?;
        writeln!(lf, "{}\t{}", lam.map(|x| format!("{x}")).unwrap_or_else(|| "NA".into()), loci.len())?;
        eprintln!(
            "[gw-catalog] single-copy baseline: {} loci -> lambda_global={} ({}.single_copy.tsv, {}.lambda_global.tsv)",
            loci.len(), lam.map(|x| format!("{x}")).unwrap_or_else(|| "NA".into()), args.out, args.out
        );
        return Ok(());
    }
```

Ensure `detect_single_copy_baseline_genome_wide` is imported — add it to the existing
`use rustle::vg_family::denovo_pipeline::{ ... }` block, and add `use std::io::Write;` if not present.

- [ ] **Step 3: Build + smoke (foreground, serial, bounded)**

Run: `cargo build --release --bin gw_family_catalog 2>&1 | grep -c '^error'`  → Expected: `0`

⚠ **Crash rule:** a genome-wide run is heavy. For the smoke, verify on a SMALL BAM or accept the full run only
if it completes foreground. Minimal check that does not need a genome-wide run — the flag parses and the emit
path compiles:

```bash
/mnt/c/Users/jfris/Desktop/Rustle/target/release/gw_family_catalog --help 2>&1 | grep -A1 single-copy-baseline
```
Expected: the flag and its doc line appear.

If a bounded real run is desired, do it FOREGROUND on the existing test BAM only, outputs under
`/home/juanfra/winloci_scratch`, and inspect `*.single_copy.tsv` + `*.lambda_global.tsv`.

- [ ] **Step 4: Commit**

```bash
git add src/bin/gw_family_catalog.rs
git commit -m "feat(gw_family_catalog): --single-copy-baseline emits single_copy.tsv + lambda_global.tsv"
```

---

### Task 4: `copy_assign --lambda-file` consumes the cached λ

**Files:**
- Modify: `src/bin/copy_assign.rs`
- Test: in-bin unit test for the parser + precedence

**Interfaces:**
- Consumes: a `lambda_global.tsv` (header `lambda_global\tn_single_copy_loci`, then one data row).
- Produces: `fn read_lambda_file(path: &str) -> Option<f64>`; the resolved λ feeds the existing `depth_cn`.

- [ ] **Step 1: Write the failing test**

`copy_assign.rs` has no `mod tests` today; add one at end of file:

```rust
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_lambda_file_parses_the_scalar() {
        let dir = std::env::temp_dir();
        let p = dir.join(format!("rustle_lam_test_{}.tsv", std::process::id()));
        std::fs::write(&p, "lambda_global\tn_single_copy_loci\n25.5\t1234\n").unwrap();
        assert_eq!(read_lambda_file(p.to_str().unwrap()), Some(25.5));
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn read_lambda_file_none_on_na_or_missing() {
        let dir = std::env::temp_dir();
        let p = dir.join(format!("rustle_lam_na_{}.tsv", std::process::id()));
        std::fs::write(&p, "lambda_global\tn_single_copy_loci\nNA\t0\n").unwrap();
        assert_eq!(read_lambda_file(p.to_str().unwrap()), None);
        assert_eq!(read_lambda_file("/nonexistent/path.tsv"), None);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn resolve_lambda_precedence_explicit_over_file() {
        // explicit scalar wins over the file
        assert_eq!(resolve_lambda(Some(30.0), Some(25.0)), Some(30.0));
        assert_eq!(resolve_lambda(None, Some(25.0)), Some(25.0));
        assert_eq!(resolve_lambda(None, None), None);
    }
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --bin copy_assign 2>&1 | tail -5`
Expected: FAIL — `cannot find function read_lambda_file` / `resolve_lambda`.

- [ ] **Step 3: Implement**

Add the flag to `Args` (next to `lambda_global`):

```rust
    /// Read λ_global from a `gw_family_catalog --single-copy-baseline` `<prefix>.lambda_global.tsv` instead of
    /// passing the scalar by hand. `--lambda-global <f64>` (if given) takes precedence. This is how the
    /// copy-number normalizer becomes an in-binary quantity rather than an external-script number.
    #[arg(long)]
    lambda_file: Option<String>,
```

Add the helpers (near the top-level fns, e.g. above `fn main`):

```rust
/// Read the `lambda_global` scalar from a `lambda_global.tsv` (header + one data row, first column). `None` if
/// missing, unreadable, or the value is `NA`.
fn read_lambda_file(path: &str) -> Option<f64> {
    let text = std::fs::read_to_string(path).ok()?;
    let data = text.lines().nth(1)?; // skip header
    data.split('\t').next()?.trim().parse::<f64>().ok()
}

/// Resolve λ: explicit scalar > file > none.
fn resolve_lambda(explicit: Option<f64>, from_file: Option<f64>) -> Option<f64> {
    explicit.or(from_file)
}
```

Wire it where λ is used. Compute the resolved value once, after `args` is parsed (near where `cfg` is built):

```rust
    let lambda = resolve_lambda(args.lambda_global, args.lambda_file.as_deref().and_then(read_lambda_file));
```

Then change the `depth_cn` call site (currently `args.lambda_global.map(|lam| depth_cn(fa.n_reads, lam)).unwrap_or(f64::NAN)`)
to use `lambda`:

```rust
                    depth_cn: lambda.map(|lam| depth_cn(fa.n_reads, lam)).unwrap_or(f64::NAN),
```

Also update the later `famcn_readonly` emit if it independently references `args.lambda_global` (grep
`lambda_global` in the file and route every copy-number use through `lambda`). The progress line that says
`depth_cn on/off` should key on `lambda.is_some()`.

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test --bin copy_assign 2>&1 | tail -5`
Expected: `test result: ok. 3 passed`
Then: `cargo build --release --bin copy_assign 2>&1 | grep -c '^error'` → `0`

- [ ] **Step 5: Commit**

```bash
git add src/bin/copy_assign.rs
git commit -m "feat(copy_assign): --lambda-file consumes the in-binary lambda_global; precedence explicit>file"
```

---

## Self-Review

- **Spec coverage:** formal transcript doc (T1) · unification via χ(H)=1 + single_copy_loci (T1) · genome-wide
  pass reusing the traversal (T2) · flag emitting `single_copy.tsv` + `lambda_global.tsv` (T3) · `--lambda-file`
  internalization with precedence (T4). Non-goals (no GTF, no isoform quant, no `--gtf` reconciliation) are
  respected — nothing in any task touches the `--gtf` path.
- **Type consistency:** `SingleCopyLocus{chrom,start,end,strand,n_reads:u32,n_exons:usize}` produced T1, consumed
  T2/T3; `single_copy_loci(&[DenovoTranscript], &[ColocatedFamily])` and `lambda_global(&[SingleCopyLocus])`
  consistent across T1→T3; `detect_single_copy_baseline_genome_wide` produced T2, consumed T3; `read_lambda_file`
  / `resolve_lambda` internal to T4.
- **DRY:** the genome-wide traversal exists once (`gw_reps_and_catalog`), used by both the catalog and the
  baseline. λ resolution is one function.
- **Known-unknown flagged for the implementer:** T2 renames the body of an existing public function; verify no
  other caller of `detect_conflict_catalog_genome_wide` exists beyond `gw_family_catalog` (grep) before relying
  on the thin-wrapper swap — if `detect_conflict_catalog_genome_wide_xchrom` or a test calls it, they keep working
  because the public signature is unchanged.
