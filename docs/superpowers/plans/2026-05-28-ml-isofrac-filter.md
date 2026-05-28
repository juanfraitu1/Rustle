# ML Isofrac Filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the `longunder` gate in `isofrac_with_summary` with a depth-3 decision tree, selectable via `--filter-mode ml`, active only in de novo mode.

**Architecture:** New `ml_filter.rs` module holds the `MlFeatures` struct, a static JSONL writer, and a placeholder `ml_predict()`. Two new bool fields in `RunConfig` (`use_ml_filter`, `guide_mode`) flow into `isofrac_with_summary` which gains a new ML branch. A Python training script reads a feature dump and prints a Rust if-else tree to paste into `ml_predict`.

**Tech Stack:** Rust (clap `ValueEnum`, `serde_json`, `OnceLock`/`Mutex`), Python (scikit-learn, pandas).

---

## File Map

| Action | Path | Responsibility |
|---|---|---|
| Create | `src/rustle/ml_filter.rs` | `MlFeatures` struct, `from_pair()`, `ml_predict()` placeholder, JSONL dump writer |
| Modify | `src/rustle/lib.rs:57–65` | Register `pub mod ml_filter;` |
| Modify | `src/rustle/types.rs:788–792` + `1197` | Add `use_ml_filter: bool` and `guide_mode: bool` to `RunConfig` and its `Default` |
| Modify | `src/rustle/main.rs` | Add `FilterMode` enum with `ValueEnum`, `--filter-mode` arg, set both new config fields |
| Modify | `src/rustle/transcript_filter.rs:2092` | Extend `isofrac_with_summary` signature with two new params; add ML branch |
| Modify | `src/rustle/transcript_filter.rs:2071` | Update `isofrac` wrapper to pass `false, false` for the two new params |
| Modify | `src/rustle/transcript_filter.rs:7913` | Update call site to pass `config.use_ml_filter, config.guide_mode` |
| Modify | `src/rustle/main.rs` (flush) | Call `rustle::ml_filter::flush_ml_dump()` at end of `run_cli` |
| Create | `bench/train_ml_filter.py` | Label features from gffcompare, train tree, print Rust if-else |
| Create | `bench/run_ml_dump.sh` | Convenience wrapper for generating training data |

---

## Task 1: Create `src/rustle/ml_filter.rs`

**Files:**
- Create: `src/rustle/ml_filter.rs`

- [ ] **Step 1: Write a failing test for `MlFeatures::from_pair` clamping behavior**

```rust
// Add at the bottom of src/rustle/ml_filter.rs after writing the file in Step 2
#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Transcript;

    fn make_tx(cov: f64, longcov: f64, bpcov_cov: f64, min_jct_mm: f64, exons: Vec<(u64, u64)>) -> Transcript {
        let mut tx = Transcript::default();
        tx.coverage = cov;
        tx.longcov = longcov;
        tx.bpcov_cov = bpcov_cov;
        tx.min_jct_mm = min_jct_mm;
        tx.exons = exons;
        tx
    }

    #[test]
    fn test_from_pair_zero_dom_cov() {
        let k = make_tx(1.0, 2.0, 100.0, 5.0, vec![(0, 100), (200, 300)]);
        let dom = make_tx(0.0, 0.0, 0.0, 0.0, vec![(0, 500)]);
        let f = MlFeatures::from_pair(&k, &dom);
        assert_eq!(f.cov_ratio, 0.0);
        assert_eq!(f.longcov_ratio, 0.0);
        assert_eq!(f.bpcov_ratio, 0.0);
        assert_eq!(f.jct_mm_ratio, 0.0);
    }

    #[test]
    fn test_from_pair_bpcov_ratio_clamped() {
        // k bpcov/tlen = 1000/100 = 10; dom bpcov/tlen = 1/500 = 0.002 → ratio = 5000 → clamped to 10
        let k = make_tx(5.0, 3.0, 1000.0, 10.0, vec![(0, 100)]);
        let dom = make_tx(100.0, 50.0, 1.0, 200.0, vec![(0, 500)]);
        let f = MlFeatures::from_pair(&k, &dom);
        assert_eq!(f.bpcov_ratio, 10.0);
    }

    #[test]
    fn test_from_pair_is_intron_subset_true() {
        // k introns: (100,200); dom introns: (100,200),(400,500)
        let k = make_tx(5.0, 2.0, 50.0, 5.0, vec![(0, 100), (200, 300)]);
        let dom = make_tx(50.0, 20.0, 500.0, 50.0, vec![(0, 100), (200, 300), (400, 500), (600, 700)]);
        let f = MlFeatures::from_pair(&k, &dom);
        assert_eq!(f.is_intron_subset, 1.0);
    }

    #[test]
    fn test_from_pair_is_intron_subset_false() {
        // k has intron (100,300) not in dom
        let k = make_tx(5.0, 2.0, 50.0, 5.0, vec![(0, 100), (300, 400)]);
        let dom = make_tx(50.0, 20.0, 500.0, 50.0, vec![(0, 100), (200, 300), (400, 500)]);
        let f = MlFeatures::from_pair(&k, &dom);
        assert_eq!(f.is_intron_subset, 0.0);
    }

    #[test]
    fn test_placeholder_predict_always_true() {
        let k = make_tx(1.0, 1.0, 10.0, 1.0, vec![(0, 100), (200, 300)]);
        let dom = make_tx(100.0, 100.0, 1000.0, 100.0, vec![(0, 100), (200, 300)]);
        let f = MlFeatures::from_pair(&k, &dom);
        assert!(ml_predict(&f), "placeholder must always return true");
    }
}
```

- [ ] **Step 2: Run the test to verify it fails (file doesn't exist yet)**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo test ml_filter 2>&1 | head -20
```

Expected: error `file not found` or `module not found` for `ml_filter`.

- [ ] **Step 3: Write `src/rustle/ml_filter.rs`**

```rust
//! ML-based isofrac filter: depth-3 decision tree trained on GGO_19 transcript features.
//! Replaces the longunder gate in isofrac_with_summary when --filter-mode ml is set.
//! Only active in de novo mode (is_guided=false); guided mode always falls back to isofrac.

use std::io::{BufWriter, Write};
use std::sync::{Mutex, OnceLock};

use crate::types::Transcript;

/// Features computed per (minority, dominant) transcript pair for the ML filter.
#[derive(Debug, Clone, Copy)]
pub struct MlFeatures {
    /// cov_k / cov_dom — current isofrac basis
    pub cov_ratio: f64,
    /// longcov_k / longcov_dom — flow read count ratio
    pub longcov_ratio: f64,
    /// (bpcov_k/tlen_k) / (bpcov_dom/tlen_dom) — per-base BAM depth ratio, clamped [0, 10]
    pub bpcov_ratio: f64,
    /// min_jct_mm_k / min_jct_mm_dom — junction read ratio, clamped [0, 10]; 0 when dom has no junctions
    pub jct_mm_ratio: f64,
    /// Number of exons in minority transcript
    pub nexons_k: f64,
    /// tlen_k / tlen_dom — relative exonic length ratio
    pub tlen_ratio: f64,
    /// Absolute longcov of minority transcript (flow read count)
    pub longcov_abs: f64,
    /// 1.0 if all k introns ⊆ dom introns, else 0.0
    pub is_intron_subset: f64,
}

impl MlFeatures {
    /// Compute features from a (minority k, dominant dom) transcript pair.
    /// Safe for zero denominators: produces 0.0 ratio rather than panicking.
    pub fn from_pair(k: &Transcript, dom: &Transcript) -> Self {
        let tlen_k = k.exons.iter().map(|(s, e)| e - s).sum::<u64>().max(1) as f64;
        let tlen_dom = dom.exons.iter().map(|(s, e)| e - s).sum::<u64>().max(1) as f64;

        let cov_ratio = if dom.coverage > 0.0 {
            k.coverage / dom.coverage
        } else {
            0.0
        };

        let longcov_ratio = if dom.longcov > 0.0 {
            k.longcov / dom.longcov
        } else {
            0.0
        };

        let bpcov_k_per_tlen = k.bpcov_cov / tlen_k;
        let bpcov_dom_per_tlen = dom.bpcov_cov / tlen_dom;
        let bpcov_ratio = if bpcov_dom_per_tlen > 0.0 {
            (bpcov_k_per_tlen / bpcov_dom_per_tlen).clamp(0.0, 10.0)
        } else {
            0.0
        };

        let jct_mm_ratio = if dom.min_jct_mm > 0.0 {
            (k.min_jct_mm / dom.min_jct_mm).clamp(0.0, 10.0)
        } else {
            0.0
        };

        let introns_k: std::collections::HashSet<(u64, u64)> =
            k.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
        let introns_dom: std::collections::HashSet<(u64, u64)> =
            dom.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
        let is_intron_subset = if introns_k.is_empty() {
            0.0
        } else if introns_k.is_subset(&introns_dom) {
            1.0
        } else {
            0.0
        };

        Self {
            cov_ratio,
            longcov_ratio,
            bpcov_ratio,
            jct_mm_ratio,
            nexons_k: k.exons.len() as f64,
            tlen_ratio: tlen_k / tlen_dom,
            longcov_abs: k.longcov,
            is_intron_subset,
        }
    }

    /// Comma-separated intron string in 1-indexed format: "donor-acceptor,..."
    /// where donor = exon[i].end + 1 and acceptor = exon[i+1].start  (1-indexed inclusive).
    /// Matches the format used in parity JSONL for cross-referencing with gffcompare output.
    pub fn intron_chain_str(tx: &Transcript) -> String {
        tx.exons
            .windows(2)
            .map(|w| format!("{}-{}", w[0].1 + 1, w[1].0))
            .collect::<Vec<_>>()
            .join(",")
    }
}

// ── JSONL dump (RUSTLE_ML_FEATURE_DUMP=1) ─────────────────────────────────────

static ML_DUMP_WRITER: OnceLock<Mutex<Option<BufWriter<std::fs::File>>>> = OnceLock::new();

fn get_ml_dump_writer() -> &'static Mutex<Option<BufWriter<std::fs::File>>> {
    ML_DUMP_WRITER.get_or_init(|| {
        if std::env::var_os("RUSTLE_ML_FEATURE_DUMP").is_none() {
            return Mutex::new(None);
        }
        let path = std::env::var("RUSTLE_ML_FEATURE_DUMP_PATH")
            .unwrap_or_else(|_| "ml_features.jsonl".to_string());
        match std::fs::File::create(&path) {
            Ok(f) => Mutex::new(Some(BufWriter::with_capacity(64 * 1024, f))),
            Err(e) => {
                eprintln!("[ml_filter] cannot create {path}: {e}");
                Mutex::new(None)
            }
        }
    })
}

/// Emit one JSONL record for a minority transcript candidate.
/// Called from isofrac_with_summary when RUSTLE_ML_FEATURE_DUMP=1.
pub fn emit_ml_candidate(f: &MlFeatures, intron_chain: &str) {
    let lk = get_ml_dump_writer();
    if let Ok(mut g) = lk.lock() {
        if let Some(w) = g.as_mut() {
            let _ = writeln!(
                w,
                r#"{{"intron_chain":"{chain}","features":{{"cov_ratio":{cr:.6},"longcov_ratio":{lr:.6},"bpcov_ratio":{br:.6},"jct_mm_ratio":{jr:.6},"nexons_k":{nk:.1},"tlen_ratio":{tr:.6},"longcov_abs":{la:.1},"is_intron_subset":{is:.1}}}}}"#,
                chain = intron_chain,
                cr = f.cov_ratio,
                lr = f.longcov_ratio,
                br = f.bpcov_ratio,
                jr = f.jct_mm_ratio,
                nk = f.nexons_k,
                tr = f.tlen_ratio,
                la = f.longcov_abs,
                is = f.is_intron_subset,
            );
        }
    }
}

/// Flush the JSONL dump writer. Call at the end of main.
pub fn flush_ml_dump() {
    if let Some(lk) = ML_DUMP_WRITER.get() {
        if let Ok(mut g) = lk.lock() {
            if let Some(w) = g.as_mut() {
                let _ = w.flush();
            }
        }
    }
}

// ── Model ─────────────────────────────────────────────────────────────────────

/// Returns `true` → keep transcript, `false` → kill (longunder).
///
/// **Placeholder**: always keeps all transcripts (equivalent to `-f 0`).
/// Replace the body with the Rust if-else tree printed by `bench/train_ml_filter.py`
/// after training on GGO_19 feature dump.
pub fn ml_predict(_f: &MlFeatures) -> bool {
    true
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::Transcript;

    fn make_tx(cov: f64, longcov: f64, bpcov_cov: f64, min_jct_mm: f64, exons: Vec<(u64, u64)>) -> Transcript {
        let mut tx = Transcript::default();
        tx.coverage = cov;
        tx.longcov = longcov;
        tx.bpcov_cov = bpcov_cov;
        tx.min_jct_mm = min_jct_mm;
        tx.exons = exons;
        tx
    }

    #[test]
    fn test_from_pair_zero_dom_cov() {
        let k = make_tx(1.0, 2.0, 100.0, 5.0, vec![(0, 100), (200, 300)]);
        let dom = make_tx(0.0, 0.0, 0.0, 0.0, vec![(0, 500)]);
        let f = MlFeatures::from_pair(&k, &dom);
        assert_eq!(f.cov_ratio, 0.0);
        assert_eq!(f.longcov_ratio, 0.0);
        assert_eq!(f.bpcov_ratio, 0.0);
        assert_eq!(f.jct_mm_ratio, 0.0);
    }

    #[test]
    fn test_from_pair_bpcov_ratio_clamped() {
        // k bpcov/tlen = 1000/100 = 10.0; dom bpcov/tlen = 1/500 = 0.002 → raw ratio 5000 → clamped to 10
        let k = make_tx(5.0, 3.0, 1000.0, 10.0, vec![(0, 100)]);
        let dom = make_tx(100.0, 50.0, 1.0, 200.0, vec![(0, 500)]);
        let f = MlFeatures::from_pair(&k, &dom);
        assert_eq!(f.bpcov_ratio, 10.0);
    }

    #[test]
    fn test_from_pair_is_intron_subset_true() {
        let k = make_tx(5.0, 2.0, 50.0, 5.0, vec![(0, 100), (200, 300)]);
        let dom = make_tx(50.0, 20.0, 500.0, 50.0, vec![(0, 100), (200, 300), (400, 500), (600, 700)]);
        let f = MlFeatures::from_pair(&k, &dom);
        assert_eq!(f.is_intron_subset, 1.0);
    }

    #[test]
    fn test_from_pair_is_intron_subset_false() {
        let k = make_tx(5.0, 2.0, 50.0, 5.0, vec![(0, 100), (300, 400)]);
        let dom = make_tx(50.0, 20.0, 500.0, 50.0, vec![(0, 100), (200, 300), (400, 500)]);
        let f = MlFeatures::from_pair(&k, &dom);
        assert_eq!(f.is_intron_subset, 0.0);
    }

    #[test]
    fn test_placeholder_predict_always_true() {
        let k = make_tx(1.0, 1.0, 10.0, 1.0, vec![(0, 100), (200, 300)]);
        let dom = make_tx(100.0, 100.0, 1000.0, 100.0, vec![(0, 100), (200, 300)]);
        let f = MlFeatures::from_pair(&k, &dom);
        assert!(ml_predict(&f));
    }
}
```

- [ ] **Step 4: Register the module in `src/rustle/lib.rs`**

In `lib.rs` at line 60 (the Stage 7 transcript filtering block), add one line after `pub mod transcript_filter;`:

```rust
pub mod ml_filter;        // ML-based isofrac alternative (depth-3 decision tree)
```

- [ ] **Step 5: Run the tests to verify they pass**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo test ml_filter 2>&1 | tail -20
```

Expected output includes:
```
test ml_filter::tests::test_from_pair_zero_dom_cov ... ok
test ml_filter::tests::test_from_pair_bpcov_ratio_clamped ... ok
test ml_filter::tests::test_from_pair_is_intron_subset_true ... ok
test ml_filter::tests::test_from_pair_is_intron_subset_false ... ok
test ml_filter::tests::test_placeholder_predict_always_true ... ok
test result: ok. 5 passed; 0 failed
```

- [ ] **Step 6: Commit**

```bash
git add src/rustle/ml_filter.rs src/rustle/lib.rs
git commit -m "feat(ml_filter): add MlFeatures struct, placeholder ml_predict, and JSONL dump writer"
```

---

## Task 2: Add `use_ml_filter` and `guide_mode` to `RunConfig`

**Files:**
- Modify: `src/rustle/types.rs`

- [ ] **Step 1: Write a failing test verifying defaults**

Add to the existing test block in `types.rs` (search for `#[cfg(test)]` at the bottom) or in a new block:

```rust
#[cfg(test)]
mod config_defaults_tests {
    use super::*;
    #[test]
    fn test_runconfig_ml_defaults() {
        let c = RunConfig::default();
        assert!(!c.use_ml_filter, "use_ml_filter must default to false");
        assert!(!c.guide_mode, "guide_mode must default to false");
    }
}
```

- [ ] **Step 2: Run to verify it fails**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo test test_runconfig_ml_defaults 2>&1 | tail -10
```

Expected: `error[E0609]: no field 'use_ml_filter' on type 'RunConfig'`

- [ ] **Step 3: Add the two fields to the `RunConfig` struct in `src/rustle/types.rs`**

Find the block around line 792 (`pub transcript_isofrac_keep_min: f64,`). Add two lines immediately after it:

```rust
    /// When true, use the ML depth-3 decision tree instead of the isofrac/longunder formula.
    /// Set by --filter-mode ml. Only applies in de novo mode (guide_mode=false).
    pub use_ml_filter: bool,
    /// True when a guide GTF was supplied (-G). Used to disable the ML filter in guided mode.
    pub guide_mode: bool,
```

- [ ] **Step 4: Add defaults in `impl Default for RunConfig`**

Find the `transcript_isofrac_keep_min: 1.0,` line (around line 1197). Add immediately after it:

```rust
            use_ml_filter: false,
            guide_mode: false,
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo test test_runconfig_ml_defaults 2>&1 | tail -10
```

Expected: `test config_defaults_tests::test_runconfig_ml_defaults ... ok`

- [ ] **Step 6: Commit**

```bash
git add src/rustle/types.rs
git commit -m "feat(config): add use_ml_filter and guide_mode fields to RunConfig"
```

---

## Task 3: Wire `--filter-mode` CLI arg in `main.rs`

**Files:**
- Modify: `src/rustle/main.rs`

- [ ] **Step 1: Add the `FilterMode` enum and the CLI arg field**

In `src/rustle/main.rs`, find the `use` imports at the top. Add after the existing `use` lines:

```rust
use clap::ValueEnum;
```

Then, find the `#[derive(Parser)]` struct block (the `Args` struct). Find the field `transcript_isofrac: f64` (around line 107). Add a new field immediately after `transcript_isofrac_keep_min`:

```rust
    /// Transcript filter mode: isofrac (default) or ml (depth-3 decision tree, de novo only)
    #[arg(long, value_enum, default_value = "isofrac")]
    filter_mode: FilterMode,
```

Then define the `FilterMode` enum just above the `#[derive(Parser)]` struct definition (before `struct Args {`):

```rust
#[derive(Clone, Copy, PartialEq, Eq, ValueEnum)]
enum FilterMode {
    Isofrac,
    Ml,
}
```

- [ ] **Step 2: Set the two new fields when building `RunConfig`**

Find the `let mut config = RunConfig {` block (around line 675). Add the two new fields inside it, after the `transcript_isofrac_keep_min: args.transcript_isofrac_keep_min,` line:

```rust
        use_ml_filter: args.filter_mode == FilterMode::Ml,
        guide_mode: args.guide.is_some(),
```

- [ ] **Step 3: Add `flush_ml_dump()` call at the end of `run_cli`**

Find the end of `run_cli()` (around line 863, just before `Ok(())`). Add:

```rust
    rustle::ml_filter::flush_ml_dump();
```

- [ ] **Step 4: Verify it compiles cleanly**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build 2>&1 | grep -E "^error|warning.*unused" | head -20
```

Expected: zero `error` lines. A few existing `unused` warnings are acceptable.

- [ ] **Step 5: Smoke-test the new flag is recognized**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
./target/debug/rustle --help 2>&1 | grep -A2 "filter-mode"
```

Expected:
```
--filter-mode <FILTER_MODE>
    Transcript filter mode: isofrac (default) or ml (depth-3 decision tree, de novo only)
    [default: isofrac] [possible values: isofrac, ml]
```

- [ ] **Step 6: Commit**

```bash
git add src/rustle/main.rs
git commit -m "feat(cli): add --filter-mode isofrac|ml flag and guide_mode wiring"
```

---

## Task 4: Wire ML branch into `isofrac_with_summary`

**Files:**
- Modify: `src/rustle/transcript_filter.rs`

- [ ] **Step 1: Add the import for `ml_filter` at the top of `transcript_filter.rs`**

Find the `use` block at the top of `transcript_filter.rs`. Add:

```rust
use crate::ml_filter::{emit_ml_candidate, ml_predict, MlFeatures};
```

- [ ] **Step 2: Extend the `isofrac_with_summary` signature**

Current signature (line 2092):
```rust
fn isofrac_with_summary(
    transcripts: Vec<Transcript>,
    isofrac: f64,
    error_perc: f64,
    drop: f64,
    _bpcov: Option<&Bpcov>,
    verbose: bool,
    keep_min_abundance: f64,
) -> (Vec<Transcript>, IsofracKillSummary) {
```

Replace with:
```rust
fn isofrac_with_summary(
    transcripts: Vec<Transcript>,
    isofrac: f64,
    error_perc: f64,
    drop: f64,
    _bpcov: Option<&Bpcov>,
    verbose: bool,
    keep_min_abundance: f64,
    use_ml_filter: bool,
    is_guided: bool,
) -> (Vec<Transcript>, IsofracKillSummary) {
```

- [ ] **Step 3: Update the `isofrac` wrapper (line 2071) to pass the two new trailing `false` params**

Current call in `isofrac` wrapper body (line 2080):
```rust
    isofrac_with_summary(
        transcripts,
        isofrac,
        error_perc,
        drop,
        _bpcov,
        verbose,
        keep_min_abundance,
    )
    .0
```

Replace with:
```rust
    isofrac_with_summary(
        transcripts,
        isofrac,
        error_perc,
        drop,
        _bpcov,
        verbose,
        keep_min_abundance,
        false,  // use_ml_filter: always isofrac in diagnostic wrapper
        false,  // is_guided
    )
    .0
```

- [ ] **Step 4: Update the main call site (line 7913) to pass the two new config fields**

Current call (lines 7913–7924):
```rust
        let (txs_after_isofrac, isofrac_summary) = isofrac_with_summary(
            txs,
            config.transcript_isofrac,
            config.pairwise_error_perc,
            config.pairwise_drop,
            bpcov,
            config.verbose,
            config.transcript_isofrac_keep_min,
        );
```

Replace with:
```rust
        let (txs_after_isofrac, isofrac_summary) = isofrac_with_summary(
            txs,
            config.transcript_isofrac,
            config.pairwise_error_perc,
            config.pairwise_drop,
            bpcov,
            config.verbose,
            config.transcript_isofrac_keep_min,
            config.use_ml_filter,
            config.guide_mode,
        );
```

- [ ] **Step 5: Add the feature extraction and ML branch inside `isofrac_with_summary`**

Find the block that computes `let mut longunder` (in the per-window loop, around where the existing `combined_factor` is used). It currently reads:

```rust
            let mut longunder = if txs[k].exons.len() > 1 {
                (cmp_multicov <= 0.0
                    && cov < isofraclong * cmp_usedcov * combined_factor
                    && cov < drop / error_perc)
                    || cov < isofraclong * cmp_multicov * combined_factor
            } else {
                cov < isofraclong * cmp_usedcov * combined_factor
            };
```

Replace the entire `let mut longunder` block with:

```rust
            let mut longunder = if use_ml_filter && !is_guided && txs[k].exons.len() > 1 {
                // ML branch: de novo multi-exon only.
                // Compute features and query the decision tree.
                let features = MlFeatures::from_pair(&txs[k], &txs[first]);
                if std::env::var_os("RUSTLE_ML_FEATURE_DUMP").is_some() {
                    let chain = MlFeatures::intron_chain_str(&txs[k]);
                    emit_ml_candidate(&features, &chain);
                }
                !ml_predict(&features)
            } else {
                // Isofrac branch: original longunder formula.
                if txs[k].exons.len() > 1 {
                    (cmp_multicov <= 0.0
                        && cov < isofraclong * cmp_usedcov * combined_factor
                        && cov < drop / error_perc)
                        || cov < isofraclong * cmp_multicov * combined_factor
                } else {
                    cov < isofraclong * cmp_usedcov * combined_factor
                }
            };
```

- [ ] **Step 6: Build in release to verify correctness**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | grep "^error" | head -20
```

Expected: zero error lines.

- [ ] **Step 7: Verify default mode produces identical output to pre-change baseline**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
./target/release/rustle -L /mnt/c/Users/jfris/Desktop/GGO_19.bam -o /tmp/ml_isofrac_test.gtf 2>/dev/null
gffcompare -RQ -r /mnt/c/Users/jfris/Desktop/GGO_19.gtf /tmp/ml_isofrac_test.gtf -o /tmp/ml_isofrac_cmp 2>/dev/null
grep "Intron chain level" /tmp/ml_isofrac_cmp.stats
```

Expected: `Intron chain level:    96.1     |    91.0` — identical to pre-change baseline. (Default mode is `--filter-mode isofrac`; no behavior change.)

- [ ] **Step 8: Verify placeholder ML mode produces low-precision (keep-all) output**

```bash
./target/release/rustle -L /mnt/c/Users/jfris/Desktop/GGO_19.bam --filter-mode ml -o /tmp/ml_placeholder_test.gtf 2>/dev/null
gffcompare -RQ -r /mnt/c/Users/jfris/Desktop/GGO_19.gtf /tmp/ml_placeholder_test.gtf -o /tmp/ml_placeholder_cmp 2>/dev/null
grep "Intron chain level\|Query mRNAs" /tmp/ml_placeholder_cmp.stats
```

Expected: precision significantly lower than 91.0% (placeholder keeps all minority transcripts → many FPs), sensitivity near 100%.

- [ ] **Step 9: Commit**

```bash
git add src/rustle/transcript_filter.rs
git commit -m "feat(filter): wire ML branch into isofrac_with_summary, de novo only"
```

---

## Task 5: Write the training script `bench/train_ml_filter.py`

**Files:**
- Create: `bench/train_ml_filter.py`
- Create: `bench/run_ml_dump.sh`

- [ ] **Step 1: Create `bench/run_ml_dump.sh`**

```bash
#!/usr/bin/env bash
# Generate ML feature dump for training.
# Usage: bash bench/run_ml_dump.sh <BAM> <REF_GTF> [output_prefix]
# Outputs: ml_features.jsonl and ml_all_transcripts.gtf in the working directory.
set -euo pipefail

BAM=${1:?Usage: run_ml_dump.sh <BAM> <REF_GTF> [prefix]}
REF_GTF=${2:?}
PREFIX=${3:-ml_training}

RUSTLE=./target/release/rustle
GTF_OUT="${PREFIX}_all_transcripts.gtf"
CMP_OUT="${PREFIX}_cmp"

echo "[run_ml_dump] Running rustle with -f 0 and feature dump..."
RUSTLE_ML_FEATURE_DUMP=1 \
RUSTLE_ML_FEATURE_DUMP_PATH="${PREFIX}_features.jsonl" \
"$RUSTLE" -L "$BAM" -f 0 -o "$GTF_OUT" 2>/dev/null

echo "[run_ml_dump] Running gffcompare..."
gffcompare -RQ -r "$REF_GTF" "$GTF_OUT" -o "$CMP_OUT" 2>/dev/null

echo "[run_ml_dump] Done."
echo "  Features : ${PREFIX}_features.jsonl"
echo "  GTF      : $GTF_OUT"
echo "  tmap     : ${CMP_OUT}.${GTF_OUT##*/}.tmap"
```

- [ ] **Step 2: Create `bench/train_ml_filter.py`**

```python
#!/usr/bin/env python3
"""
Train a depth-3 decision tree on GGO_19 ML feature dump.

Usage:
    python3 bench/train_ml_filter.py \
        ml_training_features.jsonl \
        ml_training_all_transcripts.gtf \
        ml_training_cmp.ml_training_all_transcripts.gtf.tmap

Prints the trained tree as a Rust pub fn ml_predict body to stdout.
Paste the printed function into src/rustle/ml_filter.rs replacing the
placeholder ml_predict body.
"""

import sys
import json
import re
from pathlib import Path

import pandas as pd
from sklearn.tree import DecisionTreeClassifier, export_text
from sklearn.model_selection import cross_val_score
import numpy as np


FEATURE_NAMES = [
    "cov_ratio", "longcov_ratio", "bpcov_ratio", "jct_mm_ratio",
    "nexons_k", "tlen_ratio", "longcov_abs", "is_intron_subset",
]


def parse_features(jsonl_path: str) -> dict:
    """Returns dict: intron_chain_str -> feature dict."""
    records = {}
    with open(jsonl_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            obj = json.loads(line)
            chain = obj["intron_chain"]
            records[chain] = obj["features"]
    return records


def parse_gtf_intron_chains(gtf_path: str) -> dict:
    """Returns dict: tx_id -> intron_chain_str (1-indexed start, 0-indexed end format)."""
    tx_exons: dict = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#") or "\ttranscript\t" in line:
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            start, end = int(fields[3]), int(fields[4])
            attr = fields[8]
            m = re.search(r'transcript_id "([^"]+)"', attr)
            if not m:
                continue
            tx_id = m.group(1)
            tx_exons.setdefault(tx_id, []).append((start, end))

    chains = {}
    for tx_id, exons in tx_exons.items():
        exons.sort()
        if len(exons) < 2:
            chains[tx_id] = ""
            continue
        # intron: end_of_exon_i (GTF inclusive) + 1 for 1-based start,
        # start_of_next_exon (GTF 1-based = 0-based end of intron)
        parts = [f"{exons[i][1] + 1}-{exons[i+1][0]}" for i in range(len(exons) - 1)]
        chains[tx_id] = ",".join(parts)
    return chains


def parse_tmap(tmap_path: str) -> dict:
    """Returns dict: qry_id -> class_code."""
    cc = {}
    with open(tmap_path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            # columns: ref_id, class_code, qry_id, ...
            cc[parts[2]] = parts[1]
    return cc


def tree_to_rust(clf: DecisionTreeClassifier, feature_names: list, cv_auc: float) -> str:
    """Convert a fitted DecisionTreeClassifier to a Rust if-else function body."""
    tree = clf.tree_

    def node_to_rust(node: int, indent: int) -> str:
        pad = "    " * indent
        if tree.children_left[node] == tree.children_right[node]:
            # Leaf: majority class
            values = tree.value[node][0]
            n_kill, n_keep = int(values[0]), int(values[1])
            keep = n_keep >= n_kill
            verdict = "true" if keep else "false"
            return f"{pad}return {verdict}; // keep={n_keep}, kill={n_kill}\n"

        feat = feature_names[tree.feature[node]]
        thresh = tree.threshold[node]
        left = node_to_rust(tree.children_left[node], indent + 1)
        right = node_to_rust(tree.children_right[node], indent + 1)
        return (
            f"{pad}if f.{feat} <= {thresh:.6f} {{\n"
            f"{left}"
            f"{pad}}}\n"
            f"{right}"
        )

    body = node_to_rust(0, 1)
    date = pd.Timestamp.now().strftime("%Y-%m-%d")
    return (
        f"/// Returns `true` → keep transcript, `false` → kill (longunder).\n"
        f"/// Trained on GGO_19 de novo, {date}, 5-fold CV AUC = {cv_auc:.3f}\n"
        f"pub fn ml_predict(f: &MlFeatures) -> bool {{\n"
        f"{body}"
        f"}}\n"
    )


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    jsonl_path, gtf_path, tmap_path = sys.argv[1], sys.argv[2], sys.argv[3]

    print(f"[train] Loading features from {jsonl_path}...", file=sys.stderr)
    feature_records = parse_features(jsonl_path)
    print(f"[train] {len(feature_records)} feature records loaded", file=sys.stderr)

    print(f"[train] Parsing intron chains from {gtf_path}...", file=sys.stderr)
    tx_chains = parse_gtf_intron_chains(gtf_path)

    print(f"[train] Parsing labels from {tmap_path}...", file=sys.stderr)
    tx_class = parse_tmap(tmap_path)

    # Build chain -> class_code mapping via tx_id join
    chain_to_class: dict = {}
    for tx_id, chain in tx_chains.items():
        if chain and tx_id in tx_class:
            chain_to_class[chain] = tx_class[tx_id]

    # Build training rows: join feature_records (keyed by chain) with labels
    rows = []
    for chain, feats in feature_records.items():
        if chain not in chain_to_class:
            continue  # no label — skip
        label = 1 if chain_to_class[chain] == "=" else 0
        row = {k: feats[k] for k in FEATURE_NAMES}
        row["label"] = label
        rows.append(row)

    df = pd.DataFrame(rows)
    print(f"[train] {len(df)} labeled rows ({df['label'].sum():.0f} keep, {(1-df['label']).sum():.0f} kill)",
          file=sys.stderr)

    if len(df) < 20:
        print("[train] ERROR: fewer than 20 labeled examples — check paths or re-run dump.", file=sys.stderr)
        sys.exit(1)

    X = df[FEATURE_NAMES].values
    y = df["label"].values

    clf = DecisionTreeClassifier(max_depth=3, class_weight="balanced", random_state=42)
    cv_scores = cross_val_score(clf, X, y, cv=5, scoring="roc_auc")
    print(f"[train] 5-fold CV AUC: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}", file=sys.stderr)

    clf.fit(X, y)

    text = export_text(clf, feature_names=FEATURE_NAMES)
    print("[train] Sklearn tree structure:", file=sys.stderr)
    print(text, file=sys.stderr)

    rust_fn = tree_to_rust(clf, FEATURE_NAMES, cv_scores.mean())
    print(rust_fn)


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Make `run_ml_dump.sh` executable**

```bash
chmod +x /mnt/c/Users/jfris/Desktop/Rustle/bench/run_ml_dump.sh
```

- [ ] **Step 4: Commit**

```bash
git add bench/train_ml_filter.py bench/run_ml_dump.sh
git commit -m "feat(bench): add ML filter training script and dump wrapper"
```

---

## Task 6: Generate training data, train, paste tree, benchmark

**Files:**
- Modify: `src/rustle/ml_filter.rs` (replace `ml_predict` body)

- [ ] **Step 1: Build release binary**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>&1 | grep "^error" | head -5
```

Expected: zero errors. Binary at `target/release/rustle`.

- [ ] **Step 2: Generate the feature dump**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
bash bench/run_ml_dump.sh \
    /mnt/c/Users/jfris/Desktop/GGO_19.bam \
    /mnt/c/Users/jfris/Desktop/GGO_19.gtf \
    /tmp/ml_ggo19
```

Expected output:
```
[run_ml_dump] Running rustle with -f 0 and feature dump...
[run_ml_dump] Running gffcompare...
[run_ml_dump] Done.
  Features : /tmp/ml_ggo19_features.jsonl
  GTF      : /tmp/ml_ggo19_all_transcripts.gtf
  tmap     : /tmp/ml_ggo19_cmp.ml_ggo19_all_transcripts.gtf.tmap
```

Verify the dump is non-empty:
```bash
wc -l /tmp/ml_ggo19_features.jsonl
```

Expected: ≥ 300 lines.

- [ ] **Step 3: Train the model**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
python3 bench/train_ml_filter.py \
    /tmp/ml_ggo19_features.jsonl \
    /tmp/ml_ggo19_all_transcripts.gtf \
    /tmp/ml_ggo19_cmp.ml_ggo19_all_transcripts.gtf.tmap \
    > /tmp/ml_predict_trained.rs
cat /tmp/ml_predict_trained.rs
```

Expected: prints a `pub fn ml_predict(f: &MlFeatures) -> bool { ... }` function. CV AUC should be > 0.70.

- [ ] **Step 4: Paste the trained tree into `src/rustle/ml_filter.rs`**

In `src/rustle/ml_filter.rs`, replace the entire `ml_predict` function (the placeholder at the bottom of the file):

```rust
// BEFORE (placeholder):
pub fn ml_predict(_f: &MlFeatures) -> bool {
    true
}
```

Replace with the complete function printed in Step 3 (including the doc comment with training date and CV AUC). The function signature must remain exactly:

```rust
pub fn ml_predict(f: &MlFeatures) -> bool {
```

Note: remove the leading `_` from `_f` since the trained tree uses `f`.

- [ ] **Step 5: Build and run the unit tests**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo test ml_filter 2>&1 | tail -15
```

Expected: all 5 tests pass. The `test_placeholder_predict_always_true` test will now FAIL because the trained tree no longer returns `true` for all inputs. **This is expected and correct — rename or update that test:**

Update `test_placeholder_predict_always_true` to a new test that verifies the trained tree kills an obvious noise transcript (very low cov_ratio, jct_mm_ratio, longcov_abs):

```rust
#[test]
fn test_trained_predict_kills_noise() {
    // A transcript with 0.1% coverage ratio and no junction support
    // should be killed by any reasonable trained tree.
    let k = make_tx(0.3, 1.0, 3.0, 0.0, vec![(0, 100), (200, 300)]);
    let dom = make_tx(300.0, 295.0, 3000.0, 300.0, vec![(0, 100), (200, 300)]);
    let f = MlFeatures::from_pair(&k, &dom);
    // cov_ratio=0.001, jct_mm_ratio=0.0, longcov_abs=1.0 — should be killed
    assert!(!ml_predict(&f), "noise transcript (0.1% cov, 0 jct support) should be killed");
}

#[test]
fn test_trained_predict_keeps_strong() {
    // A transcript with 15% cov_ratio, good bpcov, and longcov=10 should survive
    let k = make_tx(15.0, 10.0, 150.0, 15.0, vec![(0, 500), (700, 1000), (1200, 1500)]);
    let dom = make_tx(100.0, 80.0, 1000.0, 100.0, vec![(0, 500), (700, 1000)]);
    let f = MlFeatures::from_pair(&k, &dom);
    // cov_ratio=0.15, bpcov_ratio≈0.9, jct_mm_ratio=0.15, longcov_abs=10
    assert!(ml_predict(&f), "well-supported minority isoform should be kept");
}
```

Run again:
```bash
cargo test ml_filter 2>&1 | tail -15
```

Expected: all tests pass.

- [ ] **Step 6: Benchmark the trained ML filter on GGO_19**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
cargo build --release 2>/dev/null
./target/release/rustle -L /mnt/c/Users/jfris/Desktop/GGO_19.bam \
    --filter-mode ml -o /tmp/ml_trained_test.gtf 2>/dev/null
gffcompare -RQ -r /mnt/c/Users/jfris/Desktop/GGO_19.gtf \
    /tmp/ml_trained_test.gtf -o /tmp/ml_trained_cmp 2>/dev/null
echo "=== ML filter ==="
grep "Intron chain level\|Query mRNAs\|Matching intron" /tmp/ml_trained_cmp.stats

echo "=== Isofrac baseline ==="
grep "Intron chain level\|Query mRNAs\|Matching intron" /tmp/ggo19_default_check_cmp.stats
```

Expected: ML filter Sn ≥ 96.1% (should rescue isofrac victims), Pr as close to 91.0% as possible. Record the numbers.

- [ ] **Step 7: Verify guided mode is unaffected**

```bash
./target/release/rustle -L /mnt/c/Users/jfris/Desktop/GGO_19.bam \
    -G /mnt/c/Users/jfris/Desktop/GGO_19.gtf \
    --filter-mode ml -o /tmp/ml_guided_test.gtf 2>/dev/null
gffcompare -RQ -r /mnt/c/Users/jfris/Desktop/GGO_19.gtf \
    /tmp/ml_guided_test.gtf -o /tmp/ml_guided_cmp 2>/dev/null
echo "=== ML guided (must fall back to isofrac) ==="
grep "Intron chain level" /tmp/ml_guided_cmp.stats
```

Expected: matches the guided baseline `100.0% Sn / 99.7% Pr` — identical to running without `--filter-mode ml`.

- [ ] **Step 8: Commit**

```bash
git add src/rustle/ml_filter.rs
git commit -m "feat(ml_filter): paste trained depth-3 decision tree from GGO_19 training"
```
