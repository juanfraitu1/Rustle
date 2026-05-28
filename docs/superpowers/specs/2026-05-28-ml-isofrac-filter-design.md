# ML Isofrac Filter Design

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace Rustle's isofrac/longunder gate with a depth-3 decision tree trained on GGO_19 transcript features, selectable at runtime via `--filter-mode ml`.

**Architecture:** Feature dump → Python training → hardcoded Rust tree → CLI flag routes execution.

**Tech stack:** Rust (clap, serde_json), Python (scikit-learn, pandas).

---

## Problem Statement

Rustle's isofrac filter uses `cov_k / cov_dom < threshold` to kill minority isoforms. `cov` is derived from flow decomposition, which depletes coverage from minority paths as dominant paths are extracted first. StringTie's `pred->cov` uses a different formula (`flow × bpcov / tlen²`), so the same numeric threshold behaves differently between the two tools. The result: Rustle kills real minority isoforms (STRG.267.3, STRG.343.2, STRG.15.1) that StringTie keeps, because the relative coverage ratio is distorted by flow depletion.

A decision tree trained on multiple coverage-related signals (including BAM-depth and junction-read features that are immune to flow depletion) can learn a more stable decision boundary without being tied to one coverage formula.

---

## Scope

**In scope:**
- Replace the `longunder` boolean computation inside `isofrac_with_summary` in `transcript_filter.rs`
- New `--filter-mode isofrac|ml` CLI flag (default: `isofrac`)
- Feature dump mode (`RUSTLE_ML_FEATURE_DUMP=1`) for training data generation
- Training script `bench/train_ml_filter.py`
- Model embedded in `src/rustle/ml_filter.rs`

**Out of scope:**
- Pairwise overlap filter (`retainedintron_like`, `included_drop`) — untouched
- Cross-strand KRI filter — untouched (net +1.26pp F1, separate mechanism)
- Transcript extraction / path enumeration — untouched
- Any change to default behavior (`--filter-mode isofrac` is default)

---

## Architecture

### Data Flow

```
Step 1 — Generate training data:
  GGO_19.bam
    → rustle -f 0 RUSTLE_ML_FEATURE_DUMP=1 -o all_transcripts.gtf
    → writes ml_features.jsonl  (one line per minority transcript per window)
    → gffcompare -r GGO_19.gtf all_transcripts.gtf -o cmp
    → cmp.all_transcripts.gtf.tmap  (class_code per transcript)

Step 2 — Train:
  bench/train_ml_filter.py ml_features.jsonl cmp.all_transcripts.gtf.tmap
    → prints Rust if-else tree
    → operator pastes the printed function verbatim, replacing the placeholder body of ml_predict() in ml_filter.rs

Step 3 — Use:
  rustle --filter-mode ml [other flags] sample.bam
    → ml_predict() called instead of longunder formula
    → same pairwise filter, same KRI, same output format
```

### New Files

- `src/rustle/ml_filter.rs` — `MlFeatures` struct + `ml_predict()` (hardcoded tree)
- `bench/train_ml_filter.py` — training script
- `bench/run_ml_dump.sh` — convenience wrapper for Step 1

### Modified Files

- `src/rustle/main.rs` — `FilterMode` enum + `--filter-mode` arg
- `src/rustle/pipeline.rs` — pass `filter_mode` through `RunConfig` to the filter call
- `src/rustle/transcript_filter.rs` — feature extraction + ML branch in `isofrac_with_summary`

---

## Features

8 features computed per `(txs[k], txs[first])` pair where `first` is the dominant (highest-coverage transcript in the current isofrac window):

| Name | Formula | Notes |
|---|---|---|
| `cov_ratio` | `txs[k].coverage / txs[first].coverage` | Current isofrac basis; still informative |
| `longcov_ratio` | `txs[k].longcov / txs[first].longcov` | Flow read count; less formula-dependent |
| `bpcov_ratio` | `(bpcov_k/tlen_k) / (bpcov_dom/tlen_dom)` | Actual BAM depth per base; bypasses flow formula |
| `jct_mm_ratio` | `txs[k].min_jct_mm / txs[first].min_jct_mm` | BAM junction reads; immune to flow depletion |
| `nexons_k` | `txs[k].exons.len() as f64` | Complexity; real isoforms are multi-exon |
| `tlen_ratio` | `tlen_k / tlen_dom` | Relative exonic length |
| `longcov_abs` | `txs[k].longcov` | Absolute read support floor |
| `is_intron_subset` | 1.0 if all k introns ∈ dom introns, else 0.0 | Structural containment signal |

**Clamping:** `bpcov_ratio` and `jct_mm_ratio` clamped to `[0.0, 10.0]` to prevent outliers from distorting training. Dominants with zero `min_jct_mm` or zero `bpcov_cov` produce 0.0 for the respective ratio.

---

## Feature Dump Format

When `RUSTLE_ML_FEATURE_DUMP=1` is set, one JSON line per minority transcript per isofrac window is appended to `ml_features.jsonl` in the working directory:

```json
{
  "tx_id": "RSTL.267.3",
  "intron_chain": "38797505-38800755,38800899-38802412,...",
  "label": null,
  "features": {
    "cov_ratio": 0.00541,
    "longcov_ratio": 0.00678,
    "bpcov_ratio": 0.882,
    "jct_mm_ratio": 0.05,
    "nexons_k": 25.0,
    "tlen_ratio": 1.52,
    "longcov_abs": 2.0,
    "is_intron_subset": 0.0
  }
}
```

`intron_chain` is the comma-separated intron string (same format as parity JSONL) used to match gffcompare output.

---

## Training Script (`bench/train_ml_filter.py`)

Inputs:
- `ml_features.jsonl` — feature dump from Step 1
- `tmap` — gffcompare `.tmap` file (columns: ref_id, class_code, qry_id, ...)

Label assignment:
- Parse tmap: build dict `qry_id → class_code`
- For each JSONL record: label = 1 if class_code == `'='`, else 0
- Transcripts not in tmap → skip (shouldn't happen with `-f 0` run)

Training:
- `DecisionTreeClassifier(max_depth=3, class_weight='balanced', random_state=42)`
- Fit on all labeled examples
- Print tree structure as Rust if-else code (using `export_text` + custom formatter)
- Print 5-fold cross-validation AUC for sanity check

Output format (printed to stdout, paste into `ml_predict()`):
```rust
// Trained on GGO_19, 2026-05-28, CV-AUC=0.XX
pub fn ml_predict(f: &MlFeatures) -> bool {
    if f.bpcov_ratio < 0.15 {
        if f.jct_mm_ratio < 0.02 {
            return false;
        }
        return f.longcov_abs >= 3.0;
    }
    f.cov_ratio >= 0.005
}
```

---

## Model File (`src/rustle/ml_filter.rs`)

```rust
pub struct MlFeatures {
    pub cov_ratio: f64,
    pub longcov_ratio: f64,
    pub bpcov_ratio: f64,
    pub jct_mm_ratio: f64,
    pub nexons_k: f64,
    pub tlen_ratio: f64,
    pub longcov_abs: f64,
    pub is_intron_subset: f64,
}

/// Returns true → keep transcript, false → kill (longunder).
/// Placeholder: keep all (equivalent to -f 0). Replace body after training.
pub fn ml_predict(f: &MlFeatures) -> bool {
    let _ = f;
    true
}
```

---

## CLI Integration

In `main.rs`:

```rust
#[derive(Clone, ValueEnum)]
pub enum FilterMode {
    Isofrac,
    Ml,
}

// Inside Args struct:
#[arg(long, default_value = "isofrac", value_enum)]
pub filter_mode: FilterMode,
```

`FilterMode` flows through `RunConfig` → `print_predcluster_with_summary_multi` → `isofrac_with_summary`.

---

## Integration in `isofrac_with_summary`

The `longunder` computation is extended with a new branch:

```rust
let mut longunder = match filter_mode {
    FilterMode::Isofrac => {
        // existing formula (with combined_factor)
        if txs[k].exons.len() > 1 {
            (cmp_multicov <= 0.0 && cov < isofraclong * cmp_usedcov * combined_factor
                && cov < drop / error_perc)
                || cov < isofraclong * cmp_multicov * combined_factor
        } else {
            cov < isofraclong * cmp_usedcov * combined_factor
        }
    }
    FilterMode::Ml => {
        let features = MlFeatures::from_pair(&txs[k], &txs[first]);
        !ml_predict(&features)
    }
};
```

`MlFeatures::from_pair(k: &Transcript, dom: &Transcript) -> MlFeatures` is a constructor in `ml_filter.rs` that computes all 8 features. It uses `tx_exonic_len()` (already available in `transcript_filter.rs`) for `tlen_k` and `tlen_dom`. Zero-valued denominators produce 0.0 for the ratio (no panic). All existing rescue blocks (`RUSTLE_RELAX_LONGUNDER`, `RUSTLE_JCT_ISOFRAC_RESCUE`, etc.) remain and apply after this gate in both modes — they can still override a kill decision.

---

## Testing

- **Placeholder model test**: `rustle --filter-mode ml` with placeholder tree (`ml_predict` always returns true) should produce the same output as `rustle -f 0`. Verify with gffcompare.
- **After training**: benchmark `--filter-mode ml` on GGO_19 de novo. Target: equal or better F1 vs isofrac baseline (96.1% Sn / 91.0% Pr / 1742 chains). Report Sn, Pr, chain count.
- **Regression check**: `--filter-mode isofrac` output must be bit-identical to pre-change output (feature dump is env-gated, no behavior change in default mode).
- **Unit tests**: `MlFeatures::from_pair` with a synthetic dominant/minority pair; assert clamping behavior for zero-valued denominators.

---

## Workflow Summary

1. Implement `ml_filter.rs`, CLI flag, feature extraction, dump, integration — with placeholder `ml_predict`.
2. Run Step 1 (feature dump) on GGO_19.
3. Run `train_ml_filter.py` → paste tree into `ml_predict`.
4. Rebuild, benchmark.
5. Iterate on features or depth if result is unsatisfactory.
