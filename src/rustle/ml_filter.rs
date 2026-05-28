use std::io::{BufWriter, Write};
use std::sync::{Mutex, OnceLock};

use crate::path_extract::Transcript;

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
    /// where donor = exon[i].end + 1 and acceptor = exon[i+1].start (1-indexed inclusive).
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
    use crate::path_extract::Transcript;

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
