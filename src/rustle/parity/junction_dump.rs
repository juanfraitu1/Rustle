//! Multi-stage **junction-set** dump for StringTie vs Rustle parity.
//!
//! Enable with `RUSTLE_PARITY_JUNCTION_TSV=/path/junctions.tsv`. Each row is one bundle key
//! (`chrom`, `start`, `end`, `strand`) at one **stage** of `build_graphs`, so you can compare
//! against StringTie even when that tool’s export corresponds to an earlier internal pass
//! (e.g. before synthetic transfrags, coverage trimming, or aggressive splice snapping).
//!
//! # Stages (Rustle)
//!
//! | `stage` | Meaning |
//! |---------|---------|
//! | `read_union` | Sorted union of per-read `junctions` after `repair_reads_after_junction_quality` |
//! | `pre_canonical` | Junctions from `junction_stats_corr_final` that pass the same support gate as `filter_junctions` **before** `apply_junction_filters_and_canonicalize` |
//! | `graph_input` | `filter_junctions` output — splice graph input |
//! | `cgroup_feed` | `good_junctions_set` passed into cgroup / 3-strand bundle builder |
//!
//! The `junctions` column is `donor-acceptor` pairs sorted lexicographically, joined with `|`.
//!
//! # StringTie mirror (`PARITY_JUNCTION_TSV` in `stringtie/rlink.cpp`)
//!
//! Same TSV columns (`kind`, `source`, `stage`, `chrom`, `start`, `end`, `strand`, `junctions`).
//! Compare with `scripts/compare_junction_parity.py` using `--ref-stage` / `--cand-stage` (e.g.
//! Rustle `cgroup_feed` vs StringTie `aggregate_graph_ready_after_read_pass`). Coordinate
//! convention for `start`/`end` matches `PARITY_PARTITION_TSV` (normalized bundle key).

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::{Mutex, OnceLock};

use anyhow::Result;

use crate::types::{Junction, JunctionStats};

static WRITER: OnceLock<Mutex<Option<BufWriter<File>>>> = OnceLock::new();

fn writer_cell() -> &'static Mutex<Option<BufWriter<File>>> {
    WRITER.get_or_init(|| Mutex::new(None))
}

pub fn enabled() -> bool {
    std::env::var_os("RUSTLE_PARITY_JUNCTION_TSV").is_some()
}

pub fn init() -> Result<()> {
    let path = match std::env::var_os("RUSTLE_PARITY_JUNCTION_TSV") {
        Some(p) => p,
        None => {
            *writer_cell().lock().expect("parity junction lock poisoned") = None;
            return Ok(());
        }
    };
    let path = Path::new(&path);
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "kind\tsource\tstage\tchrom\tstart\tend\tstrand\tjunctions"
    )?;
    *writer_cell().lock().expect("parity junction lock poisoned") = Some(w);
    Ok(())
}

/// Encode junctions as sorted `donor-acceptor` tokens joined by `|`.
pub fn encode_junction_set_sorted(js: impl Iterator<Item = Junction>) -> String {
    // Parity normalization: Rustle's internal acceptor coordinate is one base before the
    // StringTie PARITY_JUNCTION_TSV acceptor convention. Shift by +1 in the dump so stage
    // comparisons are algorithmic (set membership), not coordinate-system artifacts.
    let mut v: Vec<(u64, u64)> = js
        .map(|j| (j.donor, j.acceptor.saturating_add(1)))
        .collect();
    v.sort_unstable();
    v.into_iter()
        .map(|(d, a)| format!("{}-{}", d, a))
        .collect::<Vec<_>>()
        .join("|")
}

/// Same support gate as [`crate::junctions::filter_junctions`], applied to pre-canonical stats.
pub fn encode_stats_graphlike(stats: &JunctionStats, min_reads: f64) -> String {
    encode_junction_set_sorted(
        stats
            .iter()
            .filter(|(_, s)| {
                s.strand != Some(0) && s.nreads_good >= min_reads && s.mrcount > 0.0
            })
            .map(|(j, _)| *j),
    )
}

pub fn emit_row(chrom: &str, start: u64, end: u64, strand: char, stage: &str, junctions: &str) {
    if !enabled() {
        return;
    }
    let mut guard = writer_cell().lock().expect("parity junction lock poisoned");
    let Some(w) = guard.as_mut() else {
        return;
    };
    let _ = writeln!(
        w,
        "junction_set_v1\trustle\t{}\t{}\t{}\t{}\t{}\t{}",
        stage, chrom, start, end, strand, junctions
    );
}

pub fn flush() {
    let mut guard = writer_cell().lock().expect("parity junction lock poisoned");
    if let Some(w) = guard.as_mut() {
        let _ = w.flush();
    }
}
