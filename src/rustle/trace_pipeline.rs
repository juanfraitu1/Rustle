//! Stage-level pipeline trace for diff-based bisection against a reference tool.
//!
//! Activated by setting `RUSTLE_TRACE_PIPELINE=/path/to/dir`. When set, selected
//! pipeline stages append canonical TSV rows to files under that directory:
//!
//! - `bundles.tsv`:     chrom, start, end, strand, n_reads
//! - `graph_edges.tsv`: chrom, bundle_start, bundle_end, strand, from_start, from_end, to_start, to_end
//!
//! Rows are NOT sorted as they are written. Pipe through `sort -k1,1 -k2,2n`
//! (or similar) before diffing.

use std::fs::{create_dir_all, OpenOptions};
use std::io::Write;
use std::path::PathBuf;
use std::sync::{Mutex, OnceLock};

fn trace_dir() -> Option<&'static PathBuf> {
    static DIR: OnceLock<Option<PathBuf>> = OnceLock::new();
    DIR.get_or_init(|| {
        let dir = std::env::var_os("RUSTLE_TRACE_PIPELINE")?;
        let path = PathBuf::from(dir);
        if create_dir_all(&path).is_err() {
            return None;
        }
        Some(path)
    })
    .as_ref()
}

static BUNDLES_LOCK: Mutex<()> = Mutex::new(());
static EDGES_LOCK: Mutex<()> = Mutex::new(());

fn append(filename: &str, rows: &str, lock: &Mutex<()>) {
    let Some(dir) = trace_dir() else {
        return;
    };
    let _g = lock.lock().ok();
    let path = dir.join(filename);
    if let Ok(mut f) = OpenOptions::new().create(true).append(true).open(&path) {
        let _ = f.write_all(rows.as_bytes());
    }
}

pub fn active() -> bool {
    trace_dir().is_some()
}

pub fn dump_bundle(chrom: &str, start: u64, end: u64, strand: char, n_reads: usize) {
    if !active() {
        return;
    }
    let row = format!("{}\t{}\t{}\t{}\t{}\n", chrom, start, end, strand, n_reads);
    append("bundles.tsv", &row, &BUNDLES_LOCK);
}

pub fn dump_graph_edge(
    chrom: &str,
    bundle_start: u64,
    bundle_end: u64,
    strand: char,
    from_start: u64,
    from_end: u64,
    to_start: u64,
    to_end: u64,
) {
    if !active() {
        return;
    }
    let row = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        chrom, bundle_start, bundle_end, strand, from_start, from_end, to_start, to_end
    );
    append("graph_edges.tsv", &row, &EDGES_LOCK);
}
