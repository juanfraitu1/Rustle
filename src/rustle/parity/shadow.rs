//! Shadow parity logger for layer-by-layer StringTie comparisons.
//!
//! Enable with `RUSTLE_PARITY_SHADOW=1`.
//! Optional fail-fast: `RUSTLE_PARITY_SHADOW_FAIL_FAST=1`.
//! Optional output path: `RUSTLE_PARITY_SHADOW_TSV=/path/file.tsv`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};

use anyhow::Result;

static WRITER: OnceLock<Mutex<Option<BufWriter<File>>>> = OnceLock::new();
static ENABLED: AtomicBool = AtomicBool::new(false);
static FIRST_DIVERGENCE_REPORTED: AtomicBool = AtomicBool::new(false);

fn writer_cell() -> &'static Mutex<Option<BufWriter<File>>> {
    WRITER.get_or_init(|| Mutex::new(None))
}

pub fn enabled() -> bool {
    ENABLED.load(Ordering::Relaxed)
}

pub fn init(output_gtf: &Path) -> Result<()> {
    if std::env::var_os("RUSTLE_PARITY_SHADOW").is_none() {
        ENABLED.store(false, Ordering::Relaxed);
        *writer_cell().lock().expect("parity shadow lock poisoned") = None;
        return Ok(());
    }
    let path = std::env::var("RUSTLE_PARITY_SHADOW_TSV")
        .ok()
        .unwrap_or_else(|| {
            let parent = output_gtf.parent().unwrap_or_else(|| Path::new("."));
            let stem = output_gtf
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("rustle");
            parent
                .join(format!("{}.parity_shadow.tsv", stem))
                .to_string_lossy()
                .into_owned()
        });
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);
    writeln!(
        writer,
        "layer\tchrom\tstart\tend\tstrand\trustle_a\trustle_b\trustle_c\tshadow_a\tshadow_b\tshadow_c\tdiff\tnote"
    )?;
    *writer_cell().lock().expect("parity shadow lock poisoned") = Some(writer);
    ENABLED.store(true, Ordering::Relaxed);
    FIRST_DIVERGENCE_REPORTED.store(false, Ordering::Relaxed);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn emit_counts(
    layer: &str,
    chrom: &str,
    start: u64,
    end: u64,
    strand: char,
    rustle_a: i64,
    rustle_b: i64,
    rustle_c: i64,
    shadow_a: i64,
    shadow_b: i64,
    shadow_c: i64,
    note: &str,
) {
    if !enabled() {
        return;
    }
    let mut diff = 0u8;
    if shadow_a >= 0 && shadow_a != rustle_a {
        diff = 1;
    }
    if shadow_b >= 0 && shadow_b != rustle_b {
        diff = 1;
    }
    if shadow_c >= 0 && shadow_c != rustle_c {
        diff = 1;
    }

    if diff == 1 && !FIRST_DIVERGENCE_REPORTED.swap(true, Ordering::Relaxed) {
        eprintln!(
            "[PARITY_SHADOW] first divergence at {} {}:{}-{} strand={} (rustle={}/{}/{}, shadow={}/{}/{})",
            layer, chrom, start, end, strand, rustle_a, rustle_b, rustle_c, shadow_a, shadow_b, shadow_c
        );
        if std::env::var_os("RUSTLE_PARITY_SHADOW_FAIL_FAST").is_some() {
            panic!("RUSTLE_PARITY_SHADOW_FAIL_FAST: divergence at {}", layer);
        }
    }

    let mut guard = writer_cell().lock().expect("parity shadow lock poisoned");
    let Some(writer) = guard.as_mut() else {
        return;
    };
    let _ = writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        layer,
        chrom,
        start,
        end,
        strand,
        rustle_a,
        rustle_b,
        rustle_c,
        shadow_a,
        shadow_b,
        shadow_c,
        diff,
        note
    );
}

pub fn flush() {
    let mut guard = writer_cell().lock().expect("parity shadow lock poisoned");
    if let Some(writer) = guard.as_mut() {
        let _ = writer.flush();
    }
}
