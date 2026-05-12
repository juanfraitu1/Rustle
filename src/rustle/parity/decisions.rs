//! Parity-decision log: write structured (JSONL) per-decision records during graph build,
//! transfrag construction, path extraction, and filtering. Enabled via env vars:
//!
//! - `RUSTLE_PARITY_LOG=/path/to/file.jsonl` — destination file (one event per line)
//! - `RUSTLE_PARITY_FILTER_CHROM=NC_073243.2` — restrict events to one chromosome (optional)
//! - `RUSTLE_PARITY_FILTER_RANGE=22149137-22155355` — restrict to coord range on the filtered chrom (optional)
//! - `RUSTLE_PARITY_FILTER_STEPS=junction_accept,node_create` — comma-separated whitelist (optional)
//!
//! Event format (each line is JSON):
//! ```
//! {"step":"junction_accept","tool":"rustle","chrom":"NC_073243.2","start":15400361,"end":15402062,"strand":"-","payload":{"support":12,"abundance":12.0}}
//! ```
//!
//! `step` and `tool` are mandatory. `chrom`/`start`/`end`/`strand` are conventional but optional;
//! when present they enable cheap range filtering at parse time. `payload` is a free-form JSON
//! object with step-specific fields.
//!
//! StringTie writes events to the same format (see `stringtie/parity_decisions.{h,cc}` once
//! instrumented). The diff tool at `tools/parity_decisions/diff.py` joins both files by
//! `(step, chrom, start, end, strand)` and reports divergence.

use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};

static WRITER: OnceLock<Mutex<Option<BufWriter<std::fs::File>>>> = OnceLock::new();
static ENABLED: AtomicBool = AtomicBool::new(false);
static FILTER: OnceLock<Filter> = OnceLock::new();

#[derive(Default, Clone)]
struct Filter {
    chrom: Option<String>,
    range: Option<(u64, u64)>,
    steps: Option<std::collections::HashSet<String>>,
}

#[inline]
pub fn is_enabled() -> bool {
    ENABLED.load(Ordering::Relaxed)
}

/// Initialize the parity decision logger from the standard env vars.
/// Idempotent — calling more than once is a no-op (uses OnceLock).
pub fn init_from_env() {
    let path = match std::env::var("RUSTLE_PARITY_LOG") {
        Ok(p) if !p.is_empty() => p,
        _ => return,
    };
    let file = match OpenOptions::new().create(true).write(true).truncate(true).open(&path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("[parity_decisions] failed to open {}: {}", path, e);
            return;
        }
    };
    let _ = WRITER.set(Mutex::new(Some(BufWriter::with_capacity(64 * 1024, file))));

    let mut f = Filter::default();
    if let Ok(c) = std::env::var("RUSTLE_PARITY_FILTER_CHROM") {
        if !c.is_empty() { f.chrom = Some(c); }
    }
    if let Ok(r) = std::env::var("RUSTLE_PARITY_FILTER_RANGE") {
        if let Some((s, e)) = r.split_once('-') {
            if let (Ok(s), Ok(e)) = (s.parse::<u64>(), e.parse::<u64>()) {
                f.range = Some((s, e));
            }
        }
    }
    if let Ok(s) = std::env::var("RUSTLE_PARITY_FILTER_STEPS") {
        if !s.is_empty() {
            f.steps = Some(s.split(',').map(|x| x.trim().to_string()).collect());
        }
    }
    eprintln!(
        "[parity_decisions] enabled: log={} chrom={:?} range={:?} steps={:?}",
        path, f.chrom, f.range, f.steps
    );
    let _ = FILTER.set(f);
    ENABLED.store(true, Ordering::Relaxed);

    // Header line (machine-readable: this is also a JSON object).
    let _ = emit_raw(r#"{"step":"_meta","tool":"rustle","payload":{"version":1}}"#);
}

fn passes_filter(step: &str, chrom: Option<&str>, start: u64, end: u64) -> bool {
    // _meta and other infra events always pass — they don't have meaningful coords.
    if step.starts_with('_') { return true; }
    let f = match FILTER.get() { Some(f) => f, None => return true };
    if let Some(only) = &f.steps {
        if !only.contains(step) { return false; }
    }
    // Chrom filter only applies when the event carries chrom info. Events without
    // chrom (e.g., Junction lacks chrom field in rustle) pass through — match by
    // coord range against the bundle's chrom is the user's responsibility.
    if let Some(c) = &f.chrom {
        if let Some(cc) = chrom {
            if cc != c { return false; }
        }
    }
    if let Some((rs, re_)) = f.range {
        // Skip range check if both coords are 0 (sentinel for "unknown coords").
        if !(start == 0 && end == 0) {
            if end < rs || start > re_ { return false; }
        }
    }
    true
}

fn emit_raw(line: &str) -> std::io::Result<()> {
    let lk = match WRITER.get() {
        Some(l) => l,
        None => return Ok(()),
    };
    let mut g = lk.lock().unwrap();
    if let Some(w) = g.as_mut() {
        writeln!(w, "{line}")?;
    }
    Ok(())
}

/// Emit a parity-log event. `payload_json` should be a serialized JSON object body
/// (without surrounding braces — they're added here). Use `pjson!` macro for ergonomics.
pub fn emit(step: &str, chrom: Option<&str>, start: u64, end: u64, strand: char, payload_json: &str) {
    if !is_enabled() { return; }
    if !passes_filter(step, chrom, start, end) { return; }
    let chrom_s = chrom.map(|c| format!(r#""chrom":"{}","#, c)).unwrap_or_default();
    let line = format!(
        r#"{{"step":"{step}","tool":"rustle",{chrom_s}"start":{start},"end":{end},"strand":"{strand}","payload":{{{payload_json}}}}}"#
    );
    let _ = emit_raw(&line);
}

/// Convenience macro: format a JSON object body inline.
/// Usage: `pjson!("k1": v1, "k2": v2)`
#[macro_export]
macro_rules! pjson {
    () => { String::new() };
    ($($k:literal : $v:expr),* $(,)?) => {{
        let mut s = String::new();
        let mut first = true;
        $(
            if !first { s.push(','); }
            first = false;
            s.push_str(&format!(r#""{}":{}"#, $k, serde_json::to_string(&$v).unwrap_or("null".into())));
        )*
        let _ = first;
        s
    }};
}

/// Flush + close the log on drop. Call at end of main if you want guaranteed flush.
pub fn flush() {
    if let Some(lk) = WRITER.get() {
        if let Ok(mut g) = lk.lock() {
            if let Some(w) = g.as_mut() {
                let _ = w.flush();
            }
        }
    }
}
