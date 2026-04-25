//! Per-site counters for hardstart/hardend writes. Set RUSTLE_HARDSTART_COUNT=1
//! to dump aggregated counts on program exit.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Mutex;

static SITES: Mutex<Vec<(&'static str, &'static AtomicU64)>> = Mutex::new(Vec::new());

#[inline]
pub fn bump(label: &'static str, counter: &'static AtomicU64) {
    counter.fetch_add(1, Ordering::Relaxed);
    // Register on first hit so dump knows which sites were active.
    let mut sites = SITES.lock().unwrap();
    if !sites.iter().any(|(l, _)| *l == label) {
        sites.push((label, counter));
    }
}

pub fn dump() {
    if std::env::var_os("RUSTLE_HARDSTART_COUNT").is_none() {
        return;
    }
    let sites = SITES.lock().unwrap();
    let mut rows: Vec<(&'static str, u64)> = sites
        .iter()
        .map(|(l, c)| (*l, c.load(Ordering::Relaxed)))
        .collect();
    rows.sort_by(|a, b| b.1.cmp(&a.1));
    eprintln!("HARDSTART_COUNTERS:");
    for (label, count) in rows {
        eprintln!("  {:7}  {}", count, label);
    }
}

/// Macro to declare and bump a per-site counter atomically.
#[macro_export]
macro_rules! bump_hs {
    ($label:expr) => {{
        static C: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        $crate::hard_counters::bump($label, &C);
    }};
}
