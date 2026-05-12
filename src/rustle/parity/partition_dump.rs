//! Canonical **partition geometry** dump for cross-tool parity (StringTie vs Rustle).
//!
//! Enable with `RUSTLE_PARITY_PARTITION_TSV=/path/partition.tsv`. Each bundle produces one line:
//! `partition_geometry_v1`, `rustle`, chrom, start, end, strand, signature.
//!
//! The signature is strand-agnostic: sorted `start-end` chains (segments within a chain joined
//! by `+`, chains joined by `|`). StringTie's `PARITY_PARTITION_TSV` walks **each** `CBundle`
//! linked list separately, sorts those chain strings, and joins with `|`. Use
//! [`encode_partition_signature_sub_bundle_walk`] so Rustle emits one chain per
//! [`crate::bundle_builder::SubBundleResult`] (same granularity as StringTie's per-bundle rows).

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::{Mutex, OnceLock};

use anyhow::Result;

use crate::bundle_builder::SubBundleResult;

static WRITER: OnceLock<Mutex<Option<BufWriter<File>>>> = OnceLock::new();

fn writer_cell() -> &'static Mutex<Option<BufWriter<File>>> {
    WRITER.get_or_init(|| Mutex::new(None))
}

pub fn enabled() -> bool {
    std::env::var_os("RUSTLE_PARITY_PARTITION_TSV").is_some()
}

pub fn init() -> Result<()> {
    let path = match std::env::var_os("RUSTLE_PARITY_PARTITION_TSV") {
        Some(p) => p,
        None => {
            *writer_cell().lock().expect("parity partition lock poisoned") = None;
            return Ok(());
        }
    };
    let path = Path::new(&path);
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "kind\tsource\tchrom\tstart\tend\tstrand\tsignature"
    )?;
    *writer_cell().lock().expect("parity partition lock poisoned") = Some(w);
    Ok(())
}

/// Encode color-split bundlenode geometry (intervals only) as a canonical signature string.
///
/// Segment endpoints follow StringTie's `parity_partition_emit_stringtie`: each bundlenode
/// interval is written as `(start, end - 1)` so the string matches `PARITY_PARTITION_TSV`
/// (Rustle bundlenode `end` is exclusive / one-past vs that emitted coordinate).
pub fn encode_partition_signature(components: &[Vec<(u64, u64)>]) -> String {
    let mut chains: Vec<String> = components
        .iter()
        .map(|seg| {
            seg.iter()
                .map(|(s, e)| {
                    let e_out = e.saturating_sub(1);
                    format!("{}-{}", s, e_out)
                })
                .collect::<Vec<_>>()
                .join("+")
        })
        .collect();
    chains.sort_unstable();
    chains.join("|")
}

/// One chain per sub-bundle linked list (`bnode_head` / `next`), matching StringTie's
/// `parity_partition_emit_stringtie` (one sorted chain string per `CBundle`).
pub fn encode_partition_signature_sub_bundle_walk(subs: &[SubBundleResult]) -> String {
    let mut components: Vec<Vec<(u64, u64)>> = Vec::new();
    for sb in subs {
        let mut chain: Vec<(u64, u64)> = Vec::new();
        let mut cur = sb.bnode_head.as_ref();
        while let Some(n) = cur {
            chain.push((n.start, n.end));
            cur = n.next.as_deref();
        }
        if !chain.is_empty() {
            components.push(chain);
        }
    }
    encode_partition_signature(&components)
}

pub fn emit_rustle_row(chrom: &str, start: u64, end: u64, strand: char, signature: &str) {
    if !enabled() {
        return;
    }
    let mut guard = writer_cell().lock().expect("parity partition lock poisoned");
    let Some(w) = guard.as_mut() else {
        return;
    };
    let _ = writeln!(
        w,
        "partition_geometry_v1\trustle\t{}\t{}\t{}\t{}\t{}",
        chrom, start, end, strand, signature
    );
}

pub fn flush() {
    let mut guard = writer_cell().lock().expect("parity partition lock poisoned");
    if let Some(w) = guard.as_mut() {
        let _ = w.flush();
    }
}
