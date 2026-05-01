//! Filter a BAM by alignment-score (AS) margin between primary and secondary records of the same read.
//!
//! Usage: filter_bam_by_as <in.bam> <out.bam> --mode {strict|ratio} --threshold <T>
//!   strict mode: keep secondary if AS_primary - AS_self <= T (T = absolute units, default 10)
//!   ratio mode:  keep secondary if AS_self >= T * AS_primary  (T = fraction, default 0.9)
//!
//! Always keeps:
//!   * the primary record
//!   * supplementary (chimeric/split) records of the same read
//!   * unmapped records
//!
//! Two-pass: pass 1 builds read_name_hash -> primary_AS map; pass 2 streams the BAM
//! emitting filtered records to a coordinate-sorted output BAM.

use anyhow::{anyhow, Context, Result};
use noodles_bam as bam;
use noodles_sam::alignment::io::Write as SamAlignmentWrite;
use noodles_sam::alignment::record::data::field::Tag;
use noodles_sam::alignment::RecordBuf;
use std::collections::HashMap;
use std::env;
use std::fs::File;

fn fnv1a64(s: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in s {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

#[derive(Copy, Clone, Debug)]
enum Mode {
    Strict,
    Ratio,
}

fn read_as<R: noodles_sam::alignment::Record>(rec: &R) -> Option<i32> {
    let data = rec.data();
    for entry in data.iter() {
        let (tag, value) = entry.ok()?;
        if tag == Tag::ALIGNMENT_SCORE {
            if let noodles_sam::alignment::record::data::field::Value::Int8(v) = value { return Some(v as i32); }
            if let noodles_sam::alignment::record::data::field::Value::UInt8(v) = value { return Some(v as i32); }
            if let noodles_sam::alignment::record::data::field::Value::Int16(v) = value { return Some(v as i32); }
            if let noodles_sam::alignment::record::data::field::Value::UInt16(v) = value { return Some(v as i32); }
            if let noodles_sam::alignment::record::data::field::Value::Int32(v) = value { return Some(v); }
            if let noodles_sam::alignment::record::data::field::Value::UInt32(v) = value { return Some(v as i32); }
        }
    }
    None
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "usage: filter_bam_by_as <in.bam> <out.bam> [--mode strict|ratio] [--threshold T]"
        );
        std::process::exit(2);
    }
    let in_path = &args[1];
    let out_path = &args[2];

    let mut mode = Mode::Strict;
    let mut threshold_strict: i32 = 10;
    let mut threshold_ratio: f64 = 0.9;

    let mut i = 3;
    while i < args.len() {
        match args[i].as_str() {
            "--mode" => {
                i += 1;
                match args.get(i).map(String::as_str) {
                    Some("strict") => mode = Mode::Strict,
                    Some("ratio") => mode = Mode::Ratio,
                    other => return Err(anyhow!("--mode must be strict or ratio, got {:?}", other)),
                }
            }
            "--threshold" => {
                i += 1;
                let v = args.get(i).ok_or_else(|| anyhow!("--threshold needs value"))?;
                match mode {
                    Mode::Strict => threshold_strict = v.parse()?,
                    Mode::Ratio => threshold_ratio = v.parse()?,
                }
            }
            x => return Err(anyhow!("unknown arg {x}")),
        }
        i += 1;
    }
    eprintln!(
        "[filter_bam_by_as] mode={:?} strict_t={} ratio_t={} in={} out={}",
        mode, threshold_strict, threshold_ratio, in_path, out_path
    );

    // ── PASS 1 ────────────────────────────────────────────────────────────
    // Build read_name_hash → AS_primary.
    let mut reader = bam::io::reader::Builder::default()
        .build_from_path(in_path)
        .with_context(|| format!("opening {in_path}"))?;
    let _header = reader.read_header()?;
    let mut primary_as: HashMap<u64, i32> = HashMap::new();
    let mut n_total = 0u64;
    let mut n_primary = 0u64;
    let mut n_no_as = 0u64;
    let mut rec = bam::Record::default();
    loop {
        let n = reader.read_record(&mut rec)?;
        if n == 0 { break; }
        n_total += 1;
        let flags = rec.flags();
        if flags.is_unmapped() { continue; }
        if flags.is_secondary() || flags.is_supplementary() { continue; }
        n_primary += 1;
        let name = match rec.name() { Some(n) => n, None => continue };
        let h = fnv1a64(name.as_ref());
        match read_as(&rec) {
            Some(v) => { primary_as.insert(h, v); }
            None => { n_no_as += 1; }
        }
    }
    eprintln!(
        "[filter_bam_by_as] pass1: {} records, {} primary ({} missing AS); {} primary-AS entries",
        n_total, n_primary, n_no_as, primary_as.len()
    );
    if primary_as.is_empty() {
        return Err(anyhow!("no AS tags found on any primary record — refusing to write empty BAM"));
    }

    // ── PASS 2 ────────────────────────────────────────────────────────────
    let mut reader = bam::io::reader::Builder::default().build_from_path(in_path)?;
    let header = reader.read_header()?;

    let mut writer = bam::io::writer::Builder::default()
        .build_from_writer(File::create(out_path).with_context(|| format!("creating {out_path}"))?);
    writer.write_header(&header)?;

    let mut n_kept = 0u64;
    let mut n_dropped = 0u64;
    let mut n_kept_secondary = 0u64;
    let mut rec = bam::Record::default();
    loop {
        let n = reader.read_record(&mut rec)?;
        if n == 0 { break; }
        let flags = rec.flags();

        let keep = if flags.is_unmapped() || flags.is_supplementary() || !flags.is_secondary() {
            true
        } else {
            // Secondary: apply AS-margin filter
            let name = match rec.name() { Some(n) => n, None => { n_dropped += 1; continue; } };
            let h = fnv1a64(name.as_ref());
            let prim_as = match primary_as.get(&h) { Some(&v) => v, None => { n_dropped += 1; continue; } };
            let sec_as = match read_as(&rec) { Some(v) => v, None => { n_dropped += 1; continue; } };
            let k = match mode {
                Mode::Strict => (prim_as - sec_as) <= threshold_strict,
                Mode::Ratio  => prim_as > 0 && (sec_as as f64) >= threshold_ratio * (prim_as as f64),
            };
            if k { n_kept_secondary += 1; }
            k
        };

        if keep {
            // Convert to RecordBuf so we can cap PacBio HiFi qualities (>93)
            // that noodles' encoder rejects. Cap at 93 instead of stripping —
            // preserves most quality info, only loses signal at the absolute top.
            let mut rb = RecordBuf::try_from_alignment_record(&header, &rec)?;
            for q in rb.quality_scores_mut().as_mut().iter_mut() {
                if *q > 93 { *q = 93; }
            }
            writer.write_alignment_record(&header, &rb)?;
            n_kept += 1;
        } else {
            n_dropped += 1;
        }
    }
    eprintln!(
        "[filter_bam_by_as] pass2: kept {} ({} secondary), dropped {} secondary",
        n_kept, n_kept_secondary, n_dropped
    );
    Ok(())
}
