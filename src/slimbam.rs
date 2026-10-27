//! slimbam: drop the secondaries no step reads (user 2026-10-27: a smaller BAM to iterate faster).
//!
//!   slimbam --bam IN.bam --out OUT.bam [--ratio 0.8] [--margin 10] [--threads 5]
//!
//! Keeps every primary, supplementary and unmapped record, and a secondary when, on AS OR on ms, its score >= min(ratio x P,
//! P - margin), P = the read's primary's score of the same tag. AS is what every consumer reads today (SEMISEC's cut
//! `--ratio 0.8`, O2's tie window within 10 AS of the primary, Eichler); ms is what minimap2 ranks the primary by (r1409: max ms
//! 99.3%, no intron cost), kept so a step can move to ms without the full BAM (user 2026-10-27). Exact ties stay. A secondary
//! whose read has no primary is kept; a missing tag counts 0 (as `semisec` for AS).
//! Two passes over the input (primary AS per read name, then the filter); records are copied as read, header unchanged.
//! The full BAM stays the archive: `semisec --ratio` below 0.8 needs it. Index the output with `samtools index`.
use anyhow::{Context, Result};
use clap::Parser;
use std::collections::HashMap;

#[derive(Parser)]
struct Args {
    #[arg(long)]
    bam: String,
    #[arg(long)]
    out: String,
    #[arg(long, default_value_t = 0.8)]
    ratio: f64,
    #[arg(long, default_value_t = 10)]
    margin: i32,
    #[arg(long, default_value_t = 5)]
    threads: usize,
}

/// Whether a secondary with score `sec` passes against its read's primary score `p` on one tag.
fn passes(sec: i32, p: i32, ratio: f64, margin: i32) -> bool {
    sec as f64 >= (ratio * p as f64).min((p - margin) as f64)
}

/// A secondary (AS, ms) survives against its read's primary (AS, ms) when either tag passes; no primary: kept.
fn keep_secondary(sec: (i32, i32), prim: Option<(i32, i32)>, ratio: f64, margin: i32) -> bool {
    prim.map_or(true, |p| passes(sec.0, p.0, ratio, margin) || passes(sec.1, p.1, ratio, margin))
}

const MS: noodles_sam::alignment::record::data::field::Tag = noodles_sam::alignment::record::data::field::Tag::new(b'm', b's');

fn scores(rec: &noodles_bam::Record) -> (i32, i32) {
    let ms = rustle::bam::aux_first(rec, MS, rustle::bam::value_as_i32).ok().flatten().unwrap_or(0);
    (rustle::bam::record_as(rec).unwrap_or(0), ms)
}

fn main() -> Result<()> {
    let a = Args::parse();
    let t0 = std::time::Instant::now();
    let mut prim: HashMap<Vec<u8>, (i32, i32)> = HashMap::new();
    {
        let mut r = rustle::bam::open_bam(&a.bam, a.threads)?;
        r.read_header()?;
        let mut rec = noodles_bam::Record::default();
        while r.read_record(&mut rec)? > 0 {
            let f = rec.flags();
            if f.is_unmapped() || f.is_secondary() || f.is_supplementary() {
                continue;
            }
            if let Some(n) = rec.name() {
                prim.insert(n.to_vec(), scores(&rec));
            }
        }
    }
    eprintln!("[slimbam] {} primary reads ({:.0} s)", prim.len(), t0.elapsed().as_secs_f64());
    let mut r = rustle::bam::open_bam(&a.bam, a.threads)?;
    let header = r.read_header()?;
    let file = std::fs::File::create(&a.out).with_context(|| a.out.clone())?;
    let workers = std::num::NonZeroUsize::new(a.threads.max(1)).unwrap();
    let mut w = noodles_bam::io::Writer::from(noodles_bgzf::MultithreadedWriter::with_worker_count(workers, file));
    w.write_header(&header)?;
    let (mut n_in, mut n_sec, mut n_sec_kept) = (0u64, 0u64, 0u64);
    let mut rec = noodles_bam::Record::default();
    while r.read_record(&mut rec)? > 0 {
        n_in += 1;
        if rec.flags().is_secondary() {
            n_sec += 1;
            let p = rec.name().and_then(|n| prim.get(n).copied());
            if !keep_secondary(scores(&rec), p, a.ratio, a.margin) {
                continue;
            }
            n_sec_kept += 1;
        }
        w.write_record(&header, &rec)?;
    }
    w.into_inner().finish()?;
    eprintln!(
        "[slimbam] records {n_in}; secondaries {n_sec}, kept {n_sec_kept}; written {} ({:.0} s)",
        n_in - n_sec + n_sec_kept,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::keep_secondary;

    #[test]
    fn keeps_the_union_of_the_semisec_cut_and_the_o2_tie_window_on_as_or_ms() {
        let p = Some((1000, 1000));
        assert!(keep_secondary((800, 0), p, 0.8, 10), "AS ratio 0.8 exactly");
        assert!(!keep_secondary((799, 799), p, 0.8, 10), "both tags below");
        assert!(keep_secondary((500, 800), p, 0.8, 10), "ms alone passes (a spliced copy losing on intron cost)");
        assert!(keep_secondary((1000, 1000), p, 0.8, 10), "exact tie");
        // a short read: 0.8 x 40 = 32 > 40 - 10 = 30, so the tie window decides
        assert!(keep_secondary((30, 0), Some((40, 40)), 0.8, 10));
        assert!(!keep_secondary((29, 29), Some((40, 40)), 0.8, 10));
        assert!(keep_secondary((0, 0), None, 0.8, 10), "no primary: kept");
    }
}
