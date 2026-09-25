//! `as_table` — the genome-wide best-alignment-score pre-pass behind seeding loci with GOOD secondaries.
//!
//! One streaming scan of a whole BAM (every mapped, non-supplementary record, primary and secondary
//! alike) → one row per molecule:
//!
//! ```text
//! name<TAB>best_as<TAB>second_as<TAB>n_records<TAB>primary_chrom<TAB>primary_as
//! ```
//!
//! which is exactly what `RUSTLE_GTF_SECONDARY_AS_TABLE` consumes (`global_best_as` in
//! `denovo_assemble.rs` reads the first two columns). A secondary record is then admitted into the
//! locus-seeding read pool when `AS >= ratio × best_as` with the ratio measured GENOME-WIDE, not
//! within the region — the distinction that made the per-region ratio inert (register row 850) and
//! made the genome-wide one better on the held-out substrate (row 1060: gorilla F 0.030 → 0.089 at
//! precision 1.000; row 1100: ratio 0.98 is the narrowest pool that gives the whole gain).
//!
//! Replaces the `samtools view | awk | sort | awk` pipeline (`as_scan.sh`, 113 min on a 96 GB human
//! BAM) with one pass; rows are written in name order so two tables diff cleanly. Missing `AS` counts
//! as 0, as the shell pipeline did. Memory is one entry per molecule (a 21 M-molecule human BAM ≈ 2 GB).
use anyhow::Result;
use clap::Parser;
use std::io::Write;

#[derive(Parser, Debug)]
#[command(about = "Genome-wide best AS per molecule (one BAM scan) for RUSTLE_GTF_SECONDARY_AS_TABLE.")]
struct Args {
    /// Input BAM (any order; the whole file is read once).
    #[arg(long)]
    bam: String,
    /// Output TSV (name, best_as, second_as, n_records, primary_chrom, primary_as), sorted by name.
    #[arg(long)]
    out: String,
    /// BGZF decompression threads.
    #[arg(long, default_value_t = 4)]
    threads: usize,
}

#[derive(Clone, Copy)]
struct Mol {
    best: i32,
    second: i32,
    n: u32,
    /// Reference-sequence index of the PRIMARY record (`u32::MAX` until one is seen).
    primary_ref: u32,
    primary_as: i32,
}

fn record_as(record: &noodles_bam::Record) -> Option<i32> {
    use noodles_sam::alignment::record::data::field::{Tag, Value};
    for entry in noodles_sam::alignment::Record::data(record).iter() {
        let (tag, value) = entry.ok()?;
        if tag == Tag::ALIGNMENT_SCORE {
            return match value {
                Value::Int8(v) => Some(v as i32),
                Value::UInt8(v) => Some(v as i32),
                Value::Int16(v) => Some(v as i32),
                Value::UInt16(v) => Some(v as i32),
                Value::Int32(v) => Some(v),
                Value::UInt32(v) => Some(v as i32),
                _ => None,
            };
        }
    }
    None
}

fn main() -> Result<()> {
    let args = Args::parse();
    let t0 = std::time::Instant::now();
    let mut reader = rustle::bam::open_bam(&args.bam, args.threads.max(1))?;
    let header = reader.read_header()?;
    let contigs: Vec<String> = header.reference_sequences().keys().map(|k| k.to_string()).collect();
    let mut mols: std::collections::HashMap<String, Mol> = std::collections::HashMap::new();
    let mut n_records = 0u64;
    for result in reader.records() {
        let record = result?;
        let flags = record.flags();
        // `samtools view -F 2052`: drop unmapped and supplementary, keep primary + secondary
        if flags.is_unmapped() || flags.is_supplementary() {
            continue;
        }
        n_records += 1;
        let as_ = record_as(&record).unwrap_or(0);
        let name = match record.name() {
            Some(n) => n.to_string(),
            None => continue,
        };
        let m = mols.entry(name).or_insert(Mol { best: -1, second: -1, n: 0, primary_ref: u32::MAX, primary_as: 0 });
        m.n += 1;
        if as_ > m.best {
            m.second = m.best;
            m.best = as_;
        } else if as_ > m.second {
            m.second = as_;
        }
        if !flags.is_secondary() {
            m.primary_ref = record
                .reference_sequence_id()
                .and_then(|r| r.ok())
                .map(|id| id as u32)
                .unwrap_or(u32::MAX);
            m.primary_as = as_;
        }
        if n_records % 5_000_000 == 0 {
            eprintln!("[as-table] {n_records} records, {} molecules, {:.0} s", mols.len(), t0.elapsed().as_secs_f64());
        }
    }
    let mut names: Vec<&String> = mols.keys().collect();
    names.sort_unstable();
    let mut out = std::io::BufWriter::with_capacity(1 << 20, std::fs::File::create(&args.out)?);
    // provenance header (skipped by `global_best_as`): the driver reuses a table only when `bam=` is its BAM
    let bam_abs = std::fs::canonicalize(&args.bam).map(|p| p.display().to_string()).unwrap_or_else(|_| args.bam.clone());
    writeln!(out, "#as_table\tbam={bam_abs}\trecords={n_records}\tmolecules={}", mols.len())?;
    for name in names {
        let m = &mols[name];
        let (pchrom, pas) = if m.primary_ref == u32::MAX {
            ("NA".to_string(), "NA".to_string())
        } else {
            (contigs.get(m.primary_ref as usize).cloned().unwrap_or_else(|| "NA".to_string()), m.primary_as.to_string())
        };
        writeln!(out, "{name}\t{}\t{}\t{}\t{pchrom}\t{pas}", m.best, m.second, m.n)?;
    }
    out.flush()?;
    eprintln!(
        "[as-table] {n_records} records -> {} molecules -> {} in {:.0} s",
        mols.len(),
        args.out,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}
