//! Allele-specific junction (ASJ) scanner — quantify splice junctions whose USAGE depends on the allele a
//! molecule carries. Per gene region: call a balanced het anchor SNP, partition the primary reads by allele,
//! and test each junction (per-allele PSI + Fisher 2×2). Genome-wide BH-FDR + a |ΔPSI| floor select the ASJs.
//! Writes `<out>.asj.tsv`. A `.bai` next to the BAM makes the region reads fast.

use anyhow::{Context, Result};
use clap::Parser;
use std::io::Write;

use rustle::vg_family::allele_specific_junctions::{bh_fdr, scan_gene_asj, AsjCall, AsjParams};
use rustle::vg_family::denovo_assemble::primary_aligned_reads_in_region;

#[derive(Parser, Debug)]
#[command(name = "asj", about = "Allele-specific junction scanner (per-molecule allele→junction linkage)")]
struct Args {
    /// Coordinate-sorted BAM (a `.bai` next to it enables the fast indexed region query).
    #[arg(long)]
    bam: String,
    /// A single gene region `chrom:start-end`.
    #[arg(long)]
    region: Option<String>,
    /// A regions FILE (one `chrom:start-end` per line) — the gene loci to scan.
    #[arg(long)]
    regions: Option<String>,
    /// Output prefix; writes `<out>.asj.tsv`.
    #[arg(long)]
    out: String,
    /// Effect-size floor: report a junction only if `|ΔPSI| >= dpsi`.
    #[arg(long, default_value_t = 0.30)]
    dpsi: f64,
    /// FDR threshold: report a junction only if its BH q-value `< q`.
    #[arg(long, default_value_t = 0.05)]
    q: f64,
}

fn parse_region(s: &str) -> Result<(String, u64, u64)> {
    let tok = s.split_whitespace().next().context("empty region")?;
    let (chrom, range) = tok.split_once(':').context("region must be chrom:start-end")?;
    let (lo, hi) = range.split_once('-').context("region must be chrom:start-end")?;
    Ok((chrom.to_string(), lo.replace(',', "").parse().context("bad start")?, hi.replace(',', "").parse().context("bad end")?))
}

fn main() -> Result<()> {
    let args = Args::parse();
    let regions: Vec<(String, u64, u64)> = match (&args.region, &args.regions) {
        (Some(r), None) => vec![parse_region(r)?],
        (None, Some(f)) => std::fs::read_to_string(f)
            .with_context(|| format!("reading {f}"))?
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(parse_region)
            .collect::<Result<_>>()?,
        _ => anyhow::bail!("provide exactly one of --region or --regions"),
    };

    let p = AsjParams::default();
    let mut calls: Vec<(String, AsjCall)> = Vec::new(); // all tested junctions, genome-wide
    for (chrom, lo, hi) in &regions {
        let reads = primary_aligned_reads_in_region(&args.bam, chrom, *lo, *hi)
            .with_context(|| format!("reading {chrom}:{lo}-{hi}"))?;
        for c in scan_gene_asj(&reads, *lo, *hi, &p) {
            calls.push((chrom.clone(), c));
        }
    }

    // genome-wide BH-FDR over every tested junction, then select q<thresh & |ΔPSI|>=floor.
    let q = bh_fdr(&calls.iter().map(|(_, c)| c.p).collect::<Vec<_>>());
    let mut fh = std::fs::File::create(format!("{}.asj.tsv", args.out))?;
    writeln!(fh, "chrom\tanchor_pos\tax\tay\tdonor\tacceptor\tused_x\tspan_x\tused_y\tspan_y\tpsi_x\tpsi_y\tdpsi\tp\tq")?;
    let mut n_asj = 0;
    for (i, (chrom, c)) in calls.iter().enumerate() {
        if q[i] < args.q && c.dpsi >= args.dpsi {
            n_asj += 1;
            writeln!(
                fh,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3e}\t{:.3e}",
                chrom, c.anchor_pos, c.ax as char, c.ay as char, c.donor, c.acceptor,
                c.used_x, c.span_x, c.used_y, c.span_y, c.psi_x, c.psi_y, c.dpsi, c.p, q[i]
            )?;
        }
    }
    eprintln!(
        "[asj] {} genes scanned, {} junctions tested, {} ASJ (q<{} & |dPSI|>={}) -> {}.asj.tsv",
        regions.len(), calls.len(), n_asj, args.q, args.dpsi, args.out
    );
    Ok(())
}
