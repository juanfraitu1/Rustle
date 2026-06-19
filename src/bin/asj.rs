//! Allele-specific junction (ASJ) scanner — quantify splice junctions whose USAGE depends on the allele a
//! molecule carries. Per gene region: call a balanced het anchor SNP, partition the primary reads by allele,
//! and test each junction (per-allele PSI + Fisher 2×2). Genome-wide BH-FDR + a |ΔPSI| floor select the ASJs.
//! Writes `<out>.asj.tsv`. A `.bai` next to the BAM makes the region reads fast.

use anyhow::{Context, Result};
use clap::Parser;
use std::io::Write;

use rustle::vg_family::allele_specific_junctions::{
    bh_fdr, scan_gene_asj, scan_gene_asj_multisnp, scan_gene_copy_specific_junctions, AsjParams,
};
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
    /// Phasing mode: `single` (one het anchor SNP), `multisnp` (2-means over ALL het SNPs → diploid
    /// haplotypes), or `copy` (PSVs → paralog COPIES; a call is a COPY-SPECIFIC junction).
    #[arg(long, default_value = "single")]
    mode: String,
}

/// One tested junction, flattened to (p, |ΔPSI|, output row) for genome-wide BH-FDR + emission.
struct Row {
    p: f64,
    dpsi: f64,
    line: String,
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
    let header = match args.mode.as_str() {
        "single" => "chrom\tanchor_pos\tax\tay\tdonor\tacceptor\tused_x\tspan_x\tused_y\tspan_y\tpsi_x\tpsi_y\tdpsi\tp\tq",
        "multisnp" | "copy" => "chrom\tn_snps\tphase_qual\tdonor\tacceptor\tused0\tspan0\tused1\tspan1\tpsi0\tpsi1\tdpsi\tp\tq",
        m => anyhow::bail!("unknown --mode {m} (single|multisnp|copy)"),
    };

    let mut rows: Vec<Row> = Vec::new(); // every tested junction, genome-wide
    for (chrom, lo, hi) in &regions {
        let reads = primary_aligned_reads_in_region(&args.bam, chrom, *lo, *hi)
            .with_context(|| format!("reading {chrom}:{lo}-{hi}"))?;
        match args.mode.as_str() {
            "single" => {
                for c in scan_gene_asj(&reads, *lo, *hi, &p) {
                    rows.push(Row {
                        p: c.p,
                        dpsi: c.dpsi,
                        line: format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3e}",
                            chrom, c.anchor_pos, c.ax as char, c.ay as char, c.donor, c.acceptor,
                            c.used_x, c.span_x, c.used_y, c.span_y, c.psi_x, c.psi_y, c.dpsi, c.p
                        ),
                    });
                }
            }
            _ => {
                let calls = if args.mode == "copy" {
                    scan_gene_copy_specific_junctions(&reads, *lo, *hi, &p)
                } else {
                    scan_gene_asj_multisnp(&reads, *lo, *hi, &p)
                };
                for c in calls {
                    rows.push(Row {
                        p: c.p,
                        dpsi: c.dpsi,
                        line: format!(
                            "{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3e}",
                            chrom, c.n_snps, c.phase_qual, c.donor, c.acceptor,
                            c.used0, c.span0, c.used1, c.span1, c.psi0, c.psi1, c.dpsi, c.p
                        ),
                    });
                }
            }
        }
    }

    // genome-wide BH-FDR over every tested junction, then select q<thresh & |ΔPSI|>=floor.
    let q = bh_fdr(&rows.iter().map(|r| r.p).collect::<Vec<_>>());
    let mut fh = std::fs::File::create(format!("{}.asj.tsv", args.out))?;
    writeln!(fh, "{header}")?;
    let mut n_asj = 0;
    for (i, r) in rows.iter().enumerate() {
        if q[i] < args.q && r.dpsi >= args.dpsi {
            n_asj += 1;
            writeln!(fh, "{}\t{:.3e}", r.line, q[i])?;
        }
    }
    eprintln!(
        "[asj] mode={} | {} genes scanned, {} junctions tested, {} {} (q<{} & |dPSI|>={}) -> {}.asj.tsv",
        args.mode, regions.len(), rows.len(), n_asj,
        if args.mode == "copy" { "copy-specific junctions" } else { "ASJ" },
        args.q, args.dpsi, args.out
    );
    Ok(())
}
