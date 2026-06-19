//! De-novo multi-copy family DETECTION + per-read copy ASSIGNMENT CLI.
//!
//! Annotation-free, read-coherence pipeline: scan a region of a coordinate-sorted BAM, assemble de-novo
//! transcripts, detect co-located paralog families, and assign each read — including the hard multimappers
//! minimap2 leaves at MAPQ 0 — to a specific copy via PSV bases + copy-specific junctions.
//!
//! Writes `<out>.families.tsv` (per-family roster + two-pass + silver-standard stats) and
//! `<out>.assignments.tsv` (per-read copy assignment). A `.bai` next to the BAM makes the region read fast.

use anyhow::{Context, Result};
use clap::Parser;
use std::collections::HashSet;
use std::io::Write;

use rustle::genome::GenomeIndex;
use rustle::vg_family::copy_assign::{AssignParams, AssignStatus};
use rustle::vg_family::denovo_assemble::reads_in_region;
use rustle::vg_family::denovo_pipeline::{detect_and_assign, DenovoConfig};

#[derive(Parser, Debug)]
#[command(
    name = "copy_assign",
    about = "De-novo multi-copy family detection + per-read copy assignment (PSV + copy-specific junctions)"
)]
struct Args {
    /// Coordinate-sorted BAM (a `.bai` next to it enables the fast indexed region query).
    #[arg(long)]
    bam: String,
    /// Genome FASTA (with a `.fai` so only the needed contig is loaded).
    #[arg(long)]
    fasta: String,
    /// Region to scan as `chrom:start-end` (a co-located cluster region).
    #[arg(long)]
    region: String,
    /// Output prefix; writes `<out>.families.tsv` and `<out>.assignments.tsv`.
    #[arg(long)]
    out: String,
    /// Minimum copies for a co-located family.
    #[arg(long, default_value_t = 3)]
    min_copies: usize,
    /// Co-located window (bp): copies must cluster within this span.
    #[arg(long, default_value_t = 5_000_000)]
    win: u64,
    /// BAM-reading threads.
    #[arg(long, default_value_t = 4)]
    threads: usize,
}

fn status_str(s: AssignStatus) -> &'static str {
    match s {
        AssignStatus::Assigned => "assigned",
        AssignStatus::Ambiguous => "ambiguous",
        AssignStatus::Tied => "tied",
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    let (chrom, range) = args.region.split_once(':').context("region must be chrom:start-end")?;
    let (lo_s, hi_s) = range.split_once('-').context("region must be chrom:start-end")?;
    let (lo, hi): (u64, u64) = (lo_s.parse().context("bad region start")?, hi_s.parse().context("bad region end")?);

    eprintln!("[copy_assign] region {chrom}:{lo}-{hi}");
    let (primary, bam_reads) =
        reads_in_region(&args.bam, chrom, lo, hi, args.threads).with_context(|| format!("reading {}", args.bam))?;
    eprintln!("[copy_assign] {} primary, {} mapped reads", primary.len(), bam_reads.len());

    let contigs: HashSet<String> = std::iter::once(chrom.to_string()).collect();
    let genome =
        GenomeIndex::from_fasta_contigs(&args.fasta, &contigs).with_context(|| format!("loading {}", args.fasta))?;

    let fams = detect_and_assign(
        &primary,
        &bam_reads,
        &genome,
        &DenovoConfig::default(),
        args.win,
        args.min_copies,
        &AssignParams::default(),
    );
    eprintln!("[copy_assign] {} co-located families", fams.len());

    let mut fh = std::fs::File::create(format!("{}.families.tsv", args.out))?;
    writeln!(
        fh,
        "family_id\tchrom\tn_copies\tn_reads\tpsv_cols\tresolvable_psv\tresolvable_j\tjunction_only\tassigned_j\tuniq_agree\tuniq"
    )?;
    for fa in &fams {
        writeln!(
            fh,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            fa.family_id, fa.chrom, fa.n_copies, fa.n_reads, fa.psv_cols, fa.resolvable_psv, fa.resolvable_j,
            fa.junction_only, fa.assigned_j, fa.uniq_agree, fa.uniq
        )?;
    }

    let mut ah = std::fs::File::create(format!("{}.assignments.tsv", args.out))?;
    writeln!(ah, "read_name\tfamily_id\tassigned_copy\tstatus\tn_decisive\tmargin")?;
    for fa in &fams {
        for (ri, a) in &fa.assignments {
            writeln!(
                ah,
                "{}\t{}\t{}\t{}\t{}\t{:.3}",
                bam_reads[*ri].name, fa.family_id, a.best_copy, status_str(a.status), a.n_decisive, a.log_lr_margin
            )?;
        }
    }

    let (uniq, agree): (usize, usize) = fams.iter().fold((0, 0), |(u, g), f| (u + f.uniq, g + f.uniq_agree));
    if uniq > 0 {
        eprintln!(
            "[copy_assign] silver-standard unique-mapper agreement: {agree}/{uniq} ({:.1}%)",
            100.0 * agree as f64 / uniq as f64
        );
    }
    eprintln!("[copy_assign] wrote {0}.families.tsv + {0}.assignments.tsv", args.out);
    Ok(())
}
