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
use rustle::vg_family::denovo_assemble::{reads_in_region, tied_secondary_reads_in_region};
use rustle::vg_family::denovo_pipeline::{detect_and_assign, DenovoConfig, FallbackEdge};

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
    /// A single region to scan as `chrom:start-end`.
    #[arg(long)]
    region: Option<String>,
    /// A regions FILE (one `chrom:start-end` per line; extra whitespace-separated columns ignored) for a
    /// genome-wide sweep. Regions are grouped by contig so each contig is loaded ONCE.
    #[arg(long)]
    regions: Option<String>,
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
    /// poasta memory threshold (bp) for POA homology confirmation. A candidate family pair whose larger
    /// transcript exceeds this is confirmed via the linear-memory longest-common-substring FALLBACK instead of
    /// poasta (which OOMs on long sequences); those edges are recorded in `<out>.fallback.tsv`. Lower it (e.g.
    /// 8000) on dense/large-gene regions to keep poasta off the big operands. Default 20000 matches the python.
    #[arg(long, default_value_t = 20_000)]
    max_poa_len: usize,
    /// Recover COLLAPSED copies: feed AS-tied SECONDARY reads (a copy whose reads minimap2 flagged secondary
    /// because it picked a sibling as primary) into the rescue, so the starved copy can clear the support gate.
    /// Additive to the rescue only; default OFF (primary-only, byte-identical).
    #[arg(long, default_value_t = false)]
    recover_copies: bool,
    /// AS-tie ratio for --recover-copies: a secondary counts only if its AS >= ratio * the read's best AS
    /// (1.0 = exact tie; 0.98 admits a 2% margin). Guards against homology-shadow spillover.
    #[arg(long, default_value_t = 0.98)]
    as_ratio: f64,
}

fn status_str(s: AssignStatus) -> &'static str {
    match s {
        AssignStatus::Assigned => "assigned",
        AssignStatus::Ambiguous => "ambiguous",
        AssignStatus::Tied => "tied",
    }
}

fn parse_region(s: &str) -> Result<(String, u64, u64)> {
    let tok = s.split_whitespace().next().context("empty region")?;
    let (chrom, range) = tok.split_once(':').context("region must be chrom:start-end")?;
    let (lo_s, hi_s) = range.split_once('-').context("region must be chrom:start-end")?;
    Ok((chrom.to_string(), lo_s.parse().context("bad region start")?, hi_s.parse().context("bad region end")?))
}

/// One assignment-table row (resolved while the region's reads are in scope).
struct AssignRow {
    read_name: String,
    family_id: String,
    assigned_copy: usize,
    status: &'static str,
    n_decisive: usize,
    margin: f64,
}
/// One family-table row.
struct FamilyRow {
    family_id: String,
    chrom: String,
    n_copies: usize,
    n_reads: usize,
    psv_cols: usize,
    resolvable_psv: usize,
    resolvable_j: usize,
    junction_only: usize,
    assigned_j: usize,
    uniq_agree: usize,
    uniq: usize,
    collapsed_copies: usize,
    rescued_copies: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // collect the regions, then group by contig so each contig loads ONCE (memory-bounded sweep).
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
    let mut by_contig: std::collections::BTreeMap<String, Vec<(u64, u64)>> = std::collections::BTreeMap::new();
    for (c, lo, hi) in regions {
        by_contig.entry(c).or_default().push((lo, hi));
    }
    eprintln!(
        "[copy_assign] sweeping {} region(s) over {} contig(s)",
        by_contig.values().map(|v| v.len()).sum::<usize>(),
        by_contig.len()
    );

    let mut cfg = DenovoConfig::default();
    cfg.detect.len_cap = args.max_poa_len; // poasta memory threshold: above it, the bounded LCS fallback
    let params = AssignParams::default();
    let mut family_rows: Vec<FamilyRow> = Vec::new();
    let mut assign_rows: Vec<AssignRow> = Vec::new();
    let mut fallback_all: Vec<FallbackEdge> = Vec::new(); // family edges confirmed via the LCS fallback
    let mut gfam = 0usize; // global family counter (unique ids across regions)

    for (contig, ranges) in &by_contig {
        // load this contig once; it is dropped (freed) before the next contig.
        let contigs: HashSet<String> = std::iter::once(contig.clone()).collect();
        let genome = GenomeIndex::from_fasta_contigs(&args.fasta, &contigs)
            .with_context(|| format!("loading {} for {contig}", args.fasta))?;
        for &(lo, hi) in ranges {
            let (primary, bam_reads) = reads_in_region(&args.bam, contig, lo, hi, args.threads)
                .with_context(|| format!("reading {contig}:{lo}-{hi}"))?;
            let extra = if args.recover_copies {
                tied_secondary_reads_in_region(&args.bam, contig, lo, hi, args.as_ratio).unwrap_or_default()
            } else {
                Vec::new()
            };
            let (fams, fallback) =
                detect_and_assign(&primary, &bam_reads, &genome, &cfg, args.win, args.min_copies, &params, &extra);
            fallback_all.extend(fallback);
            for fa in &fams {
                let fid = format!("CAFAM{gfam}");
                gfam += 1;
                for (ri, a) in &fa.assignments {
                    assign_rows.push(AssignRow {
                        read_name: bam_reads[*ri].name.clone(),
                        family_id: fid.clone(),
                        assigned_copy: a.best_copy,
                        status: status_str(a.status),
                        n_decisive: a.n_decisive,
                        margin: a.log_lr_margin,
                    });
                }
                family_rows.push(FamilyRow {
                    family_id: fid,
                    chrom: fa.chrom.clone(),
                    n_copies: fa.n_copies,
                    n_reads: fa.n_reads,
                    psv_cols: fa.psv_cols,
                    resolvable_psv: fa.resolvable_psv,
                    resolvable_j: fa.resolvable_j,
                    junction_only: fa.junction_only,
                    assigned_j: fa.assigned_j,
                    uniq_agree: fa.uniq_agree,
                    uniq: fa.uniq,
                    collapsed_copies: fa.collapsed_copies,
                    rescued_copies: fa.rescued_copies,
                });
            }
            eprintln!("[copy_assign]   {contig}:{lo}-{hi}: {} mapped reads -> {} families", bam_reads.len(), fams.len());
        }
    }

    let mut fh = std::fs::File::create(format!("{}.families.tsv", args.out))?;
    writeln!(
        fh,
        "family_id\tchrom\tn_copies\trescued_copies\tcollapsed_copies\tn_reads\tpsv_cols\tresolvable_psv\tresolvable_j\tjunction_only\tassigned_j\tuniq_agree\tuniq"
    )?;
    for r in &family_rows {
        writeln!(
            fh,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.family_id, r.chrom, r.n_copies, r.rescued_copies, r.collapsed_copies, r.n_reads, r.psv_cols,
            r.resolvable_psv, r.resolvable_j, r.junction_only, r.assigned_j, r.uniq_agree, r.uniq
        )?;
    }
    let mut ah = std::fs::File::create(format!("{}.assignments.tsv", args.out))?;
    writeln!(ah, "read_name\tfamily_id\tassigned_copy\tstatus\tn_decisive\tmargin")?;
    for r in &assign_rows {
        writeln!(
            ah,
            "{}\t{}\t{}\t{}\t{}\t{:.3}",
            r.read_name, r.family_id, r.assigned_copy, r.status, r.n_decisive, r.margin
        )?;
    }

    // document every family edge confirmed via the large-sequence LCS fallback (poasta memory threshold
    // exceeded), so edges resting on the approximate metric are auditable. Only written when the fallback ran.
    if !fallback_all.is_empty() {
        let mut sh = std::fs::File::create(format!("{}.fallback.tsv", args.out))?;
        writeln!(sh, "chrom\ttid_a\tstart_a\tend_a\tlen_a\ttid_b\tstart_b\tend_b\tlen_b")?;
        for s in &fallback_all {
            writeln!(
                sh,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                s.chrom, s.tid_a, s.start_a, s.end_a, s.len_a, s.tid_b, s.start_b, s.end_b, s.len_b
            )?;
        }
        eprintln!(
            "[copy_assign] {} family edge(s) confirmed via the large-seq LCS fallback (transcript > --max-poa-len={}); wrote {}.fallback.tsv",
            fallback_all.len(),
            args.max_poa_len,
            args.out
        );
    }

    let (uniq, agree): (usize, usize) = family_rows.iter().fold((0, 0), |(u, g), f| (u + f.uniq, g + f.uniq_agree));
    eprintln!(
        "[copy_assign] {} families, {} read assignments",
        family_rows.len(),
        assign_rows.len()
    );
    if uniq > 0 {
        eprintln!(
            "[copy_assign] genome-wide silver-standard unique-mapper agreement: {agree}/{uniq} ({:.1}%)",
            100.0 * agree as f64 / uniq as f64
        );
    }
    eprintln!("[copy_assign] wrote {0}.families.tsv + {0}.assignments.tsv", args.out);
    Ok(())
}
