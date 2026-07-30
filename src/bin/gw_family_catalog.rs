//! gw_family_catalog — the GENOME-WIDE multi-copy-gene-family catalog (interest I / O1).
//!
//! DEFINITION: a family is a HOMOLOGY component — a γ-quasi-clique of the transcribed-homology graph E_r
//! spanning ≥2 physically-distinct loci — enforced by `refine_families_exon_sum` (asm20 id≥0.80,
//! cov-of-shorter≥0.50 + sensitive tier + ≥2-distinct-loci gate; default ON). The read-conflict graph is
//! NOT the family definition: it SEEDS the catalog (the loci among which reads are genuinely confused = where
//! copies must be co-resolved) and is the within-family copy-number oracle (χ_H). Homology decides
//! membership; conflict counts copies. (`--homology-primary` makes E_r the membership proposer directly,
//! in-engine via minimap2 — the cleanest form; the default seeds with the read-conflict graph, then the
//! homology gate enforces the definition.) This replaces the shipped `denovo_families.tsv` (built by an
//! arbitrary `core_recip≥0.13` similarity threshold that over-merges: DNFAM0 = 728 members chr1→chrY).
//! Writes `<out>.families.tsv` + `<out>.copies.tsv`.

use anyhow::Result;
use clap::Parser;
use std::io::Write;

use rustle::vg_family::denovo_pipeline::{
    detect_single_copy_baseline_genome_wide,
    detect_conflict_catalog_genome_wide, detect_conflict_catalog_genome_wide_xchrom,
    detect_homology_catalog_genome_wide, family_protein_coheres, homology_refine_params,
    refine_families_exon_sum, DenovoConfig, RefineParams,
};
use rustle::vg_family::family_detect::DenovoTranscript;

#[derive(Parser, Debug)]
#[command(about = "Genome-wide de-tie read-conflict multi-copy-family catalog (no similarity threshold).")]
struct Args {
    /// Coordinate-sorted BAM (a `.bai` next to it enables the fast per-contig region query). Required
    /// unless `--from-genome` is given.
    #[arg(long)]
    bam: Option<String>,
    /// GENOME-ONLY mode: discover SD families from the genome alone (no reads, no annotation). Argument is
    /// a BED of search windows (e.g. `bench/soto/80_fams.chr.bed`). Reps are duplicated genomic loci found by
    /// self-alignment; grouping uses the SAME `homology_blocks` core as `--homology-primary`. Mutually
    /// exclusive with `--bam`.
    #[arg(long)]
    from_genome: Option<String>,
    /// Genome FASTA (with a `.fai`).
    #[arg(long)]
    fasta: String,
    /// Output prefix; writes `<out>.families.tsv` + `<out>.copies.tsv`.
    #[arg(long)]
    out: String,
    /// Minimum copies for a family (real paralog PAIRS are common → default 2).
    #[arg(long, default_value_t = 2)]
    min_copies: usize,
    /// Co-located window (bp): a family's copies must cluster within this span.
    #[arg(long, default_value_t = 5_000_000)]
    win: u64,
    /// BAM-reading threads.
    #[arg(long, default_value_t = 4)]
    threads: usize,
    /// Capture CROSS-CHROMOSOME paralog families (RABL2A/B, DAZ-class) — a read confused between copies on
    /// different chromosomes forms a conflict edge; families are NOT restricted to one chromosome. Without
    /// this flag, only co-located same-chrom tandem arrays are emitted (the `--win`-bounded default).
    #[arg(long, default_value_t = false)]
    cross_chrom: bool,
    /// REFINE the raw conflict-graph catalog by exon-sum (FLNC) homology (requires minimap2 on PATH). A real
    /// multi-copy family must (i) have its copies MUTUALLY HOMOLOGOUS full-length (their spliced exon-sum
    /// sequences all-vs-all align, asm20 identity>=0.80 cov-of-shorter>=0.50) and (ii) span >=2 spatially-
    /// DISTINCT loci (disjoint genomic spans). This removes the Alu-SINE-bridged cross-chrom over-merges and
    /// the same-locus isoform/fragment contamination that read-conflict alone admits. Emits one row per
    /// distinct locus. Strongly recommended with `--cross-chrom`.
    #[arg(long, default_value_t = false)]
    refine: bool,
    /// DISABLE homology refinement and emit the RAW read-conflict catalog. Refinement is ON BY DEFAULT: the raw
    /// conflict graph admits repeat-bridge (unrelated genes sharing an intronic repeat) and large-gene
    /// mis-chain false positives (~9% of families genome-wide, `bench/GW_CATALOG_FP_AUDIT.md`); the homology +
    /// distinct-locus gate removes them. Use `--no-refine` for the raw catalog (e.g. reproducing the pre-fix
    /// numbers, or when minimap2 is unavailable). The `--refine` flag is now a no-op kept for back-compat.
    #[arg(long, default_value_t = false)]
    no_refine: bool,
    /// With `--refine`, align the GENOMIC span (introns INCLUDED) instead of the exon-sum. Captures intron
    /// divergence (separates copies that are exon-identical but intron-divergent), at the cost of being
    /// stricter on older paralogs whose introns have diverged. Default off (exon-sum).
    #[arg(long, default_value_t = false)]
    refine_introns: bool,
    /// Add the PROTEIN divergent tier (longest-ORF -> mmseqs, fident>=0.50): recovers synonymous-divergent
    /// CODING paralogs (RABL2B-class) that nucleotide seeds cannot reach. Batched into one mmseqs run.
    /// Additive only. Needs `mmseqs`. Default off (the sensitive nucleotide tier is already on). Applies in
    /// TWO places: with `--refine`, as an extra within-family membership edge; with `--homology-primary`,
    /// as a THIRD genome-wide E_r DEFINITION-edge tier (alongside asm20/sensitive) that groups coding
    /// paralogs past the ~0.65 nucleotide-identity floor where reads map all-primary.
    #[arg(long, default_value_t = false)]
    protein_tail: bool,
    /// With `--refine`, DISABLE the default sensitive nucleotide tier (`-k11 -w5`), reverting to the pure
    /// asm20 (id>=0.80) baseline. For reproducing the exact asm20-only catalog.
    #[arg(long, default_value_t = false)]
    no_sensitive: bool,
    /// POA-CORE COMPLETION (cross-chrom only): after the read-conflict families are built, attach
    /// loosely-related paralogs at NEW loci whose contiguous POA core to a family member clears `t_core` —
    /// reaching divergent copies the conflict graph (confusable-only, ~>=87% identity) misses. Bounded +
    /// seeded by the conflict families. Default off.
    #[arg(long, default_value_t = false)]
    complete_core: bool,
    /// HOMOLOGY-PRIMARY (E_r) family definition: build families from exon-sum nucleotide homology
    /// (gamma-quasi-clique), not the read-conflict graph. Captures ancient paralogs + K=0 collapses the
    /// conflict path misses. Ships alongside --cross-chrom.
    #[arg(long, default_value_t = false)]
    homology_primary: bool,
    /// E_r nucleotide identity floor (sensitive tier). Default from RefineParams (~0.60). `0.98` = Soto SD98 mode.
    #[arg(long)]
    min_identity: Option<f64>,
    /// With `--homology-primary`: compute the E_r family edge on each locus's GENOMIC SPAN instead of its
    /// assembled exon-sum. Fixes family FRAGMENTATION — copies of one family whose transcripts assemble
    /// DISJOINT exon subsets share almost no assembled sequence, so no edge forms and the family splits
    /// (measured: 75% of should-link pairs had no exon-sum alignment at all). RNA still decides WHICH loci
    /// are expressed; the reference supplies complete sequence for grouping. Needs no read depth/copy number.
    /// Default off (exon-sum) — existing catalogs are byte-identical unless this is set.
    #[arg(long, default_value_t = false)]
    homology_genomic_span: bool,
    /// Enumerate genomic copy number (famCN) per family via Liftoff-style genome projection (spec §7):
    /// project each family's best-supported consensus onto the genome (minimap2) to recover K=0 collapses
    /// the RNA read-conflict/homology path merges into one locus. Writes `<out>.famcn.tsv`. Defaults ON
    /// when `--min-identity 0.98` (Soto SD98 mode), since that mode targets exactly this near-identical
    /// segdup regime.
    #[arg(long, default_value_t = false)]
    enumerate_copies: bool,
    /// Re-admit near-identical families that collapse to <2 RNA loci as K0_COLLAPSED, reporting
    /// genome-projected copy NUMBER (needs DNA parCN for per-read resolution). Writes <out>.collapsed.tsv.
    /// Default off; when off, all existing output is byte-identical.
    #[arg(long, default_value_t = false)]
    collapse_enumerate: bool,
    /// Re-admit EXON-IDENTICAL (0-PSV) but heavily-EXPRESSED families that collapse to <2 RNA loci as a
    /// K0_COLLAPSED_EXPRESSED copy-number class (projection >=2 loci, each read-supported; needs DNA parCN for
    /// per-read resolution). Writes <out>.expressed_collapsed.tsv. Requires --homology-primary (the drop
    /// point lives on that path; a silent no-op without it). Default off; byte-identical when off.
    #[arg(long, default_value_t = false)]
    collapse_expressed: bool,

    /// Emit the single-copy baseline instead of the family catalog: `<out>.single_copy.tsv` (one row per
    /// single-copy chi(H)=1 locus) + `<out>.lambda_global.tsv` (the genome-wide median n_reads = lambda_global,
    /// the copy-number normalizer). A TABLE, not a GTF -- this is copy-number calibration, not an isoform catalog.
    #[arg(long, default_value_t = false)]
    single_copy_baseline: bool,

    /// Generalized projection: project EVERY resolved copy consensus (not one per family) at id>=0.98 +
    /// >=3-read support to localize members a single-consensus projection misses. Writes <out>.allproj.tsv
    /// (DNA-localized parCN, NOT added to copies.tsv). Requires --homology-primary. Default off (byte-identical).
    #[arg(long, default_value_t = false)]
    project_all_families: bool,

    /// DNA-family fallback (DNA edge oracle): for RNA-orphan expressed loci (0 homology edges -> no RNA
    /// family), project the transcript onto the genome at a LOOSER identity floor (--dna-family-min-identity,
    /// default 0.90 vs collapse-expressed's 0.98) to recover DIVERGENT genomic paralogs whose expressed
    /// sequence is non-homologous (the DNA-family != RNA-family case). Writes <out>.dna_family.tsv (copy
    /// NUMBER only; per-read resolution needs DNA parCN). Requires --homology-primary. Default off; byte-identical when off.
    #[arg(long, default_value_t = false)]
    dna_family_fallback: bool,
    /// Projection identity floor for --dna-family-fallback (default 0.90).
    #[arg(long)]
    dna_family_min_identity: Option<f64>,
}

/// Parse a BED of search windows (`chrom\tstart\tend`, extra columns ignored) for `--from-genome`.
fn read_windows_bed(path: &str) -> Result<Vec<(String, u64, u64)>> {
    let mut out = Vec::new();
    for line in std::fs::read_to_string(path)?.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 3 {
            out.push((f[0].to_string(), f[1].parse()?, f[2].parse()?));
        }
    }
    Ok(out)
}

/// Write the family catalog (`<out>.families.tsv` + `.copies.tsv` + `.copies.fa`). Shared by the RNA
/// catalog path and the DNA `--from-genome` path: both hand it the same `Vec<Vec<DenovoTranscript>>`
/// (each inner vec = a family's copies). Sorts families by their first copy's (chrom, start) and returns
/// the sorted vec so the RNA path's downstream sections (collapse/famcn/project-all) see the same order.
fn emit_catalog(
    out: &str,
    mut fams: Vec<Vec<DenovoTranscript>>,
    refine_params: &RefineParams,
) -> Result<Vec<Vec<DenovoTranscript>>> {
    use std::collections::BTreeSet;
    let mut fh = std::fs::File::create(format!("{out}.families.tsv"))?;
    let mut ch = std::fs::File::create(format!("{out}.copies.tsv"))?;
    // The exon-sum (spliced, FLNC-derived) sequence of every copy, in transcription orientation (for the
    // RNA path); for the DNA path it is the genomic locus sequence. Annotation-free family-validation substrate.
    let mut sh = std::fs::File::create(format!("{out}.copies.fa"))?;
    writeln!(fh, "family_id\tn_copies\tn_chroms\tchroms\tcross_chrom\tavg_reads\tprotein_coheres")?;
    writeln!(ch, "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads")?;

    // deterministic order: by the first copy's (chrom, start).
    fams.sort_by(|a, b| {
        let ka = a.iter().map(|c| (c.chrom.as_str(), c.start)).min();
        let kb = b.iter().map(|c| (c.chrom.as_str(), c.start)).min();
        ka.cmp(&kb)
    });
    // orthogonal protein-coherence QC flag; "NA" when --protein-tail is off (no mmseqs call).
    let coheres = family_protein_coheres(&fams, refine_params);
    let mut n_xchrom = 0usize;
    for (fi, copies) in fams.iter().enumerate() {
        let fid = format!("GWFAM{fi}");
        let chroms: BTreeSet<&str> = copies.iter().map(|c| c.chrom.as_str()).collect();
        let cross = chroms.len() > 1;
        if cross {
            n_xchrom += 1;
        }
        let avg_reads = if copies.is_empty() {
            0.0
        } else {
            copies.iter().map(|c| c.n_reads as f64).sum::<f64>() / copies.len() as f64
        };
        let coheres_str = match coheres[fi] {
            Some(true) => "true",
            Some(false) => "false",
            None => "NA",
        };
        writeln!(
            fh,
            "{fid}\t{}\t{}\t{}\t{}\t{:.1}\t{coheres_str}",
            copies.len(),
            chroms.len(),
            chroms.iter().cloned().collect::<Vec<_>>().join(","),
            cross,
            avg_reads
        )?;
        let mut sorted = copies.clone();
        sorted.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
        for (ci, c) in sorted.iter().enumerate() {
            writeln!(
                ch,
                "{fid}\t{ci}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                c.tid,
                c.chrom,
                c.start,
                c.end,
                c.introns.len() + 1,
                c.strand,
                c.n_reads
            )?;
            if !c.seq.is_empty() {
                writeln!(
                    sh,
                    ">{fid}|{ci}|{}:{}-{}|{}|nexon={}",
                    c.chrom,
                    c.start,
                    c.end,
                    c.strand,
                    c.introns.len() + 1
                )?;
                sh.write_all(&c.seq)?;
                writeln!(sh)?;
            }
        }
    }
    eprintln!(
        "[gw-catalog] wrote {} families ({} cross-chromosome) -> {out}.families.tsv + {out}.copies.tsv",
        fams.len(),
        n_xchrom,
    );
    Ok(fams)
}

fn main() -> Result<()> {
    let args = Args::parse();
    let project_all = (args.project_all_families
        || std::env::var("RUSTLE_PROJECT_ALL_FAMILIES").ok().as_deref() == Some("1"))
        && args.homology_primary;
    let mut cfg = DenovoConfig::from_env();
    cfg.complete_poa_core = args.complete_core;
    cfg.collapse_enumerate = args.collapse_enumerate || cfg.collapse_enumerate;
    cfg.collapse_expressed = args.collapse_expressed || cfg.collapse_expressed;
    cfg.dna_family_fallback = args.dna_family_fallback || cfg.dna_family_fallback;
    if let Some(x) = args.dna_family_min_identity { cfg.dna_family_min_identity = x; }

    // --from-genome: read-free/annotation-free DNA family catalog. Discovers duplicated genomic loci by
    // self-alignment, then groups them with the SAME homology_blocks core the RNA --homology-primary path
    // uses (via families_from_reps). Returns before any BAM-consuming path.
    if let Some(win_bed) = args.from_genome.as_deref() {
        use rustle::vg_family::from_genome::{genome_reps, GenomeRepParams};
        if args.bam.is_some() {
            anyhow::bail!("--from-genome and --bam are mutually exclusive (genome-only mode takes no reads)");
        }
        let mut refine_params = homology_refine_params(args.min_identity, args.threads);
        refine_params.protein_tail = args.protein_tail;
        // DNA-mode GROUPING knobs (anti over/under-merge levers). Env-only and applied ONLY on this branch,
        // so the RNA path stays byte-identical.
        //   RUSTLE_GENOME_MIN_COVERAGE — homology-edge coverage floor (fraction of the SHORTER locus the
        //     alignment must cover). Default 0.50. This is the documented repeat/duplicon-bridge defense:
        //     raising it rejects "two families share one duplicated block" edges that fuse distinct families.
        //   RUSTLE_GENOME_GAMMA — γ-quasi-clique density for family blocks. Default 0.20. Raising it demands
        //     denser within-family connectivity, splitting blocks held together by a few bridge edges.
        if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_COVERAGE") {
            if let Ok(x) = v.parse::<f64>() { refine_params.min_coverage = x; }
        }
        let gamma: f64 = std::env::var("RUSTLE_GENOME_GAMMA")
            .ok().and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.20);
        eprintln!(
            "[gw-catalog-genome] grouping: min_identity={:.2} min_coverage={:.2} gamma={:.2}",
            refine_params.min_identity, refine_params.min_coverage, gamma
        );
        let windows = read_windows_bed(win_bed)?;
        let mut gp = GenomeRepParams::from_env();
        gp.threads = args.threads;
        gp.minimap2 = refine_params.minimap2.clone();
        // `--min-identity` sets BOTH floors: the SD locus-discovery floor here AND the family-grouping floor
        // (via refine_params above). Env RUSTLE_GENOME_MIN_IDENTITY (already applied by from_env) is overridden
        // only when the flag is given, so `--from-genome --min-identity 0.98` = SD98 discovery + grouping.
        if let Some(mi) = args.min_identity { gp.min_identity = mi; }
        let reps = genome_reps(&args.fasta, &windows, &gp)?;
        eprintln!("[gw-catalog-genome] {} windows -> {} duplicated-locus reps", windows.len(), reps.len());
        let dna_fams = rustle::vg_family::denovo_pipeline::families_from_reps(
            reps, &refine_params, gamma, args.min_copies, 0,
        )?;
        emit_catalog(&args.out, dna_fams, &refine_params)?;
        return Ok(());
    }

    // BAM is required for every path below (single-copy baseline, conflict/homology catalogs, project-all).
    let bam = args
        .bam
        .as_deref()
        .ok_or_else(|| anyhow::anyhow!("--bam is required unless --from-genome is given"))?;

    if args.single_copy_baseline {
        use rustle::vg_family::single_copy::lambda_global;
        let loci = detect_single_copy_baseline_genome_wide(
            bam, &args.fasta, args.threads, args.win, args.min_copies, !args.no_refine, &cfg,
        )?;
        let mut sc = std::fs::File::create(format!("{}.single_copy.tsv", args.out))?;
        writeln!(sc, "chrom\tstart\tend\tstrand\tn_reads\tn_exons\tchi_h\tn_psv")?;
        for l in &loci {
            writeln!(sc, "{}\t{}\t{}\t{}\t{}\t{}\t1\t0", l.chrom, l.start, l.end, l.strand, l.n_reads, l.n_exons)?;
        }
        let lam = lambda_global(&loci);
        let lam_str = lam.map(|x| format!("{x}")).unwrap_or_else(|| "NA".into());
        let mut lf = std::fs::File::create(format!("{}.lambda_global.tsv", args.out))?;
        writeln!(lf, "lambda_global\tn_single_copy_loci")?;
        writeln!(lf, "{}\t{}", lam_str, loci.len())?;
        eprintln!(
            "[gw-catalog] single-copy baseline: {} loci -> lambda_global={} ({}.single_copy.tsv, {}.lambda_global.tsv)",
            loci.len(), lam_str, args.out, args.out
        );
        return Ok(());
    }
    // unify to `Vec<Vec<DenovoTranscript>>` (each = a family's copies) for a single emit path. The
    // cross-chrom path emits same-chromosome families (incl. inverted duplications and distant paralogs)
    // and cross-chromosome families together; the default path is the tight same-strand tandem-array view.
    let mut refine_params = homology_refine_params(args.min_identity, args.threads);
    // Same knob as the --from-genome path, honoured here too: the E_r coverage floor (fraction of the
    // SHORTER exon-sum the alignment must cover, default 0.50) is the duplicon-bridge defence, and it is the
    // criterion that decides borderline members. Exposed so a member's exclusion can be attributed to
    // coverage rather than inferred -- e.g. NPIPB12, whose 13 homology edges all reach >= 0.89 identity but
    // top out at 0.25 coverage (bench/soto/merge_quality_analysis.md §24).
    if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_COVERAGE") {
        if let Ok(x) = v.parse::<f64>() { refine_params.min_coverage = x; }
    }
    // `--protein-tail` also promotes protein homology to a genome-wide E_r DEFINITION edge in the
    // homology-primary catalog (in addition to feeding the `--refine` block below): it recovers coding
    // paralogs that have diverged past the nucleotide seeds' ~0.65 identity floor.
    refine_params.protein_tail = args.protein_tail;
    // `--homology-genomic-span`: E_r edges on the genomic span of each RNA-detected locus (anti-fragmentation).
    // Needs the genome path; harmless when the flag is off (field stays false, edges stay exon-sum).
    refine_params.homology_genomic_span = args.homology_genomic_span;
    if args.homology_genomic_span {
        refine_params.intron_fasta = Some(args.fasta.clone());
        if !args.homology_primary {
            eprintln!("[gw-catalog] note: --homology-genomic-span affects the E_r edge, which only drives family membership under --homology-primary");
        }
    }
    let (raw, collapsed, expressed, dna_families): (
        Vec<Vec<DenovoTranscript>>,
        Vec<rustle::vg_family::collapse_enumerate::CollapsedFamily>,
        Vec<rustle::vg_family::collapse_enumerate::ExpressedCollapsedFamily>,
        Vec<rustle::vg_family::collapse_enumerate::ExpressedCollapsedFamily>,
    ) = if args.homology_primary {
        detect_homology_catalog_genome_wide(
            bam, &args.fasta, args.threads, args.min_copies, &cfg, &refine_params, 0.20,
        )?
    } else if args.cross_chrom {
        (
            detect_conflict_catalog_genome_wide_xchrom(
                bam, &args.fasta, args.threads, args.min_copies, &cfg,
            )?,
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
    } else {
        let catalog = detect_conflict_catalog_genome_wide(
            bam, &args.fasta, args.threads, args.win, args.min_copies, &cfg,
        )?;
        (catalog.into_iter().map(|c| c.copies).collect(), Vec::new(), Vec::new(), Vec::new())
    };
    // Exon-sum (FLNC) homology + distinct-locus refinement (the principled membership criterion). ON BY DEFAULT
    // (removes repeat-bridge + large-gene mis-chain FPs); `--no-refine` opts out to the raw conflict catalog.
    let refine = !args.no_refine;
    let fams: Vec<Vec<DenovoTranscript>> = if refine {
        let n_raw = raw.len();
        let params = RefineParams {
            threads: args.threads,
            include_introns: args.refine_introns,
            // Always supply the genome: the additive genomic-span tier (recovers near-identical segdups whose
            // partial transcript models fail the exon-sum coverage floor) needs it even without --refine-introns.
            intron_fasta: Some(args.fasta.clone()),
            nucleotide_sensitive: !args.no_sensitive,
            protein_tail: args.protein_tail,
            // Inherit the two E_r thresholds from `refine_params` rather than taking struct defaults --
            // otherwise `--min-identity` and RUSTLE_GENOME_MIN_COVERAGE are silently ignored on this path,
            // which is the code path `--refine` actually runs.
            min_identity: refine_params.min_identity,
            min_coverage: refine_params.min_coverage,
            ..Default::default()
        };
        let refined = refine_families_exon_sum(raw, &params, None, cfg.conflict.min_reads)?;
        eprintln!(
            "[gw-catalog] {}{}{} refine: {} raw families -> {} refined (homology component AND >= 2 distinct loci)",
            if args.refine_introns { "genomic(intron-inclusive)" } else { "exon-sum" },
            if args.no_sensitive { "" } else { " +sensitive" },
            if args.protein_tail { " +protein" } else { "" },
            n_raw,
            refined.len()
        );
        refined
    } else {
        raw
    };

    // Sort + write the family catalog (families.tsv/copies.tsv/copies.fa). Shared by the RNA path
    // and the DNA --from-genome path; returns the sorted families for the downstream sections.
    let fams = emit_catalog(&args.out, fams, &refine_params)?;

    // K=0-collapsed families re-admitted by Task 3 (`--collapse-enumerate`): copy-NUMBER only, kept fully
    // separate from families.tsv/copies.tsv so the OFF path (flag off or nothing collapsed) writes no file
    // and byte-identical-OFF holds.
    if cfg.collapse_enumerate && !collapsed.is_empty() {
        let mut cf = std::fs::File::create(format!("{}.collapsed.tsv", args.out))?;
        writeln!(cf, "family_id\tchrom\tstart\tend\tfamCN\tn_alt_reads\talt_frac\tstatus\tprojection_loci")?;
        for (i, fam) in collapsed.iter().enumerate() {
            writeln!(cf, "{}", rustle::vg_family::collapse_enumerate::format_collapsed_row(&format!("GWFAMc{i}"), fam))?;
        }
        eprintln!("[gw-catalog] wrote {} K=0-collapsed families -> {}.collapsed.tsv", collapsed.len(), args.out);
    }

    // K0_COLLAPSED_EXPRESSED families re-admitted by Task 3 (`--collapse-expressed`): same isolation as
    // `--collapse-enumerate` above -- OFF path (flag off or nothing expressed-collapsed) writes no file.
    if cfg.collapse_expressed && !expressed.is_empty() {
        let mut ef = std::fs::File::create(format!("{}.expressed_collapsed.tsv", args.out))?;
        writeln!(ef, "family_id\tchrom\tstart\tend\tfamCN\tmin_locus_reads\tstatus\tprojection_loci")?;
        for (i, fam) in expressed.iter().enumerate() {
            writeln!(ef, "{}", rustle::vg_family::collapse_enumerate::format_expressed_collapsed_row(&format!("GWFAMe{i}"), fam))?;
        }
        eprintln!("[gw-catalog] collapse-expressed: {} K0_COLLAPSED_EXPRESSED families -> {}.expressed_collapsed.tsv", expressed.len(), args.out);
    }

    // DNA-family fallback (`--dna-family-fallback`): RNA-orphan loci recovered by the DNA edge oracle. Same
    // isolation contract -- OFF path (flag off or nothing recovered) writes no file (byte-identical).
    if cfg.dna_family_fallback && !dna_families.is_empty() {
        let mut df = std::fs::File::create(format!("{}.dna_family.tsv", args.out))?;
        writeln!(df, "family_id\tchrom\tstart\tend\tfamCN\tmin_locus_reads\tstatus\tprojection_loci")?;
        for (i, fam) in dna_families.iter().enumerate() {
            writeln!(df, "{}", rustle::vg_family::collapse_enumerate::format_dna_family_row(&format!("GWFAMdna{i}"), fam))?;
        }
        eprintln!("[gw-catalog] dna-family-fallback: {} DNA_FAMILY_RNA_NONHOMOLOGOUS loci -> {}.dna_family.tsv", dna_families.len(), args.out);
    }

    // famCN / totalCN via genome projection (spec §7) -- the VG copy-number leg: LAND the family variation
    // graph's consensus PATH on the genome and count near-identical landing sites. A family's RNA-observed
    // copies are a LOWER bound on famCN when copies collapse onto one locus (K=0). Project each family's
    // best-supported consensus (most-read-supported copy's exon-sum) back onto the genome in ONE batched
    // minimap2 pass (one genome-index load for ALL families, not one per family -- avoids re-indexing the
    // whole genome hundreds/thousands of times) and count additional disjoint loci beyond the already-
    // known copy loci, bucketed by identity: totalCN at the >=0.80 asm20 floor, famCN at the >=0.98 Soto
    // SD98 near-identical floor.
    // totalCN cov>=0.80 (batch floor): full-length divergent copies without partial-domain fragment
    // inflation; famCN keeps cov>=0.90 (near-identical, Soto SD98). A coverage sweep on real families
    // showed cov>=0.50 inflates totalCN with partial/domain-fragment hits (GWFAM18: 16 vs 12=truth at
    // cov>=0.70), while cov>=0.90 is too strict and drops divergent full-length copies (GSTM 4 vs 19
    // truth); cov>=0.80 is the sweet spot (verified: GSTM 20 vs 19, RABL2 6 vs 5, GWFAM18 12 vs 11,
    // GWFAM21 22 vs 22).
    let enumerate = (args.enumerate_copies || args.min_identity == Some(0.98)) && args.homology_primary;
    if enumerate {
        use std::collections::HashMap;
        let consensuses: Vec<(String, Vec<u8>)> = fams
            .iter()
            .enumerate()
            .map(|(fi, copies)| {
                let cons = copies.iter().max_by_key(|c| c.n_reads).map(|c| c.seq.clone()).unwrap_or_default();
                (format!("GWFAM{fi}"), cons)
            })
            .collect();
        let known: HashMap<String, Vec<(String, u64, u64)>> = fams
            .iter()
            .enumerate()
            .map(|(fi, copies)| {
                let loci: Vec<(String, u64, u64)> = copies.iter().map(|c| (c.chrom.clone(), c.start, c.end)).collect();
                (format!("GWFAM{fi}"), loci)
            })
            .collect();
        let proj_by_fam = match rustle::vg_family::genome_projection::project_families_batch(
            &consensuses, &args.fasta, &known, 0.80, 0.80, &refine_params.minimap2, args.threads,
        ) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("[gw-catalog] batch famCN projection failed ({e}); famCN/totalCN fall back to n_rna_copies");
                HashMap::new()
            }
        };
        let mut ff = std::fs::File::create(format!("{}.famcn.tsv", args.out))?;
        writeln!(ff, "family_id\tn_rna_copies\tfamCN\ttotalCN\tprojection_loci")?;
        for (fi, copies) in fams.iter().enumerate() {
            let fid = format!("GWFAM{fi}");
            let n_rna = copies.len();
            let proj = proj_by_fam.get(&fid).cloned().unwrap_or_default();
            // totalCN: divergent copies admitted (id>=0.80; cov>=0.80 already enforced by the batch call
            // above -- full-length divergent copies without partial/domain-fragment inflation).
            let n_total_loci = proj.iter().filter(|p| p.identity >= 0.80).count();
            // famCN: Soto SD98 near-identical FULL-LENGTH copies only (id>=0.98 AND cov>=0.90) -- excludes
            // fragment/partial hits that would otherwise inflate the Soto-comparable metric. cov>=0.90 is
            // a subset of the returned cov>=0.80 loci.
            let n_fam_loci = proj.iter().filter(|p| p.identity >= 0.98 && p.cov >= 0.90).count();
            let total_cn = n_rna + n_total_loci;
            let fam_cn = n_rna + n_fam_loci;
            // projection_loci lists all totalCN-contributing loci (id>=0.80, cov>=0.80).
            let loci = proj
                .iter()
                .filter(|p| p.identity >= 0.80)
                .map(|p| format!("{}:{}-{}@{:.4}", p.chrom, p.start, p.end, p.identity))
                .collect::<Vec<_>>()
                .join(";");
            writeln!(ff, "{fid}\t{n_rna}\t{fam_cn}\t{total_cn}\t{loci}")?;
        }
        eprintln!("[gw-catalog] famCN/totalCN batch projection -> {}.famcn.tsv", args.out);
    }

    // Generalized projection (`--project-all-families`): project EVERY resolved copy's consensus (not just
    // the best-supported one per family) to localize members a single-consensus projection misses. OPTIONAL,
    // additive, DNA-localized parCN leg -- dedup'd (reciprocal-overlap >=0.50) + read-support gated (>=3
    // primary reads over the locus), written to its own file so the RNA-split catalog (families.tsv/
    // copies.tsv) and the famCN/totalCN batch projection (famcn.tsv) are untouched.
    if project_all {
        use rustle::vg_family::project_all::{CopyIn, all_copy_consensuses, known_from_fams, dedup_overlapping, overlaps_any, format_allproj_row};
        // Build (fid, copies) with the SAME fid the catalog uses.
        let fam_copies: Vec<(String, Vec<CopyIn>)> = fams.iter().enumerate().map(|(fi, copies)| {
            (format!("GWFAM{fi}"), copies.iter().map(|c| CopyIn { seq: c.seq.clone(), chrom: c.chrom.clone(), start: c.start, end: c.end }).collect())
        }).collect();
        let consensuses = all_copy_consensuses(&fam_copies);
        let known = known_from_fams(&fam_copies);
        let proj = rustle::vg_family::genome_projection::project_families_batch(
            &consensuses, &args.fasta, &known, 0.98, 0.90, &refine_params.minimap2, args.threads,
        ).unwrap_or_default();
        let spans_by_fam: std::collections::HashMap<&String, &Vec<(String,u64,u64)>> = known.iter().collect();
        let mut rows: Vec<String> = Vec::new();
        for (fid, locs) in &proj {
            for l in dedup_overlapping(locs.clone()) {
                // read-support gate: >=3 primary reads over the locus
                let n_support = rustle::vg_family::denovo_assemble::reads_in_region(bam, &l.chrom, l.start, l.end, args.threads)
                    .map(|(p, _)| p.len()).unwrap_or(0);
                if n_support < 3 { continue; }
                let overlaps = spans_by_fam.get(fid).map(|s| overlaps_any(&l.chrom, l.start, l.end, s)).unwrap_or(false);
                rows.push(format_allproj_row(fid, &l, n_support, overlaps));
            }
        }
        rows.sort(); // `proj` is a HashMap (random iteration order); sort so allproj.tsv is reproducible run-to-run
        if !rows.is_empty() {
            let mut af = std::fs::File::create(format!("{}.allproj.tsv", args.out))?;
            writeln!(af, "family_id\tchrom\tstart\tend\tidentity\tn_support_reads\toverlaps_existing_copy")?;
            for r in &rows { writeln!(af, "{r}")?; }
            eprintln!("[gw-catalog] project-all-families: {} loci (id>=0.98, >=3 reads) -> {}.allproj.tsv", rows.len(), args.out);
        }
    }
    Ok(())
}
