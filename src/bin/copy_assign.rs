//! De-novo multi-copy family DETECTION + per-read copy ASSIGNMENT CLI.
//!
//! Annotation-free, read-coherence pipeline: scan a region of a coordinate-sorted BAM, assemble de-novo
//! transcripts, detect co-located paralog families, and assign each read — including the hard multimappers
//! minimap2 leaves at MAPQ 0 — to a specific copy via PSV bases + copy-specific junctions.
//!
//! Writes `<out>.families.tsv` (per-family roster + two-pass + unique-mapper agreement stats) and
//! `<out>.assignments.tsv` (per-read copy assignment). A `.bai` next to the BAM makes the region read fast.

use anyhow::{Context, Result};
use clap::Parser;
use lru::LruCache;
use std::collections::HashSet;
use std::io::Write;
use std::num::NonZeroUsize;
use std::sync::{Arc, Mutex};

use rayon::prelude::*;
use rustle::genome::GenomeIndex;
use rustle::vg_family::absent_copy::DnaNeedsRecord;
use rustle::vg_family::copy_assign::{AssignParams, AssignStatus};
use rustle::vg_family::em_copy_assign::em_assign_family;
use rustle::vg_family::denovo_assemble::{
    assemble_gate, pass1_skeletons, reads_in_region, tied_secondary_reads_in_region, BamIndexCache, BamRead,
};
use rustle::vg_family::denovo_pipeline::{
    catalog_overlaps, detect_and_assign, DenovoConfig, FallbackEdge, FamilyAssignment, OverlapKind,
};
use rustle::vg_family::family_detect::collapse_loci_groups;
use rustle::vg_family::read_conflict::{as_evidence, AsEvidence};
use rustle::vg_family::readonly_copy_number::{chi_h_with_junctions, depth_cn};

/// One assembled isoform (FLAIR-style intron-chain collapse), kept for the optional `--gtf` emit. `gene_tid`
/// is the locus this isoform collapses into (shared-junction gene); a family copy is its own gene, so a
/// family-copy gene_tid matches a `copy_tid` in the assignment output and is tagged `multicopy` in the GTF.
struct TranscriptRec {
    tid: String,
    gene_tid: String,
    chrom: String,
    start: u64,
    end: u64,
    strand: char,
    introns: Vec<(u64, u64)>,
    n_reads: u32,
}

/// The expensive, independent per-region work (BAM read + `detect_and_assign`, which holds the dominant
/// poasta alignment). Computed possibly in parallel across regions (any contig), then drained SERIALLY in the
/// original region order so `CAFAM` ids and every output row stay byte-identical to the serial sweep.
///
/// The heavy read SEQUENCES (`BamRead`) are dropped inside the worker — the output stage needs only the read
/// NAMES (everything else lives in `fams`) — so collecting all regions' results out-of-order is lightweight.
struct RegionWork {
    contig: String,
    lo: u64,
    hi: u64,
    read_names: Vec<String>,
    /// Alignment-score evidence per read, parallel to `read_names`. Reported, never decisive — see
    /// `read_conflict::AsEvidence` for why raw AS is length-confounded and `de` makes the call.
    as_ev: Vec<AsEvidence>,
    n_mapped: usize,
    fams: Vec<FamilyAssignment>,
    fallback: Vec<FallbackEdge>,
    dna_needs: Vec<DnaNeedsRecord>,
    transcripts: Vec<TranscriptRec>, // FLAIR-style isoforms for the --gtf emit (empty unless --gtf)
}

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
    /// Minimum copies for a co-located family. Two-copy homologous families are the majority and were
    /// invisible to assignment at the old default of 3; lowering it to 2 changes default family detection
    /// on its own, independently of `--homology-primary`.
    #[arg(long, default_value_t = 2)]
    min_copies: usize,
    /// Co-located window (bp): copies must cluster within this span.
    #[arg(long, default_value_t = 5_000_000)]
    win: u64,
    /// BAM-reading threads.
    #[arg(long, default_value_t = 4)]
    threads: usize,
    /// REGION-PARALLEL sweep: process this many regions (independent families) concurrently, ACROSS contigs.
    /// Each family is an independent unit (`detect_and_assign` is pure; `BamIndexCache`/genome are read-only),
    /// so the heavy per-family poasta alignments of N families — even on different chromosomes — overlap on N
    /// cores; a bounded LRU cache of loaded contig genomes (capacity ≈ N) avoids reloads and caps resident
    /// chromosomes. Output is collected and `CAFAM` ids assigned in the SAME (serial) order afterward, so the
    /// result is BYTE-IDENTICAL to the serial run. Peak memory ≈ N regions' reads + ≈N contig genomes (raise N
    /// for speed on a many-core box, lower it if memory-bound; the documented genome-wide OOM is why this is
    /// opt-in). `1` (default) = the exact serial path. The pool is sized to N and shared with the inner
    /// per-copy alignment parallelism. Speedup ceiling = the single heaviest family (already internally
    /// parallel); on a full genome-wide sweep this approaches ~N×.
    #[arg(long, default_value_t = 1)]
    region_threads: usize,
    /// FLAIR-LIKE ASSEMBLY emit: also write `<out>.gtf` — every de-novo isoform in the swept regions as a
    /// transcript+exon GTF (IGV-loadable), grouped into genes by shared junctions. Multi-copy family copies
    /// appear as separate genes at their own loci, tagged `family_id`/`copy_index`/`multicopy "true"`; all
    /// other loci are ordinary isoforms (`multicopy "false"`). Annotation-free (intron-chain collapse + the
    /// canonical-junction gate; no short-read junction correction). Independent of the assignment outputs;
    /// default off. Pair with `bench/igv_tracks.py` for the copy-coloured reads.
    #[arg(long, default_value_t = false)]
    gtf: bool,
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
    /// Also dump the per-read PSV GENOTYPE MATRIX — `<out>.psv_reads.tsv` (each read's base at every PSV column
    /// + its assignment), `<out>.psv_copies.tsv` (each copy's PSV alleles), `<out>.psv_cols.tsv` (column →
    /// genome position). The raw per-molecule evidence behind each assignment, for the proof visualization.
    #[arg(long, default_value_t = false)]
    dump_psv: bool,
    /// Decisive-margin τ: the minimum log-likelihood-ratio over the runner-up copy to call a read ASSIGNED
    /// (else AMBIGUOUS); the single calibrated knob that replaces the vote-count (min_psv, margin) integers.
    /// τ = ln((1−p)/p) for a target per-read misassignment p (τ=6.9 default ≈ p 1e-3, the Eichler AS≥10
    /// analog; τ=2.0 ≈ p 0.12). The identifiability gate (n_decisive≥1) is independent of τ and always applied.
    #[arg(long, default_value_t = 6.9)]
    margin: f64,
    /// Per-base PSV error rate used in the likelihood when a read carries no per-base quality (HiFi ~0.003).
    #[arg(long, default_value_t = 0.003)]
    error_rate: f64,
    /// IsoCon significance level α for the DEFAULT gate: assign iff the per-read certificate p < α/(n−1)
    /// (Bonferroni over the n−1 competitors) AND the read is the strict MLE. α is the FAMILY-WIDE
    /// misassignment rate over assigned reads (1e-3 ≈ the τ=6.9 / Eichler AS≥10 precision point). Ignored
    /// when `--margin-gate` is set.
    #[arg(long, default_value_t = 1e-3)]
    alpha: f64,
    /// Use the LEGACY τ-margin gate (the `--margin` knob) instead of the IsoCon significance gate. For
    /// reproducing pre-significance-gate numbers and the gate A/B comparison.
    #[arg(long, default_value_t = false)]
    margin_gate: bool,
    /// Disable the RNA-editing filter (Clair3-RNA): by default, A↔G PSV columns showing within-copy
    /// heterogeneity are flagged as A-to-I editing sites and downweighted in the certificate so an edited
    /// base cannot fake copy-support. This reverts to trusting every PSV column uniformly.
    #[arg(long, default_value_t = false)]
    no_editing_filter: bool,
    /// εⱼ used for an editing-flagged PSV column in the certificate (the rate a base shows the other allele
    /// by editing rather than sequencing error). Default 0.2.
    #[arg(long, default_value_t = 0.2)]
    edit_rate: f64,
    /// IsoCon-style iterative copy pruning: after assignment, repeatedly merge copies that have no read with
    /// significant evidence distinguishing them from their nearest neighbor, reassigning reads until all
    /// surviving copies are defensible. Default OFF (byte-identical baseline).
    #[arg(long, default_value_t = false)]
    iterative_prune: bool,
    /// Emit FACULTATIVE long-read PHASING (dependency-free; no external phaser, no CNN). Phases reads into
    /// copy-HAPLOTYPES from their linked PSV evidence — the N-copy generalisation of read-backed phasing
    /// (min-path-cover over the PSV graph). Writes `<out>.phase_blocks.tsv` (one PHASE SET per family),
    /// `<out>.phased_haplotypes.tsv` (each haplotype's `pos:allele` variant string — the phased alleles),
    /// and `<out>.phased_reads.tsv` (read → haplotype HAPLOTAG with the decisive-margin confidence;
    /// haplotype = -1 when the read is ambiguous/tied = unphased). Block ≙ PS tag, haplotype ≙ HP tag.
    #[arg(long, default_value_t = false)]
    phase: bool,

    /// Skip the POA homology DIAGNOSTIC pass. It is the dominant per-region cost (the poasta all-pairs
    /// alignment over candidate rep pairs — ~85% of wall-clock on dense families) but is purely diagnostic:
    /// families come from the de-tie conflict graph, so the emitted families/assignments/abundance are
    /// BYTE-IDENTICAL with or without it. Only the `.fallback.tsv` report and the POA edge counts in the log
    /// are lost. STRONGLY recommended for genome-wide sweeps (measured ~6.8× faster on the heaviest family).
    /// Equivalent to setting `RUSTLE_SKIP_POA_DIAGNOSTIC=1`.
    #[arg(long, default_value_t = false)]
    skip_poa_diagnostic: bool,

    /// Discover reference-ABSENT (collapsed) copies from reads and re-thread the abstain pool against them
    /// (two-stage freeze; default OFF = byte-identical). Candidates failing the admission gate are written
    /// to `<out>.dna_needs.tsv`.
    #[arg(long, default_value_t = false)]
    absent_copies: bool,
    /// Emit `<out>.posterior.tsv`: per read, the soft per-copy POSTERIOR and the consistent ZONE (the genomic
    /// region of the copies it is compatible with) — the Bayesian complement to the hard assign/abstain, so an
    /// unassignable (Tied) read is localized to a zone with a distribution instead of a bare flag. The prior is
    /// uniform by default; set `RUSTLE_POSTERIOR_PRIOR=abundance` to weight by the EM copy abundance. Default off.
    #[arg(long, default_value_t = false)]
    posterior: bool,

    /// Run the EM soft-relaxation (Vollger 2019 PSV correlation-clustering, maximum-likelihood soft
    /// version) over each family's PSV evidence and emit `<out>.em.tsv` (per-read soft posterior +
    /// K-frontier label) and `<out>.em_abundance.tsv` (per-copy recovered abundance). This is a
    /// PSV-only reduction of the gate likelihood (editing-filtered; junctions and per-base quality
    /// are NOT used), so per-read labels may differ from the hard `.assignments.tsv` gate on
    /// junction/quality-resolvable reads; the abundance estimate is robust to these. `.em_abundance.tsv`
    /// is the convergent, gate-likelihood EM (the estimator the consistency theorem describes); it
    /// uses `error_rate` and a convergence gate, so its `pi_hat` can differ from the legacy
    /// `.quant.tsv` (`soft_quantify_em`, fixed error 0.01 / 100 iters) -- prefer `.em_abundance.tsv`
    /// for the theorem's estimator. Additive: leaves `.assignments.tsv`/`.families.tsv`/`.quant.tsv`
    /// byte-identical. Default off.
    #[arg(long, default_value_t = false)]
    em: bool,
    /// Max E/M sweeps for `--em`.
    #[arg(long, default_value_t = 500)]
    em_max_iter: usize,
    /// EM convergence tolerance (absolute+relative on the observed-data log-likelihood) for `--em`.
    #[arg(long, default_value_t = 1e-6)]
    em_eps: f64,

    /// Reference-free RNA single-copy expression floor (lambda_global): the genome-wide median
    /// n_reads over single-copy transcripts, precomputed by
    /// `bench/rna_copy_number_depth.py::global_single_copy_anchor` -- an RNA-only quantity, NOT
    /// genomic. When given, enables the `depth_cn` column of `<out>.famcn_readonly.tsv` (else
    /// `NA`); this file is ALWAYS written (additive; independent of every other output).
    #[arg(long)]
    lambda_global: Option<f64>,

    /// VG re-align supplement (Task 5, REPORT-ONLY): for every co-located family, re-align each
    /// poor-fit/candidate read (low MAPQ, heavy clipping, or high divergence — `vg_realign::is_candidate`)
    /// to the family's copy-paths and record the decision (`reassigned` / `rejected` / `novel-candidate`)
    /// to `<out>.vg_realign.tsv`. Does NOT feed corrections back into the EM/PSV assignment and does NOT
    /// admit novel candidates into the copy set — those remain a separate follow-up. Default off; leaves
    /// every other output byte-identical.
    #[arg(long, default_value_t = false)]
    vg_realign: bool,

    /// Define family MEMBERSHIP by E_r transcript homology instead of the E_c read-conflict graph. The
    /// conflict graph links two copies only when reads map ambiguously between them, so a copy whose reads
    /// all map uniquely is dropped from its family and its reads come back `tied` — not because they are
    /// unassignable, but because their true copy was never admitted. Conflict, PSVs, and chi(H) remain
    /// within-family. Admitting a dropped copy enlarges the copy set, so the Bonferroni certificate
    /// alpha/(K-1) tightens and existing assignments shift; this is why the mode is opt-in. Requires
    /// minimap2 (honors RUSTLE_MINIMAP2); aborts rather than falling back to the conflict graph.
    #[arg(long, default_value_t = false)]
    homology_primary: bool,

    /// Keep unspliced readthrough transcripts as candidate copies. A single-exon de-novo transcript that
    /// engulfs >= 5 distinct splice junctions (each with >= 2 reads) is intronic pileup / unspliced pre-mRNA,
    /// not an mRNA, and is dropped by default. Validated on 15 such transcripts (minimum 14 engulfed
    /// junctions) against 260 expressed intronless genes (maximum 4), including the EEF1A1 retrocopy whose
    /// spliced parent cross-maps onto it. Pass this to disable the filter and reproduce the old behaviour.
    #[arg(long, default_value_t = false)]
    keep_readthrough: bool,

    /// EXPERIMENTAL, OFF BY DEFAULT. Admit a single-rep locus whose reads are ambiguously placed (MAPQ 0) at a
    /// rate incompatible with a unique locus as a multi-copy family with `n_copies = chi(H)`, reads certified
    /// tied, no per-copy sequence materialised.
    ///
    /// ⚠ The instrument detects unresolvable PARALOGY, not collapse. It fires on EEF1A1, whose MAPQ-0 reads
    /// align to processed pseudogenes on other chromosomes, and reports chi(H) = 7 for a one-copy locus. A copy
    /// genuinely absent from the reference would pile reads on at HIGH mapq, giving depth excess and no
    /// ambiguity -- which is why SDA detects collapses by depth. Do not use for copy number.
    #[arg(long, default_value_t = false)]
    collapse_gate: bool,

    /// Background per-read ambiguity rate for the collapse test (fraction of PRIMARY reads at MAPQ 0
    /// genome-wide). Must be a genome-wide quantity: a region-local estimate is degenerate, since in a
    /// collapsed window the reads outside the assembled rep are precisely the ambiguous ones. Default is the
    /// value measured on GGO_mm.bam (5785 / 4404440 = 0.001313). Recompute per sample with:
    ///   `echo $(( $(samtools view -c -F 2308 b.bam) - $(samtools view -c -F 2308 -q 1 b.bam) ))`
    #[arg(long)]
    eps_amb: Option<f64>,

    /// Apply the assemble gate's `min_reads` per ISOFORM (the pre-fix behaviour) instead of per LOCUS.
    /// Diagnostic: use to isolate the effect of junction-incidence pooling.
    #[arg(long, default_value_t = false)]
    no_pool_locus_support: bool,
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

/// Alignment-score evidence for every read in the region, indexed like the region's `bam_reads`.
///
/// Each `BamRead` is ONE alignment record, so a multimapper contributes several records under one name.
/// We group by name, so `best`/`second` are that READ's top two placements anywhere in the region — the
/// familiar "AS of the best hit vs the next best hit". Purely reported; `de` decides.
fn as_evidence_per_read(bam_reads: &[BamRead]) -> Vec<AsEvidence> {
    let aligned_len = |br: &BamRead| -> u32 {
        br.read.cigar.iter().filter(|(op, _)| matches!(op, 'M' | '=' | 'X')).map(|(_, n)| *n).sum::<u64>() as u32
    };
    let mut by_name: std::collections::HashMap<&str, Vec<(i32, u32)>> = std::collections::HashMap::new();
    for br in bam_reads {
        by_name.entry(br.name.as_str()).or_default().push((br.as_score, aligned_len(br)));
    }
    bam_reads
        .iter()
        .map(|br| {
            let placements = &by_name[br.name.as_str()];
            as_evidence(placements).expect("read is its own placement, so the slice is non-empty")
        })
        .collect()
}

/// `NA` for an absent runner-up (single-placement read); otherwise the formatted value.
fn opt_i32(v: Option<i32>) -> String {
    v.map_or_else(|| "NA".to_string(), |x| x.to_string())
}
fn opt_f32(v: Option<f32>) -> String {
    v.map_or_else(|| "NA".to_string(), |x| format!("{x:.3}"))
}

/// One assignment-table row (resolved while the region's reads are in scope).
struct AssignRow {
    read_name: String,
    family_id: String,
    assigned_copy: usize,
    status: &'static str,
    n_decisive: usize,
    margin: f64,
    p_value: f64,
    min_p_value: f64,
    /// Reported alignment-score evidence (see `as_evidence_per_read`). Never feeds the decision.
    as_ev: AsEvidence,
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
/// One per-copy soft-quantification row.
struct QuantRow {
    family_id: String,
    copy_index: usize,
    copy_tid: String,
    /// Genomic span of the copy. Emitted so a catalog can be audited for the same-locus artifact:
    /// two copies of one family whose spans overlap are one locus admitted twice, not two copies.
    copy_chrom: String,
    copy_start: u64,
    copy_end: u64,
    abundance: f64,
    ci: f64,
    n_hard: usize,
}
/// One family-confirmed gene-conversion row.
struct MosaicRow {
    family_id: String,
    copy_a: usize,
    copy_b: usize,
    bp_lo: u64,
    bp_hi: u64,
    n_reads: usize,
    dispersion: u64,
    confirmed: bool,
}
/// One reference-free per-family copy-number row (Task R1: `chi_h` + `depth_cn`, no genome/assembly).
struct FamCnRow {
    family_id: String,
    chrom: String,
    n_copies: usize,
    n_reads: usize,
    chi_h: usize,
    depth_cn: f64, // NaN when --lambda-global was not given
    regime: &'static str,
}
/// One COPY-level historical gene-conversion row (a copy that is a mosaic of two others).
struct CopyConvRow {
    family_id: String,
    copy_c: String,
    copy_a: usize,
    copy_b: usize,
    bp_lo: u64,
    bp_hi: u64,
    n_decisive: usize,
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
    cfg.vg_realign = args.vg_realign; // Task 5 (report-only): off by default, byte-identical otherwise
    cfg.homology_primary = args.homology_primary; // E_r membership; off => the E_c path is untouched
    cfg.filter_readthrough = !args.keep_readthrough; // unspliced pre-mRNA spans are not copies
    cfg.gate.pool_locus_support = !args.no_pool_locus_support;
    cfg.collapse_gate = args.collapse_gate; // experimental; detects paralogy, not collapse (see module header)
    if let Some(e) = args.eps_amb {
        cfg.eps_amb = Some(e);
    }
    let params = AssignParams {
        margin: args.margin,
        error_rate: args.error_rate,
        alpha: args.alpha,
        use_margin_gate: args.margin_gate,
        rna_editing_filter: !args.no_editing_filter,
        edit_rate: args.edit_rate,
        iterative_prune: args.iterative_prune,
        ..AssignParams::default()
    };
    eprintln!("[copy_assign] decisive-margin tau={} error_rate={}", args.margin, args.error_rate);
    let mut family_rows: Vec<FamilyRow> = Vec::new();
    let mut assign_rows: Vec<AssignRow> = Vec::new();
    let mut posterior_lines: Vec<String> = Vec::new();
    // EM-abundance prior for the posterior (else uniform).
    let prior_abundance = std::env::var("RUSTLE_POSTERIOR_PRIOR").ok().as_deref() == Some("abundance");
    // locus from a de-novo tid `DN_<chrom>_<start>_<n>` (chrom may contain `_`, so split from the right).
    fn parse_locus(tid: &str) -> Option<(String, u64)> {
        let rest = tid.strip_prefix("DN_")?;
        let parts: Vec<&str> = rest.rsplitn(3, '_').collect(); // [n, start, chrom]
        Some((parts.get(2)?.to_string(), parts.get(1)?.parse().ok()?))
    }
    let mut quant_rows: Vec<QuantRow> = Vec::new();
    let mut mosaic_rows: Vec<MosaicRow> = Vec::new();
    let mut famcn_rows: Vec<FamCnRow> = Vec::new(); // reference-free chi_H + depth_cn (always emitted)
    let mut copyconv_rows: Vec<CopyConvRow> = Vec::new();
    let mut psv_read_lines: Vec<String> = Vec::new(); // --dump-psv: per-read genotype (alleles at every PSV col)
    let mut psv_copy_lines: Vec<String> = Vec::new(); // --dump-psv: per-copy PSV alleles
    let mut psv_col_lines: Vec<String> = Vec::new(); // --dump-psv: PSV column -> genome position
    let mut em_lines: Vec<String> = Vec::new();           // --em: per-read soft posterior + K-frontier label
    let mut em_abundance_lines: Vec<String> = Vec::new();  // --em: per-copy recovered abundance
    let mut phase_block_lines: Vec<String> = Vec::new();  // --phase: one phase set (PS) per family
    let mut phased_hap_lines: Vec<String> = Vec::new();   // --phase: each haplotype's PSV variant string
    let mut phased_read_lines: Vec<String> = Vec::new();  // --phase: read -> haplotype (HP) haplotag
    // --phase: a self-contained variation graph (GFA) of the phasing — PSV columns = BUBBLES
    // (one segment per allele), copies = PATHS through the bubbles. Loadable in Bandage/vg.
    let mut gfa_segs: HashSet<String> = HashSet::new();        // dedup'd S-lines (shared allele = shared node = bubble anchor)
    let mut gfa_links: HashSet<(String, String)> = HashSet::new();
    let mut gfa_paths: Vec<String> = Vec::new();
    let mut fallback_all: Vec<FallbackEdge> = Vec::new(); // family edges confirmed via the LCS fallback
    let mut dna_needs_rows: Vec<DnaNeedsRecord> = Vec::new(); // --absent-copies: candidates needing DNA validation
    let mut vg_realign_lines: Vec<String> = Vec::new(); // --vg-realign: per-read re-align decisions (report-only)
    let mut gfam = 0usize; // global family counter (unique ids across regions)
    let mut gtf_lines: Vec<String> = Vec::new(); // --gtf: FLAIR-style isoform GTF (transcript + exon rows)

    // `--skip-poa-diagnostic` is read by `detect_and_assign` via this env var (it is purely diagnostic and
    // does not change the emitted families/assignments — see the flag's help).
    if args.skip_poa_diagnostic {
        std::env::set_var("RUSTLE_SKIP_POA_DIAGNOSTIC", "1");
    }
    let timing = std::env::var_os("RUSTLE_TIMING").is_some();
    // Parse the BAM index + header ONCE and reuse across every region (the per-region path re-parses the
    // multi-MB `.bai` otherwise). None => no usable index; fall back to the per-region open (which scans).
    let bam_cache = BamIndexCache::open(&args.bam).ok();
    // Region-parallel pool (opt-in via --region-threads > 1). Sized to region_threads; the inner per-copy
    // poasta parallelism (discover_psvs) composes on the SAME pool, so total concurrency is bounded to N.
    let region_pool = if args.region_threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(args.region_threads)
                .build()
                .context("building region thread pool")?,
        )
    } else {
        None
    };
    // FLAT region list across ALL contigs, in the deterministic by_contig order (sorted contig, file-order
    // ranges). Out-of-order parallel processing over this flat list lets the globally-heaviest families —
    // which live on DIFFERENT contigs — overlap, while the serial drain below (in this same order) keeps
    // CAFAM ids + every row byte-identical to the serial sweep.
    let flat: Vec<(&String, u64, u64)> = by_contig
        .iter()
        .flat_map(|(c, ranges)| ranges.iter().map(move |&(lo, hi)| (c, lo, hi)))
        .collect();
    // Bounded LRU cache of loaded contig genomes — so a worker on any contig reuses an already-loaded genome
    // instead of reloading, and at most ~capacity contig sequences are resident (the memory bound; Arc keeps
    // a genome alive while an evicting worker still uses it). Capacity tracks the concurrency.
    let genome_cap = NonZeroUsize::new((args.region_threads + 1).max(2)).unwrap();
    let genome_cache: Arc<Mutex<LruCache<String, Arc<GenomeIndex>>>> =
        Arc::new(Mutex::new(LruCache::new(genome_cap)));
    let genome_for = |contig: &str| -> Result<Arc<GenomeIndex>> {
        if let Some(g) = genome_cache.lock().unwrap().get(contig).cloned() {
            return Ok(g);
        }
        // load OUTSIDE the lock (a chromosome load is seconds; never block other workers on it). A rare
        // double-load on a concurrent miss is harmless — the second insert just wins.
        let contigs: HashSet<String> = std::iter::once(contig.to_string()).collect();
        let g = Arc::new(
            GenomeIndex::from_fasta_contigs(&args.fasta, &contigs)
                .with_context(|| format!("loading {} for {contig}", args.fasta))?,
        );
        genome_cache.lock().unwrap().put(contig.to_string(), g.clone());
        Ok(g)
    };
    // The expensive, INDEPENDENT per-region work: BAM read + detect_and_assign (the dominant poasta alignment
    // lives here). Pure w.r.t. the read-only genome/bam_cache. The heavy read SEQUENCES are dropped here —
    // only the read NAMES + computed `fams` are returned — so collecting every region's result is lightweight.
    let compute = |contig: &String, lo: u64, hi: u64| -> Result<RegionWork> {
        let genome = genome_for(contig)?;
        let t_read = std::time::Instant::now();
        let (primary, bam_reads) = match &bam_cache {
            Some(c) => c.reads_in_region(&args.bam, contig, lo, hi),
            None => reads_in_region(&args.bam, contig, lo, hi, args.threads),
        }
        .with_context(|| format!("reading {contig}:{lo}-{hi}"))?;
        if timing {
            eprintln!(
                "[timing] reads_in_region {contig}:{lo}-{hi} ({} reads): {:.1}s",
                bam_reads.len(),
                t_read.elapsed().as_secs_f64()
            );
        }
        let extra = if args.recover_copies {
            tied_secondary_reads_in_region(&args.bam, contig, lo, hi, args.as_ratio).unwrap_or_default()
        } else {
            Vec::new()
        };
        let t_da = std::time::Instant::now();
        let (fams, fallback, dna_needs) = detect_and_assign(
            &primary, &bam_reads, &genome, &cfg, args.win, args.min_copies, &params, &extra,
            args.absent_copies, &args.fasta,
        );
        if timing {
            eprintln!("[timing] detect_and_assign {contig}:{lo}-{hi}: {:.1}s", t_da.elapsed().as_secs_f64());
        }
        // FLAIR-style isoform assembly for the optional GTF (intron-chain collapse -> gate -> gene grouping).
        // Recomputed here only under --gtf (cheap: pass1/gate are ~0s); independent of the assignment.
        let transcripts: Vec<TranscriptRec> = if args.gtf {
            let skeletons = pass1_skeletons(&primary, cfg.pass1_min_reads);
            let iso = assemble_gate(&skeletons, &genome, &cfg.gate);
            let groups = collapse_loci_groups(&iso);
            iso.iter()
                .enumerate()
                .map(|(i, t)| TranscriptRec {
                    tid: t.tid.clone(),
                    gene_tid: iso[groups[i]].tid.clone(),
                    chrom: t.chrom.clone(),
                    start: t.start,
                    end: t.end,
                    strand: t.strand,
                    introns: t.introns.clone(),
                    n_reads: t.n_reads,
                })
                .collect()
        } else {
            Vec::new()
        };
        let read_names: Vec<String> = bam_reads.iter().map(|r| r.name.clone()).collect();
        let as_ev = as_evidence_per_read(&bam_reads);
        let n_mapped = bam_reads.len();
        Ok(RegionWork { contig: contig.clone(), lo, hi, read_names, as_ev, n_mapped, fams, fallback, dna_needs, transcripts })
    };
    // Compute all regions (out-of-order across contigs when region_threads > 1), collected in the flat order.
    let works: Vec<RegionWork> = match &region_pool {
        Some(pool) => pool.install(|| {
            flat.par_iter().map(|&(c, lo, hi)| compute(c, lo, hi)).collect::<Result<Vec<_>>>()
        })?,
        None => flat.iter().map(|&(c, lo, hi)| compute(c, lo, hi)).collect::<Result<Vec<_>>>()?,
    };
    // SERIAL drain in the original region order — every row push + the `gfam` id counter is exactly the
    // serial path, so the output is byte-identical.
    {
        for work in works {
            let RegionWork { contig, lo, hi, read_names, as_ev, n_mapped, fams, fallback, dna_needs, transcripts } = work;
            let contig = &contig;
            let bam_reads = &read_names; // output stage indexes read NAMES (sequences were dropped)
            fallback_all.extend(fallback);
            dna_needs_rows.extend(dna_needs);
            // --gtf: gene_tid (a copy's own locus) -> (family id, copy index), filled as fids are assigned below.
            let mut copy_gene: std::collections::HashMap<String, (String, usize)> = std::collections::HashMap::new();
            for fa in &fams {
                let fid = format!("CAFAM{gfam}");
                gfam += 1;
                if args.gtf {
                    for (ci, tid) in fa.copy_tids.iter().enumerate() {
                        copy_gene.insert(tid.clone(), (fid.clone(), ci));
                    }
                }
                for (ri, a) in &fa.assignments {
                    assign_rows.push(AssignRow {
                        read_name: bam_reads[*ri].clone(),
                        family_id: fid.clone(),
                        assigned_copy: a.best_copy,
                        status: status_str(a.status),
                        n_decisive: a.n_decisive,
                        margin: a.log_lr_margin,
                        p_value: a.p_value,
                        min_p_value: a.min_p_value,
                        as_ev: as_ev[*ri],
                    });
                }
                // --vg-realign (report-only): the re-align supplement's per-read decisions for this family.
                // Empty unless --vg-realign was passed (cfg.vg_realign gates run_family_realign itself).
                for r in &fa.realign_records {
                    vg_realign_lines.push(format!(
                        "{}\t{}\t{}\t{}\t{:.6}\t{}",
                        r.read_name, fid, r.action, r.target_copy, r.id_best, r.linear_copy
                    ));
                }
                // soft per-copy POSTERIOR + consistent ZONE (opt-in): localize even the unassignable reads.
                if args.posterior {
                    const FLOOR: f64 = 0.01; // a copy is in the consistent ZONE if its posterior exceeds this
                    let loci: Vec<Option<(String, u64)>> =
                        fa.copy_tids.iter().map(|t| parse_locus(t)).collect();
                    for (ri, a) in &fa.assignments {
                        if a.posterior.len() != fa.copy_tids.len() {
                            continue; // posterior frame must line up with the copy roster (e.g. post-freeze)
                        }
                        // apply the prior (uniform = posterior as-is; else weight by EM abundance), renormalize.
                        let mut post: Vec<f64> = a.posterior.clone();
                        if prior_abundance {
                            for (c, x) in post.iter_mut().enumerate() {
                                *x *= fa.copy_abundance.get(c).copied().unwrap_or(0.0).max(1e-9);
                            }
                            let z: f64 = post.iter().sum();
                            if z > 0.0 {
                                for x in &mut post {
                                    *x /= z;
                                }
                            }
                        }
                        // consistent zone = copies above the floor; its genomic extent + the posterior string.
                        let mut idx: Vec<usize> = (0..post.len()).filter(|&c| post[c] > FLOOR).collect();
                        idx.sort_by(|&a2, &b2| post[b2].partial_cmp(&post[a2]).unwrap());
                        let zone: Vec<u64> = idx.iter().filter_map(|&c| loci[c].as_ref().map(|l| l.1)).collect();
                        let chrom = idx
                            .iter()
                            .find_map(|&c| loci[c].as_ref().map(|l| l.0.clone()))
                            .unwrap_or_default();
                        let (zs, ze) = (zone.iter().min().copied().unwrap_or(0), zone.iter().max().copied().unwrap_or(0));
                        let pstr = idx
                            .iter()
                            .map(|&c| format!("{}:{:.3}", c, post[c]))
                            .collect::<Vec<_>>()
                            .join(",");
                        posterior_lines.push(format!(
                            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                            bam_reads[*ri], fid, status_str(a.status), idx.len(), chrom, zs, ze, pstr
                        ));
                    }
                }
                // soft per-copy abundance (+ the hard read count for comparison)
                for (ci, tid) in fa.copy_tids.iter().enumerate() {
                    quant_rows.push(QuantRow {
                        family_id: fid.clone(),
                        copy_index: ci,
                        copy_tid: tid.clone(),
                        copy_chrom: fa.copy_spans.get(ci).map(|s| s.0.clone()).unwrap_or_default(),
                        copy_start: fa.copy_spans.get(ci).map_or(0, |s| s.1),
                        copy_end: fa.copy_spans.get(ci).map_or(0, |s| s.2),
                        abundance: fa.copy_abundance.get(ci).copied().unwrap_or(0.0),
                        ci: fa.copy_abundance_ci.get(ci).copied().unwrap_or(0.0),
                        n_hard: fa.assignments.iter().filter(|(_, a)| a.best_copy == ci).count(),
                    });
                }
                // gene-conversion: report per-read candidate switches (RT-switch-like) vs recurrence-confirmed
                // events, so the discriminator's anti-artifact gate is visible on real data.
                eprintln!(
                    "[mosaic] {fid}: {} reads showed a candidate copy-switch -> {} recurrence-confirmed conversion event(s)",
                    fa.mosaic_reads,
                    fa.conversions.iter().filter(|e| e.confirmed).count()
                );
                // gene-conversion events (per-molecule PSV-path switches that recur across reads)
                for ev in &fa.conversions {
                    mosaic_rows.push(MosaicRow {
                        family_id: fid.clone(),
                        copy_a: ev.copy_a,
                        copy_b: ev.copy_b,
                        bp_lo: ev.breakpoint_ref.0,
                        bp_hi: ev.breakpoint_ref.1,
                        n_reads: ev.n_supporting_reads,
                        dispersion: ev.breakpoint_dispersion,
                        confirmed: ev.confirmed,
                    });
                }
                // copy-level historical conversions (a copy whose PSV vector is a mosaic of two others)
                for cv in &fa.copy_conversions {
                    copyconv_rows.push(CopyConvRow {
                        family_id: fid.clone(),
                        copy_c: fa.copy_tids.get(cv.copy_c).cloned().unwrap_or_else(|| cv.copy_c.to_string()),
                        copy_a: cv.copy_a,
                        copy_b: cv.copy_b,
                        bp_lo: cv.breakpoint.0,
                        bp_hi: cv.breakpoint.1,
                        n_decisive: cv.n_decisive,
                    });
                }
                // raw per-molecule PSV genotype evidence (the assignment-proof matrix)
                if args.dump_psv {
                    let allele_str = |v: &Vec<Option<u8>>| -> String {
                        v.iter().map(|o| o.map(|b| b as char).unwrap_or('.')).collect()
                    };
                    for ((ri, a), obs) in fa.assignments.iter().zip(fa.read_psv_obs.iter()) {
                        psv_read_lines.push(format!(
                            "{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
                            bam_reads[*ri], fid, a.best_copy, status_str(a.status), a.log_lr_margin,
                            a.n_decisive, allele_str(obs)
                        ));
                    }
                    for (ci, tid) in fa.copy_tids.iter().enumerate() {
                        let alleles = fa.copy_psv_alleles.get(ci).map(allele_str).unwrap_or_default();
                        psv_copy_lines.push(format!("{}\t{}\t{}\t{}", fid, ci, tid, alleles));
                    }
                    for (col, pos) in fa.psv_col_pos.iter().enumerate() {
                        psv_col_lines.push(format!("{}\t{}\t{}", fid, col, pos.map(|x| x as i64).unwrap_or(-1)));
                    }
                }
                // EM soft-relaxation (opt-in): re-runs the family's PSV evidence through the maximum-
                // likelihood EM engine (Task 1's exact gate likelihood, Task 4's binary wiring) for a soft
                // posterior + recovered abundance, alongside (not instead of) the hard PSV+junction
                // assignment above. Fully gated: with `--em` absent this block never runs, so the hard
                // outputs are untouched.
                if args.em {
                    let em_result = em_assign_family(
                        &fa.read_psv_obs,
                        &fa.copy_psv_alleles,
                        &fa.read_junctions,
                        &fa.copy_junctions,
                        &params,
                        args.em_eps,
                        args.em_max_iter,
                    );
                    for (row_idx, (ri, _)) in fa.assignments.iter().enumerate() {
                        if row_idx >= em_result.posteriors.len() {
                            continue; // posterior frame must line up with the read roster (e.g. post-freeze)
                        }
                        let post = &em_result.posteriors[row_idx];
                        let argmax = post
                            .iter()
                            .enumerate()
                            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                            .map(|(k, _)| k)
                            .unwrap_or(0);
                        let label = match em_result.labels[row_idx] {
                            rustle::vg_family::em_copy_assign::EmLabel::Certified => "Certified",
                            rustle::vg_family::em_copy_assign::EmLabel::SoftZone => "SoftZone",
                        };
                        let post_str = post
                            .iter()
                            .enumerate()
                            .map(|(k, p)| format!("{}:{:.4}", k, p))
                            .collect::<Vec<_>>()
                            .join(";");
                        em_lines.push(format!(
                            "{}\t{}\t{}\t{}\t{}\t{}",
                            bam_reads[*ri], fid, argmax, label, post_str, em_result.n_iter
                        ));
                    }
                    for (ci, tid) in fa.copy_tids.iter().enumerate() {
                        let n_reads_soft: f64 = em_result.posteriors.iter().map(|p| p.get(ci).copied().unwrap_or(0.0)).sum();
                        em_abundance_lines.push(format!(
                            "{}\t{}\t{:.4}\t{:.2}",
                            fid,
                            tid,
                            em_result.abundances.get(ci).copied().unwrap_or(0.0),
                            n_reads_soft
                        ));
                    }
                }
                // FACULTATIVE phasing: phase set = family; haplotypes = its copies; read->haplotype =
                // the PSV assignment. A read is phased iff it clears the decisive-margin gate (Assigned);
                // Ambiguous/Tied reads are emitted with haplotype = -1 (unphaseable = K-frontier).
                if args.phase {
                    let n_phased = fa
                        .assignments
                        .iter()
                        .filter(|(_, a)| matches!(a.status, AssignStatus::Assigned))
                        .count();
                    phase_block_lines.push(format!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        fid, fa.chrom, fa.n_copies, fa.psv_cols, n_phased,
                        fa.assignments.len() - n_phased
                    ));
                    for (ci, tid) in fa.copy_tids.iter().enumerate() {
                        let mut vs = String::new();
                        if let Some(alleles) = fa.copy_psv_alleles.get(ci) {
                            for (col, a) in alleles.iter().enumerate() {
                                if let (Some(b), Some(Some(pos))) = (a, fa.psv_col_pos.get(col)) {
                                    if !vs.is_empty() {
                                        vs.push(';');
                                    }
                                    vs.push_str(&format!("{}:{}", pos, *b as char));
                                }
                            }
                        }
                        let n_sup = fa
                            .assignments
                            .iter()
                            .filter(|(_, a)| a.best_copy == ci && matches!(a.status, AssignStatus::Assigned))
                            .count();
                        phased_hap_lines.push(format!("{}\t{}\t{}\t{}\t{}", fid, ci, tid, n_sup, vs));

                        // VG view: this copy is a PATH through the PSV bubbles. A node id is
                        // (family, column, allele) — copies sharing an allele at a column share the
                        // node (a bubble anchor); where they differ the path forks (the bubble).
                        if let Some(alleles) = fa.copy_psv_alleles.get(ci) {
                            let mut path_nodes: Vec<String> = Vec::new();
                            for (col, a) in alleles.iter().enumerate() {
                                if let (Some(b), Some(Some(pos))) = (a, fa.psv_col_pos.get(col)) {
                                    let nid = format!("{}_c{}_{}", fid, col, *b as char);
                                    gfa_segs.insert(format!("S\t{}\t{}\tPO:i:{}", nid, *b as char, pos));
                                    if let Some(prev) = path_nodes.last() {
                                        gfa_links.insert((prev.clone(), nid.clone()));
                                    }
                                    path_nodes.push(nid);
                                }
                            }
                            if !path_nodes.is_empty() {
                                let p: Vec<String> = path_nodes.iter().map(|n| format!("{}+", n)).collect();
                                gfa_paths.push(format!("P\t{}_copy{}\t{}\t*", fid, ci, p.join(",")));
                            }
                        }
                    }
                    for (ri, a) in &fa.assignments {
                        let hap: i64 = if matches!(a.status, AssignStatus::Assigned) {
                            a.best_copy as i64
                        } else {
                            -1
                        };
                        phased_read_lines.push(format!(
                            "{}\t{}\t{}\t{}\t{:.3}\t{}",
                            bam_reads[*ri], fid, hap, a.n_decisive, a.log_lr_margin, status_str(a.status)
                        ));
                    }
                }
                // reference-free per-family copy number (Task R1): chi_H (PSV conflict-structure
                // lower bound) always computed; depth_cn only when --lambda-global was given.
                famcn_rows.push(FamCnRow {
                    family_id: fid.clone(),
                    chrom: fa.chrom.clone(),
                    n_copies: fa.n_copies,
                    n_reads: fa.n_reads,
                    chi_h: chi_h_with_junctions(&fa.copy_psv_alleles, &fa.copy_junctions),
                    depth_cn: args.lambda_global.map(|lam| depth_cn(fa.n_reads, lam)).unwrap_or(f64::NAN),
                    regime: if fa.collapsed_copies > 0 { "reference_collapsed" } else { "reference_resolved" },
                });
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
            eprintln!("[copy_assign]   {contig}:{lo}-{hi}: {} mapped reads -> {} families", n_mapped, fams.len());
            // --gtf: emit every isoform of this region (transcript + exon rows), tagging family-copy genes.
            for t in &transcripts {
                let (fam_attr, multicopy) = match copy_gene.get(&t.gene_tid) {
                    Some((fid, ci)) => (format!(" family_id \"{fid}\"; copy_index \"{ci}\";"), "true"),
                    None => (String::new(), "false"),
                };
                let gs = t.start + 1; // GTF is 1-based, end-inclusive (our coords are 0-based half-open)
                gtf_lines.push(format!(
                    "{}\trustle\ttranscript\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; reads \"{}\"; multicopy \"{}\";{}",
                    t.chrom, gs, t.end, t.strand, t.gene_tid, t.tid, t.n_reads, multicopy, fam_attr
                ));
                // exons = the gene span minus the introns (the read's spliced structure)
                let mut prev = t.start;
                let mut exons: Vec<(u64, u64)> = Vec::new();
                for &(d, a) in &t.introns {
                    exons.push((prev, d));
                    prev = a;
                }
                exons.push((prev, t.end));
                for (k, (es, ee)) in exons.iter().enumerate() {
                    gtf_lines.push(format!(
                        "{}\trustle\texon\t{}\t{}\t.\t{}\t.\tgene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";",
                        t.chrom, es + 1, ee, t.strand, t.gene_tid, t.tid, k + 1
                    ));
                }
            }
        } // for work in works (serial drain, region order)
    } // serial-drain block

    if args.gtf {
        let mut gh = std::fs::File::create(format!("{}.gtf", args.out))?;
        for line in &gtf_lines {
            writeln!(gh, "{line}")?;
        }
        eprintln!("[copy_assign] wrote {}.gtf ({} GTF rows = FLAIR-style isoforms; family copies tagged multicopy)",
            args.out, gtf_lines.len());
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
    writeln!(ah, "read_name\tfamily_id\tassigned_copy\tstatus\tn_decisive\tmargin\tp_value\tmin_p_value\tas_best\tas_second\tas_margin\tas_per_base_best\tas_per_base_2nd")?;
    for r in &assign_rows {
        writeln!(
            ah,
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3e}\t{:.3e}\t{}\t{}\t{}\t{:.3}\t{}",
            r.read_name, r.family_id, r.assigned_copy, r.status, r.n_decisive, r.margin, r.p_value, r.min_p_value,
            r.as_ev.best, opt_i32(r.as_ev.second), opt_i32(r.as_ev.margin()),
            r.as_ev.best_per_base, opt_f32(r.as_ev.second_per_base)
        )?;
    }

    // Reference-free per-family copy number (Task R1, additive; needs no flag): chi_H (PSV
    // conflict-structure lower bound, always computed) + depth_cn (read-depth leg, only when
    // --lambda-global was supplied -- else "NA"). famcn_readonly = the max of the two lower
    // bounds, so it recovers Tier-3 collapsed copies chi_H alone misses.
    let mut cnh = std::fs::File::create(format!("{}.famcn_readonly.tsv", args.out))?;
    writeln!(cnh, "family_id\tchrom\tn_copies\tn_reads\tchi_H\tdepth_cn\tregime\tfamcn_readonly")?;
    for r in &famcn_rows {
        if r.depth_cn.is_finite() {
            let famcn = (r.chi_h as f64).max(r.depth_cn);
            writeln!(
                cnh,
                "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{:.3}",
                r.family_id, r.chrom, r.n_copies, r.n_reads, r.chi_h, r.depth_cn, r.regime, famcn
            )?;
        } else {
            writeln!(
                cnh,
                "{}\t{}\t{}\t{}\t{}\tNA\t{}\t{}",
                r.family_id, r.chrom, r.n_copies, r.n_reads, r.chi_h, r.regime, r.chi_h
            )?;
        }
    }
    eprintln!("[copy_assign] wrote {}.famcn_readonly.tsv ({} families; depth_cn={})",
        args.out, famcn_rows.len(), if args.lambda_global.is_some() { "on" } else { "NA (pass --lambda-global)" });

    // per-read posterior + consistent zone (opt-in via --posterior).
    if args.posterior {
        let mut ph = std::fs::File::create(format!("{}.posterior.tsv", args.out))?;
        writeln!(ph, "read_name\tfamily_id\tstatus\tn_consistent\tzone_chrom\tzone_start\tzone_end\tposterior")?;
        for line in &posterior_lines {
            writeln!(ph, "{line}")?;
        }
        eprintln!("[copy_assign] wrote {}.posterior.tsv ({} reads, prior={})",
            args.out, posterior_lines.len(), if prior_abundance { "abundance" } else { "uniform" });
    }

    // EM soft-relaxation outputs (opt-in via --em): per-read soft posterior + K-frontier label, and the
    // recovered per-copy abundance. Only written under --em; the hard outputs above are unaffected either way.
    if args.em {
        let mut eh = std::fs::File::create(format!("{}.em.tsv", args.out))?;
        writeln!(eh, "read_name\tfamily_id\targmax_copy\tlabel\tposterior\tn_iter")?;
        for l in &em_lines {
            writeln!(eh, "{l}")?;
        }
        let mut eah = std::fs::File::create(format!("{}.em_abundance.tsv", args.out))?;
        writeln!(eah, "family_id\tcopy_id\tpi_hat\tn_reads_soft")?;
        for l in &em_abundance_lines {
            writeln!(eah, "{l}")?;
        }
        eprintln!(
            "[copy_assign] wrote {}.em.tsv ({} reads) + {}.em_abundance.tsv",
            args.out, em_lines.len(), args.out
        );
    }

    // soft per-copy quantification: family/copy, EM abundance ± 95% CI half-width, + the hard read count for
    // comparison. The EM uses partial PSV evidence (the benchmark: beats hard at sparse PSVs; uniform at K=0).
    let mut qh = std::fs::File::create(format!("{}.quant.tsv", args.out))?;
    writeln!(qh, "family_id\tcopy_index\tcopy_tid\tcopy_chrom\tcopy_start\tcopy_end\tabundance\tci95_halfwidth\tn_reads_hard")?;
    for r in &quant_rows {
        writeln!(qh, "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}", r.family_id, r.copy_index, r.copy_tid,
            r.copy_chrom, r.copy_start, r.copy_end, r.abundance, r.ci, r.n_hard)?;
    }

    // FACULTATIVE long-read phasing output (dependency-free): phase set (PS) per family, each haplotype's
    // PSV variant string, and read -> haplotype (HP) haplotag. Only written under --phase.
    if args.phase {
        let mut pb = std::fs::File::create(format!("{}.phase_blocks.tsv", args.out))?;
        writeln!(pb, "block_id\tchrom\tn_haplotypes\tn_psv_sites\tn_reads_phased\tn_unphased")?;
        for l in &phase_block_lines {
            writeln!(pb, "{}", l)?;
        }
        let mut ph = std::fs::File::create(format!("{}.phased_haplotypes.tsv", args.out))?;
        writeln!(ph, "block_id\thaplotype\tcopy_tid\tn_support_reads\tvariants")?;
        for l in &phased_hap_lines {
            writeln!(ph, "{}", l)?;
        }
        let mut pr = std::fs::File::create(format!("{}.phased_reads.tsv", args.out))?;
        writeln!(pr, "read_name\tblock_id\thaplotype\tn_psv_spanned\tmargin\tstatus")?;
        for l in &phased_read_lines {
            writeln!(pr, "{}", l)?;
        }
        // self-contained variation graph of the phasing (copies = paths, PSVs = bubbles)
        let mut gf = std::fs::File::create(format!("{}.phase.gfa", args.out))?;
        writeln!(gf, "H\tVN:Z:1.0")?;
        let mut segs: Vec<&String> = gfa_segs.iter().collect();
        segs.sort();
        for s in segs {
            writeln!(gf, "{}", s)?;
        }
        let mut links: Vec<&(String, String)> = gfa_links.iter().collect();
        links.sort();
        for (a, b) in links {
            writeln!(gf, "L\t{}\t+\t{}\t+\t0M", a, b)?;
        }
        for p in &gfa_paths {
            writeln!(gf, "{}", p)?;
        }
        let n_phased = phased_read_lines.iter().filter(|l| !l.contains("\t-1\t")).count();
        eprintln!(
            "[copy_assign] phasing: {} blocks, {} haplotypes, {}/{} reads phased -> {}.phased_*.tsv + {}.phase.gfa ({} bubble-nodes, {} copy-paths)",
            phase_block_lines.len(), phased_hap_lines.len(), n_phased, phased_read_lines.len(),
            args.out, args.out, gfa_segs.len(), gfa_paths.len()
        );
    }

    // gene-conversion events: per-molecule PSV-path switches confirmed by RECURRENCE across reads (vs one-off
    // chimeras). Only written when something was found. The enriched per-molecule multimapper signal.
    if !mosaic_rows.is_empty() {
        let mut mh = std::fs::File::create(format!("{}.mosaic.tsv", args.out))?;
        writeln!(mh, "family_id\tcopy_a\tcopy_b\tbreakpoint_lo\tbreakpoint_hi\tn_reads\tdispersion\tconfirmed")?;
        for r in &mosaic_rows {
            writeln!(
                mh,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.family_id, r.copy_a, r.copy_b, r.bp_lo, r.bp_hi, r.n_reads, r.dispersion, r.confirmed as u8
            )?;
        }
        let conf = mosaic_rows.iter().filter(|r| r.confirmed).count();
        eprintln!(
            "[copy_assign] {} read-level gene-conversion event(s) ({} confirmed by recurrence) -> {}.mosaic.tsv",
            mosaic_rows.len(), conf, args.out
        );
    }

    // copy-level historical gene conversions (a de-novo copy whose PSV-allele vector is a mosaic of two
    // others) -- the APOBEC3/RFPL signal, baked into the copy sequence. Written only when found.
    if !copyconv_rows.is_empty() {
        let mut ch = std::fs::File::create(format!("{}.copy_conversions.tsv", args.out))?;
        writeln!(ch, "family_id\tconverted_copy\tdonor_a\tdonor_b\tbreakpoint_lo\tbreakpoint_hi\tn_decisive")?;
        for r in &copyconv_rows {
            writeln!(
                ch, "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.family_id, r.copy_c, r.copy_a, r.copy_b, r.bp_lo, r.bp_hi, r.n_decisive
            )?;
        }
        eprintln!(
            "[copy_assign] {} COPY-level historical gene conversion(s) -> {}.copy_conversions.tsv",
            copyconv_rows.len(), args.out
        );
    }

    // --dump-psv: the raw per-molecule PSV genotype matrix (the assignment-proof evidence). reads × PSV columns
    // (each read's base + its assignment), the per-copy alleles, and the column→genome map — for the figure.
    if args.dump_psv {
        let mut rh = std::fs::File::create(format!("{}.psv_reads.tsv", args.out))?;
        writeln!(rh, "read_name\tfamily_id\tassigned_copy\tstatus\tmargin\tn_decisive\talleles")?;
        for l in &psv_read_lines {
            writeln!(rh, "{l}")?;
        }
        let mut ch = std::fs::File::create(format!("{}.psv_copies.tsv", args.out))?;
        writeln!(ch, "family_id\tcopy_index\tcopy_tid\talleles")?;
        for l in &psv_copy_lines {
            writeln!(ch, "{l}")?;
        }
        let mut lh = std::fs::File::create(format!("{}.psv_cols.tsv", args.out))?;
        writeln!(lh, "family_id\tcol_index\tgenome_pos")?;
        for l in &psv_col_lines {
            writeln!(lh, "{l}")?;
        }
        eprintln!(
            "[copy_assign] dumped PSV genotype matrix: {} read rows, {} copy rows -> {}.psv_reads.tsv/.psv_copies.tsv/.psv_cols.tsv",
            psv_read_lines.len(), psv_copy_lines.len(), args.out
        );
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

    // --absent-copies: surface candidates that failed the admission gate and need DNA-level validation.
    // Only written when --absent-copies is set so an OFF run produces exactly the same output files.
    if args.absent_copies {
        let mut dh = std::fs::File::create(format!("{}.dna_needs.tsv", args.out))?;
        writeln!(dh, "chrom\tstart\tend\tn_clusters\tread_count\treason")?;
        for r in &dna_needs_rows {
            writeln!(dh, "{}\t{}\t{}\t{}\t{}\t{}", r.chrom, r.start, r.end, r.n_clusters, r.read_count, r.reason)?;
        }
        eprintln!(
            "[copy_assign] {} DNA-needs candidate(s) -> {}.dna_needs.tsv",
            dna_needs_rows.len(),
            args.out
        );
    }

    // --vg-realign: the re-align supplement's per-family/per-read decisions (report-only — not fed back
    // into the assignment). Only written when --vg-realign is set so an OFF run produces exactly the same
    // output files (cfg.vg_realign OFF also means fa.realign_records is always empty, so this is belt-
    // and-suspenders with the flag check).
    if args.vg_realign {
        let mut vh = std::fs::File::create(format!("{}.vg_realign.tsv", args.out))?;
        writeln!(vh, "read_name\tfamily_id\taction\ttarget_copy\tid_best\tlinear_copy")?;
        for l in &vg_realign_lines {
            writeln!(vh, "{l}")?;
        }
        eprintln!(
            "[copy_assign] {} vg-realign decision(s) -> {}.vg_realign.tsv",
            vg_realign_lines.len(),
            args.out
        );
    }

    let (uniq, agree): (usize, usize) = family_rows.iter().fold((0, 0), |(u, g), f| (u + f.uniq, g + f.uniq_agree));
    eprintln!(
        "[copy_assign] {} families, {} read assignments",
        family_rows.len(),
        assign_rows.len()
    );

    // Same-locus artifact: two copies of ONE family whose genomic spans OVERLAP are one locus admitted
    // twice, not two copies. Such a family reports min_p == 1 for every read, so it abstains wholesale and
    // its reads masquerade as the K=0 identifiability wall. Warn loudly rather than fail — the catalog is
    // still emitted, but its abstention must not be read as biology. `bench/artifact_audit.py` audits this.
    {
        let catalog: Vec<(String, String, u64, u64)> = quant_rows
            .iter()
            .map(|r| (r.family_id.clone(), r.copy_chrom.clone(), r.copy_start, r.copy_end))
            .collect();
        let flagged = catalog_overlaps(&catalog);
        if !flagged.is_empty() {
            // Two distinct meanings, and they should not be described the same way (verified on RFPL/r4,
            // bench/CONTAINMENT_COVERAGE_FLOOR.md):
            //  - DuplicateLocus (recip ~ 1): one locus admitted twice. Every read scores min_p == 1, so the
            //    family abstains wholesale and its reads masquerade as the K=0 wall.
            //  - Containment (recip << 1): a shorter transcript nested/staggered inside a longer one. On
            //    low-coverage regions this is usually a fragment or a chimeric readthrough, so the COPY COUNT
            //    is inflated — NOT the min_p == 1 masquerade. It cannot be pruned without also deleting
            //    genuine overlapping tandem paralogs (they occupy the same feature cell), so it is reported,
            //    not removed.
            let n_dup = flagged.iter().filter(|f| f.3 == OverlapKind::DuplicateLocus).count();
            let n_contain = flagged.iter().filter(|f| f.3 == OverlapKind::Containment).count();
            let n_shared = flagged.iter().filter(|f| f.3 == OverlapKind::SharedAcrossFamilies).count();
            eprintln!(
                "[copy_assign] WARNING: {} copy pair(s) share genomic sequence \
                 ({n_dup} DuplicateLocus, {n_contain} Containment, {n_shared} SharedAcrossFamilies). \
                 DuplicateLocus = one locus twice, its reads abstain at min_p == 1 (not the K=0 wall). \
                 Containment = a fragment/readthrough nested in a real copy, inflating the copy count on \
                 low-coverage regions — reported, not pruned (it shares its feature cell with real \
                 overlapping paralogs).",
                flagged.len()
            );
            for &(i, j, recip, kind) in flagged.iter().take(10) {
                eprintln!(
                    "[copy_assign]   {kind:?} recip={recip:.2}  {}/{}:{}-{}  vs  {}/{}-{}",
                    catalog[i].0, catalog[i].1, catalog[i].2, catalog[i].3, catalog[j].0, catalog[j].2, catalog[j].3
                );
            }
        }
    }

    if uniq > 0 {
        eprintln!(
            "[copy_assign] genome-wide unique-mapper agreement: {agree}/{uniq} ({:.1}%)",
            100.0 * agree as f64 / uniq as f64
        );
    }
    eprintln!("[copy_assign] wrote {0}.families.tsv + {0}.assignments.tsv + {0}.quant.tsv", args.out);
    Ok(())
}
