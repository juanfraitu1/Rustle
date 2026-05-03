//! Rustle - Long-read transcript assembler (Rust).

use clap::{ArgAction, Parser};
use rustle::merge_mode::run_merge;
use rustle::run;
use rustle::RunConfig;

#[derive(Parser, Debug)]
#[command(
    name = "rustle",
    about = "Transcript assembler",
    long_about = "Assemble RNA-Seq alignments into potential transcripts.\n\
compatible options: -G (guide), -l (label), -o (output), -p (threads), -L (long reads), -e (eonly),\n\
-f/-c/-s/-m/-j/-g (coverage/length/junction/gap), -A (gene abundance), -B/-b (Ballgown),\n\
-v (verbose), -t (no coverage trim), -x (exclude seqs), -E (splice window), --nasc, --conservative.\n\
\n\
Benchmarking / gffcompare: compare de novo output to a reference GTF (run rustle without -G; pass the reference only to gffcompare -r).\n\
Use -G/--guide only for debugging or guided-mode compatibility, not for headline sensitivity/precision vs truth."
)]
struct Args {
    /// Input BAM file
    bam: Option<String>,

    /// Output GTF file -o)
    #[arg(short = 'o', long)]
    output: Option<String>,

    /// Long-read mode (default). -L is accepted for compatibility but always on.
    #[arg(short = 'L', long, default_value_t = true)]
    long_reads: bool,

    /// Process only this chromosome
    #[arg(long)]
    chrom: Option<String>,

    /// Reference annotation to guide assembly -G). GTF/GFF.
    /// For benchmarking against truth, omit this (de novo). Use only for debug or explicit guided/e-only compatibility.
    #[arg(short = 'G', long = "guide")]
    guide: Option<String>,

    /// Genome FASTA for splice consensus validation (GT/GC-AG, CT-AC/GC). CRAM: --cram-ref
    #[arg(long = "genome-fasta")]
    genome_fasta: Option<String>,

    /// Reference genome for CRAM input --ref/--cram-ref)
    #[arg(long = "cram-ref")]
    cram_ref: Option<String>,

    /// Merge mode: merge transcript models from GTF input(s) instead of assembling from BAM.
    #[arg(long)]
    merge: bool,

    /// Input GTF for merge mode (repeatable).
    #[arg(long = "merge-input")]
    merge_inputs: Vec<String>,

    /// File containing one GTF path per line for merge mode.
    #[arg(long = "merge-list")]
    merge_list: Option<String>,

    /// Disable graph-based merge (use structural merge only). Default: use CMTransfrag/process_merge_transfrags/printMergeResults-style merge.
    #[arg(long)]
    no_graph_merge: bool,

    /// Name prefix for output transcripts [default: RSTL]
    #[arg(short = 'l', long, default_value = "RSTL")]
    label: String,

    /// Minimum junction coverage -j) [default: 1]
    #[arg(short = 'j', long, default_value = "1.0")]
    min_junction_reads: f64,

    /// Minimum reads per bp for multi-exon transcript -c) [default: 1]
    #[arg(short = 'c', long, default_value = "1.0")]
    readthr: f64,

    /// Minimum reads per bp for single-exon transcript -s) [default: 4.75]
    #[arg(short = 's', long, default_value = "4.75")]
    singlethr: f64,

    /// Minimum assembled transcript length in bp -m) [default: 200]
    #[arg(short = 'm', long, default_value = "200")]
    min_transcript_length: u64,

    /// Minimum anchor length for junctions in bp -a) [default: 10]
    #[arg(short = 'a', long, default_value = "10")]
    min_anchor_length: u64,

    /// Minimum intron length (bp) to count as splice junction [default: 20]
    #[arg(long, default_value = "20")]
    min_intron_length: u64,

    /// Per-splice-site isofrac (%): remove junctions below this fraction of donor total [default: 0]
    #[arg(long, default_value = "0")]
    per_splice_site_isofrac: f64,

    /// Maximum gap between read mappings in bp -g). Long-read: 0 [default: 0]
    #[arg(short = 'g', long, default_value = "0")]
    bundle_merge_dist: u64,

    /// Bundle runoff (bp): new bundle when read start > current end + this [default: 200]
    #[arg(long, default_value = "200")]
    bundle_runoff_dist: u64,

    /// Minimum isoform fraction -f); used for transcript filter and pairwise (isofrac) [default: 0.01]
    #[arg(short = 'f', long, default_value = "0.01")]
    transcript_isofrac: f64,

    /// Disable coverage-based trimming -t): prediction trim where implemented, and
    /// long-read per-read terminal coverage trim before junction/graph (reduces gffcompare class `j`).
    #[arg(short = 't', long)]
    no_coverage_trim: bool,

    /// Do not assemble on these reference sequence(s), comma-separated -x)
    #[arg(short = 'x', long)]
    exclude_seqids: Option<String>,

    /// Window around splice sites for long-read error margin in bp -E) [default: 0].
    /// A wide window can snap PacBio/ONT donors and acceptors away from StringTie’s
    /// SSERROR-based bundle coordinates; on GGO_19 IsoSeq vs `stringtie -L` this default
    /// recovers a few exact transcript matches. Use `-E 30` for the legacy snap behavior.
    #[arg(short = 'E', long, default_value = "0")]
    junction_correction_window: u64,

    /// Conservative assembly --conservative): same as -t -c 1.5 -f 0.05
    #[arg(long)]
    conservative: bool,

    /// Output nascent transcripts --nasc)
    #[arg(long = "nasc")]
    isnascent: bool,

    /// Number of threads -p); 0 = auto-detect available CPUs [default: 0]
    #[arg(short = 'p', long, default_value = "0")]
    threads: usize,

    /// Write gene abundance to file -A)
    #[arg(short = 'A', long)]
    gene_abundance_out: Option<String>,

    /// Ballgown table files in this directory -B)
    #[arg(short = 'b', long)]
    ballgown: Option<String>,

    /// Compare bundle regions to the original algorithm -v log (e.g. GGO_19.log); print match/only_rust/only_log counts
    #[arg(long)]
    bundles: Option<String>,

    /// CGroup: junction support for small-exon keep rule [default: 10]
    #[arg(long, default_value = "10")]
    cgroup_junction_support: u64,

    /// Minimum exon length (bp); drop transcript if any exon shorter [default: 0 = off]
    #[arg(long, default_value = "0")]
    min_exon_length: u64,

    /// Junction support for graph building [default: 10]
    #[arg(long, default_value = "10")]
    junction_support: u64,

    /// Junction correction: demote when best_reads * tolerance > reads_j [default: 0 = off]
    #[arg(long, default_value = "0")]
    junction_correction_tolerance: f64,

    /// Junction canonicalization: merge boundary variants within this many bp [default: 0 = off]
    #[arg(long, default_value = "0")]
    junction_canonical_tolerance: u64,

    /// Remove transcripts fully contained in another with >= coverage [default: off]
    #[arg(long, action = ArgAction::Set, default_value_t = false)]
    filter_contained: bool,

    /// Keep isoforms with read abundance (longcov/cov) >= this through long-read isofrac (0=off) [default: 1]
    #[arg(long, default_value = "1")]
    transcript_isofrac_keep_min: f64,

    /// Low isofrac: keep in pairwise filter if cov >= this * container (lowisofrac; lower = more permissive) [default: 0.02]
    #[arg(long, default_value = "0.02")]
    lowisofrac: f64,

    /// Merge consecutive exons when gap <= this (bp); 0 = off
    #[arg(long, default_value = "0")]
    merge_micro_intron_max_gap: u64,

    /// Output TPM and coverage in GTF attributes
    #[arg(long)]
    output_tpm: bool,

    /// Rawreads mode (rawreads branch): emit models directly from active transfrags.
    #[arg(long)]
    rawreads: bool,

    /// Merge-mode BAM compatibility scaffold flag; currently warning-only.
    #[arg(long)]
    merge_bam_compat: bool,

    /// Emit gene feature lines in GTF (one per transcript)
    #[arg(long)]
    gtf_gene_lines: bool,

    /// Filter multi-exon transcripts with terminal exon < N bp (0 = off, 35 = Python SMALL_EXON)
    #[arg(long, default_value = "0")]
    min_terminal_exon_len: u64,

    /// Produce more transcripts (aim for 100+): relax filters and disable pairwise overlap filter.
    /// Equivalent to --relax-filters with pairwise-overlap-filter=false.
    #[arg(long)]
    more_transcripts: bool,

    /// Relax filters: skip intronic, longread-interval, and incomplete-vs-guide filters to keep more intron chains
    #[arg(long)]
    relax_filters: bool,

    /// Pairwise overlap filter: junction inclusion + containment (isofrac/error_perc/drop)
    #[arg(long, action = ArgAction::Set, default_value_t = true)]
    pairwise_overlap_filter: bool,

    /// Pairwise filter: isofrac override; if 0, -f (transcript_isofrac) is used [default: 0.01]
    #[arg(long, default_value = "0.01")]
    pairwise_isofrac: f64,

    /// Pairwise filter: error_perc for intronic containment [default: 0.1]
    #[arg(long, default_value = "0.1")]
    pairwise_error_perc: f64,

    /// Pairwise filter: drop ratio for junction-inclusion [default: 0.5]
    #[arg(long, default_value = "0.5")]
    pairwise_drop: f64,

    /// Retained intron filter: remove low-cov transcripts spanning introns of higher-cov
    #[arg(long)]
    retained_intron_filter: bool,

    /// Retained intron: error_perc [default: 0.1]
    #[arg(long, default_value = "0.1")]
    retained_intron_error_perc: f64,

    /// Polymerase run-off filter: remove isolated single-exon fragments near multi-exon ends
    #[arg(long)]
    polymerase_runoff_filter: bool,

    /// Enable polymerase run-on filter (bridge split between genes)
    #[arg(long)]
    polymerase_runon_filter: bool,

    /// Polymerase run-off: max distance (bp) to multi-exon end [default: 0]
    #[arg(long, default_value = "0")]
    polymerase_runoff_dist: u64,

    /// Polymerase run-off: singlethr (cov must be >= neighbor + this) [default: 4.75]
    #[arg(long, default_value = "4.75")]
    polymerase_singlethr: f64,

    /// Dedup: remove transcripts whose intron chain is a proper subset of another (same strand)
    #[arg(long, default_value_t = true)]
    dedup_intron_chain_subsets: bool,

    /// Enable graph-build longtrim-like node splitting at read-boundary peaks (lstart/lend).
    /// Disabled by default until this path matches the original algorithm more closely for long-read mode.
    #[arg(long)]
    enable_longtrim: bool,

    /// Longtrim: minimum boundary coverage to split a node [default: 2.0]
    #[arg(long, default_value = "2.0")]
    longtrim_min_boundary_cov: f64,

    /// Enable splice-aware path extension during transcript extraction.
    #[arg(long)]
    enable_splice_path_extension: bool,

    /// Strict mode: enable direct port of back_to_source_fast_long/fwd_to_sink_fast_long with safe fallback [default: on].
    #[arg(long, action = ArgAction::Set, default_value_t = true)]
    strict_longrec_port: bool,

    /// Max active internal graph nodes before pruning low-coverage nodes [default: 1000000]
    #[arg(long, default_value = "1000000")]
    allowed_graph_nodes: usize,

    /// Prune by removing junction edges (style) instead of removing nodes.
    #[arg(long)]
    use_prune_by_junction: bool,

    /// Estimate-only mode: only output transcripts matching guide annotations (-e flag).
    #[arg(short = 'e', long)]
    eonly: bool,

    /// Verbose
    #[arg(short, long)]
    verbose: bool,

    /// Debug: trace flow computation for specific bundle (e.g., "NC_073243.2:17159009-17178368")
    #[arg(long)]
    debug_bundle: Option<String>,

    /// Debug: dump junction sets for matching bundle(s).
    #[arg(long)]
    debug_junctions: bool,

    /// Restrict processing to bundles overlapping `--debug-bundle`.
    #[arg(long)]
    only_debug_bundle: bool,

    /// Debug: trace a specific intron chain (donor-acceptor pairs, comma-separated).
    /// Example: "17433882-17439906,17440121-17444180"
    #[arg(long)]
    trace_intron_chain: Option<String>,

    /// Emit contained isoforms by trimming transcript ends (may increase output).
    #[arg(long)]
    emit_contained_isoforms: bool,

    /// Emit additional transcripts by enumerating splice graph paths (max sensitivity; not default).
    #[arg(long, action = ArgAction::Set, default_value_t = true)]
    emit_junction_paths: bool,

    /// Max number of junction-path transcripts to emit per bundle [default: 1200]
    #[arg(long, default_value = "1200")]
    max_junction_paths: usize,

    /// Max nodes in a junction-path (including source/sink) [default: 200]
    #[arg(long, default_value = "200")]
    max_junction_path_len: usize,

    /// Zero all isofrac/coverage filters for maximum sensitivity.
    /// Default: false (set `--max-sensitivity true` to enable).
    #[arg(long, action = ArgAction::Set, default_value_t = false)]
    max_sensitivity: bool,

    /// Apply standard density rules (`-f` / pairwise floors, junction-path cap, containment).
    /// **On by default**. Skipped when `--max-sensitivity` is set.
    #[arg(long, action = ArgAction::SetTrue)]
    compat_preset: bool,

    /// Trace reference transcripts: report why each reference is missing (NotExtracted vs filter).
    #[arg(long)]
    trace_reference: Option<String>,

    /// Emit reference isoforms from junction support when graph/guide path extraction fails (compatibility).
    /// Also on when `--trace-reference` is set; disable with `RUSTLE_REF_JUNCTION_WITNESS=0`.
    #[arg(long, action = ArgAction::SetTrue)]
    ref_junction_witness: bool,

    /// Output path for trace report (defaults to stdout if omitted).
    #[arg(long)]
    trace_output: Option<String>,

    /// Output TSV path for per-bundle stage summaries (bundles, reads, junctions, seeds, longrec).
    /// If omitted, needy-locus runs still write `{output}.needy_bundle_summary.tsv` when tracing
    /// (see `RUSTLE_DISABLE_AUTO_BUNDLE_SUMMARY`, `RUSTLE_DEBUG_STAGE_TSV`).
    #[arg(long)]
    debug_stage_tsv: Option<String>,

    /// Output TSV path for keeptrf/usepath export from process_transfrags.
    #[arg(long)]
    keeptrf_usepath_tsv: Option<String>,

    /// Output JSONL path: dump per-bundle snapshots (graph + transfrags + transcripts) for 1:1 comparisons.
    #[arg(long)]
    snapshot_jsonl: Option<String>,

    // ── Variation graph (gene family) mode ───────────────────────────────────

    /// Enable variation graph mode: link gene family copies via multi-mapping reads,
    /// then jointly resolve multi-mapping weights using EM.
    #[arg(long)]
    vg: bool,

    /// Minimum shared multi-mapping reads to link two bundles into the same
    /// family group at discovery time [default: 3]. The post-discovery
    /// quality filter (--vg-family-min-shared) imposes a separate threshold
    /// on the FAMILY's total multimap_reads count.
    #[arg(long, default_value = "3")]
    vg_min_shared: usize,

    /// Family-quality filter: drop families with fewer than N total
    /// multi-mapping reads. Default 10 — empirically the noise/signal
    /// boundary on GGO.bam (32% of discovered families have <10 shared
    /// reads and are mostly random alignment artifacts).
    #[arg(long, default_value_t = 10)]
    vg_family_min_shared: usize,

    /// Family-quality filter: drop families with more than N copies.
    /// Default 30 — catches mtDNA mega-clusters (1509 "copies" on GGO.bam)
    /// and other non-paralog repetitive elements while keeping legitimate
    /// gene-family clusters (largest tested: TBC1D3 ~16, OR cluster ~5).
    #[arg(long, default_value_t = 30)]
    vg_family_max_copies: usize,

    /// Family-quality filter: drop families where shared multi-mapping
    /// reads per copy < N. Default 1.0 — a family where the average copy
    /// has fewer than 1 cross-mapping read is almost always an alignment
    /// artifact, not a real paralog cluster. Set to 0 to disable.
    #[arg(long, default_value_t = 1.0)]
    vg_family_min_shared_per_copy: f64,

    /// Family-quality filter: drop families with coefficient of variation
    /// (CV) of per-copy intron counts above N. Default 1.5 — only the most
    /// extreme cases (mixed-gene mega-clusters where some copies have 1
    /// intron and others have 20+) are dropped. Real paralog clusters
    /// typically have CV well below 1.0 even with uneven coverage.
    /// Single-exon paralog clusters skip this check. Set to 0 to disable.
    #[arg(long, default_value_t = 1.5)]
    vg_family_max_exon_cv: f64,

    /// Family-quality filter: minimum mean pairwise Jaccard of intron-length
    /// sets (binned at 50bp) across copies. THE GRAPH-STRUCTURAL SIGNAL —
    /// real paralogs share intron lengths; cross-cluster TE-bridge merges
    /// (chr19-GOLGA8 ↔ chr17-TBC1D3 etc.) have Jaccard near 0. Default 0.20
    /// — strict enough to filter most TE-bridges, lenient enough to keep
    /// real paralogs even with copy-number variation. Single-exon paralog
    /// clusters skip this check. Set to 0 to disable.
    #[arg(long, default_value_t = 0.20)]
    vg_family_min_primitive_jaccard: f64,

    /// Family-quality filter (optional, requires --genome-fasta): minimum
    /// mean pairwise k-mer Jaccard of FAMILY-GRAPH per-copy sequences.
    /// The 6th signal — catches mild TE-bridges that pass the
    /// intron-length-Jaccard signal by coincidence (similar intron sizes
    /// between unrelated genes) but have no actual sequence similarity.
    /// Requires building the family graph (uses --genome-fasta sequences).
    /// Default 0 (disabled) since it requires the genome FASTA. Try 0.05
    /// to start — real paralogs typically have >0.10 pairwise k-mer Jaccard
    /// at 15-mer resolution.
    #[arg(long, default_value_t = 0.0)]
    vg_family_min_kmer_jaccard: f64,

    /// Maximum EM iterations for multi-mapping resolution [default: 20]
    #[arg(long, default_value = "20")]
    vg_em_iter: usize,

    /// Discover novel gene copies from poorly-mapped/unmapped reads
    #[arg(long)]
    vg_discover_novel: bool,

    /// Output path for family group report TSV
    #[arg(long)]
    vg_report: Option<String>,

    /// Multi-mapping resolution method: none (discover only), em, em-hmm, or flow [default: em].
    /// `em` replaces the default 1/NH weighting on multi-mappers with evidence-based
    /// (junction-compatibility + context) weights. `em-hmm` replaces the heuristic
    /// E-step with HMM-based per-paralog forward log-likelihood scoring (sequence-aware;
    /// requires a `--genome` FASTA and triggers a one-pass BAM scan to collect read
    /// sequences for multi-mappers). Use `none` to disable redistribution and keep
    /// StringTie-equivalent 1/NH behaviour.
    #[arg(long, default_value = "em")]
    vg_solver: String,

    /// Use SNPs (from MD tag) for copy assignment in VG mode
    #[arg(long)]
    vg_snp: bool,

    /// Phased assembly: split bundles by HP tag haplotype in VG mode
    #[arg(long)]
    vg_phase: bool,

    /// Minimum reads to create a novel copy bundle [default: 3]
    #[arg(long, default_value = "3")]
    vg_min_novel_reads: usize,

    /// VG novel-copy discovery algorithm: `kmer` (legacy) or `hmm` (family RNA-HMM).
    #[arg(long = "vg-discover-novel-mode", default_value = "kmer")]
    vg_discover_novel_mode: String,

    /// Run external minimap2 verification on each rescued read (slower).
    #[arg(long = "vg-rescue-diagnostic")]
    vg_rescue_diagnostic: bool,

    /// Forward log-odds threshold for HMM rescue (nats). Default 30.0.
    #[arg(long = "vg-rescue-min-loglik", default_value_t = 30.0)]
    vg_rescue_min_loglik: f64,

    /// LOO experiment harness: mask reads whose primary alignment overlaps
    /// `chrom:start-end` so they enter the HMM rescue pool as if unaligned,
    /// instead of forming a normal bundle. Repeatable. Use with
    /// `--vg --vg-discover-novel --vg-discover-novel-mode hmm` and
    /// `--genome-fasta` to verify novel-copy recovery from sequence alone.
    #[arg(long = "vg-mask-region", value_name = "CHROM:START-END")]
    vg_mask_region: Vec<String>,

    /// Skip EM on families with more than N copies. Default 20 — bigger
    /// groups are typically noise (mtDNA, megafamilies) where EM is either
    /// uninformative or computationally infeasible.
    #[arg(long = "vg-em-max-copies", default_value_t = 20)]
    vg_em_max_copies: usize,

    /// Cap for routing a family to HMM-EM under `--vg-solver auto`.
    /// Default 10 — HMM scoring scales as O(reads × copies × seq × profile);
    /// above this, fall back to heuristic EM for the family.
    #[arg(long = "vg-em-hmm-max-copies", default_value_t = 10)]
    vg_em_hmm_max_copies: usize,

    /// Skip EM on intronless families (default true). Single-exon paralogs
    /// (olfactory receptors etc.) yield degenerate family graphs.
    #[arg(long = "vg-em-no-skip-intronless")]
    vg_em_no_skip_intronless: bool,
}

fn parse_intron_chain(raw: &str) -> anyhow::Result<Vec<(u64, u64)>> {
    let mut out = Vec::new();
    for item in raw.split(',') {
        let item = item.trim();
        if item.is_empty() {
            continue;
        }
        let Some((a, b)) = item.split_once('-') else {
            anyhow::bail!("invalid intron spec '{}': expected donor-acceptor", item);
        };
        let donor: u64 = a.trim().parse()?;
        let acceptor: u64 = b.trim().parse()?;
        if donor >= acceptor {
            anyhow::bail!("invalid intron spec '{}': donor must be < acceptor", item);
        }
        out.push((donor, acceptor));
    }
    if out.is_empty() {
        anyhow::bail!("trace intron chain is empty");
    }
    Ok(out)
}

pub fn run_cli() -> anyhow::Result<()> {
    let args = Args::parse();

    // Emit StringTie-exact banner if RUSTLE_STRINGTIE_EXACT=1.
    rustle::stringtie_parity::maybe_emit_banner();

    // Initialize parity-decision JSONL log if RUSTLE_PARITY_LOG is set.
    rustle::parity_decisions::init_from_env();

    let output = args
        .output
        .as_deref()
        .ok_or_else(|| anyhow::anyhow!("output path required (e.g. -o out.gtf)"))?;

    // Long-read mode is always on.
    let long_reads = true;
    let long_read_min_len = 0u64;

    let (readthr, transcript_isofrac, no_coverage_trim) = if args.conservative {
        (1.5, 0.05, true)
    } else {
        (args.readthr, args.transcript_isofrac, args.no_coverage_trim)
    };

    // -f (minimum isoform fraction) must be < 1
    if transcript_isofrac >= 1.0 {
        eprintln!(
            "Error: minimum isoform fraction (-f) must be less than 1, got {}",
            transcript_isofrac
        );
        std::process::exit(1);
    }

    let exclude_seqids: Vec<String> = args
        .exclude_seqids
        .as_deref()
        .map(|s| {
            s.split(',')
                .map(|x| x.trim().to_string())
                .filter(|x| !x.is_empty())
                .collect()
        })
        .unwrap_or_default();

    let (relax_filters, pairwise_overlap_filter) = if args.more_transcripts {
        (true, false)
    } else {
        (args.relax_filters, args.pairwise_overlap_filter)
    };

    let trace_intron_chain = if let Some(raw) = &args.trace_intron_chain {
        Some(parse_intron_chain(raw)?)
    } else {
        None
    };

    // the original algorithm -a is "minimum anchor length" for junction support around splice sites.
    // In this codebase, `junction_support` has historically been used as the anchor threshold
    // (including `mismatch_anchor` and the left/right anchor witnesses), while `min_anchor_length`
    // was added later but not wired through. Keep both CLI knobs for now, but treat them as
    // aliases with `-a/--min-anchor-length` taking precedence when it differs from the default.
    let junction_anchor = if args.min_anchor_length != 10 {
        args.min_anchor_length
    } else {
        args.junction_support
    };

    let mut config = RunConfig {
        long_reads,
        long_read_min_len,
        min_junction_reads: args.min_junction_reads,
        junction_correction_window: args.junction_correction_window,
        junction_correction_tolerance: args.junction_correction_tolerance,
        junction_canonical_tolerance: args.junction_canonical_tolerance,
        min_intron_length: args.min_intron_length,
        per_splice_site_isofrac: args.per_splice_site_isofrac,
        bundle_merge_dist: args.bundle_merge_dist,
        // assembler.cpp: runoffdist=200 (not overridden for long reads)
        // print_predcluster: runoffdist=0 only for prediction clustering
        bundle_runoff_dist: args.bundle_runoff_dist,
        cgroup_junction_support: args.cgroup_junction_support,
        min_exon_length: args.min_exon_length,
        junction_support: junction_anchor,
        readthr,
        singlethr: args.singlethr,
        min_transcript_length: args.min_transcript_length,
        min_anchor_length: junction_anchor,
        label: args.label,
        verbose: args.verbose,
        eonly: args.eonly,
        debug_bundle: args.debug_bundle,
        debug_junctions: args.debug_junctions,
        only_debug_bundle: args.only_debug_bundle,
        trace_intron_chain,
        emit_contained_isoforms: args.emit_contained_isoforms,
        emit_junction_paths: args.emit_junction_paths,
        emit_terminal_alt_acceptor: true,
        emit_micro_exon_rescue: true,
        max_junction_paths: args.max_junction_paths,
        max_junction_path_len: args.max_junction_path_len,
        threads: args.threads,
        filter_contained: args.filter_contained,
        transcript_isofrac,
        transcript_isofrac_keep_min: args.transcript_isofrac_keep_min,
        lowisofrac: args.lowisofrac,
        merge_micro_intron_max_gap: args.merge_micro_intron_max_gap,
        output_tpm: args.output_tpm,
        rawreads: args.rawreads,
        ballgown_dir: args.ballgown,
        gene_abundance_out: args.gene_abundance_out,
        gtf_gene_lines: args.gtf_gene_lines,
        min_terminal_exon_len: args.min_terminal_exon_len,
        relax_filters,
        pairwise_overlap_filter,
        pairwise_isofrac: args.pairwise_isofrac,
        pairwise_error_perc: args.pairwise_error_perc,
        pairwise_drop: args.pairwise_drop,
        retained_intron_filter: args.retained_intron_filter,
        retained_intron_error_perc: args.retained_intron_error_perc,
        polymerase_runoff_filter: args.polymerase_runoff_filter,
        polymerase_runon_filter: args.polymerase_runon_filter,
        polymerase_runoff_dist: args.polymerase_runoff_dist,
        polymerase_singlethr: args.polymerase_singlethr,
        dedup_intron_chain_subsets: args.dedup_intron_chain_subsets,
        enable_longtrim: true,
        longtrim_min_boundary_cov: args.longtrim_min_boundary_cov,
        enable_splice_path_extension: args.enable_splice_path_extension,
        strict_longrec_port: args.strict_longrec_port,
        isnascent: args.isnascent,
        allowed_graph_nodes: args.allowed_graph_nodes,
        use_prune_by_junction: args.use_prune_by_junction,
        longintron: rustle::constants::LONGINTRON,
        mismatchfrac: rustle::constants::MISMATCHFRAC,
        lowcov: rustle::constants::LOWCOV,
        sserror: rustle::constants::SSERROR,
        genome_fasta: args.genome_fasta.clone().or(args.cram_ref.clone()),
        merge_bam_compat: args.merge_bam_compat,
        use_graph_merge: !args.no_graph_merge,
        no_coverage_trim,
        exclude_seqids,
        max_sensitivity: args.max_sensitivity,
        debug_stage_tsv: args.debug_stage_tsv.clone(),
        keeptrf_usepath_tsv: args.keeptrf_usepath_tsv.clone(),
        snapshot_jsonl: args.snapshot_jsonl.clone(),
        trace_reference: args.trace_reference.is_some(),
        ref_junction_witness: args.ref_junction_witness,
        vg_mode: args.vg,
        vg_min_shared_reads: args.vg_min_shared,
        vg_em_max_iter: args.vg_em_iter,
        vg_discover_novel: args.vg_discover_novel,
        vg_report: args.vg_report.map(std::path::PathBuf::from),
        vg_solver: args.vg_solver.parse().unwrap_or(rustle::types::VgSolver::Em),
        vg_snp: args.vg_snp,
        vg_phase: args.vg_phase,
        vg_min_novel_reads: args.vg_min_novel_reads,
        vg_discover_novel_mode: args.vg_discover_novel_mode,
        vg_rescue_diagnostic: args.vg_rescue_diagnostic,
        vg_rescue_min_loglik: args.vg_rescue_min_loglik,
        vg_multimap_sequences: std::collections::HashMap::new(),
        vg_mask_regions: {
            let mut out: Vec<(String, u64, u64)> = Vec::new();
            for spec in &args.vg_mask_region {
                let (chrom, range) = match spec.split_once(':') {
                    Some(v) => v,
                    None => anyhow::bail!("--vg-mask-region: expected `chrom:start-end`, got `{}`", spec),
                };
                let (s, e) = match range.split_once('-') {
                    Some(v) => v,
                    None => anyhow::bail!("--vg-mask-region: expected `chrom:start-end`, got `{}`", spec),
                };
                let start: u64 = s.trim().parse()
                    .map_err(|_| anyhow::anyhow!("--vg-mask-region: bad start `{}`", s))?;
                let end: u64 = e.trim().parse()
                    .map_err(|_| anyhow::anyhow!("--vg-mask-region: bad end `{}`", e))?;
                out.push((chrom.trim().to_string(), start, end));
            }
            out
        },
        vg_em_max_copies: args.vg_em_max_copies,
        vg_em_hmm_max_copies: args.vg_em_hmm_max_copies,
        vg_em_skip_intronless: !args.vg_em_no_skip_intronless,
        vg_family_min_shared: args.vg_family_min_shared,
        vg_family_max_copies: args.vg_family_max_copies,
        vg_family_min_shared_per_copy: args.vg_family_min_shared_per_copy,
        vg_family_max_exon_cv: args.vg_family_max_exon_cv,
        vg_family_min_primitive_jaccard: args.vg_family_min_primitive_jaccard,
        vg_family_min_kmer_jaccard: args.vg_family_min_kmer_jaccard,
    };

    if !args.max_sensitivity && (args.compat_preset || config.long_reads) {
        config.apply_compat_preset();
    }
    if config.max_sensitivity {
        config.apply_max_sensitivity();
    }
    // After compat / max-sensitivity so LR junction-merge floors cannot override exact parity.
    config.apply_stringtie_exact_overrides();

    if args.merge {
        if args.merge_bam_compat {
            eprintln!(
                "rustle: --merge-bam-compat is scaffold-only; using existing structural merge mode."
            );
        }
        let mut merge_inputs = args.merge_inputs.clone();
        if let Some(list_path) = &args.merge_list {
            let content = std::fs::read_to_string(list_path)?;
            for line in content.lines() {
                let p = line.trim();
                if p.is_empty() || p.starts_with('#') {
                    continue;
                }
                merge_inputs.push(p.to_string());
            }
        }
        if let Some(positional) = &args.bam {
            merge_inputs.push(positional.clone());
        }
        merge_inputs.sort();
        merge_inputs.dedup();
        if merge_inputs.is_empty() {
            anyhow::bail!("merge mode requires input GTF(s): use --merge-input, --merge-list, or positional input");
        }
        return run_merge(&merge_inputs, output, &config);
    }

    let bam = args
        .bam
        .as_deref()
        .ok_or_else(|| anyhow::anyhow!("bam path required"))?;
    let genome_fasta_owned = config.genome_fasta.clone();
    let genome_fasta_ref = genome_fasta_owned.as_deref();
    run(
        bam,
        output,
        config,
        args.chrom.as_deref(),
        args.bundles.as_deref(),
        args.guide.as_deref(),
        genome_fasta_ref,
        args.trace_reference.as_deref(),
        args.trace_output.as_deref(),
    )?;

    // Explicit flush of parity log (BufWriter inside static OnceLock doesn't auto-flush).
    rustle::parity_decisions::flush();
    Ok(())
}

fn main() -> anyhow::Result<()> {
    run_cli()
}
