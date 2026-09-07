//! De-novo multi-copy family DETECTION + per-read copy ASSIGNMENT CLI.
//!
//! Annotation-free, read-coherence pipeline: scan a region of a coordinate-sorted BAM, assemble de-novo
//! transcripts, detect co-located paralog families, and assign each read — including the hard multimappers
//! minimap2 leaves at MAPQ 0. A family IS one VARIATION GRAPH (copies = PATHS, PSV columns = BUBBLES); each
//! read is THREADED through it and scored to its maximum-likelihood copy-path — by its PSV bases + copy-
//! specific junctions — then significance-gated (assign-or-abstain, never 1/k). `--phase` emits the
//! materialized GFA with the reads threaded through it as walks (the Canzar shared-evidence flip, visualized).
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
use rustle::vg_family::linearize::LinearizeCertificate;
use rustle::vg_family::copy_assign::{AssignParams, AssignStatus};
use rustle::vg_family::em_copy_assign::em_assign_family;
use rustle::vg_family::denovo_assemble::{
    assemble_gate, pass1_skeletons, reads_in_region, tied_secondary_reads_in_region, BamIndexCache, BamRead,
    GATE_MIN_READS,
};
use rustle::vg_family::catalog_input::{
    group_families, parse_copies_fa, parse_copies_tsv, to_colocated, CatalogFamily, SeqIndex,
};
use rustle::vg_family::denovo_pipeline::{
    catalog_overlaps, detect_and_assign, ColocatedFamily, DenovoConfig, FallbackEdge, FamilyAssignment,
    OverlapKind,
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
    /// Primary-alignment MAPQ per read, parallel to `read_names`. Kept (the read sequences are dropped) so the
    /// output stage can compute the tie-break invariance certificate: `mapq > 0` = a unique mapper whose copy
    /// support survives any primary/secondary relabeling. See `anchored_support`.
    read_mapqs: Vec<u8>,
    /// Per-RECORD placement `(ref_start, ref_end, flags)`, parallel to `read_names`, with `flags` bit0 =
    /// secondary and bit1 = supplementary. `read_names` alone names the MOLECULE; this names the
    /// alignment RECORD a family's claim rests on, which is what distinguishes a molecule that SPANS two
    /// loci (one record, claimed by both families) from one placed independently at two disjoint loci
    /// (two records) — the whole of the cross-family reconciliation rule. Kept for the same reason
    /// `read_mapqs` is: the heavy `BamRead` sequences are dropped in the worker.
    read_spans: Vec<(u64, u64, u8)>,
    /// Aligned blocks per record (0-based half-open; `N` closes a block, `D` extends it), parallel to
    /// `read_names`. Kept so the output stage can say whether a read has an aligned BASE inside a copy —
    /// a read spliced OVER a copy is no evidence for it (§6es hygiene; ledger §6cm).
    read_blocks: Vec<Vec<(u64, u64)>>,
    /// Alignment-score evidence per read, parallel to `read_names`. Reported, never decisive — see
    /// `read_conflict::AsEvidence` for why raw AS is length-confounded and `de` makes the call.
    as_ev: Vec<AsEvidence>,
    n_mapped: usize,
    fams: Vec<FamilyAssignment>,
    fallback: Vec<FallbackEdge>,
    dna_needs: Vec<DnaNeedsRecord>,
    /// Augment-and-linearize certificates (Task 4), one per Stage-2-admitted reference-absent copy:
    /// `(family_id, certificate, (chrom, start, end))`. Empty unless the opt-in is on (`--linearize` or
    /// `--linearize-gate`); written to `<out>.linearize.tsv` (Task 5) when set. Under `--linearize-gate`,
    /// a non-LINEARIZES verdict also demotes that candidate out of `admitted` in `detect_and_assign` (it
    /// appears here as a certificate row but not as an admitted copy).
    linearize_certs: Vec<(String, LinearizeCertificate, (String, u64, u64))>,
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
    /// poasta→minimap2 fallback length cap (bp) for the INTRON-RETENTION PSV discovery
    /// (`discover_intron_psvs`, opt-in via `RUSTLE_INTRON_PSV=1`): above this a copy pair's genomic span
    /// uses minimap2, not poasta (which is exact but O(n^2)). Threaded to the use site via `RUSTLE_POA_CAP`
    /// (an env var, not a struct field — the const sits many call frames below `main`; see the module doc at
    /// its use site). Default 20000 matches the prior hard-coded constant exactly, so leaving this flag
    /// unset is byte-identical to before it existed.
    #[arg(long, default_value_t = 20_000)]
    poa_cap: usize,
    /// Per-family multimapper read-pool cap (`o2_materialize::READ_CAP`, `MaterializeConfig::read_cap`).
    /// EXPOSED HERE FOR AUDITABILITY ONLY: `o2_materialize` is a Rust byte-parity port of the Python
    /// genome-wide-catalog materializer (`bench/o2_vg_visualization.py::materialize_family`) that no
    /// `src/bin/*.rs` binary — including this one — imports, so this flag is currently a NO-OP in
    /// `copy_assign` (parses so `RUSTLE_READ_CAP`/CLI usage never hard-errors; a non-default value warns at
    /// startup rather than silently doing nothing). Default 6000 matches the constant.
    #[arg(long, default_value_t = 6_000)]
    read_cap: usize,
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
    /// One-flag IGV bundle: implies `--dump-psv` (the PSV genotype matrix), so a subsequent
    /// `bench/igv_tracks.py --assignments <out>.assignments.tsv --bam <bam> --regions <regions> --out <out>`
    /// emits `<out>.tagged.bam` (reads coloured by assigned copy), `<out>.copies.bed`, and `<out>.psv.vcf`
    /// (PSVs as an IGV variant track, copies as samples) — everything IGV needs to SEE each assignment.
    #[arg(long, default_value_t = false)]
    igv: bool,
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
    /// ⭐ O2-8a (§6ej): abstain on junction/PSV CONFLICT — a read whose own splice junctions fit another copy
    /// strictly better than its PSV-best copy is `ambiguous`, never `assigned` (on NPIP's LCR16a cores 4 of 11
    /// junction-anchored assignments contradicted the junction at min_p 1e-16..1e-35). Writes
    /// `<out>.conflicts.tsv`. Likelihoods untouched; default OFF ⟹ byte-identical
    #[arg(long, default_value_t = false)]
    junction_conflict_abstain: bool,

    /// ⭐ §6eu / register 689: the read-support PSV filter (`read_supported_columns`: keep a copy-vs-copy column
    /// only if two alleles each reach two reads). It cannot tell a mis-assembled base from an UNEXPRESSED
    /// paralogue — both are monomorphic pileups — and on ZSCAN5 it kept 4 of 216 columns and produced 11
    /// confident wrong calls. **Default OFF (user decision 2026-09-05)**; on the 35-family sweep OFF raised
    /// assignments 18 % and MAPQ-60 agreement 95.3 → 96.0 %. `--psv-read-filter` turns it back on;
    /// `RUSTLE_PSV_READFILTER=0|1`, when set, overrides the flag (the pre-2026-09-05 escape).
    #[arg(long, default_value_t = false)]
    psv_read_filter: bool,

    /// ⭐ O2-9 / D3 (PREREG adj/d3, register 691): pool every BAM record of a molecule into ONE observation
    /// vector — "the read is the star" (§6fa): the molecule's sequence aligned to every copy's unit, its own
    /// columns, an origin certificate — instead of scoring each BAM record and abstaining on contradiction
    /// (which masked wrong column positions, rows 700–702). **Default ON (user decision 2026-09-05)**;
    /// `--no-molecule-observations` restores the record-level path byte-for-byte.
    #[arg(long, default_value_t = true)]
    molecule_observations: bool,
    /// Escape hatch: the record-level observation path (every catalog before 2026-09-05 §6fa).
    #[arg(long, default_value_t = false)]
    no_molecule_observations: bool,
    /// ⭐ §6fc: the origin certificate counts SUBSTITUTIONS only (indels = isoform structure, not origin). Higher
    /// yield (NPIP 175 → 442 assigned) at 1 wrong of 20 audited anchors and 98.6 % MAPQ-60 agreement; the
    /// default counts every edit and made no wrong anchor call. Opt-in.
    #[arg(long, default_value_t = false)]
    origin_substitutions_only: bool,
    /// ⭐ §6fc: splice junctions as pairwise evidence in read-star (opt-in: +4 points of assignment at −0.3 points
    /// of MAPQ-60 agreement on the paired 35 families).
    #[arg(long, default_value_t = false)]
    read_star_junctions: bool,
    /// ⭐ §6fd: read-star against each candidate's GENOMIC locus (unit extent padded by the family's longest
    /// read), splice-aware; the origin certificate counts X + I + D + unaligned read bases over the read length.
    /// **Default ON (user decision 2026-09-05; PREREG adj/gstar held: 0 wrong anchors, 100 % placement
    /// agreement on the paired 35)**. `--read-star-unit` restores the spliced-unit form of §6fc (10× faster).
    #[arg(long, default_value_t = true)]
    read_star_genomic: bool,
    /// Escape hatch: read-star against the spliced UNITS (§6fc), the default before 2026-09-05 §6fd.
    #[arg(long, default_value_t = false)]
    read_star_unit: bool,
    /// Escape hatch (L2): ignore the catalog's `locus_start`/`locus_end` columns and pad each unit's extent by
    /// the family's longest molecule (the §6fd rule). Without the columns both arms are byte-identical.
    #[arg(long, default_value_t = false)]
    read_star_pad_locus: bool,
    /// Escape hatch (L3, §6fi): report every single-candidate molecule as `tied` (the §6fa rule) instead of
    /// `assigned` with `sole_candidate = 1` when its only candidate passes the origin certificate.
    #[arg(long, default_value_t = false)]
    no_sole_candidate: bool,
    /// L6: write `<out>.star_reads.tsv` — per molecule its read-star proof: the read positions of its columns,
    /// its base there and every candidate's base (the assignment-proof figure's source). Default OFF.
    #[arg(long, default_value_t = false)]
    dump_star: bool,
    /// Escape hatch (§6fm): count every read-star hit for a candidate even when it does not overlap the
    /// candidate's unit inside the target (the §6fd behaviour).
    #[arg(long, default_value_t = false)]
    no_read_star_hit_in_unit: bool,
    /// Escape hatch (§6fp): certify against the genomic locus alone (the §6fd form) instead of the better of
    /// locus and expressed chain per candidate.
    #[arg(long, default_value_t = false)]
    read_star_genomic_only: bool,
    /// Escape hatch (§6fq): do NOT assign uncontested molecules (primary MAPQ ≥ 60) to their placement; run the
    /// certificate machinery on every molecule as before.
    #[arg(long, default_value_t = false)]
    no_placement_assign: bool,
    /// §6fq variant: the placement OVERRIDES the machinery's certified call on uncontested molecules (the pure
    /// "any assembler would" rule). Default is certified-first with the placement as the fallback.
    #[arg(long, default_value_t = false)]
    placement_first: bool,
    /// ⭐ O2-8c (§6eo): discover PSV columns on the GENOMIC alignment of the copies' spans (exons + introns,
    /// reverse-complement retry for inverted duplications) instead of their spliced sequences. Read-chain units
    /// of unequal exon composition sent the spliced star projection to min_p 3e-270 on a wrong call (register
    /// 683); SEDEF core hulls are collinear genomic segments. Default OFF ⟹ byte-identical
    #[arg(long, default_value_t = false)]
    psv_genomic: bool,
    /// εⱼ used for an editing-flagged PSV column in the certificate (the rate a base shows the other allele
    /// by editing rather than sequencing error). Default 0.2.
    #[arg(long, default_value_t = 0.2)]
    edit_rate: f64,
    /// IsoCon-style iterative copy pruning: after assignment, repeatedly merge copies that have no read with
    /// significant evidence distinguishing them from their nearest neighbor, reassigning reads until all
    /// surviving copies are defensible. Default OFF (byte-identical baseline).
    #[arg(long, default_value_t = false)]
    iterative_prune: bool,
    /// RE-EXPRESS THE COPY ASSIGNMENT IN PHASING VOCABULARY (PS/HP tags + a GFA). This runs NO phasing
    /// algorithm of its own: it is a RELABELING of the assignment this binary already computed —
    /// phase set = family, haplotype = the assigned copy, unphased (-1) = abstained (ambiguous/tied).
    /// Measured: over 6 historical runs / 119,524 read-rows the multiset of
    /// `(read_name, family, haplotype)` equals `(read_name, family, assigned_copy-if-assigned-else -1)`
    /// EXACTLY, symmetric difference 0. (Compare per read NAME and you get a spurious mismatch — read
    /// names repeat within a family; the identity is only visible as a multiset.)
    /// Writes `<out>.phase_blocks.tsv` (one PHASE SET per family), `<out>.phased_haplotypes.tsv` (each
    /// copy's `pos:allele` PSV string), `<out>.phased_reads.tsv` (read → copy, with `n_psv_spanned` =
    /// `n_decisive` and `margin` = the assignment's log-LR margin), and `<out>.phase.gfa`.
    /// Block ≙ PS tag, haplotype ≙ HP tag, so downstream phasing tooling can consume it.
    /// ⚠ HISTORICAL NOTE: this help previously claimed a "min-path-cover over the PSV graph". No such
    /// computation exists in the O2 path — see `docs/copy_assignment_definition.md` §9.4.
    #[arg(long, default_value_t = false)]
    phase: bool,

    /// Optional gene annotation (GFF3/GTF/BED) to tag in-genome copies annotated vs unannotated in the
    /// --phase copy graph. Without it, in-genome copies are tagged `annotation-unknown`.
    #[arg(long)]
    gff: Option<String>,

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

    /// EXPERIMENT-ONLY override of the reference-absent admission gate's cluster floor
    /// (`absent_copy.rs` gate 1, `n_clusters >= min_clusters`). Unset = the shipping default 3, and
    /// output is byte-identical to a build without this flag.
    ///
    /// Gate 1 is an identifiability claim, NOT a tuning knob: below three co-varying clusters a
    /// second COPY is indistinguishable from a heterozygous ALLELE without DNA copy-number data.
    /// Lowering it in production manufactures false positives. The one legitimate use is a
    /// removal-recovery ablation (V4b), where a known copy is DELETED from the assembly, so copy
    /// status holds BY CONSTRUCTION and the gate is re-asking a question the design already
    /// answered — and where it is unreachable anyway (deleting one copy of a 3-copy family leaves
    /// at most 2 clusters). Such a run is interpretable ONLY beside the identical INTACT-assembly
    /// control at the same value: if the control also admits a copy, the recovery is an artefact.
    ///
    /// A value that is not a positive integer falls back to 3 (a typo must not disable the gate).
    #[arg(long, value_name = "N")]
    absent_min_clusters: Option<usize>,

    /// Opt-in augment-and-linearize REPORT (requires --absent-copies; no-op otherwise): compute a
    /// `LinearizeCertificate` for every Stage-2-admitted reference-absent candidate and write it to
    /// `<out>.linearize.tsv`. This costs one minimap2 realign-pool subprocess per admitted candidate, so it
    /// is OFF by default — plain `--absent-copies` keeps its prior admission cost and emits no certificate
    /// file (byte-identical to the pre-feature `--absent-copies` output). `--linearize-gate` implies this.
    #[arg(long, default_value_t = false)]
    linearize: bool,

    /// Opt-in augment-and-linearize GATE (implies --linearize; requires --absent-copies; no-op otherwise): a
    /// Stage-2 candidate whose `LinearizeCertificate` verdict is NOT `Linearizes` (its MAPQ-0 read pool does
    /// not land on it distinguishably more often than on a dinucleotide-shuffled decoy) is DEMOTED — written
    /// to `<out>.dna_needs.tsv` instead of admitted as a copy. Because it turns the certificate into an
    /// admission decision it also enables the certificate computation + `<out>.linearize.tsv` report (as if
    /// `--linearize` were set). Default OFF: admission is unchanged (byte-identical to plain --absent-copies).
    #[arg(long, default_value_t = false)]
    linearize_gate: bool,

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

    /// Read lambda_global from a `gw_family_catalog --single-copy-baseline` `<prefix>.lambda_global.tsv`
    /// instead of passing the scalar by hand. `--lambda-global <f64>` (if given) takes precedence. This makes
    /// the copy-number normalizer an in-binary quantity rather than an external-script number.
    #[arg(long)]
    lambda_file: Option<String>,

    /// VG re-align supplement (opt-in): for every co-located family, re-align each poor-fit/candidate
    /// read (low MAPQ, heavy clipping, or high divergence — `vg_realign::is_candidate`) to the family's
    /// copy-paths and record the decision (`reassigned` / `rejected` / `novel-candidate`) to
    /// `<out>.vg_realign.tsv`. FEEDS BACK into the assignment: BOTH the correction leg (re-thread + reassign
    /// among existing copies) AND the admission leg (admit novel-read pools as NEW copies — genome-touching,
    /// widens the roster). Default off; when off, every output is byte-identical. Use `--vg-realign-correct`
    /// for the correction leg WITHOUT the FP-risk admission.
    #[arg(long, default_value_t = false)]
    vg_realign: bool,

    /// VG re-align CORRECTION leg ONLY (opt-in): re-thread hard reads through the family copy-paths and
    /// correct their assignments among the EXISTING copies, WITHOUT admitting novel copies (no genome touch,
    /// no roster widening). The safe VG-native assignment leg. Off by default; byte-identical when off.
    #[arg(long, default_value_t = false)]
    vg_realign_correct: bool,

    /// Define family MEMBERSHIP by E_r transcript homology instead of the E_c read-conflict graph. The
    /// conflict graph links two copies only when reads map ambiguously between them, so a copy whose reads
    /// all map uniquely is dropped from its family and its reads come back `tied` — not because they are
    /// unassignable, but because their true copy was never admitted. Conflict, PSVs, and chi(H) remain
    /// within-family. Admitting a dropped copy enlarges the copy set, so the Bonferroni certificate
    /// alpha/(K-1) tightens and existing assignments shift; this is why the mode is opt-in. Requires
    /// minimap2 (honors RUSTLE_MINIMAP2); aborts rather than falling back to the conflict graph.
    #[arg(long, default_value_t = false)]
    homology_primary: bool,

    /// DISABLE the mutual-homology + distinct-locus family gate. Refinement is ON BY DEFAULT: each co-located
    /// family must have its copies MUTUALLY HOMOLOGOUS (the shared E_r primary tier -- `-k 11 -w 5` at
    /// `sensitive_identity` 0.60 by default since X.4, `-x asm20` at 0.80 under `RUSTLE_ER_SENSITIVE_ONLY=0`
    /// -- cov-of-shorter>=0.50)
    /// across >= 2 distinct loci.
    /// ⚠ This is NOT "the same criterion `gw_family_catalog` refines by" and the two paths do NOT
    /// automatically agree — that claim was false from D1 (2026-08-09) onward and is corrected here
    /// (O-4). `gw_family_catalog`'s DEFAULT homology catalog does not call refine at all
    /// (`refine_enabled` = `refine_flag || !o1_homology`), so it applies no such gate; and where refine
    /// DOES run it additionally UNIONS a genomic-span tier in (see `additive_genomic_tier`), which the
    /// catalog's own E_r site never does. Both facts are certified per run: `additive_genomic_tier` in
    /// `<prefix>.rule.tsv` and `n_edges_genomic_tier_added` in `<prefix>.refine.params.tsv`.
    /// Without this gate the conflict oracle admits large-gene mis-chains (PBX1) and
    /// repeat-bridges as families (`bench/GW_CATALOG_FP_AUDIT.md`). `--no-refine` assigns the raw families and
    /// needs no minimap2.
    #[arg(long, default_value_t = false)]
    no_refine: bool,

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

    /// NO-OP in copy_assign: this flag parses here but has no consumer in this binary. The
    /// collapse-enumeration gate (re-admitting near-identical <2-RNA-loci families as K0_COLLAPSED
    /// copy NUMBER, `<out>.collapsed.tsv`) runs only in `gw_family_catalog` -- use `--collapse-enumerate`
    /// there instead. Kept here only so `RUSTLE_COLLAPSE_ENUMERATE`/CLI parsing doesn't hard-error;
    /// NOT wired into copy_assign's per-read assignment path (that would violate the
    /// COPY-NUMBER-ONLY contract of this feature).
    #[arg(long, default_value_t = false)]
    collapse_enumerate: bool,

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

    /// Seed candidate loci from AS-tied secondary reads that share an intron chain, even with no primary
    /// (recovers covered-but-tied K=0 copies as detected-but-unassignable). Implies fetching tied secondaries.
    #[arg(long, default_value_t = false)]
    tied_seed: bool,

    /// CONSUME the O1 catalog: a `gw_family_catalog` `<out>.copies.tsv`. Those rows ARE the copy set for the
    /// swept regions — this binary then detects NOTHING. Reads still come from `--bam`; only the copy set is
    /// supplied.
    ///
    /// # Why
    /// O1 and O2 share one node type, one edge engine and one admission primitive BY FUNCTION CALL and
    /// nothing BY FILE, so each binary re-derived its own families and the two tables had no join key
    /// (`GWFAM{i}` vs `CAFAM{i}`, assigned independently). Measured at defaults, they built different
    /// objects: the GSTM catalog has 4 copies where `copy_assign` found 0 families on 6031 reads. With
    /// `--families` the two copy sets agree BY CONSTRUCTION, and every emitted row carries the catalog's own
    /// `family_id` (so `family_id` here IS `GWFAM{i}`, not `CAFAM{i}`) plus a `<out>.family_join.tsv`
    /// mapping each assigned copy back to its `(family_id, copy_idx, tid)` catalog row.
    ///
    /// # What is switched OFF (documented, not silent)
    /// Family CONSTRUCTION, in full: pass-1 skeletons, the assemble gate, locus collapse, the membership
    /// oracle (E_c conflict graph / E_r homology), the POA diagnostic, co-location, the REFINE gate
    /// (`--no-refine` becomes moot), and the thin-locus RESCUE leg (it ADDS copies below the assembly read
    /// floor). Anything else that would change the roster is REFUSED rather than silently applied:
    /// `--absent-copies`, `--vg-realign` (admission leg), `--iterative-prune`, `--collapse-gate`,
    /// `--tied-seed`, `--recover-copies`. `--vg-realign-correct` is allowed (it only re-threads reads among
    /// the given copies). `--min-copies`/`--win` are not applied: the catalog already decided membership.
    ///
    /// # Contract (all loud, none silent)
    /// Every supplied copy must (a) be named by the `copies.tsv` header columns, (b) belong to a
    /// SAME-CHROMOSOME family (a cross-chrom family — RABL2's 5 contigs — is structurally unassignable by a
    /// region-scoped binary and is refused, never truncated), (c) fall inside exactly one `--region` /
    /// `--regions` entry, (d) have a sequence (see `--copies-fa`), and (e) have at least one overlapping
    /// read in the BAM. A violation aborts the run.
    #[arg(long)]
    families: Option<String>,

    /// The `gw_family_catalog` `<out>.copies.fa` beside `--families`. When given, each copy's spliced
    /// sequence is the catalog's OWN emitted bytes (checked against its `copies.tsv` row: chrom, span,
    /// strand, exon count) — so O2 assigns against exactly the sequence O1 defined the family with, with no
    /// reconstruction step that could differ. Without it the sequence is rebuilt from `--fasta` at the
    /// catalog's exon coordinates through the SAME `build_spliced_seq` the catalog used, and the strand it
    /// derives from the junction motifs must match the strand the catalog recorded (a mismatch means the
    /// FASTA is not the assembly the catalog was built against, and aborts). Prefer `--copies-fa`.
    #[arg(long)]
    copies_fa: Option<String>,
}

fn status_str(s: AssignStatus) -> &'static str {
    match s {
        AssignStatus::Assigned => "assigned",
        AssignStatus::Ambiguous => "ambiguous",
        AssignStatus::Tied => "tied",
    }
}

// ---- CROSS-FAMILY RECONCILIATION (RUSTLE_XFAM_RECONCILE) -------------------------------------------
//
// THE DEFECT. Every family is genotyped INDEPENDENTLY, and a molecule whose records fall in two
// families' read pools is put through both significance gates with no communication between them. On
// `mec` (12 supplied families, 2 regions, one NPIP-region BAM) that produces 79,175 assignment rows over
// 77,372 distinct molecules; 1,587 molecules are genotyped in >=2 families and 517 come back `assigned`
// in >=2 families at once, for 519 assigned-copy pairs. 210 of those pairs name copy intervals that are
// DISJOINT (median separation 16,341 bp, min 4,797 bp) — and a single molecule cannot originate from two
// disjoint loci, so at least one member of each such pair is wrong. Nothing downstream reconciled them.
// MEASURED split of the 519 (report arm, in-binary and authoritative): 309 `shared_locus`,
// 99 `readthrough_span`, 111 `cross_family_contradiction`; abstaining demotes 221/53,715 = 0.0041 of
// assigned rows, 110 molecules, and 221/1,035 = 0.2135 of CONTESTED rows — the third rate being the one
// that matters, since a rule that only fires on multi-placement molecules is selection-biased by
// construction. ⚠ ONE BAM, 2 regions, 6 family pairs, 170/210 disjoint pairs involving GWFAM111 alone
// and all 309 overlaps from 2 copy pairs: these COUNTS are not claimed to transport, only the mechanism.
//
// THE RULE, in three strata, and only the third one abstains:
//
//   `shared_locus`             the two claimed copy intervals OVERLAP. Two nominally different families
//                              claim one locus. That is an O1 PARTITION artifact (the same thing
//                              `catalog_overlaps` already warns about, below), not an O2 assignment
//                              error: the molecule really does come from that one interval, and both
//                              families name it. Reported, never demoted — demoting it would charge an
//                              O1 defect to O2's abstention rate and gut real copies of their support
//                              (MEASURED on `mec`, double-claimed/assigned: GWFAM113 copy0 141/141,
//                              GWFAM111 copy1 168/179, GWFAM112 copy12 141/153, GWFAM96 copy11
//                              168/234 — the same two flagged copy pairs `catalog_overlaps` reports).
//   `readthrough_span`         the two copies are disjoint but ONE alignment record is claimed by both —
//                              a molecule whose N-gap spans both loci. "A molecule cannot come from two
//                              disjoint loci" does not apply to a molecule that demonstrably spans both,
//                              so this stratum is reported, never demoted.
//   `cross_family_contradiction`  the two copies are disjoint AND the two claims rest on DIFFERENT
//                              alignment records. Two independent placements naming disjoint loci: the
//                              molecule has one origin, so at least one claim is false and there is no
//                              admissible way to tell which. This is exactly O2's assign-or-abstain
//                              contract at cross-family scope, so under `abstain` the molecule is
//                              demoted to `Ambiguous` in EVERY family it is assigned in.
//
// WHY NO ARBITRATION (a "keep the better one" arm was never built, deliberately). At cross-family scope
// none of the certificate fields is comparable (all MEASURED on `mec`, 519 contested pairs over 12
// families): the median assigned margin per family spans 9.9 to 10,745.2 (1,085x), `p_value` is 0.0 on
// BOTH sides of 434/519 = 0.8362 of contested pairs and is gated against a family-size-dependent
// `alpha/(n-1)`, and the margin winner is simply the larger-`n_decisive` side in 385/519 = 0.7418 (with
// 0 n_decisive ties). AS is byte-identical on both rows in 519/519, so an AS tie-break IS minimap2's
// primary flag in disguise — the defect that retired `uniq_agree` — and the primary flag itself agrees
// with the margin winner in only 61/109 = 0.5596 of the pairs where exactly one side is primary
// (chance, and underpowered). Abstention is
// therefore the whole of the intervention: it demotes ALL of a contradicting molecule's claims, because
// choosing which side to keep IS the arbitration being refused.
//
// WHAT THIS DOES NOT FIX (stated here and repeated above the `.quant.tsv` write). A demotion changes
// STATUS only. `n_reads_hard` counts `fa.assignments` by argmax `best_copy` with NO status filter, and
// `abundance`/`ci95` come from `soft_quantify_em` inside the per-family pipeline whose `obs_for_em` is
// populated regardless of status — so `.quant.tsv` is byte-identical between `report` and `abstain` and
// the row inflation there is untouched. Removing a molecule from a family's EM is a two-pass
// architecture change and a separate decision. Anyone claiming this fixes the abundance double-count is
// wrong.
//
// DEFAULT OFF. Under `off` (or unset) the detector does not run, nothing new is written beyond the
// unconditional `<out>.params.tsv`, `demote` is empty and every status emit site returns `status_str`
// verbatim — the OFF arm is byte-identical BY CONSTRUCTION, not by inspection.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum XfamMode {
    /// unset / `off`: the detector does not run and no side file is written (byte-identical arm).
    Off,
    /// `report`: write `<out>.xfam_conflicts.tsv` + a stderr summary; change no existing byte.
    Report,
    /// `abstain`: additionally demote every claim of a `cross_family_contradiction` molecule to `Ambiguous`.
    Abstain,
}

impl XfamMode {
    /// The EFFECTIVE value recorded in `<out>.params.tsv`. An unrecognized value is an error, not a
    /// silent `off`: a typo'd flag that quietly disables the pass is how an ON arm gets reported as OFF.
    fn from_env() -> Result<Self> {
        match std::env::var("RUSTLE_XFAM_RECONCILE").ok().as_deref() {
            None | Some("") | Some("off") => Ok(XfamMode::Off),
            Some("report") => Ok(XfamMode::Report),
            Some("abstain") => Ok(XfamMode::Abstain),
            Some(other) => anyhow::bail!(
                "RUSTLE_XFAM_RECONCILE={other:?} is not one of off|report|abstain"
            ),
        }
    }
    fn as_str(self) -> &'static str {
        match self {
            XfamMode::Off => "off",
            XfamMode::Report => "report",
            XfamMode::Abstain => "abstain",
        }
    }
}

/// A REGION-INDEPENDENT alignment-record key: `(contig, ref_start, ref_end, flags)` with
/// `flags` bit0 = secondary, bit1 = supplementary.
///
/// `fa.assignments` already carries the region-global `bam_reads` index (`idx_map[r.read_index]` in
/// `denovo_pipeline`), so the record a family's claim rests on is known exactly — but two REGIONS have
/// independent index spaces, so the index alone cannot be compared across them. This key can. Two records
/// at identical contig/start/end/flags are the same placement; collapsing them is correct and
/// conservative (it can only move a pair OUT of `cross_family_contradiction`).
type RecKey = (String, u64, u64, u8);

/// One family's `assigned` claim on one molecule: which (region, family) made it, which copy interval it
/// names, and which alignment record it rests on.
struct XfamClaim<'a> {
    g: usize,
    f: usize,
    fid: &'a str,
    copy: usize,
    span: (String, u64, u64),
    rec: RecKey,
}

/// One contested pair, as written to `<out>.xfam_conflicts.tsv`.
struct XfamConflict {
    read_name: String,
    stratum: &'static str,
    fid_a: String,
    copy_a: usize,
    span_a: (String, u64, u64),
    fid_b: String,
    copy_b: usize,
    span_b: (String, u64, u64),
    same_record: bool,
    sep_bp: u64,
    demoted: bool,
}

/// PASS 1: read-only sweep over every region's assignments, classifying every molecule claimed
/// `assigned` by two or more families.
///
/// Deterministic by construction: `BTreeMap`/`BTreeSet` keyed on `(read_name, region ordinal, family
/// ordinal)` with insertion-ordered inner vectors, and no `HashMap` iteration anywhere — so the output
/// does not depend on region-thread scheduling.
fn xfam_pass1(
    works: &[RegionWork],
    named_families: bool,
) -> (Vec<XfamConflict>, std::collections::BTreeSet<(String, usize, usize)>) {
    use std::collections::{BTreeMap, BTreeSet};
    // The family ids the DRAIN will mint, recomputed here with the same counter so the side file names
    // families exactly as `.assignments.tsv` does. Declared before `by_read` so it outlives the borrow.
    let mut gfam = 0usize;
    let mut fids: Vec<Vec<String>> = Vec::with_capacity(works.len());
    for work in works.iter() {
        let mut here = Vec::with_capacity(work.fams.len());
        for fa in &work.fams {
            here.push(if named_families { fa.family_id.clone() } else { format!("CAFAM{gfam}") });
            gfam += 1;
        }
        fids.push(here);
    }
    let mut by_read: BTreeMap<&str, Vec<XfamClaim>> = BTreeMap::new();
    for (g, work) in works.iter().enumerate() {
        for (f, fa) in work.fams.iter().enumerate() {
            for (ri, a) in &fa.assignments {
                if !matches!(a.status, AssignStatus::Assigned) {
                    continue;
                }
                let Some(span) = fa.copy_spans.get(a.best_copy).cloned() else { continue };
                let (rs, re, fl) = work.read_spans.get(*ri).copied().unwrap_or((0, 0, 0));
                by_read.entry(work.read_names[*ri].as_str()).or_default().push(XfamClaim {
                    g,
                    f,
                    fid: fids[g][f].as_str(),
                    copy: a.best_copy,
                    span,
                    rec: (work.contig.clone(), rs, re, fl),
                });
            }
        }
    }
    let mut conflicts: Vec<XfamConflict> = Vec::new();
    let mut demote: BTreeSet<(String, usize, usize)> = BTreeSet::new();
    for (name, claims) in by_read.iter() {
        // ">= 2 DISTINCT (region, family)" is the contest condition. There is at most one row per
        // (read, family) today, but the guard is stated on the key, not on the count.
        let distinct: BTreeSet<(usize, usize)> = claims.iter().map(|c| (c.g, c.f)).collect();
        if distinct.len() < 2 {
            continue;
        }
        let mut contradicts = false;
        let mut here: Vec<XfamConflict> = Vec::new();
        for i in 0..claims.len() {
            for j in (i + 1)..claims.len() {
                let (a, b) = (&claims[i], &claims[j]);
                if (a.g, a.f) == (b.g, b.f) {
                    continue; // same family: the intra-family reduction already owns this case
                }
                let overlaps =
                    a.span.0 == b.span.0 && a.span.1 < b.span.2 && b.span.1 < a.span.2;
                let same_record = a.rec == b.rec;
                let stratum = if overlaps {
                    "shared_locus"
                } else if same_record {
                    "readthrough_span"
                } else {
                    contradicts = true;
                    "cross_family_contradiction"
                };
                let sep_bp = if overlaps || a.span.0 != b.span.0 {
                    0
                } else {
                    a.span.1.max(b.span.1) - a.span.2.min(b.span.2)
                };
                here.push(XfamConflict {
                    read_name: (*name).to_string(),
                    stratum,
                    fid_a: a.fid.to_string(),
                    copy_a: a.copy,
                    span_a: a.span.clone(),
                    fid_b: b.fid.to_string(),
                    copy_b: b.copy,
                    span_b: b.span.clone(),
                    same_record,
                    sep_bp,
                    demoted: false,
                });
            }
        }
        if contradicts {
            // CONTRADICTION IS A PROPERTY OF THE MOLECULE, so every claim it makes abstains — not only
            // the two in the contradicting pair. Abstaining one side would be arbitration, which is
            // exactly what is refused here.
            for c in claims {
                demote.insert(((*name).to_string(), c.g, c.f));
            }
            for c in here.iter_mut() {
                c.demoted = true;
            }
        }
        conflicts.extend(here);
    }
    (conflicts, demote)
}

/// GFA W-line SampleId sanitizer (a walk id must be a whitespace-free GFA token).
fn sanitize_gfa_id(s: &str) -> String {
    s.chars().map(|c| if c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-' | ':') { c } else { '_' }).collect()
}

/// Per-copy status shared by `build_copy_graph` (v1, PSV-bubble graph) and `build_exon_graph` (v2,
/// exon presence/absence graph) — the single source of truth for the in-genome/absent axis so the two
/// graph builders can never disagree about a copy's status.
///
/// A copy is ABSENT iff at least one read assigned to it is `discovery_coupled` — i.e. the copy exists only
/// because absent-copy discovery (`--absent-copies`) admitted it. This is the EXACT v1 rule, unchanged:
/// `copy_map_identity` feeds the `MI:f:` tag ONLY and must NEVER drive the absent/ST status (a copy admitted
/// via `--absent-copies` carries a remap identity but is not absent unless a `discovery_coupled` read pins
/// it). (The `collapsed_copies`/`rescued_copies` fields are diagnostic COUNTS, NOT per-`copy_tids` absence
/// markers — `collapsed_copies` can exceed `n_copies` — so they must NOT drive absence either.) An absent
/// copy is `AbsentCollapsed` when its genomic span OVERLAPS a non-absent (in-genome) copy of the same
/// family — a hidden CO-LOCATED haplotype — else `AbsentDivergent` (dispersed, or spans unavailable).
/// Non-absent in-genome copies are tagged by `annotation_status`: `InGenomeAnnotated`/`InGenomeUnannotated`
/// when `--gff` was given, else `AnnotationUnknown` (we never claim "unannotated" unchecked).
fn copy_status(
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    ci: usize,
    ann: Option<&[(String, u64, u64)]>,
) -> rustle::vg_family::copy_graph::CopyStatus {
    use rustle::vg_family::copy_graph::CopyStatus;
    let n = fa.copy_tids.len();
    let is_coupled =
        |k: usize| -> bool { fa.assignments.iter().any(|(_, a)| a.discovery_coupled && a.best_copy == k) };
    if !is_coupled(ci) {
        return annotation_status(fa, ci, ann);
    }
    // Half-open genomic-span overlap of two copies (same chrom): distinguishes an absent copy that is a
    // hidden CO-LOCATED haplotype (AbsentCollapsed) from a dispersed one (AbsentDivergent).
    let span_overlap = |i: usize, j: usize| match (fa.copy_spans.get(i), fa.copy_spans.get(j)) {
        (Some((ci, si, ei)), Some((cj, sj, ej))) => ci == cj && si < ej && sj < ei,
        _ => false,
    };
    // AbsentCollapsed iff this absent copy's span overlaps some IN-GENOME (non-coupled) copy — a hidden
    // co-located haplotype; else it is dispersed / unlocalized => AbsentDivergent.
    let collapsed = (0..n).any(|k| k != ci && !is_coupled(k) && span_overlap(ci, k));
    if collapsed { CopyStatus::AbsentCollapsed } else { CopyStatus::AbsentDivergent }
}

/// Build one family's copy-graph (`--phase`): the REFERENCE walk + every copy as a tagged, corroborable
/// PATH + every read as a threaded WALK, over the family's usable PSV columns (those with both a genome
/// position and a fetchable reference base). Pure w.r.t. the caller except for `ref_base` (injected so this
/// is unit-testable without genome I/O — see `build_copy_graph_maps_family_to_graph`). Per-copy status is
/// `copy_status` (shared with `build_exon_graph` — see its doc for the absence rule).
/// `eff` is the EFFECTIVE per-row status, parallel to `fa.assignments` — `a.status` verbatim unless the
/// cross-family reconciliation demoted that row (see `XfamMode`). Passed in rather than read off
/// `a.status` so the graph's `Assigned` filters cannot disagree with `.assignments.tsv`.
fn build_copy_graph(
    fid: &str,
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    ref_base: impl Fn(&str, u64) -> Option<u8>,
    bam_reads: &[String],
    ann: Option<&[(String, u64, u64)]>,
    eff: &[AssignStatus],
) -> rustle::vg_family::copy_graph::CopyGraph {
    use rustle::vg_family::copy_graph::*;
    // usable columns (both a genome position and a reference base), remembering the original index.
    let mut cols: Vec<PsvColumn> = Vec::new();
    let mut keep: Vec<usize> = Vec::new();
    for (j, p) in fa.psv_col_pos.iter().enumerate() {
        if let Some(pos) = p {
            if let Some(rb) = ref_base(&fa.chrom, *pos) {
                cols.push(PsvColumn { col: j, genome_pos: Some(*pos), ref_allele: Some(rb) });
                keep.push(j);
            }
        }
    }
    let sel = |row: &Vec<Option<u8>>| keep.iter().map(|&j| row.get(j).copied().flatten()).collect::<Vec<_>>();
    let backbone = vec![b"NNNNNNNNNN".to_vec(); cols.len() + 1];

    let n = fa.copy_tids.len();
    let copies = (0..n)
        .map(|ci| {
            let status = copy_status(fa, ci, ann);
            let reads = fa
                .assignments
                .iter()
                .zip(eff.iter())
                .filter(|((_, a), s)| a.best_copy == ci && matches!(s, AssignStatus::Assigned))
                .count() as u32;
            CopyPath {
                id: format!("{}_copy{}", fid, ci),
                alleles: fa.copy_psv_alleles.get(ci).map(|r| sel(r)).unwrap_or_default(),
                status,
                corrob: Corrob { reads: Some(reads), suns: None, map_identity: fa.copy_map_identity.get(ci).copied().flatten() },
            }
        })
        .collect();

    let reads = fa
        .assignments
        .iter()
        .zip(fa.read_psv_obs.iter())
        .enumerate()
        .map(|(k, ((ri, a), obs))| {
            let st = eff.get(k).copied().unwrap_or(a.status);
            ReadWalk {
                name: sanitize_gfa_id(&format!("{}_{}", fid, bam_reads[*ri])),
                obs: sel(obs),
                assigned_copy: if matches!(st, AssignStatus::Assigned) { Some(a.best_copy) } else { None },
                cert: Some(ReadCert { p_value: a.p_value, min_p_value: a.min_p_value, status: st }),
            }
        })
        .collect();

    CopyGraph { family: fid.to_string(), columns: cols, backbone, copies, reads }
}

/// Build one family's exon presence/absence graph (`--phase` v2): reconstructs each copy's genomic exon
/// chain from `copy_introns` + `copy_spans` (donor/acceptor intron chain -> exon intervals) and folds it
/// through `ExonGraph::from_copies` — the exon-level sibling of `build_copy_graph`'s PSV-bubble graph.
/// Per-copy status reuses `copy_status` (DRY with `build_copy_graph`), threading the same `--gff`
/// annotation overlay (`ann`) so in-genome copies come back `InGenomeAnnotated`/`InGenomeUnannotated`
/// under `--gff`, exactly as v1's `.phase.gfa`. The per-exon reference sequence is fetched later, by the
/// caller, via `ExonGraph::to_gfa`'s own `exon_seq` closure at write time (this builder lays out intervals
/// only, never sequence — no genome I/O here).
/// `eff` — see `build_copy_graph`.
fn build_exon_graph(
    fid: &str,
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    ann: Option<&[(String, u64, u64)]>,
    eff: &[AssignStatus],
) -> rustle::vg_family::copy_graph::ExonGraph {
    use rustle::vg_family::copy_graph::*;
    let n = fa.copy_tids.len();
    let copies: Vec<(String, CopyStatus, Corrob, String, Vec<(u64, u64)>)> = (0..n)
        .map(|ci| {
            // Guard: `.get(ci)` never panics on a length-0/short fixture — missing span/introns => an
            // empty-exon copy (still walks zero classes; from_copies skips zero-length intervals anyway).
            let (chrom, start, end) =
                fa.copy_spans.get(ci).cloned().unwrap_or_else(|| (fa.chrom.clone(), 0, 0));
            let introns = fa.copy_introns.get(ci).cloned().unwrap_or_default();
            // genomic exons from the intron chain + outer span bounds:
            let mut exons = Vec::with_capacity(introns.len() + 1);
            let mut prev = start;
            for (d, a) in &introns {
                exons.push((prev, *d));
                prev = *a;
            }
            exons.push((prev, end));
            let reads = fa
                .assignments
                .iter()
                .zip(eff.iter())
                .filter(|((_, a), s)| a.best_copy == ci && matches!(s, AssignStatus::Assigned))
                .count() as u32;
            let status = copy_status(fa, ci, ann);
            let corrob = Corrob {
                reads: Some(reads),
                suns: None,
                map_identity: fa.copy_map_identity.get(ci).copied().flatten(),
            };
            (format!("{}_copy{}", fid, ci), status, corrob, chrom, exons)
        })
        .collect();
    ExonGraph::from_copies(fid, &copies)
}

/// In-genome annotation axis: overlap of copy `ci`'s span with any annotated interval.
/// `None` intervals => AnnotationUnknown (we never claim "unannotated" unchecked).
fn annotation_status(
    fa: &rustle::vg_family::denovo_pipeline::FamilyAssignment,
    ci: usize,
    ann: Option<&[(String, u64, u64)]>,
) -> rustle::vg_family::copy_graph::CopyStatus {
    use rustle::vg_family::copy_graph::CopyStatus;
    let Some(ann) = ann else { return CopyStatus::AnnotationUnknown };
    let Some((c, s, e)) = fa.copy_spans.get(ci) else { return CopyStatus::AnnotationUnknown };
    let hit = ann.iter().any(|(ac, as_, ae)| ac == c && *as_ < *e && *s < *ae);
    if hit { CopyStatus::InGenomeAnnotated } else { CopyStatus::InGenomeUnannotated }
}

/// Parse a gene annotation file (`--gff`) into `(chrom, start0, end)` intervals: accepts BED (0-based, cols
/// 0/1/2) and GFF3/GTF (1-based, cols 0/3/4 -> start converted to 0-based). Lines starting with `#` and blank
/// lines are skipped; malformed lines (wrong column count or non-numeric coords) are skipped, never panic.
///
/// Format is decided from the FILE EXTENSION first — a BED6/BED12 line whose NAME column (col 3) is numeric
/// (`chr1\t1000\t2000\t5\t0\t+`) would otherwise be silently misread as GFF, yielding the wrong interval.
/// `.bed` => BED; `.gff`/`.gff2`/`.gff3`/`.gtf` => GFF/GTF. An unknown extension falls back to the per-line
/// numeric heuristic (GFF cols 3/4 preferred, else BED cols 1/2) for best-effort on mislabeled files.
fn parse_annotation(path: &str) -> anyhow::Result<Vec<(String, u64, u64)>> {
    enum Fmt { Bed, Gff, Auto }
    let lower = path.to_ascii_lowercase();
    let fmt = if lower.ends_with(".bed") {
        Fmt::Bed
    } else if lower.ends_with(".gff") || lower.ends_with(".gff2") || lower.ends_with(".gff3") || lower.ends_with(".gtf") {
        Fmt::Gff
    } else {
        Fmt::Auto
    };
    let mut out = Vec::new();
    for line in std::fs::read_to_string(path)?.lines() {
        if line.starts_with('#') || line.trim().is_empty() { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        let bed = |f: &[&str], out: &mut Vec<(String, u64, u64)>| {
            if f.len() < 3 { return; }
            let (Ok(s), Ok(e)) = (f[1].parse::<u64>(), f[2].parse::<u64>()) else { return }; // BED 0-based
            out.push((f[0].to_string(), s, e));
        };
        let gff = |f: &[&str], out: &mut Vec<(String, u64, u64)>| {
            if f.len() < 5 { return; }
            let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else { return }; // GFF/GTF 1-based
            out.push((f[0].to_string(), s.saturating_sub(1), e));
        };
        match fmt {
            Fmt::Bed => bed(&f, &mut out),
            Fmt::Gff => gff(&f, &mut out),
            // unknown extension: prefer the GFF shape (5+ cols, numeric 3/4), else BED (3+ cols, numeric 1/2).
            Fmt::Auto => {
                let before = out.len();
                gff(&f, &mut out);
                if out.len() == before { bed(&f, &mut out); }
            }
        }
    }
    Ok(out)
}

fn parse_region(s: &str) -> Result<(String, u64, u64)> {
    let tok = s.split_whitespace().next().context("empty region")?;
    let (chrom, range) = tok.split_once(':').context("region must be chrom:start-end")?;
    let (lo_s, hi_s) = range.split_once('-').context("region must be chrom:start-end")?;
    Ok((chrom.to_string(), lo_s.parse().context("bad region start")?, hi_s.parse().context("bad region end")?))
}

/// A swept region, keyed exactly as the sweep iterates it.
type RegionKey = (String, u64, u64);
/// `--families`: the supplied catalog families BOUND to the region that will assign them.
type RegionFamilies = std::collections::BTreeMap<RegionKey, Vec<CatalogFamily>>;
/// Flank loaded around each supplied copy when gathering a dispersed family's reads (§6dh). Comfortably
/// above any long read, so a read reaching into a copy from outside its span is still collected.
const COPY_READ_PAD: u64 = 50_000;
type RegionWindows = std::collections::BTreeMap<RegionKey, Vec<(u64, u64)>>;
/// `--families`: catalog `tid` -> `(catalog family_id, catalog copy_idx)`. The JOIN KEY. Built from the
/// supplied table (never from the assignment output), so `<out>.family_join.tsv` reports the catalog's own
/// identity for a copy rather than an index this binary re-derived.
type CatalogIndex = std::collections::HashMap<String, (String, usize)>;

fn build_catalog_index(rf: &RegionFamilies) -> CatalogIndex {
    let mut ix = CatalogIndex::new();
    for fams in rf.values() {
        for f in fams {
            for c in &f.copies {
                ix.insert(c.tid.clone(), (c.family_id.clone(), c.copy_idx));
            }
        }
    }
    ix
}

/// Reference end (0-based, exclusive) of an aligned read. Local mirror of
/// `copy_assign_pipeline::read_ref_end` (which is `pub(crate)`), used only for the "every supplied copy has
/// reads" contract check below.
/// Aligned blocks of a read (0-based half-open): `M`/`=`/`X`/`D` extend, `N` closes.
fn aligned_blocks_local(read: &rustle::vg_family::copy_split::AlignedRead) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    let mut p = read.ref_start;
    let mut cur: Option<(u64, u64)> = None;
    for &(op, n) in &read.cigar {
        match op {
            'M' | '=' | 'X' | 'D' => {
                cur = Some((cur.map_or(p, |c| c.0), p + n));
                p += n;
            }
            'N' => {
                if let Some(c) = cur.take() {
                    out.push(c);
                }
                p += n;
            }
            _ => {}
        }
    }
    if let Some(c) = cur {
        out.push(c);
    }
    out
}

fn read_ref_end_local(read: &rustle::vg_family::copy_split::AlignedRead) -> u64 {
    read.ref_start
        + read.cigar.iter().filter(|(op, _)| matches!(op, 'M' | '=' | 'X' | 'D' | 'N')).map(|(_, n)| n).sum::<u64>()
}

/// Load, VALIDATE and region-bind the `--families` catalog (see the flag's help for the contract).
///
/// Returns `(None, None)` when `--families` was not given — the historical path, untouched.
///
/// Everything here is a hard error. The one thing this function must never do is drop a supplied copy:
/// a copy silently missing from O2's roster is indistinguishable, downstream, from a copy O2 legitimately
/// found no evidence for, and that ambiguity is the exact defect class this contract exists to remove.
fn load_supplied_families(
    args: &Args,
    by_contig: &std::collections::BTreeMap<String, Vec<(u64, u64)>>,
) -> Result<(Option<RegionFamilies>, Option<SeqIndex>, Option<RegionWindows>)> {
    let Some(path) = args.families.as_deref() else {
        if args.copies_fa.is_some() {
            anyhow::bail!("--copies-fa is only meaningful with --families (it supplies the copies' sequences)");
        }
        return Ok((None, None, None));
    };
    // Roster-CHANGING legs are refused, not silently applied: with --families the copy set must be exactly
    // the supplied one, and each of these adds or removes copies. --vg-realign-correct is deliberately NOT
    // here: it re-threads reads among the GIVEN copies and never touches the roster.
    for (on, flag, why) in [
        (args.absent_copies, "--absent-copies", "admits reference-absent copies"),
        (args.vg_realign, "--vg-realign", "admits novel read pools as new copies"),
        (args.iterative_prune, "--iterative-prune", "merges/drops copies"),
        (args.collapse_gate, "--collapse-gate", "admits collapsed loci as extra copies"),
        (args.tied_seed, "--tied-seed", "seeds additional loci as copies"),
        (args.recover_copies, "--recover-copies", "feeds tied secondaries into copy rescue"),
    ] {
        if on {
            anyhow::bail!(
                "--families is incompatible with {flag}: it {why}, so the assigned copy set would no longer \
                 be the supplied catalog. Drop {flag}, or run without --families."
            );
        }
    }
    if args.no_refine {
        eprintln!(
            "[copy_assign] NOTE: --no-refine is moot under --families (no family construction runs at all)."
        );
    }
    if args.homology_primary {
        eprintln!(
            "[copy_assign] NOTE: --homology-primary is moot under --families (the membership oracle is the \
             supplied catalog, not E_c or E_r)."
        );
    }
    let text = std::fs::read_to_string(path).with_context(|| format!("reading --families {path}"))?;
    let fams = group_families(parse_copies_tsv(&text).with_context(|| format!("parsing --families {path}"))?)?;
    let seqs = match args.copies_fa.as_deref() {
        Some(p) => {
            let t = std::fs::read_to_string(p).with_context(|| format!("reading --copies-fa {p}"))?;
            Some(parse_copies_fa(&t).with_context(|| format!("parsing --copies-fa {p}"))?)
        }
        None => None,
    };
    // Bind each family to the ONE swept region that contains it. Containment (not overlap) is required:
    // a family straddling a region boundary would be assigned against the reads of only part of its own
    // span, which is the truncation this mode exists to prevent.
    let mut bound: RegionFamilies = RegionFamilies::new();
    for f in fams {
        let hits: Vec<RegionKey> = by_contig
            .get(&f.chrom)
            .map(|rs| {
                rs.iter()
                    .filter(|&&(lo, hi)| f.start >= lo && f.end <= hi)
                    .map(|&(lo, hi)| (f.chrom.clone(), lo, hi))
                    .collect()
            })
            .unwrap_or_default();
        match hits.len() {
            0 => anyhow::bail!(
                "--families: {} ({}:{}-{}) lies outside every --region/--regions entry, so its reads would \
                 never be read. Add a region containing it, or remove it from the catalog.",
                f.family_id, f.chrom, f.start, f.end
            ),
            1 => bound.entry(hits[0].clone()).or_default().push(f),
            n => anyhow::bail!(
                "--families: {} ({}:{}-{}) is contained in {n} different swept regions, so which reads it \
                 would be assigned against is ambiguous. De-duplicate the region list.",
                f.family_id, f.chrom, f.start, f.end
            ),
        }
    }
    // ⭐ DISPERSED-FAMILY READ WINDOWS (§6dh). A family binds to the ONE region containing its whole
    // span, but a genuinely DISPERSED family (NPIP: 38 copies over 89.5 Mb) makes that region enormous
    // and loading it whole costs 254,726 primaries to assign copies occupying a few hundred kb — it OOMs.
    // A read that overlaps NO copy can never be assigned to one, so the region's reads are gathered from
    // the union of the copies' own neighbourhoods instead. The anti-truncation guarantee is preserved:
    // every copy's reads are still loaded in full.
    let mut windows: std::collections::BTreeMap<RegionKey, Vec<(u64, u64)>> =
        std::collections::BTreeMap::new();
    for (k, fs) in &bound {
        let mut w: Vec<(u64, u64)> = fs
            .iter()
            .flat_map(|f| f.copies.iter())
            .map(|c| (c.start.saturating_sub(COPY_READ_PAD), c.end + COPY_READ_PAD))
            .collect();
        w.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::with_capacity(w.len());
        for (lo, hi) in w {
            match merged.last_mut() {
                Some(last) if lo <= last.1 => last.1 = last.1.max(hi),
                _ => merged.push((lo, hi)),
            }
        }
        windows.insert(k.clone(), merged);
    }
    let n_fam: usize = bound.values().map(|v| v.len()).sum();
    let n_copy: usize = bound.values().flatten().map(|f| f.copies.len()).sum();
    eprintln!(
        "[copy_assign] --families {path}: {n_fam} famil{} / {n_copy} copies bound to {} region(s); sequences \
         from {}",
        if n_fam == 1 { "y" } else { "ies" },
        bound.len(),
        if seqs.is_some() { "--copies-fa (the catalog's own bytes)" } else { "--fasta (rebuilt at the catalog's exon coordinates)" },
    );
    Ok((Some(bound), seqs, Some(windows)))
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
    junction_conflict: bool,
    origin_rejected: bool,
    n_candidates: usize,
    /// The read has an aligned BASE inside a copy of its family (§6es hygiene): rows with `false` are reads
    /// gathered from the copies' neighbourhoods that overlap no copy; report O2 on `in_copy == true`.
    in_copy: bool,
    /// §6fq: the molecule's primary MAPQ < 60 — the aligner could not place it; the certificate machinery
    /// decided this row. `false` = assigned to its placement (the certificate only reported).
    contested: bool,
    /// The catalog `copy_idx` of `assigned_copy` under `--families` (copy_assign SORTS copies and reports its
    /// own index; `family_join.tsv` carries the same map). `NA` without a catalog.
    catalog_copy_idx: String,
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
    /// L4: Σ posterior over molecules the certificate did not reject (a sole candidate counts 1, a K = 0 tie
    /// 1/k, an uncertified pair its softmax). Under the read-star this is the abundance's numerator.
    n_soft: f64,
    /// Tie-break invariance certificate: reads assigned to this copy that map UNIQUELY (`mapq > 0`), so their
    /// support survives any primary/secondary relabeling of tied reads. See `anchored_support`.
    anchored: usize,
    /// `anchored >= GATE_MIN_READS`: the copy exists under every tie-break (adversarially invariant) via unique
    /// mappers. FALSE = not guaranteed by unique mappers alone (may still be junction-defined, e.g. DAZ2).
    tie_invariant: bool,
    /// `copy_junction_support >= GATE_MIN_READS`: the copy is pinned by >= 3 reads carrying a copy-specific
    /// JUNCTION (identifies it by splice structure regardless of the primary label) — the DAZ2-rescue mechanism.
    /// A copy is invariant overall if `tie_invariant || junction_invariant`.
    junction_invariant: bool,
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

/// Read the `lambda_global` scalar from a `lambda_global.tsv` (header + one data row, first column). `None` if
/// missing, unreadable, or the value is `NA`.
fn read_lambda_file(path: &str) -> Option<f64> {
    let text = std::fs::read_to_string(path).ok()?;
    let data = text.lines().nth(1)?; // skip header
    data.split('\t').next()?.trim().parse::<f64>().ok()
}

/// Resolve lambda: explicit scalar > file > none.
/// Anchored (tie-break-invariant) support for a copy: assigned reads that map UNIQUELY (`mapq > 0`) and whose
/// `best_copy == ci`. minimap2 sets `mapq = 0` exactly when the primary is tied/arbitrary, so a `mapq > 0` read
/// is never a candidate for primary/secondary relabeling — its support for `ci` is fixed under every tie-break.
/// A copy with `anchored >= GATE_MIN_READS` therefore exists under EVERY relabeling (adversarially invariant).
/// `best_copies` and `mapqs` are parallel over the family's assignments. Reproduces, in-binary, the experiment's
/// `samtools view -c -F 2308 -q 1 <copy_span>` bound.
fn anchored_support(best_copies: &[usize], mapqs: &[u8], ci: usize) -> usize {
    best_copies.iter().zip(mapqs.iter()).filter(|(&bc, &mq)| bc == ci && mq > 0).count()
}

fn resolve_lambda(explicit: Option<f64>, from_file: Option<f64>) -> Option<f64> {
    explicit.or(from_file)
}

fn verdict_str(v: rustle::vg_family::linearize::Verdict) -> &'static str {
    use rustle::vg_family::linearize::Verdict::*;
    match v {
        Linearizes => "LINEARIZES",
        Not => "NOT",
        Undetermined => "UNDETERMINED",
    }
}

/// One `<out>.linearize.tsv` row for a single Stage-2-admitted reference-absent candidate's
/// augment-and-linearize certificate (Task 4's `linearize_certs`). `NA` for NaN fracs/perm_p
/// (the `n_pool < min_pool` short-circuit in `linearize_certificate`).
fn linearize_tsv_row(fam: &str, loc: (&str, u64, u64), c: &LinearizeCertificate) -> String {
    let f = |x: f64| if x.is_nan() { "NA".to_string() } else { format!("{:.3}", x) };
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        fam,
        loc.0,
        loc.1,
        loc.2,
        c.n_pool,
        f(c.linearized_frac_real),
        f(c.mean_frac_decoy),
        f(c.delta),
        if c.perm_p.is_nan() { "NA".to_string() } else { format!("{:.4}", c.perm_p) },
        verdict_str(c.verdict)
    )
}

fn main() -> Result<()> {
    let mut args = Args::parse();
    // §6eu: the pipeline reads RUSTLE_PSV_READFILTER; an explicit env value wins, else the flag decides.
    if std::env::var_os("RUSTLE_PSV_READFILTER").is_none() {
        std::env::set_var("RUSTLE_PSV_READFILTER", if args.psv_read_filter { "1" } else { "0" });
    }
    if args.igv {
        args.dump_psv = true; // --igv is a bundle: the PSV matrix feeds bench/igv_tracks.py -> tagged BAM + PSV VCF
    }
    // --gff: parsed ONCE before the sweep into the annotation axis intervals (None => every in-genome copy in
    // the --phase copy graph stays AnnotationUnknown, byte-identical to the no-flag path).
    let annotation: Option<Vec<(String, u64, u64)>> =
        args.gff.as_deref().map(parse_annotation).transpose().context("parsing --gff")?;

    // Augment-and-linearize is opt-in: the report (`--linearize`) or the gate (`--linearize-gate`, which
    // implies the report) turns it on. When false, the certificate is skipped entirely inside
    // `detect_and_assign` (no minimap2 realign-pool per candidate) and no `<out>.linearize.tsv` is written,
    // so plain `--absent-copies` is byte-identical to the pre-feature `--absent-copies` path.
    let do_linearize = args.linearize || args.linearize_gate;

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

    // ---- O1 -> O2 FILE CONTRACT (`--families`) -------------------------------------------------------
    // The catalog is loaded, validated and BOUND TO REGIONS here, before a single read is touched, so a
    // malformed or unassignable roster fails in the first second rather than after an hour of alignment.
    // `region_families` is keyed by the exact `(contig, lo, hi)` triple the sweep iterates, so each region's
    // worker gets its own families with no re-derivation and no ambiguity about which region owns a family.
    let (region_families, catalog_seqs, region_windows) = load_supplied_families(&args, &by_contig)?;
    let catalog_index: Option<CatalogIndex> = region_families.as_ref().map(build_catalog_index);

    let lambda = resolve_lambda(args.lambda_global, args.lambda_file.as_deref().and_then(read_lambda_file));
    let mut cfg = DenovoConfig::from_env();
    cfg.detect.len_cap = args.max_poa_len; // poasta memory threshold: above it, the bounded LCS fallback
    cfg.vg_realign = args.vg_realign || args.vg_realign_correct; // VG correction leg (re-thread hard reads)
    cfg.vg_realign_admit = args.vg_realign; // roster-widening admission: only via the full --vg-realign
    cfg.homology_primary = args.homology_primary; // E_r membership; off => the E_c path is untouched
    cfg.refine = !args.no_refine; // mutual-homology family gate (matches gw_family_catalog); on by default
    cfg.filter_readthrough = !args.keep_readthrough; // unspliced pre-mRNA spans are not copies
    cfg.gate.pool_locus_support = !args.no_pool_locus_support;
    cfg.tied_seed = args.tied_seed; // seed from AS-tied secondaries
    cfg.collapse_gate = args.collapse_gate; // experimental; detects paralogy, not collapse (see module header)
    cfg.collapse_enumerate = args.collapse_enumerate || cfg.collapse_enumerate; // CLI OR env (RUSTLE_COLLAPSE_ENUMERATE)
    if cfg.collapse_enumerate {
        eprintln!(
            "[copy_assign] WARNING: --collapse-enumerate / RUSTLE_COLLAPSE_ENUMERATE is a no-op in copy_assign \
             (no consumer in this binary); run gw_family_catalog --collapse-enumerate instead."
        );
    }
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
        junction_conflict_abstain: args.junction_conflict_abstain,
        molecule_pool: args.molecule_observations && !args.no_molecule_observations,
        origin_subst_only: args.origin_substitutions_only,
        read_star_junctions: args.read_star_junctions,
        read_star_genomic: args.read_star_genomic && !args.read_star_unit,
        read_star_catalog_locus: !args.read_star_pad_locus,
        sole_candidate: !args.no_sole_candidate,
        dump_star: args.dump_star,
        read_star_hit_in_unit: !args.no_read_star_hit_in_unit,
        read_star_two_form: !args.read_star_genomic_only,
        ..AssignParams::default()
    };
    eprintln!("[copy_assign] decisive-margin tau={} error_rate={}", args.margin, args.error_rate);
    let mut family_rows: Vec<FamilyRow> = Vec::new();
    let mut placement_assigned_total = 0usize; // §6fq: uncontested molecules assigned to their placement
    let mut assign_rows: Vec<AssignRow> = Vec::new();
    let mut posterior_lines: Vec<String> = Vec::new();
    // EM-abundance prior for the posterior (else uniform).
    let prior_abundance = std::env::var("RUSTLE_POSTERIOR_PRIOR").ok().as_deref() == Some("abundance");
    // Cross-family reconciliation mode (see `XfamMode`). Parsed HERE, before any read is touched, so an
    // unrecognized value fails in the first second rather than silently running as `off`.
    let xfam_mode = XfamMode::from_env()?;
    // locus from a de-novo tid `DN_<chrom>_<start>_<n>` (chrom may contain `_`, so split from the right).
    fn parse_locus(tid: &str) -> Option<(String, u64)> {
        let rest = tid.strip_prefix("DN_")?;
        let parts: Vec<&str> = rest.rsplitn(3, '_').collect(); // [n, start, chrom]
        Some((parts.get(2)?.to_string(), parts.get(1)?.parse().ok()?))
    }
    let mut quant_rows: Vec<QuantRow> = Vec::new();
    // `--families`: one row per ASSIGNED copy, naming the catalog row it came from. The explicit join
    // between `<out>.quant.tsv` and the O1 `copies.tsv`, and the place a copy that failed to survive
    // assignment would be visible as a missing row.
    let mut join_rows: Vec<String> = Vec::new();
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
    let mut gfa_links: HashSet<String> = HashSet::new();        // dedup'd full "L\t..." strings (copy_graph emits complete lines)
    let mut gfa_paths: Vec<String> = Vec::new();
    // VG read-threading (the Canzar flip, materialized): each read WALKS the PSV-bubble nodes for the alleles
    // it observes, REUSING a copy's node wherever their alleles agree — so multimapping reads become shared
    // threaded evidence through the one family graph. W-lines (GFA 1.1) + a Bandage node/path colour CSV.
    let mut gfa_walks: Vec<String> = Vec::new();
    let mut gfa_colors: Vec<String> = Vec::new(); // "name,colour" for Bandage (copies distinct, reads by assigned copy)
    let mut legend_rows: Vec<String> = Vec::new(); // "status\tcolour" (de-duplicated at write time)
    // --phase v2: one exon presence/absence graph per family (built during the drain, where `fa` is in
    // scope; sequence-free — `to_gfa` fetches reference bases lazily at write time via `genome_for`).
    let mut exon_graphs: Vec<rustle::vg_family::copy_graph::ExonGraph> = Vec::new();
    let mut fallback_all: Vec<FallbackEdge> = Vec::new(); // family edges confirmed via the LCS fallback
    let mut dna_needs_rows: Vec<DnaNeedsRecord> = Vec::new(); // --absent-copies: candidates needing DNA validation
    // --absent-copies + opt-in --linearize/--linearize-gate: linearize certificates, one per Stage-2-admitted
    // candidate (Task 4), written to `<out>.linearize.tsv` below (Task 5) when `do_linearize`. Empty otherwise
    // (the cert is skipped in `detect_and_assign`). `--linearize-gate` also uses the verdict to gate admission
    // itself, so a demoted candidate shows up here but not in `fams`.
    let mut linearize_certs_all: Vec<(String, LinearizeCertificate, (String, u64, u64))> = Vec::new();
    let mut vg_realign_lines: Vec<String> = Vec::new(); // --vg-realign: per-read re-align decisions (report-only)
    let mut gfam = 0usize; // global family counter (unique ids across regions)
    let mut gtf_lines: Vec<String> = Vec::new(); // --gtf: FLAIR-style isoform GTF (transcript + exon rows)

    // `--skip-poa-diagnostic` is read by `detect_and_assign` via this env var (it is purely diagnostic and
    // does not change the emitted families/assignments — see the flag's help).
    if args.skip_poa_diagnostic {
        std::env::set_var("RUSTLE_SKIP_POA_DIAGNOSTIC", "1");
    }
    // `--poa-cap` is read by `discover_intron_psvs` (copy_assign_pipeline.rs) via this env var, the same
    // "flag -> env var -> deep read" idiom as `--skip-poa-diagnostic` above: the const it replaces lives many
    // call frames below `main` (through `assign_family`/`assign_family_detailed`/`find_weak_copies`), so a
    // signature thread-through would touch dozens of call sites for an opt-in (`RUSTLE_INTRON_PSV=1`) code
    // path. Always set (not gated on non-default) so the resolved value is unambiguous; the default 20000
    // reproduces the prior hard-coded constant exactly.
    std::env::set_var("RUSTLE_POA_CAP", args.poa_cap.to_string());
    if args.psv_genomic {
        std::env::set_var("RUSTLE_PSV_GENOMIC", "1");
    }
    // `--absent-min-clusters` reaches `absent_copy::AbsentCopyParams::from_env` through the same
    // "flag -> env var -> deep read" idiom (the params are built inside `detect_and_assign`'s
    // candidate loop, many frames below `main`). Set ONLY when the flag is given, so (a) an
    // unset run is byte-identical and (b) a caller who exported RUSTLE_ABSENT_MIN_CLUSTERS
    // directly is not silently clobbered by the flag's default.
    if let Some(n) = args.absent_min_clusters {
        std::env::set_var(rustle::vg_family::absent_copy::MIN_CLUSTERS_ENV, n.to_string());
        if !args.absent_copies {
            eprintln!(
                "[copy_assign] WARNING: --absent-min-clusters={n} has no effect without --absent-copies"
            );
        }
    }
    // `--read-cap` is a NO-OP in copy_assign (see the flag's help): `o2_materialize::READ_CAP` has no
    // consumer in any `src/bin/*.rs` binary. Warn rather than silently ignore a non-default value.
    if args.read_cap != 6_000 {
        eprintln!(
            "[copy_assign] WARNING: --read-cap={} has no consumer in this binary (o2_materialize's READ_CAP \
             is not on copy_assign's execution path — it backs a Rust byte-parity port of the Python \
             genome-wide-catalog materializer that no shipped binary imports); the value is ignored.",
            args.read_cap
        );
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
        // §6dh: on the --families path, gather from the supplied copies' own neighbourhoods rather than
        // the whole bound region — a dispersed family's hull can be tens of Mb while its copies occupy a
        // few hundred kb, and a read overlapping no copy can never be assigned to one.
        let wins: Vec<(u64, u64)> = region_windows
            .as_ref()
            .and_then(|w| w.get(&(contig.clone(), lo, hi)))
            .map(|v| v.iter().map(|&(a, b)| (a.max(lo), b.min(hi))).filter(|&(a, b)| a < b).collect())
            .unwrap_or_else(|| vec![(lo, hi)]);
        let (primary, bam_reads) = {
            let mut pr: Vec<_> = Vec::new();
            let mut br: Vec<_> = Vec::new();
            let mut seen = std::collections::HashSet::new();
            for &(wlo, whi) in &wins {
                let (p, b) = match &bam_cache {
                    Some(c) => c.reads_in_region(&args.bam, contig, wlo, whi),
                    None => reads_in_region(&args.bam, contig, wlo, whi, args.threads),
                }
                .with_context(|| format!("reading {contig}:{wlo}-{whi}"))?;
                // Windows are disjoint after merging, but a read spanning a boundary is returned by
                // both queries; key on (name, start) so one molecule is never two witnesses.
                for x in p {
                    // PrimaryRead has no name; its (chrom, span, intron chain) identifies the placement.
                    if seen.insert((x.chrom.clone(), x.ref_start, x.ref_end, x.introns.clone())) {
                        pr.push(x);
                    }
                }
                for x in b {
                    br.push(x);
                }
            }
            let mut bseen = std::collections::HashSet::new();
            br.retain(|x: &rustle::vg_family::denovo_assemble::BamRead| {
                bseen.insert((x.name.clone(), x.read.ref_start))
            });
            (pr, br)
        };
        if timing && wins.len() > 1 {
            eprintln!(
                "[timing] {contig}:{lo}-{hi} gathered from {} copy window(s) ({:.1} Mb of {:.1} Mb hull)",
                wins.len(),
                wins.iter().map(|&(a, b)| (b - a) as f64).sum::<f64>() / 1e6,
                (hi - lo) as f64 / 1e6
            );
        }
        if timing {
            eprintln!(
                "[timing] reads_in_region {contig}:{lo}-{hi} ({} reads): {:.1}s",
                bam_reads.len(),
                t_read.elapsed().as_secs_f64()
            );
        }
        let extra = if args.recover_copies || args.tied_seed {
            tied_secondary_reads_in_region(&args.bam, contig, lo, hi, args.as_ratio).unwrap_or_default()
        } else {
            Vec::new()
        };
        // `--families`: materialize THIS region's supplied catalog families into the copy set, and enforce
        // the last two contract clauses (a sequence exists for every copy; every copy has reads here).
        // Done inside the worker because it needs this region's genome and this region's BAM slice.
        let supplied: Option<Vec<ColocatedFamily>> = match &region_families {
            None => None,
            Some(rf) => {
                let mine = rf.get(&(contig.clone(), lo, hi)).map(|v| v.as_slice()).unwrap_or(&[]);
                let mut v: Vec<ColocatedFamily> = Vec::with_capacity(mine.len());
                for f in mine {
                    let (cf, _src) = to_colocated(f, catalog_seqs.as_ref(), &genome)?;
                    for c in &cf.copies {
                        // A supplied copy with no read here cannot be assigned anything, and dropping it
                        // silently would understate K and loosen the Bonferroni certificate for the rest.
                        let n = bam_reads
                            .iter()
                            .filter(|br| {
                                br.chrom == c.chrom
                                    && br.read.ref_start < c.end
                                    && read_ref_end_local(&br.read) > c.start
                            })
                            .count();
                        if n == 0 {
                            anyhow::bail!(
                                "--families: {} copy {} ({}:{}-{}) has NO reads in {contig}:{lo}-{hi} of \
                                 {}. It cannot be assigned, and silently dropping it would understate the \
                                 family's copy count. Check the BAM is the one the catalog was built from \
                                 (subset BAMs are the recurring trap here).",
                                cf.family_id, c.tid, c.chrom, c.start, c.end, args.bam
                            );
                        }
                    }
                    v.push(cf);
                }
                Some(v)
            }
        };
        let t_da = std::time::Instant::now();
        let (fams, fallback, dna_needs, linearize_certs) = detect_and_assign(
            &primary, &bam_reads, &genome, &cfg, args.win, args.min_copies, &params, &extra,
            args.absent_copies, do_linearize, args.linearize_gate, &args.fasta,
            supplied.as_deref(),
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
        let read_mapqs: Vec<u8> = bam_reads.iter().map(|r| r.mapq).collect();
        let read_spans: Vec<(u64, u64, u8)> = bam_reads
            .iter()
            .map(|r| {
                (
                    r.read.ref_start,
                    read_ref_end_local(&r.read),
                    (r.is_secondary as u8) | ((r.is_supplementary as u8) << 1),
                )
            })
            .collect();
        let read_blocks: Vec<Vec<(u64, u64)>> = bam_reads.iter().map(|r| aligned_blocks_local(&r.read)).collect();
        let as_ev = as_evidence_per_read(&bam_reads);
        let n_mapped = bam_reads.len();
        Ok(RegionWork { contig: contig.clone(), lo, hi, read_names, read_mapqs, read_spans, read_blocks, as_ev, n_mapped, fams, fallback, dna_needs, linearize_certs, transcripts })
    };
    // Compute all regions (out-of-order across contigs when region_threads > 1), collected in the flat order.
    let works: Vec<RegionWork> = match &region_pool {
        Some(pool) => pool.install(|| {
            flat.par_iter().map(|&(c, lo, hi)| compute(c, lo, hi)).collect::<Result<Vec<_>>>()
        })?,
        None => flat.iter().map(|&(c, lo, hi)| compute(c, lo, hi)).collect::<Result<Vec<_>>>()?,
    };
    // PASS 1 (read-only): cross-family reconciliation. It must run BEFORE the drain, not as a post-pass
    // over `assign_rows`, because a molecule's status is emitted from `fa.assignments` at FOUR sites
    // inside the drain (`.assignments.tsv`, `.posterior.tsv`, `.psv_reads.tsv`, `.phased_reads.tsv`) plus
    // the two `--phase` graph builders' `Assigned` filters. A post-pass would fix only the first and
    // leave the rest disagreeing with it — `status_consistency_across_outputs` is the regression that
    // catches exactly that. Under `off` this does not run at all.
    let (xfam_conflicts, xfam_demote) = if xfam_mode == XfamMode::Off {
        (Vec::new(), std::collections::BTreeSet::new())
    } else {
        xfam_pass1(&works, region_families.is_some())
    };
    // The status a row is EMITTED with. Under `off`/`report` (and for every non-demoted key) this is
    // `a.status` verbatim, so those arms are byte-identical BY CONSTRUCTION. The `Abstain` test comes
    // FIRST so the other two modes never even build the lookup key.
    let eff_astatus = |read_name: &str, g: usize, f: usize, a: &rustle::vg_family::copy_assign::Assignment| -> AssignStatus {
        if xfam_mode == XfamMode::Abstain
            && matches!(a.status, AssignStatus::Assigned)
            && xfam_demote.contains(&(read_name.to_string(), g, f))
        {
            // Ambiguous, NEVER Tied: `Tied` is reserved for the K=0 identifiability wall
            // (`min_p >= alpha/(n-1)`), and the warning below exists precisely so catalog artifacts do
            // not masquerade as that wall. The intra-family record contradiction already demotes to
            // Ambiguous; the cross-family scope fix inherits the same bucket.
            AssignStatus::Ambiguous
        } else {
            a.status
        }
    };
    let eff_status = |read_name: &str, g: usize, f: usize, a: &rustle::vg_family::copy_assign::Assignment| -> &'static str {
        status_str(eff_astatus(read_name, g, f, a))
    };
    // SERIAL drain (PASS 2) in the original region order — every row push + the `gfam` id counter is
    // exactly the serial path, so the output is byte-identical.
    {
        for (gwork, work) in works.into_iter().enumerate() {
            let RegionWork { contig, lo, hi, read_names, read_mapqs, read_spans: _, read_blocks, as_ev, n_mapped, fams, fallback, dna_needs, linearize_certs, transcripts } = work;
            let contig = &contig;
            let bam_reads = &read_names; // output stage indexes read NAMES (sequences were dropped)
            fallback_all.extend(fallback);
            dna_needs_rows.extend(dna_needs);
            linearize_certs_all.extend(linearize_certs);
            // Best mapq over each MOLECULE's records in this region (only its primary record can be >0).
            // See the `mqs` construction below: the tie-break invariance certificate is a molecule property.
            let mol_mapq: std::collections::HashMap<&str, u8> = {
                let mut m: std::collections::HashMap<&str, u8> = std::collections::HashMap::new();
                for (n, &q) in bam_reads.iter().zip(read_mapqs.iter()) {
                    let e = m.entry(n.as_str()).or_insert(0);
                    *e = (*e).max(q);
                }
                m
            };
            // ⭐ §6fq (user, 2026-09-06): the certificate machinery is for the CONTESTED molecules — the ones the
            // aligner could not place (primary MAPQ < 60). An uncontested molecule is assigned to its placement
            // (the copy its primary's blocks overlap most), as any assembler would use it; the certificate is
            // still computed for it and reported (`origin_rejected`), never applied. One sensitivity over every
            // read; abstention only among the contested. `--no-placement-assign` = the machinery on every read.
            let placement_assign = args.molecule_observations && !args.no_molecule_observations && !args.no_placement_assign;
            let mut placement_assigned = 0usize;
            let mut fams = fams;
            // the molecule's PRIMARY record (its highest-MAPQ record): the row's `ri` is the read-star
            // representative, which can be a secondary record at another copy
            let mol_primary: std::collections::HashMap<&str, usize> = {
                let mut m: std::collections::HashMap<&str, usize> = std::collections::HashMap::new();
                for (i, n) in bam_reads.iter().enumerate() {
                    let e = m.entry(n.as_str()).or_insert(i);
                    if read_mapqs[i] > read_mapqs[*e] {
                        *e = i;
                    }
                }
                m
            };
            if placement_assign {
                for fa in fams.iter_mut() {
                    for (ri, a) in fa.assignments.iter_mut() {
                        let mq = mol_mapq.get(bam_reads[*ri].as_str()).copied().unwrap_or(read_mapqs[*ri]);
                        if mq < 60 {
                            continue;
                        }
                        // certified first: a molecule the machinery already assigned keeps that call (it can
                        // correct a placement: 4 % of MAPQ-60 simulated reads sit at the wrong copy, §6fq);
                        // the placement is the fallback when the machinery abstains or ties
                        if !args.placement_first && a.status == rustle::vg_family::copy_assign::AssignStatus::Assigned {
                            continue;
                        }
                        let pri = mol_primary.get(bam_reads[*ri].as_str()).copied().unwrap_or(*ri);
                        let Some(bl) = read_blocks.get(pri) else { continue };
                        let mut best: Option<(usize, u64)> = None;
                        for (ci, (c, s0, e0)) in fa.copy_spans.iter().enumerate() {
                            if c != contig {
                                continue;
                            }
                            let o: u64 = bl.iter().map(|&(bs, be)| be.min(*e0).saturating_sub(bs.max(*s0))).sum();
                            if o > 0 && best.map_or(true, |(_, bo)| o > bo) {
                                best = Some((ci, o));
                            }
                        }
                        let Some((pc, _)) = best else { continue };
                        if std::env::var_os("RUSTLE_STAR_DEBUG").is_some() && placement_assigned < 5 {
                            eprintln!("[placement] read {} mapq {mq} blocks {:?} -> copy {pc} span {:?} (was best_copy {} status {:?}); spans {:?}", bam_reads[*ri], &bl[..bl.len().min(3)], fa.copy_spans.get(pc), a.best_copy, a.status, &fa.copy_spans[..fa.copy_spans.len().min(3)]);
                        }
                        a.status = rustle::vg_family::copy_assign::AssignStatus::Assigned;
                        a.best_copy = pc;
                        a.resolvable = true;
                        let mut one = vec![0.0f64; fa.copy_spans.len()];
                        one[pc] = 1.0;
                        a.posterior = one;
                        placement_assigned += 1;
                    }
                }
            }
            placement_assigned_total += placement_assigned;
            let fams = fams;
            // --gtf: gene_tid (a copy's own locus) -> (family id, copy index), filled as fids are assigned below.
            let mut copy_gene: std::collections::HashMap<String, (String, usize)> = std::collections::HashMap::new();
            for (fwork, fa) in fams.iter().enumerate() {
                // JOIN KEY. Without `--families` this binary invents its own id (`CAFAM{i}`), which is
                // precisely why the O1 and O2 tables could not be joined. With `--families` the family
                // KEEPS the catalog's own `GWFAM{i}` — no id is minted, so `family_id` means the same
                // thing in both tables. `gfam` still advances so the two modes cannot alias.
                let fid = if region_families.is_some() { fa.family_id.clone() } else { format!("CAFAM{gfam}") };
                gfam += 1;
                if args.gtf {
                    for (ci, tid) in fa.copy_tids.iter().enumerate() {
                        copy_gene.insert(tid.clone(), (fid.clone(), ci));
                    }
                }
                let cat_idx_of = |ci: usize| -> String {
                    match (&catalog_index, fa.copy_tids.get(ci)) {
                        (Some(ix), Some(tid)) => ix.get(tid).map(|(_, i)| i.to_string()).unwrap_or_else(|| "NA".into()),
                        _ => "NA".into(),
                    }
                };
                for (ri, a) in &fa.assignments {
                    let in_copy = read_blocks.get(*ri).map_or(false, |bl| {
                        bl.iter().any(|&(bs, be)| fa.copy_spans.iter().any(|(c, s0, e0)| c == contig && be > *s0 && bs < *e0))
                    });
                    assign_rows.push(AssignRow {
                        read_name: bam_reads[*ri].clone(),
                        family_id: fid.clone(),
                        assigned_copy: a.best_copy,
                        status: eff_status(&bam_reads[*ri], gwork, fwork, a),
                        n_decisive: a.n_decisive,
                        margin: a.log_lr_margin,
                        p_value: a.p_value,
                        min_p_value: a.min_p_value,
                        as_ev: as_ev[*ri],
                        junction_conflict: a.junction_conflict,
                        origin_rejected: a.origin_rejected,
                        n_candidates: a.n_candidates,
                        in_copy,
                        contested: mol_mapq.get(bam_reads[*ri].as_str()).copied().unwrap_or(read_mapqs[*ri]) < 60,
                        catalog_copy_idx: cat_idx_of(a.best_copy),
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
                            bam_reads[*ri], fid, eff_status(&bam_reads[*ri], gwork, fwork, a), idx.len(), chrom, zs, ze, pstr
                        ));
                    }
                }
                // soft per-copy abundance (+ the hard read count for comparison) + the tie-break invariance
                // certificate (anchored = unique-mapper support; see anchored_support).
                let bcs: Vec<usize> = fa.assignments.iter().map(|(_, a)| a.best_copy).collect();
                // MOLECULE-level mapq: only a molecule's PRIMARY record carries mapq>0, and after the
                // record->molecule reduction the surviving record need not be the primary. Max over the
                // molecule's records keeps `anchored`/`tie_invariant` counting the same molecules as before.
                let mqs: Vec<u8> = fa
                    .assignments
                    .iter()
                    .map(|(ri, _)| mol_mapq.get(bam_reads[*ri].as_str()).copied().unwrap_or(read_mapqs[*ri]))
                    .collect();
                for (ci, tid) in fa.copy_tids.iter().enumerate() {
                    let anchored = anchored_support(&bcs, &mqs, ci);
                    quant_rows.push(QuantRow {
                        family_id: fid.clone(),
                        copy_index: ci,
                        copy_tid: tid.clone(),
                        copy_chrom: fa.copy_spans.get(ci).map(|s| s.0.clone()).unwrap_or_default(),
                        copy_start: fa.copy_spans.get(ci).map_or(0, |s| s.1),
                        copy_end: fa.copy_spans.get(ci).map_or(0, |s| s.2),
                        abundance: fa.copy_abundance.get(ci).copied().unwrap_or(0.0),
                        ci: fa.copy_abundance_ci.get(ci).copied().unwrap_or(0.0),
                        // an ORPHAN (no candidate, §6fg) is nobody's hard read (it carried copy 0's index by default)
                        n_hard: fa.assignments.iter().filter(|(_, a)| a.best_copy == ci && !(a.origin_rejected && a.n_candidates == 0)).count(),
                        n_soft: fa.assignments.iter().filter(|(_, a)| !a.origin_rejected).map(|(_, a)| a.posterior.get(ci).copied().unwrap_or(0.0)).sum(),
                        anchored,
                        tie_invariant: anchored as u32 >= GATE_MIN_READS,
                        junction_invariant: fa.copy_junction_support.get(ci).copied().unwrap_or(0) as u32
                            >= GATE_MIN_READS,
                    });
                    // `--families`: name the catalog row this copy IS. Looked up by the catalog `tid`
                    // (carried through `DenovoTranscript::tid` untouched), never by position, so the join
                    // does not depend on the copy ordering surviving the assignment stage.
                    if let Some(ix) = &catalog_index {
                        let (cf, cidx) = match ix.get(tid) {
                            Some((f, i)) => (f.clone(), i.to_string()),
                            // Unreachable while the roster-changing legs are refused; reported rather than
                            // hidden, because a copy O2 invented is exactly what --families must never do.
                            None => ("NOT_IN_CATALOG".to_string(), "NA".to_string()),
                        };
                        join_rows.push(format!(
                            "{fid}\t{ci}\t{tid}\t{cf}\t{cidx}\t{}\t{}\t{}\t{}",
                            fa.copy_spans.get(ci).map(|s| s.0.clone()).unwrap_or_default(),
                            fa.copy_spans.get(ci).map_or(0, |s| s.1),
                            fa.copy_spans.get(ci).map_or(0, |s| s.2),
                            fa.assignments.iter().filter(|(_, a)| a.best_copy == ci).count(),
                        ));
                    }
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
                            bam_reads[*ri], fid, a.best_copy, eff_status(&bam_reads[*ri], gwork, fwork, a),
                            a.log_lr_margin, a.n_decisive, allele_str(obs)
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
                    // Effective per-row status, built ONCE so the phase-block counters, the per-copy
                    // support counts, the read haplotypes and the two graph builders' `Assigned` filters
                    // all agree row-for-row with `.assignments.tsv`. Element-wise equal to `a.status`
                    // under `off`/`report`.
                    let eff: Vec<AssignStatus> = fa
                        .assignments
                        .iter()
                        .map(|(ri, a)| eff_astatus(&bam_reads[*ri], gwork, fwork, a))
                        .collect();
                    let n_phased = eff.iter().filter(|s| matches!(s, AssignStatus::Assigned)).count();
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
                            .zip(eff.iter())
                            .filter(|((_, a), s)| a.best_copy == ci && matches!(s, AssignStatus::Assigned))
                            .count();
                        phased_hap_lines.push(format!("{}\t{}\t{}\t{}\t{}", fid, ci, tid, n_sup, vs));
                    }
                    for ((ri, a), es) in fa.assignments.iter().zip(eff.iter()) {
                        let hap: i64 = if matches!(es, AssignStatus::Assigned) {
                            a.best_copy as i64
                        } else {
                            -1
                        };
                        phased_read_lines.push(format!(
                            "{}\t{}\t{}\t{}\t{:.3}\t{}",
                            bam_reads[*ri], fid, hap, a.n_decisive, a.log_lr_margin, status_str(*es)
                        ));
                    }

                    // Materialize this family's copy-graph (REFERENCE walk + tagged copy paths + read
                    // walks over the shared PSV-bubble nodes) and fold its GFA lines / Bandage colours /
                    // status legend into the region's accumulators. Replaces the inline emitter above.
                    let ref_base = |chrom: &str, pos: u64| {
                        genome_for(chrom).ok()
                            .and_then(|g| g.fetch_sequence(chrom, pos, pos + 1))
                            .and_then(|v| v.first().copied())
                    };
                    let cg = build_copy_graph(&fid, fa, ref_base, bam_reads, annotation.as_deref(), &eff);
                    let gl = cg.gfa_lines();
                    for s in gl.segs { gfa_segs.insert(s); }
                    for l in gl.links { gfa_links.insert(l); }
                    gfa_paths.extend(gl.paths);
                    gfa_walks.extend(gl.walks);
                    for row in cg.colours_csv().lines() { gfa_colors.push(row.to_string()); }
                    legend_rows.extend(cg.legend_tsv().lines().map(|s| s.to_string()));

                    // v2: this family's exon presence/absence graph (copies = walks over shared exon
                    // classes; a copy-specific exon reads as a visible arm). Only built when the family
                    // carries an intron chain (no `copy_introns` => no exon structure to graph).
                    // Sequence-free at this point; folded into `<out>.exon.gfa` at write time below, where
                    // `genome_for` fetches each exon's bases. Same `--gff` annotation overlay as v1.
                    if !fa.copy_introns.is_empty() {
                        exon_graphs.push(build_exon_graph(&fid, fa, annotation.as_deref(), &eff));
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
                    depth_cn: lambda.map(|lam| depth_cn(fa.n_reads, lam)).unwrap_or(f64::NAN),
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
    writeln!(ah, "read_name\tfamily_id\tassigned_copy\tstatus\tn_decisive\tmargin\tp_value\tmin_p_value\tas_best\tas_second\tas_margin\tas_per_base_best\tas_per_base_2nd\tin_copy\tcatalog_copy_idx\torigin_rejected\tn_candidates\tsole_candidate\tcontested")?;
    for r in &assign_rows {
        // L3: a CONTESTED molecule assigned with exactly one candidate is a sole candidate (§6fi); an uncontested
        // one is assigned to its placement (§6fq) whatever its candidate count
        let sole = (r.status == "assigned" && r.n_candidates == 1 && r.contested) as u8;
        writeln!(
            ah,
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3e}\t{:.3e}\t{}\t{}\t{}\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.read_name, r.family_id, r.assigned_copy, r.status, r.n_decisive, r.margin, r.p_value, r.min_p_value,
            r.as_ev.best, opt_i32(r.as_ev.second), opt_i32(r.as_ev.margin()),
            r.as_ev.best_per_base, opt_f32(r.as_ev.second_per_base), r.in_copy, r.catalog_copy_idx, r.origin_rejected as u8, r.n_candidates, sole, r.contested as u8
        )?;
    }
    {
        // §6es hygiene: the numbers to quote are on reads with an aligned base inside a copy.
        let inc: Vec<&AssignRow> = assign_rows.iter().filter(|r| r.in_copy).collect();
        let cnt = |st: &str| inc.iter().filter(|r| r.status == st).count();
        eprintln!(
            "[copy_assign] reads with an aligned base inside a copy: {} of {} rows — assigned {} / tied {} / ambiguous {}",
            inc.len(),
            assign_rows.len(),
            cnt("assigned"),
            cnt("tied"),
            cnt("ambiguous")
        );
    }

    if args.junction_conflict_abstain {
        let mut cf = std::fs::File::create(format!("{}.conflicts.tsv", args.out))?;
        writeln!(cf, "read_name\tfamily_id\tpsv_best_copy\tstatus\tn_decisive\tmin_p_value")?;
        let mut n = 0usize;
        for r in assign_rows.iter().filter(|r| r.junction_conflict) {
            writeln!(cf, "{}\t{}\t{}\t{}\t{}\t{:.3e}", r.read_name, r.family_id, r.assigned_copy, r.status, r.n_decisive, r.min_p_value)?;
            n += 1;
        }
        eprintln!("[copy_assign] junction-conflict-abstain: {n} read(s) whose splice junctions contradict their PSV-best copy -> ambiguous ({}.conflicts.tsv)", args.out);
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
        args.out, famcn_rows.len(), if lambda.is_some() { "on" } else { "NA (pass --lambda-global or --lambda-file)" });

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
    //
    // SCOPE NOTE (`RUSTLE_XFAM_RECONCILE`): this file is DELIBERATELY unmoved by a cross-family abstention.
    // `n_reads_hard` below counts `fa.assignments` by argmax `best_copy` with NO status filter, and
    // `abundance`/`ci95` come from `soft_quantify_em` inside the per-family pipeline, whose `obs_for_em` is
    // populated regardless of status — so a demoted molecule still counts here, in both arms, and
    // `.quant.tsv` is byte-identical between `report` and `abstain` (`quant_is_unmoved_by_demotion` pins
    // it). The row inflation a double-assigned molecule causes here is therefore NOT fixed by that flag;
    // removing a molecule from a family's EM is a two-pass architecture change and a separate decision.
    let mut qh = std::fs::File::create(format!("{}.quant.tsv", args.out))?;
    writeln!(qh, "family_id\tcopy_index\tcopy_tid\tcopy_chrom\tcopy_start\tcopy_end\tabundance\tci95_halfwidth\tn_reads_hard\tanchored_reads\ttie_invariant\tjunction_invariant\tn_reads_soft")?;
    for r in &quant_rows {
        writeln!(qh, "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{:.2}", r.family_id, r.copy_index, r.copy_tid,
            r.copy_chrom, r.copy_start, r.copy_end, r.abundance, r.ci, r.n_hard, r.anchored, r.tie_invariant,
            r.junction_invariant, r.n_soft)?;
    }
    // A copy is invariant to the arbitrary primary/secondary label if it is pinned by unique mappers OR by a
    // copy-specific junction (splice structure identifies it regardless of the label). Report the OR bottom line.
    let n_inv = quant_rows.iter().filter(|r| r.tie_invariant || r.junction_invariant).count();
    eprintln!(
        "[copy_assign] tie-break invariance: {}/{} copies invariant (>= {} unique-mapper OR copy-specific-junction reads; FALSE = existence leans on the arbitrary primary label)",
        n_inv, quant_rows.len(), GATE_MIN_READS
    );

    // `--families`: the explicit O1<->O2 JOIN, plus the closing half of the no-silent-drop contract. The
    // input side was validated before any read was touched; this checks the OUTPUT side — that every
    // supplied catalog copy actually came back out as an assigned copy. It is the only place a copy lost
    // inside the assignment stage could be seen, so it is an ERROR, not a log line.
    if let Some(ix) = &catalog_index {
        let mut jh = std::fs::File::create(format!("{}.family_join.tsv", args.out))?;
        writeln!(
            jh,
            "family_id\tcopy_index\tcopy_tid\tcatalog_family_id\tcatalog_copy_idx\tchrom\tstart\tend\tn_reads_hard"
        )?;
        for l in &join_rows {
            writeln!(jh, "{l}")?;
        }
        let emitted: HashSet<&str> = join_rows
            .iter()
            .filter_map(|l| l.split('\t').nth(2))
            .collect();
        let missing: Vec<&String> = ix.keys().filter(|t| !emitted.contains(t.as_str())).collect();
        if !missing.is_empty() {
            anyhow::bail!(
                "--families: {} of {} supplied copies did not come back as assigned copies ({}). A supplied \
                 copy must never be dropped; this is a bug in the assignment path, not a filter.",
                missing.len(),
                ix.len(),
                missing.iter().take(10).map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
            );
        }
        eprintln!(
            "[copy_assign] wrote {}.family_join.tsv ({} copies, all {} supplied catalog copies present)",
            args.out,
            join_rows.len(),
            ix.len()
        );
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
        writeln!(gf, "H\tVN:Z:1.1")?; // GFA 1.1 (W-lines carry the read walks)
        let mut segs: Vec<&String> = gfa_segs.iter().collect();
        segs.sort();
        for s in segs {
            writeln!(gf, "{}", s)?;
        }
        let mut links: Vec<&String> = gfa_links.iter().collect();
        links.sort();
        for l in links {
            writeln!(gf, "{}", l)?;
        }
        for p in &gfa_paths {
            writeln!(gf, "{}", p)?;
        }
        // reads threaded through the family graph (W-lines) — the shared-evidence flip made visible.
        for w in &gfa_walks {
            writeln!(gf, "{}", w)?;
        }
        // Bandage colour CSV: copies distinct, reads coloured by assigned copy (grey = tied/K=0, unresolvable).
        {
            let mut cf = std::fs::File::create(format!("{}.phase.gfa.colours.csv", args.out))?;
            writeln!(cf, "Name,Colour")?;
            for c in &gfa_colors {
                writeln!(cf, "{}", c)?;
            }
        }
        // Legend: status -> colour, de-duplicated across every family's copy-graph (first-seen order, which
        // is deterministic given the serial region/family drain order above).
        {
            let mut lf = std::fs::File::create(format!("{}.phase.gfa.legend.tsv", args.out))?;
            writeln!(lf, "status\tcolour")?;
            let mut seen: HashSet<&String> = HashSet::new();
            for r in &legend_rows {
                if seen.insert(r) {
                    writeln!(lf, "{}", r)?;
                }
            }
        }
        let n_phased = phased_read_lines.iter().filter(|l| !l.contains("\t-1\t")).count();
        eprintln!(
            "[copy_assign] phasing: {} blocks, {} haplotypes, {}/{} reads phased -> {}.phased_*.tsv + {}.phase.gfa ({} bubble-nodes, {} copy-paths, {} read-walks; Bandage colours -> {}.phase.gfa.colours.csv, legend -> {}.phase.gfa.legend.tsv)",
            phase_block_lines.len(), phased_hap_lines.len(), n_phased, phased_read_lines.len(),
            args.out, args.out, gfa_segs.len(), gfa_paths.len(), gfa_walks.len(), args.out, args.out
        );

        // v2: fold every family's exon presence/absence graph (built during the drain above) into
        // <out>.exon.gfa. Each exon's reference bases are fetched HERE, lazily, via genome_for (the
        // builder itself only lays out intervals) — a missing/uncovered stretch falls back to an N-run,
        // counted rather than silently faked (never claim sequence we didn't fetch).
        let n_seq_fallback = std::cell::Cell::new(0usize);
        let exon_seq = |ec: &rustle::vg_family::copy_graph::ExonClass| -> Vec<u8> {
            genome_for(&ec.chrom)
                .ok()
                .and_then(|g| g.fetch_sequence(&ec.chrom, ec.start, ec.end))
                .unwrap_or_else(|| {
                    n_seq_fallback.set(n_seq_fallback.get() + 1);
                    vec![b'N'; (ec.end - ec.start) as usize]
                })
        };
        let mut eg_file = std::fs::File::create(format!("{}.exon.gfa", args.out))?;
        writeln!(eg_file, "H\tVN:Z:1.1")?; // one header; each family's S/L/P folded below (fids are unique, no dedup needed)
        for fam_eg in &exon_graphs {
            for line in fam_eg.to_gfa(&exon_seq).lines().skip(1) {
                // skip the per-family embedded "H\tVN:Z:1.1" line — the file header above already covers it
                writeln!(eg_file, "{}", line)?;
            }
        }
        {
            let mut cf = std::fs::File::create(format!("{}.exon.gfa.colours.csv", args.out))?;
            writeln!(cf, "Name,Colour")?;
            for fam_eg in &exon_graphs {
                for row in fam_eg.colours_csv().lines() {
                    writeln!(cf, "{}", row)?;
                }
            }
        }
        {
            let mut lf = std::fs::File::create(format!("{}.exon.gfa.legend.tsv", args.out))?;
            writeln!(lf, "status\tcolour")?;
            let mut seen: HashSet<String> = HashSet::new();
            for fam_eg in &exon_graphs {
                for row in fam_eg.legend_tsv().lines() {
                    if seen.insert(row.to_string()) {
                        writeln!(lf, "{}", row)?;
                    }
                }
            }
        }
        eprintln!(
            "[copy_assign] exon graph: {} families -> {}.exon.gfa ({} exon-node(s) fell back to an N-run, no fetchable reference sequence); Bandage colours -> {}.exon.gfa.colours.csv, legend -> {}.exon.gfa.legend.tsv",
            exon_graphs.len(), args.out, n_seq_fallback.get(), args.out, args.out
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

    // ⭐ L6 --dump-star: each molecule's read-star proof (its own columns, its bases, every candidate's bases).
    if args.dump_star {
        let proofs = rustle::vg_family::copy_assign_pipeline::take_star_proofs();
        let mut sh = std::fs::File::create(format!("{}.star_reads.tsv", args.out))?;
        writeln!(sh, "read_name\tfamily_id\tstatus\tassigned_copy\tcatalog_copy_idx\tn_candidates\tcandidates\tn_cols\tcolumns")?;
        let mut n = 0usize;
        for r in &assign_rows {
            let Some(pf) = proofs.get(&r.read_name) else { continue };
            let ch = |o: Option<u8>| o.map(|b| b as char).unwrap_or('.');
            let cols: Vec<String> = pf
                .cols
                .iter()
                .enumerate()
                .map(|(j, &pos)| {
                    let cands: String = pf.alleles.iter().map(|al| ch(al.get(j).copied().flatten())).collect();
                    format!("{pos}:{}:{cands}", ch(pf.obs.get(j).copied().flatten()))
                })
                .collect();
            writeln!(
                sh,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.read_name, r.family_id, r.status, r.assigned_copy, r.catalog_copy_idx, r.n_candidates,
                pf.cand.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(","), pf.cols.len(), cols.join(",")
            )?;
            n += 1;
        }
        eprintln!("[copy_assign] --dump-star: {n} molecule proofs -> {}.star_reads.tsv (columns = read_pos:read_base:candidate_bases in `candidates` order, `.` = uncovered)", args.out);
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
        if args.igv {
            eprintln!(
                "[copy_assign] --igv: now run  python bench/igv_tracks.py --assignments {0}.assignments.tsv \
--bam <bam> --regions <regions> --out {0}  -> {0}.tagged.bam + {0}.copies.bed + {0}.psv.vcf (load in IGV)",
                args.out
            );
        }
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

        // Augment-and-linearize certificates (Task 4/5): one row per Stage-2-admitted candidate, report-only
        // (`--linearize`) or gating admission too (`--linearize-gate`, which implies the report). ONLY written
        // when the opt-in is set (`do_linearize`) -- plain `--absent-copies` skips the certificate entirely
        // (no minimap2 per candidate) and emits no `.linearize.tsv`, byte-identical to the pre-feature path.
        if do_linearize {
            let mut lh = std::fs::File::create(format!("{}.linearize.tsv", args.out))?;
            writeln!(
                lh,
                "family_id\tchrom\tstart\tend\tn_pool\tlinearized_frac_real\tmean_frac_decoy\tdelta\tperm_p\tverdict"
            )?;
            for (fam, cert, (chrom, start, end)) in &linearize_certs_all {
                writeln!(lh, "{}", linearize_tsv_row(fam, (chrom, *start, *end), cert))?;
            }
            eprintln!(
                "[copy_assign] {} linearize certificate(s) -> {}.linearize.tsv{}",
                linearize_certs_all.len(),
                args.out,
                if args.linearize_gate { " (--linearize-gate: non-LINEARIZES candidates demoted to .dna_needs.tsv)" } else { "" }
            );
        }
    }

    // --vg-realign: the re-align supplement's per-family/per-read decisions (report-only — not fed back
    // into the assignment). Only written when --vg-realign is set so an OFF run produces exactly the same
    // output files (cfg.vg_realign OFF also means fa.realign_records is always empty, so this is belt-
    // and-suspenders with the flag check).
    if args.vg_realign || args.vg_realign_correct {
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

    // ---- <out>.xfam_conflicts.tsv (report/abstain only) -------------------------------------------
    // The NEW information goes to a NEW file: adding a column to any existing output would break the OFF
    // arm's byte-identity, which is the gate this change is judged on.
    if xfam_mode != XfamMode::Off {
        let mut xh = std::fs::File::create(format!("{}.xfam_conflicts.tsv", args.out))?;
        writeln!(
            xh,
            "read_name\tstratum\tfamily_a\tcopy_a\tchrom_a\tstart_a\tend_a\tfamily_b\tcopy_b\tchrom_b\tstart_b\tend_b\tsame_record\tsep_bp\tdemoted"
        )?;
        for c in &xfam_conflicts {
            writeln!(
                xh,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                c.read_name, c.stratum, c.fid_a, c.copy_a, c.span_a.0, c.span_a.1, c.span_a.2,
                c.fid_b, c.copy_b, c.span_b.0, c.span_b.1, c.span_b.2, c.same_record, c.sep_bp, c.demoted
            )?;
        }
        // Every rate with its denominator (a rate without one is how four headlines died here).
        let n_pairs = xfam_conflicts.len();
        let n_shared = xfam_conflicts.iter().filter(|c| c.stratum == "shared_locus").count();
        let n_rt = xfam_conflicts.iter().filter(|c| c.stratum == "readthrough_span").count();
        let n_contra = xfam_conflicts.iter().filter(|c| c.stratum == "cross_family_contradiction").count();
        let contested_reads: std::collections::BTreeSet<&str> =
            xfam_conflicts.iter().map(|c| c.read_name.as_str()).collect();
        let contested_rows: std::collections::BTreeSet<(&str, &str)> = xfam_conflicts
            .iter()
            .flat_map(|c| [(c.read_name.as_str(), c.fid_a.as_str()), (c.read_name.as_str(), c.fid_b.as_str())])
            .collect();
        // Denominator = the assigned rows the OFF arm would have emitted. A demotion only ever turns
        // "assigned" into "ambiguous", so that is the post-demotion count plus the demotions — quoting the
        // POST-demotion count would shrink the denominator by exactly the numerator.
        let n_assigned = assign_rows.iter().filter(|r| r.status == "assigned").count()
            + if xfam_mode == XfamMode::Abstain { xfam_demote.len() } else { 0 };
        let distinct_reads: std::collections::BTreeSet<&str> =
            assign_rows.iter().map(|r| r.read_name.as_str()).collect();
        let demoted_mols: std::collections::BTreeSet<&str> =
            xfam_demote.iter().map(|(n, _, _)| n.as_str()).collect();
        // Distinct molecules with at least one assigned row IN THE OFF ARM — the honest denominator for
        // "how many molecules did this demote": the post-demotion set is smaller by the numerator.
        let distinct_assigned: std::collections::BTreeSet<&str> = assign_rows
            .iter()
            .filter(|r| r.status == "assigned")
            .map(|r| r.read_name.as_str())
            .chain(demoted_mols.iter().copied())
            .collect();
        eprintln!(
            "[xfam] RUSTLE_XFAM_RECONCILE={}: {n_pairs} contested assigned-copy pair(s) over \
             {}/{} distinct molecules — {n_shared} shared_locus, {n_rt} readthrough_span, \
             {n_contra} cross_family_contradiction (the only demoting stratum). \
             Contested (read,family) rows: {}. -> {}.xfam_conflicts.tsv",
            xfam_mode.as_str(), contested_reads.len(), distinct_reads.len(), contested_rows.len(), args.out
        );
        if xfam_mode == XfamMode::Abstain {
            let n_dem = xfam_demote.len();
            eprintln!(
                "[xfam] abstain: demoted {n_dem}/{n_assigned} assigned rows ({:.4}), \
                 {}/{} distinct assigned molecules ({:.4}), {n_dem}/{} CONTESTED rows ({:.4}). \
                 The third rate is the one that matters: a rule that only ever fires on \
                 multi-placement molecules is selection-biased by construction, so the run-wide rate \
                 alone hides where all of it lands. Demotion is a STATUS change only: no row is added \
                 or deleted, and .quant.tsv (n_reads_hard, abundance, ci95) is unmoved by design.",
                if n_assigned > 0 { n_dem as f64 / n_assigned as f64 } else { 0.0 },
                demoted_mols.len(), distinct_assigned.len(),
                if !distinct_assigned.is_empty() { demoted_mols.len() as f64 / distinct_assigned.len() as f64 } else { 0.0 },
                contested_rows.len(),
                if !contested_rows.is_empty() { n_dem as f64 / contested_rows.len() as f64 } else { 0.0 }
            );
        }
    }

    // ---- <out>.params.tsv (ALWAYS) ----------------------------------------------------------------
    // The M2 landing spot this binary lacked: it writes 24 output files and no params certificate, so an
    // ON and an OFF arm of ANY env-driven knob were previously indistinguishable from their outputs.
    // A NEW file changes no existing byte, so writing it unconditionally is compatible with the OFF gate.
    {
        let mut ph = std::fs::File::create(format!("{}.params.tsv", args.out))?;
        writeln!(ph, "key\tvalue")?;
        let mut row = |k: &str, v: String| -> Result<()> {
            writeln!(ph, "{k}\t{v}")?;
            Ok(())
        };
        row("xfam_reconcile", xfam_mode.as_str().to_string())?;
        row("posterior_prior", if prior_abundance { "abundance".into() } else { "uniform".to_string() })?;
        row("margin", format!("{}", args.margin))?;
        row("error_rate", format!("{}", args.error_rate))?;
        row("alpha", format!("{}", args.alpha))?;
        row("margin_gate", format!("{}", args.margin_gate))?;
        row("rna_editing_filter", format!("{}", !args.no_editing_filter))?;
        row("junction_conflict_abstain", format!("{}", args.junction_conflict_abstain))?;
        row("psv_genomic", format!("{}", args.psv_genomic))?;
        row("psv_read_filter", std::env::var("RUSTLE_PSV_READFILTER").unwrap_or_else(|_| "unset".into()))?;
        row("molecule_observations", format!("{}", args.molecule_observations && !args.no_molecule_observations))?;
        row("origin_rejected", format!("{}", assign_rows.iter().filter(|r| r.origin_rejected).count()))?;
        row("orphans", format!("{}", assign_rows.iter().filter(|r| r.origin_rejected && r.n_candidates == 0).count()))?;
        row("origin_substitutions_only", format!("{}", args.origin_substitutions_only))?;
        row("read_star_junctions", format!("{}", args.read_star_junctions))?;
        row("read_star_genomic", format!("{}", args.read_star_genomic && !args.read_star_unit))?;
        row("read_star_catalog_locus", format!("{}", !args.read_star_pad_locus))?;
        row("sole_candidate", format!("{}", !args.no_sole_candidate))?;
        row("sole_candidates", format!("{}", assign_rows.iter().filter(|r| r.status == "assigned" && r.n_candidates == 1 && r.contested).count()))?;
        row("placement_assign", format!("{}", args.molecule_observations && !args.no_molecule_observations && !args.no_placement_assign))?;
        row("placement_assigned", format!("{}", placement_assigned_total))?;
        row("placement_first", format!("{}", args.placement_first))?;
        row("contested_rows", format!("{}", assign_rows.iter().filter(|r| r.contested).count()))?;
        row("dump_star", format!("{}", args.dump_star))?;
        row("read_star_hit_in_unit", format!("{}", !args.no_read_star_hit_in_unit))?;
        row("read_star_two_form", format!("{}", !args.read_star_genomic_only))?;
        row("junction_conflicts", format!("{}", assign_rows.iter().filter(|r| r.junction_conflict).count()))?;
        row("edit_rate", format!("{}", args.edit_rate))?;
        row("iterative_prune", format!("{}", args.iterative_prune))?;
        row("families", args.families.clone().unwrap_or_else(|| "NONE".to_string()))?;
        row("copies_fa", args.copies_fa.clone().unwrap_or_else(|| "NONE".to_string()))?;
        row("dump_psv", format!("{}", args.dump_psv))?;
        row("phase", format!("{}", args.phase))?;
        row("posterior", format!("{}", args.posterior))?;
        row("em", format!("{}", args.em))?;
        row("gtf", format!("{}", args.gtf))?;
        eprintln!("[copy_assign] wrote {}.params.tsv (run certificate)", args.out);
    }

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
            // QUANTIFY the structural warning. A flagged pair is a STRUCTURE ("these two copies share
            // sequence"); what makes it consequential is how many MOLECULES it actually double-claims —
            // the same molecules the cross-family reconciliation's `shared_locus` stratum sees, and the
            // reason that stratum is reported rather than demoted (demoting it would strip these copies
            // of most of their hard support and charge an O1 partition defect to O2's abstention rate).
            // `catalog` is built from `quant_rows` in order, so index i IS quant row i.
            let mut qidx: std::collections::BTreeMap<(&str, usize), usize> = std::collections::BTreeMap::new();
            for (i, r) in quant_rows.iter().enumerate() {
                qidx.insert((r.family_id.as_str(), r.copy_index), i);
            }
            let mut claims: std::collections::BTreeMap<&str, Vec<usize>> = std::collections::BTreeMap::new();
            for r in assign_rows.iter().filter(|r| r.status == "assigned") {
                if let Some(&i) = qidx.get(&(r.family_id.as_str(), r.assigned_copy)) {
                    claims.entry(r.read_name.as_str()).or_default().push(i);
                }
            }
            let mut double: std::collections::BTreeMap<(usize, usize), usize> = std::collections::BTreeMap::new();
            for v in claims.values().filter(|v| v.len() >= 2) {
                for a in 0..v.len() {
                    for b in (a + 1)..v.len() {
                        let (lo2, hi2) = (v[a].min(v[b]), v[a].max(v[b]));
                        *double.entry((lo2, hi2)).or_insert(0) += 1;
                    }
                }
            }
            let n_double: usize = flagged
                .iter()
                .map(|&(i, j, _, _)| double.get(&(i.min(j), i.max(j))).copied().unwrap_or(0))
                .sum();
            eprintln!(
                "[copy_assign]   ... those pairs double-claim {n_double} assigned molecule(s) in total \
                 (one molecule counted once per flagged pair it is assigned to both sides of)."
            );
            for &(i, j, recip, kind) in flagged.iter().take(10) {
                let n = double.get(&(i.min(j), i.max(j))).copied().unwrap_or(0);
                eprintln!(
                    "[copy_assign]   {kind:?} recip={recip:.2}  {}/{}:{}-{}  vs  {}/{}-{}  double_claimed_molecules={n}",
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linearize_tsv_row_formats() {
        use rustle::vg_family::linearize::{LinearizeCertificate, Verdict};
        let c = LinearizeCertificate { n_pool: 40, linearized_frac_real: 0.82, mean_frac_decoy: 0.01,
            delta: 0.81, perm_p: 0.05, verdict: Verdict::Linearizes };
        let row = linearize_tsv_row("GWFAM1", ("chr9", 100, 200), &c);
        assert_eq!(row, "GWFAM1\tchr9\t100\t200\t40\t0.820\t0.010\t0.810\t0.0500\tLINEARIZES");
    }

    #[test]
    fn read_lambda_file_parses_the_scalar() {
        let dir = std::env::temp_dir();
        let p = dir.join(format!("rustle_lam_test_{}.tsv", std::process::id()));
        std::fs::write(&p, "lambda_global\tn_single_copy_loci\n25.5\t1234\n").unwrap();
        assert_eq!(read_lambda_file(p.to_str().unwrap()), Some(25.5));
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn read_lambda_file_none_on_na_or_missing() {
        let dir = std::env::temp_dir();
        let p = dir.join(format!("rustle_lam_na_{}.tsv", std::process::id()));
        std::fs::write(&p, "lambda_global\tn_single_copy_loci\nNA\t0\n").unwrap();
        assert_eq!(read_lambda_file(p.to_str().unwrap()), None);
        assert_eq!(read_lambda_file("/nonexistent/path.tsv"), None);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn resolve_lambda_precedence_explicit_over_file() {
        assert_eq!(resolve_lambda(Some(30.0), Some(25.0)), Some(30.0));
        assert_eq!(resolve_lambda(None, Some(25.0)), Some(25.0));
        assert_eq!(resolve_lambda(None, None), None);
    }

    #[test]
    fn anchored_support_counts_only_unique_mappers_of_the_copy() {
        // reads: (best_copy, mapq) = (0,60) (0,0) (1,60) (0,60). Copy 0's unique-mapper support is the two
        // mapq>0 reads (the mapq=0 read is a tie-break-arbitrary primary and does not count).
        let bcs = [0usize, 0, 1, 0];
        let mqs = [60u8, 0, 60, 60];
        assert_eq!(anchored_support(&bcs, &mqs, 0), 2, "copy 0: two mapq>0 reads, the mapq=0 excluded");
        assert_eq!(anchored_support(&bcs, &mqs, 1), 1, "copy 1: one mapq>0 read");
        assert_eq!(anchored_support(&bcs, &mqs, 2), 0, "copy 2: no reads");
    }

    #[test]
    fn anchored_support_is_zero_when_every_primary_is_tied() {
        // The TSPY case: every assigned read is mapq=0 -> no copy has any anchored support -> none invariant.
        let bcs = [0usize, 1, 2, 0, 1];
        let mqs = [0u8, 0, 0, 0, 0];
        for ci in 0..3 {
            assert_eq!(anchored_support(&bcs, &mqs, ci), 0);
        }
    }

    #[test]
    fn tie_invariant_threshold_is_the_locus_gate() {
        use rustle::vg_family::denovo_assemble::GATE_MIN_READS;
        // The certificate boolean is anchored >= GATE_MIN_READS (=3): 2 -> false, 3 -> true.
        assert!(!(2u32 >= GATE_MIN_READS));
        assert!(3u32 >= GATE_MIN_READS);
    }

    #[test]
    fn build_copy_graph_maps_family_to_graph() {
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.n_copies = 2;
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.psv_col_pos = vec![Some(100), Some(200)];
        fa.copy_psv_alleles = vec![vec![Some(b'A'), Some(b'A')], vec![Some(b'A'), Some(b'G')]];
        fa.read_psv_obs = vec![];
        fa.assignments = vec![];
        // injected reference: A at both positions
        let ref_base = |_c: &str, _p: u64| Some(b'A');
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_copy_graph("CAFAM0", &fa, ref_base, &[], None, &eff);
        assert_eq!(g.columns.len(), 2);
        assert_eq!(g.copies.len(), 2);
        assert_eq!(g.columns[0].ref_allele, Some(b'A'));
        // graph renders the reference walk + a divergent copy
        let gfa = g.to_gfa();
        assert!(gfa.contains("P\tCAFAM0_REFERENCE"));
        assert!(gfa.contains("CAFAM0_c1_G+"));
    }

    // A discovery_coupled Assignment for read `ri` pinned to copy `best_copy` (status Assigned) — the ONLY
    // signal that makes a copy absent. Mirrors the default Assignment (copy_assign.rs) with the flag flipped.
    fn coupled_assignment(ri: usize, best_copy: usize) -> (usize, rustle::vg_family::copy_assign::Assignment) {
        (ri, rustle::vg_family::copy_assign::Assignment {
            best_copy,
            log_lr_margin: 10.0,
            n_decisive: 1,
            resolvable: true,
            status: AssignStatus::Assigned,
            p_value: 0.0,
            min_p_value: 0.0,
            discovery_coupled: true,
            junction_conflict: false,
            origin_rejected: false,
            n_candidates: 0,
            posterior: vec![],
        })
    }

    #[test]
    fn build_copy_graph_coupled_copy_colocated_is_absent_collapsed() {
        // Absence is driven by a discovery_coupled read, NOT the collapsed/rescued counts. copy1 has a
        // coupled read AND its span overlaps copy0's span => AbsentCollapsed (hidden co-located haplotype).
        // copy0 (no coupled read) stays AnnotationUnknown (no _ABSENT).
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.n_copies = 2;
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.psv_col_pos = vec![Some(100)];
        // copy0 = reference allele A; copy1 = divergent G (the absent, co-located copy).
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'G')]];
        // overlapping spans on the same chrom => co-located.
        fa.copy_spans = vec![("chr1".into(), 1000, 2000), ("chr1".into(), 1500, 2500)];
        fa.read_psv_obs = vec![vec![Some(b'G')]];
        fa.assignments = vec![coupled_assignment(0, 1)]; // one read, coupled to copy1
        let ref_base = |_c: &str, _p: u64| Some(b'A');
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_copy_graph("CAFAM0", &fa, ref_base, &["read0".to_string()], None, &eff);
        assert_eq!(g.copies.len(), 2);
        assert_eq!(g.copies[1].status, CopyStatus::AbsentCollapsed);
        assert!(g.copies[1].status.is_absent());
        let gfa = g.to_gfa();
        assert!(gfa.contains("_copy1_ABSENT"), "coupled co-located copy must render _ABSENT:\n{}", gfa);
        assert!(gfa.contains("ST:Z:absent-collapsed"), "expected ST:Z:absent-collapsed in:\n{}", gfa);
        // copy 0 is NOT absent (no coupled read) => AnnotationUnknown, no _ABSENT.
        assert_eq!(g.copies[0].status, CopyStatus::AnnotationUnknown);
        assert!(!g.copies[0].status.is_absent());
        assert!(!gfa.contains("_copy0_ABSENT"), "non-coupled copy must NOT be _ABSENT:\n{}", gfa);
    }

    #[test]
    fn build_copy_graph_coupled_copy_dispersed_is_absent_divergent() {
        // copy1 has a coupled read but its span is DISJOINT from copy0's (different chrom) => AbsentDivergent
        // (dispersed, no overlapping in-genome copy).
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.n_copies = 2;
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.psv_col_pos = vec![Some(100)];
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'G')]];
        // disjoint spans (different chrom) => dispersed.
        fa.copy_spans = vec![("chr1".into(), 1000, 2000), ("chr9".into(), 1000, 2000)];
        fa.read_psv_obs = vec![vec![Some(b'G')]];
        fa.assignments = vec![coupled_assignment(0, 1)];
        let ref_base = |_c: &str, _p: u64| Some(b'A');
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_copy_graph("CAFAM0", &fa, ref_base, &["read0".to_string()], None, &eff);
        assert_eq!(g.copies[1].status, CopyStatus::AbsentDivergent);
        let gfa = g.to_gfa();
        assert!(gfa.contains("_copy1_ABSENT"), "dispersed coupled copy must render _ABSENT:\n{}", gfa);
        assert!(gfa.contains("ST:Z:absent-divergent"), "expected ST:Z:absent-divergent in:\n{}", gfa);
        assert_eq!(g.copies[0].status, CopyStatus::AnnotationUnknown);
    }

    #[test]
    fn build_copy_graph_gstm_regression_no_coupled_reads_no_absent() {
        // REGRESSION for the GSTM bug: collapsed_copies (9) >> n_copies (3) with NO discovery_coupled reads.
        // The old code let absent_tail_start underflow to 0 and mislabeled ALL 3 in-genome copies _ABSENT.
        // Correct behavior: no coupled read => NO copy is absent; all are AnnotationUnknown, no _ABSENT.
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.n_copies = 3;
        fa.copy_tids = vec!["c0".into(), "c1".into(), "c2".into()];
        fa.psv_col_pos = vec![Some(100)];
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'C')], vec![Some(b'G')]];
        fa.read_psv_obs = vec![];
        fa.assignments = vec![]; // NO discovery_coupled reads (no --absent-copies)
        fa.collapsed_copies = 9;  // diagnostic count, far exceeds n_copies — must NOT force absence.
        fa.rescued_copies = 0;
        let ref_base = |_c: &str, _p: u64| Some(b'A');
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_copy_graph("CAFAM0", &fa, ref_base, &[], None, &eff);
        assert_eq!(g.copies.len(), 3);
        for ci in 0..3 {
            assert_eq!(g.copies[ci].status, CopyStatus::AnnotationUnknown, "copy {} must be in-genome", ci);
            assert!(!g.copies[ci].status.is_absent(), "copy {} must NOT be absent", ci);
        }
        let gfa = g.to_gfa();
        assert!(!gfa.contains("_ABSENT"), "no copy may render _ABSENT when no read is discovery_coupled:\n{}", gfa);
        assert!(!gfa.contains("absent-collapsed"), "no absent-collapsed tag expected:\n{}", gfa);
    }

    #[test]
    fn build_copy_graph_fills_mi_from_copy_map_identity() {
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.psv_col_pos = vec![Some(100)];
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'G')]];
        // copy1 has a remap identity but NO discovery_coupled read: copy_map_identity feeds the MI tag
        // ONLY and must NOT flip absent status — so copy1 stays IN-GENOME (AnnotationUnknown, no --gff).
        fa.copy_map_identity = vec![None, Some(0.952)];
        fa.assignments = vec![]; // no discovery_coupled reads => no absence
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_copy_graph("CAFAM0", &fa, |_c, _p| Some(b'A'), &[], None, &eff);
        // copy1 is in-genome, NOT absent (the LOCK on the copy_map_identity-drives-absence bug).
        assert_eq!(g.copies[1].status, CopyStatus::AnnotationUnknown, "copy_map_identity alone must NOT make a copy absent");
        assert!(!g.copies[1].status.is_absent());
        let gfa = g.to_gfa();
        // P-line name is EXACTLY CAFAM0_copy1 (the trailing tab excludes CAFAM0_copy1_ABSENT).
        let c1 = gfa.lines().find(|l| l.starts_with("P\tCAFAM0_copy1\t")).expect("copy1 must be P\\tCAFAM0_copy1 (NOT _ABSENT)");
        assert!(!gfa.contains("CAFAM0_copy1_ABSENT"), "copy_map_identity must not render _ABSENT:\n{}", gfa);
        assert!(c1.contains("MI:f:0.952"), "copy1 MI missing: {}", c1);
        assert!(c1.contains("ST:Z:annotation-unknown"), "copy1 must carry a non-absent ST:Z: status: {}", c1);
        let c0 = gfa.lines().find(|l| l.starts_with("P\tCAFAM0_copy0\t")).unwrap();
        assert!(!c0.contains("MI:f:"), "copy0 must omit MI: {}", c0);
    }

    #[test]
    fn build_exon_graph_makes_copy_specific_arm() {
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.copy_spans = vec![("chr1".into(), 0, 400), ("chr1".into(), 0, 400)];
        // copy0 introns skip 100-300 (exons 0-100, 300-400); copy1 has an extra exon (exons 0-100,150-250,300-400)
        fa.copy_introns = vec![ vec![(100,300)], vec![(100,150),(250,300)] ];
        fa.copy_map_identity = vec![None, Some(0.95)];
        // copy1 is reference-ABSENT via a discovery_coupled read (the v1 absence mechanism — NOT
        // copy_map_identity, which feeds only the MI tag); its span overlaps copy0 => AbsentCollapsed.
        fa.assignments = vec![coupled_assignment(0, 1)];
        let eff: Vec<AssignStatus> = fa.assignments.iter().map(|(_, a)| a.status).collect();
        let g = build_exon_graph("CAFAM0", &fa, None, &eff);
        // copy1 walks one more class than copy0
        assert!(g.copies[1].exon_nodes.len() > g.copies[0].exon_nodes.len());
        let gfa = g.to_gfa(|ec| vec![b'A'; (ec.end-ec.start) as usize]);
        assert!(gfa.contains("P\tCAFAM0_REFERENCE"));
        assert!(gfa.lines().any(|l| l.starts_with("P\tCAFAM0_copy1_ABSENT")));
    }

    #[test]
    fn annotation_axis_from_intervals() {
        use rustle::vg_family::denovo_pipeline::FamilyAssignment;
        use rustle::vg_family::copy_graph::CopyStatus;
        let mut fa = FamilyAssignment::empty();
        fa.chrom = "chr1".into();
        fa.copy_tids = vec!["c0".into(), "c1".into()];
        fa.copy_spans = vec![("chr1".into(), 1000, 2000), ("chr1".into(), 5000, 6000)];
        fa.psv_col_pos = vec![Some(1500)];
        fa.copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'A')]];
        // annotation covers only copy0's span
        let ann = vec![("chr1".to_string(), 900u64, 2100u64)];
        let s0 = annotation_status(&fa, 0, Some(&ann));
        let s1 = annotation_status(&fa, 1, Some(&ann));
        assert_eq!(s0, CopyStatus::InGenomeAnnotated);
        assert_eq!(s1, CopyStatus::InGenomeUnannotated);
        assert_eq!(annotation_status(&fa, 1, None), CopyStatus::AnnotationUnknown);
    }

    #[test]
    fn parse_annotation_bed_with_numeric_name() {
        // Regression: a BED6 line whose NAME column (col 3) is numeric must NOT be misread as GFF (which would
        // yield (4,5) from cols 3/4). Extension-first dispatch (`.bed` => BED cols 0/1/2) keeps it correct.
        let dir = std::env::temp_dir();
        let bed = dir.join(format!("rustle_task8_ann_{}.bed", std::process::id()));
        std::fs::write(&bed, "chr1\t1000\t2000\t5\t0\t+\n").unwrap();
        let got = parse_annotation(bed.to_str().unwrap()).unwrap();
        std::fs::remove_file(&bed).ok();
        assert_eq!(got, vec![("chr1".to_string(), 1000u64, 2000u64)], "BED numeric-name line must NOT parse as (4,5)");

        // A GFF/GTF file still parses cols 3/4 with the 1-based -> 0-based start conversion.
        let gff = dir.join(format!("rustle_task8_ann_{}.gff3", std::process::id()));
        std::fs::write(&gff, "##gff-version 3\nchr1\tsrc\tgene\t1000\t2000\t.\t+\t.\tID=g1\n").unwrap();
        let got = parse_annotation(gff.to_str().unwrap()).unwrap();
        std::fs::remove_file(&gff).ok();
        assert_eq!(got, vec![("chr1".to_string(), 999u64, 2000u64)], "GFF cols 3/4, 1-based start -> 0-based");
    }
}
