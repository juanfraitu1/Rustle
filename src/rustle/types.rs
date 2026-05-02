//! Core data structures: Bundle, Junction, Bundlenode, RunConfig.
//!
//! ## Transcript extraction mode (find_transcripts)
//! - **Eonly** (config.eonly, -e): Only guide flow runs (ref: only guides_pushmaxflow). No parse_trf.
//!   Output is only guide-matched transcripts.
//! - **Long-read** (default): extract_transcripts with long order only (get_trf_long).

use crate::constants::{LONGINTRON, LOWCOV, MISMATCHFRAC, SSERROR};
use std::cmp::Ordering;
use std::sync::Arc;

/// Fast deterministic hasher: FxHash is ~3-5x faster than SipHash for integer keys.
/// Deterministic because FxHash has no randomization seed.
pub type FixedBuild = fxhash::FxBuildHasher;
pub type DetHashMap<K, V> = std::collections::HashMap<K, V, FixedBuild>;
pub type DetHashSet<T> = std::collections::HashSet<T, FixedBuild>;

/// naming convention: generic vector.
pub type GVec<T> = Vec<T>;
/// naming convention: pointer vector (modeled as plain Vec in Rust).
pub type GPVec<T> = Vec<T>;
/// naming convention: sorted/unique array container (modeled as Vec in Rust).
pub type GArray<T> = Vec<T>;

/// Guide edge (`GEdge` in `header`): boundary value + opposite boundary + strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GEdge {
    /// One boundary coordinate.
    pub val: u64,
    /// Mate boundary coordinate from the same guide edge.
    pub endval: u64,
    /// -1, 0, +1 strand.
    pub strand: i8,
}

impl GEdge {
    pub fn new(val: u64, endval: u64, strand: i8) -> Self {
        Self {
            val,
            endval,
            strand,
        }
    }

    /// comment: when `val < endval`, this marks a start boundary.
    pub fn is_start(self) -> bool {
        self.val < self.endval
    }
}

impl Ord for GEdge {
    fn cmp(&self, other: &Self) -> Ordering {
        self.val
            .cmp(&other.val)
            .then_with(|| self.strand.cmp(&other.strand))
    }
}

impl PartialOrd for GEdge {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Assembly mode: long-read only.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyMode {
    /// Long-read assembly (-L). All reads treated as long; abundance only.
    LongRead,
}

impl AssemblyMode {
    pub fn is_long_read(self) -> bool {
        true
    }
    pub fn is_mixed(self) -> bool {
        false
    }
}

/// Genomic bundle (locus): chrom, range, reads, junctions.
#[derive(Debug, Clone)]
pub struct Bundle {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    /// Read alignments (reference start, reference end, exons as (start,end), weight, is_reverse)
    pub reads: Vec<BundleRead>,
    pub junction_stats: JunctionStats,
    /// Pre-computed bundlenodes (e.g., faithful implementation)
    pub bundlenodes: Option<CBundlenode>,
    /// Read-to-bundlenode mappings (indexes into bundlenodes)
    pub read_bnodes: Option<Vec<Vec<usize>>>,
    /// Color assignments for bundlenodes
    pub bnode_colors: Option<Vec<usize>>,
    /// True for bundles synthesized by VG rescue from unmapped reads
    /// (see src/rustle/vg_hmm/rescue.rs). Synthetic bundles bypass
    /// transcript_isofrac and cross-bundle pairwise-contained filters.
    pub synthetic: bool,
    /// Diagnostic classification bucket for this synthetic bundle (None for
    /// real bundles). Set by `vg_hmm::diagnostic::classify_internal` and
    /// optionally refined by `classify_external` when
    /// `config.vg_rescue_diagnostic == true`.
    pub rescue_class: Option<crate::vg_hmm::diagnostic::RescueClass>,
}

#[derive(Debug, Clone)]
pub struct BundleRead {
    /// Stable unique ID assigned at creation; survives coordinate mutations.
    pub read_uid: u64,
    /// Original BAM QNAME (for standard mate linking).
    pub read_name: Arc<str>,
    /// Cached 64-bit hash of `read_name` for frequent signature/dedup lookups.
    pub read_name_hash: u64,
    /// Reference sequence ID of this alignment.
    pub ref_id: Option<usize>,
    /// Mate reference sequence ID (if paired).
    pub mate_ref_id: Option<usize>,
    /// Mate alignment start, 0-based.
    pub mate_start: Option<u64>,
    /// HI tag (alignment hit index); used to disambiguate mate hashing.
    pub hi: u32,
    pub ref_start: u64,
    pub ref_end: u64,
    pub exons: Vec<(u64, u64)>,
    /// Active splice junction coordinates used by the pipeline.
    pub junctions: Vec<Junction>,
    /// Per-junction usability bit carried with the read through repair steps.
    pub junction_valid: Vec<bool>,
    /// Raw splice junctions from exon boundaries.
    pub junctions_raw: Vec<Junction>,
    /// Deletion-aware splice junctions extracted from CIGAR (juncsdel-adjusted boundaries).
    pub junctions_del: Vec<Junction>,
    pub weight: f64,
    pub is_reverse: bool,
    /// Transcript strand inferred from alignment tags (e.g. minimap2 `ts`) when available.
    /// Falls back to alignment orientation if absent.
    pub strand: char,
    /// PolyA/T at 5' (read) end (aligned OR soft-clip)
    pub has_poly_start: bool,
    /// PolyA/T at 3' (read) end (aligned OR soft-clip)
    pub has_poly_end: bool,
    /// PolyA/T in the aligned sequence at 5' end. aligned_polyT.
    pub has_poly_start_aligned: bool,
    /// PolyA/T in the 5' soft-clipped tail. unaligned_polyT.
    pub has_poly_start_unaligned: bool,
    /// PolyA/T in the aligned sequence at 3' end. aligned_polyA.
    pub has_poly_end_aligned: bool,
    /// PolyA/T in the 3' soft-clipped tail. unaligned_polyA.
    pub has_poly_end_unaligned: bool,
    /// unaligned_polyT support count (long-read CPAS evidence multiplicity).
    pub unaligned_poly_t: u16,
    /// unaligned_polyA support count (long-read CPAS evidence multiplicity).
    pub unaligned_poly_a: u16,
    /// check_last_exon_polyA >= 0.8: last exon's bases in read are ≥80% A.
    /// When true (and !ovlpguide), the entire last exon is removed before junction build.
    pub has_last_exon_polya: bool,
    /// check_first_exon_polyT >= 0.8: first exon's bases in read are ≥80% T.
    /// When true (and !ovlpguide), the entire first exon is removed before junction build.
    pub has_first_exon_polyt: bool,
    /// Query (read) length in bp; used in mixed mode to classify long vs short (long if len >= long_read_min_len).
    pub query_length: Option<u64>,
    /// Left soft-clip length (query side, bp). clipL; used for mismatch penalty and poly detection.
    pub clip_left: u32,
    /// Right soft-clip length (query side, bp). clipR; used for mismatch penalty and poly detection.
    pub clip_right: u32,
    /// NH alignment hit count; used for standard multi-mapper mismatch penalty (nh > 2).
    pub nh: u32,
    /// Edit distance (NM/nM tag when present). Used for per-junction mismatch accumulation.
    pub nm: u32,
    /// MD tag for mismatch/deletion anchor checks near splice sites.
    pub md: Option<String>,
    /// Reference positions of CIGAR insertions (anchor checks near splice sites).
    pub insertion_sites: Vec<u64>,
    /// Unitig read (from -U option). isunitig;YK tag for unitig coverage.
    pub unitig: bool,
    /// Unitig coverage from YK tag (unitig_cov=brec.tag_float("YK")).
    pub unitig_cov: f64,
    /// Read count from YC tag (rdcount=brec.tag_int("YC")).
    pub read_count_yc: f64,
    /// the original algorithm `countFragment` contribution to `frag_len` (sum exon_len / NH) for this read entry.
    /// Aggregated across all redundant alignments merged into this entry.
    pub countfrag_len: f64,
    /// the original algorithm `countFragment` contribution to `num_fragments` for this read entry.
    /// Aggregated across all redundant alignments merged into this entry.
    pub countfrag_num: f64,
    /// Sum of per-alignment `rdcount` where junction mismatch gate is true
    /// (`mismatch || nh>2` in processRead).
    /// Applied per junction when initializing BundleData-like junction stats.
    pub junc_mismatch_weight: f64,
    /// Pair read indices (pair_idx). Index into read list of paired reads.
    pub pair_idx: Vec<usize>,
    /// Pair read counts (pair_count). Weight of each pair connection.
    pub pair_count: Vec<f64>,
    // ── VG gene family fields (populated when --vg is active) ────────────
    /// Mapping quality (MAPQ). Used for filtering supplementary alignments in VG mode.
    pub mapq: u8,
    /// Per-base mismatches extracted from MD tag: (reference_position, query_base).
    /// Only populated when `--vg-snp` is active.
    pub mismatches: Vec<(u64, u8)>,
    /// Haplotype tag (HP:i:1 or HP:i:2) from phased BAMs. None = unphased.
    pub hp_tag: Option<u8>,
    /// Phase set ID (PS tag) for phased BAMs.
    pub ps_tag: Option<u32>,
    /// True if this record is the primary alignment (not secondary/supplementary).
    /// Populated at ingest; used in VG mode to filter secondaries out of
    /// non-family bundles after family discovery.
    pub is_primary_alignment: bool,
}

/// `CPred` transport for longtrim boundary points (predno, cov).
#[derive(Debug, Clone, Copy)]
pub struct CPred {
    pub predno: u64,
    pub cov: f64,
}

/// Read boundary for longtrim (CPred predno+cov; original algorithm ReadBoundary).
/// lstart = list of (position, coverage) for read starts; lend = read ends.
#[derive(Debug, Clone, Copy)]
pub struct ReadBoundary {
    pub pos: u64,
    pub cov: f64,
}

impl From<CPred> for ReadBoundary {
    fn from(v: CPred) -> Self {
        Self {
            pos: v.predno,
            cov: v.cov,
        }
    }
}

impl From<ReadBoundary> for CPred {
    fn from(v: ReadBoundary) -> Self {
        Self {
            predno: v.pos,
            cov: v.cov,
        }
    }
}

/// Prediction object compatible with `CPrediction` (bundle.h).
/// This is the legacy transport model; assembly still primarily uses `Transcript`.
#[derive(Debug, Clone)]
pub struct CPrediction {
    pub geneno: i32,
    /// Guide transcript id for equivalent annotation (`t_eq`), when available.
    pub t_eq: Option<String>,
    pub start: u64,
    pub end: u64,
    pub cov: f64,
    pub longcov: f64,
    pub strand: char,
    /// Transcript length sentinel/field (`tlen`; long-read uses negative conventions).
    pub tlen: i32,
    pub flag: bool,
    /// Linked nascent prediction index in the same prediction vector.
    pub linkpred: Option<usize>,
    pub exons: Vec<(u64, u64)>,
    pub exoncov: Vec<f64>,
    pub mergename: Option<String>,
    pub hardstart: bool,
    pub hardend: bool,
}

impl CPrediction {
    pub fn new(geneno: i32, start: u64, end: u64, cov: f64, strand: char, tlen: i32) -> Self {
        Self {
            geneno,
            t_eq: None,
            start,
            end,
            cov,
            longcov: 0.0,
            strand,
            tlen,
            flag: true,
            linkpred: None,
            exons: Vec::new(),
            exoncov: Vec::new(),
            mergename: None,
            hardstart: false,
            hardend: false,
        }
    }

    pub fn exonic_len(&self) -> u64 {
        self.exons
            .iter()
            .map(|(s, e)| e.saturating_sub(*s))
            .sum::<u64>()
    }
}

/// interval node payload (`CExon` used in `CMaxIntv.node`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CExon {
    pub predno: usize,
    pub exonno: usize,
    pub exoncov: f64,
}

impl CExon {
    pub fn new(predno: usize, exonno: usize, exoncov: f64) -> Self {
        Self {
            predno,
            exonno,
            exoncov,
        }
    }
}

/// Max-overlap interval list node compatible with `CMaxIntv` (print_predcluster internals).
#[derive(Debug, Clone)]
pub struct CMaxIntv {
    pub start: u64,
    pub end: u64,
    pub node: Vec<CExon>,
    pub cov: f64,
    pub next: Option<Box<CMaxIntv>>,
}

impl CMaxIntv {
    pub fn new(start: u64, end: u64) -> Self {
        Self {
            start,
            end,
            node: Vec::new(),
            cov: 0.0,
            next: None,
        }
    }

    pub fn with_node(node: Vec<CExon>, start: u64, end: u64, cov: f64) -> Self {
        Self {
            start,
            end,
            node,
            cov,
            next: None,
        }
    }
}

/// Merge component member compatible with `CTrInfo` (header).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct CTrInfo {
    pub trno: usize,
    pub abundance: f64,
    pub penalty: f64,
}

impl CTrInfo {
    pub fn new(trno: usize, abundance: f64, penalty: f64) -> Self {
        Self {
            trno,
            abundance,
            penalty,
        }
    }
}

/// Memoized component payload compatible with `CComponent` (header).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct CComponent {
    pub size: f64,
    pub set: Vec<usize>,
}

impl CComponent {
    pub fn new(size: f64, set: Vec<usize>) -> Self {
        Self { size, set }
    }
}

/// Junction: donor (exon end) and acceptor (next exon start). 0-based.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Junction {
    pub donor: u64,
    pub acceptor: u64,
}

impl Junction {
    pub fn new(donor: u64, acceptor: u64) -> Self {
        Self { donor, acceptor }
    }
}

/// compatible junction record (`CJunction` shape).
/// This mirrors the fields used by junction filtering/gating logic.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CJunction {
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub nreads: f64,
    pub nreads_good: f64,
    pub leftsupport: f64,
    pub rightsupport: f64,
    pub nm: f64,
    pub mm: f64,
    pub guide_match: bool,
    pub consleft: i8,
    pub consright: i8,
    pub rcount: u32,
}

impl CJunction {
    pub fn new(start: u64, end: u64, strand: i8) -> Self {
        Self {
            start,
            end,
            strand,
            nreads: 0.0,
            nreads_good: 0.0,
            leftsupport: 0.0,
            rightsupport: 0.0,
            nm: 0.0,
            mm: 0.0,
            guide_match: false,
            consleft: -1,
            consright: -1,
            rcount: 0,
        }
    }

    #[inline]
    pub fn key(self) -> Junction {
        Junction::new(self.start, self.end)
    }

    pub fn from_parts(j: Junction, st: &JunctionStat) -> Self {
        Self {
            start: j.donor,
            end: j.acceptor,
            strand: st.strand.unwrap_or(0),
            nreads: st.mrcount,
            nreads_good: st.nreads_good,
            leftsupport: st.leftsupport,
            rightsupport: st.rightsupport,
            nm: st.nm,
            mm: st.mm,
            guide_match: st.guide_match,
            consleft: st.consleft,
            consright: st.consright,
            rcount: st.rcount,
        }
    }
}

/// Deterministic conversion from internal Junction+JunctionStat map
/// to the original algorithm-shaped `CJunction` records.
pub fn junction_stats_to_cjunctions(junction_stats: &JunctionStats) -> Vec<CJunction> {
    let mut out: Vec<CJunction> = junction_stats
        .iter()
        .map(|(j, st)| CJunction::from_parts(*j, st))
        .collect();
    out.sort_by_key(|j| (j.start, j.end, j.strand));
    out
}

/// Convert the original algorithm-shaped `CJunction` records back to internal `JunctionStats`.
/// When duplicate (start,end) entries exist, numeric support fields are aggregated.
pub fn cjunctions_to_junction_stats(cjunctions: &[CJunction]) -> JunctionStats {
    let mut out: JunctionStats = Default::default();
    for cj in cjunctions {
        let key = Junction::new(cj.start, cj.end);
        let entry = out.entry(key).or_insert_with(JunctionStat::default);
        // IMPORTANT: In the original algorithm `higherr` mode, `CJunction.nreads` and `CJunction.nreads_good`
        // can temporarily store negative demotion pointers (indices) rather than counts.
        // Those pointer values must never be aggregated back into numeric support fields.
        entry.mrcount += cj.nreads.max(0.0);
        entry.nreads_good += cj.nreads_good.max(0.0);
        entry.leftsupport += cj.leftsupport;
        entry.rightsupport += cj.rightsupport;
        entry.nm += cj.nm;
        entry.mm += cj.mm;
        entry.rcount = entry.rcount.saturating_add(cj.rcount);
        entry.guide_match |= cj.guide_match;
        if entry.strand.is_none() || entry.strand == Some(0) {
            entry.strand = Some(cj.strand);
        }
        // Preserve non-canonical evidence when any source marks a site as non-canonical.
        if entry.consleft != 0 {
            entry.consleft = if cj.consleft == 0 {
                0
            } else if entry.consleft == -1 {
                cj.consleft
            } else {
                entry.consleft
            };
        }
        if entry.consright != 0 {
            entry.consright = if cj.consright == 0 {
                0
            } else if entry.consright == -1 {
                cj.consright
            } else {
                entry.consright
            };
        }
    }
    out
}

/// Per-junction statistics (mrcount, etc.). CJunction: nreads, nreads_good, strand (0=killed), leftsupport, rightsupport.
#[derive(Debug, Clone)]
pub struct JunctionStat {
    pub mrcount: f64,
    /// Reads with both donor+acceptor anchors passing threshold (nreads_good).
    pub nreads_good: f64,
    pub rcount: u32,
    /// Strand: None = unknown, Some(0) = killed (good_junc sets strand=0), Some(1) = +, Some(-1) = -.
    pub strand: Option<i8>,
    /// Mismatch-weighted count: reads with high mismatch rate crossing this junction (nm; original algorithm jd.nm). Used in junction correction.
    pub nm: f64,
    /// Well-anchored count: reads where BOTH flanking exons exceed LONGINTRONANCHOR (25bp).
    /// jd.mm; used by dynamic anchor escalation (D3: not all-bad if mm>0).
    pub mm: f64,
    /// Support on donor (left) side of splice (leftsupport); used by good_junc / witness left.
    pub leftsupport: f64,
    /// Support on acceptor (right) side of splice (rightsupport); used by good_junc / witness right.
    pub rightsupport: f64,
    /// True if this junction matches a guide annotation intron. jd.guide_match.
    /// Used by eonly mode: kills all non-guide junctions (good_junc:13972).
    pub guide_match: bool,
    /// Canonical left (donor) dinucleotide check (jd.consleft).
    /// -1 = unknown (no genome available, never fires kill); 0 = non-canonical; 1 = canonical (GT/GC for +, CT for -).
    /// kill if !guide_match && consleft==0 && nreads_good < DROP/ERROR_PERC.
    pub consleft: i8,
    /// Canonical right (acceptor) dinucleotide check (jd.consright).
    /// -1 = unknown; 0 = non-canonical; 1 = canonical (AG for +, AC/GC for -).
    /// kill if !guide_match && consright==0 && nreads_good < DROP/ERROR_PERC.
    pub consright: i8,
}

impl Default for JunctionStat {
    fn default() -> Self {
        Self {
            mrcount: 0.0,
            nreads_good: 0.0,
            rcount: 0,
            strand: None,
            nm: 0.0,
            mm: 0.0,
            leftsupport: 0.0,
            rightsupport: 0.0,
            guide_match: false,
            consleft: -1, // -1 = unknown; mimics char leftcons=-1 (no genome → !(-1)=false → never kills)
            consright: -1,
        }
    }
}

pub type JunctionStats = DetHashMap<Junction, JunctionStat>;

/// Bundlenode: contiguous segment from merged read coverage (CBundlenode).
#[derive(Debug, Clone)]
pub struct CBundlenode {
    pub start: u64,
    pub end: u64,
    pub cov: f64,
    pub bid: usize,
    pub next: Option<Box<CBundlenode>>,
    pub hardstart: bool,
    pub hardend: bool,
}

impl CBundlenode {
    /// Deep clone the linked list.
    pub fn deep_clone(&self) -> CBundlenode {
        CBundlenode {
            start: self.start,
            end: self.end,
            cov: self.cov,
            bid: self.bid,
            next: self.next.as_ref().map(|n| Box::new(n.deep_clone())),
            hardstart: self.hardstart,
            hardend: self.hardend,
        }
    }
}

/// Maps bundlenode id -> list of (graph_id, node_id) for graph nodes created from that bnode (CGraphinfo / bundle2graph).
pub type Bundle2Graph = Vec<Vec<(usize, usize)>>;

impl CBundlenode {
    pub fn new(start: u64, end: u64, cov: f64, bid: usize) -> Self {
        Self {
            start,
            end,
            cov,
            bid,
            next: None,
            hardstart: false,
            hardend: false,
        }
    }
}

/// Bundle aggregate compatible with `BundleData` (subset of high-impact fields for Rust pipeline handoff).
#[derive(Debug, Clone)]
pub struct BundleData {
    pub idx: usize,
    pub start: u64,
    pub end: u64,
    pub refseq: String,
    pub strand: char,
    pub readlist: Vec<BundleRead>,
    pub junction: Vec<CJunction>,
    pub pred: Vec<CPrediction>,
    pub numreads: u64,
    pub num_fragments: f64,
    pub sum_cov: f64,
    /// Guide identifiers carried with this bundle.
    pub keepguides: Vec<String>,
    /// Long-read CPAS cut positions by strand index: [0]=minus, [1]=plus.
    pub cpas_cuts: [Vec<u64>; 2],
    /// Adjacent-junction long-read witness counts (packed pair key -> support).
    pub lr_pair: DetHashMap<u64, u32>,
}

impl BundleData {
    pub fn from_bundle(bundle: &Bundle) -> Self {
        let num_fragments = bundle
            .reads
            .iter()
            .map(|r| 1.0 / (r.nh.max(1) as f64))
            .sum::<f64>();
        Self {
            idx: 0,
            start: bundle.start,
            end: bundle.end,
            refseq: bundle.chrom.clone(),
            strand: bundle.strand,
            readlist: bundle.reads.clone(),
            junction: junction_stats_to_cjunctions(&bundle.junction_stats),
            pred: Vec::new(),
            numreads: bundle.reads.len() as u64,
            num_fragments,
            sum_cov: 0.0,
            keepguides: Vec::new(),
            cpas_cuts: [Vec::new(), Vec::new()],
            lr_pair: Default::default(),
        }
    }

    pub fn into_bundle(self) -> Bundle {
        Bundle {
            chrom: self.refseq,
            start: self.start,
            end: self.end,
            strand: self.strand,
            reads: self.readlist,
            junction_stats: cjunctions_to_junction_stats(&self.junction),
            bundlenodes: None,
            read_bnodes: None,
            bnode_colors: None,
            synthetic: false,
            rescue_class: None,
        }
    }
}

/// Multi-mapping resolution method for VG mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VgSolver {
    /// No multi-mapping resolution — discovery/reporting only (default).
    None,
    /// Expectation-Maximization with junction-based compatibility (heuristic).
    Em,
    /// Expectation-Maximization with HMM-based per-path sequence likelihood
    /// scoring. Requires read sequences for multi-mappers (collected at BAM
    /// parse) and a fitted family-graph HMM per family group.
    EmHmm,
    /// Per-family dispatch: HMM for medium-divergence multi-copy families
    /// (2..=10 copies, has junctions, --genome-fasta provided); heuristic
    /// EM for the rest; skip noise families (n_copies > max, intronless).
    /// See loo_assembly cross-family results — HMM only pays off in the
    /// 30-90% pairwise-id band; cheaper for the high-similarity end and
    /// useless for intronless / extreme-divergence cases.
    Auto,
    /// Flow-based redistribution: two-pass assembly, redistribute proportional to transcript coverage.
    Flow,
}

impl Default for VgSolver {
    fn default() -> Self {
        Self::None
    }
}

impl std::str::FromStr for VgSolver {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "none" | "off" | "discover" => Ok(Self::None),
            "em" => Ok(Self::Em),
            "em-hmm" | "em_hmm" | "emhmm" | "hmm" => Ok(Self::EmHmm),
            "auto" => Ok(Self::Auto),
            "flow" => Ok(Self::Flow),
            _ => Err(format!("unknown VG solver '{}': expected none, em, em-hmm, auto, or flow", s)),
        }
    }
}

impl std::fmt::Display for VgSolver {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::None => write!(f, "none"),
            Self::Em => write!(f, "em"),
            Self::EmHmm => write!(f, "em-hmm"),
            Self::Auto => write!(f, "auto"),
            Self::Flow => write!(f, "flow"),
        }
    }
}

/// Global run configuration.
#[derive(Debug, Clone)]
pub struct RunConfig {
    /// Long-read mode (-L). When true, reads contribute to abundance and longstart/longend; when false, short-read mode (srabund, coverage trim).
    pub long_reads: bool,
    /// With -L: mixed mode if >0. Reads with query_length >= this are long, others short (mixedMode; original algorithm --long-read-min-len).
    pub long_read_min_len: u64,
    pub min_junction_reads: f64,
    pub junction_support: u64,
    /// Window (bp) for junction correction clustering (donor/acceptor within window).
    pub junction_correction_window: u64,
    /// Tolerance for demotion: demote when best_reads * tolerance > reads_j (e.g. 0.9).
    pub junction_correction_tolerance: f64,
    /// Tolerance (bp) for merging junction boundary variants after filtering (e.g. 5).
    pub junction_canonical_tolerance: u64,
    /// Per-splice-site isofrac (percent): remove junctions with < this % of donor total (0 = off).
    pub per_splice_site_isofrac: f64,
    /// Minimum intron length (bp) to count as splice junction; smaller gaps are deletions.
    pub min_intron_length: u64,
    /// Merge distance (bp) for building bundlenodes from exons (0 for long-read).
    pub bundle_merge_dist: u64,
    /// Gap (bp) to start a new bundle; original algorithm.cpp:598 — new bundle when read start > current_end + runoffdist. Default 200.
    pub bundle_runoff_dist: u64,
    /// Junction support for CGroup small-exon keep rule (e.g. 10).
    pub cgroup_junction_support: u64,
    pub readthr: f64,
    pub singlethr: f64,
    pub min_transcript_length: u64,
    /// Minimum exon length (bp); drop transcript if any exon shorter (0 = no filter).
    pub min_exon_length: u64,
    pub label: String,
    pub verbose: bool,
    /// Estimate-only mode: only output transcripts matching guide annotations (-e flag).
    /// Kills non-guide junctions, restricts keeptrf to guide transfrags, filters non-guide output.
    pub eonly: bool,
    /// Debug bundle ID (e.g., "chr1:1000-2000") to trace flow computation.
    pub debug_bundle: Option<String>,
    /// Debug: dump junction sets for matching bundle(s).
    pub debug_junctions: bool,
    /// Restrict the run to bundles overlapping `debug_bundle`.
    pub only_debug_bundle: bool,
    /// Debug: trace a specific intron chain during extraction (donor-acceptor pairs).
    pub trace_intron_chain: Option<Vec<(u64, u64)>>,
    /// Emit contained isoforms by trimming transcript ends (may increase output).
    pub emit_contained_isoforms: bool,
    /// Emit additional transcripts by enumerating splice graph paths (max sensitivity; not default).
    pub emit_junction_paths: bool,
    /// Emit terminal alt-acceptor rescue transcripts (LR-only; the original algorithm does not do this).
    pub emit_terminal_alt_acceptor: bool,
    /// Emit micro-exon insertion rescue transcripts (LR-only; the original algorithm does not do this).
    pub emit_micro_exon_rescue: bool,
    /// Max number of junction-path transcripts to emit per bundle.
    pub max_junction_paths: usize,
    /// Max nodes in a junction-path (including source/sink).
    pub max_junction_path_len: usize,
    pub threads: usize,
    /// Remove transcripts that are fully contained in another with >= coverage.
    pub filter_contained: bool,
    /// Isofrac: remove transcript if coverage < this fraction of locus max (0 = off).
    pub transcript_isofrac: f64,
    /// Min read-abundance (max of longcov/coverage) to retain through long-read interval isofrac
    /// even when below the dominant isoform fraction (0 = off). Default 1.0 ≈ one supporting read.
    pub transcript_isofrac_keep_min: f64,
    /// Low isofrac (header lowisofrac=0.02): keep transcript if cov >= this * neighbor (pairwise filter). Lower = more permissive.
    pub lowisofrac: f64,
    /// Merge consecutive exons when gap <= this (bp); 0 = off.
    pub merge_micro_intron_max_gap: u64,
    /// Compute and output TPM in GTF attributes.
    pub output_tpm: bool,
    /// Rawreads mode (ref): emit predictions directly from kept transfrags.
    pub rawreads: bool,
    /// If set, write Ballgown-style t_data.ctab and e_data.ctab into this directory.
    pub ballgown_dir: Option<String>,
    /// Optional gene-abundance TSV output path.
    pub gene_abundance_out: Option<String>,
    /// Emit a "gene" feature line per transcript in GTF (gene_id/transcript_id scope).
    pub gtf_gene_lines: bool,
    // --- Advanced post-assembly filters (reference script) ---
    /// Min length (bp) for first/last exon; drop multi-exon if either terminal < this (0 = off).
    pub min_terminal_exon_len: u64,
    /// Relax post-extraction filters to keep more intron chains (skip intronic, longread-interval, incomplete-vs-guide).
    pub relax_filters: bool,
    /// Enable pairwise overlap filter (junction inclusion + containment).
    pub pairwise_overlap_filter: bool,
    /// Pairwise: drop single-exon contained if cov < this fraction of container (isofrac; lower = keep more).
    pub pairwise_isofrac: f64,
    pub pairwise_error_perc: f64,
    pub pairwise_drop: f64,
    /// Enable retained-intron filter.
    pub retained_intron_filter: bool,
    pub retained_intron_error_perc: f64,
    /// Enable polymerase run-off filter (single-exon fragments near multi-exon ends).
    pub polymerase_runoff_filter: bool,
    /// Enable polymerase run-on filter (bridge split between genes).
    pub polymerase_runon_filter: bool,
    pub polymerase_runoff_dist: u64,
    pub polymerase_singlethr: f64,
    /// Dedup: remove transcripts whose intron chain is a proper subset of another (same strand).
    pub dedup_intron_chain_subsets: bool,
    /// Enable longtrim-like node splitting at read-boundary peaks (lstart/lend) in long/mixed modes.
    pub enable_longtrim: bool,
    /// Minimum accumulated boundary coverage required to split at a read boundary.
    pub longtrim_min_boundary_cov: f64,
    /// Enable splice-aware path extension during transcript extraction.
    pub enable_splice_path_extension: bool,
    /// Enable direct port trial of back_to_source_fast_long/fwd_to_sink_fast_long in strict long-read mode.
    /// When false, strict mode keeps current stable extension behavior.
    pub strict_longrec_port: bool,
    /// Nascent transcription mode (isnascent flag): run a second parse_trflong pass
    /// for non-guide transfrags to build nascent transcript extensions
    pub isnascent: bool,
    /// Maximum active internal graph nodes before prune_graph_nodes is applied.
    pub allowed_graph_nodes: usize,
    /// When true, use junction-edge pruning (delete_connection style) instead of node removal.
    pub use_prune_by_junction: bool,
    /// Long intron threshold (header `longintron`).
    pub longintron: u64,
    /// Mismatch fraction threshold (header `mismatchfrac`).
    pub mismatchfrac: f64,
    /// Low coverage threshold (header `lowcov`).
    pub lowcov: f64,
    /// Splice-site error tolerance (`sserror`).
    pub sserror: u64,
    /// Optional genome FASTA for splice consensus checks.
    pub genome_fasta: Option<String>,
    /// Merge-mode BAM compatibility scaffold flag; currently warning-only.
    pub merge_bam_compat: bool,
    /// When true (default in merge mode), use graph-based merge (CMTransfrag, process_merge_transfrags, printMergeResults).
    pub use_graph_merge: bool,
    /// Disable coverage-based trimming (the original algorithm -t), including long-read terminal coverage trim on reads.
    pub no_coverage_trim: bool,
    /// Reference sequence IDs to exclude from assembly, comma-separated (the original algorithm -x).
    pub exclude_seqids: Vec<String>,
    /// Minimum anchor length for junctions in bp (the original algorithm -a).
    pub min_anchor_length: u64,
    /// Zero all isofrac/coverage filters for maximum sensitivity (diagnostic mode).
    pub max_sensitivity: bool,
    /// Optional TSV path for compact live debug stage summaries.
    pub debug_stage_tsv: Option<String>,
    /// Optional TSV path for keeptrf/usepath export from process_transfrags.
    pub keeptrf_usepath_tsv: Option<String>,
    /// Optional JSONL path: dump per-bundle snapshots (graph + transfrags) for 1:1 comparisons.
    pub snapshot_jsonl: Option<String>,
    /// True when `--trace-reference` was passed (reference fate / debug tracing).
    pub trace_reference: bool,
    /// Emit `ref_chain_junction_witness` assists when all ref introns are in good junctions but
    /// graph/guide extraction fails. Also enabled when `trace_reference` is true unless
    /// `RUSTLE_REF_JUNCTION_WITNESS=0` (see [`RunConfig::ref_junction_witness_enabled`]).
    pub ref_junction_witness: bool,
    // ── Variation graph (gene family) mode ───────────────────────────────────
    /// Enable variation graph mode for gene family assembly (--vg).
    pub vg_mode: bool,
    /// Minimum shared multi-mapping reads to link bundles into a family group [default: 3].
    pub vg_min_shared_reads: usize,
    /// Maximum EM iterations for multi-mapping resolution [default: 20].
    pub vg_em_max_iter: usize,
    /// Discover novel gene copies from poorly-mapped/unmapped reads (--vg-discover-novel).
    pub vg_discover_novel: bool,
    /// Output path for family group report TSV (--vg-report).
    pub vg_report: Option<std::path::PathBuf>,
    /// Multi-mapping solver: em or flow [default: em].
    pub vg_solver: VgSolver,
    /// Use SNPs for copy assignment (--vg-snp).
    pub vg_snp: bool,
    /// Phased assembly using HP tags (--vg-phase).
    pub vg_phase: bool,
    /// Minimum reads to create a novel copy bundle [default: 3].
    pub vg_min_novel_reads: usize,
    /// `kmer` (legacy) or `hmm` (new). Default `kmer` until validation flips it.
    pub vg_discover_novel_mode: String,
    /// Enable external minimap2 verification of rescued reads.
    pub vg_rescue_diagnostic: bool,
    /// Sequences for multi-mapped reads (read_name_hash → bytes), populated
    /// at BAM-parse time when `vg_solver == VgSolver::EmHmm`. Empty for other
    /// solvers (sequence collection has memory cost). Consumed by
    /// `run_pre_assembly_em_hmm` to compute per-paralog forward log-likelihoods.
    pub vg_multimap_sequences: std::collections::HashMap<u64, Vec<u8>>,
    /// Forward log-odds threshold for HMM rescue (nats).
    pub vg_rescue_min_loglik: f64,
    /// Mask regions for the novel-copy LOO experiment. Reads whose primary
    /// alignment overlaps any region in this list are treated as **unaligned**:
    /// they're excluded from bundle building (no normal transcript) and
    /// included in the HMM rescue pool (sequence used; alignment ignored).
    /// Useful for verifying `--vg-discover-novel` actually recovers a paralog
    /// from its sequence alone. Repeatable: `--vg-mask-region chrom:start-end`.
    pub vg_mask_regions: Vec<(String, u64, u64)>,
    /// Skip families with > N copies during EM (`--vg-em-max-copies`).
    /// Default 20 — bigger groups are usually noise (mtDNA, repetitive
    /// elements, very large gene-family clusters where pairwise-ID is
    /// either ~100% or below useful threshold).
    pub vg_em_max_copies: usize,
    /// Cap for routing a family to HMM-EM under `--vg-solver auto`.
    /// Default 10 — HMM scoring scales as O(reads × copies × seq × profile);
    /// medium families fit, megafamilies don't.
    pub vg_em_hmm_max_copies: usize,
    /// Skip intronless families during EM (`--vg-em-skip-intronless`).
    /// Default true — intronless paralogs (e.g. olfactory receptors) yield
    /// degenerate single-node family graphs that are uninformative for
    /// HMM scoring. Per the loo_assembly cross-family test, OR cluster:
    /// 0/0 reads rescued.
    pub vg_em_skip_intronless: bool,
    /// Family-quality filter: minimum total multi-mapping reads per family
    /// (post-discovery). Below this, a "family" is more likely random
    /// alignment artifacts than a real multi-copy paralog cluster.
    pub vg_family_min_shared: usize,
    /// Family-quality filter: maximum copies per family. Above this, a
    /// "family" is more likely a repetitive element / mega-cluster than
    /// a real paralog group.
    pub vg_family_max_copies: usize,
    /// Family-quality filter: minimum shared multi-mapping reads per copy
    /// (multimap_reads / n_copies). Sparse mega-clusters (one cross-mapper
    /// linking many bundles) fail this; real paralogs maintain ≥0.5 reads/copy.
    pub vg_family_min_shared_per_copy: f64,
    /// Family-quality filter: maximum coefficient of variation (CV) of
    /// per-copy intron counts. Real paralogs of the same gene have similar
    /// exon structure (CV typically <0.3); mixed-gene clusters that happen
    /// to share a multi-mapper have wildly different exon counts (CV often
    /// >1.0). Single-exon paralog clusters skip this check.
    pub vg_family_max_exon_cv: f64,
}

impl RunConfig {
    /// Whether junction-witness reference assists are allowed for this run.
    ///
    /// Enabled when `--ref-junction-witness` or `--trace-reference` is set, unless
    /// `RUSTLE_REF_JUNCTION_WITNESS` is set to a falsey value (`0`, `false`, `no`, `off`).
    /// If `RUSTLE_REF_JUNCTION_WITNESS` is set to a truthy value (`1`, `true`, `yes`, `on`), enables
    /// witness assists even when both flags are unset.
    pub fn ref_junction_witness_enabled(&self) -> bool {
        match std::env::var("RUSTLE_REF_JUNCTION_WITNESS") {
            Ok(v) => {
                let s = v.trim();
                if s.eq_ignore_ascii_case("0")
                    || s.eq_ignore_ascii_case("false")
                    || s.eq_ignore_ascii_case("no")
                    || s.eq_ignore_ascii_case("off")
                {
                    return false;
                }
                if s.eq_ignore_ascii_case("1")
                    || s.eq_ignore_ascii_case("true")
                    || s.eq_ignore_ascii_case("yes")
                    || s.eq_ignore_ascii_case("on")
                {
                    return true;
                }
            }
            Err(_) => {}
        }
        self.ref_junction_witness || self.trace_reference
    }

    /// Zero all isofrac/coverage thresholds for maximum sensitivity (diagnostic mode).
    pub fn apply_max_sensitivity(&mut self) {
        self.transcript_isofrac = 0.0;
        self.lowisofrac = 0.0;
        self.pairwise_isofrac = 0.0;
        self.readthr = 0.0;
        self.singlethr = 0.0;
        self.per_splice_site_isofrac = 0.0;
        self.transcript_isofrac_keep_min = 0.0;
        // If we're explicitly asking for maximum sensitivity, also enable splice-graph
        // path enumeration so low-support alternative paths can still be emitted.
        self.emit_junction_paths = true;
        self.max_junction_paths = self.max_junction_paths.max(8000);
    }

    /// Long-read standard density: modest `-f` floors and a junction-path cap. Tuned to
    /// **recover sensitivity** vs an over-tight preset: redundancy is already reduced by
    /// `print_predcluster` (pairwise, longunder, dedup). We intentionally do **not** force
    /// `filter_contained` or per-splice graph pruning in the default path—those crushed gffcompare
    /// Sn on full BAM. For maximum precision, pass `--filter-contained` and/or
    /// `--per-splice-site-isofrac 0.02`, or set `RUSTLE_STRICT_PRESET=strict`.
    pub fn apply_compat_preset(&mut self) {
        let strict = std::env::var("RUSTLE_STRICT_PRESET")
            .map(|s| s.eq_ignore_ascii_case("strict"))
            .unwrap_or(false);

        if strict {
            const ISO_STRICT: f64 = 0.025;
            self.transcript_isofrac = self.transcript_isofrac.max(ISO_STRICT);
            self.pairwise_isofrac = self.pairwise_isofrac.max(ISO_STRICT);
            self.lowisofrac = self.lowisofrac.max(ISO_STRICT);
            if self.max_junction_paths > 400 {
                self.max_junction_paths = 400;
            }
            self.filter_contained = true;
            self.per_splice_site_isofrac = self.per_splice_site_isofrac.max(0.02);
            return;
        }

        // the original algorithm uses isofrac=0.01 (the default) for long-read mode (-L).
        // isofrac=0.05 is only for --conservative mode.
        // lowisofrac stays at 0.02 (header const, not changed for long reads).
        // The keep_min_abundance floor is set to 0 to match the original algorithm (no floor).
        const ISO: f64 = 0.01;
        const LOW_ISO: f64 = 0.02;
        self.transcript_isofrac = self.transcript_isofrac.max(ISO);
        self.pairwise_isofrac = self.pairwise_isofrac.max(ISO);
        self.lowisofrac = self.lowisofrac.max(LOW_ISO);
        self.transcript_isofrac_keep_min = 0.0;
        if self.max_junction_paths > 600 {
            self.max_junction_paths = 600;
        }
        // In long-read mode, the original algorithm only emits transcripts backed by complete read paths
        // (transfrags). It never enumerates novel junction combinations via DFS — that was
        // designed for short reads where no single read spans the full transcript. In LR mode,
        // junction-path DFS is the dominant FP source (~71% of overemission), generating
        // combinatorial A→B→C paths from reads that only individually span A→B and B→C.
        // Similarly, terminal_alt_acceptor_rescue and micro_exon_insertion_rescue are rustle
        // extensions that the original algorithm never performs; disabling them in LR debug mode removes
        // ~2,376 FP transcripts (1,868 + 508) at no sensitivity cost.
        if self.long_reads {
            self.emit_junction_paths = false;
            self.emit_terminal_alt_acceptor = false;
            self.emit_micro_exon_rescue = false;
            // Per-splice-site isofrac: remove junctions below 1% of their donor's total traffic.
            // Testing on GGO_19 (long reads, 1839-transcript ground truth) shows isofrac=1%
            // gains +2 TPs (+0.1pp sensitivity) while removing 48 j-FP transcripts (+0.9pp
            // precision), with zero net sensitivity cost. Higher thresholds (2-3%) hurt TPs.
            self.per_splice_site_isofrac = self.per_splice_site_isofrac.max(1.0);
            // LR aligners produce shifted splice site calls (same junction ±1-30bp).
            // Coalescing with 2bp tolerance leaves many near-identical junctions that
            // create tiny graph nodes (1-10bp) and fragment transfrags.  10bp tolerance
            // merges the shifted copies into the strongest representative.
            //
            // Skip under `RUSTLE_STRINGTIE_EXACT`: StringTie's cgroup / partition dump does not
            // apply this Rustle-specific merge before bundlenode construction.
            if !crate::stringtie_parity::stringtie_exact() {
                self.junction_canonical_tolerance = self.junction_canonical_tolerance.max(10);
            }
        }
        // filter_contained and retained_intron_filter are too aggressive for sensitivity;
        // they're available via CLI but not forced by the compatibility preset.
    }

    /// When `RUSTLE_STRINGTIE_EXACT=1`, undo Rustle-only tuning that diverges from StringTie's
    /// rlink inputs **before** bundle/cgroup construction (junction stats → read exons).
    ///
    /// Call **after** [`Self::apply_compat_preset`] so long-read defaults are applied first,
    /// then selectively tightened for bit-identical parity work (e.g. `partition_geometry` vs
    /// `PARITY_PARTITION_TSV`).
    pub fn apply_stringtie_exact_overrides(&mut self) {
        if !crate::stringtie_parity::stringtie_exact() {
            return;
        }
        // Compat LR preset raises `junction_canonical_tolerance` to ≥10bp to merge noisy
        // aligner splice shifts. StringTie's cgroup / CBundle bundlenode pass does not apply
        // that same merge before `post_bundle_partition`, so partition signatures drift.
        self.junction_canonical_tolerance = 0;
    }

    /// Assembly mode: always long-read.
    pub fn assembly_mode(&self) -> AssemblyMode {
        AssemblyMode::LongRead
    }
}

impl Default for RunConfig {
    fn default() -> Self {
        Self {
            long_reads: true,
            long_read_min_len: 0,
            min_junction_reads: 1.0, // the original algorithm default junctionthr=1 (original algorithm.cpp:142)
            junction_support: 10,
            junction_correction_window: 0,
            junction_correction_tolerance: 0.0,
            junction_canonical_tolerance: 2, // Enable 2bp tolerance for canonicalization
            per_splice_site_isofrac: 0.05, // Enable 5% per-splice-site isofrac
            min_intron_length: 20,
            bundle_merge_dist: 0,  // Long reads use 0; sets bundledist=0 for longreads
            bundle_runoff_dist: 200,
            cgroup_junction_support: 10,
            readthr: 1.0, // the original algorithm/reference -L: READTHR_LONGREAD=1.0 (original algorithm.cpp:143)
            singlethr: 4.75, // the original algorithm/reference -L: SINGLETHR_LONGREAD=4.75
            min_transcript_length: 200,
            min_exon_length: 0,
            label: "RUSTLE".to_string(),
            verbose: false,
            eonly: false,
            debug_bundle: None,
            debug_junctions: false,
            only_debug_bundle: false,
            trace_intron_chain: None,
            emit_contained_isoforms: false,
            emit_junction_paths: true,
            emit_terminal_alt_acceptor: true,
            emit_micro_exon_rescue: true,
            max_junction_paths: 600, // Balanced: sensitivity needs junction paths
            max_junction_path_len: 200,
            threads: 0, // 0 = auto-detect available CPUs
            filter_contained: false,
            transcript_isofrac: 0.01, // the original algorithm default -f 0.01
            transcript_isofrac_keep_min: 1.0,
            lowisofrac: 0.02, // the original algorithm default lowisofrac=0.02
            merge_micro_intron_max_gap: 0,
            output_tpm: false,
            rawreads: false,
            ballgown_dir: None,
            gene_abundance_out: None,
            gtf_gene_lines: false,
            min_terminal_exon_len: 0,
            relax_filters: false,
            pairwise_overlap_filter: true,
            pairwise_isofrac: 0.01, // the original algorithm default isofrac=0.01
            pairwise_error_perc: 0.1,
            pairwise_drop: 0.5,
            retained_intron_filter: false,
            retained_intron_error_perc: 0.1,
            polymerase_runoff_filter: false,
            polymerase_runon_filter: false,
            polymerase_runoff_dist: 0,
            polymerase_singlethr: 4.75,
            dedup_intron_chain_subsets: true, // Enabled by default for better precision
            enable_longtrim: true,
            longtrim_min_boundary_cov: 2.0,
            enable_splice_path_extension: false,
            strict_longrec_port: true,
            isnascent: false,
            allowed_graph_nodes: 1_000_000,
            use_prune_by_junction: false,
            longintron: LONGINTRON,
            mismatchfrac: MISMATCHFRAC,
            lowcov: LOWCOV,
            sserror: SSERROR,
            genome_fasta: None,
            merge_bam_compat: false,
            use_graph_merge: true,
            no_coverage_trim: false,
            exclude_seqids: Vec::new(),
            min_anchor_length: 10,
            max_sensitivity: false,
            debug_stage_tsv: None,
            keeptrf_usepath_tsv: None,
            snapshot_jsonl: None,
            trace_reference: false,
            ref_junction_witness: false,
            vg_mode: false,
            vg_min_shared_reads: 3,
            vg_em_max_iter: 20,
            vg_discover_novel: false,
            vg_report: None,
            vg_solver: VgSolver::Em,
            vg_snp: false,
            vg_phase: false,
            vg_min_novel_reads: 3,
            vg_discover_novel_mode: "kmer".to_string(),
            vg_rescue_diagnostic: false,
            vg_rescue_min_loglik: 30.0,
            vg_multimap_sequences: std::collections::HashMap::new(),
            vg_mask_regions: Vec::new(),
            vg_em_max_copies: 20,
            vg_em_hmm_max_copies: 10,
            vg_em_skip_intronless: true,
            vg_family_min_shared: 10,
            vg_family_max_copies: 30,
            vg_family_min_shared_per_copy: 1.0,
            vg_family_max_exon_cv: 1.5,
        }
    }
}
