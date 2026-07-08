//! De-novo family detection DRIVER (integration stage 3): chains the ported cores into a family roster.
//!
//! **Terminology.** "Locus" in this module is the **gene-locus** sense unless
//! explicitly noted: a set of isoforms collapsed by shared splice junctions
//! (`family_detect::collapse_loci`). The physical `(chrom, start, end)` span
//! used for the ≥2-distinct-loci certificate is `family_definition::distinct_loci`.
//! See `docs/VG_FAMILY_TERMS.md` for the canonical vocabulary.
//!
//!   primary reads ─► pass1 skeletons ─► assemble gate ─► collapse loci ─► detect edges ─► decompose families
//!
//! This is the read-coherence-way detection pipeline end to end (rescue + per-read copy assignment are the
//! next stage). `detect_families` is the testable transform over parsed reads + a loaded genome;
//! `detect_families_from_files` is the thin I/O wrapper (BAM + FASTA).

use std::collections::{BTreeSet, HashSet};

use anyhow::Result;

use super::absent_copy::{self, AbsentCopyParams, Admission, DnaNeedsRecord};
use super::copy_assign::{Assignment, AssignParams, AssignStatus};
use super::copy_assign_pipeline::{assign_family_detailed, assign_family_detailed_pruned, freeze_merge, read_ref_end};
use super::family_graph::contiguous_core_coverage_bounded;
use super::copy_split::{
    split_locus_copies, discover_locus_psvs, AlignedRead, CollapsedCandidate, CopyIsoform,
};
use super::denovo_assemble::{
    assemble_gate, pass1_skeletons, pass1_skeletons_robust, primary_reads_from_bam, reads_in_region,
    BamRead, GateParams, PrimaryRead, PASS1_MIN_READS,
};
use super::family_detect::{collapse_loci_span_aware, detect_edges, detect_edges_reporting, DenovoTranscript, DetectParams};
use super::family_rescue::{FamilyMember, RescueParams};
use super::family_split::{classify, community_stats, decompose_families, FamilyClass, SplitFamily, SplitParams};
use super::read_conflict::{as_tie_edges, conflict_edges, conflict_families, family_mapq0_support, ConflictParams, Placement, ReadPlacements};
use super::rescue_pipeline::{rescue_thin_loci_iterative, thin_loci, MemberSpan, RESCUE_MIN_SUPPORT};
use crate::genome::GenomeIndex;

/// Configuration for the de-novo detection pipeline (defaults mirror the python stages).
#[derive(Clone, Copy, Debug)]
pub struct DenovoConfig {
    pub pass1_min_reads: u32,
    /// A skeleton's transcript extent must be reached by at least this many reads, so a single runaway read
    /// (chimeric / intra-primed / mis-clipped terminal exon) cannot inflate the length. `1` = the old union
    /// min-start / max-end behaviour; `2` (default) trims a lone outlier at each end.
    pub min_terminal_support: u32,
    pub gate: GateParams,
    pub detect: DetectParams,
    pub split: SplitParams,
    pub conflict: ConflictParams,
    /// POA-CORE COMPLETION (default false = OFF, byte-identical): after building the read-conflict families,
    /// attach loosely-related paralogs at NEW loci whose contiguous POA core to a family member clears
    /// `detect.t_core` — reaching divergent copies the conflict graph (confusable-only) misses. Bounded +
    /// seeded by the conflict families (`family_detect::poa_core_completion_adds`).
    pub complete_poa_core: bool,
}

impl Default for DenovoConfig {
    fn default() -> Self {
        DenovoConfig {
            pass1_min_reads: PASS1_MIN_READS,
            min_terminal_support: 2,
            gate: GateParams::default(),
            detect: DetectParams::default(),
            split: SplitParams::default(),
            conflict: ConflictParams::from_env(),
            complete_poa_core: false,
        }
    }
}

/// One decomposed de-novo family: its member transcript ids + structural diagnostics + class.
#[derive(Clone, Debug)]
pub struct DenovoFamily {
    pub members: Vec<String>,
    pub n_chroms: usize,
    pub density: f64,
    pub avg_core_recip: f64,
    pub n_articulation: usize,
    pub class: FamilyClass,
}

/// The detection pipeline's result: the general-assembly + family-layer counts and the family roster.
#[derive(Clone, Debug)]
pub struct DenovoResult {
    pub n_skeletons: usize,
    pub n_transcripts: usize,
    pub n_loci: usize,
    pub n_edges: usize,
    pub families: Vec<DenovoFamily>,
}

/// Run the de-novo family DETECTION pipeline from parsed primary reads + a loaded genome.
pub fn detect_families(reads: &[PrimaryRead], genome: &GenomeIndex, cfg: &DenovoConfig) -> DenovoResult {
    let skeletons = pass1_skeletons(reads, cfg.pass1_min_reads);
    let transcripts = assemble_gate(&skeletons, genome, &cfg.gate);
    // collapse isoforms -> gene loci (rep per locus), then build the rep transcript slice.
    // Span-aware collapse recovers alternative-splice isoforms with no shared junction.
    let rep_idx = collapse_loci_span_aware(&transcripts, &cfg.detect);
    let reps: Vec<DenovoTranscript> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
    // POA-confirmed homology edges over rep-slice indices, then decompose into families.
    let edges = detect_edges(&reps, &cfg.detect);
    let split = decompose_families(&edges, &cfg.split);
    let families = split
        .into_iter()
        .map(|sf| {
            let members: Vec<String> = sf.members.iter().map(|&i| reps[i].tid.clone()).collect();
            let n_chroms = sf
                .members
                .iter()
                .map(|&i| reps[i].chrom.as_str())
                .collect::<BTreeSet<_>>()
                .len();
            DenovoFamily {
                members,
                n_chroms,
                density: sf.stats.density,
                avg_core_recip: sf.stats.avg_core_recip,
                n_articulation: sf.stats.n_articulation,
                class: sf.class,
            }
        })
        .collect();
    DenovoResult {
        n_skeletons: skeletons.len(),
        n_transcripts: transcripts.len(),
        n_loci: reps.len(),
        n_edges: edges.len(),
        families,
    }
}

/// A co-located family: a same-chrom tandem cluster of copy transcripts ready for read assignment (the
/// hard regime where the aligner pools the copies). Mirrors `copy_assign.py::colocated_families`.
#[derive(Clone, Debug)]
pub struct ColocatedFamily {
    pub family_id: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub copies: Vec<DenovoTranscript>,
}

/// Filter decomposed families to same-chrom co-located clusters of `>= min_copies` copies spanning `<= win`.
/// Webs are skipped. `reps` are the rep transcripts the `SplitFamily` member indices point into.
/// `detect` configures the same-locus collapse in `prune_same_locus` (span-aware isoform recovery).
pub fn colocated_families(
    reps: &[DenovoTranscript],
    families: &[SplitFamily],
    win: u64,
    min_copies: usize,
    detect: &DetectParams,
) -> Vec<ColocatedFamily> {
    let mut out = Vec::new();
    for (fi, fam) in families.iter().enumerate() {
        if fam.class != FamilyClass::Family {
            continue;
        }
        // Partition by (chrom, STRAND): copies of a gene family share a strand. Opposite-strand loci that
        // merge via inverted-repeat / antisense homology are NOT copies of each other (a read carries its
        // own strand, so there is no copy-level ambiguity to resolve), and force-aligning a `+`-strand
        // spliced sequence to a `-`-strand one manufactures spurious PSV columns. Splitting by strand
        // removes that antisense over-merge from the copy-assignment input.
        let mut by_key: std::collections::BTreeMap<(&str, char), Vec<usize>> =
            std::collections::BTreeMap::new();
        for &m in &fam.members {
            by_key
                .entry((reps[m].chrom.as_str(), reps[m].strand))
                .or_default()
                .push(m);
        }
        for ((chrom, _strand), mut idxs) in by_key {
            idxs.sort_by_key(|&m| reps[m].start);
            let copies: Vec<DenovoTranscript> = idxs.iter().map(|&m| reps[m].clone()).collect();
            // SAME-LOCUS de-dup: copies of a multi-copy family must be at DISJOINT loci. Two transcripts
            // are the SAME locus (different isoforms/assemblies of one gene, NOT two copies) if they share
            // a junction or one is an unspliced span of the other. Collapsing these prevents force-aligning
            // a 1-exon read-through against its own 12-exon spliced form (the CAFAM0 over-merge). Real
            // tandem paralogs have DISJOINT junction sets at distinct coords, so they survive untouched.
            let copies = prune_same_locus(copies, detect);
            if copies.len() < min_copies {
                continue;
            }
            let span = copies.last().unwrap().start - copies[0].start;
            if span <= win {
                let start = copies.iter().map(|c| c.start).min().unwrap();
                let end = copies.iter().map(|c| c.end).max().unwrap();
                out.push(ColocatedFamily { family_id: format!("DSFAM{fi}"), chrom: chrom.to_string(), start, end, copies });
            }
        }
    }
    out
}

/// Collapse SAME-LOCUS transcripts so a multi-copy family's `copies` are at DISJOINT loci (the precondition
/// for PSV-based copy assignment — see `colocated_families`). Two same-strand transcripts are the same
/// locus iff:
///   (a) they **share a junction** — an intron `(donor, acceptor)` coordinate pair is an exact genomic
///       fingerprint of a splice site; two transcripts using it are isoforms of ONE gene, not two copies; or
///   (b) one is **structureless and contained** — a `<= 1`-exon transcript whose span is `>= CONTAIN_FRAC`
///       inside a more-structured one is the unspliced / read-through version of that copy; or
///   (d) they are **span-overlapping same-strand isoforms with high homology** — even with disjoint
///       junctions, two transcripts that overlap and have POA contiguous-core coverage >=
///       `p.collapse_span_core` are recovered as alternative splice isoforms of one gene. This mirrors
///       `family_detect::collapse_loci_span_aware`.
/// Tandem paralog copies have DISJOINT junction sets at distinct coordinates and typically fail the
/// conservative core-coverage bar, so they are preserved. Among same-locus transcripts the most-structured
/// representative is kept (more exons, then more reads, then longer sequence). Input may be in any order;
/// output preserves input order of the kept representatives.
fn prune_same_locus(copies: Vec<DenovoTranscript>, p: &DetectParams) -> Vec<DenovoTranscript> {
    const CONTAIN_FRAC: f64 = 0.5; // a structureless span this fraction inside another = the same locus
    let exon_count = |t: &DenovoTranscript| t.introns.len() + 1;
    // process best-structured first so it becomes the kept representative of its locus.
    let mut order: Vec<usize> = (0..copies.len()).collect();
    order.sort_by(|&a, &b| {
        exon_count(&copies[b])
            .cmp(&exon_count(&copies[a]))
            .then_with(|| copies[b].n_reads.cmp(&copies[a].n_reads))
            .then_with(|| copies[b].seq.len().cmp(&copies[a].seq.len()))
    });
    let mut kept_idx: Vec<usize> = Vec::new();
    for &i in &order {
        let same_locus = kept_idx.iter().any(|&k| {
            let (a, b) = (&copies[k], &copies[i]);
            // same chrom is a precondition for "same locus" (cross-chrom copies are always distinct loci).
            if a.chrom != b.chrom {
                return false;
            }
            // (a) shared junction → isoforms of one gene.
            let shares_junction = a.introns.iter().any(|j| b.introns.contains(j));
            if shares_junction {
                return true;
            }
            let span_overlap = b.end.min(a.end).saturating_sub(b.start.max(a.start));
            if span_overlap == 0 {
                return false;
            }
            let shorter = (a.end - a.start).min(b.end - b.start).max(1);
            let contained = span_overlap as f64 >= CONTAIN_FRAC * shorter as f64;
            // (b) structureless containment → the unspliced/read-through version of a copy.
            let structureless = exon_count(a) <= 1 || exon_count(b) <= 1;
            // (c) ANTISENSE at the SAME locus → overlapping spans on OPPOSITE strands = an inverted-repeat
            // artifact (a + gene and a − gene sharing sequence), NOT two paralog copies. Distinct-loci
            // opposite-strand copies (cross-chrom or distant same-chrom) do NOT overlap, so this never fires
            // on a genuine inverted-duplication paralog — which is exactly what we want to keep.
            let antisense_overlap = contained && a.strand != b.strand;
            if (contained && structureless) || antisense_overlap {
                return true;
            }
            // (d) span-aware homology: same-strand, span-overlapping, disjoint-junction isoforms with
            //     strong sequence homology. Avoid re-checking structureless pairs (already handled above).
            if p.collapse_span_aware
                && a.strand == b.strand
                && !structureless
                && !a.introns.is_empty()
                && !b.introns.is_empty()
            {
                let au = a.seq.to_ascii_uppercase();
                let bu = b.seq.to_ascii_uppercase();
                let core = contiguous_core_coverage_bounded(&au, &bu, p.len_cap);
                if core >= p.collapse_span_core {
                    return true;
                }
            }
            false
        });
        if !same_locus {
            kept_idx.push(i);
        }
    }
    // keep the chosen representatives in the input (start-sorted) order.
    let keep: std::collections::HashSet<usize> = kept_idx.into_iter().collect();
    copies
        .into_iter()
        .enumerate()
        .filter(|(i, _)| keep.contains(i))
        .map(|(_, c)| c)
        .collect()
}


/// Per-co-located-family read-assignment summary, with the two-pass (PSV-only vs PSV+junction) breakdown
/// and the silver-standard unique-mapper agreement (the accuracy proxy on real data — `copy_assign.py`).
#[derive(Clone, Debug)]
pub struct FamilyAssignment {
    pub family_id: String,
    pub chrom: String,
    pub n_copies: usize,
    pub n_reads: usize,
    pub psv_cols: usize,
    /// PSV-only.
    pub resolvable_psv: usize,
    pub assigned_psv: usize,
    /// PSV + copy-specific junctions.
    pub resolvable_j: usize,
    pub assigned_j: usize,
    /// reads a copy-specific junction resolved that PSVs alone could not.
    pub junction_only: usize,
    /// silver-standard: assigned + uniquely mapped (mapq > 0), and of those how many agree with where the
    /// read confidently mapped (the assignment-vs-mapping accuracy proxy).
    pub uniq: usize,
    pub uniq_agree: usize,
    /// EXTRA copies recovered by within-locus PSV split (`copy_split::split_locus_copies`): reads piled on
    /// one detected copy that actually carry `>= 2` PSV haplotypes — collapsed copies the aligner merged.
    /// (Caveat: het/editing/segdup can inflate this; treat as candidates.)
    pub collapsed_copies: usize,
    /// copies recovered by iterative family-aware RESCUE (below the assembly gate but homologous to the
    /// family); these are added to the copy set and reads are assigned to them. `n_copies` includes them.
    pub rescued_copies: usize,
    /// `(index into bam_reads, PSV+junction assignment)` for each read over the family.
    pub assignments: Vec<(usize, Assignment)>,
    /// SOFT per-copy quantification: copy transcript id, EM abundance (fraction of reads, sums to 1), and the
    /// 95% CI half-width — the probabilistic estimator (uses partial PSV evidence; uniform at the K=0 floor).
    pub copy_tids: Vec<String>,
    pub copy_abundance: Vec<f64>,
    pub copy_abundance_ci: Vec<f64>,
    /// gene-conversion mosaic: reads whose PSV path switches copy mid-molecule, and the family-confirmed
    /// conversion events (breakpoint recurs across molecules) — the enriched per-molecule multimapper signal.
    pub mosaic_reads: usize,
    pub conversions: Vec<crate::vg_family::mosaic::ConversionEvent>,
    /// copy-level historical gene conversions (a copy whose PSV-allele vector is a mosaic of two others).
    pub copy_conversions: Vec<crate::vg_family::copy_assign_pipeline::CopyConversion>,
    /// PSV genotype-matrix visualization data: column → forward-genome position, per-copy alleles, and per-read
    /// PSV observations (`read_psv_obs[i]` aligns with `assignments[i]`) — the raw per-molecule evidence that
    /// proves each copy assignment (a read carries one copy's alleles across its covered PSVs).
    pub psv_col_pos: Vec<Option<u64>>,
    pub copy_psv_alleles: Vec<Vec<Option<u8>>>,
    pub read_psv_obs: Vec<Vec<Option<u8>>>,
}

/// END-TO-END pipeline: detect families, then for each co-located family assign every read overlapping it to
/// a copy (PSV + copy-specific-junction likelihood, two-pass). `bam_reads` carry chrom (for region
/// filtering) and mapq (for the silver-standard). The runnable detection + per-read copy-assignment pipeline.
#[allow(clippy::too_many_arguments)]
/// A homology-prefiltered candidate family pair whose larger transcript exceeded the poasta memory threshold
/// (`cfg.detect.len_cap`) and was therefore confirmed via the memory-bounded longest-common-substring FALLBACK
/// rather than the exact poasta path. The edge is still produced (no OOM, no loss); this record surfaces which
/// family edges rest on the approximate large-sequence metric, for audit.
#[derive(Debug, Clone)]
pub struct FallbackEdge {
    pub chrom: String,
    pub tid_a: String,
    pub start_a: u64,
    pub end_a: u64,
    pub len_a: usize,
    pub tid_b: String,
    pub start_b: u64,
    pub end_b: u64,
    pub len_b: usize,
}

/// Map `bam_reads` (primary + secondary, grouped by read name) into per-read placement lists over `reps`.
/// Each BAM record is attributed to AT MOST ONE rep locus — the one with the greatest coordinate overlap
/// with the read.  This prevents uniquely-mapped reads that happen to span two coordinate-adjacent or
/// nested rep loci from producing spurious conflict edges: a genuine conflict requires a read to appear
/// TWICE in the BAM (once at each locus), not just to overlap two loci from a single alignment.
/// Supplementary (chimeric/split) records are excluded — they are not multimapping alternative placements.
/// The de-tie criterion uses `de` (gap-compressed divergence) to discriminate tied from resolvable placements.
pub(super) fn build_read_placements(bam_reads: &[BamRead], reps: &[DenovoTranscript]) -> Vec<ReadPlacements> {
    use std::collections::BTreeMap;
    let mut by_name: BTreeMap<&str, Vec<Placement>> = BTreeMap::new();
    for br in bam_reads {
        if br.is_supplementary {
            continue; // chimeric/split read is not a multimapping alternative placement
        }
        let read_end = read_ref_end(&br.read);
        let best = reps
            .iter()
            .enumerate()
            .filter(|(_, rep)| br.chrom == rep.chrom && br.read.ref_start < rep.end && read_end > rep.start)
            .max_by_key(|(_, rep)| rep.end.min(read_end) - rep.start.max(br.read.ref_start));
        if let Some((li, _)) = best {
            let aln_len = br.read.cigar.iter()
                .filter(|(op, _)| matches!(op, 'M' | '=' | 'X')).map(|(_, n)| *n).sum::<u64>() as u32;
            by_name.entry(br.name.as_str()).or_default().push(Placement {
                locus: li, de: br.de, mapq: br.mapq, as_score: br.as_score, aln_len,
            });
        }
    }
    by_name.into_values().collect()
}

/// Convert the read-conflict kernel's component output to `SplitFamily` objects that `colocated_families`
/// consumes. Edge weights (read counts) are cast to `f64` for `community_stats`; class follows the same
/// size+density rule as the POA path so large, sparse conflict components can still be flagged as Webs.
fn conflict_to_split_families(
    families: &[Vec<usize>],
    c_edges: &[(usize, usize, usize)],
    p: &SplitParams,
) -> Vec<SplitFamily> {
    let float_edges: Vec<(usize, usize, f64)> =
        c_edges.iter().map(|&(a, b, w)| (a, b, w as f64)).collect();
    let mut out: Vec<SplitFamily> = families
        .iter()
        .map(|members| {
            let mut m = members.clone();
            m.sort_unstable();
            let stats = community_stats(&m, &float_edges);
            let class = classify(stats.n, stats.density, p);
            SplitFamily { members: m, stats, class }
        })
        .collect();
    // deterministic: size desc, then smallest member.
    out.sort_by(|a, b| {
        b.members.len().cmp(&a.members.len()).then_with(|| a.members[0].cmp(&b.members[0]))
    });
    out
}

/// Emit every non-host identifiable collapsed copy at each locus as a [`CollapsedCandidate`].
///
/// For each representative transcript in `reps`, overlapping BAM reads are collected and fed to
/// [`split_locus_copies`]. When ≥2 identifiable copies are found the most-supported copy is treated
/// as the "host" (reference-anchored) copy and skipped; the rest become candidates for the
/// downstream admission gate. `psv_pos` is returned parallel to each `iso.allele_vector`.
///
/// The phantom safeguard is intrinsic: a fully-tied locus does not split (copy-specific PSVs are
/// required), so only copies with real distinguishing evidence enter the candidate list.
fn recover_collapsed_candidates(reps: &[DenovoTranscript], bam_reads: &[BamRead]) -> Vec<CollapsedCandidate> {
    let mut out: Vec<CollapsedCandidate> = Vec::new();
    for rep in reps {
        let reads: Vec<AlignedRead> = bam_reads
            .iter()
            .filter(|br| br.chrom == rep.chrom && br.read.ref_start < rep.end && read_ref_end(&br.read) > rep.start)
            .map(|br| br.read.clone())
            .collect();
        if reads.len() < 6 {
            continue;
        }
        let psv_pos = discover_locus_psvs(&reads, 3);
        let copies = split_locus_copies(&reads, 3, 2, 3);
        if copies.len() < 2 {
            continue;
        }
        let n_clusters = copies.iter().filter(|c| c.identifiable).count();
        // Emit every identifiable copy EXCEPT the most-supported one (treated as the host/ref copy).
        // Sort descending by read_count so skip(1) drops the host; ties broken by allele_vector for
        // determinism (BTreeMap-order allele_vector is a reliable secondary key).
        let mut ids: Vec<&CopyIsoform> = copies.iter().filter(|c| c.identifiable).collect();
        ids.sort_by(|a, b| {
            b.read_count
                .cmp(&a.read_count)
                .then_with(|| a.allele_vector.cmp(&b.allele_vector))
        });
        for iso in ids.into_iter().skip(1) {
            out.push(CollapsedCandidate {
                host_tid: rep.tid.clone(),
                chrom: rep.chrom.clone(),
                start: rep.start,
                end: rep.end,
                iso: iso.clone(),
                psv_pos: psv_pos.clone(),
                n_clusters,
            });
        }
    }
    out
}

/// Collapsed-copy recovery PAST the family gate. The genuinely collapsed tandem arrays (DAZ-type) don't form
/// a co-located family to rescue into, so this runs the PSV copy-split directly on EACH rep's overlapping
/// reads (`bam_reads` already include the secondary multimappers). Returns the number of EXTRA
/// PSV-DISTINGUISHABLE copies found beyond one-per-rep. The phantom safeguard is intrinsic: a fully-tied
/// (identical) locus does not split (`split_locus_copies` requires copy-specific PSVs), so only copies with
/// real distinguishing evidence are counted — the identifiability gate the user asked for, for free.
fn recover_collapsed_copies(reps: &[DenovoTranscript], bam_reads: &[BamRead]) -> usize {
    recover_collapsed_candidates(reps, bam_reads).len()
}

pub fn detect_and_assign(
    primary_reads: &[PrimaryRead],
    bam_reads: &[BamRead],
    genome: &GenomeIndex,
    cfg: &DenovoConfig,
    win: u64,
    min_copies: usize,
    p: &AssignParams,
    rescue_extra: &[PrimaryRead],
    absent_copies: bool,
    fasta_path: &str,
) -> (Vec<FamilyAssignment>, Vec<FallbackEdge>, Vec<DnaNeedsRecord>) {
    let timing = std::env::var_os("RUSTLE_TIMING").is_some();
    let mut t_lap = std::time::Instant::now();
    macro_rules! lap {
        ($name:expr) => {
            if timing {
                eprintln!("[timing]   {}: {:.1}s", $name, t_lap.elapsed().as_secs_f64());
                t_lap = std::time::Instant::now();
            }
        };
    }
    let skeletons = pass1_skeletons(primary_reads, cfg.pass1_min_reads);
    let transcripts = assemble_gate(&skeletons, genome, &cfg.gate);
    let rep_idx = collapse_loci_span_aware(&transcripts, &cfg.detect);
    let reps: Vec<DenovoTranscript> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
    eprintln!(
        "[detect_and_assign] {} primary -> {} skeletons -> {} transcripts -> {} reps",
        primary_reads.len(),
        skeletons.len(),
        transcripts.len(),
        reps.len()
    );
    lap!("pass1+gate+collapse");
    // Conflict-graph: AUTHORITATIVE family criterion
    let placements = build_read_placements(bam_reads, &reps);
    lap!(format!("build_read_placements ({} reads x {} reps)", bam_reads.len(), reps.len()));
    let c_edges = conflict_edges(reps.len(), &placements, &cfg.conflict);
    let c_fams = conflict_families(reps.len(), &c_edges);
    // audit: the de-tie edge set must be a subset of the AS-tie edge set (the bake-off invariant).
    let as_edges = as_tie_edges(reps.len(), &placements, 0.9, cfg.conflict.min_reads);
    let de_subset = c_edges.iter().all(|&(i, j, _)| as_edges.contains(&(i, j)));
    eprintln!(
        "[detect_and_assign] conflict-graph (de-tie): {} edges -> {} families | de⊆AS={} (AS edges={})",
        c_edges.len(), c_fams.len(), de_subset, as_edges.len(),
    );
    for (fi, fam) in c_fams.iter().enumerate() {
        let members: Vec<&str> = fam.iter().map(|&i| reps[i].tid.as_str()).collect();
        let coords: Vec<String> = fam.iter()
            .map(|&i| format!("{}:{}-{}", reps[i].chrom, reps[i].start, reps[i].end)).collect();
        let edge_weights: Vec<usize> = c_edges.iter()
            .filter(|&&(a, b, _)| fam.contains(&a) && fam.contains(&b)).map(|&(_, _, w)| w).collect();
        let (mq_support, mq_both0) = family_mapq0_support(&placements, fam, &cfg.conflict);
        eprintln!(
            "  conflict-fam{fi} n={} reads_linking={:?} mapq0_frac={}/{}: {} @ {}",
            fam.len(), edge_weights, mq_both0, mq_support, members.join(","), coords.join(" | "),
        );
    }
    let split = conflict_to_split_families(&c_fams, &c_edges, &cfg.split);

    // POA homology edges — DIAGNOSTIC ONLY, no longer drives family membership (families come from the
    // de-tie conflict graph at `conflict_to_split_families` above). The POA pairwise (poasta over all
    // candidate rep pairs) is the dominant cost on dense regions, so it is skippable with no effect on
    // the emitted families/assignments — only the diagnostic edge count and the `.fallback.tsv` report
    // are lost. Opt out with RUSTLE_SKIP_POA_DIAGNOSTIC=1 (e.g. genome-wide sweeps).
    let skip_poa = std::env::var_os("RUSTLE_SKIP_POA_DIAGNOSTIC").is_some();
    let (poa_edges, fallback_pairs) = if skip_poa {
        (Vec::new(), Vec::new())
    } else {
        detect_edges_reporting(&reps, &cfg.detect)
    };
    eprintln!(
        "[detect_and_assign] POA homology (diagnostic): {} edges ({} via large-seq fallback){}",
        poa_edges.len(), fallback_pairs.len(),
        if skip_poa { " [SKIPPED via RUSTLE_SKIP_POA_DIAGNOSTIC]" } else { "" },
    );
    lap!("conflict-graph + POA homology diagnostic");
    let fallback: Vec<FallbackEdge> = fallback_pairs
        .iter()
        .map(|&(a, b)| FallbackEdge {
            chrom: reps[a].chrom.clone(),
            tid_a: reps[a].tid.clone(), start_a: reps[a].start, end_a: reps[a].end, len_a: reps[a].seq.len(),
            tid_b: reps[b].tid.clone(), start_b: reps[b].start, end_b: reps[b].end, len_b: reps[b].seq.len(),
        })
        .collect();
    if !rescue_extra.is_empty() {
        let rec = recover_collapsed_copies(&reps, bam_reads);
        eprintln!(
            "[detect_and_assign] rescue_extra (AS-tied secondaries): {} | collapsed copies recovered past family gate (PSV-distinguishable): {}",
            rescue_extra.len(),
            rec
        );
    }
    let mut out = Vec::new();
    let mut dna_needs: Vec<DnaNeedsRecord> = Vec::new();
    for cf in colocated_families(&reps, &split, win, min_copies, &cfg.detect) {
        // RESCUE: recover under-assembled copies homologous to this family (below the >=3-read assembly gate)
        // and ADD them to the copy set, so reads can be assigned to them too. Iterative (bridge-aware).
        let members: Vec<FamilyMember> = cf
            .copies
            .iter()
            .map(|c| FamilyMember::new(c.tid.clone(), cf.family_id.clone(), c.seq.clone()))
            .collect();
        let member_spans: Vec<MemberSpan> = cf
            .copies
            .iter()
            .map(|c| MemberSpan { chrom: c.chrom.clone(), start: c.start, end: c.end })
            .collect();
        // scan a ±1 Mb NEIGHBOURHOOD around the family span (the python rescue `WIN`) so under-assembled
        // copies just outside the tight cluster are reachable — not only the detected-copy span.
        const RESCUE_WIN: u64 = 1_000_000;
        let (rlo, rhi) = (cf.start.saturating_sub(RESCUE_WIN), cf.end + RESCUE_WIN);
        // collapsed-copy recovery: AS-tied secondary reads (rescue_extra, empty by default) join the rescue
        // input so a copy whose reads minimap2 flagged secondary can clear the thin-loci support gate.
        let region_primary: Vec<PrimaryRead> = primary_reads
            .iter()
            .chain(rescue_extra.iter())
            .filter(|r| r.chrom == cf.chrom && r.ref_start < rhi && r.ref_end > rlo)
            .cloned()
            .collect();
        let loci = thin_loci(&region_primary, RESCUE_MIN_SUPPORT);
        let rescued =
            rescue_thin_loci_iterative(&loci, &members, &member_spans, genome, &RescueParams::default(), 3);
        let rescued_copies = rescued.len();
        let mut all_copies: Vec<DenovoTranscript> = cf.copies.clone();
        for rc in &rescued {
            all_copies.push(DenovoTranscript {
                tid: format!("RC_{}_{}", rc.locus.chrom, rc.locus.start),
                chrom: rc.locus.chrom.clone(),
                start: rc.locus.start,
                end: rc.locus.end,
                n_reads: rc.locus.support,
                strand: rc.strand,
                introns: rc.locus.introns.clone(),
                seq: rc.seq.clone(),
            });
        }
        // reads on this family's chrom overlapping its span (assign overlaps by coord, so pre-filter chrom).
        let mut idx_map = Vec::new();
        let mut region = Vec::new();
        let mut region_mapq = Vec::new();
        for (i, br) in bam_reads.iter().enumerate() {
            if br.chrom == cf.chrom && br.read.ref_start < cf.end && read_ref_end(&br.read) > cf.start {
                idx_map.push(i);
                region.push(br.read.clone());
                region_mapq.push(br.mapq);
            }
        }
        // Stage-1 and Stage-2 run WITHOUT iterative pruning so freeze_merge and downstream bookkeeping stay
        // in the original index space. Pruning, when requested, is applied as a final post-process below.
        let p_once = AssignParams { iterative_prune: false, ..*p };
        // Stage-1: assign over the reference copies only (borrow scoped so `all_copies` stays reassignable).
        let mut detail = {
            let copies: Vec<&DenovoTranscript> = all_copies.iter().collect();
            assign_family_detailed(&copies, &region, &p_once, Some(genome))
        };
        // Task 5 (opt-in): two-stage freeze for reference-ABSENT (collapsed) copies. OFF => this whole block
        // is skipped, so the loop below is byte-for-byte the pre-Task-5 path (`all_copies`/`detail` unchanged).
        if absent_copies {
            let cands = recover_collapsed_candidates(&all_copies, bam_reads);
            let mut admitted: Vec<DenovoTranscript> = Vec::new();
            for cand in &cands {
                if let Some(host) = all_copies.iter().find(|t| t.tid == cand.host_tid) {
                    match absent_copy::admit_candidate(cand, host, genome, fasta_path, &AbsentCopyParams::default()) {
                        Admission::Copy(t) => admitted.push(t),
                        // Task 6: collect DNA-needs records for the caller to surface as <out>.dna_needs.tsv.
                        Admission::DnaNeeds(r) => dna_needs.push(r),
                    }
                }
            }
            if !admitted.is_empty() {
                let n_ref = all_copies.len();
                let mut copies2_owned = all_copies.clone();
                copies2_owned.extend(admitted);
                // Stage-2: re-assign over ref (indices 0..n_ref) + admitted absent copies, then FREEZE: a read
                // Stage-1-Assigned at a multi-ref family keeps its Stage-1 result; the rest take Stage-2 (a
                // surviving absent-copy assignment is flagged `discovery_coupled`). Matched by read_index.
                {
                    let copies2: Vec<&DenovoTranscript> = copies2_owned.iter().collect();
                    let mut d2 = assign_family_detailed(&copies2, &region, &p_once, Some(genome));
                    d2.results = freeze_merge(&detail.results, std::mem::take(&mut d2.results), n_ref);
                    detail = d2;
                }
                // The ref copies keep indices 0..n_ref, so every per-read `best_copy` stays valid; the augmented
                // set now drives copy_tids / abundance / by_copy below.
                all_copies = copies2_owned;
            }
        }
        // IsoCon-style iterative copy pruning (opt-in): reassign all reads on the final copy set, dropping
        // copies that have no read with significant evidence against their nearest neighbor. This is a global
        // post-process so the output copy roster is internally consistent.
        if p.iterative_prune && all_copies.len() >= 2 {
            let copies: Vec<&DenovoTranscript> = all_copies.iter().collect();
            detail = assign_family_detailed_pruned(&copies, &region, p, Some(genome));
            all_copies = detail.copy_indices.iter().map(|&i| all_copies[i].clone()).collect();
        }
        // Unified gene-conversion-vs-artifact discriminator: tag each confirmed event by recurrence
        // (already in `confirmed`) + microhomology at the breakpoint (RT/template-switch signature,
        // from the genome) + DNA support (catalog leg = None here; the bench harness adds it).
        detail.conversion_class = super::copy_assign_pipeline::classify_conversions(
            &detail, genome, |_ev| None,
        );
        if detail.mosaic_reads > 0 || !detail.conversions.is_empty() {
            use super::mosaic::Classification as Cl;
            let confirmed = detail.conversions.iter().filter(|e| e.confirmed).count();
            let rt = detail.conversion_class.iter().filter(|&&c| c == Cl::RtSwitchArtifact).count();
            let gc = detail.conversion_class.iter().filter(|&&c| c == Cl::GeneConversion).count();
            eprintln!(
                "[mosaic] {} {}: {} mosaic-path reads, {} conversion events ({} confirmed -> {} gene_conversion, {} rt_switch_artifact)",
                cf.family_id, cf.chrom, detail.mosaic_reads, detail.conversions.len(), confirmed, gc, rt
            );
        }
        // collapsed-copy recovery: group the reads by their mapped copy/locus and split each by within-locus
        // PSV haplotype; >= 2 identifiable copies at one locus means extra (collapsed) copies were merged.
        let mut by_copy: Vec<Vec<AlignedRead>> = vec![Vec::new(); all_copies.len()];
        for r in &detail.results {
            by_copy[r.mapped_copy].push(region[r.read_index].clone());
        }
        let collapsed_copies: usize = by_copy
            .iter()
            .map(|reads| split_locus_copies(reads, 3, 2, 3).len().saturating_sub(1))
            .sum();
        let mut fa = FamilyAssignment {
            family_id: cf.family_id,
            chrom: cf.chrom,
            n_copies: all_copies.len(),
            n_reads: detail.results.len(),
            psv_cols: detail.n_cols,
            resolvable_psv: 0,
            assigned_psv: 0,
            resolvable_j: 0,
            assigned_j: 0,
            junction_only: 0,
            uniq: 0,
            uniq_agree: 0,
            collapsed_copies,
            rescued_copies,
            assignments: Vec::with_capacity(detail.results.len()),
            copy_tids: all_copies.iter().map(|c| c.tid.clone()).collect(),
            copy_abundance: detail.copy_abundance.clone(),
            copy_abundance_ci: detail.copy_abundance_ci.clone(),
            mosaic_reads: detail.mosaic_reads,
            conversions: detail.conversions.clone(),
            copy_conversions: detail.copy_conversions.clone(),
            psv_col_pos: detail.psv_col_pos.clone(),
            copy_psv_alleles: detail.copy_psv_alleles.clone(),
            read_psv_obs: Vec::with_capacity(detail.results.len()),
        };
        for r in detail.results {
            let resolvable_psv = r.psv.n_decisive >= 1;
            let assigned_j = r.combined.status == AssignStatus::Assigned;
            fa.resolvable_psv += resolvable_psv as usize;
            fa.assigned_psv += (r.psv.status == AssignStatus::Assigned) as usize;
            fa.resolvable_j += (r.combined.n_decisive >= 1) as usize;
            fa.assigned_j += assigned_j as usize;
            fa.junction_only += (r.combined.n_decisive >= 1 && !resolvable_psv) as usize;
            if assigned_j && region_mapq[r.read_index] > 0 {
                fa.uniq += 1;
                fa.uniq_agree += (r.combined.best_copy == r.mapped_copy) as usize;
            }
            fa.read_psv_obs.push(r.psv_obs);
            fa.assignments.push((idx_map[r.read_index], r.combined));
        }
        out.push(fa);
    }
    (out, fallback, dna_needs)
}

/// I/O wrapper: load primary reads from a BAM and the genome from a FASTA (scoped via `.fai` to only the
/// contigs the reads touch), then run detection.
pub fn detect_families_from_files(
    bam_path: &str,
    fasta_path: &str,
    threads: usize,
    cfg: &DenovoConfig,
) -> Result<DenovoResult> {
    let reads = primary_reads_from_bam(bam_path, threads)?;
    let contigs: HashSet<String> = reads.iter().map(|r| r.chrom.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;
    Ok(detect_families(&reads, &genome, cfg))
}

/// GENOME-WIDE de-tie READ-CONFLICT family catalog (interest I / O1) — the principled, threshold-free
/// family definition run at scale, replacing the per-region scan and the similarity-threshold POA catalog
/// (`detect_families` uses `detect_edges` = `core_recip≥0.13`, an arbitrary bar that over-merges; THIS uses
/// the read-conflict graph: a family is a connected component of loci among which reads are genuinely
/// confused, with NO similarity threshold). It applies the same-strand + disjoint-loci fixes
/// (`colocated_families`), so the output is the clean multi-copy-family catalog.
///
/// Architecture (memory-bounded): (1) build genome-wide reps once (skeletons → gate → collapse_loci),
/// then free the primaries/transcripts; (2) per CHROM, load that chromosome's reads (primary + secondary)
/// and build the conflict edges over the chrom's reps — same-chrom only, which is exactly what
/// `colocated_families` keeps anyway (cross-chrom components are split out there); (3) union the per-chrom
/// edges into the genome-wide conflict graph → connected components → `colocated_families`.
pub fn detect_conflict_catalog_genome_wide(
    bam_path: &str,
    fasta_path: &str,
    threads: usize,
    win: u64,
    min_copies: usize,
    cfg: &DenovoConfig,
) -> Result<Vec<ColocatedFamily>> {
    // --- (1) genome-wide reps ---
    let reads = primary_reads_from_bam(bam_path, threads)?;
    let contigs: HashSet<String> = reads.iter().map(|r| r.chrom.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;
    let skeletons = pass1_skeletons_robust(&reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    drop(reads); // free the ~1.7M primaries before the per-chrom read load
    let transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
    let rep_idx = collapse_loci_span_aware(&transcripts, &cfg.detect);
    let reps: Vec<DenovoTranscript> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
    drop(transcripts);
    eprintln!(
        "[gw-catalog] {} skeletons -> {} reps over {} contigs",
        skeletons.len(),
        reps.len(),
        contigs.len()
    );

    // global rep indices grouped by chrom (for per-chrom processing + local↔global remap).
    let mut by_chrom: std::collections::BTreeMap<&str, Vec<usize>> = std::collections::BTreeMap::new();
    for (gi, rep) in reps.iter().enumerate() {
        by_chrom.entry(rep.chrom.as_str()).or_default().push(gi);
    }

    // --- (2) per-chrom conflict edges (load each chrom's reads ONCE; same-chrom edges only) ---
    let mut all_edges: Vec<(usize, usize, usize)> = Vec::new();
    for (chrom, glob) in &by_chrom {
        if glob.len() < 2 {
            continue; // a single rep on a chrom can form no conflict edge
        }
        let clen = genome.chrom_len(chrom);
        let (_primary, bam_reads) = match reads_in_region(bam_path, chrom, 0, clen, threads) {
            Ok(r) => r,
            Err(_) => continue,
        };
        let chrom_reps: Vec<DenovoTranscript> = glob.iter().map(|&g| reps[g].clone()).collect();
        let placements = build_read_placements(&bam_reads, &chrom_reps);
        let edges = conflict_edges(chrom_reps.len(), &placements, &cfg.conflict);
        let n_edges = edges.len(); // this chrom's edge count (not the cumulative total)
        for (i, j, w) in edges {
            all_edges.push((glob[i], glob[j], w)); // local → global rep index
        }
        eprintln!(
            "[gw-catalog] {chrom}: {} reps, {} reads, {} conflict edges",
            chrom_reps.len(),
            bam_reads.len(),
            n_edges
        );
    }

    // --- (3) global components → split → colocated (strand-split + disjoint-loci de-dup) ---
    let c_fams = conflict_families(reps.len(), &all_edges);
    let split = conflict_to_split_families(&c_fams, &all_edges, &cfg.split);
    let catalog = colocated_families(&reps, &split, win, min_copies, &cfg.detect);
    eprintln!(
        "[gw-catalog] {} conflict components -> {} clean families (same-strand, disjoint-loci, >={} copies)",
        c_fams.len(),
        catalog.len(),
        min_copies
    );
    Ok(catalog)
}

/// O1 family emission from a de-tie read-conflict graph: Louvain-decompose each raw connected component
/// into DENSE communities (`decompose_families`), drop the sparse `Web` bridges (repeat-driven transitive
/// over-merges), then collapse SAME-LOCUS artifacts (`prune_same_locus`), keeping communities with
/// `>= min_copies` distinct copies.
///
/// This is chrom- AND strand-AGNOSTIC on purpose: a same-chromosome INVERTED duplication (opposite-strand,
/// disjoint loci), a distant intra-chromosomal segdup, and a cross-chromosome paralog pair are ALL
/// legitimate O1 families and all survive here. Contrast `colocated_families`, which additionally partitions
/// by `(chrom, strand)` and caps the span by `win` — those two filters are correct for building the O2
/// copy-ASSIGNMENT input (a single-strand tight tandem array) but WRONG for the O1 catalog, because they
/// erase inverted duplications and distant same-chrom paralogs. Routing same-chrom recovery through
/// `colocated_families` is exactly what dropped ~half of the same-chrom families; the O1 catalog uses this
/// helper instead.
fn families_from_conflict_graph(
    reps: &[DenovoTranscript],
    edges: &[(usize, usize, usize)],
    min_copies: usize,
    cfg: &DenovoConfig,
) -> Vec<Vec<DenovoTranscript>> {
    let float_edges: Vec<(usize, usize, f64)> =
        edges.iter().map(|&(a, b, w)| (a, b, w as f64)).collect();
    let split = crate::vg_family::family_split::decompose_families(&float_edges, &cfg.split);
    let mut out: Vec<Vec<DenovoTranscript>> = Vec::new();
    for sf in &split {
        if sf.class != FamilyClass::Family {
            continue; // drop sparse over-merge webs (repeat-driven transitive blobs)
        }
        let comp_reps: Vec<DenovoTranscript> = sf.members.iter().map(|&i| reps[i].clone()).collect();
        // same-locus artifacts (incl. same-locus antisense overlap); distinct-loci copies (incl. inverted
        // duplications and cross-chrom) are kept.
        let clean = prune_same_locus(comp_reps, &cfg.detect);
        if clean.len() >= min_copies {
            out.push(clean);
        }
    }
    out
}

/// CROSS-CHROMOSOME-aware genome-wide de-tie conflict family catalog. Unlike
/// `detect_conflict_catalog_genome_wide` (same-chrom only, via `colocated_families`'s `(chrom,strand)`
/// partition + `win`), this captures paralog families whose copies are on DIFFERENT chromosomes
/// (RABL2A/B, DAZ-class) or distant on the same chromosome. The de-tie conflict edge is chrom-agnostic
/// (a read confused between a chr12 copy and a chr6 copy de-ties just the same), so the only thing that
/// forbade cross-chrom families was the post-hoc partition — removed here.
///
/// Difference from the same-chrom path: (1) per-chrom read loads accumulate placements into ONE GLOBAL
/// `name → [Placement]` map with GLOBAL rep indices, so a read's alignments on different chromosomes land
/// in the same entry → a cross-chrom conflict edge forms; (2) each connected component is a family
/// directly (no `(chrom,strand)` partition, no `win`), cleaned only of SAME-LOCUS artifacts by
/// `prune_same_locus` (which now keeps distinct-loci opposite-strand copies = genuine inverted
/// duplications). Returns each family's copies (possibly spanning chromosomes).
pub fn detect_conflict_catalog_genome_wide_xchrom(
    bam_path: &str,
    fasta_path: &str,
    threads: usize,
    min_copies: usize,
    cfg: &DenovoConfig,
) -> Result<Vec<Vec<DenovoTranscript>>> {
    use super::copy_assign_pipeline::read_ref_end;
    use super::read_conflict::Placement;
    // --- genome-wide reps (same as the same-chrom path) ---
    let reads = primary_reads_from_bam(bam_path, threads)?;
    let contigs: HashSet<String> = reads.iter().map(|r| r.chrom.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;
    let skeletons = pass1_skeletons_robust(&reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    drop(reads);
    let transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
    let rep_idx = collapse_loci_span_aware(&transcripts, &cfg.detect);
    let reps: Vec<DenovoTranscript> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
    drop(transcripts);
    let mut by_chrom: std::collections::BTreeMap<&str, Vec<usize>> = std::collections::BTreeMap::new();
    for (gi, rep) in reps.iter().enumerate() {
        by_chrom.entry(rep.chrom.as_str()).or_default().push(gi);
    }
    eprintln!("[gw-catalog-xchrom] {} skeletons -> {} reps over {} contigs", skeletons.len(), reps.len(), contigs.len());

    // --- GLOBAL placement accumulation: a read's placements (across chroms) keyed by name, GLOBAL idx ---
    let mut name_map: std::collections::HashMap<String, Vec<Placement>> = std::collections::HashMap::new();
    for (chrom, glob) in &by_chrom {
        let clen = genome.chrom_len(chrom);
        let (_primary, bam_reads) = match reads_in_region(bam_path, chrom, 0, clen, threads) {
            Ok(r) => r,
            Err(_) => continue,
        };
        // chrom-local reps (small) for the efficient best-overlap search; remap local → global.
        let local: Vec<(usize, &DenovoTranscript)> = glob.iter().map(|&g| (g, &reps[g])).collect();
        for br in &bam_reads {
            if br.is_supplementary {
                continue;
            }
            let read_end = read_ref_end(&br.read);
            let best = local
                .iter()
                .filter(|(_, rep)| br.read.ref_start < rep.end && read_end > rep.start)
                .max_by_key(|(_, rep)| rep.end.min(read_end) - rep.start.max(br.read.ref_start));
            if let Some((gi, _)) = best {
                let aln_len = br.read.cigar.iter()
                    .filter(|(op, _)| matches!(op, 'M' | '=' | 'X')).map(|(_, n)| *n).sum::<u64>() as u32;
                name_map.entry(br.name.clone()).or_default().push(Placement {
                    locus: *gi,
                    de: br.de,
                    mapq: br.mapq,
                    as_score: br.as_score,
                    aln_len,
                });
            }
        }
    }
    let placements: Vec<ReadPlacements> =
        name_map.into_values().filter(|v| v.len() >= 2).collect();

    // --- global conflict graph (cross-chrom edges now present) ---
    let edges = conflict_edges(reps.len(), &placements, &cfg.conflict);
    // Louvain community split + Web/density classify on the conflict graph. WITHOUT the same-chrom
    // path's `(chrom,strand)` partition + `win` (which we removed to allow cross-chrom families), the raw
    // connected components transitively over-merge: a read confusing A↔B and another confusing B↔C chains
    // unrelated loci A,B,C (and repeat-driven cross-mapping produces 100s-of-member blobs spanning all
    // chromosomes). `decompose_families` breaks each component into DENSE communities (Louvain) and flags
    // the residual sparse bridges as `Web`; a real paralog family is a dense conflict subgraph (every copy
    // confused with the others), a transitive bridge-chain is sparse. We keep only `Family`-class
    // communities, then prune same-locus artifacts. This is the cross-chrom analogue of the same-chrom
    // path's over-merge guard.
    // O1 family emission (chrom/strand-agnostic Louvain + same-locus prune). Same-chromosome families —
    // including INVERTED duplications (opposite-strand, disjoint loci) and distant intra-chromosomal
    // segdups — are emitted here alongside cross-chromosome families. (The former `colocated_families`
    // supplement path was removed: its `(chrom,strand)` + `win` filters, correct only for the O2
    // copy-assignment input, structurally erased those same-chrom families.)
    let mut out = families_from_conflict_graph(&reps, &edges, min_copies, cfg);
    let n_xchrom = out
        .iter()
        .filter(|fam| fam.iter().map(|c| c.chrom.as_str()).collect::<BTreeSet<_>>().len() > 1)
        .count();
    eprintln!(
        "[gw-catalog-xchrom] {} families (>= {} copies), of which {} CROSS-CHROMOSOME",
        out.len(),
        min_copies,
        n_xchrom
    );

    // POA-CORE COMPLETION (opt-in): attach loosely-related paralogs at NEW loci the read-conflict graph misses
    // (it links only copies reads confuse). Seeded by the families above ("when assignment is determined
    // needed"); bounded by the minimizer prefilter to family-adjacent candidates. Default OFF -> no-op.
    if cfg.complete_poa_core {
        let tid2idx: std::collections::HashMap<&str, usize> =
            reps.iter().enumerate().map(|(i, r)| (r.tid.as_str(), i)).collect();
        let fam_repidx: Vec<Vec<usize>> = out
            .iter()
            .map(|fam| fam.iter().filter_map(|c| tid2idx.get(c.tid.as_str()).copied()).collect())
            .collect();
        let adds = crate::vg_family::family_detect::poa_core_completion_adds(&reps, &fam_repidx, &cfg.detect);
        let mut n_added = 0usize;
        for (f, extra) in adds.iter().enumerate() {
            for &ri in extra {
                out[f].push(reps[ri].clone());
                n_added += 1;
            }
        }
        eprintln!(
            "[gw-catalog-xchrom] POA-core completion (t_core={:.2}): +{} loosely-related paralog copies attached to {} families",
            cfg.detect.t_core,
            n_added,
            adds.iter().filter(|a| !a.is_empty()).count()
        );
    }
    Ok(out)
}

/// RefineParams for the homology-primary (E_r) path. `min_identity` (if given) is the EFFECTIVE
/// nucleotide identity floor: it sets BOTH the asm20 tier and the sensitive tier, so the union's
/// floor is exactly `mi` (0.98 = Soto SD98 mode). None -> defaults (asm20 0.80 / sensitive 0.60).
pub fn homology_refine_params(min_identity: Option<f64>, threads: usize) -> RefineParams {
    let mut p = RefineParams { threads, ..Default::default() };
    if let Some(mi) = min_identity {
        p.min_identity = mi;
        p.sensitive_identity = mi; // BOTH tiers -> effective floor = mi (not .min())
    }
    p
}

/// GENOME-WIDE homology-primary (E_r) family catalog. reps → E_r edges → γ-quasi-clique blocks →
/// ≥2 distinct loci → families. Chrom/strand-agnostic; a superset of the conflict catalog.
pub fn detect_homology_catalog_genome_wide(
    bam_path: &str,
    fasta_path: &str,
    threads: usize,
    min_copies: usize,
    cfg: &DenovoConfig,
    refine: &RefineParams,
    gamma: f64,
) -> Result<Vec<Vec<DenovoTranscript>>> {
    // --- reps (identical to the conflict path's rep build) ---
    let reads = primary_reads_from_bam(bam_path, threads)?;
    let contigs: HashSet<String> = reads.iter().map(|r| r.chrom.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;
    let skeletons = pass1_skeletons_robust(&reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    drop(reads);
    let transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
    let rep_idx = collapse_loci_span_aware(&transcripts, &cfg.detect);
    let reps: Vec<DenovoTranscript> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
    drop(transcripts);
    eprintln!("[gw-catalog-homology] {} skeletons -> {} reps over {} contigs", skeletons.len(), reps.len(), contigs.len());

    // --- E_r edges + γ-quasi-clique blocks ---
    let edges2 = homology_edges_all_reps(&reps, refine)?;
    let edges3: Vec<(usize, usize, f64)> = edges2.iter().map(|&(a, b)| (a, b, 1.0)).collect();
    let blocks = crate::vg_family::family_split::gamma_quasi_clique_partition(reps.len(), &edges3, gamma);

    let mut out: Vec<Vec<DenovoTranscript>> = Vec::new();
    for block in blocks {
        if block.len() < min_copies { continue; }
        let copies: Vec<DenovoTranscript> = block.iter().map(|&i| reps[i].clone()).collect();
        let loci = distinct_locus_reps(copies); // ≥2 spatially-distinct loci certificate
        if loci.len() >= min_copies {
            out.push(loci);
        }
    }
    eprintln!("[gw-catalog-homology] {} E_r edges -> {} families (>= {} distinct loci)", edges2.len(), out.len(), min_copies);
    Ok(out)
}

/// Parameters for the exon-sum (FLNC) homology refinement. The defaults match the validated operating
/// point (`bench/validate_exon_sum.py`): minimap2 asm20, identity >= 0.80 (asm20's native divergence
/// envelope), coverage-of-shorter >= 0.50 (more than half the shorter spliced sequence aligns).
#[derive(Clone, Debug)]
pub struct RefineParams {
    pub min_identity: f64,
    pub min_coverage: f64,
    pub minimap2: String,
    pub threads: usize,
    /// Align the GENOMIC span (introns INCLUDED) instead of the exon-sum (spliced) sequence. The exon-sum
    /// drops introns, so two copies with identical exons but DIVERGENT introns look identical (the K=0
    /// identifiability frontier); including introns lets that divergence separate them, and adds the
    /// intron/flank signal the exon-sum cannot see. Requires `intron_fasta`. Tradeoff: introns diverge
    /// faster than exons, so at a fixed identity this is STRICTER — it can drop older paralogs whose
    /// introns have diverged below the gate. Default off (exon-sum is the robust homology-detection level).
    pub include_introns: bool,
    /// Genome FASTA path used to fetch genomic spans when `include_introns` is set.
    pub intron_fasta: Option<String>,
    /// Add the SENSITIVE NUCLEOTIDE divergent tier (`minimap2 -k11 -w5`, identity >= 0.70) to the asm20 core.
    /// asm20 cannot SEED a pair below ~80% identity (lowering its gate is a no-op — its ~0.82 floor is the
    /// seed-sensitivity envelope), so the smaller k=11 seed recovers nucleotide-divergent paralogs (e.g.
    /// rapidly-evolving KZNF / SAFB families) down to ~73-76% identity. Cheap (one extra minimap2/family) and
    /// validated to add real families with 0 false merges (the `min_coverage` floor is the false-merge
    /// defense). DEFAULT ON (promoted into the default refinement).
    pub nucleotide_sensitive: bool,
    /// Sensitive-tier nucleotide identity floor for the E_r homology edge (`-k11 -w5`). Lowered from the
    /// old within-refine 0.70 to reach ancient paralogs (KRAB-ZNF ~0.62). Repeat bridges are held off by
    /// `min_coverage`. Fixed by the family P/R sweep.
    pub sensitive_identity: f64,
    /// Add the PROTEIN divergent tier (longest-ORF 6-frame translate → `mmseqs` all-vs-all, fident >= 0.50,
    /// qcov/tcov >= `min_coverage`) — recovers SYNONYMOUS-divergent CODING paralogs (e.g. the RABL2B retrocopy
    /// family at 87-99% protein but only 70% genomic identity) that nucleotide seeds can NEVER anchor. Additive
    /// only (merges divergent copies into a family, never splits asm20's calls). DEFAULT OFF — needs `mmseqs`
    /// (absent → contributes no edges) and is the costlier tier. Batched into one `mmseqs` run across all
    /// families (within-family hits only). The qcov/tcov guard rejects lone-shared-domain merges.
    pub protein_tail: bool,
    pub mmseqs: String,
}

impl Default for RefineParams {
    fn default() -> Self {
        RefineParams {
            min_identity: 0.80,
            min_coverage: 0.50,
            minimap2: std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string()),
            threads: 4,
            include_introns: false,
            intron_fasta: None,
            nucleotide_sensitive: true,
            sensitive_identity: 0.60,
            protein_tail: false,
            mmseqs: std::env::var("RUSTLE_MMSEQS").unwrap_or_else(|_| "mmseqs".to_string()),
        }
    }
}

/// Exon-sum (FLNC) homology refinement of a raw conflict-graph catalog. The raw cross-chromosome conflict
/// graph OVER-MERGES: read-conflict fires on any shared region — including a 150–300 bp Alu SINE covering
/// only a few percent of a transcript — so unrelated genes get bridged into one "family". This converts a
/// family into a clean combinatorial object by requiring BOTH criteria a real multi-copy family must meet:
///   (i)  the copies are MUTUALLY HOMOLOGOUS full-length — their spliced exon-sum sequences all-vs-all align
///        (minimap2 asm20, identity >= `min_identity`, coverage-of-shorter >= `min_coverage`); a family is
///        split into the connected components of this homology graph (read-conflict can link a fragment;
///        full-length alignment with a coverage floor cannot); and
///   (ii) the component spans >= 2 spatially-DISTINCT loci — distinct paralog copies occupy DISJOINT genomic
///        spans, so two copies are the SAME locus iff their spans overlap on the same (chrom, strand) (a gene
///        plus its own nested fragment / a second isoform). This is the distinct-LOCUS guarantee that
///        homology alone cannot give (two isoforms of one gene trivially align), and it is threshold-free.
/// Each surviving component is emitted as a refined family of its per-locus representatives (widest span).
/// `minimap2` is required; if it cannot be spawned this returns an error (the caller can fall back to the
/// raw catalog). Independent of annotation → does not depend on a reference gene set.
pub fn refine_families_exon_sum(
    families: Vec<Vec<DenovoTranscript>>,
    params: &RefineParams,
) -> Result<Vec<Vec<DenovoTranscript>>> {
    // When including introns, load the genome once for the contigs present across all families.
    let genome: Option<GenomeIndex> = if params.include_introns {
        let path = params.intron_fasta.as_ref().ok_or_else(|| {
            anyhow::anyhow!("include_introns set but intron_fasta (genome path) is None")
        })?;
        let contigs: HashSet<String> =
            families.iter().flatten().map(|c| c.chrom.clone()).collect();
        Some(GenomeIndex::from_fasta_contigs(path, &contigs)?)
    } else {
        None
    };
    // PROTEIN tier (opt-in): one batched mmseqs run across ALL families' ORFs (within-family hits only),
    // instead of a subprocess per family. `fam_protein[f]` = the protein homology edges of family f.
    let fam_protein: Vec<Vec<(usize, usize)>> = if params.protein_tail {
        batch_protein_edges(&families, 0.50, params.min_coverage, params)?
    } else {
        Vec::new()
    };
    let mut refined: Vec<Vec<DenovoTranscript>> = Vec::new();
    for (fi, fam) in families.iter().enumerate() {
        if fam.len() < 2 {
            continue;
        }
        // exon-sum (spliced) by default; genomic span (introns INCLUDED) when include_introns is set.
        let seqs: Vec<Vec<u8>> = fam.iter().map(|c| refine_copy_seq(c, genome.as_ref())).collect();
        // base detector: asm20 on the configured sequence (the validated, high-precision recent core).
        let mut edge_set: BTreeSet<(usize, usize)> =
            nucleotide_edges(&seqs, &["-x", "asm20"], params.min_identity, params.min_coverage, params)?
                .into_iter()
                .collect();
        // divergent tiers are computed on the EXON-SUM (protein ORFs need the spliced sequence; the
        // sensitive nucleotide seed is also cleanest on the spliced copy). Edges are UNIONed in.
        if params.nucleotide_sensitive {
            let exon_seqs: Vec<Vec<u8>> = fam.iter().map(|c| c.seq.clone()).collect();
            for e in nucleotide_edges(&exon_seqs, &["-k", "11", "-w", "5"], 0.70, params.min_coverage, params)? {
                edge_set.insert(e);
            }
        }
        if params.protein_tail {
            for &e in fam_protein.get(fi).map(|v| v.as_slice()).unwrap_or(&[]) {
                edge_set.insert(e);
            }
        }
        let edges: Vec<(usize, usize)> = edge_set.into_iter().collect();
        for comp in homology_components(fam.len(), &edges) {
            if comp.len() < 2 {
                continue;
            }
            let comp_copies: Vec<DenovoTranscript> = comp.iter().map(|&i| fam[i].clone()).collect();
            let loci = distinct_locus_reps(comp_copies);
            if loci.len() >= 2 {
                refined.push(loci);
            }
        }
    }
    Ok(refined)
}

/// The sequence used to compare a copy in the refinement: its exon-sum (spliced) by default, or its GENOMIC
/// span (introns + flanks INCLUDED, in transcription orientation) when a genome is supplied (`include_introns`).
fn refine_copy_seq(copy: &DenovoTranscript, genome: Option<&GenomeIndex>) -> Vec<u8> {
    match genome {
        Some(g) => {
            let mut s = g.fetch_sequence(&copy.chrom, copy.start, copy.end).unwrap_or_default();
            if s.is_empty() {
                return copy.seq.clone(); // fetch miss → fall back to the exon-sum
            }
            if copy.strand == '-' {
                s = crate::vg::reverse_complement(&s);
            }
            s
        }
        None => copy.seq.clone(),
    }
}

/// All-vs-all minimap2 alignment of a family's per-copy sequences (`seqs[i]`) with the given preset/params
/// `mm_args` (e.g. `["-x","asm20"]` or `["-k","11","-w","5"]`) → the set of homology edges `(i, j)` (i < j)
/// whose alignment passes `identity >= min_id` AND `aligned-fraction-of-shorter >= min_cov`. One subprocess
/// per family (all sequences written to a single temp FASTA, `-X` all-vs-all). minimap2 tries both strands,
/// so a `+`/`-` paralog pair is still detected. Spawn failure → `Err`; alignment failure → no edges.
pub(crate) fn nucleotide_edges(
    seqs: &[Vec<u8>],
    mm_args: &[&str],
    min_id: f64,
    min_cov: f64,
    params: &RefineParams,
) -> Result<Vec<(usize, usize)>> {
    use std::io::Write;
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    let nonce = seqs.iter().map(|s| s.len()).sum::<usize>().wrapping_mul(1000003)
        ^ seqs.len()
        ^ mm_args.len().wrapping_mul(7);
    let path = dir.join(format!("rustle_refine_{pid}_{nonce}.fa"));
    struct Cleanup(std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    let _cl = Cleanup(path.clone());
    {
        let mut f = std::fs::File::create(&path)?;
        for (i, s) in seqs.iter().enumerate() {
            writeln!(f, ">{i}")?;
            f.write_all(s)?;
            writeln!(f)?;
        }
    }
    let out = std::process::Command::new(&params.minimap2)
        .args(["-c", "-X", "--no-long-join", "-t"])
        .arg(params.threads.to_string())
        .args(mm_args)
        .arg(&path)
        .arg(&path)
        .output()
        .map_err(|e| anyhow::anyhow!("failed to run minimap2 ('{}') for refinement: {e}", params.minimap2))?;
    if !out.status.success() {
        return Ok(Vec::new()); // an alignment-time failure on one family → no edges (family dissolves)
    }
    let text = String::from_utf8_lossy(&out.stdout);
    let mut edge_set: BTreeSet<(usize, usize)> = BTreeSet::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let (q, t) = (f[0].parse::<usize>().ok(), f[5].parse::<usize>().ok());
        let (q, t) = match (q, t) {
            (Some(q), Some(t)) if q != t => (q, t),
            _ => continue,
        };
        let ql = f[1].parse::<f64>().unwrap_or(0.0);
        let qs = f[2].parse::<f64>().unwrap_or(0.0);
        let qe = f[3].parse::<f64>().unwrap_or(0.0);
        let tl = f[6].parse::<f64>().unwrap_or(0.0);
        let de = f[12..].iter().find_map(|x| x.strip_prefix("de:f:").and_then(|v| v.parse::<f64>().ok()));
        let ident = match de {
            Some(d) => 1.0 - d,
            None => {
                let nmatch = f[9].parse::<f64>().unwrap_or(0.0);
                let alnlen = f[10].parse::<f64>().unwrap_or(1.0).max(1.0);
                nmatch / alnlen
            }
        };
        let shorter = ql.min(tl).max(1.0);
        let cov = (qe - qs) / shorter;
        if ident >= min_id && cov >= min_cov {
            edge_set.insert((q.min(t), q.max(t)));
        }
    }
    Ok(edge_set.into_iter().collect())
}

/// E_r homology edges over ALL reps' exon-sum sequences: asm20 (recent) ∪ sensitive -k11 -w5 (ancient),
/// both gated by `min_coverage`. One minimap2 all-vs-all per tier over the whole rep set (minimap2's index
/// is the prefilter). Protein is NOT an edge here — it is orthogonal QC (see per-family protein_coheres).
pub(crate) fn homology_edges_all_reps(
    reps: &[DenovoTranscript],
    params: &RefineParams,
) -> Result<Vec<(usize, usize)>> {
    let seqs: Vec<Vec<u8>> = reps.iter().map(|r| r.seq.clone()).collect();
    let mut set: BTreeSet<(usize, usize)> =
        nucleotide_edges(&seqs, &["-x", "asm20"], params.min_identity, params.min_coverage, params)?
            .into_iter().collect();
    if params.nucleotide_sensitive {
        for e in nucleotide_edges(&seqs, &["-k", "11", "-w", "5"], params.sensitive_identity, params.min_coverage, params)? {
            set.insert(e);
        }
    }
    Ok(set.into_iter().collect())
}

/// Protein-homology tier, BATCHED across all families into ONE `mmseqs easy-search`. Each copy's longest ORF
/// (6-frame) is translated and written once with a `f<family>_c<copy>` id; a single all-vs-all mmseqs run
/// replaces the ~N per-family subprocess calls (the only scaling cost of the protein tier). Only WITHIN-family
/// hits are kept (cross-family protein homology is irrelevant — refinement is per conflict family), and only
/// hits passing `fident >= min_fident` AND `min(qcov, tcov) >= min_cov` (homology over MOST of the protein —
/// a lone shared domain is rejected). Returns `out[f]` = the protein edges `(ci, cj)` of family `f`. This
/// reaches SYNONYMOUS-divergent CODING paralogs (RABL2B retrocopies: ~70% genomic but 87–99% protein) that
/// nucleotide seeds cannot anchor, while structurally excluding repeat-only merges. `mmseqs` absent/failed →
/// all-empty (the nucleotide tiers still apply).
fn batch_protein_edges(
    families: &[Vec<DenovoTranscript>],
    min_fident: f64,
    min_cov: f64,
    params: &RefineParams,
) -> Result<Vec<Vec<(usize, usize)>>> {
    use std::io::Write;
    let mut out: Vec<Vec<(usize, usize)>> = vec![Vec::new(); families.len()];
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    let work = dir.join(format!("rustle_prot_{pid}"));
    let qf = work.join("q.fa");
    let res = work.join("res.m8");
    let tmp = work.join("tmp");
    struct Cleanup(std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }
    let _cl = Cleanup(work.clone());
    std::fs::create_dir_all(&work)?;
    // one fasta: id = "f{family}_c{copy}", only copies whose family has >=2 ORFs of >=40 aa.
    let mut wrote = 0usize;
    {
        let mut f = std::fs::File::create(&qf)?;
        for (fi, fam) in families.iter().enumerate() {
            let prots: Vec<Vec<u8>> = fam.iter().map(|c| longest_orf_aa(&c.seq)).collect();
            if prots.iter().filter(|p| p.len() >= 40).count() < 2 {
                continue; // family can't form a protein edge
            }
            for (ci, p) in prots.iter().enumerate() {
                if p.len() >= 40 {
                    writeln!(f, ">f{fi}_c{ci}")?;
                    f.write_all(p)?;
                    writeln!(f)?;
                    wrote += 1;
                }
            }
        }
    }
    if wrote < 2 {
        return Ok(out);
    }
    let ran = std::process::Command::new(&params.mmseqs)
        .arg("easy-search")
        .arg(&qf)
        .arg(&qf)
        .arg(&res)
        .arg(&tmp)
        .args(["-s", "7.5", "--threads"])
        .arg(params.threads.to_string())
        .args(["--format-output", "query,target,fident,qcov,tcov"])
        .output();
    match ran {
        Ok(o) if o.status.success() => {}
        _ => return Ok(out), // mmseqs absent / failed → no protein edges (graceful)
    }
    let text = match std::fs::read_to_string(&res) {
        Ok(t) => t,
        Err(_) => return Ok(out),
    };
    // parse "f{fi}_c{ci}" ids; keep only within-family, threshold-passing hits.
    let parse_id = |s: &str| -> Option<(usize, usize)> {
        let rest = s.strip_prefix('f')?;
        let (fs, cs) = rest.split_once("_c")?;
        Some((fs.parse().ok()?, cs.parse().ok()?))
    };
    let mut seen: Vec<BTreeSet<(usize, usize)>> = vec![BTreeSet::new(); families.len()];
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            continue;
        }
        let (qf_, qc) = match parse_id(f[0]) {
            Some(x) => x,
            None => continue,
        };
        let (tf, tc) = match parse_id(f[1]) {
            Some(x) => x,
            None => continue,
        };
        if qf_ != tf || qc == tc {
            continue; // cross-family or self
        }
        let fident = f[2].parse::<f64>().unwrap_or(0.0);
        let qcov = f[3].parse::<f64>().unwrap_or(0.0);
        let tcov = f[4].parse::<f64>().unwrap_or(0.0);
        if fident >= min_fident && qcov.min(tcov) >= min_cov {
            seen[qf_].insert((qc.min(tc), qc.max(tc)));
        }
    }
    for (fi, s) in seen.into_iter().enumerate() {
        out[fi] = s.into_iter().collect();
    }
    Ok(out)
}

/// Longest open reading frame across all 6 frames, translated to single-letter amino acids (stop = frame
/// break). A frame-robust CDS proxy for de-novo transcripts (no annotation): the longest stop-free run is
/// taken as the protein. Returns the AA bytes (no stop codons).
fn longest_orf_aa(seq: &[u8]) -> Vec<u8> {
    // standard genetic code, indexed by (base3 code) with T,C,A,G = 0..3.
    const AA: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    let code = |b: u8| -> Option<usize> {
        match b.to_ascii_uppercase() {
            b'T' => Some(0),
            b'C' => Some(1),
            b'A' => Some(2),
            b'G' => Some(3),
            _ => None,
        }
    };
    let translate = |s: &[u8]| -> Vec<u8> {
        let mut best: Vec<u8> = Vec::new();
        for frame in 0..3 {
            let mut cur: Vec<u8> = Vec::new();
            let mut i = frame;
            while i + 3 <= s.len() {
                let aa = match (code(s[i]), code(s[i + 1]), code(s[i + 2])) {
                    (Some(a), Some(b), Some(c)) => AA[a * 16 + b * 4 + c],
                    _ => b'*', // N / ambiguous → treat as a break
                };
                if aa == b'*' {
                    if cur.len() > best.len() {
                        best = std::mem::take(&mut cur);
                    } else {
                        cur.clear();
                    }
                } else {
                    cur.push(aa);
                }
                i += 3;
            }
            if cur.len() > best.len() {
                best = cur;
            }
        }
        best
    };
    let fwd = translate(seq);
    let rc = crate::vg::reverse_complement(seq);
    let rev = translate(&rc);
    if rev.len() > fwd.len() {
        rev
    } else {
        fwd
    }
}

/// Connected components of the homology graph over `n` copies.
fn homology_components(n: usize, edges: &[(usize, usize)]) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..n).collect();
    for &(a, b) in edges {
        uf_union(&mut parent, a, b);
    }
    let mut comp: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
    for i in 0..n {
        let r = uf_find(&mut parent, i);
        comp.entry(r).or_default().push(i);
    }
    comp.into_values().collect()
}

/// Iterative union-find (mirrors `family_detect`'s, kept local to avoid widening that module's API).
fn uf_find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}
fn uf_union(parent: &mut [usize], a: usize, b: usize) {
    let (ra, rb) = (uf_find(parent, a), uf_find(parent, b));
    if ra != rb {
        parent[ra] = rb;
    }
}

/// Opposite-strand overlapping copies collapse only when the minority copy has fewer than `1/DENOM` the
/// reads of its overlapping partner — a strand artifact (a few antisense reads on a real transcript), not
/// a balanced sense/antisense gene pair.
const ANTISENSE_MINORITY_DENOM: u32 = 10;

/// Collapse copies at the SAME locus and return one representative per spatially-DISTINCT locus. Distinct
/// paralog copies occupy DISJOINT genomic spans, so two copies overlapping on the same chromosome are the
/// same locus. SAME-strand overlap is always one locus (a gene plus its own nested fragment / a second
/// isoform — which homology alignment cannot distinguish, but disjointness can). OPPOSITE-strand overlap is
/// the sense/antisense case: collapse ONLY when one copy is a clear read-minority (a few antisense reads on
/// a real transcript = a strand artifact, e.g. GWFAM99 = 666 `+` reads vs 3 `-` reads at one locus, a
/// sense/antisense mis-split); a BALANCED overlapping antisense pair is two genuine loci and is kept.
/// Threshold-free for the same-strand case; the antisense case uses an order-of-magnitude read-minority
/// guard. The representative of a locus is the MOST-supported copy (most reads — the real one, not the
/// minority artifact — then widest span).
fn distinct_locus_reps(copies: Vec<DenovoTranscript>) -> Vec<DenovoTranscript> {
    let n = copies.len();
    let mut parent: Vec<usize> = (0..n).collect();
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = (&copies[i], &copies[j]);
            let same_pos = a.chrom == b.chrom && a.end.min(b.end) > a.start.max(b.start);
            if !same_pos {
                continue;
            }
            let collapse = if a.strand == b.strand {
                true
            } else {
                let (lo, hi) = (a.n_reads.min(b.n_reads), a.n_reads.max(b.n_reads));
                lo.saturating_mul(ANTISENSE_MINORITY_DENOM) < hi // minority antisense = artifact
            };
            if collapse {
                uf_union(&mut parent, i, j);
            }
        }
    }
    // representative per locus = MOST reads (the real copy, not the minority artifact), then widest span.
    let key = |t: &DenovoTranscript| (t.n_reads, t.end - t.start);
    let mut by_locus: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new(); // root -> rep idx
    for i in 0..n {
        let r = uf_find(&mut parent, i);
        let rep = by_locus.entry(r).or_insert(i);
        if key(&copies[i]) > key(&copies[*rep]) {
            *rep = i;
        }
    }
    let mut reps: Vec<usize> = by_locus.into_values().collect();
    reps.sort_unstable();
    reps.into_iter().map(|i| copies[i].clone()).collect()
}

#[cfg(test)]
mod tests {
    use super::super::copy_split::AlignedRead;
    use super::super::family_detect::collapse_loci;
    use super::*;

    #[test]
    fn homology_refine_params_min_identity_sets_both_tiers() {
        let p = homology_refine_params(Some(0.98), 4);
        assert_eq!(p.min_identity, 0.98);
        assert_eq!(p.sensitive_identity, 0.98, "Soto mode: sensitive tier must also be 0.98, not left at 0.60");
        let d = homology_refine_params(None, 4);
        assert_eq!(d.min_identity, 0.80);
        assert_eq!(d.sensitive_identity, 0.60);
    }

    struct SplitMix64(u64);
    impl SplitMix64 {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
    }
    fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut rng = SplitMix64(seed);
        const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
        (0..n).map(|_| B[(rng.next_u64() % 4) as usize]).collect()
    }
    fn read(chrom: &str, s: u64, e: u64, introns: &[(u64, u64)]) -> PrimaryRead {
        PrimaryRead { chrom: chrom.into(), ref_start: s, ref_end: e, introns: introns.to_vec() }
    }
    fn cat(parts: &[&[u8]]) -> Vec<u8> {
        parts.iter().flat_map(|p| p.iter().copied()).collect()
    }

    /// Two paralog genes on one contig sharing a 300 bp homologous core (distinct flanks, distinct intron
    /// coords). Each gene: exon1 [base,base+200), canonical intron [base+200,base+220) GT..AG, exon2
    /// [base+220,base+420). Spliced = flank(50) + core(300) + flank(50). Three reads per gene.
    fn two_paralogs() -> (GenomeIndex, Vec<PrimaryRead>) {
        let core = rand_seq(300, 0xC0FE_CAFE);
        let mut g = vec![b'A'; 1500];
        let mut place = |g: &mut Vec<u8>, base: usize, f1: u64, f2: u64| {
            g[base..base + 50].copy_from_slice(&rand_seq(50, f1));
            g[base + 50..base + 200].copy_from_slice(&core[0..150]);
            // intron [base+200, base+220): GT..AG
            g[base + 200] = b'G';
            g[base + 201] = b'T';
            g[base + 218] = b'A';
            g[base + 219] = b'G';
            g[base + 220..base + 370].copy_from_slice(&core[150..300]);
            g[base + 370..base + 420].copy_from_slice(&rand_seq(50, f2));
        };
        place(&mut g, 0, 0xA1, 0xA2);
        place(&mut g, 1000, 0xB1, 0xB2);
        let genome = GenomeIndex::from_seqs(&[("c1", &g)]);
        let mut reads = Vec::new();
        for _ in 0..3 {
            reads.push(read("c1", 0, 420, &[(200, 220)]));
        }
        for _ in 0..3 {
            reads.push(read("c1", 1000, 1420, &[(1200, 1220)]));
        }
        (genome, reads)
    }

    #[test]
    fn detect_families_recovers_a_two_copy_family() {
        let (genome, reads) = two_paralogs();
        let r = detect_families(&reads, &genome, &DenovoConfig::default());
        assert_eq!(r.n_transcripts, 2, "two gated transcripts (one per gene)");
        assert_eq!(r.n_loci, 2, "distinct intron chains -> two loci");
        assert_eq!(r.n_edges, 1, "the two paralogs share a homologous core -> one POA edge");
        assert_eq!(r.families.len(), 1, "one decomposed family");
        let fam = &r.families[0];
        assert_eq!(fam.members.len(), 2);
        assert_eq!(fam.n_chroms, 1);
        assert_eq!(fam.class, FamilyClass::Family);
        let mut members = fam.members.clone();
        members.sort();
        assert_eq!(members, vec!["DN_c1_0_2".to_string(), "DN_c1_1000_2".to_string()]);
    }

    /// Two paralogs sharing a core but DIFFERING at 3 PSV offsets (so detection finds a family AND reads are
    /// resolvable), plus one copyB read aligned to copyA's region (a multimapper).
    fn two_paralogs_with_psvs() -> (GenomeIndex, Vec<PrimaryRead>, Vec<BamRead>) {
        let base = rand_seq(300, 0xC0FE_D1FF);
        let psv = [60usize, 150, 240];
        let (mut core_a, mut core_b) = (base.clone(), base);
        for &p in &psv {
            core_a[p] = b'A';
            core_b[p] = b'C';
        }
        let (fa1, fa2) = (rand_seq(50, 0xA1), rand_seq(50, 0xA2));
        let (fb1, fb2) = (rand_seq(50, 0xB1), rand_seq(50, 0xB2));
        let mut g = vec![b'A'; 1500];
        let mut put = |g: &mut Vec<u8>, at: usize, s: &[u8]| g[at..at + s.len()].copy_from_slice(s);
        // copyA at 0, copyB at 1000: exon1 = flank(50)+core[0:150], intron GT..AG, exon2 = core[150:300]+flank(50)
        for (base_off, f1, core, f2) in [(0usize, &fa1, &core_a, &fa2), (1000, &fb1, &core_b, &fb2)] {
            put(&mut g, base_off, f1);
            put(&mut g, base_off + 50, &core[0..150]);
            g[base_off + 200] = b'G';
            g[base_off + 201] = b'T';
            g[base_off + 218] = b'A';
            g[base_off + 219] = b'G';
            put(&mut g, base_off + 220, &core[150..300]);
            put(&mut g, base_off + 370, f2);
        }
        let genome = GenomeIndex::from_seqs(&[("c1", &g)]);
        let mut primary = Vec::new();
        for _ in 0..3 {
            primary.push(read("c1", 0, 420, &[(200, 220)]));
        }
        for _ in 0..3 {
            primary.push(read("c1", 1000, 1420, &[(1200, 1220)]));
        }
        // a copyB read aligned to copyA's genomic region: seq = copyB's spliced transcript.
        let copyb_spliced = cat(&[&fb1, &core_b, &fb2]);
        let copya_spliced = cat(&[&fa1, &core_a, &fa2]);
        // Three multimappers (readB/readC/readD): each primary at locus 0, secondary at locus 1000.
        // de 0.010/0.012 (both < 0.05, |Δ|=0.002 < 0.005) -> de-tied -> conflict edge. min_reads=3 satisfied.
        let mk = |name: &str, primary_seq: Vec<u8>, secondary_seq: Vec<u8>| {
            vec![
                BamRead { chrom: "c1".into(),
                    read: AlignedRead { ref_start: 0, cigar: vec![('M', 200), ('N', 20), ('M', 200)], seq: primary_seq, qual: vec![] },
                    mapq: 60, name: name.into(), as_score: 380, de: 0.010, is_supplementary: false },
                BamRead { chrom: "c1".into(),
                    read: AlignedRead { ref_start: 1000, cigar: vec![('M', 200), ('N', 20), ('M', 200)], seq: secondary_seq, qual: vec![] },
                    mapq: 0, name: name.into(), as_score: 379, de: 0.012, is_supplementary: false },
            ]
        };
        let mut bam = Vec::new();
        for nm in ["readB", "readC", "readD"] {
            bam.extend(mk(nm, copyb_spliced.clone(), copya_spliced.clone()));
        }
        (genome, primary, bam)
    }

    /// Real-data END-TO-END run: detect families + assign reads over a co-located family region. Ignored;
    /// run in RELEASE with a co-located cluster region (e.g. MAGEA on GGO):
    ///   RUSTLE_DENOVO_SMOKE_BAM=GGO.bam RUSTLE_DENOVO_SMOKE_FASTA=GGO.fasta \
    ///   RUSTLE_DENOVO_SMOKE_REGION=NC_073247.2:161251228-164865959 \
    ///     cargo test --release --lib -- --ignored smoke_detect_and_assign_real --nocapture
    #[test]
    #[ignore = "needs a real BAM + FASTA + REGION"]
    fn smoke_detect_and_assign_real() {
        use crate::vg_family::denovo_assemble::reads_in_region;
        let (bam, fasta, region) = match (
            std::env::var("RUSTLE_DENOVO_SMOKE_BAM"),
            std::env::var("RUSTLE_DENOVO_SMOKE_FASTA"),
            std::env::var("RUSTLE_DENOVO_SMOKE_REGION"),
        ) {
            (Ok(b), Ok(f), Ok(r)) => (b, f, r),
            _ => return,
        };
        let (chrom, range) = region.split_once(':').unwrap();
        let (lo, hi) = range.split_once('-').unwrap();
        let (lo, hi): (u64, u64) = (lo.parse().unwrap(), hi.parse().unwrap());
        let (primary, bam_reads) = reads_in_region(&bam, chrom, lo, hi, 4).expect("read region");
        let mut contigs = std::collections::HashSet::new();
        contigs.insert(chrom.to_string());
        let genome = GenomeIndex::from_fasta_contigs(&fasta, &contigs).expect("fasta");
        eprintln!("region {chrom}:{lo}-{hi}: {} primary reads, {} mapped reads", primary.len(), bam_reads.len());
        let (fas, fallback, _dna_needs) = detect_and_assign(
            &primary,
            &bam_reads,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
            &[],
            false,
            &fasta,
        );
        if !fallback.is_empty() {
            eprintln!("family edges via large-seq fallback (auditable): {}", fallback.len());
        }
        eprintln!("co-located families with assignments: {}", fas.len());
        for fa in &fas {
            let pct = |x: usize| if fa.n_reads > 0 { 100.0 * x as f64 / fa.n_reads as f64 } else { 0.0 };
            eprintln!(
                "  {} {} copies={} reads={} PSVcols={} resolv_PSV={:.0}% resolv_+J={:.0}% J_only={} assign_+J={:.0}% silver={}/{}",
                fa.family_id, fa.chrom, fa.n_copies, fa.n_reads, fa.psv_cols,
                pct(fa.resolvable_psv), pct(fa.resolvable_j), fa.junction_only,
                pct(fa.assigned_j), fa.uniq_agree, fa.uniq,
            );
        }
        let (tot_uniq, tot_agree): (usize, usize) = fas.iter().fold((0, 0), |(u, a), f| (u + f.uniq, a + f.uniq_agree));
        if tot_uniq > 0 {
            eprintln!("AGGREGATE silver-standard unique-mapper agreement: {tot_agree}/{tot_uniq} ({:.1}%)", 100.0 * tot_agree as f64 / tot_uniq as f64);
        }
    }

    #[test]
    fn detect_and_assign_resolves_multimapper_end_to_end() {
        let (genome, primary, aligned) = two_paralogs_with_psvs();
        let (fas, fallback, dna_needs) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
            &[],
            false,
            "",
        );
        assert!(fallback.is_empty(), "small paralogs use the exact poasta path, no fallback");
        assert!(dna_needs.is_empty(), "absent_copies=false must return empty dna_needs vec");
        assert_eq!(fas.len(), 1, "one co-located 2-copy family");
        let fa = &fas[0];
        assert_eq!(fa.n_copies, 2);
        // Six BamRecords: readB/readC/readD each primary+secondary; all overlap the family.
        assert_eq!(fa.n_reads, 6);
        assert_eq!(fa.assignments.len(), 6);
        // copies sorted by start: copyA=0, copyB=1.
        // Primary (ref_start=0, seq=copyb_spliced) -> assigned to copyB (best_copy=1).
        let primary_assign = fa.assignments.iter().find(|(_, a)| a.best_copy == 1)
            .expect("primary read (copyB seq at locus 0) should resolve to copyB (copy 1)");
        assert_eq!(primary_assign.1.status, super::super::copy_assign::AssignStatus::Assigned);
        // Secondary (ref_start=1000, seq=copya_spliced) -> assigned to copyA (best_copy=0).
        let secondary_assign = fa.assignments.iter().find(|(_, a)| a.best_copy == 0)
            .expect("secondary read (copyA seq at locus 1000) should resolve to copyA (copy 0)");
        assert_eq!(secondary_assign.1.status, super::super::copy_assign::AssignStatus::Assigned);
    }

    /// When `absent_copies=false`, the third element of `detect_and_assign`'s return tuple must
    /// always be empty — the admission block is entirely skipped, so no `DnaNeedsRecord`s are
    /// produced regardless of how many collapsed-copy candidates exist at the loci.
    #[test]
    fn detect_and_assign_absent_copies_off_returns_empty_dna_needs() {
        let (genome, primary, aligned) = two_paralogs_with_psvs();
        let (_, _, dna_needs) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
            &[],
            false, // OFF — the admission block must be completely skipped
            "",
        );
        assert!(
            dna_needs.is_empty(),
            "absent_copies=false must yield an empty dna_needs vec (byte-identical OFF path)"
        );
    }

    #[test]
    fn colocated_families_filters_same_chrom_clusters() {
        // Build a family with 2 same-chrom copies (via the detect chain), assert it is co-located.
        let core = rand_seq(400, 0xCAFE_01);
        let mut g = vec![b'A'; 2000];
        let mut put = |g: &mut Vec<u8>, at: usize, s: &[u8]| g[at..at + s.len()].copy_from_slice(s);
        for base_off in [0usize, 1000] {
            put(&mut g, base_off, &rand_seq(80, 0xF0 + base_off as u64));
            put(&mut g, base_off + 80, &core[0..160]);
            g[base_off + 240] = b'G';
            g[base_off + 241] = b'T';
            g[base_off + 258] = b'A';
            g[base_off + 259] = b'G';
            put(&mut g, base_off + 260, &core[160..400]);
            put(&mut g, base_off + 500, &rand_seq(80, 0xE0 + base_off as u64));
        }
        let genome = GenomeIndex::from_seqs(&[("c1", &g)]);
        let mut primary = Vec::new();
        for base_off in [0u64, 1000] {
            for _ in 0..3 {
                primary.push(read("c1", base_off, base_off + 580, &[(base_off + 240, base_off + 260)]));
            }
        }
        let cfg = DenovoConfig::default();
        let skeletons = pass1_skeletons(&primary, cfg.pass1_min_reads);
        let transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
        let reps: Vec<DenovoTranscript> =
            collapse_loci(&transcripts).iter().map(|&i| transcripts[i].clone()).collect();
        let split = decompose_families(&detect_edges(&reps, &cfg.detect), &cfg.split);
        let colo = colocated_families(&reps, &split, 5_000_000, 2, &cfg.detect);
        assert_eq!(colo.len(), 1);
        assert_eq!(colo[0].copies.len(), 2);
        assert_eq!(colo[0].chrom, "c1");
        // same family but min_copies=3 -> not co-located (only 2 copies)
        assert!(colocated_families(&reps, &split, 5_000_000, 3, &cfg.detect).is_empty());
    }

    #[test]
    fn detect_families_single_copy_yields_no_family() {
        // one gene only -> no homology edge -> no family.
        let (genome, reads) = two_paralogs();
        let single: Vec<PrimaryRead> = reads.into_iter().take(3).collect(); // gene A only
        let r = detect_families(&single, &genome, &DenovoConfig::default());
        assert_eq!(r.n_transcripts, 1);
        assert_eq!(r.n_edges, 0);
        assert!(r.families.is_empty());
    }

    /// Real-data end-to-end run. Ignored by default. Run in RELEASE (POA is heavy), optionally scoped to a
    /// locus to keep the POA-pair count small:
    ///   RUSTLE_DENOVO_SMOKE_BAM=...bam RUSTLE_DENOVO_SMOKE_FASTA=...fasta \
    ///   RUSTLE_DENOVO_SMOKE_REGION=NC_073243.2:104789647-104877901 \
    ///     cargo test --release --lib -- --ignored smoke_detect_families --nocapture
    #[test]
    #[ignore = "needs a real BAM + FASTA via RUSTLE_DENOVO_SMOKE_{BAM,FASTA}"]
    fn smoke_detect_families_on_real_data() {
        use std::collections::HashSet as Set;
        let (bam, fasta) = match (
            std::env::var("RUSTLE_DENOVO_SMOKE_BAM"),
            std::env::var("RUSTLE_DENOVO_SMOKE_FASTA"),
        ) {
            (Ok(b), Ok(f)) => (b, f),
            _ => return,
        };
        let mut reads = primary_reads_from_bam(&bam, 4).expect("read BAM");
        // optional region scope: chrom:start-end (keep reads overlapping the window).
        if let Ok(region) = std::env::var("RUSTLE_DENOVO_SMOKE_REGION") {
            let (chrom, range) = region.split_once(':').expect("chrom:start-end");
            let (s, e) = range.split_once('-').expect("start-end");
            let (s, e): (u64, u64) = (s.parse().unwrap(), e.parse().unwrap());
            reads.retain(|r| r.chrom == chrom && r.ref_start < e && r.ref_end > s);
        }
        let cfg = DenovoConfig::default();
        let contigs: Set<String> = reads.iter().map(|r| r.chrom.clone()).collect();
        let genome = GenomeIndex::from_fasta_contigs(&fasta, &contigs).expect("load FASTA contigs");

        // instrument the cost driver: how many POA pairs survive the k-mer + contiguous-span pre-filter.
        let skeletons = pass1_skeletons(&reads, cfg.pass1_min_reads);
        let transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
        let rep_idx = collapse_loci(&transcripts);
        let reps: Vec<_> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
        let maxlen = reps.iter().map(|r| r.seq.len()).max().unwrap_or(0);
        let big = reps.iter().filter(|r| r.seq.len() > 5000).count();
        let pairs = crate::vg_family::family_detect::candidate_pairs(&reps, &cfg.detect);
        eprintln!(
            "STAGED: reads={} skeletons={} transcripts={} loci={} maxlen={} big(>5kb)={} POA_pairs={}",
            reads.len(),
            skeletons.len(),
            transcripts.len(),
            reps.len(),
            maxlen,
            big,
            pairs.len()
        );
        // diagnostic: pairwise k-mer overlap + POA core coverage of the loci (why a pair did/didn't survive).
        for i in 0..reps.len().min(6) {
            for j in (i + 1)..reps.len().min(6) {
                let ki = crate::vg_family::family_detect::canonical_kmer_first_pos(&reps[i].seq);
                let kj = crate::vg_family::family_detect::canonical_kmer_first_pos(&reps[j].seq);
                let shared = ki.keys().filter(|c| kj.contains_key(*c)).count();
                if shared == 0 {
                    continue;
                }
                let cov = crate::vg_family::family_graph::contiguous_core_coverage(&reps[i].seq, &reps[j].seq);
                eprintln!(
                    "  DIAG {}({}bp) <-> {}({}bp): shared_kmers={} core_cov={:.3}",
                    reps[i].tid, reps[i].seq.len(), reps[j].tid, reps[j].seq.len(), shared, cov
                );
            }
        }

        let r = detect_families(&reads, &genome, &cfg);
        let webs = r.families.iter().filter(|f| f.class == FamilyClass::Web).count();
        eprintln!(
            "DETECT: {} loci -> {} edges -> {} families ({} discrete + {} webs)",
            r.n_loci,
            r.n_edges,
            r.families.len(),
            r.families.len() - webs,
            webs,
        );
        for fam in r.families.iter().take(12) {
            eprintln!(
                "  family n={} chroms={} density={:.3} avg_core_recip={:.3} {:?}: {}",
                fam.members.len(),
                fam.n_chroms,
                fam.density,
                fam.avg_core_recip,
                fam.class,
                fam.members.join(",")
            );
        }
    }

    /// DIAGNOSTIC: dump the homology prefilter internals for a region so we can see WHY a copy-copy pair did
    /// or didn't become a candidate. Prints, per rep pair: raw shared k-mers, shared INFORMATIVE k-mers,
    /// the contiguous span (the prefilter's `core`) vs the `t_core*minlen` bar, and the GROUND-TRUTH
    /// contiguous-core coverage (bounded, so the 228 kb hub cannot OOM). Compares to candidate_pairs output.
    ///   RUSTLE_DENOVO_SMOKE_BAM=...bam RUSTLE_DENOVO_SMOKE_FASTA=...fasta \
    ///   RUSTLE_DIAG_REGION=NC_073234.2:48446103-49179309 \
    ///     cargo test --release --lib -- --ignored diag_prefilter_homology --nocapture
    #[test]
    #[ignore = "diagnostic: needs RUSTLE_DENOVO_SMOKE_{BAM,FASTA} + RUSTLE_DIAG_REGION"]
    fn diag_prefilter_homology() {
        use crate::vg_family::family_detect::canonical_kmer_first_pos;
        use crate::vg_family::family_graph::contiguous_core_coverage_bounded;
        use std::collections::{BTreeMap, BTreeSet};
        let (bam, fasta, region) = match (
            std::env::var("RUSTLE_DENOVO_SMOKE_BAM"),
            std::env::var("RUSTLE_DENOVO_SMOKE_FASTA"),
            std::env::var("RUSTLE_DIAG_REGION"),
        ) {
            (Ok(b), Ok(f), Ok(r)) => (b, f, r),
            _ => return,
        };
        let (chrom, range) = region.split_once(':').expect("chrom:start-end");
        let (s, e) = range.split_once('-').expect("start-end");
        let (lo, hi): (u64, u64) = (s.parse().unwrap(), e.parse().unwrap());
        let mut reads = primary_reads_from_bam(&bam, 4).expect("read BAM");
        reads.retain(|r| r.chrom == chrom && r.ref_start < hi && r.ref_end > lo);
        let cfg = DenovoConfig::default();
        let p = &cfg.detect;
        let contigs: HashSet<String> = std::iter::once(chrom.to_string()).collect();
        let genome = GenomeIndex::from_fasta_contigs(&fasta, &contigs).expect("fasta");
        let skeletons = pass1_skeletons(&reads, cfg.pass1_min_reads);
        let transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
        let reps: Vec<DenovoTranscript> =
            collapse_loci(&transcripts).iter().map(|&i| transcripts[i].clone()).collect();
        eprintln!("=== {} reps (cnt_min={} cnt_max={} k_share={} t_core={}) ===", reps.len(), p.cnt_min, p.cnt_max, p.k_share, p.t_core);
        for (i, r) in reps.iter().enumerate() {
            eprintln!("  [{i}] {} len={} exons={}", r.tid, r.seq.len(), r.introns.len() + 1);
        }
        // replicate candidate_pairs internals: per-rep k-mers, ownership, informative set, informative sig.
        let rep_kmers: Vec<BTreeMap<u64, u32>> = reps.iter().map(|r| canonical_kmer_first_pos(&r.seq)).collect();
        let mut owner: BTreeMap<u64, usize> = BTreeMap::new();
        for km in &rep_kmers {
            for c in km.keys() {
                *owner.entry(*c).or_insert(0) += 1;
            }
        }
        let info: BTreeSet<u64> = owner.iter().filter(|(_, &c)| c >= p.cnt_min && c <= p.cnt_max).map(|(&c, _)| c).collect();
        let sig: Vec<BTreeMap<u64, u32>> = rep_kmers
            .iter()
            .map(|km| {
                km.iter()
                    .filter(|(c, _)| info.contains(*c))
                    .map(|(&c, &x)| (c, x))
                    .collect::<BTreeMap<u64, u32>>()
            })
            .collect();
        eprintln!("=== pairwise (raw=all shared k-mers, info=shared INFORMATIVE, span=prefilter core, bar=t_core*minlen, cov=ground-truth homology) ===");
        for i in 0..reps.len() {
            for j in (i + 1)..reps.len() {
                let raw = rep_kmers[i].keys().filter(|c| rep_kmers[j].contains_key(*c)).count();
                let (mut amin, mut amax, mut bmin, mut bmax, mut common) = (u32::MAX, 0u32, u32::MAX, 0u32, 0usize);
                for (c, &pa) in &sig[i] {
                    if let Some(&pb) = sig[j].get(c) {
                        common += 1;
                        amin = amin.min(pa);
                        amax = amax.max(pa);
                        bmin = bmin.min(pb);
                        bmax = bmax.max(pb);
                    }
                }
                let span = if common > 0 { (amax - amin).min(bmax - bmin) as usize } else { 0 };
                let minlen = reps[i].seq.len().min(reps[j].seq.len());
                let bar = p.t_core * minlen as f64;
                let au = reps[i].seq.to_ascii_uppercase();
                let bu = reps[j].seq.to_ascii_uppercase();
                let bru = crate::vg::reverse_complement(&bu);
                let cov = contiguous_core_coverage_bounded(&au, &bu, p.len_cap)
                    .max(contiguous_core_coverage_bounded(&au, &bru, p.len_cap));
                // forced LCS (robust to poasta's flank-threading collapse): if lcs >> poasta cov on a small
                // pair, poasta UNDER-counts a real shared core (false family split).
                let minl = au.len().min(bu.len());
                let lcs = (crate::vg_family::family_graph::longest_common_substring(&au, &bu)
                    .max(crate::vg_family::family_graph::longest_common_substring(&au, &bru))) as f64
                    / minl as f64;
                let why = if common < p.k_share { "FAIL:k_share" } else if (span as f64) < bar { "FAIL:span" } else { "PASS" };
                eprintln!(
                    "  [{i}]x[{j}] raw={raw} info={common} span={span} bar={bar:.0} poasta={cov:.3} lcs={lcs:.3} {why}{}",
                    if cov.max(lcs) >= p.t_core { "  <-HOMOLOG" } else { "" }
                );
            }
        }
        let cps = crate::vg_family::family_detect::candidate_pairs(&reps, p);
        eprintln!("candidate_pairs accepts: {cps:?}");
    }

    fn bam_read(chrom: &str, start: u64, end: u64, name: &str, de: f32, is_secondary: bool) -> BamRead {
        let mapq = if is_secondary { 0 } else { 60 };
        let len = (end - start) as u64;
        BamRead {
            chrom: chrom.into(),
            read: AlignedRead { ref_start: start, cigar: vec![('M', len)], seq: vec![b'A'; len as usize], qual: vec![] },
            mapq, name: name.into(), as_score: 500, de, is_supplementary: false,
        }
    }

    fn rep(chrom: &str, start: u64, end: u64) -> DenovoTranscript {
        DenovoTranscript {
            tid: format!("{chrom}:{start}-{end}"),
            chrom: chrom.into(),
            start,
            end,
            n_reads: 3,
            strand: '+',
            introns: vec![],
            seq: vec![b'A'; (end - start) as usize],
        }
    }

    /// Build a rep with explicit introns + read count (for prune_same_locus structural tests).
    fn rep_s(start: u64, end: u64, introns: Vec<(u64, u64)>, n_reads: u32) -> DenovoTranscript {
        DenovoTranscript {
            tid: format!("t_{start}_{end}_{}", introns.len()),
            chrom: "c1".into(),
            start,
            end,
            n_reads,
            strand: '-',
            introns,
            seq: vec![b'A'; (end - start) as usize],
        }
    }

    #[test]
    fn prune_same_locus_drops_unspliced_readthrough_keeps_tandems() {
        // The CAFAM0 shape: a 12-exon copy + its 1-exon unspliced read-through (contains it) + a genuine
        // tandem paralog at a disjoint locus with its OWN junctions.
        let spliced = rep_s(82594889, 82620183, (0..11).map(|k| (82595000 + k * 1000, 82595500 + k * 1000)).collect(), 15);
        let unspliced = rep_s(82594891, 82781357, vec![], 4); // 1 exon, contains `spliced`
        let tandem = rep_s(82739824, 82765535, (0..7).map(|k| (82740000 + k * 1000, 82740500 + k * 1000)).collect(), 14);
        let kept = prune_same_locus(vec![spliced.clone(), unspliced, tandem.clone()], &DetectParams::default());
        let starts: Vec<u64> = kept.iter().map(|c| c.start).collect();
        assert_eq!(starts, vec![82594889, 82739824], "unspliced read-through dropped; spliced + tandem kept");
        assert_eq!(kept[0].introns.len(), 11, "kept the 12-exon spliced representative, not the 1-exon");
    }

    #[test]
    fn prune_same_locus_dedups_shared_junction_isoforms() {
        // two isoforms sharing a junction = one gene -> keep the more-structured.
        let iso_a = rep_s(100, 900, vec![(200, 300), (400, 500)], 5);
        let iso_b = rep_s(100, 700, vec![(200, 300)], 8); // shares junction (200,300)
        let kept = prune_same_locus(vec![iso_a, iso_b], &DetectParams::default());
        assert_eq!(kept.len(), 1, "shared-junction isoforms collapse to one locus");
        assert_eq!(kept[0].introns.len(), 2, "kept the 3-exon isoform");
    }

    #[test]
    fn distinct_locus_reps_collapses_same_locus_keeps_disjoint() {
        // A homology component with: copy at locus X (wide), its nested low-read fragment at X (overlaps,
        // same strand), and a real paralog at a DISJOINT locus Y. -> 2 distinct loci (X's wide rep + Y).
        let wide = rep_s(1000, 9000, vec![(2000, 3000)], 500); // strand '-'
        let mut frag = rep_s(4000, 4600, vec![], 5); // nested in `wide`, same chrom+strand
        frag.strand = '-';
        let mut para = rep_s(50000, 58000, vec![(51000, 52000)], 40); // disjoint span = distinct locus
        para.strand = '-';
        let loci = distinct_locus_reps(vec![wide.clone(), frag, para.clone()]);
        let starts: Vec<u64> = loci.iter().map(|c| c.start).collect();
        assert_eq!(starts, vec![1000, 50000], "nested fragment collapses into its locus; disjoint paralog kept");
    }

    #[test]
    fn distinct_locus_reps_balanced_antisense_pair_kept_distinct() {
        // Overlapping OPPOSITE-strand spans with BALANCED read support = a genuine sense/antisense pair: keep both.
        let mut plus = rep_s(1000, 5000, vec![(2000, 3000)], 100);
        plus.strand = '+';
        let mut minus = rep_s(1000, 5000, vec![(2000, 3000)], 80);
        minus.strand = '-';
        let loci = distinct_locus_reps(vec![plus, minus]);
        assert_eq!(loci.len(), 2, "balanced opposite-strand overlapping spans are distinct loci");
    }

    #[test]
    fn distinct_locus_reps_minority_antisense_collapses_keeps_dominant() {
        // The GWFAM99 case: a real transcript (666 `+` reads) + a few antisense reads mis-assembled as a
        // `-` copy overlapping it (3 reads). The minority antisense copy is a strand artifact -> collapse,
        // and KEEP the dominant (most-reads) copy.
        let mut dom = rep_s(105549018, 105560000, vec![(105550000, 105551000)], 666);
        dom.strand = '+';
        let mut anti = rep_s(105549018, 105552021, vec![], 3); // overlaps dom, opposite strand, 3 reads
        anti.strand = '-';
        let loci = distinct_locus_reps(vec![dom, anti]);
        assert_eq!(loci.len(), 1, "minority antisense copy collapses into the real locus");
        assert_eq!(loci[0].n_reads, 666, "the dominant (real) copy is the representative");
    }

    #[test]
    fn longest_orf_aa_translates_and_picks_longest_frame() {
        // ATG AAA GGG TTT TGT = M K G F C (no stop) → the forward longest ORF.
        let seq = b"ATGAAAGGGTTTTGT";
        assert_eq!(longest_orf_aa(seq), b"MKGFC".to_vec());
        // a stop (TAA) breaks the run; the longer side is kept.
        let with_stop = b"ATGAAATAAGGGTTTTGTCCCAAA"; // frame0: MK * GFCPK -> longest stop-free = "GFCPK"
        let orf = longest_orf_aa(with_stop);
        assert!(orf.len() >= 5, "kept the longest stop-free run, got {:?}", String::from_utf8_lossy(&orf));
        assert!(!orf.contains(&b'*'), "ORF must not contain stop codons");
    }

    #[test]
    fn homology_components_splits_a_bridged_family() {
        // copies {0,1} homologous, {2,3} homologous, NO edge between the pair = two components (the Alu-
        // bridge case where read-conflict had linked them but full-length alignment does not).
        let comps = homology_components(4, &[(0, 1), (2, 3)]);
        assert_eq!(comps.len(), 2);
        assert!(comps.iter().any(|c| c == &vec![0, 1]));
        assert!(comps.iter().any(|c| c == &vec![2, 3]));
    }

    #[test]
    fn prune_same_locus_preserves_disjoint_tandems_with_span_overlap() {
        // two real tandem copies, DISJOINT junctions, spans overlap ~40% — must BOTH survive.
        // Use divergent sequences so the span-aware homology arm does not merge them.
        let mut cp1 = rep_s(1000, 3000, vec![(1500, 1600), (2000, 2100)], 10);
        cp1.seq = rand_seq(2000, 0xA1);
        let mut cp2 = rep_s(2200, 4200, vec![(2500, 2600), (3000, 3100)], 9); // overlaps 1000-3000 by 800bp
        cp2.seq = rand_seq(2000, 0xB2);
        let kept = prune_same_locus(vec![cp1, cp2], &DetectParams::default());
        assert_eq!(kept.len(), 2, "disjoint-junction tandem copies preserved despite span overlap");
    }

    #[test]
    fn prune_same_locus_keeps_cross_chrom_copies() {
        // RABL2A/B-class: two paralog copies on DIFFERENT chromosomes (same-ish coords). Cross-chrom is
        // ALWAYS distinct loci → both survive (the basis for cross-chrom family capture).
        let mut a = rep_s(15131653, 15147533, vec![(15132000, 15133000)], 8); // strand '-'
        a.chrom = "NC_073235.2".into();
        a.strand = '+';
        let mut b = rep_s(48818440, 48832011, vec![(48819000, 48820000)], 9);
        b.chrom = "NC_086018.1".into();
        b.strand = '-'; // opposite strand, different chrom — still a genuine paralog pair
        let kept = prune_same_locus(vec![a, b], &DetectParams::default());
        assert_eq!(kept.len(), 2, "cross-chrom (even opposite-strand) copies are distinct loci → both kept");
    }

    #[test]
    fn prune_same_locus_drops_same_locus_antisense() {
        // a + gene and a − gene with OVERLAPPING spans at the same locus = inverted-repeat artifact, NOT
        // two copies → de-dup to one. (Distinct-loci opposite-strand copies do not overlap, so survive.)
        let mut plus = rep_s(1000, 11000, vec![(2000, 3000), (4000, 5000)], 15);
        plus.strand = '+';
        let mut minus = rep_s(1100, 10900, vec![(6000, 7000), (8000, 9000)], 5); // overlaps plus, opposite strand
        minus.strand = '-';
        let kept = prune_same_locus(vec![plus, minus], &DetectParams::default());
        assert_eq!(kept.len(), 1, "same-locus antisense overlap collapses to one");
        assert_eq!(kept[0].strand, '+', "kept the more-structured/more-read (15-read 4-exon) copy");
    }

    #[test]
    fn prune_same_locus_span_aware_merges_disjoint_junction_isoforms() {
        // Two same-strand isoforms with overlapping spans, disjoint junctions, and identical sequences
        // (100% POA core). Span-aware prune should collapse them to one locus.
        let iso_a = rep_s(1000, 3000, vec![(1500, 1600), (2000, 2100)], 15);
        let iso_b = rep_s(2200, 4200, vec![(2500, 2600), (3000, 3100)], 9); // overlaps, same strand, disjoint junctions
        let kept = prune_same_locus(vec![iso_a, iso_b], &DetectParams::default());
        assert_eq!(
            kept.len(),
            1,
            "span-aware prune merges disjoint-junction isoforms with high span homology"
        );
        assert_eq!(kept[0].introns.len(), 2, "kept the more-structured representative");
    }

    #[test]
    fn prune_same_locus_span_aware_preserves_adjacent_paralogs() {
        // Same geometry as above but with low sequence homology. The conservative core threshold must NOT
        // merge the two paralogs.
        let mut para_a = rep_s(1000, 3000, vec![(1500, 1600), (2000, 2100)], 15);
        para_a.seq = rand_seq(2000, 0xA1); // random, not all-A
        let mut para_b = rep_s(2200, 4200, vec![(2500, 2600), (3000, 3100)], 9);
        para_b.seq = rand_seq(2000, 0xB2); // random, disjoint from para_a
        let kept = prune_same_locus(vec![para_a, para_b], &DetectParams::default());
        assert_eq!(
            kept.len(),
            2,
            "adjacent paralogs with low span homology are preserved as two loci"
        );
    }

    #[test]
    fn colocated_families_splits_by_strand() {
        use super::super::family_split::{CommunityStats, FamilyClass, SplitFamily};
        // 4 disjoint-locus copies on c1: two '+' and two '-' (distinct junctions, no containment), all in
        // ONE conflict family. colocated_families must split them into TWO families by strand (the
        // antisense fix): a '+' 2-copy family and a '-' 2-copy family.
        let plus = |s: u64, e: u64, j: (u64, u64)| {
            let mut t = rep_s(s, e, vec![j], 5);
            t.strand = '+';
            t
        };
        let minus = |s: u64, e: u64, j: (u64, u64)| rep_s(s, e, vec![j], 5); // rep_s defaults to '-'
        let reps = vec![
            plus(1000, 2000, (1200, 1300)),
            plus(5000, 6000, (5200, 5300)),
            minus(3000, 4000, (3200, 3300)),
            minus(7000, 8000, (7200, 7300)),
        ];
        let fam = SplitFamily {
            members: vec![0, 1, 2, 3],
            stats: CommunityStats { n: 4, n_edges: 0, density: 1.0, avg_core_recip: 0.0, n_articulation: 0 },
            class: FamilyClass::Family,
        };
        let colo = colocated_families(&reps, &[fam], 5_000_000, 2, &DetectParams::default());
        assert_eq!(colo.len(), 2, "the + and - copies form two separate same-strand families");
        let mut sizes: Vec<usize> = colo.iter().map(|c| c.copies.len()).collect();
        sizes.sort();
        assert_eq!(sizes, vec![2, 2], "each strand family keeps its 2 disjoint copies");
        // every family is single-strand
        for c in &colo {
            let strands: std::collections::HashSet<char> = c.copies.iter().map(|x| x.strand).collect();
            assert_eq!(strands.len(), 1, "a co-located family is single-strand");
        }
    }

    #[test]
    fn families_from_conflict_graph_keeps_same_chrom_inverted_duplication() {
        // The O1-catalog regression guard for same-chromosome families. An INVERTED duplication —
        // two homologous copies on the SAME chromosome at DISJOINT loci but OPPOSITE strands — de-ties
        // reads (one conflict edge) and is a genuine multi-copy family. `families_from_conflict_graph`
        // (the O1 emission) must KEEP it. This is exactly the case `colocated_families` drops via its
        // (chrom,strand) partition (see `colocated_families_splits_by_strand`); routing the O1 same-chrom
        // catalog through `colocated_families` is what lost inverted-dup + distant same-chrom paralogs.
        let mut plus = rep_s(1_000, 3_000, vec![(1_500, 1_600)], 8);
        plus.chrom = "c1".into();
        plus.strand = '+';
        let mut minus = rep_s(5_000, 7_000, vec![(5_500, 5_600)], 6); // disjoint span, opposite strand, same chrom
        minus.chrom = "c1".into();
        minus.strand = '-';
        let reps = vec![plus, minus];
        let edges = vec![(0usize, 1usize, 5usize)]; // reads de-tie the two copies
        let fams = families_from_conflict_graph(&reps, &edges, 2, &DenovoConfig::default());
        assert_eq!(fams.len(), 1, "the inverted-duplication pair forms one O1 family");
        assert_eq!(fams[0].len(), 2, "both copies retained (no strand split, no win cap)");
        let chroms: std::collections::BTreeSet<&str> = fams[0].iter().map(|c| c.chrom.as_str()).collect();
        assert_eq!(chroms.len(), 1, "it is a SAME-chromosome family");
        let strands: std::collections::BTreeSet<char> = fams[0].iter().map(|c| c.strand).collect();
        assert_eq!(strands.len(), 2, "the two copies are on opposite strands (inverted duplication)");
    }

    #[test]
    fn build_read_placements_multimapper_forms_conflict_family() {
        let reps = vec![rep("c1", 0, 200), rep("c1", 1000, 1200)];
        let bam = vec![
            bam_read("c1", 0, 200, "read_X", 0.010, false),
            bam_read("c1", 1000, 1200, "read_X", 0.012, true),
        ];
        let placements = build_read_placements(&bam, &reps);
        let edges = super::super::read_conflict::conflict_edges(
            2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        assert!(!edges.is_empty(), "de-tied cross-locus read must produce a conflict edge");
        assert_eq!(conflict_families(2, &edges), vec![vec![0, 1]]);
    }

    #[test]
    fn build_read_placements_domain_sharer_is_no_family() {
        let reps = vec![rep("c1", 0, 200), rep("c1", 1000, 1200)];
        let bam = vec![bam_read("c1", 0, 200, "read_Y", 0.010, false)];
        let placements = build_read_placements(&bam, &reps);
        let edges = super::super::read_conflict::conflict_edges(
            2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        assert!(edges.is_empty());
        assert!(conflict_families(2, &edges).is_empty());
    }

    #[test]
    fn build_read_placements_nested_loci_unique_read_no_conflict() {
        let reps = vec![rep("c1", 0, 1000), rep("c1", 400, 600)];
        let bam = vec![bam_read("c1", 450, 550, "unique_read", 0.010, false)];
        let placements = build_read_placements(&bam, &reps);
        let total: usize = placements.iter().map(|p| p.len()).sum();
        assert_eq!(total, 1, "one record -> one best-overlap placement -> no cross-locus pair");
        let edges = super::super::read_conflict::conflict_edges(
            2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        assert!(edges.is_empty());
    }

    #[test]
    fn build_read_placements_oversplit_fragments_form_no_family() {
        // Over-split robustness (airtight property): two near-identical co-positioned loci —
        // one gene the locus collapse left as two fragments (the PRNP case). Because each
        // record attributes to its single best-overlap locus, co-positioned fragments share
        // no conflicting read, so they cannot form a conflict family — even with many
        // error-floor-de reads. A sequence-similarity definition WOULD group them (~42% of
        // GGO similarity "families" were one over-split locus); the conflict definition needs
        // no over-split guard. See bench/family_definition_formal.md.
        let reps = vec![rep("c1", 1000, 6000), rep("c1", 1005, 6000)];
        let bam: Vec<BamRead> = (0..20)
            .map(|i| bam_read("c1", 1000, 6000, &format!("frag_read_{i}"), 0.004, false))
            .collect();
        let placements = build_read_placements(&bam, &reps);
        assert!(placements.iter().all(|p| p.len() == 1),
            "each record attributes to one best-overlap locus -> no cross-locus pair");
        let edges = super::super::read_conflict::conflict_edges(
            2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 3, sig: None });
        assert!(edges.is_empty(), "over-split fragments must not form a conflict edge");
        assert!(conflict_families(2, &edges).is_empty(),
            "an over-split single locus is not a multi-copy family");
    }

    #[test]
    fn build_read_placements_diverged_secondary_is_no_conflict() {
        // read fits copy 0 (de 0.001) far better than copy 1 (de 0.020): resolvable, no edge.
        let reps = vec![rep("c1", 0, 200), rep("c1", 1000, 1200)];
        let bam = vec![
            bam_read("c1", 0, 200, "read_Z", 0.001, false),
            bam_read("c1", 1000, 1200, "read_Z", 0.020, true),
        ];
        let placements = build_read_placements(&bam, &reps);
        let edges = super::super::read_conflict::conflict_edges(
            2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        assert!(edges.is_empty(), "diverged secondary (large de gap) is resolvable, not a conflict");
    }

    #[test]
    fn build_read_placements_excludes_supplementary() {
        // a supplementary (chimeric) record on the second locus must NOT create a conflict edge.
        let reps = vec![rep("c1", 0, 200), rep("c1", 1000, 1200)];
        let mut supp = bam_read("c1", 1000, 1200, "read_S", 0.010, true);
        supp.is_supplementary = true;
        let bam = vec![bam_read("c1", 0, 200, "read_S", 0.010, false), supp];
        let placements = build_read_placements(&bam, &reps);
        let total: usize = placements.iter().map(|p| p.len()).sum();
        assert_eq!(total, 1, "supplementary record excluded -> only the primary placement remains");
        let edges = super::super::read_conflict::conflict_edges(
            2, &placements, &ConflictParams { delta: 0.005, de_max: 0.05, min_reads: 1, sig: None });
        assert!(edges.is_empty());
    }

    #[test]
    fn denovoconfig_default_conflict_params_are_sane() {
        let cfg = DenovoConfig::default();
        assert!((cfg.conflict.delta - 0.005).abs() < 1e-9);
        assert!((cfg.conflict.de_max - 0.05).abs() < 1e-9);
        assert_eq!(cfg.conflict.min_reads, 3);
    }

    #[test]
    fn conflict_to_split_families_produces_one_family_per_component() {
        // Two disjoint conflict components: {0,1} linked by 5 reads, {2,3} by 3 reads.
        let families = vec![vec![0usize, 1], vec![2usize, 3]];
        let c_edges = vec![(0usize, 1usize, 5usize), (2, 3, 3)];
        let p = SplitParams::default();
        let split = conflict_to_split_families(&families, &c_edges, &p);
        assert_eq!(split.len(), 2);
        assert!(split.iter().all(|sf| sf.class == FamilyClass::Family));
        // members are sorted ascending within each family.
        assert_eq!(split[0].members, vec![0, 1]);
        assert_eq!(split[1].members, vec![2, 3]);
        // density of a 2-node clique (1 edge / 1 possible) = 1.0.
        for sf in &split {
            assert!((sf.stats.density - 1.0).abs() < 1e-9, "2-node clique density must be 1.0");
        }
    }

    #[test]
    fn homology_edges_links_two_similar_reps_not_divergent() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            eprintln!("minimap2 absent; skipping"); return;
        }
        // rep 0 and rep 1: same 900bp sequence with ~5% mismatch (an ancient paralog); rep 2: random (unrelated).
        let base = rand_seq(900, 0x11);
        let mut para = base.clone();
        for k in (0..para.len()).step_by(20) { para[k] = b"ACGT"[(para[k] as usize + 1) % 4]; } // ~5% divergence
        let reps = vec![
            DenovoTranscript { tid: "r0".into(), chrom: "c1".into(), start: 0, end: 900, n_reads: 10, strand: '+', introns: vec![], seq: base },
            DenovoTranscript { tid: "r1".into(), chrom: "c1".into(), start: 5000, end: 5900, n_reads: 8, strand: '+', introns: vec![], seq: para },
            DenovoTranscript { tid: "r2".into(), chrom: "c1".into(), start: 9000, end: 9900, n_reads: 5, strand: '+', introns: vec![], seq: rand_seq(900, 0x99) },
        ];
        let params = RefineParams::default(); // min_identity 0.80, sensitive_identity 0.60, min_coverage 0.50
        let edges = homology_edges_all_reps(&reps, &params).unwrap();
        assert!(edges.contains(&(0, 1)), "the two paralog reps must be E_r-linked, got {:?}", edges);
        assert!(!edges.contains(&(0, 2)) && !edges.contains(&(1, 2)), "the unrelated rep must not link");
    }

    #[test]
    fn homology_catalog_groups_fixture_family() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
        let fams = detect_homology_catalog_genome_wide(
            "tests/fixtures/same_chrom_supplement/reads.bam",
            "tests/fixtures/same_chrom_supplement/genome.fa",
            2, 2, &DenovoConfig::default(), &RefineParams::default(), 0.20,
        ).unwrap();
        // the fixture's two homologous loci (c1:A + c2:X) must land in one family of >= 2 distinct loci.
        assert!(fams.iter().any(|f| f.len() >= 2), "expected a >=2-copy homology family, got {:?}", fams.iter().map(|f| f.len()).collect::<Vec<_>>());
    }

}
