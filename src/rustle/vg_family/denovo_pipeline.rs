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
use super::copy_assign::{assign_read_editing, Assignment, AssignParams, AssignStatus, CopyProfile, ReadFeatures};
use super::copy_assign_pipeline::{
    assign_family_detailed, assign_family_detailed_pruned, best_overlap_copy, build_family_profiles,
    copy_boundaries, detect_editing_columns, freeze_merge, gen2off, read_ref_end, FamilyProfiles,
};
use super::family_graph::contiguous_core_coverage_bounded;
use super::copy_split::{
    split_locus_copies, discover_locus_psvs, AlignedRead, CollapsedCandidate, CopyIsoform,
};
use super::denovo_assemble::{
    assemble_gate, pass1_skeletons, pass1_skeletons_robust, primary_reads_from_bam, reads_in_region,
    BamRead, GateParams, PrimaryRead, GATE_MIN_READS, PASS1_MIN_READS,
};
use super::family_detect::{
    collapse_loci_span_aware, collapse_loci_span_aware_with_totals, detect_edges, detect_edges_reporting,
    DenovoTranscript, DetectParams,
};
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
    /// VG re-align supplement (Task 5, opt-in): default false = OFF, byte-identical. When true,
    /// `detect_and_assign` runs `vg_realign::run_family_realign` over each co-located family's reads and
    /// stores the resulting `Vec<RealignRecord>` on `FamilyAssignment::realign_records`. REPORT-ONLY —
    /// does not alter the copy set, the PSV/junction assignment, or any other emitted field.
    pub vg_realign: bool,
    /// E_r homology-primary family MEMBERSHIP (opt-in). Conflict/PSV/χ(H) remain within-family. Enlarges
    /// the copy set ⟹ stricter Bonferroni α/(K−1) ⟹ assignments shift. Requires minimap2.
    pub homology_primary: bool,
    /// Drop single-exon reps that engulf >= `READTHROUGH_MIN_DISTINCT` distinct junctions: unspliced
    /// pre-mRNA, never a copy. Default ON — see `is_unspliced_readthrough` for the validation.
    pub filter_readthrough: bool,
    /// Admit a COLLAPSED single-rep locus as a multi-copy family: `n_copies = χ(H)`, reads certified Tied, no
    /// per-copy consensus materialised. **Default OFF**: the ambiguity instrument detects unresolvable
    /// PARALOGY, not collapse — it fires on EEF1A1 (pseudogenes on other chromosomes) with χ(H) = 7. See the
    /// `collapse_gate` module header.
    pub collapse_gate: bool,
    /// Background per-read ambiguity rate for the collapse test. Must be GENOME-WIDE, never region-local: in the
    /// DAZ window every read outside DAZ1's span is DAZ2's and ambiguous, so a local background would be ~0.95.
    /// `None` ⇒ the gate abstains. Default = `GENOME_WIDE_EPS_AMB` measured on `GGO_mm.bam`.
    pub eps_amb: Option<f64>,
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
            vg_realign: false,
            homology_primary: false,
            filter_readthrough: true,
            collapse_gate: false,
            eps_amb: Some(crate::vg_family::collapse_gate::GENOME_WIDE_EPS_AMB),
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
    // Same `k` as every other path: a locus must have ONE canonical extent across objectives. The three
    // genome-wide O1 catalogs already trim at `cfg.min_terminal_support`; this path silently used k = 1.
    let skeletons = pass1_skeletons_robust(reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    let mut transcripts = assemble_gate(&skeletons, genome, &cfg.gate);
    if cfg.filter_readthrough {
        let support = read_junction_support(reads);
        retain_non_readthrough(&mut transcripts, &support, "detect_families");
        retain_non_mischain(&mut transcripts, &support, "detect_families");
    }
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
        // Partition by CHROM only — NOT by strand. An inverted duplicate (MAGEA4 `+` / MAGEA10 `-`,
        // RFPL2 `-` / RFPL3 `+`) is a genuine paralog pair at two disjoint loci on opposite strands. Keying
        // on `(chrom, strand)` split every such pair into two singletons, which `min_copies` then dropped,
        // so O2 could never assign an inverted duplicate even when the homology oracle found the edge.
        //
        // The strand split was really guarding against SAME-LOCUS antisense (a `+` gene and a `-` gene whose
        // spans overlap — an inverted-repeat artifact, not two copies). `prune_same_locus` clause (c) removes
        // exactly that, and only that: distinct-loci opposite-strand copies do not overlap, so it never fires
        // on a real inverted duplication. Running it over the whole chromosome group — rather than within one
        // strand, where it could never see an antisense pair at all — gives the same protection without
        // discarding the paralogs. Mixed-strand families are safe downstream: each copy carries its own
        // strand, `copy_assign_pipeline` reverse-complements read bases per copy, and two inverted-duplicate
        // mRNAs are both in transcription orientation, so they align forward to each other.
        let mut by_key: std::collections::BTreeMap<&str, Vec<usize>> = std::collections::BTreeMap::new();
        for &m in &fam.members {
            by_key.entry(reps[m].chrom.as_str()).or_default().push(m);
        }
        for (chrom, mut idxs) in by_key {
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
/// and the unique-mapper agreement (the accuracy proxy on real data — `copy_assign.py`).
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
    /// Per-copy junction-only support (parallel to `copy_tids`): reads a copy-specific junction pinned to each
    /// copy where PSVs could not. A copy with `>= GATE_MIN_READS` here is identifiable by splice structure alone,
    /// so its existence is invariant to the arbitrary primary/secondary label even without unique mappers (the
    /// DAZ2 case). Feeds the `junction_invariant` certificate column. See `add_junction_support`.
    pub copy_junction_support: Vec<usize>,
    /// unique-mapper agreement: assigned + uniquely mapped (mapq > 0), and of those how many agree with where the
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
    /// Genomic span `(chrom, start, end)` of each copy, parallel to `copy_tids`. Emitted so a catalog can be
    /// audited for the same-locus artifact: two copies of ONE family whose spans OVERLAP are not two loci,
    /// they are one locus admitted twice (nested/near-duplicate de-novo transcripts). Such a family reports
    /// `min_p == 1` for every read — the reads look unassignable when in fact the copy set is malformed.
    /// The check is structural and needs no annotation. See `copies_overlap` and `bench/artifact_audit.py`.
    pub copy_spans: Vec<(String, u64, u64)>,
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
    /// Task H2 (VG-harmony): per-copy copy-specific junction offsets, parallel to `copy_psv_alleles`
    /// (`copy_junctions[k]` aligns with `copy_psv_alleles[k]`) — threads O3's junction evidence
    /// into the same per-copy frame the EM engine (`em_assign_family`) consumes. Empty when the
    /// underlying `FamilyDetail` carries no junction data.
    pub copy_junctions: Vec<Vec<i64>>,
    /// per-read intron-boundary offsets in the assigned copy's spliced space, parallel to
    /// `read_psv_obs`/`assignments` (`read_junctions[i]` aligns with `read_psv_obs[i]`/`assignments[i]`).
    pub read_junctions: Vec<Vec<i64>>,
    /// Task 5 (opt-in via `DenovoConfig::vg_realign`, default OFF): the VG re-align supplement's per-read
    /// decisions for this family (`vg_realign::run_family_realign`) — empty unless `vg_realign` is set.
    /// REPORT-ONLY: not consumed by `assignments`/`copy_tids`/anything else above.
    pub realign_records: Vec<crate::vg_family::vg_realign::RealignRecord>,
}

impl FamilyAssignment {
    /// All-zero / all-empty. Only for constructing a family that never went through the assignment pipeline
    /// (see `gated_family`); a normal family is built field-by-field and must not silently default anything.
    pub fn empty() -> Self {
        FamilyAssignment {
            family_id: String::new(),
            chrom: String::new(),
            n_copies: 0,
            n_reads: 0,
            psv_cols: 0,
            resolvable_psv: 0,
            assigned_psv: 0,
            resolvable_j: 0,
            assigned_j: 0,
            junction_only: 0,
            copy_junction_support: Vec::new(),
            uniq: 0,
            uniq_agree: 0,
            collapsed_copies: 0,
            rescued_copies: 0,
            assignments: Vec::new(),
            copy_tids: Vec::new(),
            copy_spans: Vec::new(),
            copy_abundance: Vec::new(),
            copy_abundance_ci: Vec::new(),
            mosaic_reads: 0,
            conversions: Vec::new(),
            copy_conversions: Vec::new(),
            psv_col_pos: Vec::new(),
            copy_psv_alleles: Vec::new(),
            read_psv_obs: Vec::new(),
            copy_junctions: Vec::new(),
            read_junctions: Vec::new(),
            realign_records: Vec::new(),
        }
    }
}

/// END-TO-END pipeline: detect families, then for each co-located family assign every read overlapping it to
/// a copy (PSV + copy-specific-junction likelihood, two-pass). `bam_reads` carry chrom (for region
/// filtering) and mapq (for the unique-mapper agreement). The runnable detection + per-read copy-assignment pipeline.
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
    to_split_families(families, &float_edges, p)
}

/// Core components->`SplitFamily` conversion shared by every membership oracle (read-conflict today,
/// homology-primary E_r later): edges are already `f64` weights (read counts for the conflict path, a
/// uniform `1.0` for the homology path). Class follows the same size+density rule regardless of the
/// oracle that produced `families`/`edges`.
fn to_split_families(
    families: &[Vec<usize>],
    edges: &[(usize, usize, f64)],
    p: &SplitParams,
) -> Vec<SplitFamily> {
    let mut out: Vec<SplitFamily> = families
        .iter()
        .map(|members| {
            let mut m = members.clone();
            m.sort_unstable();
            let stats = community_stats(&m, edges);
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

/// Pairs of copy indices whose genomic spans OVERLAP, with the RECIPROCAL overlap fraction.
///
/// Two copies of one family are, by definition, two distinct loci — distinct loci occupy disjoint genomic
/// intervals. Any intersection therefore signals a defect, and the reciprocal fraction
/// `overlap / max(len_i, len_j)` says which one:
///
/// * **≈ 1.0 — one locus admitted twice.** Two de-novo transcripts on top of each other, differing only by
///   boundary wobble (observed: `164381222-164384848` vs `164381237-164384845`, 15 bp offset). They are
///   sequence-identical, so every read scores `min_p == 1` against the pair, the family abstains wholesale,
///   and its reads masquerade as the K=0 identifiability wall. `collapse_loci_groups` misses these because
///   it unions only transcripts sharing an EXACT intron `(chrom, donor, acceptor)` and never consults
///   positional overlap.
/// * **≪ 1.0 — containment.** A long readthrough transcript enclosing a short one (observed: a 188 kb span
///   containing a 3.6 kb one). Merging these would let the readthrough absorb genuinely distinct tandem
///   copies, so they must NOT be merged — it is a separate defect in transcript construction.
///
/// Purely structural: no annotation, no reference truth. Returns `(i, j, reciprocal)` with `i < j`.
pub fn copies_overlap(spans: &[(String, u64, u64)]) -> Vec<(usize, usize, f64)> {
    let mut out = Vec::new();
    for i in 0..spans.len() {
        for j in (i + 1)..spans.len() {
            let (ca, sa, ea) = (&spans[i].0, spans[i].1, spans[i].2);
            let (cb, sb, eb) = (&spans[j].0, spans[j].1, spans[j].2);
            if ca == cb && sa.max(sb) < ea.min(eb) {
                let overlap = (ea.min(eb) - sa.max(sb)) as f64;
                let longest = (ea - sa).max(eb - sb) as f64;
                out.push((i, j, if longest > 0.0 { overlap / longest } else { 1.0 }));
            }
        }
    }
    out
}

use crate::vg_family::collapse_gate::{collapse_verdict, Ambiguity, CollapseVerdict};

/// Ambiguously-placed primary reads over a rep's span.
///
/// A read is AMBIGUOUS iff `mapq == 0`: the aligner found no reason to prefer this placement over another.
/// Supplementary records are excluded for the same reason they are excluded from conflict edges — a chimeric
/// segment is adjacency, not ambiguity. Secondary records never enter `bam_reads` as separate entries here.
fn locus_ambiguity(rep: &DenovoTranscript, bam_reads: &[BamRead]) -> Ambiguity {
    let (mut n, mut k) = (0usize, 0usize);
    for br in bam_reads {
        if br.is_supplementary || br.chrom != rep.chrom {
            continue;
        }
        if br.read.ref_start < rep.end && read_ref_end(&br.read) > rep.start {
            n += 1;
            if br.mapq == 0 {
                k += 1;
            }
        }
    }
    Ambiguity { n, k }
}

/// A family whose copies were never assembled: `n_copies` is χ(H), every read is certified `Tied`, and NO
/// per-copy consensus exists — so no assignment can be wrong. `min_p_value = 1.0` is the honest certificate:
/// no distinguishing column was available to this read, because no copy sequences were materialised.
fn gated_family(rep: &DenovoTranscript, bam_reads: &[BamRead], chi: usize, family_id: String) -> FamilyAssignment {
    let assignments: Vec<(usize, Assignment)> = bam_reads
        .iter()
        .enumerate()
        .filter(|(_, br)| {
            !br.is_supplementary
                && br.chrom == rep.chrom
                && br.read.ref_start < rep.end
                && read_ref_end(&br.read) > rep.start
        })
        .map(|(i, _)| {
            (
                i,
                Assignment {
                    best_copy: 0,
                    log_lr_margin: 0.0,
                    n_decisive: 0,
                    resolvable: false,
                    status: AssignStatus::Tied,
                    p_value: 1.0,
                    min_p_value: 1.0,
                    discovery_coupled: false,
                    posterior: vec![1.0 / chi as f64; chi],
                },
            )
        })
        .collect();
    FamilyAssignment {
        family_id,
        chrom: rep.chrom.clone(),
        n_copies: chi,
        n_reads: assignments.len(),
        collapsed_copies: chi, // records HOW the count was obtained: never mistake this for an assembled family
        assignments,
        copy_tids: vec![rep.tid.clone()], // the one assembled rep; the other copies have no sequence
        copy_spans: vec![(rep.chrom.clone(), rep.start, rep.end)],
        ..FamilyAssignment::empty()
    }
}

/// Reads required before a splice junction counts as real rather than alignment noise.
pub const READTHROUGH_MIN_SUPPORT: usize = 2;
/// Distinct junctions inside a single-exon transcript's span that mark it as unspliced pre-mRNA.
///
/// Read off the data, not tuned: across 30 regions every single-exon de-novo transcript (n = 15, including
/// the GSTM, RFPL, DAZ and BPY2 readthroughs) contains **≥ 14** distinct junctions, while **no** expressed
/// annotated intronless gene out of 260 exceeds **4** — including the EEF1A1 retrocopy `LOC109023808`, whose
/// spliced parent contributes 3516 reads that cross-map onto it. `bench/READTHROUGH_RULE_VALIDATION.md`.
pub const READTHROUGH_MIN_DISTINCT: usize = 5;

/// An intron larger than this (bp) is "giant" and is checked for mis-chain support. Set above the largest intron
/// seen in a real recovered paralog on this substrate (POTE's 48 kb), so genuine large-gene introns are never
/// examined. The large-gene mis-chains (`bench/GW_CATALOG_FP_AUDIT.md`) carry giant introns of 48 kb–1.2 Mb.
pub const MISCHAIN_GIANT_INTRON_BP: u64 = 50_000;

/// A giant intron carried by fewer than this many reads (that exact junction) is a spurious splice below the
/// locus gate — it cannot anchor a real transcript. Equal to the locus gate (`GATE_MIN_READS`), so this only
/// removes giant introns that are already sub-threshold noise; it deliberately does NOT try to catch
/// well-supported gene-splits (PBX1's 115 kb spurious intron has 6 reads, indistinguishable from a real
/// low-expression intron without homology — that is `--refine`'s job).
pub const MISCHAIN_MIN_JUNCTION_READS: usize = GATE_MIN_READS as usize;

/// Splice junctions observed in the PRIMARY reads, with the number of primary reads supporting each.
/// Key is `(chrom, donor, acceptor)` in genomic coordinates.
///
/// Primary-only, i.e. `-F 2308`. `PrimaryRead` is secondary-free by construction
/// (`primary_read_from_record` returns `None` for secondary and supplementary records) and its `introns` are
/// already the CIGAR-`N` gaps. The predecessor took `&[BamRead]` and skipped only `is_supplementary`, so it
/// **counted secondary alignments** -- on this substrate (`minimap2 -N 50`, ~63% secondary) that inflated the
/// statistic from **56 to 154** distinct junctions at the DAZ readthrough span. The rule
/// ([`is_unspliced_readthrough`]) was validated on primaries alone, so its input must be too.
/// Accumulate the per-copy JUNCTION-ONLY support that certifies a copy as invariant to the arbitrary
/// primary/secondary label. A read is junction-only support for its assigned copy when a copy-specific junction
/// made it decisive (`combined_decisive`) but PSVs alone could not (`!psv_decisive`) — it carries a junction
/// distinctive to that copy, identifying it by splice structure regardless of which locus holds its primary
/// alignment. Bucketed per `best_copy`; an out-of-range copy index is ignored. Mirrors the family-level
/// `junction_only` tally, per copy — the DAZ2-rescue signal.
pub fn add_junction_support(support: &mut [usize], best_copy: usize, combined_decisive: bool, psv_decisive: bool) {
    if combined_decisive && !psv_decisive {
        if let Some(s) = support.get_mut(best_copy) {
            *s += 1;
        }
    }
}

pub fn read_junction_support(reads: &[PrimaryRead]) -> std::collections::HashMap<(String, u64, u64), usize> {
    let mut sup = std::collections::HashMap::new();
    for pr in reads {
        for &(d, a) in &pr.introns {
            *sup.entry((pr.chrom.clone(), d, a)).or_insert(0) += 1;
        }
    }
    sup
}

/// Drop unspliced readthrough transcripts in place, logging each. Shared by every path that assembles de-novo
/// transcripts, so O1 (the family catalogs) and O2 (assignment) agree about what a locus IS.
///
/// Must run on TRANSCRIPTS, before `collapse_loci_span_aware`: a readthrough is the longest object in its span,
/// so after collapse it becomes the representative of every locus it covers, and dropping that rep deletes the
/// real copies with it (the DAZ 298 kb span absorbed DAZ2 into DAZ1's group).
pub fn retain_non_readthrough(
    transcripts: &mut Vec<DenovoTranscript>,
    support: &std::collections::HashMap<(String, u64, u64), usize>,
    tag: &str,
) {
    let mut dropped = Vec::new();
    transcripts.retain(|t| {
        let rt = is_unspliced_readthrough(t, support, READTHROUGH_MIN_SUPPORT, READTHROUGH_MIN_DISTINCT);
        if rt {
            dropped.push(format!("{}:{}-{} ({} reads)", t.chrom, t.start, t.end, t.n_reads));
        }
        !rt
    });
    if !dropped.is_empty() {
        eprintln!(
            "[{tag}] readthrough filter: dropped {} single-exon transcript(s) engulfing >= {} distinct \
             junctions (unspliced pre-mRNA, not copies) -> {} transcripts",
            dropped.len(),
            READTHROUGH_MIN_DISTINCT,
            transcripts.len()
        );
        for d in &dropped {
            eprintln!("[{tag}]   readthrough {d}");
        }
    }
}

/// Is this transcript the UNSPLICED form of a locus (intronic pileup / readthrough) rather than an mRNA?
///
/// True iff it is **single-exon** and at least [`READTHROUGH_MIN_DISTINCT`] distinct junctions, each with
/// [`READTHROUGH_MIN_SUPPORT`] supporting reads, lie **entirely inside** its span.
///
/// A readthrough engulfs whole gene structures, so it contains many distinct junctions. Three properties make
/// the rule safe, and each was validated against a rule that failed without it:
/// * **Distinct**, not total — TSPYL1 has 51 junction *observations* inside its span but only 4 distinct ones.
/// * **Read-level**, not transcript-level — a rule phrased over assembled spliced transcripts misses the RFPL
///   readthrough entirely, because that window is too sparsely expressed to assemble one.
/// * **Entirely inside** — a gene nested in another gene's intron sees the host's junctions *flanking* it,
///   never within it. A retrocopy is likewise safe: it has no introns, so its spliced parent's cross-mapping
///   reads align contiguously and cannot deposit a junction inside it.
///
/// Rejected alternatives (all measured, `bench/READTHROUGH_RULE_VALIDATION.md`, `bench/YAG_CHECK.md`): any
/// single read junction inside the span (drops 21% of real intronless genes); span / longest-contained-read
/// ratio (no separation); exonic coverage breadth (the GSTM readthrough is better covered than TSPYL1).
pub fn is_unspliced_readthrough(
    rep: &DenovoTranscript,
    junction_support: &std::collections::HashMap<(String, u64, u64), usize>,
    min_support: usize,
    min_distinct: usize,
) -> bool {
    if !rep.introns.is_empty() {
        return false; // it splices; it is a transcript model, not an unspliced span
    }
    let engulfed = junction_support
        .iter()
        .filter(|((chrom, d, a), &n)| {
            n >= min_support && chrom == &rep.chrom && rep.start <= *d && *a <= rep.end
        })
        .count();
    engulfed >= min_distinct
}

/// Is this transcript held together by a SUB-GATE giant intron — a spurious splice minimap2 created by
/// mis-aligning a read's ends across a large gene's giant intron?
///
/// True iff any of the transcript's OWN introns is GIANT (`> giant_bp`) yet carried by fewer than
/// `min_junction_reads` primary reads carrying that exact junction. A junction below the locus gate is noise, so
/// a giant intron below it cannot anchor a real transcript. This removes the clearly-unsupported mis-chains
/// (467 genome-wide) that would otherwise inflate spans and seed false loci.
///
/// **Scope limit (measured, `bench/GW_CATALOG_FP_AUDIT.md`):** this does NOT catch a well-supported gene-split.
/// PBX1's 115 kb spurious intron is carried by 6 reads — above the gate — so it is intrinsically indistinguishable
/// from a real low-expression large-gene intron by any within-transcript signal. Separating those two requires
/// HOMOLOGY context (the copy shares no sequence with its supposed paralog), which is the `--refine` gate's job,
/// not the assembler's. Distinct from [`is_unspliced_readthrough`] (single-exon only).
pub fn is_giant_intron_mischain(
    rep: &DenovoTranscript,
    junction_support: &std::collections::HashMap<(String, u64, u64), usize>,
    giant_bp: u64,
    min_junction_reads: usize,
) -> bool {
    rep.introns.iter().any(|&(d, a)| {
        a.saturating_sub(d) > giant_bp
            && junction_support.get(&(rep.chrom.clone(), d, a)).copied().unwrap_or(0) < min_junction_reads
    })
}

/// Drop de-novo transcripts that are large-gene mis-chains (see [`is_giant_intron_mischain`]). Mirrors
/// [`retain_non_readthrough`]; reuses the same primary-read junction-support map. Emits a per-drop diagnostic.
pub fn retain_non_mischain(
    transcripts: &mut Vec<DenovoTranscript>,
    support: &std::collections::HashMap<(String, u64, u64), usize>,
    tag: &str,
) {
    let mut dropped = Vec::new();
    transcripts.retain(|t| {
        let mc = is_giant_intron_mischain(t, support, MISCHAIN_GIANT_INTRON_BP, MISCHAIN_MIN_JUNCTION_READS);
        if mc {
            dropped.push(format!("{}:{}-{} ({} reads)", t.chrom, t.start, t.end, t.n_reads));
        }
        !mc
    });
    if !dropped.is_empty() {
        eprintln!(
            "[{tag}] mis-chain filter: dropped {} transcript(s) with a giant intron (> {} bp) supported by < {} \
             reads (spurious splice across a large gene, not a copy) -> {} transcripts",
            dropped.len(),
            MISCHAIN_GIANT_INTRON_BP,
            MISCHAIN_MIN_JUNCTION_READS,
            transcripts.len()
        );
        for d in &dropped {
            eprintln!("[{tag}]   mis-chain {d}");
        }
    }
}

/// How two copies in an emitted catalog can wrongly share genomic sequence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OverlapKind {
    /// Same family, reciprocal overlap ≈ 1: one locus admitted twice (boundary wobble). Every read then
    /// scores `min_p == 1` against the pair, so the family abstains wholesale and its reads masquerade as
    /// the K=0 wall.
    DuplicateLocus,
    /// Same family, reciprocal overlap ≪ 1: a long readthrough transcript enclosing a shorter one. Must NOT
    /// be merged — merging would let the readthrough absorb genuinely distinct tandem copies.
    Containment,
    /// Copies of DIFFERENT families sharing sequence. A locus belongs to exactly one family, so this is
    /// always a defect: either the same family was emitted twice (overlapping input regions), or a
    /// readthrough transcript in one family spans loci that are separate copies of another. Observed at
    /// GSTM: a 30 kb single-intron transcript covering both GSTM5 and GSTM1 was admitted as a copy
    /// alongside GSTM3, while a second family held GSTM5 and GSTM1 correctly as two copies.
    SharedAcrossFamilies,
}

/// Every pair of catalog copies that wrongly shares genomic sequence, classified. Copies are
/// `(family_id, chrom, start, end)`. Sorted sweep, so it is linear in the output size, not quadratic in
/// the catalog. Structural: needs no annotation and no reference truth.
pub fn catalog_overlaps(copies: &[(String, String, u64, u64)]) -> Vec<(usize, usize, f64, OverlapKind)> {
    let mut idx: Vec<usize> = (0..copies.len()).collect();
    idx.sort_by(|&a, &b| copies[a].1.cmp(&copies[b].1).then(copies[a].2.cmp(&copies[b].2)));
    let mut out = Vec::new();
    for (pos, &i) in idx.iter().enumerate() {
        for &j in idx[pos + 1..].iter() {
            if copies[i].1 != copies[j].1 || copies[j].2 >= copies[i].3 {
                break; // sorted by start: no later copy on this contig can overlap i
            }
            let overlap = (copies[i].3.min(copies[j].3) - copies[j].2) as f64;
            let longest = (copies[i].3 - copies[i].2).max(copies[j].3 - copies[j].2) as f64;
            let recip = if longest > 0.0 { overlap / longest } else { 1.0 };
            let kind = if copies[i].0 != copies[j].0 {
                OverlapKind::SharedAcrossFamilies
            } else if recip > 0.9 {
                OverlapKind::DuplicateLocus
            } else {
                OverlapKind::Containment
            };
            out.push((i.min(j), i.max(j), recip, kind));
        }
    }
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

// ---------------------------------------------------------------------------------------------
// VG re-align END-TO-END plan, Task 3: wire `vg_realign::apply_realign`'s corrections + admissions
// into `FamilyAssignment` under `DenovoConfig::vg_realign`. Corrections + the EM abundance recompute
// are additive/must-have (`apply_realign_patch`); admission of novel-read pools into the copy set is
// best-effort/mechanism-only (`admit_novel_pools`), matching the module's own doc on `apply_realign`.
// ---------------------------------------------------------------------------------------------

/// `--em`'s CLI defaults (`copy_assign.rs`'s `em_eps`/`em_max_iter`): `detect_and_assign` has no
/// separate EM knob of its own, so the vg-realign abundance recompute reuses the same operating
/// point the `--em` mode already ships.
const VG_REALIGN_EM_EPS: f64 = 1e-6;
const VG_REALIGN_EM_MAX_ITER: usize = 500;

/// For copy `copy_idx`, the family's PSV columns in that copy's own FORWARD-GENOME coordinate frame:
/// `out[j] = Some(g)` iff copy `copy_idx` carries column `j` (i.e. isn't gapped there), at genomic
/// position `g`. Built from `FamilyProfiles::copy_gpos[copy_idx]` (a `(column, genomic_pos)` list,
/// already computed once per family by `build_family_profiles` — no re-derivation of the PSV columns
/// themselves).
fn family_col_genomic_pos(profiles: &FamilyProfiles, copy_idx: usize) -> Vec<Option<u64>> {
    let mut out = vec![None; profiles.n_cols];
    for &(col, g) in &profiles.copy_gpos[copy_idx] {
        out[col] = Some(g);
    }
    out
}

/// Invert `col_gpos` (a copy's family-PSV-column -> genomic-position map, from
/// `family_col_genomic_pos`) through `gen2off_map` (genomic position -> spliced-consensus offset,
/// `copy_assign_pipeline::gen2off`) to get the copy's family-PSV-column -> CONSENSUS-OFFSET map that
/// `vg_realign::apply_realign`/`path_obs_at` need. A column the copy doesn't carry (`None` in
/// `col_gpos`, or a genomic position `gen2off_map` doesn't cover) gets `sentinel` (one past the
/// consensus's own length) — `path_obs_at`'s exact-index search over `0..seq.len()` never matches
/// it, so it reads back as `None`, the same "gapped" observation a real gap column produces.
fn psv_positions_for(
    col_gpos: &[Option<u64>],
    gen2off_map: &std::collections::BTreeMap<u64, usize>,
    sentinel: usize,
) -> Vec<usize> {
    col_gpos
        .iter()
        .map(|g| g.and_then(|g| gen2off_map.get(&g).copied()).unwrap_or(sentinel))
        .collect()
}

/// VG re-align END-TO-END plan, Task 3 (must-have): apply `apply_realign`'s per-read CORRECTIONS
/// into `fa`.
///
/// Sets `fa.realign_records` (the decision log, same as the old report-only wiring). For each
/// `(read_idx, (new_copy, obs))` in `apply.corrected` (already in `fa.assignments`'s GLOBAL
/// `bam_reads`-index space — the caller translates `apply_realign`'s per-family-local indices before
/// calling this), finds the `fa.assignments` entry with that `read_idx` and:
///   * overwrites the PARALLEL `fa.read_psv_obs` entry (same position) to `obs`;
///   * (I3 fix) RE-DERIVES the entire `Assignment` from scratch over the corrected `obs`, via the
///     SAME gate (`copy_assign::assign_read_editing`) the one-shot pipeline uses — not just
///     `best_copy`. Before this fix, a correction only moved `best_copy`/`read_psv_obs` and left
///     `status`/`posterior`/`p_value`/`log_lr_margin`/`n_decisive` stale from the read's
///     PRE-correction call (e.g. a read reassigned off a `Tied` call still printed `status=Tied`,
///     with a stale `posterior`), making `.assignments.tsv`/`.posterior.tsv` self-inconsistent.
///
/// A corrected read with no existing `fa.assignments` entry (the one-shot gate itself couldn't
/// resolve it at all) is left alone — nothing is fabricated.
///
/// Does NOT recompute `fa.copy_abundance` — that must happen AFTER `admit_novel_pools` has finished
/// widening the copy roster (I2 fix; see [`recompute_realign_abundance`]), so an admitted copy's
/// abundance isn't computed over a stale, pre-admission `copy_psv_alleles`/`copy_junctions` frame.
///
/// Returns `apply.novel_pools` (also GLOBAL-index space) for the caller to feed to
/// `admit_novel_pools` — this function does not attempt admission itself.
pub(crate) fn apply_realign_patch(
    fa: &mut FamilyAssignment,
    apply: super::vg_realign::RealignApply,
    p: &AssignParams,
) -> Vec<Vec<usize>> {
    fa.realign_records = apply.records;

    // Snapshot the copy roster as CopyProfiles (same shape `assign_read_editing` consumes) BEFORE
    // any correction is applied -- corrections only ever retarget an EXISTING copy (`apply_realign`
    // chooses `best_copy` from the family's already-known `copy_seqs`), so this frame is stable
    // across every corrected read below.
    let copy_profiles: Vec<CopyProfile> = fa
        .copy_psv_alleles
        .iter()
        .zip(fa.copy_junctions.iter())
        .enumerate()
        .map(|(i, (alleles, junctions))| CopyProfile {
            copy_id: i,
            alleles: alleles.clone(),
            junctions: junctions.clone(),
        })
        .collect();
    // Mirrors `assign_family_detailed`'s own one-shot computation of the RNA-editing flag (from the
    // family's PRE-correction `read_psv_obs`) so a corrected read's certificate is held to the same
    // standard as every other read's.
    let editing_cols: Vec<bool> = if p.rna_editing_filter {
        detect_editing_columns(&fa.read_psv_obs, &copy_profiles)
    } else {
        Vec::new()
    };

    for (read_idx, (new_copy, obs)) in apply.corrected {
        if let Some(pos) = fa.assignments.iter().position(|(ri, _)| *ri == read_idx) {
            let rf = ReadFeatures {
                psv_obs: obs.clone(),
                psv_qual: Vec::new(),
                junctions: fa.read_junctions[pos].clone(),
            };
            fa.assignments[pos].1 = match assign_read_editing(&rf, &copy_profiles, p, &editing_cols) {
                Some(a) => a,
                None => {
                    // `copy_profiles` is empty (no copies at all) -- can't happen in practice (a
                    // correction implies >= 1 candidate copy existed), but fall back to a
                    // consistent one-hot Assignment rather than leaving anything stale.
                    let mut posterior = vec![0.0; copy_profiles.len()];
                    if let Some(slot) = posterior.get_mut(new_copy) {
                        *slot = 1.0;
                    }
                    Assignment {
                        best_copy: new_copy,
                        log_lr_margin: 0.0,
                        n_decisive: 0,
                        resolvable: false,
                        status: AssignStatus::Assigned,
                        p_value: 0.0,
                        min_p_value: 0.0,
                        discovery_coupled: false,
                        posterior,
                    }
                }
            };
            fa.read_psv_obs[pos] = obs;
        }
    }
    apply.novel_pools
}

/// VG re-align END-TO-END plan, Task 3 (must-have, I2 fix): recompute the per-copy abundance over
/// the FINAL copy roster and observation set — call this AFTER both `apply_realign_patch` (which
/// applies corrections) and `admit_novel_pools` (which may widen `fa.copy_psv_alleles`/
/// `copy_junctions`/`copy_tids`), never before. Reruns the junction-aware EM
/// (`em_copy_assign::em_assign_family`, the SAME `AssignParams` the one-shot gate used) over
/// `fa.read_psv_obs`/`fa.copy_psv_alleles`/`fa.read_junctions`/`fa.copy_junctions` and overwrites
/// `fa.copy_abundance`, then widens `fa.copy_abundance_ci` (0.0-fill) to match, so the caller's
/// invariant `copy_abundance.len() == copy_psv_alleles.len() == copy_tids.len()` holds even when an
/// admission grew the roster. Running this BEFORE admission (the pre-fix ordering) left an admitted
/// copy's abundance implicitly `0.0` and `copy_abundance` shorter than `copy_psv_alleles`.
pub(crate) fn recompute_realign_abundance(
    fa: &mut FamilyAssignment,
    p: &AssignParams,
    em_eps: f64,
    em_max_iter: usize,
) {
    let em = super::em_copy_assign::em_assign_family(
        &fa.read_psv_obs,
        &fa.copy_psv_alleles,
        &fa.read_junctions,
        &fa.copy_junctions,
        p,
        em_eps,
        em_max_iter,
    );
    fa.copy_abundance = em.abundances;
    let n = fa.copy_psv_alleles.len();
    if fa.copy_abundance_ci.len() < n {
        fa.copy_abundance_ci.resize(n, 0.0);
    }
}

/// VG re-align END-TO-END plan, Task 3 (best-effort/mechanism-only): turn each `apply_realign`
/// novel-read pool (reads fitting NO existing copy-path at all) into a `CollapsedCandidate` and run
/// it through `absent_copy::admit_candidate` — the SAME admission gate the reference-absent-copy
/// discovery path (`DenovoConfig`'s `absent_copies` flag) already uses.
///
/// `pools` holds GLOBAL `bam_reads` indices (the caller has already translated `apply_realign`'s
/// per-family-local indices). `host` = the family's most-supported copy (max `n_reads`), matching
/// `recover_collapsed_candidates`'s "most-supported copy is the host" convention.
///
/// A pool that can't be cleanly turned into an admitted copy is left alone (no panic, no fabricated
/// admission) — it stays only as the `"novel-candidate"` rows `fa.realign_records` already carries:
///   - too few PSV-distinguishing columns among the pool's own reads to attempt a candidate at all;
///   - no identifiable haplotype among the pool's reads (`split_locus_copies` finds none);
///   - the admission gate itself declines (`Admission::DnaNeeds`, e.g. remaps back to the host
///     locus at >= 98% identity, or isn't min_p-distinct from the host).
///
/// On `Admission::Copy(t)`: `t` shares the HOST's genomic coordinate frame (it's an allele-overlay of
/// the host's spliced sequence — `absent_copy::admit_candidate`'s contract), so `t`'s alleles at the
/// family's existing PSV columns are derived the SAME way a real copy's are: the host's own
/// `family_col_genomic_pos` (genomic positions), inverted through `t`'s own `gen2off` (since `t`'s
/// intron chain can differ from the host's own). Appends the new copy to
/// `copy_psv_alleles`/`copy_junctions`/`copy_tids`, bumps `n_copies`, and assigns the pool's reads to
/// it (updating an existing `fa.assignments` entry in place, or — for a pool read the one-shot gate
/// never assigned at all — appending a new `discovery_coupled` one, keeping `read_psv_obs`/
/// `read_junctions` parallel to `assignments` and `n_reads == assignments.len()`).
///
/// M4 fix: every `Assignment::posterior` vector is widened (0.0-pad) to the FINAL copy-roster width
/// (`fa.copy_tids.len()`) before returning, then renormalized. Without this, a posterior sized
/// mid-loop — either a pre-existing read's posterior from BEFORE any admission in this family, or a
/// newly-admitted-copy read's posterior sized right after ITS pool's admission but before a LATER
/// pool admits yet another copy — stays shorter than the final roster, and the
/// `.posterior.tsv` writer's frame-width guard (`a.posterior.len() != fa.copy_tids.len()`,
/// `copy_assign.rs`) silently skips that read rather than reporting it.
fn widen_posteriors_to_final_roster(fa: &mut FamilyAssignment) {
    let final_n = fa.copy_tids.len();
    for (_, a) in fa.assignments.iter_mut() {
        if a.posterior.len() < final_n {
            a.posterior.resize(final_n, 0.0);
            let z: f64 = a.posterior.iter().sum();
            if z > 0.0 {
                for x in a.posterior.iter_mut() {
                    *x /= z;
                }
            }
        }
    }
}

/// VG re-align END-TO-END plan, Task 3 (best-effort/mechanism-only): turn each `apply_realign`
/// novel-read pool (reads fitting NO existing copy-path at all) into a `CollapsedCandidate` and run
/// it through `absent_copy::admit_candidate` — the SAME admission gate the reference-absent-copy
/// discovery path (`DenovoConfig`'s `absent_copies` flag) already uses.
///
/// `pools` holds GLOBAL `bam_reads` indices (the caller has already translated `apply_realign`'s
/// per-family-local indices). `host` = the family's most-supported copy (max `n_reads`), matching
/// `recover_collapsed_candidates`'s "most-supported copy is the host" convention.
///
/// A pool that can't be cleanly turned into an admitted copy is left alone (no panic, no fabricated
/// admission) — it stays only as the `"novel-candidate"` rows `fa.realign_records` already carries:
///   - too few PSV-distinguishing columns among the pool's own reads to attempt a candidate at all;
///   - no identifiable haplotype among the pool's reads (`split_locus_copies` finds none);
///   - the admission gate itself declines (`Admission::DnaNeeds`, e.g. remaps back to the host
///     locus at >= 98% identity, or isn't min_p-distinct from the host).
///
/// On `Admission::Copy(t)`: `t` shares the HOST's genomic coordinate frame (it's an allele-overlay of
/// the host's spliced sequence — `absent_copy::admit_candidate`'s contract), so `t`'s alleles at the
/// family's existing PSV columns are derived the SAME way a real copy's are: the host's own
/// `family_col_genomic_pos` (genomic positions), inverted through `t`'s own `gen2off` (since `t`'s
/// intron chain can differ from the host's own). Appends the new copy to
/// `copy_psv_alleles`/`copy_junctions`/`copy_tids`, bumps `n_copies`, and assigns the pool's reads to
/// it (updating an existing `fa.assignments` entry in place, or — for a pool read the one-shot gate
/// never assigned at all — appending a new `discovery_coupled` one, keeping `read_psv_obs`/
/// `read_junctions` parallel to `assignments` and `n_reads == assignments.len()`).
///
/// (M4 fix, no longer a caveat) every `Assignment::posterior` — pre-existing or freshly appended — is
/// widened to the FINAL copy-roster width before this returns; see [`widen_posteriors_to_final_roster`].
///
/// Thin wrapper over [`admit_novel_pools_with_admitter`], supplying the REAL admission gate
/// (`absent_copy::admit_candidate`, which shells minimap2 for the remap-identity check). Tests drive
/// the hermetic core directly with an injected admitter, mirroring `absent_copy.rs`'s own
/// `admit_candidate`/`admit_candidate_with_remap` split.
fn admit_novel_pools(
    fa: &mut FamilyAssignment,
    pools: &[Vec<usize>],
    bam_reads: &[BamRead],
    all_copies: &[DenovoTranscript],
    genome: &GenomeIndex,
    fasta_path: &str,
    profiles: &FamilyProfiles,
) {
    admit_novel_pools_with_admitter(fa, pools, bam_reads, all_copies, profiles, |cand, host| {
        absent_copy::admit_candidate(cand, host, genome, fasta_path, &AbsentCopyParams::default())
    })
}

/// Hermetic core of [`admit_novel_pools`]: identical mechanism, but the admission gate itself is
/// injected as `admit` rather than hardcoded to the real (minimap2-shelling) `absent_copy::admit_candidate`.
fn admit_novel_pools_with_admitter<F>(
    fa: &mut FamilyAssignment,
    pools: &[Vec<usize>],
    bam_reads: &[BamRead],
    all_copies: &[DenovoTranscript],
    profiles: &FamilyProfiles,
    admit: F,
) where
    F: Fn(&CollapsedCandidate, &DenovoTranscript) -> Admission,
{
    if pools.is_empty() || all_copies.is_empty() {
        return;
    }
    let host_idx = match (0..all_copies.len()).max_by_key(|&k| all_copies[k].n_reads) {
        Some(k) => k,
        None => return,
    };
    let host = &all_copies[host_idx];
    let host_col_gpos = family_col_genomic_pos(profiles, host_idx);

    for pool in pools {
        let pool_reads: Vec<AlignedRead> = pool.iter().map(|&gi| bam_reads[gi].read.clone()).collect();
        let psv_pos = discover_locus_psvs(&pool_reads, 3);
        if psv_pos.len() < 2 {
            continue; // too few distinguishing columns -- left as a "novel-candidate" record.
        }
        let min_reads_per_copy = pool_reads.len().min(2).max(1);
        let isos = split_locus_copies(&pool_reads, 3, 2, min_reads_per_copy);
        let n_clusters = isos.len();
        let iso = match isos.into_iter().max_by_key(|c| c.read_count) {
            Some(iso) => iso,
            None => continue, // no identifiable haplotype among the pool's reads.
        };
        let cand = CollapsedCandidate {
            host_tid: host.tid.clone(),
            chrom: host.chrom.clone(),
            start: host.start,
            end: host.end,
            iso,
            psv_pos,
            n_clusters,
        };
        let admitted = admit(&cand, host);
        let t = match admitted {
            Admission::Copy(t) => t,
            Admission::DnaNeeds(_) => continue, // gate declined -- stays a "novel-candidate" record.
        };

        let t_g2o = gen2off(&t);
        let positions = psv_positions_for(&host_col_gpos, &t_g2o, t.seq.len());
        let alleles: Vec<Option<u8>> = positions.iter().map(|&off| t.seq.get(off).copied()).collect();
        let new_copy = fa.copy_psv_alleles.len();
        fa.copy_psv_alleles.push(alleles);
        fa.copy_junctions.push(copy_boundaries(&t));
        fa.copy_tids.push(t.tid.clone());
        fa.n_copies += 1;

        for &gi in pool {
            // C1: orient the read to `t`'s own (transcription-strand) frame before the traceback --
            // same fix as `vg_realign::apply_realign`'s correction path, needed here too since `t`
            // can itself be a `-`-strand copy.
            let oriented = super::vg_realign::orient_for_copy(&bam_reads[gi].read.seq, &t.seq);
            let map = super::vg_realign::align_traceback(&oriented, &t.seq);
            let obs = super::vg_realign::path_obs_at(&map, &positions, &oriented);
            if let Some(pos) = fa.assignments.iter().position(|(ri, _)| *ri == gi) {
                fa.assignments[pos].1.best_copy = new_copy;
                fa.read_psv_obs[pos] = obs;
            } else {
                let mut posterior = vec![0.0; fa.copy_psv_alleles.len()];
                posterior[new_copy] = 1.0;
                fa.assignments.push((
                    gi,
                    Assignment {
                        best_copy: new_copy,
                        log_lr_margin: 0.0,
                        n_decisive: 0,
                        resolvable: false,
                        status: AssignStatus::Assigned,
                        p_value: 0.0,
                        min_p_value: 0.0,
                        discovery_coupled: true,
                        posterior,
                    },
                ));
                fa.read_psv_obs.push(obs);
                fa.read_junctions.push(Vec::new());
            }
        }
    }
    widen_posteriors_to_final_roster(fa);
    fa.n_reads = fa.assignments.len();
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
    // Same `k` as the O1 catalogs: one canonical extent per locus across objectives (see `detect_families`).
    let skeletons = pass1_skeletons_robust(primary_reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    let mut transcripts = assemble_gate(&skeletons, genome, &cfg.gate);

    // Unspliced readthrough spans are not transcripts, so they are removed BEFORE loci are collapsed. Filter
    // the REPS instead and a readthrough becomes the representative of every locus it happens to span: at DAZ
    // its 298 kb single-exon span absorbed DAZ2's transcripts into DAZ1's locus group, and dropping that rep
    // then deleted DAZ2 along with it. Left in the copy set at all, one becomes a "copy" — at GSTM a 30 kb
    // span covering GSTM5 + GSTM1 sat beside GSTM3 — and it makes read alignment quadratic (a 128 kb
    // "transcript" hung assignment past 400 s). Never silent: every drop is logged.
    if cfg.filter_readthrough {
        let support = read_junction_support(primary_reads);
        retain_non_readthrough(&mut transcripts, &support, "detect_and_assign");
        retain_non_mischain(&mut transcripts, &support, "detect_and_assign");
    }

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
    // Family MEMBERSHIP oracle: default = conflict-graph (E_c, AUTHORITATIVE de-tie criterion); opt-in
    // (`cfg.homology_primary`) = homology (E_r) — a copy whose reads all map uniquely (high MAPQ) forms
    // no conflict edge and would otherwise be DROPPED from its family (the read-conflict oracle is an
    // AMBIGUITY oracle, not a homology oracle). Everything downstream of `families`/`edges_f64`
    // (colocated_families, PSV discovery, assignment, EM, chi_H) is unchanged by which oracle ran.
    let placements = build_read_placements(bam_reads, &reps);
    lap!(format!("build_read_placements ({} reads x {} reps)", bam_reads.len(), reps.len()));
    let (families, edges_f64): (Vec<Vec<usize>>, Vec<(usize, usize, f64)>) = if !cfg.homology_primary {
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
        let edges_f64: Vec<(usize, usize, f64)> = c_edges.iter().map(|&(a, b, w)| (a, b, w as f64)).collect();
        (c_fams, edges_f64)
    } else {
        // Loud abort on failure — NEVER silently fall back to E_c (precedent: the famCN
        // `unwrap_or_default()` silent-degradation bug). `threads`: no thread count is threaded through
        // `detect_and_assign`/`DenovoConfig` today (this is diagnostic-scale minimap2, not the BAM-reading
        // hot path), so this uses the same default (4) `RefineParams::default()`/the existing
        // `homology_refine_params` unit tests use, rather than inventing a new config field.
        let refine = homology_refine_params(None, 4);
        let e2 = homology_edges_all_reps(&reps, &refine)
            .expect("--homology-primary: homology_edges_all_reps failed — is minimap2 on PATH or RUSTLE_MINIMAP2 set?");
        let edges_f64: Vec<(usize, usize, f64)> = e2.iter().map(|&(a, b)| (a, b, 1.0)).collect();
        let families = crate::vg_family::family_split::gamma_quasi_clique_partition(reps.len(), &edges_f64, 0.20);
        eprintln!("[detect_and_assign] homology (E_r): {} edges -> {} families", e2.len(), families.len());
        (families, edges_f64)
    };
    let split = to_split_families(&families, &edges_f64, &cfg.split);

    // POA homology edges — DIAGNOSTIC ONLY, no longer drives family membership (families come from the
    // membership oracle above: the de-tie conflict graph by default, or E_r homology when
    // `cfg.homology_primary` is set). The POA pairwise (poasta over all
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
        // Task 5 / VG re-align END-TO-END plan Task 3 (opt-in, additive): OFF (default) => this whole
        // block is skipped, `vg_realign_apply` stays `None`, `fa.realign_records` stays empty, and
        // nothing else in `fa`/`detail`/`all_copies` is touched -- `fa` is byte-identical to the
        // pre-Task-3 build. ON => `vg_realign::apply_realign` runs over this family's reads (the SAME
        // `region`-filtered `BamRead`s the old report-only wiring used); its corrections + admissions
        // are applied to `fa` AFTER the per-read loop below populates `fa.assignments`/`read_psv_obs`
        // (see `apply_realign_patch`/`admit_novel_pools`).
        let mut vg_realign_apply: Option<(super::vg_realign::RealignApply, FamilyProfiles)> = None;
        if cfg.vg_realign {
            let family_bam_reads: Vec<BamRead> = idx_map.iter().map(|&i| bam_reads[i].clone()).collect();
            let copy_refs: Vec<&DenovoTranscript> = all_copies.iter().collect();
            // The family's PSV columns in EACH copy's own spliced-consensus offset frame (what
            // `apply_realign`'s `path_obs_at` re-extraction needs): `build_family_profiles` keys PSV
            // columns by FORWARD-GENOME position per copy (`FamilyProfiles::copy_gpos`), not by
            // spliced offset directly, so invert via each copy's own `gen2off`.
            let profiles = build_family_profiles(&copy_refs, Some(genome));
            let copy_seqs: Vec<Vec<u8>> = all_copies.iter().map(|c| c.seq.clone()).collect();
            let psv_pos_per_copy: Vec<Vec<usize>> = (0..all_copies.len())
                .map(|k| {
                    let col_gpos = family_col_genomic_pos(&profiles, k);
                    let g2o = gen2off(&all_copies[k]);
                    psv_positions_for(&col_gpos, &g2o, copy_seqs[k].len())
                })
                .collect();
            let linear_copy_of: Vec<Option<usize>> = family_bam_reads
                .iter()
                .map(|br| best_overlap_copy(&br.read, &copy_refs))
                .collect();
            let apply = super::vg_realign::apply_realign(
                &family_bam_reads,
                &all_copies,
                &copy_seqs,
                &psv_pos_per_copy,
                &linear_copy_of,
                &super::vg_realign::RealignParams::default(),
                p.error_rate,
                p.alpha,
            );
            // `apply`'s indices are LOCAL to `family_bam_reads`; translate to the GLOBAL `bam_reads`
            // indices `fa.assignments`/`idx_map[r.read_index]` are keyed on.
            let apply_global = super::vg_realign::RealignApply {
                corrected: apply.corrected.into_iter().map(|(l, v)| (idx_map[l], v)).collect(),
                novel_pools: apply
                    .novel_pools
                    .iter()
                    .map(|pl| pl.iter().map(|&l| idx_map[l]).collect())
                    .collect(),
                records: apply.records,
            };
            vg_realign_apply = Some((apply_global, profiles));
        }
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
            copy_junction_support: vec![0; all_copies.len()],
            uniq: 0,
            uniq_agree: 0,
            collapsed_copies,
            rescued_copies,
            assignments: Vec::with_capacity(detail.results.len()),
            copy_tids: all_copies.iter().map(|c| c.tid.clone()).collect(),
            copy_spans: all_copies.iter().map(|c| (c.chrom.clone(), c.start, c.end)).collect(),
            copy_abundance: detail.copy_abundance.clone(),
            copy_abundance_ci: detail.copy_abundance_ci.clone(),
            mosaic_reads: detail.mosaic_reads,
            conversions: detail.conversions.clone(),
            copy_conversions: detail.copy_conversions.clone(),
            psv_col_pos: detail.psv_col_pos.clone(),
            copy_psv_alleles: detail.copy_psv_alleles.clone(),
            read_psv_obs: Vec::with_capacity(detail.results.len()),
            copy_junctions: detail.copy_junctions.clone(),
            read_junctions: Vec::with_capacity(detail.results.len()),
            realign_records: Vec::new(),
        };
        for r in detail.results {
            let resolvable_psv = r.psv.n_decisive >= 1;
            let assigned_j = r.combined.status == AssignStatus::Assigned;
            fa.resolvable_psv += resolvable_psv as usize;
            fa.assigned_psv += (r.psv.status == AssignStatus::Assigned) as usize;
            fa.resolvable_j += (r.combined.n_decisive >= 1) as usize;
            fa.assigned_j += assigned_j as usize;
            fa.junction_only += (r.combined.n_decisive >= 1 && !resolvable_psv) as usize;
            add_junction_support(
                &mut fa.copy_junction_support,
                r.combined.best_copy,
                r.combined.n_decisive >= 1,
                resolvable_psv,
            );
            if assigned_j && region_mapq[r.read_index] > 0 {
                fa.uniq += 1;
                fa.uniq_agree += (r.combined.best_copy == r.mapped_copy) as usize;
            }
            fa.read_psv_obs.push(r.psv_obs);
            fa.read_junctions.push(r.junctions);
            fa.assignments.push((idx_map[r.read_index], r.combined));
        }
        // VG re-align END-TO-END plan, Task 3: apply the corrections (must-have), attempt admission
        // of any novel-read pools (best-effort) now that `fa.assignments`/`read_psv_obs`/
        // `read_junctions` are populated, THEN recompute the EM abundance (I2 fix: must run AFTER
        // admission, over the FINAL widened copy roster, or an admitted copy reports abundance 0.0
        // and a shorter `copy_abundance` than `copy_psv_alleles`). No-op (and `fa` untouched from the
        // pre-Task-3 shape) when `cfg.vg_realign` is off, since `vg_realign_apply` is `None` in that case.
        if let Some((apply_global, profiles)) = vg_realign_apply {
            let novel_pools = apply_realign_patch(&mut fa, apply_global, p);
            admit_novel_pools(&mut fa, &novel_pools, bam_reads, &all_copies, genome, fasta_path, &profiles);
            recompute_realign_abundance(&mut fa, p, VG_REALIGN_EM_EPS, VG_REALIGN_EM_MAX_ITER);
        }
        out.push(fa);
    }

    // COLLAPSE GATE. `colocated_families` needs >= min_copies assembled REPS, so a collapsed locus — one whose
    // copies are near-identical, leaving the aligner to pile every read onto one of them — is dropped. On chrY
    // that is DAZ and BPY2: both copies present in the reads, zero families emitted. Reps that no family
    // claimed are tested here. Leg 1 asks whether the locus is collapsed at all (ambiguously-placed reads at a
    // rate incompatible with a unique locus); ONLY then does leg 2 count PSV haplotypes. Running leg 2 alone
    // made the single-copy gene TSPYL1 report 12 copies against DAZ's 3 (`bench/COLLAPSED_COPY_GATE.md`).
    if cfg.collapse_gate {
        // owned, so `out` is free to grow below
        let familied: std::collections::HashSet<String> =
            out.iter().flat_map(|fa| fa.copy_tids.iter().cloned()).collect();
        for rep in reps.iter().filter(|r| !familied.contains(&r.tid)) {
            let obs = locus_ambiguity(rep, bam_reads);
            let reads: Vec<AlignedRead> = bam_reads
                .iter()
                .filter(|br| {
                    !br.is_supplementary
                        && br.chrom == rep.chrom
                        && br.read.ref_start < rep.end
                        && read_ref_end(&br.read) > rep.start
                })
                .map(|br| br.read.clone())
                .collect();
            let haplotypes: Vec<Vec<Option<u8>>> = split_locus_copies(&reads, 3, 2, 3)
                .into_iter()
                .filter(|c| c.identifiable)
                .map(|c| c.allele_vector)
                .collect();
            match collapse_verdict(obs, cfg.eps_amb, &haplotypes, p.alpha, min_copies) {
                CollapseVerdict::Fire { chi_h, p_value } => {
                    eprintln!(
                        "[detect_and_assign] collapse gate: {}:{}-{} ambiguous {}/{} p={p_value:.2e} -> \
                         chi(H)={chi_h} copies (reads certified TIED; chi(H) is a LOWER BOUND)",
                        rep.chrom, rep.start, rep.end, obs.k, obs.n
                    );
                    let fid = format!("DSFAM{}", out.len());
                    out.push(gated_family(rep, bam_reads, chi_h, fid));
                }
                CollapseVerdict::Abstain(why) => eprintln!(
                    "[detect_and_assign] collapse gate ABSTAINS at {}:{}-{}: {why}",
                    rep.chrom, rep.start, rep.end
                ),
                CollapseVerdict::NotCollapsed { .. } => {}
            }
        }
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
/// DRY core for the genome-wide conflict catalog AND the single-copy baseline: builds genome-wide reps
/// (skeletons -> gate -> readthrough filter -> span-aware collapse), per-chrom conflict edges, and the
/// co-located families, and returns BOTH the reps and the family catalog. `detect_conflict_catalog_genome_wide`
/// keeps only the families; `detect_single_copy_baseline_genome_wide` takes the reps a family does not claim.
fn gw_reps_and_catalog(
    bam_path: &str,
    fasta_path: &str,
    threads: usize,
    win: u64,
    min_copies: usize,
    cfg: &DenovoConfig,
) -> Result<(Vec<DenovoTranscript>, Vec<u32>, Vec<ColocatedFamily>)> {
    // --- (1) genome-wide reps ---
    let reads = primary_reads_from_bam(bam_path, threads)?;
    let contigs: HashSet<String> = reads.iter().map(|r| r.chrom.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;
    let skeletons = pass1_skeletons_robust(&reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    // Build the junction-support map BEFORE freeing the primaries. It is keyed by distinct junction, so it is a
    // rounding error next to the ~1.7M reads it derives from. O1 must drop readthroughs for the same reason O2
    // does, or the two objectives disagree about what a locus is.
    let rt_support = if cfg.filter_readthrough { Some(read_junction_support(&reads)) } else { None };
    drop(reads); // free the ~1.7M primaries before the per-chrom read load
    let mut transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
    if let Some(sup) = &rt_support {
        retain_non_readthrough(&mut transcripts, sup, "gw-catalog");
        retain_non_mischain(&mut transcripts, sup, "gw-catalog");
    }
    // rep_totals[k] is the LOCUS TOTAL reads for reps[k] (all isoforms summed) — the single-copy expression
    // basis for lambda_global, on the same footing as the family-total E_fam in depth_cn.
    let (rep_idx, rep_totals) = collapse_loci_span_aware_with_totals(&transcripts, &cfg.detect);
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
    Ok((reps, rep_totals, catalog))
}

/// Genome-wide multi-copy family catalog (the conflict oracle). Thin wrapper over `gw_reps_and_catalog`.
pub fn detect_conflict_catalog_genome_wide(
    bam_path: &str,
    fasta_path: &str,
    threads: usize,
    win: u64,
    min_copies: usize,
    cfg: &DenovoConfig,
) -> Result<Vec<ColocatedFamily>> {
    let (_reps, _rep_totals, catalog) = gw_reps_and_catalog(bam_path, fasta_path, threads, win, min_copies, cfg)?;
    Ok(catalog)
}

/// Genome-wide single-copy baseline: the reps that no family claims, as lightweight records. Same traversal as
/// the conflict catalog (DRY via `gw_reps_and_catalog`); a locus is single-copy iff it is not a copy of any
/// emitted family. Feeds lambda_global and the `.single_copy.tsv` baseline table.
pub fn detect_single_copy_baseline_genome_wide(
    bam_path: &str,
    fasta_path: &str,
    threads: usize,
    win: u64,
    min_copies: usize,
    cfg: &DenovoConfig,
) -> Result<Vec<crate::vg_family::single_copy::SingleCopyLocus>> {
    let (reps, rep_totals, catalog) = gw_reps_and_catalog(bam_path, fasta_path, threads, win, min_copies, cfg)?;
    Ok(crate::vg_family::single_copy::single_copy_loci(&reps, &rep_totals, &catalog))
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
    // Support map before the free; see the same-chrom catalog for why.
    let rt_support = if cfg.filter_readthrough { Some(read_junction_support(&reads)) } else { None };
    drop(reads);
    let mut transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
    if let Some(sup) = &rt_support {
        retain_non_readthrough(&mut transcripts, sup, "gw-catalog");
        retain_non_mischain(&mut transcripts, sup, "gw-catalog");
    }
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
    // Support map before the free; see the same-chrom catalog for why.
    let rt_support = if cfg.filter_readthrough { Some(read_junction_support(&reads)) } else { None };
    drop(reads);
    let mut transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
    if let Some(sup) = &rt_support {
        retain_non_readthrough(&mut transcripts, sup, "gw-catalog");
        retain_non_mischain(&mut transcripts, sup, "gw-catalog");
    }
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
    /// USED IN TWO PLACES: (1) `refine_families_exon_sum` — within an ALREADY-formed (conflict/homology)
    /// family, as an additional membership edge; (2) `homology_edges_all_reps` (the E_r homology-primary
    /// catalog) — as a THIRD genome-wide definition-edge tier alongside asm20/sensitive, promoting protein
    /// from orthogonal QC to a DEFINITION edge that groups coding paralogs past the ~0.65 nt-identity floor
    /// where both nucleotide tiers find no edge and reads map all-primary.
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

/// E_r homology edges over ALL reps' exon-sum sequences: asm20 (recent) ∪ sensitive -k11 -w5 (ancient) ∪
/// (opt-in) protein, all gated by `min_coverage`. One minimap2 all-vs-all per nucleotide tier over the
/// whole rep set (minimap2's index is the prefilter). When `params.protein_tail` is set, protein homology
/// is promoted from orthogonal QC to a THIRD definition edge: coding paralogs that have diverged past the
/// nucleotide seeds' floor (~0.65 identity, where both nt tiers find no edge and reads map all-primary)
/// still share a conserved protein (synonymous divergence) — this is the only way to recover that tail.
/// `batch_protein_edges` is family-scoped, so all reps are passed as ONE family to get within-set edges
/// over the whole rep universe (local indices == global indices since there is exactly one "family").
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
    if params.protein_tail {
        let prot = batch_protein_edges(std::slice::from_ref(&reps.to_vec()), 0.50, params.min_coverage, params)?;
        if let Some(edges) = prot.first() {
            for &e in edges {
                set.insert(e);
            }
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

/// Orthogonal protein QC (NEVER a definition edge): does each nt-defined family also cohere at the protein
/// level? `Some(true)` = its ORFs form a connected protein-homology graph, `Some(false)` = they do not,
/// `None` = mmseqs unavailable / protein tier off (no effect on membership).
pub fn family_protein_coheres(families: &[Vec<DenovoTranscript>], params: &RefineParams) -> Vec<Option<bool>> {
    if !params.protein_tail {
        return vec![None; families.len()];
    }
    let edges = match batch_protein_edges(families, 0.50, params.min_coverage, params) {
        Ok(e) => e,
        Err(_) => return vec![None; families.len()],
    };
    families.iter().enumerate().map(|(fi, fam)| {
        let fe = edges.get(fi).cloned().unwrap_or_default();
        if fam.len() < 2 { return Some(true); }
        Some(homology_components(fam.len(), &fe).iter().any(|c| c.len() == fam.len()))
    }).collect()
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
    fn add_junction_support_counts_only_pure_junction_reads_per_copy() {
        let mut sup = vec![0usize, 0, 0];
        // a copy-specific junction resolved this read (combined) where PSVs could not -> counts for copy 1.
        add_junction_support(&mut sup, 1, true, false);
        assert_eq!(sup, vec![0, 1, 0]);
        // PSVs already resolved it -> not junction-only -> no change.
        add_junction_support(&mut sup, 1, true, true);
        assert_eq!(sup, vec![0, 1, 0], "PSV-decisive read is not junction-only support");
        // combined not decisive -> no junction pinned it -> no change.
        add_junction_support(&mut sup, 2, false, false);
        assert_eq!(sup, vec![0, 1, 0]);
        // out-of-range copy index is ignored (no panic).
        add_junction_support(&mut sup, 99, true, false);
        assert_eq!(sup, vec![0, 1, 0], "out-of-range best_copy is guarded");
    }

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
                "  {} {} copies={} reads={} PSVcols={} resolv_PSV={:.0}% resolv_+J={:.0}% J_only={} assign_+J={:.0}% uniq_agree={}/{}",
                fa.family_id, fa.chrom, fa.n_copies, fa.n_reads, fa.psv_cols,
                pct(fa.resolvable_psv), pct(fa.resolvable_j), fa.junction_only,
                pct(fa.assigned_j), fa.uniq_agree, fa.uniq,
            );
        }
        let (tot_uniq, tot_agree): (usize, usize) = fas.iter().fold((0, 0), |(u, a), f| (u + f.uniq, a + f.uniq_agree));
        if tot_uniq > 0 {
            eprintln!("AGGREGATE unique-mapper agreement: {tot_agree}/{tot_uniq} ({:.1}%)", 100.0 * tot_agree as f64 / tot_uniq as f64);
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

    // -----------------------------------------------------------------------------------------
    // VG re-align END-TO-END plan, Task 3: patch FamilyAssignment under --vg-realign.
    // -----------------------------------------------------------------------------------------

    /// `cfg.vg_realign = false` (the default) must leave `FamilyAssignment` exactly as the
    /// pre-Task-3 build produced it: `realign_records` empty, and the SAME best_copy assignments
    /// `detect_and_assign_resolves_multimapper_end_to_end` already asserts (that test predates
    /// Task 3 and is unchanged by it -- this test additionally pins `realign_records` explicitly).
    #[test]
    fn vg_realign_off_is_byte_identical() {
        let (genome, primary, aligned) = two_paralogs_with_psvs();
        let cfg = DenovoConfig::default();
        assert!(!cfg.vg_realign, "default must be OFF");
        let (fas, _fallback, dna_needs) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &cfg,
            5_000_000,
            2,
            &AssignParams::default(),
            &[],
            false,
            "",
        );
        assert!(dna_needs.is_empty());
        assert_eq!(fas.len(), 1, "one co-located 2-copy family");
        let fa = &fas[0];
        assert!(fa.realign_records.is_empty(), "vg_realign OFF -> the vg_realign block must not run at all");
        assert_eq!(fa.n_copies, 2);
        assert_eq!(fa.assignments.len(), 6);
        // Same assignments the pre-Task-3 test (`detect_and_assign_resolves_multimapper_end_to_end`)
        // pins: copyB read (locus 0) -> copy 1, copyA read (locus 1000) -> copy 0, both Assigned.
        let primary_assign = fa.assignments.iter().find(|(_, a)| a.best_copy == 1)
            .expect("primary read (copyB seq at locus 0) should resolve to copyB (copy 1)");
        assert_eq!(primary_assign.1.status, AssignStatus::Assigned);
        let secondary_assign = fa.assignments.iter().find(|(_, a)| a.best_copy == 0)
            .expect("secondary read (copyA seq at locus 1000) should resolve to copyA (copy 0)");
        assert_eq!(secondary_assign.1.status, AssignStatus::Assigned);
    }

    /// `apply_realign_patch` directly (per the task's suggested factoring): a planted correction
    /// (`apply.corrected`) must move the read's `Assignment::best_copy`, overwrite its parallel
    /// `read_psv_obs` entry, and the recomputed `copy_abundance` must be the EM's actual output over
    /// the corrected evidence (not the stale pre-correction value) -- proving both the correction AND
    /// the EM recompute (steps 2 + 4) are wired, independent of the heavier full-pipeline plumbing
    /// (`apply_realign`'s PSV-position derivation, admissions) exercised elsewhere.
    #[test]
    fn correction_updates_assignment_and_recomputes_em() {
        use crate::vg_family::vg_realign::{RealignApply, RealignRecord};
        use std::collections::HashMap;

        // One PSV column distinguishing copy 0 ('A') from copy 1 ('C').
        let copy_psv_alleles = vec![vec![Some(b'A')], vec![Some(b'C')]];
        let copy_junctions = vec![Vec::new(), Vec::new()];

        // Read (global bam index) 10: the one-shot gate left it Tied (uncovered PSV, best_copy
        // defaulted to 0) -- the vg-realign significance certificate (full-sequence identity,
        // computed upstream of `apply_realign_patch`) supplies the correction to copy 1 instead.
        // Read 11 is untouched (no correction), already Assigned to copy 1.
        let mut fa = FamilyAssignment {
            family_id: "fam0".into(),
            chrom: "c1".into(),
            n_copies: 2,
            n_reads: 2,
            psv_cols: 1,
            resolvable_psv: 0,
            assigned_psv: 0,
            resolvable_j: 0,
            assigned_j: 0,
            junction_only: 0,
            copy_junction_support: Vec::new(),
            uniq: 0,
            uniq_agree: 0,
            collapsed_copies: 0,
            rescued_copies: 0,
            assignments: vec![
                (
                    10,
                    Assignment {
                        best_copy: 0,
                        log_lr_margin: 0.0,
                        n_decisive: 0,
                        resolvable: false,
                        status: AssignStatus::Tied,
                        p_value: 1.0,
                        min_p_value: 1.0,
                        discovery_coupled: false,
                        posterior: vec![0.5, 0.5],
                    },
                ),
                (
                    11,
                    Assignment {
                        best_copy: 1,
                        log_lr_margin: 5.0,
                        n_decisive: 1,
                        resolvable: true,
                        status: AssignStatus::Assigned,
                        p_value: 1e-6,
                        min_p_value: 1e-6,
                        discovery_coupled: false,
                        posterior: vec![0.0, 1.0],
                    },
                ),
            ],
            copy_tids: vec!["c0".into(), "c1".into()],
            copy_spans: vec![("c".into(), 0, 10), ("c".into(), 20, 30)],
            copy_abundance: vec![0.5, 0.5],
            copy_abundance_ci: vec![0.5, 0.5],
            mosaic_reads: 0,
            conversions: Vec::new(),
            copy_conversions: Vec::new(),
            psv_col_pos: vec![Some(100)],
            copy_psv_alleles: copy_psv_alleles.clone(),
            read_psv_obs: vec![vec![None], vec![Some(b'C')]],
            copy_junctions: copy_junctions.clone(),
            read_junctions: vec![Vec::new(), Vec::new()],
            realign_records: Vec::new(),
        };

        let mut corrected = HashMap::new();
        corrected.insert(10usize, (1usize, vec![Some(b'C')]));
        let apply = RealignApply {
            corrected,
            novel_pools: Vec::new(),
            records: vec![RealignRecord {
                read_name: "r10".into(),
                action: "reassigned".into(),
                target_copy: 1,
                id_best: 0.99,
                linear_copy: 0,
            }],
        };

        let p = AssignParams::default();
        let novel_pools = apply_realign_patch(&mut fa, apply, &p);
        assert!(novel_pools.is_empty(), "no novel pools were planted");
        assert_eq!(fa.realign_records.len(), 1, "records must be set on fa");
        assert_eq!(fa.realign_records[0].action, "reassigned");

        let pos = fa.assignments.iter().position(|(ri, _)| *ri == 10).expect("read 10 still present");
        assert_eq!(fa.assignments[pos].1.best_copy, 1, "correction must move read 10 to copy 1");
        assert_eq!(fa.read_psv_obs[pos], vec![Some(b'C')], "read_psv_obs must be overwritten in lockstep");
        // read 11 (untouched) must be unaffected.
        let pos11 = fa.assignments.iter().position(|(ri, _)| *ri == 11).unwrap();
        assert_eq!(fa.assignments[pos11].1.best_copy, 1);

        // (I2 fix) `apply_realign_patch` alone no longer recomputes `copy_abundance` -- that must
        // happen AFTER `admit_novel_pools`, via the separate `recompute_realign_abundance` step
        // (see `admission_widens_abundance_not_left_at_zero` for the admission-ordering case).
        assert_eq!(
            fa.copy_abundance,
            vec![0.5, 0.5],
            "apply_realign_patch alone must not touch copy_abundance any more"
        );
        recompute_realign_abundance(&mut fa, &p, 1e-6, 500);

        // copy_abundance must be the EM's ACTUAL recomputed output over the corrected read_psv_obs
        // (both reads now carry copy 1's allele), not the stale pre-correction 0.5/0.5.
        let expected = crate::vg_family::em_copy_assign::em_assign_family(
            &fa.read_psv_obs,
            &fa.copy_psv_alleles,
            &fa.read_junctions,
            &fa.copy_junctions,
            &p,
            1e-6,
            500,
        )
        .abundances;
        assert_eq!(fa.copy_abundance, expected, "copy_abundance must equal em_assign_family's own recompute");
        assert!(
            fa.copy_abundance[1] > fa.copy_abundance[0],
            "both reads now favor copy 1 post-correction: abundance = {:?}",
            fa.copy_abundance
        );
    }

    /// I3 fix: a correction must RE-DERIVE the read's entire `Assignment` (status, posterior,
    /// p_value, margin, n_decisive) from the corrected obs, not just move `best_copy` -- else
    /// `.assignments.tsv`/`.posterior.tsv` report a self-inconsistent story (e.g. a reassigned read
    /// still printing `status=Tied` with its stale pre-correction posterior).
    #[test]
    fn correction_rederives_full_assignment_not_just_best_copy() {
        use crate::vg_family::vg_realign::{RealignApply, RealignRecord};
        use std::collections::HashMap;

        // Two PSV columns cleanly distinguishing copy 0 ('A','A') from copy 1 ('C','C') -- well
        // clear of the alpha=1e-3 Bonferroni bound (min_p ~= (0.003/3)^2 = 1e-6 << 1e-3), so the
        // fresh certificate is unambiguously Assigned (not stuck on a single-column boundary tie).
        let copy_psv_alleles = vec![vec![Some(b'A'), Some(b'A')], vec![Some(b'C'), Some(b'C')]];
        let copy_junctions = vec![Vec::new(), Vec::new()];

        // Read (global bam index) 20: the one-shot gate never spanned either column (both `None`)
        // -> genuinely Tied, flat posterior, best_copy defaulted to 0. The vg-realign correction
        // re-aligns it against copy 1's full consensus and DOES span both columns, matching copy 1
        // exactly at both -- a clean, decisive correction.
        let mut fa = FamilyAssignment {
            family_id: "fam1".into(),
            chrom: "c1".into(),
            n_copies: 2,
            n_reads: 1,
            psv_cols: 2,
            resolvable_psv: 0,
            assigned_psv: 0,
            resolvable_j: 0,
            assigned_j: 0,
            junction_only: 0,
            copy_junction_support: Vec::new(),
            uniq: 0,
            uniq_agree: 0,
            collapsed_copies: 0,
            rescued_copies: 0,
            assignments: vec![(
                20,
                Assignment {
                    best_copy: 0,
                    log_lr_margin: 0.0,
                    n_decisive: 0,
                    resolvable: false,
                    status: AssignStatus::Tied,
                    p_value: 1.0,
                    min_p_value: 1.0,
                    discovery_coupled: false,
                    posterior: vec![0.5, 0.5],
                },
            )],
            copy_tids: vec!["c0".into(), "c1".into()],
            copy_spans: vec![("c".into(), 0, 10), ("c".into(), 20, 30)],
            copy_abundance: vec![0.5, 0.5],
            copy_abundance_ci: vec![0.5, 0.5],
            mosaic_reads: 0,
            conversions: Vec::new(),
            copy_conversions: Vec::new(),
            psv_col_pos: vec![Some(100), Some(200)],
            copy_psv_alleles: copy_psv_alleles.clone(),
            read_psv_obs: vec![vec![None, None]],
            copy_junctions: copy_junctions.clone(),
            read_junctions: vec![Vec::new()],
            realign_records: Vec::new(),
        };

        let mut corrected = HashMap::new();
        corrected.insert(20usize, (1usize, vec![Some(b'C'), Some(b'C')]));
        let apply = RealignApply {
            corrected,
            novel_pools: Vec::new(),
            records: vec![RealignRecord {
                read_name: "r20".into(),
                action: "reassigned".into(),
                target_copy: 1,
                id_best: 0.99,
                linear_copy: 0,
            }],
        };

        let p = AssignParams::default();
        apply_realign_patch(&mut fa, apply, &p);

        let pos = fa.assignments.iter().position(|(ri, _)| *ri == 20).expect("read 20 still present");
        let a = &fa.assignments[pos].1;
        assert_eq!(a.best_copy, 1, "must move to copy 1");
        assert_eq!(a.status, AssignStatus::Assigned, "must NOT keep the stale Tied status");
        assert!(a.posterior[1] > 0.9, "posterior must be one-hot-ish on copy 1, got {:?}", a.posterior);
        assert!(a.posterior[0] < 0.1, "posterior must not keep the stale flat 0.5, got {:?}", a.posterior);
        assert!(a.p_value < 1e-3, "p_value must reflect the fresh decisive certificate, not the stale 1.0, got {}", a.p_value);
        assert!(a.n_decisive >= 1, "n_decisive must reflect the fresh spanned columns, not the stale 0");
    }

    /// I2 fix: `copy_abundance` must be recomputed AFTER admission widens the copy roster, not
    /// before -- otherwise an admitted copy reports abundance 0.0 and `copy_abundance` is shorter
    /// than `copy_psv_alleles`. Drives the hermetic `admit_novel_pools_with_admitter` (an injected
    /// admitter bypasses the real minimap2-shelling gate -- `absent_copy.rs`'s own tests already
    /// cover THAT gate in isolation, per its `admit_candidate`/`admit_candidate_with_remap` split)
    /// followed by `recompute_realign_abundance`, mirroring the production `detect_and_assign`
    /// ordering.
    #[test]
    fn admission_widens_abundance_not_left_at_zero() {
        // Host copy: 80bp, '+' strand, no introns, chr1:0-80. Two positions (10, 40) will carry the
        // pool's internal PSVs; the family already tracks column 0 at genomic position 10.
        let backbone = rand_seq(80, 0xD00D);
        let mut mutated = backbone.clone();
        for &pos in &[10usize, 40usize] {
            let orig = backbone[pos];
            mutated[pos] = [b'A', b'C', b'G', b'T'].into_iter().find(|&b| b != orig).unwrap();
        }

        let host = DenovoTranscript {
            tid: "host".into(), chrom: "chr1".into(), start: 0, end: 80, n_reads: 7, strand: '+',
            introns: vec![], seq: backbone.clone(),
        };
        let all_copies = vec![host.clone()];

        let profiles = FamilyProfiles {
            profiles: vec![CopyProfile { copy_id: 0, alleles: vec![Some(backbone[10])], junctions: vec![] }],
            copy_gpos: vec![vec![(0usize, 10u64)]],
            gen2off: vec![gen2off(&host)],
            strand: vec!['+'],
            n_cols: 1,
        };

        let mut fa = FamilyAssignment {
            family_id: "famX".into(),
            chrom: "chr1".into(),
            n_copies: 1,
            n_reads: 0,
            psv_cols: 1,
            resolvable_psv: 0,
            assigned_psv: 0,
            resolvable_j: 0,
            assigned_j: 0,
            junction_only: 0,
            copy_junction_support: Vec::new(),
            uniq: 0,
            uniq_agree: 0,
            collapsed_copies: 0,
            rescued_copies: 0,
            assignments: Vec::new(),
            copy_tids: vec!["host".into()],
            copy_spans: vec![("c".into(), 0, 10)],
            copy_abundance: vec![1.0],
            copy_abundance_ci: vec![0.0],
            mosaic_reads: 0,
            conversions: Vec::new(),
            copy_conversions: Vec::new(),
            psv_col_pos: vec![Some(10)],
            copy_psv_alleles: vec![vec![Some(backbone[10])]],
            read_psv_obs: Vec::new(),
            copy_junctions: vec![Vec::new()],
            read_junctions: Vec::new(),
            realign_records: Vec::new(),
        };

        fn mk_bam(name: &str, seq: Vec<u8>) -> BamRead {
            let len = seq.len() as u64;
            BamRead {
                chrom: "chr1".into(),
                read: AlignedRead { ref_start: 0, cigar: vec![('M', len)], seq, qual: vec![] },
                mapq: 3,
                name: name.into(),
                as_score: 0,
                de: 0.0,
                is_supplementary: false,
            }
        }

        // 3 reads carrying the host's own backbone (haplotype A) + 4 reads carrying the mutated
        // haplotype (B) -- both groups clear `discover_locus_psvs`'s `min_allele_reads = 3` at
        // BOTH columns (10, 40), so the pool splits into 2 identifiable clusters; the larger (B,
        // n=4) is the one `admit_novel_pools_with_admitter` carries into the (mocked) admission gate.
        let mut bam_reads = Vec::new();
        for i in 0..3 {
            bam_reads.push(mk_bam(&format!("hapA_{i}"), backbone.clone()));
        }
        for i in 0..4 {
            bam_reads.push(mk_bam(&format!("hapB_{i}"), mutated.clone()));
        }
        let pools = vec![(0..7).collect::<Vec<usize>>()];

        let admitted_t = DenovoTranscript {
            tid: "admitted1".into(), chrom: "chr1".into(), start: 0, end: 80, n_reads: 4, strand: '+',
            introns: vec![], seq: mutated.clone(),
        };

        admit_novel_pools_with_admitter(&mut fa, &pools, &bam_reads, &all_copies, &profiles, |_c, _h| {
            Admission::Copy(admitted_t.clone())
        });

        assert_eq!(fa.n_copies, 2, "the pool must be admitted as a new copy");
        assert_eq!(fa.copy_psv_alleles.len(), 2);
        assert_eq!(fa.copy_tids, vec!["host".to_string(), "admitted1".to_string()]);
        assert_eq!(fa.assignments.len(), 7, "all 7 pool reads must be assigned (none pre-existed)");

        let p = AssignParams::default();
        recompute_realign_abundance(&mut fa, &p, 1e-6, 500);

        assert_eq!(
            fa.copy_abundance.len(),
            fa.copy_psv_alleles.len(),
            "I2: copy_abundance must be widened to the FINAL (post-admission) copy roster"
        );
        assert!(
            fa.copy_abundance[1] > 0.0,
            "the admitted copy's abundance must be > 0 -- its 4 supporting reads carry its exact \
             allele, got {:?}",
            fa.copy_abundance
        );
        assert_eq!(fa.copy_abundance_ci.len(), 2, "copy_abundance_ci must be widened too");
    }

    /// Full-pipeline smoke test with `cfg.vg_realign = true`: exercises the `psv_pos_per_copy`
    /// derivation (`build_family_profiles` -> `family_col_genomic_pos` -> `gen2off` inversion) that
    /// `apply_realign_patch`'s hand-built unit test above does NOT touch (there the PSV-position
    /// mapping is supplied directly). The fixture's low-MAPQ (0) secondary reads are re-align
    /// candidates; `apply_realign` should re-derive the SAME copy call the one-shot PSV gate already
    /// makes (both resolve via the same 3-PSV signal), so this also pins that the vg-realign path
    /// does not regress the already-correct assignment on a case where the two mechanisms agree.
    #[test]
    fn vg_realign_on_end_to_end_smoke() {
        let (genome, primary, aligned) = two_paralogs_with_psvs();
        let mut cfg = DenovoConfig::default();
        cfg.vg_realign = true;
        let (fas, _fallback, dna_needs) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &cfg,
            5_000_000,
            2,
            &AssignParams::default(),
            &[],
            false,
            "",
        );
        assert!(dna_needs.is_empty());
        assert_eq!(fas.len(), 1);
        let fa = &fas[0];
        assert_eq!(fa.n_copies, 2, "no spurious admissions on this clean fixture");
        assert_eq!(fa.assignments.len(), 6);
        assert!(!fa.realign_records.is_empty(), "the 3 low-MAPQ secondary reads must be candidates");
        assert!(
            fa.realign_records.iter().all(|r| r.action == "reassigned" || r.action == "rejected"),
            "no unfit/novel reads on this fixture: {:?}",
            fa.realign_records
        );
        // The vg-realign path must agree with (not regress) the one-shot PSV gate's already-correct
        // calls on this fixture.
        let primary_assign = fa.assignments.iter().find(|(_, a)| a.best_copy == 1)
            .expect("copyB read still resolves to copy 1");
        assert_eq!(primary_assign.1.status, AssignStatus::Assigned);
        let secondary_assign = fa.assignments.iter().find(|(_, a)| a.best_copy == 0)
            .expect("copyA read still resolves to copy 0");
        assert_eq!(secondary_assign.1.status, AssignStatus::Assigned);
        // copy_abundance was recomputed via em_assign_family -- still a valid distribution.
        assert_eq!(fa.copy_abundance.len(), 2);
        let sum: f64 = fa.copy_abundance.iter().sum();
        assert!((sum - 1.0).abs() < 1e-6, "copy_abundance must sum to 1, got {:?}", fa.copy_abundance);
    }

    /// MAGEA4 (`+`) / MAGEA10 (`-`) and RFPL2 (`-`) / RFPL3 (`+`) are real INVERTED DUPLICATE paralog pairs:
    /// homologous, on one chromosome, at DISJOINT loci, on opposite strands. Partitioning `colocated_families`
    /// by `(chrom, strand)` split each pair into two singletons, which `min_copies = 2` then dropped — so O2
    /// could never assign an inverted duplicate, even when the homology oracle found the edge. Grouping by
    /// chromosome alone keeps them; `prune_same_locus` clause (c) still removes SAME-LOCUS antisense overlap,
    /// which was the only thing the strand split was actually protecting against.
    #[test]
    fn colocated_families_keeps_inverted_duplicate_pair() {
        let mut plus = rep_s(1_000, 5_000, vec![(2_000, 3_000)], 20);
        plus.strand = '+';
        let mut minus = rep_s(50_000, 54_000, vec![(51_000, 52_000)], 18); // disjoint locus, opposite strand
        minus.strand = '-';
        let reps = vec![plus, minus];
        let stats = crate::vg_family::family_split::CommunityStats { n: 2, n_edges: 1, density: 1.0, avg_core_recip: 1.0, n_articulation: 0 };
        let fams = vec![SplitFamily { members: vec![0, 1], class: FamilyClass::Family, stats }];

        let out = colocated_families(&reps, &fams, 5_000_000, 2, &DetectParams::default());
        assert_eq!(out.len(), 1, "an inverted duplicate pair is ONE co-located family");
        assert_eq!(out[0].copies.len(), 2, "both copies survive: {:?}", out[0].copies.iter().map(|c| c.strand).collect::<Vec<_>>());
        let strands: std::collections::BTreeSet<char> = out[0].copies.iter().map(|c| c.strand).collect();
        assert_eq!(strands, ['+', '-'].into_iter().collect(), "the family spans both strands");
    }

    /// The protection the strand split was really providing must survive: a `+` gene and a `-` gene whose
    /// spans OVERLAP are an inverted-repeat/antisense artifact at one locus, not two copies.
    #[test]
    fn colocated_families_still_drops_same_locus_antisense_overlap() {
        let mut plus = rep_s(1_000, 11_000, vec![(2_000, 3_000), (4_000, 5_000)], 15);
        plus.strand = '+';
        let mut minus = rep_s(1_100, 10_900, vec![(6_000, 7_000), (8_000, 9_000)], 5); // overlaps, opposite strand
        minus.strand = '-';
        let reps = vec![plus, minus];
        let stats = crate::vg_family::family_split::CommunityStats { n: 2, n_edges: 1, density: 1.0, avg_core_recip: 1.0, n_articulation: 0 };
        let fams = vec![SplitFamily { members: vec![0, 1], class: FamilyClass::Family, stats }];

        let out = colocated_families(&reps, &fams, 5_000_000, 2, &DetectParams::default());
        assert!(out.is_empty(), "antisense overlap collapses to ONE copy, which is below min_copies=2");
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
    fn protein_coheres_is_none_without_mmseqs() {
        // With protein_tail off / mmseqs absent the flag is None (no membership effect).
        let fam = vec![vec![
            DenovoTranscript{tid:"a".into(),chrom:"c1".into(),start:0,end:300,n_reads:5,strand:'+',introns:vec![],seq:b"ATGAAAGGGTTTTGTCCCAAAGGG".to_vec()},
            DenovoTranscript{tid:"b".into(),chrom:"c1".into(),start:9000,end:9300,n_reads:4,strand:'+',introns:vec![],seq:b"ATGAAAGGGTTTTGTCCCAAAGGG".to_vec()},
        ]];
        let mut p = RefineParams::default(); p.protein_tail = false;
        let flags = family_protein_coheres(&fam, &p);
        assert_eq!(flags.len(), 1);
        assert_eq!(flags[0], None, "protein QC is None when protein tier is off");
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
    fn colocated_families_keeps_mixed_strand_disjoint_copies_in_one_family() {
        use super::super::family_split::{CommunityStats, FamilyClass, SplitFamily};
        // 4 disjoint-locus copies on c1: two '+' and two '-' (distinct junctions, no containment), all in
        // ONE conflict family. They are four copies of one inverted-duplication family, so they must stay
        // ONE family. This test previously asserted a split into two same-strand families; that split is what
        // made O2 blind to MAGEA4/MAGEA10 and RFPL2/RFPL3. Same-locus antisense is still removed — by
        // `prune_same_locus` clause (c) — and is covered by
        // `colocated_families_still_drops_same_locus_antisense_overlap`.
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
        assert_eq!(colo.len(), 1, "disjoint copies on one chromosome are ONE family regardless of strand");
        assert_eq!(colo[0].copies.len(), 4, "all four disjoint copies are kept");
        let strands: std::collections::HashSet<char> = colo[0].copies.iter().map(|x| x.strand).collect();
        assert_eq!(strands.len(), 2, "the family spans both strands (inverted duplication)");
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
    fn to_split_families_matches_conflict_wrapper() {
        // Two conflict components, edges with distinct usize weights spanning both.
        let families = vec![vec![0usize, 1, 2], vec![3usize, 4]];
        let c_edges = vec![(0usize, 1usize, 5usize), (1, 2, 3), (3, 4, 2)];
        let p = SplitParams::default();

        let via_wrapper = conflict_to_split_families(&families, &c_edges, &p);

        let float_edges: Vec<(usize, usize, f64)> =
            c_edges.iter().map(|&(a, b, w)| (a, b, w as f64)).collect();
        let via_core = to_split_families(&families, &float_edges, &p);

        assert_eq!(via_wrapper.len(), via_core.len());
        for (w, c) in via_wrapper.iter().zip(via_core.iter()) {
            assert_eq!(w.members, c.members);
            assert_eq!(w.class, c.class);
            assert_eq!(w.stats, c.stats);
        }
    }

    /// The real CAFAM69 geometry: a 65 kb de-novo transcript with a shorter one nested inside it, admitted
    /// as two "copies" of one family. Every read scored `min_p == 1` against that pair, so all 1043 reads
    /// abstained and looked like the K=0 wall. Disjoint loci (the true paralog case) must NOT be flagged.
    #[test]
    fn copies_overlap_flags_nested_transcripts_not_disjoint_paralogs() {
        let nested = vec![
            ("NC_073247.2".to_string(), 164_381_269, 164_446_025),
            ("NC_073247.2".to_string(), 164_430_240, 164_436_493),
        ];
        let got = copies_overlap(&nested);
        assert_eq!(got.len(), 1);
        assert_eq!((got[0].0, got[0].1), (0, 1));
        assert!(got[0].2 < 0.2, "a 65 kb span containing a 6 kb one is CONTAINMENT, not a duplicate: {}", got[0].2);

        // boundary wobble: the real CAFAM69 pair -- reciprocal ~1.0 => one locus admitted twice
        let wobble = vec![
            ("NC_073247.2".to_string(), 164_381_222, 164_384_848),
            ("NC_073247.2".to_string(), 164_381_237, 164_384_845),
        ];
        let w = copies_overlap(&wobble);
        assert_eq!(w.len(), 1);
        assert!(w[0].2 > 0.99, "15 bp of wobble => reciprocal overlap ~1.0, got {}", w[0].2);

        // CAFAM20: two genuinely distinct loci 43 kb apart, exonically identical (a real K=0 family).
        let disjoint = vec![
            ("NC_073229.2".to_string(), 136_502_891, 136_507_428),
            ("NC_073229.2".to_string(), 136_550_623, 136_557_225),
        ];
        assert!(copies_overlap(&disjoint).is_empty(), "disjoint paralogs are a real family");

        // Same coordinates on different contigs cannot overlap.
        let cross = vec![("c1".to_string(), 100, 200), ("c2".to_string(), 100, 200)];
        assert!(copies_overlap(&cross).is_empty());

        // Abutting spans (end == start) touch but do not overlap.
        let abut = vec![("c1".to_string(), 100, 200), ("c1".to_string(), 200, 300)];
        assert!(copies_overlap(&abut).is_empty(), "half-open spans: end == start is not an overlap");
    }

    /// The real GSTM catalog: CAFAM1 correctly holds GSTM5 and GSTM1 as two disjoint copies, while CAFAM0
    /// holds GSTM3 plus a 30 kb single-intron readthrough transcript spanning BOTH of them. The readthrough
    /// shares sequence with two copies of another family — a defect no within-family check can see.
    #[test]
    fn catalog_overlaps_flags_readthrough_shared_across_families() {
        let c = |f: &str, s: u64, e: u64| (f.to_string(), "NC_073224.2".to_string(), s, e);
        let copies = vec![
            c("CAFAM0", 129_169_623, 129_173_090), // GSTM3
            c("CAFAM0", 129_190_708, 129_220_537), // readthrough over GSTM5+GSTM1
            c("CAFAM1", 129_191_737, 129_198_214), // GSTM5
            c("CAFAM1", 129_211_328, 129_222_730), // GSTM1
        ];
        let found = catalog_overlaps(&copies);
        let shared: Vec<_> =
            found.iter().filter(|f| f.3 == OverlapKind::SharedAcrossFamilies).map(|f| (f.0, f.1)).collect();
        assert_eq!(shared, vec![(1, 2), (1, 3)], "the readthrough must be flagged against BOTH real copies");
        // GSTM3 is disjoint, and GSTM5/GSTM1 are disjoint from each other: no other pair.
        assert_eq!(found.len(), 2, "no spurious pairs: {found:?}");
    }

    #[test]
    fn catalog_overlaps_separates_duplicate_locus_from_containment() {
        let c = |f: &str, s: u64, e: u64| (f.to_string(), "c1".to_string(), s, e);
        // wobble: same interval twice, one family
        let dup = vec![c("F", 164_381_222, 164_384_848), c("F", 164_381_237, 164_384_845)];
        assert_eq!(catalog_overlaps(&dup)[0].3, OverlapKind::DuplicateLocus);
        // containment: a long readthrough enclosing a fragment, one family
        let con = vec![c("F", 20_202_981, 20_390_694), c("F", 20_243_056, 20_379_683)];
        let g = catalog_overlaps(&con);
        assert_eq!(g[0].3, OverlapKind::Containment);
        assert!(g[0].2 < 0.9);
        // disjoint copies on the same contig produce nothing
        let ok = vec![c("F", 100, 200), c("F", 300, 400)];
        assert!(catalog_overlaps(&ok).is_empty());
    }

    fn junc(list: &[(u64, u64, usize)]) -> std::collections::HashMap<(String, u64, u64), usize> {
        list.iter().map(|&(d, a, n)| (("c1".to_string(), d, a), n)).collect()
    }

    /// The GSTM / DAZ / BPY2 shape: a single-exon span engulfing many distinct junctions is unspliced pre-mRNA.
    #[test]
    fn readthrough_single_exon_engulfing_many_junctions_is_flagged() {
        let mut rt = rep_s(1_000, 100_000, vec![], 12);
        rt.introns.clear();
        let j = junc(&[(10_000, 11_000, 5), (20_000, 21_000, 4), (30_000, 31_000, 9),
                       (40_000, 41_000, 3), (50_000, 51_000, 7)]);
        assert!(is_unspliced_readthrough(&rt, &j, READTHROUGH_MIN_SUPPORT, READTHROUGH_MIN_DISTINCT));
    }

    /// A SPLICED transcript is a transcript model, never a readthrough — the rule must not touch it.
    #[test]
    fn readthrough_never_flags_a_spliced_transcript() {
        let spliced = rep_s(1_000, 100_000, vec![(10_000, 11_000), (20_000, 21_000)], 40);
        let j = junc(&[(10_000, 11_000, 5), (20_000, 21_000, 4), (30_000, 31_000, 9),
                       (40_000, 41_000, 3), (50_000, 51_000, 7)]);
        assert!(!is_unspliced_readthrough(&spliced, &j, READTHROUGH_MIN_SUPPORT, READTHROUGH_MIN_DISTINCT));
    }

    /// TSPYL1: a real intronless gene with a handful of stray junctions inside its span (4 distinct) — KEPT.
    /// This is the closest control; the margin to the artifact minimum (14) is three junctions.
    #[test]
    fn readthrough_keeps_real_intronless_gene_with_stray_junctions() {
        let mut gene = rep_s(1_000, 15_000, vec![], 2080);
        gene.introns.clear();
        let j = junc(&[(2_000, 2_500, 3), (3_000, 3_500, 2), (4_000, 4_500, 5), (5_000, 5_500, 4)]);
        assert!(!is_unspliced_readthrough(&gene, &j, READTHROUGH_MIN_SUPPORT, READTHROUGH_MIN_DISTINCT));
    }

    /// A giant 115 kb intron carried by only 2 reads (below the locus gate) is a spurious splice — flagged.
    #[test]
    fn mischain_flags_sub_gate_giant_intron() {
        let sub_gate = rep_s(92_076_242, 92_198_994, vec![(92_077_601, 92_193_251)], 5); // 115.6 kb intron
        let j = junc(&[(92_077_601, 92_193_251, 2)]); // only 2 reads cross the giant intron -> below gate
        assert!(is_giant_intron_mischain(&sub_gate, &j, MISCHAIN_GIANT_INTRON_BP, MISCHAIN_MIN_JUNCTION_READS));
    }

    /// A real large-gene paralog (POTE): its 48 kb intron is under the giant threshold AND well-supported (7
    /// reads) — never flagged. Both properties must hold for the rule to be safe.
    #[test]
    fn mischain_keeps_real_large_gene_with_well_supported_intron() {
        // 48 kb intron < 50 kb giant threshold: not even examined.
        let pote = rep_s(32_266_748, 32_352_771, vec![(32_280_000, 32_328_000)], 6);
        let j = junc(&[(32_280_000, 32_328_000, 7)]);
        assert!(!is_giant_intron_mischain(&pote, &j, MISCHAIN_GIANT_INTRON_BP, MISCHAIN_MIN_JUNCTION_READS));
        // and even a GIANT (60 kb) intron is kept when well-supported (>= 3 reads) — a real deep large gene.
        let big = rep_s(0, 200_000, vec![(20_000, 80_000)], 643);
        let jb = junc(&[(20_000, 80_000, 643)]);
        assert!(!is_giant_intron_mischain(&big, &jb, MISCHAIN_GIANT_INTRON_BP, MISCHAIN_MIN_JUNCTION_READS));
    }

    /// retain_non_mischain drops the mis-chain and keeps the real large gene.
    #[test]
    fn retain_non_mischain_drops_only_the_mischain() {
        let mut txs = vec![
            rep_s(92_076_242, 92_198_994, vec![(92_077_601, 92_193_251)], 5), // mis-chain: 115kb / 2 reads
            rep_s(0, 200_000, vec![(20_000, 80_000)], 643),                    // real: 60kb / 643 reads
        ];
        let j = junc(&[(92_077_601, 92_193_251, 2), (20_000, 80_000, 643)]);
        retain_non_mischain(&mut txs, &j, "test");
        assert_eq!(txs.len(), 1);
        assert_eq!(txs[0].start, 0, "the well-supported large gene survives; the mis-chain is dropped");
    }

    /// EEF1A1 -> LOC109023808: a retrocopy has NO introns, so the spliced parent's cross-mapping reads align
    /// contiguously and deposit no junction inside it. Measured: 0 distinct junctions, 15 reads.
    #[test]
    fn readthrough_keeps_an_intronless_retrocopy() {
        let mut retro = rep_s(97_380_144, 97_381_766, vec![], 15);
        retro.introns.clear();
        assert!(!is_unspliced_readthrough(&retro, &junc(&[]), READTHROUGH_MIN_SUPPORT, READTHROUGH_MIN_DISTINCT));
    }

    /// A junction must lie ENTIRELY inside the span. A gene nested in another gene's intron sees the host's
    /// junctions FLANKING it, so containment — not overlap — is what protects it.
    #[test]
    fn readthrough_ignores_junctions_that_merely_overlap_the_span() {
        let mut nested = rep_s(10_000, 11_000, vec![], 50);
        nested.introns.clear();
        // the host gene's intron spans the nested gene entirely: donor before, acceptor after.
        let j = junc(&[(9_000, 12_000, 90), (8_000, 13_000, 80), (7_000, 14_000, 70),
                       (6_000, 15_000, 60), (5_000, 16_000, 50), (4_000, 17_000, 40)]);
        assert!(!is_unspliced_readthrough(&nested, &j, READTHROUGH_MIN_SUPPORT, READTHROUGH_MIN_DISTINCT));
    }

    /// Junctions below the support floor are alignment noise, not splicing evidence.
    #[test]
    fn readthrough_ignores_singleton_junctions() {
        let mut rt = rep_s(1_000, 100_000, vec![], 12);
        rt.introns.clear();
        let j = junc(&[(10_000, 11_000, 1), (20_000, 21_000, 1), (30_000, 31_000, 1),
                       (40_000, 41_000, 1), (50_000, 51_000, 1), (60_000, 61_000, 1)]);
        assert!(!is_unspliced_readthrough(&rt, &j, READTHROUGH_MIN_SUPPORT, READTHROUGH_MIN_DISTINCT));
    }

    /// A locus contributes ambiguity only from PRIMARY reads that overlap it. Supplementary records are
    /// adjacency, not ambiguity; reads outside the rep belong to another locus.
    #[test]
    fn locus_ambiguity_counts_only_mapq0_primaries_overlapping_the_rep() {
        let rep = rep_s(1_000, 2_000, vec![], 10);
        let mk = |start: u64, mapq: u8, supp: bool| BamRead {
            chrom: "c1".into(),
            read: AlignedRead {
                ref_start: start,
                cigar: vec![('M', 100)],
                seq: vec![b'A'; 100],
                qual: vec![30; 100],
            },
            mapq,
            name: format!("r{start}_{mapq}"),
            as_score: 100,
            de: 0.01,
            is_supplementary: supp,
        };
        let reads = vec![
            mk(1_100, 0, false),  // ambiguous, inside  -> numerator + denominator
            mk(1_200, 60, false), // unique, inside     -> denominator only
            mk(1_300, 0, true),   // supplementary      -> ignored
            mk(9_000, 0, false),  // outside the rep    -> ignored
        ];
        assert_eq!(
            locus_ambiguity(&rep, &reads),
            crate::vg_family::collapse_gate::Ambiguity { n: 2, k: 1 }
        );
    }

    /// OFF by default: the ambiguity instrument detects unresolvable paralogy rather than collapse, and fires
    /// on the single-copy control EEF1A1 with chi(H) = 7. Enabling it is an explicit, informed choice.
    #[test]
    fn single_copy_baseline_excludes_family_copies() {
        use crate::vg_family::single_copy::single_copy_loci;
        let mk = |tid: &str, start: u64| DenovoTranscript {
            tid: tid.into(), chrom: "c1".into(), start, end: start + 100, n_reads: 12, strand: '+',
            introns: vec![], seq: vec![],
        };
        let (a, b, c) = (mk("a", 0), mk("b", 1000), mk("c", 2000));
        let families = vec![ColocatedFamily {
            family_id: "F0".into(), chrom: "c1".into(), start: 0, end: 1100, copies: vec![a.clone(), b.clone()],
        }];
        let sc = single_copy_loci(&[a, b, c], &[12, 12, 12], &families);
        assert_eq!(sc.len(), 1);
        assert_eq!(sc[0].start, 2000);
    }

    #[test]
    fn denovoconfig_default_disables_the_collapse_gate() {
        assert!(!DenovoConfig::default().collapse_gate);
    }

    /// A gated family reports chi(H) copies but only ONE copy tid (the single assembled rep), every read Tied
    /// with min_p = 1.0, and `collapsed_copies` set so it can never be mistaken for an assembled family.
    #[test]
    fn gated_family_ties_every_read_and_materialises_no_copy_sequence() {
        let rep = rep_s(1_000, 2_000, vec![], 10);
        let mk = |start: u64| BamRead {
            chrom: "c1".into(),
            read: AlignedRead { ref_start: start, cigar: vec![('M', 100)], seq: vec![b'A'; 100], qual: vec![30; 100] },
            mapq: 0,
            name: format!("r{start}"),
            as_score: 100,
            de: 0.01,
            is_supplementary: false,
        };
        let reads = vec![mk(1_100), mk(1_200), mk(9_000)];
        let fa = gated_family(&rep, &reads, 3, "DSFAM7".into());
        assert_eq!(fa.n_copies, 3, "n_copies is chi(H)");
        assert_eq!(fa.collapsed_copies, 3, "records how the count was obtained");
        assert_eq!(fa.copy_tids.len(), 1, "only the assembled rep has a sequence");
        assert_eq!(fa.n_reads, 2, "only the two overlapping reads");
        assert!(fa.copy_psv_alleles.is_empty(), "no per-copy consensus is materialised");
        for (_, a) in &fa.assignments {
            assert_eq!(a.status, AssignStatus::Tied);
            assert_eq!(a.min_p_value, 1.0, "no distinguishing column exists: that is the honest certificate");
        }
    }

    /// Junction support is PRIMARY-only by type: `PrimaryRead` cannot hold a secondary alignment
    /// (`primary_read_from_record` returns `None` for them). The predecessor took `&[BamRead]` and counted
    /// secondaries, inflating the DAZ readthrough span from 56 to 154 distinct junctions.
    #[test]
    fn read_junction_support_counts_each_primary_read_once() {
        let pr = |s: u64, introns: Vec<(u64, u64)>| PrimaryRead { chrom: "c1".into(), ref_start: s, ref_end: s + 400, introns };
        let reads = vec![pr(50, vec![(100, 200)]), pr(60, vec![(100, 200)]), pr(70, vec![(300, 400)])];
        let sup = read_junction_support(&reads);
        assert_eq!(sup.get(&("c1".to_string(), 100, 200)), Some(&2));
        assert_eq!(sup.get(&("c1".to_string(), 300, 400)), Some(&1));
        assert_eq!(sup.len(), 2);
    }

    /// The shared helper O1 and O2 now both call: a single-exon transcript engulfing >= 5 supported junctions
    /// is dropped; a spliced transcript over the same span is kept.
    #[test]
    fn retain_non_readthrough_drops_engulfing_single_exon_keeps_spliced() {
        let mut support = std::collections::HashMap::new();
        for i in 1..=6u64 {
            support.insert(("c1".to_string(), i * 100, i * 100 + 50), 2usize);
        }
        let mk = |tid: &str, introns: Vec<(u64, u64)>| DenovoTranscript {
            tid: tid.into(), chrom: "c1".into(), start: 0, end: 1000, n_reads: 5, strand: '+', introns,
            seq: vec![b'A'; 50],
        };
        let mut ts = vec![mk("readthrough", vec![]), mk("spliced", vec![(100, 150)])];
        retain_non_readthrough(&mut ts, &support, "test");
        assert_eq!(ts.len(), 1, "the single-exon engulfer is dropped");
        assert_eq!(ts[0].tid, "spliced");
    }

    /// One canonical extent per locus across objectives. `detect_families` and `detect_and_assign` used a
    /// hardcoded k = 1 while the three genome-wide O1 catalogs trimmed at `cfg.min_terminal_support` = 2, so
    /// the same locus had different boundaries in O1 and O2.
    #[test]
    fn skeleton_terminal_trim_is_uniform_k() {
        let mk = |s: u64| PrimaryRead { chrom: "c1".into(), ref_start: s, ref_end: 500, introns: vec![(200, 300)] };
        let reads = vec![mk(100), mk(100), mk(100), mk(10)]; // one runaway 5' read
        assert_eq!(pass1_skeletons(&reads, 2)[0].start, 10, "k=1 keeps the runaway");
        assert_eq!(
            pass1_skeletons_robust(&reads, 2, DenovoConfig::default().min_terminal_support)[0].start,
            100,
            "the k every path now uses trims it"
        );
    }

    #[test]
    fn denovoconfig_default_homology_primary_is_off() {
        assert!(!DenovoConfig::default().homology_primary, "default must be byte-identical to today (E_c oracle)");
    }

    /// Task 2 (O1<->O2 harmony): `detect_and_assign`'s `cfg.homology_primary == false` branch composes
    /// `conflict_edges -> conflict_families -> (map to f64) -> to_split_families`. This must be identical
    /// (per-family members/class/stats) to calling the pre-existing `conflict_to_split_families` wrapper
    /// directly on the SAME `(c_fams, c_edges)` — i.e. the routing change in `detect_and_assign` is
    /// behavior-preserving for the default (off) path.
    #[test]
    fn homology_primary_off_leaves_conflict_path_identical() {
        let n = 5;
        let c_edges: Vec<(usize, usize, usize)> = vec![(0, 1, 5), (1, 2, 3), (3, 4, 2)];
        let c_fams = conflict_families(n, &c_edges);
        let p = SplitParams::default();

        // exactly the `!cfg.homology_primary` branch body in `detect_and_assign`:
        let edges_f64: Vec<(usize, usize, f64)> = c_edges.iter().map(|&(a, b, w)| (a, b, w as f64)).collect();
        let via_off_branch = to_split_families(&c_fams, &edges_f64, &p);

        // the pre-existing path:
        let via_existing = conflict_to_split_families(&c_fams, &c_edges, &p);

        assert_eq!(via_off_branch.len(), via_existing.len());
        for (a, b) in via_off_branch.iter().zip(via_existing.iter()) {
            assert_eq!(a.members, b.members);
            assert_eq!(a.class, b.class);
            assert_eq!(a.stats, b.stats);
        }
    }

    /// Task 2's whole point: a copy pair whose reads map UNIQUELY (no conflict edge at all — the exact
    /// shape of the drop bug) is invisible to the E_c oracle (`conflict_families` drops singletons
    /// entirely, so BOTH copies vanish, not merely remain unlinked) but IS unioned by the E_r homology
    /// oracle once an edge exists between them.
    #[test]
    fn homology_oracle_unions_uniquely_mappable_pair() {
        // E_c: no conflict edge between the two loci at all -> conflict_families finds NO family (both
        // copies would vanish from the family roster, the exact drop this task fixes).
        assert!(conflict_families(2, &[]).is_empty(), "no conflict edge => E_c forms no family at all");

        // The partition step (`gamma_quasi_clique_partition`) unions a homology edge into one family. This
        // is proven with a SYNTHETIC edge first (independent of minimap2 availability), then confirmed with
        // the real minimap2-backed `homology_edges_all_reps` when the binary is present.
        let synthetic_edges_f64: Vec<(usize, usize, f64)> = vec![(0, 1, 1.0)];
        let synthetic_families =
            crate::vg_family::family_split::gamma_quasi_clique_partition(2, &synthetic_edges_f64, 0.20);
        assert_eq!(synthetic_families, vec![vec![0, 1]], "homology edge must union both loci into ONE family");

        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            eprintln!("minimap2 absent; skipping the real homology_edges_all_reps confirmation");
            return;
        }
        // Two reps, near-identical (~4% divergence, so minimap2 finds a genuine homology edge) but placed
        // on DIFFERENT chroms so no read could ever place ambiguously between them -- a uniquely-mappable
        // paralog pair, i.e. exactly the copies read-conflict (E_c) can never link.
        let base = rand_seq(900, 0x22);
        let mut para = base.clone();
        for k in (0..para.len()).step_by(25) { para[k] = b"ACGT"[(para[k] as usize + 1) % 4]; }
        let reps = vec![
            DenovoTranscript { tid: "u0".into(), chrom: "c1".into(), start: 0, end: 900, n_reads: 10, strand: '+', introns: vec![], seq: base },
            DenovoTranscript { tid: "u1".into(), chrom: "c2".into(), start: 0, end: 900, n_reads: 10, strand: '+', introns: vec![], seq: para },
        ];
        let params = RefineParams::default();
        let edges = homology_edges_all_reps(&reps, &params).unwrap();
        assert!(edges.contains(&(0, 1)), "near-identical uniquely-mapping pair must be E_r-linked, got {:?}", edges);
        let edges_f64: Vec<(usize, usize, f64)> = edges.iter().map(|&(a, b)| (a, b, 1.0)).collect();
        let families = crate::vg_family::family_split::gamma_quasi_clique_partition(2, &edges_f64, 0.20);
        assert_eq!(families, vec![vec![0, 1]], "the real homology oracle must union the uniquely-mappable pair");
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

    // --- codon-aware fixture builder for the protein-tail rescue test ------------------------------------
    // The standard genetic code as (amino acid, codon) pairs, excluding the 3 stop codons (TAA/TAG/TGA).
    // Used only to build a nt-divergent / protein-conserved pair: a codon-based "random protein" plus a
    // synonymously-recoded (+ conservative-substitution) sibling whose nucleotide sequence diverges hard
    // while the encoded protein stays near-identical.
    const CODON_TABLE: &[(u8, &[u8; 3])] = &[
        (b'F', b"TTT"), (b'F', b"TTC"),
        (b'L', b"TTA"), (b'L', b"TTG"), (b'L', b"CTT"), (b'L', b"CTC"), (b'L', b"CTA"), (b'L', b"CTG"),
        (b'I', b"ATT"), (b'I', b"ATC"), (b'I', b"ATA"),
        (b'M', b"ATG"),
        (b'V', b"GTT"), (b'V', b"GTC"), (b'V', b"GTA"), (b'V', b"GTG"),
        (b'S', b"TCT"), (b'S', b"TCC"), (b'S', b"TCA"), (b'S', b"TCG"), (b'S', b"AGT"), (b'S', b"AGC"),
        (b'P', b"CCT"), (b'P', b"CCC"), (b'P', b"CCA"), (b'P', b"CCG"),
        (b'T', b"ACT"), (b'T', b"ACC"), (b'T', b"ACA"), (b'T', b"ACG"),
        (b'A', b"GCT"), (b'A', b"GCC"), (b'A', b"GCA"), (b'A', b"GCG"),
        (b'Y', b"TAT"), (b'Y', b"TAC"),
        (b'H', b"CAT"), (b'H', b"CAC"),
        (b'Q', b"CAA"), (b'Q', b"CAG"),
        (b'N', b"AAT"), (b'N', b"AAC"),
        (b'K', b"AAA"), (b'K', b"AAG"),
        (b'D', b"GAT"), (b'D', b"GAC"),
        (b'E', b"GAA"), (b'E', b"GAG"),
        (b'C', b"TGT"), (b'C', b"TGC"),
        (b'W', b"TGG"),
        (b'R', b"CGT"), (b'R', b"CGC"), (b'R', b"CGA"), (b'R', b"CGG"), (b'R', b"AGA"), (b'R', b"AGG"),
        (b'G', b"GGT"), (b'G', b"GGC"), (b'G', b"GGA"), (b'G', b"GGG"),
    ];
    fn codons_for(aa: u8) -> Vec<&'static [u8; 3]> {
        CODON_TABLE.iter().filter(|(a, _)| *a == aa).map(|(_, c)| *c).collect()
    }
    fn all_aas() -> Vec<u8> {
        let mut v: Vec<u8> = CODON_TABLE.iter().map(|(a, _)| *a).collect();
        v.sort();
        v.dedup();
        v
    }
    fn hamming3(a: &[u8; 3], b: &[u8; 3]) -> usize {
        (0..3).filter(|&i| a[i] != b[i]).count()
    }
    /// The codon for `aa` that differs MAXIMALLY (Hamming, over the 3 nt positions) from `from`.
    fn max_divergent_codon(aa: u8, from: &[u8; 3]) -> [u8; 3] {
        let opts = codons_for(aa);
        **opts.iter().max_by_key(|c| hamming3(c, from)).expect("every table aa has >=1 codon")
    }
    /// A conservative (same biochemical class) amino-acid substitution group.
    fn conservative_group(aa: u8) -> &'static [u8] {
        match aa {
            b'A' | b'V' | b'L' | b'I' | b'M' => b"AVLIM",
            b'F' | b'Y' | b'W' => b"FYW",
            b'S' | b'T' | b'N' | b'Q' => b"STNQ",
            b'D' | b'E' => b"DE",
            b'K' | b'R' | b'H' => b"KRH",
            _ => b"GPC", // G, P, C (and any fallback)
        }
    }

    #[test]
    fn homology_edges_protein_rescues_nt_divergent_pair() {
        if std::process::Command::new("mmseqs").arg("version").output().is_err() {
            eprintln!("mmseqs absent; skipping (protein-tail rescue NOT verified in this environment)");
            return;
        }
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            eprintln!("minimap2 absent; skipping");
            return;
        }
        // Build a ~250-codon (750 nt) random coding sequence (seq_a), and a sibling (seq_b) that: (i) at
        // ~85% of codons keeps the SAME amino acid but recodes to the synonymous codon MAXIMALLY divergent
        // from seq_a's (wobble-position + codon-family divergence), and (ii) at ~15% of codons substitutes a
        // conservative (same biochemical class) DIFFERENT amino acid, again maximally-divergent-coded. This
        // pushes nucleotide identity hard down (both nt tiers must miss it) while the encoded protein stays
        // near-identical (well above mmseqs' fident>=0.50 gate) — exactly the "diverged past the nucleotide
        // floor, still protein-conserved" coding-paralog case the protein DEFINITION edge exists to rescue.
        let n_codons = 250;
        let mut rng = SplitMix64(0xC0DE_C0DE_1234_5678);
        let aas = all_aas();
        let mut seq_a: Vec<u8> = Vec::with_capacity(n_codons * 3);
        let mut seq_b: Vec<u8> = Vec::with_capacity(n_codons * 3);
        for i in 0..n_codons {
            let aa_a = aas[(rng.next_u64() as usize) % aas.len()];
            let opts_a = codons_for(aa_a);
            let codon_a = *opts_a[(rng.next_u64() as usize) % opts_a.len()];
            seq_a.extend_from_slice(&codon_a);

            let aa_b = if i % 7 == 0 {
                let group = conservative_group(aa_a);
                let alt: Vec<u8> = group.iter().copied().filter(|&a| a != aa_a).collect();
                if alt.is_empty() { aa_a } else { alt[(rng.next_u64() as usize) % alt.len()] }
            } else {
                aa_a
            };
            let codon_b = max_divergent_codon(aa_b, &codon_a);
            seq_b.extend_from_slice(&codon_b);
        }
        assert_eq!(seq_a.len(), 750);
        assert_eq!(seq_b.len(), 750);

        // Proxy identity (no indels by construction -> same length, position-wise Hamming) for reporting.
        let same = seq_a.iter().zip(seq_b.iter()).filter(|(x, y)| x == y).count();
        let proxy_identity = same as f64 / seq_a.len() as f64;
        eprintln!("[test] constructed pair proxy (position-wise) nt identity = {proxy_identity:.3}");
        assert!(proxy_identity < 0.60, "RNG did not produce enough nt divergence: proxy identity {proxy_identity:.3}");

        // Confirm via the pipeline's own aligner: BOTH nt tiers must find NO edge at this divergence.
        let p = RefineParams::default();
        let sensitive_edges = nucleotide_edges(
            &[seq_a.clone(), seq_b.clone()], &["-k", "11", "-w", "5"], 0.60, 0.50, &p,
        ).unwrap();
        assert!(sensitive_edges.is_empty(), "pair must be genuinely nt-unresolvable (< 0.60), got {:?}", sensitive_edges);

        let reps = vec![
            DenovoTranscript { tid: "protA".into(), chrom: "c1".into(), start: 1_000, end: 1_750, n_reads: 20, strand: '+', introns: vec![], seq: seq_a },
            DenovoTranscript { tid: "protB".into(), chrom: "c1".into(), start: 50_000, end: 50_750, n_reads: 18, strand: '+', introns: vec![], seq: seq_b },
        ];

        let mut params_off = RefineParams::default();
        params_off.protein_tail = false;
        let edges_off = homology_edges_all_reps(&reps, &params_off).unwrap();
        assert!(!edges_off.contains(&(0, 1)), "without protein_tail the nt-divergent pair must NOT link, got {:?}", edges_off);

        let mut params_on = RefineParams::default();
        params_on.protein_tail = true;
        let edges_on = homology_edges_all_reps(&reps, &params_on).unwrap();
        assert!(edges_on.contains(&(0, 1)), "protein_tail=true must RESCUE the nt-divergent, protein-conserved pair, got {:?}", edges_on);
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
