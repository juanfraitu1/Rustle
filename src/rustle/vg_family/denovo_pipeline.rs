//! De-novo family detection DRIVER (integration stage 3): chains the ported cores into a family roster.
//!
//! **Terminology.** "Locus" in this module is the **gene-locus** sense unless
//! explicitly noted: a set of isoforms collapsed by shared splice junctions
//! (`family_detect::collapse_loci`). The physical `(chrom, start, end)` span
//! used for the ≥2-distinct-loci certificate is `family_definition::distinct_loci`.
//! See `docs/REFERENCE.md` for the canonical vocabulary.
//!
//!   primary reads ─► pass1 skeletons ─► assemble gate ─► collapse loci ─► detect edges ─► decompose families
//!
//! This is the read-coherence-way detection pipeline end to end (rescue + per-read copy assignment are the
//! next stage). `detect_families` is the testable transform over parsed reads + a loaded genome (a
//! single-region reference path; the genome-wide catalogs are the shipped entry points).

use std::collections::{BTreeSet, HashSet, BTreeMap};

use anyhow::Result;

use super::absent_copy::{self, AbsentCopyParams, Admission, DnaNeedsRecord};
use super::copy_assign::{
    assign_read_editing, Assignment, AssignParams, AssignStatus, BubbleGraph, CopyProfile, ReadFeatures,
};
use super::copy_assign_pipeline::{
    assign_family_detailed, assign_family_detailed_pruned, best_overlap_copy, build_family_profiles,
    copy_boundaries, detect_editing_columns, freeze_merge, gen2off, read_ref_end, FamilyProfiles,
};
use super::family_graph::contiguous_core_coverage_bounded;
use super::copy_split::{
    split_locus_copies, discover_locus_psvs, AlignedRead, CollapsedCandidate, CopyIsoform,
};
use super::denovo_assemble::{
    aligned_reads_from_bam, assemble_gate, pass1_skeletons, pass1_skeletons_robust, primary_reads_from_bam,
    reads_in_region, split_mischained_reads, tied_seed_skeletons, BamRead, GateParams, PrimaryRead,
    GATE_MIN_READS, PASS1_MIN_READS,
};
use super::family_detect::{
    collapse_loci_span_aware, collapse_loci_span_aware_with_totals, detect_edges, detect_edges_reporting,
    DenovoTranscript, DetectParams,
};
use super::family_rescue::{FamilyMember, RescueParams};
use super::family_split::{classify, community_stats, decompose_families, FamilyClass, SplitFamily, SplitParams};
use super::read_conflict::{
    as_tie_edges, conflict_edges, conflict_families, family_mapq0_support, locus_unique_mapper_counts,
    reads_distinguish, ConflictParams, Placement, ReadPlacements,
};
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
    /// VG re-align supplement (opt-in): default false = OFF, byte-identical. When true,
    /// `detect_and_assign` runs `vg_realign::apply_realign` over each co-located family's reads and then
    /// FEEDS THE RESULT BACK into the emitted `FamilyAssignment`: `apply_realign_patch` corrects per-read
    /// VG re-align CORRECTION leg (opt-in): re-thread each poor-fit/candidate read through the family's
    /// copy-paths, take the best-fitting path (same ε^Δ significance certificate as the PSV gate), and CORRECT
    /// its per-read copy assignment; then `recompute_realign_abundance` re-weights the EM `copy_abundance`.
    /// Per-read decisions recorded on `realign_records`. OFF (default) leaves every emitted field
    /// byte-identical. The roster-widening novel-copy ADMISSION is the SEPARATE `vg_realign_admit` leg below.
    pub vg_realign: bool,
    /// VG re-align ADMISSION leg (opt-in, gated SEPARATELY from the correction leg because it touches the
    /// GENOME and WIDENS the copy roster on the unvalidated O4-divergent frontier = real FP-copy risk):
    /// `admit_novel_pools` clusters reads that fit NO existing copy and admits them as new copies. Requires
    /// `vg_realign` (the correction leg) to also be on. OFF (default) => no roster widening.
    pub vg_realign_admit: bool,
    /// E_r homology-primary family MEMBERSHIP (opt-in). Conflict/PSV/χ(H) remain within-family. Enlarges
    /// the copy set ⟹ stricter Bonferroni α/(K−1) ⟹ assignments shift. Requires minimap2.
    pub homology_primary: bool,
    /// Drop single-exon reps that engulf >= `READTHROUGH_MIN_DISTINCT` distinct junctions: unspliced
    /// pre-mRNA, never a copy. Default ON — see `is_unspliced_readthrough` for the validation.
    pub filter_readthrough: bool,
    /// Split mis-chained reads at spurious giant introns before seeding (opt-in). Default off = byte-identical.
    pub mischain_salvage: bool,
    /// Gate each co-located family by MUTUAL HOMOLOGY (`refine_families_exon_sum`, on the SHARED E_r
    /// primary tier since X.4 -- `-k 11 -w 5` @ 0.60 by default,
    /// cov-of-shorter>=0.50, + sensitive tier) across >= 2 distinct loci — the SAME criterion
    /// `gw_family_catalog` refines by. Default ON: without it the read-conflict oracle admits large-gene
    /// intra-gene mis-chains (PBX1) and repeat-bridges as families (`bench/GW_CATALOG_FP_AUDIT.md`). Requires
    /// minimap2 on PATH. Off ⇒ the raw conflict/homology families are assigned as-is.
    pub refine: bool,
    /// Admit a COLLAPSED single-rep locus as a multi-copy family: `n_copies = χ(H)`, reads certified Tied, no
    /// per-copy consensus materialised. **Default OFF**: the ambiguity instrument detects unresolvable
    /// PARALOGY, not collapse — it fires on EEF1A1 (pseudogenes on other chromosomes) with χ(H) = 7. See the
    /// `collapse_gate` module header.
    pub collapse_gate: bool,
    /// Re-admit near-identical families that collapse to < 2 RNA loci as K0_COLLAPSED, ENUMERATING the
    /// genome-projected copy count instead of abstaining. **Default OFF**: when off, all existing output is
    /// byte-identical. Env `RUSTLE_COLLAPSE_ENUMERATE=1`, CLI `--collapse-enumerate` on `gw_family_catalog`
    /// and `copy_assign`.
    pub collapse_enumerate: bool,
    /// Re-admit EXON-IDENTICAL (0-PSV) but heavily-EXPRESSED families that collapse to < 2 RNA loci as a
    /// K0_COLLAPSED_EXPRESSED copy-number class: no hidden-copy witness required, just >= 2 genome-projected
    /// loci that are EACH read-supported. **Default OFF**: when off, all existing output is byte-identical.
    /// Env `RUSTLE_COLLAPSE_EXPRESSED=1`, CLI `--collapse-expressed` on `gw_family_catalog`.
    pub collapse_expressed: bool,
    /// DNA-family fallback (DNA edge oracle): RNA-orphan loci (0 homology edges -> no RNA family) projected
    /// onto the genome at `dna_family_min_identity` to recover DIVERGENT genomic paralogs whose expressed
    /// transcript is non-homologous (DNA-family != RNA-family). Copy-number only. `RUSTLE_DNA_FAMILY_FALLBACK=1`,
    /// CLI `--dna-family-fallback`. Default off (byte-identical).
    pub dna_family_fallback: bool,
    /// Tied-secondary seeding (opt-in): after the primary-read pass1 skeletons are built, seed ADDITIONAL
    /// skeletons from AS-tied secondary reads (`rescue_extra`) that agree on an intron chain at loci with
    /// no primary skeleton (`tied_seed_skeletons`) — recovers starved co-located copies (K=0 members with 0
    /// primaries) as DETECTED-but-unassignable loci. Default OFF: when off, `skeletons` is unchanged and
    /// every downstream stage is byte-identical.
    pub tied_seed: bool,
    /// Projection identity floor for the DNA-family fallback (looser than collapse-expressed's 0.98, since
    /// these paralogs are divergent). CLI `--dna-family-min-identity`. Default 0.90.
    pub dna_family_min_identity: f64,
    /// Repeat gate for the DNA-family fallback: reject a re-admitted orphan locus whose reference region is
    /// >= this fraction soft-masked (RepeatMasker lowercase) — an audit showed the un-gated fallback re-admits
    /// Alus/low-complexity as spurious SD copies. `RUSTLE_DNA_FAMILY_MAX_SOFTMASK`. Default 0.30 (audit: on-member re-admissions are themselves ~70% repeat; 0.30 keeps the lone genuine recovery, cuts volume 715->81).
    pub dna_family_max_softmask: f64,
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
            vg_realign_admit: false,
            homology_primary: false,
            filter_readthrough: true,
            mischain_salvage: false,
            refine: true,
            collapse_gate: false,
            collapse_enumerate: false,
            collapse_expressed: false,
            dna_family_fallback: false,
            tied_seed: false,
            dna_family_min_identity: 0.90,
            dna_family_max_softmask: 0.30,
            eps_amb: Some(crate::vg_family::collapse_gate::GENOME_WIDE_EPS_AMB),
        }
    }
}

impl DenovoConfig {
    /// Read-count floor for "this locus is really there", shared by every path.
    ///
    /// It is STORED on `ConflictParams` because the conflict graph was the first consumer, but it is a plain
    /// scalar (default 3, `RUSTLE_CONFLICT_MIN_READS`) with no conflict-graph semantics. Reaching through
    /// `cfg.conflict.min_reads` inside `detect_homology_catalog_genome_wide` made the `--homology-primary`
    /// path *read* as if it consulted E_c, which is precisely the claim that path exists to avoid: there,
    /// family membership is decided by sequence homology alone (`homology_blocks`), and no conflict edge is
    /// ever built. Call this accessor on the homology path so the code says what is true.
    /// Enforced by `homology_catalog_never_touches_the_conflict_graph`.
    pub fn locus_min_reads(&self) -> usize {
        self.conflict.min_reads
    }

    /// Read overrides from `RUSTLE_*` env vars on top of `Default` (currently `RUSTLE_COLLAPSE_ENUMERATE`
    /// and `RUSTLE_COLLAPSE_EXPRESSED`).
    pub fn from_env() -> Self {
        DenovoConfig {
            collapse_enumerate: std::env::var("RUSTLE_COLLAPSE_ENUMERATE").ok().as_deref() == Some("1"),
            collapse_expressed: std::env::var("RUSTLE_COLLAPSE_EXPRESSED").ok().as_deref() == Some("1"),
            dna_family_fallback: std::env::var("RUSTLE_DNA_FAMILY_FALLBACK").ok().as_deref() == Some("1"),
            dna_family_max_softmask: env_num("RUSTLE_DNA_FAMILY_MAX_SOFTMASK", 0.30),
            // Readthrough/mis-chain gate SENSITIVITY toggle: RUSTLE_KEEP_READTHROUGH=1 disables both gates
            // (readthroughs that connect copies are then NOT filtered). Default (unset) = gates ON = today.
            filter_readthrough: std::env::var("RUSTLE_KEEP_READTHROUGH").ok().as_deref() != Some("1"),
            mischain_salvage: std::env::var("RUSTLE_MISCHAIN_SALVAGE").ok().as_deref() == Some("1"),
            ..Self::default()
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
                let au = crate::vg_family::family_graph::upper_cow(&a.seq);
                let bu = crate::vg_family::family_graph::upper_cow(&b.seq);
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
    /// The check is structural and needs no annotation. See `catalog_overlaps` and `bench/artifact_audit.py`.
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
    /// genomic (donor,acceptor) intron chain per copy, parallel to copy_tids. For v2 exon graph.
    pub copy_introns: Vec<Vec<(u64, u64)>>,
    /// per-copy genome remap identity, parallel to copy_tids: Some for reference-ABSENT copies
    /// (their discovery remap identity), None for in-genome copies (never remapped). For MI:f:.
    pub copy_map_identity: Vec<Option<f64>>,
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
            copy_introns: Vec::new(),
            copy_map_identity: Vec::new(),
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


/// Read-depth floor defining the core (`RUSTLE_ER_CORE_DEPTH`, default 10). Measured edge-level on
/// chr1+chr15: depth >= 10 separated best; depth >= 3 barely trims the span (core/span median 0.93) and a
/// relative floor (10% of the locus max) was slightly worse.
fn core_depth_floor() -> u32 {
    std::env::var("RUSTLE_ER_CORE_DEPTH").ok().and_then(|v| v.parse().ok()).unwrap_or(10)
}

/// Coverage floor applied against the read-supported CORE instead of the called span
/// (`RUSTLE_ER_CORE_COVERAGE`, unset = off = byte-identical). The core is a tighter target than the span,
/// so the equivalent demand is a HIGHER floor: 0.80 on the core corresponds to 0.50 on the span.
fn core_cov_floor() -> Option<f64> {
    let v = std::env::var("RUSTLE_ER_CORE_COVERAGE").ok()?;
    if v.is_empty() || v == "0" { return None; }
    v.parse::<f64>().ok().filter(|x| *x > 0.0)
}


/// SECOND coverage floor, charged on the LONGER sequence (`RUSTLE_ER_COVERAGE_LONGER_FLOOR`,
/// unset = off = byte-identical). This ADDS a clause; the shipped shorter-side floor is untouched, so
/// the rule becomes BLAST's `qcovs AND scovs`: the shorter must clear `min_coverage`, the longer must
/// clear this.
///
/// ⚠ NOT the same as `RUSTLE_ER_COVERAGE_LONGER`, which REPLACES the denominator and therefore applies
/// the SAME floor to the longer side — i.e. the SYMMETRIC variant. Offline arms (ledger §6bk/§6bl)
/// reject symmetry: at 0.50 on the longer side NPIP recall falls 14/31 -> 12/31 on `arm_f2` and the
/// genome-wide catalog loses 44% of its families, while an ASYMMETRIC 0.30 gains non-ZNF precision on
/// BOTH substrates at unchanged NPIP. The tree already predicted why: `denovo_pipeline.rs:4900-4926`
/// records that only 134/171 NPIP true pairs can reach 0.50 on the longer axis at all (NPIPB8-NPIPB2
/// caps at 0.215), because the duplicated unit is size-invariant while annotated spans are not.
///
/// Applied at the `E_r` site only. The tier-2 admission path (`RUSTLE_TIER2_ADMIT`, default off) keeps
/// its own one-sided test.
fn er_cov_longer_floor() -> Option<f64> {
    let v = std::env::var("RUSTLE_ER_COVERAGE_LONGER_FLOOR").ok()?;
    if v.is_empty() || v == "0" { return None; }
    v.parse::<f64>().ok().filter(|x| *x > 0.0)
}

/// Maximum `de` (gap-compressed divergence) for a read to define a locus BOUNDARY
/// (`RUSTLE_LOCUS_DE_EXTENT`, unset = off = byte-identical). Measured optimum 0.0005.
fn locus_de_extent() -> Option<f32> {
    let v = std::env::var("RUSTLE_LOCUS_DE_EXTENT").ok()?;
    if v.is_empty() || v == "0" { return None; }
    v.parse::<f32>().ok().filter(|x| *x > 0.0)
}

/// Locus extent defined by only the reads that CANNOT have come from a paralog.
///
/// The span of an emitted copy is inherited from whichever transcript became its representative, and that
/// span is what the E_r coverage test divides by and what every size comparison scores. Measured against
/// RefSeq gene spans (per copy, chr1+chr15): the shipped representative is in-band 0.5-2x only 25% of the
/// time with a median of 0.18 -- overwhelmingly TRUNCATED. Taking the extent of ALL primary reads instead
/// gives 36% and median 3.53 -- now overwhelmingly OVER-extended, because reads that bled in from a
/// paralog drag the boundary outward.
///
/// Filtering those out by divergence fixes both ends at once, and the effect is monotone until it turns:
///
///   primary reads, de <= 0.005    in-band 42%   median 2.05
///   primary reads, de <= 0.002    in-band 49%   median 1.29
///   primary reads, de <= 0.0005   in-band 58%   median 1.00   <- optimum
///   primary reads, de <= 0.0      in-band 49%   median 0.91   (turns over; too few reads survive)
///
/// de <= 0.0005 is roughly one mismatch per 2 kb: a read that can only have come from THIS copy. Absolute
/// thresholds beat per-locus relative trimming (58% vs 46% for dropping the most divergent half), which is
/// worth recording because relative is the more obvious design.
///
/// Returns `None` for a locus with no qualifying read (56 of 340 in the measurement) so the caller keeps
/// the representative's own span -- the change can only sharpen a boundary, never delete a locus.
///
/// The result is CLAMPED to contain the representative's own intron chain: the flanks are where
/// over-extension lives, but a junction is asserted evidence and shrinking past one would invalidate the
/// model rather than improve it.
pub(super) fn locus_confident_extent(
    bam_reads: &[BamRead],
    reps: &[DenovoTranscript],
    max_de: f32,
) -> Vec<Option<(u64, u64)>> {
    use std::collections::BTreeMap;
    let mut by_chrom: BTreeMap<&str, Vec<(u64, u64)>> = BTreeMap::new();
    for br in bam_reads {
        if br.is_supplementary || br.de > max_de {
            continue;
        }
        by_chrom
            .entry(br.chrom.as_str())
            .or_default()
            .push((br.read.ref_start, read_ref_end(&br.read)));
    }
    for v in by_chrom.values_mut() {
        v.sort_unstable();
    }
    reps.iter()
        .map(|r| {
            let v = by_chrom.get(r.chrom.as_str())?;
            let (mut lo, mut hi) = (u64::MAX, 0u64);
            for &(s, e) in v.iter() {
                if s >= r.end {
                    break;
                }
                if e > r.start {
                    lo = lo.min(s);
                    hi = hi.max(e);
                }
            }
            if lo == u64::MAX || hi <= lo {
                return None;
            }
            // never cut into the transcript's own asserted structure
            if let (Some(first), Some(last)) = (r.introns.first(), r.introns.last()) {
                lo = lo.min(first.0);
                hi = hi.max(last.1);
            }
            Some((lo, hi))
        })
        .collect()
}



/// Deny a STUB the right to CREATE family membership (`RUSTLE_ER_NO_STUB_EDGES=1`, default off).
///
/// The locus is kept; only its E_r edges are refused. Deleting stubs outright costs far more than it buys
/// (F1 0.704 -> 0.608) because the >=2-loci gate then dissolves whole families, and the ones that dissolve
/// are disproportionately small and pure. Refusing the EDGE leaves every locus in place, so a family whose
/// other members connect to each other survives intact.
///
/// Motivated by NPIP, where every false positive was a single-exon stub joining through a short perfect
/// match while the real members carried 4, 8 and 21 exons -- and where identity and coverage CANNOT reject
/// them, because the contaminant pairs score higher than the true pairs on both.
///
/// A genuinely intronless copy (`stub == false` despite one exon) is unaffected: it has no spliced evidence
/// at its locus, so it is a real intronless gene rather than a fragment of a spliced one.
fn er_no_stub_edges() -> bool {
    std::env::var("RUSTLE_ER_NO_STUB_EDGES").map(|v| v != "0" && !v.is_empty()).unwrap_or(false)
}

/// Whether reads at each rep's locus assert a splice junction carried by `>= min_reads` reads.
///
/// Combined with `n_exon == 1` this separates a STUB (single-exon representative of a spliced gene) from a
/// genuinely INTRONLESS copy. Both look identical in the copy table; only the reads distinguish them.
pub(super) fn locus_has_spliced_evidence(
    bam_reads: &[BamRead],
    reps: &[DenovoTranscript],
    min_reads: u32,
) -> Vec<bool> {
    use std::collections::BTreeMap;
    let mut by_chrom: BTreeMap<&str, BTreeMap<(u64, u64), u32>> = BTreeMap::new();
    for br in bam_reads {
        if br.is_supplementary {
            continue;
        }
        let mut p = br.read.ref_start;
        let e = by_chrom.entry(br.chrom.as_str()).or_default();
        for &(op, n) in &br.read.cigar {
            match op {
                'N' => {
                    *e.entry((p, p + n)).or_insert(0) += 1;
                    p += n;
                }
                'M' | 'D' | 'X' | '=' => p += n,
                _ => {}
            }
        }
    }
    reps.iter()
        .map(|r| {
            by_chrom.get(r.chrom.as_str()).is_some_and(|m| {
                m.range((r.start, 0)..(r.end, u64::MAX)).any(|(_, &c)| c >= min_reads)
            })
        })
        .collect()
}

/// Bases of each rep's span covered by at least `min_depth` reads: the READ-SUPPORTED CORE.
///
/// The E_r coverage test divides by the shorter copy's SPAN, so a boundary error moves the very denominator
/// that decides membership. A readthrough tail is low-depth and therefore never enters the core, so the core
/// does not inflate when the model runs past the locus. See `DenovoTranscript::core_bp` for the measurements.
///
/// Depth counts every primary read overlapping the span, matching what `samtools depth` reports there, so
/// the value is independent of how reads were attributed to loci.
pub(super) fn locus_core_bp(bam_reads: &[BamRead], reps: &[DenovoTranscript], min_depth: u32) -> Vec<u64> {
    use std::collections::BTreeMap;
    let mut by_chrom: BTreeMap<&str, Vec<(u64, u64)>> = BTreeMap::new();
    for br in bam_reads {
        if br.is_supplementary {
            continue;
        }
        by_chrom
            .entry(br.chrom.as_str())
            .or_default()
            .push((br.read.ref_start, read_ref_end(&br.read)));
    }
    for v in by_chrom.values_mut() {
        v.sort_unstable();
    }
    reps.iter()
        .map(|r| {
            let Some(v) = by_chrom.get(r.chrom.as_str()) else { return 0 };
            // Sweep only the read endpoints that fall in this rep's span; a position's depth changes only
            // there, so counting covered bases between consecutive events is exact.
            let mut events: Vec<(u64, i32)> = Vec::new();
            for &(s, e) in v.iter() {
                if s >= r.end {
                    break;
                }
                if e > r.start {
                    events.push((s.max(r.start), 1));
                    events.push((e.min(r.end), -1));
                }
            }
            if events.is_empty() {
                return 0;
            }
            events.sort_unstable();
            let (mut depth, mut prev, mut core) = (0i32, events[0].0, 0u64);
            for (pos, delta) in events {
                if depth >= min_depth as i32 && pos > prev {
                    core += pos - prev;
                }
                depth += delta;
                prev = pos;
            }
            core
        })
        .collect()
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

use crate::vg_family::collapse_gate::{collapse_verdict, Ambiguity, CollapseVerdict};

/// Ambiguously-placed primary reads over a rep's span.
///
/// A read is AMBIGUOUS iff `mapq == 0`: the aligner found no reason to prefer this placement over another.
/// Supplementary records are excluded for the same reason they are excluded from conflict edges — a chimeric
/// segment is adjacency, not ambiguity. ⚠ Secondary records DO enter `bam_reads` as separate entries
/// (`reads_in_region` collects every mapped record) — this counter deliberately includes them, because a
/// multimapper's alternative placements are exactly the ambiguity being measured.
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
        copy_introns: vec![rep.introns.clone()],
        copy_map_identity: vec![None], // rep is an in-genome assembled transcript, never remapped
        ..FamilyAssignment::empty()
    }
}

/// `--tied-seed` existence-only append (mirrors the collapse-gate `gated_family` block). Each tied rep whose
/// span overlaps NO already-emitted copy span becomes a singleton `TSFAM` family that TIES every overlapping
/// read: the copy is DETECTED (existence, O1) but every read is `AssignStatus::Tied` — assignment (O2) is
/// abstained, because a starved copy's reads are K=0 with their tie partner. Tied reps are computed OUTSIDE
/// the primary `reps`/conflict/refine/assignment pipeline, so appending them here leaves the primary copy set
/// byte-identical to a no-`--tied-seed` run (the fix for the amylase 21->6 over-merge and the os1 rep-shift).
/// `emitted_spans` = every copy span already in `out`. Deterministic (input order); tid ids run
/// `TSFAM{next_family_index + k}`. A tied rep overlapping an earlier-admitted tied rep is also skipped.
fn tied_existence_families(
    tied_reps: &[DenovoTranscript],
    emitted_spans: &[(String, u64, u64)],
    bam_reads: &[BamRead],
    next_family_index: usize,
) -> Vec<FamilyAssignment> {
    let mut out = Vec::new();
    let mut occupied: Vec<(String, u64, u64)> = emitted_spans.to_vec();
    for t in tied_reps {
        let overlaps = occupied
            .iter()
            .any(|(c, s, e)| c == &t.chrom && *s < t.end && t.start < *e);
        if overlaps {
            continue;
        }
        let fid = format!("TSFAM{}", next_family_index + out.len());
        out.push(gated_family(t, bam_reads, 1, fid));
        occupied.push((t.chrom.clone(), t.start, t.end));
    }
    out
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
/// Env override for a gate threshold (sensitivity-analysis knob). Unset => the compiled default, so a normal
/// run is byte-identical. Lets `bench/soto/gate_sensitivity_sweep.sh` vary the readthrough/mis-chain gates
/// without recompiling.
fn env_num<T: std::str::FromStr>(key: &str, default: T) -> T {
    std::env::var(key).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
}

pub fn retain_non_readthrough(
    transcripts: &mut Vec<DenovoTranscript>,
    support: &std::collections::HashMap<(String, u64, u64), usize>,
    tag: &str,
) {
    let min_distinct = env_num("RUSTLE_READTHROUGH_MIN_DISTINCT", READTHROUGH_MIN_DISTINCT);
    let mut dropped = Vec::new();
    transcripts.retain(|t| {
        let rt = is_unspliced_readthrough(t, support, READTHROUGH_MIN_SUPPORT, min_distinct);
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
            min_distinct,
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
    let giant_bp = env_num("RUSTLE_MISCHAIN_GIANT_BP", MISCHAIN_GIANT_INTRON_BP);
    let min_reads = env_num("RUSTLE_MISCHAIN_MIN_READS", MISCHAIN_MIN_JUNCTION_READS);
    let mut dropped = Vec::new();
    transcripts.retain(|t| {
        let mc = is_giant_intron_mischain(t, support, giant_bp, min_reads);
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
            giant_bp,
            min_reads,
            transcripts.len()
        );
        for d in &dropped {
            eprintln!("[{tag}]   mis-chain {d}");
        }
    }
}

/// Is this READ mis-chained — does its alignment leave the locus through a giant gap?
///
/// True iff any of the read's OWN introns exceeds `giant_bp`. **There is deliberately no support clause**,
/// and that is the entire difference from [`is_giant_intron_mischain`]:
///
/// * [`is_giant_intron_mischain`] asks *is this TRANSCRIPT MODEL spurious?* A giant intron carried by many
///   reads may be a real deep large gene, so it is kept — the documented scope limit.
/// * This asks *do THIS READ's bases belong at THIS locus?* A read chained across a giant gap lays its
///   far-end bases into a pileup that is hundreds of kb away, and it does that whether one read or a
///   hundred make the same chaining mistake. **A popular mis-chain is evidence of systematic mis-chaining,
///   not of truth.** So for any per-read statistic scoped to a locus — identity, divergence, clipping, a
///   PSV column — the support clause is not merely unhelpful, it is the wrong question.
///
/// Measured on the matched-individual substrate (`o3_mischain/fix/red/`, denominators printed there):
/// at the locus `GWFAM244:2` **150 of 150** primary reads carry an intron > `MISCHAIN_GIANT_INTRON_BP`,
/// the dominant mis-chain is **827,011 bp** and is carried by **97** primary reads — **32×**
/// [`MISCHAIN_MIN_JUNCTION_READS`] — so [`is_giant_intron_mischain`] returns `false` for every dominant
/// mis-chain at that locus while this predicate returns `true` for every one of its reads.
pub fn is_mischained_read(r: &PrimaryRead, giant_bp: u64) -> bool {
    r.introns.iter().any(|&(d, a)| a.saturating_sub(d) > giant_bp)
}

/// Drop mis-chained reads ([`is_mischained_read`]) before any locus-scoped per-read statistic. Returns the
/// number dropped so the caller can print the denominator. Threshold from the single source of truth
/// [`MISCHAIN_GIANT_INTRON_BP`], env-overridable exactly like the transcript-level filter, so a sweep moves
/// both rules together and neither can silently disagree about what "giant" means.
///
/// This is NOT [`split_mischained_reads`]: that one SALVAGES (cuts a read into its local pieces and keeps
/// both) and is the right tool when the reads are being assembled. This one REMOVES, which is the right
/// tool when the reads are being measured — a salvaged half-read is still a read whose mate half sits at
/// another copy, and a locus statistic must not average over the two.
pub fn retain_local_reads(reads: &mut Vec<PrimaryRead>) -> usize {
    let giant = env_num("RUSTLE_MISCHAIN_GIANT_BP", MISCHAIN_GIANT_INTRON_BP);
    let before = reads.len();
    reads.retain(|r| !is_mischained_read(r, giant));
    before - reads.len()
}

/// Mis-chain salvage (opt-in). `Some(split reads)` when `cfg.mischain_salvage`, else `None` (caller keeps the
/// originals — byte-identical). Reuses the gate thresholds so it splits exactly the introns the gate would drop.
fn maybe_salvage_mischain(reads: &[PrimaryRead], cfg: &DenovoConfig) -> Option<Vec<PrimaryRead>> {
    if !cfg.mischain_salvage {
        return None;
    }
    let giant = env_num("RUSTLE_MISCHAIN_GIANT_BP", MISCHAIN_GIANT_INTRON_BP);
    let min_reads = env_num("RUSTLE_MISCHAIN_MIN_READS", MISCHAIN_MIN_JUNCTION_READS);
    let support = read_junction_support(reads);
    Some(split_mischained_reads(reads, &support, giant, min_reads))
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
    // Built once from the stable pre-correction `copy_profiles` (see comment above), same as the
    // family's own one-shot pass -- every corrected read in this loop threads the same graph.
    let graph = BubbleGraph::from_copies(&copy_profiles);

    for (read_idx, (new_copy, obs)) in apply.corrected {
        if let Some(pos) = fa.assignments.iter().position(|(ri, _)| *ri == read_idx) {
            let rf = ReadFeatures {
                psv_obs: obs.clone(),
                psv_qual: Vec::new(),
                junctions: fa.read_junctions[pos].clone(),
            };
            fa.assignments[pos].1 = match assign_read_editing(&rf, &graph, &copy_profiles, p, &editing_cols) {
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
        // `from_env` == `default` unless RUSTLE_ABSENT_MIN_CLUSTERS is set (removal-ablation only).
        absent_copy::admit_candidate(cand, host, genome, fasta_path, &AbsentCopyParams::from_env())
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
        let (t, adm_id) = match admitted {
            Admission::Copy(t, id) => (t, id),
            Admission::DnaNeeds(_) => continue, // gate declined -- stays a "novel-candidate" record.
        };

        let t_g2o = gen2off(&t);
        let positions = psv_positions_for(&host_col_gpos, &t_g2o, t.seq.len());
        let alleles: Vec<Option<u8>> = positions.iter().map(|&off| t.seq.get(off).copied()).collect();
        let new_copy = fa.copy_psv_alleles.len();
        fa.copy_psv_alleles.push(alleles);
        fa.copy_junctions.push(copy_boundaries(&t));
        fa.copy_tids.push(t.tid.clone());
        fa.copy_spans.push((t.chrom.clone(), t.start, t.end));
        fa.copy_introns.push(t.introns.clone());
        fa.copy_map_identity.push(adm_id);
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

/// Rebuild a `ColocatedFamily` from a refined copy list (`refine_families_exon_sum` returns bare copy sets). The
/// span is the min/max over the copies; the id is a placeholder — the emitting binary renumbers families on output.
fn colocated_from_copies(idx: usize, mut copies: Vec<DenovoTranscript>) -> ColocatedFamily {
    // preserve the by-start copy ordering `colocated_families` guarantees (the assignment step indexes copies
    // by this order); `refine_families_exon_sum`/`distinct_locus_reps` may reorder them.
    copies.sort_by_key(|c| c.start);
    let chrom = copies.first().map(|c| c.chrom.clone()).unwrap_or_default();
    let start = copies.iter().map(|c| c.start).min().unwrap_or(0);
    let end = copies.iter().map(|c| c.end).max().unwrap_or(0);
    ColocatedFamily { family_id: format!("CAFAM{idx}"), chrom, start, end, copies }
}

/// FNV-1a of `t_seq`, reduced to a `linearize_certificate` seed. Deterministic (same candidate ->
/// same decoy shuffles -> a reproducible certificate), and distinct candidates get distinct seeds
/// without threading an RNG through `detect_and_assign`.
fn fnv1a_seed(t_seq: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in t_seq {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Package the stage-2 admission pool build + `linearize::linearize_certificate` call for one admitted
/// absent copy: `pool` is the MAPQ-0 read pool at this family's region, `copy_seqs` are the OTHER
/// already-known copies (the candidate itself is appended internally as the last contig — see
/// `linearize_certificate`). Fixed defaults (`n_decoys=20, min_pool=5, alpha=0.05`); `seed` is a
/// deterministic FNV-1a hash of `t_seq` so the certificate is reproducible per candidate.
///
/// `n_decoys=20` (not 19): a permutation test at alpha=0.05 needs at least 20 decoys, so the perm_p
/// floor `1/(n_decoys+1) = 1/21 ~= 0.048` is strictly below alpha and a perfect candidate can reach
/// `Linearizes`. (Decoys are now the 20 dinucleotide shuffles only — the reverse-complement decoy was
/// removed because minimap2 is strand-symmetric and RC(candidate) would tie `real`; see
/// `linearize_certificate`. Previously 19 shuffles + 1 RC also summed to 20 decoys, so this preserves
/// the same statistical resolution with all decoys valid.)
fn family_linearize_cert(
    t_seq: &[u8],
    copy_seqs: &[Vec<u8>],
    pool: &[Vec<u8>],
    realign: impl Fn(&[Vec<u8>], &[Vec<u8>]) -> Vec<Option<(usize, u32)>>,
) -> super::linearize::LinearizeCertificate {
    super::linearize::linearize_certificate(t_seq, copy_seqs, pool, 20, fnv1a_seed(t_seq), 5, 0.05, realign)
}

/// The augment-and-linearize OPT-IN guard: compute the certificate for one admitted candidate ONLY when
/// `do_linearize` (= --linearize or --linearize-gate) is set; otherwise return `None` and — crucially —
/// never call `realign`, so plain `--absent-copies` incurs no minimap2 realign-pool subprocess per
/// candidate. Isolating the guard here makes the opt-in deterministically testable (inject a fake
/// `realign`) without a full end-to-end absent-copy admission.
fn linearize_cert_if_enabled(
    do_linearize: bool,
    t_seq: &[u8],
    copy_seqs: &[Vec<u8>],
    pool: &[Vec<u8>],
    realign: impl Fn(&[Vec<u8>], &[Vec<u8>]) -> Vec<Option<(usize, u32)>>,
) -> Option<super::linearize::LinearizeCertificate> {
    if do_linearize {
        Some(family_linearize_cert(t_seq, copy_seqs, pool, realign))
    } else {
        None
    }
}

/// Detect co-located families in a region and assign every read to a copy.
///
/// # `supplied_families` — the O1→O2 file contract (`copy_assign --families`)
///
/// `None` (the default, and every historical caller) = the shipped behaviour: families are DERIVED here
/// from the BAM (pass1 → gate → collapse → membership oracle → co-location → refine).
///
/// `Some(f)` = the copy set is GIVEN (a `gw_family_catalog` `copies.tsv`, materialized by
/// `catalog_input::to_colocated`), and this function DOES NOT CONSTRUCT FAMILIES AT ALL. Specifically it
/// skips, and must keep skipping:
///
/// * `pass1_skeletons_robust` and everything reached through it — with no skeletons, `transcripts`/`reps`
///   are empty, so the placement pass, the membership oracle (E_c conflict graph or E_r homology), the POA
///   diagnostic, `colocated_families` and the collapse gate all degenerate to no-ops for free;
/// * the REFINE gate (`refine_families_exon_sum`), which re-clusters a copy set by its own edge tiers —
///   re-deriving membership on top of a supplied roster is exactly what this mode exists to avoid;
/// * the family RESCUE leg (`rescue_thin_loci_iterative`), which ADDS under-assembled copies to the family.
///
/// What still runs is the ASSIGNMENT of BAM reads to the given copies (PSV + junction likelihood +
/// significance gate) — reads always come from the BAM; only the copy set is supplied. The other
/// roster-changing legs (`absent_copies`, `vg_realign` admission, `iterative_prune`, `collapse_gate`,
/// `tied_seed`) are refused by the CLI when `--families` is given rather than silently widened here.
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
    do_linearize: bool,
    linearize_gate: bool,
    fasta_path: &str,
    supplied_families: Option<&[ColocatedFamily]>,
) -> (
    Vec<FamilyAssignment>,
    Vec<FallbackEdge>,
    Vec<DnaNeedsRecord>,
    Vec<(String, super::linearize::LinearizeCertificate, (String, u64, u64))>,
) {
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
    // O1→O2 FILE CONTRACT: with a supplied catalog the whole detection front end is DEAD WORK and, worse,
    // would be a second membership oracle running behind the user's back. Emptying the skeletons is enough
    // to switch it all off at once: `assemble_gate([])` = no transcripts = no reps, and every stage below
    // (`build_read_placements`, the conflict/homology oracle, the POA diagnostic, `colocated_families`, the
    // collapse gate) is a fold over `reps` that is vacuous on an empty rep set. Stated once, here, rather
    // than as five separate `if supplied` guards that could drift apart.
    let supplied = supplied_families.is_some();
    // Same `k` as the O1 catalogs: one canonical extent per locus across objectives (see `detect_families`).
    let salvaged = if supplied { None } else { maybe_salvage_mischain(primary_reads, cfg) };
    let seed_reads: &[PrimaryRead] = salvaged.as_deref().unwrap_or(primary_reads);
    let skeletons = if supplied {
        Vec::new()
    } else {
        pass1_skeletons_robust(seed_reads, cfg.pass1_min_reads, cfg.min_terminal_support)
    };
    // TIED-SEED (opt-in): assemble the tied-seed skeletons into their OWN reps, kept ENTIRELY OUT of the
    // primary `reps` / conflict / refine / assignment pipeline. K=0 tied reps mixed into `reps` add spurious
    // conflict edges that over-merge families and wreck assignment (chr1 amylase 21->6 copies; os1 rep-shift).
    // These are appended as existence-only `TSFAM` families at the END (see `tied_existence_families`), so the
    // primary copy set stays byte-identical to a no-`--tied-seed` run.
    let tied_reps: Vec<DenovoTranscript> = if cfg.tied_seed {
        let tied_sk = tied_seed_skeletons(rescue_extra, &skeletons, cfg.pass1_min_reads);
        // DIAGNOSTIC (env-gated, no behavior change when unset): dump each tied-seed skeleton with the real
        // primary-read coverage at its span. RUSTLE_TIED_SEED_DEBUG=<path>.
        if let Ok(dbg_path) = std::env::var("RUSTLE_TIED_SEED_DEBUG") {
            use std::io::Write;
            if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(&dbg_path) {
                for t in &tied_sk {
                    let ov = |p: &PrimaryRead| p.chrom == t.chrom && p.ref_start < t.end && t.start < p.ref_end;
                    let np = primary_reads.iter().filter(|p| ov(p)).count();
                    let np_uns = primary_reads.iter().filter(|p| ov(p) && p.introns.is_empty()).count();
                    let np_spl = primary_reads.iter().filter(|p| ov(p) && !p.introns.is_empty()).count();
                    let n_tied_ov = rescue_extra.iter().filter(|p| ov(p)).count();
                    let n_sk = skeletons.iter().filter(|s| s.chrom == t.chrom && s.start < t.end && t.start < s.end).count();
                    let _ = writeln!(
                        f, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        t.chrom, t.start, t.end, t.n_reads, t.introns.is_empty(),
                        np, np_uns, np_spl, n_tied_ov, n_sk
                    );
                }
            }
        }
        let mut tied_tx = assemble_gate(&tied_sk, genome, &cfg.gate);
        if cfg.filter_readthrough {
            let support = read_junction_support(primary_reads);
            retain_non_readthrough(&mut tied_tx, &support, "detect_and_assign(tied)");
            retain_non_mischain(&mut tied_tx, &support, "detect_and_assign(tied)");
        }
        let tied_rep_idx = collapse_loci_span_aware(&tied_tx, &cfg.detect);
        tied_rep_idx.iter().map(|&i| tied_tx[i].clone()).collect()
    } else {
        Vec::new()
    };
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
    let mut reps: Vec<DenovoTranscript> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
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
    // Per-copy unique-mapper support (Task 3, identifiability-merge): record each rep's count of MAPQ>0
    // (unambiguous) placements — free, `placements` is already built above — so `distinct_locus_reps`'s
    // same-strand merge guard (inside `refine_families_exon_sum` below) has real read evidence instead of
    // unconditionally collapsing every same-strand co-located pair.
    let uniq_counts = locus_unique_mapper_counts(&placements, reps.len());
    for (i, c) in reps.iter_mut().enumerate() {
        c.distinguishing_uniq = uniq_counts[i];
    }
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
    // Report-first augment-and-linearize certificate (Task 4): one entry per Stage-2-admitted reference-absent
    // copy, computed against that candidate's MAPQ-0 read pool. Admission itself is unchanged here (Task 5
    // adds the opt-in `--linearize-gate`); this is purely additive reporting threaded out to the caller.
    let mut linearize_certs: Vec<(String, super::linearize::LinearizeCertificate, (String, u64, u64))> =
        Vec::new();
    // Co-located families, then (default) the SAME mutual-homology + distinct-locus refinement gate that
    // `gw_family_catalog` applies — so the per-region path and the genome-wide catalog agree on what a family is.
    // Without it the conflict oracle admits large-gene mis-chains (PBX1) and repeat-bridges (see
    // `bench/GW_CATALOG_FP_AUDIT.md`). `cfg.refine` off ⇒ assign the raw families (no minimap2).
    //
    // ⚠ **RESIDUAL D1 (2026-08-09, KNOWN, NOT FIXED HERE).** `gw_family_catalog` made refine OPT-IN on
    // the homology (E_r) catalog, because refine appends a second clustering stage — its own core +
    // genomic-span edge SUBSTRATES, clustered by CONNECTED COMPONENTS — on top of γ-quasi-clique(E_r),
    // and none of that is in `docs/seeded_family_definition.md` §1. This site still runs refine
    // unconditionally on `cfg.refine`. That is CORRECT for copy_assign's DEFAULT (`homology_primary`
    // defaults to false ⇒ these are E_c conflict-derived co-located families, which is exactly the
    // catalog refine was written for and where it was measured to bite). It reproduces the D1 shape
    // ONLY under `copy_assign --homology-primary`. Left alone deliberately: that combination was not
    // measured, and an unmeasured behavioural change is worse than a recorded one.
    //
    // With `supplied_families` the roster is GIVEN, so neither `colocated_families` nor refine runs: this
    // mode exists precisely so O2 stops deriving its own answer to "what is a family". `reps` is empty here
    // (see the front-end note above), so `colocated_families` would return nothing anyway — the branch is
    // written explicitly so the skip is a stated decision rather than a side effect of an empty vector.
    let mut colocated: Vec<ColocatedFamily> = match supplied_families {
        Some(f) => {
            eprintln!(
                "[detect_and_assign] --families: {} supplied catalog famil{} assigned AS GIVEN \
                 (no detection, no refine, no rescue)",
                f.len(),
                if f.len() == 1 { "y" } else { "ies" }
            );
            f.to_vec()
        }
        None => colocated_families(&reps, &split, win, min_copies, &cfg.detect),
    };
    if cfg.refine && !supplied {
        let before = colocated.len();
        let copysets: Vec<Vec<DenovoTranscript>> = colocated.iter().map(|c| c.copies.clone()).collect();
        // RNA exon-sum substrate -> the forward-only E_r guard applies (DEFAULT ON 2026-08-19).
        let refine_params = RefineParams {
            intron_fasta: Some(fasta_path.to_string()),
            require_forward_alignment: true,
            substrate: Substrate::TranscriptOriented,
            ..Default::default()
        };
        let refined = refine_families_exon_sum(copysets, &refine_params, Some(genome), cfg.conflict.min_reads)
            .expect("refine (default): refine_families_exon_sum failed — is minimap2 on PATH? pass --no-refine to skip");
        // ⚠ This line used to hard-code "asm20 id>=0.80". It was already wrong for `--min-identity` runs
        // and became wrong for EVERY default run when the primary tier moved to the sensitive seed (X.4),
        // so it is derived from the same selector refine actually uses.
        let (log_seed, log_floor, _) = er_primary_tier(&refine_params);
        eprintln!(
            "[detect_and_assign] refine: {} co-located families -> {} homology-gated ({} id>={:.2}, cov>={:.2}, >= 2 distinct loci)",
            before,
            refined.len(),
            log_seed.join(" "),
            log_floor,
            refine_params.min_coverage
        );
        colocated = refined.into_iter().enumerate().map(|(i, c)| colocated_from_copies(i, c)).collect();
    }
    for cf in colocated {
        // RESCUE: recover under-assembled copies homologous to this family (below the >=3-read assembly gate)
        // and ADD them to the copy set, so reads can be assigned to them too. Iterative (bridge-aware).
        //
        // ⚠ SKIPPED under `--families`. Rescue WIDENS the roster — it is copy construction, just at a lower
        // read floor — so leaving it on would mean `copy_assign --families` assigned against a copy set that
        // is not the one it was handed, and the O1/O2 sets could not be compared. Requirement: with a
        // supplied catalog the copy set is EXACTLY the supplied catalog.
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
        let rescued = if supplied {
            Vec::new()
        } else {
            let region_primary: Vec<PrimaryRead> = primary_reads
                .iter()
                .chain(rescue_extra.iter())
                .filter(|r| r.chrom == cf.chrom && r.ref_start < rhi && r.ref_end > rlo)
                .cloned()
                .collect();
            let loci = thin_loci(&region_primary, RESCUE_MIN_SUPPORT);
            rescue_thin_loci_iterative(&loci, &members, &member_spans, genome, &RescueParams::default(), 3)
        };
        let rescued_copies = rescued.len();
        // tid -> discovery remap identity, for reference-ABSENT copies admitted below (Stage-2). Keyed on
        // tid (not index) so it survives the pruning/reassignment `all_copies` undergoes afterward; an
        // in-genome copy is simply absent from this map, which the main build below reads as `None`.
        let mut remap_id_by_tid: std::collections::HashMap<String, f64> = std::collections::HashMap::new();
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
             ..Default::default() });
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
        // O2.14: `region` holds alignment RECORDS; a multimapping molecule appears once per placement (the
        // shipped BAM is `-Y`, so its secondary records carry the full SEQ and each is independently
        // scorable). The names make the MOLECULE the unit of a result — see `assign_family_detailed_once`.
        let region_names: Vec<String> = idx_map.iter().map(|&i| bam_reads[i].name.clone()).collect();
        // Unique-mapper support is a property of the MOLECULE ("this molecule has a mapq>0 placement"),
        // not of whichever of its records represents it after the reduction — only the PRIMARY record
        // carries mapq>0, and the representative need not be the primary. Taking the max over a molecule's
        // records keeps `uniq`/`anchored` counting exactly the molecules they counted before the reduction.
        let region_mapq_mol: Vec<u8> = {
            let mut best: std::collections::HashMap<&str, u8> = std::collections::HashMap::new();
            for (n, &q) in region_names.iter().zip(region_mapq.iter()) {
                let e = best.entry(n.as_str()).or_insert(0);
                *e = (*e).max(q);
            }
            region_names.iter().map(|n| best[n.as_str()]).collect()
        };
        // Stage-1 and Stage-2 run WITHOUT iterative pruning so freeze_merge and downstream bookkeeping stay
        // in the original index space. Pruning, when requested, is applied as a final post-process below.
        let p_once = AssignParams { iterative_prune: false, ..*p };
        // Stage-1: assign over the reference copies only (borrow scoped so `all_copies` stays reassignable).
        let mut detail = {
            let copies: Vec<&DenovoTranscript> = all_copies.iter().collect();
            assign_family_detailed(&copies, &region, &p_once, Some(genome), Some(&region_names))
        };
        // Task 5 (opt-in): two-stage freeze for reference-ABSENT (collapsed) copies. OFF => this whole block
        // is skipped, so the loop below is byte-for-byte the pre-Task-5 path (`all_copies`/`detail` unchanged).
        if absent_copies {
            let cands = recover_collapsed_candidates(&all_copies, bam_reads);
            let mut admitted: Vec<DenovoTranscript> = Vec::new();
            // Augment-and-linearize inputs (opt-in via `do_linearize` = --linearize or --linearize-gate).
            // `pool` (this family's MAPQ-0 / ambiguous read pool) and `copy_seqs` (the reference copies) are
            // INVARIANT across the `for cand` loop below — `all_copies` is only reassigned AFTER the loop —
            // so build them once here rather than per candidate. Empty (no minimap2 work) when opt-out.
            let (pool, copy_seqs): (Vec<Vec<u8>>, Vec<Vec<u8>>) = if do_linearize {
                (
                    region.iter().zip(region_mapq.iter())
                        .filter(|(_, &q)| q == 0).map(|(r, _)| r.seq.clone()).collect(),
                    all_copies.iter().map(|c| c.seq.clone()).collect(),
                )
            } else {
                (Vec::new(), Vec::new())
            };
            for cand in &cands {
                if let Some(host) = all_copies.iter().find(|t| t.tid == cand.host_tid) {
                    // `from_env` == `default` unless RUSTLE_ABSENT_MIN_CLUSTERS is set (removal-ablation only).
                    match absent_copy::admit_candidate(cand, host, genome, fasta_path, &AbsentCopyParams::from_env()) {
                        Admission::Copy(t, id) => {
                            if let Some(v) = id {
                                remap_id_by_tid.insert(t.tid.clone(), v);
                            }
                            // Augment-and-linearize certificate (opt-in, report-first): re-align this family's
                            // MAPQ-0 (ambiguous) read pool against the reference copies + this candidate, and
                            // score how uniquely they land on the new copy vs a dinucleotide-shuffled decoy.
                            // SKIPPED entirely (no minimap2, no push) when `do_linearize` is false, so plain
                            // `--absent-copies` keeps its prior admission cost and emits no certificate.
                            if let Some(cert) = linearize_cert_if_enabled(
                                do_linearize, &t.seq, &copy_seqs, &pool, absent_copy::realign_pool_minimap2,
                            ) {
                                let verdict = cert.verdict;
                                linearize_certs.push((cf.family_id.clone(), cert, (t.chrom.clone(), t.start, t.end)));
                                // Opt-in gate (Task 5): a candidate that does NOT linearize is demoted to a
                                // DNA-needs record instead of admitted, when `--linearize-gate` is set. Only
                                // reachable when the cert was computed (do_linearize implied by the gate).
                                if linearize_gate && !matches!(verdict, super::linearize::Verdict::Linearizes) {
                                    // UNDETERMINED (perm_p is NaN, pool below min) is a distinct cause from a
                                    // computed non-linearizing perm_p >= alpha — report the true reason.
                                    let reason = if matches!(verdict, super::linearize::Verdict::Undetermined) {
                                        "pool too small for linearize certificate"
                                    } else {
                                        "did not linearize (perm_p >= alpha)"
                                    };
                                    dna_needs.push(DnaNeedsRecord {
                                        chrom: t.chrom.clone(),
                                        start: t.start,
                                        end: t.end,
                                        n_clusters: cand.n_clusters,
                                        reason: reason.to_string(),
                                        read_count: cand.iso.read_count,
                                    });
                                    continue;
                                }
                            }
                            admitted.push(t);
                        }
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
                    let mut d2 = assign_family_detailed(&copies2, &region, &p_once, Some(genome), Some(&region_names));
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
            detail = assign_family_detailed_pruned(&copies, &region, p, Some(genome), Some(&region_names));
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
            copy_introns: all_copies.iter().map(|c| c.introns.clone()).collect(),
            copy_map_identity: all_copies.iter().map(|c| remap_id_by_tid.get(&c.tid).copied()).collect(),
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
            if assigned_j && region_mapq_mol[r.read_index] > 0 {
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
            // Admission (roster-widening, genome-touching, O4-frontier FP risk) is gated SEPARATELY from the
            // correction leg: correction-only re-threads + reassigns among EXISTING copies, admits none.
            if cfg.vg_realign_admit {
                admit_novel_pools(&mut fa, &novel_pools, bam_reads, &all_copies, genome, fasta_path, &profiles);
            } else {
                let _ = (&novel_pools, &profiles);
            }
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

    // TIED-SEED existence-only append. `tied_reps` were assembled/collapsed separately at the top (never in
    // the primary `reps`), so `out` above is byte-identical to a no-`--tied-seed` run. Here each tied rep at a
    // locus no emitted copy occupies becomes a singleton `TSFAM` family that ties every overlapping read
    // (detected, unassignable). Additive-only: never merges or drops a primary copy.
    if cfg.tied_seed && !tied_reps.is_empty() {
        let emitted: Vec<(String, u64, u64)> =
            out.iter().flat_map(|fa| fa.copy_spans.iter().cloned()).collect();
        let tied_fams = tied_existence_families(&tied_reps, &emitted, bam_reads, out.len());
        eprintln!(
            "[detect_and_assign] tied-seed existence-only: {} tied reps -> {} appended (non-overlapping)",
            tied_reps.len(),
            tied_fams.len()
        );
        out.extend(tied_fams);
    }

    (out, fallback, dna_needs, linearize_certs)
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
    let reads = maybe_salvage_mischain(&reads, cfg).unwrap_or(reads);
    let skeletons = pass1_skeletons_robust(&reads, cfg.pass1_min_reads, cfg.min_terminal_support);
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
    let mut reps: Vec<DenovoTranscript> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
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
    // Parallelized ACROSS chromosomes: each chrom's reads_in_region + build_read_placements + conflict_edges
    // is independent and its edges are remapped to GLOBAL rep indices before merge. Results are re-assembled
    // in sorted-chrom order (rayon collect preserves the indexed order of `chrom_work`), so `all_edges` and
    // the per-chrom stderr log are BYTE-IDENTICAL to the previous serial build. A bounded pool of `threads`
    // workers caps concurrent chromosome read-loads (WSL2 memory), and each chrom decodes single-threaded
    // since the parallelism is now across chromosomes (avoids threads^2 BGZF oversubscription).
    use rayon::prelude::*;
    let chrom_work: Vec<(&str, &Vec<usize>)> =
        by_chrom.iter().filter(|(_, g)| g.len() >= 2).map(|(c, g)| (*c, g)).collect();
    // Per-copy unique-mapper support (Task 3, identifiability-merge), piggybacked on this SAME per-chrom
    // placement build: `(global_rep_idx, mapq>0_placement_count)` pairs, free — `placements` is already
    // built here for `conflict_edges`. Feeds `distinct_locus_reps`'s same-strand merge guard downstream
    // (via `refine_families_exon_sum`) so it sees real read evidence instead of unconditionally collapsing.
    type ChromEdges = (Vec<(usize, usize, usize)>, Vec<(usize, usize)>, (String, usize, usize, usize));
    let run = |&(chrom, glob): &(&str, &Vec<usize>)| -> Option<ChromEdges> {
        let clen = genome.chrom_len(chrom);
        let (_primary, bam_reads) = reads_in_region(bam_path, chrom, 0, clen, 1).ok()?;
        let chrom_reps: Vec<DenovoTranscript> = glob.iter().map(|&g| reps[g].clone()).collect();
        let placements = build_read_placements(&bam_reads, &chrom_reps);
        let edges = conflict_edges(chrom_reps.len(), &placements, &cfg.conflict);
        let uniq_local = locus_unique_mapper_counts(&placements, chrom_reps.len());
        let uniq_global: Vec<(usize, usize)> =
            uniq_local.into_iter().enumerate().map(|(li, cnt)| (glob[li], cnt)).collect();
        let n_edges = edges.len(); // this chrom's edge count (not the cumulative total)
        let remapped: Vec<(usize, usize, usize)> =
            edges.into_iter().map(|(i, j, w)| (glob[i], glob[j], w)).collect(); // local → global rep index
        Some((remapped, uniq_global, (chrom.to_string(), chrom_reps.len(), bam_reads.len(), n_edges)))
    };
    let per_chrom: Vec<ChromEdges> = match rayon::ThreadPoolBuilder::new().num_threads(threads.max(1)).build() {
        Ok(pool) => pool.install(|| chrom_work.par_iter().filter_map(&run).collect()),
        Err(_) => chrom_work.iter().filter_map(&run).collect(), // serial fallback (byte-identical)
    };
    let mut all_edges: Vec<(usize, usize, usize)> = Vec::new();
    let mut uniq_counts = vec![0usize; reps.len()];
    for (edges, uniq_global, (chrom, n_reps, n_reads, n_edges)) in per_chrom {
        all_edges.extend(edges);
        for (gi, cnt) in uniq_global {
            uniq_counts[gi] = cnt;
        }
        eprintln!("[gw-catalog] {chrom}: {n_reps} reps, {n_reads} reads, {n_edges} conflict edges");
    }
    for (i, c) in reps.iter_mut().enumerate() {
        c.distinguishing_uniq = uniq_counts[i];
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
    refine: bool,
    cfg: &DenovoConfig,
) -> Result<Vec<crate::vg_family::single_copy::SingleCopyLocus>> {
    let (reps, rep_totals, catalog) =
        gw_reps_and_catalog(bam_path, fasta_path, threads, win, min_copies, cfg)?;
    // Single-copy = a rep that no HOMOLOGY-GATED family claims — the SAME family definition as the catalog and
    // copy_assign. Without refining first, an FP conflict "family" (a large-gene mis-chain or repeat-bridge)
    // would wrongly EXCLUDE its genuinely-single-copy reps from the baseline. `--no-refine` skips it.
    let families: Vec<ColocatedFamily> = if refine {
        let copysets: Vec<Vec<DenovoTranscript>> = catalog.iter().map(|c| c.copies.clone()).collect();
        let refine_params =
            // RNA exon-sum substrate -> forward-only guard applies (DEFAULT ON 2026-08-19).
            RefineParams {
                threads,
                intron_fasta: Some(fasta_path.to_string()),
                require_forward_alignment: true,
                substrate: Substrate::TranscriptOriented,
                ..Default::default()
            };
        let refined = refine_families_exon_sum(copysets, &refine_params, None, cfg.conflict.min_reads)?;
        eprintln!(
            "[gw-catalog] single-copy: refined {} raw -> {} homology-gated families for the multi-copy exclusion",
            catalog.len(),
            refined.len()
        );
        refined.into_iter().enumerate().map(|(i, c)| colocated_from_copies(i, c)).collect()
    } else {
        catalog
    };
    Ok(crate::vg_family::single_copy::single_copy_loci(&reps, &rep_totals, &families))
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
    let reads = maybe_salvage_mischain(&reads, cfg).unwrap_or(reads);
    let skeletons = pass1_skeletons_robust(&reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    // `RUSTLE_DEBUG_LOCUS=chr:start-end` traces one locus through every stage, so a member that vanishes can
    // be attributed to the stage that dropped it instead of inferred. Off by default; pure logging.
    let dbg_locus: Option<(String, u64, u64)> = std::env::var("RUSTLE_DEBUG_LOCUS").ok().and_then(|v| {
        let (c, r) = v.split_once(':')?;
        let (a, b) = r.split_once('-')?;
        Some((c.to_string(), a.trim().parse().ok()?, b.trim().parse().ok()?))
    });
    if let Some((c, a, b)) = &dbg_locus {
        let n = reads.iter().filter(|r| &r.chrom == c && r.ref_start < *b && r.ref_end > *a).count();
        let mut sk: Vec<_> = skeletons.iter().filter(|s| &s.chrom == c && s.start < *b && s.end > *a).collect();
        sk.sort_by_key(|s| std::cmp::Reverse(s.n_reads));
        eprintln!("[dbg {c}:{a}-{b}] reads={n}  skeletons={}", sk.len());
        for s in sk.iter().take(6) {
            eprintln!("[dbg]   skeleton {}-{} introns={} reads={}", s.start, s.end, s.introns.len(), s.n_reads);
        }
    }
    // Support map before the free; see the same-chrom catalog for why.
    let rt_support = if cfg.filter_readthrough { Some(read_junction_support(&reads)) } else { None };
    drop(reads);
    let mut transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
    if let Some(sup) = &rt_support {
        retain_non_readthrough(&mut transcripts, sup, "gw-catalog");
        retain_non_mischain(&mut transcripts, sup, "gw-catalog");
    }
    if let Some((c, a, b)) = &dbg_locus {
        let t: Vec<&DenovoTranscript> = transcripts.iter().filter(|t| &t.chrom == c && t.start < *b && t.end > *a).collect();
        eprintln!("[dbg {c}:{a}-{b}] transcripts after gate+filters={}", t.len());
        let mut t = t; t.sort_by_key(|x| std::cmp::Reverse(x.n_reads));
        for x in t.iter().take(6) {
            eprintln!("[dbg]   tx {}-{} exons={} reads={} strand={} seq={}bp",
                      x.start, x.end, x.introns.len() + 1, x.n_reads, x.strand, x.seq.len());
        }
    }
    // Exon-union substrate (`RUSTLE_LOCUS_EXON_UNION=1`, default off = byte-identical): rebuild each locus
    // rep from the union of its group's exons rather than its single best chain. See `union_locus_reps`.
    let union_reps = std::env::var("RUSTLE_LOCUS_EXON_UNION").map(|v| v != "0" && !v.is_empty()).unwrap_or(false);
    let union_floor: u32 = std::env::var("RUSTLE_LOCUS_UNION_MIN_READS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3);
    let mut reps: Vec<DenovoTranscript> = if union_reps {
        union_locus_reps(&transcripts, &cfg.detect, &genome, union_floor)
    } else if let Some(fl) = cothread_rep_floor() {
        cothread_locus_reps(&transcripts, &cfg.detect, &genome, fl)
    } else {
        let rep_idx = collapse_loci_span_aware(&transcripts, &cfg.detect);
        rep_idx.iter().map(|&i| transcripts[i].clone()).collect()
    };
    if let Some((c, a, b)) = &dbg_locus {
        let r: Vec<&DenovoTranscript> = reps.iter().filter(|t| &t.chrom == c && t.start < *b && t.end > *a).collect();
        eprintln!("[dbg {c}:{a}-{b}] REPS after collapse={}", r.len());
        for x in r.iter().take(8) {
            eprintln!("[dbg]   rep {}-{} exons={} reads={} strand={}", x.start, x.end, x.introns.len() + 1, x.n_reads, x.strand);
        }
    }
    drop(transcripts);
    let mut by_chrom: std::collections::BTreeMap<&str, Vec<usize>> = std::collections::BTreeMap::new();
    for (gi, rep) in reps.iter().enumerate() {
        by_chrom.entry(rep.chrom.as_str()).or_default().push(gi);
    }
    eprintln!("[gw-catalog-xchrom] {} skeletons -> {} reps over {} contigs", skeletons.len(), reps.len(), contigs.len());

    // --- GLOBAL placement accumulation: a read's placements (across chroms) keyed by name, GLOBAL idx ---
    let mut name_map: std::collections::HashMap<String, Vec<Placement>> = std::collections::HashMap::new();
    // Per-rep read 3' ends, collected only when the TES extension is enabled (O(reads) extra memory).
    let tes_wanted = matches!(std::env::var("RUSTLE_TES_EXTEND"), Ok(v) if v != "0" && !v.is_empty());
    let mut tes_ends: std::collections::HashMap<usize, Vec<u64>> = std::collections::HashMap::new();
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
                // 3'-end evidence for the TES extension (see `extend_exon_sum_to_tes`). Strand-aware, taken
                // from the REP's transcription strand rather than the read flag, since a read's own strand
                // reflects library orientation and the rep is what the exon-sum is built in.
                if tes_wanted {
                    let three = if reps[*gi].strand == '-' { br.read.ref_start } else { read_end };
                    tes_ends.entry(*gi).or_insert_with(Vec::new).push(three);
                }
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
    // Record each rep's count of MAPQ>0 (unambiguous) placements as `distinguishing_uniq`, exactly as the
    // same-chrom (`gw_reps_and_catalog`) and homology (`detect_homology_catalog_genome_wide`) paths do.
    // WITHOUT this the field stays 0 on every rep here, `reads_distinguish(0, 0, ..)` is always false, and
    // `distinct_locus_reps` collapses ANY overlapping same-strand pair unconditionally -- i.e. the chi(H)
    // read-evidence guard added in 9e887b4 was inert on the --cross-chrom path. It only became visible when
    // 121b7ea made --refine (and hence `refine_families_exon_sum` -> `distinct_locus_reps`) default-on.
    //
    // Counted over the UNFILTERED name_map: the `v.len() >= 2` filter below keeps only multi-placement
    // (ambiguous) reads, which is precisely the complement of the unique mappers that are the evidence here.
    if tes_wanted {
        let win: u64 = std::env::var("RUSTLE_TES_WINDOW").ok().and_then(|v| v.parse().ok()).unwrap_or(400);
        let frac: f64 = std::env::var("RUSTLE_TES_FRAC").ok().and_then(|v| v.parse().ok()).unwrap_or(0.30);
        let mut n_set = 0usize;
        for (gi, ends) in &tes_ends {
            let fwd = reps[*gi].strand != '-';
            if let Some(t) = crate::vg_family::denovo_assemble::sharp_tes(ends, win, frac, fwd) {
                reps[*gi].tes = Some(t);
                n_set += 1;
            }
        }
        // Apply the extension HERE, where the genome is in scope. `refine_families_exon_sum` is called with
        // `None` on this path (it works from the exon-sum, not genomic spans), so it has no genome to fetch
        // the extension with — an earlier version put the call there and it silently never fired.
        let mut n_ext = 0usize;
        let mut ext_bp = 0usize;
        for i in 0..reps.len() {
            if let Some(newseq) = extend_exon_sum_to_tes(&reps[i], &genome) {
                ext_bp += newseq.len().saturating_sub(reps[i].seq.len());
                reps[i].seq = newseq;
                n_ext += 1;
            }
        }
        eprintln!("[gw-catalog-xchrom] TES extension: sharp 3' terminus on {n_set}/{} reps; \
                   extended {n_ext} exon-sums by {ext_bp} bp total", reps.len());
    }
    let all_placements: Vec<ReadPlacements> = name_map.values().cloned().collect();
    let uniq_counts = locus_unique_mapper_counts(&all_placements, reps.len());
    drop(all_placements);
    for (i, c) in reps.iter_mut().enumerate() {
        c.distinguishing_uniq = uniq_counts[i];
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
impl RefineParams {
    /// The ONE place the forward-only E_r guard is decided. It is meaningful only between
    /// transcript-oriented reps, so it is inert on a reference-oriented substrate no matter what
    /// `require_forward_alignment` says. Every orientation-sensitive option should route through a
    /// method like this rather than reading the flag directly.
    pub(crate) fn forward_only_active(&self) -> bool {
        self.require_forward_alignment && self.substrate == Substrate::TranscriptOriented
    }
}

pub fn homology_refine_params(min_identity: Option<f64>, threads: usize) -> RefineParams {
    let mut p = RefineParams { threads, ..Default::default() };
    if let Some(mi) = min_identity {
        p.min_identity = mi;
        p.sensitive_identity = mi; // BOTH tiers -> effective floor = mi (not .min())
    }
    p
}

/// The shared family-grouping ENGINE: E_r homology edges over rep sequences -> γ-quasi-clique blocks
/// (rep-index groups). Both the RNA genome-wide homology catalog and the DNA `--from-genome` path call
/// this, so "same engine, two substrates" is literally the same function. No read/annotation dependency:
/// it consumes only `rep.seq`.
/// Coverage floor for the within-family split. **DEFAULT OFF (0.0)**; set `RUSTLE_COVERAGE_SPLIT=0.90` to enable.
///
/// ⚠ SHIPPED DEFAULT-ON, THEN REVERTED. The "free purity gain" that justified default-on (over-merges
/// 22 -> 10 at unchanged recall) was an artifact of a broken benchmark: the scorer credited a truth member
/// as recovered by ANY overlapping copy, including copies whose family the split had just dissolved into
/// singletons. Counting only members still inside a >= 2-locus family, the real trade on the Soto benchmark
/// is:
///     identity 0.70 (baseline)          recall 80.7  over-merges 22  worst 7
///     identity 0.70 + coverage 0.90     recall 72.7  over-merges 10  worst 6
///     identity 0.98                     recall 76.0  over-merges  8  worst 2
/// i.e. **raising the identity floor dominates this split on all three axes**, and the split costs 8.0
/// recall points, more than the 0.98 floor's 5.8.
///
/// The cause is exactly the RNA failure mode the floor cannot distinguish from a real split: assembly
/// incompleteness. One copy assembles a full multi-exon transcript while its sibling yields only a
/// fragment, so their exon-sums cannot clear a high coverage floor and the family breaks apart. That is an
/// artifact of transcript reconstruction, not evidence the copies are unrelated.
///
/// Kept because the same lever DOES help on the DNA/Soto-replication side, where copies are whole
/// duplication units rather than transcripts (ARI 0.608 -> 0.698 with no WGS). Enable it there, not on RNA.
pub(crate) fn coverage_split_floor() -> f64 {
    std::env::var("RUSTLE_COVERAGE_SPLIT")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.0)
}

/// Split one homology block into sub-blocks whose members mutually align over at least `min_cov` of the
/// SHORTER exon-sum.
///
/// The E_r edge that BUILT the block already requires coverage >= `params.min_coverage` (0.50); this is a
/// second, stricter application of the same measurement, used to split a block that is already formed. It
/// can therefore only split, never merge, so recall can fall only when a sub-block drops below the
/// >= 2-distinct-loci bar.
///
/// Why coverage and not identity: on Soto's own catalog the groupings that need WGS copy number (famCN) to
/// separate are INVISIBLE to identity — median within-vs-between-group identity delta +0.0001, everything
/// sitting at 0.98-0.997 — but obvious in coverage (delta +0.524, with between-group coverage 0.000, i.e.
/// the subgroups frequently do not align to each other at all). Measured end-to-end on the Soto benchmark
/// at floor 0.90: real over-merges 22 -> 10 with member recall UNCHANGED at 81.8%, where instead raising the
/// identity floor to 0.98 bought 22 -> 8 at a cost of 5.8 recall points. Flat over 0.90-0.95.
/// ⚠It does NOT subsume the identity floor: the worst single fusion stays 6 here vs 2 at identity 0.98.
/// The two levers remove different failures — many small fusions vs one dense blob.
/// High-coverage edge set over ALL reps, computed ONCE. Pairs that mutually align over >= `min_cov` of the
/// shorter exon-sum, using the SAME tiers that build the E_r edge.
///
/// PERFORMANCE, learned the hard way: the first version aligned each block separately inside the family
/// loop. That is one minimap2 process per block -- thousands genome-wide -- and it took the
/// homology-primary integration test from seconds to **65 minutes**. Alignment cost belongs in one
/// genome-wide pass, exactly like `homology_edges_all_reps`, after which the per-block split is a pure set
/// operation.
pub(crate) fn coverage_edges_all_reps(
    reps: &[DenovoTranscript],
    min_cov: f64,
    params: &RefineParams,
) -> Result<BTreeSet<(usize, usize)>> {
    if min_cov <= 0.0 || reps.len() < 2 {
        return Ok(BTreeSet::new());
    }
    let seqs: Vec<Vec<u8>> = reps.iter().map(|r| r.seq.clone()).collect();
    if seqs.iter().all(|s| s.is_empty()) {
        return Ok(BTreeSet::new());
    }
    let mut edges: BTreeSet<(usize, usize)> =
        nucleotide_edges(&seqs, &["-x", "asm20"], params.min_identity, min_cov, None, params)?
            .into_iter()
            .collect();
    if params.nucleotide_sensitive {
        edges.extend(nucleotide_edges(&seqs, ER_SENSITIVE_SEED, params.sensitive_identity, min_cov, None, params)?);
    }
    Ok(edges)
}

/// Split one homology block into sub-blocks connected by `cov_edges` (from `coverage_edges_all_reps`).
///
/// The E_r edge that BUILT the block already requires coverage >= `params.min_coverage` (0.50); this is a
/// second, stricter application of the same measurement on an already-formed block, so it can only split,
/// never merge. Recall falls only when a sub-block drops below the >= 2-distinct-loci bar.
///
/// Why coverage and not identity: on Soto's own catalog the groupings that need WGS copy number (famCN) to
/// separate are INVISIBLE to identity -- median within-vs-between-group identity delta +0.0001, everything
/// at 0.98-0.997 -- but obvious in coverage (delta +0.524, between-group coverage 0.000, i.e. the subgroups
/// frequently do not align to each other at all). Measured end-to-end on the Soto benchmark at floor 0.90:
/// real over-merges 22 -> 10 with member recall UNCHANGED at 81.8%, where raising the identity floor to
/// 0.98 instead bought 22 -> 8 at a cost of 5.8 recall points. Flat over 0.90-0.95.
/// - It does NOT subsume the identity floor: worst single fusion stays 6 here vs 2 at identity 0.98.
///
/// NO SIGNAL => NO SPLIT: if no member pair carries a coverage edge, the block is returned intact. An E_r
/// edge already joined these reps, so silence means the aligner found nothing (short exon-sums produce no
/// asm20 alignment at all), not that the copies are unrelated. Treating silence as separateness shattered
/// the test fixture's legitimate 2-copy family into singletons.
pub(crate) fn coverage_split_block(
    block: &[usize],
    cov_edges: &BTreeSet<(usize, usize)>,
) -> Vec<Vec<usize>> {
    if block.len() < 2 {
        return vec![block.to_vec()];
    }
    let local: Vec<(usize, usize)> = (0..block.len())
        .flat_map(|i| ((i + 1)..block.len()).map(move |j| (i, j)))
        .filter(|&(i, j)| {
            let (a, b) = (block[i].min(block[j]), block[i].max(block[j]));
            cov_edges.contains(&(a, b))
        })
        .collect();
    if local.is_empty() {
        return vec![block.to_vec()];
    }
    let mut parent: Vec<usize> = (0..block.len()).collect();
    for (i, j) in local {
        uf_union(&mut parent, i, j);
    }
    let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
    for i in 0..block.len() {
        let r = uf_find(&mut parent, i);
        groups.entry(r).or_default().push(block[i]);
    }
    groups.into_values().collect()
}

pub(crate) fn homology_blocks(
    reps: &[DenovoTranscript],
    refine: &RefineParams,
    gamma: f64,
) -> Result<Vec<Vec<usize>>> {
    homology_blocks_pooled(reps, None, refine, gamma)
}

/// As `homology_blocks`, optionally pooling every isoform's exons into the shared-exon rule.
pub(crate) fn homology_blocks_pooled(
    reps: &[DenovoTranscript],
    pooled: Option<&[Vec<Vec<u8>>]>,
    refine: &RefineParams,
    gamma: f64,
) -> Result<Vec<Vec<usize>>> {
    Ok(homology_blocks_pooled_with_edges(reps, pooled, refine, gamma)?.0)
}

/// As `homology_blocks_pooled`, additionally returning the `E_r` edge set the blocks were cut from.
/// The blocks are bit-identical to `homology_blocks_pooled`'s — this only stops throwing the edges away,
/// because the λ certificate needs the graph, not just the partition of it.
/// TIER-2 ADMISSION (`RUSTLE_TIER2_ADMIT=1`, unset = OFF = byte-identical).
///
/// WHAT IT DOES. After reps are built, find read clusters that NO rep covers, build each one's
/// read-covered sequence, and admit it as a copy if it clears the SHIPPED `E_r` rule against an existing
/// rep. Admitted copies get a `T2~` tid prefix so `copies.tsv` distinguishes them: "assembled its own
/// model" and "admitted by similarity plus weak expression" are different claims and a single count
/// would hide that.
///
/// WHY (`docs/o1_ledger.md` §5r/§5s). The pipeline conflates "is this locus expressed?" with "can I
/// assemble a transcript model here?", and failing the second kills the first. At NPIP, 10 loci are
/// expressed and build no node; every one of them clears `E_r` against a node-bearing sibling — 10/10 on
/// oracle spans AND 10/10 on the read-cluster span the pipeline can actually obtain — while the same rule
/// admits only 4/200 non-NPIP control clusters (0.0200).
///
/// ⚠ THIS WORKS AROUND A DEFECT RATHER THAN FIXING IT. §5s topped every locus up to ~40 reads and
/// recovery still reached only 15/31: 13 of the 14 remaining loci PASS the read gate (pooled support up
/// to 43) and have NO overlapping rep, because `collapse_loci_span_aware` absorbs them and
/// `pick_locus_rep` keeps only the winner's span. The real repair is in locus collapse.
fn tier2_enabled() -> bool {
    matches!(std::env::var("RUSTLE_TIER2_ADMIT"), Ok(v) if v != "0" && !v.is_empty())
}

fn tier2_rescue(
    bam_path: &str,
    threads: usize,
    reps: &[DenovoTranscript],
    genome: &GenomeIndex,
    params: &RefineParams,
    min_reads: u32,
) -> Vec<DenovoTranscript> {
    use crate::vg_family::denovo_assemble::{build_footprint_seq, footprint_skeletons, Skeleton};
    use std::io::Write;
    // ⚠ The caller DROPS `reads` before the homology stage to bound peak memory (the catalog peaks near
    // 17 GB). Tier-2 is opt-in, so it re-reads the BAM itself rather than forcing every run to hold the
    // reads alive for a feature that is off by default.
    let reads = match primary_reads_from_bam(bam_path, threads) {
        Ok(r) => r,
        Err(e) => { eprintln!("[tier2] skipped: cannot re-read {bam_path}: {e}"); return Vec::new(); }
    };
    let reads = &reads[..];
    // Existing reps block their own regions: a locus that already has a rep needs no rescue.
    let blockers: Vec<Skeleton> = reps
        .iter()
        .map(|r| Skeleton {
            chrom: r.chrom.clone(), start: r.start, end: r.end, n_reads: min_reads.max(1),
            introns: r.introns.clone(), tied_seeded: false, footprint: false,
            read_strand: None, read_rev: 0, read_tot: 0,
        })
        .collect();
    let cands = footprint_skeletons(reads, &blockers, min_reads);
    if cands.is_empty() {
        return Vec::new();
    }
    // Sequence for each candidate, from its READ-COVERED blocks (never the whole genomic span).
    let mut cand_seq: Vec<(usize, Vec<u8>, char)> = Vec::new();
    for (i, sk) in cands.iter().enumerate() {
        if let Some((seq, strand)) = build_footprint_seq(
            genome, &sk.chrom, sk.start, sk.end, &sk.introns, sk.read_strand,
        ) {
            if seq.len() >= 300 {
                cand_seq.push((i, seq, strand));
            }
        }
    }
    if cand_seq.is_empty() {
        return Vec::new();
    }
    // Align candidates against the EXISTING reps with the shipped E_r invocation.
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    let (cp, rp) = (dir.join(format!("t2c_{pid}.fa")), dir.join(format!("t2r_{pid}.fa")));
    {
        let mut f = match std::fs::File::create(&cp) { Ok(f) => f, Err(_) => return Vec::new() };
        for (i, seq, _) in &cand_seq {
            let _ = writeln!(f, ">C{i}");
            let _ = f.write_all(seq);
            let _ = writeln!(f);
        }
        let mut g = match std::fs::File::create(&rp) { Ok(g) => g, Err(_) => return Vec::new() };
        for (j, r) in reps.iter().enumerate() {
            let _ = writeln!(g, ">R{j}");
            let _ = g.write_all(&r.seq);
            let _ = writeln!(g);
        }
    }
    // ANTI-DRIFT: the tier is assembled from its single definitions, never re-typed. `-X` is
    // deliberately EXCLUDED — tier-2 aligns candidates against reps (two files), not a set against
    // itself, so all-vs-all/`--dual=no` does not apply here.
    let out = {
        let mut cmd = std::process::Command::new(&params.minimap2);
        cmd.args(ER_TIER_FLAGS.iter().filter(|f| **f != "-X"))
            .args(ER_SENSITIVE_SEED)
            .arg("-t")
            .arg(params.threads.max(1).to_string())
            .arg(&rp)
            .arg(&cp);
        cmd.output()
    };
    let _ = std::fs::remove_file(&cp);
    let _ = std::fs::remove_file(&rp);
    let out = match out { Ok(o) => o, Err(_) => return Vec::new() };
    let mut admitted: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 { continue; }
        let (Ok(ql), Ok(qs), Ok(qe)) = (f[1].parse::<f64>(), f[2].parse::<f64>(), f[3].parse::<f64>())
        else { continue };
        let (Ok(tl), Ok(ts), Ok(te)) = (f[6].parse::<f64>(), f[7].parse::<f64>(), f[8].parse::<f64>())
        else { continue };
        let (Ok(nm), Ok(bl)) = (f[9].parse::<f64>(), f[10].parse::<f64>()) else { continue };
        if f[4] != "+" || bl <= 0.0 || ql <= 0.0 || tl <= 0.0 { continue; }
        // the SAME rule: identity >= floor AND coverage >= floor of the SHORTER, axis follows denominator.
        //
        // ⚠ THE ESTIMATOR MUST MATCH `nucleotide_edges_scored`, and it silently did not. E_r prefers the
        // gap-compressed `de:f:` tag and falls back to `nmatch/blocklen`; this site used `nmatch/blocklen`
        // unconditionally. They are NOT interchangeable — `de` charges a gap of any length once, `nm/bl`
        // charges every gap BASE — so against the same 0.60 constant the raw ratio fires on 1.2–6.9% of
        // records (down to 0.1478) that E_r admits. Tier-2 was therefore applying a materially STRICTER
        // identity test than the definition it admits into, while appearing to use the same threshold.
        // Measured and diagnosed in ledger §6i/§6m. Do not "simplify" this back to `nm / bl`.
        let de = f[12..]
            .iter()
            .find_map(|x| x.strip_prefix("de:f:").and_then(|v| v.parse::<f64>().ok()));
        let ident = match de { Some(d) => 1.0 - d, None => nm / bl };
        let cov = if ql <= tl { (qe - qs) / ql } else { (te - ts) / tl };
        if ident >= params.sensitive_identity && cov >= params.min_coverage {
            if let Some(i) = f[0].strip_prefix('C').and_then(|x| x.parse::<usize>().ok()) {
                admitted.insert(i);
            }
        }
    }
    cand_seq
        .into_iter()
        .filter(|(i, _, _)| admitted.contains(i))
        .map(|(i, seq, strand)| {
            let sk = &cands[i];
            DenovoTranscript {
                tid: format!("T2~{}_{}_{}", sk.chrom, sk.start, sk.end),
                chrom: sk.chrom.clone(), start: sk.start, end: sk.end, n_reads: sk.n_reads,
                strand, introns: sk.introns.clone(), seq, ..Default::default()
            }
        })
        .collect()
}

/// Flattening wrapper: identical partition, edges without their weights. Every pre-existing caller
/// keeps its exact signature and output; the weights are available from the `_weighted` core below.
pub(crate) fn homology_blocks_pooled_with_edges(
    reps: &[DenovoTranscript],
    pooled: Option<&[Vec<Vec<u8>>]>,
    refine: &RefineParams,
    gamma: f64,
) -> Result<(Vec<Vec<usize>>, Vec<(usize, usize)>)> {
    let (blocks, edges_w) = homology_blocks_pooled_with_edges_weighted(reps, pooled, refine, gamma)?;
    Ok((blocks, edges_w.iter().map(|&(a, b, _, _)| (a, b)).collect()))
}

/// As above, but RETURNS the per-edge identity and coverage instead of discarding them. They are already
/// computed here (see the note below); the only change is that they now reach the caller, so a REPORTED
/// certificate can carry them. This does NOT re-weight the partition — see the `RUSTLE_ER_WEIGHTED_PARTITION`
/// note below and ledger §5q/§6r: `induced_density` discards weights, γ is inert on 79% of families, and
/// identity weights < 1 lower every density, so weighting the TEST would split MORE.
pub(crate) fn homology_blocks_pooled_with_edges_weighted(
    reps: &[DenovoTranscript],
    pooled: Option<&[Vec<Vec<u8>>]>,
    refine: &RefineParams,
    gamma: f64,
) -> Result<(Vec<Vec<usize>>, Vec<(usize, usize, f64, f64)>)> {
    let edges_w = homology_edges_all_reps_pooled_weighted(reps, pooled, refine)?;
    let edges2: Vec<(usize, usize)> = edges_w.iter().map(|&(a, b, _, _)| (a, b)).collect();
    // Every edge carries weight 1.0, so the weighted machinery underneath runs UNWEIGHTED. This is a
    // deliberate choice, not an oversight, and it is worth stating because `de` (hence identity) IS
    // computed for every edge and then discarded here.
    //
    // Weighting was tested: using identity/coverage as an edge weight and re-clustering gives the same
    // partition to within noise on the Soto benchmark, because the edge SET is already the product of
    // hard identity+coverage floors -- the surviving edges have little dynamic range left to weight.
    // The lever that does move the result is WHICH edges exist (the floors, and the shared-exon rule),
    // not how strongly they are weighted. Revisit only if the edge rule is ever relaxed enough to admit
    // a genuinely graded population.
    // EDGE WEIGHTS (`RUSTLE_ER_WEIGHTED_PARTITION`, unset = OFF = every weight 1.0 = byte-identical).
    //
    // The partition currently discards identity and coverage: they decide whether an edge EXISTS and are
    // then flattened to 1.0, even though Louvain underneath is a WEIGHTED-modularity algorithm. Measured
    // on the 3-contig NPIP catalog (unit = E_r edge): NPIP<->NPIP edges have median identity 0.9878
    // against 0.8281 for NPIP<->non-NPIP, an AUC of 0.9491. Coverage is weaker at AUC 0.6418.
    //
    // ⚠ This is EDGE-LEVEL evidence, and this project has had three definitional changes pass edge or
    // node metrics and fail end-to-end (§4x, §5b, §5c). It must be judged on the PARTITION.
    // ⚠ P1 (seed-invariance) survives because weighting changes only the SCORE inside `split_once`; the
    // operation stays split-only from `all_components`, which is what the proof rests on.
    let w_mode = std::env::var("RUSTLE_ER_WEIGHTED_PARTITION").unwrap_or_default();
    let edges3: Vec<(usize, usize, f64)> = if w_mode.is_empty() || w_mode == "0" {
        edges2.iter().map(|&(a, b)| (a, b, 1.0)).collect()
    } else {
        eprintln!(
            "[partition] WEIGHTED by {} over {} edges",
            if w_mode == "coverage" { "coverage" } else { "identity" }, edges_w.len()
        );
        edges_w
            .iter()
            .map(|&(a, b, i, c)| (a, b, if w_mode == "coverage" { c } else { i }))
            .collect()
    };
    let blocks = crate::vg_family::family_split::gamma_quasi_clique_partition(reps.len(), &edges3, gamma);
    let _ = &edges2;
    Ok((blocks, edges_w))
}

/// Group a rep set into families with the shared engine (`homology_blocks`) and keep those spanning
/// >= `min_copies` spatially-distinct loci. Used by the DNA `--from-genome` path; the RNA genome-wide
/// homology catalog keeps its collapse-aware loop but shares `homology_blocks` for the grouping itself,
/// so both substrates are grouped by the same engine. `min_reads` is the `distinct_locus_reps`
/// same-strand-merge floor (0 for DNA — no reads; distinct genomic loci are still kept).
pub fn families_from_reps(
    reps: Vec<DenovoTranscript>,
    refine: &RefineParams,
    gamma: f64,
    min_copies: usize,
    min_reads: usize,
) -> Result<Vec<Vec<DenovoTranscript>>> {
    Ok(families_from_reps_certified(reps, refine, gamma, min_copies, min_reads)?.0)
}

/// THE STRUCTURAL CERTIFICATE of ONE emitted family, computed on the graph whose nodes are the family's
/// EMITTED copies (post `coverage_split_block`, post `distinct_locus_reps` — see
/// `distinct_locus_reps_grouped` for why the node set matters).
///
/// This is a REPORT, not a rule: nothing in the pipeline branches on it, no family is added, dropped or
/// reshaped by it, and it carries no threshold. Its purpose is to let a reader check a claim the catalog
/// would otherwise only assert — how much alignment evidence the family's connectedness actually rests on.
#[derive(Clone, Debug, PartialEq)]
pub struct FamilyCertificate {
    /// Emitted copies in the family (= the node count λ and `density` are computed over).
    pub n: usize,
    /// `E_r` edges induced on those copies, de-duplicated (a merged locus inherits its members' edges).
    pub n_edges: usize,
    /// `2|E| / (n(n-1))`; `1.0` for `n <= 1`.
    pub density: f64,
    /// Edge connectivity: the minimum number of `E_r` edges whose removal would split this family.
    /// **`lambda >= 2` certifies that no single alignment record's loss can split it.** `lambda == 1`
    /// says the family hangs on one edge — always true for a 2-copy family, which is exactly why this
    /// is reported rather than enforced.
    pub lambda: usize,
    /// REPORTED, per emitted copy, in the same order as `groups`: the maximum E_r identity on any edge
    /// incident to that copy within this family — the strongest evidence tying it to a sibling. `NaN`
    /// when the certificate was built without weights.
    ///
    /// Why it is emitted (ledger §5q, §6r): identity is computed for every edge, used once for the
    /// threshold, and was then discarded. It carries AUC 0.9491 for wanted-vs-suspect edges, and on a
    /// held-out catalog (356 families / 1,070 copies) it predicts where copy assignment meets
    /// alignment-score near-ties: Spearman ρ = 0.5804, permutation p = 0.0005, with a near-tie rate of
    /// 0.5496 above identity 0.95 against 0.0867 below.
    ///
    /// ⚠ It is a PRIOR, never a gate: 8.67% of reads below the cut ARE contested, so skipping that
    /// stratum would drop real work. Like λ, nothing branches on this value.
    pub copy_max_identity: Vec<f64>,
}

/// As `families_from_reps`, additionally returning one `FamilyCertificate` per emitted family, in the
/// SAME order. The family list is bit-identical to `families_from_reps`' — the certificate is derived
/// from the graph that was already built and is never consulted to decide membership.
pub fn families_from_reps_certified(
    reps: Vec<DenovoTranscript>,
    refine: &RefineParams,
    gamma: f64,
    min_copies: usize,
    min_reads: usize,
) -> Result<(Vec<Vec<DenovoTranscript>>, Vec<FamilyCertificate>)> {
    let (blocks, er_edges) = homology_blocks_pooled_with_edges_weighted(&reps, None, refine, gamma)?;
    let cov_split = coverage_split_floor();
    let cov_edges = coverage_edges_all_reps(&reps, cov_split, refine)?;
    let mut out = Vec::new();
    let mut certs = Vec::new();
    for block in blocks {
        for sub in coverage_split_block(&block, &cov_edges) {
            let copies: Vec<DenovoTranscript> = sub.iter().map(|&i| reps[i].clone()).collect();
            let grouped = distinct_locus_reps_grouped(copies, min_reads);
            let loci: Vec<DenovoTranscript> = grouped.iter().map(|(t, _)| t.clone()).collect();
            if sub.len() >= min_copies && loci.len() >= min_copies {
                let groups: Vec<Vec<usize>> = grouped
                    .iter()
                    .map(|(_, mem)| mem.iter().map(|&li| sub[li]).collect())
                    .collect();
                certs.push(certificate_for_weighted(&groups, &er_edges));
                out.push(loci);
            }
        }
    }
    Ok((out, certs))
}

/// Build a `FamilyCertificate` for a family whose emitted copies are `groups` — one entry per emitted
/// copy, listing the GLOBAL rep indices that merged into it — over the global `E_r` edge set.
///
/// An edge between two emitted copies exists iff ANY member of one aligns to ANY member of the other.
/// Edges INSIDE a merged copy are dropped (they became self-loops when the loci collapsed), and the
/// result is de-duplicated so several alignment records for one copy pair count as the ONE edge whose
/// loss would actually split the family.
/// Recompute λ certificates for families produced by a clustering OTHER than γ-quasi-clique(E_r) —
/// i.e. the `--refine` path, which re-clusters over a different edge set.
///
/// Before this existed the only safe option was to DROP the certificate (the column printed `NA`),
/// because carrying the pre-refine λ across would have described a partition the emitted row no longer
/// belongs to, and a wrong λ is worse than an absent one. This makes the certificate TOTAL instead:
/// the graph is rebuilt over exactly the rows that will be emitted, so what is reported describes what
/// is written.
///
/// Each family's transcripts are already one-per-locus by the time refine returns, so the certificate's
/// nodes are the emitted copies and the groups are singletons.
///
/// ⚠ COST: one alignment per family, the same order `--refine` already pays. It is not free, which is
/// why it runs only on the path that needs it.
pub fn certificates_for_families(
    fams: &[Vec<DenovoTranscript>],
    params: &RefineParams,
) -> Result<Vec<FamilyCertificate>> {
    let mut out = Vec::with_capacity(fams.len());
    for fam in fams {
        let groups: Vec<Vec<usize>> = (0..fam.len()).map(|i| vec![i]).collect();
        let edges = homology_edges_all_reps_pooled_weighted(fam, None, params)?;
        out.push(certificate_for_weighted(&groups, &edges));
    }
    Ok(out)
}

pub(crate) fn certificate_for(
    groups: &[Vec<usize>],
    er_edges: &[(usize, usize)],
) -> FamilyCertificate {
    let w: Vec<(usize, usize, f64, f64)> =
        er_edges.iter().map(|&(a, b)| (a, b, f64::NAN, f64::NAN)).collect();
    certificate_for_weighted(groups, &w)
}

/// As `certificate_for`, but also records each emitted copy's maximum incident E_r identity.
pub(crate) fn certificate_for_weighted(
    groups: &[Vec<usize>],
    er_edges: &[(usize, usize, f64, f64)],
) -> FamilyCertificate {
    use std::collections::{BTreeSet, HashMap};
    let n = groups.len();
    let mut owner: HashMap<usize, usize> = HashMap::new();
    for (gi, g) in groups.iter().enumerate() {
        for &r in g {
            owner.insert(r, gi);
        }
    }
    let mut induced: BTreeSet<(usize, usize)> = BTreeSet::new();
    // Per-copy maximum incident identity, reported only. An edge INSIDE one emitted copy (ga == gb) is
    // still evidence about that copy, so it counts here even though it induces no family edge.
    let mut cmax: Vec<f64> = vec![f64::NAN; n];
    for &(a, b, ident, _cov) in er_edges {
        if let (Some(&ga), Some(&gb)) = (owner.get(&a), owner.get(&b)) {
            if ident.is_finite() {
                for g in [ga, gb] {
                    cmax[g] = if cmax[g].is_nan() { ident } else { cmax[g].max(ident) };
                }
            }
            if ga != gb {
                induced.insert((ga.min(gb), ga.max(gb)));
            }
        }
    }
    let edges: Vec<(usize, usize)> = induced.iter().copied().collect();
    let n_edges = edges.len();
    let density = if n > 1 {
        2.0 * n_edges as f64 / (n as f64 * (n as f64 - 1.0))
    } else {
        1.0
    };
    FamilyCertificate {
        n,
        n_edges,
        density,
        lambda: crate::vg_family::family_split::edge_connectivity(n, &edges),
        copy_max_identity: cmax,
    }
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
) -> Result<(
    Vec<Vec<DenovoTranscript>>,
    Vec<FamilyCertificate>, // one per emitted family, SAME order as the families above
    Vec<crate::vg_family::collapse_enumerate::CollapsedFamily>,
    Vec<crate::vg_family::collapse_enumerate::ExpressedCollapsedFamily>,
    Vec<crate::vg_family::collapse_enumerate::ExpressedCollapsedFamily>, // DNA-family fallback (RNA-orphans)
)> {
    // --- reps (identical to the conflict path's rep build) ---
    let reads = primary_reads_from_bam(bam_path, threads)?;
    let contigs: HashSet<String> = reads.iter().map(|r| r.chrom.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;
    let reads = maybe_salvage_mischain(&reads, cfg).unwrap_or(reads);
    let skeletons = pass1_skeletons_robust(&reads, cfg.pass1_min_reads, cfg.min_terminal_support);
    // Support map before the free; see the same-chrom catalog for why.
    let rt_support = if cfg.filter_readthrough { Some(read_junction_support(&reads)) } else { None };
    // `RUSTLE_DEBUG_LOCUS` (comma-separated chr:start-end) traces each locus through every stage of THIS
    // path. The xchrom path has its own copy; `--homology-primary` runs HERE, and a first attempt to trace
    // through the xchrom hook produced ZERO output because of that. Pure logging; no behaviour change.
    let dbg_loci: Vec<(String, u64, u64)> = std::env::var("RUSTLE_DEBUG_LOCUS")
        .ok()
        .map(|v| {
            v.split(',')
                .filter_map(|one| {
                    let (c, r) = one.trim().split_once(':')?;
                    let (a, b) = r.split_once('-')?;
                    Some((c.to_string(), a.trim().parse().ok()?, b.trim().parse().ok()?))
                })
                .collect()
        })
        .unwrap_or_default();
    if !dbg_loci.is_empty() {
        let pooled = super::denovo_assemble::locus_support(&skeletons);
        for (c, a, b) in &dbg_loci {
            let n = reads.iter().filter(|r| &r.chrom == c && r.ref_start < *b && r.ref_end > *a).count();
            let mut sk: Vec<(usize, &super::denovo_assemble::Skeleton)> = skeletons
                .iter()
                .enumerate()
                .filter(|(_, s)| &s.chrom == c && s.start < *b && s.end > *a)
                .collect();
            sk.sort_by_key(|(_, s)| std::cmp::Reverse(s.n_reads));
            eprintln!("[dbg {c}:{a}-{b}] reads={n} skeletons={}", sk.len());
            for (i, sx) in sk.iter().take(8) {
                let sup = pooled.get(*i).copied().unwrap_or(sx.n_reads);
                let urs = matches!(std::env::var("RUSTLE_READ_STRAND"), Ok(v) if v != "0" && !v.is_empty());
                let marg: f64 = std::env::var("RUSTLE_READ_STRAND_MARGIN")
                    .ok()
                    .and_then(|v| v.parse().ok())
                    .unwrap_or(0.90);
                let why = super::denovo_assemble::gate_reject_reason(
                    sx, &genome, &cfg.gate, sup, urs, marg,
                )
                .unwrap_or("KEPT");
                eprintln!(
                    "[dbg]   skeleton {}-{} introns={} reads={} pooled={} -> {}",
                    sx.start, sx.end, sx.introns.len(), sx.n_reads, sup, why
                );
            }
        }
    }
    drop(reads);
    let mut transcripts = assemble_gate(&skeletons, &genome, &cfg.gate);
    if let Some(sup) = &rt_support {
        retain_non_readthrough(&mut transcripts, sup, "gw-catalog");
        retain_non_mischain(&mut transcripts, sup, "gw-catalog");
    }
    let union_reps = std::env::var("RUSTLE_LOCUS_EXON_UNION").map(|v| v != "0" && !v.is_empty()).unwrap_or(false);
    let union_floor: u32 = std::env::var("RUSTLE_LOCUS_UNION_MIN_READS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3);
    // Pooled per-locus exon sequences for `RUSTLE_SHARED_EXON_ISOFORMS` (default off). Captured HERE because
    // `transcripts` is dropped immediately below and the non-representative isoforms carry exon sequence the
    // representative may lack -- 46% of representatives covering a known member are single-exon stubs while
    // 53% of those loci have a gate-passing spliced chain. Min-bp is applied later, so slice with 0 here.
    let pool_isoforms = std::env::var("RUSTLE_SHARED_EXON_ISOFORMS")
        .map(|v| v != "0" && !v.is_empty())
        .unwrap_or(false);
    let mut pooled_exons: Vec<Vec<Vec<u8>>> = Vec::new();
    let mut reps: Vec<DenovoTranscript> = if union_reps {
        union_locus_reps(&transcripts, &cfg.detect, &genome, union_floor)
    } else if let Some(fl) = cothread_rep_floor() {
        cothread_locus_reps(&transcripts, &cfg.detect, &genome, fl)
    } else if pool_isoforms {
        let (rep_idx, members) =
            crate::vg_family::family_detect::collapse_loci_span_aware_with_members(&transcripts, &cfg.detect);
        pooled_exons = members
            .iter()
            .map(|ms| ms.iter().flat_map(|&i| exon_seqs_of(&transcripts[i], 0)).collect())
            .collect();
        let n_ex: usize = pooled_exons.iter().map(|v| v.len()).sum();
        eprintln!(
            "[shared-exon-isoforms] pooled {} exons over {} loci from {} isoforms (rep-only would give {})",
            n_ex,
            rep_idx.len(),
            members.iter().map(|m| m.len()).sum::<usize>(),
            rep_idx.iter().map(|&i| exon_seqs_of(&transcripts[i], 0).len()).sum::<usize>()
        );
        rep_idx.iter().map(|&i| transcripts[i].clone()).collect()
    } else {
        let rep_idx = collapse_loci_span_aware(&transcripts, &cfg.detect);
        rep_idx.iter().map(|&i| transcripts[i].clone()).collect()
    };
    for (c, a, b) in &dbg_loci {
        let t = transcripts.iter().filter(|t| &t.chrom == c && t.start < *b && t.end > *a).count();
        let mut r: Vec<&DenovoTranscript> =
            reps.iter().filter(|t| &t.chrom == c && t.start < *b && t.end > *a).collect();
        r.sort_by_key(|x| std::cmp::Reverse(x.n_reads));
        eprintln!("[dbg {c}:{a}-{b}] transcripts={t} REPS={}", r.len());
        for x in r.iter().take(6) {
            eprintln!(
                "[dbg]   rep {}-{} exons={} reads={} strand={}",
                x.start, x.end, x.introns.len() + 1, x.n_reads, x.strand
            );
        }
    }
    drop(transcripts);
    // TIER-2 ADMISSION (`RUSTLE_TIER2_ADMIT`, default off = byte-identical). Purely additive: it only
    // adds reps at read clusters no existing rep covers.
    if tier2_enabled() {
        let extra = tier2_rescue(bam_path, threads, &reps, &genome, refine,
            crate::vg_family::denovo_assemble::GATE_MIN_READS);
        if !extra.is_empty() {
            eprintln!("[tier2] admitted {} loci by similarity to an existing rep (tid prefix T2~)", extra.len());
        }
        reps.extend(extra);
    }
    eprintln!("[gw-catalog-homology] {} skeletons -> {} reps over {} contigs", skeletons.len(), reps.len(), contigs.len());

    // Per-copy unique-mapper support (Task 3, identifiability-merge): `distinct_locus_reps`'s same-strand
    // merge guard below needs read evidence, but `reads` above is `PrimaryRead` (no MAPQ) and is already
    // dropped. Re-scan the BAM for full alignment records (MAPQ) and attribute each read to its
    // closest-overlapping rep (`build_read_placements`, the same attribution `detect_and_assign` uses), then
    // record each rep's count of MAPQ>0 (unambiguous) placements as `distinguishing_uniq`. This is an
    // ADDITIONAL full-BAM pass (a genuine genome-wide cost — flagged for Task 4 perf validation); it changes
    // no existing catalog field except which same-strand co-located pairs the merge below collapses.
    let mapq_reads = aligned_reads_from_bam(bam_path, threads)?;
    let placements_for_uniq = build_read_placements(&mapq_reads, &reps);
    // Read-supported core, measured on the same pass rather than re-reading the BAM. Only consumed when
    // `RUSTLE_ER_CORE_COVERAGE` is set; computing it always keeps the value available to the audit dump.
    let cores = locus_core_bp(&mapq_reads, &reps, core_depth_floor());
    let spliced_ev = locus_has_spliced_evidence(&mapq_reads, &reps, 3);
    // Boundaries are derived from the SAME borrow, before the reads are dropped. Cloning them instead
    // roughly doubles peak RSS on a multimapper-rich BAM (94% of alignments secondary here) and made an
    // 18-minute arm thrash for an hour at 15.5 GB.
    let de_extent = locus_de_extent().map(|max_de| locus_confident_extent(&mapq_reads, &reps, max_de));
    drop(mapq_reads);
    let uniq_counts = locus_unique_mapper_counts(&placements_for_uniq, reps.len());
    drop(placements_for_uniq);
    for (i, c) in reps.iter_mut().enumerate() {
        c.distinguishing_uniq = uniq_counts[i];
        c.core_bp = cores[i];
        c.stub = c.introns.is_empty() && spliced_ev[i];
    }
    if let Some(ext) = de_extent {
        let (mut set, mut kept) = (0usize, 0usize);
        for (i, c) in reps.iter_mut().enumerate() {
            match ext[i] {
                Some((lo, hi)) => { c.start = lo; c.end = hi; set += 1; }
                None => kept += 1,
            }
        }
        eprintln!("[locus-extent] confident-read boundaries: {set} reps re-bounded, \
                   {kept} kept (no qualifying read)");
    }
    // AUDIT ONLY (`RUSTLE_LOCUS_AUDIT=1`, default silent → catalogs byte-identical): dump the rep table
    // that `distinct_locus_reps` (the ">= 2 spatially-distinct loci" certificate) is about to run on, so
    // the merge predicate can be re-derived offline. One line per rep, tab-separated.
    if std::env::var("RUSTLE_LOCUS_AUDIT").map(|v| v != "0" && !v.is_empty()).unwrap_or(false) {
        for (i, c) in reps.iter().enumerate() {
            eprintln!(
                "[rep-audit]\t{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                c.chrom, c.start, c.end, c.strand, c.n_reads, c.distinguishing_uniq,
                c.introns.len() + 1, c.seq.len()
            );
        }
    }

    // --- E_r edges + γ-quasi-clique blocks ---
    // `er_edges` is kept (not discarded as before) solely to compute each emitted family's λ certificate.
    let (blocks, er_edges) = homology_blocks_pooled_with_edges_weighted(
        &reps,
        if pooled_exons.is_empty() { None } else { Some(pooled_exons.as_slice()) },
        refine,
        gamma,
    )?;
    let n_blocks = blocks.len(); // captured before the loop below consumes `blocks` (for the diagnostic eprintln)
    // Within-family COVERAGE split (see `coverage_split_block`). Applied to the formed blocks, so it can
    // only split. `RUSTLE_COVERAGE_SPLIT=0` restores the pre-split catalog.
    let cov_split = coverage_split_floor();
    let cov_edges = coverage_edges_all_reps(&reps, cov_split, refine)?;
    let blocks: Vec<Vec<usize>> = {
        let mut v = Vec::with_capacity(blocks.len());
        for b in blocks {
            v.extend(coverage_split_block(&b, &cov_edges));
        }
        if v.len() != n_blocks {
            eprintln!("[gw-catalog-homology] coverage split (>= {cov_split}): {n_blocks} blocks -> {} sub-blocks", v.len());
        }
        v
    };

    let mut out: Vec<Vec<DenovoTranscript>> = Vec::new();
    // One entry per family pushed to `out`, in the same order (see `FamilyCertificate`).
    let mut certs: Vec<FamilyCertificate> = Vec::new();
    let mut collapsed: Vec<crate::vg_family::collapse_enumerate::CollapsedFamily> = Vec::new();
    // Collected here and projected in ONE batched minimap2 call after the loop (one genome index load
    // total), instead of re-indexing the genome per dropped candidate.
    let mut expressed_candidates: Vec<(String, String, u64, u64, Vec<u8>)> = Vec::new();
    // DNA-family fallback: the SAME dropped (<2 RNA-distinct-loci) orphan candidates, projected after the loop
    // at the looser `dna_family_min_identity` to reach divergent genomic paralogs the 0.98 collapse path misses.
    let mut dna_candidates: Vec<(String, String, u64, u64, Vec<u8>)> = Vec::new();
    for block in blocks {
        let copies: Vec<DenovoTranscript> = block.iter().map(|&i| reps[i].clone()).collect();
        // ≥2 spatially-distinct loci certificate. `locus_min_reads()`, not `cfg.conflict.min_reads`: this
        // path builds no conflict graph (see the accessor's doc comment).
        let grouped = distinct_locus_reps_grouped(copies.clone(), cfg.locus_min_reads());
        let loci: Vec<DenovoTranscript> = grouped.iter().map(|(t, _)| t.clone()).collect();
        if block.len() >= min_copies && loci.len() >= min_copies {
            // λ certificate on the EMITTED node set (post-merge), never on the pre-merge block.
            let groups: Vec<Vec<usize>> = grouped
                .iter()
                .map(|(_, mem)| mem.iter().map(|&li| block[li]).collect())
                .collect();
            certs.push(certificate_for_weighted(&groups, &er_edges));
            out.push(loci);
        } else if (cfg.collapse_enumerate || cfg.collapse_expressed || cfg.dna_family_fallback) && loci.len() < 2 {
            // GENUINE collapse only (< 2 RNA-distinct loci), independent of min_copies: a block with
            // >= 2 distinct loci that merely falls short of a higher --min-copies is a resolved (not
            // collapsed) candidate and must NOT get a mixed-locus union re-admission window.
            // dropped < 2-distinct-loci candidate → try re-admit as K=0-collapsed COPY-NUMBER (PSV/hidden-
            // copy witness) and/or K0_COLLAPSED_EXPRESSED (read-supported projection, no PSV requirement).
            // The two paths are independent flags; each fires only when its own config bit is set.
            let chrom = copies[0].chrom.clone();
            let lo = copies.iter().map(|c| c.start).min().unwrap_or(0);
            let hi = copies.iter().map(|c| c.end).max().unwrap_or(0);
            let consensus = copies.iter().max_by_key(|c| c.seq.len()).map(|c| c.seq.clone()).unwrap_or_default();
            if cfg.collapse_enumerate {
                if let Some(cf) = crate::vg_family::collapse_enumerate::readmit_locus(
                    bam_path, &chrom, lo, hi, &consensus, &genome, fasta_path, &refine.minimap2, threads,
                ) {
                    collapsed.push(cf);
                }
            }
            if cfg.collapse_expressed {
                let id = format!("exp{}", expressed_candidates.len());
                expressed_candidates.push((id, chrom.clone(), lo, hi, consensus.clone()));
            }
            if cfg.dna_family_fallback {
                let id = format!("dna{}", dna_candidates.len());
                dna_candidates.push((id, chrom.clone(), lo, hi, consensus.clone()));
            }
        }
    }
    let expressed = if cfg.collapse_expressed {
        crate::vg_family::collapse_enumerate::readmit_expressed_batch(
            &expressed_candidates, bam_path, fasta_path, &refine.minimap2, threads, 0.98,
        )
    } else {
        Vec::new()
    };
    let dna_families = if cfg.dna_family_fallback {
        crate::vg_family::collapse_enumerate::readmit_dna_family_batch(
            &dna_candidates, bam_path, fasta_path, &refine.minimap2, threads, cfg.dna_family_min_identity,
            cfg.dna_family_max_softmask,
        )
    } else {
        Vec::new()
    };
    eprintln!("[gw-catalog-homology] {} γ-quasi-clique blocks -> {} families (>= {} distinct loci)", n_blocks, out.len(), min_copies);
    if !collapsed.is_empty() {
        eprintln!("[gw-catalog-homology] collapse-enumerate: {} K=0-collapsed families re-admitted (copy-number only)", collapsed.len());
    }
    if !expressed.is_empty() {
        eprintln!("[gw-catalog-homology] collapse-expressed: {} K0_COLLAPSED_EXPRESSED families re-admitted (copy-number only)", expressed.len());
    }
    if !dna_families.is_empty() {
        eprintln!("[gw-catalog-homology] dna-family-fallback: {} DNA_FAMILY_RNA_NONHOMOLOGOUS loci re-admitted (copy-number only)", dna_families.len());
    }
    debug_assert_eq!(out.len(), certs.len(), "one certificate per emitted family, same order");
    Ok((out, certs, collapsed, expressed, dna_families))
}

/// Parameters for the exon-sum (FLNC) homology refinement. The defaults match the validated operating
/// point (`bench/validate_exon_sum.py`): minimap2 asm20, identity >= 0.80 (asm20's native divergence
/// envelope), coverage-of-shorter >= 0.50 (more than half the shorter spliced sequence aligns).
/// Which orientation convention a set of representative sequences follows.
///
/// This exists because `RefineParams` configures BOTH substrates, and orientation-sensitive options
/// mean opposite things on them. Between TRANSCRIPT-oriented reps a `-` alignment is
/// reverse-complement homology (an inverted repeat, not two homologous transcripts); between
/// REFERENCE-oriented reps it is an ordinary **inverted segmental duplication**, which is real and
/// must be kept.
///
/// Making the substrate explicit turns three scattered `require_forward_alignment = false`
/// force-offs into one declaration per entry point. On 2026-08-19 a default flip silently applied
/// the RNA guard to DNA and dropped an inverted duplication; only a test written for another purpose
/// caught it. `forward_only_active()` is the single place that decision now lives.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Substrate {
    /// RNA exon-sum, or an RNA-locus genomic span after `refine_copy_seq` normalisation — stored
    /// 5'→3' in TRANSCRIPTION direction.
    TranscriptOriented,
    /// Reference-oriented DNA intervals (`--from-genome`, the DNA arm of `--joint-dna-rna`).
    ReferenceOriented,
}

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
    /// Compute the genome-wide E_r homology edge on each rep's GENOMIC SPAN instead of its assembled
    /// exon-sum (`homology_edges_all_reps`). Fixes family FRAGMENTATION caused by incomplete/disjoint
    /// assembly — see the rationale in `homology_edges_all_reps`. Needs `intron_fasta`. Default OFF
    /// (exon-sum), so every existing catalog is byte-identical unless explicitly enabled.
    pub homology_genomic_span: bool,
    /// Require a forward (`+`) PAF orientation for nucleotide homology edges. This is meaningful only
    /// when every input sequence is in TRANSCRIPTION orientation (the RNA exon-sum and RNA-locus genomic
    /// span paths). A reverse-only match between such sequences supports reverse-complement sequence
    /// reuse, not a homologous transcript, and is enriched for inverted-repeat false positives.
    ///
    /// This must stay OFF for the annotation-free `--from-genome` arm: those reps are stored in reference
    /// orientation, so a real inverted genomic duplication legitimately aligns on `-`.
    pub require_forward_alignment: bool,
    /// Orientation convention of the reps these params will be applied to. Orientation-sensitive
    /// options are INERT unless this is `TranscriptOriented` — see `forward_only_active`.
    /// Defaults to `ReferenceOriented`, the conservative choice: a forgotten declaration costs
    /// precision, never an inverted duplication.
    pub substrate: Substrate,
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
            // Sensitive-tier identity floor. Overridable ALONE via RUSTLE_SENSITIVE_IDENTITY so the tier can
            // be A/B'd without moving `min_identity` (which --min-identity sets for BOTH tiers). Needed to
            // isolate the effect of any other change from this floor's own effect.
            sensitive_identity: std::env::var("RUSTLE_SENSITIVE_IDENTITY")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.60),
            protein_tail: false,
            mmseqs: std::env::var("RUSTLE_MMSEQS").unwrap_or_else(|_| "mmseqs".to_string()),
            homology_genomic_span: false,
            // ⚠ STAYS FALSE, and the 2026-08-19 default flip deliberately did NOT change it.
            // `RefineParams` is SUBSTRATE-AGNOSTIC: the same struct configures the RNA exon-sum path
            // and the reference-oriented DNA path (`--from-genome`, genome-mode grouping). A '-' record
            // is reverse-complement homology between TRANSCRIPT-oriented reps, but a REAL inverted
            // segmental duplication between REFERENCE-oriented ones. Flipping this default silently
            // applied the RNA guard to DNA and dropped an inverted duplication — caught by
            // `genome_mode_grouping_keeps_an_inverted_duplication`. The guard is therefore turned on at
            // the RNA ENTRY POINTS, never here. See docs/o1_investigations.md#false-positive-hardening-rules-that-survived-falsification.
            require_forward_alignment: false,
            // Conservative default: orientation-sensitive options stay inert until an entry
            // point declares the substrate. Forgetting costs precision, never correctness.
            substrate: Substrate::ReferenceOriented,
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
/// Do the `edges` connect all `n` nodes into a single component? (union-find). Used to short-circuit refine's
/// additive homology tiers when the asm20 core already fully connects a family.
/// Which E_r tier produced an edge. An edge can be found by more than one tier, so this is a MASK and not
/// an enum: `tier_names` reports every tier that independently supports it. Recorded rather than inferred,
/// so "this family was discovered by the sensitive tier" is an auditable fact — the union of tiers is
/// otherwise indistinguishable from a single-tier result once the edges are merged.
pub type TierMask = u8;
pub const TIER_ASM20: TierMask = 1 << 0;
pub const TIER_SENSITIVE: TierMask = 1 << 1;
pub const TIER_GENOMIC: TierMask = 1 << 2;
pub const TIER_PROTEIN: TierMask = 1 << 3;

/// Human-readable tier list for a mask, e.g. `asm20+sensitive`. Empty mask -> "none".
pub fn tier_names(m: TierMask) -> String {
    let mut v: Vec<&str> = Vec::new();
    if m & TIER_ASM20 != 0 {
        v.push("asm20");
    }
    if m & TIER_SENSITIVE != 0 {
        v.push("sensitive");
    }
    if m & TIER_GENOMIC != 0 {
        v.push("genomic-span");
    }
    if m & TIER_PROTEIN != 0 {
        v.push("protein");
    }
    if v.is_empty() {
        "none".to_string()
    } else {
        v.join("+")
    }
}

fn edges_connect_all(n: usize, edges: &BTreeSet<(usize, usize)>) -> bool {
    if n <= 1 {
        return true;
    }
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        let mut r = x;
        while p[r] != r {
            r = p[r];
        }
        let mut c = x;
        while p[c] != r {
            let nx = p[c];
            p[c] = r;
            c = nx;
        }
        r
    }
    for &(a, b) in edges {
        let ra = find(&mut parent, a);
        let rb = find(&mut parent, b);
        parent[ra] = rb;
    }
    let r0 = find(&mut parent, 0);
    (1..n).all(|i| find(&mut parent, i) == r0)
}

pub fn refine_families_exon_sum(
    families: Vec<Vec<DenovoTranscript>>,
    params: &RefineParams,
    passed_genome: Option<&GenomeIndex>,
    min_reads: usize,
) -> Result<Vec<Vec<DenovoTranscript>>> {
    // Genome for the include_introns core tier + the additive genomic-span tier. PREFER the caller's
    // already-loaded `GenomeIndex` (the per-region path caches contigs in an LRU precisely to avoid reloads); only
    // load from `intron_fasta` when none was passed. This eliminates a full-chromosome reload per region. For
    // include_introns the genome MUST be available; for the additive genomic tier a load failure is best-effort
    // (skip the tier — e.g. a synthetic test contig absent from the fasta).
    let owned_genome: Option<GenomeIndex> = if passed_genome.is_some() {
        None
    } else {
        match params.intron_fasta.as_ref() {
            Some(path) => {
                let contigs: HashSet<String> = families.iter().flatten().map(|c| c.chrom.clone()).collect();
                match GenomeIndex::from_fasta_contigs(path, &contigs) {
                    Ok(g) => Some(g),
                    Err(e) if params.include_introns => {
                        return Err(e.context("include_introns: genome load failed"));
                    }
                    Err(e) => {
                        eprintln!("[refine] genomic-span tier skipped: genome load failed ({e})");
                        None
                    }
                }
            }
            None if params.include_introns => {
                anyhow::bail!("include_introns set but intron_fasta (genome path) is None")
            }
            None => None,
        }
    };
    let genome: Option<&GenomeIndex> = passed_genome.or(owned_genome.as_ref());
    // PROTEIN tier (opt-in): one batched mmseqs run across ALL families' ORFs (within-family hits only),
    // instead of a subprocess per family. `fam_protein[f]` = the protein homology edges of family f.
    let fam_protein: Vec<Vec<(usize, usize)>> = if params.protein_tail {
        batch_protein_edges(&families, 0.50, params.min_coverage, params)?
    } else {
        Vec::new()
    };
    // ⭐ X.4 — ONE TIER. The primary tier is resolved ONCE, from the shared `er_primary_tier`, and reused
    // for every family and for the genomic-span tier below. Before this, both sites called
    // `primary_seed_args()` directly and so ignored `RUSTLE_ER_SENSITIVE_ONLY` entirely.
    let (prim_seed, prim_floor, prim_mask) = er_primary_tier(params);
    let prim_seed_ref: Vec<&str> = prim_seed.iter().map(String::as_str).collect();
    // ⭐ O-3 — THE ONE DECISION about the additive genomic-span leg, taken ONCE, here. The gate below
    // branches on it and the certificate at the bottom of this function is built from THIS value, so
    // `params.tsv` cannot claim a substrate the run did not use. It used to be a nested `if` no observer
    // could read, which is why the certificate printed `substrate = exon-sum` for an `E_x ∪ E_g` run.
    let genomic_tier = additive_genomic_tier(params, genome.is_some());
    // The substrate the CORE tier aligns. `refine_copy_seq(_, core_genome)` is what actually decides it,
    // so it is named from the same condition rather than restated.
    let core_substrate = substrate_name(params.include_introns);
    eprintln!(
        "[refine] E_r primary tier: {} @ identity {:.4} (sensitive_only={}) | core substrate: {} | additive genomic tier: {}",
        prim_seed.join(" "),
        prim_floor,
        er_sensitive_only(params),
        core_substrate,
        genomic_tier.label()
    );
    // Aggregates for the `<prefix>.refine.params.tsv` instrument (below). Data only; the RULE rows come
    // from `er_rule_rows` and are byte-comparable against the O1 catalog's.
    let mut dump_reps = 0usize;
    let mut dump_edges = 0usize;
    let mut dump_lens: Vec<usize> = Vec::new();
    // ⭐ UNION AUDIT (`RUSTLE_UNION_AUDIT=<path>`) — the instrument for "how often does the additive
    // genomic-span tier actually fire, and does it move anything?". PURE ACCOUNTING: every value below is
    // derived from data the function already computes, nothing feeds back into a decision, and the whole
    // block writes only when the env var is set. Behaviour-neutrality is asserted by the control panel
    // staying byte-identical with the var set and unset.
    let audit_on = std::env::var("RUSTLE_UNION_AUDIT").map(|v| !v.is_empty()).unwrap_or(false);
    let mut au_fam_examined = 0usize; // families with >= 2 copies (the gate's denominator)
    let mut au_gate_false = 0usize; // ... where edges_connect_all(primary) was FALSE
    let mut au_genomic_calls = 0usize; // ... where the genomic-span tier actually ran (genome present)
    let mut au_genomic_skipped_nogenome = 0usize; // gate open but no genome => tier could not run
    let mut au_fam_genomic_added = 0usize; // families where the genomic leg added >= 1 NEW edge
    let mut au_edges_genomic_new = 0usize; // NEW edges contributed by the genomic leg
    let mut au_edges_primary_total = 0usize;
    let mut au_fam_comps_changed = 0usize; // families whose COMPONENT PARTITION moved under the union
    let mut au_comps_primary_total = 0usize;
    let mut au_comps_union_total = 0usize;
    let mut au_fam_connected_by_union = 0usize; // gate-false families the union made fully connected
    let mut au_rows: Vec<String> = Vec::new();
    let mut refined: Vec<Vec<DenovoTranscript>> = Vec::new();
    for (fi, fam) in families.iter().enumerate() {
        if fam.len() < 2 {
            continue;
        }
        au_fam_examined += 1;
        // core tier: exon-sum (spliced) for a clean, intron-length-independent identity; genomic span only when
        // include_introns is explicitly set (the genomic-span tier below is the additive default-on path).
        let core_genome = if params.include_introns { genome } else { None };
        let seqs: Vec<Vec<u8>> = fam.iter().map(|c| refine_copy_seq(c, core_genome)).collect();
        dump_reps += fam.len();
        dump_lens.extend(seqs.iter().map(|s| s.len()));
        // base detector: the PRIMARY tier on the configured sequence.
        // Provenance: every edge carries the mask of tiers that independently produced it, so a family's
        // discovery path is recorded rather than reconstructed from the merged edge set.
        let mut prov: BTreeMap<(usize, usize), TierMask> = BTreeMap::new();
        {
            // The dump tag carries the CORE substrate's name (one spelling, `substrate_name`), so a
            // `<prefix>.refine.core.<substrate>.*.paf` can never be confused with the additive leg's.
            for e in nucleotide_edges_tagged(
                &seqs,
                &prim_seed_ref,
                prim_floor,
                params.min_coverage,
                None,
                params,
                Some(&format!("refine.core.{core_substrate}")),
            )? {
                *prov.entry(e).or_insert(0) |= prim_mask;
            }
        }
        let mut edge_set: BTreeSet<(usize, usize)> = prov.keys().copied().collect();
        // Snapshot of the PRIMARY (exon-sum core) edge set, before any additive tier unions into it. This is
        // the object every "union is a no-op" claim on record was actually measured on; keeping it lets the
        // audit compare the two partitions on the SHIPPED object instead of on an externally built PAF.
        let au_primary: BTreeSet<(usize, usize)> = if audit_on { edge_set.clone() } else { BTreeSet::new() };
        au_edges_primary_total += edge_set.len();
        // The additive tiers only UNION more edges. If the asm20 exon-sum core already connects every copy into a
        // single homology component, they are a provable no-op (`homology_components` on a superset of a connected
        // edge set gives the same partition, and `distinct_locus_reps` runs unchanged). So skip them — and the
        // genome fetch the genomic-span tier needs — for the common fully-connected family. This is the case for
        // most near-identical segdup families; the tiers (and the genome touch) run only when asm20 left a gap.
        if !edges_connect_all(fam.len(), &edge_set) {
            // GENOMIC-SPAN tier (recall fix; evidence + denominators in `bench/FALSE_NEGATIVES.md`, which
            // was DELETED in 9b0814f and RESTORED by O-4 — this citation dangled for a month, leaving a
            // default-ON edge rule with no justification in the tree. O-4 2026-08-13 re-measured the leg
            // on the SHIPPED object and KEPT it ON: 11/26 families' partitions moved, yet block sets vs
            // DNA-only differ 0/26 and the emitted catalogs are identical (ARI 1.0000, 0 forbidden pairs).
            // It is a recall gain with no measured precision cost BECAUSE it is gated and family-local.)
            // A near-identical segdup whose PARTIAL
            // de-novo transcript models fail the exon-sum coverage floor still covers most of its GENOMIC extent
            // at the same (gap-compressed) identity; a repeat-bridge covers < min_coverage of the genomic span
            // regardless of the repeat's identity, so this does NOT readmit bridges/gene-splits. Skipped when
            // include_introns already ran the core tier on genomic.
            au_gate_false += 1;
            // ⭐ O-3 — the gate reads `genomic_tier`, the SAME value the certificate is built from. The
            // arms below are exactly the old nested `if !params.include_introns { if let Some(g) = … }`:
            // `CoreIsGenomic` and `NotPresent` fall through silently as the outer `if` did, `NoGenome` is
            // the old `else`. Behaviour-identical by construction; what changed is that it is now READABLE.
            if genomic_tier.armed() {
                let g = genome.expect("GenomicTier::Armed is returned only when a genome is available");
                au_genomic_calls += 1;
                let gseqs: Vec<Vec<u8>> = fam.iter().map(|c| refine_copy_seq(c, Some(g))).collect();
                // SAME primary tier as the core run — a different SUBSTRATE, never a different rule. The
                // dump tag names that substrate, so the `.paf` / `.args` this leg drops can no longer be
                // mistaken for the core run's (before O-3 the two were byte-identical but for a counter).
                let mut au_new_here = 0usize;
                for e in nucleotide_edges_tagged(
                    &gseqs,
                    &prim_seed_ref,
                    prim_floor,
                    params.min_coverage,
                    None,
                    params,
                    Some(&format!("refine.additive.{SUBSTRATE_GENOMIC_SPAN}")),
                )? {
                    // `insert` returns true iff the edge was NOT already present, i.e. iff the genomic
                    // leg contributed an edge the exon-sum core did not have. That is the quantity the
                    // whole run turns on, so it is counted here rather than inferred from the tier mask
                    // (a mask can be set on an edge the primary tier already found).
                    if edge_set.insert(e) {
                        au_new_here += 1;
                    }
                    *prov.entry(e).or_insert(0) |= TIER_GENOMIC;
                }
                au_edges_genomic_new += au_new_here;
                if au_new_here > 0 {
                    au_fam_genomic_added += 1;
                }
            } else if genomic_tier == GenomicTier::NoGenome {
                au_genomic_skipped_nogenome += 1;
            }
            // divergent tiers on the EXON-SUM (protein ORFs need the spliced sequence; the sensitive nucleotide
            // seed is cleanest on the spliced copy). Edges are UNIONed in.
            //
            // ⚠ Under sensitive-only on the exon-sum substrate the PRIMARY tier already IS this run — same
            // seed, same floor, same sequences, same `cores: None` — so re-running it can only reproduce the
            // edges already in `prov`. Skipped, which removes one minimap2 subprocess per unconnected family
            // and cannot move the edge set. (With `include_introns` the core ran on GENOMIC sequence, so the
            // exon-sum run is still a distinct substrate and does run.)
            if params.nucleotide_sensitive && !(er_sensitive_only(params) && !params.include_introns) {
                let exon_seqs: Vec<Vec<u8>> = fam.iter().map(|c| c.seq.clone()).collect();
                // Identity floor comes from `params.sensitive_identity` (the SAME knob as the E_r path).
                // This was a hard-coded 0.70 that `--min-identity` could not reach, so `--min-identity 0.98`
                // silently admitted 0.70 edges here. One tier, one knob.
                for e in nucleotide_edges(
                    &exon_seqs,
                    ER_SENSITIVE_SEED,
                    params.sensitive_identity,
                    params.min_coverage,
                    None,
                    params,
                )? {
                    edge_set.insert(e);
                    *prov.entry(e).or_insert(0) |= TIER_SENSITIVE;
                }
            }
            if params.protein_tail {
                for &e in fam_protein.get(fi).map(|v| v.as_slice()).unwrap_or(&[]) {
                    edge_set.insert(e);
                    *prov.entry(e).or_insert(0) |= TIER_PROTEIN;
                }
            }
        }
        // ⭐ THE HINGE MEASUREMENT, per family: does the union change the PARTITION the pipeline goes on to
        // emit? Compared as component SIGNATURES (each component's sorted member list), because a bare
        // component COUNT can stay equal while membership moves. Read-only: `edges` below is unaffected.
        if audit_on {
            let sig = |es: &BTreeSet<(usize, usize)>| -> Vec<Vec<usize>> {
                let v: Vec<(usize, usize)> = es.iter().copied().collect();
                let mut cs = homology_components(fam.len(), &v);
                for c in cs.iter_mut() {
                    c.sort_unstable();
                }
                cs.sort();
                cs
            };
            let cp = sig(&au_primary);
            let cu = sig(&edge_set);
            au_comps_primary_total += cp.len();
            au_comps_union_total += cu.len();
            let moved = cp != cu;
            if moved {
                au_fam_comps_changed += 1;
            }
            if au_primary.len() != edge_set.len() && edges_connect_all(fam.len(), &edge_set) {
                au_fam_connected_by_union += 1;
            }
            let span = &fam[0];
            au_rows.push(format!(
                "{}\t{}:{}-{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                fi,
                span.chrom,
                span.start,
                span.end,
                fam.len(),
                au_primary.len(),
                edge_set.len(),
                edge_set.len() - au_primary.len(),
                cp.len(),
                cu.len(),
                if moved { "PARTITION_MOVED" } else { "same" }
            ));
        }
        let edges: Vec<(usize, usize)> = edge_set.into_iter().collect();
        dump_edges += edges.len();
        for comp in homology_components(fam.len(), &edges) {
            if comp.len() < 2 {
                continue;
            }
            // Provenance of THIS component: the union of its internal edges' tier masks, plus the tiers
            // without which it would fall apart. An edge found only by one tier is load-bearing for that
            // component; report it so "discovered via the sensitive tier" is checkable.
            let cset: BTreeSet<usize> = comp.iter().copied().collect();
            let mut mask: TierMask = 0;
            let mut sole: BTreeMap<TierMask, usize> = BTreeMap::new();
            for (&(a, b), &m) in prov.iter() {
                if cset.contains(&a) && cset.contains(&b) {
                    mask |= m;
                    if m.count_ones() == 1 {
                        *sole.entry(m).or_insert(0) += 1;
                    }
                }
            }
            let comp_copies: Vec<DenovoTranscript> = comp.iter().map(|&i| fam[i].clone()).collect();
            let loci = distinct_locus_reps(comp_copies, min_reads);
            if loci.len() >= 2 {
                if let Some(span) = loci.first() {
                    let uniq: Vec<String> =
                        sole.iter().map(|(&m, &n)| format!("{}={}", tier_names(m), n)).collect();
                    eprintln!(
                        "[provenance] family @ {}:{}-{} ({} loci): tiers={} sole-support[{}]",
                        span.chrom,
                        span.start,
                        span.end,
                        loci.len(),
                        tier_names(mask),
                        uniq.join(" ")
                    );
                }
                refined.push(loci);
            }
        }
    }
    // ⭐ X.4(2) — THE INSTRUMENT NOW FIRES INSIDE `copy_assign`. `write_er_edge_dump` is reachable only
    // from `homology_edges_all_reps_pooled`, which O2's refine never calls; measured before this, 25/25
    // params files were written on the O1 side and 0/25 on the O2 side, so "diff two params.tsv files"
    // was an O1-only capability. Refine now emits its own pair under the SAME `RUSTLE_ER_EDGE_DUMP`
    // prefix, with the rule rows produced by the SHARED `er_rule_rows`.
    // ⭐ UNION AUDIT emission. Appends (never truncates) so one file accumulates a whole panel, and every
    // call of this function is its own `call=` block — a driver that runs 25 regions gets 25 blocks in one
    // file. Denominators are PRE-DECLARED here: `fam_examined` is every input family with >= 2 copies, and
    // it is printed whether or not the tier ever fired.
    if audit_on {
        if let Ok(path) = std::env::var("RUSTLE_UNION_AUDIT") {
            use std::io::Write;
            static UNION_AUDIT_SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
            let call = UNION_AUDIT_SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            if let Ok(mut f) = std::fs::OpenOptions::new().create(true).append(true).open(&path) {
                let _ = writeln!(
                    f,
                    "#call\t{call}\tlabel\t{}",
                    std::env::var("RUSTLE_UNION_AUDIT_LABEL").unwrap_or_else(|_| "NA".into())
                );
                let _ = writeln!(
                    f,
                    "SUMMARY\tcall={call}\tfam_examined={au_fam_examined}\tgate_false={au_gate_false}\
                     \tgenomic_calls={au_genomic_calls}\tgenomic_skipped_nogenome={au_genomic_skipped_nogenome}\
                     \tfam_genomic_added_edges={au_fam_genomic_added}\tedges_genomic_new={au_edges_genomic_new}\
                     \tedges_primary_total={au_edges_primary_total}\tfam_partition_moved={au_fam_comps_changed}\
                     \tcomps_primary={au_comps_primary_total}\tcomps_union={au_comps_union_total}\
                     \tfam_connected_by_union={au_fam_connected_by_union}\tinclude_introns={}\tgenome_present={}",
                    params.include_introns,
                    genome.is_some()
                );
                for r in &au_rows {
                    let _ = writeln!(f, "FAM\tcall={call}\t{r}");
                }
            }
        }
    }
    if let Ok(prefix) = std::env::var("RUSTLE_ER_EDGE_DUMP") {
        if !prefix.is_empty() {
            static REFINE_DUMP_SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
            let n = REFINE_DUMP_SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let stem =
                if n == 0 { format!("{prefix}.refine") } else { format!("{prefix}.refine.call{}", n + 1) };
            // ⭐ O-3 — `genomic_tier` is THE VARIABLE THE GATE BRANCHED ON, not a second derivation of the
            // same condition. That is the whole point: a certificate re-deriving its own answer can be
            // right about the code and wrong about the run.
            let site = ErRuleSite {
                substrate_genomic: params.include_introns,
                core_lens_supplied: false,
                genomic_tier,
            };
            let rule = er_rule_rows(params, &site);
            write_kv_tsv(&format!("{stem}.rule.tsv"), &rule);
            // params.tsv = the data rows THEN the rule rows, so it stays a superset of `.rule.tsv` on both
            // sides. Only `.rule.tsv` is the parity object; these counts are input-dependent by design.
            let mut rows: Vec<(String, String)> = vec![
                ("site".into(), "refine_families_exon_sum".into()),
                ("n_families_in".into(), families.len().to_string()),
                ("n_families_out".into(), refined.len().to_string()),
                ("n_reps".into(), dump_reps.to_string()),
                ("n_edges".into(), dump_edges.to_string()),
                ("min_identity_asm20".into(), format!("{:.6}", params.min_identity)),
                ("sensitive_identity".into(), format!("{:.6}", params.sensitive_identity)),
                ("nucleotide_sensitive".into(), params.nucleotide_sensitive.to_string()),
                ("protein_tail".into(), params.protein_tail.to_string()),
                ("include_introns".into(), params.include_introns.to_string()),
                ("threads".into(), params.threads.to_string()),
                // ⭐ O-3 — DID THE ADDITIVE GENOMIC LEG ACTUALLY FIRE, AND WHAT DID IT CONTRIBUTE?
                // `additive_genomic_tier` (a RULE row, appended below) says whether it was ARMED; these
                // say what happened. Both are needed: an armed leg that never fires and a leg that
                // rebuilt the whole edge set print the same rule row and must not print the same run.
                ("n_families_examined".into(), au_fam_examined.to_string()),
                ("n_families_core_unconnected".into(), au_gate_false.to_string()),
                ("n_families_genomic_tier_ran".into(), au_genomic_calls.to_string()),
                ("n_families_genomic_tier_skipped_nogenome".into(), au_genomic_skipped_nogenome.to_string()),
                ("n_families_genomic_tier_added_edges".into(), au_fam_genomic_added.to_string()),
                // The number every "the union is a no-op" claim is about: edges present in `E_x ∪ E_g`
                // and NOT in the exon-sum core's `E_x`. 0 here is the no-op, measured on the object the
                // binary actually emits.
                ("n_edges_genomic_tier_added".into(), au_edges_genomic_new.to_string()),
                ("n_edges_core_tier".into(), au_edges_primary_total.to_string()),
                ("additive_genomic_tier_fired".into(), (au_edges_genomic_new > 0).to_string()),
                ("paf_glob_for_this_edge_set".into(), format!("{prefix}.refine.*.paf")),
                ("substrate_median_len_bp".into(), {
                    median_len(&dump_lens).map(|m| m.to_string()).unwrap_or_else(|| "NA".into())
                }),
                ("coverage_floor_median_bp_demand".into(), {
                    coverage_floor_bp_demand(params.min_coverage, &dump_lens)
                        .map(|b| b.to_string())
                        .unwrap_or_else(|| "NA".into())
                }),
            ];
            rows.extend(rule);
            write_kv_tsv(&format!("{stem}.params.tsv"), &rows);
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
                s = crate::vg_family::seq_utils::reverse_complement(&s);
            }
            s
        }
        None => copy.seq.clone(), // already carries the TES extension when RUSTLE_TES_EXTEND is on
    }
}

/// Append the terminal-exon extension (copy `end`/`start` → observed TES) to the exon-sum. Opt-in via
/// `RUSTLE_TES_EXTEND`; without it this returns the exon-sum unchanged, so every existing catalog is
/// byte-identical.
///
/// WHY the 3' end specifically. The exon-sum ends at the k-th-read boundary quantile, which is deliberately
/// conservative. For a TRUNCATED duplicate that is exactly the wrong place to stop: NOTCH2NLB's transcript
/// terminates ~58 kb from NOTCH2's, ANAPC1P2's ~59 kb from ANAPC1's, and that difference is the clearest
/// signal separating the parent from the copy. Measured over the Soto benchmark, 14/42 sibling pairs have
/// genuinely distinct 3' termini against only 3/40 for 5' termini, at a median 4.9 kb
/// (`bench/soto/bam_tie_signals.md` §9), so this adds real discriminating sequence rather than noise.
///
/// Fires only where `copy.tes` is set, i.e. where the 3'-end distribution was SHARP (`sharp_tes`); a broad
/// distribution is differential coverage and extending to it would just absorb the furthest stray read.
pub(crate) fn extend_exon_sum_to_tes(copy: &DenovoTranscript, g: &GenomeIndex) -> Option<Vec<u8>> {
    let tes = copy.tes?;
    // Extension window in genomic coordinates, in the 3' direction for this strand. An inward or flush TES
    // yields no window — this may only ever ADD sequence, never trim.
    let (lo, hi) = if copy.strand == '-' {
        if tes >= copy.start { return None; }
        (tes, copy.start)
    } else {
        if tes <= copy.end { return None; }
        (copy.end, tes)
    };
    let mut ext = g.fetch_sequence(&copy.chrom, lo, hi)?;
    if ext.is_empty() {
        return None;
    }
    if copy.strand == '-' {
        ext = crate::vg_family::seq_utils::reverse_complement(&ext);
    }
    let mut out = copy.seq.clone();
    out.extend_from_slice(&ext); // exon-sum is in transcription orientation, so the 3' extension appends
    Some(out)
}

/// Minimap2 seeding arguments for the PRIMARY homology tier.
///
/// Default `-x asm20`, which presets k=19/w=10. That seeding is the binding constraint far more often than
/// the identity threshold is: measured all-vs-all on curated gene families
/// (`bench/soto/preset_sensitivity.py`), asm20 finds 8 of 66 SIGLEC member pairs where `-k 9 -w 3` finds 60,
/// 64/153 H2BC pairs vs 153, 18/36 TUBA pairs vs 36. Critically the pairs asm20 misses are not marginal —
/// most are ALREADY above the 0.80 identity floor (SIGLEC 8 -> 31 pairs at >= 0.80, H2BC 64 -> 141). Short
/// exon-sums suffer worst, since k=19/w=10 leaves too few anchors, which is why family recovery correlated
/// with sequence length in `merge_quality_analysis.md` §19b.
///
/// `RUSTLE_MM2_SEED` overrides it, e.g. `RUSTLE_MM2_SEED="-k 9 -w 3"`.
pub(crate) fn primary_seed_args() -> Vec<String> {
    match std::env::var("RUSTLE_MM2_SEED") {
        Ok(v) if !v.trim().is_empty() => v.split_whitespace().map(str::to_string).collect(),
        _ => vec!["-x".to_string(), "asm20".to_string()],
    }
}

/// ─────────────────────────────────────────────────────────────────────────────────────────────────────
/// THE SHIPPED E_r ALIGNMENT TIER — **ONE DEFINITION, NO SECOND COPY**.
///
/// Every all-vs-all homology alignment in this crate is `minimap2 -c -X --no-long-join -t <threads>
/// <preset...>`. These three flags were previously hardcoded at FOUR sites (two `Command` builders, the
/// `RUSTLE_ER_EDGE_DUMP` `.args` line, and the `params.tsv` `mm_args_sensitive` row). That duplication is
/// exactly how the tier drifted away from the eight bench/crossspecies panel scripts, which ran
/// `-N 200 -p 0.02` with NO `-X` and therefore built a DIFFERENT graph on byte-identical FASTA (partition
/// differs on 4/14 panels, edge count on 10/14). `-N`/`-p` are INERT at this tier; `-X` is the operative
/// difference, because it implies `--dual=no`.
///
/// ⚠ `-X` means ONE orientation per pair is emitted and the query is NOT necessarily the shorter sequence
/// (measured: the query is the LONGER sequence in ~60% of records). Any per-record statistic computed here
/// must therefore choose its AXIS explicitly — see `ER_COVERAGE_FORM` and the coverage block below.
///
/// The VALUE is unchanged from the shipped binary; only the duplication is removed.
pub(crate) const ER_TIER_FLAGS: &[&str] = &["-c", "-X", "--no-long-join"];

/// Human-readable statement of the coverage form applied to every E_r record. Emitted into
/// `params.tsv` so a run self-reports its definition (diff two `params.tsv` files, never read code).
pub(crate) const ER_COVERAGE_FORM: &str = "aligned span on the SHORTER sequence / len(shorter) \
    (ql<=tl ? (qe-qs)/ql : (te-ts)/tl); axis follows the denominator";

/// Build the shipped-tier `minimap2` invocation. **The only place a `Command` for E_r is constructed.**
pub(crate) fn er_tier_command(minimap2: &str, threads: usize, mm_args: &[&str]) -> std::process::Command {
    let mut cmd = std::process::Command::new(minimap2);
    cmd.args(ER_TIER_FLAGS).arg("-t").arg(threads.max(1).to_string()).args(mm_args);
    cmd
}

/// The SENSITIVE seeding tier (`-k 11 -w 5`), previously re-typed at five sites. Paired with
/// `RefineParams::sensitive_identity`, not `min_identity`.
pub(crate) const ER_SENSITIVE_SEED: &[&str] = &["-k", "11", "-w", "5"];

/// ⭐ IS THE SHIPPED SINGLE-TIER DEFAULT ACTIVE? (`RUSTLE_ER_SENSITIVE_ONLY`, default **on**.)
///
/// **X.4 — THIS PREDICATE USED TO EXIST IN EXACTLY ONE FUNCTION.** The 2026-08-07 "sensitive-only
/// default" was written inline inside `homology_edges_all_reps_pooled`, so it governed the O1 catalog
/// edge and NOTHING ELSE. `refine_families_exon_sum` — the edge step O2 (`copy_assign`) actually runs —
/// called `primary_seed_args()` unconditionally at two sites and therefore kept running `-x asm20` at
/// `min_identity` (0.80) after the flip. Measured before this fix: O1 ran the sensitive tier on 25/25
/// panel calls while O2's refine ran asm20@0.80 on 10 of 13. "One tier" was true of one caller.
///
/// The predicate now lives here and both callers read it, so the tier cannot drift per call site again.
pub(crate) fn er_sensitive_only(params: &RefineParams) -> bool {
    params.nucleotide_sensitive
        && std::env::var("RUSTLE_ER_SENSITIVE_ONLY").map(|v| v != "0" && !v.is_empty()).unwrap_or(true)
}

/// The PRIMARY E_r tier as `(seed args, identity floor, tier bit)` — **the single definition of "what
/// alignment decides an edge"**, shared by O1's `homology_edges_all_reps_pooled` and O2's
/// `refine_families_exon_sum`.
///
/// Sensitive-only (the default): `-k 11 -w 5` at `sensitive_identity`. Otherwise the legacy primary:
/// `primary_seed_args()` (`-x asm20`, or `RUSTLE_MM2_SEED`) at `min_identity`.
pub(crate) fn er_primary_tier(params: &RefineParams) -> (Vec<String>, f64, TierMask) {
    if er_sensitive_only(params) {
        (ER_SENSITIVE_SEED.iter().map(|s| (*s).to_string()).collect(), params.sensitive_identity, TIER_SENSITIVE)
    } else {
        (primary_seed_args(), params.min_identity, TIER_ASM20)
    }
}

/// ⭐ **THE TWO SUBSTRATE NAMES, ONCE.** Spelled here and nowhere else, so the certificate row, the
/// `[refine]` log line and the `.args` / PAF filenames of a dump can never disagree about which sequence
/// an alignment was run on. (Same reason `ER_TIER_FLAGS` exists: the tier was hardcoded at four sites and
/// two of them could describe a command the binary never ran.)
pub(crate) const SUBSTRATE_EXON_SUM: &str = "exon-sum";
pub(crate) const SUBSTRATE_GENOMIC_SPAN: &str = "genomic-span";

/// The name of the substrate an E_r tier aligned. `genomic == true` ⟹ the genomic span of the locus
/// (introns + flanks), else the spliced exon-sum.
pub(crate) fn substrate_name(genomic: bool) -> &'static str {
    if genomic {
        SUBSTRATE_GENOMIC_SPAN
    } else {
        SUBSTRATE_EXON_SUM
    }
}

/// ⭐ **THE STATE OF THE ADDITIVE GENOMIC-SPAN TIER AT A CALL SITE.** `family_refine` does not compute a
/// single-substrate `E_r`: it runs the primary tier on the CORE substrate and then, for any family the
/// core leaves unconnected, UNIONS IN the same tier re-run on the genomic span. Whether that second leg is
/// armed is a property of the *rule*, not of the data, and until O-3 it was decided inline by a nested
/// `if !params.include_introns { if let Some(g) = genome { … } }` that nothing else could read — so the
/// certificate printed `substrate = exon-sum` for an edge set that is `E_x ∪ E_g`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum GenomicTier {
    /// Armed: fires on every family the core tier leaves unconnected. The shipped refine default.
    Armed,
    /// The core tier ALREADY aligns the genomic span (`include_introns`), so there is no second substrate.
    CoreIsGenomic,
    /// Armed by configuration but no genome is reachable at this call site, so the leg cannot run.
    NoGenome,
    /// This call site has no additive genomic leg at all. `homology_edges_all_reps_pooled` SWAPS its
    /// substrate (`homology_genomic_span`) rather than unioning a second one in.
    NotPresent,
}

impl GenomicTier {
    /// The one predicate the gate branches on. Reading it here rather than re-deriving the condition is
    /// what makes the certificate and the run the same decision.
    pub(crate) fn armed(self) -> bool {
        matches!(self, GenomicTier::Armed)
    }
    /// The certificate value. Says what it is AND, when off, why — an unexplained `off` cannot be told
    /// apart from a missing genome, which is the difference between a rule and an accident.
    pub(crate) fn label(self) -> String {
        match self {
            GenomicTier::Armed => format!(
                "armed ({SUBSTRATE_GENOMIC_SPAN}, same primary tier; UNIONED into any family the core tier leaves unconnected)"
            ),
            GenomicTier::CoreIsGenomic => {
                format!("off (core substrate is already {SUBSTRATE_GENOMIC_SPAN})")
            }
            GenomicTier::NoGenome => "off (no genome reachable at this call site)".into(),
            GenomicTier::NotPresent => "absent (single-substrate site)".into(),
        }
    }
}

/// ⭐ **THE ONE DECISION.** `family_refine`'s gate and `er_rule_rows` BOTH call this, so `params.tsv`
/// cannot describe a tier the run did not arm. Do NOT re-derive the condition anywhere else — that
/// duplication is exactly how the audit trail came to describe a command the binary never ran (B1).
pub(crate) fn additive_genomic_tier(params: &RefineParams, genome_available: bool) -> GenomicTier {
    if params.include_introns {
        GenomicTier::CoreIsGenomic
    } else if !genome_available {
        GenomicTier::NoGenome
    } else {
        GenomicTier::Armed
    }
}

/// Which sequence each E_r call site aligns, and whether it can offer a read-supported core denominator.
/// Only these things legitimately differ between the two sites, so they are named rather than left
/// to be inferred from a key that happens to be spelled differently on each side.
///
/// ⚠ `include_introns` (refine) and `homology_genomic_span` (the catalog) are **NOT the same substrate
/// under two names** — that earlier comment was false and is corrected here (O-4). Each SWAPS its own
/// site's core substrate, but refine ADDITIONALLY unions a genomic-span tier in (`genomic_tier`) while
/// the catalog has no additive leg at all. A swap and a gated union are different operations, and that
/// difference is the one line `diff <prefix>.rule.tsv <prefix>.refine.rule.tsv` prints — see
/// `site_divergence_policy` below for why it is intentional rather than drift.
pub(crate) struct ErRuleSite {
    /// `true` ⟹ the primary tier aligns the GENOMIC span; `false` ⟹ the spliced exon-sum. This is the
    /// CORE substrate only; a site may additionally union a second one in — see `genomic_tier`.
    pub substrate_genomic: bool,
    /// `true` ⟹ per-rep read-supported core lengths were passed, so `RUSTLE_ER_CORE_COVERAGE` can bite.
    pub core_lens_supplied: bool,
    /// The additive genomic-span leg's state, taken from `additive_genomic_tier()` at the site — never
    /// re-derived here.
    pub genomic_tier: GenomicTier,
}

/// ⭐ **THE E_r RULE AS A FILE.** Every knob that DECIDES an edge, at its EFFECTIVE value, and **nothing
/// data-dependent** — no counts, no lengths, no paths. Emitted verbatim by BOTH E_r call sites as
/// `<prefix>.rule.tsv` (O1) and `<prefix>.refine.rule.tsv` (O2), so
///
///     diff <prefix>.rule.tsv <prefix>.refine.rule.tsv
///
/// is the answer to "do O1 and O2 use the same rule?". An EMPTY diff is the certificate; any line is the
/// drift. This is what row X.2 of `docs/OBJECTIVES_AND_VERIFICATION.md` recorded as missing: the previous
/// `params.tsv` was written only by `homology_edges_all_reps_pooled`, which O2's refine never calls, so
/// the "settle it by diffing two files" claim was an O1-only capability (25/25 written on the O1 side,
/// 0/25 on the O2 side).
///
/// ⚠ Keep this free of anything a run's INPUT can move. A count in here turns a rule diff into a data
/// diff and the certificate stops meaning anything.
pub(crate) fn er_rule_rows(params: &RefineParams, site: &ErRuleSite) -> Vec<(String, String)> {
    let (seed, floor, mask) = er_primary_tier(params);
    let longer = matches!(std::env::var("RUSTLE_ER_COVERAGE_LONGER"), Ok(v) if v != "0" && !v.is_empty());
    vec![
        ("mm_flags".into(), ER_TIER_FLAGS.join(" ")),
        ("primary_tier_seed".into(), seed.join(" ")),
        ("primary_tier_identity".into(), format!("{floor:.6}")),
        ("primary_tier_name".into(), tier_names(mask)),
        ("sensitive_only".into(), er_sensitive_only(params).to_string()),
        // The ADDITIVE sensitive run. Under sensitive-only on the exon-sum it is the primary tier itself,
        // so re-running it would be a byte-identical no-op; that is stated here rather than left implicit.
        ("additive_sensitive_tier".into(), {
            if !params.nucleotide_sensitive {
                "off".into()
            } else if er_sensitive_only(params) && !site.substrate_genomic {
                "same-as-primary (not re-run)".into()
            } else {
                format!("{} @ {:.6}", ER_SENSITIVE_SEED.join(" "), params.sensitive_identity)
            }
        }),
        ("protein_tier".into(), params.protein_tail.to_string()),
        ("min_coverage".into(), format!("{:.6}", params.min_coverage)),
        // The SECOND, longer-side coverage floor. A RULE row: it decides edges, so an ON arm must
        // be distinguishable from an OFF arm by the certificate alone (defect M2).
        ("min_coverage_longer".into(), match er_cov_longer_floor() {
            Some(f) => format!("{f:.6}"),
            None => "<unset>".into(),
        }),
        // ⚠ NAMED `core_substrate`, NOT `substrate`. A run whose edge set is `E_x ∪ E_g` has no single
        // substrate, and the old key invited exactly the misreading it produced: `substrate = exon-sum`
        // printed on a run that unioned genomic-span edges in (O-3 / joint-run finding F6).
        ("core_substrate".into(), substrate_name(site.substrate_genomic).into()),
        // The ADDITIVE GENOMIC-SPAN leg. This is a RULE row: whether the leg is armed decides edges. How
        // often it FIRED and how many edges it contributed are data, and live in `params.tsv`.
        ("additive_genomic_tier".into(), site.genomic_tier.label()),
        // ⭐ O-4 — WHY THE ONE DIFFERING LINE ABOVE IS A DECISION AND NOT DRIFT. Emitted IDENTICALLY at
        // both sites (a constant), so it never widens the diff; its whole job is that a reader who diffs
        // the two certificates and finds `additive_genomic_tier` differing learns, from the file itself,
        // that the divergence was measured and kept. Constant ⟹ it cannot turn a rule diff into a data
        // diff. Evidence and denominators: `bench/FALSE_NEGATIVES.md`.
        (
            "site_divergence_policy".into(),
            "additive_genomic_tier is the ONLY axis on which the O1 and O2 E_r sites differ; \
             O-4 2026-08-13 ADJUDICATED KEEP BOTH (refine's leg is GATED+family-local: emitted \
             catalogs identical 0/26 vs DNA-only, ARI 1.0000; the O1 site's swap is global and its \
             precision sign is truth-dependent) -- see bench/FALSE_NEGATIVES.md"
                .into(),
        ),
        ("coverage_denominator".into(), match (core_cov_floor(), site.core_lens_supplied) {
            (Some(f), true) => format!("core(min)@floor={f:.6} else span"),
            // The core denominator is only reachable where core lengths were passed. Saying "core" at a
            // site that never supplies them would describe a rule the binary cannot run.
            (Some(f), false) => format!("span(min) [RUSTLE_ER_CORE_COVERAGE={f:.6} UNREACHABLE here: no core lens]"),
            (None, _) if longer => "span(max)".into(),
            (None, _) => "span(min)".into(),
        }),
        ("coverage_form".into(), ER_COVERAGE_FORM.into()),
        ("identity_metric".into(), "1-de (fallback nmatch/blocklen when de:f: absent)".into()),
        (
            "alignment_orientation".into(),
            if params.forward_only_active() {
                "forward-only (+); valid only for transcript-oriented RNA representatives".into()
            } else {
                "both (+/-)".into()
            },
        ),
        ("edge_rule".into(), "ANY single record clearing both floors".into()),
        ("summed_coverage".into(), std::env::var("RUSTLE_ER_SUM_COVERAGE").unwrap_or_else(|_| "<unset>".into())),
        ("drop_stub_edges".into(), er_no_stub_edges().to_string()),
        ("shared_exon_mode".into(), std::env::var("RUSTLE_SHARED_EXON").unwrap_or_else(|_| "<unset>".into())),
        ("repeat_hub_gate".into(), std::env::var("RUSTLE_ER_REPEAT_GATE").unwrap_or_else(|_| "<unset>".into())),
        ("junction_majority".into(), std::env::var("RUSTLE_JUNCTION_MAJORITY").unwrap_or_else(|_| "<unset>".into())),
        // Changes WHICH SKELETONS BECOME NODES, so an ON and an OFF catalog must not have byte-identical
        // params.tsv (the M2 defect).
        ("gate_min_reads".into(), super::denovo_assemble::gate_min_reads().to_string()),
        // Changes WHICH SITES ARE PROPOSED (every placement, not just the primary-flagged one), so an ON
        // and an OFF catalog must not have byte-identical params.tsv (the M2 defect).
        ("flagfree_sites".into(), std::env::var("RUSTLE_FLAGFREE_SITES").unwrap_or_else(|_| "<unset>".into())),
        ("junction_nc_max_bp".into(), std::env::var("RUSTLE_JUNCTION_NC_MAX_BP").unwrap_or_else(|_| "10000 (default)".into())),
        ("footprint_nodes".into(), std::env::var("RUSTLE_FOOTPRINT_NODES").unwrap_or_else(|_| "<unset>".into())),
        ("footprint_min_cov".into(), std::env::var("RUSTLE_FOOTPRINT_MIN_COV").unwrap_or_else(|_| "2 (default)".into())),
        ("footprint_windows".into(), std::env::var("RUSTLE_FOOTPRINT_WINDOWS").unwrap_or_else(|_| "<unset>".into())),
        ("weighted_partition".into(), std::env::var("RUSTLE_ER_WEIGHTED_PARTITION").unwrap_or_else(|_| "<unset>".into())),
        ("tier2_admit".into(), std::env::var("RUSTLE_TIER2_ADMIT").unwrap_or_else(|_| "<unset>".into())),
        ("collapse_exonic".into(), std::env::var("RUSTLE_COLLAPSE_EXONIC").unwrap_or_else(|_| "<unset>".into())),
        // Changes WHICH TRANSCRIPTS COLLAPSE INTO ONE LOCUS -- i.e. what a NODE IS -- so an ON and an OFF
        // catalog must not have byte-identical params.tsv (the M2 defect).
        ("collapse_unstranded".into(), std::env::var("RUSTLE_COLLAPSE_UNSTRANDED").unwrap_or_else(|_| "<unset>".into())),
        ("minimap2".into(), params.minimap2.clone()),
    ]
}

/// Write a `key\tvalue` TSV. Used for both `.rule.tsv` and `.params.tsv` so the two can never disagree
/// about their own format.
fn write_kv_tsv(path: &str, rows: &[(String, String)]) {
    use std::io::Write;
    match std::fs::File::create(path) {
        Ok(mut fh) => {
            let _ = writeln!(fh, "key\tvalue");
            for (k, v) in rows {
                let _ = writeln!(fh, "{k}\t{v}");
            }
            eprintln!("[er-dump] -> {path}");
        }
        Err(e) => eprintln!("[er-dump] FAILED to write {path}: {e}"),
    }
}

/// The same invocation as a string, for logs / `.args` dumps / `params.tsv`. Derived from
/// `er_tier_command`'s constants so a dump can never describe a command that was not run.
pub(crate) fn er_tier_cmdline(minimap2: &str, threads: usize, mm_args: &[&str]) -> String {
    let mut s = format!("{minimap2} {} -t {}", ER_TIER_FLAGS.join(" "), threads.max(1));
    if !mm_args.is_empty() {
        s.push(' ');
        s.push_str(&mm_args.join(" "));
    }
    s
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
    cores: Option<&[u64]>,
    params: &RefineParams,
) -> Result<Vec<(usize, usize)>> {
    nucleotide_edges_tagged(seqs, mm_args, min_id, min_cov, cores, params, None)
}

/// As `nucleotide_edges`, but LABELS the PAF/`.args` this call drops under `RUSTLE_ER_EDGE_DUMP`.
///
/// ⚠ WHY THIS EXISTS (O-3). `family_refine` runs the SAME tier twice per unconnected family — once on the
/// exon-sum core, once on the genomic span — and the two dumps differed only by a call counter, so the
/// audit trail could not say which substrate produced which alignment. The tag is built from
/// `substrate_name()`, the single spelling the certificate also uses.
pub(crate) fn nucleotide_edges_tagged(
    seqs: &[Vec<u8>],
    mm_args: &[&str],
    min_id: f64,
    min_cov: f64,
    cores: Option<&[u64]>,
    params: &RefineParams,
    dump_tag: Option<&str>,
) -> Result<Vec<(usize, usize)>> {
    Ok(nucleotide_edges_scored(seqs, mm_args, min_id, min_cov, cores, params, dump_tag, None)?
        .into_iter()
        .map(|(i, j, _ident, _cov)| (i, j))
        .collect())
}

/// ADDITIVE DISCLOSURE for one E_r edge (`docs/o1_ledger.md` §6ay/§6az), read off the SAME exemplar PAF
/// record that supplied the row's `identity`/`coverage`.
///
/// WHY THIS EXISTS. E_r charges coverage on the **SHORTER** sequence only, so `1 - coverage` is the
/// shorter member's overhang and is bounded at the 0.50 floor. **The LONGER member's overhang is measured
/// by nothing.** Measured on the shipped catalog: median `cov_longer` 0.44, and 1,035/1,726 = 59.97% of
/// directly-certified within-family pairs have `cov_longer < 0.50`; a synthetic pair built to differ by
/// 1,000 bp at EACH end reports `coverage` = 1.000. So a perfect-looking coverage number is compatible
/// with a 2 kb end difference, and nothing in the dump said so.
///
/// ⛔ THIS IS DISCLOSURE, NOT A RULE. `RUSTLE_ER_COVERAGE_LONGER` REPLACES the denominator and is
/// therefore a GATE — flipping it would fail the certifying record of up to 60% of within-family edges.
/// These three numbers gate NOTHING: they are computed, stored in a map parallel to `metrics`, and
/// printed. The edge SET and every pre-existing column are untouched by construction.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(crate) struct ErEdgeFlank {
    /// Aligned span measured on the LONGER sequence's axis, divided by the LONGER length. The mirror of
    /// the shipped `coverage`, which is the same quantity on the SHORTER axis.
    pub(crate) cov_longer: f64,
    /// Total unaligned flank (5' + 3') of rep `i` on the exemplar record, in bp.
    pub(crate) unaln_i: u64,
    /// Total unaligned flank (5' + 3') of rep `j` on the exemplar record, in bp.
    pub(crate) unaln_j: u64,
}

/// The disclosure arithmetic, isolated so it is testable without an aligner.
///
/// `q_is_i` says whether the PAF's QUERY is the lower-indexed rep of the pair (`i = q.min(t)`), which is
/// what maps `qs/qe/ql` and `ts/te/tl` onto `unaln_i`/`unaln_j`. Under the shipped `-X` (`--dual=no`)
/// exactly one orientation per pair is emitted and minimap2 does not guarantee which sequence is the
/// query, so this mapping cannot be assumed.
///
/// The longer side is chosen with `ql >= tl`, the same tie convention `RUSTLE_ER_COVERAGE_LONGER` uses
/// (`ql.max(tl)` with the axis following the denominator). On a tie the two denominators are equal, so
/// only the numerator's axis is at stake.
///
/// ⚠ `cov_longer <= coverage` is a THEOREM only when the aligned span is the same on both axes (a gapless
/// or indel-balanced record): the numerator is then shared and the denominator is larger. An indel-heavy
/// record can in principle break it, which is itself worth seeing rather than hiding.
pub(crate) fn er_edge_flank(
    ql: f64,
    qs: f64,
    qe: f64,
    tl: f64,
    ts: f64,
    te: f64,
    q_is_i: bool,
) -> ErEdgeFlank {
    let (long_aln, long_len) = if ql >= tl { (qe - qs, ql) } else { (te - ts, tl) };
    let cov_longer = long_aln / long_len.max(1.0);
    let q_unaln = (qs + (ql - qe)).max(0.0).round() as u64;
    let t_unaln = (ts + (tl - te)).max(0.0).round() as u64;
    let (unaln_i, unaln_j) = if q_is_i { (q_unaln, t_unaln) } else { (t_unaln, q_unaln) };
    ErEdgeFlank { cov_longer, unaln_i, unaln_j }
}

/// As `nucleotide_edges`, but also returns the EXEMPLAR identity and coverage of each edge — the passing
/// record with the highest coverage (ties broken by higher identity), which is deterministic.
///
/// WHY THIS SPLIT EXISTS. `homology_edges_all_reps_pooled` receives only `(usize, usize)` pairs: the two
/// numbers that actually decide an edge are computed here and were discarded one stack frame later, so the
/// Rust engine could report THAT an edge exists but never WHY. That made true parity against the canonical
/// Python mirror (`bench/soto/rustlib.py`) unverifiable — a diff could show a set difference but not
/// attribute it to the identity metric (`1 - de` vs `nmatch/blocklen`), the coverage denominator
/// (`min` vs `max` vs read-supported core), the floors, or the tier. `RUSTLE_ER_EDGE_DUMP` consumes these
/// values; see the dump block at the end of `homology_edges_all_reps_pooled`.
///
/// The exemplar is REPORTING ONLY. The edge SET is unchanged: a pair is still an edge as soon as ONE
/// record clears both floors, and coverage is still never summed across records on the default path.
///
/// `dump_tag` labels the PAF this call writes under `RUSTLE_ER_EDGE_DUMP`. `Some("er")` marks the two tier
/// calls that actually BUILD the returned E_r edge set; `None` (every other call site) marks the downstream
/// per-family refinement runs. Without the tag a single run drops several same-named PAFs in the directory
/// and a differ can silently pick one that never contributed to the dumped edge set — which would produce a
/// confident, wrong parity result, the exact failure this dump exists to prevent.
fn nucleotide_edges_scored(
    seqs: &[Vec<u8>],
    mm_args: &[&str],
    min_id: f64,
    min_cov: f64,
    cores: Option<&[u64]>,
    params: &RefineParams,
    dump_tag: Option<&str>,
    guard_exempt: Option<&dyn Fn(usize, usize) -> bool>,
) -> Result<Vec<(usize, usize, f64, f64)>> {
    nucleotide_edges_scored_disclosed(
        seqs, mm_args, min_id, min_cov, cores, params, dump_tag, guard_exempt, None,
    )
}

/// As `nucleotide_edges_scored`, but ALSO fills `disclose` with the per-edge `ErEdgeFlank` read off the
/// very same exemplar record that supplied the returned `(identity, coverage)`.
///
/// A SEPARATE map, deliberately: nothing that participates in the edge DECISION is widened or moved, so
/// the returned edge set is bit-for-bit what it was. Entries appear only for edges whose exemplar came
/// from a real PAF record — an edge admitted by the opt-in summed-coverage rule has no single certifying
/// record and therefore no entry, which the dump prints as `NA` rather than inventing.
#[allow(clippy::too_many_arguments)]
fn nucleotide_edges_scored_disclosed(
    seqs: &[Vec<u8>],
    mm_args: &[&str],
    min_id: f64,
    min_cov: f64,
    cores: Option<&[u64]>,
    params: &RefineParams,
    dump_tag: Option<&str>,
    // Pair-level exemption from the orientation guard, supplied by the caller because only it knows what
    // the sequences ARE. Passed as a predicate rather than as parallel arrays so no per-rep data has to be
    // threaded through this function's seven other call sites. `None` = the guard applies unconditionally,
    // which is the shipped behaviour.
    guard_exempt: Option<&dyn Fn(usize, usize) -> bool>,
    // Out-parameter, reporting only. `None` = not collected.
    disclose: Option<&mut BTreeMap<(usize, usize), ErEdgeFlank>>,
) -> Result<Vec<(usize, usize, f64, f64)>> {
    use std::io::Write;
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    let nonce = seqs.iter().map(|s| s.len()).sum::<usize>().wrapping_mul(1000003)
        ^ seqs.len()
        ^ mm_args.len().wrapping_mul(7);
    // The nonce is a pure function of (total length, count, arg count), so two CONCURRENT calls of the same
    // shape in one process resolve to the same path: both write it and the first `Cleanup` to drop deletes
    // it under the other. Production is sequential here, but `cargo test` runs every test in threads of a
    // single process, which made the suite flaky (2 of 3 runs). A per-call counter makes the path unique
    // without losing the nonce's diagnostic value.
    static TMP_SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let uniq = TMP_SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = dir.join(format!("rustle_refine_{pid}_{nonce}_{uniq}.fa"));
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
    // ⚠⚠ STREAMED, NOT BUFFERED (2026-08-19). This used `.output()`, which holds the ENTIRE
    // all-vs-all PAF in a `Vec<u8>` before parsing a single line. Per-family that is nothing; genome-wide
    // it is the whole 12,415-rep all-vs-all through the sensitive `-k11 -w5` tier, and it is what made the
    // genome-wide catalog unrunnable — a run reached 23.7 GB RSS + 10.2 GB swap, `D` state, ~11% CPU
    // (the parent blocked accumulating the child's stdout) and had to be killed at 1h34m.
    // Memory is now bounded by one line regardless of PAF size. Records, order and edges are unchanged.
    let mut child = er_tier_command(&params.minimap2, params.threads, mm_args)
        .arg(&path)
        .arg(&path)
        .stdout(std::process::Stdio::piped())
        .spawn()
        .map_err(|e| anyhow::anyhow!("failed to run minimap2 ('{}') for refinement: {e}", params.minimap2))?;
    let child_stdout = child
        .stdout
        .take()
        .ok_or_else(|| anyhow::anyhow!("minimap2 stdout pipe unavailable"))?;
    // Opt-in: hand the parity differ the EXACT bytes this rule was applied to. Without this the Python
    // mirror has to re-run minimap2 itself, and any difference in that invocation (the on-disk fixture PAF
    // was made with `-x asm20 ... -N 200 -p 0.02` and NO `-X`, i.e. both orientations of every pair) is
    // charged to the edge rule instead of to the alignment. The rep FASTA is copied too because the
    // `Cleanup` guard above deletes it, and its `>{i}` headers ARE the node ids used in the PAF.
    // ⚠ On a genome-wide run this file is the full all-vs-all PAF and can be large.
    // The dump is now TEED while streaming rather than written from a buffer.
    let mut dump_stem: Option<(String, String)> = None;
    let mut dump_w: Option<std::io::BufWriter<std::fs::File>> = None;
    if let Ok(prefix) = std::env::var("RUSTLE_ER_EDGE_DUMP") {
        if !prefix.is_empty() {
            let argsig: String = mm_args
                .join("")
                .chars()
                .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
                .collect();
            let tag = dump_tag.map(|t| format!("{t}.")).unwrap_or_else(|| "refine.".to_string());
            let stem = format!("{prefix}.{tag}{argsig}.{uniq}");
            if let Ok(f) = std::fs::File::create(format!("{stem}.paf")) {
                dump_w = Some(std::io::BufWriter::new(f));
                dump_stem = Some((stem, tag));
            }
        }
    }
    let mut n_records = 0usize;
    // OPTIONAL summed-coverage rule (`RUSTLE_ER_SUM_COVERAGE=1`, default off = byte-identical).
    //
    // The default rule evaluates coverage on ONE record and never sums, so two loci sharing 60% of the
    // shorter sequence across three separate blocks get no edge. That is exactly what a SHATTERED locus
    // representative produces: reads in segmental duplications fragment into many exact intron chains
    // (GTF2IP14: 287 reads over 222 distinct chains, the winner holding 13), so the rep describes a
    // fragment and fails the 0.50 floor against its own paralogs. Measured: 32 of the 61 Soto members that
    // no mode finds have a >=3-read chain available, i.e. the locus is buildable and is lost here.
    //
    // Summing raises the NUMERATOR and leaves the sequences untouched, unlike the exon-union substrate,
    // which lengthened them and inflated the denominator instead (that cost 20 recall points).
    //
    // Guards, because summing is what the per-record rule exists to prevent:
    //   - only records on the SAME strand for that pair are summed (a real fragmented gene is collinear;
    //     a repeat matches in both orientations),
    //   - only records whose query span is >= `min_block` bp count, so a swarm of short repeat hits
    //     cannot accumulate past the floor,
    //   - query intervals are UNIONed, never added, so overlapping records cannot double-count.
    let sum_cov = std::env::var("RUSTLE_ER_SUM_COVERAGE")
        .map(|v| v != "0" && !v.is_empty())
        .unwrap_or(false);
    let min_block: f64 = std::env::var("RUSTLE_ER_SUM_MIN_BLOCK")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(200.0);
    let mut sum_blocks: BTreeMap<(usize, usize, char), (Vec<(u64, u64)>, f64, f64)> = BTreeMap::new();

    // Value = the EXEMPLAR (identity, coverage) for reporting only; membership is decided exactly as before.
    let mut edge_set: BTreeMap<(usize, usize), (f64, f64)> = BTreeMap::new();
    // Parallel to `edge_set` and updated in LOCKSTEP with it, so the disclosed flank numbers always
    // describe the record that supplied this edge's exemplar identity/coverage. Reporting only.
    let mut flanks: BTreeMap<(usize, usize), ErEdgeFlank> = BTreeMap::new();
    use std::io::BufRead as _;
    for line in std::io::BufReader::new(child_stdout).lines() {
        let line = line.map_err(|e| anyhow::anyhow!("reading minimap2 stdout: {e}"))?;
        n_records += 1;
        if let Some(w) = dump_w.as_mut() {
            let _ = writeln!(w, "{line}");
        }
        let line = line.as_str();
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let (q, t) = (f[0].parse::<usize>().ok(), f[5].parse::<usize>().ok());
        let (q, t) = match (q, t) {
            (Some(q), Some(t)) if q != t => (q, t),
            _ => continue,
        };
        let strand = f[4].chars().next().unwrap_or('+');
        if params.forward_only_active() && strand != '+' && !guard_exempt.is_some_and(|ex| ex(q, t)) {
            // A MINUS-strand record means the two stored sequences disagree in orientation. That is
            // evidence of antisense only when BOTH strands were actually measured: a rep whose strand was
            // never junction-determined carries the `'+'` placeholder, so it may simply be stored
            // backwards, and its minus alignment says nothing about sense. Measured on the shipped
            // catalog: 3,951/4,009 = 0.9855 of guard-blocked pairs involve such a rep, leaving 58 genuine
            // antisense candidates genome-wide. See `guard_exempt` at the call site for the scoping.
            continue;
        }
        let ql = f[1].parse::<f64>().unwrap_or(0.0);
        let qs = f[2].parse::<f64>().unwrap_or(0.0);
        let qe = f[3].parse::<f64>().unwrap_or(0.0);
        let tl = f[6].parse::<f64>().unwrap_or(0.0);
        // M1: the TARGET-axis aligned span. Previously unparsed, which is why the numerator was always
        // taken from the query axis even when the denominator came from the target (see the coverage
        // block below). PAF target start/end are columns 8/9 (0-based fields 7/8) and are always in
        // forward-target orientation, so `te - ts` is a length regardless of the strand column.
        let ts = f[7].parse::<f64>().unwrap_or(0.0);
        let te = f[8].parse::<f64>().unwrap_or(0.0);
        let de = f[12..].iter().find_map(|x| x.strip_prefix("de:f:").and_then(|v| v.parse::<f64>().ok()));
        let ident = match de {
            Some(d) => 1.0 - d,
            None => {
                let nmatch = f[9].parse::<f64>().unwrap_or(0.0);
                let alnlen = f[10].parse::<f64>().unwrap_or(1.0).max(1.0);
                nmatch / alnlen
            }
        };
        // Denominator: the shorter SPAN by default, or the shorter READ-SUPPORTED CORE when
        // `RUSTLE_ER_CORE_COVERAGE` is set and both cores were measured. The span is produced by the
        // pipeline itself, so a boundary error moves the criterion that decides membership; the core is
        // depth-defined and a readthrough tail never enters it. Falls back to the span whenever either core
        // is 0 (not measured), so the rule can never become MORE permissive through missing data.
        // The floor must travel WITH the denominator. The core is a strictly smaller target than the span
        // (measured core/span median 0.41, quartiles 0.10-0.68), so the same fraction means a different
        // demand; applying the span's 0.50 to a core denominator silently makes the rule LOOSER.
        // `RUSTLE_ER_COVERAGE_LONGER=1` divides by the LONGER sequence instead of the shorter, i.e. it
        // demands that BOTH members be covered (equivalently BLAST's qcovs AND scovs clearing the floor).
        // Default off = byte-identical.
        //
        // WHY. Coverage-of-shorter is STRUCTURALLY BLIND to truncation: a 10% fragment that aligns fully
        // into a complete sibling scores 1.00, so the rule cannot tell a copy from a piece of one. Three
        // independent measurements point at the same fix:
        //   - as a certificate on RNA pairs it gives 129 TRUE vs 2 FALSE, precision 0.985;
        //   - it is the mechanism behind the only measured RNA-vs-DNA partition contradiction — a 2,037 bp
        //     NPIPB6 fragment reaches coverage 0.948 against a 38,653 bp chimeric read-through node while
        //     touching 5% of it, dragging EIF3CL into the NPIP family;
        //   - it is what would make component merging MONOTONE under improving evidence: with min(), a rep
        //     that gets MORE complete can BREAK its edge to a shorter sibling (the exon-union substrate lost
        //     20 recall points to exactly that), so "better reads can only merge components" is currently
        //     false, and that monotonicity is what makes a witnessed component a coarser reading of a family
        //     rather than a competing answer.
        //
        // ⚠ It is NOT free: the duplicated unit is size-invariant while annotated spans are not (NPIP genes
        // span 10.6-49.4 kb around a ~16 kb cassette, block length correlating with max(span) at only
        // r=+0.196), so dividing by the longer span penalises pairs whose ANNOTATION differs rather than
        // whose SEQUENCE does. Measured ceiling on NPIP: only 134/171 true pairs could reach 0.50 at all,
        // and NPIPB8-NPIPB2 caps at 0.215. Expect a recall cost; the question this knob exists to answer is
        // whether the precision gain pays for it end to end.
        //
        // ⚠⚠ M1 — THE COVERAGE STATISTIC MUST BE A FRACTION, AND BEFORE THIS FIX IT WAS NOT.
        //
        // The numerator was ALWAYS the query-axis aligned span `qe - qs` while the denominator was
        // `min(ql, tl)` — i.e. a query-axis numerator over a possibly-TARGET denominator. That is not an
        // aligned FRACTION of anything; it is two different sequences' measurements divided.
        //
        // It bites because the shipped tier passes `-X`, which implies `--dual=no`: exactly ONE
        // orientation per pair is emitted, and minimap2 does not guarantee the query is the shorter
        // sequence. Measured on the shipped catalog PAFs the query is the LONGER sequence in ~60% of
        // records, and 110 of 939 accepted edges (11.7%) scored coverage ABOVE 1.0, maximum 1.178 — a
        // fraction exceeding 1 is the defect making itself visible.
        //
        // THE FIX: the numerator's AXIS FOLLOWS THE DENOMINATOR. Whichever sequence supplies the
        // denominator also supplies the aligned span:
        //     default (shorter):  ql <= tl ? (qe-qs)/ql : (te-ts)/tl
        //     COVERAGE_LONGER:    ql >= tl ? (qe-qs)/ql : (te-ts)/tl
        //     CORE_COVERAGE:      the side whose CORE was selected supplies the span
        // This is NOT a null change — it removes edges whose apparent coverage was borrowed from the
        // other sequence's axis, and it can add edges where the target side is the better-covered one.
        let longer_cov = matches!(std::env::var("RUSTLE_ER_COVERAGE_LONGER"), Ok(v) if v != "0" && !v.is_empty());
        // `true` = the QUERY side supplies both numerator and denominator; `false` = the TARGET side.
        // Note this choice is orientation-independent: if a record arrived with q and t swapped, the
        // flag flips with it and the SAME physical sequence is still selected.
        let mut side_is_query = if longer_cov { ql >= tl } else { ql <= tl };
        let span_denom = if longer_cov { ql.max(tl).max(1.0) } else { ql.min(tl).max(1.0) };
        let (shorter, floor) = match core_cov_floor() {
            Some(f) => match (cores.and_then(|c| c.get(q)), cores.and_then(|c| c.get(t))) {
                (Some(&a), Some(&b)) if a > 0 && b > 0 => {
                    let d = if longer_cov { a.max(b) } else { a.min(b) };
                    // The core denominator can select the OPPOSITE side from the span, so recompute the
                    // axis against the cores rather than inheriting the span's answer.
                    side_is_query = if longer_cov { a >= b } else { a <= b };
                    ((d as f64).max(1.0), f)
                }
                // Core unmeasured for this pair: fall back to the span AND to the span's floor, so missing
                // data can never make the rule more permissive.
                _ => (span_denom, min_cov),
            },
            None => (span_denom, min_cov),
        };
        let aln_on_denom_axis = if side_is_query { qe - qs } else { te - ts };
        let cov = aln_on_denom_axis / shorter;
        // SECOND floor on the LONGER sequence, additive to the clause above (BLAST qcovs AND scovs).
        // Uses `er_edge_flank`'s formula verbatim so the gated quantity is the SAME number the dump
        // already discloses as `cov_longer` -- a reader can check the gate against the column.
        let longer_ok = match er_cov_longer_floor() {
            Some(f) => er_edge_flank(ql, qs, qe, tl, ts, te, true).cov_longer >= f,
            None => true,
        };
        if ident >= min_id && cov >= floor && longer_ok {
            let k = (q.min(t), q.max(t));
            // Keep the highest-coverage passing record as the exemplar (ties -> higher identity). Pure
            // reporting: the KEY set is identical either way, so this cannot move an edge.
            //
            // Spelled with `Entry` rather than `or_insert` + compare ONLY so the disclosure map can learn
            // whether this record actually became the exemplar. The two branches are exactly the old
            // behaviour: a vacant key takes this record; an occupied one is replaced under the identical
            // `(cov, ident) >` test (which a fresh insert could never pass against itself).
            let took_exemplar = match edge_set.entry(k) {
                std::collections::btree_map::Entry::Vacant(v) => {
                    v.insert((ident, cov));
                    true
                }
                std::collections::btree_map::Entry::Occupied(mut o) => {
                    if (cov, ident) > (o.get().1, o.get().0) {
                        o.insert((ident, cov));
                        true
                    } else {
                        false
                    }
                }
            };
            if took_exemplar {
                // §6az DISCLOSURE, from THIS record — the same one the row's identity/coverage came from.
                flanks.insert(k, er_edge_flank(ql, qs, qe, tl, ts, te, q < t));
            }
        }
        // M1 applies here too: the summed rule unions intervals into the SAME denominator, so the
        // intervals must be measured on the denominator's axis or the union is a sum of two coordinate
        // systems. `min_block` is likewise a length on that axis.
        let (aln_s, aln_e) = if side_is_query { (qs, qe) } else { (ts, te) };
        if sum_cov && ident >= min_id && aln_on_denom_axis >= min_block {
            let entry = sum_blocks
                .entry((q.min(t), q.max(t), strand))
                .or_insert_with(|| (Vec::new(), shorter, ident));
            entry.0.push((aln_s as u64, aln_e as u64));
            entry.1 = entry.1.min(shorter);
            entry.2 = entry.2.max(ident);
        }
    }
    // ---- the child is done: status check + dump finalisation ----
    // ⚠ SEMANTICS PRESERVED: a non-zero minimap2 exit is still SILENTLY an empty edge set (the family
    // dissolves). Streaming means we discover that AFTER parsing, so whatever was accumulated is
    // discarded here rather than returned.
    let status = child
        .wait()
        .map_err(|e| anyhow::anyhow!("waiting for minimap2 ('{}'): {e}", params.minimap2))?;
    if let Some(mut w) = dump_w.take() {
        let _ = w.flush();
    }
    if !status.success() {
        if std::env::var("RUSTLE_ER_EDGE_DUMP").is_ok() {
            eprintln!(
                "[er-dump] ⚠ minimap2 EXITED NON-ZERO ({}) for args {:?} on {} seqs -> ZERO edges. \
                 An empty dump here is a FAILED ALIGNMENT, not a result.",
                status,
                mm_args,
                seqs.len()
            );
        }
        return Ok(Vec::new());
    }
    if let Some((stem, tag)) = dump_stem {
        let _ = std::fs::copy(&path, format!("{stem}.reps.fa"));
        let _ = std::fs::write(
            format!("{stem}.args"),
            format!(
                // ⚠ `tier` is the O-3 fix: the exon-sum core run and the additive genomic-span run use
                // the SAME command on DIFFERENT sequence, so without this line two `.args` files from one
                // family are byte-identical and neither says what it aligned.
                "{} <reps.fa> <reps.fa>\ntier\t{}\nmin_identity\t{}\nmin_coverage\t{}\ncoverage_form\t{}\n",
                er_tier_cmdline(&params.minimap2, params.threads, mm_args),
                tag.trim_end_matches('.'),
                min_id,
                min_cov,
                ER_COVERAGE_FORM
            ),
        );
        eprintln!(
            "[er-dump] tier {:?}: {n_records} PAF records -> {stem}.paf (+ .reps.fa, .args)",
            mm_args
        );
    }
    if sum_cov {
        let mut added = 0usize;
        for ((a, b, _strand), (mut iv, shorter, best_ident)) in sum_blocks {
            if edge_set.contains_key(&(a, b)) {
                continue;
            }
            iv.sort_unstable();
            let mut union_len = 0u64;
            let mut cur: Option<(u64, u64)> = None;
            for (s0, e0) in iv {
                match cur {
                    Some((cs, ce)) if s0 <= ce => cur = Some((cs, ce.max(e0))),
                    Some((cs, ce)) => {
                        union_len += ce - cs;
                        cur = Some((s0, e0));
                    }
                    None => cur = Some((s0, e0)),
                }
            }
            if let Some((cs, ce)) = cur {
                union_len += ce - cs;
            }
            let union_cov = union_len as f64 / shorter.max(1.0);
            if union_cov >= min_cov {
                edge_set.insert((a, b), (best_ident, union_cov));
                added += 1;
            }
        }
        if added > 0 {
            eprintln!("[summed-coverage] {added} additional edge(s) from collinear blocks >= {min_block} bp");
        }
    }
    if let Some(out) = disclose {
        out.clear();
        // Summed-coverage edges (opt-in) are inserted into `edge_set` with no certifying record, so they
        // are simply absent here; the filter is the assertion that this map never outruns the edge set.
        out.extend(flanks.iter().filter(|(k, _)| edge_set.contains_key(k)).map(|(k, v)| (*k, *v)));
    }
    Ok(edge_set.into_iter().map(|((i, j), (id, cov))| (i, j, id, cov)).collect())
}

/// Run the standard nucleotide E_r tier(s) on one explicit sequence substrate. This is the small shared
/// primitive used by the joint RNA/DNA certificate below: both arms inherit the same seeds, identity and
/// coverage thresholds, while the caller deliberately controls the substrate and orientation semantics.
fn nucleotide_rule_edge_set(
    seqs: &[Vec<u8>],
    cores: Option<&[u64]>,
    params: &RefineParams,
    dump_tag: &str,
) -> Result<BTreeSet<(usize, usize)>> {
    let (seed, floor, _) = er_primary_tier(params);
    let seed_ref: Vec<&str> = seed.iter().map(String::as_str).collect();
    let mut set: BTreeSet<(usize, usize)> = nucleotide_edges_scored(
        seqs,
        &seed_ref,
        floor,
        params.min_coverage,
        cores,
        params,
        Some(dump_tag), None)?
    .into_iter()
    .map(|(i, j, _, _)| (i, j))
    .collect();

    // When sensitive-only is active, `er_primary_tier` has already run this exact tier. Otherwise it is
    // the additive sensitive tier, matching `homology_edges_all_reps_pooled` without duplicating a run.
    if params.nucleotide_sensitive && !er_sensitive_only(params) {
        set.extend(
            nucleotide_edges_scored(
                seqs,
                ER_SENSITIVE_SEED,
                params.sensitive_identity,
                params.min_coverage,
                cores,
                params,
                Some(dump_tag), None)?
            .into_iter()
            .map(|(i, j, _, _)| (i, j)),
        );
    }
    Ok(set)
}

fn joint_component_count(n: usize, edges: &[(usize, usize)]) -> usize {
    if n == 0 {
        return 0;
    }
    let mut adj = vec![Vec::new(); n];
    for &(a, b) in edges {
        if a < n && b < n && a != b {
            adj[a].push(b);
            adj[b].push(a);
        }
    }
    let mut seen = vec![false; n];
    let mut components = 0;
    for root in 0..n {
        if seen[root] {
            continue;
        }
        components += 1;
        seen[root] = true;
        let mut stack = vec![root];
        while let Some(a) = stack.pop() {
            for &b in &adj[a] {
                if !seen[b] {
                    seen[b] = true;
                    stack.push(b);
                }
            }
        }
    }
    components
}

/// Write a typed RNA/DNA evidence certificate over the RNA-detected locus universe.
///
/// This intentionally does NOT union or intersect the two graphs to alter family membership. The RNA arm
/// aligns spliced, transcript-oriented exon sums; the DNA arm aligns the same loci's complete genomic spans
/// (transcript-normalised for a stable node representation) and remains orientation-agnostic so inverted
/// structural duplications are not erased. Their union is reported edge-by-edge as `RNA_DNA`, `RNA_ONLY`,
/// or `DNA_ONLY`, including cross-family edges that can expose a possible split or a repeat bridge.
///
/// Outputs:
///   `<out>.joint_edges.tsv`    one row per edge in either arm;
///   `<out>.joint_families.tsv` within-family connectivity and concordance summaries;
///   `<out>.joint_rule.tsv`     the typed, non-membership semantics of the comparison.
pub fn write_joint_rna_dna_certificate(
    out: &str,
    fams: &[Vec<DenovoTranscript>],
    fasta_path: &str,
    params: &RefineParams,
) -> Result<()> {
    use std::io::Write;

    anyhow::ensure!(
        !matches!(std::env::var("RUSTLE_SHARED_EXON"), Ok(v) if v != "0" && !v.is_empty()),
        "--joint-dna-rna currently requires the standard nucleotide E_r tiers; RUSTLE_SHARED_EXON is a different edge definition"
    );

    #[derive(Clone)]
    struct Node {
        family: usize,
        copy: usize,
        rep: DenovoTranscript,
    }

    // Match `copies.tsv` exactly: family order is already the emitted order, and copy_idx is coordinate-sorted.
    let mut nodes = Vec::new();
    for (fi, fam) in fams.iter().enumerate() {
        let mut sorted = fam.clone();
        sorted.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
        for (ci, rep) in sorted.into_iter().enumerate() {
            nodes.push(Node { family: fi, copy: ci, rep });
        }
    }
    let contigs: HashSet<String> = nodes.iter().map(|n| n.rep.chrom.clone()).collect();
    // An empty emitted catalog must produce empty joint reports cheaply, not trigger the FASTA loader's
    // intentional empty-contig fallback to loading the entire genome.
    let genome = if nodes.is_empty() {
        GenomeIndex::default()
    } else {
        GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?
    };
    let rna_seqs: Vec<Vec<u8>> = nodes.iter().map(|n| n.rep.seq.clone()).collect();
    let dna_seqs: Vec<Vec<u8>> = nodes.iter().map(|n| refine_copy_seq(&n.rep, Some(&genome))).collect();
    let core_lens: Vec<u64> = nodes.iter().map(|n| n.rep.core_bp).collect();

    let mut rna_params = params.clone();
    rna_params.homology_genomic_span = false;
    rna_params.protein_tail = false; // this certificate compares nucleotide evidence on two substrates
    let rna_edges = if nodes.is_empty() {
        BTreeSet::new()
    } else {
        nucleotide_rule_edge_set(&rna_seqs, Some(&core_lens), &rna_params, "joint_rna")?
    };

    let mut dna_params = params.clone();
    dna_params.require_forward_alignment = false;
    dna_params.substrate = Substrate::ReferenceOriented; // genomic spans, not transcripts
    dna_params.homology_genomic_span = true;
    dna_params.protein_tail = false;
    let dna_edges = if nodes.is_empty() {
        BTreeSet::new()
    } else {
        nucleotide_rule_edge_set(&dna_seqs, None, &dna_params, "joint_dna")?
    };

    let union: BTreeSet<(usize, usize)> = rna_edges.union(&dna_edges).copied().collect();
    let mut ef = std::fs::File::create(format!("{out}.joint_edges.tsv"))?;
    writeln!(
        ef,
        "family_i\tcopy_i\tchrom_i\tstart_i\tend_i\tstrand_i\tfamily_j\tcopy_j\tchrom_j\tstart_j\tend_j\tstrand_j\trna_edge\tdna_edge\tevidence_class\tsame_emitted_family"
    )?;
    let mut cross = 0usize;
    for &(i, j) in &union {
        let (a, b) = (&nodes[i], &nodes[j]);
        let rna = rna_edges.contains(&(i, j));
        let dna = dna_edges.contains(&(i, j));
        let class = match (rna, dna) {
            (true, true) => "RNA_DNA",
            (true, false) => "RNA_ONLY",
            (false, true) => "DNA_ONLY",
            (false, false) => unreachable!(),
        };
        let same = a.family == b.family;
        cross += (!same) as usize;
        writeln!(
            ef,
            "GWFAM{}\t{}\t{}\t{}\t{}\t{}\tGWFAM{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            a.family,
            a.copy,
            a.rep.chrom,
            a.rep.start,
            a.rep.end,
            a.rep.strand,
            b.family,
            b.copy,
            b.rep.chrom,
            b.rep.start,
            b.rep.end,
            b.rep.strand,
            rna,
            dna,
            class,
            same,
        )?;
    }

    let mut ff = std::fs::File::create(format!("{out}.joint_families.tsv"))?;
    writeln!(
        ff,
        "family_id\tn_copies\trna_edges\tdna_edges\tboth_edges\trna_only_edges\tdna_only_edges\tedge_jaccard\trna_components\tdna_components\tjoint_status\tkappa"
    )?;
    for (fi, fam) in fams.iter().enumerate() {
        let mut rna_local = Vec::new();
        let mut dna_local = Vec::new();
        let (mut both, mut rna_only, mut dna_only) = (0usize, 0usize, 0usize);
        for &(i, j) in &union {
            if nodes[i].family != fi || nodes[j].family != fi {
                continue;
            }
            let pair = (nodes[i].copy, nodes[j].copy);
            let rna = rna_edges.contains(&(i, j));
            let dna = dna_edges.contains(&(i, j));
            if rna {
                rna_local.push(pair);
            }
            if dna {
                dna_local.push(pair);
            }
            match (rna, dna) {
                (true, true) => both += 1,
                (true, false) => rna_only += 1,
                (false, true) => dna_only += 1,
                (false, false) => unreachable!(),
            }
        }
        let rc = joint_component_count(fam.len(), &rna_local);
        let dc = joint_component_count(fam.len(), &dna_local);
        let status = match (rc <= 1, dc <= 1) {
            (true, true) => "BOTH_CONNECTED",
            (true, false) => "RNA_CONNECTED_DNA_FRAGMENTED",
            (false, true) => "DNA_CONNECTED_RNA_FRAGMENTED",
            (false, false) => "BOTH_FRAGMENTED",
        };
        // The pre-declared certificate in docs/o1_investigations.md#the-joint-dna-rna-family-definition-retracted: RNA must both connect the
        // emitted family and contribute no edge absent from DNA. With no RNA edge there is no held-out
        // RNA test to score, so the family is UNTESTABLE rather than silently called concordant.
        let kappa = if rna_local.is_empty() {
            "UNTESTABLE"
        } else if rc <= 1 && rna_only == 0 {
            "CONCORDANT"
        } else {
            "DISCORDANT"
        };
        let union_n = both + rna_only + dna_only;
        let jaccard = if union_n == 0 { 1.0 } else { both as f64 / union_n as f64 };
        writeln!(
            ff,
            "GWFAM{fi}\t{}\t{}\t{}\t{both}\t{rna_only}\t{dna_only}\t{jaccard:.4}\t{rc}\t{dc}\t{status}\t{kappa}",
            fam.len(),
            rna_local.len(),
            dna_local.len(),
        )?;
    }

    let mut rf = std::fs::File::create(format!("{out}.joint_rule.tsv"))?;
    writeln!(rf, "key\tvalue")?;
    writeln!(rf, "node_universe\tRNA-detected emitted loci")?;
    writeln!(rf, "membership_effect\tnone (report/certificate only)")?;
    writeln!(rf, "rna_substrate\tspliced exon-sum, transcription orientation")?;
    writeln!(rf, "rna_orientation\t{}", if rna_params.require_forward_alignment { "forward-only (+)" } else { "both (+/-)" })?;
    writeln!(rf, "dna_substrate\tcomplete genomic span of the same RNA locus, transcript-normalised")?;
    writeln!(rf, "dna_orientation\tboth (+/-); inverted structural duplication remains valid")?;
    writeln!(rf, "protein_edges\texcluded (nucleotide substrate comparison)")?;
    writeln!(rf, "primary_identity\t{:.6}", er_primary_tier(params).1)?;
    writeln!(rf, "min_coverage\t{:.6}", params.min_coverage)?;
    writeln!(rf, "edge_classes\tRNA_DNA,RNA_ONLY,DNA_ONLY")?;
    writeln!(rf, "kappa\tCONCORDANT iff RNA is connected and E_RNA is a subset of E_DNA; UNTESTABLE iff RNA has no edge; otherwise DISCORDANT")?;
    eprintln!(
        "[joint-rna-dna] {} RNA edges, {} DNA edges, {} union edges ({} cross-family) -> {out}.joint_*.tsv",
        rna_edges.len(),
        dna_edges.len(),
        union.len(),
        cross,
    );
    Ok(())
}

/// Rebuild each locus representative from the UNION of its group's exons instead of its single best chain.
///
/// Off by default (`RUSTLE_LOCUS_EXON_UNION=1` to enable) — unset is byte-identical to the shipping path.
///
/// Measured against Soto truth windows and gorilla known gene families, unioning cut cross-family edges by
/// ~70% (Soto 164 -> 50) and broke up the pathological component that fused 40 of 83 Soto families into
/// 2-family pairs, while within-family edges rose 165 -> 448. The cost is loci: groups where no chain clears
/// `min_chain_reads` produce nothing (Soto 299 -> 263 loci, gorilla 75 -> 65).
///
/// The sequence is built DIRECTLY from the merged exons, deliberately bypassing `build_spliced_seq`'s
/// canonical-junction gate: the union's "introns" are gaps between merged exons and can be chimeric
/// (donor from one chain, acceptor from another), so re-gating would discard the locus over a junction no
/// read ever asserted. Each contributing chain already passed that gate upstream.
/// Junction-support floor for the co-threaded representative; `None` (default) leaves the path untouched.
///
/// `RUSTLE_COTHREAD_REP` is the switch and carries the floor: `1` selects the measured default of 3 reads
/// per junction, any other integer sets it directly. 3 is the floor the rest of the pipeline already uses
/// for a junction, and the co-observation constraint (not a higher floor) is what keeps the chain local --
/// see `cothread_locus_geometry`.
fn cothread_rep_floor() -> Option<u32> {
    let v = std::env::var("RUSTLE_COTHREAD_REP").ok()?;
    if v.is_empty() || v == "0" {
        return None;
    }
    Some(match v.parse::<u32>() {
        Ok(1) | Err(_) => 3,
        Ok(n) => n.max(1),
    })
}

/// Rebuild each locus representative as the maximum-weight path through its READ-WITNESSED splice graph
/// (`cothread_locus_geometry`), instead of picking the single best-supported observed chain.
///
/// Off by default (`RUSTLE_COTHREAD_REP=1` to enable) — unset is byte-identical to the shipping path.
///
/// This targets the single-exon representative, which is one defect with three symptoms: 74.6% of reps are
/// single-exon; a stub cannot cover half of a spliced sibling so the pair gets NO alignment record at all
/// (median coverage 0.000 stub-vs-stub, 0.099 stub-vs-spliced, 0.636 spliced-vs-spliced); and stub-
/// represented loci come out at 0.29x their true size against 0.75x for spliced ones. Verified against
/// RefSeq gene spans rather than duplication blocks: every member with a spliced rep lands at 0.86-1.04
/// (AMY1A/B/C all 0.86, GOLGA6L10 1.00, FAM72 0.88-1.03) and every badly undersized one has a 1-exon rep
/// (AMY2A 0.52, SRGAP2/B/C 0.03-0.05).
///
/// Measured on the 52 stub-represented Soto members, taking the constructed chain where one exists at the
/// floor and keeping today's representative otherwise: correctly-sized loci (within 0.5-2x of truth) go
/// from 15 to 26 at a floor of 3 reads per junction and to 28 at 5, while over-extension falls from the
/// unconstrained variant's 35% to 27% and 13%. SRGAP2 goes 0.05 -> 0.99 and SRGAP2C 0.03 -> 0.97.
///
/// ⚠ The representative feeds the E_r substrate, so edges, families and F1 all move. `RUSTLE_SPLICED_REP`
/// regressed by touching exactly this. Judge it on family detection, not only on sizes.
fn cothread_locus_reps(
    transcripts: &[DenovoTranscript],
    detect: &crate::vg_family::family_detect::DetectParams,
    genome: &GenomeIndex,
    min_reads: u32,
) -> Vec<DenovoTranscript> {
    use crate::vg_family::family_detect::{cothread_locus_geometry, locus_groups};
    use crate::vg_family::seq_utils::reverse_complement;
    let mut out = Vec::new();
    let (mut built, mut kept) = (0usize, 0usize);
    for members in locus_groups(transcripts, detect) {
        let rep_i = *members
            .iter()
            .max_by_key(|&&i| (transcripts[i].n_reads, transcripts[i].end - transcripts[i].start))
            .expect("locus group is never empty");
        let base = &transcripts[rep_i];
        // Fall back to today's representative whenever no read-witnessed chain clears the floor, so the
        // change can only add structure and never remove a locus.
        let Some((start, end, introns)) = cothread_locus_geometry(transcripts, &members, min_reads)
        else {
            kept += 1;
            out.push(base.clone());
            continue;
        };
        let mut seq: Vec<u8> = Vec::new();
        let mut prev = start;
        let mut ok = true;
        for &(d, a) in &introns {
            match genome.fetch_sequence(&base.chrom, prev, d) {
                Some(b) => seq.extend_from_slice(&b),
                None => {
                    ok = false;
                    break;
                }
            }
            prev = a;
        }
        if ok {
            match genome.fetch_sequence(&base.chrom, prev, end) {
                Some(b) => seq.extend_from_slice(&b),
                None => ok = false,
            }
        }
        if !ok || seq.is_empty() {
            kept += 1;
            out.push(base.clone());
            continue;
        }
        seq.make_ascii_uppercase();
        if base.strand == '-' {
            seq = reverse_complement(&seq);
        }
        built += 1;
        out.push(DenovoTranscript { start, end, introns, seq, ..base.clone() });
    }
    eprintln!(
        "[cothread-rep] {} locus reps: {} rebuilt from the splice graph, {} kept as-is \
         (no read-witnessed chain >= {} reads/junction)",
        out.len(),
        built,
        kept,
        min_reads
    );
    out
}

fn union_locus_reps(
    transcripts: &[DenovoTranscript],
    detect: &crate::vg_family::family_detect::DetectParams,
    genome: &GenomeIndex,
    min_chain_reads: u32,
) -> Vec<DenovoTranscript> {
    use crate::vg_family::family_detect::{locus_groups, union_locus_geometry};
    use crate::vg_family::seq_utils::reverse_complement;
    let mut out = Vec::new();
    let (mut dropped, mut widened) = (0usize, 0usize);
    for members in locus_groups(transcripts, detect) {
        let rep_i = *members
            .iter()
            .max_by_key(|&&i| (transcripts[i].n_reads, transcripts[i].end - transcripts[i].start))
            .expect("locus group is never empty");
        let base = &transcripts[rep_i];
        let Some((start, end, introns)) = union_locus_geometry(transcripts, &members, min_chain_reads)
        else {
            dropped += 1;
            continue;
        };
        let mut seq: Vec<u8> = Vec::new();
        let mut prev = start;
        let mut ok = true;
        for &(d, a) in &introns {
            match genome.fetch_sequence(&base.chrom, prev, d) {
                Some(b) => seq.extend_from_slice(&b),
                None => {
                    ok = false;
                    break;
                }
            }
            prev = a;
        }
        if ok {
            match genome.fetch_sequence(&base.chrom, prev, end) {
                Some(b) => seq.extend_from_slice(&b),
                None => ok = false,
            }
        }
        if !ok || seq.is_empty() {
            dropped += 1;
            continue;
        }
        seq.make_ascii_uppercase();
        if base.strand == '-' {
            seq = reverse_complement(&seq);
        }
        if seq.len() > base.seq.len() {
            widened += 1;
        }
        out.push(DenovoTranscript { start, end, introns, seq, ..base.clone() });
    }
    eprintln!(
        "[exon-union] {} locus reps ({} widened vs single-chain, {} dropped: no chain >= {} reads)",
        out.len(),
        widened,
        dropped,
        min_chain_reads
    );
    out
}

/// All-vs-all over arbitrary sequences, accepting a pair on identity plus ABSOLUTE aligned length rather
/// than the coverage-of-shorter fraction `nucleotide_edges` uses. The shared-exon rule needs this because a
/// single shared exon is meaningful at, say, 300 aligned bp regardless of how long either gene is; scaling
/// by the shorter sequence would make a short exon trivially pass and a long one trivially fail.
fn nucleotide_edges_indexed(
    seqs: &[Vec<u8>],
    mm_args: &[&str],
    min_id: f64,
    min_bp: u64,
    params: &RefineParams,
) -> Result<Vec<(usize, usize)>> {
    use std::io::Write;
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    static TMP_SEQ2: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let uniq = TMP_SEQ2.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = dir.join(format!("rustle_exon_{pid}_{uniq}.fa"));
    struct Cleanup(std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    let _guard = Cleanup(path.clone());
    {
        let mut fh = std::fs::File::create(&path)?;
        for (i, sq) in seqs.iter().enumerate() {
            writeln!(fh, ">{i}")?;
            fh.write_all(sq)?;
            writeln!(fh)?;
        }
    }
    let out = er_tier_command(&params.minimap2, params.threads, mm_args)
        .arg(&path)
        .arg(&path)
        .output()
        .map_err(|e| anyhow::anyhow!("failed to run minimap2 ('{}'): {e}", params.minimap2))?;
    if !out.status.success() {
        return Ok(Vec::new());
    }
    let text = String::from_utf8_lossy(&out.stdout);
    let mut edges: BTreeSet<(usize, usize)> = BTreeSet::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let (q, t) = match (f[0].parse::<usize>().ok(), f[5].parse::<usize>().ok()) {
            (Some(q), Some(t)) if q != t => (q, t),
            _ => continue,
        };
        let qs = f[2].parse::<f64>().unwrap_or(0.0);
        let qe = f[3].parse::<f64>().unwrap_or(0.0);
        let de = f[12..].iter().find_map(|x| x.strip_prefix("de:f:").and_then(|v| v.parse::<f64>().ok()));
        let ident = match de {
            Some(d) => 1.0 - d,
            None => {
                let nm = f[9].parse::<f64>().unwrap_or(0.0);
                let al = f[10].parse::<f64>().unwrap_or(1.0).max(1.0);
                nm / al
            }
        };
        if ident >= min_id && (qe - qs) as u64 >= min_bp {
            edges.insert((q.min(t), q.max(t)));
        }
    }
    Ok(edges.into_iter().collect())
}

/// SHARED-EXON edges: Soto's clustering criterion, as an alternative to the exon-sum rule.
///
/// Soto 2025 clusters SD-98 genes into families "based on shared exons and similar famCN (MAD<1) between
/// paralogs". The famCN half needs WGS read depth, which neither IsoSeq nor the reference provides, so this
/// implements the SEQUENCE half only: two loci are linked if ANY ONE of their exons aligns to any exon of
/// the other above `min_identity`, over at least `min_bp` aligned bases.
///
/// The contrast with `homology_edges_all_reps` is the point of having both:
///   exon-sum rule  one alignment must cover >= 50% of the SHORTER WHOLE gene model
///   shared-exon    one exon pair suffices, with no whole-gene coverage requirement
/// So shared-exon is strictly more permissive on partial homology and should merge more. Whether that is
/// recall or over-merge is exactly what comparing them measures.
///
/// Exon sequences are sliced out of each rep's stored exon-sum rather than re-fetched, so this sees exactly
/// the bases the exon-sum rule sees. On the '-' strand the stored sequence is the reverse complement of the
/// concatenation, so the exon LENGTHS are reversed before slicing; each slice is then that exon in reverse
/// complement, which is immaterial to alignment since minimap2 tries both strands.
/// Slice one transcript's stored exon-sum back into its individual exon sequences.
///
/// ⚠ The stored `seq` is the CONCATENATED exons, and on the '-' strand it is the reverse complement of that
/// concatenation -- so the exon LENGTHS must be reversed before slicing. Getting this wrong silently returns
/// the right number of exons with the wrong bases. Factored out so the single-representative and pooled
/// callers cannot drift apart on it.
fn exon_seqs_of(t: &DenovoTranscript, min_bp: u64) -> Vec<Vec<u8>> {
    let mut lens: Vec<usize> = Vec::with_capacity(t.introns.len() + 1);
    let mut prev = t.start;
    for &(d, a) in &t.introns {
        lens.push(d.saturating_sub(prev) as usize);
        prev = a;
    }
    lens.push(t.end.saturating_sub(prev) as usize);
    if t.strand == '-' {
        lens.reverse();
    }
    let mut out = Vec::new();
    let mut off = 0usize;
    for l in lens {
        if l == 0 || off + l > t.seq.len() {
            off += l;
            continue;
        }
        if l as u64 >= min_bp {
            out.push(t.seq[off..off + l].to_vec());
        }
        off += l;
    }
    out
}

/// Shared-exon edges where each locus contributes the exons of EVERY isoform collapsed into it, not only
/// those of its chosen representative (`RUSTLE_SHARED_EXON_ISOFORMS`).
///
/// WHY. The pipeline keeps one representative per locus and discards the rest, and 46% of the
/// representatives that cover a known family member are single-exon stubs. Measured at those stub loci:
/// 95% still have spliced reads, and 53% have a spliced chain carried by >= 3 reads -- i.e. a gate-passing
/// spliced model exists and we threw it away. NOTCH2NLA, SRGAP2C and LIMS1 are represented by unspliced
/// stubs while 92, 65 and 124 reads respectively support a proper spliced model at the same locus. A stub
/// cannot cover half of a full transcript, so the pair never gets an edge even though both loci were found.
///
/// WHY THIS IS NOT `RUSTLE_LOCUS_EXON_UNION`, which lost 20 recall points: that concatenated the isoforms
/// into ONE LONGER sequence, and coverage is aligned-span / shorter-length, so a longer representative
/// inflated its own denominator and edges vanished. Here the exons stay SEPARATE and are matched
/// any-to-any, so no sequence gets longer and no denominator moves. It also does not touch representative
/// selection, so it cannot regress the way `RUSTLE_SPLICED_REP` did.
///
/// ⚠ Judge it on PRECISION: more exons per locus means more chances to match, including on shared repeats
/// and common domains. The question is whether the NEW edges are true co-family pairs.
pub fn shared_exon_edges_pooled(
    reps: &[DenovoTranscript],
    pooled: &[Vec<Vec<u8>>],
    min_identity: f64,
    min_bp: u64,
    params: &RefineParams,
) -> Result<Vec<(usize, usize)>> {
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut owner: Vec<usize> = Vec::new();
    for (ri, exs) in pooled.iter().enumerate() {
        for e in exs {
            if e.len() as u64 >= min_bp {
                seqs.push(e.clone());
                owner.push(ri);
            }
        }
    }
    edges_from_exon_pool(seqs, owner, reps, min_identity, min_bp, params, "shared-exon-pooled")
}

/// Shared engine for both shared-exon variants: one all-vs-all over every exon, then lift exon pairs to
/// locus pairs.
fn edges_from_exon_pool(
    seqs: Vec<Vec<u8>>,
    owner: Vec<usize>,
    reps: &[DenovoTranscript],
    min_identity: f64,
    min_bp: u64,
    params: &RefineParams,
    tag: &str,
) -> Result<Vec<(usize, usize)>> {
    if seqs.len() < 2 {
        return Ok(Vec::new());
    }
    // How many DISTINCT exon pairs must support a locus pair before it becomes an edge
    // (`RUSTLE_SHARED_EXON_MIN_COUNT`, default 1 = the original any-one-exon rule).
    //
    // One shared exon is weak evidence: a single conserved domain or an exonised repeat links two loci that
    // are not paralogs. Measured on chr1+chr15 (HUMAN) with every isoform's exons pooled, the any-one-exon
    // rule admitted 209 new locus pairs of which only 32 (15%) were true co-family pairs, and raising the
    // length floor did not fix it (25% at 600 bp, and F1 never beat the representative-only baseline).
    //
    // ⚠⚠ THE "REQUIRING SEVERAL INDEPENDENT EXONS IS THE UNTESTED AXIS" CLAIM THAT STOOD HERE IS RETRACTED
    // (2026-08-22). It is NOT untested and it does NOT distinguish structure from element:
    //
    //  * MEASURED 2026-08-03, HUMAN chr1+chr15, through the REAL BINARY (hp/iso_count.sh against the matched
    //    hp/iso_minbp.sh MB_300 control, id=0.70 min_bp=300 gamma=0.20 ISOFORMS=1), unit = EMITTED copies /
    //    EMITTED families:  MIN_COUNT 1 -> 545/121 | 2 -> 323/62 | 4 -> 243/43 | 8 -> 197/43.
    //    Monotone loss on the emitted partition. DEAD-END.
    //  * MEASURED 2026-08-22, GORILLA genome-wide, rep-only, unit = NODE PAIR: after deleting exon matches
    //    that fall on the SAME genomic interval, the rule adds 376 pairs at MIN_COUNT>=1, 9 at >=2, 0 at >=3.
    //  * ⚠ AND THE COUNTER IS NOT COUNTING WHAT THE NAME SAYS: under ISOFORMS=1 `support` accumulates exon
    //    INDEX pairs over an UN-DEDUPLICATED pool, so ONE genomic exon pair contributes n_iso_a * n_iso_b.
    //    MIN_COUNT there measures ISOFORM MULTIPLICITY, not exon independence. Every MIN_COUNT figure ever
    //    quoted under ISOFORMS=1 is void as a statement about independent exons (it remains valid as the
    //    partition COST above, which is what actually kills the axis).
    let min_shared: usize = std::env::var("RUSTLE_SHARED_EXON_MIN_COUNT")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1);
    let mut support: std::collections::BTreeMap<(usize, usize), usize> = std::collections::BTreeMap::new();
    let pairs = nucleotide_edges_indexed(&seqs, &["-x", "asm20"], min_identity, min_bp, params)?;
    let mut self_overlap = 0usize;
    for (a, b) in pairs {
        let (ra, rb) = (owner[a], owner[b]);
        if ra == rb {
            continue;
        }
        // DISTINCT-LOCUS GUARD. `ra != rb` alone is not enough: two representatives whose genomic spans
        // OVERLAP are one locus seen twice, and their exons are then the SAME PHYSICAL DNA, so they align
        // to themselves at ~100% identity and manufacture an edge that asserts no homology at all.
        //
        // Measured genome-wide on gorilla (2026-08-22, rep-only, id >= 0.98, >= 100 bp): of the 2,795 node
        // pairs this rule reported that E_r does not, 2,601 = 93.06% had OVERLAPPING SPANS and 2,419 =
        // 86.55% had every matched exon on the same physical interval -- 84.01% of them additionally on
        // OPPOSITE strands, i.e. a locus matching its own antisense annotation. The in-E_r pairs invert the
        // statistic (50/396 = 12.63% overlap), so this is a property of the pairs the rule ADDS, not of the
        // measurement. Removing them takes the genuine new payload from 2,795 to 376 node pairs.
        //
        // The same_pos test mirrors `distinct_locus_reps_grouped`, so the two agree on what "one locus" is.
        let (x, y) = (&reps[ra], &reps[rb]);
        if x.chrom == y.chrom && x.end.min(y.end) > x.start.max(y.start) {
            self_overlap += 1;
            continue;
        }
        *support.entry((ra.min(rb), ra.max(rb))).or_insert(0) += 1;
    }
    if self_overlap > 0 {
        eprintln!("[{tag}] dropped {self_overlap} same-locus exon matches (overlapping rep spans)");
    }
    let edge_set: BTreeSet<(usize, usize)> = support
        .into_iter()
        .filter(|&(_, n)| n >= min_shared)
        .map(|(k, _)| k)
        .collect();
    eprintln!(
        "[{tag}] {} exons over {} loci -> {} locus pairs linked (id >= {}, >= {} bp aligned)",
        seqs.len(),
        reps.len(),
        edge_set.len(),
        min_identity,
        min_bp
    );
    Ok(edge_set.into_iter().collect())
}

pub fn shared_exon_edges(
    reps: &[DenovoTranscript],
    min_identity: f64,
    min_bp: u64,
    params: &RefineParams,
) -> Result<Vec<(usize, usize)>> {
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut owner: Vec<usize> = Vec::new();
    for (ri, r) in reps.iter().enumerate() {
        for e in exon_seqs_of(r, min_bp) {
            seqs.push(e);
            owner.push(ri);
        }
    }
    edges_from_exon_pool(seqs, owner, reps, min_identity, min_bp, params, "shared-exon")
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
    homology_edges_all_reps_pooled(reps, None, params)
}

/// As `homology_edges_all_reps`, but when `pooled` is supplied the shared-exon rule draws on EVERY isoform's
/// exons for each locus instead of only the representative's. See `shared_exon_edges_pooled`.
/// Threshold for the genome-anchored repeat-hub veto. `RUSTLE_ER_REPEAT_GATE=<M>`; unset or 0 = OFF.
fn repeat_gate_threshold() -> Option<u32> {
    std::env::var("RUSTLE_ER_REPEAT_GATE")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .filter(|&m| m > 0)
}

/// Genome multiplicity of each representative: the number of DISTINCT, non-overlapping places in the
/// REFERENCE where that sequence occurs. Consumed only by the repeat-hub veto.
///
/// WHY GENOME-ANCHORED AND NOT CATALOG-ANCHORED. Counting a rep's partners inside the catalog (the R5
/// design) makes the gate a function of which other reps happen to be present, and that broke P1
/// seed-invariance at 94/147. Multiplicity here is a property of the rep and the reference ALONE, so
/// adding or removing other reps cannot move it and P1 holds by construction.
///
/// ⚠ `-X` is deliberately ABSENT. `-X` is minimap2's all-vs-all mode and implies `--dual=no`; this is a
/// query-vs-reference mapping, where it would be wrong. `-p` and `-N` are stated explicitly because a
/// copy count is meaningless without them (MAPKBP1 gives 1/1 at `-p 0.8` and 9/8 at `-p 0.1`).
fn genome_multiplicity(seqs: &[Vec<u8>], genome_path: &str, params: &RefineParams) -> Result<Vec<u32>> {
    use std::io::{BufRead, Write};
    const MIN_ID: f64 = 0.90;
    const MIN_COV: f64 = 0.50;
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    static TMP_SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let uniq = TMP_SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = dir.join(format!("rustle_gmult_{pid}_{uniq}.fa"));
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
    let mut cmd = std::process::Command::new(&params.minimap2);
    cmd.args(["-x", "asm20", "-c", "--eqx", "-N", "200", "-p", "0.1"])
        .arg("-t")
        .arg(params.threads.max(1).to_string())
        .arg(genome_path)
        .arg(&path)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::null());
    let mut child = cmd.spawn().map_err(|e| {
        anyhow::anyhow!("failed to run minimap2 ('{}') for the repeat gate: {e}", params.minimap2)
    })?;
    let out = child.stdout.take().ok_or_else(|| anyhow::anyhow!("minimap2 stdout unavailable"))?;
    let mut hits: Vec<Vec<(String, u64, u64)>> = vec![Vec::new(); seqs.len()];
    for line in std::io::BufReader::new(out).lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 {
            continue;
        }
        let (Ok(q), Ok(ql), Ok(qs), Ok(qe)) =
            (f[0].parse::<usize>(), f[1].parse::<f64>(), f[2].parse::<f64>(), f[3].parse::<f64>())
        else {
            continue;
        };
        let (Ok(ts), Ok(te), Ok(nm), Ok(bl)) =
            (f[7].parse::<u64>(), f[8].parse::<u64>(), f[9].parse::<f64>(), f[10].parse::<f64>())
        else {
            continue;
        };
        if q >= hits.len() || ql <= 0.0 || bl <= 0.0 {
            continue;
        }
        if nm / bl < MIN_ID || (qe - qs) / ql < MIN_COV {
            continue;
        }
        hits[q].push((f[5].to_string(), ts, te));
    }
    let _ = child.wait();
    Ok(hits.into_iter().map(count_distinct_loci).collect())
}

/// Collapse a rep's reference hits into DISTINCT occurrences: sort, then merge anything that overlaps on
/// the same contig. Split out of `genome_multiplicity` so the counting rule can be tested without
/// invoking minimap2 -- an off-by-one here silently moves every multiplicity, and with it the veto.
fn count_distinct_loci(mut v: Vec<(String, u64, u64)>) -> u32 {
    v.sort();
    let mut n = 0u32;
    let mut cur_t = String::new();
    let mut cur_e = 0u64;
    let mut started = false;
    for (t, s, e) in v {
        if !started || t != cur_t || s >= cur_e {
            n += 1;
            cur_t = t;
            cur_e = e;
            started = true;
        } else if e > cur_e {
            cur_e = e;
        }
    }
    n
}

/// `homology_edges_all_reps_pooled`, but also returning each edge's exemplar (identity, coverage).
///
/// Added because the partition DISCARDED those numbers: it flattened every edge to weight 1.0 before
/// handing them to Louvain, which is a WEIGHTED-modularity algorithm. Measured on the 3-contig NPIP
/// catalog, identity separates NPIP<->NPIP edges from NPIP<->non-NPIP at AUC 0.9491 (median 0.9878 vs
/// 0.8281), so the flattening throws away a strong signal. Consumed only under
/// `RUSTLE_ER_WEIGHTED_PARTITION`; the unweighted wrapper below keeps every existing caller byte-identical.
pub(crate) fn homology_edges_all_reps_pooled(
    reps: &[DenovoTranscript],
    pooled: Option<&[Vec<Vec<u8>>]>,
    params: &RefineParams,
) -> Result<Vec<(usize, usize)>> {
    Ok(homology_edges_all_reps_pooled_weighted(reps, pooled, params)?
        .into_iter()
        .map(|(a, b, _, _)| (a, b))
        .collect())
}

pub(crate) fn homology_edges_all_reps_pooled_weighted(
    reps: &[DenovoTranscript],
    pooled: Option<&[Vec<Vec<u8>>]>,
    params: &RefineParams,
) -> Result<Vec<(usize, usize, f64, f64)>> {
    // EDGE SUBSTRATE. By default the E_r edge is computed on the ASSEMBLED exon-sum (`rep.seq`). That makes
    // family membership hostage to assembly completeness: two copies of one family whose transcripts assemble
    // DISJOINT exon subsets (one gets exons 1-2, its sibling 7-8) share almost no sequence, so no edge forms
    // and the family fragments. Measured on the Soto benchmark: 75% of should-link pairs had NO exon-sum
    // alignment at all, and those that aligned covered a median of only 12% of the shorter sequence.
    //
    // `homology_genomic_span` (opt-in) instead computes the edge on the GENOMIC SPAN of each RNA-detected
    // locus — complete sequence regardless of what was assembled. RNA still decides WHERE (which loci are
    // expressed); the reference supplies WHAT-GROUPS-WITH-WHAT. Needs no read depth / copy number.
    // Simulated effect (863 RNA loci, id>=0.90 cov>=0.50): completeness 25%->43%, splits 45->32, with
    // homogeneity held at ~90% (the node set stays RNA-limited, so there is little repeat bridging to admit).
    let span_genome: Option<GenomeIndex> = if params.homology_genomic_span {
        match params.intron_fasta.as_ref() {
            Some(path) => {
                let contigs: HashSet<String> = reps.iter().map(|r| r.chrom.clone()).collect();
                match GenomeIndex::from_fasta_contigs(path, &contigs) {
                    Ok(g) => Some(g),
                    Err(e) => {
                        // best-effort: fall back to exon-sum edges rather than fail the catalog
                        eprintln!("[homology] genomic-span edges skipped: genome load failed ({e})");
                        None
                    }
                }
            }
            None => {
                eprintln!("[homology] genomic-span edges requested but no genome path set; using exon-sum");
                None
            }
        }
    } else {
        None
    };
    if span_genome.is_some() {
        eprintln!("[homology] E_r edge substrate: GENOMIC SPAN ({} reps)", reps.len());
    }
    let seqs: Vec<Vec<u8>> = reps.iter().map(|r| refine_copy_seq(r, span_genome.as_ref())).collect();
    // Parallel to `seqs` by construction (both are built 1:1 from `reps`), which is what lets
    // `nucleotide_edges` index cores by the same PAF sequence id.
    let core_lens: Vec<u64> = reps.iter().map(|r| r.core_bp).collect();
    // SCOPED ORIENTATION GUARD (`RUSTLE_ER_GUARD_SCOPED`, default OFF = byte-identical).
    //
    // The guard rejects every minus-strand record. That is right when both strands were MEASURED, and
    // vacuous when one was not: `denovo_assemble.rs` stamps `strand.unwrap_or('+')` on any rep with no
    // canonical junction, so all 5,928 single-exon reps carry `'+'` and none carries `'-'`. Such a rep may
    // simply be stored backwards, in which case its minus alignment is an artefact of the placeholder.
    //
    // Measured on the shipped catalog (rep-pair unit): of 4,009 guard-blocked pairs, 3,951 = 0.9855
    // involve a rep with an unmeasured strand, leaving 58 genuine antisense candidates genome-wide. And
    // among the stub loci that DO have a unique containing spliced model, 123/124 = 0.9919 of the
    // exonic-upgrade cases fail here rather than on identity or coverage, with a placeholder-strand
    // partner in 0.9535 of them.
    //
    // ⚠ THE SPAN-OVERLAP CLAUSE IS LOAD-BEARING, NOT TIDINESS. Exempting on the unmeasured strand alone
    // would admit 3,951 pairs of which 1,727 = 0.4371 have OVERLAPPING SPANS -- one locus counted twice,
    // the dominant known artefact class -- against 0.0109 in the shipped edge set, i.e. 40x enriched.
    // Requiring disjoint spans leaves 2,224 pairs whose SEDEF support is 0.2350 against the shipped set's
    // 0.3025 (0.78x): weaker than what ships, but not the artefact population.
    // ⛔⛔ MEASURED AND REJECTED 2026-08-25 — DO NOT ENABLE. The exemption is sound in principle (the
    // guard is testing a field that was never measured) and catastrophic in effect: of the 2,224 pairs it
    // admits, 1,913 join two reps that are ALREADY catalog copies, and 1,912/1,913 = 0.9995 of those join
    // two DIFFERENT families. They would fuse 112 family pairs, touching 77/627 = 0.1228 of the catalog,
    // collapsing it into 17 blocks whose largest swallows 43 families. Retained, off, because the flag is
    // the record of the experiment. See `docs/o1_ledger.md` §4t.
    //
    // ⚠ The HUMAN 150-window negative panel CANNOT adjudicate this: it offers 0 qualifying pairs (its
    // windows are single-locus, so its 25 guard-blocked pairs are all CO-LOCATED and the disjoint-span
    // clause excludes them). Its 2/150 result under this flag is uninformative, not a pass.
    let guard_exempt: Option<Box<dyn Fn(usize, usize) -> bool>> =
        if matches!(std::env::var("RUSTLE_ER_GUARD_SCOPED"), Ok(v) if v != "0" && !v.is_empty()) {
            let meta: Vec<(String, u64, u64, bool)> = reps
                .iter()
                .map(|r| (r.chrom.clone(), r.start, r.end, r.introns.is_empty()))
                .collect();
            Some(Box::new(move |a: usize, b: usize| {
                let (Some(x), Some(y)) = (meta.get(a), meta.get(b)) else { return false };
                let unmeasured = x.3 || y.3; // a junction-less rep has no measured strand
                let overlap = x.0 == y.0 && x.2.min(y.2) > x.1.max(y.1);
                unmeasured && !overlap
            }))
        } else {
            None
        };
    let drop_stub_edges = er_no_stub_edges();
    if let Some(fl) = core_cov_floor() {
        let have = core_lens.iter().filter(|&&c| c > 0).count();
        eprintln!("[homology] E_r coverage denominator: READ-SUPPORTED CORE (depth >= {}, floor {fl:.2}); \
                   measured for {have}/{} reps, rest fall back to span", core_depth_floor(), reps.len());
    }
    // SHARED-EXON mode (`RUSTLE_SHARED_EXON=1`): replace the exon-sum rule with Soto's criterion --
    // any single shared exon links two loci. Not a tier that unions in; it REPLACES the nucleotide runs,
    // so the two definitions can be compared rather than blended.
    if std::env::var("RUSTLE_SHARED_EXON").map(|v| v != "0" && !v.is_empty()).unwrap_or(false) {
        let se_id: f64 = std::env::var("RUSTLE_SHARED_EXON_IDENTITY")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.98);
        let se_bp: u64 = std::env::var("RUSTLE_SHARED_EXON_MIN_BP")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(100);
        eprintln!("[shared-exon] MODE ACTIVE: replacing the exon-sum rule (id >= {se_id}, >= {se_bp} bp)");
        let se = match pooled {
            Some(p) => shared_exon_edges_pooled(reps, p, se_id, se_bp, params)?,
            None => shared_exon_edges(reps, se_id, se_bp, params)?,
        };
        // shared-exon mode reports no per-edge exemplar, so every weight defaults to 1.0
        return Ok(se.into_iter().map(|(a, b)| (a, b, 1.0, 1.0)).collect());
    }
    let mut prov: BTreeMap<(usize, usize), TierMask> = BTreeMap::new();
    // Exemplar (identity, coverage, tier-that-reported-them) per edge, consumed by `RUSTLE_ER_EDGE_DUMP`
    // only. Populated whether or not the dump is on (at most |E| f64 pairs), so enabling the dump cannot
    // change which code paths run.
    let mut metrics: BTreeMap<(usize, usize), (f64, f64, TierMask)> = BTreeMap::new();
    // §6az ADDITIVE DISCLOSURE: `cov_longer` + both unaligned flanks, from the SAME exemplar record as
    // `metrics`. A SEPARATE map keyed identically, so the tuple that feeds the edge decision is untouched.
    // It gates nothing; `write_er_edge_dump` is its only consumer.
    let mut flanks: BTreeMap<(usize, usize), ErEdgeFlank> = BTreeMap::new();
    // ⭐ SENSITIVE-ONLY IS NOW THE DEFAULT on this path; set `RUSTLE_ER_SENSITIVE_ONLY=0` to restore the
    // asm20 run. asm20 is a SUBSET of the sensitive run, structurally and measured:
    //   * structurally: whenever `sensitive_identity <= min_identity` the sensitive run has BOTH the lower
    //     identity bar AND the denser seeding (-k11 -w5 vs asm20's -k19 -w19), so it cannot miss an edge
    //     the coarser, stricter run finds.
    //   * Soto (--homology-primary, sensitive floor 0.70): asm20 produced 2894 edges, **0 unique**; the
    //     sensitive run added 1643 on its own.
    //   * NPIP 2026-08-06 (27 discovered loci, GENOMIC): asm20 2162, sensitive 2210, union 2210,
    //     asm20-only **0**.  TBC1D3 (11 loci, genomic): asm20 55, sensitive 55, asm20-only **0**.
    // ⚠ NOT EXACTLY 0 ON THE SPLICED SUBSTRATE. Same 27 NPIP loci, read-supported exons instead of
    //   genomic span: asm20 150, sensitive 269, union **270** -- asm20 contributed ONE sole edge.
    //   The subset argument above is about seeding density and identity floor, both of which hold; the
    //   single exception comes from chaining on a fragmented (concatenated-exon) target, where the
    //   coarser -k19 -w19 seeding can chain across a junction the denser seeding splits. One edge in
    //   270 does not justify a second genome-wide pass, but the claim is "0 on genomic, 1 on spliced",
    //   not "0".
    // Skipping it removes one whole genome-wide all-vs-all from the edge step.
    // ⚠ PATH-SPECIFIC. On the deprecated `--cross-chrom` refine path asm20 runs FIRST and the sensitive run
    // is conditional on `edges_connect_all`, so THERE asm20 carried sole support in 139 families. That path
    // is retired; do not generalise this default to it.
    // X.4: the predicate moved to `er_sensitive_only` so `refine_families_exon_sum` reads the SAME one.
    // The value here is unchanged.
    let sensitive_only = er_sensitive_only(params);
    if !sensitive_only {
        let seed = primary_seed_args();
        let seed_ref: Vec<&str> = seed.iter().map(String::as_str).collect();
        let mut tier_flanks: BTreeMap<(usize, usize), ErEdgeFlank> = BTreeMap::new();
        let scored = nucleotide_edges_scored_disclosed(&seqs, &seed_ref, params.min_identity, params.min_coverage, Some(&core_lens), params, Some("er"), guard_exempt.as_ref().map(|f| f as &dyn Fn(usize, usize) -> bool), Some(&mut tier_flanks))?;
        for (i, j, ident, cov) in scored {
            *prov.entry((i, j)).or_insert(0) |= TIER_ASM20;
            if record_edge_metric(&mut metrics, (i, j), ident, cov, TIER_ASM20) {
                record_edge_flank(&mut flanks, (i, j), tier_flanks.get(&(i, j)).copied());
            }
        }
    } else {
        eprintln!("[homology] asm20 run SKIPPED (RUSTLE_ER_SENSITIVE_ONLY)");
    }
    let mut set: BTreeSet<(usize, usize)> = prov.keys().copied().collect();
    if params.nucleotide_sensitive {
        let mut tier_flanks: BTreeMap<(usize, usize), ErEdgeFlank> = BTreeMap::new();
        let scored = nucleotide_edges_scored_disclosed(
            &seqs,
            ER_SENSITIVE_SEED,
            params.sensitive_identity,
            params.min_coverage,
            Some(&core_lens),
            params,
            Some("er"), None, Some(&mut tier_flanks))?;
        for (i, j, ident, cov) in scored {
            set.insert((i, j));
            *prov.entry((i, j)).or_insert(0) |= TIER_SENSITIVE;
            if record_edge_metric(&mut metrics, (i, j), ident, cov, TIER_SENSITIVE) {
                record_edge_flank(&mut flanks, (i, j), tier_flanks.get(&(i, j)).copied());
            }
        }
    }
    if params.protein_tail {
        let prot = batch_protein_edges(std::slice::from_ref(&reps.to_vec()), 0.50, params.min_coverage, params)?;
        if let Some(edges) = prot.first() {
            for &e in edges {
                set.insert(e);
                *prov.entry(e).or_insert(0) |= TIER_PROTEIN;
            }
        }
    }
    // Tier accounting over the whole rep set: how many edges each tier found, and how many it is the SOLE
    // source of. The second number is the one that decides whether a tier earns its place -- an edge that
    // another tier also finds costs a subprocess and changes nothing.
    let mut total: BTreeMap<TierMask, usize> = BTreeMap::new();
    let mut only: BTreeMap<TierMask, usize> = BTreeMap::new();
    for &m in prov.values() {
        for bit in [TIER_ASM20, TIER_SENSITIVE, TIER_GENOMIC, TIER_PROTEIN] {
            if m & bit != 0 {
                *total.entry(bit).or_insert(0) += 1;
            }
        }
        if m.count_ones() == 1 {
            *only.entry(m).or_insert(0) += 1;
        }
    }
    let fmt: Vec<String> = total
        .iter()
        .map(|(&b, &n)| format!("{}={} (sole {})", tier_names(b), n, only.get(&b).copied().unwrap_or(0)))
        .collect();
    eprintln!("[provenance] E_r edges by tier: {} | total {}", fmt.join(" "), prov.len());
    if let Ok(path) = std::env::var("RUSTLE_EDGE_PROVENANCE") {
        use std::io::Write;
        if let Ok(mut fh) = std::fs::File::create(&path) {
            let _ = writeln!(fh, "rep_i\trep_j\tchrom_i\tstart_i\tchrom_j\tstart_j\ttiers");
            for (&(i, j), &m) in prov.iter() {
                let _ = writeln!(
                    fh,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    i, j, reps[i].chrom, reps[i].start, reps[j].chrom, reps[j].start, tier_names(m)
                );
            }
            eprintln!("[provenance] wrote per-edge tiers to {path}");
        }
    }
    if drop_stub_edges {
        let before = set.len();
        set.retain(|&(i, j)| !reps[i].stub && !reps[j].stub);
        let n_stub = reps.iter().filter(|r| r.stub).count();
        let n_intronless = reps.iter().filter(|r| r.introns.is_empty() && !r.stub).count();
        eprintln!(
            "[homology] stub-edge guard: {n_stub} stubs and {n_intronless} genuinely intronless reps; \
             edges {before} -> {} (loci all retained)",
            set.len()
        );
    }
    // GENOME-ANCHORED REPEAT-HUB VETO (`RUSTLE_ER_REPEAT_GATE=<M>`, unset = OFF = byte-identical).
    //
    // WHY THIS EXISTS. `family_define` has carried a repeat-hub gate (`min_shared_mult >= 20`) ON BY
    // DEFAULT since the annotation-driven era; this de novo catalog path has never had one. Measured on
    // the 3-contig full-depth fibroblast NPIP catalog (`docs/o1_ledger.md` §4v): among single-exon reps
    // that FORM an edge, 20/227 = 0.0881 occur at >= 20 distinct genomic loci, against 0/1007 = 0.0000
    // among the single-exon reps that form none. High genome multiplicity is therefore close to
    // diagnostic of a stub edge -- and the >= 20 cut was reached from that distribution BEFORE the
    // shipped gate's own threshold was known, i.e. two independent routes to the same number.
    //
    // ⚠ A VETO, NEVER AN ADMISSION CRITERION -- it only ever REMOVES edges. Parity with gamma would need
    // M ~ 3, which discards 48% of the shipped edge set.
    // ⚠ Cuts only when BOTH endpoints are hubs, matching the shipped `min_shared_mult >= M` semantic
    // (the MINIMUM over the pair clears M). A one-sided cut would remove edges from a hub to a
    // single-copy locus, which is the signature of a real dispersed family, not of a repeat.
    if let Some(m) = repeat_gate_threshold() {
        match params.intron_fasta.as_deref() {
            Some(g) => match genome_multiplicity(&seqs, g, params) {
                Ok(mult) => {
                    let before = set.len();
                    set.retain(|&(i, j)| {
                        !(mult.get(i).copied().unwrap_or(0) >= m
                            && mult.get(j).copied().unwrap_or(0) >= m)
                    });
                    let hubs = mult.iter().filter(|&&x| x >= m).count();
                    eprintln!(
                        "[homology] repeat-hub veto (genome-anchored, both endpoints >= {m} loci): \
                         {hubs}/{} reps are hubs; edges {before} -> {} (loci all retained)",
                        mult.len(),
                        set.len()
                    );
                }
                Err(e) => eprintln!("[homology] repeat-hub veto SKIPPED (kept all edges): {e}"),
            },
            None => eprintln!("[homology] repeat-hub veto SKIPPED (kept all edges): no genome path"),
        }
    }
    // Opt-in parity dump. Placed HERE, after the stub guard, so what it writes is exactly what this
    // function RETURNS -- the pre-existing `RUSTLE_EDGE_PROVENANCE` writer above dumps `prov`, which
    // `set.retain` never filters, so under `RUSTLE_ER_NO_STUB_EDGES=1` that file is a strict superset.
    if let Ok(prefix) = std::env::var("RUSTLE_ER_EDGE_DUMP") {
        if !prefix.is_empty() {
            write_er_edge_dump(&prefix, reps, &seqs, &set, &prov, &metrics, &flanks, params, span_genome.is_some());
        }
    }
    Ok(set
        .into_iter()
        .map(|(a, b)| {
            let (i, c, _) = metrics.get(&(a, b)).copied().unwrap_or((1.0, 1.0, TierMask::default()));
            (a, b, i, c)
        })
        .collect())
}

/// Keep the highest-coverage exemplar per edge across tiers (ties -> higher identity, then lower tier bit).
/// Deterministic and reporting-only; see `nucleotide_edges_scored`.
///
/// Returns whether THIS call's record is now the stored exemplar, so the §6az disclosure map can follow
/// the identical choice instead of re-deriving it (a second copy of this predicate is exactly how the
/// disclosed flanks would come to describe a different record than the row's identity/coverage).
fn record_edge_metric(
    metrics: &mut BTreeMap<(usize, usize), (f64, f64, TierMask)>,
    key: (usize, usize),
    ident: f64,
    cov: f64,
    tier: TierMask,
) -> bool {
    match metrics.entry(key) {
        std::collections::btree_map::Entry::Vacant(v) => {
            v.insert((ident, cov, tier));
            true
        }
        std::collections::btree_map::Entry::Occupied(mut o) => {
            if (cov, ident) > (o.get().1, o.get().0) {
                o.insert((ident, cov, tier));
                true
            } else {
                false
            }
        }
    }
}

/// Move the §6az disclosure for `key` to whatever `record_edge_metric` just accepted as the exemplar.
///
/// `None` means the winning tier reported no certifying record for this pair (the opt-in summed-coverage
/// rule). The stale entry is then REMOVED rather than left behind: an inherited flank would describe some
/// other tier's record while the printed identity/coverage described this one, which is precisely the
/// inconsistency this disclosure exists to avoid.
fn record_edge_flank(
    flanks: &mut BTreeMap<(usize, usize), ErEdgeFlank>,
    key: (usize, usize),
    flank: Option<ErEdgeFlank>,
) {
    match flank {
        Some(f) => {
            flanks.insert(key, f);
        }
        None => {
            flanks.remove(&key);
        }
    }
}

/// OPT-IN E_r EDGE DUMP (`RUSTLE_ER_EDGE_DUMP=<prefix>`, unset = not one byte changes).
///
/// WHY THIS EXISTS. `gw_family_catalog` emits `copies.tsv` and `families.tsv` but never its EDGE set, and
/// the edge set is NOT recoverable from them: `distinct_locus_reps` merges/drops reps after blocks are
/// formed and `emit_catalog` re-sorts families and assigns a fresh per-family `copy_idx`, so the node
/// numbering the edges are expressed in does not survive into any shipped output. That made parity between
/// this engine and the canonical Python mirror (`bench/soto/rustlib.py`) UNVERIFIABLE: a test could only
/// compare Python against another Python harness while claiming to check the Rust. This dump is the missing
/// artefact — it is the Rust's own answer, in a form Python can diff against.
///
/// Three files, all tab-separated, all with a header, all sorted (edges by `(i, j)`, nodes by index,
/// params in a fixed order) so a `diff` is stable across runs:
///   * `<prefix>.nodes.tsv`  — one row per rep INCLUDING isolated ones. Without it "no edge" and "not a
///     node" are indistinguishable. `aln_len` is `len(refine_copy_seq(...))`, i.e. the PAF `qlen`/`tlen`
///     that became the coverage DENOMINATOR — under `--homology-genomic-span` that is the fetched genomic
///     span, not `rep.seq.len()`, and using the wrong one silently reproduces a different rule.
///   * `<prefix>.edges.tsv` — one row per RETURNED edge, with both node keys, both intervals, the exemplar
///     identity and coverage, the tier those two numbers came from, and the full contributing TierMask.
///   * `<prefix>.params.tsv` — every knob that can silently move the edge rule, at its EFFECTIVE value.
///     A clean diff means nothing if both sides were quietly mis-set; this block is what makes it mean
///     something.
///
/// `node_key` is `L~{chrom}_{start}_{end}`, matching `rustlib.canon_node`. `DenovoTranscript.start/end`
/// are already 0-based half-open, so NO arithmetic is applied here — the `-1` inside `canon_node` exists
/// only to convert samtools' 1-based inclusive region strings.
///
/// ⚠ Silent in `RUSTLE_SHARED_EXON=1` mode, which returns from `homology_edges_all_reps_pooled` before any
/// of this. An empty dump there is that, not a zero-edge finding.
///
/// ⚠ This function is called once per invocation of `homology_edges_all_reps_pooled`. If a run calls it
/// more than once the later calls get a `.call<N>` infix rather than overwriting, and each path is echoed
/// to stderr.
/// Lower median of a length list (deterministic for even n). `None` on an empty list.
fn median_len(lens: &[usize]) -> Option<usize> {
    if lens.is_empty() {
        return None;
    }
    let mut v = lens.to_vec();
    v.sort_unstable();
    Some(v[(v.len() - 1) / 2])
}

/// **D2 — the coverage clause `c` is SCALE-DEPENDENT.** `c` is a fraction of `min(|u|,|v|)`, so the
/// same 0.50 demands a completely different number of BASES depending on what sequence each node
/// contributes. Measured on one fixed 61-node set with the identical rule and aligner, the median
/// `0.50·min(len)` demand is **8,810 bp on the genomic span, 4,658 bp on pooled exons and 895 bp on the
/// shipped rep transcript** — and the pair-bite tracks it (8.30% / 28.31% / 30.81%).
///
/// The rule was deliberately LEFT ALONE (see `docs/seeded_family_definition.md` §1a: the absolute-bp
/// alternative was swept and refuted, and genomic-span-by-default moves the residual into locus
/// boundaries). What is fixed instead is that the substrate is no longer an UNDECLARED free parameter:
/// this returns the run's own absolute number so `<prefix>.params.tsv` can carry it.
///
/// `lens` must be the lengths of the sequences actually aligned, not `rep.seq.len()` — under
/// `--homology-genomic-span` those differ by ~7.5× (median spliced `seqlen/span` 0.1292).
fn coverage_floor_bp_demand(min_coverage: f64, lens: &[usize]) -> Option<u64> {
    let m = median_len(lens)?;
    Some((min_coverage.max(0.0) * m as f64).round() as u64)
}

#[allow(clippy::too_many_arguments)]
fn write_er_edge_dump(
    prefix: &str,
    reps: &[DenovoTranscript],
    seqs: &[Vec<u8>],
    set: &BTreeSet<(usize, usize)>,
    prov: &BTreeMap<(usize, usize), TierMask>,
    metrics: &BTreeMap<(usize, usize), (f64, f64, TierMask)>,
    // §6az additive disclosure, keyed exactly like `metrics`. Consumed here and nowhere else.
    flanks: &BTreeMap<(usize, usize), ErEdgeFlank>,
    params: &RefineParams,
    genomic_span_active: bool,
) {
    use std::io::Write;
    static DUMP_SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let n = DUMP_SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let stem = if n == 0 { prefix.to_string() } else { format!("{prefix}.call{}", n + 1) };
    let key = |r: &DenovoTranscript| format!("L~{}_{}_{}", r.chrom, r.start, r.end);

    // --- nodes ---------------------------------------------------------------------------------
    match std::fs::File::create(format!("{stem}.nodes.tsv")) {
        Ok(mut fh) => {
            let _ = writeln!(
                fh,
                "idx\tnode_key\tchrom\tstart\tend\tstrand\tn_exon\tn_reads\tdistinguishing_uniq\t\
                 core_bp\tstub\taln_len\texon_sum_len\tdegree\texons\texon_bp"
            );
            for (i, r) in reps.iter().enumerate() {
                let deg = set.iter().filter(|&&(a, b)| a == i || b == i).count();
                // The rep's EXON blocks, not just their count. Without them every locational query
                // against this table has to use the rep's genomic SPAN -- and spans are 90.83% intron
                // by bp (614,224,265 bp of span vs 56,345,279 bp of exon-sum across the 17,924 reps of
                // the 2026-08-20 catalog; median per-rep exonic fraction 0.2232). That inflates
                // "an interval lies inside a rep" ~11x over the transcribed footprint and makes
                // "correctly rejected because it is intronic" indistinguishable from "real transcript
                // sequence absent from every node" -- the ambiguity that invalidated the 2026-08-21
                // SD-gap audit. `copies.tsv` has carried exons since it was written; the 15,905 reps
                // that join NO family are exactly the ones that had no other way to be seen.
                let exons = crate::vg_family::catalog_input::exon_blocks_str(r.start, r.end, &r.introns);
                // Emitted alongside `exon_sum_len` (= the spliced seq length) ON PURPOSE: the two are
                // computed by different paths and must agree. A mismatch means the malformed-chain
                // guard dropped a block, which is a data defect the consumer should see rather than a
                // silent truncation of the locus.
                let exon_bp: u64 = exons
                    .split(',')
                    .filter(|b| !b.is_empty())
                    .filter_map(|b| b.split_once('-'))
                    .filter_map(|(a, z)| Some(z.parse::<u64>().ok()?.saturating_sub(a.parse::<u64>().ok()?)))
                    .sum();
                let _ = writeln!(
                    fh,
                    "{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{deg}\t{exons}\t{exon_bp}",
                    key(r),
                    r.chrom,
                    r.start,
                    r.end,
                    r.strand,
                    r.introns.len() + 1,
                    r.n_reads,
                    r.distinguishing_uniq,
                    r.core_bp,
                    r.stub,
                    seqs.get(i).map(|s| s.len()).unwrap_or(0),
                    r.seq.len(),
                );
            }
            eprintln!("[er-dump] {} nodes -> {stem}.nodes.tsv", reps.len());
        }
        Err(e) => eprintln!("[er-dump] FAILED to write {stem}.nodes.tsv: {e}"),
    }

    // --- edges ---------------------------------------------------------------------------------
    match std::fs::File::create(format!("{stem}.edges.tsv")) {
        Ok(mut fh) => {
            let _ = writeln!(
                fh,
                "rep_i\trep_j\tnode_key_i\tnode_key_j\tchrom_i\tstart_i\tend_i\tstrand_i\t\
                 chrom_j\tstart_j\tend_j\tstrand_j\tidentity\tcoverage\tcov_longer\tunaln_i\tunaln_j\t\
                 metric_tier\ttiers"
            );
            for &(i, j) in set.iter() {
                let (ri, rj) = (&reps[i], &reps[j]);
                let (ident, cov, mtier) = match metrics.get(&(i, j)) {
                    // {:.6} on both sides of the diff; the raw f64 is compared by both engines with `>=`
                    // against the same floors, so six places is well past what can differ.
                    Some(&(id, cv, t)) => (format!("{id:.6}"), format!("{cv:.6}"), tier_names(t)),
                    // Protein-tier edges carry no nucleotide identity/coverage at all.
                    None => ("NA".to_string(), "NA".to_string(), "NA".to_string()),
                };
                // §6az DISCLOSURE. `coverage` above is the aligned fraction of the SHORTER sequence;
                // `cov_longer` is the same quantity on the LONGER one, and `unaln_i`/`unaln_j` are each
                // rep's total 5'+3' clipped flank in bp — all three off the SAME exemplar record.
                // `NA` wherever no single record certified the edge (protein tier; opt-in summed
                // coverage), never a substituted value.
                let (cov_longer, unaln_i, unaln_j) = match flanks.get(&(i, j)) {
                    Some(f) => (
                        format!("{:.6}", f.cov_longer),
                        f.unaln_i.to_string(),
                        f.unaln_j.to_string(),
                    ),
                    None => ("NA".to_string(), "NA".to_string(), "NA".to_string()),
                };
                let _ = writeln!(
                    fh,
                    "{i}\t{j}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{ident}\t{cov}\t\
                     {cov_longer}\t{unaln_i}\t{unaln_j}\t{mtier}\t{}",
                    key(ri),
                    key(rj),
                    ri.chrom,
                    ri.start,
                    ri.end,
                    ri.strand,
                    rj.chrom,
                    rj.start,
                    rj.end,
                    rj.strand,
                    tier_names(prov.get(&(i, j)).copied().unwrap_or(0)),
                );
            }
            eprintln!("[er-dump] {} edges -> {stem}.edges.tsv", set.len());
        }
        Err(e) => eprintln!("[er-dump] FAILED to write {stem}.edges.tsv: {e}"),
    }

    // --- params --------------------------------------------------------------------------------
    // The PAF stems are built from the RAW prefix (they are written by `nucleotide_edges_scored`, which
    // does not know about this function's `.call<N>` de-collision suffix).
    let stem_hint = prefix;
    let env = |k: &str| std::env::var(k).unwrap_or_else(|_| "<unset>".to_string());
    // ⭐ X.4 — THE RULE ROWS ARE NO LONGER TYPED HERE. They come from the shared `er_rule_rows`, the same
    // function `refine_families_exon_sum` calls, and are ALSO written on their own as `<stem>.rule.tsv`
    // so the O1/O2 parity question is a `diff` of two data-free files. Everything below `rows` is
    // input-dependent (counts, lengths, paths) and must NOT be used for that diff.
    let rule = er_rule_rows(
        params,
        // ⚠ O-3 — `NotPresent`, and this is a REAL ASYMMETRY, not a formality: this site SWAPS its
        // substrate (`homology_genomic_span`, default OFF) where refine UNIONS a second one in (default
        // ON). Recording it is what makes `diff <prefix>.rule.tsv <prefix>.refine.rule.tsv` able to see
        // the difference at all; it printed an empty diff across it before.
        &ErRuleSite {
            substrate_genomic: genomic_span_active,
            core_lens_supplied: true,
            genomic_tier: GenomicTier::NotPresent,
        },
    );
    write_kv_tsv(&format!("{stem}.rule.tsv"), &rule);
    let mut rows: Vec<(String, String)> = vec![
        ("site".into(), "homology_edges_all_reps_pooled".into()),
        ("n_reps".into(), reps.len().to_string()),
        ("n_edges".into(), set.len().to_string()),
        ("min_identity_asm20".into(), format!("{:.6}", params.min_identity)),
        ("sensitive_identity".into(), format!("{:.6}", params.sensitive_identity)),
        ("nucleotide_sensitive".into(), params.nucleotide_sensitive.to_string()),
        ("protein_tail".into(), params.protein_tail.to_string()),
        ("homology_genomic_span".into(), params.homology_genomic_span.to_string()),
        ("genomic_span_active".into(), genomic_span_active.to_string()),
        ("threads".into(), params.threads.to_string()),
        // Derived from ER_TIER_FLAGS, never re-typed: this row used to be a fourth hardcoded copy of the
        // tier and could therefore describe a command the binary did not run.
        (
            "mm_args_sensitive".into(),
            er_tier_cmdline(&params.minimap2, params.threads, ER_SENSITIVE_SEED),
        ),
        // THE ONLY PAFs THAT BUILT THIS EDGE SET. A run also drops `<prefix>.refine.*.paf` from the
        // downstream per-family refinement; diffing against one of those compares nothing.
        ("paf_glob_for_this_edge_set".into(), format!("{stem_hint}.er.*.paf")),
        // D2: the coverage clause is a FRACTION, so the substrate silently sets its absolute demand.
        // Record both so no catalog can be quoted without its scale (docs/seeded_family_definition.md §1a).
        ("substrate_median_len_bp".into(), {
            let lens: Vec<usize> = seqs.iter().map(|s| s.len()).collect();
            median_len(&lens).map(|m| m.to_string()).unwrap_or_else(|| "NA".into())
        }),
        ("coverage_floor_median_bp_demand".into(), {
            let lens: Vec<usize> = seqs.iter().map(|s| s.len()).collect();
            coverage_floor_bp_demand(params.min_coverage, &lens)
                .map(|b| b.to_string())
                .unwrap_or_else(|| "NA".into())
        }),
        ("env.RUSTLE_ER_SENSITIVE_ONLY".into(), env("RUSTLE_ER_SENSITIVE_ONLY")),
        ("env.RUSTLE_ER_SUM_COVERAGE".into(), env("RUSTLE_ER_SUM_COVERAGE")),
        ("env.RUSTLE_ER_SUM_MIN_BLOCK".into(), env("RUSTLE_ER_SUM_MIN_BLOCK")),
        ("env.RUSTLE_ER_CORE_COVERAGE".into(), env("RUSTLE_ER_CORE_COVERAGE")),
        ("env.RUSTLE_ER_CORE_DEPTH".into(), env("RUSTLE_ER_CORE_DEPTH")),
        ("env.RUSTLE_ER_COVERAGE_LONGER".into(), env("RUSTLE_ER_COVERAGE_LONGER")),
        ("env.RUSTLE_ER_NO_STUB_EDGES".into(), env("RUSTLE_ER_NO_STUB_EDGES")),
        ("env.RUSTLE_SHARED_EXON".into(), env("RUSTLE_SHARED_EXON")),
        ("env.RUSTLE_SENSITIVE_IDENTITY".into(), env("RUSTLE_SENSITIVE_IDENTITY")),
        ("env.RUSTLE_GENOME_MIN_COVERAGE".into(), env("RUSTLE_GENOME_MIN_COVERAGE")),
        ("env.RUSTLE_GENOME_GAMMA".into(), env("RUSTLE_GENOME_GAMMA")),
        ("env.RUSTLE_LOCUS_EXON_UNION".into(), env("RUSTLE_LOCUS_EXON_UNION")),
        ("env.RUSTLE_COTHREAD_REP".into(), env("RUSTLE_COTHREAD_REP")),
        // Third member of the set above: like EXON_UNION and COTHREAD_REP it changes what a rep IS --
        // here the `strand` column AND the stored sequence bytes of every unspliced rep. Without this row
        // an OFF and an ON catalog have byte-identical params.tsv, defeating this file's stated purpose.
        ("env.RUSTLE_READ_STRAND".into(), env("RUSTLE_READ_STRAND")),
        // PRE-EXISTING M2 HOLES, closed here. `RUSTLE_READ_STRAND_MARGIN` is READ_STRAND's abstention
        // threshold (`denovo_assemble.rs`), so two READ_STRAND arms at different margins were previously
        // indistinguishable from their params.tsv. `RUSTLE_SPLICED_REP` changes rep selection
        // (`family_detect::locus_reps`) AND the unspliced-strand collapse clause (`collapse_parent`), yet
        // an ON and an OFF pair of SPLICED_REP arms had byte-identical params.tsv until this row existed.
        ("env.RUSTLE_READ_STRAND_MARGIN".into(), env("RUSTLE_READ_STRAND_MARGIN")),
        ("env.RUSTLE_SPLICED_REP".into(), env("RUSTLE_SPLICED_REP")),
        // Also an opt-in that changes the EDGE RULE, so a catalog must be able to name it. Omitting
        // this row is the M2 defect from the read-strand work, repeated: without it an ON and an OFF
        // arm have byte-identical params.tsv and a null result cannot be distinguished from a flag
        // that never reached the process.
        ("env.RUSTLE_ER_GUARD_SCOPED".into(), env("RUSTLE_ER_GUARD_SCOPED")),
    ];
    // The rule rows land LAST, so `params.tsv` remains a strict superset of `rule.tsv` on both sides.
    rows.extend(rule);
    write_kv_tsv(&format!("{stem}.params.tsv"), &rows);
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
    let rc = crate::vg_family::seq_utils::reverse_complement(seq);
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
/// paralog copies occupy DISJOINT genomic spans, so two copies overlapping on the same chromosome are
/// candidates for the same locus. SAME-strand overlap collapses iff NO read distinguishes the pair — the
/// χ(H) read-conflict criterion (`read_conflict::reads_distinguish`), restricted to co-located copies: a
/// copy whose `distinguishing_uniq` (unique-mapper support) clears `min_reads` is read-evidenced as its own
/// locus (a gene plus its own nested fragment / a second isoform is the K=0 case homology alignment cannot
/// tell apart from a genuine collapsed paralog — reads can). This replaces the unconditional coordinate
/// collapse that produced the advisor-flagged "distinguishable-but-merged" cases. OPPOSITE-strand overlap is
/// the sense/antisense case: collapse ONLY when one copy is a clear read-minority (a few antisense reads on
/// a real transcript = a strand artifact, e.g. GWFAM99 = 666 `+` reads vs 3 `-` reads at one locus, a
/// sense/antisense mis-split); a BALANCED overlapping antisense pair is two genuine loci and is kept —
/// UNCHANGED by this fix. `min_reads` is the conflict-graph noise floor (`cfg.conflict.min_reads`), not a
/// new threshold. The representative of a locus is the MOST-supported copy (most reads — the real one, not
/// the minority artifact — then widest span).
fn distinct_locus_reps(copies: Vec<DenovoTranscript>, min_reads: usize) -> Vec<DenovoTranscript> {
    distinct_locus_reps_grouped(copies, min_reads).into_iter().map(|(t, _)| t).collect()
}

/// As `distinct_locus_reps`, but ALSO returns, for each emitted locus, the indices into the input
/// `copies` that merged into it. Identical merge behaviour — `distinct_locus_reps` is a thin wrapper.
///
/// EXISTS FOR THE λ CERTIFICATE. λ has to be computed on the graph whose nodes ARE THE EMITTED COPIES,
/// not on the pre-merge block: this function is exactly where the node set changes (co-located copies
/// collapse into one locus), so a λ measured before it would describe a different object than the one
/// the catalog prints. Computing a structural statistic on the wrong node set is the documented trap
/// "never judge a change to what a NODE IS on node-level metrics".
fn distinct_locus_reps_grouped(
    copies: Vec<DenovoTranscript>,
    min_reads: usize,
) -> Vec<(DenovoTranscript, Vec<usize>)> {
    let n = copies.len();
    // AUDIT ONLY (`RUSTLE_LOCUS_AUDIT=1`, default silent): log every co-located pair this merge examines
    // and its verdict, so "how many merges does the >=2-distinct-loci certificate perform, and at what
    // overlap fraction" is measurable without re-deriving the block structure offline.
    let audit = std::env::var("RUSTLE_LOCUS_AUDIT").map(|v| v != "0" && !v.is_empty()).unwrap_or(false);
    // A/B KNOB (`RUSTLE_LOCUS_MERGE_MIN_COVER`, default 0.0 = today's "any overlap"): additionally require
    // the overlap to cover >= this fraction of the SHORTER copy before the co-located merge may fire — i.e.
    // the same containment bar `collapse_parent` (the candidate->locus merge) uses. Exists to MEASURE the
    // advisor's proposed harmonisation of the two overlap rules; not a shipped default.
    let merge_min_cover: f64 = std::env::var("RUSTLE_LOCUS_MERGE_MIN_COVER")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let mut parent: Vec<usize> = (0..n).collect();
    // Counts pairs whose MERGE/KEEP verdict came from `reads_distinguish` — O2's χ(H) predicate — i.e.
    // the sole place O1's node set consults read evidence. Reported unconditionally below so the
    // O1 ⊥ O2 exception can never be silently load-bearing. See the block comment at the branch.
    let mut read_leg_decisions: usize = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = (&copies[i], &copies[j]);
            let same_pos = a.chrom == b.chrom && a.end.min(b.end) > a.start.max(b.start);
            if !same_pos {
                continue;
            }
            let cover = {
                let ov = a.end.min(b.end) - a.start.max(b.start);
                let minlen = (a.end - a.start).min(b.end - b.start).max(1);
                ov as f64 / minlen as f64
            };
            if cover < merge_min_cover {
                if audit {
                    eprintln!(
                        "[locus-audit]\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\tBELOW_COVER",
                        a.chrom, a.start, a.end, a.strand, a.n_reads, a.distinguishing_uniq,
                        b.start, b.end, b.strand, b.n_reads, b.distinguishing_uniq,
                        a.end.min(b.end) - a.start.max(b.start), cover
                    );
                }
                continue;
            }
            // A JUNCTION-LESS rep carries no strand evidence and no independent splice structure:
            // `DenovoTranscript::strand` is assigned "from the gate's canonical-junction classification",
            // and a copy with no introns has no junctions to classify, so its `'+'` is a DEFAULT LABEL, not
            // a call. When such a rep's span is additionally CONTAINED in an overlapping copy's span, it is
            // that copy's own unspliced read fraction — not a paralog — and it must collapse.
            //
            // Measured (2026-08-08, 7 strict single-copy control genes): every spurious "2-copy family" at
            // a single-copy gene was exactly this shape — (spliced multi-exon rep) + (single-exon rep at the
            // same start), edge identity 1.0000 because it is literally the same sequence, and the
            // single-exon rep's read count equal to the count of reads with no `N` in CIGAR. Those reads sit
            // on the SAME BAM strand, so nothing about the pair is antisense; the fake strand difference
            // merely routed it to the `ANTISENSE_MINORITY_DENOM` branch, where ratios of 0.145-0.571 cleared
            // the 1/10 artifact bar and the pair survived as two "loci". 4 of 7 controls emitted families
            // this way, and the same shape inflated real families' copy counts (MAGEA 13 records at 11 loci,
            // GSTM 5 at 4).
            //
            // THE RULE: exactly ONE side junction-less, the other spliced. A rep WITH junctions carries
            // positive evidence of being an independent transcript unit (its splice structure); an
            // overlapping rep WITHOUT junctions carries none, and cannot be told apart from that
            // transcript's own unspliced reads. When BOTH lack junctions neither can claim priority, so
            // the read-conflict criterion below decides — which is what preserves the advisor-flagged
            // "distinguishable-but-merged" cases.
            //
            // ⚠ Span CONTAINMENT was tried first and is the WRONG SHAPE: it fixed only 1 of the 4 measured
            // artifacts, because the junction-less rep's ends are read-derived and routinely overhang by a
            // base or two — HMBS by 1 bp (127034654-127036246 vs ...245), TFRC by 2 bp
            // (207807896 vs 207807898), and TBP by 5 kb. Overlap + the junction asymmetry is the invariant;
            // containment is not.
            //
            // PARAMETER-FREE (a junction count of zero, and the overlap already required by `same_pos`), so
            // this restores the threshold-free distinct-LOCUS guarantee this function's doc comment claims.
            // It cannot touch a genuine paralog pair (disjoint spans never reach here), a balanced
            // sense/antisense pair (both spliced), or two junction-less copies (read rule still decides).
            let unspliced_fragment = a.introns.is_empty() != b.introns.is_empty();
            // ⚠ O1 ⊥ O2 — THE ONE PLACE O1's NODE SET DEPENDS ON READS, STATED RATHER THAN ASSUMED.
            //
            // `reads_distinguish` is O2's χ(H) edge predicate (`read_conflict.rs:253`), so the branch
            // below is a point where READ EVIDENCE decides what a NODE IS, and therefore where O1's
            // "membership by SEQUENCE alone" does not hold literally. It reaches here from the DEFAULT
            // catalog path (`detect_homology_catalog_genome_wide`), not from behind `--refine`.
            //
            // IT IS NOT REMOVABLE, AND THAT IS THE POINT. The branch is entered only by a co-located,
            // SAME-STRAND pair that the junction-asymmetry rule above did not settle — i.e. one whose
            // two sides are both spliced or both junction-less. In the junction-less case there is no
            // splice structure to compare and the genomic sequence is shared by construction (the spans
            // overlap), so NOTHING in the sequence separates the pair: this is the true K=0 case, and
            // the unique-mapper count is the only signal that exists. Both advisor-flagged
            // "distinguishable-but-merged" tests below (`..._keeps_distinguishable_colocated_copies_
            // separate`, `..._two_junctionless_copies_still_decided_by_reads`) are exactly this shape —
            // their `rep()` helper builds `introns: vec![]`.
            //
            // MEASURED 2026-08-13 (`RUSTLE_LOCUS_AUDIT=1`, gorilla): this branch decided **0 of 109**
            // co-located pairs — 15 on the 25-region control panel and 94 on the 19-region family panel.
            // Every pair was settled by the junction-asymmetry rule above (104) or the antisense-ratio
            // rule below (5). So the dependency is real in the code and inert on the data measured.
            // ⚠ Do NOT "fix" this by comparing intron chains on the both-spliced side: two overlapping
            // same-strand reps with different chains are usually ISOFORMS OF ONE GENE, and splitting
            // them would over-split loci. There is no substrate on record where the branch fires, so
            // such a change cannot be validated — it would be judged on node-level metrics, which has
            // failed end-to-end three times in this project.
            //
            // The counter below makes the exception VISIBLE instead of silent: if the read leg ever
            // decides anything, the run says so on stderr, unconditionally.
            let collapse = if unspliced_fragment {
                true
            } else if a.strand == b.strand {
                read_leg_decisions += 1;
                // χ(H): collapse only when NO read distinguishes them (true K=0). The PSV/junction leg is
                // not wired here (no per-copy PSV/junction evidence reaches this merge point) — the
                // unique-mapper leg alone is the fix for the flagged over-merge cases. `min_reads.max(1)`:
                // distinguishability requires at least 1 unique read, so a misconfigured
                // `RUSTLE_CONFLICT_MIN_READS=0` can never flip `reads_distinguish(0, 0, false, 0)` to `true`
                // and invert this guard into "nothing ever merges" (see `ConflictParams::from_env`).
                !reads_distinguish(a.distinguishing_uniq, b.distinguishing_uniq, /*shared_psv_or_junction=*/false, min_reads.max(1))
            } else {
                let (lo, hi) = (a.n_reads.min(b.n_reads), a.n_reads.max(b.n_reads));
                lo.saturating_mul(ANTISENSE_MINORITY_DENOM) < hi // minority antisense = artifact
            };
            if audit {
                let ov = a.end.min(b.end) - a.start.max(b.start);
                let minlen = (a.end - a.start).min(b.end - b.start).max(1);
                eprintln!(
                    "[locus-audit]\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
                    a.chrom, a.start, a.end, a.strand, a.n_reads, a.distinguishing_uniq,
                    b.start, b.end, b.strand, b.n_reads, b.distinguishing_uniq,
                    ov, ov as f64 / minlen as f64,
                    // distinguish the unspliced-fragment collapse so "how many merges came from the
                    // junction-less rule" is countable straight out of the audit log.
                    match (collapse, unspliced_fragment) {
                        (true, true) => "MERGE_UNSPLICED",
                        (true, false) => "MERGE",
                        (false, _) => "KEEP",
                    }
                );
            }
            if collapse {
                uf_union(&mut parent, i, j);
            }
        }
    }
    // ⚠ O1 ⊥ O2 EXCEPTION REPORT — unconditional, not behind `RUSTLE_LOCUS_AUDIT`. Silence here is the
    // claim "this catalog's node set is a function of SEQUENCE alone"; a non-zero count is the claim
    // "N co-located same-strand pairs were decided by O2's χ(H) read predicate", which must then be
    // disclosed wherever the catalog's numbers are quoted. Measured 0 on every substrate to date.
    if read_leg_decisions > 0 {
        eprintln!(
            "[o1-perp-o2] WARNING: {read_leg_decisions} co-located same-strand pair(s) were decided by \
             reads_distinguish (O2's chi(H) predicate), not by sequence. O1's node set is NOT a \
             function of sequence alone for this run — disclose it with any number derived from it."
        );
    }
    // representative per locus = MOST reads (the real copy, not the minority artifact), then widest span.
    let key = |t: &DenovoTranscript| (t.n_reads, t.end - t.start);
    let roots: Vec<usize> = (0..n).map(|i| uf_find(&mut parent, i)).collect();
    let mut by_locus: std::collections::BTreeMap<usize, usize> = std::collections::BTreeMap::new(); // root -> rep idx
    let mut members: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
    for i in 0..n {
        let r = roots[i];
        let rep = by_locus.entry(r).or_insert(i);
        if key(&copies[i]) > key(&copies[*rep]) {
            *rep = i;
        }
        members.entry(r).or_default().push(i);
    }
    // Emitted order is by representative index, exactly as before the grouped variant existed.
    let mut out: Vec<(usize, Vec<usize>)> = by_locus
        .into_iter()
        .map(|(root, rep)| (rep, members.remove(&root).unwrap_or_default()))
        .collect();
    out.sort_by_key(|&(rep, _)| rep);
    out.into_iter().map(|(rep, mem)| (copies[rep].clone(), mem)).collect()
}

#[cfg(test)]
mod tests {
    use super::super::copy_split::AlignedRead;
    use super::super::family_detect::collapse_loci;
    use super::*;

    // ── B1: THE SHIPPED TIER IS CENTRALISED, AND STAYING CENTRALISED IS ENFORCED ─────────────────
    //
    // The value below is the SHIPPED value, transcribed from the pre-fix binary. If this assertion
    // ever needs editing, the tier changed and every O1 number in
    // docs/seeded_family_definition.md must be recomputed and re-labelled.
    #[test]
    fn er_tier_value_is_the_shipped_tier_and_is_unchanged() {
        assert_eq!(ER_TIER_FLAGS, &["-c", "-X", "--no-long-join"]);
        assert_eq!(ER_SENSITIVE_SEED, &["-k", "11", "-w", "5"]);
        // The exact argv the shipped binary ran on the sensitive tier, threads=4.
        assert_eq!(
            er_tier_cmdline("minimap2", 4, ER_SENSITIVE_SEED),
            "minimap2 -c -X --no-long-join -t 4 -k 11 -w 5"
        );
        // ...and on the asm20 seeding tier used by refine.
        let seed: Vec<&str> = vec!["-x", "asm20"];
        assert_eq!(
            er_tier_cmdline("minimap2", 4, &seed),
            "minimap2 -c -X --no-long-join -t 4 -x asm20"
        );
        // `-t 0` is not a legal minimap2 thread count; one of the two former hardcoded sites clamped
        // and the other did not. The centralised builder clamps, always.
        assert!(er_tier_cmdline("minimap2", 0, &[]).contains("-t 1"));
        // The COMMAND and the CMDLINE must not be able to disagree -- a .args dump that describes a
        // command the binary did not run is worse than no dump.
        let cmd = er_tier_command("minimap2", 4, ER_SENSITIVE_SEED);
        let argv: Vec<String> =
            cmd.get_args().map(|a| a.to_string_lossy().into_owned()).collect();
        assert_eq!(
            format!("{} {}", cmd.get_program().to_string_lossy(), argv.join(" ")),
            er_tier_cmdline("minimap2", 4, ER_SENSITIVE_SEED)
        );
    }

    /// ANTI-DRIFT GUARD. The tier was hardcoded at four sites and that is precisely how it drifted
    /// away from the panel scripts. This test fails the moment a fifth copy appears: it re-reads this
    /// source file and requires the flag to exist as a Rust string literal EXACTLY ONCE.
    #[test]
    fn er_tier_flags_appear_as_a_literal_exactly_once_in_the_source() {
        let whole = include_str!("denovo_pipeline.rs");
        // Count only SHIPPING code. The test module legitimately re-types the tier once, to PIN its
        // value (`er_tier_value_is_the_shipped_tier_and_is_unchanged`); that copy is the point.
        let src = whole.split("\n#[cfg(test)]\nmod tests {").next().expect("test module marker");
        assert!(src.len() < whole.len(), "test-module split failed; the guard would scan itself");
        // Assembled at runtime so this test's own text is not the thing it counts.
        let needle = format!("\"--no-{}\"", "long-join");
        let n = src.matches(needle.as_str()).count();
        assert_eq!(
            n, 1,
            "the E_r tier flag appears as a string literal {n} times; it must appear ONCE, inside \
             ER_TIER_FLAGS. A second copy is how B1 happened. Use er_tier_command()."
        );
        let k11 = src.matches("\"-k\", \"11\", \"-w\", \"5\"").count();
        assert_eq!(k11, 1, "the sensitive seed is re-typed {k11} times; use ER_SENSITIVE_SEED.");
    }

    // ── M1: THE COVERAGE STATISTIC IS A FRACTION ────────────────────────────────────────────────
    //
    // Pure-arithmetic guard on the FORM. The rule is: the numerator's axis follows the denominator.
    // Written against the same PAF fields the parser reads, so it pins the formula, not a wrapper.
    #[test]
    fn coverage_numerator_axis_follows_the_denominator() {
        // (ql, qs, qe, tl, ts, te)
        fn cov(ql: f64, qs: f64, qe: f64, tl: f64, ts: f64, te: f64, longer: bool) -> f64 {
            let side_is_query = if longer { ql >= tl } else { ql <= tl };
            let denom = if longer { ql.max(tl).max(1.0) } else { ql.min(tl).max(1.0) };
            let aln = if side_is_query { qe - qs } else { te - ts };
            aln / denom
        }
        // THE DEFECT, reproduced. -X/--dual=no put the LONGER sequence on the query axis: a 10,000 bp
        // query aligning 9,000 bp against a 2,000 bp target. Old form = 9000/2000 = 4.5, a "fraction"
        // of 450%. New form reads the target axis: 1800/2000 = 0.90.
        let old = (9000.0 - 0.0) / 2000.0_f64.min(10000.0);
        assert!(old > 1.0, "the pre-fix form must be reproducible as >1 to prove it was not a fraction");
        let new = cov(10000.0, 0.0, 9000.0, 2000.0, 100.0, 1900.0, false);
        assert!((new - 0.90).abs() < 1e-9, "got {new}");
        // A fraction can never exceed 1: the aligned span on a sequence is bounded by its length.
        assert!(new <= 1.0);
        // Query IS the shorter sequence -> the query axis is used, i.e. the historical behaviour is
        // preserved on exactly the records where it was already correct.
        let q_shorter = cov(2000.0, 100.0, 1900.0, 10000.0, 0.0, 9000.0, false);
        assert!((q_shorter - 0.90).abs() < 1e-9, "got {q_shorter}");
        // COVERAGE_LONGER demands BOTH members be covered, so it takes the longer sequence's axis.
        let longer = cov(10000.0, 0.0, 9000.0, 2000.0, 100.0, 1900.0, true);
        assert!((longer - 0.90).abs() < 1e-9, "got {longer}");
        // ...and the same pair with the roles swapped gives the SAME number under both settings.
        // Symmetry under q/t exchange is the property the old form lacked.
        for lg in [false, true] {
            let a = cov(10000.0, 0.0, 9000.0, 2000.0, 100.0, 1900.0, lg);
            let b = cov(2000.0, 100.0, 1900.0, 10000.0, 0.0, 9000.0, lg);
            assert!((a - b).abs() < 1e-9, "coverage must be symmetric under q/t swap (longer={lg})");
        }
    }

    /// End-to-end on the real parser: two sequences of DIFFERENT length, so the denominator axis is
    /// actually exercised. Every other aligner test in this file uses EQUAL-length sequences, which
    /// makes the denominator untestable by construction.
    #[test]
    fn nucleotide_edges_coverage_is_bounded_by_one_on_unequal_lengths() {
        // A 6 kb sequence and a 2 kb exact SUBSEQUENCE of it. The short one is fully covered (1.00);
        // the long one is one third covered. Under the shipped -X the pair is emitted once, in an
        // orientation we do not control -- which is the whole point.
        let mut rng: u64 = 0x5eed_1234;
        let long: Vec<u8> = (0..6000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[((rng >> 33) % 4) as usize]
            })
            .collect();
        let short = long[2000..4000].to_vec();
        let p = RefineParams { threads: 1, ..Default::default() };
        let scored =
            nucleotide_edges_scored(&[long, short], ER_SENSITIVE_SEED, 0.60, 0.50, None, &p, None, None)
                .expect("alignment");
        assert_eq!(scored.len(), 1, "the pair must produce exactly one edge");
        let cov = scored[0].3;
        assert!(
            cov <= 1.0 + 1e-9,
            "coverage {cov} exceeds 1.0 -- the statistic is not a fraction (M1 has regressed)"
        );
        // Coverage-of-the-SHORTER of an exact subsequence is ~1.0 whichever axis minimap2 chose.
        assert!(cov > 0.90, "coverage {cov} -- an exact 2 kb subsequence must cover the shorter fully");
    }

    /// §6az DISCLOSURE ARITHMETIC. Pins `cov_longer` and the two flank lengths on a synthetic record
    /// whose overhang is known exactly, and pins the inequality that makes the disclosure worth printing.
    ///
    /// The record is §6ay's synthetic control made concrete: a 5,000 bp sequence that aligns end-to-end
    /// into a 7,000 bp one which carries 1,000 bp of extra sequence at EACH end. Coverage-of-the-shorter
    /// reports 1.000 — a perfect number — while 2,000 bp of the longer member is outside the alignment
    /// and no shipped column said so.
    #[test]
    fn er_edge_flank_discloses_the_longer_members_overhang() {
        // query = the SHORTER member and rep i; target = the LONGER member and rep j.
        let (ql, qs, qe) = (5000.0, 0.0, 5000.0);
        let (tl, ts, te) = (7000.0, 1000.0, 6000.0);
        let f = er_edge_flank(ql, qs, qe, tl, ts, te, /* q_is_i */ true);

        // The shipped statistic, recomputed here rather than assumed: aligned span on the SHORTER axis
        // over the SHORTER length. This is the number that looks perfect.
        let coverage = (qe - qs) / ql;
        assert!((coverage - 1.0).abs() < 1e-12, "coverage {coverage} — the control must look perfect");

        // 5,000 aligned bases on the LONGER axis over 7,000 bp.
        assert!(
            (f.cov_longer - 5000.0 / 7000.0).abs() < 1e-12,
            "cov_longer {} != 5000/7000",
            f.cov_longer
        );
        assert_eq!(f.unaln_i, 0, "the shorter member aligns end to end");
        assert_eq!(f.unaln_j, 2000, "1,000 bp hangs off EACH end of the longer member");

        // Whenever the lengths differ and the aligned span is the same on both axes, the larger
        // denominator makes `cov_longer` STRICTLY smaller. That is the whole content of the concession.
        assert!(
            f.cov_longer < coverage,
            "cov_longer {} must be < coverage {coverage} when the lengths differ",
            f.cov_longer
        );

        // The i/j mapping is not cosmetic: under the shipped `-X` minimap2 chooses which sequence is the
        // query, so the same record with q and t swapped must move the flanks to the other column and
        // leave `cov_longer` — a property of the PAIR — exactly where it was.
        let g = er_edge_flank(tl, ts, te, ql, qs, qe, /* q_is_i */ true);
        assert!((g.cov_longer - f.cov_longer).abs() < 1e-12, "cov_longer must not depend on the axis order");
        assert_eq!((g.unaln_i, g.unaln_j), (2000, 0), "swapping q/t must swap the flanks");

        // Equal lengths: no member is "the longer", the denominators coincide, and the inequality above
        // becomes an equality rather than reversing.
        let h = er_edge_flank(4000.0, 100.0, 3900.0, 4000.0, 50.0, 3850.0, true);
        assert!((h.cov_longer - 3800.0 / 4000.0).abs() < 1e-12, "cov_longer {}", h.cov_longer);
        assert_eq!((h.unaln_i, h.unaln_j), (200, 200));
    }

    /// The orientation guard is typed to transcript-oriented RNA representatives. Two representatives
    /// that are exact reverse complements pass the historical strand-agnostic rule, but cannot describe
    /// the same transcript orientation and must not form an RNA E_r edge when the guard is armed.
    #[test]
    fn rna_forward_only_rejects_a_reverse_complement_only_edge() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            return;
        }
        let mut rng: u64 = 0x0a11_ce55_2026;
        let seq: Vec<u8> = (0..5000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[((rng >> 33) % 4) as usize]
            })
            .collect();
        let rev = crate::vg_family::seq_utils::reverse_complement(&seq);

        let permissive = RefineParams { threads: 1, ..Default::default() };
        let old = nucleotide_edges_scored(
            &[seq.clone(), rev.clone()],
            ER_SENSITIVE_SEED,
            0.60,
            0.50,
            None,
            &permissive,
            None,
            None,
        )
        .expect("strand-agnostic alignment");
        assert_eq!(old.len(), 1, "the historical rule must admit the reverse-complement pair");

        // ⚠ The SUBSTRATE arms the guard, not the flag alone — this is the RNA half of the typed rule,
        // so it must declare transcript orientation. (`forward_only_guard_is_inert_on_a_reference_\
        // oriented_substrate` covers the DNA half.)
        let guarded = RefineParams {
            threads: 1,
            require_forward_alignment: true,
            substrate: Substrate::TranscriptOriented,
            ..Default::default()
        };
        let new = nucleotide_edges_scored(
            &[seq, rev],
            ER_SENSITIVE_SEED,
            0.60,
            0.50,
            None,
            &guarded,
            None,
            None, // the scoped exemption is opt-in; this test pins the UNSCOPED guard
        )
        .expect("forward-only alignment");
        assert!(new.is_empty(), "a reverse-only match must not become an RNA family edge");
    }

    /// The other half of the typed rule: reference-oriented DNA reps must continue to admit a real
    /// inverted duplication. This exercises the same `families_from_reps_certified` grouping core used by
    /// `gw_family_catalog --from-genome`, not only the PAF parser in isolation.
    #[test]
    /// `--refine` used to emit `NA` for λ because its partition comes from a different edge set.
    /// `certificates_for_families` rebuilds the graph over exactly the rows that will be written, so
    /// the certificate describes the emitted family rather than a superseded one.
    fn certificates_are_recomputed_for_a_refined_partition() {
        // two families: a connected pair, and a single copy (degenerate)
        let mk = |name: &str, start: u64, seq: &str| DenovoTranscript {
            tid: name.to_string(),
            chrom: "chr1".to_string(),
            start,
            end: start + seq.len() as u64,
            n_reads: 5,
            seq: seq.as_bytes().to_vec(),
            ..Default::default()
        };
        let body: String = std::iter::repeat("ACGTTGCAAGGCTTACGGATCCTTAGGCAT").take(40).collect();
        let fams = vec![
            vec![mk("a", 1_000, &body), mk("b", 900_000, &body)],
            vec![mk("c", 2_000_000, &body)],
        ];
        let params = RefineParams { threads: 1, ..Default::default() };
        let certs = match certificates_for_families(&fams, &params) {
            Ok(c) => c,
            // minimap2 absent in this environment -> the wiring is what matters, skip silently
            Err(_) => return,
        };
        assert_eq!(certs.len(), fams.len(), "one certificate per emitted family, same order");
        assert_eq!(certs[0].n, 2, "the certificate counts the EMITTED rows");
        assert_eq!(certs[1].n, 1);
        assert_eq!(certs[1].lambda, 0, "a single-copy family has no cut to certify");
    }

    #[test]
    /// LOCK-IN, 2026-08-19. `RefineParams` configures BOTH substrates, so orientation-sensitive
    /// options must be inert until an entry point DECLARES the substrate. Flipping this default to
    /// `TranscriptOriented` would silently arm the RNA guard on reference-oriented DNA reps and drop
    /// inverted duplications — the mistake `genome_mode_grouping_keeps_an_inverted_duplication` caught.
    /// **If you are here to flip this default, that is the bug.**
    fn refine_params_default_is_orientation_agnostic() {
        assert!(
            !RefineParams::default().require_forward_alignment,
            "RefineParams::default() must stay orientation-agnostic — enable the forward-only guard at \
             the RNA entry points instead, or reference-oriented DNA reps lose inverted duplications"
        );
        assert_eq!(
            RefineParams::default().substrate,
            Substrate::ReferenceOriented,
            "the conservative default: a forgotten substrate declaration must cost PRECISION, never an \
             inverted duplication"
        );
    }

    #[test]
    /// The substrate, not the flag, is what arms an orientation-sensitive option. This is the whole
    /// point of the type: setting the flag on a reference-oriented substrate must be a no-op rather
    /// than a silent correctness bug.
    fn forward_only_guard_is_inert_on_a_reference_oriented_substrate() {
        let dna = RefineParams {
            require_forward_alignment: true,
            substrate: Substrate::ReferenceOriented,
            ..Default::default()
        };
        assert!(
            !dna.forward_only_active(),
            "a '-' record between REFERENCE-oriented reps is a real inverted segmental duplication; \
             the forward-only guard must not fire even when the flag is set"
        );
        let rna = RefineParams {
            require_forward_alignment: true,
            substrate: Substrate::TranscriptOriented,
            ..Default::default()
        };
        assert!(rna.forward_only_active(), "declared RNA substrate + flag set must arm the guard");
        let rna_off = RefineParams {
            require_forward_alignment: false,
            substrate: Substrate::TranscriptOriented,
            ..Default::default()
        };
        assert!(!rna_off.forward_only_active(), "--no-rna-forward-only must still disarm it");
    }

    #[test]
    fn genome_mode_grouping_keeps_an_inverted_duplication() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            return;
        }
        let mut rng: u64 = 0xd1a0_1a57;
        let seq: Vec<u8> = (0..5000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[((rng >> 33) % 4) as usize]
            })
            .collect();
        let rev = crate::vg_family::seq_utils::reverse_complement(&seq);
        let reps = vec![
            DenovoTranscript {
                tid: "dna_a".into(),
                chrom: "a".into(),
                start: 0,
                end: 5000,
                n_reads: 1,
                strand: '+',
                seq,
                ..Default::default()
            },
            DenovoTranscript {
                tid: "dna_b_inverted".into(),
                chrom: "b".into(),
                start: 0,
                end: 5000,
                n_reads: 1,
                strand: '+',
                seq: rev,
                ..Default::default()
            },
        ];
        let params = RefineParams { threads: 1, ..Default::default() };
        let (fams, _) = families_from_reps_certified(reps, &params, 0.20, 2, 0)
            .expect("DNA grouping");
        assert_eq!(fams.len(), 1, "the inverted DNA duplication must remain one family");
        assert_eq!(fams[0].len(), 2);
    }

    #[test]
    fn orientation_rule_is_explicit_in_the_certificate() {
        let site = ErRuleSite {
            substrate_genomic: false,
            core_lens_supplied: true,
            genomic_tier: GenomicTier::NotPresent,
        };
        let off = er_rule_rows(&RefineParams::default(), &site);
        // ⚠ The SUBSTRATE is what arms the guard, not the flag — setting `require_forward_alignment`
        // alone on the default (reference-oriented) substrate is deliberately a no-op, so the
        // certificate must be built on a declared RNA substrate to report "forward-only".
        let on = er_rule_rows(
            &RefineParams {
                require_forward_alignment: true,
                substrate: Substrate::TranscriptOriented,
                ..Default::default()
            },
            &site,
        );
        // and the same flag on the DNA substrate must still report "both"
        let dna = er_rule_rows(
            &RefineParams {
                require_forward_alignment: true,
                substrate: Substrate::ReferenceOriented,
                ..Default::default()
            },
            &site,
        );
        let value = |rows: &[(String, String)]| -> String {
            rows.iter()
                .find(|(k, _)| k == "alignment_orientation")
                .map(|(_, v)| v.clone())
                .expect("alignment_orientation rule row")
        };
        assert_eq!(value(&off), "both (+/-)");
        assert!(value(&on).starts_with("forward-only (+)"));
        assert_eq!(
            value(&off), value(&dna),
            "the flag must be inert on a reference-oriented substrate — same certificate as off"
        );
    }

    #[test]
    fn joint_certificate_reports_corroborated_edges_without_changing_membership() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            return;
        }
        let mut rng: u64 = 0xd1a0_2026;
        let genomic: Vec<u8> = (0..5000)
            .map(|_| {
                rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[((rng >> 33) % 4) as usize]
            })
            .collect();
        let exon_sum = genomic[500..3500].to_vec();
        let fam = vec![
            DenovoTranscript {
                tid: "rna_a".into(),
                chrom: "a".into(),
                start: 0,
                end: 5000,
                n_reads: 10,
                strand: '+',
                introns: vec![(1000, 3000)],
                seq: exon_sum.clone(),
                ..Default::default()
            },
            DenovoTranscript {
                tid: "rna_b".into(),
                chrom: "b".into(),
                start: 0,
                end: 5000,
                n_reads: 9,
                strand: '+',
                introns: vec![(1000, 3000)],
                seq: exon_sum,
                ..Default::default()
            },
        ];
        let tmp = tempfile::tempdir().expect("tempdir");
        let fasta = tmp.path().join("genome.fa");
        let seq = String::from_utf8(genomic).expect("DNA alphabet");
        std::fs::write(&fasta, format!(">a\n{seq}\n>b\n{seq}\n")).expect("write fasta");
        let out = tmp.path().join("joint");
        let params = RefineParams {
            threads: 1,
            require_forward_alignment: true,
            ..Default::default()
        };
        write_joint_rna_dna_certificate(
            out.to_str().unwrap(),
            &[fam],
            fasta.to_str().unwrap(),
            &params,
        )
        .expect("joint certificate");

        let edges = std::fs::read_to_string(out.with_extension("joint_edges.tsv")).expect("edges");
        assert!(edges.lines().skip(1).any(|r| r.contains("\ttrue\ttrue\tRNA_DNA\ttrue")), "{edges}");
        let families = std::fs::read_to_string(out.with_extension("joint_families.tsv")).expect("families");
        assert!(families.lines().skip(1).any(|r| r.ends_with("\tBOTH_CONNECTED\tCONCORDANT")), "{families}");
        let rules = std::fs::read_to_string(out.with_extension("joint_rule.tsv")).expect("rules");
        assert!(rules.contains("membership_effect\tnone (report/certificate only)"));
    }

    // ---- X.4: ONE TIER, and an instrument that can prove it ------------------------------------
    //
    // The 2026-08-07 "sensitive-only default" was written inline inside `homology_edges_all_reps_pooled`.
    // `refine_families_exon_sum` — the edge step `copy_assign` runs — called `primary_seed_args()`
    // unconditionally at two sites, so O2 kept running `-x asm20` at 0.80 (10 of 13 panel calls) while
    // O1 ran `-k 11 -w 5` at 0.60 (25 of 25). These three tests are what makes that non-reproducible.

    /// Serializes the two tests below that touch `RUSTLE_ER_SENSITIVE_ONLY`, which is process-global.
    static X4_ENV_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

    #[test]
    fn er_primary_tier_is_the_sensitive_seed_by_default() {
        let _g = X4_ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        let saved = std::env::var("RUSTLE_ER_SENSITIVE_ONLY").ok();
        std::env::remove_var("RUSTLE_ER_SENSITIVE_ONLY");
        let p = RefineParams::default();
        let (seed, floor, mask) = er_primary_tier(&p);
        assert_eq!(seed, ER_SENSITIVE_SEED.iter().map(|s| s.to_string()).collect::<Vec<_>>());
        assert_eq!(floor, p.sensitive_identity, "the floor must travel WITH the seed, not stay at min_identity");
        assert_eq!(mask, TIER_SENSITIVE);
        // Opting out restores the legacy primary EXACTLY, floor included.
        std::env::set_var("RUSTLE_ER_SENSITIVE_ONLY", "0");
        let (seed, floor, mask) = er_primary_tier(&p);
        assert_eq!(seed, primary_seed_args());
        assert_eq!(floor, p.min_identity);
        assert_eq!(mask, TIER_ASM20);
        match saved {
            Some(v) => std::env::set_var("RUSTLE_ER_SENSITIVE_ONLY", v),
            None => std::env::remove_var("RUSTLE_ER_SENSITIVE_ONLY"),
        }
        // `--no-sensitive` also falls back to asm20 without touching the env var.
        let off = RefineParams { nucleotide_sensitive: false, ..RefineParams::default() };
        assert_eq!(er_primary_tier(&off).2, TIER_ASM20);
    }

    /// THE PARITY CERTIFICATE, as an assertion — **on the axes the two sites are supposed to share**.
    /// Given the SAME site descriptor the shared `er_rule_rows` must emit byte-identical rows; that is
    /// what makes `diff <prefix>.rule.tsv <prefix>.refine.rule.tsv` a real answer instead of a claim.
    ///
    /// ⚠ O-3 — this test used to be read as *"O1 and O2 run the same rule"*, which it never checked: it
    /// hands BOTH sites a hand-written descriptor. The sites' ACTUAL descriptors differ on the additive
    /// genomic leg, and that is asserted separately in
    /// `the_o1_and_o2_certificates_differ_only_on_the_additive_genomic_tier`.
    #[test]
    fn er_rule_rows_are_identical_at_the_o1_and_o2_call_sites() {
        let _g = X4_ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        let saved = std::env::var("RUSTLE_ER_SENSITIVE_ONLY").ok();
        std::env::remove_var("RUSTLE_ER_SENSITIVE_ONLY");
        let p = RefineParams::default();
        let site = |cores: bool| ErRuleSite {
            substrate_genomic: false,
            core_lens_supplied: cores,
            genomic_tier: GenomicTier::NotPresent,
        };
        // O1: `homology_edges_all_reps_pooled` — exon-sum substrate, core lengths supplied.
        let o1 = er_rule_rows(&p, &site(true));
        // O2: `refine_families_exon_sum` — exon-sum substrate, no core lengths (it passes `cores: None`).
        let o2 = er_rule_rows(&p, &site(false));
        match saved {
            Some(v) => std::env::set_var("RUSTLE_ER_SENSITIVE_ONLY", v),
            None => std::env::remove_var("RUSTLE_ER_SENSITIVE_ONLY"),
        }
        assert_eq!(o1, o2, "O1 and O2 must emit the same E_r rule at shipped defaults");
        let m: std::collections::BTreeMap<&str, &str> =
            o1.iter().map(|(k, v)| (k.as_str(), v.as_str())).collect();
        assert_eq!(m["primary_tier_seed"], "-k 11 -w 5");
        assert_eq!(m["primary_tier_identity"], "0.600000");
        assert_eq!(m["sensitive_only"], "true");
        // No count, length or path may enter the rule file: a data row would turn the parity diff into
        // a data diff and the certificate would stop meaning anything.
        for k in m.keys() {
            assert!(
                !k.starts_with("n_") && !k.contains("median") && !k.contains("paf"),
                "`{k}` is data-dependent and must live in params.tsv, not rule.tsv"
            );
        }
    }

    /// ⭐ O-3 — **THE ASYMMETRY THE CERTIFICATE USED TO HIDE, AS AN ASSERTION.** At shipped defaults with a
    /// genome reachable, O2's refine unions an additive GENOMIC-SPAN leg into `E_r` and O1's catalog does
    /// not (it SWAPS its substrate instead, `homology_genomic_span` default OFF). So the two rule files
    /// must differ, and must differ on exactly ONE key. Before this, both printed `substrate = exon-sum`
    /// and the diff was empty — a certificate that answered "same rule?" with "yes" about a rule that
    /// differed.
    ///
    /// ⚠ This test is the reason the X.2 claim *"`diff O1.rule.tsv O2.refine.rule.tsv` is EMPTY"* must be
    /// requoted: it was empty because the file could not see the substrate, not because the rules agreed.
    #[test]
    fn the_o1_and_o2_certificates_differ_only_on_the_additive_genomic_tier() {
        let _g = X4_ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        let saved = std::env::var("RUSTLE_ER_SENSITIVE_ONLY").ok();
        std::env::remove_var("RUSTLE_ER_SENSITIVE_ONLY");
        let p = RefineParams::default();
        // The descriptors the two sites ACTUALLY build (see the two `ErRuleSite { … }` literals).
        let o1 = er_rule_rows(
            &p,
            &ErRuleSite {
                substrate_genomic: false,
                core_lens_supplied: true,
                genomic_tier: GenomicTier::NotPresent,
            },
        );
        let o2 = er_rule_rows(
            &p,
            &ErRuleSite {
                substrate_genomic: false,
                core_lens_supplied: false,
                genomic_tier: additive_genomic_tier(&p, true),
            },
        );
        match saved {
            Some(v) => std::env::set_var("RUSTLE_ER_SENSITIVE_ONLY", v),
            None => std::env::remove_var("RUSTLE_ER_SENSITIVE_ONLY"),
        }
        let diff: Vec<&str> = o1
            .iter()
            .zip(o2.iter())
            .filter(|(a, b)| a != b)
            .map(|(a, _)| a.0.as_str())
            .collect();
        assert_eq!(
            diff,
            vec!["additive_genomic_tier"],
            "the O1/O2 rule diff must show the additive genomic leg and NOTHING else; got {diff:?}"
        );
        let get = |rows: &[(String, String)], k: &str| {
            rows.iter().find(|(a, _)| a == k).map(|(_, v)| v.clone()).unwrap_or_default()
        };
        assert!(get(&o2, "additive_genomic_tier").starts_with("armed"));
        assert!(get(&o1, "additive_genomic_tier").starts_with("absent"));
        // Both still agree on the CORE substrate — the point is that the core is not the whole rule.
        assert_eq!(get(&o1, "core_substrate"), SUBSTRATE_EXON_SUM);
        assert_eq!(get(&o2, "core_substrate"), SUBSTRATE_EXON_SUM);
    }

    /// ⭐ O-4 — **THE DIVERGENCE MUST NOT BE SILENT, AND ITS EVIDENCE MUST BE IN THE TREE.**
    /// O-3 made the two sites' one differing line VISIBLE; a bare difference still reads as drift. This
    /// asserts the certificate carries, identically at both sites, a row saying the divergence was
    /// measured and kept, that the row NAMES the key that differs, and — the part that actually failed
    /// for a month — that the evidence it cites EXISTS. `bench/FALSE_NEGATIVES.md` was deleted in
    /// `9b0814f` while `family_refine` went on citing it, leaving a default-ON edge rule with no
    /// justification in the tree. A doc reference that no test reads is a doc reference that rots.
    #[test]
    fn the_site_divergence_is_declared_identically_at_both_sites_and_its_evidence_exists() {
        let _g = X4_ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        let p = RefineParams::default();
        let get = |rows: &[(String, String)], k: &str| {
            rows.iter().find(|(a, _)| a == k).map(|(_, v)| v.clone())
        };
        let o1 = er_rule_rows(
            &p,
            &ErRuleSite {
                substrate_genomic: false,
                core_lens_supplied: true,
                genomic_tier: GenomicTier::NotPresent,
            },
        );
        let o2 = er_rule_rows(
            &p,
            &ErRuleSite {
                substrate_genomic: false,
                core_lens_supplied: false,
                genomic_tier: additive_genomic_tier(&p, true),
            },
        );
        let a = get(&o1, "site_divergence_policy").expect("O1 certificate must declare the policy");
        let b = get(&o2, "site_divergence_policy").expect("O2 certificate must declare the policy");
        // CONSTANT across sites: the policy explains the diff, it must never widen it.
        assert_eq!(a, b, "the policy row is a constant; if it differs it becomes a second drift axis");
        // It must name the key that actually differs, or it explains the wrong line.
        assert!(
            a.contains("additive_genomic_tier"),
            "the policy must name the diverging key; got `{a}`"
        );
        // ⚠ THE GUARD THAT WOULD HAVE CAUGHT THE ORIGINAL DEFECT: the cited evidence must be in the tree.
        let cited = "bench/FALSE_NEGATIVES.md";
        assert!(a.contains(cited), "the policy must cite its evidence; got `{a}`");
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(cited);
        assert!(
            path.is_file(),
            "`{cited}` is cited by the E_r certificate AND by `family_refine`'s genomic-span tier, so it \
             must exist. It was deleted once (9b0814f) and the default-ON tier cited a missing file for a \
             month. Restore it rather than dropping the citation."
        );
        let body = std::fs::read_to_string(&path).expect("cited evidence must be readable");
        assert!(
            body.contains("genomic-span"),
            "the cited file must actually be about the genomic-span tier it justifies"
        );
    }

    /// ⭐ O-3 — the arming decision is a total function of `(include_introns, genome available)`, and every
    /// off-state says WHY. An unexplained `off` cannot be told apart from a missing genome, which is the
    /// difference between a rule and an accident of deployment.
    #[test]
    fn additive_genomic_tier_states_are_exhaustive_and_self_explaining() {
        let mut p = RefineParams::default();
        assert!(!p.include_introns, "shipped default: the core tier is the exon-sum");
        assert_eq!(additive_genomic_tier(&p, true), GenomicTier::Armed);
        assert_eq!(additive_genomic_tier(&p, false), GenomicTier::NoGenome);
        p.include_introns = true;
        // include_introns runs the CORE on genomic, so there is no second substrate to union in — the
        // gate's `if !params.include_introns` and this arm must never disagree.
        assert_eq!(additive_genomic_tier(&p, true), GenomicTier::CoreIsGenomic);
        assert_eq!(additive_genomic_tier(&p, false), GenomicTier::CoreIsGenomic);
        assert!(GenomicTier::Armed.armed());
        for t in [GenomicTier::CoreIsGenomic, GenomicTier::NoGenome, GenomicTier::NotPresent] {
            assert!(!t.armed(), "{t:?} must not arm the leg");
            let l = t.label();
            assert!(
                l.starts_with("off (") || l.starts_with("absent ("),
                "every non-armed state must state its reason; got `{l}`"
            );
        }
        // One spelling of each substrate, shared with the dump tags and the `[refine]` log line.
        assert_eq!(substrate_name(false), SUBSTRATE_EXON_SUM);
        assert_eq!(substrate_name(true), SUBSTRATE_GENOMIC_SPAN);
    }

    /// ⭐ O-3 SOURCE GUARD — **ONE SOURCE OF TRUTH, the `ER_TIER_FLAGS` precedent.** The certificate must
    /// be built from the SAME value the gate branched on. A second `additive_genomic_tier(...)` call in
    /// the dump block, or a hand-written `genomic_tier:` literal there, would re-create precisely the
    /// defect this fixes: a file that is right about the code and wrong about the run.
    #[test]
    fn refine_certificate_derives_the_genomic_tier_from_the_gate_not_a_second_derivation() {
        let src = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/rustle/vg_family/denovo_pipeline.rs"
        ))
        .expect("read denovo_pipeline.rs");
        let start = src.find("pub fn refine_families_exon_sum(").expect("refine_families_exon_sum not found");
        let end = src[start..].find("\n/// The sequence used to compare a copy").expect("end of refine not found");
        let body = &src[start..start + end];
        let code: String =
            body.lines().map(|l| l.split("//").next().unwrap_or("")).collect::<Vec<_>>().join("\n");
        assert_eq!(
            code.matches("additive_genomic_tier(").count(),
            1,
            "the additive genomic leg must be decided ONCE in refine; a second derivation is how a \
             certificate comes to describe a run that did not happen"
        );
        assert!(
            code.contains("genomic_tier.armed()"),
            "the gate must branch on the shared `GenomicTier`, not re-test `include_introns` / `genome`"
        );
        assert!(
            code.contains("genomic_tier,"),
            "the emitted ErRuleSite must carry the gate's own `genomic_tier` value"
        );
        // The old inline condition must not come back alongside it.
        assert!(
            !code.contains("if !params.include_introns {\n                if let Some(g) = genome"),
            "the nested include_introns/genome gate is what made the leg unobservable"
        );
    }

    /// Source-level guard: refine must not reach for `primary_seed_args()` again. That call is exactly
    /// how the tier drifted — it bypasses `RUSTLE_ER_SENSITIVE_ONLY` and pairs the seed with
    /// `min_identity` instead of `sensitive_identity`.
    #[test]
    fn refine_resolves_its_tier_only_through_er_primary_tier() {
        let src = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/rustle/vg_family/denovo_pipeline.rs"
        ))
        .expect("read denovo_pipeline.rs");
        let start = src.find("pub fn refine_families_exon_sum(").expect("refine_families_exon_sum not found");
        let end = src[start..].find("\n/// The sequence used to compare a copy").expect("end of refine not found");
        let body = &src[start..start + end];
        // CODE only — the function's own commentary NAMES the retired call in order to explain why it is
        // retired, and a scan that cannot tell those apart fails on its own documentation.
        let code: String = body
            .lines()
            .map(|l| l.split("//").next().unwrap_or(""))
            .collect::<Vec<_>>()
            .join("\n");
        assert!(
            !code.contains("primary_seed_args()"),
            "refine must resolve its seed via `er_primary_tier`; a direct `primary_seed_args()` call \
             ignores RUSTLE_ER_SENSITIVE_ONLY and silently re-splits O1 and O2 into two tiers"
        );
        assert!(code.contains("er_primary_tier(params)"), "refine must call `er_primary_tier`");
        assert!(
            code.contains("er_rule_rows(params, &site)"),
            "refine must emit the shared rule rows, or the O1/O2 params diff is unavailable inside copy_assign"
        );
    }

    /// ⭐ O-3 — **THE CERTIFICATE MUST NAME THE SUBSTRATE THAT ACTUALLY RAN.** `family_refine` runs its
    /// core tier on the exon-sum and then, gated only on `!edges_connect_all`, UNIONS IN an additive
    /// GENOMIC-SPAN tier (default on whenever a genome is reachable). Before this test the certificate
    /// wired `substrate` to `params.include_introns` alone, so a run whose edge set is `E_x ∪ E_g` wrote
    /// `substrate = exon-sum` and carried **no row at all** for the additive tier — and "diff two
    /// `params.tsv` files", the project's official way to settle *do O1 and O2 use the same rule?*, was
    /// blind on exactly that axis.
    ///
    /// This is an END-TO-END assertion on the FILES, not on `er_rule_rows` in isolation: the failure it
    /// guards against is a certificate that disagrees with the run, which only the emitted bytes can show.
    /// The fixture is the one from `genomic_span_edges_link_fragments_that_exon_sum_cannot` — two ~99%
    /// identical loci whose "assembled" exon-sums are DISJOINT, so the core tier cannot link them and the
    /// family exists ONLY because the genomic leg fired.
    #[test]
    fn refine_certificate_reports_the_additive_genomic_tier_that_actually_ran() {
        use crate::vg_family::family_detect::DenovoTranscript;
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            eprintln!("minimap2 absent; skip");
            return;
        }
        let fa = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/from_genome/subset.fa");
        if std::fs::metadata(fa).is_err() {
            eprintln!("fixture absent; skip");
            return;
        }
        // The dump prefix is process-global env, and so is the tier lock every other X.4 test takes.
        let _g = X4_ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        let dir = std::env::temp_dir().join(format!("rustle_o3_cert_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("mkdir dump dir");
        let prefix = dir.join("cert");
        let mk = |tid: &str, chrom: &str, end: u64, seq: Vec<u8>| DenovoTranscript {
            tid: tid.into(),
            chrom: chrom.into(),
            start: 0,
            end,
            n_reads: 5,
            strand: '+',
            introns: vec![],
            seq,
            distinguishing_uniq: 0,
            core_bp: 0,
            stub: false,
            tes: None,
        };
        // 431 bp is a FINGERPRINT, not a magic number: `RUSTLE_ER_EDGE_DUMP` is process-global, so any
        // other test calling refine while this one holds it drops its own dump in the same directory
        // under `.call<N>`. `substrate_median_len_bp` = 431 with 2 reps and 1 family in identifies THIS
        // call without asking the file the question the test is here to ask. (Selecting on the answer —
        // "the params file that says the tier fired" — would be the classic denominator-conditioned-on-
        // the-prediction defect; it would pass on an empty dump.)
        let fam = vec![
            mk("a", "NCF1", 15440, b"A".repeat(431)),
            mk("b", "NCF1B", 15319, b"C".repeat(431)),
        ];
        let mut p = homology_refine_params(Some(0.90), 2);
        p.intron_fasta = Some(fa.to_string());
        assert!(!p.include_introns, "the core substrate must be the exon-sum for this test to mean anything");
        std::env::set_var("RUSTLE_ER_EDGE_DUMP", prefix.to_str().unwrap());
        let out = refine_families_exon_sum(vec![fam], &p, None, 1);
        std::env::remove_var("RUSTLE_ER_EDGE_DUMP");
        let out = out.expect("refine");
        // The family exists ONLY via the additive genomic leg: the exon-sums are disjoint 400-mers.
        assert_eq!(out.len(), 1, "the genomic-span tier must union the two disjoint exon-sums into one family");

        let kv = |path: &std::path::Path| -> std::collections::BTreeMap<String, String> {
            std::fs::read_to_string(path)
                .unwrap_or_default()
                .lines()
                .skip(1)
                .filter_map(|l| l.split_once('\t'))
                .map(|(k, v)| (k.to_string(), v.to_string()))
                .collect()
        };
        // Pick THIS call's dump by fingerprint rather than by filename: `REFINE_DUMP_SEQ` renames the
        // second and later refine calls in the process to `.call<N>`, and OTHER tests' refine calls land
        // in this directory too while the env var is set (that is not hypothetical — it is what the
        // first full-suite run of this test caught).
        let mut mine: Vec<(std::path::PathBuf, std::collections::BTreeMap<String, String>)> = Vec::new();
        let mut seen: Vec<String> = Vec::new();
        for e in std::fs::read_dir(&dir).expect("read dump dir") {
            let path = e.expect("dirent").path();
            if !path.to_string_lossy().ends_with(".params.tsv") {
                continue;
            }
            let m = kv(&path);
            seen.push(format!(
                "{} [n_reps={:?} n_families_in={:?} substrate_median_len_bp={:?}]",
                path.display(),
                m.get("n_reps"),
                m.get("n_families_in"),
                m.get("substrate_median_len_bp")
            ));
            if m.get("n_reps").map(String::as_str) == Some("2")
                && m.get("n_families_in").map(String::as_str) == Some("1")
                && m.get("substrate_median_len_bp").map(String::as_str) == Some("431")
            {
                mine.push((path, m));
            }
        }
        assert_eq!(
            mine.len(),
            1,
            "the 431 bp / 2-rep / 1-family fingerprint must identify exactly one dump; saw {seen:#?}"
        );
        let (ppath, par) = mine.pop().expect("checked len == 1");
        let rpath = std::path::PathBuf::from(ppath.to_string_lossy().replace(".params.tsv", ".rule.tsv"));
        let rule = kv(&rpath);
        let shown = std::fs::read_to_string(&rpath).unwrap_or_default();

        // (1) THE RULE FILE — the parity object — must say the additive tier was armed.
        let tier = rule.get("additive_genomic_tier").map(String::as_str).unwrap_or("<MISSING>");
        assert!(
            tier.starts_with("armed"),
            "rule.tsv must record the additive GENOMIC-SPAN tier as armed on a run whose edge set is \
             E_x u E_g; got `additive_genomic_tier = {tier}`.\n--- rule.tsv as emitted ---\n{shown}"
        );
        // (2) ...and must not let `core_substrate` be read as "this is all that ran".
        assert_eq!(
            rule.get("core_substrate").map(String::as_str),
            Some("exon-sum"),
            "the core substrate row must be named `core_substrate`, not `substrate`: a run that unions a \
             second substrate in has no single `substrate`.\n--- rule.tsv as emitted ---\n{shown}"
        );
        // (3) THE PARAMS FILE — the data side — must say it FIRED and how much it contributed.
        assert_eq!(
            par.get("n_families_genomic_tier_ran").map(String::as_str),
            Some("1"),
            "params.tsv must record that the genomic leg RAN on this family; got {:?}",
            par.get("n_families_genomic_tier_ran")
        );
        let added: usize =
            par.get("n_edges_genomic_tier_added").and_then(|v| v.parse().ok()).unwrap_or(0);
        assert!(
            added >= 1,
            "params.tsv must record the edges the genomic leg CONTRIBUTED (the family only exists \
             because of them); got n_edges_genomic_tier_added = {:?}",
            par.get("n_edges_genomic_tier_added")
        );
        // `RUSTLE_O3_KEEP_DUMP=1` leaves the emitted certificate on disk: the before/after bytes of this
        // dump ARE the evidence that the fix landed, and a passing test that deletes them is unquotable.
        if std::env::var("RUSTLE_O3_KEEP_DUMP").map(|v| v.is_empty()).unwrap_or(true) {
            let _ = std::fs::remove_dir_all(&dir);
        } else {
            eprintln!("[o3] certificate kept at {}", dir.display());
        }
    }

    #[test]
    fn primary_seed_args_defaults_to_asm20_and_is_overridable() {
        // Default must stay asm20 so no existing catalog moves.
        std::env::remove_var("RUSTLE_MM2_SEED");
        assert_eq!(primary_seed_args(), vec!["-x".to_string(), "asm20".to_string()]);
        // An override is split on whitespace into separate argv entries -- passing "-k 9 -w 3" as ONE
        // argument would make minimap2 reject it.
        std::env::set_var("RUSTLE_MM2_SEED", "-k 9 -w 3");
        assert_eq!(primary_seed_args(), vec!["-k", "9", "-w", "3"]);
        // Blank / whitespace-only falls back rather than passing an empty argv.
        std::env::set_var("RUSTLE_MM2_SEED", "   ");
        assert_eq!(primary_seed_args(), vec!["-x".to_string(), "asm20".to_string()]);
        std::env::remove_var("RUSTLE_MM2_SEED");
    }

    #[test]
    fn tes_extension_appends_only_outward_and_never_trims() {
        use crate::vg_family::family_detect::DenovoTranscript;
        // Exon-sum "ACGT" for a copy at chr1:100-104; a sharp TES at 110 should append the 6 genomic bases
        // 104..110 in transcription orientation. With the env knob OFF nothing may change -- that is what
        // keeps every existing catalog byte-identical.
        let base = DenovoTranscript {
            chrom: "c1".into(), start: 100, end: 104, strand: '+',
            seq: b"ACGT".to_vec(), ..Default::default()
        };
        // A real genome so the extension can actually be fetched. An earlier version of this test passed
        // `None` for the genome, which made every case a trivial no-op and hid that the production call site
        // ALSO passed None -- the extension silently never fired on a full Soto run.
        let dir = std::env::temp_dir().join(format!("tes_ext_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("g.fa");
        // c1: 100 A's, then 10 C's at 100..110, then T's
        let seq: String = "A".repeat(100) + &"C".repeat(10) + &"T".repeat(50);
        std::fs::write(&fa, format!(">c1\n{seq}\n")).unwrap();
        let contigs: std::collections::HashSet<String> = ["c1".to_string()].into_iter().collect();
        let g = GenomeIndex::from_fasta_contigs(fa.to_str().unwrap(), &contigs).unwrap();

        // '+' strand, exon-sum ends at 104, sharp TES at 110 -> append genome[104..110] = "CCCCCC".
        let fwd = DenovoTranscript { tes: Some(110), ..base.clone() };
        assert_eq!(extend_exon_sum_to_tes(&fwd, &g), Some(b"ACGTCCCCCC".to_vec()));

        // An INWARD tes must never trim: extension may only ADD sequence.
        assert_eq!(extend_exon_sum_to_tes(&DenovoTranscript { tes: Some(102), ..base.clone() }, &g), None);
        // tes == end is a no-op, not a zero-length fetch.
        assert_eq!(extend_exon_sum_to_tes(&DenovoTranscript { tes: Some(104), ..base.clone() }, &g), None);
        // no tes at all -> nothing to do
        assert_eq!(extend_exon_sum_to_tes(&base, &g), None);

        // '-' strand extends from `start` DOWNWARD, reverse-complemented. start=100, tes=90 -> genome[90..100]
        // is "AAAAAAAAAA", revcomp "TTTTTTTTTT", appended in transcription orientation.
        let rev = DenovoTranscript { strand: '-', tes: Some(90), ..base.clone() };
        assert_eq!(extend_exon_sum_to_tes(&rev, &g), Some(b"ACGTTTTTTTTTTT".to_vec()));
        // a '-' tes ABOVE start is inward -> no-op
        assert_eq!(extend_exon_sum_to_tes(&DenovoTranscript { strand: '-', tes: Some(200), ..base.clone() }, &g), None);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn denovo_config_tied_seed_defaults_off() {
        assert!(!DenovoConfig::default().tied_seed);
    }

    #[test]
    fn mischain_salvage_defaults_off() {
        assert!(!DenovoConfig::default().mischain_salvage);
    }

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

    /// `copy_introns`/`copy_map_identity` must exist and stay length-aligned with `copy_tids` at every
    /// construction site (`FamilyAssignment::empty()` starts both empty; a built family -- see
    /// `admission_widens_abundance_not_left_at_zero` -- must keep all three the same length).
    #[test]
    fn family_assignment_has_aligned_intron_and_identity_fields() {
        let fa = FamilyAssignment::empty();
        assert_eq!(fa.copy_introns.len(), 0);
        assert_eq!(fa.copy_map_identity.len(), 0);
        assert_eq!(fa.copy_introns.len(), fa.copy_tids.len());
        assert_eq!(fa.copy_map_identity.len(), fa.copy_tids.len());
    }

    #[test]
    /// Anti-FRAGMENTATION invariant: two copies of one family whose ASSEMBLED transcripts cover DISJOINT
    /// exon subsets share almost no exon-sum sequence, so the default (exon-sum) edge does NOT link them.
    /// Computing the edge on the GENOMIC SPAN instead recovers the link. Uses the two-copy fixture whose
    /// contigs hold the same 600 bp motif twice.
    #[test]
    fn genomic_span_edges_link_fragments_that_exon_sum_cannot() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
        let fa = "tests/fixtures/from_genome/subset.fa";
        if std::fs::metadata(fa).is_err() { eprintln!("fixture absent; skip"); return; }
        // Two copies of the SAME family (NCF1 / NCF1B contigs are ~99% identical), but their "assembled"
        // exon-sums are DISJOINT slices: copy A carries the head, copy B the tail. They cannot align.
        let head: Vec<u8> = b"A".repeat(400);            // stand-in fragment A (no similarity to B)
        let tail: Vec<u8> = b"C".repeat(400);            // stand-in fragment B
        let mk = |tid: &str, chrom: &str, start: u64, end: u64, seq: Vec<u8>| DenovoTranscript {
            tid: tid.into(), chrom: chrom.into(), start, end,
            n_reads: 5, strand: '+', introns: vec![], seq, distinguishing_uniq: 0, core_bp: 0, stub: false, tes: None,
        };
        let reps = vec![
            mk("a", "NCF1", 0, 15440, head),
            mk("b", "NCF1B", 0, 15319, tail),
        ];
        let mut p = homology_refine_params(Some(0.90), 2);
        // exon-sum substrate: the disjoint fragments cannot align -> no edge
        let e_exon = homology_edges_all_reps(&reps, &p).unwrap();
        assert!(e_exon.is_empty(), "disjoint exon-sums must NOT link, got {e_exon:?}");
        // genomic-span substrate: the underlying loci are ~99% identical -> edge forms
        p.homology_genomic_span = true;
        p.intron_fasta = Some(fa.to_string());
        let e_span = homology_edges_all_reps(&reps, &p).unwrap();
        assert!(!e_span.is_empty(), "genomic spans of the two paralogous loci must link");
    }

    /// O1 INDEPENDENCE GUARD. The `--homology-primary` claim we make to the advisor is that family
    /// MEMBERSHIP there is decided by sequence homology alone: E_r proposes the family, and the read-conflict
    /// graph E_c plays no part. That claim is false of `--cross-chrom`, where the conflict graph proposes
    /// candidate families and E_r can only split them -- so the distinction is load-bearing and must not
    /// erode by someone adding one convenient `conflict_edges` call.
    ///
    /// This is a SOURCE-level guard rather than a behavioural one because "never calls X" is not observable
    /// from outputs: a conflict-derived merge and a homology-derived merge produce the same struct. It reads
    /// the body of `detect_homology_catalog_genome_wide` and fails if any conflict-graph constructor appears.
    ///
    /// Read evidence that IS allowed on this path, and why it is not E_c:
    ///   - locus discovery (primary alignments -> skeletons -> reps): decides WHERE loci are, not who is
    ///     related to whom.
    ///   - `distinguishing_uniq` (per-locus count of MAPQ>0 placements): a per-locus scalar consumed only by
    ///     `distinct_locus_reps`, which collapses OVERLAPPING same-strand copies inside one block. It can
    ///     never add a member and never links two blocks, so it cannot create a family.
    /// Neither is a pairwise ambiguity edge, which is what E_c is.
    /// The coverage split must SPLIT a block whose two halves share only a fragment, and must LEAVE ALONE a
    /// block whose members align full-length. Both halves matter: a split rule that fires on everything
    /// would destroy recall, and one that never fires is inert.
    #[test]
    fn coverage_split_separates_fragment_sharers_and_keeps_full_length_pairs() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            return;
        }
        let mk = |tid: &str, seq: &[u8]| DenovoTranscript {
            tid: tid.into(), chrom: "chr1".into(), start: 0, end: seq.len() as u64,
            n_reads: 5, strand: '+', introns: vec![], seq: seq.to_vec(),
            distinguishing_uniq: 0, core_bp: 0, stub: false, tes: None,
        };
        let core = rand_seq(600, 0xA11CE);
        let tail_a = rand_seq(1400, 0xBEEF);
        let tail_b = rand_seq(1400, 0xF00D);
        // a1/a2 are near-identical full length; b1 shares ONLY the 600 bp core with them
        let a1: Vec<u8> = core.iter().chain(tail_a.iter()).copied().collect();
        let mut a2 = a1.clone();
        a2[50] = if a2[50] == b'A' { b'C' } else { b'A' };
        let b1: Vec<u8> = core.iter().chain(tail_b.iter()).copied().collect();
        let reps = vec![mk("a1", &a1), mk("a2", &a2), mk("b1", &b1)];
        let params = homology_refine_params(Some(0.80), 2);

        let edges = coverage_edges_all_reps(&reps, 0.90, &params).unwrap();
        let groups = coverage_split_block(&[0, 1, 2], &edges);
        let of = |i: usize| groups.iter().position(|g| g.contains(&i)).unwrap();
        assert_eq!(of(0), of(1), "full-length near-identical copies must stay together");
        assert_ne!(of(0), of(2), "a copy sharing only a fragment must split off at coverage 0.90");

        // floor 0 disables the split entirely (the opt-out contract)
        let none = coverage_edges_all_reps(&reps, 0.0, &params).unwrap();
        assert!(none.is_empty(), "min_cov 0 must yield no split edges");
        assert_eq!(coverage_split_block(&[0, 1, 2], &none).len(), 1, "no edges must leave the block whole");
        // a singleton block is returned unchanged
        assert_eq!(coverage_split_block(&[0], &edges), vec![vec![0]]);
    }

    /// A block the aligner cannot speak to must survive intact. Short exon-sums (the 108-138 bp fixture)
    /// produce NO asm20 alignment at all; treating that silence as "unrelated" shattered a legitimate
    /// 2-copy family into singletons and made `homology_catalog_groups_fixture_family` return `[]`.
    #[test]
    fn coverage_split_leaves_a_block_intact_when_the_aligner_finds_nothing() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            return;
        }
        let mk = |tid: &str, seq: &[u8]| DenovoTranscript {
            tid: tid.into(), chrom: "c1".into(), start: 0, end: seq.len() as u64,
            n_reads: 9, strand: '+', introns: vec![], seq: seq.to_vec(),
            distinguishing_uniq: 0, core_bp: 0, stub: false, tes: None,
        };
        // two unrelated ~110 bp sequences: far too short for a reliable alignment
        let reps = vec![mk("s1", &rand_seq(110, 0x1234)), mk("s2", &rand_seq(110, 0x9876))];
        let params = homology_refine_params(None, 2);
        let edges = coverage_edges_all_reps(&reps, 0.90, &params).unwrap();
        let groups = coverage_split_block(&[0, 1], &edges);
        assert_eq!(
            groups.len(), 1,
            "no alignment signal must mean NO split -- silence is not evidence of separateness"
        );
    }

    /// The split is DEFAULT OFF. It was briefly shipped default-on on the strength of a benchmark whose
    /// recall metric could not see families dissolved into singletons; corrected, the split costs 8.0 recall
    /// points on RNA and is dominated by raising the identity floor. This test pins the default so that
    /// cannot silently regress.
    /// Pooling isoform exons must ADD the exons the representative lacks, and must not silently change the
    /// rep-only behaviour. The case that matters is a locus whose representative is an unspliced STUB while a
    /// spliced isoform at the same locus carries the real exons -- that is 46% of the representatives
    /// covering a known family member, and it is why a stub cannot form an edge with a full transcript.
    /// Two representatives over the SAME genomic span are one locus seen twice. Their exons are the same
    /// physical DNA, so without a distinct-locus guard they align to themselves at ~100% and fabricate an
    /// edge. Genome-wide on gorilla this was 93.06% of everything the shared-exon rule reported that E_r
    /// did not (2,601/2,795), 84.01% of it on OPPOSITE strands -- a locus matching its own antisense model.
    #[test]
    fn shared_exon_does_not_link_a_locus_to_itself_across_two_representatives() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            return;
        }
        let ex = rand_seq(400, 0xC0FFEE);
        // `end - start` MUST equal `seq.len()` here: `exon_seqs_of` derives exon LENGTHS from the genomic
        // coordinates and then slices `seq` by them, so a span longer than the sequence makes every exon
        // fail the `off + l > seq.len()` bound and the locus contributes NO exons at all -- which would
        // make this test pass its guard assertion for the wrong reason.
        let mk = |tid: &str, chrom: &str, start: u64, strand: char| DenovoTranscript {
            tid: tid.into(), chrom: chrom.into(), start, end: start + ex.len() as u64, n_reads: 9, strand,
            introns: vec![], seq: ex.clone(), distinguishing_uniq: 0, core_bp: 0, stub: false, tes: None,
        };
        let params = homology_refine_params(Some(0.90), 2);

        // Same span, opposite strands -- the dominant real-data shape.
        let overlapping = vec![mk("a", "chr1", 1_000, '+'), mk("b", "chr1", 1_200, '-')];
        assert!(
            shared_exon_edges(&overlapping, 0.90, 200, &params).unwrap().is_empty(),
            "two representatives over one genomic locus must not link to each other"
        );

        // Positive control: identical sequence at DISJOINT loci is a real homology claim and must survive,
        // so the guard is not simply suppressing every edge.
        let disjoint = vec![mk("a", "chr1", 1_000, '+'), mk("b", "chr2", 9_000_000, '+')];
        assert!(
            shared_exon_edges(&disjoint, 0.90, 200, &params).unwrap().contains(&(0, 1)),
            "control: the same exon at two DISJOINT loci must still form an edge"
        );
    }

    #[test]
    fn pooled_isoform_exons_recover_what_a_stub_representative_lacks() {
        if std::process::Command::new("minimap2").arg("--version").output().is_err() {
            return;
        }
        let mk = |tid: &str, start: u64, introns: Vec<(u64, u64)>, seq: Vec<u8>| DenovoTranscript {
            tid: tid.into(), chrom: "chr1".into(), start, end: start + 1200,
            n_reads: 9, strand: '+', introns, seq, distinguishing_uniq: 0, core_bp: 0, stub: false, tes: None,
        };
        let ex_a = rand_seq(400, 0x5EED1);
        let ex_b = rand_seq(400, 0x5EED2);
        // locus 1: representative is an unspliced STUB carrying only filler; a spliced isoform at the same
        // locus carries ex_a + ex_b.
        let stub = mk("stub", 0, vec![], rand_seq(600, 0xF11E));
        let spliced: Vec<u8> = ex_a.iter().chain(ex_b.iter()).copied().collect();
        // locus 2: a full transcript sharing exon A with locus 1's discarded isoform.
        let other: Vec<u8> = ex_a.iter().chain(rand_seq(400, 0xBEEF).iter()).copied().collect();
        let reps = vec![stub.clone(), mk("other", 50_000, vec![(50_400, 50_500)], other)];

        let params = homology_refine_params(Some(0.90), 2);
        let rep_only = shared_exon_edges(&reps, 0.90, 200, &params).unwrap();
        assert!(
            rep_only.is_empty(),
            "control: with only the stub representative there is no shared exon, so no edge"
        );

        let pooled = vec![vec![ex_a.clone(), ex_b.clone()], vec![ex_a.clone()]];
        let with_pool = shared_exon_edges_pooled(&reps, &pooled, 0.90, 200, &params).unwrap();
        assert!(
            with_pool.contains(&(0, 1)),
            "pooling the discarded spliced isoform's exons must recover the edge the stub could not form"
        );
    }

    #[test]
    fn coverage_split_is_off_by_default() {
        assert!(
            coverage_split_floor() == 0.0 || std::env::var("RUSTLE_COVERAGE_SPLIT").is_ok(),
            "the coverage split must be OFF by default -- it is dominated by the identity floor on RNA"
        );
    }

    #[test]
    fn homology_catalog_never_touches_the_conflict_graph() {
        let src = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/rustle/vg_family/denovo_pipeline.rs"
        ))
        .expect("read own source");
        let start = src
            .find("pub fn detect_homology_catalog_genome_wide(")
            .expect("homology catalog entry point not found -- was it renamed?");
        // Body ends at the NEAREST following top-level item (`pub fn`/`fn` at column 0) -- the nearer of the
        // two, not whichever is checked first, or the scan would run on into unrelated functions.
        let rest = &src[start + 10..];
        let end = [rest.find("\npub fn "), rest.find("\nfn ")]
            .into_iter()
            .flatten()
            .min()
            .map(|o| start + 10 + o)
            .unwrap_or(src.len());
        // Strip line comments: a comment EXPLAINING that this path avoids the conflict graph must not trip
        // the guard, and equally a banned call must not hide behind one.
        let body: String = src[start..end]
            .lines()
            .map(|l| l.split("//").next().unwrap_or(""))
            .collect::<Vec<_>>()
            .join("\n");

        // ⚠ THIS GUARD USED TO CHECK SPELLING, NOT SEMANTICS, AND PASSED WHILE THE VIOLATION FLOWED.
        // It banned four hand-picked strings; the actual `E_c` dependency on this path uses DIFFERENT
        // names (`locus_unique_mapper_counts`, and `reads_distinguish` one level down inside
        // `distinct_locus_reps`), so none of them tripped it. Same shape as the `--cross-chrom`
        // precedent, where a flag's NAME hid what it selected. Two fixes:
        //
        //   (1) the checked set is DERIVED from the `use super::read_conflict::{...}` import rather than
        //       hand-listed, so a new `E_c` symbol cannot be added without this test seeing it; and
        //   (2) the scan is TRANSITIVE over the in-file helpers the catalog calls, because the body-only
        //       scan was structurally blind to `reads_distinguish` — it is not in the catalog body, it
        //       is in `distinct_locus_reps`, which the catalog calls.
        //
        // DISCLOSED_EC_USES is the honest ledger of where O1's node set DOES consult read evidence. It
        // is not an exemption — every entry is a place the "membership by SEQUENCE alone" claim must be
        // qualified, and the spec qualifies it. Adding to this list is a THESIS EDIT, not a code edit.
        let ec_imports: Vec<String> = {
            let u = src.find("use super::read_conflict::{").expect("read_conflict import not found");
            let close = src[u..].find("};").expect("unterminated read_conflict import") + u;
            src[u + "use super::read_conflict::{".len()..close]
                .split(',')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty())
                .collect()
        };
        // Transitive scan: the catalog body PLUS the in-file helpers it calls that can carry E_c inward.
        // ⚠ Bound each helper by BRACE MATCHING, not by "the next top-level `fn`". `distinct_locus_reps`
        // is the last column-0 `fn` in this file (everything after it is the indented `mod tests`), so a
        // next-item scan finds nothing, falls back to end-of-file, and swallows the test module — where
        // the banned names legitimately appear in this very guard's own ban list. That false-positives
        // the guard, which is worse than the leak it looks for: a test that cries wolf gets relaxed.
        let mut scanned = body.clone();
        // ⚠ BOTH names are required. `distinct_locus_reps` became a one-line wrapper when the λ
        // certificate needed the merge groups, and the E_c call (`reads_distinguish`) moved into
        // `distinct_locus_reps_grouped`. Scanning only the wrapper would make this guard pass
        // VACUOUSLY — it would find no banned symbol because it would be reading an empty body.
        // The `used || !disclosed` assertion below is what catches that, and it is why the list is
        // asserted non-trivial here rather than trusted.
        for helper in ["fn distinct_locus_reps(", "fn distinct_locus_reps_grouped("] {
            if let Some(h) = src.find(helper) {
                let open = h + src[h..].find('{').expect("helper has no body");
                let mut depth = 0usize;
                let mut hend = src.len();
                for (off, ch) in src[open..].char_indices() {
                    match ch {
                        '{' => depth += 1,
                        '}' => {
                            depth -= 1;
                            if depth == 0 {
                                hend = open + off + 1;
                                break;
                            }
                        }
                        _ => {}
                    }
                }
                scanned.push('\n');
                scanned.push_str(
                    &src[h..hend].lines().map(|l| l.split("//").next().unwrap_or(""))
                        .collect::<Vec<_>>().join("\n"),
                );
            }
        }
        const DISCLOSED_EC_USES: [&str; 2] = [
            // per-locus MAPQ>0 count, fed to `distinguishing_uniq` (catalog body)
            "locus_unique_mapper_counts",
            // O2's chi(H) predicate, inside `distinct_locus_reps`; decides MERGE/KEEP for a co-located
            // SAME-STRAND pair that the junction-asymmetry rule did not settle. Measured to decide
            // 0 of 109 such pairs on gorilla and 0 firings over 451 chr1 loci, but reachable.
            "reads_distinguish",
        ];
        for sym in &ec_imports {
            // types carry no behaviour into the node set; only callables do
            if sym.chars().next().is_some_and(|c| c.is_uppercase()) {
                continue;
            }
            let used = scanned.contains(&format!("{sym}("));
            let disclosed = DISCLOSED_EC_USES.contains(&sym.as_str());
            assert!(
                !used || disclosed,
                "O1 must decide membership from E_r alone, but the homology catalog path now reaches \
                 `{sym}` from read_conflict.rs (the E_c module) — a NEW read dependency in O1's node \
                 set. Either remove it, or add it to DISCLOSED_EC_USES *and* qualify the \
                 membership-by-sequence-alone claim in docs/seeded_family_definition.md and in \
                 bench/soto/PIPELINE_READS_TO_GRAPH.txt. Do not relax this test on its own."
            );
            assert!(
                used || !disclosed,
                "`{sym}` is listed in DISCLOSED_EC_USES but no longer appears on the homology catalog \
                 path. If the dependency is genuinely gone, delete the entry AND strengthen the claim \
                 in the spec — O1 just became sequence-alone at the node and should say so."
            );
        }
        // Accepts `homology_blocks(` or `homology_blocks_pooled(` -- both route through
        // `homology_edges_all_reps*`, i.e. E_r alone. The pooled variant only widens which EXONS feed the
        // shared-exon rule (every isoform's, not just the representative's); it introduces no read-assignment
        // and no conflict edge, so the O1-from-sequence-alone invariant is unchanged.
        // `homology_blocks_pooled_with_edges(` is accepted for the same reason as the pooled variant: it
        // returns the SAME blocks and additionally hands back the E_r edge set the blocks were cut from,
        // which the λ certificate reports on. It introduces no new edge source — λ is computed from E_r
        // and nothing else — so the O1-from-sequence-alone invariant is unchanged.
        // `..._with_edges_weighted(` is the same function returning the identity/coverage it already
        // computed instead of discarding them (ledger §5q). Same blocks, same edge SET, same source —
        // the weights are REPORTED on the copy and never re-enter the partition, so the invariant holds.
        assert!(
            body.contains("homology_blocks(")
                || body.contains("homology_blocks_pooled(")
                || body.contains("homology_blocks_pooled_with_edges(")
                || body.contains("homology_blocks_pooled_with_edges_weighted("),
            "the homology catalog must still form families via homology_blocks* (E_r)"
        );
        // The path must not reach into ConflictParams either -- that reads as an E_c dependency even where
        // only the shared read-count scalar is wanted. `locus_min_reads()` is the path-neutral accessor.
        assert!(
            !body.contains("cfg.conflict."),
            "use cfg.locus_min_reads() on the homology path, not cfg.conflict.*"
        );
        // ⚠ AND THE ACCESSOR IS AN ALIAS, WHICH IS WHY THE ASSERT ABOVE PASSES WHILE THE DEPENDENCY
        // FLOWS. `locus_min_reads()` literally returns `self.conflict.min_reads`, i.e. O2's
        // `RUSTLE_CONFLICT_MIN_READS` (default 3) under an O1-sounding name. Banning `cfg.conflict.` at
        // the call site therefore enforces SPELLING only. Pin the alias so nobody reads the assert above
        // as evidence that the homology path has its own, independent read floor — it does not.
        let accessor = src
            .find("pub fn locus_min_reads(")
            .map(|i| src[i..i + 200].to_string())
            .expect("locus_min_reads accessor not found -- was it renamed?");
        assert!(
            accessor.contains("self.conflict.min_reads"),
            "locus_min_reads() no longer aliases conflict.min_reads. If the homology path gained an \
             INDEPENDENT read floor that is an improvement to O1 ⊥ O2 and the spec should say so; \
             update this pin deliberately rather than letting the change pass unnoticed."
        );
    }

    /// The veto's arithmetic: multiplicity is DISTINCT occurrences, so overlapping hits at one locus
    /// count once and abutting/disjoint ones count separately. Pinned because §4x retains the gate
    /// disabled -- an off-by-one here would move every multiplicity if it were ever re-enabled.
    #[test]
    fn genome_multiplicity_counts_distinct_occurrences_not_records() {
        let l = |c: &str, s: u64, e: u64| (c.to_string(), s, e);
        assert_eq!(count_distinct_loci(vec![]), 0, "no hits = multiplicity 0");
        assert_eq!(count_distinct_loci(vec![l("c1", 0, 100)]), 1);
        // three records piled on ONE locus are one occurrence, not three
        assert_eq!(
            count_distinct_loci(vec![l("c1", 0, 100), l("c1", 50, 150), l("c1", 10, 20)]),
            1,
            "overlapping records at one locus must not inflate multiplicity"
        );
        // same contig, disjoint -> two; abutting (s == prev end) is disjoint by half-open convention
        assert_eq!(count_distinct_loci(vec![l("c1", 0, 100), l("c1", 200, 300)]), 2);
        assert_eq!(count_distinct_loci(vec![l("c1", 0, 100), l("c1", 100, 200)]), 2);
        // a hit on a DIFFERENT contig is always its own occurrence, even at identical coordinates
        assert_eq!(count_distinct_loci(vec![l("c1", 0, 100), l("c2", 0, 100)]), 2);
        // NPIP-shaped: 31 disjoint occurrences trip any M=20 hub test purely by copy number, which is
        // the confound that refuted the veto (§4x).
        let npip: Vec<_> = (0..31).map(|i| l("c1", i * 1000, i * 1000 + 500)).collect();
        assert_eq!(count_distinct_loci(npip), 31);
    }

    /// `RUSTLE_ER_REPEAT_GATE` semantics: absent or 0 means OFF, which is what keeps the shipped path
    /// byte-identical (verified 2026-08-26 on copies/families md5).
    #[test]
    fn repeat_gate_is_off_unless_a_positive_threshold_is_set() {
        // Parsing is pure; exercise it directly rather than mutating process env (tests share one
        // process, so an env write here would race every other test).
        let parse = |v: Option<&str>| -> Option<u32> {
            v.and_then(|v| v.parse::<u32>().ok()).filter(|&m| m > 0)
        };
        assert_eq!(parse(None), None, "unset = OFF");
        assert_eq!(parse(Some("")), None, "empty = OFF");
        assert_eq!(parse(Some("0")), None, "explicit 0 = OFF");
        assert_eq!(parse(Some("nonsense")), None, "unparseable = OFF, never a silent default");
        assert_eq!(parse(Some("20")), Some(20), "the §4x threshold");
    }

    #[test]
    fn homology_blocks_groups_identical_reps_and_isolates_unrelated() {
        // three reps: two identical sequences (should share an E_r edge -> one block) + one unrelated.
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
        let mk = |tid: &str, chrom: &str, start: u64, seq: &[u8]| DenovoTranscript {
            tid: tid.into(), chrom: chrom.into(), start, end: start + seq.len() as u64,
            n_reads: 5, strand: '+', introns: vec![], seq: seq.to_vec(), distinguishing_uniq: 0, core_bp: 0, stub: false, tes: None,
        };
        // a 400 bp "gene" duplicated at two loci, plus an unrelated 400 bp sequence.
        let a: Vec<u8> = b"ACGT".iter().cycle().take(400).copied().collect();
        let b: Vec<u8> = b"TGCA".iter().cycle().take(400).copied().collect();
        let reps = vec![
            mk("d1", "chr1", 1000, &a),
            mk("d2", "chr9", 5000, &a),   // identical to d1 -> same family
            mk("d3", "chr1", 9000, &b),   // unrelated -> its own block
        ];
        let refine = homology_refine_params(None, 2);
        let blocks = homology_blocks(&reps, &refine, 0.20).unwrap();
        // d1 and d2 land in the same block; d3 is alone.
        let block_of = |i: usize| blocks.iter().position(|bl| bl.contains(&i)).unwrap();
        assert_eq!(block_of(0), block_of(1), "identical reps must share a block");
        assert_ne!(block_of(0), block_of(2), "unrelated rep must be a separate block");
    }

    /// `--min-identity` must reach BOTH E_r floors on the `--refine` path too, not just on the
    /// homology-primary path. gw_family_catalog.rs rebuilds a RefineParams for `refine_families_exon_sum`
    /// and used to copy only min_identity + min_coverage, so `--min-identity 0.93` gave run 1 a floor of
    /// 0.93 while run 3 silently kept the 0.60 default -- and since the runs are unioned, 0.60 was the
    /// effective floor. This pins the invariant the CLI advertises.
    #[test]
    fn refine_params_carry_the_sensitive_floor_too() {
        let src = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/src/bin/gw_family_catalog.rs"))
            .expect("read gw_family_catalog.rs");
        // Anchor on the refine gate in main(). (Was `let refine = !args.no_refine;` before the D1 fix
        // made refine opt-in on the homology catalog; see `refine_enabled`.)
        let start = src
            .find("let refine = refine_enabled(")
            .expect("refine block not found in gw_family_catalog::main");
        // ⚠ THIS USED TO BE A FIXED 2,400-BYTE WINDOW, and it failed the first time an unrelated line
        // was added between the anchor and the struct — the fields had simply slid past the cutoff.
        // A magic byte count is not a statement about the code, so it is replaced by the actual
        // construct: the `RefineParams { .. }` literal in the refine branch, brace-matched. That is
        // STRICTLY TIGHTER than the window (which could also have swept in unrelated code), so this
        // is a strengthening, not a relaxation of the guard.
        let lit = start
            + src[start..]
                .find("let params = RefineParams {")
                .expect("the --refine RefineParams literal must follow the refine gate");
        let open = lit + src[lit..].find('{').expect("RefineParams literal has no body");
        let mut depth = 0usize;
        let mut close = src.len();
        for (off, ch) in src[open..].char_indices() {
            match ch {
                '{' => depth += 1,
                '}' => {
                    depth -= 1;
                    if depth == 0 {
                        close = open + off + 1;
                        break;
                    }
                }
                _ => {}
            }
        }
        let block = &src[lit..close];
        for field in ["min_identity: refine_params.min_identity",
                      "min_coverage: refine_params.min_coverage",
                      "sensitive_identity: refine_params.sensitive_identity"] {
            assert!(
                block.contains(field),
                "the --refine RefineParams must inherit `{field}` from refine_params; dropping one lets \
                 --min-identity set only some of the E_r floors, and the union takes the loosest"
            );
        }
    }

    // ---- D2: the coverage clause is scale-dependent; the substrate must be RECORDED ---------------
    //
    // `c = 0.50 of min(qlen,tlen)` is a completeness filter, not a homology criterion, and the two
    // shipped paths feed it substrates ~7.5x apart in length. The RULE is deliberately unchanged (the
    // absolute-bp alternative was swept and refuted; genomic-span-by-default moves the residual into
    // locus boundaries -- docs/seeded_family_definition.md 1a). What these defend is that the run
    // reports its own absolute demand, so the substrate stops being an UNDECLARED free parameter.

    #[test]
    fn coverage_floor_bp_demand_is_the_documented_absolute_number() {
        // c * median(|node|), rounded. 0.50 of a 1,790 bp median transcript = 895 bp -- the number the
        // 61-node panel measured for the shipped RNA substrate.
        assert_eq!(coverage_floor_bp_demand(0.50, &[1790]), Some(895));
        // lower median for even n, so the value is deterministic across runs
        assert_eq!(coverage_floor_bp_demand(0.50, &[1000, 2000]), Some(500));
        assert_eq!(coverage_floor_bp_demand(0.50, &[100, 1790, 9000]), Some(895));
        // empty node set has no scale to report -- must be NA, never 0 (0 reads as "no floor at all").
        assert_eq!(coverage_floor_bp_demand(0.50, &[]), None);
    }

    #[test]
    fn coverage_floor_bp_demand_exposes_the_substrate_scale_dependence() {
        // THE D2 DEFECT, as an assertion: one identical coverage FRACTION, two substrates, ~10x
        // different base-pair demand. Panel numbers: 8,810 bp on genomic span vs 895 bp on the rep
        // transcript, from a median spliced seqlen/span ratio of ~0.13.
        let transcript = coverage_floor_bp_demand(0.50, &[1790]).unwrap();
        let genomic_span = coverage_floor_bp_demand(0.50, &[17620]).unwrap();
        assert_eq!((transcript, genomic_span), (895, 8810));
        assert!(
            genomic_span >= 5 * transcript,
            "the same 0.50 floor must be reported as a different absolute demand per substrate; \
             collapsing them back to one number is the defect this records"
        );
    }

    /// The E_r params dump is the only machine-readable place a catalog's SCALE is recorded. If these
    /// rows go away, `<prefix>.params.tsv` documents the coverage FRACTION while silently omitting what
    /// it costs in bases, and an E_r number becomes unquotable without re-running the alignment.
    #[test]
    fn er_params_dump_records_the_substrate_scale() {
        let src = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/rustle/vg_family/denovo_pipeline.rs"
        ))
        .expect("read denovo_pipeline.rs");
        let start = src.find("fn write_er_edge_dump(").expect("write_er_edge_dump not found");
        let body = &src[start..];
        for key in ["substrate_median_len_bp", "coverage_floor_median_bp_demand", "homology_genomic_span"] {
            assert!(
                body.contains(&format!("(\"{key}\".into()")),
                "{stem}.params.tsv must record `{key}` -- the E_r edge substrate sets the ABSOLUTE \
                 demand of the 0.50 coverage floor (8,810 bp on genomic span vs 895 bp on a rep \
                 transcript), so a catalog quoted without it is unattributable",
                stem = "<prefix>"
            );
        }
    }

    /// `nodes.tsv` is the only table that covers reps which join NO family (15,905 of 17,924 in the
    /// 2026-08-20 catalog). Its ABSENCE of an exon array is what forced the 2026-08-21 SD-gap audit to
    /// run every locational query on rep genomic SPANS -- 90.83% intron by bp -- which made "the
    /// interval is inside this rep's intron" (correct rejection) indistinguishable from "real
    /// transcript sequence no node covers" (a real miss), and invalidated that angle. Dropping the
    /// column again would silently reintroduce the same ambiguity, so it is pinned at the source level.
    #[test]
    fn er_node_dump_emits_the_exon_array_not_just_the_exon_count() {
        let src = std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/src/rustle/vg_family/denovo_pipeline.rs"
        ))
        .expect("read denovo_pipeline.rs");
        let start = src.find("fn write_er_edge_dump(").expect("write_er_edge_dump not found");
        let body = &src[start..];
        for col in ["\\texons\\texon_bp", "exon_blocks_str("] {
            assert!(
                body.contains(col),
                "the node dump must emit the exon ARRAY ({col} missing): n_exon alone cannot tell an \
                 intronic interval from an uncovered one, and rep spans are 90.83% intron by bp"
            );
        }
    }

    /// BEHAVIOURAL companion to the source lint above: the lint would still pass if `exon_bp` computed
    /// garbage, and that parse-back arithmetic is the only genuinely new logic in the emitter.
    ///
    /// Uses a MINUS-strand 3-exon rep whose span (1,000 bp) is twice its exon-sum (500 bp) -- the same
    /// shape, in miniature, that makes span-based locational queries misleading at genome scale (rep
    /// spans are 90.83% intron by bp across the shipped catalog).
    #[test]
    fn er_node_dump_exon_array_agrees_with_the_spliced_length_it_reports() {
        let r = DenovoTranscript {
            tid: "t".into(),
            chrom: "c1".into(),
            start: 1000,
            end: 2000,
            n_reads: 7,
            strand: '-',
            introns: vec![(1100, 1300), (1500, 1800)],
            seq: vec![b'A'; 500],
            ..Default::default()
        };
        let dir = std::env::temp_dir().join(format!("rustle_er_dump_{}", std::process::id()));
        let _ = std::fs::create_dir_all(&dir);
        write_er_edge_dump(
            dir.join("t").to_str().expect("utf8 path"),
            &[r],
            &[vec![b'A'; 500]],
            &BTreeSet::new(),
            &BTreeMap::new(),
            &BTreeMap::new(),
            &BTreeMap::new(),
            &RefineParams::default(),
            false,
        );
        // `DUMP_SEQ` is a process-global counter, so the stem may carry a `.callN` suffix.
        let f = std::fs::read_dir(&dir)
            .expect("readdir")
            .filter_map(|e| e.ok().map(|e| e.path()))
            .find(|p| p.to_string_lossy().ends_with(".nodes.tsv"))
            .expect("the dump must write a nodes.tsv");
        let text = std::fs::read_to_string(&f).expect("read nodes.tsv");
        let hdr: Vec<&str> = text.lines().next().expect("header").split('\t').collect();
        let row: Vec<&str> = text.lines().nth(1).expect("one data row").split('\t').collect();
        let col = |n: &str| {
            row[hdr.iter().position(|h| *h == n).unwrap_or_else(|| panic!("no column {n}"))]
        };
        assert_eq!(col("exons"), "1000-1100,1300-1500,1800-2000", "genomic-ascending half-open blocks");
        assert_eq!(col("n_exon"), "3");
        // The invariant the two columns exist to expose: computed by different paths, they must agree.
        // A mismatch means the malformed-chain guard dropped a block -- visible, not silent.
        assert_eq!(col("exon_bp"), "500");
        assert_eq!(col("exon_bp"), col("exon_sum_len"), "exon_bp must equal the spliced length");
        // And the reason the array is needed at all: n_exon + span cannot recover this.
        let span = col("end").parse::<u64>().expect("end") - col("start").parse::<u64>().expect("start");
        assert_eq!(span, 1000);
        assert!(span > col("exon_bp").parse::<u64>().expect("exon_bp"), "span exceeds the exon-sum");
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The scoped-guard predicate: exempt a minus-strand record ONLY when a rep's strand was never
    /// measured AND the two spans are disjoint. Both clauses are load-bearing -- dropping the span clause
    /// admits a population 40x enriched for overlapping-span pairs (one locus counted twice).
    #[test]
    fn scoped_guard_exempts_only_unmeasured_strand_at_disjoint_loci() {
        let mk = |chrom: &str, start: u64, end: u64, spliced: bool| DenovoTranscript {
            tid: "t".into(), chrom: chrom.into(), start, end, n_reads: 5, strand: '+',
            introns: if spliced { vec![(start + 10, start + 20)] } else { vec![] },
            seq: vec![b'A'; 50], distinguishing_uniq: 0, core_bp: 0, stub: false, tes: None,
        };
        // index 0 unspliced (strand NOT measured), 1 spliced elsewhere, 2 spliced OVERLAPPING 0.
        let reps = vec![mk("c1", 1_000, 2_000, false), mk("c2", 9_000, 9_500, true), mk("c1", 1_500, 2_500, true)];
        let meta: Vec<(String, u64, u64, bool)> =
            reps.iter().map(|r| (r.chrom.clone(), r.start, r.end, r.introns.is_empty())).collect();
        let ex = |a: usize, b: usize| {
            let (x, y) = (&meta[a], &meta[b]);
            (x.3 || y.3) && !(x.0 == y.0 && x.2.min(y.2) > x.1.max(y.1))
        };
        assert!(ex(0, 1), "unmeasured strand at a DISJOINT locus: exempt");
        assert!(!ex(0, 2), "unmeasured strand but OVERLAPPING spans: NOT exempt (one locus counted twice)");
        assert!(!ex(1, 2), "both strands junction-determined: NOT exempt — a real antisense candidate");
        assert!(ex(1, 0) == ex(0, 1), "the predicate must be symmetric in its arguments");
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

    // Task 4 (augment-and-linearize): `family_linearize_cert` packages the stage-2 admission pool build +
    // `linearize_certificate` call. Fake injected realign (no minimap2) — mirrors linearize.rs's own
    // `fake_realign` test double so this stays hermetic.
    #[test]
    fn family_linearize_cert_uses_mapq0_pool() {
        use crate::vg_family::linearize::Verdict;
        // Uses a longer, non-palindromic candidate (same one as linearize.rs's
        // `real_copy_linearizes_decoy_does_not`) to minimize the chance a dinucleotide shuffle returns
        // the original. Decoys are the 20 dinucleotide shuffles only (the reverse-complement decoy was
        // removed: minimap2 is strand-symmetric, so RC(candidate) would tie `real`). With the whole
        // pool matching the candidate and no decoy beating it, perm_p = 1/(20+1) = 1/21 ~= 0.048 < 0.05
        // -> Linearizes.
        let cand = b"ACGTACGTTTGGCCAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let pool: Vec<Vec<u8>> = (0..8).map(|_| cand.clone()).collect();
        let realign = |refs: &[Vec<u8>], reads: &[Vec<u8>]| {
            let ci = refs.len() - 1;
            reads.iter().map(|r| if r == &refs[ci] { Some((ci, 60u32)) } else { Some((0, 0)) }).collect::<Vec<_>>()
        };
        let cert = family_linearize_cert(&cand, &[b"TTTTGGGGCCCCAAAAAAAATTTTGGGG".to_vec()], &pool, realign);
        assert!(matches!(cert.verdict, Verdict::Linearizes), "cert = {:?}", cert);
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
        PrimaryRead { chrom: chrom.into(), ref_start: s, ref_end: e, introns: introns.to_vec(), reverse: false }
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
    /// Physical fixture: a molecule's two records carry the SAME sequence (copyB's transcript), which is
    /// what a real BAM contains -- minimap2 `-Y` copies the full SEQ onto the secondary record, and a
    /// molecule has one sequence. The primary sits at copyA's locus, so resolving it to copyB is
    /// resolution AGAINST the primary flag.
    fn two_paralogs_with_psvs() -> (GenomeIndex, Vec<PrimaryRead>, Vec<BamRead>) {
        two_paralogs_with_psvs_impl(true)
    }

    /// ADVERSARIAL fixture: the same molecule's two records carry DIFFERENT sequences (copyB's on the
    /// primary, copyA's on the secondary). No aligner can emit this -- it is here only to force the
    /// record-vs-record contradiction that the molecule reduction must answer with an abstention.
    fn two_paralogs_contradicting_records() -> (GenomeIndex, Vec<PrimaryRead>, Vec<BamRead>) {
        two_paralogs_with_psvs_impl(false)
    }

    fn two_paralogs_with_psvs_impl(same_seq: bool) -> (GenomeIndex, Vec<PrimaryRead>, Vec<BamRead>) {
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
                    mapq: 60, name: name.into(), as_score: 380, de: 0.010, is_supplementary: false, is_secondary: false },
                BamRead { chrom: "c1".into(),
                    read: AlignedRead { ref_start: 1000, cigar: vec![('M', 200), ('N', 20), ('M', 200)], seq: secondary_seq, qual: vec![] },
                    mapq: 0, name: name.into(), as_score: 379, de: 0.012, is_supplementary: false, is_secondary: false },
            ]
        };
        let mut bam = Vec::new();
        if same_seq {
            // PHYSICAL: each molecule has ONE sequence on BOTH of its records. Three molecules carry
            // copyB's transcript and three carry copyA's -- both alleles must be observed at every PSV
            // column or the read-support filter drops the column as an assembly artifact.
            for nm in ["readB", "readC", "readD"] {
                bam.extend(mk(nm, copyb_spliced.clone(), copyb_spliced.clone()));
            }
            for nm in ["readE", "readF", "readG"] {
                bam.extend(mk(nm, copya_spliced.clone(), copya_spliced.clone()));
            }
        } else {
            for nm in ["readB", "readC", "readD"] {
                bam.extend(mk(nm, copyb_spliced.clone(), copya_spliced.clone()));
            }
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
        let (fas, fallback, _dna_needs, _linearize_certs) = detect_and_assign(
            &primary,
            &bam_reads,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
            &[],
            false,
            false,
            false,
            &fasta,
            None,
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
        let (fas, fallback, dna_needs, _linearize_certs) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
            &[],
            false,
            false,
            false,
            "",
            None,
        );
        assert!(fallback.is_empty(), "small paralogs use the exact poasta path, no fallback");
        assert!(dna_needs.is_empty(), "absent_copies=false must return empty dna_needs vec");
        assert_eq!(fas.len(), 1, "one co-located 2-copy family");
        let fa = &fas[0];
        assert_eq!(fa.n_copies, 2);
        // Six BAM RECORDS -- readB/readC/readD, each primary (locus 0) + secondary (locus 1000) -- but
        // THREE MOLECULES. The unit of an assignment is the molecule, so three rows, not six.
        assert_eq!(fa.n_reads, 6, "12 records / 6 molecules -> 6 results");
        assert_eq!(fa.assignments.len(), 6);
        // copies sorted by start: copyA=0, copyB=1. EVERY molecule's primary record sits at copyA's locus,
        // so the primary flag says "copy 0" for all six. Three of them carry copyB's transcript and must
        // come back copy 1: resolution is by SEQUENCE, against the primary flag.
        use super::super::copy_assign::AssignStatus as St;
        let to_b = fa.assignments.iter().filter(|(_, a)| a.best_copy == 1 && a.status == St::Assigned).count();
        let to_a = fa.assignments.iter().filter(|(_, a)| a.best_copy == 0 && a.status == St::Assigned).count();
        assert_eq!((to_a, to_b), (3, 3), "3 copyA molecules -> copy 0, 3 copyB molecules -> copy 1");
    }

    /// O2.14. A molecule has ONE copy of origin, so two records of one molecule that are both `Assigned`
    /// and name DIFFERENT copies cannot both be right. The reduction answers with the molecule-level
    /// abstention O2's assign-or-abstain contract requires -- it does NOT pick the stronger record, because
    /// on real data the stronger record is the PRIMARY one 323/323 times, which would make the assignment
    /// the aligner's primary flag in disguise.
    #[test]
    fn detect_and_assign_contradicting_records_of_one_molecule_abstain() {
        let (genome, primary, aligned) = two_paralogs_contradicting_records();
        let (fas, _fallback, _dna_needs, _lc) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
            &[],
            false,
            false,
            false,
            "",
            None,
        );
        assert_eq!(fas.len(), 1, "one co-located 2-copy family");
        let fa = &fas[0];
        assert_eq!(fa.assignments.len(), 3, "3 molecules, 6 records");
        for (_, a) in &fa.assignments {
            assert_eq!(
                a.status,
                super::super::copy_assign::AssignStatus::Ambiguous,
                "records naming different copies -> the MOLECULE abstains"
            );
        }
    }

    /// Final-review fix #4 (opt-in coverage): the augment-and-linearize certificate is computed for an
    /// admitted candidate ONLY when `do_linearize` is set. `linearize_cert_if_enabled` is the single guard
    /// that gates it (and, crucially, never calls `realign` when off — no minimap2 for plain --absent-copies).
    /// A fake `realign` double (that would panic if called when it must not) proves both directions
    /// hermetically, without a full end-to-end absent-copy admission.
    #[test]
    fn linearize_cert_if_enabled_is_the_opt_in_guard() {
        use crate::vg_family::linearize::Verdict;
        let cand = b"ACGTACGTTTGGCCAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
        let copy_seqs = vec![b"TTTTGGGGCCCCAAAAAAAATTTTGGGG".to_vec()];
        let pool: Vec<Vec<u8>> = (0..8).map(|_| cand.clone()).collect();

        // OFF: certificate SKIPPED (None), and `realign` is NEVER invoked (this double panics if it is).
        let panic_realign = |_: &[Vec<u8>], _: &[Vec<u8>]| -> Vec<Option<(usize, u32)>> {
            panic!("realign must not be called when do_linearize is false");
        };
        let off = linearize_cert_if_enabled(false, &cand, &copy_seqs, &pool, panic_realign);
        assert!(off.is_none(), "do_linearize=false must skip the certificate entirely");

        // ON: certificate COMPUTED (Some) — populated, using a fake realign (no minimap2).
        let ok_realign = |refs: &[Vec<u8>], reads: &[Vec<u8>]| {
            let ci = refs.len() - 1;
            reads.iter().map(|r| if r == &refs[ci] { Some((ci, 60u32)) } else { Some((0, 0)) }).collect::<Vec<_>>()
        };
        let on = linearize_cert_if_enabled(true, &cand, &copy_seqs, &pool, ok_realign);
        let cert = on.expect("do_linearize=true must compute the certificate");
        assert!(matches!(cert.verdict, Verdict::Linearizes), "cert = {cert:?}");
    }

    /// End-to-end complement: under `absent_copies=true` but `do_linearize=false`, `detect_and_assign`'s
    /// returned `linearize_certs` vec is empty (the per-candidate certificate is skipped), so plain
    /// `--absent-copies` pays no linearize cost. (This fixture admits no absent copy, so the ON vec is
    /// exercised deterministically by `linearize_cert_if_enabled_is_the_opt_in_guard` above instead.)
    #[test]
    fn detect_and_assign_do_linearize_off_returns_empty_linearize_certs() {
        let (genome, primary, aligned) = two_paralogs_with_psvs();
        let (_, _, _dna_needs, linearize_certs) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
            &[],
            true,  // absent_copies ON — the admission block runs
            false, // do_linearize OFF — but the certificate is skipped
            false, // linearize_gate OFF
            "",
            None,  // families DERIVED here (the O1->O2 --families contract is the other test)
        );
        assert!(
            linearize_certs.is_empty(),
            "do_linearize=false must yield an empty linearize_certs vec even with absent_copies=true"
        );
    }

    /// When `absent_copies=false`, the third element of `detect_and_assign`'s return tuple must
    /// always be empty — the admission block is entirely skipped, so no `DnaNeedsRecord`s are
    /// produced regardless of how many collapsed-copy candidates exist at the loci.
    #[test]
    fn detect_and_assign_absent_copies_off_returns_empty_dna_needs() {
        let (genome, primary, aligned) = two_paralogs_with_psvs();
        let (_, _, dna_needs, _linearize_certs) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
            &[],
            false, // OFF — the admission block must be completely skipped
            false,
            false,
            "",
            None,
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
        let (fas, _fallback, dna_needs, _linearize_certs) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &cfg,
            5_000_000,
            2,
            &AssignParams::default(),
            &[],
            false,
            false,
            false,
            "",
            None,
        );
        assert!(dna_needs.is_empty());
        assert_eq!(fas.len(), 1, "one co-located 2-copy family");
        let fa = &fas[0];
        assert!(fa.realign_records.is_empty(), "vg_realign OFF -> the vg_realign block must not run at all");
        assert_eq!(fa.n_copies, 2);
        assert_eq!(fa.assignments.len(), 6, "6 molecules (12 records)");
        // Same assignments the pre-Task-3 test (`detect_and_assign_resolves_multimapper_end_to_end`) pins.
        let to_b = fa.assignments.iter().filter(|(_, a)| a.best_copy == 1 && a.status == AssignStatus::Assigned).count();
        let to_a = fa.assignments.iter().filter(|(_, a)| a.best_copy == 0 && a.status == AssignStatus::Assigned).count();
        assert_eq!((to_a, to_b), (3, 3));
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
            copy_introns: vec![Vec::new(), Vec::new()],
            copy_map_identity: vec![None, None],
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
            copy_introns: vec![Vec::new(), Vec::new()],
            copy_map_identity: vec![None, None],
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
         ..Default::default() };
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
            copy_introns: vec![Vec::new()],
            copy_map_identity: vec![None],
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
                is_secondary: false,
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
            introns: vec![(30, 40)], seq: mutated.clone(),
         ..Default::default() };

        // The mocked admitter reports a remap identity (as `absent_copy::admit_candidate` does for a
        // real reference-ABSENT admission) so this test also pins that `copy_map_identity` carries it
        // through, aligned 1:1 with the widened `copy_tids`/`copy_introns` roster.
        admit_novel_pools_with_admitter(&mut fa, &pools, &bam_reads, &all_copies, &profiles, |_c, _h| {
            Admission::Copy(admitted_t.clone(), Some(0.987))
        });

        assert_eq!(fa.n_copies, 2, "the pool must be admitted as a new copy");
        assert_eq!(fa.copy_psv_alleles.len(), 2);
        assert_eq!(fa.copy_tids, vec!["host".to_string(), "admitted1".to_string()]);
        assert_eq!(fa.assignments.len(), 7, "all 7 pool reads must be assigned (none pre-existed)");
        assert_eq!(
            fa.copy_introns.len(), fa.copy_tids.len(),
            "copy_introns must stay aligned with copy_tids after admission"
        );
        assert_eq!(
            fa.copy_map_identity.len(), fa.copy_tids.len(),
            "copy_map_identity must stay aligned with copy_tids after admission"
        );
        assert_eq!(fa.copy_map_identity[0], None, "host is in-genome, never remapped");
        assert!(
            fa.copy_map_identity[1].is_some(),
            "the admitted (reference-ABSENT) copy's map identity must be captured, got {:?}",
            fa.copy_map_identity[1]
        );
        assert_eq!(fa.copy_introns[1], admitted_t.introns, "admitted copy's intron chain must propagate");

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
        let (fas, _fallback, dna_needs, _linearize_certs) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &cfg,
            5_000_000,
            2,
            &AssignParams::default(),
            &[],
            false,
            false,
            false,
            "",
            None,
        );
        assert!(dna_needs.is_empty());
        assert_eq!(fas.len(), 1);
        let fa = &fas[0];
        assert_eq!(fa.n_copies, 2, "no spurious admissions on this clean fixture");
        assert_eq!(fa.assignments.len(), 6, "6 molecules (12 records)");
        assert!(!fa.realign_records.is_empty(), "the low-MAPQ secondary records must be candidates");
        assert!(
            fa.realign_records.iter().all(|r| r.action == "reassigned" || r.action == "rejected"),
            "no unfit/novel reads on this fixture: {:?}",
            fa.realign_records
        );
        // The vg-realign path must agree with (not regress) the one-shot PSV gate's already-correct
        // calls on this fixture.
        let to_b = fa.assignments.iter().filter(|(_, a)| a.best_copy == 1 && a.status == AssignStatus::Assigned).count();
        let to_a = fa.assignments.iter().filter(|(_, a)| a.best_copy == 0 && a.status == AssignStatus::Assigned).count();
        assert_eq!((to_a, to_b), (3, 3), "the vg-realign path must not regress the PSV gate's calls");
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
        let stats = crate::vg_family::family_split::CommunityStats { n: 2, n_edges: 1, density: 1.0, avg_core_recip: 1.0, n_articulation: 0, lambda: 1 };
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
        let stats = crate::vg_family::family_split::CommunityStats { n: 2, n_edges: 1, density: 1.0, avg_core_recip: 1.0, n_articulation: 0, lambda: 1 };
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
                let au = crate::vg_family::family_graph::upper_cow(&reps[i].seq);
                let bu = crate::vg_family::family_graph::upper_cow(&reps[j].seq);
                let bru = crate::vg_family::seq_utils::reverse_complement(&bu);
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
            mapq, name: name.into(), as_score: 500, de, is_supplementary: false, is_secondary,
        }
    }

    #[test]
    fn core_bp_counts_only_bases_at_or_above_the_depth_floor() {
        // Three reads stacked over 100..200, plus one lone read running out to 1000 -- the readthrough
        // tail that must NOT enter the core, since keeping it out is the whole point of the denominator.
        let reads = vec![
            bam_read("c1", 100, 200, "a", 0.0, false),
            bam_read("c1", 100, 200, "b", 0.0, false),
            bam_read("c1", 100, 200, "c", 0.0, false),
            bam_read("c1", 150, 1000, "readthrough", 0.0, false),
        ];
        let r = rep("c1", 0, 1000);
        assert_eq!(
            locus_core_bp(&reads, std::slice::from_ref(&r), 3),
            vec![100],
            "only the 3-deep block counts; the 850 bp single-read tail is excluded"
        );
        assert_eq!(
            locus_core_bp(&reads, std::slice::from_ref(&r), 5),
            vec![0],
            "nothing reaches depth 5 -> core unmeasured -> callers must fall back to the span"
        );
    }

    #[test]
    fn confident_extent_ignores_divergent_reads_and_never_cuts_a_junction() {
        // Two clean reads define the locus; a divergent read (paralog bleed) runs 3 kb further out and
        // must NOT move the boundary.
        let mut clean_a = bam_read("c1", 1_000, 2_000, "clean_a", 0.0002, false);
        let mut clean_b = bam_read("c1", 1_100, 2_100, "clean_b", 0.0003, false);
        let bleed = bam_read("c1", 1_000, 5_000, "paralog_bleed", 0.02, false);
        clean_a.de = 0.0002;
        clean_b.de = 0.0003;
        let reads = vec![clean_a, clean_b, bleed];
        let r = rep("c1", 900, 5_200);
        let got = locus_confident_extent(&reads, std::slice::from_ref(&r), 0.0005);
        assert_eq!(got, vec![Some((1_000, 2_100))], "the 0.02-divergence read must not set the boundary");

        // A locus whose reads are ALL divergent yields None, so the caller keeps today's span.
        let only_bleed = vec![bam_read("c1", 1_000, 5_000, "b", 0.02, false)];
        assert_eq!(locus_confident_extent(&only_bleed, std::slice::from_ref(&r), 0.0005), vec![None]);

        // Clamping: a junction outside the confident-read extent still has to be contained.
        let mut spliced = rep("c1", 900, 5_200);
        spliced.introns = vec![(2_500, 3_000)];
        let got = locus_confident_extent(&reads, std::slice::from_ref(&spliced), 0.0005);
        assert_eq!(got, vec![Some((1_000, 3_000))], "must not shrink past an asserted junction");
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
         ..Default::default() }
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
         ..Default::default() }
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
        let loci = distinct_locus_reps(vec![wide.clone(), frag, para.clone()], 3);
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
        let loci = distinct_locus_reps(vec![plus, minus], 3);
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
        let loci = distinct_locus_reps(vec![dom, anti], 3);
        assert_eq!(loci.len(), 1, "minority antisense copy collapses into the real locus");
        assert_eq!(loci[0].n_reads, 666, "the dominant (real) copy is the representative");
    }

    #[test]
    fn distinct_locus_reps_unspliced_fragment_collapses_despite_balanced_reads() {
        // THE MEASURED CONTROL ARTIFACT (2026-08-08). At 4 of 7 strict single-copy control genes the
        // catalog emitted a "2-copy family" that was one locus: a spliced multi-exon rep plus a
        // SINGLE-EXON rep at the same start, whose read count equalled the locus's unspliced reads (reads
        // with no `N` in CIGAR) and whose reads sat on the SAME BAM strand. The single-exon rep's `'+'` is
        // only the default label for a junction-less rep, so the pair was routed to the antisense branch,
        // where its BALANCED read ratio (0.145-0.571, far above the 1/10 artifact bar) let it survive.
        // Read support here is deliberately balanced (60 vs 100 = 0.60) so the ANTISENSE_MINORITY_DENOM
        // rule cannot be what collapses it — only the junction-less + containment rule can.
        let mut spliced = rep_s(1000, 5000, vec![(2000, 3000)], 100);
        spliced.strand = '-';
        let mut unspliced = rep_s(1000, 5200, vec![], 60); // no introns, and deliberately OVERHANGING
        unspliced.strand = '+';
        unspliced.distinguishing_uniq = 40; // would "distinguish" under the same-strand read rule too
        let loci = distinct_locus_reps(vec![spliced, unspliced], 3);
        assert_eq!(loci.len(), 1, "an overlapping junction-less rep is the locus's unspliced fraction, not a second locus");
        assert_eq!(loci[0].n_reads, 100, "the spliced (real) copy is the representative");
    }

    #[test]
    fn distinct_locus_reps_disjoint_junctionless_copies_stay_distinct() {
        // Guard against over-collapsing: the rule requires CONTAINMENT, so two junction-less copies at
        // disjoint spans (real single-exon paralogs) remain two loci.
        let a = rep("chrX", 100, 200);
        let b = rep("chrX", 5000, 5100);
        let loci = distinct_locus_reps(vec![a, b], 3);
        assert_eq!(loci.len(), 2, "disjoint junction-less copies are distinct loci");
    }

    #[test]
    fn distinct_locus_reps_two_junctionless_copies_still_decided_by_reads() {
        // BOTH sides junction-less -> the asymmetry that identifies an unspliced fragment is absent, so
        // neither can claim priority and the read-conflict rule must still decide. This is what keeps the
        // advisor-flagged "distinguishable-but-merged" behaviour intact.
        let mut a = rep("chrX", 100, 200);
        a.distinguishing_uniq = 40;
        let b = rep("chrX", 150, 250);
        let loci = distinct_locus_reps(vec![a, b], 3);
        assert_eq!(loci.len(), 2, "two junction-less copies are decided by reads, not by the unspliced rule");
    }

    #[test]
    fn distinct_locus_reps_unspliced_rule_survives_end_overhang() {
        // The measured HMBS/TFRC shape: the junction-less rep OVERHANGS the spliced one by a base or two,
        // so any containment-based rule misses it. Overlap + junction asymmetry must still collapse.
        let spliced = rep_s(127034654, 127036245, vec![(127035000, 127035500)], 6);
        let mut unspliced = rep_s(127034654, 127036246, vec![], 19); // 1 bp past the spliced end
        unspliced.strand = '+';
        let loci = distinct_locus_reps(vec![spliced, unspliced], 3);
        assert_eq!(loci.len(), 1, "a 1 bp overhang must not defeat the unspliced-fragment collapse");
    }

    #[test]
    fn distinct_locus_reps_keeps_distinguishable_colocated_copies_separate() {
        // Task 3 (identifiability-merge): two overlapping SAME-strand copies where one carries 40
        // unique-mapper reads (well above the min_reads=3 floor) -> reads DISTINGUISH them, so the pair
        // must NOT be coordinate-collapsed even though they overlap on the same strand (the advisor-flagged
        // "distinguishable-but-merged" over-merge).
        let mut a = rep("chrX", 100, 200);
        a.distinguishing_uniq = 40;
        let b = rep("chrX", 150, 250);
        let loci = distinct_locus_reps(vec![a, b], 3);
        assert_eq!(loci.len(), 2, "a copy with 40 unique reads must not be merged");
    }

    #[test]
    fn distinct_locus_reps_still_merges_indistinguishable_colocated_copies() {
        // Same shape, but NEITHER copy carries unique-mapper support -> true K=0, still one locus.
        let a = rep("chrX", 100, 200);
        let b = rep("chrX", 150, 250);
        let loci = distinct_locus_reps(vec![a, b], 3);
        assert_eq!(loci.len(), 1, "indistinguishable co-located copies still merge (K=0)");
    }

    #[test]
    fn distinct_locus_reps_min_reads_zero_still_merges_indistinguishable() {
        // Hardening: `RUSTLE_CONFLICT_MIN_READS=0` must not invert the merge guard.
        // `reads_distinguish(0, 0, false, 0)` would be `0 >= 0` = true (spuriously "distinguishable"),
        // so `distinct_locus_reps` clamps `min_reads` to `.max(1)` before calling it. Two co-located
        // same-strand copies with zero unique-mapper support must still merge to one locus, not two.
        let mut a = DenovoTranscript { chrom: "chrX".into(), start: 100, end: 200, strand: '+', n_reads: 3, ..Default::default() };
        a.distinguishing_uniq = 0;
        let mut b = DenovoTranscript { chrom: "chrX".into(), start: 150, end: 250, strand: '+', n_reads: 3, ..Default::default() };
        b.distinguishing_uniq = 0;
        let loci = distinct_locus_reps(vec![a, b], 0);
        assert_eq!(loci.len(), 1, "min_reads=0 must still merge indistinguishable co-located copies, not invert the guard");
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

    /// The summed-coverage rule's interval merge, extracted so the arithmetic is testable without
    /// spawning minimap2. Mirrors the union loop in `nucleotide_edges`.
    fn union_len(mut iv: Vec<(u64, u64)>) -> u64 {
        iv.sort_unstable();
        let mut total = 0u64;
        let mut cur: Option<(u64, u64)> = None;
        for (s, e) in iv {
            match cur {
                Some((cs, ce)) if s <= ce => cur = Some((cs, ce.max(e))),
                Some((cs, ce)) => {
                    total += ce - cs;
                    cur = Some((s, e));
                }
                None => cur = Some((s, e)),
            }
        }
        if let Some((cs, ce)) = cur {
            total += ce - cs;
        }
        total
    }

    #[test]
    fn summed_coverage_unions_disjoint_blocks() {
        // Three separate blocks of a shattered rep: 300 + 300 + 300 = 900 of a 1500 bp shorter sequence
        // = 0.60, which clears the 0.50 floor. No single block would.
        assert_eq!(union_len(vec![(0, 300), (500, 800), (1000, 1300)]), 900);
    }

    #[test]
    fn summed_coverage_never_double_counts_overlaps() {
        // Overlapping records must UNION, not add: 0-600 and 400-900 is 900 bp of query, not 1100.
        // Adding would let a repeat aligned many times to the same locus manufacture coverage.
        assert_eq!(union_len(vec![(0, 600), (400, 900)]), 900);
        assert_eq!(union_len(vec![(0, 1000), (200, 300), (400, 500)]), 1000);
    }

    #[test]
    fn summed_coverage_single_block_matches_the_per_record_rule() {
        // With one record the summed rule must reduce to the default, so enabling it cannot change
        // any pair that already had a qualifying single alignment.
        assert_eq!(union_len(vec![(120, 900)]), 780);
        assert_eq!(union_len(vec![]), 0);
    }

    #[test]
    fn shared_exon_slices_match_the_exon_sum_layout() {
        // The rule slices exons out of the stored exon-sum instead of re-fetching, so the slice offsets
        // must reproduce the layout build_spliced_seq produced. On '-' the stored seq is the reverse
        // complement of the concatenation, so exon LENGTHS run in the opposite order.
        let lens = |t: &DenovoTranscript| -> Vec<usize> {
            let mut v = Vec::new();
            let mut prev = t.start;
            for &(d, a) in &t.introns {
                v.push(d.saturating_sub(prev) as usize);
                prev = a;
            }
            v.push(t.end.saturating_sub(prev) as usize);
            if t.strand == '-' {
                v.reverse();
            }
            v
        };
        // exons 100-200 (100 bp) and 300-450 (150 bp) => 250 bp exon-sum
        let mut t = DenovoTranscript {
            chrom: "c".into(), start: 100, end: 450, introns: vec![(200, 300)],
            seq: vec![b'A'; 250], strand: '+', ..Default::default()
        };
        assert_eq!(lens(&t), vec![100, 150]);
        assert_eq!(lens(&t).iter().sum::<usize>(), t.seq.len(), "slices must tile the exon-sum exactly");
        t.strand = '-';
        assert_eq!(lens(&t), vec![150, 100], "minus strand reverses the exon order in the stored seq");
        assert_eq!(lens(&t).iter().sum::<usize>(), t.seq.len());
    }

    #[test]
    fn tier_names_reports_every_contributing_tier() {
        assert_eq!(tier_names(0), "none");
        assert_eq!(tier_names(TIER_ASM20), "asm20");
        assert_eq!(tier_names(TIER_SENSITIVE), "sensitive");
        assert_eq!(tier_names(TIER_GENOMIC), "genomic-span");
        assert_eq!(tier_names(TIER_PROTEIN), "protein");
        // An edge found by two tiers must report BOTH -- reporting only the first would make a tier look
        // load-bearing when a second tier independently supports the same edge.
        assert_eq!(tier_names(TIER_ASM20 | TIER_SENSITIVE), "asm20+sensitive");
        assert_eq!(
            tier_names(TIER_ASM20 | TIER_SENSITIVE | TIER_GENOMIC | TIER_PROTEIN),
            "asm20+sensitive+genomic-span+protein"
        );
    }

    #[test]
    fn tier_bits_are_distinct_single_bits() {
        // count_ones()==1 is the "sole support" test in the provenance accounting; it is only meaningful if
        // every tier constant is a distinct single bit.
        for b in [TIER_ASM20, TIER_SENSITIVE, TIER_GENOMIC, TIER_PROTEIN] {
            assert_eq!(b.count_ones(), 1, "tier {b} must be a single bit");
        }
        let all = [TIER_ASM20, TIER_SENSITIVE, TIER_GENOMIC, TIER_PROTEIN];
        for (i, a) in all.iter().enumerate() {
            for b in &all[i + 1..] {
                assert_eq!(a & b, 0, "tier bits must not overlap");
            }
        }
    }

    #[test]
    fn sensitive_tier_floor_follows_min_identity_on_the_refine_path() {
        // Regression: the refine path hard-coded 0.70, which --min-identity could not reach, so
        // `--min-identity 0.98` silently admitted 0.70 edges. Both tiers must now move together.
        let p = homology_refine_params(Some(0.98), 4);
        assert_eq!(p.min_identity, 0.98);
        assert_eq!(p.sensitive_identity, 0.98, "sensitive tier must track --min-identity, not stay at 0.70");
        let d = RefineParams::default();
        assert_eq!(d.min_identity, 0.80);
        assert_eq!(d.sensitive_identity, 0.60, "default sensitive floor is 0.60 -- the EFFECTIVE E_r floor");
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
            DenovoTranscript{tid:"a".into(),chrom:"c1".into(),start:0,end:300,n_reads:5,strand:'+',introns:vec![],seq:b"ATGAAAGGGTTTTGTCCCAAAGGG".to_vec(), ..Default::default() },
            DenovoTranscript{tid:"b".into(),chrom:"c1".into(),start:9000,end:9300,n_reads:4,strand:'+',introns:vec![],seq:b"ATGAAAGGGTTTTGTCCCAAAGGG".to_vec(), ..Default::default() },
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
            stats: CommunityStats { n: 4, n_edges: 0, density: 1.0, avg_core_recip: 0.0, n_articulation: 0, lambda: 0 },
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

    // ---------------------------------------------------------------------------------------------
    // READ-LEVEL mis-chain: the gap the transcript-level rule leaves open, and the guard that keeps it shut.
    // ---------------------------------------------------------------------------------------------

    fn pread(s: u64, e: u64, introns: Vec<(u64, u64)>) -> PrimaryRead {
        PrimaryRead { chrom: "c1".into(), ref_start: s, ref_end: e, introns, reverse: false }
    }

    /// THE REGRESSION. The real `GWFAM244:2` shape, measured on the matched-individual BAM
    /// (`o3_mischain/fix/red/O1_SIDE.txt`): an 827,011 bp mis-chain carried by **97** primary reads.
    ///
    /// The transcript-level rule must return FALSE on it — 97 >= `MISCHAIN_MIN_JUNCTION_READS` (3), and
    /// that is by design, not a bug. Therefore a locus-scoped PER-READ statistic that reuses the
    /// transcript rule sees **nothing**, which is exactly how 150/150 mis-chained reads reached the O3
    /// divergence statistic and manufactured two of its three surviving candidates. The read-level rule
    /// must return TRUE on the very same object. If these two ever agree here, the read-level rule has
    /// silently acquired a support clause and this test fails.
    #[test]
    fn read_level_mischain_catches_what_the_transcript_rule_is_designed_to_miss() {
        let giant = (122_919_265u64, 123_746_276u64); // 827,011 bp, the measured mis-chain
        assert_eq!(giant.1 - giant.0, 827_011);

        let tx = rep_s(122_741_457, 123_746_276, vec![giant], 97);
        let j = junc(&[(giant.0, giant.1, 97)]); // 97 reads = 32x MISCHAIN_MIN_JUNCTION_READS
        assert!(
            !is_giant_intron_mischain(&tx, &j, MISCHAIN_GIANT_INTRON_BP, MISCHAIN_MIN_JUNCTION_READS),
            "the transcript rule is SUPPOSED to keep a well-supported giant intron; if this ever flips, \
             the read-level rule below is no longer the thing standing between mis-chains and the statistic"
        );

        let r = pread(122_741_457, 123_746_276, vec![giant]);
        assert!(
            is_mischained_read(&r, MISCHAIN_GIANT_INTRON_BP),
            "READ-LEVEL rule must fire on the same 827 kb chain regardless of how many reads share it"
        );
    }

    /// `retain_local_reads` removes exactly the reads that leave the locus, keeps ordinary transcript
    /// structure, and reports the count so a caller can print its denominator. The kept read carries the
    /// largest intron measured in a surviving O3 candidate (2,137 bp at GWFAM382:2) and a 48 kb POTE-scale
    /// intron — both under the giant threshold, both untouched.
    #[test]
    fn retain_local_reads_drops_only_the_escaping_reads_and_counts_them() {
        let mut reads = vec![
            pread(122_741_457, 122_743_192, vec![(122_741_900, 122_744_037)]), // 2,137 bp — a real intron
            pread(0, 200_000, vec![(20_000, 68_000)]),                          // 48 kb POTE — under giant
            pread(122_741_457, 123_746_276, vec![(122_919_265, 123_746_276)]),  // 827 kb — escapes
            pread(19_611_985, 20_383_478, vec![(19_611_985, 20_383_478)]),      // 771 kb — escapes
            pread(1_000, 4_000, vec![]),                                        // unspliced local
        ];
        let dropped = retain_local_reads(&mut reads);
        assert_eq!(dropped, 2, "exactly the two giant-gap reads");
        assert_eq!(reads.len(), 3);
        assert!(reads.iter().all(|r| !is_mischained_read(r, MISCHAIN_GIANT_INTRON_BP)));
    }

    /// GUARD — ONE SOURCE OF TRUTH FOR THE THRESHOLD.
    ///
    /// `MISCHAIN_GIANT_INTRON_BP` must be *defined* exactly once in the Rust sources, because the O3
    /// detector's read-selection step (`o3_detector/mischain.py`) parses that literal out of this file
    /// rather than restating it. A second definition would let the two languages drift apart silently.
    ///
    /// **Bounded scan, deliberately.** The `homology_catalog_never_touches_the_conflict_graph` guard was
    /// itself false-positived by a scan that ran to EOF and matched text outside the region it meant to
    /// police, so this one (a) enumerates `src/**/*.rs` explicitly rather than the whole tree, (b) matches
    /// only a `const` DEFINITION (`const NAME`), never a use, and (c) asserts it actually visited files, so
    /// a broken path cannot make the guard pass vacuously.
    #[test]
    fn mischain_threshold_has_exactly_one_definition_site() {
        fn rs_files(dir: &std::path::Path, out: &mut Vec<std::path::PathBuf>) {
            let Ok(rd) = std::fs::read_dir(dir) else { return };
            for e in rd.flatten() {
                let p = e.path();
                if p.is_dir() {
                    rs_files(&p, out);
                } else if p.extension().and_then(|s| s.to_str()) == Some("rs") {
                    out.push(p);
                }
            }
        }
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
        let mut files = Vec::new();
        rs_files(&root, &mut files);
        assert!(files.len() > 5, "scan found {} .rs files under {root:?} — the scan is broken, not the code", files.len());

        // The needle is ASSEMBLED AT RUNTIME. Spelling it as one literal would put the needle into this
        // very file and the guard would count ITSELF as a second definition — which is exactly what it
        // did on the first run (`fix/suite.log`: found lines 943 and 8707). Requiring `pub const` at the
        // start of the line is the same predicate `mischain.py` parses with, so the two agree by
        // construction rather than by coincidence.
        let needle = format!("pub {} {}", "const", "MISCHAIN_GIANT_INTRON_BP");
        let mut defs: Vec<String> = Vec::new();
        for f in &files {
            let Ok(txt) = std::fs::read_to_string(f) else { continue };
            for (i, line) in txt.lines().enumerate() {
                let t = line.trim_start();
                if t.starts_with("//") {
                    continue; // a doc/comment mention is a reference, not a definition
                }
                if t.starts_with(&needle) {
                    defs.push(format!("{}:{}", f.display(), i + 1));
                }
            }
        }
        assert_eq!(
            defs.len(),
            1,
            "MISCHAIN_GIANT_INTRON_BP must have exactly ONE definition site; found {defs:?}"
        );
        assert_eq!(MISCHAIN_GIANT_INTRON_BP, 50_000, "the parsed-by-Python literal changed value");
    }

    /// Integration: `maybe_salvage_mischain` before `pass1_skeletons_robust` recovers LOCAL seeding at both
    /// ends of a spurious cross-locus bridge.
    ///
    /// Locus A (chr1:1000-1100) and locus B (chr1:900_000-900_100) each have 3 native unspliced reads — on
    /// their own already enough to seed a bounded skeleton (`min_reads = 3`). Two EXTRA reads mis-chain A to
    /// B across a 898.9 kb spurious intron, an exact copy of the giant/junction shape `is_giant_intron_mischain`
    /// flags: > `MISCHAIN_GIANT_INTRON_BP` (50 kb) and carried by only 2 reads — BELOW `MISCHAIN_MIN_JUNCTION_READS`
    /// (= `GATE_MIN_READS` = 3), so it is spurious by the same criterion the assembled-transcript gate uses.
    ///
    /// OFF (measured): `pass1_skeletons_robust`'s own `(chrom, intron-chain)` grouping requires `>= min_reads`
    /// (3) reads sharing the EXACT chain; the 2 bridge reads share one chain but only number 2, so that group
    /// is silently dropped by the min-reads filter (no giant/unbounded skeleton is even produced — the bridge
    /// reads simply vanish uncounted). Locus A and B skeletons are seeded ONLY by their 3 native reads each
    /// (`n_reads == 3`).
    ///
    /// ON (measured): `split_mischained_reads` cuts each bridge read at the spurious intron into two local,
    /// intron-free segments — one at A's coordinates, one at B's — BEFORE Pass-1 groups by position
    /// (`cluster_unspliced`). Those 2+2 salvaged segments join their native locus's cluster (identical span),
    /// so both A's and B's skeleton gain real read support: `n_reads` goes from 3 to 5 at EACH locus. This is
    /// the load-bearing ON-vs-OFF distinction (not merely "a skeleton exists", which is already true OFF) —
    /// the test fails if the split is a no-op.
    #[test]
    fn salvage_seeds_local_locus_that_a_mischain_bridge_would_lose() {
        use crate::vg_family::denovo_assemble::{pass1_skeletons_robust, PrimaryRead, Skeleton};
        let mk = |s: u64, e: u64, introns: Vec<(u64, u64)>| PrimaryRead { chrom: "chr1".into(), ref_start: s, ref_end: e, introns, reverse: false };
        let mut reads = vec![];
        for _ in 0..3 { reads.push(mk(1000, 1100, vec![])); } // locus A natives
        for _ in 0..3 { reads.push(mk(900_000, 900_100, vec![])); } // locus B natives
        // 2 (not 3!) spurious bridge reads: their shared junction (1100, 900_000) is then supported by
        // exactly 2 reads -- BELOW MISCHAIN_MIN_JUNCTION_READS (3), so `is_giant_intron_mischain`'s own
        // criterion (and split_mischained_reads', which reuses it) flags it as spurious. At 3 reads it would
        // sit AT the gate (not "fewer than") and never be cut -- see `mischain_keeps_real_large_gene_with_well_supported_intron`.
        for _ in 0..2 { reads.push(mk(1000, 900_100, vec![(1100, 900_000)])); }

        let cfg_off = DenovoConfig::default();
        let cfg_on = DenovoConfig { mischain_salvage: true, ..DenovoConfig::default() };

        // OFF must not transform reads at all.
        assert!(maybe_salvage_mischain(&reads, &cfg_off).is_none(), "OFF must not transform reads");

        let reads_off = maybe_salvage_mischain(&reads, &cfg_off).unwrap_or_else(|| reads.clone());
        let reads_on = maybe_salvage_mischain(&reads, &cfg_on).unwrap_or_else(|| reads.clone());
        let sk_off = pass1_skeletons_robust(&reads_off, 3, 1);
        let sk_on = pass1_skeletons_robust(&reads_on, 3, 1);

        let bounded_at = |sk: &[Skeleton], lo: u64, hi: u64| -> Option<u32> {
            sk.iter()
                .find(|s| s.chrom == "chr1" && s.start <= hi && s.end >= lo && s.end - s.start < 300_000)
                .map(|s| s.n_reads)
        };

        // OFF: each locus seeded by its 3 native reads only; the bridge reads contribute nothing (their
        // 2-read chain never reaches the skeleton min-reads gate).
        assert_eq!(bounded_at(&sk_off, 1000, 1100), Some(3), "OFF: locus A seeded only by native reads");
        assert_eq!(bounded_at(&sk_off, 900_000, 900_100), Some(3), "OFF: locus B seeded only by native reads");
        assert_eq!(sk_off.len(), 2, "OFF: no phantom/giant skeleton from the sub-gate bridge chain");

        // ON: the salvaged local segments join each locus's native cluster -- strictly MORE read support at
        // BOTH ends, the real invariant (not just "a skeleton exists", which is already true OFF).
        assert!(bounded_at(&sk_on, 900_000, 900_100).unwrap_or(0) >= 1);
        assert_eq!(bounded_at(&sk_on, 1000, 1100), Some(5), "ON: locus A gains the 2 salvaged local segments");
        assert_eq!(bounded_at(&sk_on, 900_000, 900_100), Some(5), "ON: locus B gains the 2 salvaged local segments");
        assert!(
            bounded_at(&sk_on, 1000, 1100).unwrap() > bounded_at(&sk_off, 1000, 1100).unwrap(),
            "salvage must add strictly more bounded local support at A, not merely re-derive OFF's count"
        );
        assert!(
            bounded_at(&sk_on, 900_000, 900_100).unwrap() > bounded_at(&sk_off, 900_000, 900_100).unwrap(),
            "salvage must add strictly more bounded local support at B, not merely re-derive OFF's count"
        );
        assert_eq!(sk_on.len(), 2, "ON: still exactly the two real loci, no phantom bridge skeleton");
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
            is_secondary: false,
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
         ..Default::default() };
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

    // NOTE: mutates the process-global RUSTLE_COLLAPSE_ENUMERATE env var. `serial_test` is not a
    // dev-dependency of this crate (checked Cargo.toml), so this cannot be marked `#[serial]` without
    // adding a new dependency; no other test in this crate reads/writes this var, so the race window
    // is currently theoretical, but a future env-mutating test on the same var would need coordinating.
    #[test]
    fn collapse_enumerate_defaults_off_and_reads_env() {
        let d = DenovoConfig::default();
        assert!(!d.collapse_enumerate, "must default OFF for byte-identical behavior");
        std::env::set_var("RUSTLE_COLLAPSE_ENUMERATE", "1");
        assert!(DenovoConfig::from_env().collapse_enumerate);
        std::env::remove_var("RUSTLE_COLLAPSE_ENUMERATE");
        assert!(!DenovoConfig::from_env().collapse_enumerate);
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
            is_secondary: false,
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

    #[test]
    fn tied_existence_families_appends_only_non_overlapping_as_tied() {
        // The amylase/os1 fix: tied reps become existence-only TSFAM families ONLY at loci no emitted copy
        // occupies, and every read on them abstains (Tied). Overlapping tied reps are dropped.
        let mk = |start: u64| BamRead {
            chrom: "c1".into(),
            read: AlignedRead { ref_start: start, cigar: vec![('M', 100)], seq: vec![b'A'; 100], qual: vec![30; 100] },
            mapq: 0,
            name: format!("r{start}"),
            as_score: 100,
            de: 0.01,
            is_supplementary: false,
            is_secondary: false,
        };
        let tied = vec![rep_s(1_000, 1_500, vec![], 5), rep_s(5_000, 5_500, vec![], 5)];
        let emitted = vec![("c1".to_string(), 900u64, 1_100u64)]; // overlaps the first tied rep only
        let reads = vec![mk(1_050), mk(5_100)];
        let fams = tied_existence_families(&tied, &emitted, &reads, 3);
        assert_eq!(fams.len(), 1, "only the tied rep at the free locus is appended");
        assert_eq!(fams[0].family_id, "TSFAM3", "id continues from next_family_index");
        assert_eq!(fams[0].n_copies, 1, "singleton existence copy");
        assert_eq!(fams[0].copy_spans, vec![("c1".to_string(), 5_000, 5_500)]);
        assert!(!fams[0].assignments.is_empty(), "the overlapping read is tied to it");
        for (_, a) in &fams[0].assignments {
            assert_eq!(a.status, AssignStatus::Tied, "existence-only: every read abstains (K=0)");
        }
    }

    /// Junction support is PRIMARY-only by type: `PrimaryRead` cannot hold a secondary alignment
    /// (`primary_read_from_record` returns `None` for them). The predecessor took `&[BamRead]` and counted
    /// secondaries, inflating the DAZ readthrough span from 56 to 154 distinct junctions.
    #[test]
    fn read_junction_support_counts_each_primary_read_once() {
        let pr = |s: u64, introns: Vec<(u64, u64)>| PrimaryRead { chrom: "c1".into(), ref_start: s, ref_end: s + 400, introns, reverse: false };
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
         ..Default::default() };
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
        let mk = |s: u64| PrimaryRead { chrom: "c1".into(), ref_start: s, ref_end: 500, introns: vec![(200, 300)], reverse: false };
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
            DenovoTranscript { tid: "u0".into(), chrom: "c1".into(), start: 0, end: 900, n_reads: 10, strand: '+', introns: vec![], seq: base, ..Default::default() },
            DenovoTranscript { tid: "u1".into(), chrom: "c2".into(), start: 0, end: 900, n_reads: 10, strand: '+', introns: vec![], seq: para, ..Default::default() },
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
            DenovoTranscript { tid: "r0".into(), chrom: "c1".into(), start: 0, end: 900, n_reads: 10, strand: '+', introns: vec![], seq: base, ..Default::default() },
            DenovoTranscript { tid: "r1".into(), chrom: "c1".into(), start: 5000, end: 5900, n_reads: 8, strand: '+', introns: vec![], seq: para, ..Default::default() },
            DenovoTranscript { tid: "r2".into(), chrom: "c1".into(), start: 9000, end: 9900, n_reads: 5, strand: '+', introns: vec![], seq: rand_seq(900, 0x99), ..Default::default() },
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
            &[seq_a.clone(), seq_b.clone()], ER_SENSITIVE_SEED, 0.60, 0.50, None, &p,
        ).unwrap();
        assert!(sensitive_edges.is_empty(), "pair must be genuinely nt-unresolvable (< 0.60), got {:?}", sensitive_edges);

        let reps = vec![
            DenovoTranscript { tid: "protA".into(), chrom: "c1".into(), start: 1_000, end: 1_750, n_reads: 20, strand: '+', introns: vec![], seq: seq_a, ..Default::default() },
            DenovoTranscript { tid: "protB".into(), chrom: "c1".into(), start: 50_000, end: 50_750, n_reads: 18, strand: '+', introns: vec![], seq: seq_b, ..Default::default() },
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
        let (fams, certs, _collapsed, _expressed, _dna) = detect_homology_catalog_genome_wide(
            "tests/fixtures/same_chrom_supplement/reads.bam",
            "tests/fixtures/same_chrom_supplement/genome.fa",
            2, 2, &DenovoConfig::default(), &RefineParams::default(), 0.20,
        ).unwrap();
        // the fixture's two homologous loci (c1:A + c2:X) must land in one family of >= 2 distinct loci.
        assert!(fams.iter().any(|f| f.len() >= 2), "expected a >=2-copy homology family, got {:?}", fams.iter().map(|f| f.len()).collect::<Vec<_>>());
        // ONE CERTIFICATE PER EMITTED FAMILY, SAME ORDER -- the column would otherwise attach each
        // family's lambda to a different family's row.
        assert_eq!(certs.len(), fams.len(), "one certificate per emitted family");
        for (c, f) in certs.iter().zip(fams.iter()) {
            assert_eq!(c.n, f.len(), "certificate node count must be the EMITTED copy count");
        }
    }

    // ── THE λ CERTIFICATE IS COMPUTED ON THE EMITTED NODE SET ────────────────────────────────────
    //
    // `distinct_locus_reps` MERGES co-located copies, so the emitted family has fewer nodes than the
    // block it came from. A λ measured before that merge describes a different object than the row the
    // catalog prints. This is the documented trap "never judge a change to what a NODE IS on node-level
    // metrics", and the test below is built so the two answers genuinely DISAGREE -- if `certificate_for`
    // were ever changed to take the pre-merge block, this test fails rather than drifting silently.

    #[test]
    fn certificate_is_measured_after_the_locus_merge_not_before() {
        // Pre-merge: reps {0,1,2} form a TRIANGLE -> lambda = 2 (2-edge-connected).
        let er_edges = [(0usize, 1usize), (0, 2), (1, 2)];
        assert_eq!(
            crate::vg_family::family_split::edge_connectivity(3, &er_edges),
            2,
            "the pre-merge block is a triangle"
        );
        // Emitted: reps 0 and 1 collapsed into ONE locus, so the family has 2 copies joined by the
        // single edge "something in {0,1} aligns to 2" -> lambda = 1, NOT 2.
        let groups = vec![vec![0usize, 1usize], vec![2usize]];
        let cert = certificate_for(&groups, &er_edges);
        assert_eq!(cert.n, 2, "two EMITTED copies, not three reps");
        assert_eq!(cert.n_edges, 1, "(0,2) and (1,2) are the SAME emitted edge and must not double-count");
        assert_eq!(cert.lambda, 1, "measured on the emitted node set, the family hangs on one edge");
        assert!((cert.density - 1.0).abs() < 1e-9, "2 nodes, 1 edge => density 1.0");
    }

    #[test]
    fn certificate_drops_edges_internal_to_a_merged_copy() {
        // The (0,1) edge became a SELF-LOOP when those reps collapsed into one locus: it cannot hold two
        // emitted copies together, so it must not be counted as evidence that it does.
        let er_edges = [(0usize, 1usize)];
        let cert = certificate_for(&[vec![0, 1]], &er_edges);
        assert_eq!(cert.n, 1);
        assert_eq!(cert.n_edges, 0, "an intra-copy edge is not a family edge");
        assert_eq!(cert.lambda, 0, "a 1-node family has no cut to pay");
    }

    #[test]
    fn certificate_two_copy_family_is_always_lambda_one() {
        // THE REASON λ IS REPORTED AND NOT ENFORCED. A 2-copy family cannot clear `lambda >= 2`, so a
        // membership rule built on it would delete the single most common family size in the catalog.
        let cert = certificate_for(&[vec![0], vec![1]], &[(0usize, 1usize)]);
        assert_eq!(cert.n, 2);
        assert_eq!(cert.lambda, 1);
        assert!(cert.lambda < 2, "cut_certified=false here is NOT a defect flag");
    }
}
