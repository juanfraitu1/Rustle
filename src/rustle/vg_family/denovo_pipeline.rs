//! De-novo family detection DRIVER (integration stage 3): chains the ported cores into a family roster.
//!
//!   primary reads ─► pass1 skeletons ─► assemble gate ─► collapse loci ─► detect edges ─► decompose families
//!
//! This is the read-coherence-way detection pipeline end to end (rescue + per-read copy assignment are the
//! next stage). `detect_families` is the testable transform over parsed reads + a loaded genome;
//! `detect_families_from_files` is the thin I/O wrapper (BAM + FASTA).

use std::collections::{BTreeSet, HashSet};

use anyhow::Result;

use super::copy_assign::{Assignment, AssignParams, AssignStatus};
use super::copy_assign_pipeline::{assign_family_detailed, read_ref_end};
use super::copy_split::{split_locus_copies, AlignedRead};
use super::denovo_assemble::{
    assemble_gate, pass1_skeletons, primary_reads_from_bam, BamRead, GateParams, PrimaryRead, PASS1_MIN_READS,
};
use super::family_detect::{collapse_loci, detect_edges, detect_edges_reporting, DenovoTranscript, DetectParams};
use super::family_rescue::{FamilyMember, RescueParams};
use super::family_split::{decompose_families, FamilyClass, SplitFamily, SplitParams};
use super::rescue_pipeline::{rescue_thin_loci_iterative, thin_loci, MemberSpan, RESCUE_MIN_SUPPORT};
use crate::genome::GenomeIndex;

/// Configuration for the de-novo detection pipeline (defaults mirror the python stages).
#[derive(Clone, Copy, Debug)]
pub struct DenovoConfig {
    pub pass1_min_reads: u32,
    pub gate: GateParams,
    pub detect: DetectParams,
    pub split: SplitParams,
}

impl Default for DenovoConfig {
    fn default() -> Self {
        DenovoConfig {
            pass1_min_reads: PASS1_MIN_READS,
            gate: GateParams::default(),
            detect: DetectParams::default(),
            split: SplitParams::default(),
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
    let rep_idx = collapse_loci(&transcripts);
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
pub fn colocated_families(
    reps: &[DenovoTranscript],
    families: &[SplitFamily],
    win: u64,
    min_copies: usize,
) -> Vec<ColocatedFamily> {
    let mut out = Vec::new();
    for (fi, fam) in families.iter().enumerate() {
        if fam.class != FamilyClass::Family {
            continue;
        }
        let mut by_chrom: std::collections::BTreeMap<&str, Vec<usize>> = std::collections::BTreeMap::new();
        for &m in &fam.members {
            by_chrom.entry(reps[m].chrom.as_str()).or_default().push(m);
        }
        for (chrom, mut idxs) in by_chrom {
            if idxs.len() < min_copies {
                continue;
            }
            idxs.sort_by_key(|&m| reps[m].start);
            let span = reps[*idxs.last().unwrap()].start - reps[idxs[0]].start;
            if span <= win {
                let copies: Vec<DenovoTranscript> = idxs.iter().map(|&m| reps[m].clone()).collect();
                let start = copies.iter().map(|c| c.start).min().unwrap();
                let end = copies.iter().map(|c| c.end).max().unwrap();
                out.push(ColocatedFamily { family_id: format!("DSFAM{fi}"), chrom: chrom.to_string(), start, end, copies });
            }
        }
    }
    out
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

/// Collapsed-copy recovery PAST the family gate. The genuinely collapsed tandem arrays (DAZ-type) don't form
/// a co-located family to rescue into, so this runs the PSV copy-split directly on EACH rep's overlapping
/// reads (`bam_reads` already include the secondary multimappers). Returns the number of EXTRA
/// PSV-DISTINGUISHABLE copies found beyond one-per-rep. The phantom safeguard is intrinsic: a fully-tied
/// (identical) locus does not split (`split_locus_copies` requires copy-specific PSVs), so only copies with
/// real distinguishing evidence are counted — the identifiability gate the user asked for, for free.
fn recover_collapsed_copies(reps: &[DenovoTranscript], bam_reads: &[BamRead]) -> usize {
    let mut recovered = 0usize;
    for rep in reps {
        let reads: Vec<AlignedRead> = bam_reads
            .iter()
            .filter(|br| br.chrom == rep.chrom && br.read.ref_start < rep.end && read_ref_end(&br.read) > rep.start)
            .map(|br| br.read.clone())
            .collect();
        if reads.len() < 6 {
            continue;
        }
        let copies = split_locus_copies(&reads, 3, 2, 3);
        if copies.len() >= 2 {
            recovered += copies.len() - 1;
        }
    }
    recovered
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
) -> (Vec<FamilyAssignment>, Vec<FallbackEdge>) {
    let skeletons = pass1_skeletons(primary_reads, cfg.pass1_min_reads);
    let transcripts = assemble_gate(&skeletons, genome, &cfg.gate);
    let rep_idx = collapse_loci(&transcripts);
    let reps: Vec<DenovoTranscript> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
    eprintln!(
        "[detect_and_assign] {} primary -> {} skeletons -> {} transcripts -> {} reps",
        primary_reads.len(),
        skeletons.len(),
        transcripts.len(),
        reps.len()
    );
    if !rescue_extra.is_empty() {
        let rec = recover_collapsed_copies(&reps, bam_reads);
        eprintln!(
            "[detect_and_assign] rescue_extra (AS-tied secondaries): {} | collapsed copies recovered past family gate (PSV-distinguishable): {}",
            rescue_extra.len(),
            rec
        );
    }
    let (edges, fallback_pairs) = detect_edges_reporting(&reps, &cfg.detect);
    eprintln!(
        "[detect_and_assign] {} homology edges ({} via large-seq fallback)",
        edges.len(),
        fallback_pairs.len()
    );
    let fallback: Vec<FallbackEdge> = fallback_pairs
        .iter()
        .map(|&(a, b)| FallbackEdge {
            chrom: reps[a].chrom.clone(),
            tid_a: reps[a].tid.clone(),
            start_a: reps[a].start,
            end_a: reps[a].end,
            len_a: reps[a].seq.len(),
            tid_b: reps[b].tid.clone(),
            start_b: reps[b].start,
            end_b: reps[b].end,
            len_b: reps[b].seq.len(),
        })
        .collect();
    let split = decompose_families(&edges, &cfg.split);
    let mut out = Vec::new();
    for cf in colocated_families(&reps, &split, win, min_copies) {
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
        let copies: Vec<&DenovoTranscript> = all_copies.iter().collect();
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
        let detail = assign_family_detailed(&copies, &region, p);
        if detail.mosaic_reads > 0 || !detail.conversions.is_empty() {
            eprintln!(
                "[mosaic] {} {}: {} mosaic-path reads, {} conversion events",
                cf.family_id, cf.chrom, detail.mosaic_reads, detail.conversions.len()
            );
        }
        // collapsed-copy recovery: group the reads by their mapped copy/locus and split each by within-locus
        // PSV haplotype; >= 2 identifiable copies at one locus means extra (collapsed) copies were merged.
        let mut by_copy: Vec<Vec<AlignedRead>> = vec![Vec::new(); copies.len()];
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
            fa.assignments.push((idx_map[r.read_index], r.combined));
        }
        out.push(fa);
    }
    (out, fallback)
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

#[cfg(test)]
mod tests {
    use super::super::copy_split::AlignedRead;
    use super::*;

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
        let ar = AlignedRead { ref_start: 0, cigar: vec![('M', 200), ('N', 20), ('M', 200)], seq: copyb_spliced };
        let bam = vec![BamRead { chrom: "c1".into(), read: ar, mapq: 60, name: "readB".into() }];
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
        let (fas, fallback) = detect_and_assign(
            &primary,
            &bam_reads,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
            &[],
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
        let (fas, fallback) = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
            &[],
        );
        assert!(fallback.is_empty(), "small paralogs use the exact poasta path, no fallback");
        assert_eq!(fas.len(), 1, "one co-located 2-copy family");
        let fa = &fas[0];
        assert_eq!(fa.n_copies, 2);
        assert_eq!(fa.n_reads, 1);
        assert_eq!(fa.assignments.len(), 1);
        let (_, a) = &fa.assignments[0];
        // copies sorted by start: copyA=0, copyB=1. The copyB read (aligned to copyA's region) -> copyB.
        assert_eq!(a.best_copy, 1, "multimapper resolved to its true copy (copyB)");
        assert_eq!(a.status, super::super::copy_assign::AssignStatus::Assigned);
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
        let colo = colocated_families(&reps, &split, 5_000_000, 2);
        assert_eq!(colo.len(), 1);
        assert_eq!(colo[0].copies.len(), 2);
        assert_eq!(colo[0].chrom, "c1");
        // same family but min_copies=3 -> not co-located (only 2 copies)
        assert!(colocated_families(&reps, &split, 5_000_000, 3).is_empty());
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
}
