//! De-novo family detection DRIVER (integration stage 3): chains the ported cores into a family roster.
//!
//!   primary reads ─► pass1 skeletons ─► assemble gate ─► collapse loci ─► detect edges ─► decompose families
//!
//! This is the read-coherence-way detection pipeline end to end (rescue + per-read copy assignment are the
//! next stage). `detect_families` is the testable transform over parsed reads + a loaded genome;
//! `detect_families_from_files` is the thin I/O wrapper (BAM + FASTA).

use std::collections::{BTreeSet, HashSet};

use anyhow::Result;

use super::copy_assign::{AssignParams, AssignStatus};
use super::copy_assign_pipeline::{assign_family, read_ref_end};
use super::copy_split::AlignedRead;
use super::denovo_assemble::{
    assemble_gate, pass1_skeletons, primary_reads_from_bam, GateParams, PrimaryRead, PASS1_MIN_READS,
};
use super::family_detect::{collapse_loci, detect_edges, DenovoTranscript, DetectParams};
use super::family_split::{decompose_families, FamilyClass, SplitFamily, SplitParams};
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

/// Per-co-located-family read-assignment summary.
#[derive(Clone, Debug)]
pub struct FamilyAssignment {
    pub family_id: String,
    pub chrom: String,
    pub n_copies: usize,
    pub n_reads: usize,
    pub n_resolvable: usize,
    pub n_assigned: usize,
    pub n_tied: usize,
    /// `(index into aligned_reads, assignment)` for each read over the family.
    pub assignments: Vec<(usize, super::copy_assign::Assignment)>,
}

/// END-TO-END pipeline: detect families, then for each co-located family assign every read overlapping it to
/// a copy. `aligned_reads` are `(chrom, read)` (chrom for region filtering). This is the runnable
/// detection + per-read copy-assignment pipeline.
#[allow(clippy::too_many_arguments)]
pub fn detect_and_assign(
    primary_reads: &[PrimaryRead],
    aligned_reads: &[(String, AlignedRead)],
    genome: &GenomeIndex,
    cfg: &DenovoConfig,
    win: u64,
    min_copies: usize,
    p: &AssignParams,
) -> Vec<FamilyAssignment> {
    let skeletons = pass1_skeletons(primary_reads, cfg.pass1_min_reads);
    let transcripts = assemble_gate(&skeletons, genome, &cfg.gate);
    let rep_idx = collapse_loci(&transcripts);
    let reps: Vec<DenovoTranscript> = rep_idx.iter().map(|&i| transcripts[i].clone()).collect();
    let edges = detect_edges(&reps, &cfg.detect);
    let split = decompose_families(&edges, &cfg.split);
    let mut out = Vec::new();
    for cf in colocated_families(&reps, &split, win, min_copies) {
        let copies: Vec<&DenovoTranscript> = cf.copies.iter().collect();
        // reads on this family's chrom overlapping its span (assign_family overlaps by coord, so pre-filter chrom).
        let mut idx_map = Vec::new();
        let mut region = Vec::new();
        for (i, (chrom, read)) in aligned_reads.iter().enumerate() {
            if chrom == &cf.chrom && read.ref_start < cf.end && read_ref_end(read) > cf.start {
                idx_map.push(i);
                region.push(read.clone());
            }
        }
        let assigns = assign_family(&copies, &region, p);
        let (mut n_resolvable, mut n_assigned, mut n_tied) = (0usize, 0usize, 0usize);
        let mut assignments = Vec::with_capacity(assigns.len());
        for (local_ri, a) in assigns {
            match a.status {
                AssignStatus::Tied => n_tied += 1,
                AssignStatus::Assigned => {
                    n_resolvable += 1;
                    n_assigned += 1;
                }
                AssignStatus::Ambiguous => n_resolvable += 1,
            }
            assignments.push((idx_map[local_ri], a));
        }
        out.push(FamilyAssignment {
            family_id: cf.family_id,
            chrom: cf.chrom,
            n_copies: cf.copies.len(),
            n_reads: assignments.len(),
            n_resolvable,
            n_assigned,
            n_tied,
            assignments,
        });
    }
    out
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
    fn two_paralogs_with_psvs() -> (GenomeIndex, Vec<PrimaryRead>, Vec<(String, AlignedRead)>) {
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
        (genome, primary, vec![("c1".to_string(), ar)])
    }

    #[test]
    fn detect_and_assign_resolves_multimapper_end_to_end() {
        let (genome, primary, aligned) = two_paralogs_with_psvs();
        let fas = detect_and_assign(
            &primary,
            &aligned,
            &genome,
            &DenovoConfig::default(),
            5_000_000,
            2,
            &super::super::copy_assign::AssignParams::default(),
        );
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
}
