//! Leave-one-out HMM rescue experiment for GOLGA6L7.
//!
//! For each held-out paralog c in the family:
//!   1. Build a family primitive (FamilyGraph + per-exon ProfileHmm) from
//!      the OTHER paralogs only.
//!   2. Sample reads aligned within the held-out paralog's genomic region
//!      from GGO_19.bam (positive set — these are the reads that would
//!      be unmapped if the held-out paralog were missing from the genome).
//!   3. Sample random non-family reads from a distant region (negative).
//!   4. Score both populations with `forward_against_family`.
//!   5. Report median, 90th percentile, max for each population, and the
//!      separation positive_max − random_max.
//!
//! Outcome interpretation:
//!   - If positive scores >> random scores, the family primitive can
//!     identify the held-out paralog's reads from sequence alone — i.e.,
//!     the rescue mechanism (vg_hmm/rescue.rs) would find them in the
//!     unmapped pool of a real sample.
//!   - The separation as a function of paralog identity to in-primitive
//!     paralogs quantifies the divergence ceiling of HMM rescue.

use std::collections::HashMap;
use std::sync::Arc;

use rustle::genome::GenomeIndex;
use rustle::types::{Bundle, BundleRead, Junction, JunctionStat, JunctionStats};
use rustle::vg::FamilyGroup;
use rustle::vg_hmm::family_graph::{build_family_graph, fit_profiles_in_place};
use rustle::vg_hmm::scorer::forward_against_family;

const CHROM: &str = "NC_073243.2";

const COPY0_EXONS: [(u64, u64); 9] = [
    (104789646, 104789918), (104791223, 104791352), (104792166, 104792196),
    (104792465, 104792516), (104792604, 104792698), (104792780, 104792887),
    (104794129, 104794188), (104794517, 104794659), (104794953, 104796276),
];
const COPY1_EXONS: [(u64, u64); 9] = [
    (104830535, 104830737), (104832042, 104832171), (104832985, 104833015),
    (104833284, 104833335), (104833423, 104833517), (104833599, 104833706),
    (104834952, 104835011), (104835340, 104835482), (104835776, 104837094),
];
const COPY2_EXONS: [(u64, u64); 9] = [
    (104871355, 104871545), (104872850, 104872979), (104873793, 104873823),
    (104874092, 104874143), (104874231, 104874325), (104874407, 104874514),
    (104875758, 104875817), (104876146, 104876288), (104876582, 104877901),
];

const PARALOG_LABELS: [&str; 3] = ["COPY0/LOC115931294", "COPY1/LOC134757625", "COPY2/LOC101137218"];

fn paralog_exons(c: usize) -> &'static [(u64, u64)] {
    match c { 0 => &COPY0_EXONS, 1 => &COPY1_EXONS, 2 => &COPY2_EXONS, _ => unreachable!() }
}

fn build_junction_stats(exons: &[(u64, u64)]) -> JunctionStats {
    let mut s: JunctionStats = HashMap::with_hasher(fxhash::FxBuildHasher::default());
    for i in 0..exons.len().saturating_sub(1) {
        s.insert(Junction { donor: exons[i].1, acceptor: exons[i + 1].0 }, JunctionStat::default());
    }
    s
}

fn make_bundle(exons: &[(u64, u64)]) -> Bundle {
    let start = exons.iter().map(|e| e.0).min().unwrap();
    let end = exons.iter().map(|e| e.1).max().unwrap();
    let read = BundleRead {
        read_uid: 1, read_name: Arc::from("synthetic"), read_name_hash: 0,
        ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
        ref_start: start, ref_end: end, exons: exons.to_vec(),
        junctions: Vec::new(), junction_valid: Vec::new(),
        junctions_raw: Vec::new(), junctions_del: Vec::new(),
        weight: 1.0, is_reverse: false, strand: '+',
        has_poly_start: false, has_poly_end: false,
        has_poly_start_aligned: false, has_poly_start_unaligned: false,
        has_poly_end_aligned: false, has_poly_end_unaligned: false,
        unaligned_poly_t: 0, unaligned_poly_a: 0,
        has_last_exon_polya: false, has_first_exon_polyt: false,
        query_length: None, clip_left: 0, clip_right: 0,
        nh: 1, nm: 0, md: None, insertion_sites: Vec::new(),
        unitig: false, unitig_cov: 0.0, read_count_yc: 0.0,
        countfrag_len: 0.0, countfrag_num: 0.0, junc_mismatch_weight: 0.0,
        pair_idx: Vec::new(), pair_count: Vec::new(),
        mapq: 0, mismatches: Vec::new(),
        hp_tag: None, ps_tag: None, is_primary_alignment: true,
    };
    Bundle {
        chrom: CHROM.to_string(), start, end, strand: '+',
        reads: vec![read],
        junction_stats: build_junction_stats(exons),
        bundlenodes: None, read_bnodes: None, bnode_colors: None,
        synthetic: false, rescue_class: None,
    }
}

fn decode_seq(seq_obj: &noodles_bam::record::Sequence) -> Vec<u8> {
    let raw = seq_obj.as_ref();
    let len = seq_obj.len();
    let mut out: Vec<u8> = Vec::with_capacity(len);
    for (i, &byte) in raw.iter().enumerate() {
        let nibs = [(byte >> 4) & 0xF, byte & 0xF];
        for (j, n) in nibs.iter().enumerate() {
            if i * 2 + j >= len { break; }
            out.push(match n { 1 => b'A', 2 => b'C', 4 => b'G', 8 => b'T', _ => b'N' });
        }
    }
    out
}

/// Macro-style closure: caller passes a BAM reader and we collect sequences
/// in the given region. Inlined to avoid the noodles type-bound dance for a
/// standalone fn.
macro_rules! collect_reads {
    ($reader:expr, $header:expr, $region_str:expr, $cap:expr, $min_len:expr) => {{
        let region = $region_str.parse()?;
        let q = $reader.query(&$header, &region)?;
        let mut out: Vec<Vec<u8>> = Vec::new();
        for r in q {
            let rec = r?;
            let f = rec.flags();
            if !f.is_unmapped() && !f.is_secondary() && !f.is_supplementary() {
                let seq = decode_seq(&rec.sequence());
                if seq.len() >= $min_len {
                    out.push(seq);
                    if out.len() >= $cap { break; }
                }
            }
        }
        out
    }}
}

fn percentile(sorted: &[f64], pct: f64) -> f64 {
    if sorted.is_empty() { return f64::NEG_INFINITY; }
    let idx = ((sorted.len() as f64 * pct / 100.0).ceil() as usize).saturating_sub(1);
    sorted[idx.min(sorted.len() - 1)]
}


fn run_loo_for_paralog(
    held_out: usize,
    genome: &GenomeIndex,
    pos_reads: &[Vec<u8>],
    rnd_reads: &[Vec<u8>],
) -> anyhow::Result<()> {
    let other_ids: Vec<usize> = (0..3).filter(|&i| i != held_out).collect();
    let bundles: Vec<Bundle> = other_ids.iter().map(|&i| make_bundle(paralog_exons(i))).collect();
    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: (0..bundles.len()).collect(),
        multimap_reads: HashMap::new(),
    };
    let mut fg = build_family_graph(&family, &bundles, Some(genome), 0.30, 0.30)?;
    let n_shared = fg.nodes.iter().filter(|n| !n.copy_specific).count();
    fit_profiles_in_place(&mut fg)?;

    let _ = genome;  // genome used by build_family_graph internally
    // Score read populations.
    let mut pos_scores: Vec<f64> = pos_reads.iter().map(|r| forward_against_family(&fg, r)).collect();
    let mut rnd_scores: Vec<f64> = rnd_reads.iter().map(|r| forward_against_family(&fg, r)).collect();
    pos_scores.sort_by(|a, b| a.partial_cmp(b).unwrap());
    rnd_scores.sort_by(|a, b| a.partial_cmp(b).unwrap());

    println!();
    println!("──── HOLD OUT  {}  (built primitive from {} + {}) ────",
        PARALOG_LABELS[held_out],
        PARALOG_LABELS[other_ids[0]], PARALOG_LABELS[other_ids[1]]);
    println!("  primitive: {} nodes ({} shared, {} copy-specific), {} edges",
        fg.n_nodes(), n_shared, fg.n_nodes() - n_shared, fg.n_edges());

    println!("  positive (held-out paralog region, n={}):", pos_scores.len());
    println!("    median={:.2}   90th={:.2}   max={:.2}",
        percentile(&pos_scores, 50.0),
        percentile(&pos_scores, 90.0),
        pos_scores.last().copied().unwrap_or(f64::NEG_INFINITY));

    println!("  random (non-family region, n={}):", rnd_scores.len());
    println!("    median={:.2}   90th={:.2}   max={:.2}",
        percentile(&rnd_scores, 50.0),
        percentile(&rnd_scores, 90.0),
        rnd_scores.last().copied().unwrap_or(f64::NEG_INFINITY));

    let pos_max = pos_scores.last().copied().unwrap_or(f64::NEG_INFINITY);
    let rnd_max = rnd_scores.last().copied().unwrap_or(f64::NEG_INFINITY);
    let sep = pos_max - rnd_max;
    println!("  separation (pos_max − rnd_max): {:.2} nats", sep);

    // AUROC over all (pos, rnd) pairs.
    let mut wins = 0usize; let mut ties = 0usize; let mut tot = 0usize;
    for &p in &pos_scores {
        for &r in &rnd_scores {
            tot += 1;
            if p > r { wins += 1; } else if p == r { ties += 1; }
        }
    }
    let auroc = if tot == 0 { 0.5 } else { (wins as f64 + 0.5 * ties as f64) / tot as f64 };
    println!("  AUROC (positive > random): {:.3}", auroc);
    Ok(())
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let fasta = args.iter().position(|a| a == "--fasta")
        .and_then(|i| args.get(i + 1)).map(|s| s.as_str())
        .unwrap_or("/scratch/jxi21/Assembler/GGO.fasta");
    let bam_path = args.iter().position(|a| a == "--bam")
        .and_then(|i| args.get(i + 1)).map(|s| s.as_str())
        .unwrap_or("/scratch/jxi21/Assembler/GGO_19.bam");

    eprintln!("[loo] loading genome from {} ...", fasta);
    let genome = GenomeIndex::from_fasta(fasta)?;
    eprintln!("[loo] genome loaded; opening BAM {}", bam_path);
    let mut reader = noodles_bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    println!("════════════════════════════════════════════════════════════════════");
    println!("  Leave-one-out HMM rescue experiment — GOLGA6L7");
    println!("════════════════════════════════════════════════════════════════════");

    // Sample reads per paralog region.
    let mut pos_reads: Vec<Vec<Vec<u8>>> = vec![Vec::new(); 3];
    for c in 0..3 {
        let exons = paralog_exons(c);
        let s = exons.iter().map(|e| e.0).min().unwrap() + 1;
        let e = exons.iter().map(|e| e.1).max().unwrap();
        let region = format!("{}:{}-{}", CHROM, s, e);
        let reads: Vec<Vec<u8>> = collect_reads!(reader, header, &region, 200usize, 200usize);
        println!("  reads in {} region ({}-{}): {}", PARALOG_LABELS[c], s, e, reads.len());
        pos_reads[c] = reads;
    }

    // Random reads from a non-family region elsewhere on the same chromosome.
    let rnd_region = format!("{}:1-100000000", CHROM);
    let rnd_reads: Vec<Vec<u8>> = collect_reads!(reader, header, &rnd_region, 200usize, 200usize);
    println!("  random reads (chr {}:1-100M, excl. family region not strict): {}",
        CHROM, rnd_reads.len());

    // For each held-out paralog, run the experiment.
    for c in 0..3 {
        if pos_reads[c].is_empty() {
            println!();
            println!("──── HOLD OUT  {}  (skipped — 0 reads in region) ────", PARALOG_LABELS[c]);
            continue;
        }
        run_loo_for_paralog(c, &genome, &pos_reads[c], &rnd_reads)?;
    }

    println!();
    println!("════════════════════════════════════════════════════════════════════");
    println!("Notes:");
    println!("  - GGO_19.bam contains reads only for chr19 (NC_073243.2), so this");
    println!("    experiment runs on GOLGA6L7 only. Other families (GOLGA8, AMY2A/B,");
    println!("    TBC1D3, OR cluster, NBPF, KRAB-ZNF) require their respective");
    println!("    chromosome BAMs to test.");
    println!("  - 'Positive' reads here are mapped to the held-out paralog's region.");
    println!("    In a real 'novel copy' scenario, those reads would be unmapped");
    println!("    (because the genome lacks that paralog). The HMM scoring is");
    println!("    sequence-only, so the score behavior is identical either way.");
    println!("  - The 'random' set is sampled from the same chromosome and may");
    println!("    include reads from elsewhere in the GOLGA6L7 region (we don't");
    println!("    explicitly exclude). For a strict negative control, sample from");
    println!("    a different chromosome — this BAM only has chr19.");
    Ok(())
}
