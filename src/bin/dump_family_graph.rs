/// Smoke checkpoint: Task 3.2.5 of the vg-hmm novel-copy plan.
///
/// Extends the Task 1.6.5 smoke to:
///   1. Load the genome FASTA and build the family graph WITH sequences.
///   2. Fit per-exon profile HMMs on every node.
///   3. Sample two read populations from GGO_19.bam via indexed query:
///        - Positive:  primary alignments in NC_073243.2:104789647-104877901
///        - Random:    primary alignments in NC_073243.2:1-100000000
///   4. Score each population with forward_against_family and print statistics.
///
/// Delete after Phase 7.

use std::collections::HashMap;
use std::sync::Arc;

use rustle::genome::GenomeIndex;
use rustle::types::{Bundle, BundleRead, Junction, JunctionStat, JunctionStats};
use rustle::vg::FamilyGroup;
use rustle::vg_hmm::family_graph::{build_family_graph, fit_profiles_in_place};
use rustle::vg_hmm::scorer::forward_against_family;

// ── Hard-coded GOLGA6L7 exon coordinates (0-based, half-open) ──────────────
// Extracted from GGO_genomic.gff: $4-1 gives 0-based start, $5 is end.
//
// LOC115931294 / XM_055367547.2 (copy 0)
const COPY0_EXONS: [(u64, u64); 9] = [
    (104789646, 104789918),
    (104791223, 104791352),
    (104792166, 104792196),
    (104792465, 104792516),
    (104792604, 104792698),
    (104792780, 104792887),
    (104794129, 104794188),
    (104794517, 104794659),
    (104794953, 104796276),
];

// LOC134757625 / XM_063700849.1 (copy 1)
const COPY1_EXONS: [(u64, u64); 9] = [
    (104830535, 104830737),
    (104832042, 104832171),
    (104832985, 104833015),
    (104833284, 104833335),
    (104833423, 104833517),
    (104833599, 104833706),
    (104834952, 104835011),
    (104835340, 104835482),
    (104835776, 104837094),
];

// LOC101137218 / XM_055367548.2 (copy 2)
const COPY2_EXONS: [(u64, u64); 9] = [
    (104871355, 104871545),
    (104872850, 104872979),
    (104873793, 104873823),
    (104874092, 104874143),
    (104874231, 104874325),
    (104874407, 104874514),
    (104875758, 104875817),
    (104876146, 104876288),
    (104876582, 104877901),
];

fn build_junction_stats(exons: &[(u64, u64)]) -> JunctionStats {
    let mut stats: JunctionStats = HashMap::with_hasher(fxhash::FxBuildHasher::default());
    // 8 junctions for 9 exons: donor = exon[i].end, acceptor = exon[i+1].start
    for i in 0..exons.len() - 1 {
        let junc = Junction { donor: exons[i].1, acceptor: exons[i + 1].0 };
        stats.insert(junc, JunctionStat::default());
    }
    stats
}

fn make_bundle(exons: &[(u64, u64)]) -> Bundle {
    let start = exons.iter().map(|e| e.0).min().unwrap();
    let end = exons.iter().map(|e| e.1).max().unwrap();

    let read = BundleRead {
        read_uid: 1,
        read_name: Arc::from("synthetic_r1"),
        read_name_hash: 0,
        ref_id: None,
        mate_ref_id: None,
        mate_start: None,
        hi: 0,
        ref_start: start,
        ref_end: end,
        exons: exons.to_vec(),
        junctions: Vec::new(),
        junction_valid: Vec::new(),
        junctions_raw: Vec::new(),
        junctions_del: Vec::new(),
        weight: 1.0,
        is_reverse: false,
        strand: '+',
        has_poly_start: false,
        has_poly_end: false,
        has_poly_start_aligned: false,
        has_poly_start_unaligned: false,
        has_poly_end_aligned: false,
        has_poly_end_unaligned: false,
        unaligned_poly_t: 0,
        unaligned_poly_a: 0,
        has_last_exon_polya: false,
        has_first_exon_polyt: false,
        query_length: None,
        clip_left: 0,
        clip_right: 0,
        nh: 1,
        nm: 0,
        md: None,
        insertion_sites: Vec::new(),
        unitig: false,
        unitig_cov: 0.0,
        read_count_yc: 0.0,
        countfrag_len: 0.0,
        countfrag_num: 0.0,
        junc_mismatch_weight: 0.0,
        pair_idx: Vec::new(),
        pair_count: Vec::new(),
        mapq: 0,
        mismatches: Vec::new(),
        hp_tag: None,
        ps_tag: None,
        is_primary_alignment: true,
    };

    Bundle {
        chrom: "NC_073243.2".to_string(),
        start,
        end,
        strand: '+',
        reads: vec![read],
        junction_stats: build_junction_stats(exons),
        bundlenodes: None,
        read_bnodes: None,
        bnode_colors: None,
    }
}

/// Decode a noodles BAM 4-bit-encoded sequence to ASCII bytes.
fn decode_seq(seq_obj: &noodles_bam::record::Sequence) -> Vec<u8> {
    let seq_raw = seq_obj.as_ref();
    let seq_len = seq_obj.len();
    let mut seq_bytes: Vec<u8> = Vec::with_capacity(seq_len);
    for (i, &byte) in seq_raw.iter().enumerate() {
        let nibbles = [(byte >> 4) & 0xF, byte & 0xF];
        for (j, nib) in nibbles.iter().enumerate() {
            if i * 2 + j >= seq_len {
                break;
            }
            let b = match nib {
                1 => b'A',
                2 => b'C',
                4 => b'G',
                8 => b'T',
                _ => b'N',
            };
            seq_bytes.push(b);
        }
    }
    seq_bytes
}

fn percentile(sorted: &[f64], pct: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NEG_INFINITY;
    }
    let idx = ((sorted.len() as f64 * pct / 100.0).ceil() as usize).saturating_sub(1);
    sorted[idx.min(sorted.len() - 1)]
}

fn main() -> anyhow::Result<()> {
    // ── CLI args ───────────────────────────────────────────────────────────────
    let args: Vec<String> = std::env::args().collect();
    let fasta_path = args.iter().position(|a| a == "--fasta")
        .and_then(|i| args.get(i + 1))
        .map(|s| s.as_str())
        .unwrap_or("/mnt/c/Users/jfris/Desktop/GGO.fasta");
    let bam_path = args.iter().position(|a| a == "--bam")
        .and_then(|i| args.get(i + 1))
        .map(|s| s.as_str())
        .unwrap_or("/mnt/c/Users/jfris/Desktop/GGO_19.bam");

    // ── Build bundles ──────────────────────────────────────────────────────────
    let bundles = vec![
        make_bundle(&COPY0_EXONS),
        make_bundle(&COPY1_EXONS),
        make_bundle(&COPY2_EXONS),
    ];

    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: vec![0, 1, 2],
        multimap_reads: HashMap::new(),
    };

    // ── Load genome FASTA ─────────────────────────────────────────────────────
    eprintln!("[smoke] loading genome from {} ...", fasta_path);
    let genome = GenomeIndex::from_fasta(fasta_path)?;
    eprintln!("[smoke] genome loaded");

    // ── Build family graph WITH sequences ─────────────────────────────────────
    let mut fg = build_family_graph(&family, &bundles, Some(&genome), 0.30, 0.30)?;
    println!(
        "[smoke] family graph: {} nodes, {} edges",
        fg.n_nodes(),
        fg.n_edges()
    );
    for n in &fg.nodes {
        println!(
            "  node {:>2}: {} {}-{} ({} copies, copy_specific={})",
            n.idx.0,
            n.chrom,
            n.span.0,
            n.span.1,
            n.per_copy_sequences.len(),
            n.copy_specific
        );
    }
    for e in &fg.edges {
        println!(
            "  edge: {:>2} -> {:>2} (support={}, strand={})",
            e.from.0, e.to.0, e.family_support, e.strand
        );
    }

    // ── Fit profiles ──────────────────────────────────────────────────────────
    eprintln!("[smoke] fitting profiles ...");
    fit_profiles_in_place(&mut fg)?;
    let n_with_profile = fg.nodes.iter().filter(|n| n.profile.is_some()).count();
    println!("[smoke] profiles fit: {}/{} nodes have a profile", n_with_profile, fg.n_nodes());

    // ── Open indexed BAM ──────────────────────────────────────────────────────
    eprintln!("[smoke] opening indexed BAM: {}", bam_path);
    let mut reader = noodles_bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    const CAP: usize = 50;
    const POS_REGION: &str = "NC_073243.2:104789647-104877901";
    const RND_REGION: &str = "NC_073243.2:1-100000000";

    // ── Collect positive reads ────────────────────────────────────────────────
    eprintln!("[smoke] collecting positive reads from {} ...", POS_REGION);
    let pos_region = POS_REGION.parse()?;
    let mut pos_seqs: Vec<Vec<u8>> = Vec::new();
    for result in reader.query(&header, &pos_region)? {
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let seq = decode_seq(&record.sequence());
        if seq.len() < 50 {
            continue;
        }
        pos_seqs.push(seq);
        if pos_seqs.len() >= CAP {
            break;
        }
    }
    println!("[smoke] positive reads collected: {}", pos_seqs.len());

    // ── Collect random reads ──────────────────────────────────────────────────
    // Query NC_073243.2:1-100000000 but skip reads whose start overlaps the GOLGA6L7 region.
    // We use start position only (no alignment_end on raw Record); start < 104789647 is sufficient.
    eprintln!("[smoke] collecting random reads from {} ...", RND_REGION);
    let rnd_region = RND_REGION.parse()?;
    const POS_START_1B: u64 = 104_789_647; // 1-based start of GOLGA6L7 region
    let mut rnd_seqs: Vec<Vec<u8>> = Vec::new();
    for result in reader.query(&header, &rnd_region)? {
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        // Exclude reads whose alignment start falls within the GOLGA6L7 region.
        // alignment_start() returns Option<io::Result<Position>> (1-based).
        let aln_start_1b: u64 = match record.alignment_start() {
            Some(Ok(p)) => p.get() as u64,
            _ => continue,
        };
        if aln_start_1b >= POS_START_1B {
            continue; // in or past GOLGA6L7 region — exclude
        }
        let seq = decode_seq(&record.sequence());
        if seq.len() < 50 {
            continue;
        }
        rnd_seqs.push(seq);
        if rnd_seqs.len() >= CAP {
            break;
        }
    }
    println!("[smoke] random reads collected: {}", rnd_seqs.len());

    // ── Score both populations via forward_against_family (full-length reads) ──
    // Single-DP boundary threading makes this O(N × L × profile_len) per read —
    // tractable for 3-10 kb long reads across 27 profiles.
    eprintln!("[smoke] scoring positive reads via forward_against_family (full-length) ...");
    let mut pos_scores: Vec<f64> = pos_seqs.iter()
        .map(|s| forward_against_family(&fg, s))
        .collect();
    pos_scores.sort_by(|a, b| a.partial_cmp(b).unwrap());

    eprintln!("[smoke] scoring random reads via forward_against_family (full-length) ...");
    let mut rnd_scores: Vec<f64> = rnd_seqs.iter()
        .map(|s| forward_against_family(&fg, s))
        .collect();
    rnd_scores.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // ── Print statistics ──────────────────────────────────────────────────────
    println!();
    println!("[smoke] ── Score distribution (log-likelihood, nats) ──");
    let pos_median = percentile(&pos_scores, 50.0);
    let pos_90 = percentile(&pos_scores, 90.0);
    let pos_max = pos_scores.last().copied().unwrap_or(f64::NEG_INFINITY);
    println!(
        "  positive  n={:<3}  median={:.2}  90th={:.2}  max={:.2}",
        pos_scores.len(), pos_median, pos_90, pos_max
    );

    let rnd_median = percentile(&rnd_scores, 50.0);
    let rnd_90 = percentile(&rnd_scores, 90.0);
    let rnd_max = rnd_scores.last().copied().unwrap_or(f64::NEG_INFINITY);
    println!(
        "  random    n={:<3}  median={:.2}  90th={:.2}  max={:.2}",
        rnd_scores.len(), rnd_median, rnd_90, rnd_max
    );

    let separation = pos_max - rnd_max;
    println!("  separation (positive_max - random_max): {:.2} nats", separation);

    if separation < 5.0 {
        println!("WARN: separation lower than expected (< 5 nats) — check profile fitting and BAM region");
    }

    Ok(())
}
