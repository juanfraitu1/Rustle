//! End-to-end leave-one-out HMM rescue + assembly across multiple families.
//!
//! For each family and each held-out paralog c:
//!   1. Build family primitive from the OTHER paralogs.
//!   2. Score reads from the held-out paralog's region.
//!   3. Compute Viterbi path (with per-node read spans) for each rescued read.
//!   4. Cluster by path and synthesize one Bundle per cluster.
//!   5. Compare three reconstructions to truth:
//!        [A] RAW           — family-graph node coords (production floor)
//!        [B] VG-ONLY       — Viterbi-derived per-exon length consensus
//!                            (no CIGAR or external aligner needed)
//!        [C] CIGAR-REFINED — StringTie-style consensus over read alignments
//!                            (works only when reads have CIGAR; upper bound)
//!
//! Use `--bam GGO_19.bam` (default) for chr19 families or `--bam GGO.bam`
//! for the full genome BAM.
//! Use `--family <name>` to run a single family; default runs all.

use std::collections::HashMap;
use std::sync::Arc;

use rustle::genome::GenomeIndex;
use rustle::types::{Bundle, BundleRead, Junction, JunctionStat, JunctionStats};
use rustle::vg::FamilyGroup;
use rustle::vg_hmm::family_graph::{build_family_graph, fit_profiles_in_place};
use rustle::vg_hmm::rescue::synthesize_bundles;
use rustle::vg_hmm::scorer::{forward_against_family, viterbi_path};
use rustle::tss_tts::cluster_positions_with_counts;

// ── Families ────────────────────────────────────────────────────────────────

struct Family {
    name: &'static str,
    chrom: &'static str,
    strand: char,
    paralogs: Vec<(&'static str, Vec<(u64, u64)>)>,
}

fn families() -> Vec<Family> {
    vec![
        Family {
            name: "GOLGA6L7",
            chrom: "NC_073243.2",
            strand: '+',
            paralogs: vec![
                ("COPY0/LOC115931294", vec![
                    (104789646, 104789918), (104791223, 104791352), (104792166, 104792196),
                    (104792465, 104792516), (104792604, 104792698), (104792780, 104792887),
                    (104794129, 104794188), (104794517, 104794659), (104794953, 104796276),
                ]),
                ("COPY1/LOC134757625", vec![
                    (104830535, 104830737), (104832042, 104832171), (104832985, 104833015),
                    (104833284, 104833335), (104833423, 104833517), (104833599, 104833706),
                    (104834952, 104835011), (104835340, 104835482), (104835776, 104837094),
                ]),
                ("COPY2/LOC101137218", vec![
                    (104871355, 104871545), (104872850, 104872979), (104873793, 104873823),
                    (104874092, 104874143), (104874231, 104874325), (104874407, 104874514),
                    (104875758, 104875817), (104876146, 104876288), (104876582, 104877901),
                ]),
            ],
        },
        Family {
            name: "GOLGA8",
            chrom: "NC_073240.2",
            strand: '+',
            paralogs: vec![
                ("LOC129527175", vec![
                    (31359988, 31360796), (31362458, 31362578), (31363465, 31363525),
                    (31363706, 31363787), (31364640, 31364679), (31364941, 31364989),
                    (31365077, 31365162), (31365257, 31365367), (31366414, 31366501),
                    (31366607, 31366715), (31366910, 31366998), (31367374, 31367631),
                    (31367884, 31367964), (31369216, 31369293), (31369724, 31369816),
                    (31370390, 31370491), (31370573, 31370671), (31370772, 31370928),
                    (31371013, 31374377),
                ]),
                ("LOC134757307", vec![
                    (33166121, 33166933), (33168595, 33168715), (33169603, 33169663),
                    (33169844, 33169925), (33170778, 33170817), (33171079, 33171127),
                    (33171215, 33171300), (33171395, 33171505), (33172550, 33172637),
                    (33172743, 33172851), (33173046, 33173134), (33173510, 33173767),
                    (33174020, 33174100), (33175350, 33175427), (33175859, 33175951),
                    (33176525, 33176626), (33176708, 33176806), (33176907, 33177063),
                    (33177148, 33177521),
                ]),
                ("LOC129527176", vec![
                    (33336735, 33337343), (33339005, 33339125), (33340013, 33340073),
                    (33340254, 33340335), (33341188, 33341227), (33341489, 33341537),
                    (33341625, 33341710), (33341805, 33341915), (33342960, 33343047),
                    (33343153, 33343261), (33343456, 33343544), (33343920, 33344177),
                    (33344431, 33344511), (33345762, 33345839), (33346271, 33346363),
                    (33346937, 33347038), (33347120, 33347218), (33347319, 33347475),
                    (33347560, 33347903),
                ]),
            ],
        },
        Family {
            name: "AMY2A_vs_AMY2B",
            chrom: "NC_073224.2",
            strand: '-',
            paralogs: vec![
                ("AMY2A", vec![
                    (136262907, 136263131), (136264460, 136264586), (136264696, 136264815),
                    (136266785, 136266885), (136266979, 136267102), (136267914, 136268048),
                    (136268814, 136269045), (136269490, 136269688), (136270497, 136270644),
                    (136270989, 136271203), (136312931, 136313018),
                ]),
                ("AMY2B", vec![
                    (136309082, 136309306), (136310625, 136310751), (136310862, 136310981),
                    (136312931, 136313031), (136313125, 136313248), (136314058, 136314192),
                    (136314924, 136315155), (136315604, 136315802), (136316608, 136316755),
                    (136317102, 136317316), (136334001, 136334077),
                ]),
            ],
        },
        Family {
            name: "TBC1D3",
            chrom: "NC_073228.2",
            strand: '-',
            paralogs: vec![
                ("LOC101144080", vec![
                    (44678473, 44678918), (44678961, 44679397), (44680378, 44680531),
                    (44680901, 44681001), (44681464, 44681530), (44681914, 44682009),
                    (44682722, 44682843), (44683389, 44683438), (44684689, 44684799),
                    (44685129, 44685237), (44685680, 44685761), (44686243, 44686283),
                    (44687569, 44687655), (44687817, 44687890), (44688167, 44689131),
                ]),
                ("LOC101125558", vec![
                    (54400032, 54400477), (54400520, 54400959), (54401940, 54402093),
                    (54402463, 54402563), (54403026, 54403092), (54403476, 54403571),
                    (54404284, 54404405), (54404930, 54404979), (54406231, 54406341),
                    (54406671, 54406779), (54407222, 54407303), (54407785, 54407825),
                    (54409111, 54409194), (54409356, 54409429), (54410876, 54410996),
                ]),
                ("LOC129533806", vec![
                    (54930302, 54930747), (54930790, 54931229), (54932210, 54932363),
                    (54932733, 54932833), (54933296, 54933362), (54933746, 54933841),
                    (54934554, 54934675), (54935200, 54935249), (54936501, 54936611),
                    (54936941, 54937049), (54937492, 54937573), (54938055, 54938095),
                    (54939380, 54939463), (54939625, 54939698), (54941138, 54941251),
                ]),
            ],
        },
        Family {
            name: "OR_chr1_cluster",
            chrom: "NC_073224.2",
            strand: '-',
            paralogs: vec![
                ("LOC101140082", vec![(6032726, 6033683)]),
                ("LOC115931420", vec![(6097321, 6098296)]),
                ("LOC101141577", vec![(6109018, 6109972)]),
                ("LOC101138400", vec![(6189021, 6189984)]),
                ("LOC101136561", vec![(6222733, 6223672)]),
            ],
        },
        Family {
            name: "NBPF",
            chrom: "NC_073224.2",
            strand: '-',
            paralogs: vec![
                ("NBPF1", vec![
                    (124449205, 124451073), (124451691, 124451800), (124452516, 124452689),
                    (124453343, 124453395), (124454104, 124454277), (124454875, 124454927),
                    (124455646, 124455819), (124456427, 124456479), (124457521, 124457685),
                    (124461585, 124461637), (124462914, 124463120), (124463584, 124463657),
                    (124464698, 124464913), (124465747, 124465850), (124467622, 124467832),
                    (124469184, 124469396), (124469859, 124469932), (124470983, 124471198),
                    (124472033, 124472136), (124473935, 124473983), (124473985, 124474123),
                    (124474223, 124474484),
                ]),
                ("NBPF3", vec![
                    (218361365, 218363123), (218363740, 218363792), (218364506, 218364679),
                    (218365284, 218365336), (218366041, 218366214), (218366846, 218366898),
                    (218367940, 218368122), (218368660, 218368738), (218371280, 218371332),
                    (218372626, 218372832), (218373296, 218373369), (218374419, 218374634),
                    (218375470, 218375573), (218377308, 218377518), (218378566, 218378778),
                    (218379241, 218379314), (218380364, 218380579), (218381415, 218381518),
                    (218383253, 218383463), (218385989, 218386059), (218386328, 218386411),
                    (218405598, 218405843),
                ]),
            ],
        },
        Family {
            name: "KRAB_ZNF",
            chrom: "NC_073244.2",
            strand: '+',
            paralogs: vec![
                ("ZNF765", vec![
                    (65303554, 65303625), (65328903, 65328991), (65333026, 65333153),
                    (65338655, 65340075), (65340076, 65340240), (65340242, 65342879),
                ]),
                ("ZNF813-like", vec![
                    (65299948, 65300029), (65303498, 65303625), (65308461, 65310310),
                    (65310363, 65311256), (65311321, 65311645), (65311647, 65312451),
                ]),
            ],
        },
    ]
}

// ── Bundle helpers ──────────────────────────────────────────────────────────

fn build_junction_stats(exons: &[(u64, u64)]) -> JunctionStats {
    let mut s: JunctionStats = HashMap::with_hasher(fxhash::FxBuildHasher::default());
    for i in 0..exons.len().saturating_sub(1) {
        s.insert(Junction { donor: exons[i].1, acceptor: exons[i + 1].0 }, JunctionStat::default());
    }
    s
}

fn make_bundle(exons: &[(u64, u64)], chrom: &str, strand: char) -> Bundle {
    let start = exons.iter().map(|e| e.0).min().unwrap();
    let end = exons.iter().map(|e| e.1).max().unwrap();
    let read = BundleRead {
        read_uid: 1, read_name: Arc::from("synthetic"), read_name_hash: 0,
        ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
        ref_start: start, ref_end: end, exons: exons.to_vec(),
        junctions: Vec::new(), junction_valid: Vec::new(),
        junctions_raw: Vec::new(), junctions_del: Vec::new(),
        weight: 1.0, is_reverse: strand == '-', strand,
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
        chrom: chrom.to_string(), start, end, strand,
        reads: vec![read],
        junction_stats: build_junction_stats(exons),
        bundlenodes: None, read_bnodes: None, bnode_colors: None,
        synthetic: false, rescue_class: None,
    }
}

// ── BAM-derived reads ──────────────────────────────────────────────────────

#[derive(Clone, Debug)]
struct ReadRec {
    name: String,
    seq: Vec<u8>,
    aln_exons: Vec<(u64, u64)>,
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

fn exons_from_record(rec: &noodles_bam::Record) -> anyhow::Result<Vec<(u64, u64)>> {
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record::Cigar as _;
    let aln_start_1b: u64 = match rec.alignment_start() {
        Some(Ok(p)) => p.get() as u64,
        _ => return Ok(Vec::new()),
    };
    let mut pos: u64 = aln_start_1b - 1;
    let mut block_start: u64 = pos;
    let mut exons: Vec<(u64, u64)> = Vec::new();
    let cigar = rec.cigar();
    for op_res in cigar.iter() {
        let op = op_res?;
        let len = op.len() as u64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => pos += len,
            Kind::Skip => {
                if pos > block_start { exons.push((block_start, pos)); }
                pos += len;
                block_start = pos;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    if pos > block_start { exons.push((block_start, pos)); }
    Ok(exons)
}

macro_rules! collect_reads_named {
    ($reader:expr, $header:expr, $region_str:expr, $cap:expr, $min_len:expr) => {{
        let region = $region_str.parse()?;
        let q = $reader.query(&$header, &region)?;
        let mut out: Vec<ReadRec> = Vec::new();
        for r in q {
            let rec = r?;
            let f = rec.flags();
            if !f.is_unmapped() && !f.is_secondary() && !f.is_supplementary() {
                let seq = decode_seq(&rec.sequence());
                if seq.len() >= $min_len {
                    let name = rec.name().map(|n| String::from_utf8_lossy(n.as_ref()).to_string())
                        .unwrap_or_else(|| format!("read{}", out.len()));
                    let aln_exons = exons_from_record(&rec).unwrap_or_default();
                    out.push(ReadRec { name, seq, aln_exons });
                    if out.len() >= $cap { break; }
                }
            }
        }
        out
    }}
}

// ── Helpers ─────────────────────────────────────────────────────────────────

fn intron_chain_lengths(exons: &[(u64, u64)]) -> (Vec<u64>, Vec<u64>) {
    let mut e_lens: Vec<u64> = exons.iter().map(|&(s, e)| e.saturating_sub(s)).collect();
    let mut i_lens: Vec<u64> = Vec::with_capacity(exons.len().saturating_sub(1));
    for i in 0..exons.len().saturating_sub(1) {
        i_lens.push(exons[i + 1].0.saturating_sub(exons[i].1));
    }
    (std::mem::take(&mut e_lens), std::mem::take(&mut i_lens))
}

fn chain_match(truth: &[(u64, u64)], bundle: &[(u64, u64)], tol: u64) -> (bool, bool, usize, usize) {
    let count_match = truth.len() == bundle.len();
    let (te, ti) = intron_chain_lengths(truth);
    let (be, bi) = intron_chain_lengths(bundle);
    let n_e = te.len().min(be.len());
    let n_i = ti.len().min(bi.len());
    let n_match_e = (0..n_e).filter(|&i| (te[i] as i64 - be[i] as i64).abs() as u64 <= tol).count();
    let n_match_i = (0..n_i).filter(|&i| (ti[i] as i64 - bi[i] as i64).abs() as u64 <= tol).count();
    (count_match, count_match && n_match_e == n_e && n_match_i == n_i, n_match_e, n_match_i)
}

fn yn(b: bool) -> &'static str { if b { "PASS" } else { "FAIL" } }

fn median(xs: &[u64]) -> u64 {
    if xs.is_empty() { return 0; }
    let mut v = xs.to_vec(); v.sort_unstable(); v[v.len() / 2]
}

fn refine_bundle_boundaries(bundle_exons: &[(u64, u64)], reads: &[&ReadRec]) -> Vec<(u64, u64)> {
    let n = bundle_exons.len();
    if n == 0 { return Vec::new(); }
    let full: Vec<&&ReadRec> = reads.iter().filter(|r| r.aln_exons.len() == n).collect();
    if full.is_empty() { return bundle_exons.to_vec(); }
    let mut out: Vec<(u64, u64)> = Vec::with_capacity(n);
    for i in 0..n {
        let starts: Vec<u64> = full.iter().map(|r| r.aln_exons[i].0).collect();
        let ends: Vec<u64> = full.iter().map(|r| r.aln_exons[i].1).collect();
        out.push((median(&starts), median(&ends)));
    }
    out
}

// ── Per-paralog routine ────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct Outcome { raw: bool, vg: bool, hybrid: bool, clustered: bool, cigar: bool, n_reads: usize, n_rescued: usize }

fn run_for_paralog(
    family: &Family,
    held_out: usize,
    pos_reads: &[ReadRec],
    genome: &GenomeIndex,
    score_threshold: f64,
    min_reads_per_bundle: usize,
    hard_min_reads: usize,
    cluster_window: u64,
) -> anyhow::Result<Outcome> {
    println!();
    println!("──── HOLD OUT  {}/{}  ────", family.name, family.paralogs[held_out].0);

    let other_ids: Vec<usize> = (0..family.paralogs.len()).filter(|&i| i != held_out).collect();
    let bundles: Vec<Bundle> = other_ids.iter()
        .map(|&i| make_bundle(&family.paralogs[i].1, family.chrom, family.strand))
        .collect();
    let group = FamilyGroup {
        family_id: 0,
        bundle_indices: (0..bundles.len()).collect(),
        multimap_reads: HashMap::new(),
    };
    let mut fg = build_family_graph(&group, &bundles, Some(genome), 0.30, 0.30)?;
    fit_profiles_in_place(&mut fg)?;
    println!("  primitive: {} nodes, {} edges; {} reads in region", fg.n_nodes(), fg.n_edges(), pos_reads.len());

    if pos_reads.is_empty() {
        return Ok(Outcome { raw: false, vg: false, hybrid: false, clustered: false, cigar: false, n_reads: 0, n_rescued: 0 });
    }

    // Length-normalized rescue: score per read base. Forward log-prob scales
    // linearly with read length, so a fixed nats threshold doesn't generalize
    // across families (GOLGA6L7 ≈ 2.2 kb transcripts vs GOLGA8 ≈ 5.8 kb).
    // We use score / read_length (nats per base). Empirical threshold:
    // matches typically yield about -1.0 to -2.0 nats/base; junk yields more
    // negative. The CLI --threshold is interpreted as nats/base.
    let mut scored: Vec<(usize, String, f64, Vec<rustle::vg_hmm::family_graph::NodeIdx>, Vec<(usize, usize)>)> = Vec::new();
    for (idx, r) in pos_reads.iter().enumerate() {
        let s_total = forward_against_family(&fg, &r.seq);
        let s_per_base = if r.seq.is_empty() { f64::NEG_INFINITY } else { s_total / r.seq.len() as f64 };
        if s_per_base <= score_threshold { continue; }
        if let Some(vp) = viterbi_path(&fg, &r.seq) {
            if !vp.nodes.is_empty() {
                scored.push((idx, r.name.clone(), s_total, vp.nodes, vp.read_spans));
            }
        }
    }
    println!("  rescued reads (score/base > {:.2}): {}", score_threshold, scored.len());

    if scored.is_empty() {
        return Ok(Outcome { raw: false, vg: false, hybrid: false, clustered: false, cigar: false, n_reads: pos_reads.len(), n_rescued: 0 });
    }

    let rescued_with_paths: Vec<(String, Vec<rustle::vg_hmm::family_graph::NodeIdx>)> =
        scored.iter().map(|(_, n, _, p, _)| (n.clone(), p.clone())).collect();
    let synth = synthesize_bundles(&fg, &rescued_with_paths, min_reads_per_bundle, 15, false);
    println!("  synthetic bundles (≥{} reads): {}", min_reads_per_bundle, synth.len());
    if synth.is_empty() {
        return Ok(Outcome { raw: false, vg: false, hybrid: false, clustered: false, cigar: false, n_reads: pos_reads.len(), n_rescued: scored.len() });
    }

    let mut cluster_of: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
    let mut spans_by_path: HashMap<Vec<usize>, Vec<Vec<u64>>> = HashMap::new();
    for (orig_idx, _, _, path, spans) in &scored {
        let key: Vec<usize> = path.iter().map(|n| n.0).collect();
        cluster_of.entry(key.clone()).or_default().push(*orig_idx);
        let entry = spans_by_path.entry(key).or_insert_with(|| vec![Vec::new(); path.len()]);
        let n_nodes = spans.len();
        let read_len = pos_reads[*orig_idx].seq.len() as u64;
        for (i, &(s, e)) in spans.iter().enumerate() {
            if i >= entry.len() { continue; }
            // First exon footprint = everything before the first profile's exit
            // (includes pre-profile boundary-state bases — the read's 5' UTR
            // tip past the splice donor).
            // Last exon footprint = everything after the last profile's entry
            // (includes post-profile boundary-state bases — TES extension).
            // Internal exons = within-profile span (exit - entry).
            let footprint = if i == 0 && n_nodes > 1 {
                e as u64
            } else if i + 1 == n_nodes && n_nodes > 1 {
                read_len.saturating_sub(s as u64)
            } else {
                e.saturating_sub(s) as u64
            };
            entry[i].push(footprint);
        }
    }

    let truth = &family.paralogs[held_out].1;

    let mut any_raw = false; let mut any_vg = false; let mut any_hybrid = false;
    let mut any_clustered = false; let mut any_cigar = false;
    for (bi, b) in synth.iter().enumerate() {
        let raw_exons = b.reads[0].exons.clone();
        let n_reads = b.reads.len();
        let mut contributors: Vec<&ReadRec> = Vec::new();
        for (_k, members) in &cluster_of {
            if members.len() == n_reads {
                contributors = members.iter().map(|&i| &pos_reads[i]).collect();
                break;
            }
        }

        let path_key: Vec<usize> = (0..raw_exons.len())
            .filter_map(|i| fg.nodes.iter().position(|n| n.span == raw_exons[i])).collect();

        // Helper: compute consensus exon list given an aggregator for boundary
        // exons (first / last). Internal exons keep the family-graph node span
        // (splice junctions conserved). Closure picks length for boundary exon.
        let build_consensus = |boundary_agg: &dyn Fn(&[u64], u64) -> u64| -> Vec<(u64, u64)> {
            if path_key.len() != raw_exons.len() { return raw_exons.clone(); }
            let spans = match spans_by_path.get(&path_key) {
                Some(s) => s, None => return raw_exons.clone(),
            };
            let n_ex = raw_exons.len();
            let mut out = Vec::with_capacity(n_ex);
            for (i, &(node_s, node_e)) in raw_exons.iter().enumerate() {
                let node_len = node_e - node_s;
                let len = if i == 0 || i + 1 == n_ex {
                    let footprints: &[u64] = if i < spans.len() { &spans[i] } else { &[] };
                    if footprints.is_empty() { node_len } else { boundary_agg(footprints, node_len) }
                } else {
                    // Internal: family-graph node span is the splice-junction-bounded length.
                    node_len
                };
                if i == 0 && n_ex > 1 {
                    out.push((node_e.saturating_sub(len), node_e));
                } else if i + 1 == n_ex && n_ex > 1 {
                    out.push((node_s, node_s + len));
                } else {
                    out.push((node_s, node_e));
                }
            }
            out
        };

        // [B] VG-only: median consensus on read footprints, applied to both boundary
        // exons. Symmetric.
        let vg_exons = build_consensus(&|f, _node| median(f));
        // [D] HYBRID: strand-aware asymmetric anchoring.
        // The TSS exon is where reads have their 5' end (often 5'-truncated by
        // PacBio IsoSeq), so we use max(node_len, max-of-read-footprints) to
        // recover truncation. The TES exon is where reads have their 3' end
        // (poly-A capture; max can over-extend into genomic flanking), so we
        // use median for robustness. Internal exons keep family-graph node
        // span (both junctions conserved).
        //   + strand: TSS = genomic-first exon, TES = genomic-last exon.
        //   - strand: TSS = genomic-last exon (after rev-comp alignment),
        //             TES = genomic-first exon.
        let n_ex_local = raw_exons.len();
        let strand_plus = family.strand == '+';
        let tss_is_first = strand_plus;
        let hybrid_exons = if path_key.len() == raw_exons.len() {
            if let Some(spans) = spans_by_path.get(&path_key) {
                let mut out = Vec::with_capacity(n_ex_local);
                for (i, &(node_s, node_e)) in raw_exons.iter().enumerate() {
                    let node_len = node_e - node_s;
                    let is_first = i == 0;
                    let is_last = i + 1 == n_ex_local;
                    let is_tss = (is_first && tss_is_first) || (is_last && !tss_is_first);
                    let is_tes = (is_last && tss_is_first) || (is_first && !tss_is_first);
                    let footprints: &[u64] = if i < spans.len() { &spans[i] } else { &[] };
                    let len = if (is_tss || is_tes) && n_ex_local > 1 {
                        if footprints.is_empty() { node_len }
                        else if is_tss {
                            (*footprints.iter().max().unwrap_or(&node_len)).max(node_len)
                        } else {
                            median(footprints)
                        }
                    } else {
                        node_len
                    };
                    if is_first && n_ex_local > 1 {
                        out.push((node_e.saturating_sub(len), node_e));
                    } else if is_last && n_ex_local > 1 {
                        out.push((node_s, node_s + len));
                    } else {
                        out.push((node_s, node_e));
                    }
                }
                out
            } else { raw_exons.clone() }
        } else { raw_exons.clone() };

        let cigar_exons = refine_bundle_boundaries(&raw_exons, &contributors);

        // [E] CLUSTERED: StringTie-style cluster_positions_with_counts on the
        // per-read footprint distribution. Pick the cluster with the highest
        // read count as the consensus boundary (filters singleton outliers like
        // AMY2A's 1057-bp 5'UTR-extended read). Hard-boundary fallback: if no
        // cluster has ≥hard_min_reads support, fall back to RAW (family-graph
        // node length) — i.e., refuse to invent a TSS/TES from too few reads.
        let clustered_exons: Vec<(u64, u64)> = if path_key.len() == raw_exons.len() {
            if let Some(spans) = spans_by_path.get(&path_key) {
                let mut out = Vec::with_capacity(n_ex_local);
                for (i, &(node_s, node_e)) in raw_exons.iter().enumerate() {
                    let node_len = node_e - node_s;
                    let is_first = i == 0;
                    let is_last = i + 1 == n_ex_local;
                    let is_boundary = (is_first || is_last) && n_ex_local > 1;
                    let footprints: &[u64] = if i < spans.len() { &spans[i] } else { &[] };
                    let len = if !is_boundary {
                        // Internal: keep family-graph span (junctions conserved).
                        node_len
                    } else if footprints.is_empty() {
                        node_len
                    } else {
                        // Cluster footprints (each read = weight 1) within window.
                        let pos: Vec<(u64, f64)> = footprints.iter().map(|&f| (f, 1.0)).collect();
                        let clusters = cluster_positions_with_counts(&pos, cluster_window, 1);
                        // Pick cluster with highest count.
                        let best = clusters.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
                        match best {
                            Some(&(center, w)) if (w as usize) >= hard_min_reads => center,
                            // Hard-boundary fallback: too few reads in any cluster.
                            _ => node_len,
                        }
                    };
                    if is_first && n_ex_local > 1 {
                        out.push((node_e.saturating_sub(len), node_e));
                    } else if is_last && n_ex_local > 1 {
                        out.push((node_s, node_s + len));
                    } else {
                        out.push((node_s, node_e));
                    }
                }
                out
            } else { raw_exons.clone() }
        } else { raw_exons.clone() };

        let (raw_e, raw_i) = intron_chain_lengths(&raw_exons);
        let (vg_e, vg_i)   = intron_chain_lengths(&vg_exons);
        let (hyb_e, hyb_i) = intron_chain_lengths(&hybrid_exons);
        let (clu_e, clu_i) = intron_chain_lengths(&clustered_exons);
        let (cig_e, cig_i) = intron_chain_lengths(&cigar_exons);
        let (_, raw_ok, raw_em, raw_im) = chain_match(truth, &raw_exons,       10);
        let (_, vg_ok, vg_em, vg_im)    = chain_match(truth, &vg_exons,        10);
        let (_, hyb_ok, hyb_em, hyb_im) = chain_match(truth, &hybrid_exons,    10);
        let (_, clu_ok, clu_em, clu_im) = chain_match(truth, &clustered_exons, 10);
        let (_, cig_ok, cig_em, cig_im) = chain_match(truth, &cigar_exons,     10);
        any_raw |= raw_ok; any_vg |= vg_ok; any_hybrid |= hyb_ok;
        any_clustered |= clu_ok; any_cigar |= cig_ok;

        let truth_e = intron_chain_lengths(truth).0;
        let truth_i = intron_chain_lengths(truth).1;
        println!("  bundle #{} ({} reads, {} exons; {} truth exons):", bi, n_reads, raw_exons.len(), truth.len());
        println!("    truth exon lens : {:?}", truth_e);
        println!("    truth intron    : {:?}", truth_i);
        println!("    [A] RAW        e={:?} i={:?}  ok=({}/{} e, {}/{} i, chain={})",
            raw_e, raw_i, raw_em, raw_e.len().min(truth_e.len()),
            raw_im, raw_i.len().min(truth_i.len()), yn(raw_ok));
        println!("    [B] VG-ONLY    e={:?} i={:?}  ok=({}/{} e, {}/{} i, chain={})",
            vg_e, vg_i, vg_em, vg_e.len().min(truth_e.len()),
            vg_im, vg_i.len().min(truth_i.len()), yn(vg_ok));
        println!("    [D] HYBRID     e={:?} i={:?}  ok=({}/{} e, {}/{} i, chain={})",
            hyb_e, hyb_i, hyb_em, hyb_e.len().min(truth_e.len()),
            hyb_im, hyb_i.len().min(truth_i.len()), yn(hyb_ok));
        println!("    [E] CLUSTERED  e={:?} i={:?}  ok=({}/{} e, {}/{} i, chain={})",
            clu_e, clu_i, clu_em, clu_e.len().min(truth_e.len()),
            clu_im, clu_i.len().min(truth_i.len()), yn(clu_ok));
        println!("    [C] CIGAR-REF  e={:?} i={:?}  ok=({}/{} e, {}/{} i, chain={})",
            cig_e, cig_i, cig_em, cig_e.len().min(truth_e.len()),
            cig_im, cig_i.len().min(truth_i.len()), yn(cig_ok));
    }
    Ok(Outcome { raw: any_raw, vg: any_vg, hybrid: any_hybrid, clustered: any_clustered,
        cigar: any_cigar, n_reads: pos_reads.len(), n_rescued: scored.len() })
}

// ── Main ───────────────────────────────────────────────────────────────────

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let fasta = args.iter().position(|a| a == "--fasta")
        .and_then(|i| args.get(i + 1)).map(|s| s.as_str())
        .unwrap_or("/scratch/jxi21/Assembler/GGO.fasta");
    let bam_path = args.iter().position(|a| a == "--bam")
        .and_then(|i| args.get(i + 1)).map(|s| s.as_str())
        .unwrap_or("/scratch/jxi21/Assembler/GGO.bam");
    let only_family: Option<&str> = args.iter().position(|a| a == "--family")
        .and_then(|i| args.get(i + 1)).map(|s| s.as_str());
    // Length-normalized: nats per read base. Empirical: family-matching reads
    // typically give -1.5 to -2.5 nats/base; junk reads -3.0 to -10.0 nats/base.
    // -3.0 is permissive enough for low-identity rescue (~17% identity in KRAB-ZNF).
    let score_threshold: f64 = args.iter().position(|a| a == "--threshold")
        .and_then(|i| args.get(i + 1)).and_then(|s| s.parse().ok())
        .unwrap_or(-3.0);
    let min_reads: usize = args.iter().position(|a| a == "--min-reads")
        .and_then(|i| args.get(i + 1)).and_then(|s| s.parse().ok())
        .unwrap_or(2);
    // Hard-boundary requirement (StringTie-style): a TSS/TES cluster must have
    // ≥hard_min_reads supporting reads to be reported. Below this, fall back
    // to the family-graph node length (RAW). Default 3 — works for low-coverage
    // (2-5 reads/paralog) families by leaning on the family-graph prior.
    let hard_min_reads: usize = args.iter().position(|a| a == "--hard-min-reads")
        .and_then(|i| args.get(i + 1)).and_then(|s| s.parse().ok())
        .unwrap_or(3);
    // Position-cluster window for footprint clustering (bp). 20 bp matches
    // StringTie's long-read TSS jitter tolerance.
    let cluster_window: u64 = args.iter().position(|a| a == "--cluster-window")
        .and_then(|i| args.get(i + 1)).and_then(|s| s.parse().ok())
        .unwrap_or(20);

    eprintln!("[loo] loading genome from {} ...", fasta);
    let genome = GenomeIndex::from_fasta(fasta)?;
    eprintln!("[loo] opening BAM {}", bam_path);
    let mut reader = noodles_bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)?;
    let header = reader.read_header()?;

    let fams = families();
    let mut summary: Vec<(String, String, Outcome)> = Vec::new();
    for fam in &fams {
        if let Some(only) = only_family { if fam.name != only { continue; } }
        println!();
        println!("════════════════════════════════════════════════════════════════════");
        println!("  FAMILY {}    chrom {}    strand {}    paralogs={}", fam.name, fam.chrom, fam.strand, fam.paralogs.len());
        println!("════════════════════════════════════════════════════════════════════");

        let mut pos_reads: Vec<Vec<ReadRec>> = vec![Vec::new(); fam.paralogs.len()];
        for (c, (_label, exons)) in fam.paralogs.iter().enumerate() {
            let s = exons.iter().map(|e| e.0).min().unwrap() + 1;
            let e = exons.iter().map(|e| e.1).max().unwrap();
            let region = format!("{}:{}-{}", fam.chrom, s, e);
            pos_reads[c] = collect_reads_named!(reader, header, &region, 30usize, 200usize);
            println!("  reads in {} region: {}", fam.paralogs[c].0, pos_reads[c].len());
        }

        for c in 0..fam.paralogs.len() {
            let outcome = run_for_paralog(fam, c, &pos_reads[c], &genome, score_threshold, min_reads, hard_min_reads, cluster_window)?;
            summary.push((fam.name.to_string(), fam.paralogs[c].0.to_string(), outcome));
        }
    }

    println!();
    println!("════════════════════════════════════════════════════════════════════");
    println!("  SUMMARY");
    println!("════════════════════════════════════════════════════════════════════");
    println!("  {:<22} {:<24} {:>6} {:>9} {:>6} {:>6} {:>6} {:>9} {:>6}",
        "family", "held-out", "reads", "rescued", "RAW", "VG", "HYBRID", "CLUSTERED", "CIGAR");
    for (f, p, o) in &summary {
        println!("  {:<22} {:<24} {:>6} {:>9} {:>6} {:>6} {:>6} {:>9} {:>6}",
            f, p, o.n_reads, o.n_rescued,
            yn(o.raw), yn(o.vg), yn(o.hybrid), yn(o.clustered), yn(o.cigar));
    }
    Ok(())
}
