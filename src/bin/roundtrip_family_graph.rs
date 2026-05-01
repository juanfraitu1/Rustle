//! Round-trip lossless-recovery experiment for multi-copy gene families.
//!
//! For each family, runs two experiments:
//!
//! Experiment 1 — current rustle FamilyGraph (position-clustering)
//!   Documents the limitation: position-overlap clustering can never merge
//!   paralogs that don't coordinate-overlap, so widely-spaced tandems
//!   produce a disjoint-union graph (trivially lossless but uninteresting).
//!   For overlapping paralogs (AMY2A/AMY2B), this DOES merge.
//!
//! Experiment 2 — POA-based primitive (vg + hmm)
//!   Per exon position, runs partial-order alignment via poasta across all
//!   paralogs' exon sequences. The resulting POA graph is a sequence
//!   variation graph: nodes carry base content, paralog sequences become
//!   paths through it. Round-trip walks each paralog's path through every
//!   exon's POA and compares to genomic ground truth (byte-exact). HMM
//!   layer fitted via ProfileHmm::from_msa.
//!
//! Families tested:
//!   - GOLGA6L7   3 paralogs,  9 exons, chr19 NC_073243.2 +
//!   - GOLGA8     3 paralogs, 19 exons, chr20 NC_073240.2 +
//!   - AMY2A/B    2 paralogs, 11 exons, chr01 NC_073224.2 - (overlapping)

use std::collections::HashMap;
use std::sync::Arc;

use rustle::genome::GenomeIndex;
use rustle::types::{Bundle, BundleRead, Junction, JunctionStat, JunctionStats};
use rustle::vg::FamilyGroup;
use rustle::vg_hmm::family_graph::{build_family_graph, CopyId, FamilyGraph};
use rustle::vg_hmm::profile::{poa_msa, ProfileHmm};

// ── Families ────────────────────────────────────────────────────────────────

struct Family {
    name: &'static str,
    chrom: &'static str,
    strand: char,
    /// Each paralog: (label, exons in genomic ascending order).
    paralogs: Vec<(&'static str, Vec<(u64, u64)>)>,
}

fn families() -> Vec<Family> {
    vec![
        // GOLGA6L7 — chr19 tandem triple, 40 kb apart, 9 exons each, + strand.
        Family {
            name: "GOLGA6L7",
            chrom: "NC_073243.2",
            strand: '+',
            paralogs: vec![
                ("LOC115931294 / XM_055367547.2", vec![
                    (104789646, 104789918), (104791223, 104791352), (104792166, 104792196),
                    (104792465, 104792516), (104792604, 104792698), (104792780, 104792887),
                    (104794129, 104794188), (104794517, 104794659), (104794953, 104796276),
                ]),
                ("LOC134757625 / XM_063700849.1", vec![
                    (104830535, 104830737), (104832042, 104832171), (104832985, 104833015),
                    (104833284, 104833335), (104833423, 104833517), (104833599, 104833706),
                    (104834952, 104835011), (104835340, 104835482), (104835776, 104837094),
                ]),
                ("LOC101137218 / XM_055367548.2", vec![
                    (104871355, 104871545), (104872850, 104872979), (104873793, 104873823),
                    (104874092, 104874143), (104874231, 104874325), (104874407, 104874514),
                    (104875758, 104875817), (104876146, 104876288), (104876582, 104877901),
                ]),
            ],
        },
        // GOLGA8 — chr20 tandem triple, 19 exons each, + strand. golgin A8N-like.
        Family {
            name: "GOLGA8",
            chrom: "NC_073240.2",
            strand: '+',
            paralogs: vec![
                ("LOC129527175 / XM_063698442.1", vec![
                    (31359988, 31360796), (31362458, 31362578), (31363465, 31363525),
                    (31363706, 31363787), (31364640, 31364679), (31364941, 31364989),
                    (31365077, 31365162), (31365257, 31365367), (31366414, 31366501),
                    (31366607, 31366715), (31366910, 31366998), (31367374, 31367631),
                    (31367884, 31367964), (31369216, 31369293), (31369724, 31369816),
                    (31370390, 31370491), (31370573, 31370671), (31370772, 31370928),
                    (31371013, 31374377),
                ]),
                ("LOC134757307 / XM_063698931.1", vec![
                    (33166121, 33166933), (33168595, 33168715), (33169603, 33169663),
                    (33169844, 33169925), (33170778, 33170817), (33171079, 33171127),
                    (33171215, 33171300), (33171395, 33171505), (33172550, 33172637),
                    (33172743, 33172851), (33173046, 33173134), (33173510, 33173767),
                    (33174020, 33174100), (33175350, 33175427), (33175859, 33175951),
                    (33176525, 33176626), (33176708, 33176806), (33176907, 33177063),
                    (33177148, 33177521),
                ]),
                ("LOC129527176 / XM_055363990.2", vec![
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
        // AMY2A / AMY2B — 11 exons each, - strand, OVERLAPPING coordinates.
        // AMY2A spans 136262907-136313018; AMY2B spans 136309082-136334077.
        // The two genes share the splice acceptor at 136312931. Sequences are
        // genomic-strand bytes (reverse-complement of mRNA) — paralog homology
        // is preserved because both are reverse-complemented the same way.
        Family {
            name: "AMY2A_vs_AMY2B",
            chrom: "NC_073224.2",
            strand: '-',
            paralogs: vec![
                ("AMY2A / XM_031007661.3", vec![
                    (136262907, 136263131), (136264460, 136264586), (136264696, 136264815),
                    (136266785, 136266885), (136266979, 136267102), (136267914, 136268048),
                    (136268814, 136269045), (136269490, 136269688), (136270497, 136270644),
                    (136270989, 136271203), (136312931, 136313018),
                ]),
                ("AMY2B / XM_031007502.3", vec![
                    (136309082, 136309306), (136310625, 136310751), (136310862, 136310981),
                    (136312931, 136313031), (136313125, 136313248), (136314058, 136314192),
                    (136314924, 136315155), (136315604, 136315802), (136316608, 136316755),
                    (136317102, 136317316), (136334001, 136334077),
                ]),
            ],
        },
        // TBC1D3 — chr5 hominoid-specific TBC1D3G/3G-like cluster, 3 paralogs,
        // 15 exons each, - strand. All on NC_073228.2; spread ~10 Mb apart.
        Family {
            name: "TBC1D3",
            chrom: "NC_073228.2",
            strand: '-',
            paralogs: vec![
                ("LOC101144080 / XM_063706317.1", vec![
                    (44678473, 44678918), (44678961, 44679397), (44680378, 44680531),
                    (44680901, 44681001), (44681464, 44681530), (44681914, 44682009),
                    (44682722, 44682843), (44683389, 44683438), (44684689, 44684799),
                    (44685129, 44685237), (44685680, 44685761), (44686243, 44686283),
                    (44687569, 44687655), (44687817, 44687890), (44688167, 44689131),
                ]),
                ("LOC101125558 / XM_055388331.2", vec![
                    (54400032, 54400477), (54400520, 54400959), (54401940, 54402093),
                    (54402463, 54402563), (54403026, 54403092), (54403476, 54403571),
                    (54404284, 54404405), (54404930, 54404979), (54406231, 54406341),
                    (54406671, 54406779), (54407222, 54407303), (54407785, 54407825),
                    (54409111, 54409194), (54409356, 54409429), (54410876, 54410996),
                ]),
                ("LOC129533806 / XM_055388408.2", vec![
                    (54930302, 54930747), (54930790, 54931229), (54932210, 54932363),
                    (54932733, 54932833), (54933296, 54933362), (54933746, 54933841),
                    (54934554, 54934675), (54935200, 54935249), (54936501, 54936611),
                    (54936941, 54937049), (54937492, 54937573), (54938055, 54938095),
                    (54939380, 54939463), (54939625, 54939698), (54941138, 54941251),
                ]),
            ],
        },
        // Olfactory receptor cluster on chr1, 5 paralogs, single-exon each,
        // - strand. Tightly clustered within ~200 kb. ORs are intronless
        // 7TM receptors — the largest gene family in mammals (~500 in
        // gorilla here). Per-exon POA collapses to whole-transcript POA.
        Family {
            name: "OR_chr1_cluster",
            chrom: "NC_073224.2",
            strand: '-',
            paralogs: vec![
                ("LOC101140082 / XM_004028742.4", vec![(6032726, 6033683)]),
                ("LOC115931420 / XM_031003600.3", vec![(6097321, 6098296)]),
                ("LOC101141577 / XM_004028745.5", vec![(6109018, 6109972)]),
                ("LOC101138400 / XM_019039372.4", vec![(6189021, 6189984)]),
                ("LOC101136561 / XM_004028735.4", vec![(6222733, 6223672)]),
            ],
        },
        // NBPF — Neuroblastoma BreakPoint Family, hominoid-specific expansion,
        // chr1 - strand. Famous for variable-number tandem DUF1220/Olduvai
        // domain repeats; paralogs vary widely in exon count (3 to 46+).
        // We pick two with the same exon count (22) — NBPF1 and NBPF3.
        Family {
            name: "NBPF",
            chrom: "NC_073224.2",
            strand: '-',
            paralogs: vec![
                ("NBPF1  LOC115934004 / XM_063698357.1", vec![
                    (124449205, 124451073), (124451691, 124451800), (124452516, 124452689),
                    (124453343, 124453395), (124454104, 124454277), (124454875, 124454927),
                    (124455646, 124455819), (124456427, 124456479), (124457521, 124457685),
                    (124461585, 124461637), (124462914, 124463120), (124463584, 124463657),
                    (124464698, 124464913), (124465747, 124465850), (124467622, 124467832),
                    (124469184, 124469396), (124469859, 124469932), (124470983, 124471198),
                    (124472033, 124472136), (124473935, 124473983), (124473985, 124474123),
                    (124474223, 124474484),
                ]),
                ("NBPF3  LOC101131274 / XM_063702015.1", vec![
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
        // KRAB-ZNF — chr19 zinc-finger pair, 6 exons each, + strand.
        // ZNF765 and ZNF813-like, ~5 kb apart on the same chromosome.
        // KRAB-ZNFs are the largest TF family in primates (~400 in human);
        // chr19 is their main expansion locus. Pair tests a KRAB-domain
        // family with shared 6-exon architecture.
        Family {
            name: "KRAB_ZNF",
            chrom: "NC_073244.2",
            strand: '+',
            paralogs: vec![
                ("ZNF765  LOC101133546 / XM_055370414.2", vec![
                    (65303554, 65303625), (65328903, 65328991), (65333026, 65333153),
                    (65338655, 65340075), (65340076, 65340240), (65340242, 65342879),
                ]),
                ("ZNF813-like  LOC109024216 / XM_063701888.1", vec![
                    (65299948, 65300029), (65303498, 65303625), (65308461, 65310310),
                    (65310363, 65311256), (65311321, 65311645), (65311647, 65312451),
                ]),
            ],
        },
    ]
}

// ── Bundle helpers ──────────────────────────────────────────────────────────

fn build_junction_stats(exons: &[(u64, u64)]) -> JunctionStats {
    let mut stats: JunctionStats = HashMap::with_hasher(fxhash::FxBuildHasher::default());
    for i in 0..exons.len().saturating_sub(1) {
        let junc = Junction { donor: exons[i].1, acceptor: exons[i + 1].0 };
        stats.insert(junc, JunctionStat::default());
    }
    stats
}

fn make_bundle(exons: &[(u64, u64)], chrom: &str, strand: char) -> Bundle {
    let start = exons.iter().map(|e| e.0).min().unwrap();
    let end = exons.iter().map(|e| e.1).max().unwrap();
    let read = BundleRead {
        read_uid: 1,
        read_name: Arc::from("synthetic_r1"),
        read_name_hash: 0,
        ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
        ref_start: start, ref_end: end,
        exons: exons.to_vec(),
        junctions: Vec::new(), junction_valid: Vec::new(),
        junctions_raw: Vec::new(), junctions_del: Vec::new(),
        weight: 1.0, is_reverse: strand == '-', strand,
        has_poly_start: false, has_poly_end: false,
        has_poly_start_aligned: false, has_poly_start_unaligned: false,
        has_poly_end_aligned: false, has_poly_end_unaligned: false,
        unaligned_poly_t: 0, unaligned_poly_a: 0,
        has_last_exon_polya: false, has_first_exon_polyt: false,
        query_length: None, clip_left: 0, clip_right: 0,
        nh: 1, nm: 0, md: None,
        insertion_sites: Vec::new(),
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

fn fetch_concat(genome: &GenomeIndex, chrom: &str, exons: &[(u64, u64)]) -> Vec<u8> {
    let mut out = Vec::new();
    for &(s, e) in exons {
        out.extend_from_slice(&genome.fetch_sequence(chrom, s, e).unwrap_or_default());
    }
    out
}

fn introns_from(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    (0..exons.len().saturating_sub(1))
        .map(|i| (exons[i].1, exons[i + 1].0))
        .collect()
}

fn yn(b: bool) -> &'static str { if b { "PASS" } else { "FAIL" } }
fn pct(num: usize, den: usize) -> f64 {
    if den == 0 { 0.0 } else { 100.0 * num as f64 / den as f64 }
}

// ─── Experiment 1 — rustle FamilyGraph ─────────────────────────────────────

fn recover_copy_from_family_graph(
    fg: &FamilyGraph,
    cid: CopyId,
) -> (Vec<(u64, u64)>, Vec<u8>) {
    let mut entries: Vec<((u64, u64), Vec<u8>)> = Vec::new();
    for node in &fg.nodes {
        let span = node.per_copy_spans.iter().find(|(c, _)| *c == cid).map(|(_, s)| *s);
        let seq = node.per_copy_sequences.iter().find(|(c, _)| *c == cid).map(|(_, s)| s.clone());
        if let (Some(sp), Some(sq)) = (span, seq) {
            entries.push((sp, sq));
        }
    }
    entries.sort_by_key(|(span, _)| span.0);
    let exons: Vec<(u64, u64)> = entries.iter().map(|(s, _)| *s).collect();
    let mut concat = Vec::new();
    for (_, s) in &entries { concat.extend_from_slice(s); }
    (exons, concat)
}

fn run_native_experiment(family: &Family, genome: &GenomeIndex) -> bool {
    println!();
    println!("──── EXPERIMENT 1 — rustle FamilyGraph (position-overlap clustering) ────");
    let bundles: Vec<Bundle> = family.paralogs.iter()
        .map(|(_, exons)| make_bundle(exons, family.chrom, family.strand)).collect();
    let family_group = FamilyGroup {
        family_id: 0,
        bundle_indices: (0..bundles.len()).collect(),
        multimap_reads: HashMap::new(),
    };
    let fg = match build_family_graph(&family_group, &bundles, Some(genome), 0.30, 0.30) {
        Ok(g) => g,
        Err(e) => { eprintln!("[exp 1] build_family_graph failed: {}", e); return false; }
    };
    let shared = fg.nodes.iter().filter(|n| !n.copy_specific).count();
    println!(
        "  family graph: {} nodes ({} shared, {} copy-specific), {} edges",
        fg.n_nodes(), shared, fg.n_nodes() - shared, fg.n_edges()
    );
    if shared == 0 {
        println!("  NOTE: 0 shared nodes — paralogs do not coordinate-overlap; merge degenerates");
        println!("        to disjoint union. Round-trip below is trivially lossless.");
    } else {
        println!("  NOTE: {} shared nodes — position-clustering DID merge across paralogs.", shared);
    }
    let mut all_ok = true;
    for (cid, (label, truth_exons)) in family.paralogs.iter().enumerate() {
        let (rec_exons, rec_seq) = recover_copy_from_family_graph(&fg, cid);
        let truth_seq = fetch_concat(genome, family.chrom, truth_exons);
        let chain_ok = introns_from(truth_exons) == introns_from(&rec_exons);
        let coords_ok = truth_exons.as_slice() == rec_exons.as_slice();
        let seq_ok = truth_seq == rec_seq;
        let pass = chain_ok && coords_ok && seq_ok;
        all_ok &= pass;
        println!(
            "  {} : chain={} coords={} sequence={} ({} bp)",
            label, yn(chain_ok), yn(coords_ok), yn(seq_ok), truth_seq.len()
        );
    }
    all_ok
}

// ─── Experiment 2 — POA primitive ──────────────────────────────────────────

fn ungap(row: &[u8]) -> Vec<u8> {
    row.iter().copied().filter(|&b| b != b'-').collect()
}

#[derive(Default, Clone, Copy)]
struct ColumnStats {
    n_columns: usize,
    fully_conserved: usize,   // all paralogs agree on a non-gap base
    majority_agree: usize,    // > n/2 paralogs agree on a non-gap base, but not all
    polymorphic: usize,       // all paralogs non-gap, no majority
    indel: usize,             // ≥ 1 paralog has a gap
}

#[derive(Default, Clone)]
struct FamilyStats {
    columns: ColumnStats,
    /// Sum of per-column Shannon entropy over A/C/G/T/gap (in bits, log base 2).
    /// Mean = sum / n_columns.
    entropy_sum_bits: f64,
    /// Number of columns at each integer entropy bin (0=fully conserved, 1=2 alleles 50/50, 2=4 alleles equal, etc.).
    /// Histogram bins: [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25].
    entropy_hist: [usize; 9],
    /// For each indel column, count how many paralogs are gapped.
    /// indel_breakdown[k] = number of columns with exactly k paralogs gapped.
    indel_breakdown: Vec<usize>,
    /// Pairwise Hamming distance: hamming[i][j] = # columns where paralog i and j differ
    /// (gaps count as their own symbol). Symmetric.
    hamming: Vec<Vec<usize>>,
    /// Number of distinct alleles per column, summed over columns by allele count.
    /// alleles_per_col[k] = # columns with exactly k distinct symbols (1=conserved).
    alleles_per_col: [usize; 6],  // 1..=5 (5 = all four bases + gap)
}

fn shannon_entropy(probs: impl IntoIterator<Item = f64>) -> f64 {
    let mut h = 0.0;
    for p in probs {
        if p > 0.0 { h -= p * p.log2(); }
    }
    h
}

fn entropy_bin(h: f64) -> usize {
    let bin = (h / 0.25).floor() as i64;
    bin.clamp(0, 8) as usize
}

fn family_stats(msa: &[Vec<u8>]) -> FamilyStats {
    let n_par = msa.len();
    let n_col = msa[0].len();
    let mut fs = FamilyStats {
        indel_breakdown: vec![0; n_par + 1],
        hamming: vec![vec![0; n_par]; n_par],
        ..Default::default()
    };
    fs.columns.n_columns = n_col;

    for c in 0..n_col {
        let col: Vec<u8> = msa.iter().map(|r| r[c]).collect();

        // Hamming: count differences pairwise.
        for i in 0..n_par {
            for j in (i + 1)..n_par {
                if col[i] != col[j] {
                    fs.hamming[i][j] += 1;
                    fs.hamming[j][i] += 1;
                }
            }
        }

        // Allele set (incl. gap).
        let mut symbols: Vec<u8> = col.clone();
        symbols.sort_unstable();
        symbols.dedup();
        let n_alleles = symbols.len().min(5);
        fs.alleles_per_col[n_alleles] += 1;

        // Indel handling.
        let n_gaps = col.iter().filter(|&&b| b == b'-').count();
        if n_gaps > 0 {
            fs.columns.indel += 1;
            fs.indel_breakdown[n_gaps] += 1;
        } else {
            // Categorize by majority.
            let mut counts: HashMap<u8, usize> = HashMap::new();
            for &b in &col { *counts.entry(b).or_insert(0) += 1; }
            let max_count = counts.values().copied().max().unwrap_or(0);
            if max_count == n_par {
                fs.columns.fully_conserved += 1;
            } else if max_count > n_par / 2 {
                fs.columns.majority_agree += 1;
            } else {
                fs.columns.polymorphic += 1;
            }
        }

        // Shannon entropy over the column (5-symbol alphabet incl. gap).
        let mut counts: HashMap<u8, usize> = HashMap::new();
        for &b in &col { *counts.entry(b).or_insert(0) += 1; }
        let total = n_par as f64;
        let h = shannon_entropy(counts.values().map(|&c| c as f64 / total));
        fs.entropy_sum_bits += h;
        fs.entropy_hist[entropy_bin(h)] += 1;
    }
    fs
}

fn run_poa_experiment(family: &Family, genome: &GenomeIndex) -> bool {
    println!();
    println!("──── EXPERIMENT 2 — POA-based primitive (vg + hmm) ────");
    let n_exons = family.paralogs[0].1.len();
    if !family.paralogs.iter().all(|(_, e)| e.len() == n_exons) {
        println!("  paralogs have varying exon counts — per-exon POA not applicable; skipping.");
        return true;
    }
    let n_par = family.paralogs.len();
    let truths: Vec<Vec<u8>> = family.paralogs.iter()
        .map(|(_, exons)| fetch_concat(genome, family.chrom, exons))
        .collect();

    let mut totals = FamilyStats {
        indel_breakdown: vec![0; n_par + 1],
        hamming: vec![vec![0; n_par]; n_par],
        ..Default::default()
    };
    let mut total_match_cols = 0usize;
    let mut total_insert_cols = 0usize;
    let mut recovered: Vec<Vec<u8>> = vec![Vec::new(); n_par];

    println!("  Per-exon POA conservation (n_paralogs = {}):", n_par);
    println!("  ┌─────┬─────────────┬──────────┬────────────┬──────────┬──────────┬─────────┐");
    println!("  │ ex  │ lens (bp)   │ cols_msa │ fully_cons │ majority │ polymor. │ indel   │");
    println!("  ├─────┼─────────────┼──────────┼────────────┼──────────┼──────────┼─────────┤");

    for i in 0..n_exons {
        let seqs: Vec<Vec<u8>> = family.paralogs.iter()
            .map(|(_, exons)| genome.fetch_sequence(family.chrom, exons[i].0, exons[i].1).unwrap_or_default())
            .collect();
        let lens_str = seqs.iter().map(|s| s.len().to_string())
            .collect::<Vec<_>>().join("/");

        let msa = match poa_msa(&seqs) {
            Ok(m) => m,
            Err(e) => { eprintln!("[exp 2] {} exon {} poa_msa failed: {}", family.name, i, e); return false; }
        };
        let fs = family_stats(&msa);
        // Accumulate.
        totals.columns.n_columns += fs.columns.n_columns;
        totals.columns.fully_conserved += fs.columns.fully_conserved;
        totals.columns.majority_agree += fs.columns.majority_agree;
        totals.columns.polymorphic += fs.columns.polymorphic;
        totals.columns.indel += fs.columns.indel;
        totals.entropy_sum_bits += fs.entropy_sum_bits;
        for k in 0..9 { totals.entropy_hist[k] += fs.entropy_hist[k]; }
        for k in 0..(n_par + 1) { totals.indel_breakdown[k] += fs.indel_breakdown[k]; }
        for k in 0..6 { totals.alleles_per_col[k] += fs.alleles_per_col[k]; }
        for r in 0..n_par {
            for c in 0..n_par { totals.hamming[r][c] += fs.hamming[r][c]; }
        }

        for c in 0..n_par {
            let r = ungap(&msa[c]);
            if r != seqs[c] {
                println!("    ! exon {} paralog {}: ungap(POA row) != input ({} vs {} bp)",
                    i, c, r.len(), seqs[c].len());
                return false;
            }
            recovered[c].extend_from_slice(&r);
        }

        let hmm = ProfileHmm::from_msa(&msa).unwrap_or_else(|_| ProfileHmm::empty(0));
        total_match_cols += hmm.n_columns;
        total_insert_cols += fs.columns.n_columns - hmm.n_columns;

        println!(
            "  │ {:>3} │ {:>11} │ {:>8} │ {:>10} │ {:>8} │ {:>8} │ {:>7} │",
            i, lens_str, fs.columns.n_columns, fs.columns.fully_conserved,
            fs.columns.majority_agree, fs.columns.polymorphic, fs.columns.indel
        );
    }
    println!("  └─────┴─────────────┴──────────┴────────────┴──────────┴──────────┴─────────┘");

    println!();
    println!("  Round-trip (paralog == ungap(POA path)):");
    let mut all_ok = true;
    for (cid, (label, _)) in family.paralogs.iter().enumerate() {
        let ok = recovered[cid] == truths[cid];
        all_ok &= ok;
        println!(
            "    {} : {} ({} bp recovered, {} bp truth)",
            label, yn(ok), recovered[cid].len(), truths[cid].len()
        );
    }

    println!();
    println!("  Family-wide MSA conservation:");
    let n_total = totals.columns.n_columns;
    println!(
        "    fully conserved        : {:>5} / {:>5} ({:>5.1}%)",
        totals.columns.fully_conserved, n_total, pct(totals.columns.fully_conserved, n_total)
    );
    println!(
        "    majority agree         : {:>5} / {:>5} ({:>5.1}%)",
        totals.columns.majority_agree, n_total, pct(totals.columns.majority_agree, n_total)
    );
    println!(
        "    polymorphic            : {:>5} / {:>5} ({:>5.1}%)",
        totals.columns.polymorphic, n_total, pct(totals.columns.polymorphic, n_total)
    );
    println!(
        "    indel column           : {:>5} / {:>5} ({:>5.1}%)",
        totals.columns.indel, n_total, pct(totals.columns.indel, n_total)
    );

    println!();
    println!("  Per-column Shannon entropy (5-symbol alphabet incl. gap, log2 bits):");
    let max_h = (n_par as f64).log2().min(5f64.log2());
    let mean_h = if n_total == 0 { 0.0 } else { totals.entropy_sum_bits / n_total as f64 };
    println!("    mean entropy           : {:.4} bits  (max possible ≈ {:.2} bits)", mean_h, max_h);
    println!("    distribution (cols per bin):");
    let bin_labels = ["[0.00,0.25)", "[0.25,0.50)", "[0.50,0.75)", "[0.75,1.00)",
                       "[1.00,1.25)", "[1.25,1.50)", "[1.50,1.75)", "[1.75,2.00)", "[2.00,∞)"];
    for k in 0..9 {
        if totals.entropy_hist[k] > 0 {
            println!("      H ∈ {}  : {:>6} ({:>5.1}%)",
                bin_labels[k], totals.entropy_hist[k], pct(totals.entropy_hist[k], n_total));
        }
    }

    println!();
    println!("  Distinct alleles per column (incl. gap as 5th symbol):");
    for k in 1..=5 {
        if totals.alleles_per_col[k] > 0 {
            println!("    {} alleles : {:>6} ({:>5.1}%)",
                k, totals.alleles_per_col[k], pct(totals.alleles_per_col[k], n_total));
        }
    }

    if totals.columns.indel > 0 {
        println!();
        println!("  Indel breakdown (cols by # gapped paralogs out of {}):", n_par);
        for k in 1..(n_par + 1) {
            if totals.indel_breakdown[k] > 0 {
                println!("    {} gapped : {:>6} ({:>5.1}%  of {} indel cols)",
                    k, totals.indel_breakdown[k],
                    pct(totals.indel_breakdown[k], totals.columns.indel),
                    totals.columns.indel);
            }
        }
    }

    if n_par >= 2 {
        println!();
        println!("  Pairwise Hamming distance (cols where paralogs differ; gap = 5th symbol):");
        // Header.
        print!("    {:>3} │", "");
        for j in 0..n_par { print!(" {:>7} ", format!("p{}", j)); }
        println!();
        print!("    ────┼");
        for _ in 0..n_par { print!("─────────"); }
        println!();
        for i in 0..n_par {
            print!("    p{:<2} │", i);
            for j in 0..n_par {
                if i == j { print!(" {:>7} ", "—"); }
                else { print!(" {:>7} ", totals.hamming[i][j]); }
            }
            println!();
        }
        // Mean pairwise % identity.
        let mut tot_pairs = 0usize; let mut tot_diffs = 0usize;
        for i in 0..n_par {
            for j in (i + 1)..n_par {
                tot_pairs += 1; tot_diffs += totals.hamming[i][j];
            }
        }
        if tot_pairs > 0 {
            let mean_diffs = tot_diffs as f64 / tot_pairs as f64;
            let mean_id = if n_total == 0 { 0.0 } else { 100.0 * (1.0 - mean_diffs / n_total as f64) };
            println!("    mean pairwise identity : {:.2}%  ({:.0} avg differing cols / pair)",
                mean_id, mean_diffs);
        }
    }

    println!();
    println!("  HMM layer:");
    println!("    match columns (M)  : {}", total_match_cols);
    println!("    insert columns (I) : {}", total_insert_cols);

    all_ok
}

// ─── Main ──────────────────────────────────────────────────────────────────

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let fasta_path = args.iter().position(|a| a == "--fasta")
        .and_then(|i| args.get(i + 1))
        .map(|s| s.as_str())
        .unwrap_or("/scratch/jxi21/Assembler/GGO.fasta");
    let only_family: Option<&str> = args.iter().position(|a| a == "--family")
        .and_then(|i| args.get(i + 1)).map(|s| s.as_str());

    eprintln!("[round-trip] loading genome from {} ...", fasta_path);
    let genome = GenomeIndex::from_fasta(fasta_path)?;
    eprintln!("[round-trip] genome loaded");

    let fams = families();
    let mut summary: Vec<(String, bool, bool)> = Vec::new();
    for fam in &fams {
        if let Some(only) = only_family {
            if fam.name != only { continue; }
        }
        println!();
        println!("════════════════════════════════════════════════════════════════════");
        println!("  FAMILY — {}    chrom {}    strand {}    paralogs={}",
            fam.name, fam.chrom, fam.strand, fam.paralogs.len());
        println!("════════════════════════════════════════════════════════════════════");
        let ok1 = run_native_experiment(fam, &genome);
        let ok2 = run_poa_experiment(fam, &genome);
        summary.push((fam.name.to_string(), ok1, ok2));
    }

    println!();
    println!("════════════════════════════════════════════════════════════════════");
    println!("  SUMMARY");
    println!("════════════════════════════════════════════════════════════════════");
    println!("  {:<22} {:<14} {:<14}", "family", "exp1 (FG)", "exp2 (POA)");
    let mut all_pass = true;
    for (name, ok1, ok2) in &summary {
        println!("  {:<22} {:<14} {:<14}", name, yn(*ok1), yn(*ok2));
        all_pass &= *ok1 && *ok2;
    }
    if !all_pass {
        std::process::exit(1);
    }
    Ok(())
}
