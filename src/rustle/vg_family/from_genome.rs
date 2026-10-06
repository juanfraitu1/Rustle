//! DNA front-end: discover duplicated genomic loci by self-alignment and emit them as reps for the
//! shared homology-grouping core (`denovo_pipeline::homology_blocks`). Read-free and annotation-free —
//! the genome-only counterpart of the RNA read front-end. The reps differ from RNA reps in exactly one
//! way: genomic `seq` (introns included) and an EMPTY intron chain. That single difference is the
//! scientific claim (splicing discards the intron/flank sequence that separates near-identical copies).
//!
//! **STATUS:** OPT-IN — `--from-genome <BED>` (src/bin/gw_family_catalog.rs:38-39, `#[arg(long)] from_genome: Option<String>`, default None)  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)
use std::collections::{HashMap, HashSet};
use anyhow::Result;
use crate::vg_family::family_detect::DenovoTranscript;
use crate::genome::GenomeIndex;
use crate::vg_family::genome_projection::project_families_batch;

/// Parameters for the DNA-mode SD detector. `min_identity` is the LOCUS-DISCOVERY floor ("is this segment
/// duplicated enough to be a candidate locus"); it is distinct from the family-GROUPING identity floor
/// applied later inside `homology_blocks`. Env overrides let the SD floor be swept (e.g. 0.98 for SD98).
pub struct GenomeRepParams {
    pub min_identity: f64,
    pub min_block: u64,
    pub max_locus_span: u64,
    pub minimap2: String,
    pub threads: usize,
}
impl Default for GenomeRepParams {
    fn default() -> Self {
        Self { min_identity: 0.90, min_block: 1000, max_locus_span: 3_000_000,
                minimap2: "minimap2".into(), threads: 4 }
    }
}
impl GenomeRepParams {
    pub fn from_env() -> Self {
        let mut p = Self::default();
        if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_IDENTITY") { if let Ok(x) = v.parse() { p.min_identity = x; } }
        if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_BLOCK") { if let Ok(x) = v.parse() { p.min_block = x; } }
        if let Ok(v) = std::env::var("RUSTLE_MINIMAP2") { p.minimap2 = v; }
        p
    }
}

/// Discover duplicated genomic loci within `windows` and return them as reps (genomic sequence, empty
/// intron chain) for `homology_blocks`. Steps: (1) self-align each window's sequence against the genome
/// FASTA (`project_families_batch`, the same minimap2 primitive the famCN projection uses) to find every
/// locus it recurs at; (2) collect the hit loci across all windows, keep blocks in `[min_block,
/// max_locus_span]`, merge overlapping loci; (3) emit one rep per merged locus with its genomic sequence.
pub fn genome_reps(
    fasta_path: &str,
    windows: &[(String, u64, u64)],
    p: &GenomeRepParams,
) -> Result<Vec<DenovoTranscript>> {
    let contigs: HashSet<String> = windows.iter().map(|(c, _, _)| c.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;

    // (1) SD detector: each window's sequence is a query; project_families_batch returns every genome
    // locus it recurs at (identity >= min_identity). One batched minimap2 pass. `known` empty = keep all
    // hits (incl. the self locus — a window is itself a candidate locus; grouping decides families).
    let consensuses: Vec<(String, Vec<u8>)> = windows.iter().enumerate().filter_map(|(i, (c, s, e))| {
        genome.fetch_sequence(c, *s, *e).map(|seq| (format!("w{i}"), seq))
    }).collect();
    let known: HashMap<String, Vec<(String, u64, u64)>> = HashMap::new();
    let cov = 0.0_f64; // block length is gated by min_block below, not by fractional window coverage
    let hits = project_families_batch(&consensuses, fasta_path, &known, p.min_identity, cov,
                                      &p.minimap2, p.threads)?;

    // (2) rep construction. Each WINDOW is a locus to group — emit it directly so every member locus is a
    // node the quasi-clique can place (this is Soto's setup: given the loci, group them by shared sequence).
    // The self-alignment hits then ADD paralog loci OUTSIDE the windows, giving windowed singletons a sibling.
    // Windows are NOT merged (they are distinct member loci) and NOT min_block-filtered (a short member is
    // still a locus); only the extra discovered loci are size-gated and deduped, so dense paralog regions
    // can't fuse into giant blocks that miss a member's exact coordinate.
    let mut reps = Vec::new();
    let mut window_spans: Vec<(String, u64, u64)> = Vec::with_capacity(windows.len());
    for (c, s, e) in windows {
        if let Some(seq) = genome.fetch_sequence(c, *s, *e) {
            reps.push(DenovoTranscript {
                tid: format!("DN_{c}_{s}_1"),
                chrom: c.clone(), start: *s, end: *e, n_reads: 1, strand: '+',
                introns: vec![], seq, distinguishing_uniq: 0,
                core_bp: 0,
                stub: false, tes: None,
            });
            window_spans.push((c.clone(), *s, *e));
        }
    }
    // discovered paralog loci that fall OUTSIDE every window (a paralog at a locus no member covers).
    let mut extra: Vec<(String, u64, u64)> = Vec::new();
    for hs in hits.into_values() {
        for h in hs {
            let len = h.end.saturating_sub(h.start);
            if len < p.min_block || len > p.max_locus_span { continue; }
            let inside_window = window_spans.iter().any(|(wc, ws, we)| *wc == h.chrom && h.start < *we && *ws < h.end);
            if inside_window { continue; }
            extra.push((h.chrom, h.start, h.end));
        }
    }
    extra.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    for (chrom, start, end) in merge_overlapping(&extra) {
        if let Some(seq) = genome.fetch_sequence(&chrom, start, end) {
            reps.push(DenovoTranscript {
                tid: format!("DN_{chrom}_{start}_1"),
                chrom, start, end, n_reads: 1, strand: '+',
                introns: vec![], seq, distinguishing_uniq: 0,
                core_bp: 0,
                stub: false, tes: None,
            });
        }
    }
    Ok(reps)
}

/// Derive `--from-genome` search windows from an existing SD-pairs BED (same 6-column format
/// `annotation_families::SdPairs::from_bed_str` reads: chrom1,start1,end1,chrom2,start2,end2, extra
/// columns ignored), instead of requiring the caller to hand-supply a windows BED (proposal #2,
/// docs/o1_ledger.md §6j5). Every interval on EITHER side of every SD pair becomes a candidate window
/// (genome-wide self-similarity has already flagged it as duplicated -- no new arbitrary threshold, no
/// tile-size/overlap constant to invent); overlapping intervals on the same contig are merged so a dense
/// SD region doesn't emit many redundant, near-identical windows.
pub fn windows_from_sd_bed(path: &str) -> Result<Vec<(String, u64, u64)>> {
    let mut raw: Vec<(String, u64, u64)> = Vec::new();
    for line in std::fs::read_to_string(path)?.lines() {
        if line.is_empty() || line.starts_with('#') { continue; }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 { continue; }
        let (Ok(s1), Ok(e1), Ok(s2), Ok(e2)) =
            (f[1].parse::<u64>(), f[2].parse::<u64>(), f[4].parse::<u64>(), f[5].parse::<u64>())
        else { continue };
        if e1 > s1 { raw.push((f[0].to_string(), s1, e1)); }
        if e2 > s2 { raw.push((f[3].to_string(), s2, e2)); }
    }
    raw.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    Ok(merge_overlapping(&raw))
}

/// Single-linkage merge of overlapping genomic intervals (input sorted by (chrom, start)).
fn merge_overlapping(loci: &[(String, u64, u64)]) -> Vec<(String, u64, u64)> {
    let mut out: Vec<(String, u64, u64)> = Vec::new();
    for (c, s, e) in loci.iter().cloned() {
        match out.last_mut() {
            Some((pc, _ps, pe)) if *pc == c && s <= *pe => { *pe = (*pe).max(e); }
            _ => out.push((c, s, e)),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn windows_from_sd_bed_extracts_both_sides_and_merges_overlaps() {
        let dir = std::env::temp_dir().join(format!("rustle_sdbed_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let bed = dir.join("sd.bed");
        std::fs::write(
            &bed,
            "chr1\t100\t200\tchr2\t500\t600\t95.0\t+\t+\n\
             chr1\t150\t250\tchr3\t10\t20\t93.0\t+\t-\n\
             # a comment line, ignored\n\
             chr9\t0\t5\tchr9\t0\t5\t100.0\t+\t+\n", // degenerate zero-length-safe pair, both sides valid
        ).unwrap();
        let windows = windows_from_sd_bed(bed.to_str().unwrap()).unwrap();
        std::fs::remove_dir_all(&dir).ok();
        // chr1: [100,200) and [150,250) overlap -> merged to [100,250).
        assert!(windows.contains(&("chr1".to_string(), 100, 250)), "{windows:?}");
        assert!(windows.contains(&("chr2".to_string(), 500, 600)), "{windows:?}");
        assert!(windows.contains(&("chr3".to_string(), 10, 20)), "{windows:?}");
        assert!(windows.contains(&("chr9".to_string(), 0, 5)), "{windows:?}");
        assert_eq!(windows.len(), 4, "chr1's two overlapping sides must merge into one window: {windows:?}");
    }

    #[test]
    fn merge_overlapping_joins_adjacent_and_keeps_disjoint() {
        let loci = vec![
            ("chr1".to_string(), 10, 100),
            ("chr1".to_string(), 90, 200),   // overlaps previous -> merge to 10..200
            ("chr1".to_string(), 500, 600),  // disjoint -> separate
            ("chr2".to_string(), 0, 50),     // different contig -> separate
        ];
        let m = merge_overlapping(&loci);
        assert_eq!(m, vec![
            ("chr1".to_string(), 10, 200),
            ("chr1".to_string(), 500, 600),
            ("chr2".to_string(), 0, 50),
        ]);
    }

    /// Proposal #2 Stage A ceiling test (docs/o1_ledger.md §6j5): if we already know where all 31 true
    /// gorilla NPIP loci are (the oracle, not the annotation -- this uses ONLY genomic coordinates and
    /// self-alignment, no gene model), does raw DNA self-alignment alone (the exact `project_families_batch`
    /// primitive `genome_reps` calls) recover the true NPIP family structure -- i.e. does each locus's own
    /// self-alignment hit set include at least one OTHER oracle locus? This bounds the achievable ceiling for
    /// proposal #2 before any window-generation-without-prior-knowledge design work (Stage B) is attempted:
    /// if DNA self-alignment can't even connect these loci when told exactly where they are, no amount of
    /// clever window generation will make the downstream signal appear.
    /// `#[ignore]`d: needs the real ~3.6GB gorilla genome FASTA on disk, not a checked-in fixture.
    #[test]
    #[ignore]
    fn stage_a_ceiling_dna_self_alignment_recovers_npip_structure_from_true_coords() {
        use crate::vg_family::genome_projection::project_families_batch;
        let fa = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta";
        if std::fs::metadata(fa).is_err() { eprintln!("gorilla genome fasta absent; skip"); return; }
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }

        // The 31-locus NPIP oracle (o1_oracle/npip31.regions), "chrom:1based_start-1based_end" per line,
        // converted to 0-based half-open.
        let oracle_txt = std::fs::read_to_string("/mnt/linuxdisk/home/juanfraitu/o1_oracle/npip31.regions")
            .expect("oracle regions file");
        let mut windows: Vec<(String, u64, u64)> = Vec::new();
        for line in oracle_txt.lines() {
            let line = line.trim();
            if line.is_empty() { continue; }
            let (chrom, rest) = line.split_once(':').expect("chrom:start-end");
            let (s, e) = rest.split_once('-').expect("start-end");
            let s1: u64 = s.parse().unwrap();
            let e1: u64 = e.parse().unwrap();
            windows.push((chrom.to_string(), s1 - 1, e1)); // 1-based inclusive -> 0-based half-open
        }
        assert_eq!(windows.len(), 31, "expected all 31 oracle loci");

        let contigs: HashSet<String> = windows.iter().map(|(c, _, _)| c.clone()).collect();
        let genome = GenomeIndex::from_fasta_contigs(fa, &contigs).expect("genome index");
        let consensuses: Vec<(String, Vec<u8>)> = windows.iter().enumerate().filter_map(|(i, (c, s, e))| {
            genome.fetch_sequence(c, *s, *e).map(|seq| (format!("w{i}"), seq))
        }).collect();
        assert_eq!(consensuses.len(), 31, "every oracle window must fetch real sequence");

        let known: HashMap<String, Vec<(String, u64, u64)>> = HashMap::new();
        let p = GenomeRepParams::default(); // min_identity 0.90, the same floor genome_reps() uses
        let hits = project_families_batch(&consensuses, fa, &known, p.min_identity, 0.0, &p.minimap2, p.threads)
            .expect("project_families_batch");

        let mut connected_to_another_locus = 0usize;
        let mut connected_absent_only = 0usize; // connects, and the query window is one of the 26 "absent" ones
        // (absent = has zero de novo RNA node today; computed once, offline, against ggo.nodes.tsv --
        // hardcoded here since this is a one-shot ceiling measurement, not a maintained pipeline path)
        let absent_idx: HashSet<usize> = [0,1,2,3,4,6,7,8,9,10,11,12,13,15,18,19,20,21,22,23,24,26,27,28,29,30]
            .into_iter().collect();
        for (i, _) in windows.iter().enumerate() {
            let qid = format!("w{i}");
            let Some(hs) = hits.get(&qid) else { continue };
            let hit_other_locus = hs.iter().any(|h| {
                windows.iter().enumerate().any(|(j, (oc, os, oe))| {
                    j != i && &h.chrom == oc && h.end > *os && h.start < *oe
                })
            });
            if hit_other_locus {
                connected_to_another_locus += 1;
                if absent_idx.contains(&i) { connected_absent_only += 1; }
            }
        }
        eprintln!(
            "[stage_a] {}/31 oracle loci have a self-alignment hit landing on >=1 OTHER oracle locus \
             ({}/{} of the 26 currently-RNA-absent ones)",
            connected_to_another_locus, connected_absent_only, absent_idx.len()
        );
    }

    /// Proposal #2 Stage B test (docs/o1_ledger.md §6j5): unlike Stage A (which used the oracle's OWN true
    /// coordinates as windows -- a ceiling test, not a real mechanism), this feeds `genome_reps` windows
    /// derived ONLY from real SD-pair evidence (`windows_from_sd_bed`-style extraction, restricted here to
    /// SD pairs touching the 31-locus NPIP oracle region, to keep the batch small and safe -- a genome-wide
    /// run of the same mechanism triggered a real, disclosed memory blow-up, see the ledger) -- no oracle
    /// coordinates are used as INPUT, only as the scoring truth afterward. Tests whether SD-derived (not
    /// hand-known) windows recover the same NPIP structure Stage A showed is theoretically reachable.
    /// `#[ignore]`d: needs the real gorilla genome FASTA + SD file on disk.
    #[test]
    #[ignore]
    fn stage_b_sd_derived_windows_recover_npip_structure() {
        let fa = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta";
        if std::fs::metadata(fa).is_err() { eprintln!("gorilla genome fasta absent; skip"); return; }
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
        let win_bed = "/mnt/linuxdisk/home/juanfraitu/o1_fromgenome_sd/npip_seeded_windows2.bed";
        if std::fs::metadata(win_bed).is_err() { eprintln!("seeded windows file absent; skip"); return; }

        let mut windows: Vec<(String, u64, u64)> = Vec::new();
        for line in std::fs::read_to_string(win_bed).unwrap().lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() >= 3 { windows.push((f[0].to_string(), f[1].parse().unwrap(), f[2].parse().unwrap())); }
        }
        eprintln!("[stage_b] {} SD-derived windows (no oracle coordinates used as input)", windows.len());

        let p = GenomeRepParams::default();
        let reps = genome_reps(fa, &windows, &p).expect("genome_reps");
        eprintln!("[stage_b] {} reps produced", reps.len());

        // Score AFTER the fact against the 31-locus oracle -- truth used only for scoring, never as input.
        let oracle_txt = std::fs::read_to_string("/mnt/linuxdisk/home/juanfraitu/o1_oracle/npip31.regions")
            .expect("oracle regions file");
        let mut oracle: Vec<(String, u64, u64)> = Vec::new();
        for line in oracle_txt.lines() {
            let line = line.trim();
            if line.is_empty() { continue; }
            let (chrom, rest) = line.split_once(':').unwrap();
            let (s, e) = rest.split_once('-').unwrap();
            let s1: u64 = s.parse().unwrap();
            let e1: u64 = e.parse().unwrap();
            oracle.push((chrom.to_string(), s1 - 1, e1));
        }
        let absent_idx: HashSet<usize> = [0,1,2,3,4,6,7,8,9,10,11,12,13,15,18,19,20,21,22,23,24,26,27,28,29,30]
            .into_iter().collect();
        let mut covered = 0usize;
        let mut covered_absent = 0usize;
        for (i, (oc, os_, oe)) in oracle.iter().enumerate() {
            let hit = reps.iter().any(|r| &r.chrom == oc && r.end > *os_ && r.start < *oe);
            if hit {
                covered += 1;
                if absent_idx.contains(&i) { covered_absent += 1; }
            }
        }
        eprintln!(
            "[stage_b] {}/31 oracle loci covered by an SD-derived genome_reps rep ({}/{} of the 26 \
             currently-RNA-absent ones)",
            covered, covered_absent, absent_idx.len()
        );
    }

    /// Proposal #1 + #2 COMBINED (docs/o1_ledger.md §6j6): wires both shipped-but-unconnected
    /// primitives -- `SdPairs::single_span_core` (§6j4, corrects an existing RNA rep's span) and
    /// `genome_reps`/`windows_from_sd_bed` (§6j5, adds DNA-only reps for RNA-absent loci) -- into ONE
    /// combined rep set, then runs it through the REAL, UNCHANGED downstream pipeline
    /// (`family_detect::detect_edges`, `T_CORE=0.13` untouched, then `family_split::decompose_families`,
    /// the shipped gamma-quasi-clique partition) to measure whether this actually changes real family
    /// recovery on the 31-locus NPIP oracle, AND at what real precision cost (does any oracle-containing
    /// family also swallow a non-oracle gene from the same real 392-copy corpus). Baseline reps are the
    /// REAL, already-computed footprint-off run (§6j3's control) on the real 3-contig NPIP-bearing BAM
    /// subset -- not synthetic, not re-derived. `#[ignore]`d: needs real gorilla data on disk.
    ///
    /// Uses the default edge core (LCS since docs/o1_ledger.md §6jd; under the old POA default this hung, §6j6).
    #[test]
    #[ignore]
    fn stage_c_combining_proposals_1_and_2_measures_real_family_recovery() {
        stage_c_body("default edge core");
    }

    fn stage_c_body(label_suffix: &str) {
        use crate::vg_family::annotation_families::SdPairs;
        use crate::vg_family::family_detect::{detect_edges, DetectParams};
        use crate::vg_family::family_detect::family_split::{decompose_families, SplitParams};

        let fa = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta";
        if std::fs::metadata(fa).is_err() { eprintln!("gorilla genome fasta absent; skip"); return; }
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }

        let copies_tsv = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.tsv";
        let copies_fa = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.fa";
        let sedef_bed = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed";
        let sd_windows_bed = "/mnt/linuxdisk/home/juanfraitu/o1_fromgenome_sd/npip_seeded_windows2.bed";
        for p in [copies_tsv, copies_fa, sedef_bed, sd_windows_bed] {
            if std::fs::metadata(p).is_err() { eprintln!("required real-data file {p} absent; skip"); return; }
        }

        // --- load the real, already-computed baseline RNA reps (footprint OFF, the shipped default) ---
        let tsv_text = std::fs::read_to_string(copies_tsv).unwrap();
        let fa_text = std::fs::read_to_string(copies_fa).unwrap();
        // FASTA records are in the SAME row order as the TSV (both written by the same run, one pass).
        let seqs: Vec<Vec<u8>> = fa_text.lines().filter(|l| !l.starts_with('>')).map(|l| l.as_bytes().to_vec()).collect();
        let mut baseline: Vec<DenovoTranscript> = Vec::new();
        for (i, line) in tsv_text.lines().skip(1).enumerate() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 11 { continue; }
            let chrom = f[3].to_string();
            let start: u64 = f[4].parse().unwrap();
            let end: u64 = f[5].parse().unwrap();
            let strand = f[7].chars().next().unwrap_or('+');
            let n_reads: u32 = f[8].parse().unwrap_or(1);
            let exons: Vec<(u64, u64)> = f[9].split(',').filter_map(|e| {
                let (s, en) = e.split_once('-')?;
                Some((s.parse().ok()?, en.parse().ok()?))
            }).collect();
            let introns: Vec<(u64, u64)> = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
            baseline.push(DenovoTranscript {
                tid: f[2].to_string(), chrom, start, end, n_reads, strand,
                introns, seq: seqs.get(i).cloned().unwrap_or_default(), distinguishing_uniq: 0,
                core_bp: 0, stub: false, tes: None,
            });
        }
        eprintln!("[stage_c] {} baseline RNA reps loaded from the real footprint-OFF run", baseline.len());

        // Bound the compute: the full 391-rep corpus's real pairwise POA edge-confirmation is UNTRACTABLE
        // (a live run on the unfiltered set was killed after ~5+ min CPU-heavy with no result -- an
        // independent, real reconfirmation of proposal #3's own finding, §6j0, that `confirm_edge` has
        // severe, uncapped per-pair cost variance on real production-scale sequences). Restrict to reps
        // overlapping the SAME SD-seeded window region proposal #2 already validated as safe (62 windows,
        // 5.44Mb) -- this is not a shortcut around the precision question: it is the single densest real
        // multi-copy neighborhood on this substrate (all 31 oracle loci plus their real genomic neighbors),
        // so it is if anything a HARDER, more relevant false-merge test than the full corpus, not a weaker one.
        let mut sd_windows: Vec<(String, u64, u64)> = Vec::new();
        for line in std::fs::read_to_string(sd_windows_bed).unwrap().lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() >= 3 { sd_windows.push((f[0].to_string(), f[1].parse().unwrap(), f[2].parse().unwrap())); }
        }
        let n_before_bound = baseline.len();
        baseline.retain(|r| sd_windows.iter().any(|(c, s, e)| &r.chrom == c && r.end > *s && r.start < *e));
        eprintln!(
            "[stage_c] bounded baseline to the SD-seeded-window region for tractable real edge-confirmation: \
             {n_before_bound} -> {} reps",
            baseline.len()
        );

        // --- oracle ---
        let oracle_txt = std::fs::read_to_string("/mnt/linuxdisk/home/juanfraitu/o1_oracle/npip31.regions").unwrap();
        let oracle: Vec<(String, u64, u64)> = oracle_txt.lines().filter(|l| !l.trim().is_empty()).map(|l| {
            let (c, rest) = l.trim().split_once(':').unwrap();
            let (s, e) = rest.split_once('-').unwrap();
            (c.to_string(), s.parse::<u64>().unwrap() - 1, e.parse::<u64>().unwrap())
        }).collect();
        assert_eq!(oracle.len(), 31);

        // --- score helper: real detect_edges + decompose_families, unchanged; report how many distinct
        //     families the 31 oracle loci fall into, and whether any oracle-containing family also
        //     contains a non-oracle rep (a real, on-substrate false-merge check, not a synthetic one). ---
        let score = |reps: &[DenovoTranscript], label: &str| -> (usize, usize, usize) {
            let dp = DetectParams::default();
            let t0 = std::time::Instant::now();
            let edges = detect_edges(reps, &dp);
            eprintln!("[stage_c/e:{label}] detect_edges took {:?}", t0.elapsed());
            let families = decompose_families(&edges, &SplitParams::default());
            let mut rep_family: Vec<Option<usize>> = vec![None; reps.len()];
            for (fi, fam) in families.iter().enumerate() {
                for &m in &fam.members { rep_family[m] = Some(fi); }
            }
            let mut oracle_families: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
            let mut oracle_covered = 0usize;
            for (oc, os_, oe) in &oracle {
                if let Some((ri, _)) = reps.iter().enumerate()
                    .find(|(_, r)| &r.chrom == oc && r.end > *os_ && r.start < *oe) {
                    oracle_covered += 1;
                    if let Some(fi) = rep_family[ri] { oracle_families.insert(fi); }
                }
            }
            let mut false_merge_families = 0usize;
            for &fi in &oracle_families {
                let has_foreign = families[fi].members.iter().any(|&m| {
                    !oracle.iter().any(|(oc, os_, oe)|
                        &reps[m].chrom == oc && reps[m].end > *os_ && reps[m].start < *oe)
                });
                if has_foreign { false_merge_families += 1; }
            }
            eprintln!(
                "[stage_c:{label}] {} reps, {} edges, {} families total; oracle: {}/31 covered, spread \
                 across {} families; {}/{} oracle-containing families also contain a non-oracle member",
                reps.len(), edges.len(), families.len(), oracle_covered, oracle_families.len(),
                false_merge_families, oracle_families.len()
            );
            (oracle_covered, oracle_families.len(), false_merge_families)
        };

        // --- BASELINE: the real, unmodified footprint-off RNA reps alone ---
        let baseline_result = score(&baseline, "baseline (RNA only, uncorrected)");

        // --- PROPOSAL #1: correct spans of baseline reps overlapping an oracle locus ---
        let sedef_text = std::fs::read_to_string(sedef_bed).unwrap();
        let sd_pairs = SdPairs::from_bed_str(&sedef_text);
        let contigs: HashSet<String> = baseline.iter().map(|r| r.chrom.clone())
            .chain(oracle.iter().map(|(c, _, _)| c.clone())).collect();
        let genome = GenomeIndex::from_fasta_contigs(fa, &contigs).expect("genome index");
        let mut corrected = baseline.clone();
        let mut n_corrected = 0usize;
        for r in corrected.iter_mut() {
            let overlaps_oracle = oracle.iter().any(|(oc, os_, oe)| &r.chrom == oc && r.end > *os_ && r.start < *oe);
            if !overlaps_oracle { continue; }
            if let Some((cs, ce)) = sd_pairs.single_span_core(&r.chrom, r.start, r.end, 15_000, 2) {
                if let Some(seq) = genome.fetch_sequence(&r.chrom, cs, ce) {
                    r.start = cs; r.end = ce; r.introns = vec![]; r.seq = seq;
                    n_corrected += 1;
                }
            }
        }
        eprintln!("[stage_c] proposal #1: corrected {n_corrected} of the baseline reps overlapping an oracle locus");
        let p1_result = score(&corrected, "proposal #1 only (spans corrected)");

        // --- PROPOSAL #2: add SD-seeded DNA-only reps (already-validated §6j5 Stage B windows, reused
        //     from `sd_windows` loaded above -- same file, same region used to bound the baseline) ---
        let t_genome_reps = std::time::Instant::now();
        let dna_reps = genome_reps(fa, &sd_windows, &GenomeRepParams::default()).expect("genome_reps");
        eprintln!("[stage_c/e] genome_reps took {:?}", t_genome_reps.elapsed());
        eprintln!("[stage_c] proposal #2: {} SD-seeded DNA-only reps generated", dna_reps.len());

        // --- COMBINED: corrected RNA reps + DNA-only reps, deduped where they land on the same locus ---
        let mut combined = corrected.clone();
        let mut n_dup_skipped = 0usize;
        for d in dna_reps {
            let dup = combined.iter().any(|r| r.chrom == d.chrom &&
                r.start.max(d.start) < r.end.min(d.end) &&
                (r.end.min(d.end) - r.start.max(d.start)) as f64 / (d.end - d.start).max(1) as f64 > 0.8);
            if dup { n_dup_skipped += 1; continue; }
            combined.push(d);
        }
        eprintln!(
            "[stage_c] combined: {} total reps ({} DNA-only reps skipped as >80% overlap-duplicates)",
            combined.len(), n_dup_skipped
        );
        let combined_result = score(&combined, "COMBINED (proposal #1 + #2)");

        eprintln!(
            "[stage_c/e] SUMMARY [{label_suffix}] (oracle_covered/oracle_families/false_merge_families): \
             baseline {:?}, prop#1-only {:?}, COMBINED {:?}",
            baseline_result, p1_result, combined_result
        );
    }

    /// DIAGNOSTIC (false-merge evidence dump for §6j8's "false merge" families, docs/o1_ledger.md §6j9). Runs
    /// ONLY the `baseline` and `proposal #1` cells of `stage_c_body` (same loading, bounding and proposal-#1
    /// correction; `genome_reps` is never called), scoring each with the real `detect_edges` (default edge core)
    /// + `decompose_families`, and writes `$RUSTLE_FM_OUT/<cell>/` (reps.tsv, reps.fa, candidates.tsv, edges.tsv,
    /// families.tsv, summary.tsv). Other phases (`RUSTLE_FM_PHASE`): `bridge` computes the exact POA core for a
    /// pair list (`$RUSTLE_FM_PAIRS`, one process per pair recommended), `pairs` runs `candidate_pairs` + LCS on
    /// an arbitrary rep set, `decompose` partitions an edge list. No threads are spawned by any phase.
    /// `#[ignore]`d: needs real data; run under a shell `timeout`.
    #[test]
    #[ignore]
    fn stage_f_falsemerge_evidence_dump() {
        let phase = std::env::var("RUSTLE_FM_PHASE").unwrap_or_else(|_| "dump".into());
        let root = std::env::var("RUSTLE_FM_OUT")
            .unwrap_or_else(|_| "/mnt/linuxdisk/home/juanfraitu/o1_falsemerge".into());
        match phase.as_str() {
            "dump" => fm_dump(&root),
            "bridge" => fm_bridge(&root),
            "pairs" => fm_pairs(&root),
            "decompose" => fm_decompose(),
            "sd_partition" => fm_sd_partition(),
            other => panic!("unknown RUSTLE_FM_PHASE {other}"),
        }
    }

    /// SD-atom phase (prereg Addendum J2): `families_from_edges` on atoms (`$RUSTLE_FM_NODES`: `idx chrom start end`,
    /// header) and weighted edges (`$RUSTLE_FM_EDGES`: `i j identity coverage`, header), floors 0.80 / 0.50,
    /// gamma 0.20, min_copies 2, min_reads 0. Writes a `copies.tsv`-shaped table to `$RUSTLE_FM_FAMILIES_OUT`.
    fn fm_sd_partition() {
        use crate::vg_family::denovo_pipeline::families_from_edges;
        use std::fmt::Write as _;
        let nodes_path = std::env::var("RUSTLE_FM_NODES").expect("RUSTLE_FM_NODES");
        let edges_path = std::env::var("RUSTLE_FM_EDGES").expect("RUSTLE_FM_EDGES");
        let out_path = std::env::var("RUSTLE_FM_FAMILIES_OUT").expect("RUSTLE_FM_FAMILIES_OUT");
        let mut reps = Vec::new();
        for line in std::fs::read_to_string(&nodes_path).unwrap().lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            let (Ok(s), Ok(e)) = (f[2].parse::<u64>(), f[3].parse::<u64>()) else { panic!("bad node row {line}") };
            assert_eq!(f[0].parse::<usize>().unwrap(), reps.len(), "node indices must be 0..n in order");
            reps.push(DenovoTranscript {
                tid: format!("DNA_{}_{s}_{e}", f[1]),
                chrom: f[1].to_string(), start: s, end: e, n_reads: 1, strand: '+',
                introns: vec![], seq: vec![], distinguishing_uniq: 0, core_bp: 0, stub: false, tes: None,
            });
        }
        let mut edges = Vec::new();
        for line in std::fs::read_to_string(&edges_path).unwrap().lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            edges.push((f[0].parse::<usize>().unwrap(), f[1].parse::<usize>().unwrap(),
                        f[2].parse::<f64>().unwrap(), f[3].parse::<f64>().unwrap()));
        }
        let n_nodes = reps.len();
        let mut fams = families_from_edges(reps, &edges, 0.80, 0.50, 0.20, 2, 0);
        fams.sort_by(|a, b| b.len().cmp(&a.len()).then((&a[0].chrom, a[0].start).cmp(&(&b[0].chrom, b[0].start))));
        let mut out = String::from("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\n");
        for (fi, fam) in fams.iter().enumerate() {
            for (ci, t) in fam.iter().enumerate() {
                writeln!(out, "SDFAM{fi}\t{ci}\t{}\t{}\t{}\t{}\t1\t+\t1\t{}-{}\tNA", t.tid, t.chrom, t.start, t.end, t.start, t.end).unwrap();
            }
        }
        std::fs::write(&out_path, out).unwrap();
        eprintln!("[stage_f:sd_partition] {n_nodes} atoms, {} edges -> {} families", edges.len(), fams.len());
    }

    fn fm_dump(root: &str) {
        use crate::vg_family::annotation_families::SdPairs;
        use crate::vg_family::family_detect::{detect_edges, DetectParams};

        let fa = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta";
        let copies_tsv = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.tsv";
        let copies_fa = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.fa";
        let sedef_bed = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed";
        let sd_windows_bed = "/mnt/linuxdisk/home/juanfraitu/o1_fromgenome_sd/npip_seeded_windows2.bed";
        let oracle_path = "/mnt/linuxdisk/home/juanfraitu/o1_oracle/npip31.regions";
        for p in [fa, copies_tsv, copies_fa, sedef_bed, sd_windows_bed, oracle_path] {
            if std::fs::metadata(p).is_err() { eprintln!("required real-data file {p} absent; skip"); return; }
        }
        eprintln!("[stage_f] root={root}");

        // --- loading: identical to stage_c_body (exon strings kept alongside for the dump) ---
        let tsv_text = std::fs::read_to_string(copies_tsv).unwrap();
        let fa_text = std::fs::read_to_string(copies_fa).unwrap();
        let seqs: Vec<Vec<u8>> = fa_text.lines().filter(|l| !l.starts_with('>')).map(|l| l.as_bytes().to_vec()).collect();
        let mut loaded: Vec<(DenovoTranscript, String, usize)> = Vec::new();
        for (i, line) in tsv_text.lines().skip(1).enumerate() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 11 { continue; }
            let chrom = f[3].to_string();
            let start: u64 = f[4].parse().unwrap();
            let end: u64 = f[5].parse().unwrap();
            let strand = f[7].chars().next().unwrap_or('+');
            let n_reads: u32 = f[8].parse().unwrap_or(1);
            let exons: Vec<(u64, u64)> = f[9].split(',').filter_map(|e| {
                let (s, en) = e.split_once('-')?;
                Some((s.parse().ok()?, en.parse().ok()?))
            }).collect();
            let introns: Vec<(u64, u64)> = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
            loaded.push((DenovoTranscript {
                tid: f[2].to_string(), chrom, start, end, n_reads, strand,
                introns, seq: seqs.get(i).cloned().unwrap_or_default(), distinguishing_uniq: 0,
                core_bp: 0, stub: false, tes: None,
            }, f[9].to_string(), exons.len()));
        }
        let mut sd_windows: Vec<(String, u64, u64)> = Vec::new();
        for line in std::fs::read_to_string(sd_windows_bed).unwrap().lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() >= 3 { sd_windows.push((f[0].to_string(), f[1].parse().unwrap(), f[2].parse().unwrap())); }
        }
        let n_before_bound = loaded.len();
        loaded.retain(|(r, _, _)| sd_windows.iter().any(|(c, s, e)| &r.chrom == c && r.end > *s && r.start < *e));
        let exon_strs: Vec<String> = loaded.iter().map(|x| x.1.clone()).collect();
        let n_exons: Vec<usize> = loaded.iter().map(|x| x.2).collect();
        let baseline: Vec<DenovoTranscript> = loaded.into_iter().map(|x| x.0).collect();
        eprintln!("[stage_f] bounded {n_before_bound} -> {} baseline reps", baseline.len());

        let oracle_txt = std::fs::read_to_string(oracle_path).unwrap();
        let oracle: Vec<(String, u64, u64)> = oracle_txt.lines().filter(|l| !l.trim().is_empty()).map(|l| {
            let (c, rest) = l.trim().split_once(':').unwrap();
            let (s, e) = rest.split_once('-').unwrap();
            (c.to_string(), s.parse::<u64>().unwrap() - 1, e.parse::<u64>().unwrap())
        }).collect();
        assert_eq!(oracle.len(), 31);

        let dp = DetectParams::default();
        let base_edges = detect_edges(&baseline, &dp);
        eprintln!("[stage_f:baseline] detect_edges ({:?} core): {} edges", dp.edge_core, base_edges.len());

        // --- PROPOSAL #1: identical to stage_c_body ---
        let sedef_text = std::fs::read_to_string(sedef_bed).unwrap();
        let sd_pairs = SdPairs::from_bed_str(&sedef_text);
        let contigs: HashSet<String> = baseline.iter().map(|r| r.chrom.clone())
            .chain(oracle.iter().map(|(c, _, _)| c.clone())).collect();
        let genome = GenomeIndex::from_fasta_contigs(fa, &contigs).expect("genome index");
        let mut corrected = baseline.clone();
        let mut is_corrected = vec![false; corrected.len()];
        for (ri, r) in corrected.iter_mut().enumerate() {
            let overlaps_oracle = oracle.iter().any(|(oc, os_, oe)| &r.chrom == oc && r.end > *os_ && r.start < *oe);
            if !overlaps_oracle { continue; }
            if let Some((cs, ce)) = sd_pairs.single_span_core(&r.chrom, r.start, r.end, 15_000, 2) {
                if let Some(seq) = genome.fetch_sequence(&r.chrom, cs, ce) {
                    r.start = cs; r.end = ce; r.introns = vec![]; r.seq = seq;
                    is_corrected[ri] = true;
                }
            }
        }
        eprintln!("[stage_f] proposal #1: corrected {} reps", is_corrected.iter().filter(|x| **x).count());
        let p1_edges = detect_edges(&corrected, &dp);
        eprintln!("[stage_f:proposal1] detect_edges ({:?} core): {} edges", dp.edge_core, p1_edges.len());
        let base_kinds: Vec<&'static str> = n_exons.iter().map(|&n| if n > 1 { "rna_spliced" } else { "rna_single_exon" }).collect();
        let p1_kinds: Vec<&'static str> = (0..corrected.len())
            .map(|i| if is_corrected[i] { "genomic_span_p1corrected" } else { base_kinds[i] }).collect();
        fm_write_cell(root, "baseline", &baseline, &baseline, &exon_strs, &base_kinds, &is_corrected.iter().map(|_| false).collect::<Vec<_>>(),
            &oracle, base_edges, &dp);
        fm_write_cell(root, "proposal1", &corrected, &baseline, &exon_strs, &p1_kinds, &is_corrected,
            &oracle, p1_edges, &dp);
    }

    #[allow(clippy::too_many_arguments)]
    fn fm_write_cell(
        root: &str,
        cell: &str,
        reps: &[DenovoTranscript],
        orig: &[DenovoTranscript],
        exon_strs: &[String],
        kinds: &[&'static str],
        is_corrected: &[bool],
        oracle: &[(String, u64, u64)],
        edges: Vec<(usize, usize, f64)>,
        dp: &crate::vg_family::family_detect::DetectParams,
    ) {
        use crate::vg_family::family_detect::candidate_pairs;
        use crate::vg_family::family_detect::family_graph::{longest_common_substring, upper_cow};
        use crate::vg_family::family_detect::family_split::{connected_components, decompose_families, SplitParams};
        use crate::vg_family::seq_utils::reverse_complement;
        use std::fmt::Write as _;

        let dir = format!("{root}/{cell}");
        std::fs::create_dir_all(&dir).unwrap();
        let n = reps.len();
        let pairs = candidate_pairs(reps, dp);

        // exact-substring (LCS) values for both orientations of every candidate pair (cheap, deterministic).
        let lcs: Vec<(usize, usize)> = pairs.iter().map(|&(a, b)| {
            let au = upper_cow(&reps[a].seq);
            let bu = upper_cow(&reps[b].seq);
            (longest_common_substring(&au, &bu), longest_common_substring(&au, &reverse_complement(&bu)))
        }).collect();

        let edges_source = "real_detect_edges";
        let edge_val: HashMap<(usize, usize), f64> = edges.iter().map(|&(a, b, v)| ((a, b), v)).collect();

        let families = decompose_families(&edges, &SplitParams::default());
        let comps = connected_components(&edges, 2);
        let mut rep_comp: Vec<Option<usize>> = vec![None; n];
        for (ci, c) in comps.iter().enumerate() { for &m in c { rep_comp[m] = Some(ci); } }
        let mut rep_family: Vec<Option<usize>> = vec![None; n];
        for (fi, fam) in families.iter().enumerate() { for &m in &fam.members { rep_family[m] = Some(fi); } }
        let overlaps = |r: &DenovoTranscript| -> Vec<usize> {
            oracle.iter().enumerate().filter(|(_, (oc, os_, oe))| &r.chrom == oc && r.end > *os_ && r.start < *oe)
                .map(|(oi, _)| oi).collect()
        };
        let rep_oracle: Vec<Vec<usize>> = reps.iter().map(|r| overlaps(r)).collect();
        // stage_c's oracle-family rule: for each oracle locus, the FIRST rep (by index) overlapping it.
        let mut first_for: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut oracle_families: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
        let mut oracle_covered = 0usize;
        for (oi, (oc, os_, oe)) in oracle.iter().enumerate() {
            if let Some((ri, _)) = reps.iter().enumerate().find(|(_, r)| &r.chrom == oc && r.end > *os_ && r.start < *oe) {
                oracle_covered += 1;
                first_for[ri].push(oi);
                if let Some(fi) = rep_family[ri] { oracle_families.insert(fi); }
            }
        }
        let join = |v: &[usize]| v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",");

        // reps.tsv + reps.fa
        let mut rt = String::from("idx\ttid\tchrom\tstart\tend\tstrand\tn_reads\tseq_len\tkind\tcorrected\torig_start\torig_end\tn_exons_orig\texons_orig\toracle_idx\toracle_first_rep_for\tfamily_id\tfamily_class\tcomponent_id\n");
        let mut rf = String::new();
        for (i, r) in reps.iter().enumerate() {
            let fam = rep_family[i];
            writeln!(rt, "{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.tid, r.chrom, r.start, r.end, r.strand, r.n_reads, r.seq.len(), kinds[i], is_corrected[i] as u8,
                orig[i].start, orig[i].end, exon_strs[i].split(',').count(), exon_strs[i],
                join(&rep_oracle[i]), join(&first_for[i]),
                fam.map(|f| f.to_string()).unwrap_or_default(),
                fam.map(|f| format!("{:?}", families[f].class)).unwrap_or_default(),
                rep_comp[i].map(|c| c.to_string()).unwrap_or_default()).unwrap();
            writeln!(rf, ">{i}|{}|{}:{}-{}|{}\n{}", r.tid, r.chrom, r.start, r.end, kinds[i], String::from_utf8_lossy(&r.seq)).unwrap();
        }
        std::fs::write(format!("{dir}/reps.tsv"), rt).unwrap();
        std::fs::write(format!("{dir}/reps.fa"), rf).unwrap();

        // candidates.tsv (every candidate pair) and edges.tsv (confirmed subset), same columns.
        let header = "i\tj\ttid_i\ttid_j\tkind_i\tkind_j\tlen_i\tlen_j\tmin_len\tmax_len\tover_len_cap\tconfirmed\tcore_frac_run\tlcs_fwd_bp\tlcs_rc_bp\tlcs_orientation_emulated\tlcs_bp_emulated\tlcs_over_min_len\tlcs_over_max_len\tlcs_emulated_core_frac\tlcs_emulated_pass\trun_value_equals_lcs_emulation\tfamily_i\tfamily_j\tsame_family\tcomponent_i\tcomponent_j\n";
        let mut ct = String::from(header);
        let mut et = String::from(header);
        let mut n_eq_lcs = 0usize;
        for (k, &(a, b)) in pairs.iter().enumerate() {
            let (la, lb) = (reps[a].seq.len(), reps[b].seq.len());
            let (mn, mx) = (la.min(lb), la.max(lb));
            let (lf, lr) = lcs[k];
            let fwd_frac = lf as f64 / mn as f64;
            let (orient, lbp) = if fwd_frac < dp.t_core && (lr as f64 / mn as f64) > fwd_frac { ("rc", lr) } else { ("fwd", lf) };
            let emu = lbp as f64 / mn as f64;
            let conf = edge_val.get(&(a, b)).copied();
            let eq = conf.map(|v| (v - emu).abs() < 1e-12);
            if eq == Some(true) { n_eq_lcs += 1; }
            let (fa_, fb_) = (rep_family[a], rep_family[b]);
            let row = format!("{a}\t{b}\t{}\t{}\t{}\t{}\t{la}\t{lb}\t{mn}\t{mx}\t{}\t{}\t{}\t{lf}\t{lr}\t{orient}\t{lbp}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                reps[a].tid, reps[b].tid, kinds[a], kinds[b], (mx > dp.len_cap) as u8, conf.is_some() as u8,
                conf.map(|v| format!("{v:.6}")).unwrap_or_default(),
                lbp as f64 / mn as f64, lbp as f64 / mx as f64, emu, (emu >= dp.t_core) as u8,
                eq.map(|x| (x as u8).to_string()).unwrap_or_default(),
                fa_.map(|f| f.to_string()).unwrap_or_default(), fb_.map(|f| f.to_string()).unwrap_or_default(),
                (fa_.is_some() && fa_ == fb_) as u8,
                rep_comp[a].map(|c| c.to_string()).unwrap_or_default(), rep_comp[b].map(|c| c.to_string()).unwrap_or_default());
            ct.push_str(&row);
            if conf.is_some() { et.push_str(&row); }
        }
        std::fs::write(format!("{dir}/candidates.tsv"), ct).unwrap();
        std::fs::write(format!("{dir}/edges.tsv"), et).unwrap();

        // families.tsv
        let mut ft = String::from("family_id\tclass\tn_members\tn_edges_induced\tdensity\tavg_core_recip\tcomponent_id\tcomponent_size\tmembers\tmember_tids\toracle_overlapping_members\toracle_first_rep_members\toracle_loci\tforeign_members\tis_oracle_family_stage_c\tfalse_merge_flag_stage_c\n");
        let mut false_merge_families = 0usize;
        for (fi, fam) in families.iter().enumerate() {
            let ovl: Vec<usize> = fam.members.iter().copied().filter(|&m| !rep_oracle[m].is_empty()).collect();
            let firsts: Vec<usize> = fam.members.iter().copied().filter(|&m| !first_for[m].is_empty()).collect();
            let foreign: Vec<usize> = fam.members.iter().copied().filter(|&m| rep_oracle[m].is_empty()).collect();
            let mut loci: Vec<usize> = fam.members.iter().flat_map(|&m| rep_oracle[m].iter().copied()).collect();
            loci.sort_unstable(); loci.dedup();
            let is_of = oracle_families.contains(&fi);
            let fm = is_of && !foreign.is_empty();
            if fm { false_merge_families += 1; }
            let ci = rep_comp[fam.members[0]];
            writeln!(ft, "{fi}\t{:?}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                fam.class, fam.members.len(), fam.stats.n_edges, fam.stats.density, fam.stats.avg_core_recip,
                ci.map(|c| c.to_string()).unwrap_or_default(), ci.map(|c| comps[c].len().to_string()).unwrap_or_default(),
                join(&fam.members), fam.members.iter().map(|&m| reps[m].tid.as_str()).collect::<Vec<_>>().join(","),
                join(&ovl), join(&firsts), join(&loci), join(&foreign), is_of as u8, fm as u8).unwrap();
        }
        std::fs::write(format!("{dir}/families.tsv"), ft).unwrap();

        let mut st = String::new();
        writeln!(st, "key\tvalue").unwrap();
        writeln!(st, "cell\t{cell}").unwrap();
        writeln!(st, "edge_core\t{:?}", dp.edge_core).unwrap();
        writeln!(st, "edges_source\t{edges_source}").unwrap();
        writeln!(st, "n_reps\t{n}").unwrap();
        writeln!(st, "n_candidate_pairs\t{}", pairs.len()).unwrap();
        writeln!(st, "n_candidate_pairs_over_len_cap\t{}", pairs.iter().filter(|&&(a, b)| reps[a].seq.len().max(reps[b].seq.len()) > dp.len_cap).count()).unwrap();
        writeln!(st, "n_edges\t{}", edges.len()).unwrap();
        writeln!(st, "n_families\t{}", families.len()).unwrap();
        writeln!(st, "oracle_covered\t{oracle_covered}/31").unwrap();
        writeln!(st, "oracle_families\t{}", oracle_families.len()).unwrap();
        writeln!(st, "false_merge_families\t{false_merge_families}/{}", oracle_families.len()).unwrap();
        writeln!(st, "n_edges_value_equal_lcs_emulation\t{n_eq_lcs}").unwrap();
        writeln!(st, "n_candidates_lcs_emulated_pass\t{}", pairs.iter().enumerate().filter(|(k, &(a, b))| {
            let mn = reps[a].seq.len().min(reps[b].seq.len()) as f64;
            let (lf, lr) = lcs[*k];
            (lf.max(lr) as f64 / mn) >= dp.t_core
        }).count()).unwrap();
        std::fs::write(format!("{dir}/summary.tsv"), &st).unwrap();
        eprintln!("[stage_f:{cell}] {} reps, {} candidates, {} edges ({edges_source}), {} families; oracle {}/31 in {} families; false-merge {}/{}",
            n, pairs.len(), edges.len(), families.len(), oracle_covered, oracle_families.len(), false_merge_families, oracle_families.len());
    }

    /// Pairs phase: production `candidate_pairs` over an arbitrary rep set (copies.tsv columns + FASTA in the
    /// same row order, `$RUSTLE_FM_REPS_TSV` / `$RUSTLE_FM_REPS_FA`), with serial LCS for both orientations.
    /// Writes `<root>/<$RUSTLE_FM_CELL>/{reps.fa,candidates.tsv}` so the `bridge` phase can run on it. No threads.
    fn fm_pairs(root: &str) {
        use crate::vg_family::family_detect::{candidate_pairs, DetectParams};
        use crate::vg_family::family_detect::family_graph::{longest_common_substring, upper_cow};
        use crate::vg_family::seq_utils::reverse_complement;
        use std::fmt::Write as _;
        let tsv = std::env::var("RUSTLE_FM_REPS_TSV").expect("RUSTLE_FM_REPS_TSV");
        let fa = std::env::var("RUSTLE_FM_REPS_FA").expect("RUSTLE_FM_REPS_FA");
        let cell = std::env::var("RUSTLE_FM_CELL").expect("RUSTLE_FM_CELL");
        let seqs: Vec<Vec<u8>> = std::fs::read_to_string(&fa).unwrap().lines()
            .filter(|l| !l.starts_with('>')).map(|l| l.as_bytes().to_vec()).collect();
        let mut reps: Vec<DenovoTranscript> = Vec::new();
        for (i, line) in std::fs::read_to_string(&tsv).unwrap().lines().skip(1).enumerate() {
            let f: Vec<&str> = line.split('\t').collect();
            let exons: Vec<(u64, u64)> = f[9].split(',').filter_map(|e| {
                let (s, en) = e.split_once('-')?;
                Some((s.parse().ok()?, en.parse().ok()?))
            }).collect();
            reps.push(DenovoTranscript {
                tid: f[2].to_string(), chrom: f[3].to_string(), start: f[4].parse().unwrap(), end: f[5].parse().unwrap(),
                n_reads: f[8].parse().unwrap_or(1), strand: f[7].chars().next().unwrap_or('+'),
                introns: exons.windows(2).map(|w| (w[0].1, w[1].0)).collect(),
                seq: seqs[i].clone(), distinguishing_uniq: 0, core_bp: 0, stub: false, tes: None,
            });
        }
        let dp = DetectParams::default();
        let pairs = candidate_pairs(&reps, &dp);
        eprintln!("[stage_f:pairs] {} reps -> {} candidate pairs", reps.len(), pairs.len());
        let dir = format!("{root}/{cell}");
        std::fs::create_dir_all(&dir).unwrap();
        let mut rf = String::new();
        for (i, r) in reps.iter().enumerate() {
            writeln!(rf, ">{i}|{}\n{}", r.tid, String::from_utf8_lossy(&r.seq)).unwrap();
        }
        std::fs::write(format!("{dir}/reps.fa"), rf).unwrap();
        let mut ct = String::from("i\tj\ttid_i\ttid_j\tlen_i\tlen_j\tmin_len\tmax_len\tlcs_fwd_bp\tlcs_rc_bp\n");
        for &(a, b) in &pairs {
            let au = upper_cow(&reps[a].seq);
            let bu = upper_cow(&reps[b].seq);
            let (la, lb) = (reps[a].seq.len(), reps[b].seq.len());
            let lf = longest_common_substring(&au, &bu);
            let lr = longest_common_substring(&au, &reverse_complement(&bu));
            writeln!(ct, "{a}\t{b}\t{}\t{}\t{la}\t{lb}\t{}\t{}\t{lf}\t{lr}", reps[a].tid, reps[b].tid, la.min(lb), la.max(lb)).unwrap();
        }
        std::fs::write(format!("{dir}/candidates.tsv"), ct).unwrap();
    }

    /// Decompose phase: the shipped `decompose_families(edges, SplitParams::default())` on an edge list
    /// (`$RUSTLE_FM_EDGES`: `i<TAB>j<TAB>core` rows, header optional), writing `family_id<TAB>class<TAB>members`
    /// to `$RUSTLE_FM_FAMILIES_OUT`. Lets two edge definitions be partitioned by the identical code.
    fn fm_decompose() {
        use crate::vg_family::family_detect::family_split::{decompose_families, SplitParams};
        use std::fmt::Write as _;
        let edges_path = std::env::var("RUSTLE_FM_EDGES").expect("RUSTLE_FM_EDGES");
        let out_path = std::env::var("RUSTLE_FM_FAMILIES_OUT").expect("RUSTLE_FM_FAMILIES_OUT");
        let mut edges: Vec<(usize, usize, f64)> = Vec::new();
        for line in std::fs::read_to_string(&edges_path).unwrap().lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 3 { continue; }
            let (Ok(i), Ok(j), Ok(v)) = (f[0].parse::<usize>(), f[1].parse::<usize>(), f[2].parse::<f64>()) else { continue };
            edges.push((i.min(j), i.max(j), v));
        }
        edges.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
        let families = decompose_families(&edges, &SplitParams::default());
        let mut out = String::from("family_id\tclass\tmembers\n");
        for (fi, fam) in families.iter().enumerate() {
            writeln!(out, "{fi}\t{:?}\t{}", fam.class,
                fam.members.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(",")).unwrap();
        }
        std::fs::write(&out_path, out).unwrap();
        eprintln!("[stage_f:decompose] {} edges -> {} families", edges.len(), families.len());
    }

    /// Bridge phase: exact POA-core `confirm_edge` values (EdgeCore::Poa, no budget) for a list of pairs, serial.
    fn fm_bridge(root: &str) {
        use crate::vg_family::family_detect::{confirm_edge, DetectParams, LEN_CAP, T_CORE};
        use crate::vg_family::family_detect::family_graph::{contiguous_core_coverage_bounded_with, longest_common_substring, upper_cow, EDGE_CONFIRM_ASTAR};
        use crate::vg_family::seq_utils::reverse_complement;
        use std::io::Write;
        let pairs_path = std::env::var("RUSTLE_FM_PAIRS").unwrap_or_else(|_| format!("{root}/bridging_pairs.list"));
        let out_path = std::env::var("RUSTLE_FM_BRIDGE_OUT").unwrap_or_else(|_| format!("{root}/bridging_production.tsv"));
        let mut cache: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
        let mut out = std::fs::File::create(&out_path).unwrap();
        writeln!(out, "cell\ti\tj\tlen_i\tlen_j\tmin_len\tmax_len\tproduction_path\tfwd_value\tfwd_secs\trc_value\trc_secs\tproduction_core_frac\tproduction_core_bp\tcore_over_max_len\tpasses_tcore\tconfirm_edge_agrees\tlcs_fwd_bp\tlcs_rc_bp").unwrap();
        let dp = DetectParams { edge_core: crate::vg_family::family_detect::EdgeCore::Poa, ..DetectParams::default() };
        for line in std::fs::read_to_string(&pairs_path).unwrap().lines() {
            if line.trim().is_empty() || line.starts_with('#') { continue; }
            let f: Vec<&str> = line.split('\t').collect();
            let (cell, i, j): (&str, usize, usize) = (f[0], f[1].parse().unwrap(), f[2].parse().unwrap());
            let seqs = cache.entry(cell.to_string()).or_insert_with(|| {
                let txt = std::fs::read_to_string(format!("{root}/{cell}/reps.fa")).unwrap();
                txt.lines().filter(|l| !l.starts_with('>')).map(|l| l.as_bytes().to_vec()).collect()
            });
            let (a, b) = (seqs[i].clone(), seqs[j].clone());
            let (mn, mx) = (a.len().min(b.len()), a.len().max(b.len()));
            let au = upper_cow(&a).into_owned();
            let bu = upper_cow(&b).into_owned();
            let path = if mx > LEN_CAP { "LCS (len_cap)" } else { "poasta exact" };
            eprintln!("[stage_f:bridge] START {cell} {i} {j} len {} x {} path={path}", a.len(), b.len());
            let t = std::time::Instant::now();
            let fwd = contiguous_core_coverage_bounded_with(&au, &bu, LEN_CAP, EDGE_CONFIRM_ASTAR);
            let fwd_secs = t.elapsed().as_secs_f64();
            let mut cr = fwd;
            let (mut rcv, mut rcs) = (String::from("not_run"), String::new());
            if fwd < T_CORE {
                let t2 = std::time::Instant::now();
                let rc = contiguous_core_coverage_bounded_with(&au, &reverse_complement(&bu), LEN_CAP, EDGE_CONFIRM_ASTAR);
                rcs = format!("{:.3}", t2.elapsed().as_secs_f64());
                rcv = format!("{rc:.6}");
                if rc > cr { cr = rc; }
            }
            let ce = confirm_edge(&a, &b, &dp);
            let agrees = match ce { Some(v) => (v - cr).abs() < 1e-12 && cr >= T_CORE, None => cr < T_CORE };
            let core_bp = (cr * mn as f64).round() as usize;
            let lf = longest_common_substring(&au, &bu);
            let lr = longest_common_substring(&au, &reverse_complement(&bu));
            writeln!(out, "{cell}\t{i}\t{j}\t{}\t{}\t{mn}\t{mx}\t{path}\t{fwd:.6}\t{fwd_secs:.3}\t{rcv}\t{rcs}\t{cr:.6}\t{core_bp}\t{:.6}\t{}\t{}\t{lf}\t{lr}",
                a.len(), b.len(), core_bp as f64 / mx as f64, (cr >= T_CORE) as u8, agrees as u8).unwrap();
            out.flush().unwrap();
            eprintln!("[stage_f:bridge] DONE  {cell} {i} {j} cr={cr:.4} in {:.1}s agrees={agrees}", t.elapsed().as_secs_f64());
        }
    }

    /// DIAGNOSTIC (docs/o1_ledger.md §6j7): §6j6's stage_c hung on the SAME bounded (391->42) rep set this
    /// re-loads. Instead of running `detect_edges`'s rayon-parallel batch (which hides WHICH pair is slow
    /// behind aggregate CPU%), this calls `candidate_pairs` once then `confirm_edge` SEQUENTIALLY, one pair
    /// at a time, printing a start-line BEFORE and a finish-line AFTER each call so a hang is visible as an
    /// unmatched start-line in the log even if the whole process is later killed by an external `timeout`.
    /// `#[ignore]`d: needs real gorilla data; run under a shell `timeout`, never bare.
    #[test]
    #[ignore]
    fn stage_d_profile_confirm_edge_per_pair_sequentially() {
        use crate::vg_family::family_detect::{candidate_pairs, confirm_edge, DetectParams};
        use std::io::Write;
        use std::time::Instant;

        let copies_tsv = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.tsv";
        let copies_fa = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.fa";
        let sd_windows_bed = "/mnt/linuxdisk/home/juanfraitu/o1_fromgenome_sd/npip_seeded_windows2.bed";
        for p in [copies_tsv, copies_fa, sd_windows_bed] {
            if std::fs::metadata(p).is_err() { eprintln!("required real-data file {p} absent; skip"); return; }
        }

        // --- identical loading + bounding to §6j6's stage_c, so this reproduces the EXACT same hang ---
        let tsv_text = std::fs::read_to_string(copies_tsv).unwrap();
        let fa_text = std::fs::read_to_string(copies_fa).unwrap();
        let seqs: Vec<Vec<u8>> = fa_text.lines().filter(|l| !l.starts_with('>')).map(|l| l.as_bytes().to_vec()).collect();
        let mut baseline: Vec<DenovoTranscript> = Vec::new();
        for (i, line) in tsv_text.lines().skip(1).enumerate() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 11 { continue; }
            let chrom = f[3].to_string();
            let start: u64 = f[4].parse().unwrap();
            let end: u64 = f[5].parse().unwrap();
            let strand = f[7].chars().next().unwrap_or('+');
            let n_reads: u32 = f[8].parse().unwrap_or(1);
            let exons: Vec<(u64, u64)> = f[9].split(',').filter_map(|e| {
                let (s, en) = e.split_once('-')?;
                Some((s.parse().ok()?, en.parse().ok()?))
            }).collect();
            let introns: Vec<(u64, u64)> = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
            baseline.push(DenovoTranscript {
                tid: f[2].to_string(), chrom, start, end, n_reads, strand,
                introns, seq: seqs.get(i).cloned().unwrap_or_default(), distinguishing_uniq: 0,
                core_bp: 0, stub: false, tes: None,
            });
        }
        let mut sd_windows: Vec<(String, u64, u64)> = Vec::new();
        for line in std::fs::read_to_string(sd_windows_bed).unwrap().lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() >= 3 { sd_windows.push((f[0].to_string(), f[1].parse().unwrap(), f[2].parse().unwrap())); }
        }
        baseline.retain(|r| sd_windows.iter().any(|(c, s, e)| &r.chrom == c && r.end > *s && r.start < *e));
        eprintln!("[stage_d] {} bounded reps (should match §6j6's 42)", baseline.len());
        for (i, r) in baseline.iter().enumerate() {
            eprintln!("[stage_d]   rep[{i}] {}:{}-{} len={} n_reads={}", r.chrom, r.start, r.end, r.seq.len(), r.n_reads);
        }
        std::io::stderr().flush().ok();

        let dp = DetectParams { edge_core: crate::vg_family::family_detect::EdgeCore::Poa, ..DetectParams::default() };
        let t0 = Instant::now();
        let pairs = candidate_pairs(&baseline, &dp);
        eprintln!("[stage_d] candidate_pairs: {} pairs in {:?}", pairs.len(), t0.elapsed());
        std::io::stderr().flush().ok();

        let mut n_confirmed = 0usize;
        for (k, &(a, b)) in pairs.iter().enumerate() {
            let la = baseline[a].seq.len();
            let lb = baseline[b].seq.len();
            eprintln!(
                "[stage_d] START pair {k}/{} : rep[{a}] len={la} x rep[{b}] len={lb} (max={})",
                pairs.len(), la.max(lb)
            );
            std::io::stderr().flush().ok();
            let t1 = Instant::now();
            let cr = confirm_edge(&baseline[a].seq, &baseline[b].seq, &dp);
            let dt = t1.elapsed();
            if cr.is_some() { n_confirmed += 1; }
            eprintln!("[stage_d] DONE  pair {k}/{} in {:?} -> {:?}", pairs.len(), dt, cr);
            std::io::stderr().flush().ok();
        }
        eprintln!("[stage_d] ALL {} pairs done, {} confirmed, total {:?}", pairs.len(), n_confirmed, t0.elapsed());
    }

    #[test]
    fn genome_reps_finds_family_copies_with_genomic_seq_and_no_introns() {
        // Real subset fixture: 3 near-identical NCF1 copies + 2 unrelated decoys, each as its own contig.
        // genome_reps must surface the duplicated NCF1 copies as reps (with genomic seq, empty introns).
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
        let fa = "tests/fixtures/from_genome/subset.fa";
        if std::fs::metadata(fa).is_err() { eprintln!("fixture absent; skip"); return; }
        // windows = each contig full-length (as in windows.bed).
        let windows: Vec<(String, u64, u64)> = [
            ("NCF1", 15440u64), ("NCF1B", 15319), ("NCF1C", 15406), ("DECOY1", 10001), ("DECOY2", 10001),
        ].iter().map(|(c, l)| (c.to_string(), 0u64, *l)).collect();
        let p = GenomeRepParams { min_identity: 0.90, min_block: 400, ..Default::default() };
        let reps = genome_reps(fa, &windows, &p).unwrap();
        // the three NCF1 copies must all be discovered as duplicated loci.
        for want in ["NCF1", "NCF1B", "NCF1C"] {
            assert!(reps.iter().any(|r| r.chrom == want), "missing duplicated locus {want}; got {:?}",
                    reps.iter().map(|r| r.chrom.as_str()).collect::<Vec<_>>());
        }
        // every rep is genomic: empty intron chain and seq length == span.
        for r in &reps {
            assert!(r.introns.is_empty(), "DNA rep {} must have empty intron chain", r.chrom);
            assert_eq!(r.seq.len() as u64, r.end - r.start, "DNA rep {} seq must be the full genomic span", r.chrom);
        }
    }
}

// ---- merged 2026-10-05: was `vg_family/project_all.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod project_all {
//! Generalized projection (`--project-all-families`): project EVERY resolved copy's consensus (not one per
//! family) to localize members a single family-consensus projection misses. OPTIONAL, additive, DNA-localized
//! parCN leg — never changes the RNA-split catalog or the family definition. Pure extraction + dedup + row
//! formatting here; the minimap2 projection + read-support gate are wired in gw_family_catalog.
//!
//! **STATUS:** OPT-IN — --project-all-families, `#[arg(long, default_value_t = false)]` at src/bin/gw_family_catalog.rs:195-196; equivalently env RUSTLE_PROJECT_ALL_FAMILIES=  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use std::collections::HashMap;

use crate::vg_family::genome_projection::CopyLocus;

fn recip_overlap(a: &CopyLocus, b: &CopyLocus) -> f64 {
    if a.chrom != b.chrom { return 0.0; }
    let (lo, hi) = (a.start.max(b.start), a.end.min(b.end));
    if hi <= lo { return 0.0; }
    let ov = (hi - lo) as f64;
    (ov / (a.end - a.start).max(1) as f64).min(ov / (b.end - b.start).max(1) as f64)
}

/// Collapse reciprocal-overlap ≥ 0.50 loci (from different sibling consensuses hitting one genomic locus)
/// into one, keeping the highest-identity survivor.
pub fn dedup_overlapping(mut loci: Vec<CopyLocus>) -> Vec<CopyLocus> {
    loci.sort_by(|a, b| b.identity.partial_cmp(&a.identity).unwrap_or(std::cmp::Ordering::Equal));
    let mut kept: Vec<CopyLocus> = Vec::new();
    for l in loci {
        if !kept.iter().any(|k| recip_overlap(k, &l) >= 0.50) { kept.push(l); }
    }
    kept
}

pub fn overlaps_any(chrom: &str, start: u64, end: u64, spans: &[(String, u64, u64)]) -> bool {
    spans.iter().any(|(c, s, e)| c == chrom && !(*s > end || *e < start))
}

pub fn format_allproj_row(family_id: &str, l: &CopyLocus, n_support: usize, overlaps_existing: bool) -> String {
    format!("{family_id}\t{}\t{}\t{}\t{:.3}\t{}\t{}", l.chrom, l.start, l.end, l.identity, n_support, overlaps_existing)
}

#[derive(Clone, Debug)]
pub struct CopyIn { pub seq: Vec<u8>, pub chrom: String, pub start: u64, pub end: u64 }

/// One `(family_id, consensus)` entry PER COPY, with the family_id repeated across its copies so
/// `project_families_batch` unions all copies' hits under that one family key.
pub fn all_copy_consensuses(fams: &[(String, Vec<CopyIn>)]) -> Vec<(String, Vec<u8>)> {
    fams.iter().flat_map(|(fid, copies)| copies.iter().map(move |c| (fid.clone(), c.seq.clone()))).collect()
}

/// Per-family copy spans, for the projection's `known` self-exclusion (a copy projecting back onto its own
/// catalogued locus is not a new localization).
pub fn known_from_fams(fams: &[(String, Vec<CopyIn>)]) -> HashMap<String, Vec<(String, u64, u64)>> {
    fams.iter().map(|(fid, copies)| (fid.clone(), copies.iter().map(|c| (c.chrom.clone(), c.start, c.end)).collect())).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dedup_overlap_and_row_format() {
        use crate::vg_family::genome_projection::CopyLocus;
        let mk = |s: u64, e: u64, id: f64| CopyLocus { chrom: "chr7".into(), start: s, end: e, identity: id, cov: 0.95 };
        // two overlapping hits (id .994 vs .982) + one disjoint -> 2 survivors, higher id kept
        let out = dedup_overlapping(vec![mk(75976253,75991692,0.994), mk(75976300,75991600,0.982), mk(76360590,76375995,0.990)]);
        assert_eq!(out.len(), 2);
        assert!(out.iter().any(|l| (l.identity - 0.994).abs() < 1e-9 && l.start == 75976253));
        assert!(overlaps_any("chr7", 75976300, 75991600, &[("chr7".into(), 75976253, 75991692)]));
        assert!(!overlaps_any("chr7", 76360590, 76375995, &[("chr7".into(), 75976253, 75991692)]));
        let row = format_allproj_row("GWFAM7", &mk(75976253,75991692,0.994), 41, false);
        assert_eq!(row, "GWFAM7\tchr7\t75976253\t75991692\t0.994\t41\tfalse");
    }

    #[test]
    fn extraction_repeats_fid_and_collects_spans() {
        let fams = vec![
            ("GWFAM0".to_string(), vec![
                CopyIn { seq: b"ACGT".to_vec(), chrom: "chr1".into(), start: 100, end: 200 },
                CopyIn { seq: b"ACGA".to_vec(), chrom: "chr1".into(), start: 500, end: 600 },
            ]),
            ("GWFAM1".to_string(), vec![
                CopyIn { seq: b"TTTT".to_vec(), chrom: "chr2".into(), start: 10, end: 20 },
            ]),
        ];
        let cons = all_copy_consensuses(&fams);
        // one entry per copy, family_id repeated per copy (so project_families_batch unions per family)
        assert_eq!(cons, vec![
            ("GWFAM0".to_string(), b"ACGT".to_vec()),
            ("GWFAM0".to_string(), b"ACGA".to_vec()),
            ("GWFAM1".to_string(), b"TTTT".to_vec()),
        ]);
        let known = known_from_fams(&fams);
        assert_eq!(known["GWFAM0"], vec![("chr1".to_string(),100,200), ("chr1".to_string(),500,600)]);
        assert_eq!(known["GWFAM1"], vec![("chr2".to_string(),10,20)]);
    }
}
}

// ---- merged 2026-10-05: was `vg_family/seed_projection.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod seed_projection {
//! `--seed`: a QUERY over the emitted catalog. **Not** a term in the family definition.
//!
//! ## Why this is a projection and not a pipeline
//!
//! The seeded probe in `bench/crossspecies/seed_family.sh` builds its node set *from the seed*
//! (`V(s)` = intervals receiving ≥ β aligned bp from `s`). That construction is measurably
//! seed-dependent and it is **not** what the shipped binary does:
//!
//! * strict seed-invariance FAILS on it — seeding independently from each of the 19 annotated human
//!   NPIP genes gives only **4 distinct `F(s)` as sets of loci**, agreeing on **64 of 171** seed
//!   pairs, because `V(s)` inherits the *seed's own length* (seed NPIPB8, 10.6 kb → 27 loci of
//!   ~9.9–10.6 kb; seed NPIPB11, 25.2 kb → loci of 8.8–25.6 kb). Membership at the level of
//!   annotated genes is invariant (19/19 seeds → one gene set); locus EXTENT is not;
//! * on the same probe a gorilla single-copy CONTROL (`SDHA`) returns a **14-locus "family"** at
//!   induced density 0.747 — legal at both γ = 0.20 and γ = 0.40 — where the shipped binary
//!   correctly emits nothing. 52 of its 68 accepted edges have "coverage" > 1.0 (up to 2.019),
//!   because `(qe-qs)/min(|u|,|v|)` is unbounded above;
//! * the seeded probe and the shipped catalog disagree wholesale on gorilla (MAGEA4 seed → 1 locus
//!   vs 9 shipped copies; GSTM3 seed, 4,057 bp < β → 0 loci vs 4 shipped copies; HERC2 → 164 vs 2).
//!
//! The shipped binary's node set is built from the reads (or, under `--from-genome`, from genome
//! self-alignment) and never consults a seed. `gamma_quasi_clique_partition` starts from
//! `all_components` INCLUDING singletons, so the blocks **partition** the whole node set. Over a
//! seed-free node set, "the block containing `s`" is therefore a fact about a partition, not a
//! parameterised computation: every node lands in exactly one block, and any two nodes in the same
//! block return the same block. That is what makes `--seed` a legitimate query — and it is why the
//! query must read the emitted catalog rather than re-run anything with `s` in scope.
//!
//! ## Scope limit, stated up front
//!
//! The projection can only see blocks the run actually EMITTED. At the default `--min-copies 2` a
//! singleton block is not written, so a seed landing on one is reported as `ABSTAIN_NO_OVERLAP` and
//! is indistinguishable here from a seed landing on no transcribed locus at all. Re-run with
//! `--min-copies 1` to separate the two cases. Abstention is the O1-safe answer: "no emitted family
//! at this seed" is a legitimate result and is never guessed around.
//!
//! **STATUS:** OPT-IN — --seed (src/bin/gw_family_catalog.rs:223-224, `#[arg(long)] seed: Vec<String>`, default = empty vec; project_seeds returns Ok(()) immediately on an em  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use anyhow::{anyhow, Result};

/// A seed interval, stored 0-based half-open to match `copies.tsv`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    /// The spec exactly as the user typed it, echoed into the output so a row is traceable.
    pub spec: String,
}

/// Outcome of projecting one seed onto the emitted catalog.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeedHit {
    /// The seed overlaps `copy_idx` of family `family_idx` by `overlap_bp` bases.
    Hit {
        family_idx: usize,
        copy_idx: usize,
        overlap_bp: u64,
    },
    /// No emitted copy overlaps the seed. Honest answer; not an error.
    Abstain,
}

/// Parse `chrom:start-end`.
///
/// Coordinates are **1-based inclusive** (the GFF/samtools region convention, which is what a seed
/// pasted out of an annotation carries) and are converted to 0-based half-open on the way in, so
/// they line up with `copies.tsv`. `,` and `_` are accepted as digit separators. The chromosome
/// name may itself contain `:` (none of ours do, but RefSeq-style names contain `.` and `_`), so the
/// split is on the LAST `:`.
pub fn parse_seed(spec: &str) -> Result<SeedLocus> {
    let s = spec.trim();
    let (chrom, range) = s
        .rsplit_once(':')
        .ok_or_else(|| anyhow!("--seed {spec:?}: expected chrom:start-end (1-based inclusive)"))?;
    if chrom.is_empty() {
        return Err(anyhow!("--seed {spec:?}: empty chromosome"));
    }
    let (a, b) = range
        .rsplit_once('-')
        .ok_or_else(|| anyhow!("--seed {spec:?}: expected chrom:start-end (1-based inclusive)"))?;
    let clean = |x: &str| x.replace([',', '_'], "");
    let a: u64 = clean(a)
        .parse()
        .map_err(|_| anyhow!("--seed {spec:?}: start is not an integer"))?;
    let b: u64 = clean(b)
        .parse()
        .map_err(|_| anyhow!("--seed {spec:?}: end is not an integer"))?;
    if a == 0 {
        return Err(anyhow!("--seed {spec:?}: coordinates are 1-based, so start must be >= 1"));
    }
    if b < a {
        return Err(anyhow!("--seed {spec:?}: end < start"));
    }
    Ok(SeedLocus {
        chrom: chrom.to_string(),
        start: a - 1, // 1-based inclusive -> 0-based half-open
        end: b,
        spec: s.to_string(),
    })
}

/// Overlap of two half-open intervals, in bp (0 when disjoint or touching).
pub fn overlap_bp(a0: u64, a1: u64, b0: u64, b1: u64) -> u64 {
    let lo = a0.max(b0);
    let hi = a1.min(b1);
    hi.saturating_sub(lo)
}

/// Project a seed onto the emitted catalog.
///
/// `fams[fi][ci]` must be `(chrom, start, end)` in **exactly the order `emit_catalog` writes them** —
/// families in emitted order (so `fi` ↔ `GWFAM{fi}`), copies sorted by `(chrom, start)` (so `ci` ↔
/// `copy_idx`) — otherwise the reported ids do not name the rows the user can look up.
///
/// Rule: maximum overlapping bases wins. This is deliberate and sufficient rather than clever — on
/// the gorilla RABL2 check the shipped copy boundaries sat within 3–46 bp of the independently
/// seeded interval, so no seed is anywhere near two copies at once. Ties (including the degenerate
/// zero-length seed) resolve to the lowest `(family_idx, copy_idx)`, which is deterministic because
/// the caller's order is deterministic.
pub fn project_seed(seed: &SeedLocus, fams: &[Vec<(String, u64, u64)>]) -> SeedHit {
    let mut best: Option<(u64, usize, usize)> = None;
    for (fi, copies) in fams.iter().enumerate() {
        for (ci, (chrom, start, end)) in copies.iter().enumerate() {
            if chrom != &seed.chrom {
                continue;
            }
            let ov = overlap_bp(seed.start, seed.end, *start, *end);
            if ov == 0 {
                continue;
            }
            if best.map_or(true, |(b, _, _)| ov > b) {
                best = Some((ov, fi, ci));
            }
        }
    }
    match best {
        Some((overlap_bp, family_idx, copy_idx)) => SeedHit::Hit {
            family_idx,
            copy_idx,
            overlap_bp,
        },
        None => SeedHit::Abstain,
    }
}

/// Header of `<out>.seed.tsv`.
pub const SEED_TSV_HEADER: &str =
    "seed\tstatus\tfamily_id\tn_copies\tmember_idx\tchrom\tstart\tend\tis_seed_locus\tseed_overlap_bp";

/// Render the projection of one seed as the rows of `<out>.seed.tsv`: one row per MEMBER of the
/// component containing the seed (that is the object the query returns — the component, not the
/// single locus the seed happened to land on), or exactly one `ABSTAIN_NO_OVERLAP` row.
pub fn format_seed_rows(seed: &SeedLocus, fams: &[Vec<(String, u64, u64)>], hit: &SeedHit) -> Vec<String> {
    match hit {
        SeedHit::Abstain => vec![format!(
            "{}\tABSTAIN_NO_OVERLAP\t.\t0\t.\t.\t.\t.\t.\t0",
            seed.spec
        )],
        SeedHit::Hit {
            family_idx,
            copy_idx,
            overlap_bp,
        } => {
            let copies = &fams[*family_idx];
            copies
                .iter()
                .enumerate()
                .map(|(ci, (chrom, start, end))| {
                    format!(
                        "{}\tHIT\tGWFAM{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        seed.spec,
                        family_idx,
                        copies.len(),
                        ci,
                        chrom,
                        start,
                        end,
                        if ci == *copy_idx { "true" } else { "false" },
                        if ci == *copy_idx { *overlap_bp } else { 0 },
                    )
                })
                .collect()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fams() -> Vec<Vec<(String, u64, u64)>> {
        vec![
            // GWFAM0: two copies on chr1
            vec![
                ("chr1".to_string(), 1000, 2000),
                ("chr1".to_string(), 5000, 6000),
            ],
            // GWFAM1: two copies, one of them adjacent-but-disjoint to GWFAM0's second copy
            vec![
                ("chr1".to_string(), 6000, 7000),
                ("chr2".to_string(), 100, 900),
            ],
        ]
    }

    // ---- parse_seed ----

    #[test]
    fn parse_converts_one_based_inclusive_to_half_open() {
        let s = parse_seed("chr1:1001-2000").unwrap();
        assert_eq!(s.chrom, "chr1");
        assert_eq!((s.start, s.end), (1000, 2000));
        assert_eq!(s.spec, "chr1:1001-2000");
    }

    #[test]
    fn parse_accepts_refseq_style_names_and_digit_separators() {
        let s = parse_seed("NC_073236.2:140,592,916-140_608_794").unwrap();
        assert_eq!(s.chrom, "NC_073236.2");
        assert_eq!((s.start, s.end), (140_592_915, 140_608_794));
    }

    #[test]
    fn parse_rejects_zero_start_because_coordinates_are_one_based() {
        assert!(parse_seed("chr1:0-100").is_err());
    }

    #[test]
    fn parse_rejects_malformed_specs() {
        assert!(parse_seed("chr1").is_err());
        assert!(parse_seed("chr1:100").is_err());
        assert!(parse_seed("chr1:abc-100").is_err());
        assert!(parse_seed("chr1:200-100").is_err());
        assert!(parse_seed(":100-200").is_err());
    }

    #[test]
    fn single_base_seed_is_legal_and_has_length_one() {
        let s = parse_seed("chr1:1500-1500").unwrap();
        assert_eq!((s.start, s.end), (1499, 1500));
    }

    // ---- overlap ----

    #[test]
    fn touching_intervals_do_not_overlap() {
        assert_eq!(overlap_bp(0, 100, 100, 200), 0);
        assert_eq!(overlap_bp(0, 101, 100, 200), 1);
    }

    // ---- project_seed ----

    #[test]
    fn seed_inside_a_copy_returns_that_copy() {
        let f = fams();
        let s = parse_seed("chr1:1201-1300").unwrap();
        assert_eq!(
            project_seed(&s, &f),
            SeedHit::Hit { family_idx: 0, copy_idx: 0, overlap_bp: 100 }
        );
    }

    #[test]
    fn seed_off_every_copy_abstains_rather_than_guessing() {
        let f = fams();
        let s = parse_seed("chr1:3001-4000").unwrap();
        assert_eq!(project_seed(&s, &f), SeedHit::Abstain);
    }

    #[test]
    fn seed_on_an_absent_chromosome_abstains() {
        let f = fams();
        let s = parse_seed("chrX:1001-2000").unwrap();
        assert_eq!(project_seed(&s, &f), SeedHit::Abstain);
    }

    #[test]
    fn straddling_seed_goes_to_the_copy_it_overlaps_most_not_the_first_one() {
        // 5900..6600 overlaps GWFAM0 copy1 (5000..6000) by 100 bp and GWFAM1 copy0 (6000..7000) by
        // 600 bp. Max-overlap must win over first-seen, or the answer would depend on emit order.
        let f = fams();
        let s = parse_seed("chr1:5901-6600").unwrap();
        assert_eq!(
            project_seed(&s, &f),
            SeedHit::Hit { family_idx: 1, copy_idx: 0, overlap_bp: 600 }
        );
    }

    #[test]
    fn exact_ties_resolve_to_the_lowest_family_then_copy_index() {
        let f = vec![
            vec![("chr1".to_string(), 0, 100)],
            vec![("chr1".to_string(), 0, 100)],
        ];
        let s = parse_seed("chr1:1-100").unwrap();
        assert_eq!(
            project_seed(&s, &f),
            SeedHit::Hit { family_idx: 0, copy_idx: 0, overlap_bp: 100 }
        );
    }

    #[test]
    fn empty_catalog_abstains() {
        let s = parse_seed("chr1:1-100").unwrap();
        assert_eq!(project_seed(&s, &[]), SeedHit::Abstain);
    }

    // ---- the returned object is the COMPONENT, not the locus ----

    #[test]
    fn hit_rows_carry_every_member_of_the_component() {
        let f = fams();
        let s = parse_seed("chr1:1201-1300").unwrap();
        let hit = project_seed(&s, &f);
        let rows = format_seed_rows(&s, &f, &hit);
        assert_eq!(rows.len(), 2, "the query returns the whole component: {rows:?}");
        assert!(rows.iter().all(|r| r.contains("\tGWFAM0\t")));
        assert_eq!(rows.iter().filter(|r| r.ends_with("\ttrue\t100")).count(), 1);
        assert!(rows[1].contains("\tfalse\t0"), "non-seed members carry overlap 0: {}", rows[1]);
    }

    #[test]
    fn abstain_emits_exactly_one_row_and_names_the_reason() {
        let f = fams();
        let s = parse_seed("chr9:1-100").unwrap();
        let rows = format_seed_rows(&s, &f, &project_seed(&s, &f));
        assert_eq!(rows.len(), 1);
        assert!(rows[0].contains("ABSTAIN_NO_OVERLAP"), "{}", rows[0]);
        assert!(rows[0].starts_with("chr9:1-100\t"));
    }

    #[test]
    fn every_row_has_the_same_field_count_as_the_header() {
        let f = fams();
        let ncol = SEED_TSV_HEADER.split('\t').count();
        for spec in ["chr1:1201-1300", "chr9:1-100"] {
            let s = parse_seed(spec).unwrap();
            for r in format_seed_rows(&s, &f, &project_seed(&s, &f)) {
                assert_eq!(r.split('\t').count(), ncol, "{r}");
            }
        }
    }

    // ---- the property the projection inherits from the partition ----

    #[test]
    fn any_member_of_a_component_returns_the_same_component() {
        // P1 is a THEOREM over a seed-free node set: blocks partition the nodes, so seeding from any
        // member returns the identical block. Asserted here so a future change that makes the answer
        // seed-dependent fails a test rather than a re-measurement.
        let f = fams();
        let seeds = ["chr1:1001-2000", "chr1:5001-6000"];
        let rows: Vec<Vec<String>> = seeds
            .iter()
            .map(|spec| {
                let s = parse_seed(spec).unwrap();
                // strip the echoed seed spec + per-seed overlap columns; compare the COMPONENT.
                format_seed_rows(&s, &f, &project_seed(&s, &f))
                    .iter()
                    .map(|r| {
                        let c: Vec<&str> = r.split('\t').collect();
                        c[2..8].join("\t")
                    })
                    .collect()
            })
            .collect();
        assert_eq!(rows[0], rows[1]);
    }
}
}
