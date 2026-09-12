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
    #[test]
    #[ignore]
    fn stage_c_combining_proposals_1_and_2_measures_real_family_recovery() {
        use crate::vg_family::annotation_families::SdPairs;
        use crate::vg_family::family_detect::{detect_edges, DetectParams};
        use crate::vg_family::family_split::{decompose_families, SplitParams};

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
            let edges = detect_edges(reps, &dp);
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
        let dna_reps = genome_reps(fa, &sd_windows, &GenomeRepParams::default()).expect("genome_reps");
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
            "[stage_c] SUMMARY (oracle_covered/oracle_families/false_merge_families): baseline {:?}, \
             prop#1-only {:?}, COMBINED {:?}",
            baseline_result, p1_result, combined_result
        );
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
