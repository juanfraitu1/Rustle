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
    /// Body factored into [`stage_c_or_e_body`] so the historical unbudgeted reproduction (this test,
    /// kept for the record -- it WILL hang if run bare, that is the documented §6j6 finding) and the
    /// budgeted re-attempt (`stage_e_...`, §6j8) share one real implementation instead of two copies
    /// that could silently drift apart.
    #[test]
    #[ignore]
    fn stage_c_combining_proposals_1_and_2_measures_real_family_recovery() {
        stage_c_or_e_body(None, "UNBUDGETED (historical §6j6 hang reproduction -- do not run without `timeout`)");
    }

    /// §6j8: the SAME real #1+#2 integration as `stage_c_...` above, with `confirm_edge`'s new per-pair
    /// time budget enabled (500ms -- justified in docs/o1_ledger.md §6j8 from the real per-pair
    /// distribution measured in §6j7: 192 real pairs at 2.51s-19.75s each, so 500ms is comfortably below
    /// the observed MINIMUM -- a 5x margin, not a hair's-width cut -- forcing the faithful fallback for
    /// every pair in this documented hard neighborhood while remaining orders of magnitude above what any
    /// ordinary/fast pair elsewhere needs. Live wall-clock runs at 2s and 1s budgets on this same
    /// substrate (§6j8) measured total time scaling roughly with the budget (164s at 2s, 95s at 1s, for
    /// the baseline cell alone) -- tightening further here trades zero additional fidelity risk (every
    /// pair in this neighborhood already exceeds even the smallest of these budgets by a wide margin, so
    /// which ones fall back does not change) for wall-clock headroom to let all three cells complete.
    /// This is the measurement §6j6 was blocked from getting.
    /// `#[ignore]`d: needs real gorilla data; run under a shell `timeout`, never bare (should now finish
    /// well inside it).
    #[test]
    #[ignore]
    fn stage_e_combining_proposals_with_time_budget_measures_real_family_recovery() {
        stage_c_or_e_body(
            Some(std::time::Duration::from_millis(500)),
            "BUDGETED 500ms (§6j8 fix -- the measurement §6j6 was blocked from getting)",
        );
    }

    fn stage_c_or_e_body(time_budget: Option<std::time::Duration>, label_suffix: &str) {
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
            let dp = DetectParams { time_budget, ..DetectParams::default() };
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

    /// DIAGNOSTIC (false-merge evidence dump for §6j8's "false merge" families). Additive only: runs
    /// ONLY the `baseline` and `proposal #1` cells of `stage_c_or_e_body` (same loading, same bounding,
    /// same proposal-#1 correction, same `detect_edges` + `decompose_families` scoring, same order --
    /// `genome_reps` is never called) and writes a self-describing evidence dump per cell under
    /// `$RUSTLE_FM_OUT/<cell>/` (reps.tsv, reps.fa, candidates.tsv, edges.tsv, families.tsv, summary.tsv).
    ///
    /// Env: `RUSTLE_FM_PHASE` = `dump` (default) | `bridge`; `RUSTLE_FM_OUT` (default
    /// /mnt/linuxdisk/home/juanfraitu/o1_falsemerge); `RUSTLE_FM_BUDGET_MS` (default 500);
    /// `RUSTLE_FM_REAL` (default 1: run the real `detect_edges`); `RUSTLE_FM_REPLICA` (default 1: also run
    /// an instrumented replica of `confirm_edge` that records, per orientation, whether the budget fired).
    /// `bridge` phase: reads `$RUSTLE_FM_PAIRS` (cell\ti\tj) + `$RUSTLE_FM_OUT/<cell>/reps.fa` and computes
    /// the PRODUCTION (budget None) `confirm_edge` value one pair at a time into `$RUSTLE_FM_BRIDGE_OUT`.
    /// `#[ignore]`d: needs real gorilla data; run under a shell `timeout`, never bare.
    /// ⚠ `RUSTLE_FM_REAL=1`/`RUSTLE_FM_REPLICA=1` use the time budget, which spawns one uncancellable poasta
    /// thread per timed-out pair; with both on this test was OOM-killed at 25.7 GB (docs/o1_ledger.md §6j9).
    /// Leave them off (the default): edges then come from a serial LCS emulation in ~2 s / ~0.6 GB.
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
            other => panic!("unknown RUSTLE_FM_PHASE {other}"),
        }
    }

    fn fm_env_flag(name: &str, default: bool) -> bool {
        match std::env::var(name) {
            Ok(v) => v != "0",
            Err(_) => default,
        }
    }

    /// Per-orientation instrumented replica of `family_graph::contiguous_core_coverage_bounded_budgeted`
    /// (identical control flow) that also reports which branch produced the value and whether the memo
    /// already held the exact key when the call started.
    fn fm_budgeted_instrumented(
        a: &[u8],
        b: &[u8],
        cap: usize,
        astar: bool,
        budget: Option<std::time::Duration>,
    ) -> (f64, &'static str, f64, bool) {
        use crate::vg_family::family_graph::{
            contiguous_core_coverage_bounded_with, core_memo_contains, longest_common_substring,
        };
        let t = std::time::Instant::now();
        let over = a.len().max(b.len()) > cap;
        let memo_had = core_memo_contains(a, b, cap, astar);
        let Some(d) = budget else {
            let v = contiguous_core_coverage_bounded_with(a, b, cap, astar);
            return (v, if over { "lcs_len_cap" } else { "poasta_exact" }, t.elapsed().as_secs_f64(), memo_had);
        };
        let (tx, rx) = std::sync::mpsc::channel();
        let (ao, bo) = (a.to_vec(), b.to_vec());
        std::thread::spawn(move || {
            let v = contiguous_core_coverage_bounded_with(&ao, &bo, cap, astar);
            let _ = tx.send(v);
        });
        match rx.recv_timeout(d) {
            Ok(v) => (v, if over { "lcs_len_cap" } else { "poasta_within_budget" }, t.elapsed().as_secs_f64(), memo_had),
            Err(_) => {
                let minlen = a.len().min(b.len());
                let v = if minlen == 0 { 0.0 } else { longest_common_substring(a, b) as f64 / minlen as f64 };
                (v, if over { "lcs_len_cap_timeout" } else { "lcs_budget_timeout" }, t.elapsed().as_secs_f64(), memo_had)
            }
        }
    }

    struct FmReplica {
        fwd: (f64, &'static str, f64, bool),
        rc: Option<(f64, &'static str, f64, bool)>,
        cr: f64,
    }

    fn fm_replica_confirm(
        a: &[u8],
        b: &[u8],
        dp: &crate::vg_family::family_detect::DetectParams,
    ) -> FmReplica {
        use crate::vg_family::family_graph::{upper_cow, EDGE_CONFIRM_ASTAR};
        use crate::vg_family::seq_utils::reverse_complement;
        let au = upper_cow(a);
        let bu = upper_cow(b);
        let fwd = fm_budgeted_instrumented(&au, &bu, dp.len_cap, EDGE_CONFIRM_ASTAR, dp.time_budget);
        let mut cr = fwd.0;
        let mut rc = None;
        if cr < dp.t_core {
            let r = fm_budgeted_instrumented(&au, &reverse_complement(&bu), dp.len_cap, EDGE_CONFIRM_ASTAR, dp.time_budget);
            if r.0 > cr { cr = r.0; }
            rc = Some(r);
        }
        FmReplica { fwd, rc, cr }
    }

    fn fm_dump(root: &str) {
        use crate::vg_family::annotation_families::SdPairs;
        use crate::vg_family::family_detect::{detect_edges, DetectParams};
        use crate::vg_family::family_graph::core_memo_stats;

        let fa = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta";
        let copies_tsv = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.tsv";
        let copies_fa = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.fa";
        let sedef_bed = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed";
        let sd_windows_bed = "/mnt/linuxdisk/home/juanfraitu/o1_fromgenome_sd/npip_seeded_windows2.bed";
        let oracle_path = "/mnt/linuxdisk/home/juanfraitu/o1_oracle/npip31.regions";
        for p in [fa, copies_tsv, copies_fa, sedef_bed, sd_windows_bed, oracle_path] {
            if std::fs::metadata(p).is_err() { eprintln!("required real-data file {p} absent; skip"); return; }
        }
        let budget_ms: u64 = std::env::var("RUSTLE_FM_BUDGET_MS").ok().and_then(|v| v.parse().ok()).unwrap_or(500);
        // Both default OFF: with a time budget they spawn one uncancellable poasta thread per timed-out pair,
        // which OOM-killed this test at 25.7 GB. Off => edges come from a serial LCS emulation, no threads.
        let run_real = fm_env_flag("RUSTLE_FM_REAL", false);
        let run_replica = fm_env_flag("RUSTLE_FM_REPLICA", false);
        let memo_env = std::env::var("RUSTLE_POA_MEMO").unwrap_or_else(|_| "<unset>".into());
        eprintln!("[stage_f] root={root} budget_ms={budget_ms} real={run_real} replica={run_replica} RUSTLE_POA_MEMO={memo_env}");

        // --- loading: identical to stage_c_or_e_body (exon strings kept alongside for the dump) ---
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

        let time_budget = (run_real || run_replica).then(|| std::time::Duration::from_millis(budget_ms));
        let dp = DetectParams { time_budget, ..DetectParams::default() };
        let memo_start = core_memo_stats();

        // --- BASELINE real detect_edges (same position in the order as stage_c_or_e_body) ---
        let (base_edges, base_secs) = if run_real {
            let t0 = std::time::Instant::now();
            let e = detect_edges(&baseline, &dp);
            let s = t0.elapsed().as_secs_f64();
            eprintln!("[stage_f:baseline] real detect_edges: {} edges in {s:.1}s", e.len());
            (Some(e), s)
        } else { (None, 0.0) };
        let memo_after_base = core_memo_stats();

        // --- PROPOSAL #1: identical to stage_c_or_e_body ---
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
        let (p1_edges, p1_secs) = if run_real {
            let t0 = std::time::Instant::now();
            let e = detect_edges(&corrected, &dp);
            let s = t0.elapsed().as_secs_f64();
            eprintln!("[stage_f:proposal1] real detect_edges: {} edges in {s:.1}s", e.len());
            (Some(e), s)
        } else { (None, 0.0) };
        let memo_after_p1 = core_memo_stats();

        let memo_info = format!(
            "memo_stats_start\t{:?}\nmemo_stats_after_baseline_real\t{:?}\nmemo_stats_after_proposal1_real\t{:?}\n",
            memo_start, memo_after_base, memo_after_p1
        );
        let base_kinds: Vec<&'static str> = n_exons.iter().map(|&n| if n > 1 { "rna_spliced" } else { "rna_single_exon" }).collect();
        let p1_kinds: Vec<&'static str> = (0..corrected.len())
            .map(|i| if is_corrected[i] { "genomic_span_p1corrected" } else { base_kinds[i] }).collect();
        fm_write_cell(root, "baseline", &baseline, &baseline, &exon_strs, &base_kinds, &is_corrected.iter().map(|_| false).collect::<Vec<_>>(),
            &oracle, base_edges, base_secs, run_replica, &dp, budget_ms, &memo_env, &memo_info);
        fm_write_cell(root, "proposal1", &corrected, &baseline, &exon_strs, &p1_kinds, &is_corrected,
            &oracle, p1_edges, p1_secs, run_replica, &dp, budget_ms, &memo_env, &memo_info);
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
        real_edges: Option<Vec<(usize, usize, f64)>>,
        real_secs: f64,
        run_replica: bool,
        dp: &crate::vg_family::family_detect::DetectParams,
        budget_ms: u64,
        memo_env: &str,
        memo_info: &str,
    ) {
        use crate::vg_family::family_detect::candidate_pairs;
        use crate::vg_family::family_graph::{core_memo_stats, longest_common_substring, upper_cow};
        use crate::vg_family::family_split::{connected_components, decompose_families, SplitParams};
        use crate::vg_family::seq_utils::reverse_complement;
        use rayon::prelude::*;
        use std::fmt::Write as _;

        let dir = format!("{root}/{cell}");
        std::fs::create_dir_all(&dir).unwrap();
        let n = reps.len();
        let pairs = candidate_pairs(reps, dp);

        // replica (instrumented confirm_edge) over every candidate pair, rayon-parallel like detect_edges.
        let t_rep = std::time::Instant::now();
        let replica: Option<Vec<FmReplica>> = if run_replica {
            Some(pairs.par_iter().map(|&(a, b)| fm_replica_confirm(&reps[a].seq, &reps[b].seq, dp)).collect())
        } else { None };
        let replica_secs = t_rep.elapsed().as_secs_f64();
        let memo_after_replica = core_memo_stats();

        // exact-substring (LCS) values for both orientations of every candidate pair (cheap, deterministic).
        let lcs: Vec<(usize, usize)> = pairs.iter().map(|&(a, b)| {
            let au = upper_cow(&reps[a].seq);
            let bu = upper_cow(&reps[b].seq);
            (longest_common_substring(&au, &bu), longest_common_substring(&au, &reverse_complement(&bu)))
        }).collect();

        // LCS emulation of confirm_edge: forward, then reverse complement only if forward misses T_CORE.
        let lcs_emulated: Vec<(usize, usize, f64)> = pairs.iter().zip(lcs.iter()).filter_map(|(&(a, b), &(lf, lr))| {
            let mn = reps[a].seq.len().min(reps[b].seq.len());
            if mn == 0 { return None; }
            let fwd = lf as f64 / mn as f64;
            let cr = if fwd < dp.t_core { fwd.max(lr as f64 / mn as f64) } else { fwd };
            (cr >= dp.t_core).then_some((a, b, cr))
        }).collect();
        let (edges_source, edges): (&str, Vec<(usize, usize, f64)>) = match (real_edges, replica.as_ref()) {
            (Some(e), _) => ("real_detect_edges", e),
            (None, Some(rv)) => ("replica", pairs.iter().zip(rv.iter())
                .filter(|(_, r)| r.cr >= dp.t_core).map(|(&(a, b), r)| (a, b, r.cr)).collect()),
            (None, None) => ("lcs_emulation", lcs_emulated),
        };
        let edge_val: HashMap<(usize, usize), f64> = edges.iter().map(|&(a, b, v)| ((a, b), v)).collect();
        let replica_edges = replica.as_ref().map(|rv| rv.iter().filter(|r| r.cr >= dp.t_core).count());

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
        let header = "i\tj\ttid_i\ttid_j\tkind_i\tkind_j\tlen_i\tlen_j\tmin_len\tmax_len\tover_len_cap\tconfirmed\tcore_frac_run\tlcs_fwd_bp\tlcs_rc_bp\tlcs_orientation_emulated\tlcs_bp_emulated\tlcs_over_min_len\tlcs_over_max_len\tlcs_emulated_core_frac\tlcs_emulated_pass\trun_value_equals_lcs_emulation\treplica_fwd_path\treplica_fwd_value\treplica_fwd_secs\treplica_fwd_memo_had_key\treplica_rc_path\treplica_rc_value\treplica_rc_secs\treplica_rc_memo_had_key\treplica_core_frac\treplica_pass\tfamily_i\tfamily_j\tsame_family\tcomponent_i\tcomponent_j\n";
        let mut ct = String::from(header);
        let mut et = String::from(header);
        let mut n_eq_lcs = 0usize;
        let mut path_counts: std::collections::BTreeMap<String, usize> = std::collections::BTreeMap::new();
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
            let (rfp, rfv, rfs, rfm, rrp, rrv, rrs, rrm, rcr, rpass) = match replica.as_ref() {
                Some(rv) => {
                    let r = &rv[k];
                    *path_counts.entry(format!("fwd:{}", r.fwd.1)).or_insert(0) += 1;
                    if let Some(rc) = r.rc { *path_counts.entry(format!("rc:{}", rc.1)).or_insert(0) += 1; }
                    (r.fwd.1.to_string(), format!("{:.6}", r.fwd.0), format!("{:.3}", r.fwd.2), (r.fwd.3 as u8).to_string(),
                     r.rc.map(|x| x.1.to_string()).unwrap_or_else(|| "not_run".into()),
                     r.rc.map(|x| format!("{:.6}", x.0)).unwrap_or_default(),
                     r.rc.map(|x| format!("{:.3}", x.2)).unwrap_or_default(),
                     r.rc.map(|x| (x.3 as u8).to_string()).unwrap_or_default(),
                     format!("{:.6}", r.cr), ((r.cr >= dp.t_core) as u8).to_string())
                }
                None => ("not_run".into(), String::new(), String::new(), String::new(), "not_run".into(), String::new(), String::new(), String::new(), String::new(), String::new()),
            };
            let (fa_, fb_) = (rep_family[a], rep_family[b]);
            let row = format!("{a}\t{b}\t{}\t{}\t{}\t{}\t{la}\t{lb}\t{mn}\t{mx}\t{}\t{}\t{}\t{lf}\t{lr}\t{orient}\t{lbp}\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{rfp}\t{rfv}\t{rfs}\t{rfm}\t{rrp}\t{rrv}\t{rrs}\t{rrm}\t{rcr}\t{rpass}\t{}\t{}\t{}\t{}\t{}\n",
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
        writeln!(st, "budget_ms\t{budget_ms}").unwrap();
        writeln!(st, "RUSTLE_POA_MEMO\t{memo_env}").unwrap();
        writeln!(st, "edges_source\t{edges_source}").unwrap();
        writeln!(st, "n_reps\t{n}").unwrap();
        writeln!(st, "n_candidate_pairs\t{}", pairs.len()).unwrap();
        writeln!(st, "n_candidate_pairs_over_len_cap\t{}", pairs.iter().filter(|&&(a, b)| reps[a].seq.len().max(reps[b].seq.len()) > dp.len_cap).count()).unwrap();
        writeln!(st, "n_edges\t{}", edges.len()).unwrap();
        writeln!(st, "n_families\t{}", families.len()).unwrap();
        writeln!(st, "oracle_covered\t{oracle_covered}/31").unwrap();
        writeln!(st, "oracle_families\t{}", oracle_families.len()).unwrap();
        writeln!(st, "false_merge_families\t{false_merge_families}/{}", oracle_families.len()).unwrap();
        writeln!(st, "real_detect_edges_secs\t{real_secs:.1}").unwrap();
        writeln!(st, "n_edges_value_equal_lcs_emulation\t{n_eq_lcs}").unwrap();
        writeln!(st, "n_candidates_lcs_emulated_pass\t{}", pairs.iter().enumerate().filter(|(k, &(a, b))| {
            let mn = reps[a].seq.len().min(reps[b].seq.len()) as f64;
            let (lf, lr) = lcs[*k];
            (lf.max(lr) as f64 / mn) >= dp.t_core
        }).count()).unwrap();
        writeln!(st, "replica_run\t{run_replica}").unwrap();
        writeln!(st, "replica_secs\t{replica_secs:.1}").unwrap();
        writeln!(st, "replica_n_edges\t{}", replica_edges.map(|x| x.to_string()).unwrap_or_else(|| "not_run".into())).unwrap();
        for (k, v) in &path_counts { writeln!(st, "replica_path_count_{k}\t{v}").unwrap(); }
        st.push_str(memo_info);
        writeln!(st, "memo_stats_after_{cell}_replica\t{:?}", memo_after_replica).unwrap();
        std::fs::write(format!("{dir}/summary.tsv"), &st).unwrap();
        eprintln!("[stage_f:{cell}] {} reps, {} candidates, {} edges ({edges_source}), {} families; oracle {}/31 in {} families; false-merge {}/{}; replica edges {:?}",
            n, pairs.len(), edges.len(), families.len(), oracle_covered, oracle_families.len(), false_merge_families, oracle_families.len(), replica_edges);
    }

    /// Pairs phase: production `candidate_pairs` over an arbitrary rep set (copies.tsv columns + FASTA in the
    /// same row order, `$RUSTLE_FM_REPS_TSV` / `$RUSTLE_FM_REPS_FA`), with serial LCS for both orientations.
    /// Writes `<root>/<$RUSTLE_FM_CELL>/{reps.fa,candidates.tsv}` so the `bridge` phase can run on it. No threads.
    fn fm_pairs(root: &str) {
        use crate::vg_family::family_detect::{candidate_pairs, DetectParams};
        use crate::vg_family::family_graph::{longest_common_substring, upper_cow};
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
        assert!(dp.time_budget.is_none());
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

    /// Bridge phase: PRODUCTION-default (`time_budget: None`) edge values for a list of pairs, serial.
    fn fm_bridge(root: &str) {
        use crate::vg_family::family_detect::{confirm_edge, DetectParams, LEN_CAP, T_CORE};
        use crate::vg_family::family_graph::{contiguous_core_coverage_bounded_with, longest_common_substring, upper_cow, EDGE_CONFIRM_ASTAR};
        use crate::vg_family::seq_utils::reverse_complement;
        use std::io::Write;
        let pairs_path = std::env::var("RUSTLE_FM_PAIRS").unwrap_or_else(|_| format!("{root}/bridging_pairs.list"));
        let out_path = std::env::var("RUSTLE_FM_BRIDGE_OUT").unwrap_or_else(|_| format!("{root}/bridging_production.tsv"));
        let mut cache: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
        let mut out = std::fs::File::create(&out_path).unwrap();
        writeln!(out, "cell\ti\tj\tlen_i\tlen_j\tmin_len\tmax_len\tproduction_path\tfwd_value\tfwd_secs\trc_value\trc_secs\tproduction_core_frac\tproduction_core_bp\tcore_over_max_len\tpasses_tcore\tconfirm_edge_agrees\tlcs_fwd_bp\tlcs_rc_bp").unwrap();
        let dp = DetectParams::default();
        assert!(dp.time_budget.is_none());
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

        let dp = DetectParams::default();
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
