//! Bulk-merged family arms: genome projection, copy graph, split, utilities, etc.
//!
//! **STATUS:** INFRASTRUCTURE  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

pub mod from_genome {
    //! DNA front-end: discover duplicated genomic loci by self-alignment and emit them as reps for the
    //! shared homology-grouping core (`denovo_pipeline::homology_blocks`). Read-free and annotation-free —
    //! the genome-only counterpart of the RNA read front-end. The reps differ from RNA reps in exactly one
    //! way: genomic `seq` (introns included) and an EMPTY intron chain. That single difference is the
    //! scientific claim (splicing discards the intron/flank sequence that separates near-identical copies).
    //!
    //! **STATUS:** OPT-IN — `--from-genome <BED>` (src/bin/gw_family_catalog.rs:38-39, `#[arg(long)] from_genome: Option<String>`, default None)  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)
    use crate::family::family_detect::DenovoTranscript;
    use crate::family::genome_projection::project_families_batch;
    use crate::genome::GenomeIndex;
    use anyhow::Result;
    use std::collections::{HashMap, HashSet};

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
            Self {
                min_identity: 0.90,
                min_block: 1000,
                max_locus_span: 3_000_000,
                minimap2: "minimap2".into(),
                threads: 4,
            }
        }
    }
    impl GenomeRepParams {
        pub fn from_env() -> Self {
            let mut p = Self::default();
            if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_IDENTITY") {
                if let Ok(x) = v.parse() {
                    p.min_identity = x;
                }
            }
            if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_BLOCK") {
                if let Ok(x) = v.parse() {
                    p.min_block = x;
                }
            }
            if let Ok(v) = std::env::var("RUSTLE_MINIMAP2") {
                p.minimap2 = v;
            }
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
        let contigs: crate::types::DetHashSet<String> =
            windows.iter().map(|(c, _, _)| c.clone()).collect();
        let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;

        // (1) SD detector: each window's sequence is a query; project_families_batch returns every genome
        // locus it recurs at (identity >= min_identity). One batched minimap2 pass. `known` empty = keep all
        // hits (incl. the self locus — a window is itself a candidate locus; grouping decides families).
        let consensuses: Vec<(String, Vec<u8>)> = windows
            .iter()
            .enumerate()
            .filter_map(|(i, (c, s, e))| {
                genome
                    .fetch_sequence(c, *s, *e)
                    .map(|seq| (format!("w{i}"), seq))
            })
            .collect();
        let known: HashMap<String, Vec<(String, u64, u64)>> = HashMap::new();
        let cov = 0.0_f64; // block length is gated by min_block below, not by fractional window coverage
        let hits = project_families_batch(
            &consensuses,
            fasta_path,
            &known,
            p.min_identity,
            cov,
            &p.minimap2,
            p.threads,
        )?;

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
                    chrom: c.clone(),
                    start: *s,
                    end: *e,
                    n_reads: 1,
                    strand: '+',
                    introns: vec![],
                    seq,
                    distinguishing_uniq: 0,
                    core_bp: 0,
                    stub: false,
                    tes: None,
                });
                window_spans.push((c.clone(), *s, *e));
            }
        }
        // discovered paralog loci that fall OUTSIDE every window (a paralog at a locus no member covers).
        let mut extra: Vec<(String, u64, u64)> = Vec::new();
        for hs in hits.into_values() {
            for h in hs {
                let len = h.end.saturating_sub(h.start);
                if len < p.min_block || len > p.max_locus_span {
                    continue;
                }
                let inside_window = window_spans
                    .iter()
                    .any(|(wc, ws, we)| *wc == h.chrom && h.start < *we && *ws < h.end);
                if inside_window {
                    continue;
                }
                extra.push((h.chrom, h.start, h.end));
            }
        }
        extra.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        for (chrom, start, end) in merge_overlapping(&extra) {
            if let Some(seq) = genome.fetch_sequence(&chrom, start, end) {
                reps.push(DenovoTranscript {
                    tid: format!("DN_{chrom}_{start}_1"),
                    chrom,
                    start,
                    end,
                    n_reads: 1,
                    strand: '+',
                    introns: vec![],
                    seq,
                    distinguishing_uniq: 0,
                    core_bp: 0,
                    stub: false,
                    tes: None,
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
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 6 {
                continue;
            }
            let (Ok(s1), Ok(e1), Ok(s2), Ok(e2)) = (
                f[1].parse::<u64>(),
                f[2].parse::<u64>(),
                f[4].parse::<u64>(),
                f[5].parse::<u64>(),
            ) else {
                continue;
            };
            if e1 > s1 {
                raw.push((f[0].to_string(), s1, e1));
            }
            if e2 > s2 {
                raw.push((f[3].to_string(), s2, e2));
            }
        }
        raw.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        Ok(merge_overlapping(&raw))
    }

    /// Single-linkage merge of overlapping genomic intervals (input sorted by (chrom, start)).
    fn merge_overlapping(loci: &[(String, u64, u64)]) -> Vec<(String, u64, u64)> {
        let mut out: Vec<(String, u64, u64)> = Vec::new();
        for (c, s, e) in loci.iter().cloned() {
            match out.last_mut() {
                Some((pc, _ps, pe)) if *pc == c && s <= *pe => {
                    *pe = (*pe).max(e);
                }
                _ => out.push((c, s, e)),
            }
        }
        out
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::types::{DetHashMap, DetHashSet};

        #[test]
        fn windows_from_sd_bed_extracts_both_sides_and_merges_overlaps() {
            let dir =
                std::env::temp_dir().join(format!("rustle_sdbed_test_{}", std::process::id()));
            std::fs::create_dir_all(&dir).unwrap();
            let bed = dir.join("sd.bed");
            std::fs::write(
                &bed,
                "chr1\t100\t200\tchr2\t500\t600\t95.0\t+\t+\n\
                 chr1\t150\t250\tchr3\t10\t20\t93.0\t+\t-\n\
                 # a comment line, ignored\n\
                 chr9\t0\t5\tchr9\t0\t5\t100.0\t+\t+\n", // degenerate zero-length-safe pair, both sides valid
            )
            .unwrap();
            let windows = windows_from_sd_bed(bed.to_str().unwrap()).unwrap();
            std::fs::remove_dir_all(&dir).ok();
            // chr1: [100,200) and [150,250) overlap -> merged to [100,250).
            assert!(
                windows.contains(&("chr1".to_string(), 100, 250)),
                "{windows:?}"
            );
            assert!(
                windows.contains(&("chr2".to_string(), 500, 600)),
                "{windows:?}"
            );
            assert!(
                windows.contains(&("chr3".to_string(), 10, 20)),
                "{windows:?}"
            );
            assert!(windows.contains(&("chr9".to_string(), 0, 5)), "{windows:?}");
            assert_eq!(
                windows.len(),
                4,
                "chr1's two overlapping sides must merge into one window: {windows:?}"
            );
        }

        #[test]
        fn merge_overlapping_joins_adjacent_and_keeps_disjoint() {
            let loci = vec![
                ("chr1".to_string(), 10, 100),
                ("chr1".to_string(), 90, 200), // overlaps previous -> merge to 10..200
                ("chr1".to_string(), 500, 600), // disjoint -> separate
                ("chr2".to_string(), 0, 50),   // different contig -> separate
            ];
            let m = merge_overlapping(&loci);
            assert_eq!(
                m,
                vec![
                    ("chr1".to_string(), 10, 200),
                    ("chr1".to_string(), 500, 600),
                    ("chr2".to_string(), 0, 50),
                ]
            );
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
            use crate::family::genome_projection::project_families_batch;
            let fa = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta";
            if std::fs::metadata(fa).is_err() {
                eprintln!("gorilla genome fasta absent; skip");
                return;
            }
            if std::process::Command::new("minimap2")
                .arg("--version")
                .output()
                .is_err()
            {
                return;
            }

            // The 31-locus NPIP oracle (o1_oracle/npip31.regions), "chrom:1based_start-1based_end" per line,
            // converted to 0-based half-open.
            let oracle_txt =
                std::fs::read_to_string("/mnt/linuxdisk/home/juanfraitu/o1_oracle/npip31.regions")
                    .expect("oracle regions file");
            let mut windows: Vec<(String, u64, u64)> = Vec::new();
            for line in oracle_txt.lines() {
                let line = line.trim();
                if line.is_empty() {
                    continue;
                }
                let (chrom, rest) = line.split_once(':').expect("chrom:start-end");
                let (s, e) = rest.split_once('-').expect("start-end");
                let s1: u64 = s.parse().unwrap();
                let e1: u64 = e.parse().unwrap();
                windows.push((chrom.to_string(), s1 - 1, e1)); // 1-based inclusive -> 0-based half-open
            }
            assert_eq!(windows.len(), 31, "expected all 31 oracle loci");

            let contigs: crate::types::DetHashSet<String> =
                windows.iter().map(|(c, _, _)| c.clone()).collect();
            let genome = GenomeIndex::from_fasta_contigs(fa, &contigs).expect("genome index");
            let consensuses: Vec<(String, Vec<u8>)> = windows
                .iter()
                .enumerate()
                .filter_map(|(i, (c, s, e))| {
                    genome
                        .fetch_sequence(c, *s, *e)
                        .map(|seq| (format!("w{i}"), seq))
                })
                .collect();
            assert_eq!(
                consensuses.len(),
                31,
                "every oracle window must fetch real sequence"
            );

            let known: HashMap<String, Vec<(String, u64, u64)>> = HashMap::new();
            let p = GenomeRepParams::default(); // min_identity 0.90, the same floor genome_reps() uses
            let hits = project_families_batch(
                &consensuses,
                fa,
                &known,
                p.min_identity,
                0.0,
                &p.minimap2,
                p.threads,
            )
            .expect("project_families_batch");

            let mut connected_to_another_locus = 0usize;
            let mut connected_absent_only = 0usize; // connects, and the query window is one of the 26 "absent" ones
                                                    // (absent = has zero de novo RNA node today; computed once, offline, against ggo.nodes.tsv --
                                                    // hardcoded here since this is a one-shot ceiling measurement, not a maintained pipeline path)
            let absent_idx: HashSet<usize> = [
                0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 15, 18, 19, 20, 21, 22, 23, 24, 26, 27,
                28, 29, 30,
            ]
            .into_iter()
            .collect();
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
                    if absent_idx.contains(&i) {
                        connected_absent_only += 1;
                    }
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
            if std::fs::metadata(fa).is_err() {
                eprintln!("gorilla genome fasta absent; skip");
                return;
            }
            if std::process::Command::new("minimap2")
                .arg("--version")
                .output()
                .is_err()
            {
                return;
            }
            let win_bed =
                "/mnt/linuxdisk/home/juanfraitu/o1_fromgenome_sd/npip_seeded_windows2.bed";
            if std::fs::metadata(win_bed).is_err() {
                eprintln!("seeded windows file absent; skip");
                return;
            }

            let mut windows: Vec<(String, u64, u64)> = Vec::new();
            for line in std::fs::read_to_string(win_bed).unwrap().lines() {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() >= 3 {
                    windows.push((
                        f[0].to_string(),
                        f[1].parse().unwrap(),
                        f[2].parse().unwrap(),
                    ));
                }
            }
            eprintln!(
                "[stage_b] {} SD-derived windows (no oracle coordinates used as input)",
                windows.len()
            );

            let p = GenomeRepParams::default();
            let reps = genome_reps(fa, &windows, &p).expect("genome_reps");
            eprintln!("[stage_b] {} reps produced", reps.len());

            // Score AFTER the fact against the 31-locus oracle -- truth used only for scoring, never as input.
            let oracle_txt =
                std::fs::read_to_string("/mnt/linuxdisk/home/juanfraitu/o1_oracle/npip31.regions")
                    .expect("oracle regions file");
            let mut oracle: Vec<(String, u64, u64)> = Vec::new();
            for line in oracle_txt.lines() {
                let line = line.trim();
                if line.is_empty() {
                    continue;
                }
                let (chrom, rest) = line.split_once(':').unwrap();
                let (s, e) = rest.split_once('-').unwrap();
                let s1: u64 = s.parse().unwrap();
                let e1: u64 = e.parse().unwrap();
                oracle.push((chrom.to_string(), s1 - 1, e1));
            }
            let absent_idx: HashSet<usize> = [
                0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 15, 18, 19, 20, 21, 22, 23, 24, 26, 27,
                28, 29, 30,
            ]
            .into_iter()
            .collect();
            let mut covered = 0usize;
            let mut covered_absent = 0usize;
            for (i, (oc, os_, oe)) in oracle.iter().enumerate() {
                let hit = reps
                    .iter()
                    .any(|r| &r.chrom == oc && r.end > *os_ && r.start < *oe);
                if hit {
                    covered += 1;
                    if absent_idx.contains(&i) {
                        covered_absent += 1;
                    }
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
            use crate::family::annotation_families::SdPairs;
            use crate::family::family_detect::family_split::{decompose_families, SplitParams};
            use crate::family::family_detect::{detect_edges, DetectParams};

            let fa = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta";
            if std::fs::metadata(fa).is_err() {
                eprintln!("gorilla genome fasta absent; skip");
                return;
            }
            if std::process::Command::new("minimap2")
                .arg("--version")
                .output()
                .is_err()
            {
                return;
            }

            let copies_tsv = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.tsv";
            let copies_fa = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.fa";
            let sedef_bed = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed";
            let sd_windows_bed =
                "/mnt/linuxdisk/home/juanfraitu/o1_fromgenome_sd/npip_seeded_windows2.bed";
            for p in [copies_tsv, copies_fa, sedef_bed, sd_windows_bed] {
                if std::fs::metadata(p).is_err() {
                    eprintln!("required real-data file {p} absent; skip");
                    return;
                }
            }

            // --- load the real, already-computed baseline RNA reps (footprint OFF, the shipped default) ---
            let tsv_text = std::fs::read_to_string(copies_tsv).unwrap();
            let fa_text = std::fs::read_to_string(copies_fa).unwrap();
            // FASTA records are in the SAME row order as the TSV (both written by the same run, one pass).
            let seqs: Vec<Vec<u8>> = fa_text
                .lines()
                .filter(|l| !l.starts_with('>'))
                .map(|l| l.as_bytes().to_vec())
                .collect();
            let mut baseline: Vec<DenovoTranscript> = Vec::new();
            for (i, line) in tsv_text.lines().skip(1).enumerate() {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 11 {
                    continue;
                }
                let chrom = f[3].to_string();
                let start: u64 = f[4].parse().unwrap();
                let end: u64 = f[5].parse().unwrap();
                let strand = f[7].chars().next().unwrap_or('+');
                let n_reads: u32 = f[8].parse().unwrap_or(1);
                let exons: Vec<(u64, u64)> = f[9]
                    .split(',')
                    .filter_map(|e| {
                        let (s, en) = e.split_once('-')?;
                        Some((s.parse().ok()?, en.parse().ok()?))
                    })
                    .collect();
                let introns: Vec<(u64, u64)> = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
                baseline.push(DenovoTranscript {
                    tid: f[2].to_string(),
                    chrom,
                    start,
                    end,
                    n_reads,
                    strand,
                    introns,
                    seq: seqs.get(i).cloned().unwrap_or_default(),
                    distinguishing_uniq: 0,
                    core_bp: 0,
                    stub: false,
                    tes: None,
                });
            }
            eprintln!(
                "[stage_c] {} baseline RNA reps loaded from the real footprint-OFF run",
                baseline.len()
            );

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
                if f.len() >= 3 {
                    sd_windows.push((
                        f[0].to_string(),
                        f[1].parse().unwrap(),
                        f[2].parse().unwrap(),
                    ));
                }
            }
            let n_before_bound = baseline.len();
            baseline.retain(|r| {
                sd_windows
                    .iter()
                    .any(|(c, s, e)| &r.chrom == c && r.end > *s && r.start < *e)
            });
            eprintln!(
                "[stage_c] bounded baseline to the SD-seeded-window region for tractable real edge-confirmation: \
                 {n_before_bound} -> {} reps",
                baseline.len()
            );

            // --- oracle ---
            let oracle_txt =
                std::fs::read_to_string("/mnt/linuxdisk/home/juanfraitu/o1_oracle/npip31.regions")
                    .unwrap();
            let oracle: Vec<(String, u64, u64)> = oracle_txt
                .lines()
                .filter(|l| !l.trim().is_empty())
                .map(|l| {
                    let (c, rest) = l.trim().split_once(':').unwrap();
                    let (s, e) = rest.split_once('-').unwrap();
                    (
                        c.to_string(),
                        s.parse::<u64>().unwrap() - 1,
                        e.parse::<u64>().unwrap(),
                    )
                })
                .collect();
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
                    for &m in &fam.members {
                        rep_family[m] = Some(fi);
                    }
                }
                let mut oracle_families: std::collections::BTreeSet<usize> =
                    std::collections::BTreeSet::new();
                let mut oracle_covered = 0usize;
                for (oc, os_, oe) in &oracle {
                    if let Some((ri, _)) = reps
                        .iter()
                        .enumerate()
                        .find(|(_, r)| &r.chrom == oc && r.end > *os_ && r.start < *oe)
                    {
                        oracle_covered += 1;
                        if let Some(fi) = rep_family[ri] {
                            oracle_families.insert(fi);
                        }
                    }
                }
                let mut false_merge_families = 0usize;
                for &fi in &oracle_families {
                    let has_foreign = families[fi].members.iter().any(|&m| {
                        !oracle.iter().any(|(oc, os_, oe)| {
                            &reps[m].chrom == oc && reps[m].end > *os_ && reps[m].start < *oe
                        })
                    });
                    if has_foreign {
                        false_merge_families += 1;
                    }
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
            let contigs: DetHashSet<String> = baseline
                .iter()
                .map(|r| r.chrom.clone())
                .chain(oracle.iter().map(|(c, _, _)| c.clone()))
                .collect();
            let genome = GenomeIndex::from_fasta_contigs(fa, &contigs).expect("genome index");
            let mut corrected = baseline.clone();
            let mut n_corrected = 0usize;
            for r in corrected.iter_mut() {
                let overlaps_oracle = oracle
                    .iter()
                    .any(|(oc, os_, oe)| &r.chrom == oc && r.end > *os_ && r.start < *oe);
                if !overlaps_oracle {
                    continue;
                }
                if let Some((cs, ce)) =
                    sd_pairs.single_span_core(&r.chrom, r.start, r.end, 15_000, 2)
                {
                    if let Some(seq) = genome.fetch_sequence(&r.chrom, cs, ce) {
                        r.start = cs;
                        r.end = ce;
                        r.introns = vec![];
                        r.seq = seq;
                        n_corrected += 1;
                    }
                }
            }
            eprintln!("[stage_c] proposal #1: corrected {n_corrected} of the baseline reps overlapping an oracle locus");
            let p1_result = score(&corrected, "proposal #1 only (spans corrected)");

            // --- PROPOSAL #2: add SD-seeded DNA-only reps (already-validated §6j5 Stage B windows, reused
            //     from `sd_windows` loaded above -- same file, same region used to bound the baseline) ---
            let t_genome_reps = std::time::Instant::now();
            let dna_reps =
                genome_reps(fa, &sd_windows, &GenomeRepParams::default()).expect("genome_reps");
            eprintln!("[stage_c/e] genome_reps took {:?}", t_genome_reps.elapsed());
            eprintln!(
                "[stage_c] proposal #2: {} SD-seeded DNA-only reps generated",
                dna_reps.len()
            );

            // --- COMBINED: corrected RNA reps + DNA-only reps, deduped where they land on the same locus ---
            let mut combined = corrected.clone();
            let mut n_dup_skipped = 0usize;
            for d in dna_reps {
                let dup = combined.iter().any(|r| {
                    r.chrom == d.chrom
                        && r.start.max(d.start) < r.end.min(d.end)
                        && (r.end.min(d.end) - r.start.max(d.start)) as f64
                            / (d.end - d.start).max(1) as f64
                            > 0.8
                });
                if dup {
                    n_dup_skipped += 1;
                    continue;
                }
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
            use crate::family::denovo_pipeline::families_from_edges;
            use std::fmt::Write as _;
            let nodes_path = std::env::var("RUSTLE_FM_NODES").expect("RUSTLE_FM_NODES");
            let edges_path = std::env::var("RUSTLE_FM_EDGES").expect("RUSTLE_FM_EDGES");
            let out_path = std::env::var("RUSTLE_FM_FAMILIES_OUT").expect("RUSTLE_FM_FAMILIES_OUT");
            let mut reps = Vec::new();
            for line in std::fs::read_to_string(&nodes_path)
                .unwrap()
                .lines()
                .skip(1)
            {
                let f: Vec<&str> = line.split('\t').collect();
                let (Ok(s), Ok(e)) = (f[2].parse::<u64>(), f[3].parse::<u64>()) else {
                    panic!("bad node row {line}")
                };
                assert_eq!(
                    f[0].parse::<usize>().unwrap(),
                    reps.len(),
                    "node indices must be 0..n in order"
                );
                reps.push(DenovoTranscript {
                    tid: format!("DNA_{}_{s}_{e}", f[1]),
                    chrom: f[1].to_string(),
                    start: s,
                    end: e,
                    n_reads: 1,
                    strand: '+',
                    introns: vec![],
                    seq: vec![],
                    distinguishing_uniq: 0,
                    core_bp: 0,
                    stub: false,
                    tes: None,
                });
            }
            let mut edges = Vec::new();
            for line in std::fs::read_to_string(&edges_path)
                .unwrap()
                .lines()
                .skip(1)
            {
                let f: Vec<&str> = line.split('\t').collect();
                edges.push((
                    f[0].parse::<usize>().unwrap(),
                    f[1].parse::<usize>().unwrap(),
                    f[2].parse::<f64>().unwrap(),
                    f[3].parse::<f64>().unwrap(),
                ));
            }
            let n_nodes = reps.len();
            let mut fams = families_from_edges(reps, &edges, 0.80, 0.50, 0.20, 2, 0);
            fams.sort_by(|a, b| {
                b.len()
                    .cmp(&a.len())
                    .then((&a[0].chrom, a[0].start).cmp(&(&b[0].chrom, b[0].start)))
            });
            let mut out = String::from("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tmax_family_identity\n");
            for (fi, fam) in fams.iter().enumerate() {
                for (ci, t) in fam.iter().enumerate() {
                    writeln!(
                        out,
                        "SDFAM{fi}\t{ci}\t{}\t{}\t{}\t{}\t1\t+\t1\t{}-{}\tNA",
                        t.tid, t.chrom, t.start, t.end, t.start, t.end
                    )
                    .unwrap();
                }
            }
            std::fs::write(&out_path, out).unwrap();
            eprintln!(
                "[stage_f:sd_partition] {n_nodes} atoms, {} edges -> {} families",
                edges.len(),
                fams.len()
            );
        }

        fn fm_dump(root: &str) {
            use crate::family::annotation_families::SdPairs;
            use crate::family::family_detect::{detect_edges, DetectParams};

            let fa = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta";
            let copies_tsv = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.tsv";
            let copies_fa = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.fa";
            let sedef_bed = "/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_sedef_final.bed";
            let sd_windows_bed =
                "/mnt/linuxdisk/home/juanfraitu/o1_fromgenome_sd/npip_seeded_windows2.bed";
            let oracle_path = "/mnt/linuxdisk/home/juanfraitu/o1_oracle/npip31.regions";
            for p in [
                fa,
                copies_tsv,
                copies_fa,
                sedef_bed,
                sd_windows_bed,
                oracle_path,
            ] {
                if std::fs::metadata(p).is_err() {
                    eprintln!("required real-data file {p} absent; skip");
                    return;
                }
            }
            eprintln!("[stage_f] root={root}");

            // --- loading: identical to stage_c_body (exon strings kept alongside for the dump) ---
            let tsv_text = std::fs::read_to_string(copies_tsv).unwrap();
            let fa_text = std::fs::read_to_string(copies_fa).unwrap();
            let seqs: Vec<Vec<u8>> = fa_text
                .lines()
                .filter(|l| !l.starts_with('>'))
                .map(|l| l.as_bytes().to_vec())
                .collect();
            let mut loaded: Vec<(DenovoTranscript, String, usize)> = Vec::new();
            for (i, line) in tsv_text.lines().skip(1).enumerate() {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 11 {
                    continue;
                }
                let chrom = f[3].to_string();
                let start: u64 = f[4].parse().unwrap();
                let end: u64 = f[5].parse().unwrap();
                let strand = f[7].chars().next().unwrap_or('+');
                let n_reads: u32 = f[8].parse().unwrap_or(1);
                let exons: Vec<(u64, u64)> = f[9]
                    .split(',')
                    .filter_map(|e| {
                        let (s, en) = e.split_once('-')?;
                        Some((s.parse().ok()?, en.parse().ok()?))
                    })
                    .collect();
                let introns: Vec<(u64, u64)> = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
                loaded.push((
                    DenovoTranscript {
                        tid: f[2].to_string(),
                        chrom,
                        start,
                        end,
                        n_reads,
                        strand,
                        introns,
                        seq: seqs.get(i).cloned().unwrap_or_default(),
                        distinguishing_uniq: 0,
                        core_bp: 0,
                        stub: false,
                        tes: None,
                    },
                    f[9].to_string(),
                    exons.len(),
                ));
            }
            let mut sd_windows: Vec<(String, u64, u64)> = Vec::new();
            for line in std::fs::read_to_string(sd_windows_bed).unwrap().lines() {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() >= 3 {
                    sd_windows.push((
                        f[0].to_string(),
                        f[1].parse().unwrap(),
                        f[2].parse().unwrap(),
                    ));
                }
            }
            let n_before_bound = loaded.len();
            loaded.retain(|(r, _, _)| {
                sd_windows
                    .iter()
                    .any(|(c, s, e)| &r.chrom == c && r.end > *s && r.start < *e)
            });
            let exon_strs: Vec<String> = loaded.iter().map(|x| x.1.clone()).collect();
            let n_exons: Vec<usize> = loaded.iter().map(|x| x.2).collect();
            let baseline: Vec<DenovoTranscript> = loaded.into_iter().map(|x| x.0).collect();
            eprintln!(
                "[stage_f] bounded {n_before_bound} -> {} baseline reps",
                baseline.len()
            );

            let oracle_txt = std::fs::read_to_string(oracle_path).unwrap();
            let oracle: Vec<(String, u64, u64)> = oracle_txt
                .lines()
                .filter(|l| !l.trim().is_empty())
                .map(|l| {
                    let (c, rest) = l.trim().split_once(':').unwrap();
                    let (s, e) = rest.split_once('-').unwrap();
                    (
                        c.to_string(),
                        s.parse::<u64>().unwrap() - 1,
                        e.parse::<u64>().unwrap(),
                    )
                })
                .collect();
            assert_eq!(oracle.len(), 31);

            let dp = DetectParams::default();
            let base_edges = detect_edges(&baseline, &dp);
            eprintln!(
                "[stage_f:baseline] detect_edges ({:?} core): {} edges",
                dp.edge_core,
                base_edges.len()
            );

            // --- PROPOSAL #1: identical to stage_c_body ---
            let sedef_text = std::fs::read_to_string(sedef_bed).unwrap();
            let sd_pairs = SdPairs::from_bed_str(&sedef_text);
            let contigs: DetHashSet<String> = baseline
                .iter()
                .map(|r| r.chrom.clone())
                .chain(oracle.iter().map(|(c, _, _)| c.clone()))
                .collect();
            let genome = GenomeIndex::from_fasta_contigs(fa, &contigs).expect("genome index");
            let mut corrected = baseline.clone();
            let mut is_corrected = vec![false; corrected.len()];
            for (ri, r) in corrected.iter_mut().enumerate() {
                let overlaps_oracle = oracle
                    .iter()
                    .any(|(oc, os_, oe)| &r.chrom == oc && r.end > *os_ && r.start < *oe);
                if !overlaps_oracle {
                    continue;
                }
                if let Some((cs, ce)) =
                    sd_pairs.single_span_core(&r.chrom, r.start, r.end, 15_000, 2)
                {
                    if let Some(seq) = genome.fetch_sequence(&r.chrom, cs, ce) {
                        r.start = cs;
                        r.end = ce;
                        r.introns = vec![];
                        r.seq = seq;
                        is_corrected[ri] = true;
                    }
                }
            }
            eprintln!(
                "[stage_f] proposal #1: corrected {} reps",
                is_corrected.iter().filter(|x| **x).count()
            );
            let p1_edges = detect_edges(&corrected, &dp);
            eprintln!(
                "[stage_f:proposal1] detect_edges ({:?} core): {} edges",
                dp.edge_core,
                p1_edges.len()
            );
            let base_kinds: Vec<&'static str> = n_exons
                .iter()
                .map(|&n| {
                    if n > 1 {
                        "rna_spliced"
                    } else {
                        "rna_single_exon"
                    }
                })
                .collect();
            let p1_kinds: Vec<&'static str> = (0..corrected.len())
                .map(|i| {
                    if is_corrected[i] {
                        "genomic_span_p1corrected"
                    } else {
                        base_kinds[i]
                    }
                })
                .collect();
            fm_write_cell(
                root,
                "baseline",
                &baseline,
                &baseline,
                &exon_strs,
                &base_kinds,
                &is_corrected.iter().map(|_| false).collect::<Vec<_>>(),
                &oracle,
                base_edges,
                &dp,
            );
            fm_write_cell(
                root,
                "proposal1",
                &corrected,
                &baseline,
                &exon_strs,
                &p1_kinds,
                &is_corrected,
                &oracle,
                p1_edges,
                &dp,
            );
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
            dp: &crate::family::family_detect::DetectParams,
        ) {
            use crate::family::family_detect::candidate_pairs;
            use crate::family::family_detect::family_graph::{longest_common_substring, upper_cow};
            use crate::family::family_detect::family_split::{
                connected_components, decompose_families, SplitParams,
            };
            use crate::family::seq_utils::reverse_complement;
            use std::fmt::Write as _;

            let dir = format!("{root}/{cell}");
            std::fs::create_dir_all(&dir).unwrap();
            let n = reps.len();
            let pairs = candidate_pairs(reps, dp);

            // exact-substring (LCS) values for both orientations of every candidate pair (cheap, deterministic).
            let lcs: Vec<(usize, usize)> = pairs
                .iter()
                .map(|&(a, b)| {
                    let au = upper_cow(&reps[a].seq);
                    let bu = upper_cow(&reps[b].seq);
                    (
                        longest_common_substring(&au, &bu),
                        longest_common_substring(&au, &reverse_complement(&bu)),
                    )
                })
                .collect();

            let edges_source = "real_detect_edges";
            let edge_val: HashMap<(usize, usize), f64> =
                edges.iter().map(|&(a, b, v)| ((a, b), v)).collect();

            let families = decompose_families(&edges, &SplitParams::default());
            let comps = connected_components(&edges, 2);
            let mut rep_comp: Vec<Option<usize>> = vec![None; n];
            for (ci, c) in comps.iter().enumerate() {
                for &m in c {
                    rep_comp[m] = Some(ci);
                }
            }
            let mut rep_family: Vec<Option<usize>> = vec![None; n];
            for (fi, fam) in families.iter().enumerate() {
                for &m in &fam.members {
                    rep_family[m] = Some(fi);
                }
            }
            let overlaps = |r: &DenovoTranscript| -> Vec<usize> {
                oracle
                    .iter()
                    .enumerate()
                    .filter(|(_, (oc, os_, oe))| &r.chrom == oc && r.end > *os_ && r.start < *oe)
                    .map(|(oi, _)| oi)
                    .collect()
            };
            let rep_oracle: Vec<Vec<usize>> = reps.iter().map(|r| overlaps(r)).collect();
            // stage_c's oracle-family rule: for each oracle locus, the FIRST rep (by index) overlapping it.
            let mut first_for: Vec<Vec<usize>> = vec![Vec::new(); n];
            let mut oracle_families: std::collections::BTreeSet<usize> =
                std::collections::BTreeSet::new();
            let mut oracle_covered = 0usize;
            for (oi, (oc, os_, oe)) in oracle.iter().enumerate() {
                if let Some((ri, _)) = reps
                    .iter()
                    .enumerate()
                    .find(|(_, r)| &r.chrom == oc && r.end > *os_ && r.start < *oe)
                {
                    oracle_covered += 1;
                    first_for[ri].push(oi);
                    if let Some(fi) = rep_family[ri] {
                        oracle_families.insert(fi);
                    }
                }
            }
            let join = |v: &[usize]| {
                v.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            };

            // reps.tsv + reps.fa
            let mut rt = String::from("idx\ttid\tchrom\tstart\tend\tstrand\tn_reads\tseq_len\tkind\tcorrected\torig_start\torig_end\tn_exons_orig\texons_orig\toracle_idx\toracle_first_rep_for\tfamily_id\tfamily_class\tcomponent_id\n");
            let mut rf = String::new();
            for (i, r) in reps.iter().enumerate() {
                let fam = rep_family[i];
                writeln!(
                    rt,
                    "{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    r.tid,
                    r.chrom,
                    r.start,
                    r.end,
                    r.strand,
                    r.n_reads,
                    r.seq.len(),
                    kinds[i],
                    is_corrected[i] as u8,
                    orig[i].start,
                    orig[i].end,
                    exon_strs[i].split(',').count(),
                    exon_strs[i],
                    join(&rep_oracle[i]),
                    join(&first_for[i]),
                    fam.map(|f| f.to_string()).unwrap_or_default(),
                    fam.map(|f| format!("{:?}", families[f].class))
                        .unwrap_or_default(),
                    rep_comp[i].map(|c| c.to_string()).unwrap_or_default()
                )
                .unwrap();
                writeln!(
                    rf,
                    ">{i}|{}|{}:{}-{}|{}\n{}",
                    r.tid,
                    r.chrom,
                    r.start,
                    r.end,
                    kinds[i],
                    String::from_utf8_lossy(&r.seq)
                )
                .unwrap();
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
                let (orient, lbp) = if fwd_frac < dp.t_core && (lr as f64 / mn as f64) > fwd_frac {
                    ("rc", lr)
                } else {
                    ("fwd", lf)
                };
                let emu = lbp as f64 / mn as f64;
                let conf = edge_val.get(&(a, b)).copied();
                let eq = conf.map(|v| (v - emu).abs() < 1e-12);
                if eq == Some(true) {
                    n_eq_lcs += 1;
                }
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
                if conf.is_some() {
                    et.push_str(&row);
                }
            }
            std::fs::write(format!("{dir}/candidates.tsv"), ct).unwrap();
            std::fs::write(format!("{dir}/edges.tsv"), et).unwrap();

            // families.tsv
            let mut ft = String::from("family_id\tclass\tn_members\tn_edges_induced\tdensity\tavg_core_recip\tcomponent_id\tcomponent_size\tmembers\tmember_tids\toracle_overlapping_members\toracle_first_rep_members\toracle_loci\tforeign_members\tis_oracle_family_stage_c\tfalse_merge_flag_stage_c\n");
            let mut false_merge_families = 0usize;
            for (fi, fam) in families.iter().enumerate() {
                let ovl: Vec<usize> = fam
                    .members
                    .iter()
                    .copied()
                    .filter(|&m| !rep_oracle[m].is_empty())
                    .collect();
                let firsts: Vec<usize> = fam
                    .members
                    .iter()
                    .copied()
                    .filter(|&m| !first_for[m].is_empty())
                    .collect();
                let foreign: Vec<usize> = fam
                    .members
                    .iter()
                    .copied()
                    .filter(|&m| rep_oracle[m].is_empty())
                    .collect();
                let mut loci: Vec<usize> = fam
                    .members
                    .iter()
                    .flat_map(|&m| rep_oracle[m].iter().copied())
                    .collect();
                loci.sort_unstable();
                loci.dedup();
                let is_of = oracle_families.contains(&fi);
                let fm = is_of && !foreign.is_empty();
                if fm {
                    false_merge_families += 1;
                }
                let ci = rep_comp[fam.members[0]];
                writeln!(
                    ft,
                    "{fi}\t{:?}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    fam.class,
                    fam.members.len(),
                    fam.stats.n_edges,
                    fam.stats.density,
                    fam.stats.avg_core_recip,
                    ci.map(|c| c.to_string()).unwrap_or_default(),
                    ci.map(|c| comps[c].len().to_string()).unwrap_or_default(),
                    join(&fam.members),
                    fam.members
                        .iter()
                        .map(|&m| reps[m].tid.as_str())
                        .collect::<Vec<_>>()
                        .join(","),
                    join(&ovl),
                    join(&firsts),
                    join(&loci),
                    join(&foreign),
                    is_of as u8,
                    fm as u8
                )
                .unwrap();
            }
            std::fs::write(format!("{dir}/families.tsv"), ft).unwrap();

            let mut st = String::new();
            writeln!(st, "key\tvalue").unwrap();
            writeln!(st, "cell\t{cell}").unwrap();
            writeln!(st, "edge_core\t{:?}", dp.edge_core).unwrap();
            writeln!(st, "edges_source\t{edges_source}").unwrap();
            writeln!(st, "n_reps\t{n}").unwrap();
            writeln!(st, "n_candidate_pairs\t{}", pairs.len()).unwrap();
            writeln!(
                st,
                "n_candidate_pairs_over_len_cap\t{}",
                pairs
                    .iter()
                    .filter(|&&(a, b)| reps[a].seq.len().max(reps[b].seq.len()) > dp.len_cap)
                    .count()
            )
            .unwrap();
            writeln!(st, "n_edges\t{}", edges.len()).unwrap();
            writeln!(st, "n_families\t{}", families.len()).unwrap();
            writeln!(st, "oracle_covered\t{oracle_covered}/31").unwrap();
            writeln!(st, "oracle_families\t{}", oracle_families.len()).unwrap();
            writeln!(
                st,
                "false_merge_families\t{false_merge_families}/{}",
                oracle_families.len()
            )
            .unwrap();
            writeln!(st, "n_edges_value_equal_lcs_emulation\t{n_eq_lcs}").unwrap();
            writeln!(
                st,
                "n_candidates_lcs_emulated_pass\t{}",
                pairs
                    .iter()
                    .enumerate()
                    .filter(|(k, &(a, b))| {
                        let mn = reps[a].seq.len().min(reps[b].seq.len()) as f64;
                        let (lf, lr) = lcs[*k];
                        (lf.max(lr) as f64 / mn) >= dp.t_core
                    })
                    .count()
            )
            .unwrap();
            std::fs::write(format!("{dir}/summary.tsv"), &st).unwrap();
            eprintln!("[stage_f:{cell}] {} reps, {} candidates, {} edges ({edges_source}), {} families; oracle {}/31 in {} families; false-merge {}/{}",
                n, pairs.len(), edges.len(), families.len(), oracle_covered, oracle_families.len(), false_merge_families, oracle_families.len());
        }

        /// Pairs phase: production `candidate_pairs` over an arbitrary rep set (copies.tsv columns + FASTA in the
        /// same row order, `$RUSTLE_FM_REPS_TSV` / `$RUSTLE_FM_REPS_FA`), with serial LCS for both orientations.
        /// Writes `<root>/<$RUSTLE_FM_CELL>/{reps.fa,candidates.tsv}` so the `bridge` phase can run on it. No threads.
        fn fm_pairs(root: &str) {
            use crate::family::family_detect::family_graph::{longest_common_substring, upper_cow};
            use crate::family::family_detect::{candidate_pairs, DetectParams};
            use crate::family::seq_utils::reverse_complement;
            use std::fmt::Write as _;
            let tsv = std::env::var("RUSTLE_FM_REPS_TSV").expect("RUSTLE_FM_REPS_TSV");
            let fa = std::env::var("RUSTLE_FM_REPS_FA").expect("RUSTLE_FM_REPS_FA");
            let cell = std::env::var("RUSTLE_FM_CELL").expect("RUSTLE_FM_CELL");
            let seqs: Vec<Vec<u8>> = std::fs::read_to_string(&fa)
                .unwrap()
                .lines()
                .filter(|l| !l.starts_with('>'))
                .map(|l| l.as_bytes().to_vec())
                .collect();
            let mut reps: Vec<DenovoTranscript> = Vec::new();
            for (i, line) in std::fs::read_to_string(&tsv)
                .unwrap()
                .lines()
                .skip(1)
                .enumerate()
            {
                let f: Vec<&str> = line.split('\t').collect();
                let exons: Vec<(u64, u64)> = f[9]
                    .split(',')
                    .filter_map(|e| {
                        let (s, en) = e.split_once('-')?;
                        Some((s.parse().ok()?, en.parse().ok()?))
                    })
                    .collect();
                reps.push(DenovoTranscript {
                    tid: f[2].to_string(),
                    chrom: f[3].to_string(),
                    start: f[4].parse().unwrap(),
                    end: f[5].parse().unwrap(),
                    n_reads: f[8].parse().unwrap_or(1),
                    strand: f[7].chars().next().unwrap_or('+'),
                    introns: exons.windows(2).map(|w| (w[0].1, w[1].0)).collect(),
                    seq: seqs[i].clone(),
                    distinguishing_uniq: 0,
                    core_bp: 0,
                    stub: false,
                    tes: None,
                });
            }
            let dp = DetectParams::default();
            let pairs = candidate_pairs(&reps, &dp);
            eprintln!(
                "[stage_f:pairs] {} reps -> {} candidate pairs",
                reps.len(),
                pairs.len()
            );
            let dir = format!("{root}/{cell}");
            std::fs::create_dir_all(&dir).unwrap();
            let mut rf = String::new();
            for (i, r) in reps.iter().enumerate() {
                writeln!(rf, ">{i}|{}\n{}", r.tid, String::from_utf8_lossy(&r.seq)).unwrap();
            }
            std::fs::write(format!("{dir}/reps.fa"), rf).unwrap();
            let mut ct = String::from(
                "i\tj\ttid_i\ttid_j\tlen_i\tlen_j\tmin_len\tmax_len\tlcs_fwd_bp\tlcs_rc_bp\n",
            );
            for &(a, b) in &pairs {
                let au = upper_cow(&reps[a].seq);
                let bu = upper_cow(&reps[b].seq);
                let (la, lb) = (reps[a].seq.len(), reps[b].seq.len());
                let lf = longest_common_substring(&au, &bu);
                let lr = longest_common_substring(&au, &reverse_complement(&bu));
                writeln!(
                    ct,
                    "{a}\t{b}\t{}\t{}\t{la}\t{lb}\t{}\t{}\t{lf}\t{lr}",
                    reps[a].tid,
                    reps[b].tid,
                    la.min(lb),
                    la.max(lb)
                )
                .unwrap();
            }
            std::fs::write(format!("{dir}/candidates.tsv"), ct).unwrap();
        }

        /// Decompose phase: the shipped `decompose_families(edges, SplitParams::default())` on an edge list
        /// (`$RUSTLE_FM_EDGES`: `i<TAB>j<TAB>core` rows, header optional), writing `family_id<TAB>class<TAB>members`
        /// to `$RUSTLE_FM_FAMILIES_OUT`. Lets two edge definitions be partitioned by the identical code.
        fn fm_decompose() {
            use crate::family::family_detect::family_split::{decompose_families, SplitParams};
            use std::fmt::Write as _;
            let edges_path = std::env::var("RUSTLE_FM_EDGES").expect("RUSTLE_FM_EDGES");
            let out_path = std::env::var("RUSTLE_FM_FAMILIES_OUT").expect("RUSTLE_FM_FAMILIES_OUT");
            let mut edges: Vec<(usize, usize, f64)> = Vec::new();
            for line in std::fs::read_to_string(&edges_path).unwrap().lines() {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 3 {
                    continue;
                }
                let (Ok(i), Ok(j), Ok(v)) = (
                    f[0].parse::<usize>(),
                    f[1].parse::<usize>(),
                    f[2].parse::<f64>(),
                ) else {
                    continue;
                };
                edges.push((i.min(j), i.max(j), v));
            }
            edges.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
            let families = decompose_families(&edges, &SplitParams::default());
            let mut out = String::from("family_id\tclass\tmembers\n");
            for (fi, fam) in families.iter().enumerate() {
                writeln!(
                    out,
                    "{fi}\t{:?}\t{}",
                    fam.class,
                    fam.members
                        .iter()
                        .map(|m| m.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                )
                .unwrap();
            }
            std::fs::write(&out_path, out).unwrap();
            eprintln!(
                "[stage_f:decompose] {} edges -> {} families",
                edges.len(),
                families.len()
            );
        }

        /// Bridge phase: exact POA-core `confirm_edge` values (EdgeCore::Poa, no budget) for a list of pairs, serial.
        fn fm_bridge(root: &str) {
            use crate::family::family_detect::family_graph::{
                contiguous_core_coverage_bounded_with, longest_common_substring, upper_cow,
                EDGE_CONFIRM_ASTAR,
            };
            use crate::family::family_detect::{confirm_edge, DetectParams, LEN_CAP, T_CORE};
            use crate::family::seq_utils::reverse_complement;
            use std::io::Write;
            let pairs_path = std::env::var("RUSTLE_FM_PAIRS")
                .unwrap_or_else(|_| format!("{root}/bridging_pairs.list"));
            let out_path = std::env::var("RUSTLE_FM_BRIDGE_OUT")
                .unwrap_or_else(|_| format!("{root}/bridging_production.tsv"));
            let mut cache: HashMap<String, Vec<Vec<u8>>> = HashMap::new();
            let mut out = std::fs::File::create(&out_path).unwrap();
            writeln!(out, "cell\ti\tj\tlen_i\tlen_j\tmin_len\tmax_len\tproduction_path\tfwd_value\tfwd_secs\trc_value\trc_secs\tproduction_core_frac\tproduction_core_bp\tcore_over_max_len\tpasses_tcore\tconfirm_edge_agrees\tlcs_fwd_bp\tlcs_rc_bp").unwrap();
            let dp = DetectParams {
                edge_core: crate::family::family_detect::EdgeCore::Poa,
                ..DetectParams::default()
            };
            for line in std::fs::read_to_string(&pairs_path).unwrap().lines() {
                if line.trim().is_empty() || line.starts_with('#') {
                    continue;
                }
                let f: Vec<&str> = line.split('\t').collect();
                let (cell, i, j): (&str, usize, usize) =
                    (f[0], f[1].parse().unwrap(), f[2].parse().unwrap());
                let seqs = cache.entry(cell.to_string()).or_insert_with(|| {
                    let txt = std::fs::read_to_string(format!("{root}/{cell}/reps.fa")).unwrap();
                    txt.lines()
                        .filter(|l| !l.starts_with('>'))
                        .map(|l| l.as_bytes().to_vec())
                        .collect()
                });
                let (a, b) = (seqs[i].clone(), seqs[j].clone());
                let (mn, mx) = (a.len().min(b.len()), a.len().max(b.len()));
                let au = upper_cow(&a).into_owned();
                let bu = upper_cow(&b).into_owned();
                let path = if mx > LEN_CAP {
                    "LCS (len_cap)"
                } else {
                    "poasta exact"
                };
                eprintln!(
                    "[stage_f:bridge] START {cell} {i} {j} len {} x {} path={path}",
                    a.len(),
                    b.len()
                );
                let t = std::time::Instant::now();
                let fwd =
                    contiguous_core_coverage_bounded_with(&au, &bu, LEN_CAP, EDGE_CONFIRM_ASTAR);
                let fwd_secs = t.elapsed().as_secs_f64();
                let mut cr = fwd;
                let (mut rcv, mut rcs) = (String::from("not_run"), String::new());
                if fwd < T_CORE {
                    let t2 = std::time::Instant::now();
                    let rc = contiguous_core_coverage_bounded_with(
                        &au,
                        &reverse_complement(&bu),
                        LEN_CAP,
                        EDGE_CONFIRM_ASTAR,
                    );
                    rcs = format!("{:.3}", t2.elapsed().as_secs_f64());
                    rcv = format!("{rc:.6}");
                    if rc > cr {
                        cr = rc;
                    }
                }
                let ce = confirm_edge(&a, &b, &dp);
                let agrees = match ce {
                    Some(v) => (v - cr).abs() < 1e-12 && cr >= T_CORE,
                    None => cr < T_CORE,
                };
                let core_bp = (cr * mn as f64).round() as usize;
                let lf = longest_common_substring(&au, &bu);
                let lr = longest_common_substring(&au, &reverse_complement(&bu));
                writeln!(out, "{cell}\t{i}\t{j}\t{}\t{}\t{mn}\t{mx}\t{path}\t{fwd:.6}\t{fwd_secs:.3}\t{rcv}\t{rcs}\t{cr:.6}\t{core_bp}\t{:.6}\t{}\t{}\t{lf}\t{lr}",
                    a.len(), b.len(), core_bp as f64 / mx as f64, (cr >= T_CORE) as u8, agrees as u8).unwrap();
                out.flush().unwrap();
                eprintln!(
                    "[stage_f:bridge] DONE  {cell} {i} {j} cr={cr:.4} in {:.1}s agrees={agrees}",
                    t.elapsed().as_secs_f64()
                );
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
            use crate::family::family_detect::{candidate_pairs, confirm_edge, DetectParams};
            use std::io::Write;
            use std::time::Instant;

            let copies_tsv = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.tsv";
            let copies_fa = "/mnt/linuxdisk/home/juanfraitu/o1_bundle6/ggo_off.copies.fa";
            let sd_windows_bed =
                "/mnt/linuxdisk/home/juanfraitu/o1_fromgenome_sd/npip_seeded_windows2.bed";
            for p in [copies_tsv, copies_fa, sd_windows_bed] {
                if std::fs::metadata(p).is_err() {
                    eprintln!("required real-data file {p} absent; skip");
                    return;
                }
            }

            // --- identical loading + bounding to §6j6's stage_c, so this reproduces the EXACT same hang ---
            let tsv_text = std::fs::read_to_string(copies_tsv).unwrap();
            let fa_text = std::fs::read_to_string(copies_fa).unwrap();
            let seqs: Vec<Vec<u8>> = fa_text
                .lines()
                .filter(|l| !l.starts_with('>'))
                .map(|l| l.as_bytes().to_vec())
                .collect();
            let mut baseline: Vec<DenovoTranscript> = Vec::new();
            for (i, line) in tsv_text.lines().skip(1).enumerate() {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 11 {
                    continue;
                }
                let chrom = f[3].to_string();
                let start: u64 = f[4].parse().unwrap();
                let end: u64 = f[5].parse().unwrap();
                let strand = f[7].chars().next().unwrap_or('+');
                let n_reads: u32 = f[8].parse().unwrap_or(1);
                let exons: Vec<(u64, u64)> = f[9]
                    .split(',')
                    .filter_map(|e| {
                        let (s, en) = e.split_once('-')?;
                        Some((s.parse().ok()?, en.parse().ok()?))
                    })
                    .collect();
                let introns: Vec<(u64, u64)> = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
                baseline.push(DenovoTranscript {
                    tid: f[2].to_string(),
                    chrom,
                    start,
                    end,
                    n_reads,
                    strand,
                    introns,
                    seq: seqs.get(i).cloned().unwrap_or_default(),
                    distinguishing_uniq: 0,
                    core_bp: 0,
                    stub: false,
                    tes: None,
                });
            }
            let mut sd_windows: Vec<(String, u64, u64)> = Vec::new();
            for line in std::fs::read_to_string(sd_windows_bed).unwrap().lines() {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() >= 3 {
                    sd_windows.push((
                        f[0].to_string(),
                        f[1].parse().unwrap(),
                        f[2].parse().unwrap(),
                    ));
                }
            }
            baseline.retain(|r| {
                sd_windows
                    .iter()
                    .any(|(c, s, e)| &r.chrom == c && r.end > *s && r.start < *e)
            });
            eprintln!(
                "[stage_d] {} bounded reps (should match §6j6's 42)",
                baseline.len()
            );
            for (i, r) in baseline.iter().enumerate() {
                eprintln!(
                    "[stage_d]   rep[{i}] {}:{}-{} len={} n_reads={}",
                    r.chrom,
                    r.start,
                    r.end,
                    r.seq.len(),
                    r.n_reads
                );
            }
            std::io::stderr().flush().ok();

            let dp = DetectParams {
                edge_core: crate::family::family_detect::EdgeCore::Poa,
                ..DetectParams::default()
            };
            let t0 = Instant::now();
            let pairs = candidate_pairs(&baseline, &dp);
            eprintln!(
                "[stage_d] candidate_pairs: {} pairs in {:?}",
                pairs.len(),
                t0.elapsed()
            );
            std::io::stderr().flush().ok();

            let mut n_confirmed = 0usize;
            for (k, &(a, b)) in pairs.iter().enumerate() {
                let la = baseline[a].seq.len();
                let lb = baseline[b].seq.len();
                eprintln!(
                    "[stage_d] START pair {k}/{} : rep[{a}] len={la} x rep[{b}] len={lb} (max={})",
                    pairs.len(),
                    la.max(lb)
                );
                std::io::stderr().flush().ok();
                let t1 = Instant::now();
                let cr = confirm_edge(&baseline[a].seq, &baseline[b].seq, &dp);
                let dt = t1.elapsed();
                if cr.is_some() {
                    n_confirmed += 1;
                }
                eprintln!(
                    "[stage_d] DONE  pair {k}/{} in {:?} -> {:?}",
                    pairs.len(),
                    dt,
                    cr
                );
                std::io::stderr().flush().ok();
            }
            eprintln!(
                "[stage_d] ALL {} pairs done, {} confirmed, total {:?}",
                pairs.len(),
                n_confirmed,
                t0.elapsed()
            );
        }

        #[test]
        fn genome_reps_finds_family_copies_with_genomic_seq_and_no_introns() {
            // Real subset fixture: 3 near-identical NCF1 copies + 2 unrelated decoys, each as its own contig.
            // genome_reps must surface the duplicated NCF1 copies as reps (with genomic seq, empty introns).
            if std::process::Command::new("minimap2")
                .arg("--version")
                .output()
                .is_err()
            {
                return;
            }
            let fa = "tests/fixtures/from_genome/subset.fa";
            if std::fs::metadata(fa).is_err() {
                eprintln!("fixture absent; skip");
                return;
            }
            // windows = each contig full-length (as in windows.bed).
            let windows: Vec<(String, u64, u64)> = [
                ("NCF1", 15440u64),
                ("NCF1B", 15319),
                ("NCF1C", 15406),
                ("DECOY1", 10001),
                ("DECOY2", 10001),
            ]
            .iter()
            .map(|(c, l)| (c.to_string(), 0u64, *l))
            .collect();
            let p = GenomeRepParams {
                min_identity: 0.90,
                min_block: 400,
                ..Default::default()
            };
            let reps = genome_reps(fa, &windows, &p).unwrap();
            // the three NCF1 copies must all be discovered as duplicated loci.
            for want in ["NCF1", "NCF1B", "NCF1C"] {
                assert!(
                    reps.iter().any(|r| r.chrom == want),
                    "missing duplicated locus {want}; got {:?}",
                    reps.iter().map(|r| r.chrom.as_str()).collect::<Vec<_>>()
                );
            }
            // every rep is genomic: empty intron chain and seq length == span.
            for r in &reps {
                assert!(
                    r.introns.is_empty(),
                    "DNA rep {} must have empty intron chain",
                    r.chrom
                );
                assert_eq!(
                    r.seq.len() as u64,
                    r.end - r.start,
                    "DNA rep {} seq must be the full genomic span",
                    r.chrom
                );
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

        use crate::family::genome_projection::CopyLocus;

        fn recip_overlap(a: &CopyLocus, b: &CopyLocus) -> f64 {
            if a.chrom != b.chrom {
                return 0.0;
            }
            let (lo, hi) = (a.start.max(b.start), a.end.min(b.end));
            if hi <= lo {
                return 0.0;
            }
            let ov = (hi - lo) as f64;
            (ov / (a.end - a.start).max(1) as f64).min(ov / (b.end - b.start).max(1) as f64)
        }

        /// Collapse reciprocal-overlap ≥ 0.50 loci (from different sibling consensuses hitting one genomic locus)
        /// into one, keeping the highest-identity survivor.
        pub fn dedup_overlapping(mut loci: Vec<CopyLocus>) -> Vec<CopyLocus> {
            loci.sort_by(|a, b| {
                b.identity
                    .partial_cmp(&a.identity)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            let mut kept: Vec<CopyLocus> = Vec::new();
            for l in loci {
                if !kept.iter().any(|k| recip_overlap(k, &l) >= 0.50) {
                    kept.push(l);
                }
            }
            kept
        }

        pub fn overlaps_any(
            chrom: &str,
            start: u64,
            end: u64,
            spans: &[(String, u64, u64)],
        ) -> bool {
            spans
                .iter()
                .any(|(c, s, e)| c == chrom && !(*s > end || *e < start))
        }

        pub fn format_allproj_row(
            family_id: &str,
            l: &CopyLocus,
            n_support: usize,
            overlaps_existing: bool,
        ) -> String {
            format!(
                "{family_id}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
                l.chrom, l.start, l.end, l.identity, n_support, overlaps_existing
            )
        }

        #[derive(Clone, Debug)]
        pub struct CopyIn {
            pub seq: Vec<u8>,
            pub chrom: String,
            pub start: u64,
            pub end: u64,
        }

        /// One `(family_id, consensus)` entry PER COPY, with the family_id repeated across its copies so
        /// `project_families_batch` unions all copies' hits under that one family key.
        pub fn all_copy_consensuses(fams: &[(String, Vec<CopyIn>)]) -> Vec<(String, Vec<u8>)> {
            fams.iter()
                .flat_map(|(fid, copies)| copies.iter().map(move |c| (fid.clone(), c.seq.clone())))
                .collect()
        }

        /// Per-family copy spans, for the projection's `known` self-exclusion (a copy projecting back onto its own
        /// catalogued locus is not a new localization).
        pub fn known_from_fams(
            fams: &[(String, Vec<CopyIn>)],
        ) -> HashMap<String, Vec<(String, u64, u64)>> {
            fams.iter()
                .map(|(fid, copies)| {
                    (
                        fid.clone(),
                        copies
                            .iter()
                            .map(|c| (c.chrom.clone(), c.start, c.end))
                            .collect(),
                    )
                })
                .collect()
        }

        #[cfg(test)]
        mod tests {
            use super::*;
            use crate::types::{DetHashMap, DetHashSet};

            #[test]
            fn dedup_overlap_and_row_format() {
                use crate::family::genome_projection::CopyLocus;
                let mk = |s: u64, e: u64, id: f64| CopyLocus {
                    chrom: "chr7".into(),
                    start: s,
                    end: e,
                    identity: id,
                    cov: 0.95,
                };
                // two overlapping hits (id .994 vs .982) + one disjoint -> 2 survivors, higher id kept
                let out = dedup_overlapping(vec![
                    mk(75976253, 75991692, 0.994),
                    mk(75976300, 75991600, 0.982),
                    mk(76360590, 76375995, 0.990),
                ]);
                assert_eq!(out.len(), 2);
                assert!(out
                    .iter()
                    .any(|l| (l.identity - 0.994).abs() < 1e-9 && l.start == 75976253));
                assert!(overlaps_any(
                    "chr7",
                    75976300,
                    75991600,
                    &[("chr7".into(), 75976253, 75991692)]
                ));
                assert!(!overlaps_any(
                    "chr7",
                    76360590,
                    76375995,
                    &[("chr7".into(), 75976253, 75991692)]
                ));
                let row = format_allproj_row("GWFAM7", &mk(75976253, 75991692, 0.994), 41, false);
                assert_eq!(row, "GWFAM7\tchr7\t75976253\t75991692\t0.994\t41\tfalse");
            }

            #[test]
            fn extraction_repeats_fid_and_collects_spans() {
                let fams = vec![
                    (
                        "GWFAM0".to_string(),
                        vec![
                            CopyIn {
                                seq: b"ACGT".to_vec(),
                                chrom: "chr1".into(),
                                start: 100,
                                end: 200,
                            },
                            CopyIn {
                                seq: b"ACGA".to_vec(),
                                chrom: "chr1".into(),
                                start: 500,
                                end: 600,
                            },
                        ],
                    ),
                    (
                        "GWFAM1".to_string(),
                        vec![CopyIn {
                            seq: b"TTTT".to_vec(),
                            chrom: "chr2".into(),
                            start: 10,
                            end: 20,
                        }],
                    ),
                ];
                let cons = all_copy_consensuses(&fams);
                // one entry per copy, family_id repeated per copy (so project_families_batch unions per family)
                assert_eq!(
                    cons,
                    vec![
                        ("GWFAM0".to_string(), b"ACGT".to_vec()),
                        ("GWFAM0".to_string(), b"ACGA".to_vec()),
                        ("GWFAM1".to_string(), b"TTTT".to_vec()),
                    ]
                );
                let known = known_from_fams(&fams);
                assert_eq!(
                    known["GWFAM0"],
                    vec![
                        ("chr1".to_string(), 100, 200),
                        ("chr1".to_string(), 500, 600)
                    ]
                );
                assert_eq!(known["GWFAM1"], vec![("chr2".to_string(), 10, 20)]);
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
            let (chrom, range) = s.rsplit_once(':').ok_or_else(|| {
                anyhow!("--seed {spec:?}: expected chrom:start-end (1-based inclusive)")
            })?;
            if chrom.is_empty() {
                return Err(anyhow!("--seed {spec:?}: empty chromosome"));
            }
            let (a, b) = range.rsplit_once('-').ok_or_else(|| {
                anyhow!("--seed {spec:?}: expected chrom:start-end (1-based inclusive)")
            })?;
            let clean = |x: &str| x.replace([',', '_'], "");
            let a: u64 = clean(a)
                .parse()
                .map_err(|_| anyhow!("--seed {spec:?}: start is not an integer"))?;
            let b: u64 = clean(b)
                .parse()
                .map_err(|_| anyhow!("--seed {spec:?}: end is not an integer"))?;
            if a == 0 {
                return Err(anyhow!(
                    "--seed {spec:?}: coordinates are 1-based, so start must be >= 1"
                ));
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
        pub fn format_seed_rows(
            seed: &SeedLocus,
            fams: &[Vec<(String, u64, u64)>],
            hit: &SeedHit,
        ) -> Vec<String> {
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
            use crate::types::{DetHashMap, DetHashSet};

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
                    SeedHit::Hit {
                        family_idx: 0,
                        copy_idx: 0,
                        overlap_bp: 100
                    }
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
                    SeedHit::Hit {
                        family_idx: 1,
                        copy_idx: 0,
                        overlap_bp: 600
                    }
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
                    SeedHit::Hit {
                        family_idx: 0,
                        copy_idx: 0,
                        overlap_bp: 100
                    }
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
                assert_eq!(
                    rows.len(),
                    2,
                    "the query returns the whole component: {rows:?}"
                );
                assert!(rows.iter().all(|r| r.contains("\tGWFAM0\t")));
                assert_eq!(
                    rows.iter().filter(|r| r.ends_with("\ttrue\t100")).count(),
                    1
                );
                assert!(
                    rows[1].contains("\tfalse\t0"),
                    "non-seed members carry overlap 0: {}",
                    rows[1]
                );
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
}

pub mod genome_projection {
    //! Genome-projection copy enumeration (spec §7) — the VG copy-number leg: LAND the family variation graph's
    //! consensus PATH onto the genome (minimap2) and count the near-identical landing sites. famCN = the number
    //! of genomic sites the family graph maps to, recovering K=0 collapses the RNA read graph merges into one
    //! locus. (minimap2 does the alignment; the VG framing is that famCN counts where the family graph lands.) In-engine
    //! minimap2 (no Liftoff dependency); seeded by our own consensus, so no reference-annotation circularity.
    //!
    //! **STATUS:** OPT-IN — `--enumerate-copies` (src/bin/gw_family_catalog.rs:172-173, default false) or `--min-identity 0.98` (gw_family_catalog.rs:138-139, Option, default Non  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)
    use anyhow::Result;
    use std::collections::HashMap;
    use std::io::Write;

    #[derive(Clone, Debug)]
    pub struct CopyLocus {
        pub chrom: String,
        pub start: u64,
        pub end: u64,
        pub identity: f64,
        pub cov: f64,
    }

    /// RAII temp-file cleanup helper shared by the single and batch projection entry points.
    struct TempFile(std::path::PathBuf);
    impl Drop for TempFile {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    /// Resolve the minimap2 TARGET for a genome-projection call, optionally reusing a pre-built splice `.mmi`
    /// instead of re-indexing the multi-GB FASTA on every invocation. A gw_family_catalog run projects up to 3x
    /// (famCN, `--project-all-families`, `--collapse-enumerate`) and minimap2 rebuilds the whole splice minimizer
    /// index from scratch on each projection call while only a few hundred short consensuses are aligned — index
    /// build dominates each call.
    ///
    /// Point `RUSTLE_PROJECT_MMI` at a `minimap2 -x splice -d genome.splice.mmi genome.fasta` index and every
    /// projection call targets it, collapsing per-call indexing to a ONE-TIME, cross-run cost (build it once, it
    /// persists and is reused by every future run). BYTE-IDENTICAL: the `.mmi` carries the same `-x splice` k/w
    /// and the map-time `-N/-p/-c/--cs` options are unaffected by indexing (the exact guarantee `RUSTLE_ABSENT_MMI`
    /// in `absent_copy.rs` already relies on — the SAME index serves both). Unset = the raw FASTA path, exactly as
    /// before. No temp index is auto-built: a splice `.mmi` of a ~3 Gb genome is >13 GB, so the caller owns the
    /// (persistent) file rather than have the binary leak one into the temp dir per run.
    fn projection_target(genome: &str) -> String {
        std::env::var("RUSTLE_PROJECT_MMI")
            .ok()
            .filter(|m| !m.is_empty())
            .unwrap_or_else(|| genome.to_string())
    }

    /// Run `minimap2 -c -x splice -N 50 -p 0.01` (query FASTA at `query_path`) against `genome_fasta` and
    /// return the raw PAF stdout. `-p 0.01`: report divergent secondaries too (default -x splice -p suppresses
    /// them, hiding all but near-identical copies); the id/cov filter downstream decides which to keep. Returns
    /// `Ok(None)` (not an error) if minimap2 exits non-zero, matching the existing graceful-degradation contract.
    fn run_minimap2_paf(
        query_path: &std::path::Path,
        genome_fasta: &str,
        minimap2: &str,
        threads: usize,
    ) -> Result<Option<String>> {
        let out = std::process::Command::new(minimap2)
            .args(["-c", "-x", "splice", "-N", "50", "-p", "0.01", "-t"])
            .arg(threads.to_string())
            .arg(projection_target(genome_fasta))
            .arg(query_path)
            .output()
            .map_err(|e| anyhow::anyhow!("minimap2 ('{minimap2}') projection failed: {e}"))?;
        if !out.status.success() {
            return Ok(None);
        }
        Ok(Some(String::from_utf8_lossy(&out.stdout).into_owned()))
    }

    /// Parse PAF text into per-query-name hit lists, filtered by identity/coverage. `qlen` (PAF field 1) is
    /// read directly per-record, so this works uniformly whether the PAF came from a single-query or a
    /// multi-query (batch) minimap2 run — no external query-length bookkeeping needed.
    fn parse_paf_hits(
        paf: &str,
        min_identity: f64,
        min_cov: f64,
    ) -> HashMap<String, Vec<CopyLocus>> {
        let mut by_query: HashMap<String, Vec<CopyLocus>> = HashMap::new();
        for line in paf.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 12 {
                continue;
            }
            let qname = f[0].to_string();
            let qlen = f[1].parse::<f64>().unwrap_or(1.0).max(1.0);
            let tname = f[5].to_string();
            let ts = f[7].parse::<u64>().unwrap_or(0);
            let te = f[8].parse::<u64>().unwrap_or(0);
            let qs = f[2].parse::<f64>().unwrap_or(0.0);
            let qe = f[3].parse::<f64>().unwrap_or(0.0);
            let de = f[12..]
                .iter()
                .find_map(|x| x.strip_prefix("de:f:").and_then(|v| v.parse::<f64>().ok()));
            let ident = de.map(|d| 1.0 - d).unwrap_or_else(|| {
                f[9].parse::<f64>().unwrap_or(0.0) / f[10].parse::<f64>().unwrap_or(1.0).max(1.0)
            });
            let cov = (qe - qs) / qlen;
            if ident >= min_identity && cov >= min_cov {
                by_query.entry(qname).or_default().push(CopyLocus {
                    chrom: tname,
                    start: ts,
                    end: te,
                    identity: ident,
                    cov,
                });
            }
        }
        by_query
    }

    /// Disjoint filter: sort by (chrom, start), drop hits overlapping an already-kept hit or a `known` locus.
    fn disjoint_filter(mut hits: Vec<CopyLocus>, known: &[(String, u64, u64)]) -> Vec<CopyLocus> {
        hits.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
        let mut kept: Vec<CopyLocus> = Vec::new();
        let overlaps =
            |c: &CopyLocus, k: &(String, u64, u64)| c.chrom == k.0 && c.start < k.2 && k.1 < c.end;
        for h in hits {
            if kept
                .iter()
                .any(|k| k.chrom == h.chrom && k.start < h.end && h.start < k.end)
            {
                continue;
            }
            if known.iter().any(|k| overlaps(&h, k)) {
                continue;
            }
            kept.push(h);
        }
        kept
    }

    /// minimap2 the consensus against the genome; keep hits with identity ≥ `min_identity`, aligned-fraction
    /// of the consensus ≥ `min_cov` (structure-preserving), disjoint from each other and from `known` loci.
    pub fn project_family_copies(
        consensus: &[u8],
        genome_fasta: &str,
        known: &[(String, u64, u64)],
        min_identity: f64,
        min_cov: f64,
        minimap2: &str,
        threads: usize,
    ) -> Result<Vec<CopyLocus>> {
        let dir = std::env::temp_dir();
        let q = dir.join(format!(
            "rustle_proj_q_{}_{}.fa",
            std::process::id(),
            consensus.len()
        ));
        let _c = TempFile(q.clone());
        {
            let mut f = std::fs::File::create(&q)?;
            writeln!(f, ">cons")?;
            f.write_all(consensus)?;
            writeln!(f)?;
        }
        let paf = match run_minimap2_paf(&q, genome_fasta, minimap2, threads)? {
            Some(p) => p,
            None => return Ok(Vec::new()),
        };
        let mut by_query = parse_paf_hits(&paf, min_identity, min_cov);
        let hits = by_query.remove("cons").unwrap_or_default();
        Ok(disjoint_filter(hits, known))
    }

    /// Batch variant: ONE minimap2 invocation (one genome index load) for MANY families' consensuses, instead
    /// of one invocation per family. Writes all consensuses to a single multi-record query FASTA (header =
    /// `family_id`), runs minimap2 once, groups PAF hits by query name, then applies the same per-family
    /// identity/coverage filter and disjoint+known-exclusion filter as `project_family_copies`.
    pub fn project_families_batch(
        consensuses: &[(String, Vec<u8>)],
        genome_fasta: &str,
        known: &HashMap<String, Vec<(String, u64, u64)>>,
        min_identity: f64,
        min_cov: f64,
        minimap2: &str,
        threads: usize,
    ) -> Result<HashMap<String, Vec<CopyLocus>>> {
        let dir = std::env::temp_dir();
        let q = dir.join(format!(
            "rustle_proj_batch_{}_{}.fa",
            std::process::id(),
            consensuses.len()
        ));
        let _c = TempFile(q.clone());
        {
            let mut f = std::fs::File::create(&q)?;
            for (fam_id, seq) in consensuses {
                if seq.is_empty() {
                    continue;
                }
                writeln!(f, ">{fam_id}")?;
                f.write_all(seq)?;
                writeln!(f)?;
            }
        }
        let empty: HashMap<String, Vec<CopyLocus>> = consensuses
            .iter()
            .map(|(id, _)| (id.clone(), Vec::new()))
            .collect();
        let paf = match run_minimap2_paf(&q, genome_fasta, minimap2, threads)? {
            Some(p) => p,
            None => return Ok(empty),
        };
        let by_query = parse_paf_hits(&paf, min_identity, min_cov);
        let no_known: Vec<(String, u64, u64)> = Vec::new();
        let mut result: HashMap<String, Vec<CopyLocus>> = HashMap::new();
        for (fam_id, _) in consensuses {
            let hits = by_query.get(fam_id).cloned().unwrap_or_default();
            let k = known
                .get(fam_id)
                .map(|v| v.as_slice())
                .unwrap_or(no_known.as_slice());
            result.insert(fam_id.clone(), disjoint_filter(hits, k));
        }
        Ok(result)
    }

    #[derive(Clone, Debug)]
    pub struct ProjHit {
        pub qname: String,
        pub chrom: String,
        pub start: u64,
        pub end: u64,
        pub identity: f64,
        pub cov: f64,
        pub cs: String,
        /// PAF query-start/end (`[qs,qe)`, forward-query coordinates) and strand (`+`/`-`) of THIS hit's aligned
        /// segment. The `cs` tag only describes `[qs,qe)` in target-forward order (query-reversed if `strand=='-'`),
        /// so any per-query-offset reader of `cs` (e.g. parcn's SUN confirm) needs these to map a forward-query
        /// offset to the right position in the cs walk.
        pub qs: u64,
        pub qe: u64,
        pub strand: char,
    }

    /// Same as `run_minimap2_paf` but with `--cs` so the PAF carries a `cs:Z:` tag per hit -- needed by parcn
    /// to read the assembly base at a copy's private (PSV) positions. Kept as a sibling runner (not a flag on
    /// `run_minimap2_paf`) so the existing callers/tests are untouched.
    fn run_minimap2_paf_cs(
        query_path: &std::path::Path,
        target: &str,
        minimap2: &str,
        threads: usize,
    ) -> Result<Option<String>> {
        let out = match std::process::Command::new(minimap2)
            .args(["-c", "--cs", "-x", "splice", "-N", "50", "-p", "0.01", "-t"])
            .arg(threads.to_string())
            .arg(projection_target(target))
            .arg(query_path)
            .output()
        {
            Ok(o) => o,
            Err(_) => return Ok(None), // minimap2 missing/not spawnable -> graceful empty (contract)
        };
        if !out.status.success() {
            return Ok(None);
        }
        Ok(Some(String::from_utf8_lossy(&out.stdout).into_owned()))
    }

    /// CIGAR/cs-retaining sibling of `project_families_batch`: same single-pass batch minimap2 invocation, but
    /// returns per-hit `ProjHit` (with the `cs:Z:` tag) instead of coordinate-only `CopyLocus`, so parcn can
    /// read the assembly base at a copy's private positions.
    pub fn project_with_cs(
        consensuses: &[(String, Vec<u8>)],
        target: &str,
        min_identity: f64,
        min_cov: f64,
        minimap2: &str,
        threads: usize,
    ) -> Result<Vec<ProjHit>> {
        // `pid + consensuses.len()` alone can collide when two tests in the same process (default parallel
        // `cargo test`) both call this with equal-length consensus lists -- a monotonic counter makes the
        // query-FASTA path unique per call, avoiding a write/delete race on a shared temp file.
        static CALL_COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let call_id = CALL_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let dir = std::env::temp_dir();
        let q = dir.join(format!(
            "rustle_parcn_q_{}_{}_{}.fa",
            std::process::id(),
            consensuses.len(),
            call_id
        ));
        let _c = TempFile(q.clone());
        {
            let mut f = std::fs::File::create(&q)?;
            for (name, seq) in consensuses {
                if seq.is_empty() {
                    continue;
                }
                writeln!(f, ">{name}")?;
                f.write_all(seq)?;
                writeln!(f)?;
            }
        }
        let qlen: std::collections::HashMap<&str, usize> = consensuses
            .iter()
            .map(|(n, s)| (n.as_str(), s.len()))
            .collect();
        let paf = match run_minimap2_paf_cs(&q, target, minimap2, threads)? {
            Some(p) => p,
            None => return Ok(Vec::new()),
        };
        let mut hits = Vec::new();
        for line in paf.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 11 {
                continue;
            }
            let qname = f[0].to_string();
            let (qs, qe): (u64, u64) = (f[2].parse().unwrap_or(0), f[3].parse().unwrap_or(0));
            let strand = f.get(4).and_then(|s| s.chars().next()).unwrap_or('+');
            let (ts, te): (u64, u64) = (f[7].parse().unwrap_or(0), f[8].parse().unwrap_or(0));
            let (matches, blk): (f64, f64) =
                (f[9].parse().unwrap_or(0.0), f[10].parse().unwrap_or(1.0));
            let identity = if blk > 0.0 { matches / blk } else { 0.0 };
            let ql = match qlen.get(f[0]) {
                Some(&l) if l > 0 => l as f64,
                _ => continue,
            };
            let cov = (qe.saturating_sub(qs)) as f64 / ql;
            if identity < min_identity || cov < min_cov {
                continue;
            }
            let cs = f
                .iter()
                .find_map(|t| t.strip_prefix("cs:Z:"))
                .unwrap_or("")
                .to_string();
            hits.push(ProjHit {
                qname,
                chrom: f[5].to_string(),
                start: ts,
                end: te,
                identity,
                cov,
                cs,
                qs,
                qe,
                strand,
            });
        }
        Ok(hits)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn projection_enumerates_disjoint_copies() {
            if std::process::Command::new("minimap2")
                .arg("--version")
                .output()
                .is_err()
            {
                return;
            }
            // Build a tiny 2-contig genome file with the SAME 1kb sequence at two loci = 2 genomic copies.
            let dir = std::env::temp_dir().join(format!("rustle_proj_{}", std::process::id()));
            std::fs::create_dir_all(&dir).unwrap();
            // NOTE: a naive affine index formula like `(i*7+3)%4` is periodic with period 4 (any `i*k+c mod 4`
            // with k odd cycles every 4 bases) -- i.e. literally "TGCA" repeated. minimap2 then finds a match at
            // every 4bp-phase-aligned offset across the whole c1 contig, bridging the two real tandem copies into
            // one contiguous disjoint-kept interval (verified: 51 raw PAF hits -> 2 kept, not 3). Use a splitmix64
            // hash per index instead so the sequence is non-degenerate, deterministic, and dependency-free.
            let splitmix = |i: u64| -> u64 {
                let mut z = i.wrapping_add(0x9E3779B97F4A7C15);
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
                z ^ (z >> 31)
            };
            let seq: String = (0..1000u64)
                .map(|i| "ACGT".as_bytes()[(splitmix(i) % 4) as usize] as char)
                .collect();
            let fa = dir.join("g.fa");
            std::fs::write(&fa, format!(">c1\n{seq}{seq}\n>c2\n{seq}\n")).unwrap(); // c1 has 2 tandem copies, c2 has 1
            let hits = project_family_copies(
                seq.as_bytes(),
                fa.to_str().unwrap(),
                &[],
                0.98,
                0.90,
                "minimap2",
                2,
            )
            .unwrap();
            assert!(
                hits.len() >= 3,
                "3 genomic copies (2 on c1, 1 on c2) expected, got {}",
                hits.len()
            );
            let _ = std::fs::remove_dir_all(&dir);
        }

        /// Regression for the `-p` (secondary-to-primary score ratio) bug: `minimap2 -x splice`'s DEFAULT `-p`
        /// suppresses divergent secondary alignments whose chain overlaps the (higher-scoring) near-identical
        /// hit in query space -- so a family with one 0%-divergent copy and two divergent copies (~8%, ~15%)
        /// projects to only the identical copy, silently dropping the divergent ones. `-p 0.01` makes minimap2
        /// report them; the existing identity/coverage filter then decides which to keep.
        #[test]
        fn projection_finds_divergent_copies_not_just_identical() {
            if std::process::Command::new("minimap2")
                .arg("--version")
                .output()
                .is_err()
            {
                return;
            }

            // Deterministic splitmix64-based non-degenerate sequence generator (see note on the test above:
            // avoid periodic index formulas that alias into a repeated motif and confuse minimap2).
            let splitmix = |i: u64| -> u64 {
                let mut z = i.wrapping_add(0x9E3779B97F4A7C15);
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
                z ^ (z >> 31)
            };
            let bases = [b'A', b'C', b'G', b'T'];
            let gen_seq = |seed: u64, len: u64| -> Vec<u8> {
                (0..len)
                    .map(|i| {
                        bases[(splitmix(seed.wrapping_mul(0x2545_F491_4F6C_DD1D).wrapping_add(i))
                            % 4) as usize]
                    })
                    .collect::<Vec<u8>>()
            };
            // Mutate an EXACT `frac` of positions to a different base (rank-selected by hash, not a per-base
            // coin flip) so the realized divergence matches `frac` precisely regardless of sequence length.
            let mutate = |seq: &[u8], frac: f64, seed: u64| -> Vec<u8> {
                let n = seq.len();
                let mut idx: Vec<usize> = (0..n).collect();
                idx.sort_by_key(|&i| splitmix(seed.wrapping_add(i as u64)));
                let n_mut = (frac * n as f64).round() as usize;
                let mut out = seq.to_vec();
                for &i in idx.iter().take(n_mut) {
                    let orig = out[i];
                    let h = splitmix(seed ^ 0xABCDEF ^ i as u64);
                    let mut alt = bases[(h % 4) as usize];
                    if alt == orig {
                        alt = bases[((h % 4) as usize + 1) % 4];
                    }
                    out[i] = alt;
                }
                out
            };

            let copy_len = 1500u64;
            let copy0 = gen_seq(1, copy_len); // 0% divergent: the consensus itself
            let copy8 = mutate(&copy0, 0.08, 100); // ~8% divergent
            let copy15 = mutate(&copy0, 0.15, 200); // ~15% divergent

            // Assemble one contig: bg + copy0 + ~20kb bg + copy8 + ~20kb bg + copy15 + bg, so the three copies
            // sit at distinct loci ~20kb apart, not adjacent/tandem.
            let dir = std::env::temp_dir().join(format!("rustle_proj_div_{}", std::process::id()));
            std::fs::create_dir_all(&dir).unwrap();
            let mut genome: Vec<u8> = Vec::new();
            genome.extend_from_slice(&gen_seq(9001, 2_000));
            genome.extend_from_slice(&copy0);
            genome.extend_from_slice(&gen_seq(9002, 20_000));
            genome.extend_from_slice(&copy8);
            genome.extend_from_slice(&gen_seq(9003, 20_000));
            genome.extend_from_slice(&copy15);
            genome.extend_from_slice(&gen_seq(9004, 2_000));

            let fa = dir.join("g.fa");
            std::fs::write(
                &fa,
                format!(">chr1\n{}\n", String::from_utf8(genome).unwrap()),
            )
            .unwrap();

            let hits =
                project_family_copies(&copy0, fa.to_str().unwrap(), &[], 0.80, 0.50, "minimap2", 2)
                    .unwrap();
            assert!(
                hits.len() >= 3,
                "expected >=3 disjoint copies (identical + ~8% + ~15% divergent), got {}: {:?}",
                hits.len(),
                hits
            );
            let _ = std::fs::remove_dir_all(&dir);
        }

        /// TDD for `project_families_batch`: ONE minimap2 pass must correctly enumerate + bucket-by-identity
        /// loci for MULTIPLE families sharing one query FASTA, without cross-family leakage.
        #[test]
        fn project_families_batch_one_pass_buckets_by_identity() {
            if std::process::Command::new("minimap2")
                .arg("--version")
                .output()
                .is_err()
            {
                return;
            }

            let splitmix = |i: u64| -> u64 {
                let mut z = i.wrapping_add(0x9E3779B97F4A7C15);
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
                z ^ (z >> 31)
            };
            let bases = [b'A', b'C', b'G', b'T'];
            let gen_seq = |seed: u64, len: u64| -> Vec<u8> {
                (0..len)
                    .map(|i| {
                        bases[(splitmix(seed.wrapping_mul(0x2545_F491_4F6C_DD1D).wrapping_add(i))
                            % 4) as usize]
                    })
                    .collect::<Vec<u8>>()
            };
            let mutate = |seq: &[u8], frac: f64, seed: u64| -> Vec<u8> {
                let n = seq.len();
                let mut idx: Vec<usize> = (0..n).collect();
                idx.sort_by_key(|&i| splitmix(seed.wrapping_add(i as u64)));
                let n_mut = (frac * n as f64).round() as usize;
                let mut out = seq.to_vec();
                for &i in idx.iter().take(n_mut) {
                    let orig = out[i];
                    let h = splitmix(seed ^ 0xABCDEF ^ i as u64);
                    let mut alt = bases[(h % 4) as usize];
                    if alt == orig {
                        alt = bases[((h % 4) as usize + 1) % 4];
                    }
                    out[i] = alt;
                }
                out
            };

            let copy_len = 1500u64;
            // Family F1: sequence "A" at 0% / ~8% / ~15% divergence.
            let a0 = gen_seq(1, copy_len);
            let a8 = mutate(&a0, 0.08, 100);
            let a15 = mutate(&a0, 0.15, 200);
            // A SHORT ~50%-length near-identical fragment of A (first half, ~1% mutated): pins the fragment
            // concern -- it clears the totalCN floor (cov>=0.50) but must be EXCLUDED from famCN (cov<0.90).
            let a_frag = mutate(&a0[..(copy_len as usize / 2)], 0.01, 400);
            // Family F2: a DIFFERENT sequence "B" (distinct seed) at 0% / ~2% divergence.
            let b0 = gen_seq(500_001, copy_len);
            let b2 = mutate(&b0, 0.02, 300);

            let dir =
                std::env::temp_dir().join(format!("rustle_proj_batch_test_{}", std::process::id()));
            std::fs::create_dir_all(&dir).unwrap();
            let mut genome: Vec<u8> = Vec::new();
            genome.extend_from_slice(&gen_seq(9001, 2_000));
            genome.extend_from_slice(&a0);
            genome.extend_from_slice(&gen_seq(9002, 20_000));
            genome.extend_from_slice(&a8);
            genome.extend_from_slice(&gen_seq(9003, 20_000));
            genome.extend_from_slice(&a15);
            genome.extend_from_slice(&gen_seq(9007, 20_000));
            genome.extend_from_slice(&a_frag);
            genome.extend_from_slice(&gen_seq(9004, 20_000));
            genome.extend_from_slice(&b0);
            genome.extend_from_slice(&gen_seq(9005, 20_000));
            genome.extend_from_slice(&b2);
            genome.extend_from_slice(&gen_seq(9006, 2_000));

            let fa = dir.join("g.fa");
            std::fs::write(
                &fa,
                format!(">chr1\n{}\n", String::from_utf8(genome).unwrap()),
            )
            .unwrap();

            let consensuses: Vec<(String, Vec<u8>)> = vec![
                ("F1".to_string(), a0.clone()),
                ("F2".to_string(), b0.clone()),
            ];
            let known: std::collections::HashMap<String, Vec<(String, u64, u64)>> =
                std::collections::HashMap::new();
            let result = project_families_batch(
                &consensuses,
                fa.to_str().unwrap(),
                &known,
                0.80,
                0.50,
                "minimap2",
                2,
            )
            .unwrap();

            let f1 = result
                .get("F1")
                .expect("F1 must be present in batch result");
            let f2 = result
                .get("F2")
                .expect("F2 must be present in batch result");

            assert!(
                f1.len() >= 3,
                "F1 expected >=3 loci (0%/~8%/~15% divergent), got {}: {:?}",
                f1.len(),
                f1
            );
            assert!(
                f2.len() >= 2,
                "F2 expected >=2 loci (0%/~2% divergent), got {}: {:?}",
                f2.len(),
                f2
            );

            // Identity spans: F1 should carry a near-1.0, a ~0.92, and a ~0.85 hit.
            let f1_has = |lo: f64, hi: f64| f1.iter().any(|c| c.identity >= lo && c.identity <= hi);
            assert!(f1_has(0.97, 1.0), "F1 missing ~1.0-identity hit: {:?}", f1);
            assert!(
                f1_has(0.88, 0.96),
                "F1 missing ~0.92-identity hit: {:?}",
                f1
            );
            assert!(
                f1_has(0.80, 0.88),
                "F1 missing ~0.85-identity hit: {:?}",
                f1
            );

            // DISJOINT identity bands for F2, so the two asserts pin two distinct hits (not one): the
            // identical copy sits in [0.99,1.0], the ~2%-divergent copy in [0.96,0.99).
            let f2_has = |lo: f64, hi: f64| f2.iter().any(|c| c.identity >= lo && c.identity < hi);
            assert!(
                f2_has(0.99, 1.0001),
                "F2 missing ~1.0-identity hit: {:?}",
                f2
            );
            assert!(
                f2_has(0.96, 0.99),
                "F2 missing ~0.98-identity hit: {:?}",
                f2
            );

            // CopyLocus.cov must be populated: F1's full-length copies align (nearly) the whole query.
            let f1_full: Vec<&CopyLocus> = f1.iter().filter(|c| c.cov >= 0.90).collect();
            assert!(
                f1_full.len() >= 3,
                "F1 expected >=3 FULL-LENGTH copies (cov>=0.90), got {}: {:?}",
                f1_full.len(),
                f1
            );
            assert!(
                f1.iter().any(|c| c.identity >= 0.97 && c.cov >= 0.95),
                "F1's identical copy should have cov close to 1.0: {:?}",
                f1
            );

            // The planted ~50% fragment (cov~0.5, id~0.99) must appear (clears totalCN floor cov>=0.50) but
            // must be EXCLUDED from the famCN bucket (cov<0.90), even though its identity is >=0.98.
            let f1_famcn = f1
                .iter()
                .filter(|c| c.identity >= 0.98 && c.cov >= 0.90)
                .count();
            let f1_frag_hits = f1.iter().filter(|c| c.cov >= 0.40 && c.cov < 0.75).count();
            assert!(
                f1_frag_hits >= 1,
                "F1 half-length fragment (cov~0.5) should be present: {:?}",
                f1
            );
            // famCN counts only the full-length near-identical copy (a0); the fragment must NOT inflate it.
            assert!(f1_famcn <= f1_full.len(),
                "famCN bucket (id>=0.98,cov>=0.90) must exclude the half-length fragment: famCN_loci={f1_famcn}, full={} {:?}",
                f1_full.len(), f1);
            assert!(f1.iter().any(|c| c.identity >= 0.98 && c.cov < 0.90),
                "expected the planted fragment to be a >=0.98-identity but <0.90-cov hit (famCN-excluded): {:?}", f1);

            // No cross-family leakage: F1's loci must sit in the A-region of the contig (a0/a8/a15/a_frag
            // offsets, well before the B region begins), and F2's loci must sit in the B region.
            // Layout: 2000 bg | a0 | 20k bg | a8 | 20k bg | a15 | 20k bg | a_frag | 20k bg | b0 | ...
            let b_region_start = 2_000
                + copy_len
                + 20_000
                + copy_len
                + 20_000
                + copy_len
                + 20_000
                + (copy_len / 2)
                + 20_000;
            assert!(
                f1.iter().all(|c| c.start < b_region_start),
                "F1 locus leaked into the B region (query-name grouping bug): {:?}",
                f1
            );
            assert!(
                f2.iter().all(|c| c.start >= b_region_start - 100),
                "F2 locus leaked into the A region (query-name grouping bug): {:?}",
                f2
            );

            // Bucketing: F1 must have >=2 loci that are >=0.80 identity but <0.98 (the divergent 8%/15%
            // copies) -- proving totalCN (>=0.80) > famCN (>=0.98) is derivable from one projection result.
            let f1_divergent_but_above_floor = f1
                .iter()
                .filter(|c| c.identity >= 0.80 && c.identity < 0.98)
                .count();
            assert!(f1_divergent_but_above_floor >= 2,
                "F1 expected >=2 loci with 0.80<=identity<0.98 (totalCN>famCN bucketing), got {}: {:?}",
                f1_divergent_but_above_floor, f1);

            let _ = std::fs::remove_dir_all(&dir);
        }

        /// TDD for the `gw_family_catalog --enumerate-copies` totalCN cov floor: a coverage sweep on real
        /// families showed `min_cov=0.50` INFLATES totalCN with partial/domain-fragment hits (GWFAM18: 16 at
        /// cov0.50 vs 12=truth at cov>=0.70), while `min_cov=0.90` is too strict and drops divergent
        /// full-length copies (GSTM 4 vs 19 truth). `min_cov=0.80` is the sweet spot. This test proves that
        /// policy at the `project_families_batch` call level: with `min_cov=0.80`, the planted ~50%-length
        /// near-identical fragment (cov~0.5) must be EXCLUDED from the batch result entirely, while the
        /// full-length divergent copies (0%/~8%/~15% divergence, cov>=0.80) are retained.
        #[test]
        fn project_families_batch_cov80_excludes_fragments() {
            if std::process::Command::new("minimap2")
                .arg("--version")
                .output()
                .is_err()
            {
                return;
            }

            let splitmix = |i: u64| -> u64 {
                let mut z = i.wrapping_add(0x9E3779B97F4A7C15);
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
                z ^ (z >> 31)
            };
            let bases = [b'A', b'C', b'G', b'T'];
            let gen_seq = |seed: u64, len: u64| -> Vec<u8> {
                (0..len)
                    .map(|i| {
                        bases[(splitmix(seed.wrapping_mul(0x2545_F491_4F6C_DD1D).wrapping_add(i))
                            % 4) as usize]
                    })
                    .collect::<Vec<u8>>()
            };
            let mutate = |seq: &[u8], frac: f64, seed: u64| -> Vec<u8> {
                let n = seq.len();
                let mut idx: Vec<usize> = (0..n).collect();
                idx.sort_by_key(|&i| splitmix(seed.wrapping_add(i as u64)));
                let n_mut = (frac * n as f64).round() as usize;
                let mut out = seq.to_vec();
                for &i in idx.iter().take(n_mut) {
                    let orig = out[i];
                    let h = splitmix(seed ^ 0xABCDEF ^ i as u64);
                    let mut alt = bases[(h % 4) as usize];
                    if alt == orig {
                        alt = bases[((h % 4) as usize + 1) % 4];
                    }
                    out[i] = alt;
                }
                out
            };

            let copy_len = 1500u64;
            let a0 = gen_seq(1, copy_len);
            let a8 = mutate(&a0, 0.08, 100);
            let a15 = mutate(&a0, 0.15, 200);
            // Same ~50%-length near-identical fragment as the sibling test: clears the old cov>=0.50 floor
            // but must NOT clear the new cov>=0.80 floor.
            let a_frag = mutate(&a0[..(copy_len as usize / 2)], 0.01, 400);

            let dir =
                std::env::temp_dir().join(format!("rustle_proj_cov80_test_{}", std::process::id()));
            std::fs::create_dir_all(&dir).unwrap();
            let mut genome: Vec<u8> = Vec::new();
            genome.extend_from_slice(&gen_seq(9001, 2_000));
            genome.extend_from_slice(&a0);
            genome.extend_from_slice(&gen_seq(9002, 20_000));
            genome.extend_from_slice(&a8);
            genome.extend_from_slice(&gen_seq(9003, 20_000));
            genome.extend_from_slice(&a15);
            genome.extend_from_slice(&gen_seq(9007, 20_000));
            genome.extend_from_slice(&a_frag);
            genome.extend_from_slice(&gen_seq(9004, 2_000));

            let fa = dir.join("g.fa");
            std::fs::write(
                &fa,
                format!(">chr1\n{}\n", String::from_utf8(genome).unwrap()),
            )
            .unwrap();

            let consensuses: Vec<(String, Vec<u8>)> = vec![("F1".to_string(), a0.clone())];
            let known: std::collections::HashMap<String, Vec<(String, u64, u64)>> =
                std::collections::HashMap::new();
            let result = project_families_batch(
                &consensuses,
                fa.to_str().unwrap(),
                &known,
                0.80,
                0.80,
                "minimap2",
                2,
            )
            .unwrap();

            let f1 = result
                .get("F1")
                .expect("F1 must be present in batch result");

            // The full-length divergent copies (0%/~8%/~15%) all clear cov>=0.80 and must be present.
            assert!(f1.len() >= 3,
                "expected >=3 full-length loci (0%/~8%/~15% divergent) at min_cov=0.80, got {}: {:?}", f1.len(), f1);
            assert!(
                f1.iter().all(|c| c.cov >= 0.80),
                "min_cov=0.80 call must not return any sub-0.80-coverage hit: {:?}",
                f1
            );
            // The ~50%-length fragment must be gone entirely -- no hit with cov in the fragment's ~0.5 band.
            assert!(
                f1.iter().all(|c| c.cov < 0.40 || c.cov >= 0.75),
                "min_cov=0.80 must exclude the half-length (cov~0.5) fragment hit: {:?}",
                f1
            );

            let _ = std::fs::remove_dir_all(&dir);
        }

        #[test]
        fn project_with_cs_returns_cs_tag() {
            if std::process::Command::new("minimap2")
                .arg("--version")
                .output()
                .is_err()
            {
                return;
            }
            // Build a tiny "genome" with two near-identical copies of a query, one carrying a SNV.
            let dir = std::env::temp_dir();
            let g = dir.join(format!("parcn_g_{}.fa", std::process::id()));
            // 300bp non-degenerate query; genome = query at 0, query+SNV at 500 (padded with Ns).
            let q: Vec<u8> = (0..300)
                .map(|i| b"ACGT"[((i * 2654435761usize) >> 13) & 3])
                .collect();
            let mut snv = q.clone();
            snv[150] ^= 0b100; // flip a base
            let pad = vec![b'A'; 200];
            let mut gseq = q.clone();
            gseq.extend(&pad);
            gseq.extend(&snv);
            gseq.extend(&pad);
            std::fs::write(&g, format!(">c1\n{}\n", String::from_utf8_lossy(&gseq))).unwrap();
            let hits = project_with_cs(
                &[("F|0".into(), q.clone())],
                g.to_str().unwrap(),
                0.90,
                0.80,
                "minimap2",
                2,
            )
            .unwrap();
            std::fs::remove_file(&g).ok();
            assert!(hits.len() >= 2, "should find both genomic copies");
            assert!(
                hits.iter().all(|h| !h.cs.is_empty()),
                "each hit carries a cs tag"
            );
            assert!(hits.iter().all(|h| h.qname == "F|0"));
        }

        #[test]
        fn project_with_cs_missing_minimap2_is_graceful() {
            // a binary that cannot be spawned must yield an empty result, not an Err.
            let got = project_with_cs(
                &[("F|0".into(), b"ACGTACGTACGT".to_vec())],
                "/nonexistent/genome.fa",
                0.9,
                0.8,
                "definitely_not_a_real_minimap2_xyz",
                1,
            );
            assert!(got.is_ok(), "spawn failure must not error");
            assert!(
                got.unwrap().is_empty(),
                "spawn failure must yield empty hits"
            );
        }
    }
}

pub mod copy_graph {
    //! Copy-graph objects (v1): every copy of a family is a tagged, corroborable PATH in one GFA 1.1
    //! variation graph. A REFERENCE walk makes a reference-absent copy visibly an arm the reference does
    //! not take. Pure builder — no I/O; the caller fills the parallel vectors and writes the strings.
    //!
    //! **STATUS:** OPT-IN — --phase (src/bin/copy_assign.rs:235-236, default_value_t = false)  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    use std::collections::{BTreeMap, BTreeSet};

    /// Neutral faint colour for allele nodes observed ONLY in reads (carried by no CopyPath and unequal to
    /// the reference allele). Distinct from the backbone light-grey so read-only arms stay legible in Bandage.
    const READ_ONLY_COLOUR: &str = "#e8eaed";

    /// Per-copy status across the (in-genome / annotated) axes and the absent subtypes.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum CopyStatus {
        Reference,
        InGenomeAnnotated,
        InGenomeUnannotated,
        AnnotationUnknown,
        AbsentCollapsed,
        AbsentDivergent,
    }

    impl CopyStatus {
        /// `ST:Z:` tag value.
        pub fn tag(&self) -> &'static str {
            match self {
                CopyStatus::Reference => "reference",
                CopyStatus::InGenomeAnnotated => "in-genome-annotated",
                CopyStatus::InGenomeUnannotated => "in-genome-unannotated",
                CopyStatus::AnnotationUnknown => "annotation-unknown",
                CopyStatus::AbsentCollapsed => "absent-collapsed",
                CopyStatus::AbsentDivergent => "absent-divergent",
            }
        }
        pub fn is_absent(&self) -> bool {
            matches!(
                self,
                CopyStatus::AbsentCollapsed | CopyStatus::AbsentDivergent
            )
        }
        /// Bandage colour for arms unique to this status.
        pub fn colour(&self) -> &'static str {
            match self {
                CopyStatus::Reference => "#9aa0a6",
                CopyStatus::AbsentCollapsed | CopyStatus::AbsentDivergent => "#d93025",
                CopyStatus::InGenomeUnannotated => "#188038",
                CopyStatus::InGenomeAnnotated => "#1a73e8",
                CopyStatus::AnnotationUnknown => "#a142f4",
            }
        }
    }

    /// Corroboration evidence carried as GFA tags. `None` => tag omitted (never faked).
    #[derive(Clone, Debug, Default)]
    pub struct Corrob {
        pub reads: Option<u32>,        // RC:i:
        pub suns: Option<u32>,         // SU:i: (filled by the builder if left None)
        pub map_identity: Option<f64>, // MI:f:
    }

    /// One PSV column, already known to be usable (genome_pos + ref_allele both Some).
    #[derive(Clone, Debug)]
    pub struct PsvColumn {
        pub col: usize, // original column index (for provenance only)
        pub genome_pos: Option<u64>,
        pub ref_allele: Option<u8>,
    }

    /// One copy as a path: its allele per column (None = gap => routes through the reference allele node).
    #[derive(Clone, Debug)]
    pub struct CopyPath {
        pub id: String,
        pub alleles: Vec<Option<u8>>,
        pub status: CopyStatus,
        pub corrob: Corrob,
    }

    /// Per-read significance certificate carried onto the audit W-line.
    #[derive(Clone, Debug)]
    pub struct ReadCert {
        pub p_value: f64,
        pub min_p_value: f64,
        pub status: crate::family::copy_assign::AssignStatus,
    }

    /// One read as a walk over the columns it observed (None = unobserved).
    #[derive(Clone, Debug)]
    pub struct ReadWalk {
        pub name: String,
        pub obs: Vec<Option<u8>>,
        pub assigned_copy: Option<usize>, // index into CopyGraph.copies; None = tied/K=0 (grey)
        /// Significance certificate from copy_assign (p_value/min_p_value/status). `None` => the old,
        /// untagged W-line (backward-compatible / opt-in); `Some` appends CP/PV/MP/ST tags.
        pub cert: Option<ReadCert>,
    }

    /// A shared exon node in the exon presence/absence graph — one genomic exon interval.
    #[derive(Clone, Debug)]
    pub struct ExonClass {
        pub chrom: String,
        pub start: u64,
        pub end: u64,
    }

    /// A copy as an ordered walk over exon-class indices.
    #[derive(Clone, Debug)]
    pub struct CopyExonPath {
        pub id: String,
        pub exon_nodes: Vec<usize>,
        pub status: CopyStatus,
        pub corrob: Corrob,
    }

    /// Family exon presence/absence graph. `nodes` sorted by genomic start; each copy walks a subset.
    #[derive(Clone, Debug)]
    pub struct ExonGraph {
        pub family: String,
        pub nodes: Vec<ExonClass>,
        pub copies: Vec<CopyExonPath>,
    }

    /// A whole family's variation graph. columns, every copy.alleles, every read.obs are length M and
    /// share column order; backbone is length M+1.
    #[derive(Clone, Debug)]
    pub struct CopyGraph {
        pub family: String,
        pub columns: Vec<PsvColumn>,
        pub backbone: Vec<Vec<u8>>,
        pub copies: Vec<CopyPath>,
        pub reads: Vec<ReadWalk>,
    }

    /// Assembled GFA line groups (dedup + ordering handled by the caller/writer).
    #[derive(Default, Debug)]
    pub struct GfaLines {
        pub header: String,
        pub segs: Vec<String>,
        pub links: Vec<String>,
        pub paths: Vec<String>,
        pub walks: Vec<String>,
    }

    impl CopyGraph {
        fn m(&self) -> usize {
            self.columns.len()
        }
        fn bb(&self, i: usize) -> String {
            format!("{}_bb{}", self.family, i)
        }
        fn allele_node(&self, ci: usize, b: u8) -> String {
            format!("{}_c{}_{}", self.family, ci, b as char)
        }

        /// Number of columns where copy `c`'s allele is BOTH unique among the family's copies AND differs
        /// from the reference allele — a private *divergent* marker (SUN). A copy identical to the
        /// reference scores 0. (This reference-exclusion filter deviates from the plan's illustrative
        /// snippet but matches the acceptance test's SUN semantics.)
        fn private_columns(&self, c: usize) -> u32 {
            let mut n = 0u32;
            for ci in 0..self.m() {
                let Some(Some(b)) = self.copies[c].alleles.get(ci) else {
                    continue;
                };
                // Only count columns where this allele differs from the reference
                if self.columns[ci].ref_allele == Some(*b) {
                    continue;
                }
                // Check if no other copy has this same allele
                let unique =
                    self.copies.iter().enumerate().all(|(k, other)| {
                        k == c || other.alleles.get(ci).and_then(|o| *o) != Some(*b)
                    });
                if unique {
                    n += 1;
                }
            }
            n
        }

        /// The set of alleles present at column `ci` (reference ∪ all copies ∪ all reads), sorted.
        fn alleles_at(&self, ci: usize) -> BTreeSet<u8> {
            let mut set = BTreeSet::new();
            if let Some(b) = self.columns[ci].ref_allele {
                set.insert(b);
            }
            for c in &self.copies {
                if let Some(Some(b)) = c.alleles.get(ci) {
                    set.insert(*b);
                }
            }
            for r in &self.reads {
                if let Some(Some(b)) = r.obs.get(ci) {
                    set.insert(*b);
                }
            }
            set
        }

        /// Walk string "bb0+,c0_x+,bb1+,...,bbM+" given the allele taken at each column (`taken[ci]`).
        /// A `None` in `taken` routes through the reference allele node at that column.
        fn walk_tokens(&self, taken: &[Option<u8>]) -> Vec<String> {
            let m = self.m();
            let mut toks = Vec::with_capacity(2 * m + 1);
            for ci in 0..m {
                toks.push(format!("{}+", self.bb(ci)));
                let b = taken
                    .get(ci)
                    .and_then(|o| *o)
                    .or(self.columns[ci].ref_allele);
                if let Some(b) = b {
                    toks.push(format!("{}+", self.allele_node(ci, b)));
                }
            }
            toks.push(format!("{}+", self.bb(m)));
            toks
        }

        pub fn gfa_lines(&self) -> GfaLines {
            let mut out = GfaLines {
                header: "H\tVN:Z:1.1".into(),
                ..Default::default()
            };
            let m = self.m();
            // backbone spacer S-nodes bb0..=bbM
            for i in 0..=m {
                let seq = String::from_utf8_lossy(&self.backbone[i]).to_string();
                out.segs
                    .push(format!("S\t{}\t{}\tSN:Z:spacer", self.bb(i), seq));
            }
            // allele S-nodes + L-lines bb{ci} -> allele -> bb{ci+1}
            for ci in 0..m {
                let pos = self.columns[ci].genome_pos.unwrap_or(0);
                for b in self.alleles_at(ci) {
                    let nid = self.allele_node(ci, b);
                    out.segs
                        .push(format!("S\t{}\t{}\tPO:i:{}", nid, b as char, pos));
                    out.links
                        .push(format!("L\t{}\t+\t{}\t+\t0M", self.bb(ci), nid));
                    out.links
                        .push(format!("L\t{}\t+\t{}\t+\t0M", nid, self.bb(ci + 1)));
                }
            }
            // REFERENCE walk: the genome's own allele at each column.
            let ref_taken: Vec<Option<u8>> = self.columns.iter().map(|c| c.ref_allele).collect();
            let ref_walk = self.walk_tokens(&ref_taken).join(",");
            out.paths.push(format!(
                "P\t{}_REFERENCE\t{}\t*\tST:Z:reference",
                self.family, ref_walk
            ));

            // copy P-lines with corroboration tags
            for (copy_idx, cp) in self.copies.iter().enumerate() {
                let walk = self.walk_tokens(&cp.alleles).join(",");
                let name = if cp.status.is_absent() {
                    format!("{}_copy{}_ABSENT", self.family, copy_idx)
                } else {
                    format!("{}_copy{}", self.family, copy_idx)
                };
                let mut tags = String::new();
                if let Some(rc) = cp.corrob.reads {
                    tags.push_str(&format!("\tRC:i:{}", rc));
                }
                let su = cp
                    .corrob
                    .suns
                    .unwrap_or_else(|| self.private_columns(copy_idx));
                tags.push_str(&format!("\tSU:i:{}", su));
                if let Some(mi) = cp.corrob.map_identity {
                    tags.push_str(&format!("\tMI:f:{:.3}", mi));
                }
                tags.push_str(&format!("\tST:Z:{}", cp.status.tag()));
                out.paths.push(format!("P\t{}\t{}\t*{}", name, walk, tags));
            }

            // read W-lines over each read's observed span; gaps within span route through the reference node.
            for r in &self.reads {
                let first = r.obs.iter().position(|o| o.is_some());
                let last = r.obs.iter().rposition(|o| o.is_some());
                let (Some(first), Some(last)) = (first, last) else {
                    continue;
                };
                let mut toks: Vec<String> = Vec::new();
                for ci in first..=last {
                    toks.push(format!(">{}", self.bb(ci)));
                    let b = r.obs[ci].or(self.columns[ci].ref_allele);
                    if let Some(b) = b {
                        toks.push(format!(">{}", self.allele_node(ci, b)));
                    }
                }
                toks.push(format!(">{}", self.bb(last + 1)));
                let w = toks.join("");
                let hap = r.assigned_copy.map(|c| c as i64).unwrap_or(-1).max(0);
                let mut line = format!(
                    "W\t{}\t{}\t{}\t0\t{}\t{}",
                    r.name,
                    hap,
                    self.family,
                    toks.len(),
                    w
                );
                if let Some(c) = &r.cert {
                    use crate::family::copy_assign::AssignStatus::*;
                    let st = match c.status {
                        Assigned => "Assigned",
                        Ambiguous => "Ambiguous",
                        Tied => "Tied",
                    };
                    let cp = r
                        .assigned_copy
                        .map(|c| format!("copy{c}"))
                        .unwrap_or_else(|| "none".into());
                    line.push_str(&format!(
                        "\tCP:Z:{cp}\tPV:f:{}\tMP:f:{}\tST:Z:{st}",
                        c.p_value, c.min_p_value
                    ));
                }
                out.walks.push(line);
            }

            out
        }

        /// One self-contained GFA string (header + this family's lines). Convenience for tests / single-family use.
        pub fn to_gfa(&self) -> String {
            let g = self.gfa_lines();
            let mut s = String::new();
            s.push_str(&g.header);
            s.push('\n');
            for l in g
                .segs
                .iter()
                .chain(g.links.iter())
                .chain(g.paths.iter())
                .chain(g.walks.iter())
            {
                s.push_str(l);
                s.push('\n');
            }
            s
        }

        /// Bandage node colours (keyed on SEGMENT names): reference-walk nodes grey, absent-only divergent
        /// nodes red, other copy-divergent nodes their copy's status colour, backbone light grey.
        pub fn colours_csv(&self) -> String {
            use std::collections::BTreeMap;
            let mut colour: BTreeMap<String, &'static str> = BTreeMap::new();
            // backbone
            for i in 0..=self.m() {
                colour.insert(self.bb(i), "#dadce0");
            }
            for ci in 0..self.m() {
                let refb = self.columns[ci].ref_allele;
                for b in self.alleles_at(ci) {
                    let nid = self.allele_node(ci, b);
                    if Some(b) == refb {
                        colour.insert(nid, CopyStatus::Reference.colour());
                        continue;
                    }
                    // walked by any absent copy? (and not the reference allele) — absent wins over non-absent.
                    if let Some(c) = self.copies.iter().find(|c| {
                        c.status.is_absent() && c.alleles.get(ci).and_then(|o| *o) == Some(b)
                    }) {
                        colour.insert(nid, c.status.colour());
                    } else if let Some(c) = self
                        .copies
                        .iter()
                        .find(|c| c.alleles.get(ci).and_then(|o| *o) == Some(b))
                    {
                        colour.insert(nid, c.status.colour());
                    } else {
                        // observed only in reads (no copy carries it, not the reference) — neutral read-only.
                        colour.insert(nid, READ_ONLY_COLOUR);
                    }
                }
            }
            let mut s = String::new();
            for (k, v) in colour {
                s.push_str(&format!("{},{}\n", k, v));
            }
            s
        }

        /// Legend: each status actually present (plus reference) → its colour.
        pub fn legend_tsv(&self) -> String {
            use std::collections::BTreeSet;
            let mut statuses: BTreeSet<&'static str> = BTreeSet::new();
            statuses.insert("reference");
            let mut rows: Vec<(&'static str, &'static str)> =
                vec![("reference", CopyStatus::Reference.colour())];
            for c in &self.copies {
                if statuses.insert(c.status.tag()) {
                    rows.push((c.status.tag(), c.status.colour()));
                }
            }
            let mut s = String::new();
            for (st, col) in rows {
                s.push_str(&format!("{}\t{}\n", st, col));
            }
            s
        }
    }

    impl ExonGraph {
        fn node(&self, k: usize) -> String {
            format!("{}_E{}", self.family, k)
        }

        /// Bandage node colours (keyed on exon node names): a class walked by ≥1 non-absent copy → grey;
        /// else if only absent copies walk it → the owner's colour (red). Skip classes with no walkers.
        pub fn colours_csv(&self) -> String {
            let n = self.nodes.len();
            let mut colour: BTreeMap<String, &'static str> = BTreeMap::new();
            for k in 0..n {
                let walkers: Vec<&CopyExonPath> = self
                    .copies
                    .iter()
                    .filter(|c| c.exon_nodes.contains(&k))
                    .collect();
                let on_ref = walkers.iter().any(|c| !c.status.is_absent());
                let col = if on_ref {
                    CopyStatus::Reference.colour() // grey shared/reference exon
                } else if let Some(c) = walkers.first() {
                    c.status.colour() // copy-specific arm -> owner colour (red for absent)
                } else {
                    continue;
                };
                colour.insert(self.node(k), col);
            }
            let mut s = String::new();
            for (kk, v) in colour {
                s.push_str(&format!("{},{}\n", kk, v));
            }
            s
        }

        /// Legend: each status actually present (plus reference) → its colour.
        pub fn legend_tsv(&self) -> String {
            let mut seen: BTreeSet<&'static str> = BTreeSet::new();
            let mut s = format!("reference\t{}\n", CopyStatus::Reference.colour());
            seen.insert("reference");
            for c in &self.copies {
                if seen.insert(c.status.tag()) {
                    s.push_str(&format!("{}\t{}\n", c.status.tag(), c.status.colour()));
                }
            }
            s
        }

        /// Reciprocal overlap = min(inter/len_a, inter/len_b); 0 if disjoint or different chrom.
        fn recip_overlap(a: (&str, u64, u64), b: (&str, u64, u64)) -> f64 {
            if a.0 != b.0 {
                return 0.0;
            }
            let lo = a.1.max(b.1);
            let hi = a.2.min(b.2);
            if hi <= lo {
                return 0.0;
            }
            let inter = (hi - lo) as f64;
            (inter / (a.2 - a.1) as f64).min(inter / (b.2 - b.1) as f64)
        }

        pub fn to_gfa(&self, exon_seq: impl Fn(&ExonClass) -> Vec<u8>) -> String {
            // reference = classes present in >=1 non-absent copy; fallback = present in all copies
            let n = self.nodes.len();
            let non_absent: Vec<&CopyExonPath> = self
                .copies
                .iter()
                .filter(|c| !c.status.is_absent())
                .collect();
            let mut ref_nodes: Vec<usize> = (0..n)
                .filter(|&k| {
                    if non_absent.is_empty() {
                        self.copies.iter().all(|c| c.exon_nodes.contains(&k))
                    } else {
                        non_absent.iter().any(|c| c.exon_nodes.contains(&k))
                    }
                })
                .collect();
            ref_nodes.sort();

            // per-class RC = sum of reads over copies walking it
            let rc = |k: usize| -> u32 {
                self.copies
                    .iter()
                    .filter(|c| c.exon_nodes.contains(&k))
                    .filter_map(|c| c.corrob.reads)
                    .sum()
            };

            let mut s = String::from("H\tVN:Z:1.1\n");
            for (k, ec) in self.nodes.iter().enumerate() {
                let seq = String::from_utf8_lossy(&exon_seq(ec)).to_string();
                s.push_str(&format!(
                    "S\t{}\t{}\tPO:i:{}\tRC:i:{}\n",
                    self.node(k),
                    seq,
                    ec.start,
                    rc(k)
                ));
            }
            // L-lines = union of consecutive adjacencies across reference + all copy walks
            use std::collections::BTreeSet;
            let mut links: BTreeSet<(usize, usize)> = BTreeSet::new();
            let mut add_walk = |walk: &[usize], set: &mut BTreeSet<(usize, usize)>| {
                for w in walk.windows(2) {
                    set.insert((w[0], w[1]));
                }
            };
            add_walk(&ref_nodes, &mut links);
            for c in &self.copies {
                add_walk(&c.exon_nodes, &mut links);
            }
            for (a, b) in &links {
                s.push_str(&format!(
                    "L\t{}\t+\t{}\t+\t0M\n",
                    self.node(*a),
                    self.node(*b)
                ));
            }
            // REFERENCE P-line
            let rwalk: String = ref_nodes
                .iter()
                .map(|k| format!("{}+", self.node(*k)))
                .collect::<Vec<_>>()
                .join(",");
            s.push_str(&format!(
                "P\t{}_REFERENCE\t{}\t*\tST:Z:reference\n",
                self.family, rwalk
            ));
            // copy P-lines
            for (ci, c) in self.copies.iter().enumerate() {
                let name = if c.status.is_absent() {
                    format!("{}_copy{}_ABSENT", self.family, ci)
                } else {
                    format!("{}_copy{}", self.family, ci)
                };
                let walk: String = c
                    .exon_nodes
                    .iter()
                    .map(|k| format!("{}+", self.node(*k)))
                    .collect::<Vec<_>>()
                    .join(",");
                let mut tags = String::new();
                if let Some(r) = c.corrob.reads {
                    tags.push_str(&format!("\tRC:i:{}", r));
                }
                if let Some(mi) = c.corrob.map_identity {
                    tags.push_str(&format!("\tMI:f:{:.3}", mi));
                }
                tags.push_str(&format!("\tST:Z:{}", c.status.tag()));
                s.push_str(&format!("P\t{}\t{}\t*{}\n", name, walk, tags));
            }
            s
        }

        pub fn from_copies(
            family: &str,
            copies: &[(String, CopyStatus, Corrob, String, Vec<(u64, u64)>)],
        ) -> ExonGraph {
            // flatten every exon as (copy_idx, chrom, start, end); skip malformed/zero-length (end <= start)
            // so the `(a.2 - a.1)` length divisor never underflows.
            let mut flat: Vec<(usize, String, u64, u64)> = Vec::new();
            for (ci, (_, _, _, chrom, exons)) in copies.iter().enumerate() {
                for &(s, e) in exons {
                    if e > s {
                        flat.push((ci, chrom.clone(), s, e));
                    }
                }
            }

            // union-find over flat exon items
            let n = flat.len();
            let mut parent: Vec<usize> = (0..n).collect();

            fn find(p: &mut Vec<usize>, x: usize) -> usize {
                if p[x] != x {
                    let r = find(p, p[x]);
                    p[x] = r;
                }
                p[x]
            }

            for i in 0..n {
                for j in (i + 1)..n {
                    if Self::recip_overlap(
                        (&flat[i].1, flat[i].2, flat[i].3),
                        (&flat[j].1, flat[j].2, flat[j].3),
                    ) >= 0.30
                    {
                        let (a, b) = (find(&mut parent, i), find(&mut parent, j));
                        if a != b {
                            parent[a] = b;
                        }
                    }
                }
            }

            // group items by root, build one ExonClass per group (min start, max end)
            let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for i in 0..n {
                let r = find(&mut parent, i);
                groups.entry(r).or_default().push(i);
            }

            let mut classes: Vec<(String, u64, u64, Vec<usize>)> = groups
                .values()
                .map(|idxs| {
                    let start = idxs.iter().map(|&i| flat[i].2).min().unwrap();
                    let end = idxs.iter().map(|&i| flat[i].3).max().unwrap();
                    let members: Vec<usize> = idxs.iter().map(|&i| flat[i].0).collect();
                    (flat[idxs[0]].1.clone(), start, end, members)
                })
                .collect();

            classes.sort_by_key(|c| (c.1, c.2)); // by genomic (start, end) -> canonical E0..En order

            let nodes: Vec<ExonClass> = classes
                .iter()
                .map(|(chrom, s, e, _)| ExonClass {
                    chrom: chrom.clone(),
                    start: *s,
                    end: *e,
                })
                .collect();

            let copy_paths: Vec<CopyExonPath> = copies
                .iter()
                .enumerate()
                .map(|(ci, (id, status, corrob, _chrom, _exons))| {
                    let mut exon_nodes: Vec<usize> = (0..classes.len())
                        .filter(|&k| classes[k].3.contains(&ci))
                        .collect();
                    exon_nodes.sort();
                    CopyExonPath {
                        id: id.clone(),
                        exon_nodes,
                        status: *status,
                        corrob: corrob.clone(),
                    }
                })
                .collect();

            ExonGraph {
                family: family.to_string(),
                nodes,
                copies: copy_paths,
            }
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn tiny_graph() -> CopyGraph {
            // 2 columns, backbone spacers of len 3, one reference + one copy
            CopyGraph {
                family: "FAM1".into(),
                columns: vec![
                    PsvColumn {
                        col: 0,
                        genome_pos: Some(100),
                        ref_allele: Some(b'A'),
                    },
                    PsvColumn {
                        col: 1,
                        genome_pos: Some(200),
                        ref_allele: Some(b'C'),
                    },
                ],
                backbone: vec![b"NNN".to_vec(), b"NNN".to_vec(), b"NNN".to_vec()],
                copies: vec![CopyPath {
                    id: "FAM1_copy0".into(),
                    alleles: vec![Some(b'A'), Some(b'G')],
                    status: CopyStatus::InGenomeAnnotated,
                    corrob: Corrob {
                        reads: Some(5),
                        suns: None,
                        map_identity: Some(0.99),
                    },
                }],
                reads: vec![],
            }
        }

        #[test]
        fn constructs_and_reports_shape() {
            let g = tiny_graph();
            assert_eq!(g.columns.len(), 2);
            assert_eq!(g.backbone.len(), 3);
            assert_eq!(g.copies[0].status.tag(), "in-genome-annotated");
            assert!(!g.copies[0].status.is_absent());
            assert!(CopyStatus::AbsentCollapsed.is_absent());
        }

        fn parse_line_prefixes(gfa: &str) -> (usize, usize, usize) {
            let (mut s, mut l, mut h) = (0, 0, 0);
            for line in gfa.lines() {
                match line.chars().next() {
                    Some('S') => s += 1,
                    Some('L') => l += 1,
                    Some('H') => h += 1,
                    _ => {}
                }
            }
            (h, s, l)
        }

        #[test]
        fn skeleton_has_backbone_alleles_and_links() {
            let g = tiny_graph(); // col0 alleles {A(ref), A(copy)} => {A}; col1 alleles {C(ref), G(copy)} => {C,G}
            let gfa = g.to_gfa();
            let (h, s, _l) = parse_line_prefixes(&gfa);
            assert_eq!(h, 1, "one header");
            // backbone: bb0,bb1,bb2 (3) + allele nodes: col0 {A}=1, col1 {C,G}=2 => 3 => total S = 6
            assert_eq!(s, 6, "3 backbone + 3 allele S-nodes");
            assert!(gfa.contains("S\tFAM1_bb0\tNNN"));
            assert!(gfa.contains("S\tFAM1_c0_A\tA\tPO:i:100"));
            assert!(gfa.contains("S\tFAM1_c1_G\tG\tPO:i:200"));
            // every allele node linked to its flanking backbone (no dangling by construction)
            assert!(gfa.contains("L\tFAM1_bb0\t+\tFAM1_c0_A\t+\t0M"));
            assert!(gfa.contains("L\tFAM1_c0_A\t+\tFAM1_bb1\t+\t0M"));
            assert!(gfa.contains("L\tFAM1_bb1\t+\tFAM1_c1_G\t+\t0M"));
        }

        #[test]
        fn reference_walk_threads_reference_alleles() {
            let g = tiny_graph();
            let gfa = g.to_gfa();
            // reference alleles are A (col0) and C (col1)
            assert!(gfa.contains(
                "P\tFAM1_REFERENCE\tFAM1_bb0+,FAM1_c0_A+,FAM1_bb1+,FAM1_c1_C+,FAM1_bb2+\t*\tST:Z:reference"
            ), "reference P-line missing or wrong:\n{}", gfa);
        }

        #[test]
        fn copy_paths_carry_tags_and_absent_diverges() {
            // 3 columns; reference = A,A,A. copy0 in-genome matches ref. copy1 ABSENT diverges at col1 & col2.
            let g = CopyGraph {
                family: "FAM2".into(),
                columns: (0..3)
                    .map(|i| PsvColumn {
                        col: i,
                        genome_pos: Some(100 + i as u64),
                        ref_allele: Some(b'A'),
                    })
                    .collect(),
                backbone: vec![b"NN".to_vec(); 4],
                copies: vec![
                    CopyPath {
                        id: "FAM2_copy0".into(),
                        alleles: vec![Some(b'A'), Some(b'A'), Some(b'A')],
                        status: CopyStatus::InGenomeAnnotated,
                        corrob: Corrob {
                            reads: Some(8),
                            suns: None,
                            map_identity: Some(0.998),
                        },
                    },
                    CopyPath {
                        id: "FAM2_copy1".into(),
                        alleles: vec![Some(b'A'), Some(b'G'), Some(b'T')],
                        status: CopyStatus::AbsentDivergent,
                        corrob: Corrob {
                            reads: Some(12),
                            suns: None,
                            map_identity: Some(0.952),
                        },
                    },
                ],
                reads: vec![],
            };
            let gfa = g.to_gfa();
            // absent copy P-line named with _ABSENT and tagged
            let absent = gfa
                .lines()
                .find(|l| l.starts_with("P\tFAM2_copy1_ABSENT"))
                .expect("absent P-line");
            assert!(absent.contains("RC:i:12"));
            assert!(absent.contains("MI:f:0.952"));
            assert!(absent.contains("ST:Z:absent-divergent"));
            // SU (private columns): copy1's allele is unique (vs copy0) at col1(G) and col2(T) => SU:i:2
            assert!(absent.contains("SU:i:2"), "expected SU:i:2 in: {}", absent);
            // it walks the divergent nodes the reference walk does NOT (c1_G, c2_T)
            assert!(absent.contains("FAM2_c1_G+"));
            assert!(absent.contains("FAM2_c2_T+"));
            // in-genome copy0 has SU:i:0 (never unique) and is not _ABSENT
            let c0 = gfa
                .lines()
                .find(|l| l.starts_with("P\tFAM2_copy0\t"))
                .expect("copy0 P-line");
            assert!(c0.contains("SU:i:0"));
            assert!(c0.contains("ST:Z:in-genome-annotated"));
        }

        #[test]
        fn omits_unknown_corrob_tags() {
            // Honesty rule NEGATIVE path: when reads/map_identity are None, RC:i: and MI:f: are OMITTED,
            // while SU:i: (always computed) and ST:Z: (always emitted) remain.
            let g = CopyGraph {
                family: "FAM3".into(),
                columns: vec![
                    PsvColumn {
                        col: 0,
                        genome_pos: Some(100),
                        ref_allele: Some(b'A'),
                    },
                    PsvColumn {
                        col: 1,
                        genome_pos: Some(200),
                        ref_allele: Some(b'C'),
                    },
                ],
                backbone: vec![b"NN".to_vec(); 3],
                copies: vec![CopyPath {
                    id: "FAM3_copy0".into(),
                    alleles: vec![Some(b'A'), Some(b'G')],
                    status: CopyStatus::AnnotationUnknown,
                    corrob: Corrob {
                        reads: None,
                        suns: None,
                        map_identity: None,
                    },
                }],
                reads: vec![],
            };
            let gfa = g.to_gfa();
            let cp = gfa
                .lines()
                .find(|l| l.starts_with("P\tFAM3_copy0\t"))
                .expect("copy0 P-line");
            // unknown values => tags omitted (never faked)
            assert!(
                !cp.contains("RC:i:"),
                "RC:i: must be omitted when reads is None: {}",
                cp
            );
            assert!(
                !cp.contains("MI:f:"),
                "MI:f: must be omitted when map_identity is None: {}",
                cp
            );
            // always-present tags
            assert!(cp.contains("SU:i:"), "SU:i: must always be present: {}", cp);
            assert!(
                cp.contains("ST:Z:annotation-unknown"),
                "ST:Z: must always be present: {}",
                cp
            );
        }

        // Assert every P-line and W-line step is backed by an L-line (parses walks, checks adjacency set).
        fn assert_no_dangling(gfa: &str) {
            use std::collections::HashSet;
            let mut links: HashSet<(String, String)> = HashSet::new();
            for l in gfa.lines().filter(|l| l.starts_with("L\t")) {
                let f: Vec<&str> = l.split('\t').collect(); // L from + to + 0M
                links.insert((f[1].to_string(), f[3].to_string()));
            }
            let node = |tok: &str| {
                tok.trim_start_matches(['>', '<'])
                    .trim_end_matches(['>', '<', '+', '-'])
                    .to_string()
            };
            for l in gfa.lines() {
                let seq: Vec<String> = if l.starts_with("P\t") {
                    l.split('\t').nth(2).unwrap().split(',').map(node).collect()
                } else if l.starts_with("W\t") {
                    let w = l.split('\t').nth(6).unwrap();
                    w.split_inclusive(['>', '<'])
                        .filter(|s| s.len() > 1)
                        .map(node)
                        .collect()
                } else {
                    continue;
                };
                for pair in seq.windows(2) {
                    assert!(
                        links.contains(&(pair[0].clone(), pair[1].clone())),
                        "dangling walk edge {}->{} in line: {}",
                        pair[0],
                        pair[1],
                        l
                    );
                }
            }
        }

        #[test]
        fn reads_walk_with_backing_links() {
            let mut g = tiny_graph(); // 2 cols, ref A,C
            g.reads = vec![
                ReadWalk {
                    name: "readX".into(),
                    obs: vec![Some(b'A'), Some(b'C')],
                    assigned_copy: Some(0),
                    cert: None,
                },
                ReadWalk {
                    name: "readY".into(),
                    obs: vec![None, Some(b'C')],
                    assigned_copy: None,
                    cert: None,
                },
            ];
            let gfa = g.to_gfa();
            assert!(
                gfa.lines().any(|l| l.starts_with("W\treadX")),
                "readX walk missing"
            );
            assert!(
                gfa.lines().any(|l| l.starts_with("W\treadY")),
                "readY walk missing"
            );
            assert_no_dangling(&gfa);
        }

        #[test]
        fn read_walk_cert_tags_emitted_when_present() {
            use crate::family::copy_assign::AssignStatus;
            let mut g = tiny_graph(); // 2 cols, ref A,C; copy0 alleles A,G
            g.reads = vec![
                ReadWalk {
                    name: "r1".into(),
                    obs: vec![Some(b'A'), Some(b'C')],
                    assigned_copy: Some(0),
                    cert: Some(ReadCert {
                        p_value: 0.001,
                        min_p_value: 0.0005,
                        status: AssignStatus::Assigned,
                    }),
                },
                ReadWalk {
                    name: "r2".into(),
                    obs: vec![Some(b'A'), Some(b'C')],
                    assigned_copy: None,
                    cert: None,
                },
            ];
            let gfa = g.to_gfa();
            let r1 = gfa
                .lines()
                .find(|l| l.starts_with("W\tr1"))
                .expect("r1 walk");
            assert!(
                r1.contains("CP:Z:copy0")
                    && r1.contains("PV:f:0.001")
                    && r1.contains("MP:f:0.0005")
                    && r1.contains("ST:Z:Assigned")
            );
            let r2 = gfa
                .lines()
                .find(|l| l.starts_with("W\tr2"))
                .expect("r2 walk");
            assert!(
                !r2.contains("CP:Z") && !r2.contains("PV:f"),
                "no cert -> no tags (backward-compatible)"
            );
        }

        #[test]
        fn read_walk_internal_gap_routes_through_reference() {
            // 3 columns, reference A,A,A. A read observes col0 (A) and col2 (T) but NOT col1 (internal gap):
            // the gap at col1 must route through the reference allele node c1_A. The walk spans
            // bb(first_obs=0) .. bb(last_obs+1=3).
            let g = CopyGraph {
                family: "FAM4".into(),
                columns: (0..3)
                    .map(|i| PsvColumn {
                        col: i,
                        genome_pos: Some(100 + i as u64),
                        ref_allele: Some(b'A'),
                    })
                    .collect(),
                backbone: vec![b"NN".to_vec(); 4],
                copies: vec![],
                reads: vec![ReadWalk {
                    name: "readG".into(),
                    obs: vec![Some(b'A'), None, Some(b'T')],
                    assigned_copy: None,
                    cert: None,
                }],
            };
            let gfa = g.to_gfa();
            // exact W-line: gap col1 routes through reference node c1_A; starts bb0, ends bb3; 7 tokens.
            assert!(gfa.lines().any(|l| l ==
                "W\treadG\t0\tFAM4\t0\t7\t>FAM4_bb0>FAM4_c0_A>FAM4_bb1>FAM4_c1_A>FAM4_bb2>FAM4_c2_T>FAM4_bb3"
            ), "internal-gap W-line missing or wrong:\n{}", gfa);
            assert_no_dangling(&gfa);
        }

        #[test]
        fn read_observing_zero_columns_emits_no_walk() {
            // A read with no observations (all None) must emit NO W-line (early continue).
            let mut g = tiny_graph(); // 2 cols
            g.reads = vec![ReadWalk {
                name: "readEmpty".into(),
                obs: vec![None, None],
                assigned_copy: None,
                cert: None,
            }];
            let gfa = g.to_gfa();
            assert!(
                !gfa.lines().any(|l| l.starts_with("W\t")),
                "no W-line expected for a read with zero observations:\n{}",
                gfa
            );
        }

        #[test]
        fn colours_mark_absent_red_reference_grey() {
            // reuse the 3-column absent-copy graph
            let g = CopyGraph {
                family: "FAM3".into(),
                columns: (0..3)
                    .map(|i| PsvColumn {
                        col: i,
                        genome_pos: Some(10 + i as u64),
                        ref_allele: Some(b'A'),
                    })
                    .collect(),
                backbone: vec![b"NN".to_vec(); 4],
                copies: vec![
                    CopyPath {
                        id: "FAM3_copy0".into(),
                        alleles: vec![Some(b'A'), Some(b'A'), Some(b'A')],
                        status: CopyStatus::InGenomeAnnotated,
                        corrob: Corrob::default(),
                    },
                    CopyPath {
                        id: "FAM3_copy1".into(),
                        alleles: vec![Some(b'A'), Some(b'G'), Some(b'T')],
                        status: CopyStatus::AbsentDivergent,
                        corrob: Corrob::default(),
                    },
                ],
                reads: vec![],
            };
            let csv = g.colours_csv();
            // reference allele node grey
            assert!(
                csv.lines().any(|l| l == "FAM3_c0_A,#9aa0a6"),
                "ref node not grey:\n{}",
                csv
            );
            // absent-only divergent nodes red
            assert!(
                csv.lines().any(|l| l == "FAM3_c1_G,#d93025"),
                "absent node not red:\n{}",
                csv
            );
            assert!(csv.lines().any(|l| l == "FAM3_c2_T,#d93025"));
            // legend lists the two statuses in use
            let legend = g.legend_tsv();
            assert!(legend.contains("reference\t#9aa0a6"));
            assert!(legend.contains("absent-divergent\t#d93025"));
        }

        #[test]
        fn colours_read_only_allele_gets_neutral() {
            // A base observed ONLY in a read (differs from ref_allele AND carried by no CopyPath) still
            // gets a GFA allele segment via alleles_at — it must receive the neutral read-only colour, not
            // be silently dropped from colours.csv (which would render uncoloured in Bandage).
            let g = CopyGraph {
                family: "FAM5".into(),
                columns: vec![PsvColumn {
                    col: 0,
                    genome_pos: Some(100),
                    ref_allele: Some(b'A'),
                }],
                backbone: vec![b"NN".to_vec(); 2],
                copies: vec![CopyPath {
                    id: "FAM5_copy0".into(),
                    alleles: vec![Some(b'A')],
                    status: CopyStatus::InGenomeAnnotated,
                    corrob: Corrob::default(),
                }],
                // read observes 'G' at col0 — neither the reference (A) nor any copy (A) carries it.
                reads: vec![ReadWalk {
                    name: "readR".into(),
                    obs: vec![Some(b'G')],
                    assigned_copy: None,
                    cert: None,
                }],
            };
            let csv = g.colours_csv();
            // the read-only node exists in the GFA (alleles_at folds it in)…
            assert!(
                g.to_gfa().contains("S\tFAM5_c0_G\tG"),
                "read-only allele node missing from GFA:\n{}",
                g.to_gfa()
            );
            // …and it must be deliberately coloured neutral, distinct from the backbone light-grey.
            assert!(
                csv.lines().any(|l| l == "FAM5_c0_G,#e8eaed"),
                "read-only node not neutral:\n{}",
                csv
            );
            // reference allele still grey
            assert!(
                csv.lines().any(|l| l == "FAM5_c0_A,#9aa0a6"),
                "ref node not grey:\n{}",
                csv
            );
        }

        #[test]
        fn colours_absent_wins_over_non_absent_at_shared_node() {
            // A single non-reference allele node walked by BOTH an absent copy and a non-absent copy at the
            // same column/base must come out RED (absent precedence), not the non-absent status colour.
            let g = CopyGraph {
                family: "FAM6".into(),
                columns: vec![PsvColumn {
                    col: 0,
                    genome_pos: Some(100),
                    ref_allele: Some(b'A'),
                }],
                backbone: vec![b"NN".to_vec(); 2],
                copies: vec![
                    // non-absent copy walks G at col0…
                    CopyPath {
                        id: "FAM6_copy0".into(),
                        alleles: vec![Some(b'G')],
                        status: CopyStatus::InGenomeAnnotated,
                        corrob: Corrob::default(),
                    },
                    // …and an absent copy walks the SAME G at col0.
                    CopyPath {
                        id: "FAM6_copy1".into(),
                        alleles: vec![Some(b'G')],
                        status: CopyStatus::AbsentDivergent,
                        corrob: Corrob::default(),
                    },
                ],
                reads: vec![],
            };
            let csv = g.colours_csv();
            // absent wins: the shared node is red, NOT the in-genome blue (#1a73e8).
            assert!(
                csv.lines().any(|l| l == "FAM6_c0_G,#d93025"),
                "shared absent/non-absent node must be red (absent wins):\n{}",
                csv
            );
            assert!(
                !csv.lines().any(|l| l == "FAM6_c0_G,#1a73e8"),
                "shared node must NOT take the non-absent colour:\n{}",
                csv
            );
        }

        #[test]
        fn exon_graph_constructs() {
            let g = ExonGraph {
                family: "F".into(),
                nodes: vec![ExonClass {
                    chrom: "c".into(),
                    start: 0,
                    end: 100,
                }],
                copies: vec![CopyExonPath {
                    id: "F_copy0".into(),
                    exon_nodes: vec![0],
                    status: CopyStatus::InGenomeAnnotated,
                    corrob: Corrob::default(),
                }],
            };
            assert_eq!(g.nodes.len(), 1);
            assert_eq!(g.copies[0].exon_nodes, vec![0]);
        }

        #[test]
        fn from_copies_clusters_and_flags_copy_specific_exon() {
            // copy0 exons E1,E3 ; copy1 exons E1,E2(extra),E3 — E2 is copy1-specific.
            let copies = vec![
                (
                    "F_copy0".to_string(),
                    CopyStatus::InGenomeAnnotated,
                    Corrob::default(),
                    "chr1".to_string(),
                    vec![(0u64, 100u64), (300, 400)],
                ),
                (
                    "F_copy1".to_string(),
                    CopyStatus::AbsentDivergent,
                    Corrob::default(),
                    "chr1".to_string(),
                    vec![(0, 100), (150, 250), (300, 400)],
                ),
            ];
            let g = ExonGraph::from_copies("F", &copies);
            assert_eq!(g.nodes.len(), 3, "E1,E2,E3");
            // find the class only copy1 walks (the extra exon ~150-250)
            let owners: Vec<Vec<usize>> = (0..g.nodes.len())
                .map(|k| {
                    g.copies
                        .iter()
                        .enumerate()
                        .filter(|(_, c)| c.exon_nodes.contains(&k))
                        .map(|(i, _)| i)
                        .collect()
                })
                .collect();
            let copy_specific: Vec<usize> = (0..g.nodes.len())
                .filter(|&k| owners[k] == vec![1])
                .collect();
            assert_eq!(copy_specific.len(), 1, "exactly one copy1-specific exon");
            // copy0 walks 2 classes, copy1 walks 3
            assert_eq!(g.copies[0].exon_nodes.len(), 2);
            assert_eq!(g.copies[1].exon_nodes.len(), 3);
        }

        #[test]
        fn from_copies_respects_overlap_threshold() {
            // A=(0,100), B=(70,170): inter=30 over len 100 => recip=0.30 => AT threshold => MERGE (1 class).
            let merge = vec![
                (
                    "F_copy0".to_string(),
                    CopyStatus::InGenomeAnnotated,
                    Corrob::default(),
                    "chr1".to_string(),
                    vec![(0u64, 100u64)],
                ),
                (
                    "F_copy1".to_string(),
                    CopyStatus::InGenomeAnnotated,
                    Corrob::default(),
                    "chr1".to_string(),
                    vec![(70u64, 170u64)],
                ),
            ];
            let g = ExonGraph::from_copies("F", &merge);
            assert_eq!(
                g.nodes.len(),
                1,
                "recip overlap exactly 0.30 merges into one class"
            );
            assert_eq!(g.copies[0].exon_nodes, vec![0]);
            assert_eq!(g.copies[1].exon_nodes, vec![0]);

            // A=(0,100), B'=(71,171): inter=29 over len 100 => recip=0.29 => below threshold => SEPARATE (2 classes).
            let split = vec![
                (
                    "F_copy0".to_string(),
                    CopyStatus::InGenomeAnnotated,
                    Corrob::default(),
                    "chr1".to_string(),
                    vec![(0u64, 100u64)],
                ),
                (
                    "F_copy1".to_string(),
                    CopyStatus::InGenomeAnnotated,
                    Corrob::default(),
                    "chr1".to_string(),
                    vec![(71u64, 171u64)],
                ),
            ];
            let g = ExonGraph::from_copies("F", &split);
            assert_eq!(
                g.nodes.len(),
                2,
                "recip overlap 0.29 (< 0.30) stays two separate classes"
            );
            // classes sorted by start: E0=(0,100) is copy0's, E1=(71,171) is copy1's.
            assert_eq!(g.copies[0].exon_nodes, vec![0]);
            assert_eq!(g.copies[1].exon_nodes, vec![1]);
        }

        #[test]
        fn exon_gfa_has_reference_skip_and_arm_no_dangling() {
            let copies = vec![
                (
                    "F_copy0".to_string(),
                    CopyStatus::InGenomeAnnotated,
                    Corrob {
                        reads: Some(10),
                        suns: None,
                        map_identity: None,
                    },
                    "chr1".to_string(),
                    vec![(0u64, 100u64), (300, 400)],
                ),
                (
                    "F_copy1".to_string(),
                    CopyStatus::AbsentDivergent,
                    Corrob {
                        reads: Some(5),
                        suns: None,
                        map_identity: Some(0.95),
                    },
                    "chr1".to_string(),
                    vec![(0, 100), (150, 250), (300, 400)],
                ),
            ];
            let g = ExonGraph::from_copies("F", &copies);
            let gfa = g.to_gfa(|ec| vec![b'A'; (ec.end - ec.start) as usize]);
            // reference exists and is the shared backbone (2 classes), copy1 absent walks 3
            assert!(gfa.contains("P\tF_REFERENCE"));
            let c1 = gfa
                .lines()
                .find(|l| l.starts_with("P\tF_copy1_ABSENT"))
                .unwrap();
            assert!(c1.contains("MI:f:0.950"));
            assert!(c1.contains("ST:Z:absent-divergent"));
            // the copy1-specific exon node exists with RC:i:5 (only copy1, 5 reads)
            let arm = g.copies[1]
                .exon_nodes
                .iter()
                .find(|&&k| !g.copies[0].exon_nodes.contains(&k))
                .copied()
                .unwrap();
            assert!(gfa.contains(&format!("RC:i:5")));
            assert!(gfa
                .lines()
                .any(|l| l.starts_with(&format!("S\tF_E{}", arm))));
            // no dangling: every P-line step is backed by an L-line
            assert_no_dangling(&gfa);
        }

        #[test]
        fn exon_colours_arm_red_shared_grey() {
            let copies = vec![
                (
                    "F_copy0".to_string(),
                    CopyStatus::InGenomeAnnotated,
                    Corrob::default(),
                    "chr1".to_string(),
                    vec![(0u64, 100u64), (300, 400)],
                ),
                (
                    "F_copy1".to_string(),
                    CopyStatus::AbsentDivergent,
                    Corrob::default(),
                    "chr1".to_string(),
                    vec![(0, 100), (150, 250), (300, 400)],
                ),
            ];
            let g = ExonGraph::from_copies("F", &copies);
            let csv = g.colours_csv();
            let arm = g.copies[1]
                .exon_nodes
                .iter()
                .find(|&&k| !g.copies[0].exon_nodes.contains(&k))
                .copied()
                .unwrap();
            assert!(
                csv.lines().any(|l| l == format!("F_E{},#d93025", arm)),
                "arm not red:\n{}",
                csv
            );
            // a shared class (walked by the in-genome copy0) is grey
            let shared = g.copies[0].exon_nodes[0];
            assert!(csv.lines().any(|l| l == format!("F_E{},#9aa0a6", shared)));
            assert!(g.legend_tsv().contains("absent-divergent\t#d93025"));
        }
    }

    // ---- merged 2026-10-05: was `vg_family/copy_discovery.rs`, now the inline module below (one component) ----
    #[allow(clippy::all)]
    pub mod copy_discovery {
        //! Discovery of candidate gene-family copies from read alignment ties.
        //!
        //! **STATUS:** OPT-IN
        //!
        //! Reached from `copy_assign` behind `--discover-copies` (`src/bin/copy_assign.rs`, `default_value_t =
        //! false`): `discover_copies_for_family` calls [`cluster_tie_partners`] once per family inside that
        //! flag's own `if args.discover_copies {` block. REPORT ONLY -- the result is written to
        //! `<out>.discovered_copies.tsv` and never mutates the input catalog or this run's own assignments.
        //!
        //! ⚠ A [`DiscoveredCopy`] row is a candidate ALIGNED-BLOCK CLUSTER, not necessarily a whole candidate
        //! copy: clustering is per-block (see [`aligned_blocks`]), so a genuine multi-exon copy supported by only
        //! a few reads can fragment into several rows, one per exon, each independently only needing its own
        //! `TIE_PARTNER_MIN_SUPPORT` reads -- more chances for a coincidental cluster than a per-placement scheme
        //! would give. Not detected or merged here (docs/o1_ledger.md §6l6 addendum, "Named limitation"); a
        //! follow-up would regroup blocks sharing a supporting placement into one candidate with its own
        //! `exon_blocks` column, mirroring `catalog_input::exon_blocks_str`.

        use crate::family::copy_split::AlignedRead;
        use crate::family::denovo_assemble::BamRead;
        use std::collections::{HashMap, HashSet};

        pub const TIE_PARTNER_MERGE_DISTANCE_BP: u64 = 500;
        pub const TIE_PARTNER_MIN_SUPPORT: usize = 2;

        #[derive(Clone, Debug, PartialEq)]
        pub struct DiscoveredCopy {
            pub family_id: String,
            pub chrom: String,
            pub start: u64,
            pub end: u64,
            /// MAJORITY strand across the placements supporting this cluster (`false` -> `+`, `true` -> `-` from
            /// each record's own SAM FLAG 0x10), an exact tie defaulting to `+` -- the rule the design doc's own
            /// "Open Question" section resolved, matching `majority_read_strand`'s existing convention
            /// (`denovo_pipeline.rs`) and `build_footprint_seq`'s documented `+`-on-tie placeholder caveat.
            pub strand: char,
            pub n_supporting_reads: usize,
            pub read_names: Vec<String>,
            pub nearest_copy_tid: String,
            /// `None` when this family has NO copy on the candidate's own chromosome (reachable for a
            /// cross-chromosome family) -- written as `NA`, never as a `u64::MAX` sentinel.
            pub nearest_copy_distance: Option<u64>,
        }

        /// One max-AS placement of one AS-tied read, carried as ALIGNED BLOCKS rather than a bounding span.
        ///
        /// ⚠ The blocks, not a `(ref_start, ref_end)` pair, are the load-bearing part. A bounding span counts a
        /// spliced-out intron (`N`) as if the read covered it, which made `inside_any_copy` call a placement
        /// "inside" a catalog copy that only an intron spans (no aligned base ever landing there) and inflated
        /// reported cluster widths by whole intron lengths. This is the same bug class `block_overlap`
        /// (`src/bin/copy_assign.rs`) was introduced to fix at the O3 truth gate.
        #[derive(Clone, Debug, PartialEq)]
        pub struct TiePlacement {
            pub chrom: String,
            /// One `(start, end)` per `M`/`=`/`X` run, in reference order. `D`/`N` advance the reference position
            /// without producing a block (exactly `block_overlap`'s own accounting).
            pub blocks: Vec<(u64, u64)>,
            /// SAM FLAG 0x10 of the record this placement came from -- the strand majority vote's input.
            pub reverse: bool,
        }

        /// Every `M`/`=`/`X` run of an alignment as its own reference interval `[pos, pos + n)`. The copy_assign
        /// binary's `block_overlap` (pysam `get_blocks()` semantics) is computed from these blocks.
        /// `D` and `N` advance the reference cursor without emitting a block, so a deletion splits a run here too
        /// (matching `block_overlap`) -- usually immaterial since alignment `D` runs are short, but a `D` longer
        /// than `TIE_PARTNER_MERGE_DISTANCE_BP` would split one placement into two reported clusters, same as a
        /// real gap between two placements. Not observed on the real substrate this module was built against.
        pub fn aligned_blocks(read: &AlignedRead) -> Vec<(u64, u64)> {
            let mut pos = read.ref_start;
            let mut out = Vec::new();
            for &(op, n) in &read.cigar {
                match op {
                    'M' | '=' | 'X' => {
                        out.push((pos, pos + n));
                        pos += n;
                    }
                    'D' | 'N' => pos += n,
                    _ => {}
                }
            }
            out
        }

        /// Does ANY aligned block of this placement overlap a copy already in the catalog?
        ///
        /// Per-block, not a bounding-box test: a read whose intron merely SPANS a catalog copy has no aligned
        /// base inside it and is NOT "inside" that copy.
        fn inside_any_copy(
            existing_copies: &[(String, u64, u64, String)],
            p: &TiePlacement,
        ) -> bool {
            existing_copies.iter().any(|(c_chrom, c_start, c_end, _)| {
                *c_chrom == p.chrom
                    && p.blocks
                        .iter()
                        .any(|(b_start, b_end)| b_start < c_end && b_end > c_start)
            })
        }

        /// Nearest catalog copy of this family ON THE SAME CHROMOSOME, and its distance in bp (`0` when the
        /// candidate overlaps it). `("NA", None)` when the family has no copy on that chromosome at all.
        fn nearest_copy(
            existing_copies: &[(String, u64, u64, String)],
            chrom: &str,
            start: u64,
            end: u64,
        ) -> (String, Option<u64>) {
            existing_copies
                .iter()
                .filter(|(c_chrom, ..)| c_chrom == chrom)
                .map(|(_, c_start, c_end, tid)| {
                    let d = if end <= *c_start {
                        c_start - end
                    } else if start >= *c_end {
                        start - c_end
                    } else {
                        0
                    };
                    (tid.clone(), Some(d))
                })
                .min_by_key(|(_, d)| *d)
                .unwrap_or(("NA".to_string(), None))
        }

        /// Cluster out-of-catalog AS-tie-partner positions into candidate copies for one family.
        ///
        /// `tied_reads`: `(read_name, that read's max-AS placements)` -- ALREADY RESTRICTED BY THE CALLER to the
        /// reads this family actually considered (`discover_copies_for_family`, `src/bin/copy_assign.rs`). Passing
        /// a whole region's tied reads to every family is the cross-family pooling bug the final whole-branch
        /// review caught: it attributed one identical site, with an identical read list, to 3-8 different
        /// `family_id`s.
        /// `existing_copies`: `(chrom, start, end, tid)` for every copy already catalogued in this family (the
        /// caller zips `FamilyAssignment::copy_spans` with `copy_tids`).
        /// Positions already inside an existing copy span are defensively re-excluded here (never surfaced), even
        /// though the caller is expected to have filtered them out already.
        ///
        /// Pure: no BAM, no catalog types, no I/O.
        pub fn cluster_tie_partners(
            tied_reads: &[(String, Vec<TiePlacement>)],
            family_id: &str,
            existing_copies: &[(String, u64, u64, String)],
            merge_distance: u64,
            min_support: usize,
        ) -> Vec<DiscoveredCopy> {
            // Flatten to one site per ALIGNED BLOCK of every placement that lands outside every existing copy.
            // `pid` identifies the contributing placement so the strand vote counts each placement once, not once
            // per exon block. Exclusion is per-PLACEMENT (a placement with any block inside a catalog copy is
            // dropped whole); clustering is per-BLOCK, so a spliced read's intron never chains two real loci into
            // one candidate.
            struct Site {
                pid: usize,
                name: String,
                chrom: String,
                start: u64,
                end: u64,
                reverse: bool,
            }
            let mut sites: Vec<Site> = Vec::new();
            let mut pid = 0usize;
            for (name, placements) in tied_reads {
                for p in placements {
                    if !inside_any_copy(existing_copies, p) {
                        for &(start, end) in &p.blocks {
                            sites.push(Site {
                                pid,
                                name: name.clone(),
                                chrom: p.chrom.clone(),
                                start,
                                end,
                                reverse: p.reverse,
                            });
                        }
                    }
                    pid += 1;
                }
            }
            // Sort by (chrom, start, end) so overlap/proximity clustering is a single linear pass. `sort_by` is
            // stable, so equal keys keep `tied_reads`' own (already deterministic) order.
            sites.sort_by(|a, b| {
                (a.chrom.as_str(), a.start, a.end).cmp(&(b.chrom.as_str(), b.start, b.end))
            });

            struct Cluster {
                chrom: String,
                start: u64,
                end: u64,
                read_names: Vec<String>,
                seen_names: HashSet<String>,
                voted: HashSet<usize>,
                n_fwd: usize,
                n_rev: usize,
            }
            let mut clusters: Vec<Cluster> = Vec::new();
            for s in sites {
                let merge = match clusters.last() {
                    Some(last) => last.chrom == s.chrom && s.start <= last.end + merge_distance,
                    None => false,
                };
                if !merge {
                    clusters.push(Cluster {
                        chrom: s.chrom.clone(),
                        start: s.start,
                        end: s.end,
                        read_names: Vec::new(),
                        seen_names: HashSet::new(),
                        voted: HashSet::new(),
                        n_fwd: 0,
                        n_rev: 0,
                    });
                }
                let c = clusters.last_mut().expect("just pushed or matched");
                c.end = c.end.max(s.end);
                // `seen_names` keeps membership O(1) while `read_names` keeps insertion order for the report.
                if c.seen_names.insert(s.name.clone()) {
                    c.read_names.push(s.name);
                }
                if c.voted.insert(s.pid) {
                    if s.reverse {
                        c.n_rev += 1;
                    } else {
                        c.n_fwd += 1;
                    }
                }
            }

            clusters
                .into_iter()
                .filter(|c| c.read_names.len() >= min_support)
                .map(|c| {
                    let (nearest_copy_tid, nearest_copy_distance) =
                        nearest_copy(existing_copies, &c.chrom, c.start, c.end);
                    DiscoveredCopy {
                        family_id: family_id.to_string(),
                        chrom: c.chrom,
                        start: c.start,
                        end: c.end,
                        // Majority vote; an exact tie defaults to `+` (design doc's Open Question, resolved).
                        strand: if c.n_rev > c.n_fwd { '-' } else { '+' },
                        n_supporting_reads: c.read_names.len(),
                        read_names: c.read_names,
                        nearest_copy_tid,
                        nearest_copy_distance,
                    }
                })
                .collect()
        }

        /// Every AS-tied read in the region, with all of its own max-AS placements as aligned blocks.
        ///
        /// A read is AS-tied here iff at least two of its non-supplementary records share its maximum `AS`.
        pub fn tie_partner_placements(bam_reads: &[BamRead]) -> Vec<(String, Vec<TiePlacement>)> {
            let mut by_name: HashMap<&str, Vec<&BamRead>> = HashMap::new();
            for br in bam_reads.iter().filter(|b| !b.is_supplementary) {
                by_name.entry(br.name.as_str()).or_default().push(br);
            }
            let mut out = Vec::new();
            for (name, placements) in by_name {
                if placements.len() < 2 {
                    continue;
                }
                let max_as = placements.iter().map(|b| b.as_score).max().unwrap();
                let tied: Vec<&BamRead> = placements
                    .iter()
                    .copied()
                    .filter(|b| b.as_score == max_as)
                    .collect();
                if tied.len() >= 2 {
                    out.push((
                        name.to_string(),
                        tied.into_iter()
                            .map(|b| TiePlacement {
                                chrom: b.chrom.clone(),
                                blocks: aligned_blocks(&b.read),
                                reverse: b.reverse,
                            })
                            .collect(),
                    ));
                }
            }
            // Sort by read name for deterministic output order (HashMap iteration is randomized).
            out.sort_by(|a, b| a.0.cmp(&b.0));
            out
        }

        #[cfg(test)]
        mod tests {
            use super::*;

            fn one_copy(
                tid: &str,
                chrom: &str,
                start: u64,
                end: u64,
            ) -> Vec<(String, u64, u64, String)> {
                vec![(chrom.to_string(), start, end, tid.to_string())]
            }

            /// A single-block (unspliced) placement on the forward strand.
            fn pl(chrom: &str, start: u64, end: u64) -> TiePlacement {
                TiePlacement {
                    chrom: chrom.to_string(),
                    blocks: vec![(start, end)],
                    reverse: false,
                }
            }

            fn pl_rev(chrom: &str, start: u64, end: u64) -> TiePlacement {
                TiePlacement {
                    chrom: chrom.to_string(),
                    blocks: vec![(start, end)],
                    reverse: true,
                }
            }

            #[test]
            fn defensively_excludes_positions_inside_a_catalog_copy() {
                let existing = one_copy("c0", "chr1", 1000, 2000);
                // this "tied" position sits INSIDE c0's span -- must never surface as a discovery
                let tied = vec![("read1".to_string(), vec![pl("chr1", 1200, 1300)])];
                let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
                assert!(
                    out.is_empty(),
                    "a position already inside a catalog copy must never be reported"
                );
            }

            #[test]
            fn merges_positions_within_merge_distance_and_respects_min_support() {
                let existing = one_copy("c0", "chr1", 1000, 2000);
                let tied = vec![
                    ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
                    ("read2".to_string(), vec![pl("chr1", 5050, 5150)]), // within 500bp of read1's site
                ];
                let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
                assert_eq!(
                    out.len(),
                    1,
                    "two nearby out-of-catalog positions with 2 supporting reads = 1 cluster"
                );
                assert_eq!(out[0].n_supporting_reads, 2);
                assert_eq!(out[0].nearest_copy_tid, "c0");
                assert_eq!(out[0].nearest_copy_distance, Some(3000)); // 5000 - 2000
            }

            #[test]
            fn keeps_clusters_separate_beyond_merge_distance() {
                let existing = one_copy("c0", "chr1", 1000, 2000);
                let tied = vec![
                    ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
                    ("read2".to_string(), vec![pl("chr1", 5100, 5100)]),
                    ("read3".to_string(), vec![pl("chr1", 9000, 9100)]),
                    ("read4".to_string(), vec![pl("chr1", 9050, 9150)]),
                ];
                let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
                assert_eq!(
                    out.len(),
                    2,
                    "two far-apart pairs must stay two separate clusters"
                );
            }

            #[test]
            fn drops_clusters_below_min_support() {
                let existing = one_copy("c0", "chr1", 1000, 2000);
                let tied = vec![("read1".to_string(), vec![pl("chr1", 5000, 5100)])];
                let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
                assert!(
                    out.is_empty(),
                    "a single supporting read must not clear min_support=2"
                );
            }

            #[test]
            fn a_copy_inside_a_spliced_out_intron_does_not_exclude_the_placement() {
                // Mirrors `block_overlap_ignores_an_intron_that_merely_spans_the_window`
                // (src/bin/copy_assign.rs) at the containment layer. ref_start=1000, "10M5000N10M": aligned blocks
                // are [1000,1010) and [6010,6020); the intron covers [1010,6010) with no aligned base in it. A
                // catalog copy sitting entirely inside that intron must NOT swallow the placement -- with the old
                // bounding-span test (`ref_start..ref_end` = 1000..6020) it did, silently killing the candidate.
                let spliced = TiePlacement {
                    chrom: "chr1".to_string(),
                    blocks: aligned_blocks(&AlignedRead {
                        ref_start: 1000,
                        cigar: vec![('M', 10), ('N', 5000), ('M', 10)],
                        seq: vec![],
                        qual: vec![],
                    }),
                    reverse: false,
                };
                assert_eq!(spliced.blocks, vec![(1000, 1010), (6010, 6020)]);
                let copy_in_the_intron = one_copy("c0", "chr1", 3000, 3100);
                let tied = vec![
                    ("read1".to_string(), vec![spliced.clone()]),
                    ("read2".to_string(), vec![spliced.clone()]),
                ];
                let out = cluster_tie_partners(&tied, "FAM0", &copy_in_the_intron, 500, 2);
                assert_eq!(
                    out.len(),
                    2,
                    "both aligned blocks survive as their own candidates; neither is 'inside' c0"
                );
                assert_eq!((out[0].start, out[0].end), (1000, 1010));
                assert_eq!(
                    (out[1].start, out[1].end),
                    (6010, 6020),
                    "the 5kb intron never chains the two blocks"
                );

                // Contrast: a copy that genuinely overlaps one of the ALIGNED blocks does exclude the placement.
                let copy_on_the_block = one_copy("c0", "chr1", 1005, 1100);
                let out2 = cluster_tie_partners(&tied, "FAM0", &copy_on_the_block, 500, 2);
                assert!(
                    out2.is_empty(),
                    "a real aligned-base overlap must still exclude the whole placement"
                );
            }

            #[test]
            fn strand_is_the_majority_vote_and_a_tie_defaults_to_plus() {
                let existing = one_copy("c0", "chr1", 1000, 2000);
                // 2 forward + 1 reverse -> '+'
                let fwd_majority = vec![
                    ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
                    ("read2".to_string(), vec![pl("chr1", 5050, 5150)]),
                    ("read3".to_string(), vec![pl_rev("chr1", 5060, 5160)]),
                ];
                let out = cluster_tie_partners(&fwd_majority, "FAM0", &existing, 500, 2);
                assert_eq!(out.len(), 1);
                assert_eq!(out[0].strand, '+', "2 forward vs 1 reverse");

                // 2 reverse + 1 forward -> '-'
                let rev_majority = vec![
                    ("read1".to_string(), vec![pl_rev("chr1", 5000, 5100)]),
                    ("read2".to_string(), vec![pl_rev("chr1", 5050, 5150)]),
                    ("read3".to_string(), vec![pl("chr1", 5060, 5160)]),
                ];
                let out = cluster_tie_partners(&rev_majority, "FAM0", &existing, 500, 2);
                assert_eq!(out.len(), 1);
                assert_eq!(out[0].strand, '-', "2 reverse vs 1 forward");

                // exact tie -> '+' by the documented default
                let tie = vec![
                    ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
                    ("read2".to_string(), vec![pl_rev("chr1", 5050, 5150)]),
                ];
                let out = cluster_tie_partners(&tie, "FAM0", &existing, 500, 2);
                assert_eq!(out.len(), 1);
                assert_eq!(out[0].strand, '+', "an exact strand tie defaults to '+'");
            }

            #[test]
            fn a_multi_exon_placement_votes_once_for_strand() {
                // One reverse placement with 3 exon blocks close enough to land in ONE cluster must not outvote
                // two forward single-block reads: the vote is per PLACEMENT, not per block.
                let existing = one_copy("c0", "chr1", 1000, 2000);
                let three_exons = TiePlacement {
                    chrom: "chr1".to_string(),
                    blocks: vec![(5000, 5050), (5100, 5150), (5200, 5250)],
                    reverse: true,
                };
                let tied = vec![
                    ("read1".to_string(), vec![pl("chr1", 5000, 5100)]),
                    ("read2".to_string(), vec![pl("chr1", 5050, 5150)]),
                    ("read3".to_string(), vec![three_exons]),
                ];
                let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
                assert_eq!(out.len(), 1);
                assert_eq!(
                    out[0].n_supporting_reads, 3,
                    "distinct read NAMES, not blocks"
                );
                assert_eq!(
                    out[0].strand, '+',
                    "2 forward placements outvote 1 reverse placement's 3 blocks"
                );
            }

            #[test]
            fn no_copy_on_this_chromosome_yields_na_not_a_sentinel() {
                // A cross-chromosome family: the candidate lands on chr2, every catalog copy is on chr1.
                let existing = one_copy("c0", "chr1", 1000, 2000);
                let tied = vec![
                    ("read1".to_string(), vec![pl("chr2", 5000, 5100)]),
                    ("read2".to_string(), vec![pl("chr2", 5050, 5150)]),
                ];
                let out = cluster_tie_partners(&tied, "FAM0", &existing, 500, 2);
                assert_eq!(out.len(), 1);
                assert_eq!(out[0].nearest_copy_tid, "NA");
                assert_eq!(
                    out[0].nearest_copy_distance, None,
                    "no copy on chr2 -> NA, never u64::MAX"
                );
            }

            #[test]
            fn tie_partner_placements_finds_reads_tied_at_their_own_max_as() {
                use crate::family::copy_split::AlignedRead;
                use crate::family::denovo_assemble::BamRead;
                let mk = |name: &str, chrom: &str, start: u64, as_score: i32| BamRead {
                    chrom: chrom.into(),
                    read: AlignedRead {
                        ref_start: start,
                        cigar: vec![('M', 100)],
                        seq: vec![],
                        qual: vec![],
                    },
                    mapq: 0,
                    name: name.into(),
                    as_score,
                    de: 0.0,
                    is_supplementary: false,
                    is_secondary: as_score != 200,
                    reverse: false,
                    ts: None,
                };
                let reads = vec![
                    mk("tied_read", "chr1", 1000, 200),   // best
                    mk("tied_read", "chr1", 5000, 200),   // tied with the above
                    mk("tied_read", "chr1", 9000, 150),   // worse, not part of the tie
                    mk("unique_read", "chr1", 2000, 300), // only one placement, never tied
                ];
                let out = tie_partner_placements(&reads);
                assert_eq!(out.len(), 1);
                assert_eq!(out[0].0, "tied_read");
                assert_eq!(
                    out[0].1.len(),
                    2,
                    "only the 2 max-scoring placements, not the 150-scoring one"
                );
                assert_eq!(
                    out[0].1[0].blocks,
                    vec![(1000, 1100)],
                    "placements carry aligned blocks, not a span"
                );
            }

            #[test]
            fn tie_partner_placements_carries_the_strand_flag_and_splits_on_introns() {
                use crate::family::copy_split::AlignedRead;
                use crate::family::denovo_assemble::BamRead;
                let mk = |start: u64, reverse: bool| BamRead {
                    chrom: "chr1".into(),
                    read: AlignedRead {
                        ref_start: start,
                        cigar: vec![('M', 10), ('N', 500), ('M', 10)],
                        seq: vec![],
                        qual: vec![],
                    },
                    mapq: 0,
                    name: "r".into(),
                    as_score: 200,
                    de: 0.0,
                    is_supplementary: false,
                    is_secondary: false,
                    reverse,
                    ts: None,
                };
                let out = tie_partner_placements(&[mk(100, false), mk(9000, true)]);
                assert_eq!(out.len(), 1);
                let p = &out[0].1;
                assert_eq!(p.len(), 2);
                assert_eq!(
                    p[0].blocks,
                    vec![(100, 110), (610, 620)],
                    "the intron is not part of any block"
                );
                assert!(!p[0].reverse);
                assert!(p[1].reverse, "FLAG 0x10 is carried through per placement");
            }
        }
    }
}

pub mod collapse_enumerate {
    //! K=0-collapsed family re-admission (behind `--collapse-enumerate`). A near-identical family that
    //! collapses to <2 RNA-distinct loci is re-admitted as copy NUMBER iff it shows a LOCAL collapse:
    //! a `hidden_copy` second-haplotype witness that is BALANCED (co-equal depth) AND projects to >=2 genomic loci.
    //!
    //! **STATUS:** OPT-IN — --collapse-enumerate (src/bin/gw_family_catalog.rs:177-178, default_value_t = false) or env RUSTLE_COLLAPSE_ENUMERATE=1 (denovo_pipeline.rs:180); sibl  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)
    use crate::family::collapse_enumerate::hidden_copy::HiddenCopyEvidence;
    use crate::family::genome_projection::CopyLocus;
    use std::collections::HashMap;

    use crate::family::collapse_enumerate::hidden_copy::{
        detect_hidden_copy, HiddenCopyParams, ReadObs,
    };
    use crate::family::denovo_assemble::{reads_in_region, BamRead};
    use crate::family::genome_projection::{project_families_batch, project_family_copies};
    use crate::genome::GenomeIndex;

    /// Min balanced-alt fraction: a co-equal collapsed 2nd copy (~0.5), not a minor het/edit (~<=0.1).
    pub const MIN_ALT_FRAC: f64 = 0.30;

    /// famCN from a genome projection's loci count. `project_family_copies` is called with the candidate's
    /// own span as the sole `known` entry, so its returned loci are the OTHER genomic copies -- the seed
    /// locus itself is excluded by construction. Total genomic copy number = the seed copy + the projected
    /// others, hence the `+ 1`.
    pub fn famcn_from_projection(n_projection_loci: usize) -> usize {
        n_projection_loci + 1
    }

    /// The three-signal gate. ALL must hold (see spec / Global Constraints).
    pub fn admit_collapse(ev: &HiddenCopyEvidence, n_projection_loci: usize) -> bool {
        ev.flagged && ev.alt_read_fraction >= MIN_ALT_FRAC && n_projection_loci >= 2
    }

    #[derive(Debug, Clone)]
    pub struct CollapsedFamily {
        pub chrom: String,
        pub start: u64,
        pub end: u64,
        pub famcn: usize, // total genomic copy number = seed locus + projected other loci
        pub n_alt_reads: usize, // hidden 2nd-haplotype depth
        pub alt_read_fraction: f64,
        pub projection: Vec<CopyLocus>,
    }

    /// Per-read mismatch positions vs the reference window `[lo,hi)` — the `alts` a `detect_hidden_copy`
    /// column analysis consumes. PRIMARY alignments only (project invariant, `-F 2308`): secondary and
    /// supplementary records are additional placements of the SAME physical molecule, and `detect_hidden_copy`
    /// counts every `ReadObs` as an independent molecule — including them here would double-count exactly the
    /// multimapping reads this feature targets, corrupting `alt_read_fraction`. Only `M/=/X` CIGAR ops advance
    /// both ref and query; `I/S` advance query, `D/N` advance ref. A position counts as an alt iff the read
    /// base differs from the reference base. `refwin` can be shorter than `hi - lo` when `[lo,hi)` runs past
    /// the contig end (subtelomeric candidates); `refwin.get(..)` skips positions past the fetched window
    /// instead of panicking.
    fn read_obs_from_bam_reads(
        reads: &[BamRead],
        chrom: &str,
        lo: u64,
        hi: u64,
        genome: &GenomeIndex,
    ) -> Vec<ReadObs> {
        let refwin = match genome.fetch_sequence(chrom, lo, hi) {
            Some(s) => s,
            None => return Vec::new(),
        };
        let mut out = Vec::with_capacity(reads.len());
        for br in reads {
            if br.is_secondary || br.is_supplementary {
                continue;
            }
            let r = &br.read;
            let mut ref_pos = r.ref_start;
            let mut q = 0usize;
            let mut alts = Vec::new();
            for &(op, len) in &r.cigar {
                match op {
                    'M' | '=' | 'X' => {
                        for k in 0..len {
                            let rp = ref_pos + k;
                            if rp >= lo && rp < hi {
                                if let Some(&rb) = refwin.get((rp - lo) as usize) {
                                    if let Some(&qb) = r.seq.get(q + k as usize) {
                                        if qb.to_ascii_uppercase() != rb.to_ascii_uppercase() {
                                            alts.push(rp);
                                        }
                                    }
                                }
                            }
                        }
                        ref_pos += len;
                        q += len as usize;
                    }
                    'I' | 'S' => {
                        q += len as usize;
                    }
                    'D' | 'N' => {
                        ref_pos += len;
                    }
                    _ => {}
                }
            }
            out.push(ReadObs {
                start: r.ref_start,
                end: ref_pos,
                alts,
            });
        }
        out
    }

    /// Re-admission driver for one dropped collapsed candidate locus. Three-signal gate:
    /// local hidden-copy witness (balanced 2nd haplotype) + >=2 genome-projected loci.
    /// Returns `Some(CollapsedFamily)` on admit, `None` otherwise (including any I/O failure — a
    /// dropped candidate that cannot be evaluated stays dropped, exactly as today).
    pub fn readmit_locus(
        bam_path: &str,
        chrom: &str,
        lo: u64,
        hi: u64,
        consensus: &[u8],
        genome: &GenomeIndex,
        fasta_path: &str,
        minimap2: &str,
        threads: usize,
    ) -> Option<CollapsedFamily> {
        let (_p, bam_reads) = reads_in_region(bam_path, chrom, lo, hi, threads).ok()?;
        let obs = read_obs_from_bam_reads(&bam_reads, chrom, lo, hi, genome);
        let ev = detect_hidden_copy(&obs, &HiddenCopyParams::default());
        if !ev.flagged {
            return None;
        } // short-circuit before the expensive projection
        let known = vec![(chrom.to_string(), lo, hi)];
        let loci =
            project_family_copies(consensus, fasta_path, &known, 0.98, 0.90, minimap2, threads)
                .ok()?;
        if !admit_collapse(&ev, loci.len()) {
            return None;
        }
        Some(CollapsedFamily {
            chrom: chrom.to_string(),
            start: lo,
            end: hi,
            famcn: famcn_from_projection(loci.len()),
            n_alt_reads: ev.n_alt_reads,
            alt_read_fraction: ev.alt_read_fraction,
            projection: loci,
        })
    }

    /// One `<out>.collapsed.tsv` data row for a re-admitted K=0-collapsed family. Columns:
    /// family_id, chrom, start, end, famCN, n_alt_reads, alt_frac(3dp), status, projection_loci
    /// (`chrom:start-end@identity` joined by `;`). Copy-NUMBER only — these families never appear in copies.tsv.
    pub fn format_collapsed_row(family_id: &str, f: &CollapsedFamily) -> String {
        let proj = f
            .projection
            .iter()
            .map(|c| format!("{}:{}-{}@{:.3}", c.chrom, c.start, c.end, c.identity))
            .collect::<Vec<_>>()
            .join(";");
        format!(
            "{family_id}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
            f.chrom,
            f.start,
            f.end,
            f.famcn,
            f.n_alt_reads,
            f.alt_read_fraction,
            "K0_COLLAPSED",
            proj
        )
    }

    /// Minimum PRIMARY reads at a projected locus for it to count as an EXPRESSED copy (the EEF1A1 guard:
    /// silent pseudogenes fall below this).
    pub const MIN_LOCUS_READS: usize = 3;

    /// Expressed-collapsed admission (PSV-free): a dropped candidate is a real multi-copy family iff it
    /// projects to >= 2 genomic loci that are EACH read-supported. No hidden-copy witness required — these
    /// families are exon-identical (0 PSVs) but their copies are all transcribed.
    pub fn admit_expressed_collapse(read_supported_loci: usize) -> bool {
        read_supported_loci >= 2
    }

    #[derive(Debug, Clone)]
    pub struct ExpressedCollapsedFamily {
        pub chrom: String,
        pub start: u64,
        pub end: u64,
        pub famcn: usize, // seed-inclusive: read-supported projected loci + the seed
        pub min_locus_reads: usize, // weakest admitted locus's support (transparency)
        pub projection: Vec<CopyLocus>,
    }

    /// One `<out>.expressed_collapsed.tsv` row: family_id, chrom, start, end, famCN, min_locus_reads,
    /// status (`K0_COLLAPSED_EXPRESSED`), projection_loci (`chrom:start-end@identity` joined by `;`).
    pub fn format_expressed_collapsed_row(family_id: &str, f: &ExpressedCollapsedFamily) -> String {
        let proj = f
            .projection
            .iter()
            .map(|c| format!("{}:{}-{}@{:.3}", c.chrom, c.start, c.end, c.identity))
            .collect::<Vec<_>>()
            .join(";");
        format!(
            "{family_id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            f.chrom, f.start, f.end, f.famcn, f.min_locus_reads, "K0_COLLAPSED_EXPRESSED", proj
        )
    }

    /// One `<out>.dna_family.tsv` row for an RNA-orphan locus recovered by the DNA edge oracle (`--dna-family-
    /// fallback`): same schema as `expressed_collapsed`, status `DNA_FAMILY_RNA_NONHOMOLOGOUS`. The locus's
    /// EXPRESSED transcript is non-homologous to its paralogs (no RNA family forms), yet it projects to >= 2
    /// DIVERGENT genomic copies. Copy NUMBER only; per-read resolution needs DNA parCN.
    pub fn format_dna_family_row(family_id: &str, f: &ExpressedCollapsedFamily) -> String {
        let proj = f
            .projection
            .iter()
            .map(|c| format!("{}:{}-{}@{:.3}", c.chrom, c.start, c.end, c.identity))
            .collect::<Vec<_>>()
            .join(";");
        format!(
            "{family_id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            f.chrom,
            f.start,
            f.end,
            f.famcn,
            f.min_locus_reads,
            "DNA_FAMILY_RNA_NONHOMOLOGOUS",
            proj
        )
    }

    /// Pure: assemble an `ExpressedCollapsedFamily` from already-read-supported projection loci (paired with
    /// their support counts). Returns `None` unless `>= 2` loci are supported. Factored out so the admit/famCN
    /// logic is unit-testable without minimap2/BAM I/O.
    fn build_expressed_family(
        chrom: &str,
        lo: u64,
        hi: u64,
        loci: Vec<CopyLocus>,
        supports: &[usize],
    ) -> Option<ExpressedCollapsedFamily> {
        if !admit_expressed_collapse(loci.len()) {
            return None;
        }
        let min_locus_reads = supports.iter().copied().min().unwrap_or(0);
        Some(ExpressedCollapsedFamily {
            chrom: chrom.to_string(),
            start: lo,
            end: hi,
            famcn: famcn_from_projection(loci.len()),
            min_locus_reads,
            projection: loci,
        })
    }

    /// Batched re-admission for ALL dropped expressed candidates in ONE minimap2 index load (vs re-indexing the
    /// genome per candidate). `candidates` = `(unique_id, chrom, lo, hi, consensus)`. Projects every consensus
    /// at `min_identity`/cov>=0.90 (each candidate's own span excluded via `known`), then per-candidate keeps
    /// projection loci with >= MIN_LOCUS_READS primary reads (the EEF1A1 guard) and admits if >= 2 remain.
    /// `min_identity` = 0.98 for the exon-identical collapse-expressed path; the DNA-family fallback passes a
    /// looser floor (~0.90) to reach DIVERGENT genomic paralogs whose expressed transcript is non-homologous.
    pub fn readmit_expressed_batch(
        candidates: &[(String, String, u64, u64, Vec<u8>)],
        bam_path: &str,
        fasta_path: &str,
        minimap2: &str,
        threads: usize,
        min_identity: f64,
    ) -> Vec<ExpressedCollapsedFamily> {
        if candidates.is_empty() {
            return Vec::new();
        }
        let consensuses: Vec<(String, Vec<u8>)> = candidates
            .iter()
            .map(|(id, _, _, _, seq)| (id.clone(), seq.clone()))
            .collect();
        let known: HashMap<String, Vec<(String, u64, u64)>> = candidates
            .iter()
            .map(|(id, ch, lo, hi, _)| (id.clone(), vec![(ch.clone(), *lo, *hi)]))
            .collect();
        let proj = project_families_batch(
            &consensuses,
            fasta_path,
            &known,
            min_identity,
            0.90,
            minimap2,
            threads,
        )
        .unwrap_or_default();
        let mut out = Vec::new();
        for (id, chrom, lo, hi, _seq) in candidates {
            let loci = match proj.get(id) {
                Some(l) => l.clone(),
                None => continue,
            };
            let mut supported = Vec::new();
            let mut supports = Vec::new();
            for l in loci {
                let n = reads_in_region(bam_path, &l.chrom, l.start, l.end, threads)
                    .map(|(p, _)| p.len())
                    .unwrap_or(0);
                if n >= MIN_LOCUS_READS {
                    supports.push(n);
                    supported.push(l);
                }
            }
            if let Some(f) = build_expressed_family(chrom, *lo, *hi, supported, &supports) {
                out.push(f);
            }
        }
        out
    }

    /// Soft-mask (RepeatMasker lowercase) fraction of a reference sequence: lowercase bases over total. The repeat
    /// signal for the dna-family gate — Alu SINEs and low-complexity, which project to ~99%-identical dispersed
    /// copies genome-wide and masquerade as SD paralogs, are soft-masked in the reference. An empty slice is 0.0.
    pub fn softmask_frac(seq: &[u8]) -> f64 {
        if seq.is_empty() {
            return 0.0;
        }
        let lc = seq.iter().filter(|b| b.is_ascii_lowercase()).count();
        lc as f64 / seq.len() as f64
    }

    /// DNA-family re-admission for RNA-orphan loci (the DNA edge oracle, `--dna-family-fallback`). Projects each
    /// orphan consensus at `min_identity`/cov>=0.90; the seed's own span is excluded via `known`, so the returned
    /// loci are the OTHER genomic paralog copies. Admits as a DNA-family iff >= 1 other genomic copy exists (total
    /// famCN >= 2 = the >= 2-genomic-copy definition). UNLIKE `readmit_expressed_batch`, the projected loci need
    /// NOT be read-supported: a DNA-family's paralogs are frequently SILENT — that is precisely WHY the orphan's
    /// expressed transcript is non-homologous to them (DNA-family != RNA-family). Per-locus read counts are still
    /// computed and reported (min over paralogs; typically 0 for the silent copies) but do NOT gate admission.
    /// Copy-NUMBER only; per-read resolution of the silent copies requires DNA parCN.
    ///
    /// REPEAT GATE (`max_softmask`): a candidate whose re-admitted locus is >= `max_softmask` soft-masked
    /// (repeat/transposon-derived) is REJECTED — an adversarial audit showed the un-gated fallback projects orphan
    /// loci onto Alus/low-complexity (11/12 off-benchmark, both controls) rather than real gene families. This
    /// mirrors Soto's own RepeatMasker exclusion. Soft-mask is read verbatim via `IndexedFasta` (the loaded
    /// `GenomeIndex` upper-cases and cannot feed it).
    pub fn readmit_dna_family_batch(
        candidates: &[(String, String, u64, u64, Vec<u8>)],
        bam_path: &str,
        fasta_path: &str,
        minimap2: &str,
        threads: usize,
        min_identity: f64,
        max_softmask: f64,
    ) -> Vec<ExpressedCollapsedFamily> {
        if candidates.is_empty() {
            return Vec::new();
        }
        let consensuses: Vec<(String, Vec<u8>)> = candidates
            .iter()
            .map(|(id, _, _, _, seq)| (id.clone(), seq.clone()))
            .collect();
        let known: HashMap<String, Vec<(String, u64, u64)>> = candidates
            .iter()
            .map(|(id, ch, lo, hi, _)| (id.clone(), vec![(ch.clone(), *lo, *hi)]))
            .collect();
        let proj = project_families_batch(
            &consensuses,
            fasta_path,
            &known,
            min_identity,
            0.90,
            minimap2,
            threads,
        )
        .unwrap_or_default();
        let fa = crate::genome::IndexedFasta::open(fasta_path).ok();
        let mut out = Vec::new();
        for (id, chrom, lo, hi, _seq) in candidates {
            // repeat gate: reject re-admitted loci that are repeat-dominated (soft-mask >= max_softmask).
            if let Some(fa) = &fa {
                if let Some(bytes) = fa.fetch(chrom, *lo as i64, *hi as i64) {
                    if softmask_frac(&bytes) >= max_softmask {
                        continue;
                    }
                }
            }
            let loci = match proj.get(id) {
                Some(l) => l.clone(),
                None => continue,
            };
            if loci.is_empty() {
                continue;
            } // >= 1 other genomic copy => famCN >= 2 => DNA-family
            let min_reads = loci
                .iter()
                .map(|l| {
                    reads_in_region(bam_path, &l.chrom, l.start, l.end, threads)
                        .map(|(p, _)| p.len())
                        .unwrap_or(0)
                })
                .min()
                .unwrap_or(0);
            out.push(ExpressedCollapsedFamily {
                chrom: chrom.clone(),
                start: *lo,
                end: *hi,
                famcn: famcn_from_projection(loci.len()),
                min_locus_reads: min_reads,
                projection: loci,
            });
        }
        out
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        fn ev(flagged: bool, frac: f64) -> HiddenCopyEvidence {
            HiddenCopyEvidence {
                n_primary_reads: 300,
                n_alt_positions: 40,
                n_alt_reads: (300.0 * frac) as usize,
                alt_read_fraction: frac,
                flagged,
            }
        }
        #[test]
        fn softmask_frac_counts_lowercase_over_total() {
            assert_eq!(softmask_frac(b"ACGT"), 0.0);
            assert_eq!(softmask_frac(b"acgt"), 1.0);
            assert_eq!(softmask_frac(b"ACac"), 0.5);
            assert_eq!(softmask_frac(b""), 0.0);
            assert_eq!(
                softmask_frac(b"ACGTn"),
                0.2,
                "lowercase n (masked) counts as soft-masked"
            );
        }

        #[test]
        fn admits_only_when_all_three_signals_hold() {
            assert!(
                admit_collapse(&ev(true, 0.50), 2),
                "flagged + balanced + >=2 loci -> admit"
            );
            assert!(
                !admit_collapse(&ev(false, 0.50), 2),
                "not flagged -> reject"
            );
            assert!(
                !admit_collapse(&ev(true, 0.10), 2),
                "minor 2nd haplotype (het/edit-like) -> reject"
            );
            assert!(
                !admit_collapse(&ev(true, 0.50), 1),
                "single projection locus -> reject"
            );
        }

        /// The gate is `n_projection_loci >= 2`: exactly `MIN_ALT_FRAC` and exactly 2 projection loci
        /// must ADMIT (not a strict `>` gate on either signal).
        #[test]
        fn admit_collapse_boundary_at_min_alt_frac_and_two_loci() {
            assert!(
                admit_collapse(&ev(true, MIN_ALT_FRAC), 2),
                "alt_read_fraction == MIN_ALT_FRAC exactly -> admit (>=)"
            );
        }

        /// `famcn` is seed-inclusive: `project_family_copies` excludes the seed's own locus (it is the sole
        /// `known` entry), so the projection count is copies OTHER than the seed. Total famCN = seed + others.
        #[test]
        fn famcn_is_seed_inclusive() {
            assert_eq!(
                famcn_from_projection(0),
                1,
                "no other projected loci -> just the seed copy"
            );
            assert_eq!(
                famcn_from_projection(3),
                4,
                "3 other projected loci + the seed copy"
            );
        }

        #[test]
        fn readmit_decision_from_readobs_balanced_vs_het() {
            use crate::family::collapse_enumerate::hidden_copy::{
                detect_hidden_copy, HiddenCopyParams, ReadObs,
            };
            // 20 candidate columns; ~half the reads carry every alt (a co-equal collapsed 2nd copy)
            let cols: Vec<u64> = (0..20).map(|i| 1000 + i * 10).collect();
            let mk = |carry: bool| ReadObs {
                start: 1000,
                end: 1200,
                alts: if carry { cols.clone() } else { vec![] },
            };
            let mut collapse: Vec<ReadObs> = (0..150).map(|_| mk(true)).collect();
            collapse.extend((0..150).map(|_| mk(false))); // 0.50 balanced
            let ev = detect_hidden_copy(&collapse, &HiddenCopyParams::default());
            assert!(ev.flagged && ev.alt_read_fraction >= MIN_ALT_FRAC);
            assert!(admit_collapse(&ev, 2), "balanced collapse + 2 loci admits");
            // minor het: only 8% carry the alts
            let mut het: Vec<ReadObs> = (0..24).map(|_| mk(true)).collect();
            het.extend((0..276).map(|_| mk(false))); // 0.08
            let ev2 = detect_hidden_copy(&het, &HiddenCopyParams::default());
            assert!(!admit_collapse(&ev2, 2), "minor het does not admit");
        }

        use crate::family::copy_split::AlignedRead;

        /// Build a `BamRead` with a single mismatch (alt) baked into an otherwise-reference-matching
        /// `M`-only alignment at `mismatch_offset` (relative to `ref_start`), so `read_obs_from_bam_reads`
        /// always produces exactly one alt -- at `ref_start + mismatch_offset` -- per surviving read. A
        /// distinct offset per read gives each read an identifiable alt position, so tests can confirm
        /// WHICH reads survived a filter, not just how many.
        fn mk_bam_read(
            ref_start: u64,
            len: u64,
            mismatch_offset: u64,
            is_secondary: bool,
            is_supplementary: bool,
            name: &str,
        ) -> BamRead {
            let mut seq = vec![b'A'; len as usize];
            seq[mismatch_offset as usize] = b'C'; // mismatch vs an all-'A' reference at relative offset `mismatch_offset`
            BamRead {
                chrom: "c1".to_string(),
                read: AlignedRead {
                    ref_start,
                    cigar: vec![('M', len)],
                    seq,
                    qual: Vec::new(),
                },
                mapq: 60,
                name: name.to_string(),
                as_score: 0,
                de: 0.0,
                is_supplementary,
                is_secondary,
                reverse: false,
                ts: None,
            }
        }

        #[test]
        fn read_obs_skips_secondary_and_supplementary() {
            let seq = vec![b'A'; 200];
            let genome = GenomeIndex::from_seqs(&[("c1", &seq[..])]);
            // Each read carries a distinct, identifiable alt offset so the surviving ReadObs can be
            // matched back to the exact primary reads that produced them (not merely counted).
            let reads = vec![
                mk_bam_read(10, 50, 5, false, false, "primary1"), // alt at 15
                mk_bam_read(10, 50, 10, true, false, "secondary"), // alt at 20 (must NOT survive)
                mk_bam_read(10, 50, 15, false, true, "supplementary"), // alt at 25 (must NOT survive)
                mk_bam_read(10, 50, 20, false, false, "primary2"),     // alt at 30
            ];
            let obs = read_obs_from_bam_reads(&reads, "c1", 0, 200, &genome);
            assert_eq!(
                obs.len(),
                2,
                "only the two primary (non-secondary, non-supplementary) reads survive"
            );
            let alts: Vec<Vec<u64>> = obs.iter().map(|o| o.alts.clone()).collect();
            assert_eq!(
                alts,
                vec![vec![15], vec![30]],
                "survivors are exactly primary1 (alt@15) and primary2 (alt@30), in read order -- not the \
                 secondary (alt@20) or supplementary (alt@25)"
            );
        }

        #[test]
        fn read_obs_window_past_contig_end_does_not_panic() {
            let seq = vec![b'A'; 100];
            let genome = GenomeIndex::from_seqs(&[("c1", &seq[..])]);
            // Read near the contig end; window requested [0, 100_000) runs far past the 100bp contig,
            // so the fetched reference window is truncated to 100bp while reads may extend past it.
            let reads = vec![mk_bam_read(90, 20, 5, false, false, "primary_near_end")];
            let obs = read_obs_from_bam_reads(&reads, "c1", 0, 100_000, &genome);
            assert_eq!(
                obs.len(),
                1,
                "the read is kept; only its out-of-window positions are dropped"
            );
            assert_eq!(
                obs[0].alts,
                vec![95],
                "the read's in-window alt (ref_start 90 + offset 5 = 95, within the truncated 100bp \
                 reference window) is present -- only positions past the fetched window are dropped"
            );
        }

        #[test]
        fn collapsed_tsv_row_format() {
            use crate::family::genome_projection::CopyLocus;
            let f = CollapsedFamily {
                chrom: "chr2".into(),
                start: 108994973,
                end: 109147842,
                famcn: 2,
                n_alt_reads: 600,
                alt_read_fraction: 0.49,
                projection: vec![
                    CopyLocus {
                        chrom: "chr2".into(),
                        start: 108994973,
                        end: 109147842,
                        identity: 0.99,
                        cov: 0.95,
                    },
                    CopyLocus {
                        chrom: "chr2".into(),
                        start: 110869109,
                        end: 110895544,
                        identity: 0.993,
                        cov: 0.92,
                    },
                ],
            };
            let row = format_collapsed_row("GWFAMc0", &f);
            assert_eq!(row, "GWFAMc0\tchr2\t108994973\t109147842\t2\t600\t0.490\tK0_COLLAPSED\tchr2:108994973-109147842@0.990;chr2:110869109-110895544@0.993");
        }

        #[test]
        fn admit_expressed_needs_two_read_supported_loci() {
            assert!(!admit_expressed_collapse(0));
            assert!(!admit_expressed_collapse(1));
            assert!(admit_expressed_collapse(2));
            assert!(admit_expressed_collapse(5));
        }

        #[test]
        fn expressed_collapsed_row_format() {
            let f = ExpressedCollapsedFamily {
                chrom: "chr2".into(),
                start: 97950885,
                end: 98048181,
                famcn: 3,
                min_locus_reads: 32,
                projection: vec![
                    CopyLocus {
                        chrom: "chr2".into(),
                        start: 97950885,
                        end: 98048181,
                        identity: 0.998,
                        cov: 0.95,
                    },
                    CopyLocus {
                        chrom: "chr2".into(),
                        start: 99100000,
                        end: 99198000,
                        identity: 0.994,
                        cov: 0.93,
                    },
                ],
            };
            assert_eq!(format_expressed_collapsed_row("GWFAMe0", &f),
                "GWFAMe0\tchr2\t97950885\t98048181\t3\t32\tK0_COLLAPSED_EXPRESSED\tchr2:97950885-98048181@0.998;chr2:99100000-99198000@0.994");
        }

        #[test]
        fn expressed_driver_builds_family_from_supported_loci() {
            // exercises the post-projection assembly logic: given read-supported loci, famCN is seed-inclusive and
            // min_locus_reads is the weakest. (The minimap2/BAM I/O path is covered by the Soto/EEF1A1 live runs.)
            let loci = vec![
                CopyLocus {
                    chrom: "c".into(),
                    start: 0,
                    end: 9,
                    identity: 0.99,
                    cov: 0.95,
                },
                CopyLocus {
                    chrom: "c".into(),
                    start: 50,
                    end: 59,
                    identity: 0.995,
                    cov: 0.95,
                },
            ];
            let supports = vec![32usize, 4usize];
            let fam = build_expressed_family("c", 0, 9, loci.clone(), &supports);
            assert!(fam.is_some());
            let fam = fam.unwrap();
            assert_eq!(fam.famcn, 3); // 2 supported loci + seed
            assert_eq!(fam.min_locus_reads, 4); // weakest
                                                // fewer than 2 supported -> None
            assert!(
                build_expressed_family("c", 0, 9, loci[..1].to_vec(), &supports[..1]).is_none()
            );
        }
    }

    // ---- merged 2026-10-05: was `vg_family/hidden_copy.rs`, now the inline module below (one component) ----
    #[allow(clippy::all)]
    pub mod hidden_copy {
        //! Detect gene-family copies PRESENT in the reads but ABSENT from the reference genome
        //! (collapsed segdup / assembly gap / CNV / private duplication).
        //!
        //! A hidden copy's reads have no correct home, so they mismap to the closest sibling reference
        //! copy carrying their PRIVATE SNPs — a COHERENT second haplotype the reference (one copy at
        //! that locus) cannot explain. This detector finds that second haplotype among the locus's
        //! PRIMARY alignments and FLAGS the discrepancy. Per the DAZ3 discipline it DETECTS and reports
        //! evidence ("the reads imply ≥2 copies; the reference models 1") and ABSTAINS from placing or
        //! fabricating the missing copy — it never manufactures a copy.
        //!
        //! Design = synthesis of an independent design panel (statistical / algorithmic / honesty
        //! lenses), which converged on: PRIMARY-alignments-only matrix (the paralog-bleed firewall, since an
        //! in-reference paralog's reads are primary at THEIR locus), candidate columns where the
        //! non-reference allele frequency sits in a balanced band (excludes 0.5% sequencing error and
        //! fixed differences), and a co-segregation/block test (a hidden copy's alt columns co-occur on
        //! ONE read subset — distinguishing it from scattered heterozygous SNPs by requiring many).
        //!
        //! **STATUS:** OPT-IN — `--collapse-enumerate` (src/bin/gw_family_catalog.rs:177-178, `#[arg(long, default_value_t = false)]`)  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

        /// One primary alignment's observation at the locus: its covered span [start, end) and the positions
        /// where it carries a non-reference allele (a mismatch).
        #[derive(Debug, Clone)]
        pub struct ReadObs {
            pub start: u64,
            pub end: u64,
            pub alts: Vec<u64>,
        }

        #[derive(Debug, Clone)]
        pub struct HiddenCopyParams {
            pub balanced_lo: f64, // min alt-allele fraction for a candidate column (≫ error rate)
            pub balanced_hi: f64, // max alt-allele fraction (above = fixed diff / ref error)
            pub min_depth: usize, // min coverage at a candidate column
            pub min_alt_positions: usize, // min candidate columns to call a hidden copy (≫ a few hets)
            pub min_alt_reads: usize, // min reads in the alt haplotype (the hidden copy's depth)
            pub share_hi: f64, // a read joins H if alt at ≥ this fraction of candidate cols it covers
        }

        impl Default for HiddenCopyParams {
            fn default() -> Self {
                HiddenCopyParams {
                    balanced_lo: 0.20,
                    balanced_hi: 0.60,
                    min_depth: 8,
                    min_alt_positions: 12, // het firewall: a diploid het is 1-2 columns; a copy is dozens
                    min_alt_reads: 5,
                    share_hi: 0.60,
                }
            }
        }

        impl HiddenCopyParams {
            pub fn from_env() -> Self {
                let mut p = HiddenCopyParams::default();
                let getf = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<f64>().ok());
                let getu = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<usize>().ok());
                if let Some(v) = getf("RUSTLE_VG_HIDDEN_ALT_LO") {
                    p.balanced_lo = v;
                }
                if let Some(v) = getf("RUSTLE_VG_HIDDEN_ALT_HI") {
                    p.balanced_hi = v;
                }
                if let Some(v) = getu("RUSTLE_VG_HIDDEN_MIN_DEPTH") {
                    p.min_depth = v;
                }
                if let Some(v) = getu("RUSTLE_VG_HIDDEN_MIN_POSITIONS") {
                    p.min_alt_positions = v;
                }
                if let Some(v) = getu("RUSTLE_VG_HIDDEN_MIN_READS") {
                    p.min_alt_reads = v;
                }
                p
            }
        }

        /// Evidence for a copy not in the reference. DETECT + FLAG only — no placement, no sequence.
        #[derive(Debug, Clone, PartialEq)]
        pub struct HiddenCopyEvidence {
            pub n_primary_reads: usize,
            pub n_alt_positions: usize, // coherent second-haplotype columns
            pub n_alt_reads: usize, // reads in the alt haplotype (the hidden copy's apparent depth)
            pub alt_read_fraction: f64, // n_alt_reads / n_primary_reads
            pub flagged: bool,      // evidence of an unmodeled copy at this locus
        }

        /// Pure detector over PRIMARY alignments at one reference-copy locus. The caller MUST pass primary
        /// reads only (the paralog-bleed firewall). Deterministic; no I/O.
        pub fn detect_hidden_copy(reads: &[ReadObs], p: &HiddenCopyParams) -> HiddenCopyEvidence {
            let n = reads.len();
            let none = HiddenCopyEvidence {
                n_primary_reads: n,
                n_alt_positions: 0,
                n_alt_reads: 0,
                alt_read_fraction: 0.0,
                flagged: false,
            };
            if n < p.min_depth {
                return none;
            }

            // Per-position alt-read count (only positions some read calls alt are candidates).
            let mut alt_count: crate::types::DetHashMap<u64, usize> = Default::default();
            for r in reads {
                for &pos in &r.alts {
                    *alt_count.entry(pos).or_insert(0) += 1;
                }
            }

            // Candidate columns: balanced alt fraction over sufficient depth. Error (~0.5%) never
            // reaches balanced_lo; a fixed difference / reference error sits above balanced_hi.
            let mut candidates: Vec<u64> = Vec::new();
            for (&pos, &ac) in &alt_count {
                if ac < 2 {
                    continue; // a singleton is error, not a haplotype
                }
                let cov = reads
                    .iter()
                    .filter(|r| r.start <= pos && pos < r.end)
                    .count();
                if cov < p.min_depth {
                    continue;
                }
                let frac = ac as f64 / cov as f64;
                if frac >= p.balanced_lo && frac <= p.balanced_hi {
                    candidates.push(pos);
                }
            }
            candidates.sort_unstable();
            let n_alt_positions = candidates.len();

            // Co-segregation: a hidden COPY's candidate columns co-occur on ONE read subset (H). Partition
            // reads by their alt-share over the candidate columns they cover; H = the alt haplotype.
            let cand_set: crate::types::DetHashSet<u64> = candidates.iter().copied().collect();
            let n_alt_reads = reads
                .iter()
                .filter(|r| {
                    let covered = candidates
                        .iter()
                        .filter(|&&pos| r.start <= pos && pos < r.end)
                        .count();
                    if covered == 0 {
                        return false;
                    }
                    let alt_at = r.alts.iter().filter(|pos| cand_set.contains(pos)).count();
                    (alt_at as f64 / covered as f64) >= p.share_hi
                })
                .count();

            // Flag only with MANY co-segregating positions (≫ a het) AND a real alt-haplotype read group.
            let flagged = n_alt_positions >= p.min_alt_positions && n_alt_reads >= p.min_alt_reads;

            HiddenCopyEvidence {
                n_primary_reads: n,
                n_alt_positions,
                n_alt_reads,
                alt_read_fraction: if n > 0 {
                    n_alt_reads as f64 / n as f64
                } else {
                    0.0
                },
                flagged,
            }
        }

        #[cfg(test)]
        mod tests {
            use super::*;

            fn p() -> HiddenCopyParams {
                HiddenCopyParams::default()
            }

            // n reads all spanning [0, span); `hap` reads carry alt at every position in `shared`,
            // the rest carry alt at `noise` random-but-distinct positions each (sequencing error).
            fn reads(
                n: usize,
                span: u64,
                hap: usize,
                shared: &[u64],
                noise_per_read: u64,
            ) -> Vec<ReadObs> {
                (0..n)
                    .map(|r| {
                        let mut alts: Vec<u64> = if r < hap { shared.to_vec() } else { Vec::new() };
                        // distinct error positions per read (no cross-read coherence)
                        for k in 0..noise_per_read {
                            alts.push(span - 1 - (r as u64 * 17 + k)); // deterministic, scattered, unique-ish
                        }
                        ReadObs {
                            start: 0,
                            end: span,
                            alts,
                        }
                    })
                    .collect()
            }

            #[test]
            fn hidden_copy_is_flagged() {
                // 60 reads, 30 carry a 20-position shared alt haplotype → a hidden copy.
                let shared: Vec<u64> = (100..120).collect();
                let rs = reads(60, 2000, 30, &shared, 3);
                let e = detect_hidden_copy(&rs, &p());
                assert!(e.flagged);
                assert_eq!(e.n_alt_positions, 20);
                assert_eq!(e.n_alt_reads, 30);
                assert!((e.alt_read_fraction - 0.5).abs() < 1e-9);
            }

            #[test]
            fn sequencing_error_is_not_flagged() {
                // No shared haplotype, just ~9 random errors per read → no candidate columns.
                let rs = reads(60, 2000, 0, &[], 9);
                let e = detect_hidden_copy(&rs, &p());
                assert!(!e.flagged);
                assert_eq!(e.n_alt_positions, 0);
            }

            #[test]
            fn heterozygous_snp_is_not_flagged() {
                // 30 of 60 reads share just 2 alt positions (a diploid het) → below min_alt_positions.
                let rs = reads(60, 2000, 30, &[500, 900], 3);
                let e = detect_hidden_copy(&rs, &p());
                assert!(
                    !e.flagged,
                    "a het (2 positions) must not be called a hidden copy"
                );
                assert!(e.n_alt_positions < p().min_alt_positions);
            }

            #[test]
            fn fixed_difference_above_band_is_not_a_candidate() {
                // A position alt in ~all reads (fixed diff / reference error) is above balanced_hi.
                let rs = reads(40, 2000, 40, &[700], 0); // all 40 alt at 700 → frac 1.0 > 0.60
                let e = detect_hidden_copy(&rs, &p());
                assert_eq!(e.n_alt_positions, 0);
                assert!(!e.flagged);
            }

            #[test]
            fn low_depth_abstains() {
                let shared: Vec<u64> = (100..120).collect();
                let rs = reads(4, 2000, 2, &shared, 0); // below min_depth
                assert!(!detect_hidden_copy(&rs, &p()).flagged);
            }
        }
    }

    // ---- merged 2026-10-05: was `vg_family/collapse_gate.rs`, now the inline module below (one component) ----
    #[allow(clippy::all)]
    pub mod collapse_gate {
        //! Collapse gate: is a single-rep locus actually several copies the aligner could not separate?
        //!
        //! SDA (Vollger et al., Nat Methods 2019) detects a collapse by read-depth excess and only THEN defines PSVs,
        //! "requiring sequence coverages consistent with a single-copy locus in order to distinguish PSVs from allelic
        //! variants". We ran that second stage alone, and a single-copy gene (TSPYL1) reported 12 collapsed copies
        //! against DAZ's 3 — see `bench/COLLAPSED_COPY_GATE.md`. **Collapse first, haplotypes second.**
        //!
        //! In DNA a collapse shows as excess depth. In RNA depth is copy number × expression (Clair3-RNA: "the coverage
        //! is uneven across genomic regions in RNA-seq"), and allele-specific expression destroys the allelic balance
        //! that would otherwise separate a het allele from a PSV ("zygosity flipping can happen"). So we keep SDA's
        //! structure and change its instrument: a collapse shows instead as **reads the aligner cannot place uniquely**.
        //!
        //! Measured on GGO Iso-Seq (primary records only): **0 MAPQ-0 primaries across 9449 reads** at five single-copy
        //! loci (TSPYL1, DERPC, ATXN7L3B, GSPT2, EEF1A1), against 19/20 at DAZ2 and 30/34 at TSPY. The statistic is
        //! expression-invariant — TSPYL1 has 2151 reads and no ambiguity; DAZ2 has 20 reads and 95%.
        //!
        //! ⚠⚠ **DEFAULT OFF. The instrument is not what this module's name claims, and a control proved it.**
        //!
        //! MAPQ 0 means "this read maps equally well somewhere else". It does NOT mean "this locus is collapsed". Run
        //! genome-wide, the gate fires on **EEF1A1** — whose MAPQ-0 reads align to its processed pseudogenes on
        //! NC_073224.2 and NC_073227.2, other chromosomes entirely — and reports `chi(H) = 7` for a locus with one copy.
        //!
        //! The logic actually inverts. If a copy were truly ABSENT from the reference, its reads would pile onto the
        //! present copy at HIGH mapping quality, giving depth excess and *no ambiguity at all*. That is precisely why
        //! SDA detects collapses by read depth and not by mapping quality. Ambiguity detects **unresolvable paralogy**
        //! (which is what `read_conflict`'s E_c oracle already does), not collapse.
        //!
        //! On DAZ the gate emits `chi(H) = 2`, matching the annotation — but DAZ2 is present in the reference 20 kb
        //! away, so even there the signal is paralogy, not collapse. What actually failed at DAZ is assembly: DAZ2 has
        //! 20 primary reads and never becomes a rep.
        //!
        //! Kept, off by default, because the machinery and its tests are sound and the p-value is correct for the
        //! question it truly answers. Do not enable it until `chi(H)` is shown to bound copies rather than haplotypes.
        //! See `bench/COLLAPSE_GATE_VALIDATION.md`.
        //!
        //! **STATUS:** REFUTED  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

        use crate::family::copy_assign::poisson_binomial_upper_tail;
        use crate::family::readonly_copy_number::chi_h;

        /// Ambiguously-placed primary reads (`k`) out of primary reads (`n`) at a locus.
        #[derive(Clone, Copy, Debug, PartialEq, Eq)]
        pub struct Ambiguity {
            pub n: usize,
            pub k: usize,
        }

        /// Background per-read ambiguity rate, under a Jeffreys prior so it is never exactly zero.
        ///
        /// `bg` is pooled over the region's uniquely-mappable reps (those that are not gate candidates). The controls
        /// observe `k = 0`, whose MLE is 0 — and a zero rate would make a single stray MAPQ-0 read infinitely
        /// significant. `None` when there is no background to estimate from: the caller must then ABSTAIN, never fire.
        pub fn estimate_eps_amb(bg: Ambiguity) -> Option<f64> {
            if bg.n == 0 {
                return None;
            }
            Some((bg.k as f64 + 0.5) / (bg.n as f64 + 1.0))
        }

        /// Upper-tail probability of seeing `obs.k` or more ambiguous reads among `obs.n`, if the locus were unique.
        ///
        /// Under the null `k ~ Binomial(n, eps_amb)`. Reuses the shipped Poisson-binomial tail with a constant
        /// probability vector rather than adding a second distribution implementation. That routine is O(n²) in the
        /// number of trials, and here a trial is a READ rather than a distinguishing position — but `k == 0` returns
        /// immediately, and a clean locus is exactly the `k == 0` case, so uniquely-mapping loci cost nothing.
        pub fn collapse_pvalue(obs: Ambiguity, eps_amb: f64) -> f64 {
            if obs.k == 0 || obs.n == 0 {
                return 1.0;
            }
            let probs = vec![eps_amb; obs.n];
            poisson_binomial_upper_tail(obs.k, &probs)
        }

        /// What the gate decided about a locus that has only one assembled rep.
        #[derive(Clone, Debug, PartialEq)]
        pub enum CollapseVerdict {
            /// Collapsed, and its reads resolve into `chi_h >= min_copies` conflicting haplotypes.
            Fire { chi_h: usize, p_value: f64 },
            /// The locus places its reads unambiguously, or its haplotypes do not reach `min_copies`.
            NotCollapsed { p_value: f64 },
            /// Cannot decide — never fire on an unbounded statistic.
            Abstain(&'static str),
        }

        /// Two legs, in SDA's order. Leg 2 is NOT consulted unless leg 1 fires: a single-copy gene reports plenty of
        /// haplotypes (a het allele *is* a haplotype), and gating on them alone made TSPYL1 report 12 copies.
        ///
        /// `eps_amb` is the background per-read ambiguity rate. It must be a **genome-wide** quantity, not a
        /// region-local one: SDA estimates its background from "unique regions" of the genome, and for good reason —
        /// in the DAZ window the only reads that fall outside DAZ1's span are DAZ2's, every one of them ambiguous, so a
        /// region-local background would be ~0.95 and the gate could never fire. Measured on `GGO_mm.bam`:
        /// 5785 MAPQ-0 primaries in 4,404,440, i.e. `eps_amb = 0.0013`, itself conservative because it includes the
        /// genuinely collapsed loci. `None` ⇒ ABSTAIN; never fire on an unbounded statistic.
        ///
        /// `haplotypes` are the `allele_vector`s of the **identifiable** copies at the locus.
        ///
        /// χ(H) is a **lower bound** on copy number, never an estimate: two copies × two alleles also yields four
        /// haplotypes. That is why the caller emits a copy NUMBER with reads certified tied, not an assignment.
        pub fn collapse_verdict(
            obs: Ambiguity,
            eps_amb: Option<f64>,
            haplotypes: &[Vec<Option<u8>>],
            alpha: f64,
            min_copies: usize,
        ) -> CollapseVerdict {
            // leg 1 — is this locus collapsed at all?
            let Some(eps_amb) = eps_amb else {
                return CollapseVerdict::Abstain("no background ambiguity rate (pass --eps-amb)");
            };
            let p_value = collapse_pvalue(obs, eps_amb);
            if p_value >= alpha {
                return CollapseVerdict::NotCollapsed { p_value };
            }
            // leg 2 — and how many copies collapsed?
            let chi = chi_h(haplotypes);
            if chi < min_copies {
                return CollapseVerdict::NotCollapsed { p_value };
            }
            CollapseVerdict::Fire {
                chi_h: chi,
                p_value,
            }
        }

        /// Measured background on `GGO_mm.bam`: 5785 MAPQ-0 primaries out of 4,404,440. Conservative — it includes the
        /// genuinely collapsed loci, so the true unique-region rate is lower and the gate is harder to fire, not easier.
        /// Recompute per sample:
        /// `echo $(( $(samtools view -c -F 2308 b.bam) - $(samtools view -c -F 2308 -q 1 b.bam) ))`
        pub const GENOME_WIDE_EPS_AMB: f64 = 0.001313;

        #[cfg(test)]
        mod tests {
            use super::*;

            /// DAZ1 itself: 22 ambiguous of 200. Against the genome-wide background this is overwhelming.
            #[test]
            fn verdict_fires_on_daz1_against_the_genome_wide_background() {
                let v = collapse_verdict(
                    Ambiguity { n: 200, k: 22 },
                    Some(GENOME_WIDE_EPS_AMB),
                    &haps(&[b"ACGT", b"ACGA"]),
                    1e-3,
                    2,
                );
                assert!(
                    matches!(v, CollapseVerdict::Fire { chi_h: 2, .. }),
                    "DAZ1 must fire, got {v:?}"
                );
            }

            /// Three stray ambiguous reads in a well-covered locus must NOT fire at alpha = 1e-3.
            #[test]
            fn verdict_does_not_fire_on_a_few_stray_ambiguous_reads() {
                let v = collapse_verdict(
                    Ambiguity { n: 500, k: 3 },
                    Some(GENOME_WIDE_EPS_AMB),
                    &haps(&[b"ACGT", b"ACGA"]),
                    1e-3,
                    2,
                );
                assert!(
                    matches!(v, CollapseVerdict::NotCollapsed { .. }),
                    "3 strays must not fire, got {v:?}"
                );
            }

            #[test]
            fn eps_amb_is_never_zero_even_when_no_background_read_is_ambiguous() {
                // The five single-copy controls give 0 ambiguous reads in 9449. The MLE is 0, under which ONE stray
                // MAPQ-0 read would be infinitely significant. Jeffreys keeps it strictly positive.
                let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
                assert!(eps > 0.0, "eps_amb must be strictly positive, got {eps}");
                assert!(
                    (eps - 0.5 / 9450.0).abs() < 1e-12,
                    "Jeffreys: (k + 1/2) / (n + 1)"
                );
            }

            #[test]
            fn eps_amb_abstains_without_background_reads() {
                assert_eq!(
                    estimate_eps_amb(Ambiguity { n: 0, k: 0 }),
                    None,
                    "no background => cannot estimate => abstain"
                );
            }

            #[test]
            fn collapse_pvalue_is_significant_for_daz2_and_not_for_a_clean_locus() {
                let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
                let p_daz2 = collapse_pvalue(Ambiguity { n: 20, k: 19 }, eps); // DAZ2: 19 of 20 ambiguous
                assert!(
                    p_daz2 < 1e-6,
                    "DAZ2 must be overwhelmingly significant, got {p_daz2}"
                );
                let p_clean = collapse_pvalue(Ambiguity { n: 2151, k: 0 }, eps); // TSPYL1
                assert!(
                    (p_clean - 1.0).abs() < 1e-12,
                    "k = 0 => p = 1, got {p_clean}"
                );
            }

            #[test]
            fn collapse_pvalue_of_a_single_stray_read_is_not_significant_at_alpha() {
                let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
                let p = collapse_pvalue(Ambiguity { n: 500, k: 1 }, eps);
                assert!(
                    p > 1e-3,
                    "a single stray MAPQ-0 read must not fire the gate, got {p}"
                );
            }

            /// Allele vectors: haplotypes differing at a shared column conflict, so `chi_h` counts them separately.
            fn haps(rows: &[&[u8]]) -> Vec<Vec<Option<u8>>> {
                rows.iter()
                    .map(|r| {
                        r.iter()
                            .map(|&b| if b == b'.' { None } else { Some(b) })
                            .collect()
                    })
                    .collect()
            }

            /// DAZ: 19/20 ambiguous, background clean, two distinguishable haplotypes.
            #[test]
            fn verdict_fires_on_a_collapsed_locus_with_two_haplotypes() {
                let v = collapse_verdict(
                    Ambiguity { n: 20, k: 19 },
                    Some(GENOME_WIDE_EPS_AMB),
                    &haps(&[b"ACGT", b"ACGA"]),
                    1e-3,
                    2,
                );
                match v {
                    CollapseVerdict::Fire { chi_h, .. } => assert_eq!(chi_h, 2),
                    other => panic!("expected Fire, got {other:?}"),
                }
            }

            /// TSPYL1: a single-copy gene whose reads are ALL uniquely placed. Leg 2 would report many haplotypes
            /// (het alleles, editing, isoform noise) — leg 1 must stop it before leg 2 is ever consulted.
            #[test]
            fn verdict_rejects_a_unique_locus_however_many_haplotypes_it_reports() {
                let twelve: Vec<Vec<Option<u8>>> = (0..12u8)
                    .map(|i| vec![Some(b'A' + i), Some(b'C')])
                    .collect();
                let v = collapse_verdict(
                    Ambiguity { n: 2151, k: 0 },
                    Some(GENOME_WIDE_EPS_AMB),
                    &twelve,
                    1e-3,
                    2,
                );
                assert!(
                    matches!(v, CollapseVerdict::NotCollapsed { .. }),
                    "unique locus must never fire, got {v:?}"
                );
            }

            #[test]
            fn verdict_abstains_without_a_background_estimate() {
                let v = collapse_verdict(
                    Ambiguity { n: 20, k: 19 },
                    None,
                    &haps(&[b"AC", b"AG"]),
                    1e-3,
                    2,
                );
                assert!(
                    matches!(v, CollapseVerdict::Abstain(_)),
                    "no background => abstain, got {v:?}"
                );
            }

            /// `min_copies` applies to χ(H), not to the rep count: one haplotype is not a family.
            #[test]
            fn verdict_rejects_a_collapse_that_resolves_to_one_haplotype() {
                let v = collapse_verdict(
                    Ambiguity { n: 20, k: 19 },
                    Some(GENOME_WIDE_EPS_AMB),
                    &haps(&[b"ACGT"]),
                    1e-3,
                    2,
                );
                assert!(
                    matches!(v, CollapseVerdict::NotCollapsed { .. }),
                    "chi_h = 1 < min_copies => no family, got {v:?}"
                );
            }
        }
    }
}

pub mod catalog_input {
    //! O1 CATALOG → O2 INPUT: parse a `gw_family_catalog` `<out>.copies.tsv` (+ the optional
    //! `<out>.copies.fa`) back into the copy set `copy_assign` assigns reads to.
    //!
    //! # Why this module exists
    //!
    //! O1 (`gw_family_catalog`) and O2 (`copy_assign`) already share ONE node type
    //! (`family_detect::DenovoTranscript` = catalog row = O2 copy), ONE rep-build front end, ONE edge engine
    //! and ONE admission primitive — but they shared all of it BY FUNCTION CALL and NOTHING BY FILE. Each
    //! binary re-derived its own families from the BAM, so at defaults they built DIFFERENT objects (measured:
    //! GSTM catalog 4 copies vs `copy_assign` 0 families on 6031 reads) and their family ids (`GWFAM{i}` vs
    //! `CAFAM{i}`) were assigned independently, leaving NO JOIN KEY between the two tables.
    //!
    //! This module is the file-level contract that closes that gap: `copy_assign --families <copies.tsv>`
    //! CONSUMES the O1 roster instead of re-deriving one, so the two objects agree BY CONSTRUCTION rather
    //! than by coincidence, and every emitted row carries the catalog's own `family_id`/`tid`.
    //!
    //! # Contract (fail loudly, never silently drop)
    //!
    //! Every check below is an ERROR, not a filter. A silently dropped copy is exactly the defect class this
    //! audit keeps finding (the `unwrap_or_default()` famCN silent degradation; the subset-BAM traps), and it
    //! would be undetectable here — a copy quietly missing from O2's roster looks identical to a copy O2
    //! legitimately could not assign reads to.
    //!
    //! * the TSV must carry the `gw_family_catalog` header and every required column (parsed BY NAME, so a
    //!   future appended column cannot shift the meaning of a field);
    //! * a CROSS-CHROMOSOME family (RABL2's 5 contigs) is not truncated to whichever copies happen to fall in
    //!   one swept region: `copy_assign` gathers reads for such a family from every one of its copies' own
    //!   (chromosome, span) windows directly (2026-09-15; see `catalog_input::group_families`'s doc and
    //!   `copy_assign.rs`'s cross-chromosome pass) instead of binding it to a single region;
    //! * the exon blocks must be well formed and reconstruct the copy's own `start`/`end` and `n_exon`;
    //! * with `--copies-fa`, EVERY supplied copy must have a FASTA record whose header coordinates match its
    //!   TSV row (the header is `>{family_id}|{copy_idx}|{chrom}:{start}-{end}|{strand}|nexon={n}`);
    //! * without `--copies-fa`, the spliced sequence is rebuilt from the genome at the catalog's own exon
    //!   coordinates via the SAME `build_spliced_seq` the catalog used, and the strand it derives from the
    //!   junction motifs must agree with the strand the catalog recorded (a disagreement means the FASTA is
    //!   not the assembly the catalog was built against).
    //!
    //! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    use std::collections::BTreeMap;

    use anyhow::{bail, Context, Result};

    use crate::family::denovo_pipeline::ColocatedFamily;
    use crate::family::family_detect::DenovoTranscript;
    use crate::genome::GenomeIndex;

    /// One `<out>.copies.tsv` row: a catalog COPY, with its catalog identity (`family_id`, `copy_idx`, `tid`)
    /// preserved so the O2 output can be joined back to the O1 table it came from.
    #[derive(Clone, Debug, PartialEq, Eq)]
    pub struct CatalogCopy {
        pub family_id: String,
        pub copy_idx: usize,
        pub tid: String,
        pub chrom: String,
        pub start: u64,
        pub end: u64,
        pub strand: char,
        pub n_reads: u32,
        /// Half-open genomic exon blocks, ascending — the `exons` column, verbatim.
        pub exons: Vec<(u64, u64)>,
        /// ⭐ O2-8c′ (§6ep): the copy's SEDEF core hull (`core_hull` column, 0-based half-open `s-e`, or absent /
        /// `NA`). Under `copy_assign --psv-genomic` the PSV alignment span is this hull — the segment homologous
        /// to ≥ half the family by construction — not the copy's extent.
        pub core_hull: Option<(u64, u64)>,
        /// ⭐ L2 (`docs/O1_O2_LOOSE_ENDS.md`): the read-supported locus extent (`locus_start`/`locus_end`, 0-based
        /// half-open) written by `mcl_families`; under the genomic read-star this is the copy's alignment target
        /// (replacing the §6fd padding rule). Absent on older catalogs.
        pub locus: Option<(u64, u64)>,
        /// §6ft: `member_status = partner` — a neighbouring unit of another family, aligned to explain read-through
        /// tails, never a candidate for assignment.
        pub partner: bool,
        /// `member_status = candidate` (spec 2026-10-02 §7): a reference-absent copy proposed by `o3_candidates`, written by
        /// the driver's augmentation on its own contig `cand_<family>_<k>`. A full assignment target (unlike a partner);
        /// like a partner it is exempt from the every-copy-has-a-read contract check ([`CatalogCopy::may_have_no_reads`]).
        pub candidate: bool,
    }

    impl CatalogCopy {
        /// True when the `--families` contract lets this copy have no overlapping read: the catalog marks it unexpressed
        /// (`n_reads 0`, an annotated model kept as the unit), it is a partner (§6ft), or it is an O3 candidate (whose
        /// contig only the patch realignment's reads can reach). Every other copy without a read aborts the run.
        pub fn may_have_no_reads(&self) -> bool {
            self.n_reads == 0 || self.partner || self.candidate
        }
    }

    /// A `copy_assign --only-families` / `--skip-families` list: one family id per line, trimmed, blank lines ignored.
    pub fn parse_family_list(text: &str) -> std::collections::BTreeSet<String> {
        text.lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .map(String::from)
            .collect()
    }

    /// Keep the rows of the families in `only` (when given) that are not in `skip` (when given), in file order — the
    /// driver's two-run O2 split (spec 2026-10-02 §7: the candidate families on the augmented inputs, the rest on the
    /// originals). Applied to the parsed rows BEFORE every contract check, so the checks see exactly the families assigned.
    /// Loud, like the rest of the contract: an id the table does not hold is an error (a typo would otherwise drop or keep a
    /// family silently), and so is a selection that keeps no family.
    pub fn select_families(
        copies: Vec<CatalogCopy>,
        only: Option<&std::collections::BTreeSet<String>>,
        skip: Option<&std::collections::BTreeSet<String>>,
    ) -> Result<Vec<CatalogCopy>> {
        if only.is_none() && skip.is_none() {
            return Ok(copies);
        }
        let have: std::collections::BTreeSet<&str> =
            copies.iter().map(|c| c.family_id.as_str()).collect();
        for (flag, list) in [("--only-families", only), ("--skip-families", skip)] {
            let unknown: Vec<&str> = list
                .into_iter()
                .flatten()
                .map(String::as_str)
                .filter(|id| !have.contains(id))
                .collect();
            if !unknown.is_empty() {
                bail!(
                    "{flag} names {} not in the --families table: {}",
                    if unknown.len() == 1 {
                        "a family"
                    } else {
                        "families"
                    },
                    unknown.join(", ")
                );
            }
        }
        let kept: Vec<CatalogCopy> = copies
            .into_iter()
            .filter(|c| {
                only.map_or(true, |o| o.contains(&c.family_id))
                    && skip.map_or(true, |s| !s.contains(&c.family_id))
            })
            .collect();
        if kept.is_empty() {
            bail!("--only-families / --skip-families leave no family of the --families table to assign");
        }
        Ok(kept)
    }

    /// A catalog FAMILY: its rows grouped by `family_id`, in first-seen (file) order.
    ///
    /// `chrom`/`start`/`end` are meaningful only when [`CatalogFamily::is_cross_chrom`] is false — they are the
    /// single-chromosome span every catalog family had until 2026-09-15. A cross-chromosome family's real
    /// extent is a SET of per-chromosome spans, not one triple; use [`CatalogFamily::chrom_spans`] for that
    /// family instead of reading these three fields.
    #[derive(Clone, Debug)]
    pub struct CatalogFamily {
        pub family_id: String,
        pub chrom: String,
        pub start: u64,
        pub end: u64,
        pub copies: Vec<CatalogCopy>,
    }

    impl CatalogFamily {
        /// True when this family's copies do not all share one chromosome.
        pub fn is_cross_chrom(&self) -> bool {
            self.copies.windows(2).any(|w| w[0].chrom != w[1].chrom)
        }

        /// Per-chromosome `(min start, max end)` span across this family's own copies, one entry per distinct
        /// chromosome. For a same-chromosome family this is a single-entry map equal to `(chrom, (start, end))`.
        pub fn chrom_spans(&self) -> BTreeMap<String, (u64, u64)> {
            let mut m: BTreeMap<String, (u64, u64)> = BTreeMap::new();
            for c in &self.copies {
                let e = m.entry(c.chrom.clone()).or_insert((c.start, c.end));
                e.0 = e.0.min(c.start);
                e.1 = e.1.max(c.end);
            }
            m
        }
    }

    /// A `<out>.copies.fa` record: the sequence plus the header fields it is keyed by, so the header can be
    /// CHECKED against the TSV row rather than trusted.
    #[derive(Clone, Debug)]
    pub struct CatalogSeq {
        pub chrom: String,
        pub start: u64,
        pub end: u64,
        pub strand: char,
        pub n_exon: usize,
        pub seq: Vec<u8>,
    }

    /// `(family_id, copy_idx)` → the catalog's own emitted sequence.
    pub type SeqIndex = BTreeMap<(String, usize), CatalogSeq>;

    /// Format an ascending half-open exon block list as the canonical `start-end,start-end,...` string —
    /// the ONE representation `copies.tsv`, `nodes.tsv` and `parse_exons` all speak. Inverse of `introns_of`.
    ///
    /// Takes the fields rather than a `DenovoTranscript` so the emitters can share it without this module
    /// depending on the detector's types; `gw_family_catalog::exon_blocks` is the thin typed wrapper.
    ///
    /// An unspliced copy (no introns) yields the single block `start-end`.
    ///
    /// Defensive: a malformed chain (donor before the running cursor, or acceptor < donor) is SKIPPED
    /// rather than emitted, because a reversed block would read as a valid interval to every consumer.
    /// The dropped bases are then visible as `exon_bp < exon_sum_len` instead of silently corrupting a span.
    pub fn exon_blocks_str(start: u64, end: u64, introns: &[(u64, u64)]) -> String {
        exon_blocks(start, end, introns)
            .into_iter()
            .map(|(a, b)| format!("{a}-{b}"))
            .collect::<Vec<_>>()
            .join(",")
    }

    /// The exon blocks themselves, as ascending half-open `[start, end)` intervals — the complement of
    /// `introns` within `[start, end)`.
    ///
    /// `exon_blocks_str` is a rendering of exactly this list, so the two can never drift: any consumer that
    /// needs the intervals (e.g. fetching per-exon sequence to build a repeat mask positionally parallel to a
    /// rep's own sequence) gets the same walk that produced the `exons` column of the node dump.
    pub fn exon_blocks(start: u64, end: u64, introns: &[(u64, u64)]) -> Vec<(u64, u64)> {
        let mut out = Vec::new();
        let mut prev = start;
        for &(d, a) in introns {
            if d > prev && a >= d {
                out.push((prev, d));
                prev = a;
            }
        }
        if end > prev {
            out.push((prev, end));
        }
        out
    }

    /// Intron `(donor, acceptor)` chain implied by an ascending half-open exon block list. Inverse of
    /// `exon_blocks_str`.
    pub fn introns_of(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
        exons.windows(2).map(|w| (w[0].1, w[1].0)).collect()
    }

    fn parse_exons(s: &str) -> Result<Vec<(u64, u64)>> {
        let mut out = Vec::new();
        for blk in s.split(',') {
            let blk = blk.trim();
            if blk.is_empty() {
                continue;
            }
            let (a, b) = blk
                .split_once('-')
                .with_context(|| format!("malformed exon block {blk:?}"))?;
            let a: u64 = a
                .parse()
                .with_context(|| format!("malformed exon block {blk:?}"))?;
            let b: u64 = b
                .parse()
                .with_context(|| format!("malformed exon block {blk:?}"))?;
            out.push((a, b));
        }
        Ok(out)
    }

    /// Parse a `gw_family_catalog` `<out>.copies.tsv`. Columns are located BY HEADER NAME.
    pub fn parse_copies_tsv(text: &str) -> Result<Vec<CatalogCopy>> {
        let mut lines = text.lines();
        let header = lines
            .next()
            .context("--families file is empty (expected a gw_family_catalog copies.tsv)")?;
        let cols: Vec<&str> = header.split('\t').collect();
        let idx = |name: &str| -> Result<usize> {
            cols.iter().position(|c| *c == name).with_context(|| {
                format!("--families: copies.tsv has no `{name}` column (header was: {header:?})")
            })
        };
        let (i_fam, i_ci, i_tid, i_chrom, i_start, i_end, i_nexon, i_strand, i_reads, i_exons) = (
            idx("family_id")?,
            idx("copy_idx")?,
            idx("tid")?,
            idx("chrom")?,
            idx("start")?,
            idx("end")?,
            idx("n_exon")?,
            idx("strand")?,
            idx("n_reads")?,
            idx("exons")?,
        );
        let i_hull: Option<usize> = cols.iter().position(|c| *c == "core_hull"); // optional (§6ep)
        let i_status: Option<usize> = cols.iter().position(|c| *c == "member_status"); // optional (§6ft partners)
        let i_locus: Option<(usize, usize)> = cols // optional (L2): both columns or neither
            .iter()
            .position(|c| *c == "locus_start")
            .zip(cols.iter().position(|c| *c == "locus_end"));
        let need = cols.len();
        let mut out = Vec::new();
        for (ln, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < need {
                bail!(
                    "--families: copies.tsv line {} has {} fields, expected {need}",
                    ln + 2,
                    f.len()
                );
            }
            let at = |i: usize| f[i];
            let strand_s = at(i_strand);
            let strand = match strand_s {
                "+" => '+',
                "-" => '-',
                other => bail!(
                    "--families: copies.tsv line {}: strand {other:?} is neither + nor -",
                    ln + 2
                ),
            };
            let exons = parse_exons(at(i_exons))
                .with_context(|| format!("--families: copies.tsv line {}", ln + 2))?;
            let n_exon: usize = at(i_nexon)
                .parse()
                .with_context(|| format!("--families: copies.tsv line {}: bad n_exon", ln + 2))?;
            let c = CatalogCopy {
                family_id: at(i_fam).to_string(),
                copy_idx: at(i_ci)
                    .parse()
                    .with_context(|| format!("line {}: bad copy_idx", ln + 2))?,
                tid: at(i_tid).to_string(),
                chrom: at(i_chrom).to_string(),
                start: at(i_start)
                    .parse()
                    .with_context(|| format!("line {}: bad start", ln + 2))?,
                end: at(i_end)
                    .parse()
                    .with_context(|| format!("line {}: bad end", ln + 2))?,
                strand,
                n_reads: at(i_reads)
                    .parse()
                    .with_context(|| format!("line {}: bad n_reads", ln + 2))?,
                exons,
                core_hull: match i_hull.map(at) {
                    None | Some("NA") | Some("") => None,
                    Some(h) => {
                        let (a, b) = h.split_once('-').with_context(|| {
                            format!(
                                "--families: copies.tsv line {}: bad core_hull {h:?}",
                                ln + 2
                            )
                        })?;
                        Some((
                            a.parse()
                                .with_context(|| format!("line {}: bad core_hull", ln + 2))?,
                            b.parse()
                                .with_context(|| format!("line {}: bad core_hull", ln + 2))?,
                        ))
                    }
                },
                partner: i_status.map(at) == Some("partner"),
                candidate: i_status.map(at) == Some("candidate"),
                locus: match i_locus.map(|(a, b)| (at(a), at(b))) {
                    None | Some(("NA", _)) | Some(("", _)) => None,
                    Some((a, b)) => Some((
                        a.parse().with_context(|| {
                            format!(
                                "--families: copies.tsv line {}: bad locus_start {a:?}",
                                ln + 2
                            )
                        })?,
                        b.parse().with_context(|| {
                            format!(
                                "--families: copies.tsv line {}: bad locus_end {b:?}",
                                ln + 2
                            )
                        })?,
                    )),
                },
            };
            if let Some((a, b)) = c.locus {
                if a > c.start || b < c.end {
                    bail!(
                        "--families: {} copy {} ({}:{}-{}): locus extent {a}-{b} does not contain the copy's span",
                        c.family_id, c.copy_idx, c.chrom, c.start, c.end
                    );
                }
            }
            // Structural checks on the row itself: an exon chain that does not reconstruct the copy's own
            // span/exon count means the row is not the object the catalog wrote.
            if c.exons.is_empty() {
                bail!(
                    "--families: {} copy {} has no exon blocks",
                    c.family_id,
                    c.copy_idx
                );
            }
            if c.exons.len() != n_exon {
                bail!(
                    "--families: {} copy {}: n_exon={n_exon} but the exons column has {} blocks",
                    c.family_id,
                    c.copy_idx,
                    c.exons.len()
                );
            }
            for w in c.exons.windows(2) {
                if w[0].1 > w[1].0 {
                    bail!(
                        "--families: {} copy {}: exon blocks are not ascending/disjoint",
                        c.family_id,
                        c.copy_idx
                    );
                }
            }
            if c.exons.iter().any(|&(s, e)| s >= e) {
                bail!(
                    "--families: {} copy {}: an exon block is empty or reversed",
                    c.family_id,
                    c.copy_idx
                );
            }
            if c.exons[0].0 != c.start || c.exons[c.exons.len() - 1].1 != c.end {
                bail!(
                    "--families: {} copy {}: exon blocks span {}-{} but the row says {}-{}",
                    c.family_id,
                    c.copy_idx,
                    c.exons[0].0,
                    c.exons[c.exons.len() - 1].1,
                    c.start,
                    c.end
                );
            }
            out.push(c);
        }
        if out.is_empty() {
            bail!("--families: copies.tsv has a header but no copy rows");
        }
        Ok(out)
    }

    /// Parse a `gw_family_catalog` `<out>.copies.fa`, keyed by `(family_id, copy_idx)` from the header
    /// `>{family_id}|{copy_idx}|{chrom}:{start}-{end}|{strand}|nexon={n}`.
    pub fn parse_copies_fa(text: &str) -> Result<SeqIndex> {
        let mut out: SeqIndex = BTreeMap::new();
        let mut cur: Option<((String, usize), CatalogSeq)> = None;
        for line in text.lines() {
            if let Some(h) = line.strip_prefix('>') {
                if let Some((k, v)) = cur.take() {
                    out.insert(k, v);
                }
                let parts: Vec<&str> = h.split('|').collect();
                if parts.len() < 5 {
                    bail!("--copies-fa: malformed header {h:?} (expected fam|idx|chrom:start-end|strand|nexon=N)");
                }
                let copy_idx: usize = parts[1]
                    .parse()
                    .with_context(|| format!("--copies-fa: bad copy index in {h:?}"))?;
                let (chrom, span) = parts[2]
                    .rsplit_once(':')
                    .with_context(|| format!("--copies-fa: bad locus field in {h:?}"))?;
                let (s, e) = span
                    .split_once('-')
                    .with_context(|| format!("--copies-fa: bad span in {h:?}"))?;
                let strand = match parts[3] {
                    "+" => '+',
                    "-" => '-',
                    other => bail!("--copies-fa: strand {other:?} in {h:?} is neither + nor -"),
                };
                let n_exon: usize = parts[4]
                    .strip_prefix("nexon=")
                    .with_context(|| format!("--copies-fa: missing nexon= in {h:?}"))?
                    .parse()
                    .with_context(|| format!("--copies-fa: bad nexon in {h:?}"))?;
                cur = Some((
                    (parts[0].to_string(), copy_idx),
                    CatalogSeq {
                        chrom: chrom.to_string(),
                        start: s
                            .parse()
                            .with_context(|| format!("--copies-fa: bad start in {h:?}"))?,
                        end: e
                            .parse()
                            .with_context(|| format!("--copies-fa: bad end in {h:?}"))?,
                        strand,
                        n_exon,
                        seq: Vec::new(),
                    },
                ));
            } else if let Some((_, v)) = cur.as_mut() {
                v.seq
                    .extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
            } else if !line.trim().is_empty() {
                bail!("--copies-fa: sequence line before any `>` header");
            }
        }
        if let Some((k, v)) = cur.take() {
            out.insert(k, v);
        }
        Ok(out)
    }

    /// Group parsed rows into families, in first-seen order.
    ///
    /// A cross-chromosome family (RABL2's 5 contigs) used to be refused outright: `copy_assign`'s region sweep
    /// bound a family to the ONE region containing its whole span, so honouring only the copies that fall in
    /// one region would have silently assigned reads against a TRUNCATED roster. As of 2026-09-15 `copy_assign`
    /// instead gathers such a family's reads directly from every one of its copies' own windows, across however
    /// many chromosomes they sit on, and pools them before assignment — so the roster is never truncated and
    /// this function no longer needs to reject the family. `chrom`/`start`/`end` on the returned
    /// [`CatalogFamily`] are only meaningful for a same-chromosome family; call [`CatalogFamily::is_cross_chrom`]
    /// before trusting them.
    pub fn group_families(copies: Vec<CatalogCopy>) -> Result<Vec<CatalogFamily>> {
        let mut order: Vec<String> = Vec::new();
        let mut by_id: BTreeMap<String, Vec<CatalogCopy>> = BTreeMap::new();
        for c in copies {
            if !by_id.contains_key(&c.family_id) {
                order.push(c.family_id.clone());
            }
            by_id.entry(c.family_id.clone()).or_default().push(c);
        }
        let mut out = Vec::new();
        for fid in order {
            let mut cs = by_id
                .remove(&fid)
                .expect("family id was recorded in `order`");
            // Same ordering guarantee `colocated_families`/`colocated_from_copies` give the assignment step;
            // sorting by `(chrom, start)` rather than bare `start` is identical to the old order for every
            // same-chromosome family (the only case that existed before cross-chrom support) and gives a
            // well-defined, chromosome-grouped order for a cross-chrom one.
            cs.sort_by(|a, b| (a.chrom.as_str(), a.start).cmp(&(b.chrom.as_str(), b.start)));
            let chrom = cs[0].chrom.clone();
            let start = cs.iter().map(|c| c.start).min().unwrap_or(0);
            let end = cs.iter().map(|c| c.end).max().unwrap_or(0);
            out.push(CatalogFamily {
                family_id: fid,
                chrom,
                start,
                end,
                copies: cs,
            });
        }
        Ok(out)
    }

    /// Where a copy's spliced sequence came from — reported so a run is never ambiguous about which substrate
    /// its copy set was materialized from.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum SeqSource {
        /// `--copies-fa`: the catalog's OWN emitted bytes. Agreement with O1 then holds by construction.
        CopiesFa,
        /// Rebuilt from `--fasta` at the catalog's exon coordinates via `build_spliced_seq` (the same function
        /// `assemble_gate` used to write the catalog), with the derived strand checked against the TSV.
        Genome,
    }

    /// Materialize one catalog family as the `ColocatedFamily` the assignment stage consumes, KEEPING the
    /// catalog's `family_id` and per-copy `tid` so every emitted row joins back to `copies.tsv`.
    ///
    /// `seqs` = `Some` ⟹ `--copies-fa` (exact catalog bytes, checked against the TSV row); `None` ⟹ rebuild
    /// from `genome`. Every failure is an error — nothing is dropped.
    pub fn to_colocated(
        fam: &CatalogFamily,
        seqs: Option<&SeqIndex>,
        genome: &GenomeIndex,
    ) -> Result<(ColocatedFamily, SeqSource)> {
        let mut copies: Vec<DenovoTranscript> = Vec::with_capacity(fam.copies.len());
        let source = if seqs.is_some() {
            SeqSource::CopiesFa
        } else {
            SeqSource::Genome
        };
        for c in &fam.copies {
            let introns = introns_of(&c.exons);
            let seq: Vec<u8> = match seqs {
                Some(ix) => {
                    let rec = ix.get(&(c.family_id.clone(), c.copy_idx)).with_context(|| {
                        format!(
                            "--copies-fa has no record for {} copy {} ({}:{}-{}); the FASTA does not match the \
                             --families table",
                            c.family_id, c.copy_idx, c.chrom, c.start, c.end
                        )
                    })?;
                    if rec.chrom != c.chrom
                        || rec.start != c.start
                        || rec.end != c.end
                        || rec.strand != c.strand
                        || rec.n_exon != c.exons.len()
                    {
                        bail!(
                            "--copies-fa record for {} copy {} says {}:{}-{} {} nexon={} but copies.tsv says \
                             {}:{}-{} {} nexon={}",
                            c.family_id, c.copy_idx, rec.chrom, rec.start, rec.end, rec.strand, rec.n_exon,
                            c.chrom, c.start, c.end, c.strand, c.exons.len()
                        );
                    }
                    if rec.seq.is_empty() {
                        bail!(
                            "--copies-fa record for {} copy {} is empty",
                            c.family_id,
                            c.copy_idx
                        );
                    }
                    rec.seq.clone()
                }
                None => {
                    let (seq, derived) =
                        crate::family::denovo_assemble::build_spliced_seq(genome, &c.chrom, c.start, c.end, &introns, Some(c.strand))
                            .with_context(|| {
                                format!(
                                    "--families: could not rebuild {} copy {} ({}:{}-{}) from --fasta — the \
                                     junction motifs are not canonical in this assembly, or the contig/coords \
                                     are absent. Pass --copies-fa to use the catalog's own sequences.",
                                    c.family_id, c.copy_idx, c.chrom, c.start, c.end
                                )
                            })?;
                    if derived != c.strand {
                        bail!(
                            "--families: {} copy {} ({}:{}-{}): --fasta gives strand {derived} but copies.tsv \
                             says {} — the FASTA is not the assembly this catalog was built against",
                            c.family_id, c.copy_idx, c.chrom, c.start, c.end, c.strand
                        );
                    }
                    seq
                }
            };
            if let Some(h) = c.core_hull {
                crate::family::copy_assign::copy_assign_pipeline::register_core_hull(&c.tid, h);
            }
            if let Some(l) = c.locus {
                crate::family::copy_assign::copy_assign_pipeline::register_locus_extent(&c.tid, l);
            }
            if c.partner {
                crate::family::copy_assign::copy_assign_pipeline::register_partner(&c.tid);
            }
            copies.push(DenovoTranscript {
                tid: c.tid.clone(),
                chrom: c.chrom.clone(),
                start: c.start,
                end: c.end,
                n_reads: c.n_reads,
                strand: c.strand,
                introns,
                seq,
                ..Default::default()
            });
        }
        Ok((
            ColocatedFamily {
                family_id: fam.family_id.clone(),
                chrom: fam.chrom.clone(),
                start: fam.start,
                end: fam.end,
                copies,
            },
            source,
        ))
    }

    #[cfg(test)]
    mod tests {
        /// `exon_blocks_str` must be exactly a rendering of `exon_blocks` — they are used as the same walk
        /// (the node dump's `exons` column, and the repeat mask that must be positionally parallel to a rep's
        /// own sequence). If they drift, a mask silently misaligns against the sequence it masks.
        #[test]
        fn exon_blocks_and_its_string_rendering_are_the_same_walk() {
            let cases: &[(u64, u64, &[(u64, u64)])] = &[
                (100, 200, &[]),
                (100, 500, &[(200, 300)]),
                (100, 900, &[(200, 300), (400, 550), (700, 800)]),
                (100, 200, &[(50, 80)]),   // intron entirely before the span
                (100, 200, &[(150, 150)]), // degenerate zero-length intron
            ];
            for &(start, end, introns) in cases {
                let blocks = exon_blocks(start, end, introns);
                let rendered: String = blocks
                    .iter()
                    .map(|(a, b)| format!("{a}-{b}"))
                    .collect::<Vec<_>>()
                    .join(",");
                assert_eq!(
                    rendered,
                    exon_blocks_str(start, end, introns),
                    "walk drift at {start}-{end}"
                );
                // The blocks must be ascending, disjoint and inside the span — the properties the mask relies on.
                let mut prev = start;
                for &(a, b) in &blocks {
                    assert!(
                        a >= prev && b > a && b <= end,
                        "bad block {a}-{b} in {start}-{end}"
                    );
                    prev = b;
                }
            }
        }

        /// The exon blocks must sum to the same length the node dump reports as `exon_bp`, because the repeat
        /// mask is built by concatenating per-block fetches and is then length-checked against the rep sequence.
        #[test]
        fn exon_blocks_sum_to_the_spliced_length() {
            let introns = [(200u64, 300u64), (400, 550)];
            let blocks = exon_blocks(100, 700, &introns);
            let total: u64 = blocks.iter().map(|(a, b)| b - a).sum();
            assert_eq!(total, (700 - 100) - (300 - 200) - (550 - 400));
        }

        use super::*;

        const HDR: &str =
            "family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons";

        fn row(
            fam: &str,
            ci: usize,
            chrom: &str,
            s: u64,
            e: u64,
            exons: &str,
            n_exon: usize,
        ) -> String {
            format!(
                "{fam}\t{ci}\tDN_{chrom}_{s}_{n_exon}\t{chrom}\t{s}\t{e}\t{n_exon}\t+\t7\t{exons}"
            )
        }

        #[test]
        fn parses_a_catalog_row_and_keeps_the_catalog_identity() {
            let t = format!(
                "{HDR}\n{}\n",
                row("GWFAM3", 1, "c1", 100, 400, "100-200,300-400", 2)
            );
            let cs = parse_copies_tsv(&t).unwrap();
            assert_eq!(cs.len(), 1);
            assert_eq!(cs[0].family_id, "GWFAM3");
            assert_eq!(cs[0].copy_idx, 1);
            assert_eq!(cs[0].tid, "DN_c1_100_2");
            assert_eq!(cs[0].exons, vec![(100, 200), (300, 400)]);
            assert_eq!(introns_of(&cs[0].exons), vec![(200, 300)]);
        }

        /// L2: the optional `locus_start`/`locus_end` columns are read by name, absent on older catalogs, and must
        /// contain the copy's own span (the extent is a superset by construction).
        #[test]
        fn locus_extent_columns_are_optional_and_must_contain_the_copy() {
            let plain = format!(
                "{HDR}\n{}\n",
                row("F", 0, "c1", 100, 400, "100-200,300-400", 2)
            );
            assert_eq!(parse_copies_tsv(&plain).unwrap()[0].locus, None);
            let with = format!(
                "{HDR}\tmember_status\tlocus_start\tlocus_end\n{}\tdropped\t50\t900\n",
                row("F", 0, "c1", 100, 400, "100-200,300-400", 2)
            );
            assert_eq!(parse_copies_tsv(&with).unwrap()[0].locus, Some((50, 900)));
            let na = format!(
                "{HDR}\tlocus_start\tlocus_end\n{}\tNA\tNA\n",
                row("F", 0, "c1", 100, 400, "100-200,300-400", 2)
            );
            assert_eq!(parse_copies_tsv(&na).unwrap()[0].locus, None);
            let bad = format!(
                "{HDR}\tlocus_start\tlocus_end\n{}\t150\t900\n",
                row("F", 0, "c1", 100, 400, "100-200,300-400", 2)
            );
            let err = parse_copies_tsv(&bad).unwrap_err().to_string();
            assert!(err.contains("does not contain the copy's span"), "{err}");
        }

        #[test]
        fn columns_are_located_by_name_not_position() {
            // an EXTRA leading column must not shift the meaning of any field
            let hdr = format!("extra\t{HDR}");
            let t = format!("{hdr}\nZ\t{}\n", row("GWFAM0", 0, "c1", 0, 60, "0-60", 1));
            let cs = parse_copies_tsv(&t).unwrap();
            assert_eq!(cs[0].chrom, "c1");
            assert_eq!(cs[0].end, 60);
        }

        #[test]
        fn a_missing_required_column_is_an_error_not_a_default() {
            let t = "family_id\tcopy_idx\ttid\nGWFAM0\t0\tDN_x\n";
            let e = parse_copies_tsv(t).unwrap_err().to_string();
            assert!(e.contains("chrom"), "{e}");
        }

        #[test]
        fn exon_blocks_must_reconstruct_the_rows_own_span() {
            let t = format!(
                "{HDR}\n{}\n",
                row("GWFAM0", 0, "c1", 100, 400, "100-200,300-390", 2)
            );
            let e = parse_copies_tsv(&t).unwrap_err().to_string();
            assert!(e.contains("exon blocks span"), "{e}");
        }

        #[test]
        fn n_exon_must_agree_with_the_exon_column() {
            let t = format!(
                "{HDR}\n{}\n",
                row("GWFAM0", 0, "c1", 100, 400, "100-200,300-400", 3)
            );
            let e = parse_copies_tsv(&t).unwrap_err().to_string();
            assert!(e.contains("n_exon=3"), "{e}");
        }

        #[test]
        fn a_header_only_file_is_an_error_rather_than_an_empty_roster() {
            let e = parse_copies_tsv(&format!("{HDR}\n"))
                .unwrap_err()
                .to_string();
            assert!(e.contains("no copy rows"), "{e}");
        }

        /// `copy_assign --only-families` / `--skip-families` (the O2 split over candidate families, spec 2026-10-02 §7): a
        /// list is one family id per line (blank lines ignored, ids trimmed); the selection keeps the rows' file order; an id
        /// the table does not hold is an error (a typo would otherwise silently drop or keep a family), and so is a selection
        /// that keeps nothing.
        #[test]
        fn family_lists_select_rows_by_family_id_and_fail_loudly() {
            let t = format!(
                "{HDR}\n{}\n{}\n{}\n{}\n",
                row("F1", 0, "c1", 100, 200, "100-200", 1),
                row("F2", 0, "c1", 300, 400, "300-400", 1),
                row("F1", 1, "c2", 100, 200, "100-200", 1),
                row("F3", 0, "c3", 500, 600, "500-600", 1),
            );
            let all = parse_copies_tsv(&t).unwrap();
            let ids = |v: &[CatalogCopy]| {
                v.iter()
                    .map(|c| format!("{}:{}", c.family_id, c.copy_idx))
                    .collect::<Vec<_>>()
            };
            let list = parse_family_list("F1\n\n  F3 \n");
            assert_eq!(
                list.iter().map(String::as_str).collect::<Vec<_>>(),
                vec!["F1", "F3"]
            );
            assert_eq!(
                ids(&select_families(all.clone(), Some(&list), None).unwrap()),
                vec!["F1:0", "F1:1", "F3:0"]
            );
            assert_eq!(
                ids(&select_families(all.clone(), None, Some(&list)).unwrap()),
                vec!["F2:0"]
            );
            assert_eq!(
                select_families(all.clone(), None, None).unwrap(),
                all,
                "no list: the table unchanged"
            );
            // both lists compose: in `only` and not in `skip`
            let f1 = parse_family_list("F1\n");
            assert_eq!(
                ids(&select_families(all.clone(), Some(&list), Some(&f1)).unwrap()),
                vec!["F3:0"]
            );
            // an id the table does not hold, in either list
            let typo = parse_family_list("F1\nF9\n");
            let e = select_families(all.clone(), Some(&typo), None)
                .unwrap_err()
                .to_string();
            assert!(e.contains("--only-families") && e.contains("F9"), "{e}");
            let e = select_families(all.clone(), None, Some(&typo))
                .unwrap_err()
                .to_string();
            assert!(e.contains("--skip-families") && e.contains("F9"), "{e}");
            // a selection that keeps nothing (every family skipped, or an empty only-list)
            let every = parse_family_list("F1\nF2\nF3\n");
            let e = select_families(all.clone(), None, Some(&every))
                .unwrap_err()
                .to_string();
            assert!(e.contains("no family"), "{e}");
            let e = select_families(all.clone(), Some(&parse_family_list("\n")), None)
                .unwrap_err()
                .to_string();
            assert!(e.contains("no family"), "{e}");
        }

        /// `member_status = candidate` (an `o3_candidates` copy added by the driver's augmentation, spec 2026-10-02 §7) is a
        /// full assignment target, NOT a partner (a partner never receives a molecule), but it is exempt from the "every copy
        /// has a read" contract check exactly as a partner or an `n_reads = 0` copy is.
        #[test]
        fn a_candidate_copy_is_a_target_exempt_from_the_read_check_like_a_partner() {
            let t = format!(
                "{HDR}\tmember_status\n{}\tcandidate\n{}\tpartner\n{}\tkept\n",
                row("F", 0, "cand_F_0", 0, 900, "0-900", 1),
                row("F", 1, "c1", 100, 200, "100-200", 1),
                row("F", 2, "c1", 300, 400, "300-400", 1),
            );
            let cs = parse_copies_tsv(&t).unwrap();
            assert!(
                cs[0].candidate && !cs[0].partner,
                "a candidate is not a partner"
            );
            assert!(cs[1].partner && !cs[1].candidate);
            assert!(!cs[2].partner && !cs[2].candidate);
            // the row helper writes n_reads 7: only the status exempts here
            assert!(
                cs[0].may_have_no_reads()
                    && cs[1].may_have_no_reads()
                    && !cs[2].may_have_no_reads()
            );
            let zero = format!(
                "{HDR}\n{}\n",
                row("F", 0, "c1", 100, 200, "100-200", 1).replace("\t7\t", "\t0\t")
            );
            assert!(
                parse_copies_tsv(&zero).unwrap()[0].may_have_no_reads(),
                "n_reads 0 keeps its exemption"
            );
        }

        #[test]
        fn families_group_in_first_seen_order_and_sort_copies_by_start() {
            let t = format!(
                "{HDR}\n{}\n{}\n{}\n",
                row("GWFAM9", 0, "c1", 500, 600, "500-600", 1),
                row("GWFAM9", 1, "c1", 100, 200, "100-200", 1),
                row("GWFAM2", 0, "c1", 900, 950, "900-950", 1),
            );
            let fams = group_families(parse_copies_tsv(&t).unwrap()).unwrap();
            assert_eq!(
                fams.iter()
                    .map(|f| f.family_id.as_str())
                    .collect::<Vec<_>>(),
                vec!["GWFAM9", "GWFAM2"]
            );
            assert_eq!(
                fams[0].copies.iter().map(|c| c.start).collect::<Vec<_>>(),
                vec![100, 500]
            );
            assert_eq!((fams[0].start, fams[0].end), (100, 600));
        }

        /// 2026-09-15: a cross-chromosome family is no longer refused — `copy_assign` gathers its reads
        /// directly from every one of its copies' own chromosomes instead of binding it to one region.
        #[test]
        fn a_cross_chrom_family_is_grouped_not_refused() {
            let t = format!(
                "{HDR}\n{}\n{}\n",
                row("GWFAM0", 0, "c1", 0, 60, "0-60", 1),
                row("GWFAM0", 1, "c2", 0, 60, "0-60", 1),
            );
            let fams = group_families(parse_copies_tsv(&t).unwrap()).unwrap();
            assert_eq!(fams.len(), 1);
            assert!(fams[0].is_cross_chrom());
            assert_eq!(
                fams[0]
                    .copies
                    .iter()
                    .map(|c| c.chrom.as_str())
                    .collect::<Vec<_>>(),
                vec!["c1", "c2"]
            );
        }

        #[test]
        fn chrom_spans_gives_one_min_max_span_per_chromosome() {
            let t = format!(
                "{HDR}\n{}\n{}\n{}\n",
                row("GWFAM0", 0, "c1", 100, 200, "100-200", 1),
                row("GWFAM0", 1, "c1", 300, 400, "300-400", 1),
                row("GWFAM0", 2, "c2", 50, 60, "50-60", 1),
            );
            let fams = group_families(parse_copies_tsv(&t).unwrap()).unwrap();
            assert!(fams[0].is_cross_chrom());
            let spans = fams[0].chrom_spans();
            assert_eq!(spans.get("c1"), Some(&(100, 400)));
            assert_eq!(spans.get("c2"), Some(&(50, 60)));
        }

        #[test]
        fn is_cross_chrom_is_false_for_a_single_chromosome_family() {
            let t = format!(
                "{HDR}\n{}\n{}\n",
                row("GWFAM0", 0, "c1", 0, 60, "0-60", 1),
                row("GWFAM0", 1, "c1", 100, 160, "100-160", 1),
            );
            let fams = group_families(parse_copies_tsv(&t).unwrap()).unwrap();
            assert!(!fams[0].is_cross_chrom());
        }

        #[test]
        fn copies_fa_is_keyed_by_family_and_copy_index_with_checkable_coordinates() {
            let fa =
                ">GWFAM0|1|c1:100-400|-|nexon=2\nacgt\nACGT\n>GWFAM1|0|c2:0-50|+|nexon=1\nTTTT\n";
            let ix = parse_copies_fa(fa).unwrap();
            assert_eq!(ix.len(), 2);
            let r = &ix[&("GWFAM0".to_string(), 1)];
            assert_eq!(
                r.seq,
                b"ACGTACGT".to_vec(),
                "multi-line records concatenate, uppercased"
            );
            assert_eq!(
                (r.chrom.as_str(), r.start, r.end, r.strand, r.n_exon),
                ("c1", 100, 400, '-', 2)
            );
        }

        #[test]
        fn a_malformed_copies_fa_header_is_an_error() {
            let e = parse_copies_fa(">GWFAM0|0|c1:0-10\nAC\n")
                .unwrap_err()
                .to_string();
            assert!(e.contains("malformed header"), "{e}");
        }

        /// The emitted string and the parser are a round trip, and `introns_of` inverts it. If this ever
        /// fails, `nodes.tsv` and `copies.tsv` stop being joinable on exons -- which is the only way to ask
        /// "is this interval inside a rep's EXON or inside its INTRON", a distinction rep spans cannot make
        /// (spans are 90.83% intron by bp).
        #[test]
        fn exon_blocks_str_round_trips_through_parse_and_introns_of() {
            let introns = vec![(1100u64, 1300u64), (1500, 1800)];
            let s = exon_blocks_str(1000, 2000, &introns);
            assert_eq!(s, "1000-1100,1300-1500,1800-2000");
            let parsed = parse_exons(&s).expect("parse");
            assert_eq!(parsed, vec![(1000, 1100), (1300, 1500), (1800, 2000)]);
            assert_eq!(
                introns_of(&parsed),
                introns,
                "introns_of must invert exon_blocks_str exactly"
            );
        }

        /// The sum of the emitted blocks IS the spliced length. `write_er_edge_dump` emits both
        /// (`exon_bp` and `exon_sum_len`) from different code paths precisely so a consumer can check this
        /// invariant per row instead of trusting it.
        #[test]
        fn exon_blocks_str_lengths_sum_to_the_spliced_length() {
            let introns = vec![(1100u64, 1300u64), (1500, 1800)];
            let total: u64 = exon_blocks_str(1000, 2000, &introns)
                .split(',')
                .map(|b| {
                    let (a, z) = b.split_once('-').unwrap();
                    z.parse::<u64>().unwrap() - a.parse::<u64>().unwrap()
                })
                .sum();
            assert_eq!(
                total,
                100 + 200 + 200,
                "exon_bp must equal the exon-sum, not the span"
            );
            assert_eq!(
                exon_blocks_str(500, 900, &[]),
                "500-900",
                "an unspliced rep is one block"
            );
        }

        /// A minus-strand rep must still emit GENOMIC-ASCENDING blocks. The transcription orientation lives
        /// in `strand` and in `seq`; the coordinates never flip. A sibling Python tool (`mkreps.py`) shipped
        /// the opposite convention and scrambled exon order on minus-strand genes, so this is pinned here.
        #[test]
        fn exon_blocks_str_is_genomic_ascending_regardless_of_transcription_strand() {
            let s = exon_blocks_str(1000, 2000, &[(1100, 1300), (1500, 1800)]);
            let starts: Vec<u64> = s
                .split(',')
                .map(|b| b.split_once('-').unwrap().0.parse().unwrap())
                .collect();
            let mut sorted = starts.clone();
            sorted.sort_unstable();
            assert_eq!(
                starts, sorted,
                "blocks must be ascending in GENOMIC coordinates"
            );
        }
    }
}

pub mod copy_split {
    //! Joint read-coherence + PSV "copy split" decomposition.
    //!
    //! Two orthogonal axes group a family's reads into emitted (copy, isoform) units:
    //!   1. STRUCTURAL axis (read-coherence): the exact ordered intron chain a read traverses.
    //!      Different chains are different transcripts regardless of allele content.
    //!   2. COPY axis (PSV haplotype): within one chain, reads that span the family's PSV columns
    //!      carry the allele base distinguishing the paralog copy they came from. When >= 2 distinct
    //!      allele groups differ at >= K columns, the chain SPLITS into multiple identifiable copies.
    //!
    //! The split is deliberately conservative against the indistinguishability wall: a lone sporadic
    //! single-column mismatch never spawns a copy, copies separated by < K columns merge (identifiable
    //! = false), and fully-ambiguous (all-None) reads never force a copy of their own.
    //!
    //! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    use std::collections::BTreeMap;

    /// One read's observation for joint copy+isoform grouping.
    #[derive(Clone, Debug)]
    pub struct ReadObs {
        /// read-coherence STRUCTURAL key: ordered intron chain (donor, acceptor).
        pub intron_chain: Vec<(u64, u64)>,
        /// COPY axis: allele base at each family PSV column (index = column; None if the read
        /// does not span that column). length == n_psv_columns.
        pub psv_alleles: Vec<Option<u8>>,
    }

    /// One emitted (copy, isoform): an intron chain + the PSV haplotype that distinguishes the copy.
    #[derive(Clone, Debug, PartialEq, Eq)]
    pub struct CopyIsoform {
        pub intron_chain: Vec<(u64, u64)>,
        pub allele_vector: Vec<Option<u8>>, // consensus PSV haplotype of this copy (None = unobserved/merged)
        pub read_count: usize,
        pub identifiable: bool, // true if produced by a PSV split; false if a merged (non-identifiable) chain group
    }

    /// A candidate collapsed copy discovered at a locus, ready for the downstream admission gate.
    /// `psv_pos` is parallel to `iso.allele_vector` (same indexing by the discovered PSV columns at
    /// `min_allele_reads=3`).  `n_clusters` is the total number of identifiable copies found at this
    /// locus (including the most-supported host copy that is NOT emitted as a candidate).
    #[derive(Clone, Debug)]
    pub struct CollapsedCandidate {
        pub host_tid: String,
        pub chrom: String,
        pub start: u64,
        pub end: u64,
        pub iso: CopyIsoform,
        pub psv_pos: Vec<u64>, // parallel to iso.allele_vector (genome coords)
        pub n_clusters: usize, // # identifiable copies discovered at this locus
    }

    /// Joint read-coherence + PSV decomposition.
    /// 1. group reads by EXACT intron_chain (read-coherence).
    /// 2. within a chain-group: among reads that span the PSVs, form candidate copies = distinct
    ///    allele vectors with >= min_reads_per_copy supporting reads. If >= 2 candidate copies that
    ///    PAIRWISE differ at >= min_psv_k columns exist, SPLIT: emit one CopyIsoform per such copy
    ///    (identifiable=true); reads not matching a candidate copy (sub-threshold vectors, sporadic
    ///    single-column disagreements, or fully-ambiguous all-None reads) are apportioned to the
    ///    consistent copy if unique, else left shared (counted toward the group but not forcing a copy).
    /// 3. otherwise emit ONE merged CopyIsoform for the chain-group (identifiable=false), read_count=all.
    /// Deterministic: sort output by (intron_chain, allele_vector). Use BTreeMap, no HashMap-into-output.
    pub fn split_readchain_by_psv(
        reads: &[ReadObs],
        n_psv_columns: usize,
        min_psv_k: usize, // K: identifiability threshold (distinct PSV columns)
        min_reads_per_copy: usize,
    ) -> Vec<CopyIsoform> {
        // 1. Group reads by EXACT intron chain (read-coherence / structural axis).
        //    BTreeMap keeps the structural axis deterministically ordered.
        let mut by_chain: BTreeMap<Vec<(u64, u64)>, Vec<&ReadObs>> = BTreeMap::new();
        for r in reads {
            by_chain.entry(r.intron_chain.clone()).or_default().push(r);
        }

        let mut out: Vec<CopyIsoform> = Vec::new();

        for (chain, group) in by_chain {
            // 2. Within a chain-group, tally distinct fully-observed-at-their-columns allele
            //    vectors among reads that span at least one PSV column. A candidate copy is a
            //    distinct allele vector with >= min_reads_per_copy supporting reads. BTreeMap
            //    keeps copy enumeration deterministic.
            let mut vector_counts: BTreeMap<Vec<Option<u8>>, usize> = BTreeMap::new();
            for r in &group {
                // Skip fully-ambiguous reads: they span no column and can never define a copy.
                if r.psv_alleles.iter().all(|a| a.is_none()) {
                    continue;
                }
                *vector_counts.entry(r.psv_alleles.clone()).or_insert(0) += 1;
            }

            // Candidate copies: distinct vectors with enough support.
            let candidates: Vec<Vec<Option<u8>>> = vector_counts
                .iter()
                .filter(|(_, &c)| c >= min_reads_per_copy)
                .map(|(v, _)| v.clone())
                .collect();

            // A valid identifiable split needs the PSV axis to actually carry >= K columns,
            // and >= 2 candidate copies that PAIRWISE differ at >= K columns.
            let split_copies: Vec<Vec<Option<u8>>> = if n_psv_columns >= min_psv_k {
                select_identifiable_copies(&candidates, min_psv_k)
            } else {
                Vec::new()
            };

            if split_copies.len() >= 2 {
                // 2a. SPLIT. Apportion every read in the group: a read joins a copy iff that copy
                //     is the UNIQUE copy consistent with the read's observed (non-None) columns.
                //     Reads consistent with >1 copy (e.g. all-None, or sub-threshold ambiguous
                //     vectors) are left shared and do not force or inflate any copy's count.
                let mut counts = vec![0usize; split_copies.len()];
                for r in &group {
                    let mut consistent: Option<usize> = None;
                    let mut unique = true;
                    for (i, copy) in split_copies.iter().enumerate() {
                        if read_consistent_with(&r.psv_alleles, copy) {
                            if consistent.is_some() {
                                unique = false;
                                break;
                            }
                            consistent = Some(i);
                        }
                    }
                    if unique {
                        if let Some(i) = consistent {
                            counts[i] += 1;
                        }
                    }
                }

                for (copy, count) in split_copies.into_iter().zip(counts) {
                    out.push(CopyIsoform {
                        intron_chain: chain.clone(),
                        allele_vector: copy,
                        read_count: count,
                        identifiable: true,
                    });
                }
            } else {
                // 3. MERGED chain-group: not identifiable. read_count = all reads in the group.
                //    allele_vector = consensus haplotype (per-column majority; None where unobserved
                //    or no majority emerges), purely informational since the copy isn't split out.
                let allele_vector = consensus_haplotype(&group, n_psv_columns);
                out.push(CopyIsoform {
                    intron_chain: chain.clone(),
                    allele_vector,
                    read_count: group.len(),
                    identifiable: false,
                });
            }
        }

        // Deterministic ordering: (intron_chain, allele_vector). by_chain already ordered chains;
        // within a chain copies were enumerated from a BTreeMap (allele-vector ordered). A final
        // stable sort makes the contract explicit regardless of construction order.
        out.sort_by(|a, b| {
            a.intron_chain
                .cmp(&b.intron_chain)
                .then_with(|| a.allele_vector.cmp(&b.allele_vector))
        });
        out
    }

    // ===========================================================================
    // ReadObs BRIDGE: build a ReadObs from a read's spliced alignment + the family's
    // PSV genomic positions. The load-bearing primitive is `allele_at`, which reads a
    // read's base at a REFERENCE position by walking the CIGAR (the error-prone part).
    // ===========================================================================

    /// One read's spliced alignment (minimal model): ref_start (0-based), CIGAR ops as (op,len)
    /// with op in {'M','I','D','N','S'} (match/ins/del/intron/softclip; '=','X' treated as M),
    /// the read sequence (no hard-clipped bases), and the per-base Phred qualities parallel to
    /// `seq`. `qual` is empty when the BAM carried no quality string — callers then fall back to
    /// a flat per-base error; a populated `qual` lets the PSV likelihood weight each base by its
    /// own quality (a HiFi read's distal per-base signal).
    #[derive(Clone, Debug)]
    pub struct AlignedRead {
        pub ref_start: u64,
        pub cigar: Vec<(char, u64)>,
        pub seq: Vec<u8>,
        pub qual: Vec<u8>,
    }

    impl AlignedRead {
        /// The read's aligned blocks as EXONS, 0-based half-open: `M`/`=`/`X`/`D` extend the current block, `N`
        /// closes it (a deletion lies within an exon; only a spliced-out intron separates two). Clips and
        /// insertions consume no reference. Compare `copy_discovery::aligned_blocks`, which also splits on `D`,
        /// and `copy_assign_pipeline::read_ref_end`, which counts `N` toward the span.
        pub fn exon_blocks(&self) -> Vec<(u64, u64)> {
            let mut pos = self.ref_start;
            let mut cur: Option<(u64, u64)> = None;
            let mut out = Vec::new();
            for &(op, n) in &self.cigar {
                match op {
                    'M' | '=' | 'X' | 'D' => {
                        cur = Some((cur.map_or(pos, |c| c.0), pos + n));
                        pos += n;
                    }
                    'N' => {
                        if let Some(c) = cur.take() {
                            out.push(c);
                        }
                        pos += n;
                    }
                    _ => {}
                }
            }
            if let Some(c) = cur {
                out.push(c);
            }
            out
        }
    }

    /// Phred quality `q` -> per-base error probability `10^(-q/10)`, clamped to `[1e-4, 0.25]`
    /// (HiFi QVs run very high; the floor avoids `ln(0)` and the cap avoids over-trusting a
    /// pathologically low QV). A missing/zero QV maps to the cap, i.e. maximally uninformative.
    pub fn phred_err(q: u8) -> f64 {
        if q == 0 {
            return 0.25;
        }
        (10f64.powf(-(q as f64) / 10.0)).clamp(1e-4, 0.25)
    }

    /// Read base aligned to reference position ref_pos (0-based), or None if ref_pos is not a
    /// matched position in this read (inside an intron N, a deletion D, or outside the read span).
    pub fn allele_at(read: &AlignedRead, ref_pos: u64) -> Option<u8> {
        // Walk the CIGAR, tracking the current reference coordinate and read (seq) coordinate.
        // Only M/=/X ops consume BOTH and align a seq base to a ref position; we return that
        // base when the ref coordinate hits ref_pos. N (intron) and D (deletion) advance the
        // ref coordinate but consume no read base -> any ref_pos inside them is unmatched (None).
        // I (insertion) and S (softclip) advance the read coordinate but consume no ref.
        let mut ref_cur = read.ref_start;
        let mut seq_cur: u64 = 0;
        for &(op, len) in &read.cigar {
            match op {
                'M' | '=' | 'X' => {
                    // ref_pos in [ref_cur, ref_cur+len) maps to seq[seq_cur + (ref_pos-ref_cur)].
                    if ref_pos >= ref_cur && ref_pos < ref_cur + len {
                        let off = ref_pos - ref_cur;
                        return read.seq.get((seq_cur + off) as usize).copied();
                    }
                    ref_cur += len;
                    seq_cur += len;
                }
                'N' | 'D' => {
                    // Consumes reference only; no read base aligns here.
                    if ref_pos >= ref_cur && ref_pos < ref_cur + len {
                        return None;
                    }
                    ref_cur += len;
                }
                'I' | 'S' => {
                    // Consumes read only; reference coordinate unchanged.
                    seq_cur += len;
                }
                _ => {
                    // Unknown op (e.g. 'H' hard-clip, 'P' pad): consume nothing here. Hard-clipped
                    // bases are not present in seq, and pads touch neither coordinate.
                }
            }
        }
        None
    }

    /// Intron chain = the (ref_end_of_block, ref_start_of_next_block) gaps produced by N ops.
    pub fn intron_chain_of(read: &AlignedRead) -> Vec<(u64, u64)> {
        let mut out = Vec::new();
        let mut ref_cur = read.ref_start;
        for &(op, len) in &read.cigar {
            match op {
                'M' | '=' | 'X' | 'D' => {
                    // Reference-consuming, alignment-contiguous ops (a D does not break the chain).
                    ref_cur += len;
                }
                'N' => {
                    // Intron: gap from current ref position to current+len.
                    out.push((ref_cur, ref_cur + len));
                    ref_cur += len;
                }
                'I' | 'S' => {
                    // Read-only ops do not move the reference coordinate.
                }
                _ => {}
            }
        }
        out
    }

    /// Build a ReadObs: intron chain from N ops; psv_alleles[i] = allele_at(read, psv_positions[i]).
    pub fn build_read_obs(read: &AlignedRead, psv_positions: &[u64]) -> ReadObs {
        ReadObs {
            intron_chain: intron_chain_of(read),
            psv_alleles: psv_positions.iter().map(|&p| allele_at(read, p)).collect(),
        }
    }

    fn base_index(b: u8) -> Option<usize> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }

    /// Discover within-locus PSV positions from a read PILEUP: genomic positions where `>= 2` distinct ACGT
    /// alleles each have `>= min_allele_reads` support — the signature of `>= 2` collapsed copies at one locus
    /// (also het sites / RNA-editing, which the `>= min_psv_k` requirement in `split_readchain_by_psv` then
    /// filters by demanding a copy span MULTIPLE such columns). Returns the positions sorted ascending.
    pub fn discover_locus_psvs(reads: &[AlignedRead], min_allele_reads: usize) -> Vec<u64> {
        // pileup: genome position -> per-base [A,C,G,T] counts, accumulated in a single walk per read.
        let mut pileup: BTreeMap<u64, [usize; 4]> = BTreeMap::new();
        for read in reads {
            let mut ref_cur = read.ref_start;
            let mut seq_cur = 0u64;
            for &(op, len) in &read.cigar {
                match op {
                    'M' | '=' | 'X' => {
                        for k in 0..len {
                            if let Some(&b) = read.seq.get((seq_cur + k) as usize) {
                                if let Some(bi) = base_index(b) {
                                    pileup.entry(ref_cur + k).or_insert([0; 4])[bi] += 1;
                                }
                            }
                        }
                        ref_cur += len;
                        seq_cur += len;
                    }
                    'N' | 'D' => ref_cur += len,
                    'I' | 'S' => seq_cur += len,
                    _ => {}
                }
            }
        }
        pileup
            .into_iter()
            .filter(|(_, counts)| counts.iter().filter(|&&c| c >= min_allele_reads).count() >= 2)
            .map(|(g, _)| g)
            .collect()
    }

    /// Recover COLLAPSED copies at a single locus: discover within-locus PSVs, then split the reads by PSV
    /// haplotype (`split_readchain_by_psv`). Returns the IDENTIFIABLE (PSV-split) copies — `>= 2` means the
    /// locus is actually multiple collapsed copies the aligner piled onto one place.
    ///
    /// Caveat (measured on GGO): het sites / RNA editing / segdup spillover can mimic a collapsed copy; the
    /// `>= min_psv_k` gate filters single-column noise but not diploid haplotypes, so the real-copy headroom
    /// here was ~0 (the apparent collapses being het / domain-sharer confounds). Genuine wins are in the
    /// collapsed-tandem regime (DAZ/RBMY-like) where multiple real copies share one locus.
    pub fn split_locus_copies(
        reads: &[AlignedRead],
        min_allele_reads: usize,
        min_psv_k: usize,
        min_reads_per_copy: usize,
    ) -> Vec<CopyIsoform> {
        let psv = discover_locus_psvs(reads, min_allele_reads);
        if psv.len() < min_psv_k {
            return Vec::new(); // not enough variant columns to identify a collapsed copy
        }
        let obs: Vec<ReadObs> = reads.iter().map(|r| build_read_obs(r, &psv)).collect();
        split_readchain_by_psv(&obs, psv.len(), min_psv_k, min_reads_per_copy)
            .into_iter()
            .filter(|c| c.identifiable)
            .collect()
    }

    /// True if a read's observed (non-None) PSV columns all agree with `copy`'s alleles.
    /// A None in the read means "did not span" -> imposes no constraint at that column.
    fn read_consistent_with(read: &[Option<u8>], copy: &[Option<u8>]) -> bool {
        read.iter().zip(copy.iter()).all(|(r, c)| match (r, c) {
            (Some(rb), Some(cb)) => rb == cb,
            (Some(_), None) => false, // read observes a base where the copy haplotype is unobserved
            (None, _) => true,        // read did not span this column: no constraint
        })
    }

    /// Greedily select a maximal set of candidate copies that PAIRWISE differ at >= K observed
    /// columns. Candidates are pre-sorted (BTreeMap order); a new candidate is admitted only if it
    /// is >= K-distinct from every already-admitted copy, guaranteeing the emitted copies are
    /// mutually identifiable. Sub-K-distinct near-duplicates fall to the indistinguishability wall.
    fn select_identifiable_copies(
        candidates: &[Vec<Option<u8>>],
        min_psv_k: usize,
    ) -> Vec<Vec<Option<u8>>> {
        let mut chosen: Vec<Vec<Option<u8>>> = Vec::new();
        for cand in candidates {
            if chosen
                .iter()
                .all(|c| distinct_columns(cand, c) >= min_psv_k)
            {
                chosen.push(cand.clone());
            }
        }
        chosen
    }

    /// Count columns where both vectors observe a base and the bases differ.
    fn distinct_columns(a: &[Option<u8>], b: &[Option<u8>]) -> usize {
        a.iter()
            .zip(b.iter())
            .filter(|(x, y)| match (x, y) {
                (Some(xb), Some(yb)) => xb != yb,
                _ => false,
            })
            .count()
    }

    /// Overlay `iso.allele_vector` (parallel to `psv_pos`, genome coords) onto an ALREADY-spliced host sequence
    /// `host` (with its `exon_map`), returning a synthetic transcript carrying the collapsed copy's distinguishing
    /// bases. Substitution-only (v1): a PSV whose genome position is not in the host's exon map (intron/indel) is
    /// SKIPPED and the copy is flagged via `None` return when ANY allele cannot be placed (caller routes to
    /// DNA-needs). Forward-genome coords; the host's own strand/RC is already baked into `host.seq`.
    pub(crate) fn collapsed_copy_to_transcript_from_host_seq(
        iso: &CopyIsoform,
        psv_pos: &[u64],
        host: &crate::family::family_detect::DenovoTranscript,
    ) -> Option<crate::family::family_detect::DenovoTranscript> {
        use crate::family::copy_assign::copy_assign_pipeline::gen2off;
        if iso.allele_vector.len() != psv_pos.len() {
            return None; // parallel-vector invariant violated
        }
        // forward-genome coord -> spliced offset (inverse of exon_map), computed from host.introns so it
        // MUST match the chain that built host.seq (the wrapper passes the copy's chain in both).
        let g2o = gen2off(host);
        let mut seq = host.seq.clone();
        let mut placed = 0usize;
        for (k, &pos) in psv_pos.iter().enumerate() {
            if let Some(base) = iso.allele_vector[k] {
                match g2o.get(&pos) {
                    Some(&off) if off < seq.len() => {
                        seq[off] = base.to_ascii_uppercase();
                        placed += 1;
                    }
                    _ => return None, // PSV not placeable in host exon frame -> DNA-needs (indel/intron)
                }
            }
        }
        if placed == 0 {
            return None; // no distinguishing base placed -> not a usable synthetic copy
        }
        Some(crate::family::family_detect::DenovoTranscript {
            tid: format!("AC_{}_{}", host.chrom, host.start),
            chrom: host.chrom.clone(),
            start: host.start,
            end: host.end,
            n_reads: iso.read_count as u32,
            strand: host.strand,
            introns: iso.intron_chain.clone(),
            seq,
            ..Default::default()
        })
    }

    /// Public wrapper: fetch the host's spliced sequence from the genome (using the COPY's intron chain), then
    /// overlay the discovered alleles. Returns None if the host sequence can't be built or any allele can't be placed.
    pub fn collapsed_copy_to_transcript(
        iso: &CopyIsoform,
        psv_pos: &[u64],
        host: &crate::family::family_detect::DenovoTranscript,
        genome: &crate::genome::GenomeIndex,
    ) -> Option<crate::family::family_detect::DenovoTranscript> {
        use crate::family::denovo_assemble::build_spliced_seq;
        let (seq, strand) = build_spliced_seq(
            genome,
            &host.chrom,
            host.start,
            host.end,
            &iso.intron_chain,
            Some(host.strand),
        )?;
        // `seq` follows the COPY's intron chain, so `host_spliced.introns` MUST be that same chain for
        // `exon_map`/`gen2off` to agree with the seq bytes (the host's own chain may differ — private junction).
        let host_spliced = crate::family::family_detect::DenovoTranscript {
            seq,
            strand,
            introns: iso.intron_chain.clone(),
            ..host.clone()
        };
        collapsed_copy_to_transcript_from_host_seq(iso, psv_pos, &host_spliced)
    }

    /// Per-column majority allele over a merged chain-group; None where no read observes the column
    /// or no strict majority exists. Deterministic (BTreeMap tally, lowest base wins ties).
    fn consensus_haplotype(group: &[&ReadObs], n_psv_columns: usize) -> Vec<Option<u8>> {
        let mut out = Vec::with_capacity(n_psv_columns);
        for col in 0..n_psv_columns {
            let mut tally: BTreeMap<u8, usize> = BTreeMap::new();
            for r in group {
                if let Some(Some(base)) = r.psv_alleles.get(col) {
                    *tally.entry(*base).or_insert(0) += 1;
                }
            }
            // pick the base with the highest count; ties broken by lowest base (BTreeMap order).
            let best = tally
                .iter()
                .max_by(|x, y| x.1.cmp(y.1).then_with(|| y.0.cmp(x.0)))
                .map(|(b, _)| *b);
            out.push(best);
        }
        out
    }

    /// `min_p` identifiability bound (the same construction the assignment gate uses, copy_assign.rs:320):
    /// Π over distinguishing columns of (error_rate/3). `< alpha` ⇒ certifiably distinct.
    pub fn min_p_distinct(
        cand: &[Option<u8>],
        reference: &[Option<u8>],
        error_rate: f64,
        alpha: f64,
    ) -> bool {
        let eps = (error_rate / 3.0).clamp(0.0, 1.0);
        let mut prod = 1.0f64;
        let mut any = false;
        for (a, b) in cand.iter().zip(reference.iter()) {
            if let (Some(x), Some(y)) = (a, b) {
                if x != y {
                    prod *= eps;
                    any = true;
                }
            }
        }
        any && prod < alpha
    }

    /// Reject candidates whose differing columns are ALL A->G (plus-strand) — the RNA-editing signature
    /// (Clair3-RNA). PSV alleles are in transcription orientation, so editing shows as A->G uniformly.
    /// Returns true (=keep) iff at least one differing column is NOT A->G.
    pub fn strand_symmetric_spectrum(
        host_alleles: &[Option<u8>],
        cand_alleles: &[Option<u8>],
    ) -> bool {
        let mut diffs = 0usize;
        let mut non_ag = 0usize;
        for (h, c) in host_alleles.iter().zip(cand_alleles.iter()) {
            if let (Some(hb), Some(cb)) = (h, c) {
                if hb != cb {
                    diffs += 1;
                    if !(*hb == b'A' && *cb == b'G') {
                        non_ag += 1;
                    }
                }
            }
        }
        diffs > 0 && non_ag > 0
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::family::copy_assign::copy_assign_pipeline::exon_map;
        use crate::family::family_detect::DenovoTranscript;

        #[test]
        fn collapsed_copy_to_transcript_overlays_alleles_at_psv_positions() {
            // host: single-exon transcript on '+', chrom "c1", spliced seq fetched from a tiny genome.
            // Build a GenomeIndex stub via the existing test helper (see how other tests build it);
            // here we exercise the OVERLAY arithmetic with a host whose exon_map is identity.
            let host = DenovoTranscript {
                tid: "H".into(),
                chrom: "c1".into(),
                start: 100,
                end: 110,
                n_reads: 9,
                strand: '+',
                introns: vec![],
                seq: b"AAAAAAAAAA".to_vec(),
                ..Default::default()
            };
            // PSV at genome positions 102 and 107 → spliced offsets 2 and 7 (identity exon_map for a single exon).
            let psv_pos = vec![102u64, 107u64];
            let iso = CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'G')],
                read_count: 5,
                identifiable: true,
            };
            // Overlay directly against the host seq (no genome fetch needed for a single-exon identity map):
            let t = collapsed_copy_to_transcript_from_host_seq(&iso, &psv_pos, &host)
                .expect("transcript built");
            assert_eq!(
                t.seq,
                b"AACAAAAGAA".to_vec(),
                "C at offset 2, G at offset 7, rest = host"
            );
            assert_eq!(
                t.seq.len(),
                exon_map(&t).len(),
                "seq/exon_map length invariant"
            );
            assert_eq!(t.chrom, "c1");
            assert_eq!(t.strand, '+');
        }

        #[test]
        fn collapsed_copy_to_transcript_none_when_allele_unplaceable() {
            let host = DenovoTranscript {
                tid: "H".into(),
                chrom: "c1".into(),
                start: 100,
                end: 105,
                n_reads: 9,
                strand: '+',
                introns: vec![],
                seq: b"AAAAA".to_vec(),
                ..Default::default()
            };
            let psv_pos = vec![999u64]; // not in host exon frame
            let iso = CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C')],
                read_count: 5,
                identifiable: true,
            };
            assert!(collapsed_copy_to_transcript_from_host_seq(&iso, &psv_pos, &host).is_none());
        }

        /// REGRESSION: the public wrapper must compute `exon_map` from the COPY's intron chain (the chain that
        /// built the spliced seq), NOT the host's. When the collapsed copy carries a PRIVATE junction differing
        /// from the host's, a PSV exonic under the copy chain falls inside the host's intron — placing it via the
        /// host chain would wrongly return None (or misplace the base). Host intron [105,115); copy intron
        /// [108,112); PSV at genome 106 is exonic (offset 6) under the copy chain but inside the host intron.
        #[test]
        fn collapsed_copy_to_transcript_uses_copy_intron_chain_for_exon_map() {
            use crate::genome::GenomeIndex;
            // Genome c1: all 'A' except the COPY intron [108,112) carrying a canonical GT..AG '+'-junction.
            let mut s = vec![b'A'; 130];
            s[108] = b'G';
            s[109] = b'T'; // donor GT
            s[110] = b'A';
            s[111] = b'G'; // acceptor AG
            let genome = GenomeIndex::from_seqs(&[("c1", &s)]);
            // Host's own chain differs from the collapsed copy's (private-junction case).
            let host = DenovoTranscript {
                tid: "H".into(),
                chrom: "c1".into(),
                start: 100,
                end: 120,
                n_reads: 9,
                strand: '+',
                introns: vec![(105, 115)],
                seq: vec![],
                ..Default::default()
            };
            // Copy chain: exon1 [100,108), exon2 [112,120). PSV genome 106 -> copy offset 6.
            let iso = CopyIsoform {
                intron_chain: vec![(108, 112)],
                allele_vector: vec![Some(b'C')],
                read_count: 5,
                identifiable: true,
            };
            let psv_pos = vec![106u64];
            let t = collapsed_copy_to_transcript(&iso, &psv_pos, &host, &genome)
                .expect("copy-chain exon_map places the exonic PSV");
            assert_eq!(
                t.seq.len(),
                exon_map(&t).len(),
                "seq/exon_map length invariant"
            );
            let mut expect = vec![b'A'; 16]; // 8 + 8 exon bases
            expect[6] = b'C'; // overlaid allele at the COPY-chain offset
            assert_eq!(
                t.seq, expect,
                "C placed at the copy-chain offset 6, rest = reference"
            );
            assert_eq!(
                t.introns,
                vec![(108, 112)],
                "emitted transcript carries the copy intron chain"
            );
            assert_eq!(t.strand, '+');
        }

        /// The length-mismatch invariant branch: `allele_vector` and `psv_pos` of unequal length -> None.
        #[test]
        fn collapsed_copy_to_transcript_none_on_length_mismatch() {
            let host = DenovoTranscript {
                tid: "H".into(),
                chrom: "c1".into(),
                start: 100,
                end: 105,
                n_reads: 9,
                strand: '+',
                introns: vec![],
                seq: b"AAAAA".to_vec(),
                ..Default::default()
            };
            let iso = CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![Some(b'C'), Some(b'G')],
                read_count: 5,
                identifiable: true,
            };
            let psv_pos = vec![102u64]; // 1 position vs 2 alleles
            assert!(collapsed_copy_to_transcript_from_host_seq(&iso, &psv_pos, &host).is_none());
        }

        /// The all-None branch: every allele is unobserved -> nothing placed -> not a usable copy -> None.
        #[test]
        fn collapsed_copy_to_transcript_none_when_all_alleles_none() {
            let host = DenovoTranscript {
                tid: "H".into(),
                chrom: "c1".into(),
                start: 100,
                end: 105,
                n_reads: 9,
                strand: '+',
                introns: vec![],
                seq: b"AAAAA".to_vec(),
                ..Default::default()
            };
            let psv_pos = vec![102u64, 103u64];
            let iso = CopyIsoform {
                intron_chain: vec![],
                allele_vector: vec![None, None],
                read_count: 5,
                identifiable: true,
            };
            assert!(collapsed_copy_to_transcript_from_host_seq(&iso, &psv_pos, &host).is_none());
        }

        fn obs(chain: &[(u64, u64)], alleles: &[Option<u8>]) -> ReadObs {
            ReadObs {
                intron_chain: chain.to_vec(),
                psv_alleles: alleles.to_vec(),
            }
        }

        // Allele shorthands.
        const A: Option<u8> = Some(b'A');
        const C: Option<u8> = Some(b'C');
        const G: Option<u8> = Some(b'G');
        const T: Option<u8> = Some(b'T');
        const N: Option<u8> = None;

        /// Two reference intron chains to play the structural axis.
        fn chain_x() -> Vec<(u64, u64)> {
            vec![(100, 200), (300, 400)]
        }
        fn chain_y() -> Vec<(u64, u64)> {
            vec![(100, 200), (300, 450)]
        }

        // ---- discover_locus_psvs + split_locus_copies (collapsed-copy recovery from raw reads) ----

        /// A single-exon read over [start, start+len) carrying `seq`.
        fn pile_read(start: u64, seq: &[u8]) -> AlignedRead {
            aligned(start, &[('M', seq.len() as u64)], seq)
        }

        #[test]
        fn discover_locus_psvs_finds_split_positions() {
            // background all-G; half the reads carry A and half C at genome positions 130 and 160.
            let mut reads = Vec::new();
            for base in [b'A', b'C'] {
                for _ in 0..5 {
                    let mut s = vec![b'G'; 100];
                    s[30] = base; // genome 130
                    s[60] = base; // genome 160
                    reads.push(pile_read(100, &s));
                }
            }
            assert_eq!(discover_locus_psvs(&reads, 3), vec![130, 160]);
            // raising the support bar above the per-allele depth (5) drops them
            assert!(discover_locus_psvs(&reads, 6).is_empty());
        }

        #[test]
        fn split_locus_recovers_two_collapsed_copies() {
            // reads pile on ONE locus (same single-exon chain) but split into two PSV haplotypes -> 2 copies.
            let mut reads = Vec::new();
            for base in [b'A', b'C'] {
                for _ in 0..5 {
                    let mut s = vec![b'G'; 100];
                    s[30] = base;
                    s[60] = base;
                    reads.push(pile_read(100, &s));
                }
            }
            let copies = split_locus_copies(&reads, 3, 2, 3);
            assert_eq!(copies.len(), 2, "two collapsed copies recovered");
            assert!(copies.iter().all(|c| c.identifiable));
            assert_eq!(
                copies.iter().map(|c| c.read_count).collect::<Vec<_>>(),
                vec![5, 5]
            );
        }

        #[test]
        fn split_locus_no_variant_no_split() {
            let reads: Vec<AlignedRead> =
                (0..10).map(|_| pile_read(100, &vec![b'G'; 100])).collect();
            assert!(
                split_locus_copies(&reads, 3, 2, 3).is_empty(),
                "identical reads -> no collapsed copy"
            );
        }

        #[test]
        fn split_locus_single_variant_is_not_a_copy() {
            // reads differ at ONE position only (a het site / single PSV) -> < min_psv_k=2 -> no split.
            let mut reads = Vec::new();
            for base in [b'A', b'C'] {
                for _ in 0..5 {
                    let mut s = vec![b'G'; 100];
                    s[30] = base;
                    reads.push(pile_read(100, &s));
                }
            }
            assert!(
                split_locus_copies(&reads, 3, 2, 3).is_empty(),
                "a single variant is not a collapsed copy"
            );
        }

        #[test]
        fn split_two_copies_same_chain() {
            let cx = chain_x();
            let mut reads = Vec::new();
            for _ in 0..3 {
                reads.push(obs(&cx, &[A, C, T]));
            }
            for _ in 0..3 {
                reads.push(obs(&cx, &[G, T, A]));
            }
            let out = split_readchain_by_psv(&reads, 3, 2, 2);
            assert_eq!(out.len(), 2, "expected exactly 2 copies");
            assert!(
                out.iter().all(|c| c.identifiable),
                "both must be identifiable"
            );
            let counts: Vec<usize> = out.iter().map(|c| c.read_count).collect();
            assert_eq!(counts, vec![3, 3], "read counts 3 and 3");
            assert_ne!(
                out[0].allele_vector, out[1].allele_vector,
                "distinct allele vectors"
            );
            assert!(out.iter().all(|c| c.intron_chain == cx));
        }

        #[test]
        fn no_oversplit_identical_alleles() {
            let cx = chain_x();
            let reads: Vec<ReadObs> = (0..6).map(|_| obs(&cx, &[A, C, T])).collect();
            let out = split_readchain_by_psv(&reads, 3, 2, 2);
            assert_eq!(out.len(), 1, "single copy, no over-split");
            assert_eq!(out[0].read_count, 6);
            assert!(!out[0].identifiable, "single copy is not identifiable");
        }

        #[test]
        fn different_chains_are_distinct() {
            let cx = chain_x();
            let cy = chain_y();
            let mut reads = Vec::new();
            for _ in 0..3 {
                reads.push(obs(&cx, &[A, C, T]));
            }
            for _ in 0..3 {
                reads.push(obs(&cy, &[A, C, T]));
            }
            let out = split_readchain_by_psv(&reads, 3, 2, 2);
            assert_eq!(out.len(), 2, "two distinct transcripts keyed by chain");
            let chains: Vec<&Vec<(u64, u64)>> = out.iter().map(|c| &c.intron_chain).collect();
            assert!(chains.contains(&&cx));
            assert!(chains.contains(&&cy));
            // identical alleles must NOT merge across chains.
            assert_ne!(out[0].intron_chain, out[1].intron_chain);
        }

        #[test]
        fn sporadic_error_does_not_split() {
            let cx = chain_x();
            let mut reads: Vec<ReadObs> = (0..5).map(|_| obs(&cx, &[A, C, T])).collect();
            reads.push(obs(&cx, &[A, C, A])); // lone variant differs at last col only
            let out = split_readchain_by_psv(&reads, 3, 2, 2);
            assert_eq!(
                out.len(),
                1,
                "lone variant neither >=K-distinct nor >=min_reads"
            );
            assert_eq!(out[0].read_count, 6);
        }

        #[test]
        fn below_k_columns_merges() {
            let cx = chain_x();
            let mut reads: Vec<ReadObs> = (0..3).map(|_| obs(&cx, &[A, C, T])).collect();
            reads.extend((0..3).map(|_| obs(&cx, &[A, C, A]))); // differ at exactly 1 column
            let out = split_readchain_by_psv(&reads, 3, 2, 2);
            assert_eq!(out.len(), 1, "differ at < K columns -> merge");
            assert!(!out[0].identifiable, "the indistinguishability wall");
        }

        #[test]
        fn ambiguous_read_not_forced() {
            let cx = chain_x();
            let mut reads = Vec::new();
            for _ in 0..3 {
                reads.push(obs(&cx, &[A, C, T]));
            }
            for _ in 0..3 {
                reads.push(obs(&cx, &[G, T, A]));
            }
            for _ in 0..2 {
                reads.push(obs(&cx, &[N, N, N])); // fully ambiguous
            }
            let out = split_readchain_by_psv(&reads, 3, 2, 2);
            assert_eq!(out.len(), 2, "all-None reads create no third copy");
            assert!(out.iter().all(|c| c.identifiable));
        }

        // ---- Adversarial cases added by review ----

        /// Indistinguishability wall: zero PSV columns => no allele information at all =>
        /// exactly one merged, non-identifiable chain group regardless of read count.
        #[test]
        fn adv_no_psv_columns_merges_unidentifiable() {
            let cx = chain_x();
            // Reads carry empty allele vectors (n_psv_columns == 0).
            let reads: Vec<ReadObs> = (0..8).map(|_| obs(&cx, &[])).collect();
            let out = split_readchain_by_psv(&reads, 0, 2, 2);
            assert_eq!(out.len(), 1, "no PSV info => single merged copy");
            assert!(!out[0].identifiable, "n_psv=0 can never be identifiable");
            assert_eq!(
                out[0].read_count, 8,
                "all reads counted in the merged group"
            );
            assert!(
                out[0].allele_vector.is_empty(),
                "no columns => empty haplotype"
            );
        }

        /// A minority allele group with < min_reads_per_copy support must NOT become its own
        /// copy. Here [G,T,A] has only 1 supporting read (< min_reads=2) while [A,C,T] has 3.
        /// Result must be a single merged group, and the lone read folded into its count.
        #[test]
        fn adv_minority_below_min_reads_not_a_copy() {
            let cx = chain_x();
            let mut reads: Vec<ReadObs> = (0..3).map(|_| obs(&cx, &[A, C, T])).collect();
            reads.push(obs(&cx, &[G, T, A])); // only 1 read for a >=K-distinct vector
            let out = split_readchain_by_psv(&reads, 3, 2, 2);
            assert_eq!(
                out.len(),
                1,
                "lone strong-but-undersupported variant is not a copy"
            );
            assert!(
                !out[0].identifiable,
                "only one candidate copy => not identifiable"
            );
            assert_eq!(
                out[0].read_count, 4,
                "minority read folded into the merged count"
            );
        }

        /// Apportionment: a partially-observed read consistent with exactly one copy is
        /// assigned to it; a fully-ambiguous read is left shared (no copy inflation).
        #[test]
        fn adv_partial_read_apportioned_uniquely() {
            let cx = chain_x();
            let mut reads = Vec::new();
            for _ in 0..2 {
                reads.push(obs(&cx, &[A, C, T]));
            }
            for _ in 0..2 {
                reads.push(obs(&cx, &[G, T, A]));
            }
            // [A,N,N] disagrees with [G,T,A] at col0 -> consistent ONLY with [A,C,T].
            reads.push(obs(&cx, &[A, N, N]));
            // [N,N,N] consistent with both -> shared, counts toward neither copy.
            reads.push(obs(&cx, &[N, N, N]));
            let out = split_readchain_by_psv(&reads, 3, 2, 2);
            assert_eq!(out.len(), 2, "two identifiable copies");
            assert!(out.iter().all(|c| c.identifiable));
            let total: usize = out.iter().map(|c| c.read_count).sum();
            // 2 + 2 + 1 apportioned = 5; the all-None read is shared and excluded.
            assert_eq!(total, 5, "ambiguous read must not inflate any copy");
            // The copy whose col0==A must have picked up the [A,N,N] read (3), the other 2.
            let mut counts: Vec<usize> = out.iter().map(|c| c.read_count).collect();
            counts.sort();
            assert_eq!(
                counts,
                vec![2, 3],
                "partial read joined the unique consistent copy"
            );
        }

        /// Determinism / order-independence: shuffling the input read order yields an
        /// identical output vector (same order, same counts).
        #[test]
        fn adv_deterministic_under_reorder() {
            let cx = chain_x();
            let mut reads = Vec::new();
            for _ in 0..3 {
                reads.push(obs(&cx, &[A, A, A]));
            }
            for _ in 0..3 {
                reads.push(obs(&cx, &[C, C, C]));
            }
            for _ in 0..3 {
                reads.push(obs(&cx, &[G, G, G]));
            }
            let out1 = split_readchain_by_psv(&reads, 3, 2, 2);
            reads.reverse();
            let out2 = split_readchain_by_psv(&reads, 3, 2, 2);
            assert_eq!(out1, out2, "output independent of input read order");
        }

        #[test]
        fn three_copies() {
            let cx = chain_x();
            let mut reads = Vec::new();
            for _ in 0..3 {
                reads.push(obs(&cx, &[A, A, A]));
            }
            for _ in 0..3 {
                reads.push(obs(&cx, &[C, C, C]));
            }
            for _ in 0..3 {
                reads.push(obs(&cx, &[G, G, G]));
            }
            let out = split_readchain_by_psv(&reads, 3, 2, 2);
            assert_eq!(out.len(), 3, "three pairwise >=K-distinct copies");
            assert!(out.iter().all(|c| c.identifiable));
        }

        // ---- ReadObs BRIDGE tests (CIGAR-walking allele_at is the load-bearing primitive) ----

        fn aligned(ref_start: u64, cigar: &[(char, u64)], seq: &[u8]) -> AlignedRead {
            AlignedRead {
                ref_start,
                cigar: cigar.to_vec(),
                seq: seq.to_vec(),
                qual: vec![],
            }
        }

        /// "10M": straight match. Interior ref position returns the right base; before the
        /// start and past the end return None.
        #[test]
        fn bridge_allele_at_simple_match() {
            // ref 100..110 <-> seq[0..10]. seq = ACGTACGTAC
            let read = aligned(100, &[('M', 10)], b"ACGTACGTAC");
            assert_eq!(
                allele_at(&read, 100),
                Some(b'A'),
                "first matched ref pos -> seq[0]"
            );
            assert_eq!(
                allele_at(&read, 103),
                Some(b'T'),
                "interior pos ref 103 -> seq[3]"
            );
            assert_eq!(
                allele_at(&read, 109),
                Some(b'C'),
                "last matched pos -> seq[9]"
            );
            assert_eq!(allele_at(&read, 99), None, "before ref_start -> None");
            assert_eq!(allele_at(&read, 110), None, "past the end -> None");
            assert_eq!(allele_at(&read, 200), None, "far past the end -> None");
        }

        /// "5M100N5M": a spliced read. Positions inside the N intron return None; positions in
        /// the second block return the correct base (ref advanced past the intron, read did NOT).
        #[test]
        fn bridge_allele_at_spliced_intron() {
            // block1: ref 100..105 <-> seq[0..5] = "AAAAA"
            // intron: ref 105..205 (N100), consumes no read
            // block2: ref 205..210 <-> seq[5..10] = "CCCGT"
            let read = aligned(100, &[('M', 5), ('N', 100), ('M', 5)], b"AAAAACCCGT");
            // exon1
            assert_eq!(allele_at(&read, 100), Some(b'A'));
            assert_eq!(
                allele_at(&read, 104),
                Some(b'A'),
                "last base of block1 -> seq[4]"
            );
            // inside intron
            assert_eq!(allele_at(&read, 105), None, "intron start -> None");
            assert_eq!(allele_at(&read, 150), None, "deep in intron -> None");
            assert_eq!(allele_at(&read, 204), None, "last intron pos -> None");
            // exon2: ref 205 -> seq[5], NOT seq[105]; read coord did not advance over the intron.
            assert_eq!(allele_at(&read, 205), Some(b'C'), "block2 start -> seq[5]");
            assert_eq!(
                allele_at(&read, 208),
                Some(b'G'),
                "block2 interior ref 208 -> seq[8]"
            );
            assert_eq!(allele_at(&read, 209), Some(b'T'), "block2 end -> seq[9]");
            assert_eq!(allele_at(&read, 210), None, "past read -> None");
        }

        /// "3M2D3M": a deletion. Positions inside the D return None; positions after map correctly
        /// (ref advanced over the deletion, read did NOT).
        #[test]
        fn bridge_allele_at_deletion() {
            // block1: ref 50..53 <-> seq[0..3] = "GGG"
            // del:    ref 53..55 (D2), consumes no read
            // block2: ref 55..58 <-> seq[3..6] = "TAC"
            let read = aligned(50, &[('M', 3), ('D', 2), ('M', 3)], b"GGGTAC");
            assert_eq!(allele_at(&read, 50), Some(b'G'));
            assert_eq!(
                allele_at(&read, 52),
                Some(b'G'),
                "last base before del -> seq[2]"
            );
            assert_eq!(allele_at(&read, 53), None, "first deleted ref pos -> None");
            assert_eq!(allele_at(&read, 54), None, "second deleted ref pos -> None");
            assert_eq!(
                allele_at(&read, 55),
                Some(b'T'),
                "after del ref 55 -> seq[3]"
            );
            assert_eq!(
                allele_at(&read, 57),
                Some(b'C'),
                "after del ref 57 -> seq[5]"
            );
            assert_eq!(allele_at(&read, 58), None, "past read -> None");
        }

        /// "4S10M": softclip consumes read seq but not ref. The first matched ref position
        /// returns seq[4] (the 4 soft-clipped bases are skipped on the read axis).
        #[test]
        fn bridge_allele_at_softclip() {
            // softclip: seq[0..4] = "XXXX" (not aligned)
            // block:    ref 1000..1010 <-> seq[4..14] = "ACGTACGTAC"
            let read = aligned(1000, &[('S', 4), ('M', 10)], b"XXXXACGTACGTAC");
            assert_eq!(
                allele_at(&read, 1000),
                Some(b'A'),
                "first matched ref -> seq[4]"
            );
            assert_eq!(allele_at(&read, 1001), Some(b'C'), "-> seq[5]");
            assert_eq!(
                allele_at(&read, 1009),
                Some(b'C'),
                "last matched -> seq[13]"
            );
            assert_eq!(allele_at(&read, 999), None, "before ref_start -> None");
            assert_eq!(allele_at(&read, 1010), None, "past the end -> None");
        }

        /// "5M2I5M": an insertion. Ref positions are continuous across the insertion; the base
        /// after the insertion uses an advanced read coordinate (read += inserted bases).
        #[test]
        fn bridge_allele_at_insertion() {
            // block1: ref 10..15 <-> seq[0..5] = "AAAAA"
            // ins:    seq[5..7] = "II" (consumes read, no ref)
            // block2: ref 15..20 <-> seq[7..12] = "CCCCG"
            let read = aligned(10, &[('M', 5), ('I', 2), ('M', 5)], b"AAAAAIICCCCG");
            assert_eq!(
                allele_at(&read, 14),
                Some(b'A'),
                "last base of block1 -> seq[4]"
            );
            // ref is continuous: 15 is the first base of block2, read coord jumped over the 2 ins.
            assert_eq!(
                allele_at(&read, 15),
                Some(b'C'),
                "after insertion ref 15 -> seq[7]"
            );
            assert_eq!(
                allele_at(&read, 19),
                Some(b'G'),
                "block2 end ref 19 -> seq[11]"
            );
            assert_eq!(allele_at(&read, 20), None, "past read -> None");
        }

        /// intron_chain_of("5M100N5M") == [(ref_start+5, ref_start+105)].
        #[test]
        fn bridge_intron_chain_single() {
            let read = aligned(100, &[('M', 5), ('N', 100), ('M', 5)], b"AAAAACCCCC");
            assert_eq!(intron_chain_of(&read), vec![(105, 205)]);
        }

        // ---- Task-2: parallel-vector invariant (PSV positions || allele_vector) ----

        #[test]
        fn split_and_positions_are_parallel_vectors() {
            // Build N AlignedReads at one locus with two co-varying alleles at 2 positions; assert the
            // discovered PSV positions length == each emitted CopyIsoform.allele_vector length.
            let reads = make_two_copy_locus_reads(); // helper: 6+ reads, A/C split at 2 cols
            let pos = discover_locus_psvs(&reads, 3);
            let copies = split_locus_copies(&reads, 3, 2, 3);
            assert!(copies.len() >= 2, "two identifiable copies");
            for c in &copies {
                assert_eq!(
                    c.allele_vector.len(),
                    pos.len(),
                    "allele_vector parallel to discovered positions"
                );
            }
        }

        /// Two interleaved haplotypes: A-reads and C-reads, each carrying allele at genome 130 and 160,
        /// with a constant background of G elsewhere.  8 reads total (4+4), ≥ min_allele_reads=3 each.
        fn make_two_copy_locus_reads() -> Vec<AlignedRead> {
            let mut reads = Vec::new();
            for base in [b'A', b'C'] {
                for _ in 0..4 {
                    let mut s = vec![b'G'; 100];
                    s[30] = base; // genome pos 100+30 = 130
                    s[60] = base; // genome pos 100+60 = 160
                    reads.push(pile_read(100, &s));
                }
            }
            reads
        }

        /// intron_chain over two introns, with a deletion that must NOT split the chain.
        #[test]
        fn bridge_intron_chain_multi_with_deletion() {
            // ref 100: 3M (100..103), 50N (103..153), 2M(153..155), 1D(155..156),
            //          2M(156..158), 20N(158..178), 4M(178..182)
            let read = aligned(
                100,
                &[
                    ('M', 3),
                    ('N', 50),
                    ('M', 2),
                    ('D', 1),
                    ('M', 2),
                    ('N', 20),
                    ('M', 4),
                ],
                b"AAACCGGTTTT",
            );
            assert_eq!(
                intron_chain_of(&read),
                vec![(103, 153), (158, 178)],
                "two introns; the deletion does not break the chain"
            );
        }

        /// End-to-end build_read_obs: a spliced read + 3 PSV positions (exon1 / intron / exon2)
        /// -> ReadObs with the right intron_chain and psv_alleles = [Some, None, Some].
        #[test]
        fn bridge_build_read_obs_end_to_end() {
            // block1: ref 100..105 <-> seq[0..5] = "ACGTA"
            // intron: ref 105..205 (N100)
            // block2: ref 205..210 <-> seq[5..10] = "TTGCA"
            let read = aligned(100, &[('M', 5), ('N', 100), ('M', 5)], b"ACGTATTGCA");
            // PSV positions: one in exon1 (102 -> seq[2]='G'), one in the intron (150 -> None),
            // one in exon2 (207 -> seq[7]='G').
            let psv = [102u64, 150u64, 207u64];
            let obs = build_read_obs(&read, &psv);
            assert_eq!(
                obs.intron_chain,
                vec![(105, 205)],
                "intron chain from the N op"
            );
            assert_eq!(
                obs.psv_alleles,
                vec![Some(b'G'), None, Some(b'G')],
                "[exon1 base, intron None, exon2 base]"
            );
        }

        // ---- Task 3: Discovery discriminators (min_p_distinct, strand_symmetric_spectrum) ----

        #[test]
        fn min_p_distinct_requires_enough_columns() {
            // 1 distinguishing column: min_p = e/3 ≈ 1e-3 -> NOT < alpha=1e-3 -> false
            let cand = vec![Some(b'C'), Some(b'A')];
            let reference = vec![Some(b'A'), Some(b'A')];
            assert!(
                !min_p_distinct(&cand, &reference, 0.003, 1e-3),
                "1 column insufficient at alpha=1e-3"
            );
            // 3 distinguishing columns: min_p ≈ (1e-3)^3 << alpha -> true
            let cand3 = vec![Some(b'C'), Some(b'C'), Some(b'C')];
            let ref3 = vec![Some(b'A'), Some(b'A'), Some(b'A')];
            assert!(
                min_p_distinct(&cand3, &ref3, 0.003, 1e-3),
                "3 columns sufficient"
            );
        }

        #[test]
        fn strand_symmetric_rejects_pure_a_to_g() {
            // every difference is A->G: editing-like -> reject (false)
            let host = vec![Some(b'A'), Some(b'A')];
            let cand = vec![Some(b'G'), Some(b'G')];
            assert!(
                !strand_symmetric_spectrum(&host, &cand),
                "pure A->G is editing-like"
            );
            // mixed spectrum (A->C, C->T): real divergence -> accept (true)
            let host2 = vec![Some(b'A'), Some(b'C')];
            let cand2 = vec![Some(b'C'), Some(b'T')];
            assert!(
                strand_symmetric_spectrum(&host2, &cand2),
                "mixed spectrum is divergence-like"
            );
        }
    }
}

pub mod single_copy {
    //! Single-copy baseline: the χ(H)=1 loci that calibrate copy number.
    //!
    //! A TRANSCRIPT is the shipped `DenovoTranscript` object, defined by the `assemble_gate` predicate
    //! (`denovo_assemble.rs`): an exact-intron-chain read cluster whose LOCUS (junction-incidence component,
    //! `locus_support`) carries >= GATE_MIN_READS reads, whose junctions are all canonical and consistent-strand,
    //! and whose spliced length is in [MIN_SPLICED, MAX_SPLICED]. A locus is a copy family of size >= 1:
    //! single-copy is the degenerate χ(H)=1 case (0 PSVs); its transcripts are the λ_global baseline that
    //! calibrates `depth_cn = E_fam / λ_global`. This module is the copy-number baseline, NOT an isoform catalog.
    //!
    //! **STATUS:** OPT-IN — --single-copy-baseline (src/bin/gw_family_catalog.rs:189-190, `#[arg(long, default_value_t = false)] single_copy_baseline: bool`).  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    use crate::family::denovo_pipeline::ColocatedFamily;
    use crate::family::family_detect::DenovoTranscript;

    /// One single-copy locus (χ(H)=1). Carries NO `seq` — this is a lightweight baseline record, so a genome-wide
    /// accumulator of ~100k of these is cheap.
    #[derive(Clone, Debug, PartialEq)]
    pub struct SingleCopyLocus {
        pub chrom: String,
        pub start: u64,
        pub end: u64,
        pub strand: char,
        pub n_reads: u32,
        pub n_exons: usize,
    }

    /// Reps that NO family claims, as lightweight `SingleCopyLocus` records. Membership is by `tid`.
    ///
    /// `rep_totals[i]` is the LOCUS TOTAL read count for `reps[i]` — the summed reads over every isoform that
    /// collapsed into that rep's locus (from `collapse_loci_span_aware_with_totals`), NOT the rep isoform's own
    /// count. This is the expression basis λ_global needs: a gene's total, so it is on the same footing as the
    /// family-total `E_fam` in `depth_cn = E_fam / λ_global`.
    pub fn single_copy_loci(
        reps: &[DenovoTranscript],
        rep_totals: &[u32],
        families: &[ColocatedFamily],
    ) -> Vec<SingleCopyLocus> {
        assert_eq!(
            reps.len(),
            rep_totals.len(),
            "rep_totals must be parallel to reps"
        );
        let claimed: std::collections::HashSet<&str> = families
            .iter()
            .flat_map(|f| f.copies.iter().map(|c| c.tid.as_str()))
            .collect();
        reps.iter()
            .zip(rep_totals.iter())
            .filter(|(r, _)| !claimed.contains(r.tid.as_str()))
            .map(|(r, &total)| SingleCopyLocus {
                chrom: r.chrom.clone(),
                start: r.start,
                end: r.end,
                strand: r.strand,
                n_reads: total,
                n_exons: r.introns.len() + 1,
            })
            .collect()
    }

    /// Median read count over the single-copy loci — the genome-wide single-copy expression floor `λ_global`.
    /// `None` on an empty set.
    pub fn lambda_global(loci: &[SingleCopyLocus]) -> Option<f64> {
        if loci.is_empty() {
            return None;
        }
        let mut v: Vec<u32> = loci.iter().map(|l| l.n_reads).collect();
        v.sort_unstable();
        let n = v.len();
        Some(if n % 2 == 1 {
            v[n / 2] as f64
        } else {
            (v[n / 2 - 1] as f64 + v[n / 2] as f64) / 2.0
        })
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn tx(
            tid: &str,
            chrom: &str,
            start: u64,
            end: u64,
            n_reads: u32,
            strand: char,
            introns: Vec<(u64, u64)>,
        ) -> DenovoTranscript {
            DenovoTranscript {
                tid: tid.into(),
                chrom: chrom.into(),
                start,
                end,
                n_reads,
                strand,
                introns,
                seq: vec![],
                ..Default::default()
            }
        }
        fn fam(id: &str, copies: Vec<DenovoTranscript>) -> ColocatedFamily {
            let chrom = copies[0].chrom.clone();
            ColocatedFamily {
                family_id: id.into(),
                chrom,
                start: 0,
                end: 0,
                copies,
            }
        }

        #[test]
        fn single_copy_loci_are_the_reps_no_family_claims() {
            let a = tx("a", "c1", 100, 200, 10, '+', vec![(120, 150)]);
            let b = tx("b", "c1", 300, 400, 20, '-', vec![]);
            let c = tx("c", "c1", 500, 600, 30, '+', vec![(520, 540), (560, 580)]);
            let families = vec![fam("F0", vec![a.clone(), b.clone()])];
            let sc = single_copy_loci(&[a, b, c], &[10, 20, 30], &families);
            assert_eq!(sc.len(), 1);
            assert_eq!(sc[0].chrom, "c1");
            assert_eq!(sc[0].start, 500);
            assert_eq!(sc[0].end, 600);
            assert_eq!(sc[0].strand, '+');
            assert_eq!(sc[0].n_reads, 30);
            assert_eq!(sc[0].n_exons, 3, "n_exons = introns.len() + 1");
        }

        #[test]
        fn single_copy_locus_n_reads_is_the_locus_total_not_the_rep_isoform() {
            // A single-copy locus whose rep isoform holds 30 reads but whose locus (all collapsed isoforms) totals
            // 62. n_reads must be the total, matching the E_fam basis in depth_cn = E_fam / lambda_global.
            let c = tx("c", "c1", 500, 600, 30, '+', vec![(520, 540)]);
            let sc = single_copy_loci(&[c], &[62], &[]);
            assert_eq!(sc.len(), 1);
            assert_eq!(
                sc[0].n_reads, 62,
                "n_reads = locus total (summed isoforms), not the rep's 30"
            );
        }

        #[test]
        fn single_copy_loci_membership_is_by_tid() {
            let a = tx("a", "c1", 100, 200, 10, '+', vec![]);
            let claimed = tx("a", "c1", 100, 200, 10, '+', vec![]);
            let other = tx("z", "c1", 100, 200, 10, '+', vec![]);
            let families = vec![fam("F0", vec![claimed, other])];
            assert!(
                single_copy_loci(&[a], &[10], &families).is_empty(),
                "a's tid is claimed -> not single-copy"
            );
        }

        #[test]
        fn single_copy_locus_carries_no_seq() {
            let a = tx("a", "c1", 1, 2, 5, '+', vec![]);
            let sc = single_copy_loci(&[a], &[5], &[]);
            let _ = SingleCopyLocus {
                chrom: sc[0].chrom.clone(),
                start: sc[0].start,
                end: sc[0].end,
                strand: sc[0].strand,
                n_reads: sc[0].n_reads,
                n_exons: sc[0].n_exons,
            };
        }

        #[test]
        fn lambda_global_is_the_median_n_reads() {
            let mk = |n: u32| SingleCopyLocus {
                chrom: "c".into(),
                start: 0,
                end: 1,
                strand: '+',
                n_reads: n,
                n_exons: 1,
            };
            assert_eq!(
                lambda_global(&[mk(10), mk(20), mk(30), mk(40)]),
                Some(25.0),
                "even -> mean of middle two"
            );
            assert_eq!(
                lambda_global(&[mk(30), mk(10), mk(20)]),
                Some(20.0),
                "odd, unsorted -> middle after sort"
            );
            assert_eq!(lambda_global(&[]), None, "empty -> None");
            assert_eq!(lambda_global(&[mk(7)]), Some(7.0));
        }
    }
}

pub mod readonly_copy_number {
    //! Reference-free per-family copy number (Task R1 of the reference-free copy-number plan).
    //!
    //! The advisor's key idea: even reads that can't be ASSIGNED to a specific copy still INDICATE
    //! copy number. Two reference-free legs (no genome/assembly needed):
    //!   * `chi_h`    -- distinguishable copies from the PSV conflict structure (Rust port of
    //!                   `bench/family_copy_number.py`'s `copyonly_K`: the number of distinct
    //!                   pairwise-conflicting hap-vectors). A LOWER bound -- it collapses identical/
    //!                   compatible copies into one group.
    //!                   families e.g. `chi_H=1` on a locus whose true copy number is ~11.
    //!   * `depth_cn` -- read-depth estimate (E_fam / lambda_global) that recovers identical/
    //!                   collapsed copies `chi_h` misses, because it counts ALL family reads
    //!                   (including the unassignable ones) rather than distinct hap-vectors.
    //!
    //! Both legs are genome-free: `chi_h` only needs the per-copy PSV allele vectors already
    //! computed for assignment, and `depth_cn` only needs a read count and an externally-supplied
    //! single-copy expression floor (`lambda_global`, itself an RNA-only quantity -- see
    //! `bench/rna_copy_number_depth.py::global_single_copy_anchor`).
    //!
    //! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    /// Number of DISTINCT pairwise-conflicting copy hap-vectors (chi(H) / MCC, THEORY.md Lemma 1).
    ///
    /// Two copies CONFLICT iff they differ at >= 1 PSV column where BOTH have an observed allele
    /// (`Some(_)`); columns where either copy is `None` (no evidence) never contribute a conflict.
    /// Copies are grouped so that within a group no two members conflict; `chi_h` is the resulting
    /// group count. Implemented greedily (matches `copyonly_K`'s semantics: a family's colors are the
    /// distinct pairwise-conflicting hap-vectors) -- place each copy in the first existing group none
    /// of whose members it conflicts with, else start a new group.
    pub fn chi_h(copy_alleles: &[Vec<Option<u8>>]) -> usize {
        fn conflicts(a: &[Option<u8>], b: &[Option<u8>]) -> bool {
            a.iter()
                .zip(b.iter())
                .any(|(x, y)| matches!((x, y), (Some(xa), Some(yb)) if xa != yb))
        }

        let mut groups: Vec<Vec<&Vec<Option<u8>>>> = Vec::new();
        for alleles in copy_alleles {
            let home = groups
                .iter()
                .position(|g| g.iter().all(|m| !conflicts(m, alleles)));
            match home {
                Some(gi) => groups[gi].push(alleles),
                None => groups.push(vec![alleles]),
            }
        }
        groups.len()
    }

    /// `chi_h`, but two copies also conflict when BOTH carry a non-empty, DIFFERING copy-specific junction set.
    ///
    /// The reference-free copy COUNT must see the same evidence as the copy ASSIGNMENT it certifies. The gate and
    /// the EM assign on PSVs **and** copy-specific junctions; `chi_h` sees only PSVs. On the recovered DAZ family
    /// that contradiction is visible in the output: `n_copies = 2`, O2 assigns 2213 of 2353 placements across the
    /// two copies, and yet `chi_H = 1`, so `famcn_readonly` reports **1 copy** for a two-copy family. DAZ1 and DAZ2
    /// are near-identical exonically (one PSV column) and are separated by their junction structure.
    ///
    /// The `both non-empty` guard mirrors the PSV rule's `both Some`: absence of junction evidence is not evidence
    /// of sameness, so `chi_H` remains a valid LOWER bound (`max(n_loci, chi_H) <= true`). With no junction
    /// evidence at all this is exactly [`chi_h`].
    pub fn chi_h_with_junctions(
        copy_alleles: &[Vec<Option<u8>>],
        copy_junctions: &[Vec<i64>],
    ) -> usize {
        fn psv_conflict(a: &[Option<u8>], b: &[Option<u8>]) -> bool {
            a.iter()
                .zip(b.iter())
                .any(|(x, y)| matches!((x, y), (Some(xa), Some(yb)) if xa != yb))
        }
        fn junction_conflict(a: &[i64], b: &[i64]) -> bool {
            if a.is_empty() || b.is_empty() {
                return false; // no evidence is not evidence of difference
            }
            let (sa, sb): (std::collections::BTreeSet<_>, std::collections::BTreeSet<_>) =
                (a.iter().collect(), b.iter().collect());
            sa != sb
        }
        let empty: Vec<i64> = Vec::new();
        let junc = |i: usize| copy_junctions.get(i).unwrap_or(&empty);

        let mut groups: Vec<Vec<usize>> = Vec::new();
        for i in 0..copy_alleles.len() {
            let home = groups.iter().position(|g| {
                g.iter().all(|&m| {
                    !psv_conflict(&copy_alleles[m], &copy_alleles[i])
                        && !junction_conflict(junc(m), junc(i))
                })
            });
            match home {
                Some(gi) => groups[gi].push(i),
                None => groups.push(vec![i]),
            }
        }
        groups.len()
    }

    /// Reference-free read-depth copy-number estimate: `E_fam / lambda_global`.
    ///
    /// `e_fam` = total reads over the family (E_fam, includes unassignable reads); `lambda_global` =
    /// the genome-wide RNA single-copy expression floor (median n_reads over single-copy transcripts,
    /// precomputed by `bench/rna_copy_number_depth.py::global_single_copy_anchor` -- a genome-wide RNA
    /// quantity, NOT genomic). Returns `NaN` when `lambda_global <= 0` (undefined / not supplied).
    pub fn depth_cn(e_fam: usize, lambda_global: f64) -> f64 {
        if lambda_global <= 0.0 {
            return f64::NAN;
        }
        e_fam as f64 / lambda_global
    }

    #[cfg(test)]
    mod chi_h_junction_tests {
        use super::*;

        /// The recovered DAZ family: two copies that are exonically near-identical (one PSV column, PSV-compatible)
        /// but carry different junction structures. chi_h alone reports 1 for a family O2 assigns across 2 copies.
        #[test]
        fn chi_h_with_junctions_separates_psv_compatible_junction_distinct_copies() {
            let alleles = vec![vec![None], vec![None]];
            assert_eq!(
                chi_h(&alleles),
                1,
                "PSV-only evidence cannot see the second copy"
            );
            assert_eq!(
                chi_h_with_junctions(&alleles, &[vec![0i64], vec![7i64]]),
                2,
                "junctions separate them"
            );
        }

        /// Absence of junction evidence is not evidence of sameness -- chi_H stays a LOWER bound.
        #[test]
        fn chi_h_with_junctions_without_junction_evidence_is_exactly_chi_h() {
            let alleles = vec![vec![Some(b'A')], vec![Some(b'C')], vec![None]];
            assert_eq!(
                chi_h_with_junctions(&alleles, &[vec![], vec![], vec![]]),
                chi_h(&alleles)
            );
            // one copy has junctions, the other does not => no junction conflict is inferred
            assert_eq!(
                chi_h_with_junctions(&[vec![None], vec![None]], &[vec![3i64], vec![]]),
                1
            );
        }

        /// Identical junction sets do not manufacture a conflict, whatever the ordering.
        #[test]
        fn chi_h_with_junctions_ignores_order_and_identical_sets() {
            let alleles = vec![vec![None], vec![None]];
            assert_eq!(
                chi_h_with_junctions(&alleles, &[vec![5i64, 1], vec![1i64, 5]]),
                1
            );
        }

        /// A PSV conflict still separates copies even when their junctions agree.
        #[test]
        fn chi_h_with_junctions_still_honours_psv_conflicts() {
            let alleles = vec![vec![Some(b'A')], vec![Some(b'G')]];
            assert_eq!(chi_h_with_junctions(&alleles, &[vec![2i64], vec![2i64]]), 2);
        }

        /// It can never exceed the copy count, so it remains a lower bound on true copy number.
        #[test]
        fn chi_h_with_junctions_never_exceeds_n_copies() {
            let alleles = vec![vec![Some(b'A')], vec![Some(b'C')], vec![Some(b'G')]];
            let j = vec![vec![1i64], vec![2i64], vec![3i64]];
            assert!(chi_h_with_junctions(&alleles, &j) <= alleles.len());
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::family::copy_assign::em_copy_assign::{em_assign_family, EmLabel};
        use crate::family::copy_assign::AssignParams;

        /// Task R2 (O1<->O2 harmony pin): `chi_h` (this module, O1's conflict-graph copy COUNT) and
        /// `em_assign_family` (`em_copy_assign`, O2's EM ASSIGNMENT) both consume the SAME per-copy
        /// hap-vector (`copy_alleles` / `CopyProfile.alleles`) -- there is exactly one copy object, not
        /// two independently-tunable ones. A better O1 (every copy pairwise-conflicts, no de-tie needed)
        /// must raise `chi_h` AND let the EM certify every read against that K; collapsing two copies to
        /// an identical hap-vector (a worse/over-merged O1) must drop `chi_h` by exactly 1 AND make the EM
        /// unable to separate reads from that pair -- both legs hit the K-frontier (`SoftZone`) together.
        /// This is a REGRESSION PIN: no new behavior, just nailing down that the two already compose.
        #[test]
        fn o1_o2_share_one_copy_object() {
            let argmax = |row: &Vec<f64>| {
                row.iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .map(|(k, _)| k)
                    .unwrap()
            };
            let params = AssignParams::for_alpha(1e-3);

            // --- Regime 1: distinguishable (good O1) -- 3 copies, each privately marked at its own
            // column against a shared background, so every pair conflicts at 2 of the 3 columns (well
            // clear of the K-frontier).
            let distinguishable: Vec<Vec<Option<u8>>> = vec![
                vec![Some(b'C'), Some(b'A'), Some(b'A')],
                vec![Some(b'A'), Some(b'C'), Some(b'A')],
                vec![Some(b'A'), Some(b'A'), Some(b'C')],
            ];
            assert_eq!(
                chi_h(&distinguishable),
                3,
                "O1: 3 pairwise-conflicting copies -> chi_h = 3"
            );

            // one read per true copy, carrying that copy's exact alleles.
            let reads = distinguishable.clone();
            let result = em_assign_family(&reads, &distinguishable, &[], &[], &params, 1e-6, 500);
            assert_eq!(
                result.abundances.len(),
                distinguishable.len(),
                "O2 (EM) must sum abundance over the SAME K copies O1's chi_h counts"
            );
            for (true_copy, row) in result.posteriors.iter().enumerate() {
                assert_eq!(
                    argmax(row),
                    true_copy,
                    "read {true_copy} must be assigned to its true copy"
                );
            }
            assert!(
                result
                    .labels
                    .iter()
                    .all(|l| matches!(l, EmLabel::Certified)),
                "good O1 (all copies pairwise-conflict) -> the EM certifies every read: {:?}",
                result.labels
            );

            // --- Regime 2: collapsed (worse O1) -- copy 1 and copy 2 forced to the IDENTICAL hap-vector
            // (a de-tie / over-merge failure upstream would produce exactly this).
            let mut collapsed = distinguishable.clone();
            collapsed[2] = collapsed[1].clone();
            assert_eq!(
                chi_h(&collapsed),
                chi_h(&distinguishable) - 1,
                "collapsing 2 copies to one hap-vector must drop chi_h by exactly 1 (3 -> 2)"
            );

            // reads carrying the POST-collapse hap-vectors: one from the untouched copy 0 (control),
            // and one from each of the two now-identical copies.
            let collapsed_reads: Vec<Vec<Option<u8>>> = vec![
                distinguishable[0].clone(),
                collapsed[1].clone(),
                collapsed[2].clone(),
            ];
            let collapsed_result =
                em_assign_family(&collapsed_reads, &collapsed, &[], &[], &params, 1e-6, 500);

            assert!(
                matches!(collapsed_result.labels[0], EmLabel::Certified),
                "the un-collapsed copy must remain identifiable: {:?}",
                collapsed_result.labels
            );
            // reads from EITHER of the two now-identical copies can no longer be separated: the EM hits
            // the K-frontier and abstains (SoftZone) on both, in lockstep with chi_h's drop.
            assert!(
                matches!(collapsed_result.labels[1], EmLabel::SoftZone),
                "read from collapsed copy 1 must go SoftZone: {:?}",
                collapsed_result.labels
            );
            assert!(
                matches!(collapsed_result.labels[2], EmLabel::SoftZone),
                "read from collapsed copy 2 must go SoftZone: {:?}",
                collapsed_result.labels
            );
        }

        /// Task H1 (VG-harmony pin): a copy NOT in the reference genome — admitted by O4's
        /// `absent_copy::admit_candidate` gate as a synthetic `DenovoTranscript` — must be a
        /// first-class copy to BOTH downstream consumers, not just threaded through pipeline
        /// plumbing that happens to compile.
        ///
        /// Wiring trace (verified by inspection; this test simulates the RESULT of that trace rather
        /// than re-running it, since the genome/BAM/minimap2 machinery is exercised by
        /// `absent_copy.rs`'s own hermetic tests):
        /// `admit_candidate_with_remap` (`absent_copy.rs`) returns `Admission::Copy(t)` once a
        /// candidate clears all five gates (cluster count, `min_p_distinct` from host, strand-
        /// symmetric spectrum, placeable overlay, remap identity `< 98%`) — `t` is a synthetic
        /// `DenovoTranscript` whose sequence bakes in the copy's distinguishing alleles.
        /// `denovo_pipeline.rs` (~line 626-654, behind `--absent-copies`) collects every
        /// `Admission::Copy(t)` into `admitted`, then extends the family's copy list (`all_copies`)
        /// with it and re-runs `assign_family_detailed` (Stage-2) over the augmented set — so the
        /// admitted copy sits at a real index alongside the reference copies, indistinguishable in
        /// type from them. From there `assign_family_detailed` -> `build_family_profiles`
        /// (`copy_assign_pipeline.rs`) extracts each copy's per-PSV-column allele vector
        /// (`CopyProfile.alleles`, i.e. `copy_psv_alleles`) for every copy in the set — there is no
        /// branch for "reference" vs "admitted-absent" once a `DenovoTranscript` exists. That shared
        /// `copy_alleles` vector is exactly what both `chi_h` (this module — O1's reference-free
        /// copy-number COUNT) and `em_assign_family` (`em_copy_assign` — O2's read ASSIGNMENT)
        /// consume next. This test starts from that shared post-admission state directly: a
        /// `copy_alleles` with 2 reference copies + 1 admitted-absent copy, and reads carrying the
        /// absent copy's alleles.
        #[test]
        fn absent_copy_is_assigned_and_counted() {
            let argmax = |row: &Vec<f64>| {
                row.iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .map(|(k, _)| k)
                    .unwrap()
            };

            // copy 0, copy 1 = "reference" copies (present in the linear genome); copy 2 = the
            // O4-admitted absent copy (its alleles came from `admit_candidate`'s synthetic
            // DenovoTranscript, not from any genome coordinate). Same private-column-per-copy layout
            // as `o1_o2_share_one_copy_object`: every pair differs at 2 of the 3 columns, well clear
            // of the alpha=1e-3, K=3 Bonferroni bound.
            let copy_alleles: Vec<Vec<Option<u8>>> = vec![
                vec![Some(b'C'), Some(b'A'), Some(b'A')], // ref copy 0
                vec![Some(b'A'), Some(b'C'), Some(b'A')], // ref copy 1
                vec![Some(b'A'), Some(b'A'), Some(b'C')], // absent copy 2 (O4-admitted)
            ];
            let absent_idx = 2;

            // (a) chi_h counts the absent copy as its own color: the reference-free copy number
            // rises to 3, not 2 -- the admitted copy is not silently dropped or merged.
            assert_eq!(
                chi_h(&copy_alleles),
                3,
                "reference-free copy number must count the O4-admitted absent copy as a distinct color"
            );

            // A handful of reads carrying the absent copy's exact alleles (as if minimap2/discovery
            // had routed them to this locus and copy_psv_alleles had extracted this vector for them).
            let absent_reads: Vec<Vec<Option<u8>>> = vec![
                copy_alleles[absent_idx].clone(),
                copy_alleles[absent_idx].clone(),
                copy_alleles[absent_idx].clone(),
            ];
            let params = AssignParams::for_alpha(1e-3);
            let result =
                em_assign_family(&absent_reads, &copy_alleles, &[], &[], &params, 1e-6, 500);

            // (b) the EM assigns those reads to the absent copy, confidently.
            for (i, row) in result.posteriors.iter().enumerate() {
                assert_eq!(
                    argmax(row),
                    absent_idx,
                    "read {i} carrying the absent copy's alleles must be assigned to the absent copy"
                );
                assert!(
                    matches!(result.labels[i], EmLabel::Certified),
                    "read {i} must be Certified, not stuck in the K-frontier soft zone: {:?}",
                    result.labels[i]
                );
            }
            assert!(
                result.abundances[absent_idx] > 0.0,
                "the EM's recovered abundance for the absent copy must be > 0: {:?}",
                result.abundances
            );
        }

        #[test]
        fn chi_h_three_pairwise_distinct_private_alleles() {
            // each copy has a private allele at its own column -> all three pairwise conflict.
            let copies = vec![
                vec![Some(b'A'), Some(b'C'), Some(b'C')],
                vec![Some(b'C'), Some(b'A'), Some(b'C')],
                vec![Some(b'C'), Some(b'C'), Some(b'A')],
            ];
            assert_eq!(chi_h(&copies), 3);
        }

        #[test]
        fn chi_h_three_identical_copies_collapse_to_one() {
            let copies = vec![
                vec![Some(b'A'), Some(b'C')],
                vec![Some(b'A'), Some(b'C')],
                vec![Some(b'A'), Some(b'C')],
            ];
            assert_eq!(chi_h(&copies), 1);
        }

        #[test]
        fn chi_h_mixed_two_groups() {
            // [A,A] vs [A,C]: conflict at column 1 (A vs C). The two [A,C] copies agree everywhere ->
            // share a group. So groups = {[A,A]}, {[A,C], [A,C]} = 2.
            let copies = vec![
                vec![Some(b'A'), Some(b'A')],
                vec![Some(b'A'), Some(b'C')],
                vec![Some(b'A'), Some(b'C')],
            ];
            assert_eq!(chi_h(&copies), 2);
        }

        #[test]
        fn chi_h_empty_is_zero() {
            let copies: Vec<Vec<Option<u8>>> = Vec::new();
            assert_eq!(chi_h(&copies), 0);
        }

        #[test]
        fn chi_h_missing_evidence_never_conflicts() {
            // Columns with None on either side never contribute a conflict; these two copies never
            // disagree where both are observed, so they share a group.
            let copies = vec![
                vec![Some(b'A'), None, Some(b'C')],
                vec![None, Some(b'G'), Some(b'C')],
            ];
            assert_eq!(chi_h(&copies), 1);
        }

        #[test]
        fn depth_cn_basic_ratio() {
            assert!((depth_cn(100, 25.0) - 4.0).abs() < 1e-9);
        }

        #[test]
        fn depth_cn_zero_e_fam_is_zero() {
            assert_eq!(depth_cn(0, 10.0), 0.0);
        }

        #[test]
        fn depth_cn_nonpositive_lambda_is_nan() {
            assert!(depth_cn(50, 0.0).is_nan());
            assert!(depth_cn(50, -1.0).is_nan());
        }
    }
}

pub mod linearize {
    //! Augment-and-linearize certificate for reference-absent copies (v2.1). Pure statistics + a deterministic
    //! dinucleotide-preserving decoy; the minimap2 re-alignment is injected as a closure so this module is testable
    //! without any subprocess.
    //!
    //! **STATUS:** OPT-IN — --linearize / --linearize-gate, both `#[arg(long, default_value_t = false)]` at src/bin/copy_assign.rs:280-281 and :289-290; both additionally require  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    /// Deterministic LCG for reproducible shuffles (no external rng crate; Date/rand not needed).
    struct Lcg(u64);
    impl Lcg {
        fn new(seed: u64) -> Self {
            Lcg(seed ^ 0x9E37_79B9_7F4A_7C15)
        }
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0
        }
        fn below(&mut self, n: usize) -> usize {
            if n <= 1 {
                0
            } else {
                (self.next() >> 33) as usize % n
            }
        }
    }
    fn fisher_yates<T>(v: &mut [T], rng: &mut Lcg) {
        for i in (1..v.len()).rev() {
            let j = rng.below(i + 1);
            v.swap(i, j);
        }
    }

    /// Altschul-Erikson (1985) dinucleotide-preserving shuffle via a random Eulerian path through the
    /// nucleotide-transition multigraph. Same length, first/last base, and exact dinucleotide counts; scrambled
    /// order. Deterministic in `seed`. Degenerate/short sequences (<3) or no valid Euler path -> a copy (a decoy
    /// equal to the real candidate is conservative: it inflates the decoy fraction, never the real-minus-decoy gap).
    pub fn dinucleotide_shuffle(seq: &[u8], seed: u64) -> Vec<u8> {
        use std::collections::BTreeMap;
        let n = seq.len();
        if n < 3 {
            return seq.to_vec();
        }
        let last = seq[n - 1];
        let mut base_edges: BTreeMap<u8, Vec<u8>> = BTreeMap::new();
        for w in seq.windows(2) {
            base_edges.entry(w[0]).or_default().push(w[1]);
        }
        let nodes: Vec<u8> = base_edges.keys().copied().collect();
        let mut rng = Lcg::new(seed);
        for _attempt in 0..1000 {
            // shuffle each node's outgoing edges; the LAST edge of each node (!= terminal) is its "tree" edge.
            let mut edges = base_edges.clone();
            for v in edges.values_mut() {
                fisher_yates(v, &mut rng);
            }
            let tree: BTreeMap<u8, u8> = edges
                .iter()
                .filter(|(&b, _)| b != last)
                .map(|(&b, v)| (b, *v.last().unwrap()))
                .collect();
            // valid iff following tree edges from every node reaches `last` without cycling (arborescence into last)
            let ok = nodes.iter().all(|&start| {
                if start == last {
                    return true;
                }
                let (mut cur, mut steps) = (start, 0usize);
                loop {
                    match tree.get(&cur) {
                        Some(&nx) => {
                            cur = nx;
                            if cur == last {
                                break true;
                            }
                            steps += 1;
                            if steps > nodes.len() {
                                break false;
                            }
                        }
                        None => break false,
                    }
                }
            });
            if !ok {
                continue;
            }
            // traverse: consume each node's edges in order (tree edge, at the end, is used last)
            let mut cursor: BTreeMap<u8, usize> = nodes.iter().map(|&b| (b, 0)).collect();
            let mut out = Vec::with_capacity(n);
            out.push(seq[0]);
            let mut cur = seq[0];
            while out.len() < n {
                let v = match edges.get(&cur) {
                    Some(v) => v,
                    None => break,
                };
                let i = cursor[&cur];
                if i >= v.len() {
                    break;
                }
                let nx = v[i];
                *cursor.get_mut(&cur).unwrap() += 1;
                out.push(nx);
                cur = nx;
            }
            if out.len() == n {
                return out;
            }
        }
        seq.to_vec()
    }

    /// Verdict from linearize_certificate: the candidate linearizes relative to null/decoy shuffles.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub enum Verdict {
        Linearizes,
        Not,
        Undetermined,
    }

    /// Result of the linearize_certificate test: candidate fraction, decoy statistics, and verdict.
    #[derive(Clone, Debug)]
    pub struct LinearizeCertificate {
        pub n_pool: usize,
        pub linearized_frac_real: f64,
        pub mean_frac_decoy: f64,
        pub delta: f64,
        pub perm_p: f64,
        pub verdict: Verdict,
    }

    /// Fraction of reads whose primary hit is the candidate contig with MAPQ > 0.
    fn frac_on_candidate(hits: &[Option<(usize, u32)>], cand_idx: usize) -> f64 {
        if hits.is_empty() {
            return 0.0;
        }
        let k = hits
            .iter()
            .filter(|h| matches!(h, Some((i, q)) if *i == cand_idx && *q > 0))
            .count();
        k as f64 / hits.len() as f64
    }

    /// Test whether a sequence (candidate) linearizes: its primary-with-MAPQ>0 fraction significantly
    /// exceeds the mean of N dinucleotide-shuffled decoys. Pure function; the realign closure is
    /// injected for testability.
    ///
    /// The decoys are the N dinucleotide shuffles ONLY. The reverse complement is NOT a valid decoy
    /// for an alignment-based test — minimap2 is strand-symmetric, so revcomp(candidate) attracts the
    /// same reads as the candidate (opposite strand) and would tie `real`, forcing a false NOT. The
    /// dinucleotide shuffles are the composition-matched control (a shuffle matches neither strand of a
    /// read).
    ///
    /// # Arguments
    /// - `candidate_seq`: the candidate contig (appended as the last element to family_copy_seqs).
    /// - `family_copy_seqs`: the background copies (NOT including the candidate).
    /// - `pool_reads`: the read pool for evaluation.
    /// - `n_decoys`: count of dinucleotide-shuffled decoys.
    /// - `seed`: RNG seed for deterministic shuffles.
    /// - `min_pool`: minimum pool size; if `n_pool < min_pool`, returns Undetermined with NaN fields.
    /// - `alpha`: permutation p-value threshold (e.g., 0.05) for the Linearizes verdict.
    /// - `realign`: injected closure `Fn(refs, reads) -> Vec<Option<(contig_idx, mapq)>>` per read.
    pub fn linearize_certificate(
        candidate_seq: &[u8],
        family_copy_seqs: &[Vec<u8>],
        pool_reads: &[Vec<u8>],
        n_decoys: usize,
        seed: u64,
        min_pool: usize,
        alpha: f64,
        realign: impl Fn(&[Vec<u8>], &[Vec<u8>]) -> Vec<Option<(usize, u32)>>,
    ) -> LinearizeCertificate {
        let n_pool = pool_reads.len();
        if n_pool < min_pool {
            return LinearizeCertificate {
                n_pool,
                linearized_frac_real: f64::NAN,
                mean_frac_decoy: f64::NAN,
                delta: f64::NAN,
                perm_p: f64::NAN,
                verdict: Verdict::Undetermined,
            };
        }

        let cand_idx = family_copy_seqs.len(); // candidate is the LAST contig
        let build = |extra: &[u8]| -> Vec<Vec<u8>> {
            let mut v = family_copy_seqs.to_vec();
            v.push(extra.to_vec());
            v
        };

        // Compute the real candidate's linearized fraction.
        let real = frac_on_candidate(&realign(&build(candidate_seq), pool_reads), cand_idx);

        // Generate decoys: N dinucleotide shuffles (distinct seeds) ONLY. No reverse-complement decoy:
        // minimap2 is strand-symmetric, so revcomp(candidate) would attract the same reads on the
        // opposite strand and always tie `real` -> false NOT. (See fn-level doc.)
        let mut decoy_fracs: Vec<f64> = Vec::with_capacity(n_decoys);
        for d in 0..n_decoys {
            let decoy = dinucleotide_shuffle(candidate_seq, seed.wrapping_add(d as u64 + 1));
            decoy_fracs.push(frac_on_candidate(
                &realign(&build(&decoy), pool_reads),
                cand_idx,
            ));
        }

        // Compute statistics.
        let nd = decoy_fracs.len();
        let mean_decoy = decoy_fracs.iter().sum::<f64>() / nd as f64;
        let n_ge = decoy_fracs.iter().filter(|&&d| d >= real).count();
        let perm_p = (n_ge as f64 + 1.0) / (nd as f64 + 1.0);
        let verdict = if perm_p < alpha {
            Verdict::Linearizes
        } else {
            Verdict::Not
        };

        LinearizeCertificate {
            n_pool,
            linearized_frac_real: real,
            mean_frac_decoy: mean_decoy,
            delta: real - mean_decoy,
            perm_p,
            verdict,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        fn dinuc_counts(s: &[u8]) -> std::collections::BTreeMap<(u8, u8), usize> {
            let mut m = std::collections::BTreeMap::new();
            for w in s.windows(2) {
                *m.entry((w[0], w[1])).or_insert(0) += 1;
            }
            m
        }
        #[test]
        fn dinucleotide_shuffle_preserves_composition_and_is_deterministic() {
            let seq = b"ACGTACGTTTGGCCAAACGTACGTGGGCCCAAATTT";
            let a = dinucleotide_shuffle(seq, 42);
            let b = dinucleotide_shuffle(seq, 42);
            assert_eq!(a, b, "deterministic for a given seed");
            assert_eq!(a.len(), seq.len(), "length preserved");
            assert_eq!(a[0], seq[0], "first base preserved");
            assert_eq!(
                *a.last().unwrap(),
                *seq.last().unwrap(),
                "last base preserved"
            );
            assert_eq!(
                dinuc_counts(&a),
                dinuc_counts(seq),
                "exact dinucleotide counts preserved"
            );
            let c = dinucleotide_shuffle(seq, 43);
            assert_ne!(
                a, c,
                "different seed -> different shuffle (for a shufflable seq)"
            );
        }

        // A fake realign: a read "belongs" to the candidate iff its bytes equal the candidate contig's bytes;
        // then it maps uniquely (mapq 60) to the candidate index; otherwise it maps to copy 0 with mapq 0 (tied).
        fn fake_realign(refs: &[Vec<u8>], reads: &[Vec<u8>]) -> Vec<Option<(usize, u32)>> {
            let cand_idx = refs.len() - 1;
            reads
                .iter()
                .map(|r| {
                    if r == &refs[cand_idx] {
                        Some((cand_idx, 60))
                    } else {
                        Some((0, 0))
                    }
                })
                .collect()
        }

        #[test]
        fn real_copy_linearizes_decoy_does_not() {
            // Use a longer sequence to minimize chance of shuffle returning the original
            let cand = b"ACGTACGTTTGGCCAAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
            let copies = vec![b"TTTTGGGGCCCCAAAAAAAATTTTGGGG".to_vec()];
            // 8 pool reads that ARE the candidate (they linearize on it), 2 that are noise
            let mut pool: Vec<Vec<u8>> = (0..8).map(|_| cand.clone()).collect();
            pool.push(b"NNNNNN".to_vec());
            pool.push(b"NNNNNN".to_vec());
            // n_decoys=20: a permutation test at alpha=0.05 needs >=20 decoys so the perm_p floor
            // 1/(n_decoys+1) is strictly below alpha. (The now-removed RC decoy had been silently
            // supplying this 20th decoy; without it, n_decoys=19 would floor perm_p at 1/20 = 0.05,
            // never strictly < 0.05.)
            let cert = linearize_certificate(&cand, &copies, &pool, 20, 7, 5, 0.05, fake_realign);
            assert!(
                (cert.linearized_frac_real - 0.8).abs() < 1e-9,
                "8/10 land on candidate"
            );
            assert!(
                cert.mean_frac_decoy == 0.0,
                "decoys != candidate bytes -> no read lands on them"
            );
            assert!(cert.delta > 0.5);
            // Decoys are the N=20 dinucleotide shuffles ONLY (no RC decoy). No decoy beats real, so
            // perm_p = (0 + 1) / (n_decoys + 1) = 1/21 < 0.05.
            assert!(
                cert.perm_p <= 1.0 / (20.0 + 1.0) + 1e-9,
                "no decoy beats real -> perm_p = 1/(n_decoys+1)"
            );
            assert!(matches!(cert.verdict, Verdict::Linearizes));
        }

        #[test]
        fn null_candidate_is_not_linearized() {
            let cand = b"ACGTACGTTTGGCCAAACGTACGT".to_vec();
            let copies = vec![b"TTTTGGGGCCCCAAAA".to_vec()];
            let pool: Vec<Vec<u8>> = (0..10).map(|_| b"NNNNNN".to_vec()).collect(); // nothing matches candidate
            let cert = linearize_certificate(&cand, &copies, &pool, 19, 7, 5, 0.05, fake_realign);
            assert_eq!(cert.linearized_frac_real, 0.0);
            assert!(
                cert.perm_p > 0.05,
                "real == decoys (both 0) -> perm_p large"
            );
            assert!(matches!(cert.verdict, Verdict::Not));
        }

        #[test]
        fn small_pool_is_undetermined() {
            let cand = b"ACGT".to_vec();
            let pool = vec![cand.clone(), cand.clone()];
            let cert = linearize_certificate(&cand, &[], &pool, 19, 7, 5, 0.05, fake_realign);
            assert!(matches!(cert.verdict, Verdict::Undetermined));
        }
    }
}

pub mod vg_realign {
    //! VG re-align supplement -- re-align poor-fit/unmapped reads to O1's family copy-paths,
    //! significance-gated (correct + discover). Task 1: candidate selection. Task 3: re-align a
    //! candidate read to the family's copy-paths (identity-based, DRY on `seq_utils::aln_id`).
    //! Task 4: gate the re-align correction behind a min_p significance certificate (same `epsilon^delta` form as
    //! `copy_assign::read_copy_evidence`), and greedily pool reads that fit no existing copy into
    //! candidate novel-copy clusters. Task 5: wire Tasks 1/3/4, behind `DenovoConfig::vg_realign`, into the
    //! pipeline where it FEEDS BACK:
    //! `apply_realign_patch` CORRECTS per-read copy assignments (re-thread the hard read through the copy-paths,
    //! take the best-fitting path, same epsilon^delta significance certificate as the PSV gate); `admit_novel_pools`
    //! may ADMIT novel-read clusters as new copies (widening the roster — the O4-frontier leg); then the EM copy
    //! abundance is recomputed (denovo_pipeline.rs ~1664). Default OFF => every output byte-identical. (The
    //! `<out>.vg_realign.tsv` dump is a separate, additive report; "report-only" referred only to that file.)
    //!
    //! VG re-align END-TO-END plan, Task 1: `align_traceback` + `path_obs_at`. There is no `edlib`
    //! crate; `seq_utils::hw_distance` is a hand-rolled 2-row DP that gives the HW/infix edit
    //! DISTANCE only, no alignment path. To re-extract a read's base at a copy-path's PSV columns
    //! (follow-up (c) in `bench/VG_REALIGN.md`) we need the actual traceback, so this keeps a full DP
    //! + backtrack matrix (not the rolling 2-row form) and reconstructs the aligned columns.
    //!
    //! **STATUS:** OPT-IN — --vg-realign or --vg-realign-correct (src/bin/copy_assign.rs:340-341 and :346-347, both `default_value_t = false`; combined into cfg.vg_realign at cop  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    use crate::family::denovo_assemble::BamRead;
    use crate::family::family_detect::DenovoTranscript;
    use crate::family::seq_utils::{aln_id, hw_distance, revcomp_keep_case};

    /// Backtrack pointer for one DP cell of `align_traceback`.
    #[derive(Clone, Copy, PartialEq, Eq, Debug)]
    enum Trace {
        /// Row 0 (free leading gap on `target`) -- backtrack terminates here.
        Start,
        /// From `dp[i-1][j-1]`: consumes `query[i-1]` AND `target[j-1]` (match or mismatch).
        Diag,
        /// From `dp[i-1][j]`: consumes `query[i-1]` only -- a gap in `target`.
        Up,
        /// From `dp[i][j-1]`: consumes `target[j-1]` only -- a gap in `query`.
        Left,
    }

    /// HW/infix alignment of `query` against `target`, WITH the traceback path (unlike
    /// `seq_utils::hw_distance`, which only returns the distance via a rolling 2-row DP).
    ///
    /// Same semantics as `hw_distance`: `query` is rows, `target` is columns; row 0 is `0` across
    /// every target column (free leading gap on `target`); the alignment ends at the min-cost cell in
    /// the LAST query row (free trailing gap on `target`) and is backtracked to query row 0. Match
    /// cost 0, substitution/indel cost 1.
    ///
    /// Returns the aligned columns in order from the start of the alignment to its end: `(Some(qi),
    /// Some(ti))` for a match/mismatch, `(Some(qi), None)` for a gap in `target`, `(None, Some(ti))`
    /// for a gap in `query`. Every `query` index `0..query.len()` appears exactly once (query has no
    /// free end-gaps under HW); `target` indices outside the aligned span (the free leading/trailing
    /// gap) never appear at all.
    pub fn align_traceback(query: &[u8], target: &[u8]) -> Vec<(Option<usize>, Option<usize>)> {
        let lq = query.len();
        let lt = target.len();
        let cols = lt + 1;

        let mut dp: Vec<usize> = vec![0; (lq + 1) * cols];
        let mut back: Vec<Trace> = vec![Trace::Start; (lq + 1) * cols];

        // Row 0: free leading gap on target -> cost 0 everywhere; Start marks the backtrack terminus.
        // (dp[0][j] is already 0 from the `vec!` initializer; `back` is already `Trace::Start`.)

        for i in 1..=lq {
            dp[i * cols] = i;
            back[i * cols] = Trace::Up; // dp[i][0] = dp[i-1][0] + 1 (gap in target, column stays 0).
            let qi = query[i - 1];
            for j in 1..=lt {
                let sub = dp[(i - 1) * cols + (j - 1)] + if qi == target[j - 1] { 0 } else { 1 };
                let up = dp[(i - 1) * cols + j] + 1;
                let left = dp[i * cols + (j - 1)] + 1;

                let (best, tr) = if sub <= up && sub <= left {
                    (sub, Trace::Diag)
                } else if up <= left {
                    (up, Trace::Up)
                } else {
                    (left, Trace::Left)
                };
                dp[i * cols + j] = best;
                back[i * cols + j] = tr;
            }
        }

        // Free trailing gap on target: the alignment ends at the min-cost cell of the last query row.
        let last_row = lq * cols;
        let mut best_j = 0usize;
        let mut best_cost = dp[last_row];
        for j in 1..=lt {
            let c = dp[last_row + j];
            if c < best_cost {
                best_cost = c;
                best_j = j;
            }
        }

        let mut pairs: Vec<(Option<usize>, Option<usize>)> = Vec::new();
        let mut i = lq;
        let mut j = best_j;
        while i > 0 {
            match back[i * cols + j] {
                Trace::Diag => {
                    pairs.push((Some(i - 1), Some(j - 1)));
                    i -= 1;
                    j -= 1;
                }
                Trace::Up => {
                    pairs.push((Some(i - 1), None));
                    i -= 1;
                }
                Trace::Left => {
                    pairs.push((None, Some(j - 1)));
                    j -= 1;
                }
                Trace::Start => break, // unreachable for i > 0, but avoid looping forever if it were.
            }
        }
        pairs.reverse();
        pairs
    }

    /// For each PSV column `t` (a position in `target`, the copy-path/consensus that `align_map` was
    /// computed against), the read's base observed there: `Some(query[qi])` when `align_map` has an
    /// aligned pair `(Some(qi), Some(t))`, or `None` when the read gaps at that column or the
    /// alignment doesn't span it at all (`t` falls in the free leading/trailing target gap, or lands
    /// on a `(None, Some(t))` query-gap column). Output length matches
    /// `psv_positions_in_consensus.len()`, in the same order.
    pub fn path_obs_at(
        align_map: &[(Option<usize>, Option<usize>)],
        psv_positions_in_consensus: &[usize],
        query: &[u8],
    ) -> Vec<Option<u8>> {
        psv_positions_in_consensus
            .iter()
            .map(|&t| {
                align_map
                    .iter()
                    .find(|&&(_, ti)| ti == Some(t))
                    .and_then(|&(qi, _)| qi.map(|q| query[q]))
            })
            .collect()
    }

    /// C1 fix: orient `read_seq` to match `copy_seq`'s coordinate frame before `align_traceback`.
    ///
    /// Invariant (see `copy_assign_pipeline::fill_psv_obs`/`FamilyProfiles::strand`): a copy's
    /// `DenovoTranscript::seq` (and hence `copy_seqs[k]` here) is in TRANSCRIPTION-strand orientation,
    /// while a `BamRead`'s `read.seq` is FORWARD-GENOME orientation (as SAM/BAM already store it). For
    /// a `+`-strand copy the two frames coincide; for a `-`-strand copy the copy's consensus is the
    /// REVERSE COMPLEMENT of the forward-genome sequence at that locus, so aligning `read.seq` against
    /// it literally (as `align_traceback` does, byte-for-byte, no orientation search of its own --
    /// unlike `aln_id`, which tries both orientations for its identity SCORE) produces a ~0.5-identity
    /// garbage alignment: the wrong or `None` `path_obs` bases that then feed the EM as bogus PSV
    /// evidence.
    ///
    /// Fix: try both `read_seq` and `revcomp_keep_case(read_seq)` against `copy_seq` via `hw_distance` (the same
    /// distance `aln_id` uses internally) and return whichever orientation fits better (ties keep the
    /// forward/as-given orientation). Callers must align (and extract `path_obs_at` from) THIS
    /// returned, oriented sequence -- not the raw `read_seq` -- so the observed bases end up in the
    /// copy's own transcription-strand frame, directly comparable to `copy_psv_alleles`.
    pub fn orient_for_copy(read_seq: &[u8], copy_seq: &[u8]) -> Vec<u8> {
        let d_fwd = hw_distance(read_seq, copy_seq);
        let d_rev = hw_distance(&revcomp_keep_case(read_seq), copy_seq);
        if d_rev < d_fwd {
            revcomp_keep_case(read_seq)
        } else {
            read_seq.to_vec()
        }
    }

    /// Thresholds for flagging a read as a re-align CANDIDATE (poor-fit or unmapped).
    ///
    /// A read is a candidate if it is low-MAPQ (ambiguous/multi-mapped), heavily soft/hard-clipped
    /// (partial alignment, suggesting the reference copy it landed on isn't its true source), or
    /// highly divergent (many mismatches relative to the copy it aligned to). Divergence and clip
    /// fraction are computed by the CALLER from the CIGAR/NM at wiring time (later task); this
    /// struct only carries the thresholds `is_candidate` applies.
    pub struct RealignParams {
        pub max_mapq: u8,
        pub min_clip_frac: f64,
        pub min_div: f64,
        pub min_reads: usize,
    }

    impl Default for RealignParams {
        fn default() -> Self {
            RealignParams {
                max_mapq: 20,
                min_clip_frac: 0.20,
                min_div: 0.05,
                min_reads: 3,
            }
        }
    }

    /// True iff `(mapq, div, clip_frac)` indicates a poor primary-alignment fit under `p`'s
    /// thresholds: low MAPQ (`<= max_mapq`), OR heavy clipping (`>= min_clip_frac`), OR high
    /// divergence (`>= min_div`). Pure function of the three scalars + params -- `div` and
    /// `clip_frac` are computed by the caller from the CIGAR/NM, not derived here.
    pub fn is_candidate(mapq: u8, div: f64, clip_frac: f64, p: &RealignParams) -> bool {
        mapq <= p.max_mapq || clip_frac >= p.min_clip_frac || div >= p.min_div
    }

    /// Minimal identity floor for `realign_to_paths`: a read scoring below this against every
    /// copy-path in the family fits none of them well enough to be a re-align candidate at all.
    /// (Genuine novel copies pooled from low-fit reads are handled separately, by Task 4 --
    /// this floor only gates obvious non-members out of the per-family re-align step.)
    pub const MIN_ALN_ID: f64 = 0.5;

    /// A candidate read's best fit among a family's copy-paths (Task 3), used by Task 4 to decide
    /// whether re-aligning to the family beats the read's existing linear-locus alignment.
    pub struct RealignHit {
        /// Index into `copy_seqs` of the best-fitting copy-path (earliest index on ties).
        pub best_copy: usize,
        /// `aln_id(read_seq, copy_seqs[best_copy])` -- the best copy-path identity.
        pub id_best: f64,
        /// `aln_id(read_seq, copy_seqs[linear_copy])` when `linear_copy` is `Some` and in range,
        /// else `0.0` -- the read's fit to the copy it would otherwise be attributed to linearly.
        pub id_linear: f64,
    }

    /// Re-align `read_seq` to every copy-path in `copy_seqs`, reusing `seq_utils::aln_id` as
    /// the fit score (best infix identity, tries both orientations). Returns the best-fitting copy
    /// plus (for Task 4's accept comparison) the fit to `linear_copy`'s sequence, if given.
    ///
    /// Returns `None` when `copy_seqs` is empty, or when the best fit is below `MIN_ALN_ID`: a read
    /// that fits no copy-path at all is not a re-align candidate for this family.
    /// `aln_id` with the SHORTER sequence as the query — the idiom the retired
    /// `bridge_detector::exon_match_tensor` used ("align the SHORTER as the query").
    ///
    /// ⛔ WHY THIS EXISTS (ledger §6df). `aln_id` is `1 - hw_distance(q, t) / len(q)` with free end-gaps on the
    /// TARGET only, so a target SHORTER than the query cannot contain it and the identity is CAPPED at
    /// `len(t) / len(q)` — a pure length artefact with no homology content. Measured on MCLFAM2: a 721 bp
    /// fragment copy against 2,685 bp median reads is capped at 721/2,685 = **0.269**, so it could never win
    /// and **362/362 = 100%** of the reads whose linear copy it was were reassigned away from it
    /// unconditionally. Any family with copies of unequal length is affected, and fragments are common.
    fn aln_id_len_safe(read: &[u8], copy: &[u8]) -> f64 {
        if read.len() <= copy.len() {
            aln_id(read, copy)
        } else {
            aln_id(copy, read)
        }
    }

    pub fn realign_to_paths(
        read_seq: &[u8],
        copy_seqs: &[Vec<u8>],
        linear_copy: Option<usize>,
    ) -> Option<RealignHit> {
        if copy_seqs.is_empty() {
            return None;
        }

        let mut best_copy = 0usize;
        let mut id_best = f64::NEG_INFINITY;
        for (k, seq) in copy_seqs.iter().enumerate() {
            let id = aln_id_len_safe(read_seq, seq);
            // Strict `>` keeps the earliest index on ties.
            if id > id_best {
                id_best = id;
                best_copy = k;
            }
        }

        if id_best < MIN_ALN_ID {
            return None;
        }

        let id_linear = match linear_copy {
            Some(lc) if lc < copy_seqs.len() => aln_id_len_safe(read_seq, &copy_seqs[lc]),
            _ => 0.0,
        };

        Some(RealignHit {
            best_copy,
            id_best,
            id_linear,
        })
    }

    /// Task 4's verdict on a candidate re-alignment: either correct the read's copy attribution to
    /// `Reassign(best_copy)`, or `Reject` the correction (keep whatever attribution the caller
    /// already had -- linear, or none). Admission of genuinely novel copies (no existing attribution
    /// at all) is a separate concern, handled by `pool_novel` here and `absent_copy::admit_candidate`
    /// at wiring time -- this enum only covers correcting an EXISTING attribution.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum RealignAction {
        Reassign(usize),
        Reject,
    }

    /// Decide whether `hit` (a candidate read's re-alignment result from `realign_to_paths`) beats
    /// its existing linear-locus attribution `linear_copy` significantly enough to correct it.
    ///
    /// Mirrors `copy_assign::read_copy_evidence`'s `min_p` certificate: `n_decisive` is the number of
    /// read positions (out of `read_len`) that support the best copy-path over the linear locus
    /// (`(id_best - id_linear) * read_len`, rounded), and `min_p = (error_rate / 3)^n_decisive` is the
    /// probability that all of those decisive differences arose by sequencing error alone (an
    /// `epsilon^delta` bound: each independent error has probability `error_rate / 3` of landing on
    /// the specific alternate base that agrees with the best copy-path). `min_p < alpha` certifies the
    /// correction; otherwise the evidence isn't strong enough to overturn the existing attribution.
    ///
    /// No correction is needed (and none is offered) when the read's best copy-path already IS its
    /// linear attribution, or when there's no decisive evidence (`n_decisive < 1`) at all.
    ///
    /// When `linear_copy` is `None` this always `Reject`s. The `min_p` certificate here certifies a
    /// *correction* -- that the best copy-path beats an existing linear attribution significantly
    /// enough to overturn it. `realign_to_paths` fills `id_linear` with a `0.0` sentinel when there is
    /// no linear copy to compare against, which is not a real identity and cannot serve as a
    /// baseline: certifying against it would accept any moderately-well-fitting `id_best` (even a
    /// near-random ~0.5 identity) as a "correction" of nothing. A read with no linear attribution at
    /// all isn't a correction case -- it's handled by the separate novel-copy path (`pool_novel`).
    /// Count the PSV columns at which the read's own observed allele SUPPORTS `best` and CONTRADICTS
    /// `linear` — the quantity `min_p = (error_rate/3)^n` actually assumes.
    ///
    /// WHY THIS REPLACES `(id_best - id_linear) * read_len` (ledger §6df). That expression is a whole-read
    /// EDIT-DISTANCE difference, so one indel or one length mismatch manufactures hundreds of "independent
    /// allele observations" and `min_p` underflows to 0, making the reassignment unconditional. Measured on
    /// MCLFAM2 it gave n_decisive ~ 754 for reads whose linear copy was a 721 bp fragment, and 362/362 =
    /// 100% of those reads were reassigned without a single rejection.
    ///
    /// Columns are compared BY INDEX: `psv_best[j]` and `psv_linear[j]` are the same family PSV column j in
    /// each copy's own offset frame (`psv_positions_for` builds both from `family_col_genomic_pos`). A column
    /// counts only when BOTH copies yield a called base and the read's base matches one and not the other.
    fn psv_decisive_count(
        read_seq: &[u8],
        copy_best: &[u8],
        copy_linear: &[u8],
        psv_best: &[usize],
        psv_linear: &[usize],
    ) -> usize {
        // The read's own base at each PSV column, read through its alignment to the BEST copy.
        // (`path_obs_at` yields the READ's base, not the copy's — comparing two such vectors to each
        // other would compare the read to itself and always count 0.)
        let oriented = orient_for_copy(read_seq, copy_best);
        let map = align_traceback(&oriented, copy_best);
        let obs = path_obs_at(&map, psv_best, &oriented);

        (0..obs.len().min(psv_linear.len()))
            .filter(|&j| {
                let Some(r) = obs[j] else { return false };
                let (Some(&cb), Some(&cl)) =
                    (copy_best.get(psv_best[j]), copy_linear.get(psv_linear[j]))
                else {
                    return false;
                };
                // A column is decisive only where the two COPIES actually differ and the read picks best.
                cb != cl && r == cb
            })
            .count()
    }

    pub fn accept_realignment(
        hit: &RealignHit,
        linear_copy: Option<usize>,
        n_decisive: usize,
        error_rate: f64,
        alpha: f64,
    ) -> RealignAction {
        let Some(linear_copy) = linear_copy else {
            return RealignAction::Reject;
        };

        if linear_copy == hit.best_copy {
            return RealignAction::Reject;
        }

        if n_decisive < 1 {
            return RealignAction::Reject;
        }
        let n_decisive = n_decisive as i64;

        let min_p = (error_rate / 3.0).powi(n_decisive as i32);
        if min_p < alpha {
            RealignAction::Reassign(hit.best_copy)
        } else {
            RealignAction::Reject
        }
    }

    /// Greedily single-linkage cluster `unfit` reads (those `realign_to_paths` matched to NO existing
    /// copy -- candidate reference-absent/novel-copy material) by pairwise `aln_id >= min_id`.
    ///
    /// Each read joins the first existing cluster whose FIRST member (the cluster's representative)
    /// it matches at `>= min_id`; if it matches no cluster's representative, it starts a new
    /// singleton cluster. This is a cheap O(n * clusters) pass, not full correlation clustering --
    /// good enough to pool obviously-related novel-copy candidates for the Task-5 wiring, which is
    /// where the actual `absent_copy::admit_candidate` admission gate (needing the genome + remap)
    /// runs. Returns only clusters with `>= min_reads` members, as index vectors into `unfit`.
    pub fn pool_novel(
        unfit: &[(String, Vec<u8>)],
        min_id: f64,
        min_reads: usize,
    ) -> Vec<Vec<usize>> {
        let mut clusters: Vec<Vec<usize>> = Vec::new();

        for (i, (_, seq)) in unfit.iter().enumerate() {
            let mut joined = false;
            for cluster in clusters.iter_mut() {
                let rep = cluster[0];
                if aln_id(seq, &unfit[rep].1) >= min_id {
                    cluster.push(i);
                    joined = true;
                    break;
                }
            }
            if !joined {
                clusters.push(vec![i]);
            }
        }

        clusters.retain(|c| c.len() >= min_reads);
        clusters
    }

    /// One per-read decision from the per-family re-align supplement (Task 5), emitted verbatim as a row of
    /// `<out>.vg_realign.tsv` by the `copy_assign` binary.
    #[derive(Debug, Clone, PartialEq)]
    pub struct RealignRecord {
        pub read_name: String,
        /// `"reassigned"` (a significant correction to a different copy), `"rejected"` (a candidate that
        /// re-aligned but didn't clear Task 4's significance certificate), or `"novel-candidate"` (fits no
        /// existing copy-path at all -- `Task 4`'s `pool_novel`/`absent_copy` admission gate is the separate,
        /// out-of-scope next step for these).
        pub action: String,
        /// The copy index the record concerns: `best_copy` when `action == "reassigned"`, else `-1` (no
        /// correction target -- `"rejected"` keeps the read's existing linear attribution, `"novel-candidate"`
        /// has no copy-path fit at all).
        pub target_copy: i64,
        /// `RealignHit::id_best` (0.0 for `"novel-candidate"`, where `realign_to_paths` found no fit at all so
        /// no `id_best` was computed).
        pub id_best: f64,
        /// The copy index the read's own linear (BAM-coordinate) alignment placed it on within this family, or
        /// `-1` if none of the family's copies on the read's chromosome overlap its aligned span
        /// (`best_overlap_copy_on` returned `None`).
        pub linear_copy: i64,
    }

    /// VG re-align END-TO-END plan, Task 2: `apply_realign`'s output.
    ///
    /// `admitted` reference-absent copies are NOT produced here -- see the field doc on
    /// `novel_pools` and `apply_realign`'s doc for why (this stays genome-free/testable; admission
    /// is the follow-up wiring task's job).
    pub struct RealignApply {
        /// `read_index -> (new_copy_idx, path_obs)`: a significant correction (Task 4's
        /// `RealignAction::Reassign`) to a DIFFERENT copy than the read's existing linear
        /// attribution, plus the read's base at each of the new copy's PSV columns (from
        /// `align_traceback` + `path_obs_at` against the new copy's consensus).
        pub corrected: crate::types::DetHashMap<usize, (usize, Vec<Option<u8>>)>,
        /// Clusters (each `>= rp.min_reads` members) of read INDICES (into `bam_reads`) that fit no
        /// existing copy-path at all (`realign_to_paths` returned `None`) but are mutually similar
        /// enough (`pool_novel`, `min_id ~= 0.9`) to be candidate novel/reference-absent copies. Not
        /// yet admitted -- the wiring task turns each pool into a `CollapsedCandidate` and runs it
        /// through `absent_copy::admit_candidate` with the real genome + remap.
        pub novel_pools: Vec<Vec<usize>>,
        /// One record per candidate read processed (`"reassigned"`, `"rejected"`, or
        /// `"novel-candidate"`).
        pub records: Vec<RealignRecord>,
    }

    /// VG re-align END-TO-END plan, Task 2: apply the per-read decisions (Tasks 1/3/4) into a
    /// ready-to-consume correction map + novel-copy candidate pools, over reads spanning potentially
    /// SEVERAL families at once.
    ///
    /// `copies`/`copy_seqs` are parallel (one spliced consensus per copy, `copy_seqs[k] ==
    /// copies[k].seq` is the expected caller invariant but only `copy_seqs` is actually read here --
    /// `copies` is carried for callers/future use, e.g. locus metadata alongside the correction map).
    /// `psv_pos_per_copy[k]` are the family's PSV positions in `copy_seqs[k]`'s consensus coordinates.
    /// `linear_copy_of[i]` is read `i`'s existing linear-locus copy attribution (or `None`), parallel
    /// to `bam_reads`.
    ///
    /// CORRECTIONS (must-have): a candidate read (`is_candidate`) that re-aligns to a copy-path
    /// (`realign_to_paths`) and clears `accept_realignment`'s significance certificate
    /// (`RealignAction::Reassign`) is entered into `corrected[read_index] = (new_copy, path_obs)`,
    /// where `path_obs` is the read's base at each of the new copy's PSV columns
    /// (`align_traceback` + `path_obs_at` against `copy_seqs[new_copy]`).
    ///
    /// ADMISSIONS (best-effort/mechanism-only): a candidate read that fits NO copy-path at all
    /// (`realign_to_paths` returns `None`) is pooled with other such "unfit" reads by
    /// `pool_novel`; `novel_pools` holds the resulting clusters (mapped back to indices into
    /// `bam_reads`) with `>= rp.min_reads` members. This is genome-free and does NOT run
    /// `absent_copy::admit_candidate` -- turning a pool into an admitted reference-absent copy needs
    /// the real genome + remap, which is the follow-up wiring task's job, not this one's. Real yield
    /// here is data-limited (the O4 divergent frontier): most families will produce zero pools.
    ///
    /// Reads failing `is_candidate` are skipped entirely (no record, no correction, no pooling) --
    /// a clean primary fit has nothing to reconsider. Supplementary alignments (`is_supplementary`)
    /// are also skipped.
    pub fn apply_realign(
        bam_reads: &[BamRead],
        copies: &[DenovoTranscript],
        copy_seqs: &[Vec<u8>],
        psv_pos_per_copy: &[Vec<usize>],
        linear_copy_of: &[Option<usize>],
        rp: &RealignParams,
        error_rate: f64,
        alpha: f64,
    ) -> RealignApply {
        let _ = copies; // parallel to copy_seqs; not read directly here (see doc).

        let mut corrected = crate::types::DetHashMap::default();
        let mut records = Vec::new();
        let mut unfit: Vec<(String, Vec<u8>)> = Vec::new();
        let mut unfit_idx: Vec<usize> = Vec::new();

        // Per-candidate realignment is independent (reads shared &copy_seqs/&psv_pos_per_copy, no shared mut),
        // so parallelize it across reads and merge the outcomes back in READ-INDEX order — byte-identical to the
        // former serial loop (`corrected` is a HashMap so order-free; `records`/`unfit`/`unfit_idx` are rebuilt in
        // ascending read index, exactly as the serial loop pushed them). This is the correction leg's cost centre.
        use rayon::prelude::*;
        enum Outcome {
            Skip,
            Corrected {
                i: usize,
                copy: usize,
                obs: Vec<Option<u8>>,
                rec: RealignRecord,
            },
            Novel {
                i: usize,
                name: String,
                seq: Vec<u8>,
                rec: RealignRecord,
            },
            Rejected {
                rec: RealignRecord,
            },
        }
        let outcomes: Vec<Outcome> = bam_reads
            .par_iter()
            .enumerate()
            .map(|(i, br)| {
                if br.is_supplementary {
                    return Outcome::Skip;
                }
                let read_len = br.read.seq.len();
                let clip: u64 = br
                    .read
                    .cigar
                    .iter()
                    .filter(|&&(op, _)| op == 'S')
                    .map(|&(_, n)| n)
                    .sum();
                let clip_frac = if read_len > 0 {
                    clip as f64 / read_len as f64
                } else {
                    0.0
                };
                let div = br.de as f64;
                if !is_candidate(br.mapq, div, clip_frac, rp) {
                    return Outcome::Skip;
                }
                let linear_copy = linear_copy_of.get(i).copied().flatten();
                let linear_copy_i64 = linear_copy.map(|c| c as i64).unwrap_or(-1);
                match realign_to_paths(&br.read.seq, copy_seqs, linear_copy) {
                    None => Outcome::Novel {
                        i,
                        name: br.name.clone(),
                        seq: br.read.seq.clone(),
                        rec: RealignRecord {
                            read_name: br.name.clone(),
                            action: "novel-candidate".to_string(),
                            target_copy: -1,
                            id_best: 0.0,
                            linear_copy: linear_copy_i64,
                        },
                    },
                    Some(hit) => {
                        let id_best = hit.id_best;
                        // §6df: the decisive count is a PSV-COLUMN count, never a scaled edit distance.
                        let n_dec = match linear_copy {
                            Some(lc) if lc != hit.best_copy => psv_decisive_count(
                                &br.read.seq,
                                &copy_seqs[hit.best_copy],
                                &copy_seqs[lc],
                                &psv_pos_per_copy[hit.best_copy],
                                &psv_pos_per_copy[lc],
                            ),
                            _ => 0,
                        };
                        match accept_realignment(&hit, linear_copy, n_dec, error_rate, alpha) {
                            RealignAction::Reassign(best_copy) => {
                                // C1: orient the read to the copy's own (transcription-strand) frame before the
                                // traceback -- a `-`-strand copy's consensus is the reverse complement of
                                // forward-genome, so a literal forward alignment here would produce garbage `obs`.
                                let oriented = orient_for_copy(&br.read.seq, &copy_seqs[best_copy]);
                                let map = align_traceback(&oriented, &copy_seqs[best_copy]);
                                let obs =
                                    path_obs_at(&map, &psv_pos_per_copy[best_copy], &oriented);
                                Outcome::Corrected {
                                    i,
                                    copy: best_copy,
                                    obs,
                                    rec: RealignRecord {
                                        read_name: br.name.clone(),
                                        action: "reassigned".to_string(),
                                        target_copy: best_copy as i64,
                                        id_best,
                                        linear_copy: linear_copy_i64,
                                    },
                                }
                            }
                            RealignAction::Reject => Outcome::Rejected {
                                rec: RealignRecord {
                                    read_name: br.name.clone(),
                                    action: "rejected".to_string(),
                                    target_copy: -1,
                                    id_best,
                                    linear_copy: linear_copy_i64,
                                },
                            },
                        }
                    }
                }
            })
            .collect();
        for o in outcomes {
            match o {
                Outcome::Skip => {}
                Outcome::Corrected { i, copy, obs, rec } => {
                    corrected.insert(i, (copy, obs));
                    records.push(rec);
                }
                Outcome::Novel { i, name, seq, rec } => {
                    unfit.push((name, seq));
                    unfit_idx.push(i);
                    records.push(rec);
                }
                Outcome::Rejected { rec } => records.push(rec),
            }
        }

        let clusters = pool_novel(&unfit, 0.9, rp.min_reads);
        let novel_pools: Vec<Vec<usize>> = clusters
            .into_iter()
            .map(|c| c.into_iter().map(|j| unfit_idx[j]).collect())
            .collect();

        RealignApply {
            corrected,
            novel_pools,
            records,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::family::copy_split::AlignedRead;
        use crate::family::seq_utils::hw_distance;

        #[test]
        fn is_candidate_flags_poor_fit() {
            let p = RealignParams::default();

            // low MAPQ alone -> true
            assert!(is_candidate(5, 0.0, 0.0, &p));
            // high clip alone -> true
            assert!(is_candidate(60, 0.0, 0.30, &p));
            // high divergence alone -> true
            assert!(is_candidate(60, 0.08, 0.0, &p));
            // clean read: high MAPQ, no clip, low div -> false
            assert!(!is_candidate(60, 0.0, 0.0, &p));

            // boundary: mapq == max_mapq is still <= -> true
            assert!(is_candidate(20, 0.0, 0.0, &p));
            // boundary: just under both clip and div thresholds, and mapq above max -> false
            assert!(!is_candidate(21, 0.049, 0.19, &p));
        }

        /// Deterministic pseudo-random ACGT sequence generator (xorshift64), so test sequences are
        /// reproducible without hand-typing long strings or depending on real fixture data.
        fn pseudo_seq(seed: u64, len: usize) -> Vec<u8> {
            let bases = [b'A', b'C', b'G', b'T'];
            // splitmix64-style mix so nearby/small seeds (1, 2, 3, ...) don't collapse to the same
            // or correlated xorshift states (plain `seed | 1` made seeds 2 and 3 identical).
            let mut state = seed
                .wrapping_mul(0x9E3779B97F4A7C15)
                .wrapping_add(0x2545F4914F6CDD1D)
                | 1;
            (0..len)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    bases[(state % 4) as usize]
                })
                .collect()
        }

        #[test]
        fn realign_picks_best_copy_path() {
            let copy0 = pseudo_seq(1, 60);
            let copy1 = pseudo_seq(2, 60);
            let copy2 = pseudo_seq(3, 60);
            let copy_seqs = vec![copy0.clone(), copy1.clone(), copy2.clone()];

            // Read == copy 1 exactly -> best_copy == 1, id_best ~1.0.
            let hit = realign_to_paths(&copy1, &copy_seqs, Some(0))
                .expect("exact copy-path match must be a candidate");
            assert_eq!(hit.best_copy, 1);
            assert!(hit.id_best > 0.99, "id_best = {}", hit.id_best);
            // id_linear (fit to copy 0, the "linear locus" copy) must be strictly lower than the
            // true best-copy fit -- the read really belongs to copy 1, not copy 0.
            assert!(
                hit.id_linear < hit.id_best,
                "id_linear = {} should be < id_best = {}",
                hit.id_linear,
                hit.id_best
            );

            // Read == copy 2 with a couple of substitutions -> still best_copy == 2, id_best < 1.0.
            let mut mutated = copy2.clone();
            // Flip two bases to a base guaranteed different from the original.
            for &pos in &[10usize, 40usize] {
                let orig = mutated[pos];
                mutated[pos] = [b'A', b'C', b'G', b'T']
                    .into_iter()
                    .find(|&b| b != orig)
                    .unwrap();
            }
            let hit2 = realign_to_paths(&mutated, &copy_seqs, None)
                .expect("near-exact copy-path match must be a candidate");
            assert_eq!(hit2.best_copy, 2);
            assert!(hit2.id_best < 1.0, "id_best = {}", hit2.id_best);
            assert!(
                hit2.id_best > 0.9,
                "id_best = {} should still be a strong fit",
                hit2.id_best
            );
            // No linear_copy given -> id_linear is the documented 0.0 sentinel.
            assert_eq!(hit2.id_linear, 0.0);
        }

        #[test]
        fn realign_returns_none_for_nonmember() {
            let copy_seqs = vec![pseudo_seq(1, 60), pseudo_seq(2, 60), pseudo_seq(3, 60)];
            // A read unrelated to any copy-path (seed chosen so its best infix identity to all
            // three copies lands below MIN_ALN_ID; random same-length sequences under free-end-gap
            // edit distance land around ~0.5 identity by chance, so this isn't a free lunch).
            let read = pseudo_seq(3130, 60);
            let hit = realign_to_paths(&read, &copy_seqs, None);
            match hit {
                None => {}
                Some(h) => panic!(
                    "expected None for a non-member read, got id_best = {}",
                    h.id_best
                ),
            }

            // Empty copy_seqs -> always None regardless of the read.
            assert!(realign_to_paths(&read, &[], None).is_none());
        }

        #[test]
        fn realign_params_defaults() {
            let p = RealignParams::default();
            assert_eq!(p.max_mapq, 20);
            assert_eq!(p.min_clip_frac, 0.20);
            assert_eq!(p.min_div, 0.05);
            assert_eq!(p.min_reads, 3);
        }

        #[test]
        fn accept_significant_reassigns() {
            // id diff 0.14 * read_len 1000 -> n_decisive = 140 decisive positions favoring copy 2
            // over the read's current linear attribution (copy 0). min_p = (0.003/3)^140 is
            // astronomically small (<< alpha = 1e-3) -- certifies the correction.
            let hit = RealignHit {
                best_copy: 2,
                id_best: 0.99,
                id_linear: 0.85,
            };
            let n_decisive = ((hit.id_best - hit.id_linear) * 1000.0).round() as i64;
            assert_eq!(n_decisive, 140);

            let action = accept_realignment(&hit, Some(0), 140, 0.003, 1e-3);
            assert_eq!(action, RealignAction::Reassign(2));
        }

        #[test]
        fn accept_marginal_rejects() {
            // ONE decisive PSV column. With a deliberately high error_rate = 0.05,
            // min_p = (0.05/3)^1 = 0.01666... >= alpha = 1e-3 -- a single decisive column isn't enough
            // to overturn the read's existing linear attribution (copy 0) under this noisy an error model.
            // §6df: `n_decisive` is now a PSV-COLUMN COUNT supplied by the caller, not a scaled edit
            // distance; passing 1 here is what one distinguishing column means.
            let hit = RealignHit {
                best_copy: 2,
                id_best: 0.90,
                id_linear: 0.899,
            };
            let min_p = (0.05_f64 / 3.0).powi(1);
            assert!(min_p >= 1e-3, "min_p = {min_p} should be >= alpha");

            let action = accept_realignment(&hit, Some(0), 1, 0.05, 1e-3);
            assert_eq!(action, RealignAction::Reject);
        }

        #[test]
        fn accept_zero_decisive_rejects() {
            // No PSV column separates the two copies for this read -> n_decisive = 0 -> no decisive
            // evidence at all -> Reject, regardless of how permissive alpha/error_rate are.
            let hit = RealignHit {
                best_copy: 2,
                id_best: 0.95,
                id_linear: 0.95,
            };
            let action = accept_realignment(&hit, Some(0), 0, 0.003, 1.0);
            assert_eq!(action, RealignAction::Reject);
        }

        #[test]
        fn accept_best_equals_linear_rejects() {
            // best_copy already IS the read's linear attribution -- no correction needed even
            // though id_best/id_linear here would otherwise look decisive.
            let hit = RealignHit {
                best_copy: 2,
                id_best: 0.99,
                id_linear: 0.10,
            };
            let action = accept_realignment(&hit, Some(2), 890, 0.003, 1e-3);
            assert_eq!(action, RealignAction::Reject);
        }

        #[test]
        fn accept_none_linear_rejects() {
            // No existing linear attribution at all (unmapped read routed by Task 3) -- id_linear's
            // 0.0 sentinel is not a real baseline, so there is nothing to "correct" against. Even a
            // high id_best (0.99) must Reject here, not Reassign off the meaningless zero baseline;
            // genuinely unattributed reads are handled by the separate novel-copy path.
            let hit = RealignHit {
                best_copy: 1,
                id_best: 0.99,
                id_linear: 0.0,
            };
            let action = accept_realignment(&hit, None, 140, 0.003, 1e-3);
            assert_eq!(action, RealignAction::Reject);
        }

        #[test]
        fn pool_novel_clusters_unfit() {
            // 3 mutually near-identical reads (one exact copy + two 1-base mutants of it) plus 1
            // unrelated random read. min_id = 0.9, min_reads = 3.
            let base = pseudo_seq(100, 80);
            let mut mut1 = base.clone();
            mut1[5] = [b'A', b'C', b'G', b'T']
                .into_iter()
                .find(|&b| b != mut1[5])
                .unwrap();
            let mut mut2 = base.clone();
            mut2[60] = [b'A', b'C', b'G', b'T']
                .into_iter()
                .find(|&b| b != mut2[60])
                .unwrap();
            // Unrelated: different seed, same length, chosen so its identity to `base` lands well
            // under 0.9 (random same-length sequences under free-end-gap identity hover ~0.5).
            let unrelated = pseudo_seq(9999, 80);
            assert!(
                aln_id(&base, &unrelated) < 0.9,
                "fixture assumption broken: unrelated read too similar to base"
            );

            let unfit: Vec<(String, Vec<u8>)> = vec![
                ("r0".to_string(), base),
                ("r1".to_string(), mut1),
                ("r2".to_string(), mut2),
                ("r3".to_string(), unrelated),
            ];

            let clusters = pool_novel(&unfit, 0.9, 3);
            assert_eq!(
                clusters.len(),
                1,
                "expected exactly one cluster to survive min_reads, got {clusters:?}"
            );
            let mut got = clusters[0].clone();
            got.sort_unstable();
            assert_eq!(
                got,
                vec![0, 1, 2],
                "expected the 3 mutually-similar reads clustered together"
            );
        }

        #[test]
        fn pool_novel_below_min_reads_dropped() {
            let a = pseudo_seq(1, 60);
            let b = a.clone();
            let unfit: Vec<(String, Vec<u8>)> = vec![("a".to_string(), a), ("b".to_string(), b)];
            // Only 2 mutually-identical reads, but min_reads = 3 -> nothing survives.
            let clusters = pool_novel(&unfit, 0.9, 3);
            assert!(
                clusters.is_empty(),
                "expected no clusters to meet min_reads = 3, got {clusters:?}"
            );
        }

        /// Minimal `DenovoTranscript` builder for the `apply_realign` tests -- unspliced (no introns),
        /// mirroring how `copy_assign_pipeline`'s own tests construct copies.
        fn transcript(tid: &str, chrom: &str, start: u64, seq: Vec<u8>) -> DenovoTranscript {
            let end = start + seq.len() as u64;
            DenovoTranscript {
                tid: tid.to_string(),
                chrom: chrom.to_string(),
                start,
                end,
                n_reads: 10,
                strand: '+',
                introns: vec![],
                seq,
                ..Default::default()
            }
        }

        fn bam_read(
            name: &str,
            chrom: &str,
            ref_start: u64,
            seq: Vec<u8>,
            mapq: u8,
            de: f32,
        ) -> BamRead {
            let len = seq.len() as u64;
            BamRead {
                chrom: chrom.to_string(),
                read: AlignedRead {
                    ref_start,
                    cigar: vec![('M', len)],
                    seq,
                    qual: vec![],
                },
                mapq,
                name: name.to_string(),
                as_score: 0,
                de,
                is_supplementary: false,
                is_secondary: false,
                reverse: false,
                ts: None,
            }
        }

        // -----------------------------------------------------------------------------------------
        // Task 1 (end-to-end plan): align_traceback + path_obs_at
        // -----------------------------------------------------------------------------------------

        /// Sum of mismatched-diagonal columns + gap columns (either side) in an alignment map --
        /// this must equal `hw_distance`'s edit distance for any valid traceback of that DP.
        fn edits_in_map(
            align_map: &[(Option<usize>, Option<usize>)],
            query: &[u8],
            target: &[u8],
        ) -> usize {
            align_map
                .iter()
                .filter(|&&(qi, ti)| match (qi, ti) {
                    (Some(q), Some(t)) => query[q] != target[t],
                    (Some(_), None) | (None, Some(_)) => true,
                    (None, None) => panic!("align_traceback must never emit a (None, None) column"),
                })
                .count()
        }

        #[test]
        fn traceback_edit_distance_matches_hw() {
            let cases: Vec<(&[u8], &[u8])> = vec![
                // exact infix: ACGT occurs verbatim inside TTACGTGG -> dist 0.
                (b"ACGT", b"TTACGTGG"),
                // 1 substitution: ACGT vs an infix that differs at one base (ACCT inside TT-ACCT-GG).
                (b"ACGT", b"TTACCTGG"),
                // 1 insertion (relative to target): query has an extra base not in any target infix.
                (b"ACGGT", b"TTACGTGG"),
                // 1 deletion (relative to target): query is missing a base present in the target infix.
                (b"ACT", b"TTACGTGG"),
                // longer query with 2 edits scattered through it.
                (b"ACGTACGTAC", b"TTACGTACCTAGGG"),
            ];

            for (query, target) in cases {
                let expected = hw_distance(query, target);
                let map = align_traceback(query, target);
                let got = edits_in_map(&map, query, target);
                assert_eq!(
                    got,
                    expected,
                    "query={:?} target={:?}: traceback edits {got} != hw_distance {expected}",
                    std::str::from_utf8(query).unwrap(),
                    std::str::from_utf8(target).unwrap()
                );

                // Sanity: the aligned columns must walk query positions 0..query.len() in order (every
                // query base consumed exactly once, monotonically), since HW gives free end-gaps only on
                // the TARGET, not the query.
                let q_positions: Vec<usize> = map.iter().filter_map(|&(qi, _)| qi).collect();
                let expected_q_positions: Vec<usize> = (0..query.len()).collect();
                assert_eq!(
                    q_positions, expected_q_positions,
                    "every query position must appear exactly once, in order"
                );
            }
        }

        // -----------------------------------------------------------------------------------------
        // VG re-align END-TO-END plan, Task 2: apply_realign orchestrator
        // -----------------------------------------------------------------------------------------

        #[test]
        fn apply_corrects_mismapped_read() {
            // Two distinct-consensus copies; copy 1 has a distinguishing base at consensus position
            // `p` relative to copy 0 (guaranteed distinct at that column by construction).
            let copy0_seq = pseudo_seq(1, 200);
            let mut copy1_seq = pseudo_seq(2, 200);
            let p = 50usize;
            let p2 = 120usize; // §6df: one PSV gives min_p = (0.003/3)^1 = 1e-3, which does NOT clear
                               // alpha = 1e-3. Two decisive columns is the same bar `absent_copy` uses.
                               // Force column p to differ between the two copies (pseudo_seq draws from independent
                               // seeds so this already usually holds, but force it so the test isn't seed-lucky).
            if copy1_seq[p] == copy0_seq[p] {
                copy1_seq[p] = [b'A', b'C', b'G', b'T']
                    .into_iter()
                    .find(|&b| b != copy0_seq[p])
                    .unwrap();
            }
            if copy1_seq[p2] == copy0_seq[p2] {
                copy1_seq[p2] = [b'A', b'C', b'G', b'T']
                    .into_iter()
                    .find(|&b| b != copy0_seq[p2])
                    .unwrap();
            }
            let copy0 = transcript("copy0", "chr1", 0, copy0_seq.clone());
            let copy1 = transcript("copy1", "chr1", 5000, copy1_seq.clone());
            let copies = vec![copy0, copy1];
            let copy_seqs = vec![copy0_seq, copy1_seq.clone()];
            let psv_pos_per_copy = vec![vec![p, p2], vec![p, p2]];

            // Read == copy 1's consensus exactly, but low MAPQ and linearly attributed to copy 0
            // (its BAM placement) -- a Task-1 candidate whose true source is copy 1.
            let br = bam_read("readA", "chr1", 0, copy1_seq.clone(), 3, 0.0);
            let linear_copy_of = vec![Some(0usize)];

            let out = apply_realign(
                &[br],
                &copies,
                &copy_seqs,
                &psv_pos_per_copy,
                &linear_copy_of,
                &RealignParams::default(),
                0.003,
                1e-3,
            );

            assert_eq!(
                out.corrected.len(),
                1,
                "expected exactly one correction, got {:?}",
                out.corrected
            );
            let (new_copy, obs) = out.corrected.get(&0).expect("read 0 must be corrected");
            assert_eq!(*new_copy, 1, "read's true best-fit copy is copy 1");
            assert_eq!(
                obs,
                &vec![Some(copy1_seq[p]), Some(copy1_seq[p2])],
                "obs at the PSV column must equal copy 1's base"
            );
            assert!(out.novel_pools.is_empty());
            assert_eq!(
                out.records
                    .iter()
                    .filter(|r| r.action == "reassigned")
                    .count(),
                1,
                "expected one reassigned record, got {:?}",
                out.records
            );
        }

        /// C1: a `-`-strand copy's `copy_seqs[k]` is TRANSCRIPTION strand, i.e. the reverse complement
        /// of the forward-genome sequence at that locus. A read is always FORWARD-GENOME (`read.seq`,
        /// per BAM convention), so a read that is truly this copy's source material arrives as
        /// `revcomp_keep_case(copy_seq)`, not `copy_seq` itself. `realign_to_paths`/`aln_id` already handle this
        /// (they try both orientations for the identity SCORE), so the correction is still detected --
        /// but before the C1 fix, `apply_realign`'s `align_traceback`/`path_obs_at` aligned the raw
        /// forward `read.seq` literally against `copy_seq`, producing a ~0.5-identity garbage
        /// alignment and a wrong/`None` `path_obs` at the PSV column. This test pins that the
        /// corrected `path_obs` instead equals copy 1's own TRANSCRIPTION-strand allele -- exactly what
        /// `copy_assign_pipeline::fill_psv_obs`'s per-base `rc_base` for `-` copies would produce.
        #[test]
        fn apply_realign_strand_orients_path_obs_for_minus_copy() {
            let copy0_seq = pseudo_seq(1, 200);
            let mut copy1_seq = pseudo_seq(2, 200); // TRANSCRIPTION-strand consensus of a '-'-strand copy
            let p = 50usize;
            let p2 = 120usize; // §6df: one PSV gives min_p = (0.003/3)^1 = 1e-3, which does NOT clear
                               // alpha = 1e-3. Two decisive columns is the same bar `absent_copy` uses.
            if copy1_seq[p] == copy0_seq[p] {
                copy1_seq[p] = [b'A', b'C', b'G', b'T']
                    .into_iter()
                    .find(|&b| b != copy0_seq[p])
                    .unwrap();
            }
            if copy1_seq[p2] == copy0_seq[p2] {
                copy1_seq[p2] = [b'A', b'C', b'G', b'T']
                    .into_iter()
                    .find(|&b| b != copy0_seq[p2])
                    .unwrap();
            }
            let copy0 = transcript("copy0", "chr1", 0, copy0_seq.clone());
            let mut copy1 = transcript("copy1", "chr1", 5000, copy1_seq.clone());
            copy1.strand = '-';
            let copies = vec![copy0, copy1];
            let copy_seqs = vec![copy0_seq, copy1_seq.clone()];
            let psv_pos_per_copy = vec![vec![p, p2], vec![p, p2]];

            // The read is FORWARD-GENOME: the reverse complement of copy 1's transcription-strand
            // consensus. Low MAPQ (a Task-1 candidate) and linearly misattributed to copy 0.
            let read_seq = revcomp_keep_case(&copy1_seq);
            let br = bam_read("readA", "chr1", 0, read_seq, 3, 0.0);
            let linear_copy_of = vec![Some(0usize)];

            let out = apply_realign(
                &[br],
                &copies,
                &copy_seqs,
                &psv_pos_per_copy,
                &linear_copy_of,
                &RealignParams::default(),
                0.003,
                1e-3,
            );

            assert_eq!(
                out.corrected.len(),
                1,
                "expected exactly one correction, got {:?}",
                out.corrected
            );
            let (new_copy, obs) = out.corrected.get(&0).expect("read 0 must be corrected");
            assert_eq!(
                *new_copy, 1,
                "read's true best-fit copy is copy 1 (the '-'-strand copy)"
            );
            assert_eq!(
                obs,
                &vec![Some(copy1_seq[p]), Some(copy1_seq[p2])],
                "path_obs at the PSV column must equal copy 1's TRANSCRIPTION-strand allele, not a \
                 garbage/wrong-orientation base"
            );
        }

        #[test]
        fn apply_pools_novel_unfit() {
            // 3 mutually near-identical reads that match NEITHER copy (id_best < MIN_ALN_ID so
            // realign_to_paths returns None for all of them) -- candidate novel-copy material.
            let copy0_seq = pseudo_seq(1, 80);
            let copy1_seq = pseudo_seq(2, 80);
            let copy0 = transcript("copy0", "chr1", 0, copy0_seq.clone());
            let copy1 = transcript("copy1", "chr1", 5000, copy1_seq.clone());
            let copies = vec![copy0, copy1];
            let copy_seqs = vec![copy0_seq.clone(), copy1_seq.clone()];
            let psv_pos_per_copy = vec![vec![], vec![]];

            let novel_base = pseudo_seq(69, 80);
            assert!(
                aln_id(&novel_base, &copy0_seq) < MIN_ALN_ID
                    && aln_id(&novel_base, &copy1_seq) < MIN_ALN_ID,
                "fixture assumption broken: novel_base must fit neither existing copy"
            );
            // Exact clones of `novel_base` (not single-base mutants): with a random ~80bp sequence
            // sitting right at the ~0.5 "no better than chance" identity floor against the existing
            // copies, even a 1-base mutation can nudge `aln_id` across the `MIN_ALN_ID` boundary in
            // either direction (edit-distance realignment isn't strictly monotone in Hamming
            // distance). Cloning keeps every read's fit to the existing copies IDENTICAL to the
            // already-asserted `novel_base` fit, while still being "mutually near-identical"
            // (identity 1.0) for `pool_novel`'s `>= min_id` clustering.
            let novel1 = novel_base.clone();
            let novel2 = novel_base.clone();

            // A clean read on its correct copy (high MAPQ, no clip, low div) -- not a candidate at
            // all, must produce no record and not enter the pool.
            let clean = bam_read("clean", "chr1", 0, copy0_seq, 60, 0.0);
            let n0 = bam_read("n0", "chr1", 0, novel_base, 3, 0.0);
            let n1 = bam_read("n1", "chr1", 0, novel1, 3, 0.0);
            let n2 = bam_read("n2", "chr1", 0, novel2, 3, 0.0);

            let bam_reads = vec![clean, n0, n1, n2];
            let linear_copy_of = vec![Some(0usize), None, None, None];

            let out = apply_realign(
                &bam_reads,
                &copies,
                &copy_seqs,
                &psv_pos_per_copy,
                &linear_copy_of,
                &RealignParams::default(),
                0.003,
                1e-3,
            );

            assert!(
                out.corrected.is_empty(),
                "no corrections expected, got {:?}",
                out.corrected
            );
            assert_eq!(
                out.novel_pools.len(),
                1,
                "expected exactly one novel pool, got {:?}",
                out.novel_pools
            );
            let mut got = out.novel_pools[0].clone();
            got.sort_unstable();
            assert_eq!(
                got,
                vec![1, 2, 3],
                "expected read indices 1,2,3 (the 3 novel reads) pooled"
            );

            let novel_records: Vec<_> = out
                .records
                .iter()
                .filter(|r| r.action == "novel-candidate")
                .collect();
            assert_eq!(
                novel_records.len(),
                3,
                "expected one novel-candidate record per unfit read"
            );
            assert!(
                out.records.iter().all(|r| r.read_name != "clean"),
                "the clean read must produce no record at all"
            );
        }

        #[test]
        fn apply_clean_family_noop() {
            let copy0_seq = pseudo_seq(1, 200);
            let copy1_seq = pseudo_seq(2, 200);
            let copy0 = transcript("copy0", "chr1", 0, copy0_seq.clone());
            let copy1 = transcript("copy1", "chr1", 5000, copy1_seq.clone());
            let copies = vec![copy0, copy1];
            let copy_seqs = vec![copy0_seq.clone(), copy1_seq.clone()];
            let psv_pos_per_copy = vec![vec![], vec![]];

            let r0 = bam_read("r0", "chr1", 0, copy0_seq, 60, 0.0);
            let r1 = bam_read("r1", "chr1", 5000, copy1_seq, 60, 0.0);
            let bam_reads = vec![r0, r1];
            let linear_copy_of = vec![Some(0usize), Some(1usize)];

            let out = apply_realign(
                &bam_reads,
                &copies,
                &copy_seqs,
                &psv_pos_per_copy,
                &linear_copy_of,
                &RealignParams::default(),
                0.003,
                1e-3,
            );

            assert!(
                out.corrected.is_empty(),
                "expected no corrections, got {:?}",
                out.corrected
            );
            assert!(
                out.novel_pools.is_empty(),
                "expected no novel pools, got {:?}",
                out.novel_pools
            );
            assert!(
                out.records.is_empty(),
                "expected no records at all, got {:?}",
                out.records
            );
        }

        #[test]
        fn path_obs_reads_bases_at_psv_columns() {
            // T-idx 5 = 'C', T-idx 10 = 'G' (0-based).
            let target: Vec<u8> = b"AAAAACAAAAGAAAAA".to_vec();
            assert_eq!(target[5], b'C');
            assert_eq!(target[10], b'G');

            // Exact full-length match -> both PSV columns are spanned and read verbatim.
            let query = target.clone();
            let map = align_traceback(&query, &target);
            let obs = path_obs_at(&map, &[5, 10], &query);
            assert_eq!(obs, vec![Some(b'C'), Some(b'G')]);

            // A short read that only covers target[..8] -- doesn't reach column 12 at all.
            let short_query: Vec<u8> = target[..8].to_vec();
            let map2 = align_traceback(&short_query, &target);
            let obs2 = path_obs_at(&map2, &[12], &short_query);
            assert_eq!(obs2, vec![None]);
        }
    }
}

pub mod shared_definition {
    //! The shared multi-copy family definition (`docs/seeded_family_definition.md` §0★★), instantiated on the
    //! de novo homology catalog. OPT-IN: `RUSTLE_SHARED_DEFINITION=1`; unset leaves the catalog byte-identical.
    //!
    //! Port of the pre-registered Python prototype `bench/denovo_shared_def.py` (ledger §6jz-§6kd, prereg Addenda
    //! AB-AF). Construction, in order:
    //! 1. NODES — the catalog's read-supported reps, consolidated to gene level: exon chains cut at introns longer
    //!    than [`MAX_INTRON`] (P99.9 of annotated GGO introns), pieces under [`MIN_PIECE`] bp dropped, same-strand
    //!    pieces with overlapping exons merged (Addendum AB).
    //! 2. READ-LOCUS NODES — loci the pipeline dropped: same-strand primary MAPQ >= 1 reads grouped by exon
    //!    overlap, exons = bases at read depth >= 2, added when they overlap no existing node (Addendum AC). With
    //!    `RUSTLE_SD_READ_LOCUS_SPLIT=1`, a chained read group is first split into sub-loci linked by >= 2 reads
    //!    (Addendum AF-3).
    //! 3. EDGES — the guided finders on genomic DNA: a node's spliced transcript (`minimap2 -x splice -uf`) hitting
    //!    another node's exons at identity >= 0.80 over >= 0.50 of the transcript, or its gene body (`-x asm20`)
    //!    chaining onto another node's exons at identity >= 0.80 over >= 0.50 of the shorter body.
    //! 4. FAMILIES — triangle-supported leader neighbourhoods (confirmed on a fresh gorilla substrate, §6kd).
    //!
    //! Every tie-break mirrors the prototype so the two produce the same families on the same inputs (AF-1 gate).
    //!
    //! **STATUS:** OPT-IN  (docs/MODULE_STATUS.md; reached only when `RUSTLE_SHARED_DEFINITION` is set)

    use std::collections::{BTreeMap, BTreeSet, HashMap};
    use std::io::Write;

    use anyhow::{Context, Result};

    use crate::family::denovo_assemble::BamRead;
    use crate::family::family_detect::DenovoTranscript;
    use crate::genome::GenomeIndex;

    /// Introns longer than this cut a node's exon chain (P99.9 of 1,092,233 annotated GGO_genomic.gff introns).
    pub const MAX_INTRON: u64 = 271_359;
    /// Pieces and read loci with fewer exonic bases than this are not nodes.
    pub const MIN_PIECE: u64 = 100;
    const MIN_ID: f64 = 0.80;
    const MIN_COV: f64 = 0.50;
    const MIN_LOCUS_READS: usize = 3;
    const FLAGS: [&str; 5] = ["-c", "-N", "50", "-p", "0.1"];

    fn env_on(name: &str) -> bool {
        matches!(std::env::var(name), Ok(v) if !v.is_empty() && v != "0")
    }

    /// `RUSTLE_SHARED_DEFINITION=1`: build the homology catalog with the shared definition.
    pub fn enabled() -> bool {
        env_on("RUSTLE_SHARED_DEFINITION")
    }

    /// `RUSTLE_SD_READ_LOCUS_SPLIT=1`: split chained read groups before adding read-locus nodes (Addendum AF-3).
    pub fn split_enabled() -> bool {
        env_on("RUSTLE_SD_READ_LOCUS_SPLIT")
    }

    /// Default junction-support floor for read-isoform widening (§6m0: k = 5 keeps FAMILY R at the annotated
    /// ceiling 0.963 and lifts FAMILY F strict 0.450 -> 0.515 on the development substrate; lower k buys more
    /// isoforms at the cost of precision, higher k the reverse).
    pub const ISOFORM_MIN_READS: u64 = 5;

    /// Read-isoform widening is ON by default; `RUSTLE_SD_READ_ISOFORM=0` restores the single-representative node.
    pub fn isoform_enabled() -> bool {
        !matches!(std::env::var("RUSTLE_SD_READ_ISOFORM"), Ok(v) if v == "0")
    }

    /// `RUSTLE_SD_ISOFORM_K`: junction-support floor k, default [`ISOFORM_MIN_READS`].
    pub fn isoform_k() -> u64 {
        std::env::var("RUSTLE_SD_ISOFORM_K")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(ISOFORM_MIN_READS)
    }

    /// One node of the copy graph: a gene-level locus with its exons and the exons of its representative transcript.
    #[derive(Clone, Debug, PartialEq, Eq)]
    pub struct SdNode {
        pub chrom: String,
        pub strand: char,
        pub n_reads: u64,
        pub exons: Vec<(u64, u64)>,
        pub rep_exons: Vec<(u64, u64)>,
        /// Spliced query chains for this node: the representative chain, plus every read-supported isoform
        /// admitted by [`widen_with_read_isoforms`] (§6m0). Always non-empty; `[rep_exons]` when widening is off.
        pub tx_chains: Vec<Vec<(u64, u64)>>,
    }

    impl SdNode {
        pub fn start(&self) -> u64 {
            self.exons.first().map(|e| e.0).unwrap_or(0)
        }
        pub fn end(&self) -> u64 {
            self.exons.last().map(|e| e.1).unwrap_or(0)
        }
    }

    /// Sort and merge intervals; touching intervals (`s <= last.end`) merge, as the prototype's `merge`.
    pub fn merge(mut iv: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
        iv.sort();
        let mut out: Vec<(u64, u64)> = Vec::with_capacity(iv.len());
        for (s, e) in iv {
            match out.last_mut() {
                Some(last) if s <= last.1 => last.1 = last.1.max(e),
                _ => out.push((s, e)),
            }
        }
        out
    }

    fn exon_sum(ex: &[(u64, u64)]) -> u64 {
        ex.iter().map(|(s, e)| e - s).sum()
    }

    struct UnionFind(Vec<usize>);

    impl UnionFind {
        fn new(n: usize) -> Self {
            UnionFind((0..n).collect())
        }
        fn find(&mut self, mut x: usize) -> usize {
            while self.0[x] != x {
                self.0[x] = self.0[self.0[x]];
                x = self.0[x];
            }
            x
        }
        /// `parent[find(a)] = find(b)`, the prototype's orientation.
        fn attach(&mut self, a: usize, b: usize) {
            let (ra, rb) = (self.find(a), self.find(b));
            self.0[ra] = rb;
        }
    }

    /// Connected components of intervals that overlap by >= 1 bp, per key, by the prototype's sweep; groups are
    /// returned in order of their first member index.
    fn overlap_groups<K: Ord + Clone>(items: &[(K, Vec<(u64, u64)>)]) -> Vec<Vec<usize>> {
        let mut uf = UnionFind::new(items.len());
        let mut by: BTreeMap<K, Vec<(u64, u64, usize)>> = BTreeMap::new();
        for (i, (k, blocks)) in items.iter().enumerate() {
            for &(s, e) in blocks {
                by.entry(k.clone()).or_default().push((s, e, i));
            }
        }
        for iv in by.values_mut() {
            iv.sort();
            let (mut end, mut owner): (i64, Option<usize>) = (-1, None);
            for &(s, e, i) in iv.iter() {
                if let Some(o) = owner {
                    if (s as i64) < end {
                        uf.attach(i, o);
                    }
                }
                if (e as i64) > end {
                    end = e as i64;
                    owner = Some(i);
                }
            }
        }
        let mut order: Vec<usize> = Vec::new();
        let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..items.len() {
            let r = uf.find(i);
            groups
                .entry(r)
                .or_insert_with(|| {
                    order.push(r);
                    Vec::new()
                })
                .push(i);
        }
        order
            .into_iter()
            .map(|r| groups.remove(&r).unwrap_or_default())
            .collect()
    }

    /// Nodes as the catalog emits them: one per rep, exons from its intron chain (the node dump's `exons`).
    pub fn nodes_from_reps(reps: &[DenovoTranscript]) -> Vec<SdNode> {
        reps.iter()
            .map(|r| {
                let mut exons =
                    crate::family::catalog_input::exon_blocks(r.start, r.end, &r.introns);
                exons.sort();
                SdNode {
                    chrom: r.chrom.clone(),
                    strand: r.strand,
                    n_reads: r.n_reads as u64,
                    tx_chains: vec![exons.clone()],
                    rep_exons: exons.clone(),
                    exons,
                }
            })
            .collect()
    }

    /// Addendum AB consolidation to gene-level loci.
    pub fn consolidate(ab1: &[SdNode]) -> Vec<SdNode> {
        let mut pieces: Vec<SdNode> = Vec::new();
        for n in ab1 {
            if n.exons.is_empty() {
                continue;
            }
            let mut cur = vec![n.exons[0]];
            for &b in &n.exons[1..] {
                if b.0 > cur.last().unwrap().1 && b.0 - cur.last().unwrap().1 > MAX_INTRON {
                    pieces.push(SdNode {
                        exons: cur.clone(),
                        rep_exons: cur.clone(),
                        tx_chains: vec![cur.clone()],
                        ..n.clone()
                    });
                    cur = vec![b];
                } else {
                    cur.push(b);
                }
            }
            pieces.push(SdNode {
                rep_exons: cur.clone(),
                tx_chains: vec![cur.clone()],
                exons: cur,
                ..n.clone()
            });
        }
        pieces.retain(|p| exon_sum(&p.exons) >= MIN_PIECE);
        let items: Vec<((String, char), Vec<(u64, u64)>)> = pieces
            .iter()
            .map(|p| ((p.chrom.clone(), p.strand), p.exons.clone()))
            .collect();
        let mut out: Vec<SdNode> = overlap_groups(&items)
            .into_iter()
            .map(|g| {
                let mut rep = &pieces[g[0]];
                for &i in &g[1..] {
                    let p = &pieces[i];
                    if (p.n_reads, exon_sum(&p.exons)) > (rep.n_reads, exon_sum(&rep.exons)) {
                        rep = p;
                    }
                }
                SdNode {
                    chrom: rep.chrom.clone(),
                    strand: rep.strand,
                    n_reads: g.iter().map(|&i| pieces[i].n_reads).sum(),
                    exons: merge(
                        g.iter()
                            .flat_map(|&i| pieces[i].exons.iter().copied())
                            .collect(),
                    ),
                    tx_chains: vec![rep.exons.clone()],
                    rep_exons: rep.exons.clone(),
                }
            })
            .collect();
        out.sort_by(|a, b| {
            (a.chrom.as_str(), a.start(), a.strand).cmp(&(b.chrom.as_str(), b.start(), b.strand))
        });
        out
    }

    /// One read's exon blocks (CIGAR segments between `N`), chromosome and transcript strand.
    #[derive(Clone, Debug)]
    pub struct ReadBlocks {
        pub chrom: String,
        pub strand: char,
        pub blocks: Vec<(u64, u64)>,
    }

    /// Primary, MAPQ >= 1 reads (the expression gate's reads); strand = read orientation, flipped by `ts:A:-`.
    pub fn read_blocks(reads: &[BamRead]) -> Vec<ReadBlocks> {
        reads
            .iter()
            .filter(|r| !r.is_secondary && !r.is_supplementary && r.mapq >= 1)
            .map(|r| {
                let mut strand = if r.reverse { '-' } else { '+' };
                if r.ts == Some('-') {
                    strand = if strand == '-' { '+' } else { '-' };
                }
                let (mut blocks, mut pos, mut cur) =
                    (Vec::new(), r.read.ref_start, r.read.ref_start);
                for &(op, n) in &r.read.cigar {
                    match op {
                        'M' | 'D' | '=' | 'X' => pos += n,
                        'N' => {
                            if pos > cur {
                                blocks.push((cur, pos));
                            }
                            pos += n;
                            cur = pos;
                        }
                        _ => {}
                    }
                }
                if pos > cur {
                    blocks.push((cur, pos));
                }
                ReadBlocks {
                    chrom: r.chrom.clone(),
                    strand,
                    blocks,
                }
            })
            .collect()
    }

    /// Bases covered by >= 2 reads' blocks, merged (at equal coordinates block ends sort before starts).
    pub fn depth2_exons(block_lists: &[&[(u64, u64)]]) -> Vec<(u64, u64)> {
        let mut ev: Vec<(u64, i8)> = Vec::new();
        for blocks in block_lists {
            for &(s, e) in blocks.iter() {
                ev.push((s, 1));
                ev.push((e, -1));
            }
        }
        ev.sort();
        let (mut depth, mut start, mut ex) = (0i64, None::<u64>, Vec::new());
        for (x, d) in ev {
            let prev = depth;
            depth += d as i64;
            if prev < 2 && depth >= 2 {
                start = Some(x);
            } else if prev >= 2 && depth < 2 {
                if let Some(s) = start {
                    if x > s {
                        ex.push((s, x));
                    }
                }
            }
        }
        merge(ex)
    }

    /// Addendum AF-3: components of segments linked by >= 2 reads overlapping both, with their supporting read counts,
    /// in order of each component's first segment.
    pub fn split_linked(
        segments: &[(u64, u64)],
        block_lists: &[&[(u64, u64)]],
    ) -> Vec<(Vec<usize>, usize)> {
        let touch: Vec<Vec<usize>> = block_lists
            .iter()
            .map(|blocks| {
                let set: BTreeSet<usize> = segments
                    .iter()
                    .enumerate()
                    .filter(|(_, &(s, e))| blocks.iter().any(|&(b0, b1)| b0 < e && s < b1))
                    .map(|(k, _)| k)
                    .collect();
                set.into_iter().collect()
            })
            .collect();
        let mut link: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        for t in &touch {
            for i in 0..t.len() {
                for j in (i + 1)..t.len() {
                    *link.entry((t[i], t[j])).or_insert(0) += 1;
                }
            }
        }
        let mut uf = UnionFind::new(segments.len());
        for (&(i, j), &n) in &link {
            if n >= 2 {
                uf.attach(j, i);
            }
        }
        let mut comps: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for k in 0..segments.len() {
            let r = uf.find(k);
            comps.entry(r).or_default().push(k);
        }
        let mut out: Vec<Vec<usize>> = comps.into_values().collect();
        out.sort_by_key(|v| v[0]);
        out.into_iter()
            .map(|ks| {
                let sup = touch
                    .iter()
                    .filter(|t| t.iter().any(|k| ks.contains(k)))
                    .count();
                (ks, sup)
            })
            .collect()
    }

    /// Exon intervals of a node set, for overlap queries.
    pub struct ExonIndex {
        by: HashMap<String, Vec<(u64, u64, usize)>>,
    }

    impl ExonIndex {
        pub fn new(nodes: &[SdNode]) -> Self {
            let mut by: HashMap<String, Vec<(u64, u64, usize)>> = HashMap::new();
            for (i, n) in nodes.iter().enumerate() {
                for &(s, e) in &n.exons {
                    by.entry(n.chrom.clone()).or_default().push((s, e, i));
                }
            }
            for v in by.values_mut() {
                v.sort();
            }
            ExonIndex { by }
        }
        /// Nodes with an exon overlapping `[s, e)` by >= 1 bp.
        pub fn hits(&self, chrom: &str, s: u64, e: u64) -> BTreeSet<usize> {
            let mut out = BTreeSet::new();
            if let Some(v) = self.by.get(chrom) {
                let hi = v.partition_point(|x| x.0 < e);
                for &(x0, x1, i) in &v[..hi] {
                    if x1 > s && x0 < e {
                        out.insert(i);
                    }
                }
            }
            out
        }
    }

    /// Addendum AC read-locus nodes (AF-3 split when `split`), merged into `base` and sorted as the prototype.
    pub fn with_read_locus_nodes(
        base: &[SdNode],
        reads: &[ReadBlocks],
        split: bool,
    ) -> (Vec<SdNode>, usize) {
        let items: Vec<((String, char), Vec<(u64, u64)>)> = reads
            .iter()
            .map(|r| ((r.chrom.clone(), r.strand), r.blocks.clone()))
            .collect();
        let bidx = ExonIndex::new(base);
        let mut added: Vec<SdNode> = Vec::new();
        for members in overlap_groups(&items) {
            if members.len() < MIN_LOCUS_READS {
                continue;
            }
            let (chrom, strand) = (reads[members[0]].chrom.clone(), reads[members[0]].strand);
            let lists: Vec<&[(u64, u64)]> = members
                .iter()
                .map(|&i| reads[i].blocks.as_slice())
                .collect();
            let ex = depth2_exons(&lists);
            let consider = |sub: Vec<(u64, u64)>, n: usize, added: &mut Vec<SdNode>| {
                if n < MIN_LOCUS_READS || exon_sum(&sub) < MIN_PIECE {
                    return;
                }
                if sub
                    .iter()
                    .any(|&(s, e)| !bidx.hits(&chrom, s, e).is_empty())
                {
                    return;
                }
                added.push(SdNode {
                    chrom: chrom.clone(),
                    strand,
                    n_reads: n as u64,
                    tx_chains: vec![sub.clone()],
                    rep_exons: sub.clone(),
                    exons: sub,
                });
            };
            if split {
                for (ks, sup) in split_linked(&ex, &lists) {
                    consider(ks.iter().map(|&k| ex[k]).collect(), sup, &mut added);
                }
            } else {
                consider(ex, members.len(), &mut added);
            }
        }
        let n_added = added.len();
        let mut nodes: Vec<SdNode> = base.to_vec();
        nodes.extend(added);
        nodes.sort_by(|a, b| {
            (a.chrom.as_str(), a.start(), a.strand).cmp(&(b.chrom.as_str(), b.start(), b.strand))
        });
        (nodes, n_added)
    }

    /// Read-isoform widening (§6m0, ledger 2026-09-17). For each node, every intron chain observed in its own
    /// same-strand primary reads whose junctions are EACH carried by >= `k` reads becomes a spliced query, and its
    /// blocks are merged into the node's exon union. A node can only widen: its previous exons and representative
    /// chain are always kept, and no node is created, removed or merged here.
    ///
    /// Why: the catalog emits one representative per locus before [`consolidate`] runs, so a node's exon set is one
    /// isoform (705/710 nodes on the human ideal substrate). Widening lifts full-length NPIP copies from 5/27 to
    /// 17/27 at k = 5 with FAMILY R held at the annotated arm's own 0.963.
    pub fn widen_with_read_isoforms(
        nodes: &[SdNode],
        reads: &[ReadBlocks],
        k: u64,
    ) -> (Vec<SdNode>, usize) {
        let idx = ExonIndex::new(nodes);
        // Per node: chain (intron vector) -> (support, min start, max end).
        let mut per_node: Vec<BTreeMap<Vec<(u64, u64)>, (u64, u64, u64)>> =
            vec![BTreeMap::new(); nodes.len()];
        for r in reads {
            if r.blocks.len() < 2 {
                continue;
            }
            let introns: Vec<(u64, u64)> = r.blocks.windows(2).map(|w| (w[0].1, w[1].0)).collect();
            let (first, last) = (r.blocks[0].0, r.blocks[r.blocks.len() - 1].1);
            let mut seen: BTreeSet<usize> = BTreeSet::new();
            for &(s, e) in &r.blocks {
                for v in idx.hits(&r.chrom, s, e) {
                    if nodes[v].strand == r.strand {
                        seen.insert(v);
                    }
                }
            }
            for v in seen {
                let ent = per_node[v]
                    .entry(introns.clone())
                    .or_insert((0, u64::MAX, 0));
                ent.0 += 1;
                ent.1 = ent.1.min(first);
                ent.2 = ent.2.max(last);
            }
        }
        let mut widened = 0usize;
        let out: Vec<SdNode> = nodes
            .iter()
            .enumerate()
            .map(|(v, n)| {
                // A junction is supported when >= k reads of this node carry it, counted over all of the node's chains.
                let mut junc: BTreeMap<(u64, u64), u64> = BTreeMap::new();
                for (chain, &(sup, _, _)) in &per_node[v] {
                    for &j in chain {
                        *junc.entry(j).or_insert(0) += sup;
                    }
                }
                let mut chains: Vec<Vec<(u64, u64)>> = n.tx_chains.clone();
                let mut exons = n.exons.clone();
                let mut added = false;
                for (chain, &(_, first, last)) in &per_node[v] {
                    if chain.iter().any(|j| junc.get(j).copied().unwrap_or(0) < k) {
                        continue;
                    }
                    let blocks = crate::family::catalog_input::exon_blocks(first, last, chain);
                    if blocks.is_empty() || chains.iter().any(|c| *c == blocks) {
                        continue;
                    }
                    exons.extend(blocks.iter().copied());
                    chains.push(blocks);
                    added = true;
                }
                if !added {
                    return n.clone();
                }
                widened += 1;
                chains.sort();
                chains.dedup();
                SdNode {
                    exons: merge(exons),
                    tx_chains: chains,
                    ..n.clone()
                }
            })
            .collect();
        (out, widened)
    }

    /// One PAF record (the fields the finders use).
    #[derive(Clone, Debug)]
    pub struct PafRec {
        pub q: String,
        pub qlen: u64,
        pub qs: u64,
        pub qe: u64,
        pub strand: char,
        pub chrom: String,
        pub clen: u64,
        pub ts: u64,
        pub te: u64,
        pub nm: u64,
        pub bl: u64,
        pub cg: String,
    }

    pub fn parse_paf(text: &str) -> Vec<PafRec> {
        text.lines()
            .filter_map(|line| {
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() < 12 {
                    return None;
                }
                Some(PafRec {
                    q: f[0].to_string(),
                    qlen: f[1].parse().ok()?,
                    qs: f[2].parse().ok()?,
                    qe: f[3].parse().ok()?,
                    strand: f[4].chars().next()?,
                    chrom: f[5].to_string(),
                    clen: f[6].parse().ok()?,
                    ts: f[7].parse().ok()?,
                    te: f[8].parse().ok()?,
                    nm: f[9].parse().ok()?,
                    bl: f[10].parse().ok()?,
                    cg: f[12..]
                        .iter()
                        .find_map(|x| x.strip_prefix("cg:Z:"))
                        .unwrap_or("")
                        .to_string(),
                })
            })
            .collect()
    }

    /// A PAF `cg:Z:` CIGAR string as `(length, op)` runs, e.g. "88=3D510=1I436=" -> [(88,'='),(3,'D'),(510,'='),
    /// (1,'I'),(436,'=')]. Also used by `missing_copy_flag_pass`.
    pub(crate) fn cigar_ops(cg: &str) -> Vec<(u64, char)> {
        let mut out = Vec::new();
        let mut n: u64 = 0;
        for c in cg.chars() {
            if let Some(d) = c.to_digit(10) {
                n = n * 10 + d as u64;
            } else {
                out.push((n, c));
                n = 0;
            }
        }
        out
    }

    /// Spliced hits passing identity >= 0.80 (matches / block) and query coverage >= 0.50.
    pub fn transcript_hits(recs: &[PafRec]) -> Vec<PafRec> {
        recs.iter()
            .filter(|h| {
                h.bl > 0
                    && h.qlen > 0
                    && h.nm as f64 / h.bl as f64 >= MIN_ID
                    && (h.qe - h.qs) as f64 / h.qlen as f64 >= MIN_COV
            })
            .cloned()
            .collect()
    }

    /// Genomic exon blocks of a spliced alignment (split at `N`).
    pub fn tx_exon_blocks(h: &PafRec) -> Vec<(u64, u64)> {
        let (mut blocks, mut pos, mut cur) = (Vec::new(), h.ts, h.ts);
        for (n, op) in cigar_ops(&h.cg) {
            match op {
                'M' | '=' | 'X' | 'D' => pos += n,
                'N' => {
                    if pos > cur {
                        blocks.push((cur, pos));
                    }
                    pos += n;
                    cur = pos;
                }
                _ => {}
            }
        }
        if pos > cur {
            blocks.push((cur, pos));
        }
        blocks
    }

    /// A gene-body chain: query, chromosome, strand and its target interval.
    #[derive(Clone, Debug, PartialEq, Eq)]
    pub struct Chain {
        pub q: String,
        pub chrom: String,
        pub strand: char,
        pub s: u64,
        pub e: u64,
    }

    /// The prototype's `gene_body_chains`: records of one (query, chromosome, strand) chained in target order while
    /// the gap stays <= the query length, the span <= twice it and the query advances; a chain is kept at identity
    /// >= 0.80 with aligned query >= 0.50 of min(query length, extrapolated target span).
    pub fn gene_body_chains(recs: &[PafRec]) -> Vec<Chain> {
        let mut order: Vec<(String, String, char)> = Vec::new();
        let mut by: HashMap<(String, String, char), Vec<&PafRec>> = HashMap::new();
        for r in recs {
            let k = (r.q.clone(), r.chrom.clone(), r.strand);
            by.entry(k.clone())
                .or_insert_with(|| {
                    order.push(k);
                    Vec::new()
                })
                .push(r);
        }
        let mut chains = Vec::new();
        for k in order {
            let mut rs = by.remove(&k).unwrap_or_default();
            rs.sort_by_key(|r| (r.ts, r.te));
            let lq = rs[0].qlen as i64;
            let strand = k.2;
            let mut cur: Vec<&PafRec> = Vec::new();
            let emit = |cur: &Vec<&PafRec>, chains: &mut Vec<Chain>| {
                let nm: u64 = cur.iter().map(|x| x.nm).sum();
                let bl: u64 = cur.iter().map(|x| x.bl).sum();
                let qiv = merge(cur.iter().map(|x| (x.qs, x.qe)).collect());
                let aligned: i64 = qiv.iter().map(|(s, e)| (e - s) as i64).sum();
                let ts = cur.iter().map(|x| x.ts).min().unwrap() as i64;
                let te = cur.iter().map(|x| x.te).max().unwrap() as i64;
                let (q0, q1) = (qiv[0].0 as i64, qiv.last().unwrap().1 as i64);
                let (xs, xe) = if strand == '+' {
                    (ts - q0, te + (lq - q1))
                } else {
                    (ts - (lq - q1), te + q0)
                };
                let (xs, xe) = (xs.max(0), xe.min(cur[0].clen as i64));
                if bl > 0
                    && nm as f64 / bl as f64 >= MIN_ID
                    && aligned as f64 >= MIN_COV * (lq.min(xe - xs)) as f64
                {
                    chains.push(Chain {
                        q: k.0.clone(),
                        chrom: k.1.clone(),
                        strand,
                        s: ts as u64,
                        e: te as u64,
                    });
                }
            };
            for r in rs {
                if !cur.is_empty() {
                    let gap = r.ts as i64 - cur.iter().map(|x| x.te).max().unwrap() as i64;
                    let span = r.te as i64 - cur.iter().map(|x| x.ts).min().unwrap() as i64;
                    let last_qs = cur.last().unwrap().qs;
                    let ordered = if strand == '+' {
                        r.qs >= last_qs
                    } else {
                        r.qs <= last_qs
                    };
                    if !(gap <= lq && span <= 2 * lq && ordered) {
                        emit(&cur, &mut chains);
                        cur.clear();
                    }
                }
                cur.push(r);
            }
            if !cur.is_empty() {
                emit(&cur, &mut chains);
            }
        }
        chains
    }

    /// Query keys shared by nodes with the same representative transcript / gene body (the prototype's dedupe keys).
    pub fn tx_key(n: &SdNode) -> String {
        tx_key_for(n, &n.rep_exons)
    }

    /// The query key of one spliced chain of a node (the same form as [`tx_key`], which is this for `rep_exons`).
    pub fn tx_key_for(n: &SdNode, chain: &[(u64, u64)]) -> String {
        let ex: Vec<String> = chain.iter().map(|(a, b)| format!("{a}-{b}")).collect();
        format!("{}|{}|{}", n.chrom, n.strand, ex.join(","))
    }

    pub fn body_key(n: &SdNode) -> String {
        format!("{}|{}|{}", n.chrom, n.start(), n.end())
    }

    /// Symmetric node pairs joined by an exon edge or a gene-body edge, with counts per kind.
    pub fn edges(
        nodes: &[SdNode],
        tx_by_key: &HashMap<String, Vec<PafRec>>,
        chains_by_key: &HashMap<String, Vec<Chain>>,
    ) -> (BTreeSet<(usize, usize)>, usize, usize) {
        let idx = ExonIndex::new(nodes);
        let spliced: Vec<bool> = nodes.iter().map(|n| n.exons.len() >= 2).collect();
        let mut kinds: BTreeSet<((usize, usize), u8)> = BTreeSet::new();
        for (u, nu) in nodes.iter().enumerate() {
            let (us, ue) = (nu.start(), nu.end());
            // (kind, chrom, s, e, hit_s, hit_e, orientation on the target)
            let mut cand: Vec<(u8, &str, u64, u64, u64, u64, char)> = Vec::new();
            for chain in &nu.tx_chains {
                for h in tx_by_key
                    .get(&tx_key_for(nu, chain))
                    .map(|v| v.as_slice())
                    .unwrap_or(&[])
                {
                    for (s, e) in tx_exon_blocks(h) {
                        cand.push((0, &h.chrom, s, e, h.ts, h.te, h.strand));
                    }
                }
            }
            for c in chains_by_key
                .get(&body_key(nu))
                .map(|v| v.as_slice())
                .unwrap_or(&[])
            {
                let orient = if c.strand == '+' {
                    nu.strand
                } else if nu.strand == '+' {
                    '-'
                } else if nu.strand == '-' {
                    '+'
                } else {
                    nu.strand
                };
                cand.push((1, &c.chrom, c.s, c.e, c.s, c.e, orient));
            }
            for (kind, chrom, s, e, hs, he, orient) in cand {
                if chrom == nu.chrom && hs < ue && us < he {
                    continue;
                }
                for v in idx.hits(chrom, s, e) {
                    if v == u || (spliced[u] && spliced[v] && nodes[v].strand != orient) {
                        continue;
                    }
                    kinds.insert(((u.min(v), u.max(v)), kind));
                }
            }
        }
        let n_exon = kinds.iter().filter(|(_, k)| *k == 0).count();
        let n_body = kinds.len() - n_exon;
        (kinds.into_iter().map(|(p, _)| p).collect(), n_exon, n_body)
    }

    /// Triangle-supported leader neighbourhoods: nodes in order of reads (desc), degree (desc), chromosome, start and
    /// index; an unassigned node with unassigned neighbours leads a family of itself, those neighbours, and unassigned
    /// nodes adjacent to >= 2 members of that star.
    pub fn triangle_leaders(nodes: &[SdNode], pairs: &BTreeSet<(usize, usize)>) -> Vec<Vec<usize>> {
        let mut adj: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); nodes.len()];
        for &(a, b) in pairs {
            adj[a].insert(b);
            adj[b].insert(a);
        }
        let mut order: Vec<usize> = (0..nodes.len()).collect();
        order.sort_by(|&a, &b| {
            (
                std::cmp::Reverse(nodes[a].n_reads),
                std::cmp::Reverse(adj[a].len()),
                nodes[a].chrom.as_str(),
                nodes[a].start(),
                a,
            )
                .cmp(&(
                    std::cmp::Reverse(nodes[b].n_reads),
                    std::cmp::Reverse(adj[b].len()),
                    nodes[b].chrom.as_str(),
                    nodes[b].start(),
                    b,
                ))
        });
        let mut seen = vec![false; nodes.len()];
        let mut fams = Vec::new();
        for i in order {
            if seen[i] {
                continue;
            }
            let free: Vec<usize> = adj[i].iter().copied().filter(|&y| !seen[y]).collect();
            if free.is_empty() {
                continue;
            }
            let mut star: BTreeSet<usize> = free.iter().copied().collect();
            star.insert(i);
            let second: BTreeSet<usize> = star
                .iter()
                .flat_map(|&y| adj[y].iter().copied())
                .filter(|&z| {
                    !seen[z] && !star.contains(&z) && adj[z].intersection(&star).count() >= 2
                })
                .collect();
            let fam: BTreeSet<usize> = star.union(&second).copied().collect();
            for &x in &fam {
                seen[x] = true;
            }
            fams.push(fam.into_iter().collect::<Vec<usize>>());
        }
        fams.retain(|f| f.len() >= 2);
        fams.sort_by(|a, b| {
            (std::cmp::Reverse(a.len()), a[0]).cmp(&(std::cmp::Reverse(b.len()), b[0]))
        });
        fams
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|b| match b.to_ascii_uppercase() {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect()
    }

    fn spliced_seq(
        genome: &GenomeIndex,
        chrom: &str,
        exons: &[(u64, u64)],
        strand: char,
    ) -> Vec<u8> {
        let mut s: Vec<u8> = Vec::new();
        for &(a, b) in exons {
            if let Some(x) = genome.fetch_sequence(chrom, a, b) {
                s.extend(x.iter().map(|c| c.to_ascii_uppercase()));
            }
        }
        if strand == '-' {
            revcomp(&s)
        } else {
            s
        }
    }

    fn run_minimap2(
        minimap2: &str,
        preset: &[&str],
        threads: usize,
        target: &std::path::Path,
        query: &std::path::Path,
    ) -> Result<String> {
        let out = std::process::Command::new(minimap2)
            .args(FLAGS)
            .args(preset)
            .arg("-t")
            .arg(threads.max(1).to_string())
            .arg(target)
            .arg(query)
            .output()
            .with_context(|| format!("running {minimap2}"))?;
        anyhow::ensure!(
            out.status.success(),
            "minimap2 {:?} failed: {}",
            preset,
            String::from_utf8_lossy(&out.stderr)
        );
        Ok(String::from_utf8_lossy(&out.stdout).into_owned())
    }

    /// The full construction on a catalog's reps and reads. Returns the nodes, the families (node indices) and the
    /// edge pairs.
    pub fn build(
        reps: &[DenovoTranscript],
        reads: &[ReadBlocks],
        genome: &GenomeIndex,
        minimap2: &str,
        threads: usize,
    ) -> Result<(Vec<SdNode>, Vec<Vec<usize>>, BTreeSet<(usize, usize)>)> {
        let ab2 = consolidate(&nodes_from_reps(reps));
        let (nodes, n_added) = with_read_locus_nodes(&ab2, reads, split_enabled());
        eprintln!(
            "[shared-definition] {} reps -> {} gene-level loci + {} read-locus nodes{} = {} nodes",
            reps.len(),
            ab2.len(),
            n_added,
            if split_enabled() { " (split)" } else { "" },
            nodes.len()
        );
        let nodes = if isoform_enabled() {
            let k = isoform_k();
            let (w, n_widened) = widen_with_read_isoforms(&nodes, reads, k);
            let queries: usize = w.iter().map(|n| n.tx_chains.len()).sum();
            eprintln!(
                "[shared-definition] read-isoform widening k={k}: {n_widened} of {} nodes widened, {queries} spliced queries",
                w.len()
            );
            w
        } else {
            eprintln!("[shared-definition] read-isoform widening OFF (RUSTLE_SD_READ_ISOFORM=0)");
            nodes
        };
        let dir =
            std::env::temp_dir().join(format!("rustle_sd_{}_{}", std::process::id(), reps.len()));
        std::fs::create_dir_all(&dir)?;
        let (target, txfa, bodyfa) = (
            dir.join("target.fa"),
            dir.join("tx.fa"),
            dir.join("body.fa"),
        );
        {
            // Target contig order: `RUSTLE_SD_TARGET_ORDER` (comma list) when set, otherwise sorted names.
            let mut names: Vec<String> = genome.chroms().map(|(c, _)| c.to_string()).collect();
            names.sort();
            if let Ok(v) = std::env::var("RUSTLE_SD_TARGET_ORDER") {
                let want: Vec<String> = v
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect();
                let mut rest: Vec<String> = names
                    .iter()
                    .filter(|n| !want.contains(n))
                    .cloned()
                    .collect();
                names = want
                    .into_iter()
                    .filter(|w| genome.chrom_len(w) > 0)
                    .collect();
                names.append(&mut rest);
            }
            let mut fh = std::io::BufWriter::new(std::fs::File::create(&target)?);
            for c in &names {
                if let Some(seq) = genome.fetch_sequence(c, 0, genome.chrom_len(c)) {
                    writeln!(fh, ">{c}")?;
                    fh.write_all(&seq)?;
                    writeln!(fh)?;
                }
            }
        }
        let (mut tx_seen, mut body_seen) = (BTreeSet::new(), BTreeSet::new());
        {
            let mut ft = std::io::BufWriter::new(std::fs::File::create(&txfa)?);
            let mut fb = std::io::BufWriter::new(std::fs::File::create(&bodyfa)?);
            for n in &nodes {
                let bk = body_key(n);
                for chain in &n.tx_chains {
                    let tk = tx_key_for(n, chain);
                    if tx_seen.insert(tk.clone()) {
                        writeln!(ft, ">{tk}")?;
                        ft.write_all(&spliced_seq(genome, &n.chrom, chain, n.strand))?;
                        writeln!(ft)?;
                    }
                }
                if body_seen.insert(bk.clone()) {
                    writeln!(fb, ">{bk}")?;
                    let body = genome
                        .fetch_sequence(&n.chrom, n.start(), n.end())
                        .unwrap_or_default();
                    fb.write_all(
                        &body
                            .iter()
                            .map(|c| c.to_ascii_uppercase())
                            .collect::<Vec<u8>>(),
                    )?;
                    writeln!(fb)?;
                }
            }
        }
        let tx_paf = run_minimap2(minimap2, &["-x", "splice", "-uf"], threads, &target, &txfa)?;
        let body_paf = run_minimap2(minimap2, &["-x", "asm20"], threads, &target, &bodyfa)?;
        let _ = std::fs::remove_dir_all(&dir);
        let mut tx_by_key: HashMap<String, Vec<PafRec>> = HashMap::new();
        for h in transcript_hits(&parse_paf(&tx_paf)) {
            tx_by_key.entry(h.q.clone()).or_default().push(h);
        }
        let mut chains_by_key: HashMap<String, Vec<Chain>> = HashMap::new();
        for c in gene_body_chains(&parse_paf(&body_paf)) {
            chains_by_key.entry(c.q.clone()).or_default().push(c);
        }
        let (pairs, n_exon, n_body) = edges(&nodes, &tx_by_key, &chains_by_key);
        let fams = triangle_leaders(&nodes, &pairs);
        eprintln!(
            "[shared-definition] edges: {n_exon} exon + {n_body} gene-body ({} pairs); {} triangle-supported families holding {} loci",
            pairs.len(), fams.len(), fams.iter().map(|f| f.len()).sum::<usize>()
        );
        Ok((nodes, fams, pairs))
    }

    /// A node as a catalog copy (sequence = its exons in transcript orientation).
    pub fn node_transcript(n: &SdNode, genome: &GenomeIndex) -> DenovoTranscript {
        DenovoTranscript {
            tid: format!("SD~{}_{}_{}", n.chrom, n.start(), n.end()),
            chrom: n.chrom.clone(),
            start: n.start(),
            end: n.end(),
            n_reads: n.n_reads.min(u32::MAX as u64) as u32,
            strand: n.strand,
            introns: n.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect(),
            seq: spliced_seq(genome, &n.chrom, &n.exons, n.strand),
            distinguishing_uniq: 0,
            core_bp: 0,
            stub: false,
            tes: None,
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn node(chrom: &str, strand: char, reads: u64, exons: &[(u64, u64)]) -> SdNode {
            SdNode {
                chrom: chrom.into(),
                strand,
                n_reads: reads,
                exons: exons.to_vec(),
                rep_exons: exons.to_vec(),
                tx_chains: vec![exons.to_vec()],
            }
        }

        #[test]
        fn read_isoform_widening_admits_supported_chains_and_only_widens() {
            // One node with a 2-exon representative; reads carry a second isoform with an extra exon.
            let n = node("chr1", '+', 10, &[(100, 200), (500, 600)]);
            let rb = |blocks: &[(u64, u64)]| ReadBlocks {
                chrom: "chr1".into(),
                strand: '+',
                blocks: blocks.to_vec(),
            };
            let mut reads: Vec<ReadBlocks> = Vec::new();
            for _ in 0..5 {
                reads.push(rb(&[(100, 200), (300, 350), (500, 600)]));
            }
            // A third chain seen twice only: below k, so it must not be admitted.
            for _ in 0..2 {
                reads.push(rb(&[(100, 200), (700, 800)]));
            }
            let (out, widened) = widen_with_read_isoforms(&[n.clone()], &reads, 5);
            assert_eq!(widened, 1);
            let w = &out[0];
            // Only widens: every original exon block survives, and the rep chain stays a query.
            for e in &n.exons {
                assert!(
                    w.exons.iter().any(|x| x.0 <= e.0 && e.1 <= x.1),
                    "{e:?} lost"
                );
            }
            assert!(w.tx_chains.contains(&n.rep_exons));
            assert_eq!(w.rep_exons, n.rep_exons);
            // The supported isoform is now a query and its exon is in the union.
            assert_eq!(w.tx_chains.len(), 2);
            assert!(w.exons.iter().any(|&(a, b)| a <= 300 && 350 <= b));
            // The 2-read chain contributed nothing.
            assert!(!w.exons.iter().any(|&(a, b)| a <= 700 && 800 <= b));
        }

        #[test]
        fn read_isoform_widening_is_inert_without_support_and_on_unspliced_reads() {
            let n = node("chr1", '-', 4, &[(100, 200), (500, 600)]);
            let unspliced = ReadBlocks {
                chrom: "chr1".into(),
                strand: '-',
                blocks: vec![(100, 600)],
            };
            let wrong_strand = ReadBlocks {
                chrom: "chr1".into(),
                strand: '+',
                blocks: vec![(100, 200), (300, 350), (500, 600)],
            };
            let reads: Vec<ReadBlocks> = std::iter::repeat(unspliced)
                .take(9)
                .chain(std::iter::repeat(wrong_strand).take(9))
                .collect();
            let (out, widened) = widen_with_read_isoforms(&[n.clone()], &reads, 5);
            assert_eq!(widened, 0);
            assert_eq!(out[0], n);
        }

        #[test]
        fn isoform_knobs_default_to_on_at_k_five() {
            if std::env::var("RUSTLE_SD_READ_ISOFORM").is_err() {
                assert!(isoform_enabled());
            }
            if std::env::var("RUSTLE_SD_ISOFORM_K").is_err() {
                assert_eq!(isoform_k(), 5);
                assert_eq!(ISOFORM_MIN_READS, 5);
            }
        }

        #[test]
        fn shared_definition_is_off_by_default() {
            if std::env::var("RUSTLE_SHARED_DEFINITION").is_err() {
                assert!(!enabled());
            }
        }

        #[test]
        fn merge_joins_touching_intervals() {
            assert_eq!(
                merge(vec![(5, 10), (0, 5), (20, 30), (25, 26)]),
                vec![(0, 10), (20, 30)]
            );
        }

        #[test]
        fn consolidate_cuts_long_introns_drops_small_pieces_and_merges_overlaps() {
            let a = node(
                "c",
                '+',
                5,
                &[(0, 200), (200 + MAX_INTRON + 1, 200 + MAX_INTRON + 50)],
            );
            let b = node("c", '+', 9, &[(150, 400)]);
            let c = node("c", '-', 1, &[(150, 400)]);
            let out = consolidate(&[a, b, c]);
            // the 49 bp tail piece is dropped; a's first piece and b merge on '+'; c stays on '-'
            assert_eq!(out.len(), 2);
            assert_eq!(out[0].exons, vec![(0, 400)]);
            assert_eq!(out[0].n_reads, 14);
            assert_eq!(out[0].rep_exons, vec![(150, 400)]);
            assert_eq!(out[1].strand, '-');
        }

        #[test]
        fn depth_and_split_match_the_prototype_toy() {
            let a: Vec<(u64, u64)> = vec![(0, 100), (200, 300)];
            let b: Vec<(u64, u64)> = vec![(1000, 1100)];
            let rt: Vec<(u64, u64)> = vec![(250, 300), (1000, 1050)];
            let lists: Vec<&[(u64, u64)]> = vec![&a, &a, &a, &b, &b, &b, &rt];
            let ex = depth2_exons(&lists);
            assert_eq!(ex, vec![(0, 100), (200, 300), (1000, 1100)]);
            assert_eq!(
                split_linked(&ex, &lists),
                vec![(vec![0, 1], 4), (vec![2], 4)]
            );
        }

        #[test]
        fn read_locus_nodes_skip_groups_touching_existing_nodes_unless_split() {
            let base = vec![node("c", '+', 10, &[(1000, 1100)])];
            let mk = |b: Vec<(u64, u64)>| ReadBlocks {
                chrom: "c".into(),
                strand: '+',
                blocks: b,
            };
            let mut reads = vec![mk(vec![(0, 100), (200, 300)]); 3];
            reads.extend(vec![mk(vec![(1000, 1100)]); 3]);
            reads.push(mk(vec![(250, 300), (1000, 1050)]));
            let (whole, n_whole) = with_read_locus_nodes(&base, &reads, false);
            assert_eq!(n_whole, 0);
            assert_eq!(whole.len(), 1);
            let (split, n_split) = with_read_locus_nodes(&base, &reads, true);
            assert_eq!(n_split, 1);
            assert_eq!(split[0].exons, vec![(0, 100), (200, 300)]);
            assert_eq!(split[0].n_reads, 4);
        }

        #[test]
        fn transcript_blocks_and_hits() {
            let h = PafRec {
                q: "q".into(),
                qlen: 830,
                qs: 0,
                qe: 830,
                strand: '+',
                chrom: "c".into(),
                clen: 5000,
                ts: 2000,
                te: 4230,
                nm: 830,
                bl: 830,
                cg: "300M600N250M800N280M".into(),
            };
            assert_eq!(
                tx_exon_blocks(&h),
                vec![(2000, 2300), (2900, 3150), (3950, 4230)]
            );
            assert_eq!(transcript_hits(&[h.clone()]).len(), 1);
            let low = PafRec { nm: 600, ..h };
            assert!(transcript_hits(&[low]).is_empty());
        }

        #[test]
        fn gene_body_chain_joins_ordered_records_and_applies_floors() {
            let r = |qs, qe, ts, te| PafRec {
                q: "b".into(),
                qlen: 1000,
                qs,
                qe,
                strand: '+',
                chrom: "c".into(),
                clen: 100_000,
                ts,
                te,
                nm: 95,
                bl: 100,
                cg: String::new(),
            };
            let chains =
                gene_body_chains(&[r(500, 1000, 10_600, 11_100), r(0, 400, 10_000, 10_400)]);
            assert_eq!(
                chains,
                vec![Chain {
                    q: "b".into(),
                    chrom: "c".into(),
                    strand: '+',
                    s: 10_000,
                    e: 11_100
                }]
            );
            // a lone 200 bp record covers < 0.50 of the 1,000 bp body
            assert!(gene_body_chains(&[r(0, 200, 50_000, 50_200)]).is_empty());
        }

        #[test]
        fn triangle_leaders_stop_at_one_hop_unless_two_star_members_support() {
            // leader 0 (most reads) - {1,2}; 3 touches 1 and 2 (joins); 4 touches only 2 (does not); 4-5 pair
            let nodes: Vec<SdNode> = (0..6)
                .map(|i| {
                    node(
                        "c",
                        '+',
                        if i == 0 { 50 } else { 5 },
                        &[(i * 1000, i * 1000 + 500)],
                    )
                })
                .collect();
            let pairs: BTreeSet<(usize, usize)> = [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (4, 5)]
                .into_iter()
                .collect();
            assert_eq!(
                triangle_leaders(&nodes, &pairs),
                vec![vec![0, 1, 2, 3], vec![4, 5]]
            );
        }

        #[test]
        fn edges_ignore_own_locus_and_wrong_orientation() {
            let nodes = vec![
                node("c", '+', 5, &[(0, 100), (200, 300)]),
                node("c", '+', 5, &[(10_000, 10_100), (10_200, 10_300)]),
                node("c", '-', 5, &[(20_000, 20_100), (20_200, 20_300)]),
            ];
            let hit = |ts: u64, te: u64, strand: char| PafRec {
                q: tx_key(&nodes[0]),
                qlen: 200,
                qs: 0,
                qe: 200,
                strand,
                chrom: "c".into(),
                clen: 100_000,
                ts,
                te,
                nm: 200,
                bl: 200,
                cg: "100M100N100M".into(),
            };
            let mut tx: HashMap<String, Vec<PafRec>> = HashMap::new();
            tx.insert(
                tx_key(&nodes[0]),
                vec![
                    hit(0, 300, '+'),
                    hit(10_000, 10_300, '+'),
                    hit(20_000, 20_300, '+'),
                ],
            );
            let (pairs, n_exon, n_body) = edges(&nodes, &tx, &HashMap::new());
            assert_eq!(pairs.into_iter().collect::<Vec<_>>(), vec![(0, 1)]);
            assert_eq!((n_exon, n_body), (1, 0));
        }
    }
}

pub mod seq_utils {
    //! Small sequence utilities the family-analysis modules depend on.
    //!
    //! Two reverse complements with DIFFERENT semantics live here, named by what they do:
    //! - `reverse_complement` — uppercase ACGT only; every other byte (lowercase included) maps to `N`.
    //!   Relocated verbatim from the retired assembler `vg.rs` (`docs/RETIREMENT_AND_MIGRATION.md`).
    //! - `revcomp_keep_case` — complements both cases (`N`/`n` kept), passes every other byte through
    //!   unchanged (Python `str.translate` semantics).
    //!
    //! plus `hw_distance` / `aln_id`, the HW (infix) edit distance and identity that equal edlib's
    //! `mode="HW"` bit for bit. `revcomp_keep_case`, `hw_distance` and `aln_id` were moved here from
    //! `bridge_detector.rs` (a port of `bench/recombination_bridge_detector.py`) when the rest of that module
    //! was removed as dead code (2026-09-24, tag `notebook-2026-09-24`).
    //!
    //! **STATUS:** INFRASTRUCTURE  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    /// Reverse-complement a nucleotide sequence. Non-ACGT bytes map to `b'N'`.
    pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
        seq.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'N',
            })
            .collect()
    }

    /// One-base complement, exactly the Python `_COMP = str.maketrans("ACGTacgtNn",
    /// "TGCAtgcaNn")`: A<->T, C<->G, G<->C, T<->A (both cases), N->N, n->n; EVERY OTHER
    /// byte is left unchanged (Python `str.translate` passes unmapped chars through).
    #[inline]
    fn comp_base_keep_case(b: u8) -> u8 {
        match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            b'N' => b'N',
            b'n' => b'n',
            other => other,
        }
    }

    /// Case-preserving reverse complement: complement every base via `comp_base_keep_case`, then reverse
    /// (`s.translate(_COMP)[::-1]`).
    pub fn revcomp_keep_case(s: &[u8]) -> Vec<u8> {
        s.iter().rev().map(|&b| comp_base_keep_case(b)).collect()
    }

    /// HW (infix) edit distance of `q` against `t`: the minimum edit distance of `q` to ANY
    /// substring of `t` (free gaps at both ends of `t`). Equals `edlib.align(q, t,
    /// mode="HW", task="distance")["editDistance"]`. Two rolling DP rows over `t` columns
    /// (`dp[0][*]=0`, `dp[i][0]=i`, answer `= min_j dp[|q|][j]`).
    pub fn hw_distance(q: &[u8], t: &[u8]) -> usize {
        let lq = q.len();
        let lt = t.len();
        // dp row 0 = 0 across all t columns (a free start position in t).
        let mut prev: Vec<usize> = vec![0; lt + 1];
        let mut cur: Vec<usize> = vec![0; lt + 1];
        for i in 1..=lq {
            cur[0] = i; // dp[i][0] = i
            let qi = q[i - 1];
            for j in 1..=lt {
                let sub = prev[j - 1] + if qi == t[j - 1] { 0 } else { 1 };
                let del = prev[j] + 1;
                let ins = cur[j - 1] + 1;
                cur[j] = sub.min(del).min(ins);
            }
            std::mem::swap(&mut prev, &mut cur);
        }
        // answer = min over j of dp[|q|][j]  (free end position in t). After the final swap
        // the last computed row is in `prev`. For lq == 0, prev is the all-zero row 0.
        *prev.iter().min().expect("row has >= 1 column")
    }

    /// Best HW (infix) identity of `q` inside `t`, trying `q` AND `revcomp_keep_case(q)`.
    /// `id = 1 - min_dist / len(q)`; `0.0` if either string is empty
    /// (`recombination_bridge_detector.py:70`).
    pub fn aln_id(q: &[u8], t: &[u8]) -> f64 {
        if q.is_empty() || t.is_empty() {
            return 0.0;
        }
        let d_fwd = hw_distance(q, t);
        let d_rev = hw_distance(&revcomp_keep_case(q), t);
        let best = d_fwd.min(d_rev);
        // 1.0 - best / len(q)  (single f64 division, matches Python exactly).
        1.0 - best as f64 / q.len() as f64
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use std::collections::BTreeSet;

        #[test]
        fn reverse_complement_matches_legacy_semantics() {
            assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
            assert_eq!(reverse_complement(b"AACG"), b"CGTT");
            // non-ACGT (incl. lowercase) -> N; reverse THEN complement: [A,N,a] -> [a,N,A] -> [N,N,T]
            assert_eq!(reverse_complement(b"ANa"), b"NNT");
        }

        #[test]
        fn revcomp_table_semantics() {
            // A<->T, C<->G, N->N (both cases), everything else passes through, then reverse.
            assert_eq!(revcomp_keep_case(b"ACGT"), b"ACGT".to_vec()); // palindrome
            assert_eq!(revcomp_keep_case(b"AAAA"), b"TTTT".to_vec());
            assert_eq!(revcomp_keep_case(b"ACGTN"), b"NACGT".to_vec()); // N->N, reversed
            assert_eq!(revcomp_keep_case(b"acgt"), b"acgt".to_vec()); // lowercase palindrome
            assert_eq!(revcomp_keep_case(b"acg"), b"cgt".to_vec()); // a->t,c->g,g->c reversed
            assert_eq!(revcomp_keep_case(b""), Vec::<u8>::new());
            // unmapped byte passes through unchanged (Python str.translate leaves it)
            assert_eq!(revcomp_keep_case(b"AXT"), b"AXT".to_vec()); // T->A, X->X, A->T ; reversed = A X T
        }

        #[test]
        fn hw_distance_corners() {
            assert_eq!(hw_distance(b"ACGT", b"TTACGTGG"), 0); // exact infix
            assert_eq!(hw_distance(b"ACGT", b"ACGT"), 0);
            assert_eq!(hw_distance(b"ACGTACGTACGT", b"ACGT"), 8); // q longer -> 8 deletions
            assert_eq!(hw_distance(b"ACGT", b"ACAT"), 1); // one substitution
            assert_eq!(hw_distance(b"A", b"TTTT"), 1); // best infix "" or one mismatch
        }

        /// PROVES the pure-Rust HW-DP `aln_id` == edlib on >= 80 adversarial pairs (exact
        /// bit-for-bit float equality). Reports the first mismatch. The fixture is the `aln_id` section of
        /// the retired `bridge_detector_fixture.json`.
        #[test]
        fn aln_id_parity_vs_edlib() {
            let fx: serde_json::Value =
                serde_json::from_str(include_str!("testdata/aln_id_fixture.json"))
                    .expect("parse aln_id fixture json");
            let cases = fx["aln_id"].as_array().expect("aln_id array");
            assert!(
                cases.len() >= 80,
                "need >= 80 aln_id cases, got {}",
                cases.len()
            );
            let mut classes: BTreeSet<String> = BTreeSet::new();
            let mut n_empty = 0;
            let mut n_n = 0;
            let mut n_revcomp = 0;
            for c in cases {
                let q = c["q"].as_str().unwrap();
                let t = c["t"].as_str().unwrap();
                // EXACT bits (serde_json's decimal float parser is 1-ULP imprecise).
                let want = f64::from_bits(c["bits"].as_u64().unwrap());
                let cls = c["cls"].as_str().unwrap_or("");
                classes.insert(cls.to_string());
                if q.is_empty() || t.is_empty() {
                    n_empty += 1;
                }
                if q.contains('N') || t.contains('N') {
                    n_n += 1;
                }
                if cls.contains("revcomp") {
                    n_revcomp += 1;
                }
                let got = aln_id(q.as_bytes(), t.as_bytes());
                assert_eq!(
                    got.to_bits(),
                    want.to_bits(),
                    "aln_id MISMATCH cls='{cls}' q='{q}' t='{t}': rust {got} (bits {:#x}) != python {want} (bits {:#x})",
                    got.to_bits(),
                    want.to_bits()
                );
            }
            // edge-class coverage guard
            for required in [
                "both_empty",
                "q_empty",
                "t_empty",
                "q_longer",
                "exact_infix",
                "tandem",
            ] {
                assert!(
                    classes.contains(required),
                    "aln_id fixture missing class '{required}'"
                );
            }
            assert!(
                n_empty >= 3 && n_n >= 3 && n_revcomp >= 3,
                "aln_id fixture lacks edge coverage"
            );
        }
    }
}

pub mod parcn {
    //! Assembly-based paralog-specific copy number (parCN). OPTIONAL assembly/DNA-side supplement:
    //! projects catalog copy consensuses onto phased haplotype assemblies and counts per-copy genomic loci,
    //! disambiguated by deterministic SUN witnesses. Consumes only copies.fa + the assemblies; never wires
    //! into the RNA-exclusive core. See docs/superpowers/specs/2026-07-14-assembly-parcn-design.md.
    //!
    //! **STATUS:** OTHER-BINARY — binary `parcn` (Cargo.toml:116-118, path src/bin/parcn.rs); no flag inside it gates the module — the whole binary is the opt-in  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    use std::collections::BTreeMap;

    use crate::family::copy_assign::copy_assign_pipeline::banded_msa_pair;

    #[derive(Clone, Debug)]
    pub struct Copy {
        pub family_id: String,
        pub copy_id: String,
        pub seq: Vec<u8>,
    }

    #[derive(Clone, Debug, PartialEq)]
    pub enum Tier {
        T1,
        T2,
        T3,
        NA,
    }

    #[derive(Clone, Debug)]
    pub struct CopySun {
        pub copy_id: String,
        pub tier: Tier,
        pub private: Vec<(usize, u8)>,
    }

    /// For copy `b` vs sibling `s`: the SET of offsets in `b` (non-gap) whose aligned base in `s` differs
    /// (substitution or gap), plus `(matches, aligned_cols)` for identity. Uses the banded 2-row MSA; if the
    /// pair can't be aligned in-band, a difference from `s` cannot be CONFIRMED for any offset, so this
    /// returns an EMPTY diff set. Callers intersect diff sets across siblings to build a copy's private-position
    /// set (`private = ∩ diff_offsets(b, s_j)`); an empty set is absorbing under intersection (`∅ ∩ X = ∅`), so a
    /// failed comparison correctly removes ALL candidates from privateness (conservative: nothing private).
    /// Returning "every offset differs" instead would be a no-op under intersection (`U ∩ X = X`) and could
    /// fabricate a spurious Tier-1 (SUN) call from zero evidence — the bug this fallback fixes.
    fn diff_offsets(
        b: &[u8],
        s: &[u8],
        band: usize,
    ) -> (std::collections::HashSet<usize>, usize, usize) {
        let msa = match banded_msa_pair(b, s, band) {
            Some(m) => m,
            None => return (std::collections::HashSet::new(), 0, b.len().max(1)),
        };
        let (ab, asb) = (&msa[0], &msa[1]);
        let (mut boff, mut diff, mut matches, mut cols) =
            (0usize, std::collections::HashSet::new(), 0usize, 0usize);
        for k in 0..ab.len() {
            let (cb, cs) = (ab[k], asb[k]);
            if cb != b'-' {
                cols += 1;
                if cs == cb {
                    matches += 1;
                } else {
                    diff.insert(boff);
                }
                boff += 1;
            }
        }
        (diff, matches, cols.max(1))
    }

    /// Per-copy private positions (a SUN = an offset in copy B whose base differs from EVERY sibling) + tier.
    /// Band scales with the family's copy-length spread. Threshold-free: T1 iff ≥1 private position; T3 iff a
    /// sibling is ≥99.9% identical (indistinguishable); T2 otherwise; NA for a single-copy family.
    pub fn sun_positions(copies: &[Copy], band: usize) -> Vec<CopySun> {
        let mut out = Vec::with_capacity(copies.len());
        for (i, b) in copies.iter().enumerate() {
            if copies.len() == 1 {
                out.push(CopySun {
                    copy_id: b.copy_id.clone(),
                    tier: Tier::NA,
                    private: Vec::new(),
                });
                continue;
            }
            // private = offsets differing from ALL siblings; max_id = closest sibling identity.
            let mut private_set: Option<std::collections::HashSet<usize>> = None;
            let mut max_id = 0.0f64;
            for (j, s) in copies.iter().enumerate() {
                if i == j {
                    continue;
                }
                let (diff, matches, cols) = diff_offsets(&b.seq, &s.seq, band);
                max_id = max_id.max(matches as f64 / cols as f64);
                private_set = Some(match private_set {
                    None => diff,
                    Some(acc) => acc.intersection(&diff).copied().collect(),
                });
            }
            let private_set = private_set.unwrap_or_default();
            let mut private: Vec<(usize, u8)> =
                private_set.iter().map(|&p| (p, b.seq[p])).collect();
            private.sort_unstable();
            let tier = if !private.is_empty() {
                Tier::T1
            } else if max_id >= 0.999 {
                Tier::T3
            } else {
                Tier::T2
            };
            out.push(CopySun {
                copy_id: b.copy_id.clone(),
                tier,
                private,
            });
        }
        out
    }

    /// Parse a `gw_family_catalog` `copies.fa` (`>{family}|{copy_idx}|{chrom}:{s}-{e}|{strand}|nexon={n}`) into
    /// families → copies, preserving file order. Sequence lines are concatenated and upper-cased.
    pub fn parse_copies_fa(path: &str) -> anyhow::Result<BTreeMap<String, Vec<Copy>>> {
        let text = std::fs::read_to_string(path)?;
        let mut fams: BTreeMap<String, Vec<Copy>> = BTreeMap::new();
        let (mut fam, mut cid, mut seq) = (String::new(), String::new(), Vec::<u8>::new());
        let mut have = false;
        let flush =
            |fams: &mut BTreeMap<String, Vec<Copy>>, fam: &str, cid: &str, seq: &mut Vec<u8>| {
                if !fam.is_empty() {
                    fams.entry(fam.to_string()).or_default().push(Copy {
                        family_id: fam.to_string(),
                        copy_id: cid.to_string(),
                        seq: std::mem::take(seq),
                    });
                }
            };
        for line in text.lines() {
            if let Some(h) = line.strip_prefix('>') {
                if have {
                    flush(&mut fams, &fam, &cid, &mut seq);
                }
                let mut it = h.split('|');
                fam = it.next().unwrap_or("").to_string();
                cid = it.next().unwrap_or("0").to_string();
                have = true;
            } else {
                seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
            }
        }
        if have {
            flush(&mut fams, &fam, &cid, &mut seq);
        }
        Ok(fams)
    }

    /// Walk a minimap2 `cs:Z:` short string and return, per requested QUERY offset, the aligned TARGET
    /// (assembly) base. cs ops: `:N` = N matches (target base == query base); `*xy` = substitution, x = target
    /// base, y = query base (advance query 1); `+seq` = insertion, query-only (target base None); `-seq` =
    /// deletion, target-only (no query advance); `~` splice = target-only. Bases are returned upper-cased.
    pub fn cs_bases_at(cs: &str, query: &[u8], positions: &[usize]) -> Vec<Option<u8>> {
        let want: std::collections::HashSet<usize> = positions.iter().copied().collect();
        let mut base_at: std::collections::HashMap<usize, Option<u8>> =
            std::collections::HashMap::new();
        let bytes = cs.as_bytes();
        let (mut k, mut qoff) = (0usize, 0usize);
        while k < bytes.len() {
            match bytes[k] {
                b':' => {
                    let mut j = k + 1;
                    let mut n = 0usize;
                    while j < bytes.len() && bytes[j].is_ascii_digit() {
                        n = n * 10 + (bytes[j] - b'0') as usize;
                        j += 1;
                    }
                    for _ in 0..n {
                        if want.contains(&qoff) {
                            base_at.insert(qoff, query.get(qoff).map(|b| b.to_ascii_uppercase()));
                        }
                        qoff += 1;
                    }
                    k = j;
                }
                b'*' => {
                    // *<target><query>
                    let tgt = bytes.get(k + 1).map(|b| b.to_ascii_uppercase());
                    if want.contains(&qoff) {
                        base_at.insert(qoff, tgt);
                    }
                    qoff += 1;
                    k += 3;
                }
                b'+' => {
                    let mut j = k + 1;
                    while j < bytes.len() && bytes[j].is_ascii_alphabetic() {
                        if want.contains(&qoff) {
                            base_at.insert(qoff, None);
                        }
                        qoff += 1;
                        j += 1;
                    }
                    k = j;
                }
                b'-' | b'~' => {
                    // target-only: skip the following letters/coords, no query advance
                    let mut j = k + 1;
                    while j < bytes.len()
                        && bytes[j] != b':'
                        && bytes[j] != b'*'
                        && bytes[j] != b'+'
                        && bytes[j] != b'-'
                        && bytes[j] != b'~'
                    {
                        j += 1;
                    }
                    k = j;
                }
                _ => {
                    k += 1;
                }
            }
        }
        positions
            .iter()
            .map(|p| base_at.get(p).copied().flatten())
            .collect()
    }

    #[derive(Clone, Debug, PartialEq)]
    pub enum Method {
        Sun,
        AlignFallback,
        Unresolved,
        SingleCopy,
    }

    #[derive(Clone, Debug)]
    pub struct Locus {
        pub chrom: String,
        pub start: u64,
        pub end: u64,
        pub best_copy: String,
        pub identity: f64,
        pub runner_up_identity: f64,
        pub cs: String,
        /// PAF query-start/end (forward-query coordinates of the aligned segment `[qs,qe)`) and strand of this
        /// hit, needed to map a copy's private forward-query offset to the right column in `cs` (see
        /// `cs_match_states`): the cs tag only covers `[qs,qe)`, walked in TARGET-forward order (i.e.
        /// query-reversed when `strand=='-'`).
        pub qs: u64,
        pub qe: u64,
        pub strand: char,
    }

    #[derive(Clone, Debug)]
    pub struct Assignment {
        pub copy_id: Option<String>,
        pub method: Method,
    }

    const ALIGN_MARGIN: f64 = 0.002;

    /// Whether the projection alignment carries a MATCH at each forward-query offset. `Some(true)` = the aligned
    /// assembly base equals the copy's base there (so the copy's private allele IS present); `Some(false)` =
    /// substitution; `None` = indel or outside the aligned segment `[qs,qe)`. Strand-symmetric: a cs `:` match
    /// means identical aligned bases regardless of orientation, so only the POSITION mapping uses strand
    /// (`+`: a = p-qs walking the query forward; `-`: a = qe-1-p, since cs walks the target forward = query
    /// reversed). This is what the SUN gate needs: the private base is confirmed iff its column is a match.
    pub fn cs_match_states(
        cs: &str,
        qs: u64,
        qe: u64,
        strand: char,
        positions: &[usize],
    ) -> Vec<Option<bool>> {
        // walk cs once: aligned-segment query-offset `a` -> is_match
        let mut st: std::collections::HashMap<u64, bool> = std::collections::HashMap::new();
        let bytes = cs.as_bytes();
        let (mut k, mut a) = (0usize, 0u64);
        while k < bytes.len() {
            match bytes[k] {
                b':' => {
                    let mut j = k + 1;
                    let mut n = 0u64;
                    while j < bytes.len() && bytes[j].is_ascii_digit() {
                        n = n * 10 + (bytes[j] - b'0') as u64;
                        j += 1;
                    }
                    for _ in 0..n {
                        st.insert(a, true);
                        a += 1;
                    }
                    k = j;
                }
                b'*' => {
                    st.insert(a, false);
                    a += 1;
                    k += 3;
                }
                b'+' => {
                    let mut j = k + 1;
                    while j < bytes.len() && bytes[j].is_ascii_alphabetic() {
                        a += 1;
                        j += 1;
                    }
                    k = j;
                } // insertion: query-only -> leave unset (None)
                b'-' | b'~' => {
                    let mut j = k + 1;
                    while j < bytes.len() && !b":*+-~".contains(&bytes[j]) {
                        j += 1;
                    }
                    k = j;
                } // target-only
                _ => {
                    k += 1;
                }
            }
        }
        positions
            .iter()
            .map(|&p| {
                let p = p as u64;
                if p < qs || p >= qe {
                    return None;
                } // outside the aligned segment -> unconfirmable
                let a = if strand == '-' { qe - 1 - p } else { p - qs };
                st.get(&a).copied()
            })
            .collect()
    }

    /// Hybrid assignment of a projected locus to its best copy. Tier-1: confirm the assembly carries the best
    /// copy's private base (a cs MATCH at that private offset, mapped through `qs`/`qe`/`strand`) at ≥1 private
    /// position → deterministic SUN. Tier-2: assign to the best copy iff its identity beats the runner-up by
    /// ≥ ALIGN_MARGIN (flagged fallback). Tier-3 / private-not-confirmed / near-tie → UNRESOLVED. NA (single
    /// copy) → single_copy.
    pub fn assign_locus(locus: &Locus, sun: &CopySun) -> Assignment {
        match sun.tier {
            Tier::NA => Assignment {
                copy_id: Some(locus.best_copy.clone()),
                method: Method::SingleCopy,
            },
            Tier::T1 => {
                let positions: Vec<usize> = sun.private.iter().map(|&(p, _)| p).collect();
                let states =
                    cs_match_states(&locus.cs, locus.qs, locus.qe, locus.strand, &positions);
                let confirmed = states.iter().any(|s| *s == Some(true)); // assembly carries the copy's private base
                if confirmed {
                    Assignment {
                        copy_id: Some(locus.best_copy.clone()),
                        method: Method::Sun,
                    }
                } else {
                    Assignment {
                        copy_id: None,
                        method: Method::Unresolved,
                    }
                }
            }
            Tier::T2 => {
                if locus.identity - locus.runner_up_identity >= ALIGN_MARGIN {
                    Assignment {
                        copy_id: Some(locus.best_copy.clone()),
                        method: Method::AlignFallback,
                    }
                } else {
                    Assignment {
                        copy_id: None,
                        method: Method::Unresolved,
                    }
                }
            }
            Tier::T3 => Assignment {
                copy_id: None,
                method: Method::Unresolved,
            },
        }
    }

    #[derive(Clone, Debug)]
    pub struct ParcnRow {
        pub family_id: String,
        pub copy_id: String,
        pub tier: Tier,
        pub loci_mat: usize,
        pub loci_pat: usize,
        pub method: Method,
    }

    fn tier_str(t: &Tier) -> &'static str {
        match t {
            Tier::T1 => "T1",
            Tier::T2 => "T2",
            Tier::T3 => "T3",
            Tier::NA => "NA",
        }
    }
    fn method_str(m: &Method) -> &'static str {
        match m {
            Method::Sun => "SUN",
            Method::AlignFallback => "align_fallback",
            Method::Unresolved => "UNRESOLVED",
            Method::SingleCopy => "single_copy",
        }
    }

    /// Count per-copy assigned loci (mat/pat), pick each copy's dominant assignment method, and total the
    /// unresolved loci across both haplotypes. Copies with no assigned locus still get a row (parCN 0).
    pub fn tabulate(
        family_id: &str,
        copies: &[Copy],
        suns: &[CopySun],
        mat: &[Assignment],
        pat: &[Assignment],
    ) -> (Vec<ParcnRow>, usize) {
        use std::collections::HashMap;
        let tier_of: HashMap<&str, &Tier> =
            suns.iter().map(|s| (s.copy_id.as_str(), &s.tier)).collect();
        let mut mat_c: HashMap<String, usize> = HashMap::new();
        let mut pat_c: HashMap<String, usize> = HashMap::new();
        let mut method_of: HashMap<String, Method> = HashMap::new();
        let mut n_unres = 0usize;
        for (side, counts) in [(mat, &mut mat_c), (pat, &mut pat_c)] {
            for a in side {
                match &a.copy_id {
                    Some(cp) => {
                        *counts.entry(cp.clone()).or_insert(0) += 1;
                        method_of
                            .entry(cp.clone())
                            .or_insert_with(|| a.method.clone());
                    }
                    None => n_unres += 1,
                }
            }
        }
        let mut rows = Vec::with_capacity(copies.len());
        for c in copies {
            let tier = (*tier_of.get(c.copy_id.as_str()).unwrap_or(&&Tier::NA)).clone();
            let method = method_of
                .get(&c.copy_id)
                .cloned()
                .unwrap_or(Method::Unresolved);
            rows.push(ParcnRow {
                family_id: family_id.to_string(),
                copy_id: c.copy_id.clone(),
                tier,
                loci_mat: *mat_c.get(&c.copy_id).unwrap_or(&0),
                loci_pat: *pat_c.get(&c.copy_id).unwrap_or(&0),
                method,
            });
        }
        (rows, n_unres)
    }

    pub fn format_parcn_row(r: &ParcnRow) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.family_id,
            r.copy_id,
            tier_str(&r.tier),
            r.loci_mat,
            r.loci_pat,
            r.loci_mat + r.loci_pat,
            method_str(&r.method)
        )
    }
    pub fn format_family_row(family_id: &str, rows: &[ParcnRow], n_unresolved: usize) -> String {
        let famcn: usize = rows.iter().map(|r| r.loci_mat + r.loci_pat).sum();
        format!("{}\t{}\t{}\t{}", family_id, rows.len(), famcn, n_unresolved)
    }

    fn recip_overlap(a: &Locus, b: &Locus) -> f64 {
        if a.chrom != b.chrom {
            return 0.0;
        }
        let (lo, hi) = (a.start.max(b.start), a.end.min(b.end));
        if hi <= lo {
            return 0.0;
        }
        let ov = (hi - lo) as f64;
        let la = (a.end - a.start).max(1) as f64;
        let lb = (b.end - b.start).max(1) as f64;
        (ov / la).min(ov / lb)
    }

    /// Collapse reciprocal-overlap ≥ 0.50 loci into one, keeping the highest-identity member (its best_copy)
    /// and recording the next-highest overlapping identity as runner_up_identity (for the Tier-2 margin gate).
    pub fn dedup_loci(mut loci: Vec<Locus>) -> Vec<Locus> {
        loci.sort_by(|a, b| {
            b.identity
                .partial_cmp(&a.identity)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut kept: Vec<Locus> = Vec::new();
        for l in loci {
            if let Some(k) = kept.iter_mut().find(|k| recip_overlap(k, &l) >= 0.50) {
                // loci are sorted DESC by identity, so the kept member `k` already has the higher identity and
                // `l` is a runner-up candidate for that overlap group; record the highest runner-up seen.
                k.runner_up_identity = k.runner_up_identity.max(l.identity);
            } else {
                kept.push(l);
            }
        }
        kept
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn parse_copies_fa_groups_by_family() {
            let dir = std::env::temp_dir();
            let p = dir.join(format!("parcn_copies_{}.fa", std::process::id()));
            std::fs::write(&p, ">RBMY|0|chrY:1-9|+|nexon=3\nACGTACGTA\n>RBMY|1|chrY:20-28|+|nexon=3\nACGTACGTT\n>DAZ|0|chrY:99-104|-|nexon=1\nGGGCCC\n").unwrap();
            let fams = parse_copies_fa(p.to_str().unwrap()).unwrap();
            std::fs::remove_file(&p).ok();
            assert_eq!(fams.len(), 2);
            assert_eq!(fams["RBMY"].len(), 2);
            assert_eq!(fams["RBMY"][0].copy_id, "0");
            assert_eq!(fams["RBMY"][1].seq, b"ACGTACGTT");
            assert_eq!(fams["DAZ"][0].seq, b"GGGCCC");
        }

        #[test]
        fn sun_positions_finds_private_snv_and_tiers() {
            // copy0 vs copy1 differ ONLY at offset 4 (A vs T) -> each has a private position there (Tier-1).
            // copy2 is identical to copy0 -> Tier-3 (indistinguishable), and offset 4 is no longer private to copy0.
            let copies = vec![
                Copy {
                    family_id: "F".into(),
                    copy_id: "0".into(),
                    seq: b"ACGTAGGTCA".to_vec(),
                },
                Copy {
                    family_id: "F".into(),
                    copy_id: "1".into(),
                    seq: b"ACGTTGGTCA".to_vec(),
                },
                Copy {
                    family_id: "F".into(),
                    copy_id: "2".into(),
                    seq: b"ACGTAGGTCA".to_vec(),
                },
            ];
            let suns = sun_positions(&copies, 8);
            let s0 = suns.iter().find(|s| s.copy_id == "0").unwrap();
            let s1 = suns.iter().find(|s| s.copy_id == "1").unwrap();
            let s2 = suns.iter().find(|s| s.copy_id == "2").unwrap();
            // copy1's 'T' at offset 4 is unique among the three -> private -> Tier-1.
            assert_eq!(s1.tier, Tier::T1);
            assert!(s1.private.iter().any(|&(p, b)| p == 4 && b == b'T'));
            // copy0 and copy2 are identical -> neither can have a private position -> Tier-3.
            assert_eq!(s0.tier, Tier::T3);
            assert_eq!(s2.tier, Tier::T3);
            assert!(s0.private.is_empty());
        }

        #[test]
        fn sun_positions_single_copy_is_na() {
            let copies = vec![Copy {
                family_id: "F".into(),
                copy_id: "0".into(),
                seq: b"ACGTACGT".to_vec(),
            }];
            let suns = sun_positions(&copies, 8);
            assert_eq!(suns[0].tier, Tier::NA);
        }

        #[test]
        fn cs_bases_reads_match_and_substitution_and_insertion() {
            // query ACGTACGT (len 8). cs: 3 matches, sub (target g / query t) at q=3, 2 matches,
            // insertion of "AA" at q=6..8, tail is target-only del (does not advance query).
            // cs grammar: :N match run; *<tgt><qry> substitution; +<seq> insertion (query-only); -<seq> deletion (target-only).
            let cs = ":3*gt:2+aa-cc";
            let q = b"ACGTACGT";
            // q0 match -> assembly base 'A'(=query); q3 substitution -> assembly base 'G'(target, upper); q5 match 'C'; q6 insertion -> None.
            let got = cs_bases_at(cs, q, &[0, 3, 5, 6]);
            assert_eq!(got, vec![Some(b'A'), Some(b'G'), Some(b'C'), None]);
        }

        #[test]
        fn sun_positions_band_edge_yields_no_private() {
            // length difference (8 vs 12) exceeds a tiny band -> banded_msa_pair returns None ->
            // the conservative fallback must yield NO private positions (never a spurious Tier-1).
            let copies = vec![
                Copy {
                    family_id: "F".into(),
                    copy_id: "0".into(),
                    seq: b"ACGTACGT".to_vec(),
                },
                Copy {
                    family_id: "F".into(),
                    copy_id: "1".into(),
                    seq: b"ACGTACGTACGT".to_vec(),
                },
            ];
            let suns = sun_positions(&copies, 1); // band=1 << |8-12|
            for s in &suns {
                assert!(
                    s.private.is_empty(),
                    "band-edge failure must not fabricate private positions"
                );
                assert_ne!(
                    s.tier,
                    super::Tier::T1,
                    "a copy that could not be compared must not be Tier-1"
                );
            }
        }

        #[test]
        fn assign_locus_hybrid_tiers() {
            let mk = |cs: &str, id: f64, ru: f64| Locus {
                chrom: "c".into(),
                start: 0,
                end: 9,
                best_copy: "0".into(),
                identity: id,
                runner_up_identity: ru,
                cs: cs.into(),
                qs: 0,
                qe: 10,
                strand: '+',
            };
            let sun_t1 = CopySun {
                copy_id: "0".into(),
                tier: Tier::T1,
                private: vec![(4, b'A')],
            };
            // cs shows a MATCH across offset 4 -> assembly carries 'A' -> SUN confirmed.
            assert_eq!(
                assign_locus(&mk(":10", 0.99, 0.90), &sun_t1).method,
                Method::Sun
            );
            // cs shows a substitution at offset 4 (target g) -> assembly does NOT carry 'A' -> UNRESOLVED.
            assert_eq!(
                assign_locus(&mk(":4*ga:5", 0.99, 0.90), &sun_t1).method,
                Method::Unresolved
            );
            // Tier-2 with a clear identity margin -> align_fallback.
            let sun_t2 = CopySun {
                copy_id: "0".into(),
                tier: Tier::T2,
                private: vec![],
            };
            assert_eq!(
                assign_locus(&mk(":10", 0.99, 0.90), &sun_t2).method,
                Method::AlignFallback
            );
            // Tier-2 near-tie -> UNRESOLVED.
            assert_eq!(
                assign_locus(&mk(":10", 0.991, 0.990), &sun_t2).method,
                Method::Unresolved
            );
            // Tier-3 -> UNRESOLVED. NA -> single_copy.
            let sun_t3 = CopySun {
                copy_id: "0".into(),
                tier: Tier::T3,
                private: vec![],
            };
            assert_eq!(
                assign_locus(&mk(":10", 0.99, 0.0), &sun_t3).method,
                Method::Unresolved
            );
            let sun_na = CopySun {
                copy_id: "0".into(),
                tier: Tier::NA,
                private: vec![],
            };
            assert_eq!(
                assign_locus(&mk(":10", 0.99, 0.0), &sun_na).method,
                Method::SingleCopy
            );
        }

        #[test]
        fn tabulate_counts_and_formats() {
            let copies = vec![
                Copy {
                    family_id: "RBMY".into(),
                    copy_id: "0".into(),
                    seq: b"AAAA".to_vec(),
                },
                Copy {
                    family_id: "RBMY".into(),
                    copy_id: "1".into(),
                    seq: b"AAAT".to_vec(),
                },
            ];
            let suns = vec![
                CopySun {
                    copy_id: "0".into(),
                    tier: Tier::T1,
                    private: vec![(3, b'A')],
                },
                CopySun {
                    copy_id: "1".into(),
                    tier: Tier::T2,
                    private: vec![],
                },
            ];
            let a = |cp: &str, m: Method| Assignment {
                copy_id: Some(cp.into()),
                method: m,
            };
            let un = || Assignment {
                copy_id: None,
                method: Method::Unresolved,
            };
            // mat: copy0 SUN once, one unresolved. pat: copy0 SUN once, copy1 fallback once.
            let mat = vec![a("0", Method::Sun), un()];
            let pat = vec![a("0", Method::Sun), a("1", Method::AlignFallback)];
            let (rows, n_unres) = tabulate("RBMY", &copies, &suns, &mat, &pat);
            let r0 = rows.iter().find(|r| r.copy_id == "0").unwrap();
            assert_eq!((r0.loci_mat, r0.loci_pat), (1, 1)); // parCN 2
            let r1 = rows.iter().find(|r| r.copy_id == "1").unwrap();
            assert_eq!((r1.loci_mat, r1.loci_pat), (0, 1)); // parCN 1
            assert_eq!(n_unres, 1);
            assert_eq!(format_parcn_row(r0), "RBMY\t0\tT1\t1\t1\t2\tSUN");
            assert_eq!(format_family_row("RBMY", &rows, n_unres), "RBMY\t2\t3\t1");
        }

        #[test]
        fn dedup_collapses_overlapping_keeps_best() {
            let mk = |s: u64, e: u64, cp: &str, id: f64| Locus {
                chrom: "c1".into(),
                start: s,
                end: e,
                best_copy: cp.into(),
                identity: id,
                runner_up_identity: 0.0,
                cs: ":1".into(),
                qs: 0,
                qe: 1,
                strand: '+',
            };
            // Two heavily-overlapping hits (copy0 id .99, copy1 id .97) + one disjoint locus.
            let loci = vec![
                mk(1000, 2000, "0", 0.99),
                mk(1010, 1990, "1", 0.97),
                mk(50000, 51000, "3", 0.98),
            ];
            let mut out = dedup_loci(loci);
            out.sort_by_key(|l| l.start);
            assert_eq!(out.len(), 2);
            assert_eq!(out[0].best_copy, "0"); // highest identity wins
            assert!((out[0].runner_up_identity - 0.97).abs() < 1e-9); // runner-up recorded
            assert_eq!(out[1].best_copy, "3");
        }

        #[test]
        fn sun_confirm_respects_qs_offset() {
            // aligned segment starts at qs=5; the copy's private base is at forward-query offset 7.
            // cs (segment-local): match up to seg-offset 1, then a substitution at seg-offset 2 (= query offset 7).
            // So the assembly does NOT carry the private base at offset 7 -> must be UNRESOLVED, not SUN.
            let sun = CopySun {
                copy_id: "0".into(),
                tier: Tier::T1,
                private: vec![(7, b'A')],
            };
            let locus = Locus {
                chrom: "c".into(),
                start: 0,
                end: 9,
                best_copy: "0".into(),
                identity: 0.99,
                runner_up_identity: 0.0,
                qs: 5,
                qe: 20,
                strand: '+',
                cs: ":2*gc:10".into(),
            }; // seg-offset 2 (query 7) is a substitution
            assert_eq!(assign_locus(&locus, &sun).method, Method::Unresolved);
            // sanity: if that same column were a match, it confirms.
            let locus_ok = Locus {
                cs: ":15".into(),
                ..locus.clone()
            };
            assert_eq!(assign_locus(&locus_ok, &sun).method, Method::Sun);
        }

        #[test]
        fn sun_confirm_respects_minus_strand() {
            // minus strand: aligned segment [qs=0, qe=10); private position at forward-query offset 8 maps to
            // segment-offset qe-1-8 = 1. cs has a substitution at segment-offset 1 -> assembly differs -> UNRESOLVED.
            let sun = CopySun {
                copy_id: "0".into(),
                tier: Tier::T1,
                private: vec![(8, b'A')],
            };
            let locus = Locus {
                chrom: "c".into(),
                start: 0,
                end: 10,
                best_copy: "0".into(),
                identity: 0.99,
                runner_up_identity: 0.0,
                qs: 0,
                qe: 10,
                strand: '-',
                cs: ":1*gc:8".into(),
            };
            assert_eq!(assign_locus(&locus, &sun).method, Method::Unresolved);
            // if segment-offset 1 is a match instead -> confirmed.
            let locus_ok = Locus {
                cs: ":10".into(),
                ..locus.clone()
            };
            assert_eq!(assign_locus(&locus_ok, &sun).method, Method::Sun);
        }
    }
}

pub mod run_cache {
    //! On-disk cache of the catalog builder's expensive intermediates, for fast re-runs and for analysis.
    //!
    //! **STATUS:** OPT-IN  (docs/MODULE_STATUS.md; `RUSTLE_CACHE_DIR`, set by default by `tools/rustle_pipeline.sh`)
    //!
    //! Enabled by `RUSTLE_CACHE_DIR=<dir>` (the pipeline driver sets `PREFIX.cache`); unset = nothing is read or
    //! written and every output is byte-identical to a build without this module.
    //!
    //! These objects are cached, each at a boundary where everything downstream reads only the cached object:
    //!
    //! * **`reps/<key>/`** — the collapsed locus representatives with their read statistics, i.e. the state of
    //!   `detect_homology_catalog_genome_wide` after both BAM passes and the locus collapse (human chr16: ~260 of
    //!   ~370 s). Files: `reps.tsv` (one row per representative, in representative-index order), `reps.fa` (the
    //!   exact sequence bytes, one line each), `key.tsv` (the full key material), `DONE`.
    //!   Under `gw_family_catalog --piecewise` the same kind holds one entry PER CONTIG (key line
    //!   `rustle catalog reps v1 contig=<name>`) and one for their merge (`rustle catalog reps v1 piecewise-merge`),
    //!   same files; see `denovo_pipeline::detect_homology_catalog_piecewise`. A contig split by `--piece-records` has
    //!   one entry per read-free sub-range (`contig=<name>:<lo>-<hi>`) and its cut plan in **`plan/<key>/pieces.tsv`**.
    //! * **`paf/<key>/`** — one all-vs-all minimap2 PAF per (query bytes, command line, minimap2 version):
    //!   `out.paf`, `key.tsv`, `DONE`.
    //!   The families stage (`mcl_families --from-gtf`, key line `rustle families paf v2`) keys its entry on a
    //!   [`ContentHash`] of the loci FASTA taken WHILE the FASTA is written (no re-read), and replays a hit by
    //!   HARD-LINKING `out.paf` to `<out>.loci.paf` ([`Entry::replay`]; a copy where a link is impossible) instead of
    //!   copying it. Its entry is PINNED ([`Entry::pinned`]): `DONE` also records the payload's mtime (ns), inode, a
    //!   sampled-content fingerprint and its full content hash, so an in-place write to the shared inode (through the
    //!   linked product) is a miss, never a stale replay. The one theoretical stale replay left is the `.asbin` sidecar's
    //!   edge: a same-size in-place rewrite that restores the mtime to the nanosecond and changes no sampled block;
    //!   `RUSTLE_CACHE_VERIFY=1` re-hashes every pinned payload on a hit and closes it (audit mode). Every writer of a
    //!   replayed product unlinks it first, never truncates it, and a payload is linked to at most one product (a
    //!   second prefix sharing the cache directory gets a copy), so two products never share an inode.
    //!   `o3_candidates::minimap2` (key line `rustle o3 minimap2 v1`) keeps each of its minimap2 calls here the same way, keyed on the
    //!   content hashes of the target and of the query file.
    //! * **`cand/<key>/`** — the result of the `o3_candidates` stage (spec `docs/superpowers/specs/2026-10-02-o3-candidates-design.md`
    //!   §5.8): `candidates.tsv` and `contigs.fa` are required (an empty `contigs.fa`, a run that flagged nothing, is a complete payload).
    //!
    //! A hit requires `DONE` and a `key.tsv` byte-identical to the key the current run computes, so a hash
    //! collision in the directory name can only cause a miss, never a wrong hit. Writes go to a temporary
    //! directory renamed into place, so an interrupted run leaves no partial entry.
    //!
    //! The representatives key covers: the executable itself (path, size, mtime — any rebuild invalidates), the
    //! BAM and its index, the FASTA and its index (path, size, mtime), the `DenovoConfig`, and every `RUSTLE_*`
    //! variable EXCEPT [`DOWNSTREAM_ONLY_ENV`], the settings read only after the boundary (each verified by grep
    //! on 2026-09-24: the E_r edge rule, gamma, the coverage split, the edge dump and logging switches). An
    //! unknown variable therefore over-invalidates; it can never produce a stale hit.
    use crate::family::family_detect::DenovoTranscript;
    use anyhow::{Context, Result};
    use std::io::{BufRead, Write};
    use std::path::{Path, PathBuf};

    /// `RUSTLE_*` settings read only downstream of the representatives boundary (or output-neutral there). Diagnostics
    /// that print UPSTREAM (`RUSTLE_COLLAPSE_STATS`, `RUSTLE_LOCUS_AUDIT`, `RUSTLE_DEBUG_LOCUS`) are deliberately NOT
    /// listed: setting one changes the key, so the representatives are recomputed and the diagnostics print. Varying them
    /// re-uses the cached representatives — which is the point: edge-rule sweeps then skip both BAM passes and
    /// the collapse. Anything else in the environment is part of the key.
    pub const DOWNSTREAM_ONLY_ENV: &[&str] = &[
        "RUSTLE_CACHE_DIR",
        "RUSTLE_ER_MIN_COVERAGE",
        "RUSTLE_ER_COVERAGE_LONGER_FLOOR",
        "RUSTLE_ER_SENSITIVE_ONLY",
        "RUSTLE_ER_EDGE_DUMP",
        "RUSTLE_GENOME_GAMMA",
        "RUSTLE_COVERAGE_SPLIT",
        "RUSTLE_POA_MEMO",
        "RUSTLE_CACHE_VERIFY",
    ];

    /// The cache root, or `None` when caching is off.
    pub fn cache_root() -> Option<PathBuf> {
        std::env::var("RUSTLE_CACHE_DIR")
            .ok()
            .filter(|v| !v.is_empty())
            .map(PathBuf::from)
    }

    /// FNV-1a 64 — stable across Rust releases (unlike `DefaultHasher`), used only to NAME entries; the key text
    /// itself is compared in full on every hit.
    pub fn fnv1a64(bytes: &[u8]) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for &b in bytes {
            h ^= b as u64;
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
        h
    }

    /// Incremental FNV-1a 64 (same values as [`fnv1a64`] over the concatenated input).
    #[derive(Clone, Copy)]
    pub struct Fnv(u64);
    impl Default for Fnv {
        fn default() -> Self {
            Fnv(0xcbf2_9ce4_8422_2325)
        }
    }
    impl Fnv {
        pub fn update(&mut self, bytes: &[u8]) {
            for &b in bytes {
                self.0 ^= b as u64;
                self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
        pub fn finish(self) -> u64 {
            self.0
        }
    }

    /// A stable, word-at-a-time 128-bit content hash: two multiply-rotate lanes over the input's little-endian 8-byte
    /// words (the last partial word zero-padded), then the byte length, then a murmur3 finaliser per lane. The value
    /// depends only on the byte stream, never on how [`ContentHash::update`] calls split it, and not on the Rust release
    /// or the platform. About 8x faster than byte-wise [`Fnv`]; used where a cache key must cover every byte of a large
    /// file this process writes itself (hashed as it is written, so a hit never re-reads it).
    #[derive(Clone, Debug)]
    pub struct ContentHash {
        a: u64,
        b: u64,
        len: u64,
        tail: [u8; 8],
        ntail: usize,
    }
    impl Default for ContentHash {
        fn default() -> Self {
            ContentHash {
                a: 0x243f_6a88_85a3_08d3,
                b: 0x1319_8a2e_0370_7344,
                len: 0,
                tail: [0; 8],
                ntail: 0,
            }
        }
    }
    impl ContentHash {
        #[inline]
        fn word(&mut self, w: u64) {
            self.a = (self.a ^ w)
                .wrapping_mul(0x9e37_79b9_7f4a_7c15)
                .rotate_left(31);
            self.b = (self.b ^ w.rotate_left(23))
                .wrapping_mul(0xc2b2_ae3d_27d4_eb4f)
                .rotate_left(27);
        }
        pub fn update(&mut self, mut bytes: &[u8]) {
            self.len += bytes.len() as u64;
            if self.ntail > 0 {
                let take = (8 - self.ntail).min(bytes.len());
                self.tail[self.ntail..self.ntail + take].copy_from_slice(&bytes[..take]);
                self.ntail += take;
                bytes = &bytes[take..];
                if self.ntail < 8 {
                    return;
                }
                let w = u64::from_le_bytes(self.tail);
                self.word(w);
                self.ntail = 0;
            }
            let mut chunks = bytes.chunks_exact(8);
            for c in &mut chunks {
                self.word(u64::from_le_bytes([
                    c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7],
                ]));
            }
            let rem = chunks.remainder();
            self.tail[..rem.len()].copy_from_slice(rem);
            self.ntail = rem.len();
        }
        /// Bytes hashed so far.
        pub fn len(&self) -> u64 {
            self.len
        }
        pub fn is_empty(&self) -> bool {
            self.len == 0
        }
        /// 32 hex digits (the two lanes).
        pub fn hex(&self) -> String {
            fn fmix(mut k: u64) -> u64 {
                k ^= k >> 33;
                k = k.wrapping_mul(0xff51_afd7_ed55_8ccd);
                k ^= k >> 33;
                k = k.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
                k ^ (k >> 33)
            }
            let mut s = self.clone();
            if s.ntail > 0 {
                let mut t = [0u8; 8];
                t[..s.ntail].copy_from_slice(&s.tail[..s.ntail]);
                s.word(u64::from_le_bytes(t));
            }
            s.word(s.len);
            format!("{:016x}{:016x}", fmix(s.a), fmix(s.b ^ s.a.rotate_left(32)))
        }
        /// The hash of a whole file, streamed in 4 MiB pieces.
        pub fn of_file(path: &Path) -> std::io::Result<ContentHash> {
            use std::io::Read;
            let mut f = std::fs::File::open(path)?;
            let mut h = ContentHash::default();
            let mut buf = vec![0u8; 4 << 20];
            loop {
                let n = f.read(&mut buf)?;
                if n == 0 {
                    return Ok(h);
                }
                h.update(&buf[..n]);
            }
        }
    }

    /// A writer that hashes exactly the bytes it passes on (a [`ContentHash`] of the file as written), or passes them
    /// on untouched when built with `enabled = false` (`hash` is then `None`).
    pub struct HashingWriter<W: Write> {
        pub inner: W,
        pub hash: Option<ContentHash>,
    }
    impl<W: Write> HashingWriter<W> {
        pub fn new(inner: W, enabled: bool) -> Self {
            HashingWriter {
                inner,
                hash: enabled.then(ContentHash::default),
            }
        }
    }
    impl<W: Write> Write for HashingWriter<W> {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            let n = self.inner.write(buf)?;
            if let Some(h) = self.hash.as_mut() {
                h.update(&buf[..n]);
            }
            Ok(n)
        }
        fn flush(&mut self) -> std::io::Result<()> {
            self.inner.flush()
        }
    }

    /// The pins of a file, `mtime_ns<TAB>inode<TAB>sample`: what an in-place write through a hard link changes (a link
    /// does not). `sample` is an FNV-1a of 16 evenly spaced 4 KiB blocks (the last included) with their offsets, the
    /// `.asbin` sidecar's sampling (`denovo_assemble::AsTsvIdentity`), so even a same-size rewrite with the mtime put
    /// back is caught when it touches a sampled block.
    fn file_pins(path: &Path) -> Option<String> {
        use std::io::{Read, Seek, SeekFrom};
        let m = std::fs::metadata(path).ok()?;
        #[cfg(unix)]
        let (t, ino) = {
            use std::os::unix::fs::MetadataExt;
            (
                m.mtime() as i128 * 1_000_000_000 + m.mtime_nsec() as i128,
                m.ino(),
            )
        };
        #[cfg(not(unix))]
        let (t, ino) = (
            m.modified()
                .ok()?
                .duration_since(std::time::UNIX_EPOCH)
                .ok()?
                .as_nanos() as i128,
            0u64,
        );
        let size = m.len();
        let mut f = std::fs::File::open(path).ok()?;
        let mut fnv = Fnv::default();
        let block = 4096u64;
        let mut buf: Vec<u8> = Vec::with_capacity(block as usize);
        for i in 0..16u64 {
            let at = if size <= block {
                0
            } else {
                (size - block) * i / 15
            };
            f.seek(SeekFrom::Start(at)).ok()?;
            buf.clear();
            (&mut f).take(block).read_to_end(&mut buf).ok()?;
            fnv.update(&at.to_le_bytes());
            fnv.update(&buf);
        }
        Some(format!("{t}\t{ino}\t{:016x}", fnv.finish()))
    }

    /// `RUSTLE_CACHE_VERIFY=1`: re-hash every pinned payload on a hit (audit mode; see the module doc).
    pub fn verify_mode() -> bool {
        std::env::var("RUSTLE_CACHE_VERIFY").map_or(false, |v| v == "1")
    }

    /// Link `src` to `dest` (a new name for the same inode) after unlinking `dest`, or copy when the two are on
    /// different file systems, links are unsupported, or (`sole`) `src` already has another name besides its cache
    /// entry. `dest` is never truncated in place, so an older inode that `dest` named (a cache payload it was linked to)
    /// is left intact. Returns whether a link was made.
    ///
    /// `sole` keeps a payload linked to at most ONE product: when several output prefixes share one cache directory,
    /// the first replay links and every other prefix gets a copy, so an in-place write to one product can change the
    /// cache entry (which its pins then turn into a miss) but never another prefix's product. The driver gives each
    /// prefix its own `PREFIX.cache`, so a re-run there always links (it unlinks its own old product first).
    pub fn link_or_copy(src: &Path, dest: &Path, sole: bool) -> std::io::Result<bool> {
        match std::fs::remove_file(dest) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e),
        }
        #[cfg(unix)]
        let unshared = {
            use std::os::unix::fs::MetadataExt;
            !sole || std::fs::metadata(src)?.nlink() == 1
        };
        #[cfg(not(unix))]
        let unshared = !sole;
        if unshared && std::fs::hard_link(src, dest).is_ok() {
            return Ok(true);
        }
        std::fs::copy(src, dest)?;
        Ok(false)
    }

    /// `path<TAB>size<TAB>mtime_ns` of a file (canonical path), or `path<TAB>absent`.
    pub fn file_fingerprint(path: &str) -> String {
        let canon = std::fs::canonicalize(path)
            .map(|p| p.display().to_string())
            .unwrap_or_else(|_| path.to_string());
        match std::fs::metadata(path) {
            Ok(m) => {
                let mtime = m
                    .modified()
                    .ok()
                    .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
                    .map(|d| d.as_nanos())
                    .unwrap_or(0);
                format!("{canon}\t{}\t{mtime}", m.len())
            }
            Err(_) => format!("{canon}\tabsent"),
        }
    }

    /// The running executable's fingerprint: a rebuild changes size or mtime, so cached intermediates never
    /// survive a code change.
    pub fn exe_fingerprint() -> String {
        std::env::current_exe()
            .map(|p| file_fingerprint(&p.display().to_string()))
            .unwrap_or_else(|_| "unknown".into())
    }

    /// Every `RUSTLE_*` variable, sorted, one `k=v` per line, minus `exclude`.
    pub fn env_fingerprint(exclude: &[&str]) -> String {
        let mut v: Vec<(String, String)> = std::env::vars()
            .filter(|(k, _)| k.starts_with("RUSTLE_") && !exclude.contains(&k.as_str()))
            .collect();
        v.sort();
        v.into_iter()
            .map(|(k, val)| format!("env\t{k}={val}\n"))
            .collect()
    }

    /// One keyed cache entry: `<root>/<kind>/<fnv(key)>/`.
    pub struct Entry {
        pub dir: PathBuf,
        pub key: String,
        /// Payload files a complete entry of this kind must list in `DONE` (`reps`: reps.tsv + reps.fa, `paf`: out.paf, `cand`: candidates.tsv + contigs.fa).
        pub required: &'static [&'static str],
        /// A PINNED entry's payloads may be hard-linked out ([`Entry::replay`], [`Entry::stage_link`]): its `DONE` lines
        /// are `name<TAB>bytes<TAB>mtime_ns<TAB>inode<TAB>sample<TAB>content_hash` and a hit needs all of them to match
        /// (the full content hash only under [`verify_mode`]), so a write through a linked product invalidates the entry.
        pub pin: bool,
    }

    impl Entry {
        pub fn new(root: &Path, kind: &str, key: String) -> Entry {
            let dir = root
                .join(kind)
                .join(format!("{:016x}", fnv1a64(key.as_bytes())));
            let required: &'static [&'static str] = match kind {
                "reps" => &["reps.tsv", "reps.fa"],
                "paf" => &["out.paf"],
                "plan" => &["pieces.tsv"],
                "cand" => &["candidates.tsv", "contigs.fa"],
                _ => &[],
            };
            Entry {
                dir,
                key,
                required,
                pin: false,
            }
        }
        /// This entry, pinned (see [`Entry::pin`]).
        pub fn pinned(mut self) -> Entry {
            self.pin = true;
            self
        }
        /// A complete entry whose recorded key equals this one and whose files still have the sizes recorded at
        /// commit (`DONE` lists `name<TAB>bytes`), so a truncated file is a miss, not a silent partial replay. A pinned
        /// entry's files must also keep the mtime and inode recorded at commit ([`Entry::pin`]).
        pub fn is_hit(&self) -> bool {
            self.is_hit_verify(self.pin && verify_mode())
        }
        /// [`Entry::is_hit`], with the payload re-hash of a pinned entry forced on or off (`verify`).
        pub fn is_hit_verify(&self, verify: bool) -> bool {
            let Ok(done) = std::fs::read_to_string(self.dir.join("DONE")) else {
                return false;
            };
            // an empty or partial DONE (a crash between rename and write-back) is a miss, never a vacuous hit
            let listed: Vec<&str> = done
                .lines()
                .filter_map(|l| l.split_once('\t').map(|(n, _)| n))
                .collect();
            if listed.is_empty() || !self.required.iter().all(|r| listed.contains(r)) {
                return false;
            }
            let sizes_ok = done.lines().all(|l| {
                let c: Vec<&str> = l.split('\t').collect();
                let path = self.dir.join(c[0]);
                let size_ok = c.len() >= 2
                    && std::fs::metadata(&path)
                        .map(|m| m.len().to_string() == c[1])
                        .unwrap_or(false);
                if !size_ok {
                    return false;
                }
                if c.len() == 2 {
                    return !self.pin; // a pinned entry must carry its pins
                }
                // pinned: `name bytes mtime_ns inode sample hash` (the pins of the file as committed)
                c.len() == 6
                    && file_pins(&path).map_or(false, |p| p == c[2..5].join("\t"))
                    && (!verify || ContentHash::of_file(&path).map_or(false, |h| h.hex() == c[5]))
            });
            sizes_ok
                && std::fs::read_to_string(self.dir.join("key.tsv"))
                    .map(|k| k == self.key)
                    .unwrap_or(false)
        }
        /// Replay payload `name` of a hit as `dest` without copying it: `dest` becomes a hard link to the cached file,
        /// after unlinking whatever `dest` was; a copy where a link is impossible or the payload is already linked to
        /// another product ([`link_or_copy`] `sole`). Only for a pinned entry, whose `DONE` pins turn any later in-place
        /// write through `dest` into a miss.
        pub fn replay(&self, name: &str, dest: &Path) -> std::io::Result<bool> {
            assert!(self.pin, "only a pinned entry's payload may be linked out");
            link_or_copy(&self.dir.join(name), dest, true)
        }
        /// Put the product `src` (just written) into a staging directory as payload `name` by hard link (a copy where
        /// impossible); see [`Entry::replay`].
        pub fn stage_link(&self, staging: &Path, name: &str, src: &Path) -> std::io::Result<bool> {
            assert!(self.pin, "only a pinned entry's payload may be linked in");
            link_or_copy(src, &staging.join(name), true)
        }
        /// A fresh temporary directory next to the entry; [`Entry::commit`] renames it into place.
        pub fn staging(&self) -> Result<PathBuf> {
            static SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
            let n = SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let tmp = self
                .dir
                .with_extension(format!("tmp{}_{n}", std::process::id()));
            let _ = std::fs::remove_dir_all(&tmp);
            std::fs::create_dir_all(&tmp).with_context(|| format!("creating {}", tmp.display()))?;
            std::fs::write(tmp.join("key.tsv"), &self.key)?;
            Ok(tmp)
        }
        /// Durably publish a staging directory: every payload is fsynced, `DONE` (name + size of each file) is
        /// written and fsynced, the directory is renamed into place and the parent fsynced. If another run already
        /// published a complete entry for the same key, the staging copy is discarded instead.
        pub fn commit(&self, staging: &Path) -> Result<()> {
            if self.is_hit() {
                let _ = std::fs::remove_dir_all(staging);
                return Ok(());
            }
            let mut done = String::new();
            let mut names: Vec<String> = std::fs::read_dir(staging)?
                .filter_map(|e| e.ok())
                .map(|e| e.file_name().to_string_lossy().to_string())
                .filter(|n| n != "DONE")
                .collect();
            names.sort();
            for n in names {
                let path = staging.join(&n);
                std::fs::File::open(&path)?.sync_all()?;
                let len = std::fs::metadata(&path)?.len();
                if self.pin {
                    let pins = file_pins(&path).context("pinning a cache payload")?;
                    let h = ContentHash::of_file(&path)?.hex();
                    done.push_str(&format!("{n}\t{len}\t{pins}\t{h}\n"));
                } else {
                    done.push_str(&format!("{n}\t{len}\n"));
                }
            }
            {
                let mut f = std::fs::File::create(staging.join("DONE"))?;
                f.write_all(done.as_bytes())?;
                f.sync_all()?;
            }
            let _ = std::fs::remove_dir_all(&self.dir); // an incomplete or stale entry under this name
            if let Some(p) = self.dir.parent() {
                std::fs::create_dir_all(p)?;
            }
            std::fs::rename(staging, &self.dir)
                .with_context(|| format!("committing {}", self.dir.display()))?;
            if let Some(p) = self.dir.parent() {
                if let Ok(d) = std::fs::File::open(p) {
                    let _ = d.sync_all();
                }
            }
            Ok(())
        }
    }

    /// Write the representatives: `reps.tsv` + `reps.fa`, exact and in index order.
    pub fn write_reps(dir: &Path, reps: &[DenovoTranscript]) -> Result<()> {
        let mut t = std::io::BufWriter::new(std::fs::File::create(dir.join("reps.tsv"))?);
        writeln!(t, "idx\ttid\tchrom\tstart\tend\tn_reads\tstrand\tintrons\tdistinguishing_uniq\tcore_bp\tstub\ttes\tseq_len")?;
        let mut f = std::io::BufWriter::new(std::fs::File::create(dir.join("reps.fa"))?);
        for (i, r) in reps.iter().enumerate() {
            let introns: Vec<String> = r.introns.iter().map(|(d, a)| format!("{d}-{a}")).collect();
            writeln!(
                t,
                "{i}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                r.tid,
                r.chrom,
                r.start,
                r.end,
                r.n_reads,
                r.strand,
                if introns.is_empty() {
                    "-".to_string()
                } else {
                    introns.join(",")
                },
                r.distinguishing_uniq,
                r.core_bp,
                r.stub as u8,
                r.tes.map_or_else(|| "-".to_string(), |v| v.to_string()),
                r.seq.len()
            )?;
            writeln!(f, ">{i}")?;
            f.write_all(&r.seq)?;
            writeln!(f)?;
        }
        t.flush()?;
        f.flush()?;
        Ok(())
    }

    /// Read back what [`write_reps`] wrote.
    pub fn read_reps(dir: &Path) -> Result<Vec<DenovoTranscript>> {
        let mut seqs: Vec<Vec<u8>> = Vec::new();
        for (k, line) in std::io::BufReader::new(std::fs::File::open(dir.join("reps.fa"))?)
            .split(b'\n')
            .enumerate()
        {
            let line = line?;
            if k % 2 == 1 {
                seqs.push(line);
            }
        }
        let mut out = Vec::new();
        for (k, line) in std::io::BufReader::new(std::fs::File::open(dir.join("reps.tsv"))?)
            .lines()
            .enumerate()
        {
            let line = line?;
            if k == 0 {
                continue;
            }
            let c: Vec<&str> = line.split('\t').collect();
            anyhow::ensure!(c.len() == 13, "reps.tsv: bad row {k}");
            let i: usize = c[0].parse()?;
            let introns: Vec<(u64, u64)> = if c[7] == "-" {
                Vec::new()
            } else {
                c[7].split(',')
                    .map(|p| {
                        let (d, a) = p.split_once('-').context("intron")?;
                        Ok((d.parse()?, a.parse()?))
                    })
                    .collect::<Result<_>>()?
            };
            let seq = seqs
                .get(i)
                .cloned()
                .context("reps.fa shorter than reps.tsv")?;
            anyhow::ensure!(
                seq.len() == c[12].parse::<usize>()?,
                "reps.fa/reps.tsv length mismatch at {i}"
            );
            out.push(DenovoTranscript {
                tid: c[1].to_string(),
                chrom: c[2].to_string(),
                start: c[3].parse()?,
                end: c[4].parse()?,
                n_reads: c[5].parse()?,
                strand: c[6].chars().next().context("strand")?,
                introns,
                seq,
                distinguishing_uniq: c[8].parse()?,
                core_bp: c[9].parse()?,
                stub: c[10] == "1",
                tes: if c[11] == "-" {
                    None
                } else {
                    Some(c[11].parse()?)
                },
            });
        }
        anyhow::ensure!(
            out.len() == seqs.len(),
            "reps.tsv/reps.fa row count mismatch"
        );
        Ok(out)
    }

    /// `minimap2 --version`, once per process (part of every PAF key).
    pub fn minimap2_version(minimap2: &str) -> String {
        static V: std::sync::OnceLock<std::sync::Mutex<std::collections::HashMap<String, String>>> =
            std::sync::OnceLock::new();
        let m = V.get_or_init(Default::default);
        if let Some(v) = m.lock().unwrap().get(minimap2) {
            return v.clone();
        }
        let v = std::process::Command::new(minimap2)
            .arg("--version")
            .output()
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
            .unwrap_or_else(|_| "unknown".into());
        m.lock().unwrap().insert(minimap2.to_string(), v.clone());
        v
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn tx(i: u64, introns: Vec<(u64, u64)>, seq: &[u8]) -> DenovoTranscript {
            DenovoTranscript {
                tid: format!("DN_c_{i}"),
                chrom: "c1".into(),
                start: 100 * i,
                end: 100 * i + 50,
                n_reads: i as u32 + 2,
                strand: if i % 2 == 0 { '+' } else { '-' },
                introns,
                seq: seq.to_vec(),
                distinguishing_uniq: i as usize,
                core_bp: 7 * i,
                stub: i % 3 == 0,
                tes: if i % 2 == 0 { Some(100 * i + 60) } else { None },
            }
        }

        #[test]
        fn reps_round_trip_exactly() {
            let dir =
                std::env::temp_dir().join(format!("rustle_run_cache_test_{}", std::process::id()));
            let _ = std::fs::remove_dir_all(&dir);
            std::fs::create_dir_all(&dir).unwrap();
            let reps = vec![
                tx(0, vec![], b"ACGTacgtNN"),
                tx(1, vec![(110, 120), (130, 140)], b"GGGccc"),
                tx(3, vec![], b""),
            ];
            write_reps(&dir, &reps).unwrap();
            let back = read_reps(&dir).unwrap();
            assert_eq!(back.len(), reps.len());
            for (a, b) in reps.iter().zip(&back) {
                assert_eq!(
                    (
                        &a.tid,
                        &a.chrom,
                        a.start,
                        a.end,
                        a.n_reads,
                        a.strand,
                        &a.introns,
                        &a.seq,
                        a.distinguishing_uniq,
                        a.core_bp,
                        a.stub,
                        a.tes
                    ),
                    (
                        &b.tid,
                        &b.chrom,
                        b.start,
                        b.end,
                        b.n_reads,
                        b.strand,
                        &b.introns,
                        &b.seq,
                        b.distinguishing_uniq,
                        b.core_bp,
                        b.stub,
                        b.tes
                    )
                );
            }
            let _ = std::fs::remove_dir_all(&dir);
        }

        #[test]
        fn a_hit_needs_done_and_the_identical_key() {
            let root =
                std::env::temp_dir().join(format!("rustle_run_cache_hit_{}", std::process::id()));
            let _ = std::fs::remove_dir_all(&root);
            let e = Entry::new(&root, "reps", "key one\n".into());
            assert!(!e.is_hit());
            let st = e.staging().unwrap();
            assert!(!e.is_hit(), "staging is not a hit");
            write_reps(&st, &[]).unwrap(); // a reps entry needs its payload (reps.tsv + reps.fa) to be complete
            e.commit(&st).unwrap();
            assert!(e.is_hit());
            // same directory name forced, different key text: must miss
            let other = Entry {
                dir: e.dir.clone(),
                key: "key two\n".into(),
                required: e.required,
                pin: false,
            };
            assert!(!other.is_hit());
            // a file truncated after commit: must miss
            let e2 = Entry::new(&root, "paf", "k\n".into());
            let st2 = e2.staging().unwrap();
            std::fs::write(st2.join("out.paf"), b"line one\nline two\n").unwrap();
            e2.commit(&st2).unwrap();
            assert!(e2.is_hit());
            std::fs::write(e2.dir.join("out.paf"), b"line one\n").unwrap();
            assert!(!e2.is_hit(), "a truncated file must not replay");
            // an empty DONE (crash between rename and write-back) is a miss, not a vacuous hit
            std::fs::write(e2.dir.join("out.paf"), b"line one\nline two\n").unwrap();
            assert!(e2.is_hit());
            std::fs::write(e2.dir.join("DONE"), b"").unwrap();
            assert!(!e2.is_hit(), "an empty DONE must not be a hit");
            // a reps entry whose DONE lacks reps.fa is a miss
            let e3 = Entry::new(&root, "reps", "r\n".into());
            let st3 = e3.staging().unwrap();
            std::fs::write(st3.join("reps.tsv"), b"idx\n").unwrap();
            e3.commit(&st3).unwrap();
            assert!(
                !e3.is_hit(),
                "a reps entry without reps.fa must not be a hit"
            );
            let _ = std::fs::remove_dir_all(&root);
        }

        /// The `cand` kind (the `o3_candidates` result) requires `candidates.tsv` and `contigs.fa`: a half-written entry is a miss, and an
        /// EMPTY `contigs.fa` (a run that flagged no candidate) is a complete payload, not a missing one.
        #[test]
        fn a_cand_entry_is_a_hit_only_with_both_its_candidates_table_and_its_contigs() {
            let dir = tempfile::tempdir().unwrap();
            let e = Entry::new(dir.path(), "cand", "rustle o3 candidates v1\n".into());
            assert_eq!(e.required, ["candidates.tsv", "contigs.fa"]);
            assert!(
                e.dir.starts_with(dir.path().join("cand")),
                "{}",
                e.dir.display()
            );
            let st = e.staging().unwrap();
            std::fs::write(st.join("candidates.tsv"), b"family\n").unwrap();
            e.commit(&st).unwrap();
            assert!(
                !e.is_hit(),
                "a cand entry without contigs.fa must not be a hit"
            );
            let st = e.staging().unwrap();
            std::fs::write(st.join("candidates.tsv"), b"family\n").unwrap();
            std::fs::write(st.join("contigs.fa"), b"").unwrap();
            e.commit(&st).unwrap();
            assert!(
                e.is_hit(),
                "both payloads present (one of them empty): a hit"
            );
        }

        #[test]
        fn content_hash_is_the_byte_stream_whatever_the_split() {
            let data: Vec<u8> = (0..1000u32)
                .map(|i| (i.wrapping_mul(2_654_435_761) >> 13) as u8)
                .collect();
            let mut one = ContentHash::default();
            one.update(&data);
            for split in [
                &[1usize, 7, 8, 9, 100][..],
                &[3, 3, 3, 3],
                &[999],
                &[0, 0, 5, 0],
            ] {
                let mut h = ContentHash::default();
                let mut at = 0;
                for &n in split {
                    let n = n.min(data.len() - at);
                    h.update(&data[at..at + n]);
                    at += n;
                }
                h.update(&data[at..]);
                assert_eq!((h.hex(), h.len()), (one.hex(), 1000));
            }
            // zero padding of the last word is not ambiguous (the length is hashed), every byte counts
            let hx = |b: &[u8]| {
                let mut h = ContentHash::default();
                h.update(b);
                h.hex()
            };
            assert_ne!(hx(b"ab"), hx(b"ab\0"));
            assert_ne!(hx(b""), hx(b"\0"));
            let mut flipped = data.clone();
            for i in [0usize, 7, 8, 500, 999] {
                flipped[i] ^= 1;
                assert_ne!(hx(&flipped), one.hex(), "byte {i}");
                flipped[i] ^= 1;
            }
            // a HashingWriter hashes exactly what it writes; disabled, it only passes the bytes on
            let mut w = HashingWriter::new(Vec::new(), true);
            w.write_all(&data[..10]).unwrap();
            w.write_all(&data[10..]).unwrap();
            assert_eq!(
                (w.inner.as_slice(), w.hash.unwrap().hex()),
                (&data[..], one.hex())
            );
            let mut off = HashingWriter::new(Vec::new(), false);
            off.write_all(&data).unwrap();
            assert!(off.hash.is_none() && off.inner == data);
        }

        /// A pinned entry's payload is replayed as a hard link; `dest` is unlinked, never truncated; a write through the
        /// link is a miss (mtime), a same-size rewrite with the mtime put back is a miss when it touches a sampled block,
        /// and the one edge left (outside every sampled block) is caught by the verify (full re-hash) mode.
        #[test]
        fn a_pinned_entry_replays_by_link_and_a_write_through_the_link_is_a_miss() {
            let root =
                std::env::temp_dir().join(format!("rustle_run_cache_pin_{}", std::process::id()));
            let _ = std::fs::remove_dir_all(&root);
            std::fs::create_dir_all(&root).unwrap();
            let data: Vec<u8> = (0..(1u32 << 20))
                .map(|i| b"ACGT\t\n"[(i % 6) as usize])
                .collect();
            let product = root.join("run.loci.paf");
            std::fs::write(&product, &data).unwrap();
            let e = Entry::new(&root.join("cache"), "paf", "families k\n".into()).pinned();
            let st = e.staging().unwrap();
            assert!(
                e.stage_link(&st, "out.paf", &product).unwrap(),
                "same file system: a link, not a copy"
            );
            e.commit(&st).unwrap();
            assert!(e.is_hit_verify(false) && e.is_hit_verify(true));
            let done = std::fs::read_to_string(e.dir.join("DONE")).unwrap();
            assert!(
                done.lines().all(|l| l.split('\t').count() == 6),
                "pinned DONE rows: {done}"
            );
            // the same entry seen unpinned-style (a DONE without pins) must not satisfy a pinned lookup
            let plain: String = done
                .lines()
                .map(|l| l.split('\t').take(2).collect::<Vec<_>>().join("\t") + "\n")
                .collect();
            std::fs::write(e.dir.join("DONE"), &plain).unwrap();
            assert!(
                !e.is_hit_verify(false),
                "a pinned entry without pins is a miss"
            );
            std::fs::write(e.dir.join("DONE"), &done).unwrap();
            assert!(e.is_hit_verify(false));
            #[cfg(unix)]
            let ino = |p: &Path| {
                use std::os::unix::fs::MetadataExt;
                std::fs::metadata(p).unwrap().ino()
            };
            #[cfg(unix)]
            assert_eq!(
                ino(&product),
                ino(&e.dir.join("out.paf")),
                "the product and the payload are one inode"
            );
            // a SECOND prefix sharing the cache, whose old product names another inode: that inode is left intact, and
            // the payload (already linked to `product`) is copied, not linked, so the two products never alias
            let other = root.join("again.loci.paf");
            std::fs::write(&other, b"old product").unwrap();
            let keep = root.join("keep");
            std::fs::hard_link(&other, &keep).unwrap();
            assert!(
                !e.replay("out.paf", &other).unwrap(),
                "already linked to another product: a copy"
            );
            assert_eq!(std::fs::read(&keep).unwrap(), b"old product");
            assert_eq!(std::fs::read(&other).unwrap(), data);
            #[cfg(unix)]
            assert_ne!(ino(&other), ino(&e.dir.join("out.paf")));
            // the same prefix re-run: its own product is unlinked first, so the payload is unshared again and linked
            let dest = product.clone();
            assert!(
                e.replay("out.paf", &dest).unwrap(),
                "a re-run of the linked prefix links again"
            );
            assert_eq!(std::fs::read(&dest).unwrap(), data);
            #[cfg(unix)]
            assert_eq!(ino(&dest), ino(&e.dir.join("out.paf")));
            assert!(
                e.is_hit_verify(true),
                "a link changes neither mtime nor content"
            );
            // an in-place same-size write through the link (a later run of an older binary, a shell redirect)
            let mtime = std::fs::metadata(&dest).unwrap().modified().unwrap();
            std::thread::sleep(std::time::Duration::from_millis(30));
            let mut w = data.clone();
            w[0] = b'X';
            std::fs::write(&dest, &w).unwrap();
            assert!(!e.is_hit_verify(false), "the mtime moved");
            // the same write with the mtime put back: caught by the sampled block at offset 0
            std::fs::File::options()
                .write(true)
                .open(&dest)
                .unwrap()
                .set_modified(mtime)
                .unwrap();
            assert!(!e.is_hit_verify(false), "a sampled block changed");
            // the documented edge: a change outside every sampled block, mtime restored -> only the audit mode sees it
            let mut w2 = data.clone();
            w2[100_000] = b'X';
            std::fs::write(&dest, &w2).unwrap();
            std::fs::File::options()
                .write(true)
                .open(&dest)
                .unwrap()
                .set_modified(mtime)
                .unwrap();
            assert!(
                e.is_hit_verify(false),
                "the documented edge (same size, mtime restored, no sampled block touched)"
            );
            assert!(
                !e.is_hit_verify(true),
                "RUSTLE_CACHE_VERIFY=1 re-hashes the payload and rejects it"
            );
            let _ = std::fs::remove_dir_all(&root);
        }

        #[test]
        fn fnv_incremental_matches_one_shot() {
            let mut f = Fnv::default();
            f.update(b"abc");
            f.update(b"def");
            assert_eq!(f.finish(), fnv1a64(b"abcdef"));
        }
    }
}
