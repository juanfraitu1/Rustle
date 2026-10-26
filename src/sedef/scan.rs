//! `utilities sd-scan`: the four sedef-rs stages in one call, through files, on a thread pool of its own
//! (docs/superpowers/specs/2026-10-13-o1-layered-hierarchy-design.md, section 1). For `--out PREFIX`:
//!
//! - `PREFIX.sd.candidates.tsv`  stage 1 (seed): candidate window pairs
//! - `PREFIX.sd.blocks.tsv`      stage 2 (chain): SD blocks
//! - `PREFIX.sd.pairs.tsv`       stage 3 (refine): bp-resolution SD pairs, strand, identity (percent), status ok/self/weak
//! - `PREFIX.sd.anchors.bin`     stage 3: the anchors stage 4 reads
//! - `PREFIX.sd.segments.tsv`    stage 4 (decompose): elementary segment pairs
//! - `PREFIX.sd.families.tsv`    stage 4: duplication families; `n_copies` counts segment endpoints before duplicates are
//!   removed (two segments sharing an endpoint count it twice), so it can exceed the family's rows in duplicons.tsv
//! - `PREFIX.sd.duplicons.tsv`   one row per distinct family copy, `chrom start end F<id>`: the `--hierarchy-duplicons`
//!   format (`family_hierarchy::DupIndex::from_text`)
//! - `PREFIX.sd.params.tsv`      the version key, the inputs, the counts and the seconds per stage
//!
//! With the default parameters (the sandbox CLI's) on a genome of ONE sequence, the candidates, blocks, pairs and
//! segments equal the sandbox's `seed`, `chain`, `refine --anchors-out` and `decompose --anchors` outputs byte for byte,
//! except that a sequence is named by the first word of its FASTA header (the sandbox kept the whole header line). On
//! several sequences the vendored chain and dedup fixes (src/sedef/mod.rs) change the output.

use anyhow::{anyhow, bail, Context, Result};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use crate::sedef::{anchors, chain, decompose, fasta, refine, seed};

/// Bumped whenever an output of `run` can change for the same genome and parameters; part of the version key.
pub const SCAN_VERSION: u32 = 2;

/// The sandbox CLI's default parameters (`sedef_rs seed|chain|refine|decompose --help`), plus the output-neutral
/// shard count of the seed stage's pair counting. Scan v2 (2026-10-25, LEDGER_O1 SDSCAN2): seeding windows every
/// `seed_stride` = half a window, and `max_postings` doubled with it (a k-mer of one copy now sits in two windows).
#[derive(Clone, Debug, PartialEq)]
pub struct ScanParams {
    pub k: usize,
    pub window: usize,
    pub seed_stride: usize,
    /// Scan v2 (BISER's `mask`): seed, chain and refine on the genome without its soft-masked (lowercase) bases, and
    /// map pairs, segments and duplicons back to genome coordinates. Common repeats then seed nothing, and a repeat
    /// inside a duplication no longer splits it.
    pub masked_space: bool,
    /// Scan v2 (SEDEF): a candidate needs a shared all-unmasked seed k-mer (`seed::SeedParams::uppercase_seed`).
    pub seed_uppercase: bool,
    /// Scan v2 (SEDEF's `--uppercase`, default 100): unmasked bases per side and unmasked matches an ok pair needs.
    pub min_uppercase: u64,
    pub sketch: usize,
    pub min_kmers: usize,
    pub max_postings: usize,
    pub min_shared: usize,
    pub chain_gap: u32,
    pub chain_min_score: u32,
    pub chain_merge_frac: f64,
    pub min_cell_id: f64,
    pub min_block_id: f64,
    pub refine_max_gap: u64,
    pub refine_diag_bin: i64,
    pub refine_min_span: u64,
    pub decompose_max_gap: u64,
    pub decompose_diag_bin: i64,
    pub min_seg: u64,
    pub max_pairs: u64,
    pub count_shards: usize,
}

impl Default for ScanParams {
    fn default() -> Self {
        ScanParams {
            k: 16,
            window: 2000,
            seed_stride: 1000,
            masked_space: false,
            seed_uppercase: true,
            min_uppercase: 100,
            sketch: 256,
            min_kmers: 32,
            max_postings: 200,
            min_shared: 40,
            chain_gap: 2,
            chain_min_score: 20,
            chain_merge_frac: 0.5,
            min_cell_id: 88.0,
            min_block_id: 90.0,
            refine_max_gap: 400,
            refine_diag_bin: 50,
            refine_min_span: 500,
            decompose_max_gap: 200,
            decompose_diag_bin: 50,
            min_seg: 1000,
            max_pairs: 1024,
            count_shards: 256,
        }
    }
}

impl ScanParams {
    /// Everything that can change an output, on one line; the `sd` stage hashes it into its cache key.
    /// `count_shards` is left out: it changes no output.
    pub fn version_key(&self) -> String {
        format!(
            "sedef-rs 0.1.0 vendored 2026-10-13 scan v{SCAN_VERSION}; masked={} upper_seeds={} min_upper={}; seed k={} window={} stride={} sketch={} min_kmers={} \
             max_postings={} min_shared={}; chain gap={} min_score={} merge_frac={} min_cell_id={} min_block_id={}; \
             refine max_gap={} diag_bin={} min_span={}; decompose max_gap={} diag_bin={} min_seg={}; max_pairs={}; \
             names=first-word",
            self.masked_space as u8,
            self.seed_uppercase as u8,
            self.min_uppercase,
            self.k,
            self.window,
            self.seed_stride,
            self.sketch,
            self.min_kmers,
            self.max_postings,
            self.min_shared,
            self.chain_gap,
            self.chain_min_score,
            self.chain_merge_frac,
            self.min_cell_id,
            self.min_block_id,
            self.refine_max_gap,
            self.refine_diag_bin,
            self.refine_min_span,
            self.decompose_max_gap,
            self.decompose_diag_bin,
            self.min_seg,
            self.max_pairs
        )
    }
}

/// What `run` wrote, counted.
#[derive(Clone, Debug, PartialEq)]
pub struct ScanSummary {
    pub sequences: usize,
    pub bases: u64,
    pub candidates: usize,
    pub blocks: usize,
    pub pairs: usize,
    pub pairs_ok: usize,
    pub segments: usize,
    pub families: usize,
    pub duplicon_rows: usize,
}

/// Per sequence, the kept runs of the masked projection as (projected start, genome start), ascending; a sequence
/// without an entry is unprojected.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Projection {
    pub runs: HashMap<String, Vec<(u64, u64)>>,
}

impl Projection {
    /// Genome position of projected position `p`.
    pub fn start(&self, chr: &str, p: u64) -> u64 {
        match self.runs.get(chr) {
            None => p,
            Some(r) => {
                let i = r.partition_point(|&(ps, _)| ps <= p).saturating_sub(1);
                r.get(i).map_or(p, |&(ps, gs)| gs + (p - ps))
            }
        }
    }

    /// Genome end of a projected interval ending (exclusive) at `p`.
    pub fn end(&self, chr: &str, p: u64) -> u64 {
        if p == 0 {
            0
        } else {
            self.start(chr, p - 1) + 1
        }
    }
}

/// A sequence without its lowercase bases, and its runs (see `Projection`).
pub fn project(seq: &[u8]) -> (Vec<u8>, Vec<(u64, u64)>) {
    let mut out = Vec::with_capacity(seq.len());
    let mut runs = Vec::new();
    let mut prev_kept = false;
    for (i, &c) in seq.iter().enumerate() {
        let keep = !c.is_ascii_lowercase();
        if keep {
            if !prev_kept {
                runs.push((out.len() as u64, i as u64));
            }
            out.push(c);
        }
        prev_kept = keep;
    }
    (out, runs)
}

/// Project every sequence in place when `masked` (each shrunk to its kept bases); the identity projection otherwise.
pub fn project_genome(seqs: &mut [(String, Vec<u8>)], masked: bool) -> Projection {
    let mut proj = Projection::default();
    if masked {
        for (name, seq) in seqs.iter_mut() {
            let (p, runs) = project(seq);
            *seq = p;
            seq.shrink_to_fit();
            proj.runs.insert(name.clone(), runs);
        }
    }
    proj
}

/// Every FASTA record (plain, gzip or bgzip), named by the first word of its header line, each sequence shrunk to its
/// length. A record without a name, two records with one name, sequence data before the first '>' (or a file that is
/// not FASTA), or a file without records is refused.
pub fn read_genome(path: &Path) -> Result<Vec<(String, Vec<u8>)>> {
    let recs = fasta::read_fasta(path, None).map_err(|e| anyhow!("reading {}: {e}", path.display()))?;
    let mut seen = HashSet::new();
    let mut out = Vec::with_capacity(recs.len());
    for (header, seq) in recs {
        let name = header.split_whitespace().next().unwrap_or("").to_string();
        if name.is_empty() {
            bail!("{}: a FASTA record has no name", path.display());
        }
        if !seen.insert(name.clone()) {
            bail!("{}: two FASTA records are named {name}", path.display());
        }
        let mut seq = seq;
        seq.shrink_to_fit(); // the reader grows each sequence by doubling: up to 2x its length
        out.push((name, seq));
    }
    if out.is_empty() {
        bail!("{}: no FASTA records", path.display());
    }
    Ok(out)
}

/// The four stages on a pool of `threads` threads of its own (the sandbox's `build_global` can run once per process).
pub fn run(genome: &Path, out: &str, threads: usize, p: &ScanParams) -> Result<ScanSummary> {
    let threads = threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("building the sd-scan thread pool")?;
    pool.install(|| run_stages(genome, out, threads, p, true))
}

/// The first three stages only (candidates, blocks, pairs, anchors and params tables; no segments, families or
/// duplicons): for scoring the pairs without paying for decompose.
pub fn run_pairs_only(genome: &Path, out: &str, threads: usize, p: &ScanParams) -> Result<ScanSummary> {
    let threads = threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("building the sd-scan thread pool")?;
    pool.install(|| run_stages(genome, out, threads, p, false))
}

fn run_stages(genome: &Path, out: &str, threads: usize, p: &ScanParams, decompose: bool) -> Result<ScanSummary> {
    let f = |ext: &str| format!("{out}.sd.{ext}");
    let mut secs: Vec<(&str, f64)> = Vec::new();

    let t = Instant::now();
    let (sequences, bases, candidates, proj) = seed_stage(genome, &f("candidates.tsv"), p)?;
    secs.push(("seed", t.elapsed().as_secs_f64()));
    eprintln!("[sd-scan] seed: {sequences} sequences, {bases} bp, {candidates} candidate window pairs ({:.1} s)", secs[0].1);

    let t = Instant::now();
    let blocks = chain_stage(&f("candidates.tsv"), &f("blocks.tsv"), p)?;
    secs.push(("chain", t.elapsed().as_secs_f64()));
    eprintln!("[sd-scan] chain: {blocks} blocks ({:.1} s)", secs[1].1);

    let t = Instant::now();
    let (pairs, pairs_ok) = refine_stage(genome, &f("blocks.tsv"), &f("pairs.tsv"), &f("anchors.bin"), p, &proj)?;
    secs.push(("refine", t.elapsed().as_secs_f64()));
    eprintln!("[sd-scan] refine: {pairs} pairs, {pairs_ok} ok ({:.1} s)", secs[2].1);

    let (mut segments, mut families, mut duplicon_rows) = (0, 0, 0);
    if decompose {
        let t = Instant::now();
        (segments, families, duplicon_rows) =
            decompose_stage(&f("anchors.bin"), &f("segments.tsv"), &f("families.tsv"), &f("duplicons.tsv"), p, &proj)?;
        secs.push(("decompose", t.elapsed().as_secs_f64()));
        eprintln!(
            "[sd-scan] decompose: {segments} segments, {families} families, {duplicon_rows} duplicon rows ({:.1} s)",
            secs[3].1
        );
    }

    let s = ScanSummary { sequences, bases, candidates, blocks, pairs, pairs_ok, segments, families, duplicon_rows };
    write_params(&f("params.tsv"), genome, threads, p, &s, &secs)?;
    Ok(s)
}

fn create(path: &str) -> Result<std::io::BufWriter<std::fs::File>> {
    Ok(std::io::BufWriter::new(std::fs::File::create(path).with_context(|| format!("creating {path}"))?))
}

/// Stage 1, as the sandbox's `run_seed`; the genome is released before the pairs are counted.
/// Returns (sequences, bases, candidates).
fn seed_stage(genome: &Path, out: &str, p: &ScanParams) -> Result<(usize, u64, usize, Projection)> {
    let mut seqs = read_genome(genome)?;
    let sequences = seqs.len();
    let bases: u64 = seqs.iter().map(|(_, s)| s.len() as u64).sum();
    let proj = project_genome(&mut seqs, p.masked_space);
    let sp = seed::SeedParams {
        k: p.k,
        window: p.window,
        stride: p.seed_stride,
        uppercase_seed: p.seed_uppercase,
        sketch: p.sketch,
        min_kmers: p.min_kmers,
        max_postings: p.max_postings,
        min_shared: p.min_shared,
        count_shards: p.count_shards,
    };
    let (mut wins, infos) = seed::build_windows(&seqs, &sp);
    drop(seqs);
    let mut cands = seed::find_candidates(&mut wins, &sp);
    sort_candidates(&mut cands, &wins);
    let mut w = create(out)?;
    writeln!(w, "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tshared\test_jaccard\tupper_seed")?;
    for c in &cands {
        let (wa, wb) = (&wins[c.a], &wins[c.b]);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{}",
            infos[wa.seq as usize].name,
            wa.start,
            wa.end,
            infos[wb.seq as usize].name,
            wb.start,
            wb.end,
            c.shared,
            seed::est_jaccard(c.shared, sp.sketch),
            c.upper_seed as u8,
        )?;
    }
    w.flush()?;
    Ok((sequences, bases, cands.len(), proj))
}

/// The sandbox's candidate order (most shared first, then window a's sequence and start, then window b's start),
/// with window b's sequence as the last key: the sandbox left ties between sequences to its single map's order,
/// which the shards change. On a one-sequence genome there are no such ties.
pub fn sort_candidates(cands: &mut [seed::Candidate], wins: &[seed::WindowSketch]) {
    cands.sort_by(|x, y| {
        y.shared.cmp(&x.shared).then_with(|| {
            (wins[x.a].seq, wins[x.a].start, wins[x.b].start, wins[x.b].seq)
                .cmp(&(wins[y.a].seq, wins[y.a].start, wins[y.b].start, wins[y.b].seq))
        })
    });
}

/// Stage 2, as the sandbox's `run_chain`. Returns the number of blocks.
fn chain_stage(cands: &str, out: &str, p: &ScanParams) -> Result<usize> {
    let min_shared_cell = chain::shared_for_identity(p.min_cell_id / 100.0, p.sketch, p.k);
    let pairs = chain::load_pairs(Path::new(cands), min_shared_cell).map_err(|e| anyhow!("reading {cands}: {e}"))?;
    let cp = chain::ChainParams {
        window: p.window as u64,
        stride: p.seed_stride as u64,
        gap: p.chain_gap,
        min_score: p.chain_min_score,
        merge_frac: p.chain_merge_frac,
        sketch: p.sketch,
        k: p.k,
        min_cell_id: p.min_cell_id,
        min_block_id: p.min_block_id,
        min_shared_cell,
        require_upper_seed: p.seed_uppercase,
    };
    let (mut blocks, _) = chain::chain_and_merge(&pairs, &cp);
    blocks.sort_by(|x, y| {
        y.score
            .cmp(&x.score)
            .then_with(|| (&x.chr_a, x.a0, &x.chr_b, x.b0).cmp(&(&y.chr_a, y.a0, &y.chr_b, y.b0)))
    });
    let mut w = create(out)?;
    writeln!(w, "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tcells\tscore\test_id")?;
    for b in &blocks {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
            b.chr_a, b.a0, b.a1, b.chr_b, b.b0, b.b1, b.cells, b.score, b.est_id
        )?;
    }
    w.flush()?;
    Ok(blocks.len())
}

/// Stage 3, as the sandbox's `run_refine --anchors-out`. Returns (pairs, ok pairs).
fn refine_stage(
    genome: &Path,
    blocks_path: &str,
    out: &str,
    anchors_out: &str,
    p: &ScanParams,
    proj: &Projection,
) -> Result<(usize, usize)> {
    let blocks = chain::load_blocks(Path::new(blocks_path)).map_err(|e| anyhow!("reading {blocks_path}: {e}"))?;
    let mut genome_seqs = read_genome(genome)?;
    project_genome(&mut genome_seqs, p.masked_space);
    let seqs: HashMap<String, Vec<u8>> = genome_seqs.into_iter().collect();
    let rp = refine::RefineParams {
        k: p.k,
        max_gap: p.refine_max_gap,
        diag_bin: p.refine_diag_bin,
        min_span: p.refine_min_span,
        max_pairs: p.max_pairs,
        min_uppercase: p.min_uppercase,
        min_identity: p.min_block_id,
    };
    let refined_pairs = refine::refine_blocks(&blocks, &seqs, &rp);
    drop(seqs);
    let (refined, anchored): (Vec<refine::Refined>, Vec<anchors::BlockAnchors>) = refined_pairs.into_iter().unzip();
    anchors::write(Path::new(anchors_out), p.k, &anchored).map_err(|e| anyhow!("writing {anchors_out}: {e}"))?;
    drop(anchored);
    let mut w = create(out)?;
    writeln!(w, "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tstrand\tanchors\test_id\tstage2_id\tstatus")?;
    let mut ok = 0;
    for r in &refined {
        let status = if !r.ok {
            "weak"
        } else if r.self_pair {
            "self"
        } else if r.repeat {
            "repeat"
        } else {
            ok += 1;
            "ok"
        };
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{}",
            r.chr_a,
            proj.start(&r.chr_a, r.ra0),
            proj.end(&r.chr_a, r.ra1),
            r.chr_b,
            proj.start(&r.chr_b, r.rb0),
            proj.end(&r.chr_b, r.rb1),
            r.strand,
            r.anchors,
            r.est_id,
            r.stage2_id,
            status
        )?;
    }
    w.flush()?;
    Ok((refined.len(), ok))
}

/// Stage 4, as the sandbox's `decompose --anchors`, plus the duplicon table. Returns (segments, families, rows).
fn decompose_stage(
    anchors_path: &str,
    seg_out: &str,
    fam_out: &str,
    dup_out: &str,
    p: &ScanParams,
    proj: &Projection,
) -> Result<(usize, usize, usize)> {
    let (k, recs) = anchors::read(Path::new(anchors_path)).map_err(|e| anyhow!("reading {anchors_path}: {e}"))?;
    let dp = decompose::DecomposeParams {
        k,
        max_gap: p.decompose_max_gap,
        diag_bin: p.decompose_diag_bin,
        min_seg: p.min_seg,
        max_pairs: p.max_pairs,
    };
    let (mut segs, mut fams, _) = decompose::decompose(recs, &dp);
    if !proj.runs.is_empty() {
        // segments back to genome coordinates; families are built on those (their overlap and span are genome bases)
        for sg in segs.iter_mut() {
            (sg.a0, sg.a1) = (proj.start(&sg.chr_a, sg.a0), proj.end(&sg.chr_a, sg.a1));
            (sg.b0, sg.b1) = (proj.start(&sg.chr_b, sg.b0), proj.end(&sg.chr_b, sg.b1));
        }
        fams = decompose::families(&segs);
    }
    let mut w = create(seg_out)?;
    writeln!(w, "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tstrand\tanchors\test_id")?;
    for s in &segs {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}",
            s.chr_a, s.a0, s.a1, s.chr_b, s.b0, s.b1, s.strand, s.anchors, s.est_id
        )?;
    }
    w.flush()?;
    let mut w = create(fam_out)?;
    writeln!(w, "#family_id\tn_copies\tspan_bp\tchromosomes")?;
    for f in &fams {
        writeln!(w, "F{}\t{}\t{}\t{}", f.id, f.n_copies, f.span_bp, f.chrs.join(";"))?;
    }
    w.flush()?;
    let rows = write_duplicons(dup_out, &fams)?;
    Ok((segs.len(), fams.len(), rows))
}

/// One row per distinct copy of every family, `chrom start end F<id>`, sorted by chromosome, start, end and family:
/// the `--hierarchy-duplicons` input (`family_hierarchy::DupIndex::from_text` skips the `#` header).
pub fn write_duplicons(path: &str, fams: &[decompose::Family]) -> Result<usize> {
    let mut rows: Vec<(&str, u64, u64, usize)> = fams
        .iter()
        .flat_map(|f| f.copies.iter().map(move |(c, s, e)| (c.as_str(), *s, *e, f.id)))
        .collect();
    rows.sort();
    let mut w = create(path)?;
    writeln!(w, "#chrom\tstart\tend\tfamily_id")?;
    for (c, s, e, id) in &rows {
        writeln!(w, "{c}\t{s}\t{e}\tF{id}")?;
    }
    w.flush()?;
    Ok(rows.len())
}

fn write_params(path: &str, genome: &Path, threads: usize, p: &ScanParams, s: &ScanSummary, secs: &[(&str, f64)]) -> Result<()> {
    let mut w = create(path)?;
    writeln!(w, "version_key\t{}", p.version_key())?;
    writeln!(w, "genome\t{}", genome.display())?;
    let bytes = std::fs::metadata(genome).with_context(|| format!("reading {}", genome.display()))?.len();
    writeln!(w, "genome_bytes\t{bytes}")?;
    writeln!(w, "threads\t{threads}")?;
    writeln!(w, "count_shards\t{}", p.count_shards)?;
    for (k, v) in [
        ("sequences", s.sequences as u64),
        ("bases", s.bases),
        ("candidates", s.candidates as u64),
        ("blocks", s.blocks as u64),
        ("pairs", s.pairs as u64),
        ("pairs_ok", s.pairs_ok as u64),
        ("segments", s.segments as u64),
        ("families", s.families as u64),
        ("duplicon_rows", s.duplicon_rows as u64),
    ] {
        writeln!(w, "{k}\t{v}")?;
    }
    for (k, v) in secs {
        writeln!(w, "seconds_{k}\t{v:.1}")?;
    }
    w.flush()?;
    Ok(())
}

/// The test genome of this module's tests and of `tests/sd_scan.rs`. chrA and chrB, 200,000 bp each, from one LCG
/// (x = x * 6364136223846793005 + 1442695040888963407, base = "ACGT"[x >> 62], seed 7). With `planted`:
/// chrA:50,000-58,000, every 200th base changed (A->C->G->T->A), is copied inverted over chrB:120,000-128,000; then
/// chrB:10,000-15,000, every 100th base changed, is copied direct over chrB:170,000-175,000. Headers carry a
/// description after the name. Checked while writing the plan against the sandbox binary (plan Facts 3-4).
#[doc(hidden)]
pub fn toy_genome_fasta(planted: bool) -> String {
    fn seq(n: usize, x: &mut u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                *x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(*x >> 62) as usize]
            })
            .collect()
    }
    fn mutate(s: &[u8], every: usize) -> Vec<u8> {
        s.iter()
            .enumerate()
            .map(|(i, &c)| {
                if i % every == every - 1 {
                    match c {
                        b'A' => b'C',
                        b'C' => b'G',
                        b'G' => b'T',
                        _ => b'A',
                    }
                } else {
                    c
                }
            })
            .collect()
    }
    fn rc(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&c| match c {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            })
            .collect()
    }
    let mut x = 7u64;
    let a = seq(200_000, &mut x);
    let mut b = seq(200_000, &mut x);
    if planted {
        let inv = rc(&mutate(&a[50_000..58_000], 200));
        b.splice(120_000..128_000, inv);
        let dir = mutate(&b[10_000..15_000], 100);
        b.splice(170_000..175_000, dir);
    }
    let mut out = String::new();
    for (name, s) in [("chrA", &a), ("chrB", &b)] {
        out.push_str(&format!(">{name} toy genome\n"));
        for line in s.chunks(80) {
            out.push_str(std::str::from_utf8(line).unwrap());
            out.push('\n');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::Path;

    // The golden output of SCAN_VERSION 2 on the toy genome (v1 split the inverted 8 kb copy into four 2 kb pairs); tests/sd_scan.rs reads the same tables back through the
    // hierarchy's readers. Change these constants only together with a SCAN_VERSION bump: a change here means sd-scan's
    // output for the same genome and parameters changed, and with it the driver's cache key.
    const PLANTED_PAIRS: &str = "#chr_a\tstart_a\tend_a\tchr_b\tstart_b\tend_b\tstrand\tanchors\test_id\tstage2_id\tstatus\n\
chrA\t50000\t57999\tchrB\t120001\t128000\t-\t7360\t99.5\t97.7\tok\n\
chrB\t10000\t14999\tchrB\t170000\t174999\t+\t4200\t99.0\t96.6\tok\n";

    const PLANTED_FAMILIES: &str = "#family_id\tn_copies\tspan_bp\tchromosomes\n\
F0\t2\t15998\tchrA;chrB\n\
F1\t2\t9998\tchrB\n";

    const PLANTED_DUPLICONS: &str = "#chrom\tstart\tend\tfamily_id\n\
chrA\t50000\t57999\tF0\n\
chrB\t10000\t14999\tF1\n\
chrB\t120001\t128000\tF0\n\
chrB\t170000\t174999\tF1\n";

    fn tmp(name: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("sd_scan_{name}_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn the_projection_drops_lowercase_and_maps_back() {
        let (proj, runs) = project(b"ACgtNNaT");
        assert_eq!(proj, b"ACNNT".to_vec());
        let mut m = Projection::default();
        m.runs.insert("c".into(), runs);
        // projected 0..5 = original A0 C1 N4 N5 T7
        assert_eq!([0, 1, 2, 3, 4].map(|p| m.start("c", p)), [0, 1, 4, 5, 7]);
        assert_eq!(m.end("c", 2), 2); // [0, 2) ends after C1
        assert_eq!(m.end("c", 3), 5); // [0, 3) ends after N4
        assert_eq!(m.end("c", 5), 8);
        assert_eq!(m.start("other", 7), 7); // a sequence without runs is unprojected
    }

    /// 30 identical soft-masked 2 kb repeat copies, and a 4 kb unique duplication whose second copy carries a 1.5 kb
    /// soft-masked insertion.
    fn masked_toy_fasta() -> String {
        let mut x = 11u64;
        let mut rnd = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    b"ACGT"[(x >> 62) as usize]
                })
                .collect()
        };
        let rep: Vec<u8> = rnd(2_000).to_ascii_lowercase();
        let ins: Vec<u8> = rnd(1_500).to_ascii_lowercase();
        let sd = rnd(4_000);
        let mut g = rnd(300_000);
        for i in 0..30 {
            let at = 100_000 + i * 6_000;
            g[at..at + 2_000].copy_from_slice(&rep);
        }
        g[20_000..24_000].copy_from_slice(&sd);
        let mut copy = sd[..2_000].to_vec();
        copy.extend_from_slice(&ins);
        copy.extend_from_slice(&sd[2_000..]);
        g[60_000..65_500].copy_from_slice(&copy);
        let mut out = String::from(">chrM1\n");
        for line in g.chunks(80) {
            out.push_str(std::str::from_utf8(line).unwrap());
            out.push('\n');
        }
        out
    }

    fn ok_pairs(text: &str) -> Vec<(u64, u64, u64, u64)> {
        text.lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').collect::<Vec<_>>())
            .filter(|f| f[10] == "ok")
            .map(|f| (f[1].parse().unwrap(), f[2].parse().unwrap(), f[4].parse().unwrap(), f[5].parse().unwrap()))
            .collect()
    }

    #[test]
    fn soft_masked_repeats_pair_nothing_and_a_masked_insertion_stays_inside_the_pair() {
        let d = tmp("masked");
        let fa = d.join("m.fa");
        std::fs::write(&fa, masked_toy_fasta()).unwrap();
        let in_repeats = |s: u64| s >= 99_000;
        // (masked_space, seed_uppercase, min_uppercase): repeat-blind, BISER's projection, SEDEF's unmasked seeds + filter
        for (arm, masked, upper, min_up) in [("blind", false, false, 0), ("biser", true, false, 0), ("sedef", false, true, 100)] {
            let out = d.join(arm).to_str().unwrap().to_string();
            let p = ScanParams { masked_space: masked, seed_uppercase: upper, min_uppercase: min_up, ..ScanParams::default() };
            run_pairs_only(&fa, &out, 2, &p).unwrap();
            let pairs = ok_pairs(&read(&out, "pairs.tsv"));
            let rep = pairs.iter().filter(|q| in_repeats(q.0) && in_repeats(q.2)).count();
            if arm == "blind" {
                assert!(rep > 0, "repeat-blind scan pairs the repeat copies (the reason for the other arms)");
                continue;
            }
            assert_eq!(rep, 0, "{arm}: no repeat pairs");
            assert_eq!(pairs.len(), 1, "{arm}: {pairs:?}");
            let (a0, a1, b0, b1) = pairs[0];
            assert!(a0.abs_diff(20_000) <= 50 && a1.abs_diff(24_000) <= 50, "{arm}: {pairs:?}");
            assert!(b0.abs_diff(60_000) <= 50 && b1.abs_diff(65_500) <= 50, "{arm}: {pairs:?}");
        }
    }

    fn scan_toy(dir: &Path, name: &str, planted: bool) -> (ScanSummary, String) {
        let fa = dir.join(format!("{name}.fa"));
        std::fs::write(&fa, toy_genome_fasta(planted)).unwrap();
        let out = dir.join(name).to_str().unwrap().to_string();
        let s = run(&fa, &out, 2, &ScanParams::default()).expect("sd-scan");
        (s, out)
    }

    fn read(out: &str, ext: &str) -> String {
        std::fs::read_to_string(format!("{out}.sd.{ext}")).unwrap()
    }

    #[test]
    fn the_planted_duplications_come_out_as_pairs_and_two_copy_families() {
        let d = tmp("planted");
        let (s, out) = scan_toy(&d, "toy", true);
        assert_eq!(
            (s.sequences, s.bases, s.candidates, s.blocks, s.pairs, s.pairs_ok, s.segments, s.families, s.duplicon_rows),
            (2, 400_000, 41, 2, 2, 2, 2, 2, 4)
        );
        assert_eq!(read(&out, "pairs.tsv"), PLANTED_PAIRS);
        assert_eq!(read(&out, "families.tsv"), PLANTED_FAMILIES);
        assert_eq!(read(&out, "duplicons.tsv"), PLANTED_DUPLICONS);
        assert!(read(&out, "params.tsv").starts_with(&format!("version_key\t{}\n", ScanParams::default().version_key())));
    }

    #[test]
    fn two_runs_in_one_process_write_identical_tables() {
        let d = tmp("twice");
        let (_, a) = scan_toy(&d, "one", true);
        let (_, b) = scan_toy(&d, "two", true);
        for ext in ["candidates.tsv", "blocks.tsv", "pairs.tsv", "anchors.bin", "segments.tsv", "families.tsv", "duplicons.tsv"] {
            assert_eq!(
                std::fs::read(format!("{a}.sd.{ext}")).unwrap(),
                std::fs::read(format!("{b}.sd.{ext}")).unwrap(),
                "{ext}"
            );
        }
    }

    #[test]
    fn a_gzipped_genome_gives_the_same_tables() {
        let d = tmp("gz");
        let (_, plain) = scan_toy(&d, "plain", true);
        let gz = d.join("toy.fa.gz");
        let mut enc = flate2::write::GzEncoder::new(std::fs::File::create(&gz).unwrap(), flate2::Compression::default());
        enc.write_all(toy_genome_fasta(true).as_bytes()).unwrap();
        enc.finish().unwrap();
        let out = d.join("gz").to_str().unwrap().to_string();
        run(&gz, &out, 2, &ScanParams::default()).expect("sd-scan on .fa.gz");
        for ext in ["candidates.tsv", "pairs.tsv", "segments.tsv", "duplicons.tsv"] {
            assert_eq!(read(&plain, ext), read(&out, ext), "{ext}");
        }
    }

    #[test]
    fn a_genome_without_duplications_writes_header_only_tables() {
        let d = tmp("none");
        let (s, out) = scan_toy(&d, "none", false);
        assert_eq!((s.candidates, s.blocks, s.pairs, s.segments, s.families, s.duplicon_rows), (0, 0, 0, 0, 0, 0));
        for ext in ["candidates.tsv", "blocks.tsv", "pairs.tsv", "segments.tsv", "families.tsv", "duplicons.tsv"] {
            assert_eq!(read(&out, ext).lines().count(), 1, "{ext}: header only");
        }
    }

    #[test]
    fn records_are_named_by_the_first_word_and_a_repeated_name_is_refused() {
        let d = tmp("names");
        let fa = d.join("g.fa");
        std::fs::write(&fa, ">chr1 described here\nACGT\n>chr2\nAC\n").unwrap();
        let g = read_genome(&fa).unwrap();
        assert_eq!(g.iter().map(|(n, _)| n.as_str()).collect::<Vec<_>>(), vec!["chr1", "chr2"]);
        assert!(g.iter().all(|(_, s)| s.capacity() == s.len()), "sequences are shrunk to their length");
        std::fs::write(&fa, ">chr1 first\nACGT\n>chr1 second\nACGT\n").unwrap();
        let e = read_genome(&fa).unwrap_err().to_string();
        assert!(e.contains("two FASTA records are named chr1"), "{e}");
    }

    #[test]
    fn sequence_before_the_first_header_is_refused() {
        // vendored change (2026-10-13): the sandbox's reader glued data before the first '>' onto record 1 (and a file
        // that is not FASTA, e.g. a bgzip file read as text, gave a genome of one garbage record)
        let d = tmp("headerless");
        let fa = d.join("g.fa");
        for text in ["ACGT\n>chr1\nAC\n", "ACGT\n", "\n \nAC\n>chr1\nAC\n"] {
            std::fs::write(&fa, text).unwrap();
            let e = read_genome(&fa).unwrap_err().to_string();
            assert!(e.contains("sequence data before the first '>'"), "{text:?}: {e}");
        }
        // blank lines before the first '>' hold no sequence: the genome is unchanged
        std::fs::write(&fa, "\n\r\n>chr1\nAC\nGT\n>chr2\nA\n").unwrap();
        let g = read_genome(&fa).unwrap();
        assert_eq!(g, vec![("chr1".to_string(), b"ACGT".to_vec()), ("chr2".to_string(), b"A".to_vec())]);
    }

    #[test]
    fn the_version_key_names_every_parameter_but_the_shard_count() {
        let p = ScanParams::default();
        let k = p.version_key();
        assert!(k.starts_with("sedef-rs 0.1.0 vendored 2026-10-13 scan v2;"), "{k}");
        assert_eq!(ScanParams { count_shards: 7, ..p.clone() }.version_key(), k);
        assert_ne!(ScanParams { min_seg: 999, ..p }.version_key(), k);
    }
}
