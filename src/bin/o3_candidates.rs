//! `o3_candidates` — O3 candidate copies from each family's read net (spec `docs/superpowers/specs/2026-10-02-o3-candidates-design.md`
//! §5; plan `docs/superpowers/plans/2026-10-02-o3-candidates.md` task 8): the reads of a family -> clusters at delta -> one consensus per
//! cluster -> one refinement pass -> the significance merge -> flag / link against the genome -> components of the new-copy clusters
//! (candidate copies), each represented by the exon union of its clusters. The algorithm is `rustle::vg_family::o3_candidates`; this file is
//! the glue: arguments, the two BAM passes, the batched minimap2 calls, the per-family loop, the outputs and the stage cache.
//!
//! usage: o3_candidates --bam B --fasta G --copies P.fam.copies.tsv --copies-fa P.fam.copies.fa --index G.splice.mmi --out P.cand
//!        [--delta 0.00958] [--max-reads 1000] [--min-cluster 3] [--min-support 6] [--threads 4] [--families F1,F2]
//!
//! Products (`<out>.`): `candidates.tsv`, `clusters.tsv`, `contigs.fa`, `nets.fa` (`write_outputs`), `families.tsv` (one row per family, ruling
//! R8), `reads.tsv` (read -> cluster, spec §4) and `clusters.fa` (the consensus of every cluster in `clusters.tsv`; with `reads.tsv` the input
//! of Amendment 12's A12-2 measure). Temporary files live in `<out>.tmp/` and are removed on success.
//!
//! Environment: `RUSTLE_MINIMAP2` (the minimap2 binary, default `minimap2`), `RUSTLE_CACHE_DIR` (run_cache: the whole result is one `cand`
//! entry replayed on a hit; every minimap2 PAF is a `paf` entry). Exit status 2: a usage error, `--copies` missing / empty / unreadable, an
//! input file that does not exist, minimap2 that cannot be started (spec §8). Deterministic: seeded sampling over sorted names, stable orders;
//! threads only inside minimap2 and the BAM decompression.

use anyhow::{Context, Result};
use noodles_core::{Position, Region};
use rustle::vg_family::catalog_input::{group_families, parse_copies_fa, parse_copies_tsv, CatalogFamily};
use rustle::vg_family::copy_assign::AssignParams;
use rustle::vg_family::o3_candidates::{
    best_by_id_cov, best_by_matches, candidate_id, classify, cluster_reads, components, consensus_from_template, distinguishing_columns, is_flagged,
    minimap2, minimap2_binary, minimap2_keyed, minimizer_sketch, parse_cs, parse_paf, refine_cluster, sample_net, sketch_share,
    union_sequence_with_note, variant_is_real, write_cluster_members, write_family_table, write_outputs, Candidate, ClusterSeq, FamilyCounts,
    FamilyKmerIndex, Fate, PafHit, ATTRIB_MAX_FAMILIES, KMER_K, MIN_UNMAPPED_LEN, MM2_AVA, MM2_GENOME, MM2_MEMBERS, MM2_UNION, SKETCH_W,
};
use rustle::vg_family::run_cache as rc;
use rustle::vg_family::seq_utils::reverse_complement;
use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

const USAGE: &str = "usage: o3_candidates --bam B --fasta G --copies P.fam.copies.tsv --copies-fa P.fam.copies.fa --index G.splice.mmi \
--out P.cand [--delta 0.00958] [--max-reads 1000] [--min-cluster 3] [--min-support 6] [--threads 4] [--families F1,F2]";

/// spec §5.5: the per-column error proxy eps of the significance merge (`read_conflict.rs:77`, the default of `RUSTLE_CONFLICT_EPS`).
const MERGE_EPS: f64 = 0.001;
/// spec §5.5: two cluster consensus sequences are tested for a merge only when their minimizer sketches share at least this fraction.
const MERGE_MIN_SKETCH_SHARE: f64 = 0.5;
/// Every product, stored in and replayed from the run_cache `cand` entry (`candidates.tsv` and `contigs.fa` are the kind's required files).
const PRODUCTS: [&str; 7] = ["candidates.tsv", "contigs.fa", "clusters.tsv", "nets.fa", "families.tsv", "reads.tsv", "clusters.fa"];

/// A usage or input error: exit status 2 (spec §8), with a message naming the argument or the file.
#[derive(Debug)]
struct ExitTwo(String);
impl std::fmt::Display for ExitTwo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { f.write_str(&self.0) }
}
impl std::error::Error for ExitTwo {}
fn exit_two(msg: impl Into<String>) -> anyhow::Error { anyhow::Error::new(ExitTwo(msg.into())) }

fn main() {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    if raw.iter().any(|a| a == "-h" || a == "--help") {
        println!("{USAGE}");
        return;
    }
    if let Err(e) = run(&raw) {
        eprintln!("o3_candidates: {e:#}");
        std::process::exit(if e.downcast_ref::<ExitTwo>().is_some() { 2 } else { 1 });
    }
}

struct Args {
    bam: String,
    fasta: String,
    copies: String,
    copies_fa: String,
    index: String,
    out: String,
    delta: f64,
    max_reads: usize,
    min_cluster: usize,
    min_support: usize,
    threads: usize,
    families: Option<Vec<String>>,
    /// The command line without `--threads` and its value: the stage key's `cmd` (the thread count never changes the result).
    key_cmd: String,
}

/// Hand-rolled like `missing_copy_flag`'s, but strict: every argument is a known flag followed by its value, none given twice.
fn parse_args(raw: &[String]) -> Result<Args> {
    const FLAGS: [&str; 12] =
        ["--bam", "--fasta", "--copies", "--copies-fa", "--index", "--out", "--delta", "--max-reads", "--min-cluster", "--min-support", "--threads", "--families"];
    let usage = |msg: String| exit_two(format!("{msg}\n{USAGE}"));
    let mut kv: BTreeMap<&str, &str> = BTreeMap::new();
    let mut key_cmd: Vec<&str> = Vec::new();
    let mut i = 0;
    while i < raw.len() {
        let flag = raw[i].as_str();
        if !FLAGS.contains(&flag) {
            return Err(usage(format!("unknown argument `{flag}`")));
        }
        let value = raw.get(i + 1).ok_or_else(|| usage(format!("{flag} needs a value")))?;
        if kv.insert(flag, value).is_some() {
            return Err(usage(format!("{flag} is given twice")));
        }
        if flag != "--threads" {
            key_cmd.extend([flag, value.as_str()]);
        }
        i += 2;
    }
    let need = |k: &str| kv.get(k).map(|v| v.to_string()).ok_or_else(|| usage(format!("missing {k}")));
    fn num<T: std::str::FromStr>(kv: &BTreeMap<&str, &str>, k: &str, default: T) -> Result<T> {
        match kv.get(k) {
            None => Ok(default),
            Some(v) => v.parse().map_err(|_| exit_two(format!("{k}: `{v}` is not a valid value\n{USAGE}"))),
        }
    }
    let args = Args {
        bam: need("--bam")?,
        fasta: need("--fasta")?,
        copies: need("--copies")?,
        copies_fa: need("--copies-fa")?,
        index: need("--index")?,
        out: need("--out")?,
        delta: num(&kv, "--delta", 0.00958)?,
        max_reads: num(&kv, "--max-reads", 1000)?,
        min_cluster: num(&kv, "--min-cluster", 3)?,
        min_support: num(&kv, "--min-support", 6)?,
        threads: num(&kv, "--threads", 4)?,
        families: kv.get("--families").map(|v| v.split(',').filter(|s| !s.is_empty()).map(str::to_string).collect()),
        key_cmd: key_cmd.join(" "),
    };
    if !(args.delta.is_finite() && args.delta >= 0.0) {
        return Err(usage(format!("--delta must be a finite number >= 0, not {}", args.delta)));
    }
    for (k, v) in [("--max-reads", args.max_reads), ("--min-cluster", args.min_cluster), ("--threads", args.threads)] {
        if v == 0 {
            return Err(usage(format!("{k} must be at least 1")));
        }
    }
    Ok(args)
}

/// Run-wide settings of the minimap2 calls, and their count (reported).
struct Mm2 {
    cache: Option<PathBuf>,
    threads: usize,
    calls: Cell<usize>,
}
impl Mm2 {
    fn run(&self, args: &[&str], target: &Path, query: &Path, out: &Path) -> Result<Vec<PafHit>> {
        self.calls.set(self.calls.get() + 1);
        minimap2(args, target, query, out, self.cache.as_deref(), self.threads)?;
        read_paf(out)
    }
}

fn read_paf(path: &Path) -> Result<Vec<PafHit>> {
    Ok(parse_paf(&std::fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?))
}

/// A FASTA of `(name, sequence)` records, one line each; any old file is unlinked first (never written through a link).
fn write_fasta<'a>(path: &Path, records: impl IntoIterator<Item = (String, &'a [u8])>) -> Result<()> {
    match std::fs::remove_file(path) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
        Err(e) => return Err(e).with_context(|| format!("removing {}", path.display())),
    }
    let mut w = std::io::BufWriter::new(std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?);
    for (name, seq) in records {
        writeln!(w, ">{name}")?;
        w.write_all(seq)?;
        writeln!(w)?;
    }
    w.flush().with_context(|| format!("writing {}", path.display()))
}

fn fresh_dir(dir: &Path) -> Result<()> {
    if dir.exists() {
        std::fs::remove_dir_all(dir).with_context(|| format!("removing {}", dir.display()))?;
    }
    std::fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))
}

fn run(raw: &[String]) -> Result<()> {
    let t0 = Instant::now();
    let args = parse_args(raw)?;
    let mm2 = minimap2_binary();
    std::process::Command::new(&mm2)
        .arg("--version")
        .output()
        .map_err(|e| exit_two(format!("minimap2 could not be started: `{mm2}` ({e}); RUSTLE_MINIMAP2 names the minimap2 binary")))?;
    for (flag, path) in [("--bam", args.bam.clone()), ("--bam index", format!("{}.bai", args.bam)), ("--fasta", args.fasta.clone()), ("--index", args.index.clone())] {
        if !Path::new(&path).is_file() {
            return Err(exit_two(format!("{flag}: {path} does not exist")));
        }
    }
    let (families, kmer_input) = load_copies(&args)?;
    eprintln!("[o3_candidates] {} families from {} ({} copy sequences index the unmapped reads)", families.len(), args.copies, kmer_input.len());

    // the stage cache: one `cand` entry holds the whole result (spec §5.8)
    let cache = rc::cache_root();
    let stage = cache.as_ref().map(|root| rc::Entry::new(root, "cand", stage_key(&args, &mm2)).pinned());
    if let Some(e) = stage.as_ref() {
        if replay_stage(e, &args.out) {
            eprintln!("[cache] o3_candidates: stage result replayed from {} (BAM passes and minimap2 skipped)", e.dir.display());
            return Ok(());
        }
    }
    // a run that fails must not leave an earlier run's products in place, looking current
    for p in PRODUCTS {
        let path = format!("{}.{p}", args.out);
        match std::fs::remove_file(&path) {
            Ok(()) => {}
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
            Err(e) => return Err(e).with_context(|| format!("removing the earlier product {path}")),
        }
    }
    let tmp = PathBuf::from(format!("{}.tmp", args.out));
    fresh_dir(&tmp)?;
    let mm = Mm2 { cache, threads: args.threads, calls: Cell::new(0) };

    // spec §5.1-5.2: the nets, then the cap
    let nets = collect_nets(&args, &families, &kmer_input)?;
    let t_nets = t0.elapsed().as_secs_f64();

    // phase 1, per family: clusters, consensus, refinement, significance merge (spec §5.3-5.5)
    let alpha = AssignParams::default().alpha;
    let mut work: Vec<FamilyWork> = Vec::with_capacity(families.len());
    for (fi, fam) in families.iter().enumerate() {
        let used = sample_net(&nets.names[fi], args.max_reads);
        let dir = tmp.join(format!("f{fi}"));
        let mut w = FamilyWork { family: fam.family_id.clone(), n_net: nets.names[fi].len(), n_used: used.len(), clusters: Vec::new(), members: Vec::new() };
        if used.len() >= args.min_cluster {
            fresh_dir(&dir)?;
            let net = Net::write(&dir, used, &nets.seqs)?;
            let (clusters, log) = family_clusters(&net, &args, &dir, &mm, alpha)?;
            eprintln!("[o3_candidates] {}: net {} reads ({} used) | {}", w.family, w.n_net, w.n_used, log);
            for (k, c) in clusters.into_iter().enumerate() {
                let id = format!("{}:c{k}", w.family);
                w.members.push(c.members.iter().map(|&m| net.names[m].clone()).collect());
                w.clusters.push(ClusterSeq { family: w.family.clone(), id, n_reads: c.members.len(), seq: c.consensus });
            }
            std::fs::remove_dir_all(&dir).with_context(|| format!("removing {}", dir.display()))?;
        } else {
            eprintln!("[o3_candidates] {}: net {} reads ({} used) | fewer than --min-cluster {}: no cluster", w.family, w.n_net, w.n_used, args.min_cluster);
        }
        work.push(w);
    }
    let t_clusters = t0.elapsed().as_secs_f64();

    // phase 2: every consensus of every family against the genome, one call (spec §5.6.1); ruling R4: the hit judged is the best by identity x coverage
    let dir = tmp.join("genome");
    fresh_dir(&dir)?;
    let n_consensus: usize = work.iter().map(|w| w.clusters.len()).sum();
    let best = if n_consensus == 0 {
        HashMap::new()
    } else {
        let fa = dir.join("consensus.fa");
        write_fasta(&fa, work.iter().flat_map(|w| w.clusters.iter().map(|c| (c.id.clone(), c.seq.as_slice()))))?;
        let paf = dir.join("genome.paf");
        mm.calls.set(mm.calls.get() + 1);
        // R7: the multi-GB index is keyed by its file fingerprint, never read for the key
        minimap2_keyed(MM2_GENOME, Path::new(&args.index), &rc::file_fingerprint(&args.index), &fa, &paf, mm.cache.as_deref(), mm.threads)?;
        best_by_id_cov(&read_paf(&paf)?)
    };
    let t_genome = t0.elapsed().as_secs_f64();

    // phase 3, per family: fates, components of the new-copy clusters, flag floor (R1), exon unions (spec §5.6-5.7)
    let mut cands: Vec<Candidate> = Vec::new();
    let mut linked: Vec<(ClusterSeq, String, f64)> = Vec::new();
    let mut counts: Vec<FamilyCounts> = Vec::new();
    let mut members_out: Vec<(ClusterSeq, Vec<String>)> = Vec::new();
    let mut flagged_families: Vec<usize> = Vec::new();
    for (fi, w) in work.into_iter().enumerate() {
        let mut fc = FamilyCounts { family: w.family.clone(), n_net: w.n_net, n_used: w.n_used, n_clusters: w.clusters.len(), ..Default::default() };
        let mut new: Vec<(usize, String, f64)> = Vec::new();
        for (k, c) in w.clusters.iter().enumerate() {
            match classify(c, best.get(&c.id), args.delta) {
                Fate::InReference => fc.n_in_reference += 1,
                Fate::Linked { locus, d } => {
                    fc.n_linked += 1;
                    linked.push((c.clone(), locus, d));
                    members_out.push((c.clone(), w.members[k].clone()));
                }
                Fate::NewCopy { nearest, d } => {
                    fc.n_new += 1;
                    new.push((k, nearest, d));
                    members_out.push((c.clone(), w.members[k].clone()));
                }
            }
        }
        let fam_cands = if new.is_empty() { Vec::new() } else { family_candidates(&w, &new, &args, &tmp.join(format!("u{fi}")), &mm)? };
        fc.n_candidates = fam_cands.len();
        fc.n_flagged = fam_cands.iter().filter(|c| c.flagged).count();
        if fc.n_flagged > 0 {
            flagged_families.push(fi);
        }
        if fc.n_clusters > 0 {
            eprintln!(
                "[o3_candidates] {}: {} clusters: {} in the reference, {} linked, {} new copy -> {} candidates ({} flagged)",
                fc.family, fc.n_clusters, fc.n_in_reference, fc.n_linked, fc.n_new, fc.n_candidates, fc.n_flagged
            );
        }
        cands.extend(fam_cands);
        counts.push(fc);
    }

    // the nets of the families with a flagged candidate, whole (before the cap: the patch realignment and O2 need every read), each read once
    // (R9: the first family in --copies order keeps a read two nets share)
    let mut written: HashSet<&str> = HashSet::new();
    let nets_for_patch: Vec<(String, Vec<(String, Vec<u8>)>)> = flagged_families
        .iter()
        .map(|&fi| {
            let reads = nets.names[fi].iter().filter(|n| written.insert(n.as_str())).map(|n| (n.clone(), nets.seqs[n].clone())).collect();
            (families[fi].family_id.clone(), reads)
        })
        .collect();
    write_outputs(&args.out, &cands, &linked, &nets_for_patch)?;
    write_family_table(&args.out, &counts)?;
    let member_refs: Vec<(&ClusterSeq, &[String])> = members_out.iter().map(|(c, m)| (c, m.as_slice())).collect();
    write_cluster_members(&args.out, &member_refs)?;
    if let Some(e) = stage.as_ref() {
        store_stage(e, &args.out);
    }
    std::fs::remove_dir_all(&tmp).with_context(|| format!("removing {}", tmp.display()))?;
    let n_flagged: usize = counts.iter().map(|c| c.n_flagged).sum();
    eprintln!(
        "[o3_candidates] done: {} families, {} clusters, {} candidates, {} flagged in {} families; {} minimap2 calls; {:.1} s (nets {:.1} s, clusters {:.1} s, genome {:.1} s)",
        counts.len(), n_consensus, cands.len(), n_flagged, flagged_families.len(), mm.calls.get(), t0.elapsed().as_secs_f64(),
        t_nets, t_clusters - t_nets, t_genome - t_clusters
    );
    Ok(())
}

/// The families of `--copies` the stage runs on (`--families`, else all), and the input of the unmapped-read index: the sequence of every
/// copy of EVERY family (a read is judged against all of them, spec §5.2) from `--copies-fa`, checked against its row. `partner` rows (another
/// family's unit, `catalog_input::CatalogCopy::partner`) are no copy of the family: they bring no reads into its net and no k-mers into its index.
fn load_copies(args: &Args) -> Result<(Vec<CatalogFamily>, Vec<(String, Vec<u8>)>)> {
    let text = std::fs::read_to_string(&args.copies).map_err(|e| exit_two(format!("--copies {}: {e}", args.copies)))?;
    if text.trim().is_empty() {
        return Err(exit_two(format!("--copies {} is empty", args.copies)));
    }
    let rows = parse_copies_tsv(&text).map_err(|e| exit_two(format!("--copies {}: {e:#}", args.copies)))?;
    let all = group_families(rows)?;
    if let Some(f) = all.iter().find(|f| f.family_id.is_empty() || f.family_id.contains(char::is_whitespace)) {
        return Err(exit_two(format!("--copies {}: family id {:?} is empty or holds whitespace (it names FASTA records)", args.copies, f.family_id)));
    }
    let fa = std::fs::read_to_string(&args.copies_fa).map_err(|e| exit_two(format!("--copies-fa {}: {e}", args.copies_fa)))?;
    let seqs = parse_copies_fa(&fa).map_err(|e| exit_two(format!("--copies-fa {}: {e:#}", args.copies_fa)))?;
    let mut kmer_input = Vec::new();
    for c in all.iter().flat_map(|f| f.copies.iter()).filter(|c| !c.partner) {
        let s = seqs
            .get(&(c.family_id.clone(), c.copy_idx))
            .filter(|s| (s.chrom.as_str(), s.start, s.end) == (c.chrom.as_str(), c.start, c.end))
            .ok_or_else(|| {
                exit_two(format!(
                    "--copies-fa {}: no record {}|{}|{}:{}-{} for copy {} of family {} in --copies",
                    args.copies_fa, c.family_id, c.copy_idx, c.chrom, c.start, c.end, c.copy_idx, c.family_id
                ))
            })?;
        kmer_input.push((c.family_id.clone(), s.seq.clone()));
    }
    let families = match &args.families {
        None => all,
        Some(want) => {
            let known: HashSet<&str> = all.iter().map(|f| f.family_id.as_str()).collect();
            let unknown: Vec<&String> = want.iter().filter(|w| !known.contains(w.as_str())).collect();
            if !unknown.is_empty() {
                return Err(exit_two(format!("--families: {unknown:?} not in --copies {}", args.copies)));
            }
            all.into_iter().filter(|f| want.contains(&f.family_id)).collect()
        }
    };
    Ok((families, kmer_input))
}

/// The stage key (spec §5.8): the command line without `--threads`, the executable, the minimap2 build, the fingerprints (path, size, mtime) of
/// every input file and of the BAM and FASTA indexes, and every `RUSTLE_*` variable.
fn stage_key(args: &Args, mm2: &str) -> String {
    let fp = |p: &str| rc::file_fingerprint(p);
    format!(
        "rustle o3 candidates v1\ncmd\t{}\nexe\t{}\nminimap2\t{}\nbam\t{}\nbai\t{}\nfasta\t{}\nfai\t{}\ncopies\t{}\ncopies_fa\t{}\nindex\t{}\n{}",
        args.key_cmd,
        rc::exe_fingerprint(),
        rc::minimap2_version(mm2),
        fp(&args.bam),
        fp(&format!("{}.bai", args.bam)),
        fp(&args.fasta),
        fp(&format!("{}.fai", args.fasta)),
        fp(&args.copies),
        fp(&args.copies_fa),
        fp(&args.index),
        rc::env_fingerprint(&[])
    )
}

/// Replays every product of a complete entry to `<out>.<product>` (hard links: the products are unlinked before any rewrite); `false` (the
/// stage is computed) when the entry is no hit, lacks a product, or a product cannot be replayed.
fn replay_stage(e: &rc::Entry, out: &str) -> bool {
    e.is_hit()
        && PRODUCTS.iter().all(|p| e.dir.join(p).is_file())
        && PRODUCTS.iter().all(|p| e.replay(p, Path::new(&format!("{out}.{p}"))).is_ok())
}

fn store_stage(e: &rc::Entry, out: &str) {
    let stored = e.staging().and_then(|st| {
        for p in PRODUCTS {
            e.stage_link(&st, p, Path::new(&format!("{out}.{p}")))?;
        }
        e.commit(&st)
    });
    match stored {
        Ok(()) => eprintln!("[cache] o3_candidates: stage result stored in {}", e.dir.display()),
        Err(err) => eprintln!("[cache] could not store the stage result ({err:#}); continuing"),
    }
}

/// The nets of the families the stage runs on (spec §5.1-5.2): per family the sorted names of its reads (every read with a sequence), and
/// each read's sequence as sequenced.
struct Nets {
    names: Vec<Vec<String>>,
    seqs: HashMap<String, Vec<u8>>,
}

/// The read's sequence as sequenced: the record's, reverse-complemented when the record is reverse.
fn oriented_sequence(record: &noodles_bam::Record) -> Vec<u8> {
    let seq: Vec<u8> = record.sequence().iter().collect();
    if record.flags().is_reverse_complemented() {
        reverse_complement(&seq)
    } else {
        seq
    }
}

/// Pass A (indexed, per copy interval): every primary or secondary record overlapping a copy names its read for that family's net
/// (supplementary records do not); a primary record gives the read's sequence. Pass B (one sequential sweep of the whole BAM): the sequence
/// of every read that only secondary records named, from its primary record wherever it lies, and every unmapped record of
/// >= `MIN_UNMAPPED_LEN` bases that `FamilyKmerIndex::attribute` gives to a family the stage runs on. A read whose sequence never appears
/// (a secondary-only name without a primary record in the BAM) stays out of every net.
fn collect_nets(args: &Args, families: &[CatalogFamily], kmer_input: &[(String, Vec<u8>)]) -> Result<Nets> {
    let mut seqs: HashMap<String, Vec<u8>> = HashMap::new();
    let mut names: Vec<BTreeSet<String>> = vec![BTreeSet::new(); families.len()];
    let (mut n_primary, mut n_secondary) = (0usize, 0usize);
    {
        let mut reader = rustle::bam::open_bam(&args.bam, args.threads).with_context(|| format!("opening --bam {}", args.bam))?;
        let header = reader.read_header().with_context(|| format!("reading the header of {}", args.bam))?;
        let bai = format!("{}.bai", args.bam);
        let index = noodles_bam::bai::read(&bai).with_context(|| format!("reading the BAM index {bai}"))?;
        for (fi, fam) in families.iter().enumerate() {
            for c in fam.copies.iter().filter(|c| !c.partner) {
                anyhow::ensure!(
                    header.reference_sequences().contains_key(c.chrom.as_bytes()),
                    "copy {} of family {} lies on {}, which the header of {} does not name",
                    c.copy_idx, c.family_id, c.chrom, args.bam
                );
                // spec §5.1: the read-supported locus extent when the row has one, else the copy's span (0-based half-open -> 1-based closed)
                let (lo, hi) = c.locus.unwrap_or((c.start, c.end));
                let region = Region::new(c.chrom.as_str(), Position::try_from(lo as usize + 1)?..=Position::try_from(hi.max(lo + 1) as usize)?);
                let query = reader.query(&header, &index, &region).with_context(|| format!("querying {} at {}:{lo}-{hi}", args.bam, c.chrom))?;
                for result in query {
                    let record = result.with_context(|| format!("reading {} at {}:{lo}-{hi}", args.bam, c.chrom))?;
                    let flags = record.flags();
                    if flags.is_unmapped() || flags.is_supplementary() {
                        continue;
                    }
                    let Some(name) = record.name() else { continue };
                    let name = name.to_string();
                    if flags.is_secondary() {
                        n_secondary += 1;
                    } else if !seqs.contains_key(&name) {
                        n_primary += 1;
                        let seq = oriented_sequence(&record);
                        if !seq.is_empty() {
                            seqs.insert(name.clone(), seq);
                        }
                    }
                    names[fi].insert(name);
                }
            }
        }
    }
    let need: HashSet<String> = names.iter().flatten().filter(|n| !seqs.contains_key(*n)).cloned().collect();
    let kmer_index = FamilyKmerIndex::build(kmer_input, KMER_K, ATTRIB_MAX_FAMILIES);
    let fam_pos: HashMap<&str, usize> = families.iter().enumerate().map(|(i, f)| (f.family_id.as_str(), i)).collect();
    let (mut n_found, mut n_long_unmapped, mut n_attributed) = (0usize, 0usize, 0usize);
    let mut reader = rustle::bam::open_bam(&args.bam, args.threads).with_context(|| format!("opening --bam {}", args.bam))?;
    reader.read_header().with_context(|| format!("reading the header of {}", args.bam))?;
    for result in reader.records() {
        let record = result.with_context(|| format!("reading {} (the sequential pass)", args.bam))?;
        let flags = record.flags();
        if flags.is_unmapped() {
            // Review Focus 3: a short unmapped read is never decoded, attributed or counted
            if record.sequence().len() < MIN_UNMAPPED_LEN {
                continue;
            }
            n_long_unmapped += 1;
            let seq = oriented_sequence(&record);
            let Some(&fi) = kmer_index.attribute(&seq).and_then(|f| fam_pos.get(f.as_str())) else { continue };
            let Some(name) = record.name() else { continue };
            let name = name.to_string();
            n_attributed += 1;
            names[fi].insert(name.clone());
            seqs.entry(name).or_insert(seq);
            continue;
        }
        if need.is_empty() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let Some(name) = record.name() else { continue };
        let Ok(name) = std::str::from_utf8(name.as_ref()) else { continue };
        if need.contains(name) && !seqs.contains_key(name) {
            let seq = oriented_sequence(&record);
            if !seq.is_empty() {
                seqs.insert(name.to_string(), seq);
                n_found += 1;
            }
        }
    }
    let n_lost = need.iter().filter(|n| !seqs.contains_key(n.as_str())).count();
    eprintln!(
        "[o3_candidates] BAM: pass A {n_primary} reads by a primary record and {n_secondary} secondary records on the copies; pass B {n_found} of {} \
         secondary-only reads found by their primary record ({n_lost} without one: left out), {n_attributed} of {n_long_unmapped} unmapped reads >= \
         {MIN_UNMAPPED_LEN} bp attributed",
        need.len()
    );
    let names = names.into_iter().map(|set| set.into_iter().filter(|n| seqs.contains_key(n)).collect()).collect();
    Ok(Nets { names, seqs })
}

/// One family's used net, written to `net.fa` (sorted names) in its temporary directory: the query of every members alignment.
struct Net<'a> {
    names: Vec<String>,
    seqs: Vec<&'a [u8]>,
    fasta: PathBuf,
}
impl<'a> Net<'a> {
    fn write(dir: &Path, names: Vec<String>, all: &'a HashMap<String, Vec<u8>>) -> Result<Net<'a>> {
        let seqs: Vec<&[u8]> = names.iter().map(|n| all[n].as_slice()).collect();
        let fasta = dir.join("net.fa");
        write_fasta(&fasta, names.iter().cloned().zip(seqs.iter().copied()))?;
        Ok(Net { names, seqs, fasta })
    }
    /// The index of the longest member, the first in net order on a tie: a cluster's template (spec §5.3).
    fn longest(&self, members: &[usize]) -> usize {
        *members.iter().max_by(|&&a, &&b| self.seqs[a].len().cmp(&self.seqs[b].len()).then(b.cmp(&a))).expect("a cluster has members")
    }
}

/// A cluster of net reads: its members (indices into the net, ascending), the read it is polished on, and its consensus.
struct Cluster {
    members: Vec<usize>,
    template: usize,
    consensus: Vec<u8>,
}

/// One family after phase 1: its final clusters (named `<family>:c<k>`) with their members' names, and the net sizes.
struct FamilyWork {
    family: String,
    n_net: usize,
    n_used: usize,
    clusters: Vec<ClusterSeq>,
    members: Vec<Vec<String>>,
}

/// What phase 1 did to one family's clusters (the log line).
#[derive(Default)]
struct ClusterLog {
    read_clusters: usize,
    empty: usize,
    split_off: usize,
    split_clusters: usize,
    dropped_small: usize,
    rounds: usize,
    absorbed: usize,
    fin: usize,
}
impl std::fmt::Display for ClusterLog {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} read clusters >= --min-cluster; {} empty consensus dropped; refinement split off {} reads ({} new clusters, {} clusters fell under \
             --min-cluster); significance merge absorbed {} clusters in {} rounds; {} clusters",
            self.read_clusters, self.empty, self.split_off, self.split_clusters, self.dropped_small, self.absorbed, self.rounds, self.fin
        )
    }
}

/// Each member's best hit (most matches; the first on a tie) against its OWN group's target, from one `MM2_MEMBERS` run with the net as the
/// query and the targets, named `<tag><k>`, as the target. A member without such a hit is left out of its group's list (it abstains).
fn hits_on_own_target(net: &Net, dir: &Path, tag: &str, targets: &[&[u8]], groups: &[&[usize]], mm: &Mm2) -> Result<Vec<Vec<(usize, PafHit)>>> {
    let (fa, paf) = (dir.join(format!("{tag}.fa")), dir.join(format!("{tag}.paf")));
    write_fasta(&fa, targets.iter().enumerate().map(|(k, s)| (format!("{tag}{k}"), *s)))?;
    let hits = mm.run(MM2_MEMBERS, &fa, &net.fasta, &paf)?;
    let owner: HashMap<&str, (usize, usize)> =
        groups.iter().enumerate().flat_map(|(k, g)| g.iter().map(move |&m| (net.names[m].as_str(), (k, m)))).collect();
    let target_of: HashMap<String, usize> = (0..targets.len()).map(|k| (format!("{tag}{k}"), k)).collect();
    let mut best: HashMap<usize, PafHit> = HashMap::new();
    for h in hits {
        let Some(&(k, m)) = owner.get(h.q.as_str()) else { continue };
        if target_of.get(&h.t) != Some(&k) {
            continue;
        }
        if best.get(&m).map_or(true, |b| h.matches > b.matches) {
            best.insert(m, h);
        }
    }
    Ok(groups.iter().map(|g| g.iter().filter_map(|m| best.remove(m).map(|h| (*m, h))).collect()).collect())
}

/// The template-and-vote consensus of each `(members, template)` group (spec §5.4, rulings R2/R5): one `MM2_MEMBERS` run of the net against
/// the templates; the template read is one of the members, so its self hit covers and votes. Returns each consensus with the members' hits
/// on the template (a refined cluster is re-polished on the same template from these).
fn polish(net: &Net, dir: &Path, tag: &str, groups: &[(Vec<usize>, usize)], mm: &Mm2) -> Result<Vec<(Vec<u8>, Vec<(usize, PafHit)>)>> {
    let targets: Vec<&[u8]> = groups.iter().map(|(_, t)| net.seqs[*t]).collect();
    let members: Vec<&[usize]> = groups.iter().map(|(m, _)| m.as_slice()).collect();
    let hits = hits_on_own_target(net, dir, tag, &targets, &members, mm)?;
    groups
        .iter()
        .zip(hits)
        .map(|((_, t), h)| {
            let votes: Vec<(&[u8], &PafHit)> = h.iter().map(|(m, hit)| (net.seqs[*m], hit)).collect();
            Ok((consensus_from_template(net.seqs[*t], &votes)?, h))
        })
        .collect()
}

/// Polishes new groups (each on its longest member) and keeps those with a non-empty consensus (R5).
fn seeded_clusters(net: &Net, dir: &Path, tag: &str, groups: Vec<Vec<usize>>, mm: &Mm2, log: &mut ClusterLog) -> Result<(Vec<Cluster>, Vec<Vec<(usize, PafHit)>>)> {
    let seeds: Vec<(Vec<usize>, usize)> = groups.into_iter().map(|g| { let t = net.longest(&g); (g, t) }).collect();
    let polished = polish(net, dir, tag, &seeds, mm)?;
    let (mut clusters, mut hits) = (Vec::new(), Vec::new());
    for ((members, template), (consensus, h)) in seeds.into_iter().zip(polished) {
        if consensus.is_empty() {
            log.empty += 1;
            continue;
        }
        clusters.push(Cluster { members, template, consensus });
        hits.push(h);
    }
    Ok((clusters, hits))
}

/// Phase 1 of one family (spec §5.3-5.5 with the minimap2 engine of §9b, rulings R2/R5): the reads clustered at delta on their all-vs-all,
/// clusters under `--min-cluster` dropped, a consensus per cluster on its longest read, one refinement pass, the significance merge. The
/// final clusters come back ordered by size (descending), then first member.
fn family_clusters(net: &Net, args: &Args, dir: &Path, mm: &Mm2, alpha: f64) -> Result<(Vec<Cluster>, ClusterLog)> {
    let mut log = ClusterLog::default();
    let ava = mm.run(MM2_AVA, &net.fasta, &net.fasta, &dir.join("ava.paf"))?;
    let groups: Vec<Vec<usize>> = cluster_reads(&net.names, &ava, args.delta).into_iter().filter(|g| g.len() >= args.min_cluster).collect();
    drop(ava);
    log.read_clusters = groups.len();
    if groups.is_empty() {
        return Ok((Vec::new(), log));
    }
    let (clusters, template_hits) = seeded_clusters(net, dir, "T", groups, mm, &mut log)?;
    let clusters = if clusters.is_empty() { clusters } else { refine(net, args, dir, mm, clusters, template_hits, &mut log)? };
    let mut clusters = merge(net, dir, mm, clusters, alpha, &mut log)?;
    clusters.sort_by(|a, b| b.members.len().cmp(&a.members.len()).then(a.members[0].cmp(&b.members[0])));
    log.fin = clusters.len();
    Ok((clusters, log))
}

/// Ruling R5: one refinement pass. The members are aligned to their cluster's consensus (`MM2_MEMBERS`, the consensus as the target) and
/// `refine_cluster` keeps those that fit (`de <= delta`, half of the shorter covered); a member with no hit on the consensus does not fit
/// either. When members were split off, the kept set is re-polished on the same template (from the hits it already has) and the split-off
/// members become one new cluster, polished on its longest read, when they are at least `--min-cluster` (else they are dropped); a kept set
/// that falls under `--min-cluster` is dropped like any cluster under the floor.
fn refine(net: &Net, args: &Args, dir: &Path, mm: &Mm2, clusters: Vec<Cluster>, template_hits: Vec<Vec<(usize, PafHit)>>, log: &mut ClusterLog) -> Result<Vec<Cluster>> {
    let targets: Vec<&[u8]> = clusters.iter().map(|c| c.consensus.as_slice()).collect();
    let groups: Vec<&[usize]> = clusters.iter().map(|c| c.members.as_slice()).collect();
    let on_consensus = hits_on_own_target(net, dir, "C", &targets, &groups, mm)?;
    let (mut out, mut split): (Vec<Cluster>, Vec<Vec<usize>>) = (Vec::new(), Vec::new());
    for ((c, on_template), on_cons) in clusters.into_iter().zip(template_hits).zip(on_consensus) {
        let pairs: Vec<(&[u8], &PafHit)> = on_cons.iter().map(|(m, h)| (net.seqs[*m], h)).collect();
        let fit: HashSet<usize> = refine_cluster(&c.consensus, &pairs, args.delta).0.into_iter().map(|i| on_cons[i].0).collect();
        if fit.len() == c.members.len() {
            out.push(c);
            continue;
        }
        let (kept, rest): (Vec<usize>, Vec<usize>) = c.members.iter().copied().partition(|m| fit.contains(m));
        log.split_off += rest.len();
        if rest.len() >= args.min_cluster {
            split.push(rest);
        }
        if kept.len() < args.min_cluster {
            log.dropped_small += 1;
            continue;
        }
        let votes: Vec<(&[u8], &PafHit)> = on_template.iter().filter(|(m, _)| fit.contains(m)).map(|(m, h)| (net.seqs[*m], h)).collect();
        let consensus = consensus_from_template(net.seqs[c.template], &votes)?;
        if consensus.is_empty() {
            log.empty += 1;
            continue;
        }
        out.push(Cluster { members: kept, template: c.template, consensus });
    }
    if !split.is_empty() {
        log.split_clusters = split.len();
        out.extend(seeded_clusters(net, dir, "S", split, mm, log)?.0);
    }
    Ok(out)
}

/// Per unordered pair of `names`, the hit with the most matches (the first on a tie), whichever of the two is the query: the pair rule of
/// `cluster_reads`. Ordered by pair.
fn best_pairs<'h>(names: &[String], hits: &'h [PafHit]) -> BTreeMap<(usize, usize), &'h PafHit> {
    let idx: HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut best: BTreeMap<(usize, usize), &PafHit> = BTreeMap::new();
    for h in hits {
        let (Some(&a), Some(&b)) = (idx.get(h.q.as_str()), idx.get(h.t.as_str())) else { continue };
        if a == b {
            continue;
        }
        let key = (a.min(b), a.max(b));
        if best.get(&key).map_or(true, |o| h.matches > o.matches) {
            best.insert(key, h);
        }
    }
    best
}

/// The significance merge (spec §5.5, minimap2 engine): in each round the consensus sequences are aligned all-vs-all (`MM2_AVA`); a pair whose
/// minimizer sketches share >= 50% merges when the smaller cluster is not a real variant of the larger (`variant_is_real` with eps 0.001, the
/// assignment gate's alpha and k = the substitutions of the pair's best hit; k = 0 merges). Merges are applied by absorption: clusters in
/// order of size (descending, then index) each absorb every smaller, not yet absorbed cluster they merge with, and an absorbed cluster
/// absorbs nothing in that round, so every merge rests on a test between the absorber and the absorbed (a chain A-B-C is A+B, then (A+B)
/// against C in the next round). Each merged cluster is re-polished on the absorber's (the larger cluster's) template. Rounds repeat until
/// none merges.
fn merge(net: &Net, dir: &Path, mm: &Mm2, mut clusters: Vec<Cluster>, alpha: f64, log: &mut ClusterLog) -> Result<Vec<Cluster>> {
    while clusters.len() >= 2 {
        log.rounds += 1;
        let names: Vec<String> = (0..clusters.len()).map(|k| format!("M{k}")).collect();
        let fa = dir.join("merge.fa");
        write_fasta(&fa, names.iter().cloned().zip(clusters.iter().map(|c| c.consensus.as_slice())))?;
        let hits = mm.run(MM2_AVA, &fa, &fa, &dir.join("merge.paf"))?;
        let sketches: Vec<Vec<u64>> = clusters.iter().map(|c| minimizer_sketch(&c.consensus, KMER_K, SKETCH_W)).collect();
        let mut joins: HashSet<(usize, usize)> = HashSet::new();
        for (&(a, b), h) in &best_pairs(&names, &hits) {
            if sketch_share(&sketches[a], &sketches[b]) < MERGE_MIN_SKETCH_SHARE {
                continue;
            }
            let Some(cs) = h.cs.as_deref() else { continue };
            let k = distinguishing_columns(&parse_cs(cs).with_context(|| format!("the significance merge: {} against {}", h.q, h.t))?);
            let (na, nb) = (clusters[a].members.len(), clusters[b].members.len());
            if !variant_is_real(na.min(nb), na.max(nb), k, MERGE_EPS, alpha) {
                joins.insert((a, b));
            }
        }
        let mut order: Vec<usize> = (0..clusters.len()).collect();
        order.sort_by(|&x, &y| clusters[y].members.len().cmp(&clusters[x].members.len()).then(x.cmp(&y)));
        let mut absorbed_by: Vec<Option<usize>> = vec![None; clusters.len()];
        for (i, &x) in order.iter().enumerate() {
            if absorbed_by[x].is_some() {
                continue;
            }
            for &y in &order[i + 1..] {
                if absorbed_by[y].is_none() && joins.contains(&(x.min(y), x.max(y))) {
                    absorbed_by[y] = Some(x);
                }
            }
        }
        let n_absorbed = absorbed_by.iter().filter(|a| a.is_some()).count();
        if n_absorbed == 0 {
            break;
        }
        log.absorbed += n_absorbed;
        // the merged groups (absorber -> its members and those of every cluster it absorbed), re-polished on the absorber's template
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (y, a) in absorbed_by.iter().enumerate() {
            if let Some(x) = a {
                groups.entry(*x).or_insert_with(|| clusters[*x].members.clone()).extend(clusters[y].members.iter().copied());
            }
        }
        let seeds: Vec<(Vec<usize>, usize)> = groups
            .iter()
            .map(|(&x, m)| {
                let mut m = m.clone();
                m.sort_unstable();
                (m, clusters[x].template)
            })
            .collect();
        let polished = polish(net, dir, "R", &seeds, mm)?;
        let mut merged: HashMap<usize, Cluster> = HashMap::new();
        for ((&x, _), ((members, template), (consensus, _))) in groups.iter().zip(seeds.into_iter().zip(polished)) {
            if consensus.is_empty() {
                log.empty += 1;
                continue;
            }
            merged.insert(x, Cluster { members, template, consensus });
        }
        clusters = clusters
            .into_iter()
            .enumerate()
            .filter(|(k, _)| absorbed_by[*k].is_none())
            .filter_map(|(k, c)| if groups.contains_key(&k) { merged.remove(&k) } else { Some(c) })
            .collect();
    }
    Ok(clusters)
}

/// One member against the CURRENT union (ruling R6: `MM2_UNION`, the member as the query and the union as the target): its best hit by
/// matches, `None` without one.
fn union_hit(member: &[u8], union: &[u8], dir: &Path, mm: &Mm2) -> Result<Option<PafHit>> {
    let (u, m) = (dir.join("union.fa"), dir.join("member.fa"));
    write_fasta(&u, [("union".to_string(), union)])?;
    write_fasta(&m, [("member".to_string(), member)])?;
    let hits = mm.run(MM2_UNION, &u, &m, &dir.join("union.paf"))?;
    Ok(best_by_matches(&hits).remove("member"))
}

/// Phase 3 of one family with new-copy clusters (spec §5.6.4-5.7, ruling R1): the components of its new-copy consensus sequences (their
/// all-vs-all, `components`), each a candidate `cand_<family>_<k>` in component order, flagged when its clusters hold >= `--min-support`
/// reads, represented by the exon union of its clusters (longest first). The candidate's nearest locus and d are those of its cluster
/// closest to the reference (the smallest d; the first on a tie).
fn family_candidates(w: &FamilyWork, new: &[(usize, String, f64)], args: &Args, dir: &Path, mm: &Mm2) -> Result<Vec<Candidate>> {
    fresh_dir(dir)?;
    let comps: Vec<Vec<usize>> = if new.len() == 1 {
        vec![vec![0]]
    } else {
        let ids: Vec<String> = new.iter().map(|(k, _, _)| w.clusters[*k].id.clone()).collect();
        let fa = dir.join("new.fa");
        write_fasta(&fa, new.iter().map(|(k, _, _)| (w.clusters[*k].id.clone(), w.clusters[*k].seq.as_slice())))?;
        let hits = mm.run(MM2_AVA, &fa, &fa, &dir.join("new.paf"))?;
        components(&ids, &hits, args.delta)
    };
    let mut out = Vec::with_capacity(comps.len());
    for (k, comp) in comps.iter().enumerate() {
        let clusters: Vec<&ClusterSeq> = comp.iter().map(|&i| &w.clusters[new[i].0]).collect();
        let (nearest, d) = comp.iter().map(|&i| (&new[i].1, new[i].2)).fold(None::<(&String, f64)>, |best, x| match best {
            Some(b) if b.1 <= x.1 => Some(b),
            _ => Some(x),
        }).expect("a component has clusters");
        let mut members: Vec<Vec<u8>> = clusters.iter().map(|c| c.seq.clone()).collect();
        members.sort_by(|a, b| b.len().cmp(&a.len()));          // longest first; stable: component order on a tie
        // R6: a minimap2 failure inside the closure is kept here and returned after the union
        let mut failure: Option<anyhow::Error> = None;
        let (union, note) = union_sequence_with_note(&members, |m, u| {
            if failure.is_some() {
                return None;
            }
            union_hit(m, u, dir, mm).unwrap_or_else(|e| {
                failure = Some(e);
                None
            })
        })?;
        if let Some(e) = failure {
            return Err(e.context(format!("the exon union of {}", candidate_id(&w.family, k))));
        }
        if note.no_hit + note.skipped_minus > 0 {
            eprintln!(
                "[o3_candidates] {}: union of {} members: {} without a hit on the union, {} on the - strand (skipped)",
                candidate_id(&w.family, k), members.len(), note.no_hit, note.skipped_minus
            );
        }
        out.push(Candidate {
            family: w.family.clone(),
            id: candidate_id(&w.family, k),
            flagged: is_flagged(&clusters, args.min_support),
            clusters: clusters.into_iter().cloned().collect(),
            union,
            nearest: nearest.clone(),
            d,
            n_net: w.n_net,
            n_used: w.n_used,
        });
    }
    std::fs::remove_dir_all(dir).with_context(|| format!("removing {}", dir.display()))?;
    Ok(out)
}
