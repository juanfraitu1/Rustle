//! `missing_copy_flag` — flag expressed copies the reference does not contain, from RNA alone
//! (thesis objective O3; §6ze, `docs/PREREG_o3_rna_only_2026-09-23.md`).
//!
//! For every locus: the per-read `de` divergence mixture (S2 statistic), PSV consistency of the divergent
//! sub-pile, a spliced patched consensus, a whole-genome home search, the hypermutation / contamination /
//! editing screens, the pre-registered verdict, what DNA would have to show, and — when `--confirm` genomes
//! are given — whether the consensus has a near-perfect home there (the DNA confirmation, here the parental
//! haplotypes of the assembly's own animal).
//!
//! usage: missing_copy_flag --bam B --fasta PRIMARY.fa --loci LOCI.{gtf,gff,bed} --index PRIMARY.mmi --out PREFIX
//!        [--gff ANNOTATION.gff] [--confirm NAME=GENOME.mmi ...] [--m-min 0.10] [--delta-min 0.01]
//!        [--min-reads 10] [--min-sub 3] [--max-reads 2000] [--pi 0.002] [--threads 2] [--contigs c1,c2]
//!        [--foreign NAME=GENOME.mmi ...] [--scan-only] [--from-scan PREFIX1,PREFIX2,...]
//!        [--candidates P.cand.candidates.tsv] [--assignments P.assignments.tsv]
//!
//! Outputs `<PREFIX>.missing_copy.tsv` (one row per locus with >= --min-reads reads) and `<PREFIX>.consensus.fa`.
//! With `--candidates` (the `o3_candidates` stage's table, spec 2026-10-02 §7) the table gains a last column
//! `o3_candidate`: the flagged candidate copies whose nearest reference locus (`nearest_locus`, the consensus's best
//! genome hit) lies on the row's chromosome and overlaps its span, comma-separated, else `-` — the two O3 sources
//! corroborating each other. Without it the table is unchanged.
//! With `--assignments` (a `copy_assign` `<prefix>.assignments.tsv`) AS-tied multimappers are not trusted to
//! their primary overlap. Instead, reads with status `assigned` or `tied` and `in_copy == 1` are attributed only
//! to a locus whose label encodes the same `(family_id, copy_index)`. Loci that do not encode a clean copy key
//! fall back to primary-overlap attribution, so without `--assignments` output is byte-identical.
//!
//! Heavy work (one minimap2 run per genome index) happens once at the end, never per locus. A genome-wide run
//! on a 5-core laptop is split: `--scan-only` (BAM scan + mixture + consistency + consensus for a contig batch,
//! writes `<PREFIX>.scan.tsv` + `<PREFIX>.consensus.fa`), then one `--from-scan A,B,C` call that aligns every
//! batch's consensus sequences once per genome and writes the final table.

use anyhow::{Context, Result};
use noodles_sam::alignment::record_buf::RecordBuf;
use rustle::family::denovo_assemble::aligned_read_from_record;
use rustle::family::missing_copy::*;
use std::collections::HashMap;
use std::io::Write;

struct Args {
    bam: String,
    fasta: String,
    loci: String,
    index: String,
    out: String,
    gff: Option<String>,
    confirm: Vec<(String, String)>,
    foreign: Vec<(String, String)>,
    m_min: f64,
    delta_min: f64,
    min_reads: usize,
    min_sub: usize,
    max_reads: usize,
    pi: f64,
    threads: usize,
    contigs: Option<Vec<String>>,
    scan_only: bool,
    from_scan: Option<Vec<String>>,
    candidates: Option<String>,
    assignments: Option<String>,
}

fn parse_args() -> Result<Args> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str| -> Option<String> {
        a.iter()
            .position(|x| x == k)
            .and_then(|i| a.get(i + 1).cloned())
    };
    let need = |k: &str| -> Result<String> { get(k).with_context(|| format!("missing {k}")) };
    let mut confirm = Vec::new();
    let mut foreign = Vec::new();
    let mut i = 0;
    while i < a.len() {
        if a[i] == "--confirm" || a[i] == "--foreign" {
            let v = a.get(i + 1).context("--confirm/--foreign NAME=INDEX")?;
            let (n, p) = v
                .split_once('=')
                .context("--confirm/--foreign NAME=INDEX")?;
            if a[i] == "--confirm" {
                confirm.push((n.to_string(), p.to_string()))
            } else {
                foreign.push((n.to_string(), p.to_string()))
            }
            i += 1;
        }
        i += 1;
    }
    Ok(Args {
        bam: need("--bam")?,
        fasta: need("--fasta")?,
        loci: need("--loci")?,
        index: need("--index")?,
        out: need("--out")?,
        gff: get("--gff"),
        confirm,
        foreign,
        m_min: get("--m-min")
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(0.10),
        delta_min: get("--delta-min")
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(0.01),
        min_reads: get("--min-reads")
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(10),
        min_sub: get("--min-sub")
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(3),
        max_reads: get("--max-reads")
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(2000),
        pi: get("--pi").map(|v| v.parse()).transpose()?.unwrap_or(0.002),
        threads: get("--threads")
            .map(|v| v.parse())
            .transpose()?
            .unwrap_or(2),
        contigs: get("--contigs").map(|v| v.split(',').map(|s| s.to_string()).collect()),
        scan_only: a.iter().any(|x| x == "--scan-only"),
        from_scan: get("--from-scan").map(|v| v.split(',').map(|s| s.to_string()).collect()),
        candidates: get("--candidates"),
        assignments: get("--assignments"),
    })
}

/// A flagged `o3_candidates` candidate copy and its nearest reference locus (0-based half-open, as `nearest_locus`).
#[derive(Clone, Debug, PartialEq)]
struct FlaggedCandidate {
    id: String,
    chrom: String,
    start: u64,
    end: u64,
}

/// The flagged rows of an `o3_candidates` `<prefix>.candidates.tsv` (columns by name) that carry a nearest locus
/// (`chrom:start-end`; `none` = no genome hit, nothing to corroborate), in table order.
fn parse_flagged_candidates(text: &str) -> Result<Vec<FlaggedCandidate>> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .context("--candidates: empty file (expected an o3_candidates candidates.tsv)")?;
    let cols: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| {
        cols.iter()
            .position(|c| *c == name)
            .with_context(|| format!("--candidates: no `{name}` column (header {header:?})"))
    };
    let (i_id, i_flag, i_near) = (idx("candidate")?, idx("flagged")?, idx("nearest_locus")?);
    let mut out = Vec::new();
    for (ln, line) in lines.enumerate().filter(|(_, l)| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        let at = |i: usize| {
            f.get(i)
                .copied()
                .with_context(|| format!("--candidates line {}: too few fields", ln + 2))
        };
        if at(i_flag)? != "1" || at(i_near)? == "none" {
            continue;
        }
        let near = at(i_near)?;
        let parsed = near.rsplit_once(':').and_then(|(c, span)| {
            let (a, b) = span.split_once('-')?;
            Some((
                c.to_string(),
                a.parse::<u64>().ok()?,
                b.parse::<u64>().ok()?,
            ))
        });
        let Some((chrom, start, end)) = parsed else {
            anyhow::bail!(
                "--candidates line {}: nearest_locus {near:?} is not chrom:start-end",
                ln + 2
            );
        };
        out.push(FlaggedCandidate {
            id: at(i_id)?.to_string(),
            chrom,
            start,
            end,
        });
    }
    Ok(out)
}

/// The `o3_candidate` cell of a verdict row: the flagged candidates whose nearest locus lies on the row's chromosome and
/// overlaps its (0-based half-open) span, comma-separated in table order, else `-`.
fn o3_candidate_column(cands: &[FlaggedCandidate], r: &Row) -> String {
    let on: Vec<&str> = cands
        .iter()
        .filter(|c| c.chrom == r.chrom && c.start < r.end && r.start < c.end)
        .map(|c| c.id.as_str())
        .collect();
    if on.is_empty() {
        "-".to_string()
    } else {
        on.join(",")
    }
}

/// A PSV-resolved copy assignment for one molecule, as read from a `copy_assign` `<out>.assignments.tsv`.
#[derive(Clone, Debug, PartialEq)]
struct Assignment {
    family_id: String,
    assigned_copy: usize,
    status: String,
    in_copy: bool,
}

/// Parse `<path>.assignments.tsv` into `read_name -> Assignment`. Only rows whose `status` is `assigned`
/// or `tied` and whose `in_copy` flag is `1` are kept as usable assignments; `abstained`/`ambiguous`
/// rows are skipped entirely. The header line is required and used for column positions.
fn load_assignments(path: &str) -> Result<HashMap<String, Assignment>> {
    let text =
        std::fs::read_to_string(path).with_context(|| format!("reading --assignments {path}"))?;
    parse_assignments(&text)
}

fn parse_assignments(text: &str) -> Result<HashMap<String, Assignment>> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .context("--assignments: empty file (expected a header line)")?;
    let cols: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| {
        cols.iter()
            .position(|c| *c == name)
            .with_context(|| format!("--assignments: no `{name}` column (header {header:?})"))
    };
    let (i_read, i_fam, i_copy, i_status, i_incopy) = (
        idx("read_name")?,
        idx("family_id")?,
        idx("assigned_copy")?,
        idx("status")?,
        idx("in_copy")?,
    );
    let mut out = HashMap::new();
    for (ln, line) in lines.enumerate().filter(|(_, l)| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        let at = |i: usize| {
            f.get(i)
                .copied()
                .with_context(|| format!("--assignments line {}: too few fields", ln + 2))
        };
        let status = at(i_status)?;
        if status != "assigned" && status != "tied" {
            continue;
        }
        let in_copy = at(i_incopy)? == "1";
        let read_name = at(i_read)?.to_string();
        let family_id = at(i_fam)?.to_string();
        let assigned_copy = at(i_copy)?.parse::<usize>().with_context(|| {
            format!(
                "--assignments line {}: assigned_copy is not an integer",
                ln + 2
            )
        })?;
        out.insert(
            read_name,
            Assignment {
                family_id,
                assigned_copy,
                status: status.to_string(),
                in_copy,
            },
        );
    }
    Ok(out)
}

/// Try to extract `(family_id, copy_index)` from a locus label using an explicit `COPY` marker.
///
/// Supported conventions are case-insensitive: `FAMILY_COPY_1`, `FAMILY_0_COPY_1`,
/// `FAMILY-COPY-1`, etc. The digits immediately following the marker are taken as the copy index.
///
/// Returns `None` for labels that do not contain a clean copy marker (e.g. a plain family name).
/// Callers must fall back to primary-overlap attribution in that case.
fn locus_family_copy(label: &str) -> Option<(String, usize)> {
    let lower = label.to_lowercase();
    for sep in ["_copy_", "-copy-"] {
        if let Some(pos) = lower.find(sep) {
            let family = label[..pos].to_string();
            let rest = &label[pos + sep.len()..];
            let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
            if let Ok(n) = digits.parse::<usize>() {
                return Some((family, n));
            }
        }
    }
    None
}

fn locus_key_for_locus(id: &str, name: &str) -> Option<(String, usize)> {
    locus_family_copy(id).or_else(|| locus_family_copy(name))
}

/// Decide whether an overlapping primary read should be attributed to a locus when `--assignments` is in use.
///
/// - If the locus label does not encode `family_id/copy_index`, or no assignment file was supplied,
///   the read is attributed by primary overlap (returns `true`).
/// - If the read has no usable assignment, or its assignment is not confidently `assigned`/`tied`,
///   or its `in_copy` flag is `0`, primary overlap is also used.
/// - Otherwise the read is attributed only when the locus's `(family_id, copy_index)` exactly matches
///   the assignment.
fn read_may_attribute_to_locus(
    name: &str,
    locus_key: Option<&(String, usize)>,
    assignments: Option<&HashMap<String, Assignment>>,
) -> bool {
    let Some(lk) = locus_key else { return true };
    let Some(am) = assignments else { return true };
    let Some(a) = am.get(name) else { return true };
    if (a.status != "assigned" && a.status != "tied") || !a.in_copy {
        return true;
    }
    a.family_id == lk.0 && a.assigned_copy == lk.1
}

/// Primary reads (`-F 2308`) overlapping a region, capped by name order. With `only`, records whose name is not
/// in the set are skipped BEFORE the RecordBuf decode (the decode is the cost: a structural-only locus needs its
/// few insertion-carrying reads, not the whole pile).
fn pile(
    reader: &mut noodles_bam::io::Reader<
        noodles_bgzf::MultithreadedReader<std::io::BufReader<std::fs::File>>,
    >,
    header: &noodles_sam::Header,
    index: &noodles_bam::bai::Index,
    chrom: &str,
    lo: u64,
    hi: u64,
    cap: usize,
    only: Option<&std::collections::HashSet<&str>>,
    locus_key: Option<&(String, usize)>,
    assignments: Option<&HashMap<String, Assignment>>,
) -> Result<Vec<PileRead>> {
    let region: noodles_core::Region = format!("{chrom}:{}-{}", lo + 1, hi.max(lo + 1)).parse()?;
    let mut out = Vec::new();
    for result in reader.query(header, index, &region)? {
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        if let Some(set) = only {
            let keep = record.name().map_or(false, |n| {
                set.contains(std::str::from_utf8(n.as_ref()).unwrap_or(""))
            });
            if !keep {
                continue;
            }
        }
        let rb = RecordBuf::try_from_alignment_record(header, &record)?;
        let Some((read, _mapq, name, _as, de, _sup, _sec)) = aligned_read_from_record(&rb) else {
            continue;
        };
        if !read_may_attribute_to_locus(&name, locus_key, assignments) {
            continue;
        }
        out.push(PileRead {
            name,
            de: de as f64,
            ref_start: read.ref_start,
            ops: read.cigar,
            seq: read.seq,
        });
    }
    out.sort_by(|a, b| a.name.cmp(&b.name));
    out.truncate(cap);
    Ok(out)
}

fn minimap2_batch(index: &str, fasta: &std::path::Path, threads: usize) -> Result<String> {
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let out = std::process::Command::new(&mm2)
        .args([
            "-x",
            "splice:hq",
            "-c",
            "--eqx",
            "-N",
            "20",
            "-t",
            &threads.to_string(),
        ])
        .arg(index)
        .arg(fasta)
        .output()
        .with_context(|| format!("running {mm2}"))?;
    anyhow::ensure!(
        out.status.success(),
        "minimap2 failed on {index}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    Ok(String::from_utf8_lossy(&out.stdout).into_owned())
}

/// Everything the scan phase knows about a locus; the align phase adds home/verdict/confirmation.
#[derive(Clone, Debug)]
struct Row {
    id: String,
    name: String,
    chrom: String,
    start: u64,
    end: u64,
    n_reads: usize,
    fired: bool,
    m: f64,
    d_high: f64,
    delta: f64,
    n_sub: usize,
    n_host: usize,
    n_psv: usize,
    shared_frac: f64,
    editing_frac: f64,
    run_p: f64,
    run_top: String,
    run_top_frac: f64,
    n_runs: usize,
    is_ig_tr: bool,
    cons_len: usize,
    cons_blocks: usize,
    /// hash of the sorted sub-pile read names: overlapping loci flagged by the SAME reads share it
    sub_hash: String,
    /// fraction of the consensus blocks (the template's exons) that carry >= 1 PSV site — a real copy's
    /// PSVs spread over the transcript, an alignment artefact's cluster in one block (informational)
    psv_blocks_frac: f64,
    /// addendum 2: reads with an insertion >= 50 bp; reads whose insertion is a skipped exon (rearranged);
    /// reads whose insertion is an exon they also align (duplicated); the largest rearrangement cluster
    n_bigins: usize,
    n_rearr: usize,
    n_dup: usize,
    struct_cluster: usize,
    rearr_exon: String,
    rearr_site: u64,
    /// "fired" (mixture) / "fired_structural" / "fired_both" / "no_mixture"
    status: String,
}

const SCAN_HEAD: &str = "locus\tname\tchrom\tstart\tend\tn_reads\tstatus\tm\tde_high\tdelta\tn_sub\tn_host\tn_psv\tshared_frac\tediting_frac\trun_p\trun_top\trun_top_frac\tn_runs\tig_tr\tcons_len\tcons_blocks\tsub_hash\tpsv_blocks_frac\tn_bigins\tn_rearr\tn_dup\tstruct_cluster\trearr_exon\trearr_site";

impl Row {
    fn to_scan_line(&self) -> String {
        [
            self.id.clone(),
            self.name.clone(),
            self.chrom.clone(),
            self.start.to_string(),
            self.end.to_string(),
            self.n_reads.to_string(),
            self.status.clone(),
            format!("{:.6}", self.m),
            format!("{:.6}", self.d_high),
            format!("{:.6}", self.delta),
            self.n_sub.to_string(),
            self.n_host.to_string(),
            self.n_psv.to_string(),
            format!("{:.6}", self.shared_frac),
            format!("{:.6}", self.editing_frac),
            format!("{:.3e}", self.run_p),
            self.run_top.clone(),
            format!("{:.6}", self.run_top_frac),
            self.n_runs.to_string(),
            (self.is_ig_tr as u8).to_string(),
            self.cons_len.to_string(),
            self.cons_blocks.to_string(),
            self.sub_hash.clone(),
            format!("{:.3}", self.psv_blocks_frac),
            self.n_bigins.to_string(),
            self.n_rearr.to_string(),
            self.n_dup.to_string(),
            self.struct_cluster.to_string(),
            self.rearr_exon.clone(),
            self.rearr_site.to_string(),
        ]
        .join("\t")
    }
    fn from_scan_line(line: &str) -> Option<Row> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 30 {
            return None;
        }
        Some(Row {
            id: f[0].into(),
            name: f[1].into(),
            chrom: f[2].into(),
            start: f[3].parse().ok()?,
            end: f[4].parse().ok()?,
            n_reads: f[5].parse().ok()?,
            fired: f[6].starts_with("fired"),
            status: f[6].to_string(),
            m: f[7].parse().ok()?,
            d_high: f[8].parse().ok()?,
            delta: f[9].parse().ok()?,
            n_sub: f[10].parse().ok()?,
            n_host: f[11].parse().ok()?,
            n_psv: f[12].parse().ok()?,
            shared_frac: f[13].parse().ok()?,
            editing_frac: f[14].parse().ok()?,
            run_p: f[15].parse().ok()?,
            run_top: f[16].into(),
            run_top_frac: f[17].parse().ok()?,
            n_runs: f[18].parse().ok()?,
            is_ig_tr: f[19] == "1",
            cons_len: f[20].parse().ok()?,
            cons_blocks: f[21].parse().ok()?,
            sub_hash: f[22].into(),
            psv_blocks_frac: f[23].parse().ok()?,
            n_bigins: f[24].parse().ok()?,
            n_rearr: f[25].parse().ok()?,
            n_dup: f[26].parse().ok()?,
            struct_cluster: f[27].parse().ok()?,
            rearr_exon: f[28].into(),
            rearr_site: f[29].parse().ok()?,
        })
    }
}

fn scan(args: &Args) -> Result<Vec<Row>> {
    let loci = load_loci(&args.loci)?;
    let loci: Vec<_> = match &args.contigs {
        Some(cs) => loci.into_iter().filter(|l| cs.contains(&l.1)).collect(),
        None => loci,
    };
    eprintln!("[missing_copy_flag] {} loci", loci.len());
    let assignments: Option<HashMap<String, Assignment>> = match &args.assignments {
        Some(p) => Some(load_assignments(p)?),
        None => None,
    };
    let locus_keys: Vec<Option<(String, usize)>> = loci
        .iter()
        .map(|(id, _c, _s, _e, name)| locus_key_for_locus(id, name))
        .collect();
    if assignments.is_some() {
        let n_keyed = locus_keys.iter().filter(|k| k.is_some()).count();
        eprintln!("[missing_copy_flag] --assignments active; {n_keyed}/{} loci encode family_id/copy_index", loci.len());
    }
    let ig: rustle::types::DetHashMap<String, Vec<(u64, u64)>> = match &args.gff {
        Some(g) => load_ig_tr(g)?,
        None => Default::default(),
    };
    let contigs: rustle::types::DetHashSet<String> = loci.iter().map(|l| l.1.clone()).collect();
    let genome = rustle::genome::GenomeIndex::from_fasta_contigs(&args.fasta, &contigs)?;
    let bai = format!("{}.bai", args.bam);
    let file = std::fs::File::open(&args.bam)?;
    let worker = std::num::NonZeroUsize::new(args.threads.max(1)).unwrap();
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(
        worker,
        std::io::BufReader::with_capacity(1 << 20, file),
    );
    let mut reader = noodles_bam::io::Reader::from(bgzf);
    let header = reader.read_header()?;
    let index = noodles_bam::bai::read(&bai)?;
    let cons_path = format!("{}.consensus.fa", args.out);
    let mut cons_fa = std::fs::File::create(&cons_path)?;
    let mut rows: Vec<Row> = Vec::new();
    let (mut n_scanned, mut n_enough, mut n_fired) = (0usize, 0usize, 0usize);
    // PASS 1 (per contig, one sequential sweep): (name, de) per locus, so the mixture test never pays a
    // per-locus index query or a RecordBuf decode. PASS 2 decodes only the loci that fired.
    let mut by_contig: Vec<(String, Vec<usize>)> = Vec::new();
    for (i, l) in loci.iter().enumerate() {
        match by_contig.iter_mut().find(|(c, _)| *c == l.1) {
            Some((_, v)) => v.push(i),
            None => by_contig.push((l.1.clone(), vec![i])),
        }
    }
    for (chrom, idxs) in &by_contig {
        let Some(len) = header
            .reference_sequences()
            .get(chrom.as_bytes())
            .map(|r| usize::from(r.length()))
        else {
            continue;
        };
        let mut order: Vec<usize> = idxs.clone();
        order.sort_by_key(|&i| (loci[i].2, loci[i].3));
        let mut piles: Vec<Vec<(String, f64, u32)>> = vec![Vec::new(); order.len()];
        // insertion-carrying reads (>= 50 bp), decoded ONCE here while the record is in hand and attached to every
        // locus they overlap; the structural branch reads them from here instead of re-querying the BAM per locus
        let mut ins_reads: Vec<Vec<std::sync::Arc<PileRead>>> = vec![Vec::new(); order.len()];
        const INS_CAP: usize = 300;
        let region: noodles_core::Region = format!("{chrom}:1-{len}").parse()?;
        let (mut ptr, mut active): (usize, Vec<usize>) = (0, Vec::new());
        for result in reader.query(&header, &index, &region)? {
            let record = result?;
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let Some(start) = record.alignment_start() else {
                continue;
            };
            let rs = (usize::from(start?) as u64).saturating_sub(1);
            let mut span = 0u64;
            let mut max_ins = 0u32;
            for op in record.cigar().iter() {
                let op = op?;
                use noodles_sam::alignment::record::cigar::op::Kind;
                match op.kind() {
                    Kind::Match
                    | Kind::SequenceMatch
                    | Kind::SequenceMismatch
                    | Kind::Deletion
                    | Kind::Skip => span += op.len() as u64,
                    Kind::Insertion => max_ins = max_ins.max(op.len() as u32),
                    _ => {}
                }
            }
            let re = rs + span;
            while ptr < order.len() && loci[order[ptr]].2 < re {
                active.push(ptr);
                ptr += 1;
            }
            active.retain(|&k| loci[order[k]].3 > rs);
            if active.is_empty() {
                continue;
            }
            let de = rustle::bam::record_de(&record).map_or(0.0, |v| v as f64);
            let name = record.name().map(|n| n.to_string()).unwrap_or_default();
            let decoded: Option<std::sync::Arc<PileRead>> = if max_ins >= 50 {
                let rb = RecordBuf::try_from_alignment_record(&header, &record)?;
                aligned_read_from_record(&rb).map(|(read, _mapq, name, _as, de, _sup, _sec)| {
                    std::sync::Arc::new(PileRead {
                        name,
                        de: de as f64,
                        ref_start: read.ref_start,
                        ops: read.cigar,
                        seq: read.seq,
                    })
                })
            } else {
                None
            };
            for &k in &active {
                let l = &loci[order[k]];
                if l.2 < re && l.3 > rs {
                    if !read_may_attribute_to_locus(
                        &name,
                        locus_keys[order[k]].as_ref(),
                        assignments.as_ref(),
                    ) {
                        continue;
                    }
                    piles[k].push((name.clone(), de, max_ins));
                    if let Some(d) = &decoded {
                        if ins_reads[k].len() < INS_CAP {
                            ins_reads[k].push(d.clone());
                        }
                    }
                }
            }
        }
        let mut ins_reads = ins_reads;
        for (k, mut pile1) in piles.into_iter().enumerate() {
            let ins_here: Vec<PileRead> = std::mem::take(&mut ins_reads[k])
                .into_iter()
                .map(|a| (*a).clone())
                .collect();
            let (id, chrom, start, end, name) = &loci[order[k]];
            n_scanned += 1;
            if n_scanned % 2000 == 0 {
                eprintln!("[missing_copy_flag] {n_scanned}/{} loci, {n_enough} with reads, {n_fired} fired", loci.len());
            }
            if pile1.len() < args.min_reads {
                continue;
            }
            n_enough += 1;
            pile1.sort_by(|a, b| a.0.cmp(&b.0));
            pile1.truncate(args.max_reads);
            let de: Vec<f64> = pile1.iter().map(|x| x.1).collect();
            let fires_mixture = match two_means(&de) {
                Some((m, _, delta, mid)) => {
                    m >= args.m_min
                        && m <= 0.5
                        && delta >= args.delta_min
                        && de.iter().filter(|&&d| d > mid).count() >= args.min_sub
                }
                None => false,
            };
            let bigins: std::collections::HashSet<&str> = pile1
                .iter()
                .filter(|x| x.2 >= 50)
                .map(|x| x.0.as_str())
                .collect();
            let fires = fires_mixture || bigins.len() >= 3;
            let is_ig_tr = ig
                .get(chrom)
                .map_or(false, |v| v.iter().any(|(s, e)| *s < *end && *e > *start));
            let mut row = Row {
                id: id.clone(),
                name: name.clone(),
                chrom: chrom.clone(),
                start: *start,
                end: *end,
                n_reads: pile1.len(),
                fired: false,
                m: 0.0,
                d_high: 0.0,
                delta: 0.0,
                n_sub: 0,
                n_host: 0,
                n_psv: 0,
                shared_frac: 0.0,
                editing_frac: 0.0,
                run_p: 1.0,
                run_top: String::new(),
                run_top_frac: 0.0,
                n_runs: 0,
                is_ig_tr,
                cons_len: 0,
                cons_blocks: 0,
                sub_hash: String::new(),
                psv_blocks_frac: 0.0,
                n_bigins: bigins.len(),
                n_rearr: 0,
                n_dup: 0,
                struct_cluster: 0,
                rearr_exon: String::new(),
                rearr_site: 0,
                status: "no_mixture".to_string(),
            };
            if !fires {
                rows.push(row);
                continue;
            }
            // PASS 2: decode this locus's reads (same cap, same name order => the same pile as pass 1). A locus that
            // fires only on the structural trigger decodes just its insertion-carrying reads; the mixture branch
            // needs every read (host and sub-pile), so it decodes the whole pile.
            let reads = if fires_mixture {
                pile(
                    &mut reader,
                    &header,
                    &index,
                    chrom,
                    *start,
                    *end,
                    args.max_reads,
                    None,
                    locus_keys[order[k]].as_ref(),
                    assignments.as_ref(),
                )?
            } else {
                ins_here // structural-only: the insertion reads captured in pass 1, no second BAM query
            };
            let span_of = |r: &PileRead| {
                r.ref_start
                    + r.ops
                        .iter()
                        .filter(|(o, _)| matches!(o, '=' | 'X' | 'M' | 'D' | 'N'))
                        .map(|(_, n)| *n)
                        .sum::<u64>()
            };
            let mut wrote_consensus = false;
            let mut mixture_fired = false;
            if fires_mixture {
                if let Some(split) = split_pile(&reads, args.m_min, args.delta_min, args.min_sub) {
                    mixture_fired = true;
                    let sub: Vec<&PileRead> = split.sub.iter().map(|&i| &reads[i]).collect();
                    let host: Vec<&PileRead> = split.host.iter().map(|&i| &reads[i]).collect();
                    let lo = sub
                        .iter()
                        .map(|r| r.ref_start)
                        .min()
                        .unwrap_or(*start)
                        .min(*start);
                    let hi = sub
                        .iter()
                        .map(|r| span_of(r))
                        .max()
                        .unwrap_or(*end)
                        .max(*end);
                    let ref_seq = genome.fetch_sequence(chrom, lo, hi).unwrap_or_default();
                    let c = consistency(&sub, &host, &ref_seq, lo);
                    let (run_p, run_top, run_top_frac, n_runs) = run_screen(&sub, &host);
                    if let Some(t) = template_read(&sub) {
                        let (seq, blocks) = patched_consensus(t, &c.sites, &ref_seq, lo);
                        row.cons_len = seq.len();
                        row.cons_blocks = blocks.len();
                        if !blocks.is_empty() {
                            let with = blocks
                                .iter()
                                .filter(|(a, b)| c.sites.iter().any(|s| s.pos >= *a && s.pos < *b))
                                .count();
                            row.psv_blocks_frac = with as f64 / blocks.len() as f64;
                        }
                        if !seq.is_empty() {
                            writeln!(cons_fa, ">{id}")?;
                            cons_fa.write_all(&seq)?;
                            writeln!(cons_fa)?;
                            wrote_consensus = true;
                        }
                    }
                    let mut names: Vec<&str> = sub.iter().map(|r| r.name.as_str()).collect();
                    names.sort();
                    use std::hash::{Hash, Hasher};
                    let mut h = std::collections::hash_map::DefaultHasher::new();
                    names.hash(&mut h);
                    row.m = split.m;
                    row.d_high = split.d_high;
                    row.delta = split.delta;
                    row.n_sub = sub.len();
                    row.n_host = host.len();
                    row.n_psv = c.sites.len();
                    row.shared_frac = c.shared_frac;
                    row.editing_frac = c.editing_frac;
                    row.run_p = run_p;
                    row.run_top = run_top;
                    row.run_top_frac = run_top_frac;
                    row.n_runs = n_runs;
                    row.sub_hash = format!("{:016x}", h.finish());
                }
            }
            // STRUCTURAL branch (addendum 2): reads with an insertion >= 50 bp whose inserted sequence is an exon
            // they skip; a cluster of >= 3 sharing the same exon and insertion site fires the locus
            let mut structural_fired = false;
            if bigins.len() >= 3 {
                let big: Vec<&PileRead> = reads
                    .iter()
                    .filter(|r| bigins.contains(r.name.as_str()))
                    .collect();
                let lo = big
                    .iter()
                    .map(|r| r.ref_start)
                    .min()
                    .unwrap_or(*start)
                    .min(*start);
                let hi = big
                    .iter()
                    .map(|r| span_of(r))
                    .max()
                    .unwrap_or(*end)
                    .max(*end);
                let ref_seq = genome.fetch_sequence(chrom, lo, hi).unwrap_or_default();
                let rs = find_rearrangements(&big, &ref_seq, lo, 50);
                row.n_rearr = rs.iter().filter(|r| !r.duplicated).count();
                row.n_dup = rs.iter().filter(|r| r.duplicated).count();
                let clusters = rearrangement_clusters(&rs, 20, 3);
                if let Some(cl) = clusters.first() {
                    structural_fired = true;
                    row.struct_cluster = cl.len();
                    row.rearr_exon = format!("{}:{}-{}", chrom, cl[0].exon.0, cl[0].exon.1);
                    row.rearr_site = cl[0].ins_ref_pos;
                    let names: std::collections::HashSet<&str> =
                        cl.iter().map(|r| r.name.as_str()).collect();
                    let sub: Vec<&PileRead> = reads
                        .iter()
                        .filter(|r| names.contains(r.name.as_str()))
                        .collect();
                    let host: Vec<&PileRead> = reads
                        .iter()
                        .filter(|r| !names.contains(r.name.as_str()))
                        .collect();
                    // structural-only loci decoded just the insertion reads: n_host below counts those, not the pile;
                    // the pile size is n_reads and the run screen compares the cluster against the other decoded reads
                    if !mixture_fired {
                        let (run_p, run_top, run_top_frac, n_runs) = run_screen(&sub, &host);
                        row.run_p = run_p;
                        row.run_top = run_top;
                        row.run_top_frac = run_top_frac;
                        row.n_runs = n_runs;
                        row.n_sub = sub.len();
                        row.n_host = host.len();
                        let mut nm: Vec<&str> = sub.iter().map(|r| r.name.as_str()).collect();
                        nm.sort();
                        use std::hash::{Hash, Hasher};
                        let mut h = std::collections::hash_map::DefaultHasher::new();
                        nm.hash(&mut h);
                        row.sub_hash = format!("{:016x}", h.finish());
                    }
                    if !wrote_consensus {
                        // the copy's transcript in READ order: the longest rearranged read's own sequence
                        if let Some(t) = sub
                            .iter()
                            .max_by_key(|r| (r.seq.len(), std::cmp::Reverse(r.name.clone())))
                        {
                            if !t.seq.is_empty() {
                                writeln!(cons_fa, ">{id}")?;
                                cons_fa.write_all(&t.seq)?;
                                writeln!(cons_fa)?;
                                row.cons_len = t.seq.len();
                                row.cons_blocks = 0;
                                wrote_consensus = true;
                            }
                        }
                    }
                }
            }
            row.fired = mixture_fired || structural_fired;
            row.status = match (mixture_fired, structural_fired) {
                (true, true) => "fired_both",
                (true, false) => "fired",
                (false, true) => "fired_structural",
                (false, false) => "no_mixture",
            }
            .to_string();
            if row.fired {
                n_fired += 1;
            }
            rows.push(row);
        }
    }
    eprintln!(
        "[missing_copy_flag] scanned {n_scanned}, with >= {} reads {n_enough}, fired {n_fired}",
        args.min_reads
    );
    let scan_path = format!("{}.scan.tsv", args.out);
    let mut f = std::fs::File::create(&scan_path)?;
    writeln!(f, "{SCAN_HEAD}")?;
    for r in &rows {
        writeln!(f, "{}", r.to_scan_line())?;
    }
    eprintln!("[missing_copy_flag] wrote {scan_path} and {cons_path}");
    Ok(rows)
}

fn main() -> Result<()> {
    let args = parse_args()?;
    // read before the scan, so a bad table fails in the first second (the scan-only phase writes no verdict table)
    let cands: Option<Vec<FlaggedCandidate>> = match (&args.candidates, args.scan_only) {
        (Some(p), false) => Some(parse_flagged_candidates(
            &std::fs::read_to_string(p).with_context(|| format!("reading --candidates {p}"))?,
        )?),
        _ => None,
    };
    // phase 1: scan (this call's batch), unless resuming from earlier scans
    let (rows, cons_fasta): (Vec<Row>, std::path::PathBuf) = match &args.from_scan {
        None => {
            let rows = scan(&args)?;
            if args.scan_only {
                return Ok(());
            }
            (
                rows,
                std::path::PathBuf::from(format!("{}.consensus.fa", args.out)),
            )
        }
        Some(prefixes) => {
            use std::io::BufRead;
            let mut rows = Vec::new();
            let merged = std::path::PathBuf::from(format!("{}.consensus.fa", args.out));
            let mut out_fa = std::fs::File::create(&merged)?;
            for p in prefixes {
                let f = std::fs::File::open(format!("{p}.scan.tsv"))
                    .with_context(|| format!("{p}.scan.tsv"))?;
                for line in std::io::BufReader::new(f).lines().skip(1) {
                    if let Some(r) = Row::from_scan_line(&line?) {
                        rows.push(r);
                    }
                }
                let fa = std::fs::read(format!("{p}.consensus.fa"))
                    .with_context(|| format!("{p}.consensus.fa"))?;
                out_fa.write_all(&fa)?;
            }
            eprintln!(
                "[missing_copy_flag] resumed {} rows from {} scan batches",
                rows.len(),
                prefixes.len()
            );
            (rows, merged)
        }
    };
    let n_fired = rows.iter().filter(|r| r.fired).count();

    // phase 2: home search + confirmation, one minimap2 run per genome
    let mut hits_by: Vec<(String, HashMap<String, Vec<Hit>>)> = Vec::new();
    let mut genomes: Vec<(String, String)> = vec![("primary".to_string(), args.index.clone())];
    genomes.extend(args.confirm.iter().cloned());
    genomes.extend(args.foreign.iter().cloned());
    let n_conf = args.confirm.len();
    for (gname, gidx) in &genomes {
        let mut by: HashMap<String, Vec<Hit>> = HashMap::new();
        if n_fired > 0 {
            eprintln!(
                "[missing_copy_flag] aligning {n_fired} consensus sequences to {gname} ({gidx})"
            );
            for h in parse_paf_hits(&minimap2_batch(gidx, &cons_fasta, args.threads.max(2))?) {
                by.entry(h.query.clone()).or_default().push(h);
            }
        }
        hits_by.push((gname.clone(), by));
    }

    let tsv_path = format!("{}.missing_copy.tsv", args.out);
    let mut tsv = std::fs::File::create(&tsv_path)?;
    let mut head: Vec<String> = SCAN_HEAD.split('\t').map(String::from).collect();
    head.extend(
        [
            "class",
            "delta_over_pi",
            "consistency",
            "host_identity",
            "other_identity",
            "other_locus",
            "verdict",
            "expected_dna_depth_ratio",
        ]
        .into_iter()
        .map(String::from),
    );
    for (g, _) in hits_by.iter().skip(1).take(n_conf) {
        head.push(format!("conf_{g}_identity"));
        head.push(format!("conf_{g}_locus"));
        head.push(format!("conf_{g}"));
    }
    for (g, _) in hits_by.iter().skip(1 + n_conf) {
        head.push(format!("foreign_{g}_identity"));
        head.push(format!("foreign_{g}_locus"));
    }
    let n_extra = head.len() - SCAN_HEAD.split('\t').count();
    if cands.is_some() {
        head.push("o3_candidate".to_string());
    }
    writeln!(tsv, "{}", head.join("\t"))?;
    let mut counts: HashMap<&'static str, usize> = HashMap::new();
    let (mut n_confirmed, mut n_candidates) = (0usize, 0usize);
    let mut n_corroborated = 0usize;
    for r in &rows {
        let mut f: Vec<String> = vec![r.to_scan_line()];
        if r.fired {
            let (at, other) = home(
                hits_by[0].1.get(&r.id).map(|v| v.as_slice()).unwrap_or(&[]),
                &r.chrom,
                r.start,
                r.end,
                0.8,
            );
            let host_id = at.as_ref().map(|h| h.identity);
            let other_id = other.as_ref().map(|h| h.identity);
            let best_in = |by: &HashMap<String, Vec<Hit>>| {
                by.get(&r.id).and_then(|hs| {
                    hs.iter()
                        .filter(|h| h.qcov >= 0.8)
                        .max_by(|a, b| a.identity.partial_cmp(&b.identity).unwrap())
                        .cloned()
                })
            };
            let foreign_best: Vec<Option<Hit>> = hits_by
                .iter()
                .skip(1 + n_conf)
                .map(|(_, by)| best_in(by))
                .collect();
            let foreign_id = foreign_best
                .iter()
                .filter_map(|h| h.as_ref().map(|h| h.identity))
                .fold(None, |m: Option<f64>, x| Some(m.map_or(x, |m| m.max(x))));
            let structural_only = r.status == "fired_structural";
            let class = match r.status.as_str() {
                "fired_both" => "both",
                "fired_structural" => "structural",
                _ => "divergent",
            };
            let v = verdict(&VerdictInput {
                run_p: r.run_p,
                run_top_frac: r.run_top_frac,
                is_ig_tr: r.is_ig_tr,
                n_psv: r.n_psv,
                shared_frac: r.shared_frac,
                editing_frac: r.editing_frac,
                host_identity: host_id,
                other_identity: other_id,
                foreign_identity: foreign_id,
                delta: r.delta.max(if structural_only { 0.02 } else { 0.0 }),
                structural_only,
            });
            *counts.entry(v.as_str()).or_insert(0) += 1;
            let consistent = if r.n_psv >= 3 && r.shared_frac >= 0.5 {
                "copy_consistent"
            } else {
                "scattered"
            };
            f.extend([
                class.to_string(),
                format!("{:.1}", r.delta / args.pi),
                consistent.to_string(),
                host_id
                    .map(|x| format!("{x:.4}"))
                    .unwrap_or_else(|| "NA".into()),
                other_id
                    .map(|x| format!("{x:.4}"))
                    .unwrap_or_else(|| "NA".into()),
                other
                    .as_ref()
                    .map(|h| format!("{}:{}-{}", h.chrom, h.start, h.end))
                    .unwrap_or_else(|| "NA".into()),
                v.as_str().to_string(),
                format!("{:.2}", 1.0 / (1.0 - r.m)),
            ]);
            let is_cand = v == Verdict::ReferenceAbsentCandidate;
            n_candidates += is_cand as usize;
            let mut any_conf = false;
            for (_, by) in hits_by.iter().skip(1).take(n_conf) {
                let best = best_in(by);
                let ok = confirmed(best.as_ref().map(|h| h.identity), host_id, r.delta);
                any_conf |= ok;
                f.push(
                    best.as_ref()
                        .map(|h| format!("{:.4}", h.identity))
                        .unwrap_or_else(|| "NA".into()),
                );
                f.push(
                    best.as_ref()
                        .map(|h| format!("{}:{}-{}", h.chrom, h.start, h.end))
                        .unwrap_or_else(|| "NA".into()),
                );
                f.push((ok as u8).to_string());
            }
            for best in &foreign_best {
                f.push(
                    best.as_ref()
                        .map(|h| format!("{:.4}", h.identity))
                        .unwrap_or_else(|| "NA".into()),
                );
                f.push(
                    best.as_ref()
                        .map(|h| format!("{}:{}-{}", h.chrom, h.start, h.end))
                        .unwrap_or_else(|| "NA".into()),
                );
            }
            n_confirmed += (is_cand && any_conf) as usize;
        } else {
            f.extend(std::iter::repeat("NA".to_string()).take(n_extra));
        }
        if let Some(c) = &cands {
            let cell = o3_candidate_column(c, r);
            n_corroborated += (cell != "-") as usize;
            f.push(cell);
        }
        writeln!(tsv, "{}", f.join("\t"))?;
    }
    let mut vc: Vec<_> = counts.iter().collect();
    vc.sort();
    eprintln!("[missing_copy_flag] verdicts: {vc:?}");
    if n_conf > 0 {
        eprintln!("[missing_copy_flag] reference_absent_candidate {n_candidates}, confirmed by a confirm genome {n_confirmed}");
    }
    if let (Some(c), Some(p)) = (&cands, &args.candidates) {
        eprintln!("[missing_copy_flag] --candidates {p}: {} flagged candidates with a nearest locus; {n_corroborated} rows name one (o3_candidate)", c.len());
    }
    eprintln!("[missing_copy_flag] wrote {tsv_path}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A scan-phase row with the locus fields set and every statistic 0 (only the locus matters to `--candidates`).
    fn row(id: &str, chrom: &str, start: u64, end: u64) -> Row {
        let mut f: Vec<String> = vec![
            id.into(),
            id.into(),
            chrom.into(),
            start.to_string(),
            end.to_string(),
            "40".into(),
            "fired".into(),
        ];
        f.resize(SCAN_HEAD.split('\t').count(), "0".into());
        Row::from_scan_line(&f.join("\t")).expect("a well-formed scan line")
    }

    const CAND_HEAD: &str = "family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used";

    /// `--candidates P.cand.candidates.tsv` (spec 2026-10-02 §7): the `o3_candidate` column names the FLAGGED candidates
    /// whose nearest locus (0-based half-open, the consensus's best genome hit) lies on the row's chromosome and overlaps
    /// its span, else `-`; an unflagged candidate, one with no genome hit (`none`) and a span that only touches the
    /// locus do not corroborate.
    #[test]
    fn the_o3_candidate_column_names_a_flagged_candidate_whose_nearest_locus_overlaps_the_row() {
        let text = format!(
            "{CAND_HEAD}\n\
             MCL0\tcand_MCL0_0\t1\t60\t1\t895\tchrT:10005-13100\t0.03017\t120\t120\n\
             MCL0\tcand_MCL0_1\t1\t4\t0\t700\tchrT:10005-13100\t0.05000\t120\t120\n\
             MCL1\tcand_MCL1_0\t2\t9\t1\t600\tnone\t1.00000\t30\t30\n"
        );
        let cands = parse_flagged_candidates(&text).unwrap();
        assert_eq!(
            cands.len(),
            1,
            "only flagged candidates with a nearest locus: {cands:?}"
        );
        // the two synthetic rows: the deleted copy's surviving sibling, and a locus elsewhere on the chromosome
        assert_eq!(
            o3_candidate_column(&cands, &row("geneA", "chrT", 10000, 13100)),
            "cand_MCL0_0"
        );
        assert_eq!(
            o3_candidate_column(&cands, &row("geneB", "chrT", 20000, 21000)),
            "-"
        );
        // same coordinates on another chromosome, and a span that only touches the hit (half-open), do not corroborate
        assert_eq!(
            o3_candidate_column(&cands, &row("geneC", "chrU", 10000, 13100)),
            "-"
        );
        assert_eq!(
            o3_candidate_column(&cands, &row("geneD", "chrT", 13100, 14000)),
            "-"
        );
        // two flagged candidates on one locus: both, in table order
        let two =
            format!("{text}MCL2\tcand_MCL2_0\t1\t7\t1\t500\tchrT:12000-12500\t0.02000\t10\t10\n");
        let cands = parse_flagged_candidates(&two).unwrap();
        assert_eq!(
            o3_candidate_column(&cands, &row("geneA", "chrT", 10000, 13100)),
            "cand_MCL0_0,cand_MCL2_0"
        );
    }

    /// The table is read by column NAME, and a malformed one is an error, not an empty corroboration.
    #[test]
    fn a_candidates_table_without_its_columns_or_with_a_bad_locus_is_an_error() {
        let e = parse_flagged_candidates("family\tcandidate\tflagged\nMCL0\tcand_MCL0_0\t1\n")
            .unwrap_err()
            .to_string();
        assert!(e.contains("nearest_locus"), "{e}");
        let bad =
            format!("{CAND_HEAD}\nMCL0\tcand_MCL0_0\t1\t60\t1\t895\tchrT:10005\t0.03\t120\t120\n");
        let e = parse_flagged_candidates(&bad).unwrap_err().to_string();
        assert!(e.contains("chrT:10005"), "{e}");
        assert!(
            parse_flagged_candidates("").is_err(),
            "an empty file has no header"
        );
    }

    const ASSIGN_HEAD: &str = "read_name\tfamily_id\tassigned_copy\tstatus\tn_decisive\tmargin\tp_value\tmin_p_value\tas_best\tas_second\tas_margin\tas_per_base_best\tas_per_base_2nd\tin_copy\tcatalog_copy_idx\torigin_rejected\tn_candidates\tsole_candidate\tcontested\treadthrough_into\tprimary_local";

    /// `load_assignments` keeps only assigned/tied rows, requires the expected columns, and parses in_copy.
    #[test]
    fn assignments_parser_keeps_assigned_and_tied_and_skips_abstained_and_ambiguous() {
        let text = format!(
            "{ASSIGN_HEAD}\n\
             r1\tMCL0\t0\tassigned\t5\t1.0\t1e-5\t1e-5\t100\t80\t20\t0.1\t0.08\t1\t0\t0\t1\t1\t0\t-\t1\n\
             r2\tMCL0\t1\ttied\t2\t0.0\t1e-2\t1e-2\t90\t90\t0\t0.09\t0.09\t1\t1\t0\t2\t0\t1\t-\t1\n\
             r3\tMCL0\t0\tabstained\t0\t0.0\t1.0\t1.0\t0\t0\t0\t0.0\t0.0\t0\tNA\t0\t0\t0\t0\t-\t0\n\
             r4\tMCL0\t1\tambiguous\t1\t0.0\t1e-1\t1e-1\t70\t70\t0\t0.07\t0.07\t1\t1\t0\t2\t0\t1\t-\t1\n"
        );
        let am = parse_assignments(&text).unwrap();
        assert_eq!(am.len(), 2, "only assigned and tied are kept");
        let a1 = am.get("r1").unwrap();
        assert_eq!(a1.family_id, "MCL0");
        assert_eq!(a1.assigned_copy, 0);
        assert_eq!(a1.status, "assigned");
        assert!(a1.in_copy);
        let a2 = am.get("r2").unwrap();
        assert_eq!(a2.assigned_copy, 1);
        assert_eq!(a2.status, "tied");
        assert!(a2.in_copy);
        assert!(parse_assignments("").is_err(), "empty file has no header");
    }

    /// Locus labels with an explicit `COPY` marker parse to (family, copy); plain labels do not.
    #[test]
    fn locus_family_copy_parses_explicit_copy_markers_only() {
        assert_eq!(
            locus_family_copy("MCL0_COPY_1"),
            Some(("MCL0".to_string(), 1))
        );
        assert_eq!(
            locus_family_copy("MCL0_copy_1"),
            Some(("MCL0".to_string(), 1))
        );
        assert_eq!(
            locus_family_copy("FAM_0_COPY_12"),
            Some(("FAM_0".to_string(), 12))
        );
        assert_eq!(
            locus_family_copy("FAM-0-COPY-3"),
            Some(("FAM-0".to_string(), 3))
        );
        assert_eq!(locus_family_copy("MCL0"), None);
        assert_eq!(locus_family_copy("MCL0_1"), None);
        assert_eq!(locus_family_copy("geneA"), None);
        // the name/id fallback
        assert_eq!(
            locus_key_for_locus("geneA", "MCL0_COPY_2"),
            Some(("MCL0".to_string(), 2))
        );
    }

    /// Assignment-aware attribution: a confidently assigned read only attributes to its own copy locus;
    /// unassigned reads and loci without a clean copy key fall back to primary overlap.
    #[test]
    fn attribution_respects_assignment_when_locus_encodes_copy_and_falls_back_otherwise() {
        let mut am = HashMap::new();
        am.insert(
            "r1".to_string(),
            Assignment {
                family_id: "MCL0".to_string(),
                assigned_copy: 1,
                status: "assigned".to_string(),
                in_copy: true,
            },
        );
        let lk_copy1 = Some(("MCL0".to_string(), 1));
        let lk_copy2 = Some(("MCL0".to_string(), 2));
        let lk_none: Option<&(String, usize)> = None;
        let am_ref = Some(&am);

        // assigned read matches its copy
        assert!(read_may_attribute_to_locus("r1", lk_copy1.as_ref(), am_ref));
        // assigned read does NOT match a different copy
        assert!(!read_may_attribute_to_locus(
            "r1",
            lk_copy2.as_ref(),
            am_ref
        ));
        // unassigned read keeps primary-overlap behaviour regardless of locus key
        assert!(read_may_attribute_to_locus("rX", lk_copy1.as_ref(), am_ref));
        // locus with no clean copy key falls back to primary overlap
        assert!(read_may_attribute_to_locus("r1", lk_none, am_ref));
    }
}
