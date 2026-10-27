//! `semisec` — SEMISEC steps 1-4 (docs/LEDGER_O1.md, PREREG SEMISEC): silent copies of expressed genes found through
//! their reads' secondary alignments, classed from sequence only. The Rust port of `bench/semisec.py`'s `run`, with
//! identical output; `semisec.py report` scores it against the annotation.
//!
//! 1. Seeds: primary records (`-F 2308`) merged by overlap per strand; a seed = >= 2 reads.
//! 2. Hits: secondaries of seed reads with AS >= ratio x the read's primary AS, off the read's own seed; merged per
//!    strand into candidates of >= 2 distinct reads.
//! 3. Silent: < 2 primary records overlap the candidate.
//! 4. Class of a silent candidate from its best hit (highest AS): PROCESSED (primary spliced, hit not) or UNPROCESSED,
//!    x DISRUPTED / INTACT (the read's longest ORF, >= 100 codons, translated on the copy) or LNC-LIKE (no such ORF).
//!
//! One parallel scan of the BAM per contig (through the `.bai`), then indexed fetches of the best hits and their reads.
//! Writes OUT/seeds.tsv, OUT/cands.tsv, OUT/prim_iv.bin (+ .tsv index: primary intervals for the report's silent
//! denominator) and OUT/stats.tsv.
use anyhow::{Context, Result};
use clap::Parser;
use noodles_sam::alignment::record::cigar::op::Kind;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io::Write;

#[derive(Parser, Debug)]
#[command(about = "SEMISEC: silent copies of expressed genes from their reads' secondary alignments.")]
struct Args {
    /// Coordinate-sorted, indexed BAM (minimap2 with secondaries).
    #[arg(long)]
    bam: String,
    /// Reference FASTA with .fai.
    #[arg(long)]
    genome: String,
    /// Output directory.
    #[arg(long)]
    out: String,
    /// Keep a secondary when AS >= ratio x the read's primary AS (0.8 = minimap2's own -p default; 0 keeps every
    /// secondary the aligner reported with AS >= 0).
    #[arg(long, default_value_t = 0.8)]
    ratio: f64,
    /// Contigs scanned at once.
    #[arg(long, default_value_t = 5)]
    threads: usize,
    /// RepeatMasker .out (plain or .gz): drop a secondary hit when most of its aligned reference bases (> half) are repeat
    /// (SEMIREP). Off when absent.
    #[arg(long)]
    repeats: Option<String>,
}

const MIN_READS: usize = 2;
const MIN_CODONS: usize = 100;

struct Prim {
    name: Box<str>,
    pos: u32,
    end: u32,
    rev: bool,
    as_: i32,
    spliced: bool,
}

struct Sec {
    name: Box<str>,
    pos: u32,
    end: u32,
    rev: bool,
    as_: i32,
    /// most aligned reference bases lie in repeats (only with --repeats)
    repeat: bool,
}

/// RepeatMasker intervals per contig, 0-based half-open, merged and sorted.
fn load_repeats(path: &str) -> Result<HashMap<String, Vec<(u32, u32)>>> {
    use std::io::BufRead;
    let f = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
    let r: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(f)))
    } else {
        Box::new(std::io::BufReader::new(f))
    };
    let mut by: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    for line in r.lines() {
        let line = line?;
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 15 || !f[0].bytes().all(|b| b.is_ascii_digit()) {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[5].parse::<u32>(), f[6].parse::<u32>()) else { continue };
        by.entry(f[4].to_string()).or_default().push((s - 1, e));
    }
    for v in by.values_mut() {
        v.sort_unstable();
        let mut m: Vec<(u32, u32)> = Vec::with_capacity(v.len());
        for &(s, e) in v.iter() {
            match m.last_mut() {
                Some(l) if s <= l.1 => l.1 = l.1.max(e),
                _ => m.push((s, e)),
            }
        }
        *v = m;
    }
    Ok(by)
}

/// Bases of [s, e) inside the merged intervals `m`.
fn repeat_bp(m: &[(u32, u32)], s: u32, e: u32) -> u32 {
    let i = m.partition_point(|x| x.1 <= s);
    m[i..].iter().take_while(|x| x.0 < e).map(|x| x.1.min(e) - x.0.max(s)).sum()
}

/// Whether most aligned reference bases (M/=/X blocks) of an alignment starting at `pos` are repeat.
fn mostly_repeat(m: &[(u32, u32)], pos: u32, cig: &[(Kind, u32)]) -> bool {
    let (mut r, mut aligned, mut rep) = (pos, 0u32, 0u32);
    for &(k, n) in cig {
        match k {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                aligned += n;
                rep += repeat_bp(m, r, r + n);
                r += n;
            }
            Kind::Deletion | Kind::Skip => r += n,
            _ => {}
        }
    }
    2 * rep > aligned
}

fn ref_span(cig: &[(Kind, u32)]) -> u32 {
    cig.iter()
        .filter(|(k, _)| matches!(k, Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch))
        .map(|(_, n)| *n)
        .sum()
}

fn cigar_of(record: &noodles_bam::Record) -> Result<Vec<(Kind, u32)>> {
    let mut v = Vec::new();
    for op in record.cigar().iter() {
        let op = op?;
        v.push((op.kind(), op.len() as u32));
    }
    Ok(v)
}

/// Contigs are scanned in pieces of this many bases (by record start), so the largest contig does not bound the run.
const SPAN: u64 = 20_000_000;

/// Every primary and secondary record of one contig piece (records STARTING in it), in BAM order.
fn scan_span(bam: &str, span: &rustle::bam::ContigSpan, reps: Option<&[(u32, u32)]>) -> Result<(Vec<Prim>, Vec<Sec>)> {
    let mut reader = rustle::bam::open_bam(bam, 1)?;
    let header = reader.read_header()?;
    let (mut prims, mut secs) = (Vec::new(), Vec::new());
    let Some(id) = rustle::bam::seek_to_span(&mut reader, &header, bam, span)? else {
        return Ok((prims, secs));
    };
    let mut record = noodles_bam::Record::default();
    while reader.read_record(&mut record)? > 0 {
        match record.reference_sequence_id().transpose()? {
            Some(r) if r == id => {}
            Some(r) if r < id => continue,
            _ => break,
        }
        let start = record.alignment_start().transpose()?;
        match span.place(start.map(|p| usize::from(p) as u64 - 1)) {
            std::cmp::Ordering::Less => continue,
            std::cmp::Ordering::Greater => break,
            std::cmp::Ordering::Equal => {}
        }
        let f = record.flags();
        if f.is_unmapped() || f.is_supplementary() {
            continue;
        }
        let Some(start) = start else { continue };
        let Some(name) = record.name() else { continue };
        let pos = usize::from(start) as u32 - 1;
        let cig = cigar_of(&record)?;
        let end = pos + ref_span(&cig);
        let as_ = rustle::bam::record_as(&record).unwrap_or(0);
        let name: Box<str> = String::from_utf8_lossy(name).into();
        if f.is_secondary() {
            let repeat = reps.is_some_and(|m| mostly_repeat(m, pos, &cig));
            secs.push(Sec { name, pos, end, rev: f.is_reverse_complemented(), as_, repeat });
        } else {
            let spliced = cig.iter().any(|(k, _)| *k == Kind::Skip);
            prims.push(Prim { name, pos, end, rev: f.is_reverse_complemented(), as_, spliced });
        }
    }
    Ok((prims, secs))
}

/// Merge (start, end, payload) by strict overlap after sorting, as `semisec.py merge`.
fn merge<T: Clone>(mut ivs: Vec<(u32, u32, T)>, key: impl for<'a> Fn(&'a T) -> (&'a str, i32)) -> Vec<(u32, u32, Vec<T>)> {
    ivs.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)).then_with(|| key(&a.2).cmp(&key(&b.2))));
    let mut out: Vec<(u32, u32, Vec<T>)> = Vec::new();
    for (s, e, x) in ivs {
        match out.last_mut() {
            Some(last) if s < last.1 => {
                last.1 = last.1.max(e);
                last.2.push(x);
            }
            _ => out.push((s, e, vec![x])),
        }
    }
    out
}

/// Sorted starts, ends by start, longest span: count records overlapping [s, e).
struct Iv {
    v: Vec<(u32, u32)>,
    m: u32,
}
impl Iv {
    fn new(mut v: Vec<(u32, u32)>) -> Iv {
        v.sort_unstable();
        let m = v.iter().map(|x| x.1 - x.0).max().unwrap_or(0);
        Iv { v, m }
    }
    fn count(&self, s: u32, e: u32) -> usize {
        let i = self.v.partition_point(|x| x.0 < s.saturating_sub(self.m));
        self.v[i..].iter().take_while(|x| x.0 < e).filter(|x| x.1 > s).count()
    }
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'N' => b'N',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        b'n' => b'n',
        x => x,
    }
}

fn rc(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&b| comp(b)).collect()
}

fn is_stop(c: &[u8]) -> bool {
    matches!(c, b"TAA" | b"TAG" | b"TGA")
}

/// Longest ATG..stop ORF over the three forward frames; the first longest wins, as `semisec.py longest_orf`.
fn longest_orf(seq: &[u8]) -> (usize, usize) {
    let mut best = (0usize, 0usize);
    for fr in 0..3 {
        let mut start: Option<usize> = None;
        let mut i = fr;
        while i + 2 < seq.len() {
            let c = &seq[i..i + 3];
            match start {
                None if c == b"ATG" => start = Some(i),
                Some(st) if is_stop(c) => {
                    if i + 3 - st > best.1 - best.0 {
                        best = (st, i + 3);
                    }
                    start = None;
                }
                _ => {}
            }
            i += 3;
        }
    }
    best
}

/// Candidate bases in original-read order. Positions are doubled: a match or insertion at read position q is 2q, a
/// deletion run after read position q-1 is 2q-1 (`semisec.py`'s q - 0.5).
fn candidate_events(cig: &[(Kind, u32)], rpos: u32, reverse: bool, qlen: usize, genome: &rustle::genome::IndexedFasta, chrom: &str) -> Vec<(i64, Vec<u8>)> {
    let mut ev: Vec<(i64, Vec<u8>)> = Vec::new();
    let (mut q, mut r) = (0i64, rpos as i64);
    for &(k, n) in cig {
        let n = n as i64;
        match k {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                let rf = genome.fetch(chrom, r, r + n).unwrap_or_default();
                for j in 0..n {
                    let b = rf.get(j as usize).map(|b| b.to_ascii_uppercase());
                    ev.push((2 * (q + j), b.map(|b| vec![b]).unwrap_or_default()));
                }
                q += n;
                r += n;
            }
            Kind::Deletion => {
                let rf = genome.fetch(chrom, r, r + n).unwrap_or_default();
                ev.push((2 * q - 1, rf.iter().map(|b| b.to_ascii_uppercase()).collect()));
                r += n;
            }
            Kind::Skip => r += n,
            Kind::Insertion => {
                for j in 0..n {
                    ev.push((2 * (q + j), Vec::new()));
                }
                q += n;
            }
            Kind::SoftClip | Kind::HardClip => q += n,
            Kind::Pad => {}
        }
    }
    if reverse {
        let ql = 2 * (qlen as i64 - 1);
        ev = ev.into_iter().rev().map(|(p, b)| (ql - p, rc(&b))).collect();
    }
    ev
}

/// DISRUPTED on a frameshift or a premature in-frame stop of the read's ORF [o0, o1) on the copy.
fn orf_status(ev: &[(i64, Vec<u8>)], o0: usize, o1: usize) -> (&'static str, f64) {
    let (lo, hi) = (2 * o0 as i64, 2 * o1 as i64);
    let seg: Vec<&(i64, Vec<u8>)> = ev.iter().filter(|(p, _)| lo <= *p && *p < hi).collect();
    if seg.is_empty() {
        return ("UNCOVERED", 0.0);
    }
    let covered = seg.iter().filter(|(p, b)| !b.is_empty() && p % 2 == 0).count();
    let first = seg.iter().filter(|(p, _)| p % 2 == 0).map(|(p, _)| p / 2).min().unwrap_or(o0 as i64);
    let lead = ((first - o0 as i64) % 3) as usize;
    let cand: Vec<u8> = seg.iter().flat_map(|(_, b)| b.iter().copied()).collect();
    let net: i64 = seg.iter().filter(|(p, _)| p % 2 != 0).map(|(_, b)| b.len() as i64).sum::<i64>()
        - seg.iter().filter(|(_, b)| b.is_empty()).count() as i64;
    let frac = covered as f64 / (o1 - o0) as f64;
    if net % 3 != 0 {
        return ("DISRUPTED", frac);
    }
    let cs = &cand[((3 - lead) % 3).min(cand.len())..];
    let codons: Vec<&[u8]> = (0..cs.len().saturating_sub(2)).step_by(3).map(|i| &cs[i..i + 3]).collect();
    if codons.len() > 1 && codons[..codons.len() - 1].iter().any(|c| is_stop(c)) {
        return ("DISRUPTED", frac);
    }
    ("INTACT", frac)
}

/// The record named `name` starting at `pos` on `chrom` (secondary or primary as asked): its CIGAR and SEQ.
fn fetch_record(
    reader: &mut rustle::bam::BamFile,
    header: &noodles_sam::Header,
    index: &noodles_bam::bai::Index,
    chrom: &str,
    pos: u32,
    name: &str,
    secondary: bool,
) -> Result<Option<(Vec<(Kind, u32)>, Vec<u8>)>> {
    let region: noodles_core::Region = format!("{chrom}:{}-{}", pos + 1, pos + 1).parse()?;
    for result in reader.query(header, index, &region)? {
        let record = result?;
        let f = record.flags();
        if f.is_unmapped() || f.is_supplementary() || f.is_secondary() != secondary {
            continue;
        }
        let Some(Ok(start)) = record.alignment_start() else { continue };
        if usize::from(start) as u32 - 1 != pos || record.name().map(|n| n.as_ref() as &[u8]) != Some(name.as_bytes()) {
            continue;
        }
        return Ok(Some((cigar_of(&record)?, record.sequence().iter().collect())));
    }
    Ok(None)
}

fn main() -> Result<()> {
    let args = Args::parse();
    let t0 = std::time::Instant::now();
    std::fs::create_dir_all(&args.out)?;
    rayon::ThreadPoolBuilder::new().num_threads(args.threads.max(1)).build_global().ok();
    let contigs: Vec<String> = {
        let mut r = rustle::bam::open_bam(&args.bam, 1)?;
        r.read_header()?.reference_sequences().keys().map(|k| k.to_string()).collect()
    };
    // largest contigs first keeps the pool busy
    let counts: HashMap<String, u64> = rustle::bam::contig_record_counts(&args.bam)?
        .into_iter()
        .map(|(c, n)| (c, n.unwrap_or(u64::MAX)))
        .collect();
    let mut order: Vec<usize> = (0..contigs.len()).filter(|&i| counts.get(&contigs[i]).copied().unwrap_or(1) > 0).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(counts.get(&contigs[i]).copied().unwrap_or(0)));
    let lens: HashMap<String, u64> = {
        let mut r = rustle::bam::open_bam(&args.bam, 1)?;
        r.read_header()?
            .reference_sequences()
            .iter()
            .map(|(k, v)| (k.to_string(), usize::from(v.length()) as u64))
            .collect()
    };
    let reps = args.repeats.as_deref().map(load_repeats).transpose()?;
    let mut pieces: Vec<(usize, u64, rustle::bam::ContigSpan)> = Vec::new();
    for &i in &order {
        let (c, len) = (&contigs[i], lens[&contigs[i]]);
        let mut lo = 0u64;
        while lo < len {
            let hi = if lo + SPAN >= len { u64::MAX } else { lo + SPAN };
            pieces.push((i, lo, rustle::bam::ContigSpan { contig: c.clone(), lo, hi }));
            lo += SPAN;
        }
    }
    let mut parts: Vec<(usize, u64, Vec<Prim>, Vec<Sec>)> = pieces
        .par_iter()
        .map(|(i, lo, sp)| {
            let r = reps.as_ref().map(|h| h.get(&sp.contig).map_or(&[][..], |v| v.as_slice()));
            scan_span(&args.bam, sp, r).map(|(p, s)| (*i, *lo, p, s))
        })
        .collect::<Result<_>>()?;
    parts.sort_by_key(|x| (x.0, x.1));
    let mut scanned: Vec<(usize, Vec<Prim>, Vec<Sec>)> = Vec::new();
    for (i, _, p, s) in parts {
        match scanned.last_mut() {
            Some(last) if last.0 == i => {
                last.1.extend(p);
                last.2.extend(s);
            }
            _ => scanned.push((i, p, s)),
        }
    }
    eprintln!("[semisec] scanned {} contigs in {:.0} s", scanned.len(), t0.elapsed().as_secs_f64());

    // 1. seeds
    let mut prim: HashMap<&str, (usize, &Prim)> = HashMap::new();
    let mut read_seed: HashMap<&str, (usize, u32, u32)> = HashMap::new();
    let mut seeds: Vec<(usize, char, u32, u32, usize)> = Vec::new();
    let mut pidx: HashMap<usize, Iv> = HashMap::new();
    for (ci, prims, _) in &scanned {
        for p in prims {
            prim.insert(&p.name, (*ci, p));
        }
        for strand in [false, true] {
            let ivs: Vec<(u32, u32, &Prim)> = prims.iter().filter(|p| p.rev == strand).map(|p| (p.pos, p.end, p)).collect();
            for (s, e, ps) in merge(ivs, |p| (&*p.name, 0)) {
                let n: HashSet<&str> = ps.iter().map(|p| &*p.name).collect();
                if n.len() >= MIN_READS {
                    seeds.push((*ci, if strand { '-' } else { '+' }, s, e, n.len()));
                    for p in ps {
                        read_seed.insert(&p.name, (*ci, s, e));
                    }
                }
            }
        }
        pidx.insert(*ci, Iv::new(prims.iter().map(|p| (p.pos, p.end)).collect()));
    }
    // 2. candidate hits
    let (mut n_sec, mut n_kept, mut n_repeat) = (0u64, 0u64, 0u64);
    struct Cand<'a> {
        ci: usize,
        strand: char,
        s: u32,
        e: u32,
        reads: usize,
        nprim: usize,
        best: &'a Sec,
    }
    let mut cands: Vec<Cand> = Vec::new();
    for (ci, _, secs) in &scanned {
        let mut by: [Vec<(u32, u32, &Sec)>; 2] = [Vec::new(), Vec::new()];
        for x in secs {
            n_sec += 1;
            let Some(&(sc, ss, se)) = read_seed.get(&*x.name) else { continue };
            let p = prim[&*x.name].1;
            if (x.as_ as f64) < args.ratio * p.as_ as f64 {
                continue;
            }
            if x.repeat {
                n_repeat += 1;
                continue;
            }
            if sc == *ci && x.pos < se && x.end > ss {
                continue;
            }
            n_kept += 1;
            by[x.rev as usize].push((x.pos, x.end, x));
        }
        for (si, ivs) in by.into_iter().enumerate() {
            for (s, e, xs) in merge(ivs, |x| (&*x.name, x.as_)) {
                let reads: HashSet<&str> = xs.iter().map(|x| &*x.name).collect();
                if reads.len() < MIN_READS {
                    continue;
                }
                let mut best = xs[0];
                for x in &xs[1..] {
                    if x.as_ > best.as_ {
                        best = x;
                    }
                }
                let nprim = pidx.get(ci).map_or(0, |iv| iv.count(s, e));
                cands.push(Cand { ci: *ci, strand: if si == 1 { '-' } else { '+' }, s, e, reads: reads.len(), nprim, best });
            }
        }
    }
    cands.sort_by_key(|c| (c.ci, c.strand == '-', c.s));
    eprintln!("[semisec] {} seeds, {} candidates ({:.0} s)", seeds.len(), cands.len(), t0.elapsed().as_secs_f64());

    // 4. class of silent candidates: fetch each best hit's CIGAR and its read's primary SEQ
    let silent: Vec<usize> = (0..cands.len()).filter(|&i| cands[i].nprim < MIN_READS).collect();
    let chunk = silent.len().div_ceil(args.threads.max(1)).max(1);
    let classes: Vec<(usize, String, &'static str, f64, bool)> = silent
        .par_chunks(chunk)
        .map(|part| -> Result<Vec<_>> {
            let mut reader = rustle::bam::open_bam(&args.bam, 1)?;
            let header = reader.read_header()?;
            let index = noodles_bam::bai::read(format!("{}.bai", args.bam)).with_context(|| "reading .bai")?;
            let genome = rustle::genome::IndexedFasta::open(&args.genome)?;
            let mut out = Vec::new();
            for &i in part {
                let c = &cands[i];
                let (pci, p) = prim[&*c.best.name];
                let read = match fetch_record(&mut reader, &header, &index, &contigs[pci], p.pos, &c.best.name, false)? {
                    Some((_, seq)) if p.rev => rc(&seq),
                    Some((_, seq)) => seq,
                    None => Vec::new(),
                };
                let cig = fetch_record(&mut reader, &header, &index, &contigs[c.ci], c.best.pos, &c.best.name, true)?
                    .map(|x| x.0)
                    .with_context(|| format!("secondary {} at {}:{} not found", c.best.name, contigs[c.ci], c.best.pos))?;
                let processed = p.spliced && !cig.iter().any(|(k, _)| *k == Kind::Skip);
                let (o0, o1) = longest_orf(&read);
                let (orf, cov) = if (o1 - o0) / 3 >= MIN_CODONS {
                    let ev = candidate_events(&cig, c.best.pos, c.best.rev, read.len(), &genome, &contigs[c.ci]);
                    orf_status(&ev, o0, o1)
                } else {
                    ("LNC-LIKE", 0.0)
                };
                let class = format!("{}/{orf}", if processed { "PROCESSED" } else { "UNPROCESSED" });
                out.push((i, class, orf, cov, processed));
            }
            Ok(out)
        })
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .flatten()
        .collect();
    let cls: HashMap<usize, (String, &str, f64, bool)> = classes.into_iter().map(|(i, a, b, c, d)| (i, (a, b, c, d))).collect();

    let o = &args.out;
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{o}/cands.tsv"))?);
    writeln!(w, "chrom\tstrand\tstart\tend\treads\tprimaries\tsilent\tbest\tbest_as\tbest_seed\tclass\torf\torf_cov\tprocessed")?;
    for (i, c) in cands.iter().enumerate() {
        let (class, orf, cov, pr) = cls.get(&i).map(|x| (x.0.as_str(), x.1, format!("{:.6}", x.2), x.3.to_string())).unwrap_or(("", "", String::new(), String::new()));
        // the best read's own seed locus: the expressed gene the copy was found from
        let (sci, ss, se) = read_seed[&*c.best.name];
        writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}:{ss}-{se}\t{class}\t{orf}\t{cov}\t{pr}", contigs[c.ci], c.strand, c.s, c.e, c.reads, c.nprim, c.nprim < MIN_READS, c.best.name, c.best.as_, contigs[sci])?;
    }
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{o}/seeds.tsv"))?);
    for (ci, st, s, e, n) in &seeds {
        writeln!(w, "{}\t{st}\t{s}\t{e}\t{n}", contigs[*ci])?;
    }
    // primary intervals, sorted per contig, as little-endian u32 (start, end) pairs; the .tsv gives contig, offset, count
    let mut wb = std::io::BufWriter::new(std::fs::File::create(format!("{o}/prim_iv.bin"))?);
    let mut wi = std::io::BufWriter::new(std::fs::File::create(format!("{o}/prim_iv.tsv"))?);
    let mut off = 0usize;
    let mut cis: Vec<&usize> = pidx.keys().collect();
    cis.sort();
    for ci in cis {
        let iv = &pidx[ci];
        writeln!(wi, "{}\t{off}\t{}", contigs[*ci], iv.v.len())?;
        for &(s, e) in &iv.v {
            wb.write_all(&s.to_le_bytes())?;
            wb.write_all(&e.to_le_bytes())?;
        }
        off += iv.v.len();
    }
    let mut w = std::fs::File::create(format!("{o}/stats.tsv"))?;
    writeln!(w, "seeds\t{}\nsecondaries\t{n_sec}\nkept_hits\t{n_kept}\nrepeat_hits_dropped\t{n_repeat}\ncandidates\t{}\nsilent\t{}\nratio\t{}", seeds.len(), cands.len(), silent.len(), args.ratio)?;
    eprintln!("[semisec] {} silent classed; done in {:.0} s", silent.len(), t0.elapsed().as_secs_f64());
    Ok(())
}
