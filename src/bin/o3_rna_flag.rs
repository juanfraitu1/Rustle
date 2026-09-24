//! `o3_rna_flag` — O3 as far as RNA can carry it (§6ze, `docs/PREREG_o3_rna_only_2026-09-23.md`).
//!
//! For every locus: the per-read `de` divergence mixture (S2 statistic), PSV consistency of the divergent
//! sub-pile, a spliced patched consensus, a whole-genome home search, the hypermutation / contamination /
//! editing screens, the pre-registered verdict, what DNA would have to show, and — when `--confirm` genomes
//! are given — whether the consensus has a near-perfect home there (the DNA confirmation, here the parental
//! haplotypes of the assembly's own animal).
//!
//! usage: o3_rna_flag --bam B --fasta PRIMARY.fa --loci LOCI.{gtf,gff,bed} --index PRIMARY.mmi --out PREFIX
//!        [--gff ANNOTATION.gff] [--confirm NAME=GENOME.mmi ...] [--m-min 0.10] [--delta-min 0.01]
//!        [--min-reads 10] [--min-sub 3] [--max-reads 2000] [--pi 0.002] [--threads 2] [--contigs c1,c2]
//!        [--foreign NAME=GENOME.mmi ...] [--scan-only] [--from-scan PREFIX1,PREFIX2,...]
//!
//! Outputs `<PREFIX>.o3_rna.tsv` (one row per locus with >= --min-reads reads) and `<PREFIX>.consensus.fa`.
//! Heavy work (one minimap2 run per genome index) happens once at the end, never per locus. A genome-wide run
//! on a 5-core laptop is split: `--scan-only` (BAM scan + mixture + consistency + consensus for a contig batch,
//! writes `<PREFIX>.scan.tsv` + `<PREFIX>.consensus.fa`), then one `--from-scan A,B,C` call that aligns every
//! batch's consensus sequences once per genome and writes the final table.

use anyhow::{Context, Result};
use noodles_sam::alignment::record_buf::RecordBuf;
use rustle::vg_family::denovo_assemble::aligned_read_from_record;
use rustle::vg_family::o3_rna::*;
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
}

fn parse_args() -> Result<Args> {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str| -> Option<String> { a.iter().position(|x| x == k).and_then(|i| a.get(i + 1).cloned()) };
    let need = |k: &str| -> Result<String> { get(k).with_context(|| format!("missing {k}")) };
    let mut confirm = Vec::new();
    let mut foreign = Vec::new();
    let mut i = 0;
    while i < a.len() {
        if a[i] == "--confirm" || a[i] == "--foreign" {
            let v = a.get(i + 1).context("--confirm/--foreign NAME=INDEX")?;
            let (n, p) = v.split_once('=').context("--confirm/--foreign NAME=INDEX")?;
            if a[i] == "--confirm" { confirm.push((n.to_string(), p.to_string())) } else { foreign.push((n.to_string(), p.to_string())) }
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
        m_min: get("--m-min").map(|v| v.parse()).transpose()?.unwrap_or(0.10),
        delta_min: get("--delta-min").map(|v| v.parse()).transpose()?.unwrap_or(0.01),
        min_reads: get("--min-reads").map(|v| v.parse()).transpose()?.unwrap_or(10),
        min_sub: get("--min-sub").map(|v| v.parse()).transpose()?.unwrap_or(3),
        max_reads: get("--max-reads").map(|v| v.parse()).transpose()?.unwrap_or(2000),
        pi: get("--pi").map(|v| v.parse()).transpose()?.unwrap_or(0.002),
        threads: get("--threads").map(|v| v.parse()).transpose()?.unwrap_or(2),
        contigs: get("--contigs").map(|v| v.split(',').map(|s| s.to_string()).collect()),
        scan_only: a.iter().any(|x| x == "--scan-only"),
        from_scan: get("--from-scan").map(|v| v.split(',').map(|s| s.to_string()).collect()),
    })
}

/// Primary reads (`-F 2308`) overlapping a region, capped by name order.
fn pile(reader: &mut noodles_bam::io::Reader<noodles_bgzf::MultithreadedReader<std::io::BufReader<std::fs::File>>>, header: &noodles_sam::Header, index: &noodles_bam::bai::Index, chrom: &str, lo: u64, hi: u64, cap: usize) -> Result<Vec<PileRead>> {
    let region: noodles_core::Region = format!("{chrom}:{}-{}", lo + 1, hi.max(lo + 1)).parse()?;
    let mut out = Vec::new();
    for result in reader.query(header, index, &region)? {
        let record = result?;
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let rb = RecordBuf::try_from_alignment_record(header, &record)?;
        let Some((read, _mapq, name, _as, de, _sup, _sec)) = aligned_read_from_record(&rb) else { continue };
        out.push(PileRead { name, de: de as f64, ref_start: read.ref_start, ops: read.cigar, seq: read.seq });
    }
    out.sort_by(|a, b| a.name.cmp(&b.name));
    out.truncate(cap);
    Ok(out)
}

fn minimap2_batch(index: &str, fasta: &std::path::Path, threads: usize) -> Result<String> {
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let out = std::process::Command::new(&mm2)
        .args(["-x", "splice:hq", "-c", "--eqx", "-N", "20", "-t", &threads.to_string()])
        .arg(index)
        .arg(fasta)
        .output()
        .with_context(|| format!("running {mm2}"))?;
    anyhow::ensure!(out.status.success(), "minimap2 failed on {index}: {}", String::from_utf8_lossy(&out.stderr));
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
}

const SCAN_HEAD: &str = "locus\tname\tchrom\tstart\tend\tn_reads\tstatus\tm\tde_high\tdelta\tn_sub\tn_host\tn_psv\tshared_frac\tediting_frac\trun_p\trun_top\trun_top_frac\tn_runs\tig_tr\tcons_len\tcons_blocks\tsub_hash\tpsv_blocks_frac";

impl Row {
    fn to_scan_line(&self) -> String {
        [
            self.id.clone(), self.name.clone(), self.chrom.clone(), self.start.to_string(), self.end.to_string(), self.n_reads.to_string(),
            (if self.fired { "fired" } else { "no_mixture" }).to_string(),
            format!("{:.6}", self.m), format!("{:.6}", self.d_high), format!("{:.6}", self.delta), self.n_sub.to_string(), self.n_host.to_string(),
            self.n_psv.to_string(), format!("{:.6}", self.shared_frac), format!("{:.6}", self.editing_frac), format!("{:.3e}", self.run_p), self.run_top.clone(),
            format!("{:.6}", self.run_top_frac), self.n_runs.to_string(), (self.is_ig_tr as u8).to_string(), self.cons_len.to_string(), self.cons_blocks.to_string(), self.sub_hash.clone(),
            format!("{:.3}", self.psv_blocks_frac),
        ]
        .join("\t")
    }
    fn from_scan_line(line: &str) -> Option<Row> {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 24 {
            return None;
        }
        Some(Row {
            id: f[0].into(), name: f[1].into(), chrom: f[2].into(), start: f[3].parse().ok()?, end: f[4].parse().ok()?, n_reads: f[5].parse().ok()?,
            fired: f[6] == "fired", m: f[7].parse().ok()?, d_high: f[8].parse().ok()?, delta: f[9].parse().ok()?, n_sub: f[10].parse().ok()?, n_host: f[11].parse().ok()?,
            n_psv: f[12].parse().ok()?, shared_frac: f[13].parse().ok()?, editing_frac: f[14].parse().ok()?, run_p: f[15].parse().ok()?, run_top: f[16].into(),
            run_top_frac: f[17].parse().ok()?, n_runs: f[18].parse().ok()?, is_ig_tr: f[19] == "1", cons_len: f[20].parse().ok()?, cons_blocks: f[21].parse().ok()?, sub_hash: f[22].into(),
            psv_blocks_frac: f[23].parse().ok()?,
        })
    }
}

fn scan(args: &Args) -> Result<Vec<Row>> {
    let loci = load_loci(&args.loci)?;
    let loci: Vec<_> = match &args.contigs {
        Some(cs) => loci.into_iter().filter(|l| cs.contains(&l.1)).collect(),
        None => loci,
    };
    eprintln!("[o3_rna_flag] {} loci", loci.len());
    let ig: HashMap<String, Vec<(u64, u64)>> = match &args.gff {
        Some(g) => load_ig_tr(g)?,
        None => HashMap::new(),
    };
    let contigs: std::collections::HashSet<String> = loci.iter().map(|l| l.1.clone()).collect();
    let genome = rustle::genome::GenomeIndex::from_fasta_contigs(&args.fasta, &contigs)?;
    let bai = format!("{}.bai", args.bam);
    let file = std::fs::File::open(&args.bam)?;
    let worker = std::num::NonZeroUsize::new(args.threads.max(1)).unwrap();
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(worker, std::io::BufReader::with_capacity(1 << 20, file));
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
        let Some(len) = header.reference_sequences().get(chrom.as_bytes()).map(|r| usize::from(r.length())) else { continue };
        let mut order: Vec<usize> = idxs.clone();
        order.sort_by_key(|&i| (loci[i].2, loci[i].3));
        let mut piles: Vec<Vec<(String, f64)>> = vec![Vec::new(); order.len()];
        let region: noodles_core::Region = format!("{chrom}:1-{len}").parse()?;
        let de_tag = noodles_sam::alignment::record::data::field::Tag::new(b'd', b'e');
        let (mut ptr, mut active): (usize, Vec<usize>) = (0, Vec::new());
        for result in reader.query(&header, &index, &region)? {
            let record = result?;
            let flags = record.flags();
            if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
                continue;
            }
            let Some(start) = record.alignment_start() else { continue };
            let rs = (usize::from(start?) as u64).saturating_sub(1);
            let mut span = 0u64;
            for op in record.cigar().iter() {
                let op = op?;
                use noodles_sam::alignment::record::cigar::op::Kind;
                if matches!(op.kind(), Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion | Kind::Skip) {
                    span += op.len() as u64;
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
            let de = match record.data().get(&de_tag) {
                Some(Ok(noodles_sam::alignment::record::data::field::Value::Float(v))) => v as f64,
                _ => 0.0,
            };
            let name = record.name().map(|n| n.to_string()).unwrap_or_default();
            for &k in &active {
                let l = &loci[order[k]];
                if l.2 < re && l.3 > rs {
                    piles[k].push((name.clone(), de));
                }
            }
        }
        for (k, mut pile1) in piles.into_iter().enumerate() {
            let (id, chrom, start, end, name) = &loci[order[k]];
            n_scanned += 1;
            if n_scanned % 2000 == 0 {
                eprintln!("[o3_rna_flag] {n_scanned}/{} loci, {n_enough} with reads, {n_fired} fired", loci.len());
            }
            if pile1.len() < args.min_reads {
                continue;
            }
            n_enough += 1;
            pile1.sort_by(|a, b| a.0.cmp(&b.0));
            pile1.truncate(args.max_reads);
            let de: Vec<f64> = pile1.iter().map(|x| x.1).collect();
            let fires = match two_means(&de) {
                Some((m, _, delta, mid)) => m >= args.m_min && m <= 0.5 && delta >= args.delta_min && de.iter().filter(|&&d| d > mid).count() >= args.min_sub,
                None => false,
            };
            let is_ig_tr = ig.get(chrom).map_or(false, |v| v.iter().any(|(s, e)| *s < *end && *e > *start));
            let mut row = Row { id: id.clone(), name: name.clone(), chrom: chrom.clone(), start: *start, end: *end, n_reads: pile1.len(), fired: false, m: 0.0, d_high: 0.0, delta: 0.0, n_sub: 0, n_host: 0, n_psv: 0, shared_frac: 0.0, editing_frac: 0.0, run_p: 1.0, run_top: String::new(), run_top_frac: 0.0, n_runs: 0, is_ig_tr, cons_len: 0, cons_blocks: 0, sub_hash: String::new(), psv_blocks_frac: 0.0 };
            if !fires {
                rows.push(row);
                continue;
            }
            // PASS 2: decode this locus's reads (same cap, same name order => the same pile as pass 1)
            let reads = pile(&mut reader, &header, &index, chrom, *start, *end, args.max_reads)?;
            if let Some(split) = split_pile(&reads, args.m_min, args.delta_min, args.min_sub) {
                n_fired += 1;
                let sub: Vec<&PileRead> = split.sub.iter().map(|&i| &reads[i]).collect();
                let host: Vec<&PileRead> = split.host.iter().map(|&i| &reads[i]).collect();
                let lo = sub.iter().map(|r| r.ref_start).min().unwrap_or(*start).min(*start);
                let hi = sub.iter().map(|r| r.ref_start + r.ops.iter().filter(|(o, _)| matches!(o, '=' | 'X' | 'M' | 'D' | 'N')).map(|(_, n)| *n).sum::<u64>()).max().unwrap_or(*end).max(*end);
                let ref_seq = genome.fetch_sequence(chrom, lo, hi).unwrap_or_default();
                let c = consistency(&sub, &host, &ref_seq, lo);
                let (run_p, run_top, run_top_frac, n_runs) = run_screen(&sub, &host);
                if let Some(t) = template_read(&sub) {
                    let (seq, blocks) = patched_consensus(t, &c.sites, &ref_seq, lo);
                    row.cons_len = seq.len();
                    row.cons_blocks = blocks.len();
                    if !blocks.is_empty() {
                        let with = blocks.iter().filter(|(a, b)| c.sites.iter().any(|s| s.pos >= *a && s.pos < *b)).count();
                        row.psv_blocks_frac = with as f64 / blocks.len() as f64;
                    }
                    if !seq.is_empty() {
                        writeln!(cons_fa, ">{id}")?;
                        cons_fa.write_all(&seq)?;
                        writeln!(cons_fa)?;
                    }
                }
                let mut names: Vec<&str> = sub.iter().map(|r| r.name.as_str()).collect();
                names.sort();
                use std::hash::{Hash, Hasher};
                let mut h = std::collections::hash_map::DefaultHasher::new();
                names.hash(&mut h);
                row.fired = true;
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
            rows.push(row);
        }
    }
    eprintln!("[o3_rna_flag] scanned {n_scanned}, with >= {} reads {n_enough}, fired {n_fired}", args.min_reads);
    let scan_path = format!("{}.scan.tsv", args.out);
    let mut f = std::fs::File::create(&scan_path)?;
    writeln!(f, "{SCAN_HEAD}")?;
    for r in &rows {
        writeln!(f, "{}", r.to_scan_line())?;
    }
    eprintln!("[o3_rna_flag] wrote {scan_path} and {cons_path}");
    Ok(rows)
}

fn main() -> Result<()> {
    let args = parse_args()?;
    // phase 1: scan (this call's batch), unless resuming from earlier scans
    let (rows, cons_fasta): (Vec<Row>, std::path::PathBuf) = match &args.from_scan {
        None => {
            let rows = scan(&args)?;
            if args.scan_only {
                return Ok(());
            }
            (rows, std::path::PathBuf::from(format!("{}.consensus.fa", args.out)))
        }
        Some(prefixes) => {
            use std::io::BufRead;
            let mut rows = Vec::new();
            let merged = std::path::PathBuf::from(format!("{}.consensus.fa", args.out));
            let mut out_fa = std::fs::File::create(&merged)?;
            for p in prefixes {
                let f = std::fs::File::open(format!("{p}.scan.tsv")).with_context(|| format!("{p}.scan.tsv"))?;
                for line in std::io::BufReader::new(f).lines().skip(1) {
                    if let Some(r) = Row::from_scan_line(&line?) {
                        rows.push(r);
                    }
                }
                let fa = std::fs::read(format!("{p}.consensus.fa")).with_context(|| format!("{p}.consensus.fa"))?;
                out_fa.write_all(&fa)?;
            }
            eprintln!("[o3_rna_flag] resumed {} rows from {} scan batches", rows.len(), prefixes.len());
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
            eprintln!("[o3_rna_flag] aligning {n_fired} consensus sequences to {gname} ({gidx})");
            for h in parse_paf_hits(&minimap2_batch(gidx, &cons_fasta, args.threads.max(2))?) {
                by.entry(h.query.clone()).or_default().push(h);
            }
        }
        hits_by.push((gname.clone(), by));
    }

    let tsv_path = format!("{}.o3_rna.tsv", args.out);
    let mut tsv = std::fs::File::create(&tsv_path)?;
    let mut head: Vec<String> = SCAN_HEAD.split('\t').map(String::from).collect();
    head.extend(["delta_over_pi", "consistency", "host_identity", "other_identity", "other_locus", "verdict", "expected_dna_depth_ratio"].into_iter().map(String::from));
    for (g, _) in hits_by.iter().skip(1).take(n_conf) {
        head.push(format!("conf_{g}_identity"));
        head.push(format!("conf_{g}_locus"));
        head.push(format!("conf_{g}"));
    }
    for (g, _) in hits_by.iter().skip(1 + n_conf) {
        head.push(format!("foreign_{g}_identity"));
        head.push(format!("foreign_{g}_locus"));
    }
    writeln!(tsv, "{}", head.join("\t"))?;
    let mut counts: HashMap<&'static str, usize> = HashMap::new();
    let (mut n_confirmed, mut n_candidates) = (0usize, 0usize);
    let n_extra = head.len() - SCAN_HEAD.split('\t').count();
    for r in &rows {
        let mut f: Vec<String> = vec![r.to_scan_line()];
        if r.fired {
            let (at, other) = home(hits_by[0].1.get(&r.id).map(|v| v.as_slice()).unwrap_or(&[]), &r.chrom, r.start, r.end, 0.8);
            let host_id = at.as_ref().map(|h| h.identity);
            let other_id = other.as_ref().map(|h| h.identity);
            let best_in = |by: &HashMap<String, Vec<Hit>>| by.get(&r.id).and_then(|hs| hs.iter().filter(|h| h.qcov >= 0.8).max_by(|a, b| a.identity.partial_cmp(&b.identity).unwrap()).cloned());
            let foreign_best: Vec<Option<Hit>> = hits_by.iter().skip(1 + n_conf).map(|(_, by)| best_in(by)).collect();
            let foreign_id = foreign_best.iter().filter_map(|h| h.as_ref().map(|h| h.identity)).fold(None, |m: Option<f64>, x| Some(m.map_or(x, |m| m.max(x))));
            let v = verdict(&VerdictInput { run_p: r.run_p, run_top_frac: r.run_top_frac, is_ig_tr: r.is_ig_tr, n_psv: r.n_psv, shared_frac: r.shared_frac, editing_frac: r.editing_frac, host_identity: host_id, other_identity: other_id, foreign_identity: foreign_id, delta: r.delta });
            *counts.entry(v.as_str()).or_insert(0) += 1;
            let consistent = if r.n_psv >= 3 && r.shared_frac >= 0.5 { "copy_consistent" } else { "scattered" };
            f.extend([
                format!("{:.1}", r.delta / args.pi),
                consistent.to_string(),
                host_id.map(|x| format!("{x:.4}")).unwrap_or_else(|| "NA".into()),
                other_id.map(|x| format!("{x:.4}")).unwrap_or_else(|| "NA".into()),
                other.as_ref().map(|h| format!("{}:{}-{}", h.chrom, h.start, h.end)).unwrap_or_else(|| "NA".into()),
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
                f.push(best.as_ref().map(|h| format!("{:.4}", h.identity)).unwrap_or_else(|| "NA".into()));
                f.push(best.as_ref().map(|h| format!("{}:{}-{}", h.chrom, h.start, h.end)).unwrap_or_else(|| "NA".into()));
                f.push((ok as u8).to_string());
            }
            for best in &foreign_best {
                f.push(best.as_ref().map(|h| format!("{:.4}", h.identity)).unwrap_or_else(|| "NA".into()));
                f.push(best.as_ref().map(|h| format!("{}:{}-{}", h.chrom, h.start, h.end)).unwrap_or_else(|| "NA".into()));
            }
            n_confirmed += (is_cand && any_conf) as usize;
        } else {
            f.extend(std::iter::repeat("NA".to_string()).take(n_extra));
        }
        writeln!(tsv, "{}", f.join("\t"))?;
    }
    let mut vc: Vec<_> = counts.iter().collect();
    vc.sort();
    eprintln!("[o3_rna_flag] verdicts: {vc:?}");
    if n_conf > 0 {
        eprintln!("[o3_rna_flag] reference_absent_candidate {n_candidates}, confirmed by a confirm genome {n_confirmed}");
    }
    eprintln!("[o3_rna_flag] wrote {tsv_path}");
    Ok(())
}
