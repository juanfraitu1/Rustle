//! `semifam` — SEMI-GUIDED FAMILIES (docs/LEDGER_O1.md PREREG SEMIGUIDED, register 1472): the Rust port of
//! `bench/denovo_eval/semifam.py`, byte-identical output. De novo families (`mcl_families` copies.tsv + loci.gff3) plus the
//! copies that the reads' secondary alignments reveal (`semisec` cands.tsv): every candidate that overlaps no copy of the run is
//! kept when the SEED LOCUS's own spliced sequence (the de novo locus holding the best read's seed pile) aligns onto the candidate
//! span (minimap2 -x splice, template coverage >= --cov, identity >= --min-id); the copy's exons are the aligned blocks and it
//! joins the seed locus's family (source `semisec`, member_status `silent`); a seed locus outside every family forms a new
//! family `SEMI<n>` with its confirmed copies. Families are not re-clustered.
//!   semifam --copies BASE.fam.copies.tsv --loci BASE.fam.loci.gff3 --semisec SEMISEC_DIR --genome FA --out OUT
//!           [--cov 0.5] [--min-id 0.7] [--pad 2000] [--threads 8] [--no-new-families] [--silent-only]
//! Writes OUT.fam.copies.tsv and OUT.semi.tsv (one row per candidate with its seed locus, family, gate numbers and verdict).
use anyhow::{Context, Result};
use clap::Parser;
use std::collections::HashMap;
use std::io::Write;

#[derive(Parser, Debug)]
#[command(about = "Semi-guided families: de novo families + copies found through secondary alignments, confirmed by the seed locus.")]
struct Args {
    /// `mcl_families` copies table of the base run.
    #[arg(long)]
    copies: String,
    /// `mcl_families` loci GFF3 of the base run (`<out>.fam.loci.gff3`).
    #[arg(long)]
    loci: String,
    /// `semisec` output directory (cands.tsv).
    #[arg(long)]
    semisec: String,
    /// Reference FASTA with .fai.
    #[arg(long)]
    genome: String,
    /// Output prefix.
    #[arg(long)]
    out: String,
    /// Template coverage gate.
    #[arg(long, default_value_t = 0.5)]
    cov: f64,
    /// Identity gate (matches / alignment block length).
    #[arg(long, default_value_t = 0.7)]
    min_id: f64,
    /// Padding around the candidate span (bp).
    #[arg(long, default_value_t = 2000)]
    pad: i64,
    /// minimap2 threads.
    #[arg(long, default_value_t = 8)]
    threads: usize,
    /// Keep base families only: candidates seeded by a locus outside every family are reported, not attached.
    #[arg(long, default_value_t = false)]
    no_new_families: bool,
    /// Only `semisec` SILENT candidates (< 2 primaries); default = every candidate (SEMIGUIDED amendment 1).
    #[arg(long, default_value_t = false)]
    silent_only: bool,
}

struct Locus {
    chrom: String,
    strand: char,
    exons: Vec<(i64, i64)>,
}

struct Cand {
    fields: HashMap<String, String>,
}

/// Interval index: per key, intervals sorted by start with the longest length, as the Python `loci_over` / `overlaps_copy`.
struct Index {
    v: Vec<(i64, i64, String)>,
    maxlen: i64,
}

impl Index {
    fn new(mut v: Vec<(i64, i64, String)>) -> Self {
        v.sort();
        let maxlen = v.iter().map(|x| x.1 - x.0).max().unwrap_or(0);
        Index { v, maxlen }
    }
    fn over(&self, s: i64, e: i64) -> Vec<(i64, String)> {
        let mut i = self.v.partition_point(|x| x.0 < s - self.maxlen);
        let mut out = Vec::new();
        while i < self.v.len() && self.v[i].0 < e {
            if self.v[i].1 > s {
                out.push((e.min(self.v[i].1) - s.max(self.v[i].0), self.v[i].2.clone()));
            }
            i += 1;
        }
        out
    }
    fn any(&self, s: i64, e: i64) -> bool {
        let mut i = self.v.partition_point(|x| x.0 < s - 2_000_000);
        while i < self.v.len() && self.v[i].0 < e {
            if self.v[i].1 > s {
                return true;
            }
            i += 1;
        }
        false
    }
}

fn revcomp_upper(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            x => x,
        })
        .collect()
}

fn main() -> Result<()> {
    let args = Args::parse();
    // base copies: header + raw rows
    let text = std::fs::read_to_string(&args.copies).with_context(|| format!("reading {}", args.copies))?;
    let mut lines = text.lines();
    let hdr: Vec<String> = lines.next().context("empty copies table")?.split('\t').map(|s| s.to_string()).collect();
    let col = |name: &str| hdr.iter().position(|h| h == name).with_context(|| format!("copies table lacks column {name}"));
    let (c_fam, c_tid, c_chrom, c_start, c_end, c_strand, c_exons) =
        (col("family_id")?, col("tid")?, col("chrom")?, col("start")?, col("end")?, col("strand")?, col("exons")?);
    let rows: Vec<Vec<String>> = lines.filter(|l| !l.is_empty()).map(|l| l.split('\t').map(|s| s.to_string()).collect()).collect();
    let mut fam_of: HashMap<String, String> = HashMap::new();
    let mut row_of: HashMap<String, usize> = HashMap::new();
    let mut copy_iv: HashMap<(String, char), Vec<(i64, i64, String)>> = HashMap::new();
    for (i, r) in rows.iter().enumerate() {
        fam_of.insert(r[c_tid].clone(), r[c_fam].clone());
        row_of.insert(r[c_tid].clone(), i);
        let st = r[c_strand].chars().next().unwrap_or('+');
        copy_iv.entry((r[c_chrom].clone(), st)).or_default().push((r[c_start].parse()?, r[c_end].parse()?, r[c_tid].clone()));
    }
    let copy_idx: HashMap<(String, char), Index> = copy_iv.into_iter().map(|(k, v)| (k, Index::new(v))).collect();
    // de novo loci (all records, in families or not)
    let mut loci: HashMap<String, Locus> = HashMap::new();
    let mut order: Vec<String> = Vec::new();
    for line in std::fs::read_to_string(&args.loci)?.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 || f[2] != "exon" {
            continue;
        }
        let g = f[8].split("gene=").nth(1).context("exon without gene=")?.split(';').next().unwrap().to_string();
        let strand = if f[6] == "+" || f[6] == "-" { f[6].chars().next().unwrap() } else { '+' };
        let e = loci.entry(g.clone()).or_insert_with(|| {
            order.push(g.clone());
            Locus { chrom: f[0].to_string(), strand, exons: Vec::new() }
        });
        e.exons.push((f[3].parse::<i64>()? - 1, f[4].parse()?));
    }
    let mut by_chrom: HashMap<String, Vec<(i64, i64, String)>> = HashMap::new();
    for g in &order {
        let l = loci.get_mut(g).unwrap();
        l.exons.sort();
        by_chrom.entry(l.chrom.clone()).or_default().push((l.exons[0].0, l.exons.last().unwrap().1, g.clone()));
    }
    let loci_idx: HashMap<String, Index> = by_chrom.into_iter().map(|(k, v)| (k, Index::new(v))).collect();
    // candidates
    let ctext = std::fs::read_to_string(format!("{}/cands.tsv", args.semisec))?;
    let mut cl = ctext.lines();
    let chdr: Vec<&str> = cl.next().context("empty cands.tsv")?.split('\t').collect();
    let mut cands: Vec<Cand> = Vec::new();
    for line in cl.filter(|l| !l.is_empty()) {
        let fields: HashMap<String, String> = chdr.iter().zip(line.split('\t')).map(|(k, v)| (k.to_string(), v.to_string())).collect();
        if args.silent_only && fields["silent"] != "true" {
            continue;
        }
        cands.push(Cand { fields });
    }
    let n_fams = fam_of.values().collect::<std::collections::HashSet<_>>().len();
    eprintln!("{} base copies in {} families; {} de novo loci; {} semisec candidates", rows.len(), n_fams, loci.len(), cands.len());
    // seed locus per candidate
    enum Plan {
        Skip(&'static str),
        Align { seed: String, chrom: String, strand: char, exons: Vec<(i64, i64)> },
    }
    let mut plan: Vec<Plan> = Vec::new();
    for c in &cands {
        let (sc, iv) = c.fields["best_seed"].rsplit_once(':').context("best_seed")?;
        let (ss, se) = iv.split_once('-').context("best_seed")?;
        let (ss, se): (i64, i64) = (ss.parse()?, se.parse()?);
        let hits = loci_idx.get(sc).map(|ix| ix.over(ss, se)).unwrap_or_default();
        let (s, e): (i64, i64) = (c.fields["start"].parse()?, c.fields["end"].parse()?);
        let cst = c.fields["strand"].chars().next().unwrap_or('+');
        if copy_idx.get(&(c.fields["chrom"].clone(), cst)).map(|ix| ix.any(s, e)).unwrap_or(false) {
            plan.push(Plan::Skip("overlaps_copy"));
            continue;
        }
        if hits.is_empty() {
            plan.push(Plan::Skip("seed_no_locus"));
            continue;
        }
        // a copy already in a family first (overlapping loci records are one locus), else the loci record with most overlap
        let mut inc: Vec<(i64, String)> = hits.iter().filter(|(_, g)| row_of.contains_key(g)).cloned().collect();
        if inc.is_empty() {
            let best_l = hits.iter().max().unwrap().1.clone();
            let l = &loci[&best_l];
            let (ls, le) = (l.exons[0].0, l.exons.last().unwrap().1);
            if let Some(ix) = copy_idx.get(&(l.chrom.clone(), l.strand)) {
                for (cs, ce, t) in &ix.v {
                    if *ce > ls && *cs < le {
                        inc.push((le.min(*ce) - ls.max(*cs), t.clone()));
                    }
                }
            }
        }
        let seed = if inc.is_empty() { hits.iter().max().unwrap().1.clone() } else { inc.iter().max().unwrap().1.clone() };
        let (chrom, strand, exons) = if let Some(&i) = row_of.get(&seed) {
            let r = &rows[i];
            let ex: Vec<(i64, i64)> = r[c_exons]
                .split(',')
                .filter(|b| !b.is_empty())
                .map(|b| {
                    let (a, bb) = b.split_once('-').unwrap();
                    (a.parse().unwrap(), bb.parse().unwrap())
                })
                .collect();
            (r[c_chrom].clone(), r[c_strand].chars().next().unwrap_or('+'), ex)
        } else {
            let l = &loci[&seed];
            (l.chrom.clone(), l.strand, l.exons.clone())
        };
        plan.push(Plan::Align { seed, chrom, strand, exons });
    }
    // sequences
    let genome = rustle::genome::IndexedFasta::open(&args.genome)?;
    let tmp = std::env::temp_dir().join(format!("semifam_{}", std::process::id()));
    std::fs::create_dir_all(&tmp)?;
    let mut region: HashMap<usize, (i64, i64)> = HashMap::new();
    {
        let mut wt = std::io::BufWriter::new(std::fs::File::create(tmp.join("tpl.fa"))?);
        let mut wr = std::io::BufWriter::new(std::fs::File::create(tmp.join("reg.fa"))?);
        let mut chrom_len: HashMap<String, i64> = HashMap::new();
        let mut len_of = |chrom: &str| -> i64 {
            *chrom_len.entry(chrom.to_string()).or_insert_with(|| genome.contig_length(chrom).unwrap_or(0) as i64)
        };
        for (i, p) in plan.iter().enumerate() {
            if let Plan::Align { chrom, strand, exons, .. } = p {
                let mut t: Vec<u8> = Vec::new();
                for (s0, e0) in exons {
                    t.extend(genome.fetch(chrom, *s0, *e0).unwrap_or_default());
                }
                t.make_ascii_uppercase();
                if *strand == '-' {
                    t = revcomp_upper(&t);
                }
                let c = &cands[i];
                let cc = &c.fields["chrom"];
                let rs = (c.fields["start"].parse::<i64>()? - args.pad).max(0);
                let re_ = (c.fields["end"].parse::<i64>()? + args.pad).min(len_of(cc));
                region.insert(i, (rs, re_));
                let mut r = genome.fetch(cc, rs, re_).unwrap_or_default();
                r.make_ascii_uppercase();
                writeln!(wt, ">c{i}")?;
                wt.write_all(&t)?;
                writeln!(wt)?;
                writeln!(wr, ">c{i}")?;
                wr.write_all(&r)?;
                writeln!(wr)?;
            }
        }
    }
    let out = std::process::Command::new("minimap2")
        .args(["-c", "-x", "splice", "--secondary=yes", "-N", "200", "-p", "0.05", "-t", &args.threads.to_string()])
        .arg(tmp.join("reg.fa"))
        .arg(tmp.join("tpl.fa"))
        .stderr(std::process::Stdio::null())
        .output()
        .context("running minimap2")?;
    anyhow::ensure!(out.status.success(), "minimap2 failed");
    // best alignment per candidate: max cov x identity, first wins on ties
    let mut best: HashMap<usize, (f64, f64, char, Vec<(i64, i64)>)> = HashMap::new();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 || f[0] != f[5] {
            continue;
        }
        let i: usize = f[0][1..].parse()?;
        let ql: f64 = f[1].parse()?;
        let cov = (f[3].parse::<f64>()? - f[2].parse::<f64>()?) / ql;
        let bl: f64 = f[10].parse()?;
        let idn = if bl > 0.0 { f[9].parse::<f64>()? / bl } else { 0.0 };
        let cg = f[12..].iter().find(|x| x.starts_with("cg:Z:")).map(|x| &x[5..]).unwrap_or("");
        let mut blocks = Vec::new();
        let mut pos: i64 = f[7].parse()?;
        let mut bs = pos;
        let mut n = String::new();
        for ch in cg.chars() {
            if ch.is_ascii_digit() {
                n.push(ch);
                continue;
            }
            let k: i64 = n.parse().unwrap_or(0);
            n.clear();
            match ch {
                'M' | 'D' | 'X' | '=' => pos += k,
                'N' => {
                    blocks.push((bs, pos));
                    pos += k;
                    bs = pos;
                }
                _ => {}
            }
        }
        blocks.push((bs, pos));
        let better = match best.get(&i) {
            Some(b) => cov * idn > b.0 * b.1,
            None => true,
        };
        if better {
            best.insert(i, (cov, idn, f[4].chars().next().unwrap(), blocks));
        }
    }
    let _ = std::fs::remove_dir_all(&tmp);
    // attach
    let mut new_rows: Vec<Vec<String>> = Vec::new();
    let mut newfam: HashMap<String, String> = HashMap::new();
    let mut n_new = 0usize;
    let mut verdicts: HashMap<&str, usize> = HashMap::new();
    let mk_row = |fields: Vec<(&str, String)>| -> Vec<String> {
        let m: HashMap<&str, String> = fields.into_iter().collect();
        hdr.iter().map(|h| m.get(h.as_str()).cloned().unwrap_or_else(|| "NA".to_string())).collect()
    };
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{}.semi.tsv", args.out))?);
    writeln!(w, "chrom\tstrand\tstart\tend\treads\tsilent\tclass\tbest_read\tseed_locus\tfamily\ttpl_cov\tidentity\tn_exon\tverdict")?;
    for (i, p) in plan.iter().enumerate() {
        let c = &cands[i];
        let (mut cov, mut idn, mut fam, mut nex) = (0.0f64, 0.0f64, String::new(), 0usize);
        let mut seed_name = String::new();
        let v: &str = match p {
            Plan::Skip(why) => why,
            Plan::Align { seed, strand: tst, .. } => {
                seed_name = seed.clone();
                let _ = tst;
                match best.get(&i) {
                    None => "no_alignment",
                    Some((bc, bi, pst, bl)) => {
                        cov = *bc;
                        idn = *bi;
                        if cov < args.cov || idn < args.min_id {
                            "below_gate"
                        } else {
                            let rs = region[&i].0;
                            let ex: Vec<(i64, i64)> = bl.iter().map(|(s0, e0)| (rs + s0, rs + e0)).collect();
                            let strand = *pst; // the template is transcript-sense, the region + genomic: the PAF strand IS the copy's strand
                            let mut v = "attached";
                            if copy_idx.get(&(c.fields["chrom"].clone(), strand)).map(|ix| ix.any(ex[0].0, ex.last().unwrap().1)).unwrap_or(false) {
                                v = "overlaps_copy";
                            }
                            let mut f = fam_of.get(seed).cloned();
                            if v != "attached" {
                                fam = f.clone().unwrap_or_default(); // reported as the Python does
                            } else if f.is_none() {
                                if args.no_new_families {
                                    v = "seed_singleton";
                                } else {
                                    if !newfam.contains_key(seed) {
                                        n_new += 1;
                                        newfam.insert(seed.clone(), format!("SEMI{n_new}"));
                                        let l = &loci[seed];
                                        let sex = &l.exons;
                                        new_rows.push(mk_row(vec![
                                            ("family_id", format!("SEMI{n_new}")),
                                            ("copy_idx", "0".into()),
                                            ("tid", seed.clone()),
                                            ("chrom", l.chrom.clone()),
                                            ("start", sex[0].0.to_string()),
                                            ("end", sex.last().unwrap().1.to_string()),
                                            ("n_exon", sex.len().to_string()),
                                            ("strand", l.strand.to_string()),
                                            ("n_reads", "0".into()),
                                            ("exons", sex.iter().map(|(a, b)| format!("{a}-{b}")).collect::<Vec<_>>().join(",")),
                                            ("source", "locus_rep".into()),
                                            ("gene_id", seed.clone()),
                                            ("member_status", "ungated".into()),
                                            ("locus_start", sex[0].0.to_string()),
                                            ("locus_end", sex.last().unwrap().1.to_string()),
                                        ]));
                                    }
                                    f = newfam.get(seed).cloned();
                                }
                            }
                            if v == "attached" {
                                let f = f.unwrap();
                                fam = f.clone();
                                nex = ex.len();
                                new_rows.push(mk_row(vec![
                                    ("family_id", f),
                                    ("copy_idx", "0".into()),
                                    ("tid", format!("SS_{}_{}_{}", c.fields["chrom"], ex[0].0, strand)),
                                    ("chrom", c.fields["chrom"].clone()),
                                    ("start", ex[0].0.to_string()),
                                    ("end", ex.last().unwrap().1.to_string()),
                                    ("n_exon", nex.to_string()),
                                    ("strand", strand.to_string()),
                                    ("n_reads", c.fields["reads"].clone()),
                                    ("exons", ex.iter().map(|(a, b)| format!("{a}-{b}")).collect::<Vec<_>>().join(",")),
                                    ("max_family_identity", format!("{idn:.4}")),
                                    ("source", "semisec".into()),
                                    ("gene_id", seed.clone()),
                                    ("member_status", "silent".into()),
                                    ("locus_start", ex[0].0.to_string()),
                                    ("locus_end", ex.last().unwrap().1.to_string()),
                                ]));
                            }
                            v
                        }
                    }
                }
            }
        };
        *verdicts.entry(v).or_insert(0) += 1;
        let class = c.fields.get("class").cloned().unwrap_or_default();
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{cov:.2}\t{idn:.3}\t{nex}\t{v}",
            c.fields["chrom"], c.fields["strand"], c.fields["start"], c.fields["end"], c.fields["reads"], c.fields["silent"], class, c.fields["best"], seed_name, fam
        )?;
    }
    drop(w);
    let mut cnt: HashMap<String, usize> = HashMap::new();
    for r in &rows {
        *cnt.entry(r[c_fam].clone()).or_insert(0) += 1;
    }
    let c_idx = col("copy_idx")?;
    for r in new_rows.iter_mut() {
        let n = cnt.entry(r[c_fam].clone()).or_insert(0);
        r[c_idx] = n.to_string();
        *n += 1;
    }
    let mut w = std::io::BufWriter::new(std::fs::File::create(format!("{}.fam.copies.tsv", args.out))?);
    writeln!(w, "{}", hdr.join("\t"))?;
    for r in rows.iter().chain(new_rows.iter()) {
        writeln!(w, "{}", r.join("\t"))?;
    }
    let mut vs: Vec<_> = verdicts.iter().collect();
    vs.sort();
    eprintln!(
        "verdicts: {}; new families {}; copies {} -> {}",
        vs.iter().map(|(k, n)| format!("{k} {n}")).collect::<Vec<_>>().join(", "),
        n_new,
        rows.len(),
        rows.len() + new_rows.len()
    );
    Ok(())
}
