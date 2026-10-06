//! Small utility commands merged into one multi-call binary.
//!
//! Subcommands preserve the old standalone CLIs:
//!   utilities gff-to-gtf REFSEQ.gff[.gz] CHROM OUT.gtf
//!   utilities locus-bed PRED.gtf --out PREFIX [--ref REF.gtf] [--min-overlap 0.10]
//!   utilities mcl-port --graph EDGES.tsv [--inflation 2.8] [--prune 1e-9] [--max-iter 100]
//!   utilities parcn --copies-fa FA --mat FA --pat FA --out PREFIX [--minimap2 PATH] [--threads 4]
//!   utilities candidate-augment --fasta G.fa --copies P.fam.copies.tsv --copies-fa P.fam.copies.fa
//!                        --regions P.fam.copies.regions --cand P.cand --out P.aug

use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "utilities", about = "Small utility subcommands")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Extract one chromosome from a RefSeq GFF3 and write a gffread-style GTF.
    GffToGtf {
        /// RefSeq GFF3 (optionally gzipped).
        gff: String,
        /// Chromosome to extract.
        chrom: String,
        /// Output GTF path.
        out: String,
    },
    /// Collapse a GTF to loci, write BED, and optionally match against a reference.
    LocusBed {
        /// Predicted GTF.
        pred: String,
        /// Output prefix.
        #[arg(long)]
        out: String,
        /// Optional reference GTF.
        #[arg(long)]
        r#ref: Option<String>,
        /// Minimum reciprocal overlap.
        #[arg(long, default_value_t = 0.10)]
        min_overlap: f64,
    },
    /// Bit-faithful port of bench/mcl_port.py (the Python MCL comparator).
    MclPort {
        /// Edge file: u<TAB>v<TAB>w per undirected edge.
        #[arg(long)]
        graph: String,
        #[arg(long, default_value_t = 2.8)]
        inflation: f64,
        #[arg(long, default_value_t = 1e-9)]
        prune: f64,
        #[arg(long, default_value_t = 100)]
        max_iter: usize,
    },
    /// Assembly-based paralog-specific copy number (parCN) orchestrator.
    Parcn {
        #[arg(long)]
        copies_fa: String,
        #[arg(long)]
        mat: String,
        #[arg(long)]
        pat: String,
        #[arg(long)]
        out: String,
        #[arg(long, default_value = "minimap2")]
        minimap2: String,
        #[arg(long, default_value_t = 4)]
        threads: usize,
    },
    /// Augment the genome and copies table with flagged candidates.
    CandidateAugment {
        #[arg(long)]
        fasta: String,
        #[arg(long)]
        copies: String,
        #[arg(long)]
        copies_fa: String,
        #[arg(long)]
        regions: String,
        #[arg(long)]
        cand: String,
        #[arg(long)]
        out: String,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Cmd::GffToGtf { gff, chrom, out } => gff_to_gtf::run(&gff, &chrom, &out),
        Cmd::LocusBed {
            pred,
            out,
            r#ref,
            min_overlap,
        } => locus_bed::run(&pred, &out, r#ref.as_deref(), min_overlap),
        Cmd::MclPort {
            graph,
            inflation,
            prune,
            max_iter,
        } => mcl_port::run(&graph, inflation, prune, max_iter),
        Cmd::Parcn {
            copies_fa,
            mat,
            pat,
            out,
            minimap2,
            threads,
        } => parcn::run(&copies_fa, &mat, &pat, &out, &minimap2, threads),
        Cmd::CandidateAugment {
            fasta,
            copies,
            copies_fa,
            regions,
            cand,
            out,
        } => {
            if let Err(e) =
                candidate_augment::run(&fasta, &copies, &copies_fa, &regions, &cand, &out)
            {
                eprintln!("candidate_augment: {e}");
                std::process::exit(2);
            }
            Ok(())
        }
    }
}

mod gff_to_gtf {
    use anyhow::{Context, Result};
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader, Write};

    fn attrs(s: &str) -> HashMap<&str, &str> {
        let mut d = HashMap::new();
        for kv in s.trim_end_matches(';').split(';') {
            if let Some((k, v)) = kv.split_once('=') {
                d.insert(k.trim(), v.trim());
            }
        }
        d
    }

    pub fn run(src: &str, chrom: &str, out: &str) -> Result<()> {
        let file = std::fs::File::open(src).with_context(|| format!("opening {src}"))?;
        let reader: Box<dyn BufRead> = if src.ends_with(".gz") {
            Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut exons: HashMap<String, Vec<(u64, u64, String, String)>> = HashMap::new();
        let mut order: Vec<String> = Vec::new();
        let mut info: HashMap<String, (String, String)> = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 9 || f[0] != chrom {
                continue;
            }
            let a = attrs(f[8]);
            if f[2] == "exon" {
                let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else {
                    continue;
                };
                for pid in a.get("Parent").copied().unwrap_or("").split(',') {
                    if pid.is_empty() {
                        continue;
                    }
                    let v = exons.entry(pid.to_string()).or_insert_with(|| {
                        order.push(pid.to_string());
                        Vec::new()
                    });
                    v.push((s, e, f[1].to_string(), f[6].to_string()));
                }
            } else if let Some(id) = a.get("ID") {
                if !matches!(f[2], "gene" | "pseudogene" | "CDS" | "region") {
                    let gene = a
                        .get("gene")
                        .or_else(|| a.get("Name"))
                        .copied()
                        .unwrap_or("");
                    info.insert(
                        id.to_string(),
                        (
                            a.get("Parent").copied().unwrap_or("").to_string(),
                            gene.to_string(),
                        ),
                    );
                }
            }
        }

        let mut fo = std::io::BufWriter::new(
            std::fs::File::create(out).with_context(|| format!("creating {out}"))?,
        );
        let mut n = 0usize;
        for tid in &order {
            let ex = exons.get_mut(tid).expect("ordered id must have exons");
            ex.sort();
            let (gid, gname) = info.get(tid).cloned().unwrap_or_default();
            let (src_f, strand) = (ex[0].2.clone(), ex[0].3.clone());
            let at = format!("transcript_id \"{tid}\"; gene_id \"{gid}\"; gene_name \"{gname}\"");
            writeln!(
                fo,
                "{chrom}\t{src_f}\ttranscript\t{}\t{}\t.\t{strand}\t.\t{at}",
                ex[0].0,
                ex[ex.len() - 1].1
            )?;
            for (i, (s, e, _, _)) in ex.iter().enumerate() {
                writeln!(
                    fo,
                    "{chrom}\t{src_f}\texon\t{s}\t{e}\t.\t{strand}\t.\t{at}; exon_number \"{}\";",
                    i + 1
                )?;
            }
            n += 1;
        }
        fo.flush()?;
        eprintln!("{chrom}: {n} transcripts");
        Ok(())
    }
}

mod locus_bed {
    use anyhow::{Context, Result};
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader, Write};

    type Locus = (String, i64, i64, String, u64, u32);

    fn attr<'a>(s: &'a str, key: &str) -> Option<&'a str> {
        let pat = format!("{key} \"");
        let i = s.find(&pat)? + pat.len();
        let j = s[i..].find('"')? + i;
        Some(&s[i..j])
    }

    fn loci(path: &str) -> Result<(Vec<String>, HashMap<String, Locus>)> {
        let f = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
        let mut tx: HashMap<String, Vec<(String, i64, i64, String)>> = HashMap::new();
        let mut gene_of: HashMap<String, String> = HashMap::new();
        let mut reads_of: HashMap<String, u64> = HashMap::new();
        let mut tx_order: Vec<String> = Vec::new();
        for line in BufReader::new(f).lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            let fs: Vec<&str> = line.split('\t').collect();
            if fs.len() < 9 {
                continue;
            }
            let Some(t) = attr(fs[8], "transcript_id") else {
                continue;
            };
            if let Some(g) = attr(fs[8], "gene_id") {
                gene_of
                    .entry(t.to_string())
                    .or_insert_with(|| g.to_string());
            }
            if let Some(r) = attr(fs[8], "reads").and_then(|v| v.parse::<u64>().ok()) {
                let e = reads_of.entry(t.to_string()).or_insert(0);
                *e = (*e).max(r);
            }
            if fs[2] == "exon" {
                let (Ok(a), Ok(b)) = (fs[3].parse::<i64>(), fs[4].parse::<i64>()) else {
                    continue;
                };
                tx.entry(t.to_string())
                    .or_insert_with(|| {
                        tx_order.push(t.to_string());
                        Vec::new()
                    })
                    .push((fs[0].to_string(), a - 1, b, fs[6].to_string()));
            }
        }
        let mut order: Vec<String> = Vec::new();
        let mut out: HashMap<String, Locus> = HashMap::new();
        for t in &tx_order {
            let ex = tx.get_mut(t).expect("ordered tid must have exons");
            ex.sort_by_key(|x| x.1);
            let g = gene_of.get(t).cloned().unwrap_or_else(|| t.clone());
            let (c, s, e, st) = (
                ex[0].0.clone(),
                ex[0].1,
                ex[ex.len() - 1].2,
                ex[0].3.clone(),
            );
            let r = reads_of.get(t).copied().unwrap_or(0);
            match out.get_mut(&g) {
                Some(v) => {
                    v.1 = v.1.min(s);
                    v.2 = v.2.max(e);
                    v.4 += r;
                    v.5 += 1;
                }
                None => {
                    order.push(g.clone());
                    out.insert(g, (c, s, e, st, r, 1));
                }
            }
        }
        Ok((order, out))
    }

    fn write_bed(order: &[String], d: &HashMap<String, Locus>, path: &str) -> Result<()> {
        let mut rows: Vec<(&String, &Locus)> = order.iter().map(|g| (g, &d[g])).collect();
        rows.sort_by(|a, b| a.1 .0.cmp(&b.1 .0).then(a.1 .1.cmp(&b.1 .1)));
        let mut fo = std::io::BufWriter::new(std::fs::File::create(path)?);
        for (g, (c, s, e, st, r, _)) in rows {
            let strand = if st == "+" || st == "-" {
                st.as_str()
            } else {
                "."
            };
            writeln!(fo, "{c}\t{s}\t{e}\t{g}\t{}\t{strand}", (*r).min(1000))?;
        }
        Ok(())
    }

    pub fn run(pred: &str, out: &str, refgtf: Option<&str>, min_ov: f64) -> Result<()> {
        let (porder, p) = loci(pred)?;
        write_bed(&porder, &p, &format!("{out}.loci.bed"))?;
        eprintln!("  {} predicted loci -> {out}.loci.bed", p.len());
        if refgtf.is_none() {
            return Ok(());
        }
        let refgtf = refgtf.unwrap();
        let (rorder, r) = loci(refgtf)?;
        write_bed(&rorder, &r, &format!("{out}.ref_loci.bed"))?;
        eprintln!("  {} annotated loci -> {out}.ref_loci.bed", r.len());

        let mut bych: HashMap<&str, Vec<&String>> = HashMap::new();
        for g in &rorder {
            bych.entry(r[g].0.as_str()).or_default().push(g);
        }
        let mut cands: Vec<(f64, &String, &String)> = Vec::new();
        for pg in &porder {
            let pv = &p[pg];
            for rg in bych.get(pv.0.as_str()).map(|v| v.as_slice()).unwrap_or(&[]) {
                let rv = &r[*rg];
                let ov = pv.2.min(rv.2) - pv.1.max(rv.1);
                if ov <= 0 {
                    continue;
                }
                let rec = (ov as f64 / (pv.2 - pv.1).max(1) as f64)
                    .min(ov as f64 / (rv.2 - rv.1).max(1) as f64);
                if rec >= min_ov {
                    cands.push((rec, pg, rg));
                }
            }
        }
        cands.sort_by(|a, b| {
            b.0.partial_cmp(&a.0)
                .unwrap()
                .then(a.1.cmp(b.1))
                .then(a.2.cmp(b.2))
        });
        let (mut usedp, mut usedr) = (
            std::collections::HashSet::new(),
            std::collections::HashSet::new(),
        );
        let mut pairs: Vec<(&String, &String, f64)> = Vec::new();
        for (rec, pg, rg) in &cands {
            if usedp.contains(*pg) || usedr.contains(*rg) {
                continue;
            }
            usedp.insert(*pg);
            usedr.insert(*rg);
            pairs.push((pg, rg, *rec));
        }

        let mut fo =
            std::io::BufWriter::new(std::fs::File::create(format!("{out}.locus_match.tsv"))?);
        writeln!(fo, "pred_locus\tref_locus\tchrom\tpred_start\tpred_end\tref_start\tref_end\tpred_span\tref_span\tsize_ratio\trecip_overlap\tpred_reads\tpred_n_tx")?;
        pairs.sort_by(|a, b| a.0.cmp(b.0));
        let mut ratios: Vec<f64> = Vec::new();
        for (pg, rg, rec) in &pairs {
            let (pc, ps, pe, _, pr, pn) = &p[*pg];
            let (_, rs, re_, _, _, _) = &r[*rg];
            let (psp, rsp) = (pe - ps, re_ - rs);
            let ratio = if rsp != 0 {
                psp as f64 / rsp as f64
            } else {
                f64::NAN
            };
            ratios.push(ratio);
            writeln!(fo, "{pg}\t{rg}\t{pc}\t{ps}\t{pe}\t{rs}\t{re_}\t{psp}\t{rsp}\t{ratio:.4}\t{rec:.4}\t{pr}\t{pn}")?;
        }
        fo.flush()?;
        ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = ratios.len();
        let within = |lo: f64, hi: f64| ratios.iter().filter(|&&x| x >= lo && x <= hi).count();
        println!(
            "\n  GREEDY one-to-one matching (evaluation only), min reciprocal overlap {min_ov}"
        );
        println!(
            "    matched pairs            : {n}  of {} predicted / {} annotated",
            p.len(),
            r.len()
        );
        println!(
            "    predicted loci unmatched : {}    annotated loci unmatched: {}",
            p.len() - n,
            r.len() - n
        );
        if n > 0 {
            println!(
                "    size ratio pred/ref      : median {:.3}  q25 {:.3}  q75 {:.3}",
                ratios[n / 2],
                ratios[n / 4],
                ratios[3 * n / 4]
            );
            for (lo, hi) in [(0.9, 1.1), (0.8, 1.25), (0.5, 2.0)] {
                let w = within(lo, hi);
                println!(
                    "      within [{lo},{hi}]        : {w} ({:.1}%)",
                    100.0 * w as f64 / n as f64
                );
            }
        }
        println!("    -> {out}.locus_match.tsv");
        Ok(())
    }
}

mod mcl_port {
    use anyhow::{Context, Result};
    use std::collections::HashMap;
    use std::io::{BufRead, BufReader, Write};

    #[derive(Clone)]
    struct Csc {
        n: usize,
        indptr: Vec<usize>,
        indices: Vec<usize>,
        data: Vec<f64>,
    }

    impl Csc {
        fn from_coo(n: usize, rows: &[usize], cols: &[usize], vals: &[f64]) -> Self {
            let mut per_col: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
            for k in 0..rows.len() {
                per_col[cols[k]].push((rows[k], vals[k]));
            }
            let (mut indptr, mut indices, mut data) = (vec![0usize; n + 1], Vec::new(), Vec::new());
            for j in 0..n {
                let col = &mut per_col[j];
                col.sort_by_key(|e| e.0);
                let mut k = 0;
                while k < col.len() {
                    let (i, mut v) = col[k];
                    let mut m = k + 1;
                    while m < col.len() && col[m].0 == i {
                        v += col[m].1;
                        m += 1;
                    }
                    indices.push(i);
                    data.push(v);
                    k = m;
                }
                indptr[j + 1] = indices.len();
            }
            Csc {
                n,
                indptr,
                indices,
                data,
            }
        }

        fn matmul(&self, b: &Csc) -> Csc {
            let n = self.n;
            let mut sums = vec![0.0f64; n];
            let mut next = vec![usize::MAX; n];
            const HEAD_END: usize = usize::MAX - 1;
            let (mut indptr, mut indices, mut data) = (vec![0usize; n + 1], Vec::new(), Vec::new());
            for j in 0..n {
                let mut head = HEAD_END;
                let mut length = 0usize;
                for kk in b.indptr[j]..b.indptr[j + 1] {
                    let k = b.indices[kk];
                    let v = b.data[kk];
                    for ii in self.indptr[k]..self.indptr[k + 1] {
                        let i = self.indices[ii];
                        sums[i] += v * self.data[ii];
                        if next[i] == usize::MAX {
                            next[i] = head;
                            head = i;
                            length += 1;
                        }
                    }
                }
                for _ in 0..length {
                    if sums[head] != 0.0 {
                        indices.push(head);
                        data.push(sums[head]);
                    }
                    let temp = head;
                    head = next[head];
                    next[temp] = usize::MAX;
                    sums[temp] = 0.0;
                }
                indptr[j + 1] = indices.len();
            }
            Csc {
                n,
                indptr,
                indices,
                data,
            }
        }

        fn norm(&self) -> Csc {
            let n = self.n;
            let mut inv = vec![0.0f64; n];
            for j in 0..n {
                let mut s = 0.0f64;
                for kk in self.indptr[j]..self.indptr[j + 1] {
                    s += self.data[kk];
                }
                if s == 0.0 {
                    s = 1.0;
                }
                inv[j] = 1.0 / s;
            }
            let (mut indptr, mut indices, mut data) = (vec![0usize; n + 1], Vec::new(), Vec::new());
            for j in 0..n {
                let s = self.indptr[j];
                let e = self.indptr[j + 1];
                for kk in (s..e).rev() {
                    let v = inv[j] * self.data[kk];
                    if v != 0.0 {
                        indices.push(self.indices[kk]);
                        data.push(v);
                    }
                }
                indptr[j + 1] = indices.len();
            }
            Csc {
                n,
                indptr,
                indices,
                data,
            }
        }

        fn diff_stats(&self, m: &Csc) -> (bool, f64) {
            let mut any = false;
            let mut mx = 0.0f64;
            let mut row: HashMap<usize, f64> = HashMap::new();
            for j in 0..self.n {
                row.clear();
                for kk in self.indptr[j]..self.indptr[j + 1] {
                    *row.entry(self.indices[kk]).or_insert(0.0) += self.data[kk];
                }
                for kk in m.indptr[j]..m.indptr[j + 1] {
                    *row.entry(m.indices[kk]).or_insert(0.0) -= m.data[kk];
                }
                for &d in row.values() {
                    if d != 0.0 {
                        any = true;
                        let a = d.abs();
                        if a > mx {
                            mx = a;
                        }
                    }
                }
            }
            (!any, mx)
        }
    }

    fn mcl(
        nodes: &[String],
        edges: &[(usize, usize, f64)],
        inflation: f64,
        prune: f64,
        max_iter: usize,
    ) -> Vec<Vec<String>> {
        let n = nodes.len();
        if n == 0 {
            return Vec::new();
        }
        let (mut rows, mut cols, mut vals): (Vec<usize>, Vec<usize>, Vec<f64>) =
            ((0..n).collect(), (0..n).collect(), vec![1.0; n]);
        for &(a, b, w) in edges {
            if a == b {
                continue;
            }
            rows.extend([a, b]);
            cols.extend([b, a]);
            vals.extend([w, w]);
        }
        let mut m = Csc::from_coo(n, &rows, &cols, &vals).norm();
        for _ in 0..max_iter {
            let mut nn = m.matmul(&m);
            for x in nn.data.iter_mut() {
                *x = x.powf(inflation);
            }
            let mut keep_ip = vec![0usize; n + 1];
            let (mut ki, mut kd) = (
                Vec::with_capacity(nn.indices.len()),
                Vec::with_capacity(nn.data.len()),
            );
            for j in 0..n {
                for kk in nn.indptr[j]..nn.indptr[j + 1] {
                    let v = if nn.data[kk] < prune {
                        0.0
                    } else {
                        nn.data[kk]
                    };
                    if v != 0.0 {
                        ki.push(nn.indices[kk]);
                        kd.push(v);
                    }
                }
                keep_ip[j + 1] = ki.len();
            }
            nn = Csc {
                n,
                indptr: keep_ip,
                indices: ki,
                data: kd,
            }
            .norm();
            let (empty, mx) = nn.diff_stats(&m);
            let done = empty || mx < 1e-7;
            m = nn;
            if done {
                break;
            }
        }
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(p: &mut [usize], mut x: usize) -> usize {
            while p[x] != x {
                p[x] = p[p[x]];
                x = p[x];
            }
            x
        }
        for j in 0..n {
            let (s, e) = (m.indptr[j], m.indptr[j + 1]);
            if s == e {
                continue;
            }
            let mut best = s;
            for kk in s + 1..e {
                if m.data[kk] > m.data[best] {
                    best = kk;
                }
            }
            let r = m.indices[best];
            let (ra, rb) = (find(&mut parent, j), find(&mut parent, r));
            if ra != rb {
                parent[ra.max(rb)] = ra.min(rb);
            }
        }
        let mut order: Vec<usize> = Vec::new();
        let mut groups: HashMap<usize, Vec<String>> = HashMap::new();
        for i in 0..n {
            let r = find(&mut parent, i);
            if !groups.contains_key(&r) {
                order.push(r);
            }
            groups.entry(r).or_default().push(nodes[i].clone());
        }
        order
            .into_iter()
            .map(|r| groups.remove(&r).unwrap())
            .collect()
    }

    pub fn run(graph: &str, inflation: f64, prune: f64, max_iter: usize) -> Result<()> {
        let f = std::fs::File::open(graph).with_context(|| format!("opening {}", graph))?;
        let mut raw: Vec<(String, String, f64)> = Vec::new();
        for line in BufReader::new(f).lines() {
            let line = line?;
            let fs: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
            if fs.len() >= 3 {
                raw.push((fs[0].to_string(), fs[1].to_string(), fs[2].parse::<f64>()?));
            }
        }
        let mut nodes: Vec<String> = raw
            .iter()
            .flat_map(|(u, v, _)| [u.clone(), v.clone()])
            .collect();
        nodes.sort();
        nodes.dedup();
        let ix: HashMap<&str, usize> = nodes
            .iter()
            .enumerate()
            .map(|(i, s)| (s.as_str(), i))
            .collect();
        let edges: Vec<(usize, usize, f64)> = raw
            .iter()
            .map(|(u, v, w)| (ix[u.as_str()], ix[v.as_str()], *w))
            .collect();
        let out = std::io::stdout();
        let mut w = std::io::BufWriter::new(out.lock());
        for c in mcl(&nodes, &edges, inflation, prune, max_iter) {
            writeln!(w, "{}", c.join("\t"))?;
        }
        Ok(())
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn two_triangles_joined_by_one_weak_edge_split_at_the_bridge() {
            let nodes: Vec<String> = ["a", "b", "c", "d", "e", "f"]
                .iter()
                .map(|s| s.to_string())
                .collect();
            let edges = vec![
                (0, 1, 1.0),
                (1, 2, 1.0),
                (0, 2, 1.0),
                (3, 4, 1.0),
                (4, 5, 1.0),
                (3, 5, 1.0),
                (2, 3, 0.05),
            ];
            let c = mcl(&nodes, &edges, 2.8, 1e-9, 100);
            assert_eq!(c, vec![vec!["a", "b", "c"], vec!["d", "e", "f"]]);
            let n4: Vec<String> = ["p", "q", "r", "s"].iter().map(|s| s.to_string()).collect();
            let cyc = vec![(0, 1, 1.0), (1, 2, 1.0), (2, 3, 1.0), (0, 3, 1.0)];
            let c4 = mcl(&n4, &cyc, 2.8, 1e-9, 100);
            assert_eq!(c4.iter().map(|c| c.len()).sum::<usize>(), 4);
            assert!(c4.iter().all(|c| !c.is_empty()));
        }
    }
}

mod parcn {
    use anyhow::{Context, Result};
    use std::collections::{BTreeMap, HashMap};
    use std::io::Write;

    use rustle::genome_projection::project_with_cs;
    use rustle::parcn::{
        assign_locus, dedup_loci, format_family_row, format_parcn_row, parse_copies_fa,
        sun_positions, tabulate, Assignment, CopySun, Locus,
    };

    pub fn run(
        copies_fa: &str,
        mat: &str,
        pat: &str,
        out: &str,
        minimap2: &str,
        threads: usize,
    ) -> Result<()> {
        if std::process::Command::new(minimap2)
            .arg("--version")
            .output()
            .is_err()
        {
            anyhow::bail!("minimap2 ('{minimap2}') not found on PATH — parcn requires minimap2");
        }

        let fams = parse_copies_fa(copies_fa).with_context(|| format!("parsing {}", copies_fa))?;

        let mut queries: Vec<(String, Vec<u8>)> = Vec::new();
        for (fam_id, copies) in &fams {
            for c in copies {
                queries.push((format!("{fam_id}|{}", c.copy_id), c.seq.clone()));
            }
        }

        let suns_by_fam: BTreeMap<String, Vec<CopySun>> = fams
            .iter()
            .map(|(f, c)| (f.clone(), sun_positions(c, band_for(c))))
            .collect();

        let mat_by_fam = project_and_assign(&queries, mat, minimap2, threads, &suns_by_fam)
            .with_context(|| format!("projecting onto maternal haplotype {mat}"))?;
        let pat_by_fam = project_and_assign(&queries, pat, minimap2, threads, &suns_by_fam)
            .with_context(|| format!("projecting onto paternal haplotype {pat}"))?;

        let parcn_path = format!("{}.parcn.tsv", out);
        let fam_path = format!("{}.parcn_families.tsv", out);
        let mut pw =
            std::fs::File::create(&parcn_path).with_context(|| format!("creating {parcn_path}"))?;
        let mut fw =
            std::fs::File::create(&fam_path).with_context(|| format!("creating {fam_path}"))?;
        writeln!(
            pw,
            "family_id\tcopy_id\tsun_tier\tloci_mat\tloci_pat\tparCN\tassign_method"
        )?;
        writeln!(fw, "family_id\tn_copies\tfamCN_diploid\tn_unresolved_loci")?;

        let empty: Vec<Assignment> = Vec::new();
        for (fam_id, copies) in &fams {
            let suns = &suns_by_fam[fam_id];
            let mat_assign = mat_by_fam.get(fam_id).unwrap_or(&empty);
            let pat_assign = pat_by_fam.get(fam_id).unwrap_or(&empty);
            let (rows, n_unresolved) = tabulate(fam_id, copies, suns, mat_assign, pat_assign);
            for r in &rows {
                writeln!(pw, "{}", format_parcn_row(r))?;
            }
            writeln!(fw, "{}", format_family_row(fam_id, &rows, n_unresolved))?;
        }

        Ok(())
    }

    fn band_for(copies: &[rustle::parcn::Copy]) -> usize {
        let lens = copies.iter().map(|c| c.seq.len());
        let (lo, hi) = lens.fold((usize::MAX, 0usize), |(lo, hi), l| (lo.min(l), hi.max(l)));
        (if hi < lo { 64 } else { (hi - lo) + 64 }).min(8192)
    }

    fn project_and_assign(
        queries: &[(String, Vec<u8>)],
        target: &str,
        minimap2: &str,
        threads: usize,
        suns_by_fam: &BTreeMap<String, Vec<CopySun>>,
    ) -> Result<HashMap<String, Vec<Assignment>>> {
        let hits = project_with_cs(queries, target, 0.95, 0.90, minimap2, threads)?;
        if hits.is_empty() {
            eprintln!("[parcn] WARNING: no projection hits against haplotype target '{target}' (0 loci this side)");
            return Ok(HashMap::new());
        }
        let mut by_fam: HashMap<String, Vec<Locus>> = HashMap::new();
        for h in &hits {
            let Some((fam, copy)) = h.qname.split_once('|') else {
                continue;
            };
            by_fam.entry(fam.to_string()).or_default().push(Locus {
                chrom: h.chrom.clone(),
                start: h.start,
                end: h.end,
                best_copy: copy.to_string(),
                identity: h.identity,
                runner_up_identity: 0.0,
                cs: h.cs.clone(),
                qs: h.qs,
                qe: h.qe,
                strand: h.strand,
            });
        }
        let mut out: HashMap<String, Vec<Assignment>> = HashMap::new();
        for (fam, loci) in by_fam {
            let Some(suns) = suns_by_fam.get(&fam) else {
                continue;
            };
            let deduped = dedup_loci(loci);
            let mut assignments = Vec::with_capacity(deduped.len());
            for locus in &deduped {
                let Some(sun) = suns.iter().find(|s| s.copy_id == locus.best_copy) else {
                    continue;
                };
                assignments.push(assign_locus(locus, sun));
            }
            out.insert(fam, assignments);
        }
        Ok(out)
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn parcn_end_to_end_heterozygous_mat_pat_split() {
            if std::process::Command::new("minimap2")
                .arg("--version")
                .output()
                .is_err()
            {
                return;
            }
            let dir = std::env::temp_dir();
            let tag = std::process::id();
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
            let base = gen_seq(1, 400);
            let c0 = base.clone();
            let mut c1 = base.clone();
            let mut p = 4usize;
            while p < 400 {
                let cur_idx = bases.iter().position(|&b| b == c0[p]).unwrap();
                c1[p] = bases[(cur_idx + 1) % 4];
                p += 8;
            }
            let copies_fa = dir.join(format!("parcn_e2e_copies_{tag}.fa"));
            std::fs::write(
                &copies_fa,
                format!(
                    ">F|0|c:1-1|+|nexon=1\n{}\n>F|1|c:1-1|+|nexon=1\n{}\n",
                    String::from_utf8_lossy(&c0),
                    String::from_utf8_lossy(&c1)
                ),
            )
            .unwrap();
            let pad = gen_seq(99, 300);
            let mut mat_seq = c0.clone();
            mat_seq.extend(&pad);
            mat_seq.extend(&c1);
            mat_seq.extend(&pad);
            let mut pat_seq = c0.clone();
            pat_seq.extend(&pad);
            let write_hap = |name: &str, seq: &[u8]| {
                let p = dir.join(format!("parcn_e2e_{name}_{tag}.fa"));
                std::fs::write(&p, format!(">h_{name}\n{}\n", String::from_utf8_lossy(seq)))
                    .unwrap();
                p
            };
            let mat = write_hap("mat", &mat_seq);
            let pat = write_hap("pat", &pat_seq);
            let out = dir.join(format!("parcn_e2e_out_{tag}"));
            run(
                copies_fa.to_string_lossy().as_ref(),
                mat.to_string_lossy().as_ref(),
                pat.to_string_lossy().as_ref(),
                out.to_string_lossy().as_ref(),
                "minimap2",
                2,
            )
            .unwrap();
            let parcn =
                std::fs::read_to_string(format!("{}.parcn.tsv", out.to_string_lossy())).unwrap();

            let row = |cp: &str| -> Vec<String> {
                parcn
                    .lines()
                    .find(|l| {
                        let f: Vec<&str> = l.split('\t').collect();
                        f.len() > 1 && f[0] == "F" && f[1] == cp
                    })
                    .unwrap_or_else(|| panic!("no row for copy {cp} in:\n{parcn}"))
                    .split('\t')
                    .map(|s| s.to_string())
                    .collect()
            };
            let r0 = row("0");
            assert_eq!(r0[3], "1", "c0 loci_mat: {r0:?}");
            assert_eq!(r0[4], "1", "c0 loci_pat: {r0:?}");
            assert_eq!(r0[5], "2", "c0 parCN: {r0:?}");
            assert_eq!(r0[6], "SUN", "c0 method: {r0:?}");
            let r1 = row("1");
            assert_eq!(r1[3], "1", "c1 loci_mat: {r1:?}");
            assert_eq!(r1[4], "0", "c1 loci_pat: {r1:?}");
            assert_eq!(r1[5], "1", "c1 parCN: {r1:?}");
            assert_eq!(r1[6], "SUN", "c1 method: {r1:?}");
            assert_ne!(r1[3], r1[4], "c1 must show a real mat/pat split: {r1:?}");

            for p in [copies_fa, mat, pat] {
                std::fs::remove_file(p).ok();
            }
            std::fs::remove_file(format!("{}.parcn.tsv", out.to_string_lossy())).ok();
            std::fs::remove_file(format!("{}.parcn_families.tsv", out.to_string_lossy())).ok();
        }
    }
}

/// Augment the genome, copies table, copies FASTA and regions with flagged o3 candidates
/// (spec `docs/superpowers/specs/2026-10-02-o3-candidates-design.md` §4, §7).
mod candidate_augment {
    use anyhow::{anyhow, bail, Context, Result};
    use std::collections::{HashMap, HashSet};
    use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};

    const COPY_REQUIRED: &[&str] = &[
        "family_id",
        "copy_idx",
        "tid",
        "chrom",
        "start",
        "end",
        "n_exon",
        "strand",
        "n_reads",
        "exons",
    ];
    const CAND_REQUIRED: &[&str] = &["family", "candidate", "flagged", "union_len"];

    pub fn run(
        fasta: &str,
        copies: &str,
        copies_fa: &str,
        regions: &str,
        cand: &str,
        out: &str,
    ) -> Result<()> {
        let cand_path = format!("{cand}.candidates.tsv");
        let contigs_path = format!("{cand}.contigs.fa");

        let cand_table = read_table(&cand_path, CAND_REQUIRED, "candidates table")?;
        let flagged: Vec<usize> = cand_table
            .rows
            .iter()
            .enumerate()
            .filter(|(_, r)| cand_table.get_by_row(r, "flagged") == Some("1"))
            .map(|(i, _)| i)
            .collect();
        if flagged.is_empty() {
            bail!("{cand_path} has no flagged candidate: nothing to augment");
        }
        let contigs = read_fasta(&contigs_path, "contigs")?;
        let (genome_names, names_src) = genome_names(fasta)?;

        let mut seen: HashSet<&str> = HashSet::new();
        for &ri in &flagged {
            let r = &cand_table.rows[ri];
            let cid = cand_table.get_by_row(r, "candidate").unwrap();
            if genome_names.contains(cid) {
                bail!("candidate {cid} already names a sequence of {fasta} ({names_src}); refusing to write {out}.*");
            }
            if !seen.insert(cid) {
                bail!("candidate {cid} is flagged twice in {cand_path}");
            }
            let seq = contigs.get(cid).with_context(|| {
                format!("flagged candidate {cid} has no record in {contigs_path}")
            })?;
            let union_len = cand_table.get_by_row(r, "union_len").unwrap();
            if !union_len.chars().all(|c| c.is_ascii_digit())
                || union_len.parse::<usize>().unwrap_or(0) != seq.len()
                || seq.is_empty()
            {
                bail!(
                    "{cid}: contig of {} bp but union_len {union_len:?} in {cand_path}",
                    seq.len()
                );
            }
        }

        let copy_table = read_table(copies, COPY_REQUIRED, "copies table")?;
        let mut next_idx: HashMap<String, usize> = HashMap::new();
        for r in &copy_table.rows {
            let fid = copy_table.get_by_row(r, "family_id").unwrap().to_string();
            let idx_str = copy_table.get_by_row(r, "copy_idx").unwrap();
            let idx = idx_str
                .parse::<usize>()
                .with_context(|| format!("{copies}: bad copy_idx {idx_str:?} for {fid}"))?;
            next_idx
                .entry(fid)
                .and_modify(|v| *v = (*v).max(idx + 1))
                .or_insert(idx + 1);
        }

        let mut new_rows: Vec<String> = Vec::new();
        let mut new_fa: Vec<String> = Vec::new();
        let mut cand_regions: Vec<String> = Vec::new();
        let mut families: Vec<String> = Vec::new();
        let mut families_set: HashSet<&str> = HashSet::new();

        for &ri in &flagged {
            let r = &cand_table.rows[ri];
            let fid = cand_table.get_by_row(r, "family").unwrap();
            let cid = cand_table.get_by_row(r, "candidate").unwrap();
            let n = contigs[cid].len();
            let idx = *next_idx
                .get(fid)
                .ok_or_else(|| anyhow!("candidate {cid}: family {fid} has no row in {copies}"))?;
            next_idx.insert(fid.to_string(), idx + 1);

            let mut val: HashMap<String, String> = HashMap::new();
            val.insert("family_id".to_string(), fid.to_string());
            val.insert("copy_idx".to_string(), idx.to_string());
            val.insert("tid".to_string(), cid.to_string());
            val.insert("chrom".to_string(), cid.to_string());
            val.insert("start".to_string(), "0".to_string());
            val.insert("end".to_string(), n.to_string());
            val.insert("n_exon".to_string(), "1".to_string());
            val.insert("strand".to_string(), "+".to_string());
            val.insert("n_reads".to_string(), "0".to_string());
            val.insert("exons".to_string(), format!("0-{n}"));
            val.insert("max_family_identity".to_string(), "0".to_string());
            val.insert("source".to_string(), "o3_candidate".to_string());
            val.insert("gene_id".to_string(), ".".to_string());
            val.insert("core_hull".to_string(), "NA".to_string());
            val.insert("sd_depth".to_string(), "0".to_string());
            val.insert("core_bp".to_string(), "0".to_string());
            val.insert("rep_frac".to_string(), "0".to_string());
            val.insert("member_status".to_string(), "candidate".to_string());
            val.insert("locus_start".to_string(), "0".to_string());
            val.insert("locus_end".to_string(), n.to_string());

            new_rows.push(
                copy_table
                    .header
                    .iter()
                    .map(|c| val.get(c).unwrap_or(&"NA".to_string()).clone())
                    .collect::<Vec<_>>()
                    .join("\t"),
            );
            new_fa.push(format!(
                ">{fid}|{idx}|{cid}:0-{n}|+|nexon=1\n{}\n",
                contigs[cid]
            ));
            cand_regions.push(format!("{cid}:0-{n}"));
            if families_set.insert(fid) {
                families.push(fid.to_string());
            }
        }

        let fam_set: HashSet<String> = families.iter().cloned().collect();
        let mut intervals: Vec<(String, i64, i64)> = Vec::new();
        let mut with_region: HashSet<String> = HashSet::new();
        let reg_file = std::fs::File::open(regions)
            .with_context(|| format!("cannot read --regions {regions}"))?;
        for (n, line) in BufReader::new(reg_file).lines().enumerate() {
            let line = line?;
            if line.trim().is_empty() {
                continue;
            }
            let fs: Vec<&str> = line.split('\t').collect();
            if fs.len() < 2 {
                bail!(
                    "{regions} line {}: expected `family<TAB>chrom:start-end`",
                    n + 1
                );
            }
            if fam_set.contains(fs[0]) {
                intervals.push(parse_region(fs[1], &format!("{regions} line {}", n + 1))?);
                with_region.insert(fs[0].to_string());
            }
        }
        let lost: Vec<&str> = families
            .iter()
            .map(|s| s.as_str())
            .filter(|f| !with_region.contains(*f))
            .collect();
        if !lost.is_empty() {
            let suffix = if lost.len() == 1 { "y" } else { "ies" };
            bail!(
                "{regions} has no region for the candidate famil{suffix} {}",
                lost.join(", ")
            );
        }

        // Write every product to `.tmp`, then rename them all. An I/O error removes every `.tmp`
        // file and every product already renamed, then exits 2.
        let mut tmp_files: Vec<(String, String)> = Vec::new();
        let mut done_files: Vec<String> = Vec::new();

        let r: Result<()> = (|| {
            // P.aug.fa
            let fa_tmp = format!("{out}.fa.tmp");
            {
                let mut g = std::fs::File::open(fasta)
                    .with_context(|| format!("cannot read --fasta {fasta}"))?;
                let mut w = std::fs::File::create(&fa_tmp)
                    .with_context(|| format!("cannot create {fa_tmp}"))?;
                let size = std::io::copy(&mut g, &mut w)
                    .with_context(|| format!("cannot copy --fasta {fasta}"))?;
                if size > 0 {
                    g.seek(SeekFrom::End(-1))?;
                    let mut buf = [0u8; 1];
                    g.read_exact(&mut buf)?;
                    if buf[0] != b'\n' {
                        w.write_all(b"\n")?;
                    }
                }
                for &ri in &flagged {
                    let r = &cand_table.rows[ri];
                    let cid = cand_table.get_by_row(r, "candidate").unwrap();
                    let seq = &contigs[cid];
                    write!(w, ">{cid}\n{seq}\n")?;
                }
            }
            tmp_files.push(("fa".to_string(), fa_tmp));

            // P.aug.copies.tsv and P.aug.copies.fa
            for (suffix, src, extra) in [
                ("copies.tsv", copies, &new_rows as &[String]),
                ("copies.fa", copies_fa, &new_fa as &[String]),
            ] {
                let tmp = format!("{out}.{suffix}.tmp");
                let text =
                    std::fs::read_to_string(src).with_context(|| format!("cannot read {src}"))?;
                let mut w =
                    std::fs::File::create(&tmp).with_context(|| format!("cannot create {tmp}"))?;
                w.write_all(text.as_bytes())?;
                if !text.is_empty() && !text.ends_with('\n') {
                    w.write_all(b"\n")?;
                }
                for line in extra {
                    w.write_all(line.as_bytes())?;
                    if !line.ends_with('\n') {
                        w.write_all(b"\n")?;
                    }
                }
                tmp_files.push((suffix.to_string(), tmp));
            }

            // P.aug.regions.txt
            let reg_tmp = format!("{out}.regions.txt.tmp");
            {
                let mut w = std::fs::File::create(&reg_tmp)
                    .with_context(|| format!("cannot create {reg_tmp}"))?;
                for line in merge_regions(&mut intervals) {
                    writeln!(w, "{line}")?;
                }
                for reg in &cand_regions {
                    writeln!(w, "{reg}")?;
                }
            }
            tmp_files.push(("regions.txt".to_string(), reg_tmp));

            // P.aug.families.txt
            let fam_tmp = format!("{out}.families.txt.tmp");
            {
                let mut w = std::fs::File::create(&fam_tmp)
                    .with_context(|| format!("cannot create {fam_tmp}"))?;
                for fid in &families {
                    writeln!(w, "{fid}")?;
                }
            }
            tmp_files.push(("families.txt".to_string(), fam_tmp));

            // Commit: rename every `.tmp` to its final path.
            for (suffix, tmp) in &tmp_files {
                let final_path = format!("{out}.{suffix}");
                std::fs::rename(tmp, &final_path)
                    .with_context(|| format!("cannot rename {tmp} to {final_path}"))?;
                done_files.push(final_path);
            }
            Ok(())
        })();

        if let Err(e) = r {
            for p in &done_files {
                let _ = std::fs::remove_file(p);
            }
            for (_, p) in &tmp_files {
                let _ = std::fs::remove_file(p);
            }
            bail!("cannot make {out}.*: {e} (nothing written)");
        }

        let family_word = if families.len() == 1 {
            "family"
        } else {
            "families"
        };
        eprintln!(
            "candidate_augment: {} candidate contig(s) of {} {family_word} -> {out}.{{fa,copies.tsv,copies.fa,regions.txt,families.txt}}",
            flagged.len(),
            families.len()
        );
        Ok(())
    }

    struct Table {
        header: Vec<String>,
        rows: Vec<Vec<String>>,
    }

    impl Table {
        fn get_by_row<'a>(&'a self, row: &'a [String], col: &str) -> Option<&'a str> {
            self.header
                .iter()
                .position(|h| h == col)
                .map(|i| row[i].as_str())
        }
    }

    fn read_table(path: &str, required: &[&str], what: &str) -> Result<Table> {
        let content =
            std::fs::read_to_string(path).with_context(|| format!("cannot read {what} {path}"))?;
        let mut lines = content.lines();
        let header_line = lines
            .next()
            .ok_or_else(|| anyhow!("{what} {path} is empty (expected a header line)"))?;
        if header_line.is_empty() {
            bail!("{what} {path} is empty (expected a header line)");
        }
        let header: Vec<String> = header_line.split('\t').map(String::from).collect();
        let missing: Vec<_> = required
            .iter()
            .filter(|c| !header.contains(&c.to_string()))
            .copied()
            .collect();
        if !missing.is_empty() {
            bail!(
                "{what} {path} has no column {} (header: {header_line:?})",
                missing.join(", ")
            );
        }
        let mut rows = Vec::new();
        for (n, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let fields: Vec<String> = line.split('\t').map(String::from).collect();
            if fields.len() < header.len() {
                bail!(
                    "{what} {path} line {} has {} fields, the header {}",
                    n + 2,
                    fields.len(),
                    header.len()
                );
            }
            rows.push(fields);
        }
        Ok(Table { header, rows })
    }

    fn read_fasta(path: &str, what: &str) -> Result<HashMap<String, String>> {
        let content =
            std::fs::read_to_string(path).with_context(|| format!("cannot read {what} {path}"))?;
        let mut seqs: HashMap<String, String> = HashMap::new();
        let mut name: Option<String> = None;
        let mut seq = String::new();
        for line in content.lines() {
            let line = line.trim();
            if line.starts_with('>') {
                if let Some(n) = name.take() {
                    seqs.insert(n, std::mem::take(&mut seq));
                }
                let new_name = line[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                if seqs.contains_key(&new_name) {
                    bail!("{what} {path} holds {new_name} twice");
                }
                name = Some(new_name);
            } else if !line.is_empty() {
                if name.is_none() {
                    bail!("{what} {path} has sequence before its first header");
                }
                seq.push_str(line);
            }
        }
        if let Some(n) = name.take() {
            seqs.insert(n, seq);
        }
        Ok(seqs)
    }

    fn genome_names(fasta: &str) -> Result<(HashSet<String>, String)> {
        let mut g =
            std::fs::File::open(fasta).with_context(|| format!("cannot read --fasta {fasta}"))?;
        let mut magic = [0u8; 2];
        let n = g
            .read(&mut magic)
            .with_context(|| format!("cannot read --fasta {fasta}"))?;
        if n == 2 && magic == [0x1f, 0x8b] {
            bail!("{fasta} is compressed: the augmentation appends contigs to a plain FASTA");
        }
        drop(g);
        let fai = format!("{fasta}.fai");
        if std::path::Path::new(&fai).exists() {
            let f = std::fs::File::open(&fai)
                .with_context(|| format!("cannot read the sequence names of --fasta {fasta}"))?;
            let mut names = HashSet::new();
            for line in BufReader::new(f).lines() {
                let line = line?;
                if line.trim().is_empty() {
                    continue;
                }
                let name = line.split('\t').next().unwrap_or("").to_string();
                names.insert(name);
            }
            return Ok((names, fai));
        }
        let f = std::fs::File::open(fasta)
            .with_context(|| format!("cannot read the sequence names of --fasta {fasta}"))?;
        let mut names = HashSet::new();
        for line in BufReader::new(f).lines() {
            let line = line?;
            if line.starts_with('>') {
                let name = line[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string();
                names.insert(name);
            }
        }
        Ok((names, format!("the headers of {fasta}")))
    }

    fn parse_region(text: &str, where_: &str) -> Result<(String, i64, i64)> {
        let (chrom, span) = text
            .rsplit_once(':')
            .ok_or_else(|| anyhow!("bad region {text:?} in {where_} (expected chrom:start-end)"))?;
        let (lo, hi) = span
            .split_once('-')
            .ok_or_else(|| anyhow!("bad region {text:?} in {where_} (expected chrom:start-end)"))?;
        if chrom.is_empty()
            || !lo.chars().all(|c| c.is_ascii_digit())
            || !hi.chars().all(|c| c.is_ascii_digit())
        {
            bail!("bad region {text:?} in {where_} (expected chrom:start-end)");
        }
        Ok((chrom.to_string(), lo.parse()?, hi.parse()?))
    }

    fn merge_regions(intervals: &mut Vec<(String, i64, i64)>) -> Vec<String> {
        intervals.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)).then(a.2.cmp(&b.2)));
        let mut out: Vec<(String, i64, i64)> = Vec::new();
        for (c, lo, hi) in intervals.drain(..) {
            if let Some(last) = out.last_mut() {
                if last.0 == c && lo <= last.2 {
                    last.2 = last.2.max(hi);
                    continue;
                }
            }
            out.push((c, lo, hi));
        }
        out.into_iter()
            .map(|(c, lo, hi)| format!("{c}:{lo}-{hi}"))
            .collect()
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn merge_regions_merges_touching_and_preserves_sort() {
            let mut v = vec![
                ("c1".to_string(), 10, 20),
                ("c1".to_string(), 20, 30),
                ("c2".to_string(), 5, 15),
                ("c1".to_string(), 40, 50),
            ];
            assert_eq!(
                merge_regions(&mut v),
                vec!["c1:10-30", "c1:40-50", "c2:5-15"]
            );
        }

        #[test]
        fn parse_region_splits_on_last_colon() {
            assert_eq!(
                parse_region("chr1:2-3", "test").unwrap(),
                ("chr1".to_string(), 2, 3)
            );
            assert_eq!(
                parse_region("chr1:foo:2-3", "test").unwrap(),
                ("chr1:foo".to_string(), 2, 3)
            );
        }
    }
}
