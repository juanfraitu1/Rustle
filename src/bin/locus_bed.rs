//! Emit predicted loci as BED, and match them one-to-one against an annotation to show size agreement.
//!
//! Native Rust replacement for `bench/locus_bed.py` (§6r9), reproducing it byte for byte.
//!
//! Collapses a GTF to LOCI (one record per `gene_id`, spanning its transcripts), writes them as BED, and
//! greedily matches predicted against annotated loci on reciprocal overlap, reporting
//! `size_ratio = pred_span / ref_span`.
//!
//! ⚠ The matching here is EVALUATION ONLY. Loci are never built with bipartite matching — that is a
//! standing project constraint; this only scores loci that were already built without it. The greedy pass
//! is a LOWER BOUND on the optimal assignment, not the optimum; ties break on locus id so it is
//! deterministic.
//!
//! usage: locus_bed PRED.gtf --out PREFIX [--ref REF.gtf] [--min-overlap 0.10]

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};

/// chrom, start (0-based), end, strand, summed reads, transcript count.
type Locus = (String, i64, i64, String, u64, u32);

fn attr<'a>(s: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("{key} \"");
    let i = s.find(&pat)? + pat.len();
    let j = s[i..].find('"')? + i;
    Some(&s[i..j])
}

/// Collapse a GTF into loci keyed by `gene_id`, preserving first-seen order.
fn loci(path: &str) -> Result<(Vec<String>, HashMap<String, Locus>)> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
    // per transcript: exons, plus its gene and read support
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
        let Some(t) = attr(fs[8], "transcript_id") else { continue };
        if let Some(g) = attr(fs[8], "gene_id") {
            gene_of.entry(t.to_string()).or_insert_with(|| g.to_string());
        }
        if let Some(r) = attr(fs[8], "reads").and_then(|v| v.parse::<u64>().ok()) {
            let e = reads_of.entry(t.to_string()).or_insert(0);
            *e = (*e).max(r);
        }
        if fs[2] == "exon" {
            let (Ok(a), Ok(b)) = (fs[3].parse::<i64>(), fs[4].parse::<i64>()) else { continue };
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
        let (c, s, e, st) = (ex[0].0.clone(), ex[0].1, ex[ex.len() - 1].2, ex[0].3.clone());
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
        let strand = if st == "+" || st == "-" { st.as_str() } else { "." };
        writeln!(fo, "{c}\t{s}\t{e}\t{g}\t{}\t{strand}", (*r).min(1000))?;
    }
    Ok(())
}

fn main() -> Result<()> {
    let argv: Vec<String> = std::env::args().collect();
    let mut pred = String::new();
    let (mut out, mut refgtf, mut min_ov) = (String::new(), String::new(), 0.10f64);
    let mut i = 1;
    while i < argv.len() {
        match argv[i].as_str() {
            "--out" => { out = argv[i + 1].clone(); i += 2 }
            "--ref" => { refgtf = argv[i + 1].clone(); i += 2 }
            "--min-overlap" => { min_ov = argv[i + 1].parse().unwrap_or(0.10); i += 2 }
            a => { pred = a.to_string(); i += 1 }
        }
    }
    if pred.is_empty() || out.is_empty() {
        eprintln!("usage: locus_bed PRED.gtf --out PREFIX [--ref REF.gtf] [--min-overlap 0.10]");
        std::process::exit(2);
    }

    let (porder, p) = loci(&pred)?;
    write_bed(&porder, &p, &format!("{out}.loci.bed"))?;
    eprintln!("  {} predicted loci -> {out}.loci.bed", p.len());
    if refgtf.is_empty() {
        return Ok(());
    }
    let (rorder, r) = loci(&refgtf)?;
    write_bed(&rorder, &r, &format!("{out}.ref_loci.bed"))?;
    eprintln!("  {} annotated loci -> {out}.ref_loci.bed", r.len());

    // candidate pairs: same contig, reciprocal overlap >= min_ov
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
            let rec = (ov as f64 / (pv.2 - pv.1).max(1) as f64).min(ov as f64 / (rv.2 - rv.1).max(1) as f64);
            if rec >= min_ov {
                cands.push((rec, pg, rg));
            }
        }
    }
    // deterministic ties: best overlap first, then locus id
    cands.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap().then(a.1.cmp(b.1)).then(a.2.cmp(b.2)));
    let (mut usedp, mut usedr) = (std::collections::HashSet::new(), std::collections::HashSet::new());
    let mut pairs: Vec<(&String, &String, f64)> = Vec::new();
    for (rec, pg, rg) in &cands {
        if usedp.contains(*pg) || usedr.contains(*rg) {
            continue;
        }
        usedp.insert(*pg);
        usedr.insert(*rg);
        pairs.push((pg, rg, *rec));
    }

    let mut fo = std::io::BufWriter::new(std::fs::File::create(format!("{out}.locus_match.tsv"))?);
    writeln!(fo, "pred_locus\tref_locus\tchrom\tpred_start\tpred_end\tref_start\tref_end\tpred_span\tref_span\tsize_ratio\trecip_overlap\tpred_reads\tpred_n_tx")?;
    pairs.sort_by(|a, b| a.0.cmp(b.0));
    let mut ratios: Vec<f64> = Vec::new();
    for (pg, rg, rec) in &pairs {
        let (pc, ps, pe, _, pr, pn) = &p[*pg];
        let (_, rs, re_, _, _, _) = &r[*rg];
        let (psp, rsp) = (pe - ps, re_ - rs);
        let ratio = if rsp != 0 { psp as f64 / rsp as f64 } else { f64::NAN };
        ratios.push(ratio);
        writeln!(fo, "{pg}\t{rg}\t{pc}\t{ps}\t{pe}\t{rs}\t{re_}\t{psp}\t{rsp}\t{ratio:.4}\t{rec:.4}\t{pr}\t{pn}")?;
    }
    fo.flush()?;
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = ratios.len();
    let within = |lo: f64, hi: f64| ratios.iter().filter(|&&x| x >= lo && x <= hi).count();
    println!("\n  GREEDY one-to-one matching (evaluation only), min reciprocal overlap {min_ov}");
    println!("    matched pairs            : {n}  of {} predicted / {} annotated", p.len(), r.len());
    println!("    predicted loci unmatched : {}    annotated loci unmatched: {}", p.len() - n, r.len() - n);
    if n > 0 {
        println!("    size ratio pred/ref      : median {:.3}  q25 {:.3}  q75 {:.3}", ratios[n / 2], ratios[n / 4], ratios[3 * n / 4]);
        for (lo, hi) in [(0.9, 1.1), (0.8, 1.25), (0.5, 2.0)] {
            let w = within(lo, hi);
            println!("      within [{lo},{hi}]        : {w} ({:.1}%)", 100.0 * w as f64 / n as f64);
        }
    }
    println!("    -> {out}.locus_match.tsv");
    Ok(())
}
