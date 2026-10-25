//! o3port — the O3 port (`bench/o3_graph/o3port.py run`) in Rust, byte-identical products (plan
//! docs/superpowers/plans/2026-10-24-o3port-rust.md; the algorithm is documented in `rustle::o3port` and in the Python).
//!
//! usage: o3port run --bam B --copies COPIES.tsv --genome G.fa --net net_reads.tsv --out PREFIX [--families F1,F2]
//!        [--slot-order divergence|support] [--isocon PREFIX] [--clusterer cc|lp] [--coverage] [--semi] [--anchor] [--one-site]
//!        [--guard] [--prefer] [--no-reconstruct] [--no-clipping] [--correct] [--no-structure] [--v5-clique] [--absorb]
//!        [--isocon-rule no-pile|via-b] [--dump PATH] [--threads N]
//! `--threads` (default 4) never changes a product. Exit status 2: a usage error.

use anyhow::{bail, Context, Result};
use rustle::o3port::pile::{pile_records, Bam};
use rustle::o3port::{Clusterer, Copy, Params, ReadId, SlotOrder};
use std::collections::{BTreeMap, BTreeSet, HashMap};

const USAGE: &str = "usage: o3port run --bam B --copies COPIES.tsv --genome G.fa --net net_reads.tsv --out PREFIX [--families F1,F2] \
[--slot-order divergence|support] [--isocon PREFIX] [--clusterer cc|lp] [--coverage] [--semi] [--anchor] [--one-site] [--guard] \
[--prefer] [--no-reconstruct] [--no-clipping] [--correct] [--no-structure] [--v5-clique] [--absorb] [--isocon-rule no-pile|via-b] \
[--dump PATH] [--threads N]";

#[derive(Debug)]
struct Args {
    bam: String,
    copies: String,
    genome: String,
    net: String,
    out: String,
    families: Option<BTreeSet<String>>,
    isocon: Option<String>,
    isocon_via_b: bool,
    dump: Option<String>,
    threads: usize,
    params: Params,
}

fn main() {
    let raw: Vec<String> = std::env::args().skip(1).collect();
    match parse_args(&raw).and_then(|a| run(&a)) {
        Ok(()) => {}
        Err(e) if e.to_string().starts_with("usage") || e.to_string().starts_with("o3port: ") => {
            eprintln!("{e}");
            std::process::exit(2);
        }
        Err(e) => {
            eprintln!("o3port: {e:#}");
            std::process::exit(1);
        }
    }
}

fn parse_args(raw: &[String]) -> Result<Args> {
    if raw.first().map(String::as_str) != Some("run") {
        bail!("{USAGE}");
    }
    let mut kv: BTreeMap<&str, String> = BTreeMap::new();
    let mut p = Params::default();
    let mut it = raw[1..].iter();
    while let Some(a) = it.next() {
        let (flag, inline) = match a.split_once('=') {
            Some((f, v)) if f.starts_with("--") => (f, Some(v.to_string())),
            _ => (a.as_str(), None),
        };
        match flag {
            "--coverage" => p.coverage = true,
            "--semi" => p.semi = true,
            "--anchor" => p.anchor = true,
            "--one-site" => p.one_site = true,
            "--guard" => p.guard = true,
            "--prefer" => p.prefer = true,
            "--no-reconstruct" => p.reconstruct = false,
            "--no-clipping" => p.clipping = false,
            "--correct" => p.correct = true,
            "--no-structure" => p.structure = false,
            "--v5-clique" => p.floor_first = false,
            "--absorb" => p.absorb = true,
            "--bam" | "--copies" | "--genome" | "--net" | "--out" | "--families" | "--slot-order" | "--isocon" | "--clusterer"
            | "--isocon-rule" | "--dump" | "--threads" => {
                let v = match inline {
                    Some(v) => v,
                    None => it.next().cloned().with_context(|| format!("o3port: {flag} needs a value\n{USAGE}"))?,
                };
                kv.insert(flag, v);
            }
            "-h" | "--help" => bail!("{USAGE}"),
            _ => bail!("o3port: unknown argument {a:?}\n{USAGE}"),
        }
    }
    let need = |k: &str| kv.get(k).cloned().with_context(|| format!("o3port: {k} is required\n{USAGE}"));
    p.clusterer = match kv.get("--clusterer").map(String::as_str) {
        None | Some("cc") => Clusterer::Cc,
        Some("lp") => Clusterer::Lp,
        Some(v) => bail!("o3port: --clusterer {v:?} (cc or lp)\n{USAGE}"),
    };
    p.slot_order = match kv.get("--slot-order").map(String::as_str) {
        None | Some("divergence") => SlotOrder::Divergence,
        Some("support") => SlotOrder::Support,
        Some(v) => bail!("o3port: --slot-order {v:?} (divergence or support)\n{USAGE}"),
    };
    let isocon_via_b = match kv.get("--isocon-rule").map(String::as_str) {
        None | Some("no-pile") => false,
        Some("via-b") => true,
        Some(v) => bail!("o3port: --isocon-rule {v:?} (no-pile or via-b)\n{USAGE}"),
    };
    let threads = match kv.get("--threads") {
        None => 4,
        Some(v) => v.parse::<usize>().ok().filter(|&n| n > 0).with_context(|| format!("o3port: --threads {v:?}\n{USAGE}"))?,
    };
    Ok(Args {
        bam: need("--bam")?,
        copies: need("--copies")?,
        genome: need("--genome")?,
        net: need("--net")?,
        out: need("--out")?,
        families: kv.get("--families").filter(|v| !v.is_empty()).map(|v| v.split(',').map(str::to_string).collect()),
        isocon: kv.get("--isocon").filter(|v| !v.is_empty()).cloned(),
        isocon_via_b,
        dump: kv.get("--dump").filter(|v| !v.is_empty()).cloned(),
        threads,
        params: p,
    })
}

/// Rows of a TSV with a header, as {column: value}.
fn tsv(path: &str) -> Result<Vec<BTreeMap<String, String>>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("reading {path}"))?;
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().unwrap_or("").split('\t').collect();
    Ok(lines
        .filter(|l| !l.is_empty())
        .map(|l| head.iter().zip(l.split('\t')).map(|(h, v)| (h.to_string(), v.to_string())).collect())
        .collect())
}

/// `load_copies`: {family: [copy]} in file order; a copy spans its first sorted exon's start to its last sorted exon's end.
fn load_copies(path: &str) -> Result<BTreeMap<String, Vec<Copy>>> {
    let mut by_fam: BTreeMap<String, Vec<Copy>> = BTreeMap::new();
    for r in tsv(path)? {
        let mut blocks: Vec<(i64, i64)> = Vec::new();
        for b in r["exons"].split(',').filter(|b| !b.is_empty()) {
            let (s, e) = b.split_once('-').with_context(|| format!("{path}: exon {b:?}"))?;
            blocks.push((s.parse()?, e.parse()?));
        }
        blocks.sort();
        let (start, end) = (blocks[0].0, blocks[blocks.len() - 1].1);
        by_fam.entry(r["family_id"].clone()).or_default().push(Copy {
            tid: r["tid"].clone(),
            chrom: r["chrom"].clone(),
            start,
            end,
        });
    }
    Ok(by_fam)
}

/// `load_net`: ({family: used reads}, {family: rows}) of the families asked for.
#[allow(clippy::type_complexity)]
fn load_net(path: &str, families: &Option<BTreeSet<String>>) -> Result<(BTreeMap<String, BTreeSet<String>>, BTreeMap<String, usize>)> {
    let (mut net, mut n_net): (BTreeMap<String, BTreeSet<String>>, BTreeMap<String, usize>) = Default::default();
    for r in tsv(path)? {
        let f = &r["family"];
        if families.as_ref().is_some_and(|fs| !fs.contains(f)) {
            continue;
        }
        *n_net.entry(f.clone()).or_default() += 1;
        if r["used"] == "1" {
            net.entry(f.clone()).or_default().insert(r["read"].clone());
        }
    }
    Ok((net, n_net))
}

fn run(a: &Args) -> Result<()> {
    let by_fam = load_copies(&a.copies)?;
    let (net, _n_net) = load_net(&a.net, &a.families)?;
    let _ = (&a.genome, &a.out, &a.isocon, a.isocon_via_b, &a.dump, a.threads, &a.params);
    // the run's read table: every used net read, sorted, so ReadId order = name order
    let table: Vec<String> = net.values().flatten().cloned().collect::<BTreeSet<_>>().into_iter().collect();
    let names: HashMap<Vec<u8>, ReadId> = table.iter().enumerate().map(|(i, n)| (n.as_bytes().to_vec(), i as ReadId)).collect();
    let bam = Bam::open(&a.bam)?;
    let (prim, _records) = bam.read_primaries(&names)?;
    if std::env::var_os("O3PORT_DEBUG_PILES").is_some() {
        for (fam, reads) in &net {
            let ids: BTreeSet<ReadId> = reads.iter().map(|r| names[r.as_bytes()]).collect();
            let piles = pile_records(&bam, &by_fam[fam], &ids, &prim, &names)?;
            for (c, p) in &piles {
                let mut v: Vec<&str> = p.keys().map(|&i| table[i as usize].as_str()).collect();
                v.sort();
                println!("{fam}\t{c}\t{}", v.join(","));
            }
        }
    }
    Ok(())
}
