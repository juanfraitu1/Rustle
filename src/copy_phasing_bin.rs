//! copy_phasing — `bench/o3_graph/o3port.py run` in Rust, byte-identical products (plan
//! docs/superpowers/plans/2026-10-24-o3port-rust.md; the algorithm is documented in `rustle::copy_phasing` and in the Python).
//!
//! usage: copy_phasing run --bam B --copies COPIES.tsv --genome G.fa --net net_reads.tsv --out PREFIX [--families F1,F2]
//!        [--slot-order divergence|support] [--isocon PREFIX] [--clusterer cc|lp] [--coverage] [--semi] [--anchor] [--one-site]
//!        [--guard] [--prefer] [--no-reconstruct] [--no-clipping] [--correct] [--no-structure] [--v5-clique] [--absorb]
//!        [--isocon-rule no-pile|via-b] [--dump PATH] [--threads N]
//! `--threads` (default 4) never changes a product. Exit status 2: a usage error.

use anyhow::{bail, Context, Result};
use rustle::copy_phasing::pile::{pile_records, Bam};
use rustle::genome::IndexedFasta;
use rustle::copy_phasing::{Clusterer, Copy, Params, ReadId, SlotOrder};
use std::collections::{BTreeMap, BTreeSet, HashMap};

const USAGE: &str = "usage: copy_phasing run --bam B --copies COPIES.tsv --genome G.fa --net net_reads.tsv --out PREFIX [--families F1,F2] \
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
        Err(e) if e.to_string().starts_with("usage") || e.to_string().starts_with("copy_phasing: ") => {
            eprintln!("{e}");
            std::process::exit(2);
        }
        Err(e) => {
            eprintln!("copy_phasing: {e:#}");
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
                    None => it.next().cloned().with_context(|| format!("copy_phasing: {flag} needs a value\n{USAGE}"))?,
                };
                kv.insert(flag, v);
            }
            "-h" | "--help" => bail!("{USAGE}"),
            _ => bail!("copy_phasing: unknown argument {a:?}\n{USAGE}"),
        }
    }
    let need = |k: &str| kv.get(k).cloned().with_context(|| format!("copy_phasing: {k} is required\n{USAGE}"));
    p.clusterer = match kv.get("--clusterer").map(String::as_str) {
        None | Some("cc") => Clusterer::Cc,
        Some("lp") => Clusterer::Lp,
        Some(v) => bail!("copy_phasing: --clusterer {v:?} (cc or lp)\n{USAGE}"),
    };
    p.slot_order = match kv.get("--slot-order").map(String::as_str) {
        None | Some("divergence") => SlotOrder::Divergence,
        Some("support") => SlotOrder::Support,
        Some(v) => bail!("copy_phasing: --slot-order {v:?} (divergence or support)\n{USAGE}"),
    };
    let isocon_via_b = match kv.get("--isocon-rule").map(String::as_str) {
        None | Some("no-pile") => false,
        Some("via-b") => true,
        Some(v) => bail!("copy_phasing: --isocon-rule {v:?} (no-pile or via-b)\n{USAGE}"),
    };
    let threads = match kv.get("--threads") {
        None => 4,
        Some(v) => v.parse::<usize>().ok().filter(|&n| n > 0).with_context(|| format!("copy_phasing: --threads {v:?}\n{USAGE}"))?,
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
    use rustle::copy_phasing::family::{isocon_keep, isocon_side, run_family};
    use rustle::copy_phasing::pyjson::{dict_repr, fmt3, Json};
    use rustle::copy_phasing::slots::landing_loci;
    use rustle::copy_phasing::{AS_TIE, MIN_SUPPORT};
    use std::io::Write;
    let by_fam = load_copies(&a.copies)?;
    let (net, n_net) = load_net(&a.net, &a.families)?;
    // the run's read table: every used net read, sorted, so ReadId order = name order
    let table: Vec<String> = net.values().flatten().cloned().collect::<BTreeSet<_>>().into_iter().collect();
    let names: HashMap<Vec<u8>, ReadId> = table.iter().enumerate().map(|(i, n)| (n.as_bytes().to_vec(), i as ReadId)).collect();
    let bam = Bam::open(&a.bam)?;
    let (prim, records) = bam.read_primaries(&names)?;
    let genome = IndexedFasta::open(&a.genome).with_context(|| format!("opening --genome {} (and its .fai)", a.genome))?;
    let iso = match &a.isocon {
        Some(prefix) => isocon_side(prefix, &a.net)?,
        None => BTreeMap::new(),
    };
    let out_dir = std::path::Path::new(&a.out).parent().filter(|d| !d.as_os_str().is_empty()).unwrap_or(std::path::Path::new("."));
    let tmp = out_dir.join(format!(".copy_phasing.{}", std::process::id()));
    std::fs::create_dir_all(&tmp).with_context(|| format!("creating {}", tmp.display()))?;
    let open = |suffix: &str| -> Result<std::io::BufWriter<std::fs::File>> {
        let path = format!("{}.{suffix}", a.out);
        Ok(std::io::BufWriter::new(std::fs::File::create(&path).with_context(|| format!("creating {path}"))?))
    };
    let (mut ct, mut cf, mut cl, mut rd, mut ft) = (open("candidates.tsv")?, open("contigs.fa")?, open("clusters.tsv")?, open("reads.tsv")?, open("families.tsv")?);
    let (mut lh, mut lf, mut lo, mut o2f) = (open("putative.tsv")?, open("putative.fa")?, open("loci.tsv")?, open("o2.tsv")?);
    let mut dump_fh = match &a.dump {
        Some(d) => Some(std::io::BufWriter::new(std::fs::File::create(d).with_context(|| format!("creating {d}"))?)),
        None => None,
    };
    writeln!(ct, "family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used")?;
    writeln!(cl, "family\tcluster\tcandidate\tn_reads\tconsensus_len\tfate\tlinked_to\td")?;
    writeln!(rd, "read\tfamily\tcluster")?;
    writeln!(ft, "family\tn_net\tn_used\tpiles\tpile_reads\tclusters\tconflicting\thaplotypes\tcertified")?;
    writeln!(lh, "family\tcandidate\tcopy\tn_reads\tmulti\tother_holders_multi\tconsensus_len")?;
    writeln!(lo, "family\tlocus\tn_reads\tnear_best\tweaker_only")?;
    writeln!(o2f, "family\tread\tassigned\tn_piles")?;
    let nm = |r: ReadId| table[r as usize].as_str();
    let empty = Vec::new();
    for (fam, reads) in &net {
        let ids: BTreeSet<ReadId> = reads.iter().map(|r| names[r.as_bytes()]).collect();
        let copies = by_fam.get(fam).unwrap_or(&empty);
        let piles = pile_records(&bam, copies, &ids, &prim, &names)?;
        let out = run_family(fam, copies, &piles, &prim, &genome, &a.params, &tmp, dump_fh.is_some().then_some(&table[..]))?;
        let mut o2: Vec<(&str, &str, usize)> = out.o2.iter().map(|(r, c, k)| (nm(*r), c.as_str(), *k)).collect();
        o2.sort();
        for (r, c, k) in o2 {
            writeln!(o2f, "{fam}\t{r}\t{c}\t{k}")?;
        }
        for (k, (x, rs, seq, mi, mc)) in out.looks.iter().enumerate() {
            let mc: Vec<String> = mc.iter().map(|v| fmt3(*v)).collect();
            writeln!(lh, "{fam}\tputative_{fam}_{k}\t{x}\t{}\t{}\t{}\t{}", rs.len(), fmt3(*mi), mc.join(","), seq.len())?;
            writeln!(lf, ">putative_{fam}_{k}\n{}", String::from_utf8_lossy(seq))?;
        }
        let mut recs = Vec::new();
        for &n in &ids {
            let rs = records.get(&n).map_or(&[][..], |v| &v[..]);
            if let Some(b) = rs.iter().map(|x| x.3).max() {
                recs.extend(rs.iter().map(|(ch, s0, e0, asc, _)| (ch.to_string(), *s0, *e0, n, *asc >= b - AS_TIE)));
            }
        }
        let spans: Vec<(String, i64, i64)> = copies.iter().map(|c| (c.chrom.clone(), c.start, c.end)).collect();
        for (ch, s0, e0, nr, nb, so) in landing_loci(&recs, &spans, MIN_SUPPORT) {
            writeln!(lo, "{fam}\t{ch}:{s0}-{e0}\t{nr}\t{nb}\t{so}")?;
        }
        let piled: BTreeSet<String> = out.piled.iter().map(|&r| nm(r).to_string()).collect();
        let side: Vec<(&String, &BTreeSet<String>, &Vec<u8>)> = iso
            .get(fam)
            .map_or(&[][..], |v| &v[..])
            .iter()
            .filter(|(_, rs, _, vb)| isocon_keep(rs, *vb, &piled, a.isocon_via_b))
            .map(|(c, rs, seq, _)| (c, rs, seq))
            .collect();
        if let (Some(fh), Some(mut d)) = (dump_fh.as_mut(), out.dump) {
            d.push((
                "isocon".into(),
                Json::List(
                    side.iter()
                        .map(|(c, rs, _)| Json::Obj(vec![("candidate".into(), Json::str((*c).clone())), ("reads".into(), Json::List(rs.iter().map(|r| Json::str(r.clone())).collect()))]))
                        .collect(),
                ),
            ));
            writeln!(fh, "{}", Json::Obj(d).dumps())?;
        }
        let mut certified: Vec<(String, BTreeSet<String>, &[u8])> =
            out.certified.iter().map(|(h, rs, seq)| (h.clone(), rs.iter().map(|&r| nm(r).to_string()).collect(), &seq[..])).collect();
        certified.extend(side.iter().map(|(c, rs, seq)| (format!("isocon:{c}"), (*rs).clone(), &seq[..])));
        let pc = out.per_copy.values();
        let (s0, s1, s2) = pc.fold((0, 0, 0), |t, v| (t.0 + v.0, t.1 + v.1, t.2 + v.2));
        writeln!(ft, "{fam}\t{}\t{}\t{}\t{s0}\t{s1}\t{s2}\t{}\t{}", n_net[fam], reads.len(), out.per_copy.len(), out.n_h, certified.len())?;
        for (k, (home, rs, seq)) in certified.iter().enumerate() {
            let cand = format!("cand_{fam}_{k}");
            writeln!(ct, "{fam}\t{cand}\t1\t{}\t1\t{}\t{home}\tNA\t{}\t{}", rs.len(), seq.len(), n_net[fam], reads.len())?;
            writeln!(cf, ">{cand}\n{}", String::from_utf8_lossy(seq))?;
            writeln!(cl, "{fam}\t{fam}:h{k}\t{cand}\t{}\t{}\tnew_copy\t-\tNA", rs.len(), seq.len())?;
            for r in rs {
                writeln!(rd, "{r}\t{fam}\t{fam}:h{k}")?;
            }
        }
        eprintln!(
            "[copy_phasing] {fam}: reads {}, piles {}, haplotypes {} {}, certified {}, putative {}",
            reads.len(),
            out.per_copy.len(),
            out.n_h,
            dict_repr(&out.why),
            certified.len(),
            out.looks.len()
        );
    }
    for w in [&mut ct, &mut cf, &mut cl, &mut rd, &mut ft, &mut lh, &mut lf, &mut lo, &mut o2f] {
        w.flush()?;
    }
    if let Some(fh) = dump_fh.as_mut() {
        fh.flush()?;
    }
    let _ = std::fs::remove_dir_all(&tmp);
    let _ = a.threads;
    Ok(())
}
