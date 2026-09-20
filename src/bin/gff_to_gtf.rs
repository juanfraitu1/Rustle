//! Extract one chromosome from a RefSeq GFF3 and write a gffread-style GTF (transcript + exon rows).
//!
//! Native Rust replacement for `tools/refseq_gff_to_gtf.py` (§6r9), which it reproduces byte for byte.
//! Stands in for `gffread -T`, which is not installed on the work machine, and is how every
//! per-chromosome reference GTF in `docs/DATA.md` is built. Validated against gffread's own
//! `chr20_ref.gtf` at 4,574 == 4,574 transcripts.
//!
//! Transcripts are emitted in FIRST-SEEN order (Python dicts preserve insertion order, and the byte
//! comparison against the Python output depends on it), and each transcript's exons are sorted by
//! `(start, end, source, strand)` exactly as the Python tuple sort did.
//!
//! usage: gff_to_gtf REFSEQ.gff[.gz] CHROM OUT.gtf

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};

/// Parse a GFF3 attribute column into key -> value. Keys and values are trimmed; a trailing `;` is
/// ignored; entries without `=` are skipped.
fn attrs(s: &str) -> HashMap<&str, &str> {
    let mut d = HashMap::new();
    for kv in s.trim_end_matches(';').split(';') {
        if let Some((k, v)) = kv.split_once('=') {
            d.insert(k.trim(), v.trim());
        }
    }
    d
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: {} REFSEQ.gff[.gz] CHROM OUT.gtf", args[0]);
        std::process::exit(2);
    }
    let (src, chrom, out) = (&args[1], &args[2], &args[3]);

    let file = std::fs::File::open(src).with_context(|| format!("opening {src}"))?;
    let reader: Box<dyn BufRead> = if src.ends_with(".gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    // exons per transcript id, plus the first-seen order of those ids
    let mut exons: HashMap<String, Vec<(u64, u64, String, String)>> = HashMap::new();
    let mut order: Vec<String> = Vec::new();
    // transcript id -> (gene id, gene name)
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
            let (Ok(s), Ok(e)) = (f[3].parse::<u64>(), f[4].parse::<u64>()) else { continue };
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
                let gene = a.get("gene").or_else(|| a.get("Name")).copied().unwrap_or("");
                info.insert(id.to_string(), (a.get("Parent").copied().unwrap_or("").to_string(), gene.to_string()));
            }
        }
    }

    let mut fo = std::io::BufWriter::new(std::fs::File::create(out).with_context(|| format!("creating {out}"))?);
    let mut n = 0usize;
    for tid in &order {
        let ex = exons.get_mut(tid).expect("ordered id must have exons");
        ex.sort();
        let (gid, gname) = info.get(tid).cloned().unwrap_or_default();
        let (src_f, strand) = (ex[0].2.clone(), ex[0].3.clone());
        let at = format!("transcript_id \"{tid}\"; gene_id \"{gid}\"; gene_name \"{gname}\"");
        writeln!(fo, "{chrom}\t{src_f}\ttranscript\t{}\t{}\t.\t{strand}\t.\t{at}", ex[0].0, ex[ex.len() - 1].1)?;
        for (i, (s, e, _, _)) in ex.iter().enumerate() {
            writeln!(fo, "{chrom}\t{src_f}\texon\t{s}\t{e}\t.\t{strand}\t.\t{at}; exon_number \"{}\";", i + 1)?;
        }
        n += 1;
    }
    fo.flush()?;
    eprintln!("{chrom}: {n} transcripts");
    Ok(())
}
