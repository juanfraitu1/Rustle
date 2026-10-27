//! family_gfa: a graph view of O1 families.
//!
//! Re-encodes what the families stage (`mcl_families`) already computed. Every exonic base of every member locus is a
//! vertex; every aligned column of O1's own all-vs-all locus PAF (`<prefix>.loci.paf`, minimap2 `-c`) whose two bases
//! are both exonic joins two vertices; the joined sets (their transitive closure, as in seqwish) are the graph's
//! columns. Columns traversed together are compacted into nodes, cut at every exon boundary so nodes line up with
//! exons, and each copy is a path spelling its spliced transcript. Fusions between families (`--relations`, from
//! `mcl_families --emit-relations`) get joint graphs.
//!
//! ⚠ Representation only: built FROM O1's output, never used to define or detect families (register L472, L1086).
//! Say "a graph view of an O1 family", never "a family is one VG" (docs/O1_O2_COMPOSITION.md).
//!
//! usage: family_gfa --prefix <out>.fam --outdir DIR [--families MCL1,MCL4] [--relations <out>.fam.relations.tsv]
//!
//! Spec: docs/superpowers/specs/2026-09-30-family-gfa-design.md

use anyhow::{bail, Context, Result};
use clap::Parser;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

// ------------------------------------------------------------------------------------------------------------------
// Inputs
// ------------------------------------------------------------------------------------------------------------------

/// A locus key: chrom, start, end (1-based inclusive). Identical across `.clusters.tsv`, the `.loci.gff3` gene spans
/// and the `.loci.fa` headers (checked on 353/353 clustered loci of the 2026-09-25 chr16 catalog).
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Key {
    chrom: String,
    start: u64,
    end: u64,
}

impl Key {
    fn parse(s: &str) -> Result<Key> {
        let (c, r) = s.rsplit_once(':').with_context(|| format!("bad locus key {s:?}"))?;
        let (a, b) = r.split_once('-').with_context(|| format!("bad locus key {s:?}"))?;
        let k = Key {
            chrom: c.to_string(),
            start: a.parse().with_context(|| format!("bad start in locus key {s:?}"))?,
            end: b.parse().with_context(|| format!("bad end in locus key {s:?}"))?,
        };
        if k.end < k.start {
            bail!("locus key {s:?} ends before it starts");
        }
        Ok(k)
    }

    fn name(&self) -> String {
        format!("{}:{}-{}", self.chrom, self.start, self.end)
    }

    fn len(&self) -> usize {
        (self.end - self.start + 1) as usize
    }
}

/// One locus: gene id and strand from `.loci.gff3`, exons (1-based inclusive, genomic, sorted, overlaps merged) and
/// the genomic sequence of its span from `.loci.fa` (upper case).
#[derive(Clone, Debug)]
struct Locus {
    key: Key,
    gene_id: String,
    strand: char,
    exons: Vec<(u64, u64)>,
    seq: Vec<u8>,
}

impl Locus {
    /// (offset into the span, exon index) of every exonic base, in genomic order.
    fn exonic(&self) -> Vec<(usize, usize)> {
        let mut v = Vec::new();
        for (i, &(a, b)) in self.exons.iter().enumerate() {
            for g in a..=b {
                v.push(((g - self.key.start) as usize, i));
            }
        }
        v
    }

    fn exonic_bp(&self) -> usize {
        self.exons.iter().map(|&(a, b)| (b - a + 1) as usize).sum()
    }

    /// The spliced transcript, 5' to 3': the exonic bases, reverse-complemented for `-` loci.
    fn transcript(&self) -> Vec<u8> {
        let mut t: Vec<u8> = self.exonic().iter().map(|&(o, _)| self.seq[o]).collect();
        if self.strand == '-' {
            t.reverse();
            for b in t.iter_mut() {
                *b = comp(*b);
            }
        }
        t
    }

    fn path_name(&self) -> String {
        format!("{}@{}", self.gene_id, self.key.name())
    }
}

fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

fn field<'a>(f: &[&'a str], i: usize, line: &str) -> Result<&'a str> {
    f.get(i).copied().with_context(|| format!("line has no column {}: {line}", i + 1))
}

/// Family id -> member keys (sorted, unique), from `.clusters.tsv` (`cluster_id ... chrom start end`, one row per member).
fn read_clusters(path: &Path) -> Result<BTreeMap<String, Vec<Key>>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut lines = BufReader::new(f).lines();
    let header = lines.next().with_context(|| format!("{} is empty", path.display()))??;
    let cols: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| {
        cols.iter().position(|c| *c == name).with_context(|| format!("{} has no {name} column", path.display()))
    };
    let (ci, chi, si, ei) = (idx("cluster_id")?, idx("chrom")?, idx("start")?, idx("end")?);
    let mut fams: BTreeMap<String, BTreeSet<Key>> = BTreeMap::new();
    for line in lines {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let key = Key {
            chrom: field(&f, chi, &line)?.to_string(),
            start: field(&f, si, &line)?.parse().with_context(|| format!("bad start: {line}"))?,
            end: field(&f, ei, &line)?.parse().with_context(|| format!("bad end: {line}"))?,
        };
        fams.entry(field(&f, ci, &line)?.to_string()).or_default().insert(key);
    }
    Ok(fams.into_iter().map(|(k, v)| (k, v.into_iter().collect())).collect())
}

/// Gene span -> (gene id, strand, exons) from `.loci.gff3`: `gene` rows keyed by span, `exon` rows by `Parent`.
/// Exons are sorted and overlapping ones merged.
fn read_gff3(path: &Path) -> Result<HashMap<Key, (String, char, Vec<(u64, u64)>)>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut gene_key: HashMap<String, Key> = HashMap::new();
    let mut out: HashMap<Key, (String, char, Vec<(u64, u64)>)> = HashMap::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let attr = |k: &str| {
            let pat = format!("{k}=");
            f[8].split(';').find_map(|kv| kv.strip_prefix(pat.as_str()).map(|s| s.to_string()))
        };
        let (a, b): (u64, u64) = (
            f[3].parse().with_context(|| format!("bad start: {line}"))?,
            f[4].parse().with_context(|| format!("bad end: {line}"))?,
        );
        match f[2] {
            "gene" => {
                let id = attr("ID").with_context(|| format!("gene row without ID: {line}"))?;
                let gid = id.strip_prefix("gene-").unwrap_or(&id).to_string();
                let key = Key { chrom: f[0].to_string(), start: a, end: b };
                let strand = f[6].chars().next().unwrap_or('.');
                if gene_key.insert(gid.clone(), key.clone()).is_some() {
                    bail!("gene id {gid} appears twice in {}", path.display());
                }
                out.insert(key, (gid, strand, Vec::new()));
            }
            "exon" => {
                let parent = attr("Parent").with_context(|| format!("exon row without Parent: {line}"))?;
                let gid = parent.strip_prefix("gene-").unwrap_or(&parent).to_string();
                let key = gene_key.get(&gid).with_context(|| format!("exon of unknown gene {gid} in {}", path.display()))?;
                out.get_mut(key).expect("gene inserted with its key").2.push((a, b));
            }
            _ => {}
        }
    }
    for (_, (_, _, ex)) in out.iter_mut() {
        ex.sort();
        let mut merged: Vec<(u64, u64)> = Vec::new();
        for &(a, b) in ex.iter() {
            match merged.last_mut() {
                Some(last) if a <= last.1 => last.1 = last.1.max(b),
                _ => merged.push((a, b)),
            }
        }
        *ex = merged;
    }
    Ok(out)
}

/// Sequences of the wanted keys from `.loci.fa` (headers `chrom:start-end`), upper-cased (soft-masked bases count).
fn read_fasta(path: &Path, want: &HashSet<Key>) -> Result<HashMap<Key, Vec<u8>>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut out = HashMap::new();
    let mut cur: Option<Key> = None;
    let mut buf: Vec<u8> = Vec::new();
    for line in BufReader::new(f).lines() {
        let line = line?;
        if let Some(h) = line.strip_prefix('>') {
            if let Some(k) = cur.take() {
                out.insert(k, std::mem::take(&mut buf));
            }
            buf.clear();
            let k = Key::parse(h.split_whitespace().next().unwrap_or(""))?;
            cur = if want.contains(&k) { Some(k) } else { None };
        } else if cur.is_some() {
            buf.extend(line.trim_end().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(k) = cur.take() {
        out.insert(k, buf);
    }
    Ok(out)
}

/// One all-vs-all record: 0-based half-open coordinates on each locus span, strand, and the `cg:Z` CIGAR.
#[derive(Clone, Debug)]
struct Paf {
    q: Key,
    qs: usize,
    qe: usize,
    rev: bool,
    t: Key,
    ts: usize,
    te: usize,
    cigar: Vec<(u32, u8)>,
}

fn parse_cigar(s: &str) -> Result<Vec<(u32, u8)>> {
    let mut out = Vec::new();
    let mut n: u32 = 0;
    let mut have = false;
    for c in s.bytes() {
        if c.is_ascii_digit() {
            n = n
                .checked_mul(10)
                .and_then(|x| x.checked_add((c - b'0') as u32))
                .with_context(|| format!("CIGAR length overflow in {s:?}"))?;
            have = true;
        } else {
            if !have {
                bail!("bad CIGAR {s:?}");
            }
            out.push((n, c));
            n = 0;
            have = false;
        }
    }
    if have {
        bail!("bad CIGAR {s:?}");
    }
    Ok(out)
}

/// Records of `.loci.paf` whose query and target are both wanted and differ. Lines that do not start with two locus
/// keys are ignored; records without `cg:Z` are skipped and counted (the second value).
fn read_paf(path: &Path, want: &HashSet<Key>) -> Result<(Vec<Paf>, usize)> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut out = Vec::new();
    let mut no_cigar = 0;
    for line in BufReader::new(f).lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let (Ok(q), Ok(t)) = (Key::parse(f[0]), Key::parse(f[5])) else { continue };
        if q == t || !want.contains(&q) || !want.contains(&t) {
            continue;
        }
        let Some(cg) = f[12..].iter().find_map(|x| x.strip_prefix("cg:Z:")) else {
            no_cigar += 1;
            continue;
        };
        out.push(Paf {
            q,
            qs: f[2].parse().with_context(|| format!("bad query start: {line}"))?,
            qe: f[3].parse().with_context(|| format!("bad query end: {line}"))?,
            rev: f[4] == "-",
            t,
            ts: f[7].parse().with_context(|| format!("bad target start: {line}"))?,
            te: f[8].parse().with_context(|| format!("bad target end: {line}"))?,
            cigar: parse_cigar(cg)?,
        });
    }
    Ok((out, no_cigar))
}

/// Everything read from one families-stage prefix: all families (for membership), the selected ones, and the loci
/// and alignments they need.
struct Inputs {
    all: BTreeMap<String, Vec<Key>>,
    selected: Vec<String>,
    member_of: HashMap<Key, String>,
    loci: HashMap<Key, Locus>,
    pafs: Vec<Paf>,
    no_cigar: usize,
}

/// Loads the selected families (all when `select` is empty) and the `extra` keys (relation units) with the members of
/// their families. A family member missing from any file is an error naming it; an extra key missing from the loci
/// files is skipped (fusions report it).
fn load(prefix: &str, select: &[String], extra: &[Key]) -> Result<Inputs> {
    let clusters = format!("{prefix}.clusters.tsv");
    let all = read_clusters(Path::new(&clusters))?;
    let mut member_of: HashMap<Key, String> = HashMap::new();
    for (f, ks) in &all {
        for k in ks {
            if let Some(prev) = member_of.insert(k.clone(), f.clone()) {
                bail!("locus {} is in two families ({prev} and {f}) in {clusters}", k.name());
            }
        }
    }
    let selected: Vec<String> = if select.is_empty() {
        all.keys().cloned().collect()
    } else {
        for f in select {
            if !all.contains_key(f) {
                bail!("family {f} is not in {clusters}");
            }
        }
        select.to_vec()
    };
    // who needs each key, for error messages; None = a relation unit (optional)
    let mut owner: BTreeMap<Key, Option<String>> = BTreeMap::new();
    for f in &selected {
        for k in &all[f] {
            owner.entry(k.clone()).or_insert_with(|| Some(f.clone()));
        }
    }
    for k in extra {
        if let Some(f) = member_of.get(k) {
            for m in &all[f] {
                owner.entry(m.clone()).or_insert_with(|| Some(f.clone()));
            }
        }
        owner.entry(k.clone()).or_insert(None);
    }
    let want: HashSet<Key> = owner.keys().cloned().collect();
    let gff_path = format!("{prefix}.loci.gff3");
    let gff = read_gff3(Path::new(&gff_path))?;
    let fa_path = format!("{prefix}.loci.fa");
    let fa = read_fasta(Path::new(&fa_path), &want)?;
    let mut loci = HashMap::new();
    for (k, who) in &owner {
        let (Some(g), Some(seq)) = (gff.get(k), fa.get(k)) else {
            match who {
                None => continue,
                Some(f) => {
                    let file = if gff.contains_key(k) { &fa_path } else { &gff_path };
                    bail!("locus {} of family {f} is missing from {file}", k.name());
                }
            }
        };
        let (gid, strand, exons) = g.clone();
        if seq.len() != k.len() {
            bail!("locus {}: {fa_path} has {} bases for a {} bp span", k.name(), seq.len(), k.len());
        }
        if exons.is_empty() {
            bail!("locus {} has no exons in {gff_path}", k.name());
        }
        for &(a, b) in &exons {
            if a < k.start || b > k.end {
                bail!("locus {}: exon {a}-{b} lies outside its span", k.name());
            }
        }
        loci.insert(k.clone(), Locus { key: k.clone(), gene_id: gid, strand, exons, seq: seq.clone() });
    }
    let (pafs, no_cigar) = read_paf(Path::new(&format!("{prefix}.loci.paf")), &want)?;
    Ok(Inputs { all, selected, member_of, loci, pafs, no_cigar })
}

// ------------------------------------------------------------------------------------------------------------------
// Command line
// ------------------------------------------------------------------------------------------------------------------

#[derive(Parser, Debug)]
#[command(about = "Graph view of O1 families: an exon-block GFA per family (and per fusion across families)")]
struct Args {
    /// Families-stage prefix: reads PREFIX.clusters.tsv, .loci.gff3, .loci.fa and .loci.paf.
    #[arg(long)]
    prefix: String,
    /// Output directory (created if missing).
    #[arg(long)]
    outdir: PathBuf,
    /// Comma-separated family ids (default: every family in PREFIX.clusters.tsv).
    #[arg(long, value_delimiter = ',')]
    families: Vec<String>,
}

fn main() -> Result<()> {
    let a = Args::parse();
    let inp = load(&a.prefix, &a.families, &[])?;
    eprintln!(
        "family_gfa: {} families selected, {} loci, {} alignment records ({} without cg:Z skipped)",
        inp.selected.len(),
        inp.loci.len(),
        inp.pafs.len(),
        inp.no_cigar
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmpfile(dir: &tempfile::TempDir, name: &str, text: &str) -> PathBuf {
        let p = dir.path().join(name);
        std::fs::write(&p, text).unwrap();
        p
    }

    #[test]
    fn key_parses_and_names() {
        let k = Key::parse("chr16:3546-7021").unwrap();
        assert_eq!((k.chrom.as_str(), k.start, k.end, k.len()), ("chr16", 3546, 7021, 3476));
        assert_eq!(k.name(), "chr16:3546-7021");
        assert!(Key::parse("chr16:7021-3546").is_err());
        assert!(Key::parse("chr16").is_err());
    }

    #[test]
    fn cigar_parses_and_rejects_garbage() {
        assert_eq!(parse_cigar("908M319I50M").unwrap(), vec![(908, b'M'), (319, b'I'), (50, b'M')]);
        assert!(parse_cigar("M5").is_err());
        assert!(parse_cigar("5M3").is_err());
    }

    #[test]
    fn gff3_merges_overlapping_exons_and_keeps_strand() {
        let d = tempfile::tempdir().unwrap();
        let p = tmpfile(
            &d,
            "x.gff3",
            "##gff-version 3\nc1\t.\tgene\t100\t200\t.\t-\t.\tID=gene-G1;Name=G1\n\
             c1\t.\texon\t150\t170\t.\t-\t.\tParent=gene-G1;gene=G1\n\
             c1\t.\texon\t100\t120\t.\t-\t.\tParent=gene-G1;gene=G1\n\
             c1\t.\texon\t110\t130\t.\t-\t.\tParent=gene-G1;gene=G1\n",
        );
        let m = read_gff3(&p).unwrap();
        let (gid, strand, ex) = &m[&Key::parse("c1:100-200").unwrap()];
        assert_eq!((gid.as_str(), *strand), ("G1", '-'));
        assert_eq!(ex, &vec![(100, 130), (150, 170)]);
    }

    #[test]
    fn fasta_uppercases_and_keeps_only_wanted() {
        let d = tempfile::tempdir().unwrap();
        let p = tmpfile(&d, "x.fa", ">c1:1-6\nacgT\nNn\n>c1:10-12\nGGG\n");
        let want: HashSet<Key> = [Key::parse("c1:1-6").unwrap()].into_iter().collect();
        let m = read_fasta(&p, &want).unwrap();
        assert_eq!(m.len(), 1);
        assert_eq!(m[&Key::parse("c1:1-6").unwrap()], b"ACGTNN".to_vec());
    }

    #[test]
    fn transcript_is_reverse_complemented_on_minus() {
        let l = Locus {
            key: Key::parse("c1:1-8").unwrap(),
            gene_id: "G".into(),
            strand: '-',
            exons: vec![(1, 2), (7, 8)],
            seq: b"ACGTTTGG".to_vec(),
        };
        assert_eq!(l.exonic_bp(), 4);
        assert_eq!(l.transcript(), b"CCGT".to_vec()); // exonic ACGG -> reverse complement CCGT
        assert_eq!(l.path_name(), "G@c1:1-8");
    }
}
