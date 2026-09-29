//! The CONTAINER of a family member's extra pieces (its ACCESSORY exon blocks) and their relations to other
//! families: the Rust port of `bench/family_container.py` (frozen sha1 e197ccb3, the binding definition of
//! `docs/PREREG_fusion_container_sim_2026-09-28.md` §1 + Amendment 1), run by `mcl_families --from-gtf
//! --emit-container` after the families are written. It reads the families products and never changes a family.
//!
//! **STATUS:** OPT-IN  (docs/MODULE_STATUS.md; `mcl_families --emit-container`, default off; driver `RUSTLE_FAMILY_CONTAINER=1`)
//!
//! Definition (prereg §1, Amendment 1):
//! * LOCUS m = a member row of `clusters.tsv` (a graph node key `CONTIG:START-END`) plus every annotation record that
//!   `loci.tsv` folds into it; FAMILY F = its `cluster_id`. A record is one assembled `gene_id`: `loci.gff3`'s `gene`
//!   line gives its key and `Name=` its gene_id (two gene_ids with the same span share one key and are both taken).
//! * EXON BLOCKS of m = the union of the exons of ALL transcripts of ALL gene_ids of m's records, merged where they
//!   OVERLAP (share >= 1 base; abutting exons stay separate).
//! * Block b of m is CORE iff some PAF record between a record of m and a record of another member m' of F has an
//!   aligned column (CIGAR M/=/X) whose m-side base lies in b and whose m'-side base is an exon base of m' (m''s own
//!   all-transcript blocks). Every PAF record counts, whatever its identity, length or primary flag. ACCESSORY = not
//!   core.
//! * RELATION: each accessory block gets the same test against the members of every OTHER family F'; the family
//!   relation F -> F' exists iff some member of F carries such a block; `reciprocal` says whether F' -> F exists too.
//!   Core blocks are not tested for relations (their relation columns are `.`).
//! * Unclustered loci get no rows and are never partners. Records between two records of the same locus are skipped.
//!
//! Coordinates: loci.fa / PAF names are record keys `CONTIG:START-END` (GFF 1-based closed) and a record's sequence is
//! the genome's FORWARD strand from START to END, so PAF offset o is genome base START + o (0-based START - 1 + o).
//! PAF `+`: the CIGAR walks target [ts,te) and query [qs,qe) ascending; `-`: target ascending against the reverse
//! complement of query [qs,qe), i.e. forward query offsets from qe-1 downwards. M/=/X consume both sides, I the query,
//! D/N the target; any other op, a CIGAR whose lengths disagree with the PAF columns, or a projected record without
//! `cg:Z` is an error. Outputs are GFF 1-based closed, blocks numbered in GENOMIC order.
//!
//! Outputs (`write`): `<out>.container.tsv` (one row per clustered locus x exon block), `<out>.container_relations.tsv`
//! (one row per directed family relation) and `<out>.container_summary.tsv` (counts): byte for byte the frozen
//! script's `OUT.blocks.tsv`, `OUT.relations.tsv` and `OUT.summary.tsv`, including its input conventions (the first
//! `key "` occurrence of a GTF attribute, an empty `gene_id` read as the transcript id, the fold table's
//! last-write-wins value at the first-write position, Python's `int()` on the decimal fields). Only a lone `\r` line
//! break (which Python's universal newlines would split on) is not reproduced.
use anyhow::{bail, Context, Result};
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

/// A record / locus key `(CONTIG, START, END)`, GFF 1-based closed.
pub type Key = (String, i64, i64);

pub const BLOCK_HEADER: [&str; 17] = [
    "family_id",
    "locus",
    "chrom",
    "strand",
    "gene_ids",
    "n_records",
    "block",
    "n_blocks",
    "start",
    "end",
    "bp",
    "class",
    "core_bp",
    "core_partners",
    "rel_families",
    "rel_bp",
    "rel_partners",
];
pub const REL_HEADER: [&str; 7] = ["family_id", "related_family", "n_members", "n_blocks", "bp", "reciprocal", "members"];
pub const SUMMARY_KEYS: [&str; 26] = [
    "families",
    "loci",
    "folded_records",
    "loci_gff3_keys",
    "gene_key_collisions",
    "records_with_key_collision",
    "records_span_mismatch",
    "blocks",
    "core_blocks",
    "accessory_blocks",
    "core_block_bp",
    "accessory_bp",
    "accessory_blocks_related",
    "loci_all_core",
    "loci_with_accessory",
    "families_with_relation",
    "family_relations_directed",
    "family_relations_reciprocal",
    "paf_records",
    "paf_malformed",
    "paf_self",
    "paf_unclustered",
    "paf_same_locus",
    "paf_no_exon_interval",
    "paf_projected",
    "paf_exon_exon",
];

pub fn key_str(k: &Key) -> String {
    format!("{}:{}-{}", k.0, k.1, k.2)
}

/// Python's `int(s)` on a decimal field: surrounding whitespace, an optional sign, ASCII digits with single `_`
/// between them. `None` where Python raises `ValueError` (and beyond i64).
pub fn py_int(s: &str) -> Option<i64> {
    let t = s.trim();
    let (neg, body) = match t.as_bytes().first() {
        Some(b'+') => (false, &t[1..]),
        Some(b'-') => (true, &t[1..]),
        _ => (false, t),
    };
    if body.is_empty() || body.starts_with('_') || body.ends_with('_') || body.contains("__") {
        return None;
    }
    let digits: String = body.chars().filter(|&c| c != '_').collect();
    if !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let v: i64 = digits.parse().ok()?;
    Some(if neg { -v } else { v })
}

fn int_field(s: &str, what: &str) -> Result<i64> {
    py_int(s).with_context(|| format!("invalid literal for int(): {s:?} ({what})"))
}

/// `CONTIG:START-END` -> key (the LAST `:` splits the contig, the FIRST `-` after it the range); `None` if malformed.
pub fn parse_key(name: &str) -> Option<Key> {
    let (c, r) = name.rsplit_once(':')?;
    let (a, b) = r.split_once('-')?;
    Some((c.to_string(), py_int(a)?, py_int(b)?))
}

/// The value of `key "..."` in a GTF attribute column: the FIRST occurrence of `key "` (the same substring rule as
/// `mcl_families::gtf_loci`), up to the next `"`.
fn gtf_attr<'a>(s: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("{key} \"");
    let i = s.find(&pat)? + pat.len();
    let j = s[i..].find('"')? + i;
    Some(&s[i..j])
}

/// Merge 0-based half-open intervals that OVERLAP (share >= 1 base); abutting intervals stay separate.
pub fn merge_blocks(mut iv: Vec<(i64, i64)>) -> Vec<(i64, i64)> {
    iv.sort();
    let mut out: Vec<(i64, i64)> = Vec::new();
    for (s, e) in iv {
        if let Some(last) = out.last_mut() {
            if s < last.1 {
                if e > last.1 {
                    last.1 = e;
                }
                continue;
            }
        }
        out.push((s, e));
    }
    out
}

/// Bases covered by 0-based half-open intervals (overlapping or abutting).
pub fn union_len<'a, I: IntoIterator<Item = &'a (i64, i64)>>(iv: I) -> i64 {
    let mut v: Vec<(i64, i64)> = iv.into_iter().copied().collect();
    v.sort();
    let mut n = 0;
    let mut cur: Option<(i64, i64)> = None;
    for (s, e) in v {
        cur = match cur {
            None => Some((s, e)),
            Some((cs, ce)) if s > ce => {
                n += ce - cs;
                Some((s, e))
            }
            Some((cs, ce)) => Some((cs, ce.max(e))),
        };
    }
    if let Some((cs, ce)) = cur {
        n += ce - cs;
    }
    n
}

/// `(\d+)([MIDNSHP=X])` tokens whose concatenation is the whole CIGAR, else an error.
fn parse_cigar(cigar: &str) -> Result<Vec<(i64, u8)>> {
    let b = cigar.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let j = i;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
        }
        if i == j || i >= b.len() || !b"MIDNSHP=X".contains(&b[i]) {
            bail!("unparseable CIGAR {:?}", &cigar[..cigar.len().min(60)]);
        }
        let n: i64 = cigar[j..i].parse().with_context(|| format!("CIGAR length {:?}", &cigar[j..i]))?;
        out.push((n, b[i]));
        i += 1;
    }
    Ok(out)
}

/// The aligned runs of one PAF record, in the record's own offset frame: `(t0, q0, n)` = n aligned columns; column i
/// joins target offset t0 + i with query offset q0 + i on `+`, and with q0 + n - 1 - i on `-` (so [q0, q0 + n) is the
/// run's query interval on both strands).
pub fn aligned_runs(cigar: &str, strand: &str, qs: i64, qe: i64, ts: i64, te: i64) -> Result<Vec<(i64, i64, i64)>> {
    let ops = parse_cigar(cigar)?;
    if strand != "+" && strand != "-" {
        bail!("PAF strand {strand:?}");
    }
    let fwd = strand == "+";
    let (mut t, mut qc) = (ts, 0i64);
    let mut out = Vec::new();
    for (n, op) in ops {
        match op {
            b'M' | b'=' | b'X' => {
                out.push((t, if fwd { qs + qc } else { qe - qc - n }, n));
                t += n;
                qc += n;
            }
            b'I' => qc += n,
            b'D' | b'N' => t += n,
            _ => bail!("CIGAR op {} not allowed in a PAF record", op as char),
        }
    }
    if t != te || qc != qe - qs {
        bail!("CIGAR consumes target {} / query {qc} but the record spans {} / {}", t - ts, te - ts, qe - qs);
    }
    Ok(out)
}

/// `(lo, hi, block_index)` for every block of a locus (sorted, disjoint `starts`/`ends`) intersecting [a, b).
fn mask_ranges(starts: &[i64], ends: &[i64], a: i64, b: i64) -> Vec<(i64, i64, usize)> {
    let mut out = Vec::new();
    let mut i = ends.partition_point(|&e| e <= a);
    while i < starts.len() && starts[i] < b {
        let (lo, hi) = (starts[i].max(a), ends[i].min(b));
        if hi > lo {
            out.push((lo, hi, i));
        }
        i += 1;
    }
    out
}

fn overlaps_any(starts: &[i64], ends: &[i64], a: i64, b: i64) -> bool {
    let i = ends.partition_point(|&e| e <= a);
    i < starts.len() && starts[i] < b
}

/// One exon-exon stretch of a record: `(q_lo, q_hi, q_block, t_lo, t_hi, t_block)` in genome coordinates; the two
/// intervals have the same length and are joined column by column (reversed on `-`).
pub type ExonColumns = (i64, i64, usize, i64, i64, usize);

/// Aligned columns of one record joining an exon base of the query locus to an exon base of the target locus.
/// `q_off` / `t_off`: the genome 0-based coordinate of offset 0 of the query / target record; `q_blocks` / `t_blocks`:
/// `(starts, ends)` of the query / target LOCUS blocks, genome 0-based half-open.
#[allow(clippy::too_many_arguments)]
pub fn exon_columns(
    cigar: &str,
    strand: &str,
    qs: i64,
    qe: i64,
    ts: i64,
    te: i64,
    q_off: i64,
    t_off: i64,
    q_blocks: (&[i64], &[i64]),
    t_blocks: (&[i64], &[i64]),
) -> Result<Vec<ExonColumns>> {
    let (qst, qen) = q_blocks;
    let (tst, ten) = t_blocks;
    let fwd = strand == "+";
    let mut out = Vec::new();
    for (t0, q0, n) in aligned_runs(cigar, strand, qs, qe, ts, te)? {
        let (tg, qg) = (t_off + t0, q_off + q0);
        let tr: Vec<(i64, i64, usize)> =
            mask_ranges(tst, ten, tg, tg + n).into_iter().map(|(lo, hi, bi)| (lo - tg, hi - tg, bi)).collect();
        if tr.is_empty() {
            continue;
        }
        let qm = mask_ranges(qst, qen, qg, qg + n);
        if qm.is_empty() {
            continue;
        }
        let qr: Vec<(i64, i64, usize)> = if fwd {
            qm.into_iter().map(|(lo, hi, bi)| (lo - qg, hi - qg, bi)).collect()
        } else {
            // column i <-> query qg + n - 1 - i, so query [lo, hi) <-> i in [qg + n - hi, qg + n - lo)
            qm.into_iter().rev().map(|(lo, hi, bi)| (qg + n - hi, qg + n - lo, bi)).collect()
        };
        let (mut x, mut y) = (0, 0);
        while x < qr.len() && y < tr.len() {
            let (lo, hi) = (qr[x].0.max(tr[y].0), qr[x].1.min(tr[y].1));
            if hi > lo {
                let (ql, qh) = if fwd { (qg + lo, qg + hi) } else { (qg + n - hi, qg + n - lo) };
                out.push((ql, qh, qr[x].2, tg + lo, tg + hi, tr[y].2));
            }
            if qr[x].1 <= tr[y].1 {
                x += 1;
            } else {
                y += 1;
            }
        }
    }
    Ok(out)
}

// ---------------------------------------------------------------------------------------------------------- inputs

/// A text input, gzip-decoded when the path ends in `.gz` (as the script's `open_text`).
pub fn open_text(path: &str) -> Result<Box<dyn BufRead>> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {path}"))?;
    Ok(if path.ends_with(".gz") {
        Box::new(std::io::BufReader::with_capacity(1 << 20, flate2::read::MultiGzDecoder::new(f)))
    } else {
        Box::new(std::io::BufReader::with_capacity(1 << 20, f))
    })
}

/// `clusters.tsv`: the family of each member key, the members per family in file order (families in
/// first-appearance order), and the member keys in first-appearance order.
struct Clusters {
    fam_of: HashMap<Key, usize>,
    fam_ids: Vec<String>,
    members: Vec<Vec<Key>>,
    keys: Vec<Key>,
}

fn read_clusters(r: &mut dyn BufRead, path: &str) -> Result<Clusters> {
    let mut lines = r.lines();
    let header_line = lines.next().transpose()?.unwrap_or_default();
    let header: Vec<&str> = header_line.split('\t').collect();
    let mut col: HashMap<&str, usize> = HashMap::new();
    for (i, h) in header.iter().enumerate() {
        col.insert(h, i);
    }
    for need in ["cluster_id", "chrom", "start", "end"] {
        if !col.contains_key(need) {
            bail!("{path}: no {need} column (header {header:?})");
        }
    }
    let (ci, cc, cs, ce) = (col["cluster_id"], col["chrom"], col["start"], col["end"]);
    let mut out = Clusters { fam_of: HashMap::new(), fam_ids: Vec::new(), members: Vec::new(), keys: Vec::new() };
    let mut order: HashMap<String, usize> = HashMap::new();
    for line in lines {
        let line = line?;
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < header.len() {
            continue;
        }
        let fid = r[ci];
        let key: Key = (r[cc].to_string(), int_field(r[cs], "clusters start")?, int_field(r[ce], "clusters end")?);
        if let Some(&f) = out.fam_of.get(&key) {
            if out.fam_ids[f] != fid {
                bail!("{path}: {} is in {} and {fid} (not a strict partition)", key_str(&key), out.fam_ids[f]);
            }
            continue;
        }
        let f = match order.get(fid) {
            Some(&f) => f,
            None => {
                order.insert(fid.to_string(), out.fam_ids.len());
                out.fam_ids.push(fid.to_string());
                out.members.push(Vec::new());
                out.fam_ids.len() - 1
            }
        };
        out.fam_of.insert(key.clone(), f);
        out.members[f].push(key.clone());
        out.keys.push(key);
    }
    Ok(out)
}

/// `loci.tsv` (`annotation representative`) -> the folds `annotation -> representative` (a != b), in first-insertion
/// order with the last value (a Python dict's semantics).
fn read_folds(r: &mut dyn BufRead, path: &str) -> Result<Vec<(Key, Key)>> {
    let mut lines = r.lines();
    let header_line = lines.next().transpose()?.unwrap_or_default();
    let header: Vec<&str> = header_line.split('\t').collect();
    if header.len() < 2 || header[0] != "annotation" || header[1] != "representative" {
        bail!("{path}: header {header:?} is not `annotation representative`");
    }
    let mut folds: Vec<(Key, Key)> = Vec::new();
    let mut at: HashMap<Key, usize> = HashMap::new();
    for line in lines {
        let line = line?;
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < 2 {
            continue;
        }
        let (Some(a), Some(b)) = (parse_key(r[0]), parse_key(r[1])) else {
            bail!("{path}: malformed keys {:?}", &r[..2]);
        };
        if a != b {
            match at.get(&a) {
                Some(&i) => folds[i].1 = b,
                None => {
                    at.insert(a.clone(), folds.len());
                    folds.push((a, b));
                }
            }
        }
    }
    Ok(folds)
}

/// `loci.gff3` `gene` lines -> ({key: [gene_id, ...]}, {key: strand of its first gene line}).
fn read_loci_gff3(r: &mut dyn BufRead, path: &str) -> Result<(HashMap<Key, Vec<String>>, HashMap<Key, String>)> {
    let mut genes: HashMap<Key, Vec<String>> = HashMap::new();
    let mut strand: HashMap<Key, String> = HashMap::new();
    for line in r.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < 9 || r[2] != "gene" {
            continue;
        }
        let key: Key = (r[0].to_string(), int_field(r[3], "gff3 start")?, int_field(r[4], "gff3 end")?);
        let Some(name) = r[8].split(';').find_map(|kv| kv.strip_prefix("Name=")) else {
            bail!("{path}: gene line without Name=: {}", line.trim());
        };
        genes.entry(key.clone()).or_default().push(name.to_string());
        strand.entry(key).or_insert_with(|| r[6].to_string());
    }
    Ok((genes, strand))
}

/// Assembled GTF -> {gene_id: [(chrom, start1, end1), ...]} over all its transcripts, for the wanted gene_ids.
/// Transcript -> gene comes from `transcript` lines and exons join by `transcript_id`, as `mcl_families::gtf_loci`.
fn read_gtf(r: &mut dyn BufRead, wanted: &HashSet<String>) -> Result<HashMap<String, Vec<(String, i64, i64)>>> {
    let mut gene_of: HashMap<String, String> = HashMap::new();
    let mut exons: HashMap<String, Vec<(String, i64, i64)>> = HashMap::new();
    for line in r.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let r: Vec<&str> = line.split('\t').collect();
        if r.len() < 9 {
            continue;
        }
        let Some(t) = gtf_attr(r[8], "transcript_id") else { continue };
        if r[2] == "transcript" {
            let g = gtf_attr(r[8], "gene_id").filter(|g| !g.is_empty()).unwrap_or(t);
            if wanted.contains(g) {
                gene_of.insert(t.to_string(), g.to_string());
            }
        } else if r[2] == "exon" {
            let ex = (r[0].to_string(), int_field(r[3], "gtf exon start")?, int_field(r[4], "gtf exon end")?);
            exons.entry(t.to_string()).or_default().push(ex);
        }
    }
    let mut out: HashMap<String, Vec<(String, i64, i64)>> = HashMap::new();
    for (t, g) in gene_of {
        let v = out.entry(g).or_default();
        if let Some(ex) = exons.get(&t) {
            v.extend(ex.iter().cloned());
        }
    }
    Ok(out)
}

// ------------------------------------------------------------------------------------------------------------ core

/// One clustered locus: its records (the member key first, then the folded ones), gene_ids, exon blocks (genome
/// 0-based half-open, sorted, disjoint) and family.
struct Locus {
    key: Key,
    family: usize,
    records: Vec<Key>,
    gene_ids: Vec<String>,
    starts: Vec<i64>,
    ends: Vec<i64>,
}

#[derive(Default)]
struct Counts(HashMap<&'static str, i64>);
impl Counts {
    fn add(&mut self, k: &'static str, n: i64) {
        *self.0.entry(k).or_insert(0) += n;
    }
    fn inc(&mut self, k: &'static str) {
        self.add(k, 1);
    }
    fn get(&self, k: &str) -> i64 {
        self.0.get(k).copied().unwrap_or(0)
    }
}

fn build_loci(
    cl: &Clusters,
    folds: &[(Key, Key)],
    loci_genes: &HashMap<Key, Vec<String>>,
    gene_exons: &HashMap<String, Vec<(String, i64, i64)>>,
    cnt: &mut Counts,
) -> Result<(Vec<Locus>, HashMap<Key, usize>)> {
    let idx_of: HashMap<&Key, usize> = cl.keys.iter().enumerate().map(|(i, k)| (k, i)).collect();
    let mut records: Vec<Vec<Key>> = cl.keys.iter().map(|m| vec![m.clone()]).collect();
    for (ann, rep) in folds {
        if let Some(&i) = idx_of.get(rep) {
            if cl.fam_of.contains_key(ann) {
                bail!("loci.tsv folds {} into {} but it is itself a cluster member", key_str(ann), key_str(rep));
            }
            records[i].push(ann.clone());
            cnt.inc("folded_records");
        }
    }
    let mut locus_of_record: HashMap<Key, usize> = HashMap::new();
    let mut loci = Vec::with_capacity(records.len());
    for (i, recs) in records.into_iter().enumerate() {
        let m = &cl.keys[i];
        let mut gids: Vec<String> = Vec::new();
        let mut exs: Vec<(i64, i64)> = Vec::new();
        for rec in &recs {
            if locus_of_record.contains_key(rec) {
                bail!("record {} belongs to two loci", key_str(rec));
            }
            locus_of_record.insert(rec.clone(), i);
            let Some(names) = loci_genes.get(rec) else {
                bail!("{} has no gene line in loci.gff3 (was the families stage run --from-gtf?)", key_str(rec));
            };
            let mut rec_ex: Vec<&(String, i64, i64)> = Vec::new();
            for g in names {
                match gene_exons.get(g) {
                    Some(v) if !v.is_empty() => rec_ex.extend(v.iter()),
                    _ => bail!("gene_id {g} of {} has no exons in the GTF", key_str(rec)),
                }
                gids.push(g.clone());
            }
            if names.len() > 1 {
                cnt.inc("records_with_key_collision");
            }
            if rec_ex.iter().any(|(c, _, _)| *c != rec.0) {
                bail!("{} has exons on another contig", key_str(rec));
            }
            let lo = rec_ex.iter().map(|x| x.1).min().unwrap();
            let hi = rec_ex.iter().map(|x| x.2).max().unwrap();
            if (lo, hi) != (rec.1, rec.2) {
                cnt.inc("records_span_mismatch"); // the key is not the gene's all-transcript span
            }
            exs.extend(rec_ex.iter().map(|x| (x.1 - 1, x.2)));
        }
        let blocks = merge_blocks(exs);
        loci.push(Locus {
            key: m.clone(),
            family: cl.fam_of[m],
            records: recs,
            gene_ids: gids,
            starts: blocks.iter().map(|b| b.0).collect(),
            ends: blocks.iter().map(|b| b.1).collect(),
        });
    }
    Ok((loci, locus_of_record))
}

/// hits[m][block][partner] = genome intervals of m's block joined to an exon base of partner (loci by index).
type Hits = Vec<HashMap<usize, HashMap<usize, Vec<(i64, i64)>>>>;

fn add_iv(lst: &mut Vec<(i64, i64)>, lo: i64, hi: i64) {
    // coalesce with the last interval when they touch (the union is what is ever read)
    if let Some(last) = lst.last_mut() {
        if lo <= last.1 && hi >= last.0 {
            *last = (lo.min(last.0), hi.max(last.1));
            return;
        }
    }
    lst.push((lo, hi));
}

fn project_paf(
    r: &mut dyn BufRead,
    path: &str,
    loci: &[Locus],
    locus_of_record: &HashMap<Key, usize>,
    cnt: &mut Counts,
) -> Result<Hits> {
    let mut hits: Hits = (0..loci.len()).map(|_| HashMap::new()).collect();
    // record name -> (its key's contig start, its locus), memoised on the raw name (names repeat on every record)
    let mut memo: HashMap<String, Option<(i64, usize)>> = HashMap::new();
    let mut lookup = |name: &str| -> Option<(i64, usize)> {
        if let Some(v) = memo.get(name) {
            return *v;
        }
        let v = parse_key(name).and_then(|k| locus_of_record.get(&k).map(|&l| (k.1, l)));
        memo.insert(name.to_string(), v);
        v
    };
    for line in r.lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            cnt.inc("paf_malformed");
            continue;
        }
        cnt.inc("paf_records");
        if f[0] == f[5] {
            cnt.inc("paf_self");
            continue;
        }
        let (Some((q_start, lq)), Some((t_start, lt))) = (lookup(f[0]), lookup(f[5])) else {
            cnt.inc("paf_unclustered");
            continue;
        };
        if lq == lt {
            cnt.inc("paf_same_locus");
            continue;
        }
        let (qs, qe) = (int_field(f[2], "PAF qs")?, int_field(f[3], "PAF qe")?);
        let strand = f[4];
        let (ts, te) = (int_field(f[7], "PAF ts")?, int_field(f[8], "PAF te")?);
        let (q_off, t_off) = (q_start - 1, t_start - 1);
        let (lqq, ltt) = (&loci[lq], &loci[lt]);
        if !(overlaps_any(&lqq.starts, &lqq.ends, q_off + qs, q_off + qe)
            && overlaps_any(&ltt.starts, &ltt.ends, t_off + ts, t_off + te))
        {
            cnt.inc("paf_no_exon_interval");
            continue;
        }
        let Some(cg) = f[12..].iter().find_map(|x| x.strip_prefix("cg:Z:")) else {
            bail!("{path}: record {} -> {} has no cg:Z CIGAR (the families PAF is run with -c)", f[0], f[5]);
        };
        cnt.inc("paf_projected");
        let cols = exon_columns(
            cg,
            strand,
            qs,
            qe,
            ts,
            te,
            q_off,
            t_off,
            (&lqq.starts, &lqq.ends),
            (&ltt.starts, &ltt.ends),
        )?;
        for &(ql, qh, qb, tl, th, tb) in &cols {
            add_iv(hits[lq].entry(qb).or_default().entry(lt).or_default(), ql, qh);
            add_iv(hits[lt].entry(tb).or_default().entry(lq).or_default(), tl, th);
        }
        if !cols.is_empty() {
            cnt.inc("paf_exon_exon");
        }
    }
    Ok(hits)
}

/// The container of one families run: the rows of the three tables (fields as written) and the counts.
pub struct Container {
    pub blocks: Vec<Vec<String>>,
    pub relations: Vec<Vec<String>>,
    counts: Counts,
}

fn classify(cl: &Clusters, loci: &[Locus], hits: &Hits, strand_of: &HashMap<Key, String>) -> Container {
    let mut cnt = Counts::default();
    let mut rows: Vec<Vec<String>> = Vec::new();
    // (family, related family) -> (members, n_blocks, bp); a family index IS its first-appearance order
    let mut rel: HashMap<(usize, usize), (Vec<usize>, i64, i64)> = HashMap::new();
    let locus_idx: HashMap<&Key, usize> = loci.iter().enumerate().map(|(i, l)| (&l.key, i)).collect();
    let empty: HashMap<usize, Vec<(i64, i64)>> = HashMap::new();
    for (fid, fam_members) in cl.members.iter().enumerate() {
        for mk in fam_members {
            let mi = locus_idx[mk];
            let l = &loci[mi];
            let nb = l.starts.len();
            let mut has_acc = false;
            for bi in 0..nb {
                let (s, e) = (l.starts[bi], l.ends[bi]);
                let partners = hits[mi].get(&bi).unwrap_or(&empty);
                let mut same: Vec<usize> = partners.keys().copied().filter(|&p| loci[p].family == fid).collect();
                same.sort_by(|&a, &b| loci[a].key.cmp(&loci[b].key));
                let mut row: Vec<String> = vec![
                    cl.fam_ids[fid].clone(),
                    key_str(&l.key),
                    l.key.0.clone(),
                    strand_of.get(&l.key).cloned().unwrap_or_else(|| ".".into()),
                    l.gene_ids.join(","),
                    l.records.len().to_string(),
                    (bi + 1).to_string(),
                    nb.to_string(),
                    (s + 1).to_string(),
                    e.to_string(),
                    (e - s).to_string(),
                ];
                cnt.inc("blocks");
                if !same.is_empty() {
                    let core_bp = union_len(same.iter().flat_map(|p| partners[p].iter()));
                    let core_partners: Vec<String> = same
                        .iter()
                        .map(|p| format!("{}={}", key_str(&loci[*p].key), union_len(partners[p].iter())))
                        .collect();
                    row.extend([
                        "core".to_string(),
                        core_bp.to_string(),
                        core_partners.join(","),
                        ".".into(),
                        ".".into(),
                        ".".into(),
                    ]);
                    cnt.inc("core_blocks");
                    cnt.add("core_block_bp", e - s);
                } else {
                    has_acc = true;
                    let mut other: Vec<usize> = partners.keys().copied().filter(|&p| loci[p].family != fid).collect();
                    other.sort_by(|&a, &b| (loci[a].family, &loci[a].key).cmp(&(loci[b].family, &loci[b].key)));
                    row.extend(["accessory".to_string(), "0".into(), ".".into()]);
                    cnt.inc("accessory_blocks");
                    cnt.add("accessory_bp", e - s);
                    if !other.is_empty() {
                        let mut fams: Vec<usize> = other.iter().map(|&p| loci[p].family).collect();
                        fams.dedup(); // `other` is sorted by family order
                        let rel_bp = union_len(other.iter().flat_map(|p| partners[p].iter()));
                        let rel_partners: Vec<String> = other
                            .iter()
                            .map(|&p| {
                                format!(
                                    "{}|{}={}",
                                    cl.fam_ids[loci[p].family],
                                    key_str(&loci[p].key),
                                    union_len(partners[&p].iter())
                                )
                            })
                            .collect();
                        row.extend([
                            fams.iter().map(|&f| cl.fam_ids[f].as_str()).collect::<Vec<_>>().join(","),
                            rel_bp.to_string(),
                            rel_partners.join(","),
                        ]);
                        cnt.inc("accessory_blocks_related");
                        for &f2 in &fams {
                            let r = rel.entry((fid, f2)).or_default();
                            if r.0.last() != Some(&mi) {
                                r.0.push(mi);
                            }
                            r.1 += 1;
                            r.2 += union_len(
                                other.iter().filter(|&&p| loci[p].family == f2).flat_map(|p| partners[p].iter()),
                            );
                        }
                    } else {
                        row.extend([".".to_string(), "0".into(), ".".into()]);
                    }
                }
                rows.push(row);
            }
            cnt.inc("loci");
            cnt.inc(if has_acc { "loci_with_accessory" } else { "loci_all_core" });
        }
    }
    let mut keys: Vec<(usize, usize)> = rel.keys().copied().collect();
    keys.sort();
    let mut rel_rows = Vec::with_capacity(keys.len());
    for (f1, f2) in keys {
        let r = &rel[&(f1, f2)];
        let reciprocal = rel.contains_key(&(f2, f1));
        rel_rows.push(vec![
            cl.fam_ids[f1].clone(),
            cl.fam_ids[f2].clone(),
            r.0.len().to_string(),
            r.1.to_string(),
            r.2.to_string(),
            if reciprocal { "yes" } else { "no" }.to_string(),
            r.0.iter().map(|&m| key_str(&loci[m].key)).collect::<Vec<_>>().join(","),
        ]);
        if reciprocal {
            cnt.inc("family_relations_reciprocal");
        }
    }
    cnt.add("family_relations_directed", rel_rows.len() as i64);
    cnt.add("families", cl.members.len() as i64);
    cnt.add("families_with_relation", rel.keys().map(|k| k.0).collect::<HashSet<_>>().len() as i64);
    Container { blocks: rows, relations: rel_rows, counts: cnt }
}

/// The whole container on open inputs: the assembled GTF, the families' `clusters.tsv`, `loci.gff3`, the fold table
/// `loci.tsv` (None = no folded records) and the all-vs-all PAF (with `cg:Z`). `paf_name` names the PAF in errors.
pub fn run(
    gtf: &mut dyn BufRead,
    clusters: &mut dyn BufRead,
    loci_gff3: &mut dyn BufRead,
    loci_tsv: Option<&mut dyn BufRead>,
    paf: &mut dyn BufRead,
    paf_name: &str,
) -> Result<Container> {
    let cl = read_clusters(clusters, "clusters.tsv")?;
    let folds = match loci_tsv {
        Some(r) => read_folds(r, "loci.tsv")?,
        None => Vec::new(),
    };
    let (loci_genes, strand_of) = read_loci_gff3(loci_gff3, "loci.gff3")?;
    let mut wanted: HashSet<String> = HashSet::new();
    for m in &cl.keys {
        wanted.extend(loci_genes.get(m).into_iter().flatten().cloned());
    }
    for (ann, rep) in &folds {
        if cl.fam_of.contains_key(rep) {
            wanted.extend(loci_genes.get(ann).into_iter().flatten().cloned());
        }
    }
    let gene_exons = read_gtf(gtf, &wanted)?;
    let mut cnt = Counts::default();
    let (loci, locus_of_record) = build_loci(&cl, &folds, &loci_genes, &gene_exons, &mut cnt)?;
    cnt.add("loci_gff3_keys", loci_genes.len() as i64);
    cnt.add("gene_key_collisions", loci_genes.values().filter(|v| v.len() > 1).count() as i64);
    let hits = project_paf(paf, paf_name, &loci, &locus_of_record, &mut cnt)?;
    let mut c = classify(&cl, &loci, &hits, &strand_of);
    for (k, v) in cnt.0 {
        c.counts.add(k, v);
    }
    Ok(c)
}

/// [`run`] on file paths (`.gz` read through gzip, as the script); `loci_tsv` None or a missing file = no folds.
pub fn run_paths(gtf: &str, clusters: &str, loci_gff3: &str, loci_tsv: Option<&str>, paf: &str) -> Result<Container> {
    let mut lt = match loci_tsv {
        Some(p) if std::path::Path::new(p).exists() => Some(open_text(p)?),
        _ => None,
    };
    run(
        &mut *open_text(gtf)?,
        &mut *open_text(clusters)?,
        &mut *open_text(loci_gff3)?,
        lt.as_deref_mut().map(|r| r as &mut dyn BufRead),
        &mut *open_text(paf)?,
        paf,
    )
}

fn tsv(header: &[&str], rows: &[Vec<String>]) -> String {
    let mut s = header.join("\t");
    s.push('\n');
    for r in rows {
        s.push_str(&r.join("\t"));
        s.push('\n');
    }
    s
}

impl Container {
    pub fn count(&self, k: &str) -> i64 {
        self.counts.get(k)
    }
    /// `<out>.container.tsv` (the script's `OUT.blocks.tsv`).
    pub fn blocks_tsv(&self) -> String {
        tsv(&BLOCK_HEADER, &self.blocks)
    }
    /// `<out>.container_relations.tsv` (the script's `OUT.relations.tsv`).
    pub fn relations_tsv(&self) -> String {
        tsv(&REL_HEADER, &self.relations)
    }
    /// `<out>.container_summary.tsv` (the script's `OUT.summary.tsv`).
    pub fn summary_tsv(&self) -> String {
        let mut s = String::from("key\tvalue\n");
        for k in SUMMARY_KEYS {
            s.push_str(&format!("{k}\t{}\n", self.count(k)));
        }
        s
    }
    /// The script's stderr summary line.
    pub fn summary_line(&self) -> String {
        let kv: Vec<String> = SUMMARY_KEYS.iter().map(|k| format!("{k}={}", self.count(k))).collect();
        format!("[family_container] {}", kv.join(" "))
    }
    /// Write the three tables next to the families products (`<out>.container.tsv`, `<out>.container_relations.tsv`,
    /// `<out>.container_summary.tsv`).
    pub fn write(&self, out: &str) -> Result<()> {
        for (suffix, text) in [
            ("container.tsv", self.blocks_tsv()),
            ("container_relations.tsv", self.relations_tsv()),
            ("container_summary.tsv", self.summary_tsv()),
        ] {
            let path = format!("{out}.{suffix}");
            std::fs::write(&path, text).with_context(|| format!("writing {path}"))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    //! The frozen script's 20 unit tests (`bench/test_family_container.py`, sha1 a2f4a13d), ported one to one on the
    //! same hand-made fixture, plus byte-identity of all three tables against the script's own outputs on that
    //! fixture (`testdata/family_container/`, written by the frozen script).
    //!
    //! The fixture (GFF 1-based closed; PAF offsets are 0-based within each record, record key START = offset 0):
    //!
    //!   family MCL0: A = c1:1001-2000  gA  tA1 1001-1100,1301-1400,1901-2000 ; tA2 1051-1150,1901-1950
    //!                  + folded record A2 = c1:1951-2100 (gA2, one exon 1951-2100)  => blocks 1001-1150, 1301-1400, 1901-2100
    //!                B = c1:5001-6000  gB  tB1 5001-5150, 5301-5400, 5701-6000
    //!   family MCL1: D = c2:1001-2000  gD  tD1 1001-1200, 1601-2000
    //!                E = c2:5001-6000  gE  tE1 5001-5200, 5501-5600, 5901-6000
    //!   unclustered: U = c3:1001-2000  gU  tU1 1001-1100, 1901-2000
    //!
    //!   R1  A->B  +  q[0,150)    t[0,150)   50M2I48M2D50M     A.b1/B.b1 exon-exon, 148 aligned columns each side
    //!   R2  A->B  +  q[399,900)  t[399,900) 501M              touches A.b2/B.b2 by exactly 1 base; ends 0 bases before A.b3
    //!   R3  A->D  -  q[700,1000) t[100,400) 100M5D95M5I100M   reverse strand: only the first 100M is exon-exon
    //!   R4  B->A  +  q[700,760)  t[450,510) 60M               B.b3 onto A's INTRON: not evidence
    //!   R5  B->U  +  q[800,900)  t[0,100)   100M              U unclustered: skipped
    //!   R6  A2->D +  q[100,150)  t[650,700) 50M               the folded record counts for A (A.b3 genome 2051-2100)
    //!   R7  A->A2 +  q[950,1000) t[0,50)    50M               same locus: skipped
    //!   R8  D->E  +  q[600,650)  t[0,50)    50M               D.b2 / E.b1 core
    //!   R9  A->U  +  q[900,1000) t[900,1000) 100M             U unclustered: skipped
    //!   R10 E->B  +  q[500,600)  t[0,100)   100M              E.b2 accessory related to MCL0 (B.b1 is core already)
    use super::*;
    use std::collections::BTreeMap;
    use std::io::Write;

    type Genes = Vec<(&'static str, &'static str, &'static str, Vec<(&'static str, Vec<(i64, i64)>)>)>;
    fn genes() -> Genes {
        vec![
            ("gA", "c1", "+", vec![("tA1", vec![(1001, 1100), (1301, 1400), (1901, 2000)]), ("tA2", vec![(1051, 1150), (1901, 1950)])]),
            ("gA2", "c1", "+", vec![("tA2x", vec![(1951, 2100)])]),
            ("gB", "c1", "+", vec![("tB1", vec![(5001, 5150), (5301, 5400), (5701, 6000)])]),
            ("gD", "c2", "-", vec![("tD1", vec![(1001, 1200), (1601, 2000)])]),
            ("gE", "c2", "-", vec![("tE1", vec![(5001, 5200), (5501, 5600), (5901, 6000)])]),
            ("gU", "c3", "+", vec![("tU1", vec![(1001, 1100), (1901, 2000)])]),
        ]
    }
    fn key(g: &str) -> &'static str {
        match g {
            "gA" => "c1:1001-2000",
            "gA2" => "c1:1951-2100",
            "gB" => "c1:5001-6000",
            "gD" => "c2:1001-2000",
            "gE" => "c2:5001-6000",
            "gU" => "c3:1001-2000",
            _ => panic!("{g}"),
        }
    }
    const CLUSTERS: [(&str, &str); 4] = [("MCL0", "gA"), ("MCL0", "gB"), ("MCL1", "gD"), ("MCL1", "gE")];
    type PafRow = (&'static str, i64, i64, &'static str, &'static str, i64, i64, &'static str);
    const PAF: [PafRow; 10] = [
        ("gA", 0, 150, "+", "gB", 0, 150, "50M2I48M2D50M"),
        ("gA", 399, 900, "+", "gB", 399, 900, "501M"),
        ("gA", 700, 1000, "-", "gD", 100, 400, "100M5D95M5I100M"),
        ("gB", 700, 760, "+", "gA", 450, 510, "60M"),
        ("gB", 800, 900, "+", "gU", 0, 100, "100M"),
        ("gA2", 100, 150, "+", "gD", 650, 700, "50M"),
        ("gA", 950, 1000, "+", "gA2", 0, 50, "50M"),
        ("gD", 600, 650, "+", "gE", 0, 50, "50M"),
        ("gA", 900, 1000, "+", "gU", 900, 1000, "100M"),
        ("gE", 500, 600, "+", "gB", 0, 100, "100M"),
    ];
    fn key_len(g: &str) -> i64 {
        let k = parse_key(key(g)).unwrap();
        k.2 - k.1 + 1
    }

    fn tmpdir(tag: &str) -> std::path::PathBuf {
        static N: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = N.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let d = std::env::temp_dir().join(format!("rustle_family_container_{tag}_{}_{n}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    /// `write_fixture` of the Python tests, byte for byte. Returns (gtf, fam prefix).
    fn write_fixture(d: &std::path::Path, clusters: &[(&str, &str)], folds: &[(&str, &str)], with_cg: bool) -> (String, String) {
        let gtf = d.join("x.gtf").display().to_string();
        let mut fh = std::fs::File::create(&gtf).unwrap();
        for (g, c, st, txs) in genes() {
            for (t, exs) in &txs {
                writeln!(fh, "{c}\trustle\ttranscript\t{}\t{}\t.\t{st}\t.\tgene_id \"{g}\"; transcript_id \"{t}\"; reads \"5\";", exs[0].0, exs[exs.len() - 1].1).unwrap();
                for (i, (s, e)) in exs.iter().enumerate() {
                    writeln!(fh, "{c}\trustle\texon\t{s}\t{e}\t.\t{st}\t.\tgene_id \"{g}\"; transcript_id \"{t}\"; exon_number \"{}\";", i + 1).unwrap();
                }
            }
        }
        let fam = d.join("x.fam").display().to_string();
        let mut fh = std::fs::File::create(format!("{fam}.loci.gff3")).unwrap();
        writeln!(fh, "##gff-version 3").unwrap();
        for (g, c, st, txs) in genes() {
            let k = parse_key(key(g)).unwrap();
            writeln!(fh, "{c}\t.\tgene\t{}\t{}\t.\t{st}\t.\tID=gene-{g};Name={g}", k.1, k.2).unwrap();
            for (s, e) in &txs[0].1 {
                writeln!(fh, "{c}\t.\texon\t{s}\t{e}\t.\t{st}\t.\tParent=gene-{g};gene={g}").unwrap();
            }
        }
        let mut fh = std::fs::File::create(format!("{fam}.clusters.tsv")).unwrap();
        writeln!(fh, "cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend").unwrap();
        for (fid, g) in clusters {
            let k = parse_key(key(g)).unwrap();
            let n = clusters.iter().filter(|(f, _)| f == fid).count();
            writeln!(fh, "{fid}\t{n}\t1.0000\t1.0000\tNA\t{}\t{}\t{}", k.0, k.1, k.2).unwrap();
        }
        let mut fh = std::fs::File::create(format!("{fam}.loci.tsv")).unwrap();
        writeln!(fh, "annotation\trepresentative").unwrap();
        for (a, r) in folds {
            writeln!(fh, "{}\t{}", key(a), key(r)).unwrap();
        }
        let mut fh = std::fs::File::create(format!("{fam}.loci.paf")).unwrap();
        for (qg, qs, qe, st, tg, ts, te, cg) in PAF {
            let tags = if with_cg { format!("\tNM:i:0\ttp:A:P\tcg:Z:{cg}") } else { "\tNM:i:0\ttp:A:P".to_string() };
            writeln!(
                fh,
                "{}\t{}\t{qs}\t{qe}\t{st}\t{}\t{}\t{ts}\t{te}\t{}\t{}\t60{tags}",
                key(qg),
                key_len(qg),
                key(tg),
                key_len(tg),
                qe - qs,
                (qe - qs).max(te - ts)
            )
            .unwrap();
        }
        (gtf, fam)
    }

    fn run_fam(gtf: &str, fam: &str) -> Result<Container> {
        run_paths(
            gtf,
            &format!("{fam}.clusters.tsv"),
            &format!("{fam}.loci.gff3"),
            Some(&format!("{fam}.loci.tsv")),
            &format!("{fam}.loci.paf"),
        )
    }

    /// The `Container` class fixture: the main fixture, run and written, read back as the Python tests read it.
    struct Out {
        header: Vec<String>,
        rows: BTreeMap<(String, i64), BTreeMap<String, String>>,
        rel: Vec<BTreeMap<String, String>>,
        summary: BTreeMap<String, i64>,
    }
    fn container_fixture() -> Out {
        let d = tmpdir("main");
        let (gtf, fam) = write_fixture(&d, &CLUSTERS, &[("gA2", "gA")], true);
        let out = d.join("out").display().to_string();
        run_fam(&gtf, &fam).unwrap().write(&out).unwrap();
        let read = |p: String| std::fs::read_to_string(p).unwrap();
        let blocks = read(format!("{out}.container.tsv"));
        let mut lines = blocks.lines();
        let header: Vec<String> = lines.next().unwrap().split('\t').map(String::from).collect();
        let mut rows = BTreeMap::new();
        for l in lines {
            let r: BTreeMap<String, String> = header.iter().cloned().zip(l.split('\t').map(String::from)).collect();
            rows.insert((r["locus"].clone(), r["block"].parse().unwrap()), r);
        }
        let rels = read(format!("{out}.container_relations.tsv"));
        let mut lines = rels.lines();
        let rh: Vec<String> = lines.next().unwrap().split('\t').map(String::from).collect();
        let rel = lines.map(|l| rh.iter().cloned().zip(l.split('\t').map(String::from)).collect()).collect();
        let summary = read(format!("{out}.container_summary.tsv"))
            .lines()
            .skip(1)
            .map(|l| {
                let (k, v) = l.split_once('\t').unwrap();
                (k.to_string(), v.parse().unwrap())
            })
            .collect();
        let _ = std::fs::remove_dir_all(&d);
        Out { header, rows, rel, summary }
    }
    impl Out {
        fn row(&self, g: &str, b: i64) -> &BTreeMap<String, String> {
            &self.rows[&(key(g).to_string(), b)]
        }
    }
    fn f<'a>(r: &'a BTreeMap<String, String>, k: &str) -> &'a str {
        r[k].as_str()
    }

    // ---------------------------------------------------------------------------------------- Primitives (6)

    #[test]
    fn merge_blocks_overlap_merges_abutting_stays_separate() {
        assert_eq!(merge_blocks(vec![(10, 20), (0, 10), (15, 30), (40, 50), (45, 46)]), vec![(0, 10), (10, 30), (40, 50)]);
    }

    #[test]
    fn union_len_counts_overlapping_and_abutting_once() {
        assert_eq!(union_len(&[(0, 10), (10, 20), (15, 25), (30, 31)]), 26);
        assert_eq!(union_len(&[]), 0);
    }

    #[test]
    fn aligned_runs_forward_with_eq_x_and_n() {
        let runs = aligned_runs("10=2X3N4I5M1D6M", "+", 100, 127, 50, 77).unwrap();
        assert_eq!(runs, vec![(50, 100, 10), (60, 110, 2), (65, 116, 5), (71, 121, 6)]);
    }

    #[test]
    fn aligned_runs_reverse() {
        // '-': query consumed from qe downwards; run query interval [qe - consumed - n, qe - consumed)
        let runs = aligned_runs("100M5D95M5I100M", "-", 700, 1000, 100, 400).unwrap();
        assert_eq!(runs, vec![(100, 900, 100), (205, 805, 95), (300, 700, 100)]);
    }

    #[test]
    fn cigar_length_mismatch_and_bad_ops_raise() {
        assert!(aligned_runs("10M", "+", 0, 11, 0, 10).is_err());
        assert!(aligned_runs("5S10M", "+", 0, 10, 0, 10).is_err());
        assert!(aligned_runs("10M", ".", 0, 10, 0, 10).is_err());
    }

    #[test]
    fn exon_columns_reverse_maps_each_column_exactly() {
        // query record offset 0 = genome 0; target offset 0 = genome 1000. Query exon [3,5), target exon [1000,1002):
        // '-' 10M on q[0,10) / t[0,10): column i joins t 1000+i with q 9-i, so t [1000,1002) <-> q [8,10) -- not an
        // exon on the query -- and q [3,5) <-> t [1005,1007) -- not an exon on the target. No exon-exon column.
        let q: (&[i64], &[i64]) = (&[3], &[5]);
        assert!(exon_columns("10M", "-", 0, 10, 0, 10, 0, 1000, q, (&[1000], &[1002])).unwrap().is_empty());
        // a target exon at [1005,1007) is exactly the mirror of the query exon: two columns
        let t: (&[i64], &[i64]) = (&[1005], &[1007]);
        assert_eq!(exon_columns("10M", "-", 0, 10, 0, 10, 0, 1000, q, t).unwrap(), vec![(3, 5, 0, 1005, 1007, 0)]);
        // and on '+' the same exons do not meet (q [3,5) <-> t [1003,1005))
        assert!(exon_columns("10M", "+", 0, 10, 0, 10, 0, 1000, q, t).unwrap().is_empty());
    }

    // ----------------------------------------------------------------------------------------- Container (11)

    #[test]
    fn container_header() {
        assert_eq!(container_fixture().header, BLOCK_HEADER.to_vec());
    }

    #[test]
    fn multi_transcript_and_folded_record_union() {
        let o = container_fixture();
        let a: Vec<_> = (1..=3).map(|b| o.row("gA", b)).collect();
        let spans: Vec<(i64, i64)> = a.iter().map(|r| (f(r, "start").parse().unwrap(), f(r, "end").parse().unwrap())).collect();
        assert_eq!(spans, vec![(1001, 1150), (1301, 1400), (1901, 2100)]);
        assert_eq!(f(a[0], "gene_ids"), "gA,gA2");
        assert_eq!(f(a[0], "n_records"), "2");
        assert_eq!(f(a[0], "n_blocks"), "3");
        assert!(!o.rows.contains_key(&(key("gA2").to_string(), 1)), "a folded record is part of its locus, not a row of its own");
    }

    #[test]
    fn forward_cigar_with_insertion_and_deletion() {
        let o = container_fixture();
        let (a1, b1) = (o.row("gA", 1), o.row("gB", 1));
        assert_eq!((f(a1, "class"), f(a1, "core_bp"), f(a1, "core_partners")), ("core", "148", format!("{}=148", key("gB")).as_str()));
        assert_eq!((f(b1, "class"), f(b1, "core_bp")), ("core", "148"));
        assert_eq!(f(a1, "rel_families"), ".");
    }

    #[test]
    fn one_base_touch_is_core_zero_is_accessory() {
        let o = container_fixture();
        assert_eq!((f(o.row("gA", 2), "class"), f(o.row("gA", 2), "core_bp")), ("core", "1"));
        assert_eq!((f(o.row("gB", 2), "class"), f(o.row("gB", 2), "core_bp")), ("core", "1"));
        let a3 = o.row("gA", 3); // R2 ends one base before it; its partner base there IS a B exon base
        assert_eq!(f(a3, "class"), "accessory");
        assert_eq!(f(a3, "core_partners"), ".");
    }

    #[test]
    fn reverse_strand_projection_and_relation() {
        let o = container_fixture();
        let a3 = o.row("gA", 3);
        assert_eq!(f(a3, "rel_families"), "MCL1");
        assert_eq!(f(a3, "rel_bp"), "150"); // R3's 100 (genome 1901-2000) + R6's 50 via the folded record
        assert_eq!(f(a3, "rel_partners"), format!("MCL1|{}=150", key("gD")));
        let d1 = o.row("gD", 1);
        assert_eq!((f(d1, "class"), f(d1, "rel_families"), f(d1, "rel_bp")), ("accessory", "MCL0", "100"));
        assert_eq!(f(d1, "rel_partners"), format!("MCL0|{}=100", key("gA")));
    }

    #[test]
    fn partner_intron_is_not_evidence_and_unclustered_is_ignored() {
        let o = container_fixture();
        let b3 = o.row("gB", 3);
        assert_eq!((f(b3, "class"), f(b3, "rel_families"), f(b3, "rel_bp")), ("accessory", ".", "0"));
    }

    #[test]
    fn core_blocks_are_not_tested_for_relations() {
        let o = container_fixture();
        let (d2, e1, b1) = (o.row("gD", 2), o.row("gE", 1), o.row("gB", 1));
        assert_eq!((f(d2, "class"), f(d2, "core_bp"), f(d2, "rel_families")), ("core", "50", "."));
        assert_eq!((f(e1, "class"), f(e1, "core_bp")), ("core", "50"));
        assert_eq!(f(b1, "rel_families"), "."); // R10 hits B.b1 from MCL1 but B.b1 is core
    }

    #[test]
    fn accessory_without_any_alignment() {
        let o = container_fixture();
        let e3 = o.row("gE", 3);
        assert_eq!((f(e3, "class"), f(e3, "rel_families"), f(e3, "bp")), ("accessory", ".", "100"));
        let e2 = o.row("gE", 2);
        assert_eq!((f(e2, "class"), f(e2, "rel_families"), f(e2, "rel_bp")), ("accessory", "MCL0", "100"));
    }

    #[test]
    fn unclustered_locus_has_no_rows() {
        let o = container_fixture();
        assert!(!o.rows.keys().any(|k| k.0 == key("gU")));
        assert_eq!(o.rows.len(), 3 + 3 + 2 + 3);
    }

    #[test]
    fn family_relation_table() {
        let o = container_fixture();
        let got: Vec<Vec<&str>> = o
            .rel
            .iter()
            .map(|r| ["family_id", "related_family", "n_members", "n_blocks", "bp", "reciprocal", "members"].iter().map(|k| f(r, k)).collect())
            .collect();
        let both = format!("{},{}", key("gD"), key("gE"));
        assert_eq!(
            got,
            vec![vec!["MCL0", "MCL1", "1", "1", "150", "yes", key("gA")], vec!["MCL1", "MCL0", "2", "2", "200", "yes", both.as_str()]]
        );
    }

    #[test]
    fn summary_counts() {
        let s = container_fixture().summary;
        let g = |k: &str| s[k];
        assert_eq!((g("families"), g("loci"), g("folded_records"), g("blocks")), (2, 4, 1, 11));
        assert_eq!((g("core_blocks"), g("accessory_blocks"), g("accessory_blocks_related")), (6, 5, 3));
        assert_eq!((g("paf_records"), g("paf_unclustered"), g("paf_same_locus"), g("paf_no_exon_interval")), (10, 2, 1, 1));
        assert_eq!((g("paf_projected"), g("paf_exon_exon")), (6, 6));
        assert_eq!((g("loci_all_core"), g("loci_with_accessory")), (0, 4));
        assert_eq!(g("records_span_mismatch"), 0);
    }

    // --------------------------------------------------------------------------------------------- Guards (3)

    #[test]
    fn record_without_cigar_raises() {
        let d = tmpdir("nocg");
        let (gtf, fam) = write_fixture(&d, &CLUSTERS, &[("gA2", "gA")], false);
        assert!(run_fam(&gtf, &fam).is_err());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn member_in_two_families_raises() {
        let d = tmpdir("twofam");
        let mut cl = CLUSTERS.to_vec();
        cl.push(("MCL1", "gA"));
        let (gtf, fam) = write_fixture(&d, &cl, &[("gA2", "gA")], true);
        assert!(run_fam(&gtf, &fam).is_err());
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn without_the_fold_the_folded_record_is_unclustered() {
        let d = tmpdir("nofold");
        let (gtf, fam) = write_fixture(&d, &CLUSTERS, &[], true);
        let c = run_fam(&gtf, &fam).unwrap();
        let col = |name: &str| BLOCK_HEADER.iter().position(|h| *h == name).unwrap();
        let a: Vec<&Vec<String>> = c.blocks.iter().filter(|r| r[col("locus")] == key("gA")).collect();
        let spans: Vec<(&str, &str)> = a.iter().map(|r| (r[col("start")].as_str(), r[col("end")].as_str())).collect();
        assert_eq!(spans, vec![("1001", "1150"), ("1301", "1400"), ("1901", "2000")]);
        assert_eq!(a[2][col("rel_bp")], "100"); // R6 (from the no-longer-folded record) no longer counts
        assert_eq!(c.count("paf_same_locus"), 0);
        assert_eq!(c.count("paf_unclustered"), 4); // R5, R6, R7, R9
        let _ = std::fs::remove_dir_all(&d);
    }

    // ---------------------------------------------------------------- byte identity with the frozen script

    /// All three tables, byte for byte, against the frozen script's outputs on the same fixture (main and no-fold).
    #[test]
    fn tables_are_byte_identical_to_the_frozen_script() {
        for (tag, folds, blocks, rels, summary) in [
            (
                "main",
                &[("gA2", "gA")][..],
                include_str!("testdata/family_container/main.blocks.tsv"),
                include_str!("testdata/family_container/main.relations.tsv"),
                include_str!("testdata/family_container/main.summary.tsv"),
            ),
            (
                "nofold",
                &[][..],
                include_str!("testdata/family_container/nofold.blocks.tsv"),
                include_str!("testdata/family_container/nofold.relations.tsv"),
                include_str!("testdata/family_container/nofold.summary.tsv"),
            ),
        ] {
            let d = tmpdir(tag);
            let (gtf, fam) = write_fixture(&d, &CLUSTERS, folds, true);
            let c = run_fam(&gtf, &fam).unwrap();
            assert_eq!(c.blocks_tsv(), blocks, "{tag} blocks");
            assert_eq!(c.relations_tsv(), rels, "{tag} relations");
            assert_eq!(c.summary_tsv(), summary, "{tag} summary");
            let _ = std::fs::remove_dir_all(&d);
        }
    }

    #[test]
    fn python_int_and_key_parsing() {
        assert_eq!(py_int(" 12 "), Some(12));
        assert_eq!(py_int("+7"), Some(7));
        assert_eq!(py_int("-3"), Some(-3));
        assert_eq!(py_int("1_000"), Some(1000));
        for bad in ["", "1__0", "_1", "1_", "1.0", "0x1", "+-1", "1e3"] {
            assert_eq!(py_int(bad), None, "{bad:?}");
        }
        assert_eq!(parse_key("chrUn:KI270:1-20"), Some(("chrUn:KI270".to_string(), 1, 20)));
        assert_eq!(parse_key("c1:0100-200"), Some(("c1".to_string(), 100, 200)));
        assert_eq!(parse_key("c1:1-2-3"), None);
        assert_eq!(parse_key("c1"), None);
        assert_eq!(gtf_attr("gene_id \"\"; transcript_id \"t\";", "gene_id"), Some(""));
    }
}
